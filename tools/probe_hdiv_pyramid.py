#!/usr/bin/env python
"""Detect whether the INSTALLED NGSolve implements H(div) on PYRAMID elements.

This is the gate for HDiv-VIM tet/hex COUPLING (mixed meshes): a conforming tet+hex mesh needs pyramid
transition elements (4 tri faces + 1 quad face) at every hex(quad-face) <-> tet(tri-face) interface, and
those pyramids must carry an H(div) flux.  Joachim (Schoberl, NGSolve) committed to add HDiv-pyramid on the
NGSolve side; we WAIT for it (do NOT reimplement -- same "complement NGSolve" pattern as wait-for-Hlib).

As of NGSolve 6.2.2606 the state is ALLOC-BUT-UNIMPLEMENTED: HDiv(pyramid_mesh, order=1) constructs, but the
first Assemble raises `HDivHighOrderFESpace: Pyramid elements not implemented yet!`.  So the check must be
FUNCTIONAL (assemble a mass form + reproduce a constant field), not just "did HDiv() not raise".

Verdicts (also the exit code, so a skill / CI can branch):
  IMPLEMENTED       (exit 0)  -- HDiv-pyramid assembles, mass is nonzero+finite, a constant field is
                                 reproduced (L2 err < tol).  The block is LIFTED: add a Radia pyramid
                                 charge-Gram mode (mirror the wedge port) + enable mixed meshes.
  NOT_IMPLEMENTED   (exit 10) -- HDiv-pyramid raises "...not implemented..." (the expected current state).
  ALLOC_BUT_BROKEN  (exit 11) -- HDiv-pyramid assembles but is degenerate (ndof=0 / zero-or-NaN mass /
                                 constant field NOT reproduced): partially wired, still unusable.
  ERROR             (exit 20) -- the probe itself failed (mesh build / import / unexpected exception).

Usage:  python tools/probe_hdiv_pyramid.py [--json]
"""
import argparse
import json
import sys
import traceback

TOL_CONST_REPRO = 1e-6      # a working H(div) contains constants -> reproduces (1,0,0) to ~machine zero


def _build_pyramid_mesh():
    """One square-base pyramid (base quad z=0 + apex at (0.5,0.5,1)); volume = 1/3.  Raises on API drift."""
    import ngsolve as ng
    from netgen.meshing import Mesh as NgMesh, MeshPoint, Element3D, Element2D, FaceDescriptor
    from netgen.csg import Pnt
    ngm = NgMesh(dim=3)
    ngm.SetMaterial(1, "pyr")
    fd = ngm.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    pts = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0.5, 0.5, 1)]
    pid = [ngm.Add(MeshPoint(Pnt(*p))) for p in pts]
    b, a = pid[:4], pid[4]
    ngm.Add(Element3D(1, [b[0], b[1], b[2], b[3], a]))      # PYRAMID = base quad (4) + apex (1)
    ngm.Add(Element2D(fd, [b[0], b[1], b[2], b[3]]))        # quad base
    for i in range(4):                                     # 4 tri sides
        ngm.Add(Element2D(fd, [b[i], b[(i + 1) % 4], a]))
    return ng.Mesh(ngm)


def probe():
    """Return a verdict dict (never raises for the HDiv-pyramid outcome; only mesh/import failure -> ERROR)."""
    import numpy as np
    result = {"verdict": None, "ngsolve_version": None, "detail": {}}
    try:
        import ngsolve as ng
        result["ngsolve_version"] = ng.__version__
        mesh = _build_pyramid_mesh()
        nv = sorted({len(e.vertices) for e in mesh.Elements(ng.VOL)})
        vol = float(ng.Integrate(ng.CoefficientFunction(1.0), mesh))
        result["detail"]["nv_per_elem"] = nv
        result["detail"]["mesh_volume"] = vol
        # sanity: the pyramid MESH itself must be valid (H1 works on it) -- isolates the HDiv question
        with ng.TaskManager():
            h1 = ng.H1(mesh, order=1)
            gfh = ng.GridFunction(h1); gfh.Set(ng.CoefficientFunction(1.0))
            h1_ok = abs(float(ng.Integrate((gfh - 1.0) ** 2, mesh))) < 1e-12
        result["detail"]["mesh_valid_h1"] = bool(h1_ok and nv == [5] and abs(vol - 1.0 / 3.0) < 1e-6)
        if not result["detail"]["mesh_valid_h1"]:
            result["verdict"] = "ERROR"
            result["detail"]["reason"] = "pyramid mesh did not validate (nv/volume/H1) -- probe is unreliable"
            return result
    except Exception as e:
        result["verdict"] = "ERROR"
        result["detail"]["reason"] = f"mesh build / import failed: {type(e).__name__}: {e}"
        result["detail"]["traceback"] = traceback.format_exc()
        return result

    # ---- the actual HDiv-pyramid functional probe ----
    try:
        with ng.TaskManager():
            fes = ng.HDiv(mesh, order=1)
            ndof = fes.ndof
            u, v = fes.TnT()
            a = ng.BilinearForm(fes); a += u * v * ng.dx
            a.Assemble()                                   # <-- raises "Pyramid ... not implemented" today
            rows, cols, vals = a.mat.COO()
            vv = np.array(vals)
            mass_fro = float(np.sqrt(np.sum(vv * vv))) if len(vv) else 0.0
            finite = bool(np.all(np.isfinite(vv))) if len(vv) else True
            gf = ng.GridFunction(fes); gf.Set(ng.CoefficientFunction((1, 0, 0)))
            const_err = abs(float(ng.Integrate((gf[0] - 1.0) ** 2 + gf[1] ** 2 + gf[2] ** 2, mesh)))
        result["detail"].update(hdiv_ndof=ndof, mass_frobenius=mass_fro,
                                 mass_finite=finite, const_field_L2err=const_err)
        functional = (ndof > 0 and mass_fro > 0.0 and finite and const_err < TOL_CONST_REPRO)
        result["verdict"] = "IMPLEMENTED" if functional else "ALLOC_BUT_BROKEN"
        if not functional:
            result["detail"]["reason"] = ("HDiv-pyramid assembled but is degenerate "
                                          f"(ndof={ndof}, massFro={mass_fro:.2e}, finite={finite}, "
                                          f"const_reproL2={const_err:.2e})")
    except Exception as e:
        msg = str(e)
        result["detail"]["exception"] = f"{type(e).__name__}: {msg[:200]}"
        low = msg.lower()
        if ("pyramid" in low and ("not implemented" in low or "implemented yet" in low)) \
                or "hdivhighorderfespace: pyramid" in low:
            result["verdict"] = "NOT_IMPLEMENTED"
        else:
            result["verdict"] = "ERROR"
            result["detail"]["reason"] = "HDiv-pyramid raised an UNEXPECTED exception (not the known " \
                                         "'not implemented' message) -- inspect before trusting the verdict"
            result["detail"]["traceback"] = traceback.format_exc()
    return result


_EXIT = {"IMPLEMENTED": 0, "NOT_IMPLEMENTED": 10, "ALLOC_BUT_BROKEN": 11, "ERROR": 20}


def main():
    ap = argparse.ArgumentParser(description="Detect NGSolve H(div)-pyramid support (HDiv-VIM coupling gate).")
    ap.add_argument("--json", action="store_true", help="emit the full verdict dict as JSON")
    args = ap.parse_args()
    r = probe()
    v = r["verdict"]
    if args.json:
        print(json.dumps(r, indent=2))
    else:
        print(f"NGSolve {r['ngsolve_version']}: HDiv-pyramid = {v}")
        d = r["detail"]
        if v == "NOT_IMPLEMENTED":
            print("  -> still BLOCKED: HDiv-VIM tet/hex coupling (mixed meshes) waits for NGSolve.")
            print(f"     ({d.get('exception', 'HDiv-pyramid Assemble raised the not-implemented guard')})")
        elif v == "IMPLEMENTED":
            print("  -> UNBLOCKED! Add a Radia pyramid charge-Gram mode (mirror the wedge port) + enable")
            print("     mixed meshes.  First probe the pyramid div-image (L2 order) + face types.")
            print(f"     (ndof={d.get('hdiv_ndof')}, massFro={d.get('mass_frobenius'):.3e}, "
                  f"const-repro L2err={d.get('const_field_L2err'):.2e})")
        elif v == "ALLOC_BUT_BROKEN":
            print(f"  -> PARTIAL/unusable: {d.get('reason')}")
        else:
            print(f"  -> ERROR: {d.get('reason')}")
    sys.exit(_EXIT.get(v, 20))


if __name__ == "__main__":
    main()
