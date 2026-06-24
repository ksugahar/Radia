"""
Mesh evaluation: format quality assurance + p-convergence.

Called as subprocess from the Cubit Export Mesh PySide6 menu:
    python calc_mesh_eval.py --vol-base /path/to/model --max-order 5

The Cubit-side export command has ALREADY exported:
  - {base}_p1.vol ... {base}_p5.vol  (Netgen .vol at each order)
  - {base}_p1.vol.json               (CAD reference values)
  - {base}_qa_gmsh_o1.msh ...        (format QA files)
  - {base}_qa_nastran_o1.bdf ...
  - {base}_qa_vtk_o1.vtk ...

This script does NOT import cubit.  It only uses NGSolve and GMSH.

Phase 1: Read CAD values from companion .vol.json
Phase 2: Format QA via GMSH API getJacobians (vs CAD)
Phase 3: p-convergence via NGSolve Integrate (vs CAD)

Outputs JSON to stdout (last line).
"""

import argparse
import json
import os
import sys
import time

_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import setup_paths, progress, calc_main


def _log(msg):
    progress("EVAL", msg)


def evaluate_mesh(vol_base, max_order=5):
    """Evaluate mesh quality from pre-exported .vol and format QA files.

    Args:
        vol_base: Base path without extension (e.g. /path/to/model).
                  Expects {base}_p1.vol ... {base}_p5.vol etc.
        max_order: Maximum polynomial order (1-5).
    """
    setup_paths()
    from ngsolve import Mesh as NGMesh, Integrate, CF, BND, BBND
    from ngsolve import TaskManager

    # === Phase 1: CAD reference values from companion JSON ===
    # Use p=1 companion JSON (all orders share same CAD geometry)
    json_path = vol_base + "_p1.vol.json"
    if not os.path.isfile(json_path):
        return {"error": f"Companion JSON not found: {json_path}"}

    with open(json_path, "r", encoding="utf-8") as f:
        cad = json.load(f)

    cad_mats = cad.get("materials", {})
    cad_bnds = cad.get("boundaries", {})
    cad_edges = cad.get("edges", {})

    cad_vol_total = sum(cad_mats.values())
    cad_area_total = sum(cad_bnds.values())
    cad_length_total = sum(cad_edges.values())

    _log(f"CAD: V={cad_vol_total:.6e}, A={cad_area_total:.6e}, "
         f"L={cad_length_total:.6e}")

    # === Phase 2: Format QA via GMSH API (vs CAD) ===
    format_results = []
    try:
        import gmsh
        gmsh.initialize()

        def _gmsh_integrate(dim, gauss_order="Gauss5"):
            """Compute integral(1) over elements via GMSH getJacobians."""
            etypes, etags, ntags = gmsh.model.mesh.getElements(dim=dim)
            total = 0.0
            neg = 0
            for et in etypes:
                lc, w = gmsh.model.mesh.getIntegrationPoints(
                    int(et), gauss_order)
                jac, det, pts = gmsh.model.mesh.getJacobians(int(et), lc)
                npts = len(w)
                nel = len(det) // npts
                for i in range(nel):
                    for j in range(npts):
                        d = det[i * npts + j]
                        total += d * w[j]
                        if d < 0:
                            neg += 1
            return total, neg

        formats = [
            # (name, ext, max_order_for_format)
            # GMSH order 4-5 is an error (unreliable HO interior nodes)
            ("gmsh",    "msh", min(max_order, 3)),
            ("nastran", "bdf", min(max_order, 2)),
            ("vtk",     "vtk", min(max_order, 2)),
        ]
        for fmt_name, ext, fmt_max_order in formats:
            for order in range(1, fmt_max_order + 1):
                fpath = f"{vol_base}_qa_{fmt_name}_o{order}.{ext}"
                if not os.path.isfile(fpath):
                    format_results.append({
                        "format": fmt_name, "order": order,
                        "error": "file not found"})
                    continue
                try:
                    gmsh.clear()
                    gmsh.open(fpath)

                    total_vol, neg_det = _gmsh_integrate(3)
                    vol_err = ((total_vol - cad_vol_total) / cad_vol_total
                               * 100 if cad_vol_total else 0)

                    total_area, _ = _gmsh_integrate(2)
                    area_err = ((total_area - cad_area_total) / cad_area_total
                                * 100 if cad_area_total else 0)

                    format_results.append({
                        "format": fmt_name, "order": order,
                        "volume": total_vol, "vol_error_pct": vol_err,
                        "area": total_area, "area_error_pct": area_err,
                        "neg_det": neg_det,
                    })
                    _log(f"{fmt_name} o{order}: V={total_vol:.6e} "
                         f"({vol_err:+.2e}% vs CAD), "
                         f"A={total_area:.6e} ({area_err:+.2e}% vs CAD), "
                         f"neg_det={neg_det}")
                except Exception as e:
                    format_results.append({
                        "format": fmt_name, "order": order,
                        "error": str(e)})

        gmsh.finalize()
    except ImportError:
        _log("GMSH not available, skipping format QA")
    except Exception as e:
        _log(f"GMSH format QA error: {e}")

    # === Phase 3: p-convergence via NGSolve (order 1-5, vs CAD) ===
    orders = []
    for p in range(1, max_order + 1):
        vol_path = f"{vol_base}_p{p}.vol"
        if not os.path.isfile(vol_path):
            orders.append({"order": p, "error": "file not found"})
            continue

        _log(f"p={p}: loading {os.path.basename(vol_path)}...")
        t0 = time.perf_counter()

        try:
            mesh = NGMesh(vol_path)
        except Exception as e:
            orders.append({"order": p, "error": str(e)})
            continue

        with TaskManager():
            ng_vol = Integrate(CF(1), mesh)
            ng_area = Integrate(CF(1), mesh, VOL_or_BND=BND)
            try:
                ng_length = Integrate(CF(1), mesh, BBND)
            except Exception:
                ng_length = 0.0

            t_elapsed = time.perf_counter() - t0

            vol_err = ((ng_vol - cad_vol_total) / cad_vol_total * 100
                       if cad_vol_total else 0)
            area_err = ((ng_area - cad_area_total) / cad_area_total * 100
                        if cad_area_total else 0)
            len_err = ((ng_length - cad_length_total) / cad_length_total * 100
                       if cad_length_total else 0)

            _log(f"p={p}: V_err={vol_err:+.2e}%, A_err={area_err:+.2e}%, "
                 f"L_err={len_err:+.2e}%, t={t_elapsed:.2f}s")

            orders.append({
                "order": p,
                "ng_volume": ng_vol,
                "ng_area": ng_area,
                "ng_length": ng_length,
                "vol_error_pct": vol_err,
                "area_error_pct": area_err,
                "len_error_pct": len_err,
                "time": t_elapsed,
            })

    return {
        "cad_vol_total": cad_vol_total,
        "cad_area_total": cad_area_total,
        "cad_length_total": cad_length_total,
        "format_qa": format_results,
        "orders": orders,
        "max_order": max_order,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Mesh evaluation: format QA + p-convergence (no Cubit)")
    parser.add_argument("--vol-base", required=True,
                        help="Base path (expects {base}_p1.vol etc.)")
    parser.add_argument("--max-order", type=int, default=5)

    def run(args):
        return evaluate_mesh(args.vol_base, args.max_order)

    calc_main(run, parser)


if __name__ == "__main__":
    main()
