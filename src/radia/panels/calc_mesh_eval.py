"""
Mesh evaluation: export format quality assurance + p-convergence.

Called as subprocess from Cubit Mesh Evaluation:
    python calc_mesh_eval.py --cub5 temp.cub5 --max-order 5

Phase 1: Open cub5, get CAD values (volume, area, length)
Phase 2: Format QA (order 1-2): export .msh/.bdf/.vtk, verify via GMSH API
         getJacobians vs CAD. This is OUR quality guarantee.
Phase 3: p-convergence (order 1-5): export .vol, read with NGSolve Integrate
         vs CAD. Confirms curving accuracy across orders.

IMPORTANT: NGSolve must be imported BEFORE cubit.
Outputs JSON to stdout (last line).
"""

import argparse
import json
import os
import sys
import tempfile
import time

_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import setup_paths, setup_cubit, progress, calc_main


def _log(msg):
    progress("EVAL", msg)


def evaluate_mesh(cub5_file, max_order=5):
    # Import NGSolve FIRST (before GMSH) to avoid MKL DLL conflicts
    import tempfile
    from ngsolve import Mesh as NGMesh, Integrate, CF, BND, BBND

    setup_paths()
    cubit = setup_cubit(cub5_file)
    if cubit is None:
        return {"error": "Cubit not available"}

    vol_ids = list(cubit.get_entities("volume"))
    if not vol_ids:
        return {"error": "No volumes found in model"}

    # === Phase 1: CAD reference values ===
    cad_mats = {}
    for vid in vol_ids:
        name = cubit.get_entity_name("volume", vid) or f"Volume {vid}"
        cad_mats[name] = cubit.volume(vid).volume()
    cad_bnds = {}
    for sid in cubit.parse_cubit_list("surface", "all"):
        cad_bnds[f"surface_{sid}"] = cubit.surface(sid).area()
    cad_edges = {}
    for cid in cubit.parse_cubit_list("curve", "all"):
        ps = list(cubit.parse_cubit_list("surface", f"in curve {cid}"))
        if len(ps) >= 2:
            cad_edges[f"curve_{cid}"] = cubit.curve(cid).length()

    cad_vol_total = sum(cad_mats.values())
    cad_area_total = sum(cad_bnds.values())
    cad_length_total = sum(cad_edges.values())

    _log(f"CAD: V={cad_vol_total:.6e}, A={cad_area_total:.6e}, "
         f"L={cad_length_total:.6e}")

    # Check meshed
    total_elems = sum(
        len(cubit.get_volume_tets(vid)) + len(cubit.get_volume_hexes(vid))
        for vid in vol_ids)
    if total_elems == 0:
        return {"error": "Volumes are not meshed",
                "cad_vol_total": cad_vol_total,
                "cad_area_total": cad_area_total,
                "cad_length_total": cad_length_total}

    tmpdir = tempfile.mkdtemp(prefix="radia_eval_")

    # === Phase 2: Format QA via GMSH API (order 1-2, vs CAD) ===
    # This is our quality guarantee: .msh/.bdf/.vtk export correctness
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
            ("gmsh",    "msh", 'radia_export gmsh "{path}" overwrite'),
            ("nastran", "bdf", 'radia_export nastran "{path}" overwrite'),
            ("vtk",     "vtk", 'radia_export vtk "{path}" overwrite'),
        ]
        for fmt_name, ext, cmd_template in formats:
            for order in [1, 2]:
                fname = f"qa_{fmt_name}_o{order}.{ext}"
                fpath = os.path.join(tmpdir, fname).replace("\\", "/")
                cmd_str = cmd_template.format(path=fpath)
                if order == 2:
                    cmd_str = cmd_str.replace("overwrite",
                                              f"order {order} overwrite")
                try:
                    cubit.cmd(cmd_str)
                    if not os.path.exists(fpath):
                        format_results.append({
                            "format": fmt_name, "order": order,
                            "error": "export failed"})
                        continue

                    gmsh.clear()
                    gmsh.open(fpath)

                    # Volume (dim=3) vs CAD
                    total_vol, neg_det = _gmsh_integrate(3)
                    vol_err = ((total_vol - cad_vol_total) / cad_vol_total
                               * 100 if cad_vol_total else 0)

                    # Area (dim=2) vs CAD
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
        _log(f"p={p}: exporting .vol...")
        t0 = time.perf_counter()

        try:
            vol_path = tempfile.mktemp(suffix='.vol')
            cubit.cmd(f'radia_export netgen "{vol_path}" order {p} overwrite')
            mesh = NGMesh(vol_path)
        except Exception as e:
            orders.append({"order": p, "error": str(e)})
            continue

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

        _log(f"p={p}: V_err={vol_err:+.4f}%, A_err={area_err:+.4f}%, "
             f"L_err={len_err:+.4f}%, t={t_elapsed:.2f}s")

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

    # Cleanup
    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)

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
        description="Mesh evaluation: format QA + p-convergence")
    parser.add_argument("--cub5", required=True,
                        help="Cubit .cub5 model file")
    parser.add_argument("--max-order", type=int, default=5)

    def run(args):
        return evaluate_mesh(args.cub5, args.max_order)

    calc_main(run, parser)


if __name__ == "__main__":
    main()
