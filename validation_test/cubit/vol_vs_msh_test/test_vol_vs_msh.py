"""
Compare .vol (NGSolve) vs .msh (GMSH API) for volume, area, and length.

Demonstrates the difference between three integration paths:
  1. .vol + NGSolve Integrate (curvedElements -> exact)
  2. .msh + GMSH API getJacobians (isoparametric -> approximate for coarse HO)
  3. .msh + NGSolve ReadGmsh + Integrate (HO nodes ignored -> linear)

Prerequisites:
  Run export_sphere.jou in Cubit first (via cubit-run skill or manually).
  Output files: sphere_o1.vol, sphere_o2.vol, sphere_o1.msh, sphere_o2.msh, ...

Usage:
  python validation_test/cubit/vol_vs_msh_test/test_vol_vs_msh.py
"""

import sys
import os
import math
import numpy as np

_test_dir = os.path.dirname(os.path.abspath(__file__))

RADIUS = 0.05
V_EXACT = (4.0 / 3.0) * math.pi * RADIUS**3
A_EXACT = 4.0 * math.pi * RADIUS**2


# ================================================================
# NGSolve: read .vol and compute volume + area
# ================================================================
def ngsolve_vol_integrate(vol_path):
    """Read .vol with NGSolve, return (volume, area, curve_order, nv, ne)."""
    from ngsolve import Mesh, Integrate, CF, BND
    mesh = Mesh(vol_path)
    with TaskManager():
        vol = Integrate(CF(1), mesh)
        area = Integrate(CF(1), mesh, BND)
        return {
            'volume': vol,
            'area': area,
            'curve_order': mesh.GetCurveOrder(),
            'nv': mesh.nv,
            'ne': mesh.ne,
        }


# ================================================================
# NGSolve: read .msh via ReadGmsh and compute volume + area
# ================================================================
def ngsolve_readgmsh_integrate(msh_path):
    """Read .msh with ReadGmsh, return (volume, area, nv, ne)."""
    from netgen.read_gmsh import ReadGmsh
    from ngsolve import Mesh, Integrate, CF, BND
    from ngsolve import TaskManager
    ng = ReadGmsh(msh_path)
    mesh = Mesh(ng)
    vol = Integrate(CF(1), mesh)
    area = Integrate(CF(1), mesh, BND)
    return {
        'volume': vol,
        'area': area,
        'curve_order': mesh.GetCurveOrder(),
        'nv': mesh.nv,
        'ne': mesh.ne,
    }


# ================================================================
# GMSH API: read .msh and compute volume + area via getJacobians
# ================================================================
def gmsh_api_integrate(msh_path):
    """Read .msh with GMSH API, compute volume (3D) and area (2D)."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(msh_path)

    result = {'volume': 0.0, 'area': 0.0,
              'n_neg_det_3d': 0, 'n_neg_det_2d': 0,
              'n_nodes': 0, 'n_elems_3d': 0, 'n_elems_2d': 0,
              'elem_info': ''}

    node_tags, _, _ = gmsh.model.mesh.getNodes()
    result['n_nodes'] = len(node_tags)

    # Volume from 3D elements
    etypes3, etags3, _ = gmsh.model.mesh.getElements(dim=3)
    infos = []
    for idx, et in enumerate(etypes3):
        name, dim, order, nn, _, _ = gmsh.model.mesh.getElementProperties(int(et))
        n_elem = len(etags3[idx])
        result['n_elems_3d'] += n_elem

        gauss_order = max(2 * order, 1)
        lc, w = gmsh.model.mesh.getIntegrationPoints(int(et), f"Gauss{gauss_order}")
        n_gp = len(w)
        _, det, _ = gmsh.model.mesh.getJacobians(int(et), lc)
        det = np.array(det)
        neg = int((det < 0).sum())
        result['n_neg_det_3d'] += neg

        for i in range(n_elem):
            for j in range(n_gp):
                result['volume'] += det[i * n_gp + j] * w[j]

        infos.append(f"{name}={n_elem}")

    # Area from 2D elements
    etypes2, etags2, _ = gmsh.model.mesh.getElements(dim=2)
    for idx, et in enumerate(etypes2):
        name, dim, order, nn, _, _ = gmsh.model.mesh.getElementProperties(int(et))
        n_elem = len(etags2[idx])
        result['n_elems_2d'] += n_elem

        gauss_order = max(2 * order, 1)
        lc, w = gmsh.model.mesh.getIntegrationPoints(int(et), f"Gauss{gauss_order}")
        n_gp = len(w)
        _, det, _ = gmsh.model.mesh.getJacobians(int(et), lc)
        det = np.array(det)
        neg = int((det < 0).sum())
        result['n_neg_det_2d'] += neg

        for i in range(n_elem):
            for j in range(n_gp):
                result['area'] += det[i * n_gp + j] * w[j]

        infos.append(f"{name}={n_elem}")

    result['elem_info'] = ", ".join(infos)
    gmsh.finalize()
    return result


# ================================================================
# Main comparison
# ================================================================
def main():
    print("=" * 78)
    print("  .vol vs .msh: Volume / Area / Length Comparison")
    print(f"  Sphere r={RADIUS} m, V_exact={V_EXACT:.6e}, A_exact={A_EXACT:.6e}")
    print("=" * 78)

    orders = []
    for o in [1, 2, 3]:
        vol_path = os.path.join(_test_dir, f"sphere_o{o}.vol")
        msh_path = os.path.join(_test_dir, f"sphere_o{o}.msh")
        if os.path.isfile(vol_path) and os.path.isfile(msh_path):
            orders.append(o)
    if not orders:
        print("\nERROR: No sphere_o*.vol / sphere_o*.msh found.")
        print("Run export_sphere.jou in Cubit first:")
        print(f"  cd {_test_dir}")
        print("  cubit -batch -nojournal export_sphere.jou")
        return 1

    all_results = []

    for o in orders:
        vol_path = os.path.join(_test_dir, f"sphere_o{o}.vol")
        msh_path = os.path.join(_test_dir, f"sphere_o{o}.msh")

        print(f"\n{'─' * 78}")
        print(f"  Order {o}")
        print(f"{'─' * 78}")

        # 1) .vol + NGSolve
        r_vol = ngsolve_vol_integrate(vol_path)
        v_err = (r_vol['volume'] - V_EXACT) / V_EXACT * 100
        a_err = (r_vol['area'] - A_EXACT) / A_EXACT * 100
        print(f"  .vol + NGSolve:     V={r_vol['volume']:.6e} ({v_err:+.2e}%)  "
              f"A={r_vol['area']:.6e} ({a_err:+.2e}%)  "
              f"curve_order={r_vol['curve_order']}")

        # 2) .msh + GMSH API
        try:
            r_gmsh = gmsh_api_integrate(msh_path)
        except Exception as e:
            print(f"  .msh + GMSH API:    FAILED ({e})")
            r_gmsh = None

        if r_gmsh:
            v_err_g = (r_gmsh['volume'] - V_EXACT) / V_EXACT * 100
            a_err_g = (r_gmsh['area'] - A_EXACT) / A_EXACT * 100
            neg_info = ""
            if r_gmsh['n_neg_det_3d'] > 0:
                neg_info = f"  [!{r_gmsh['n_neg_det_3d']} neg det 3D]"
            if r_gmsh['n_neg_det_2d'] > 0:
                neg_info += f"  [!{r_gmsh['n_neg_det_2d']} neg det 2D]"
            print(f"  .msh + GMSH API:    V={r_gmsh['volume']:.6e} ({v_err_g:+.2e}%)  "
                  f"A={r_gmsh['area']:.6e} ({a_err_g:+.2e}%)"
                  f"{neg_info}")

        # 3) .msh + NGSolve ReadGmsh
        try:
            r_read = ngsolve_readgmsh_integrate(msh_path)
            v_err_r = (r_read['volume'] - V_EXACT) / V_EXACT * 100
            a_err_r = (r_read['area'] - A_EXACT) / A_EXACT * 100
            print(f"  .msh + ReadGmsh:    V={r_read['volume']:.6e} ({v_err_r:+.2e}%)  "
                  f"A={r_read['area']:.6e} ({a_err_r:+.2e}%)  "
                  f"curve_order={r_read['curve_order']}")
        except Exception as e:
            r_read = None
            print(f"  .msh + ReadGmsh:    FAILED ({e})")

        # Cross-check: .vol vs .msh (GMSH API)
        if r_gmsh:
            vol_diff = (r_gmsh['volume'] - r_vol['volume']) / r_vol['volume'] * 100
            area_diff = (r_gmsh['area'] - r_vol['area']) / r_vol['area'] * 100
            print(f"  vol vs msh diff:    V={vol_diff:+.4e}%  A={area_diff:+.4e}%")

        all_results.append({
            'order': o,
            'vol_ngsolve': r_vol['volume'],
            'vol_gmsh': r_gmsh['volume'] if r_gmsh else None,
            'vol_readgmsh': r_read['volume'] if r_read else None,
            'area_ngsolve': r_vol['area'],
            'area_gmsh': r_gmsh['area'] if r_gmsh else None,
            'area_readgmsh': r_read['area'] if r_read else None,
        })

    # Summary table
    print(f"\n{'=' * 78}")
    print(f"{'':5} {'--- Volume Error ---':>40}   {'--- Area Error ---':>36}")
    print(f"{'Order':>5} {'.vol+NGSolve':>14} {'.msh+GMSH':>14} {'.msh+ReadGmsh':>14}"
          f"   {'.vol+NGSolve':>12} {'.msh+GMSH':>12} {'.msh+ReadGmsh':>12}")
    print(f"{'─' * 78}")

    for r in all_results:
        ve_v = (r['vol_ngsolve'] - V_EXACT) / V_EXACT * 100
        ve_g = ((r['vol_gmsh'] - V_EXACT) / V_EXACT * 100
                if r['vol_gmsh'] is not None else float('nan'))
        ve_r = ((r['vol_readgmsh'] - V_EXACT) / V_EXACT * 100
                if r['vol_readgmsh'] is not None else float('nan'))
        ae_v = (r['area_ngsolve'] - A_EXACT) / A_EXACT * 100
        ae_g = ((r['area_gmsh'] - A_EXACT) / A_EXACT * 100
                if r['area_gmsh'] is not None else float('nan'))
        ae_r = ((r['area_readgmsh'] - A_EXACT) / A_EXACT * 100
                if r['area_readgmsh'] is not None else float('nan'))

        def fmt_pct(v): return f"{v:+13.2e}%" if not math.isnan(v) else "          N/A"
        def fmt_pct_s(v): return f"{v:+11.2e}%" if not math.isnan(v) else "        N/A"
        print(f"{r['order']:>5} {fmt_pct(ve_v)} {fmt_pct(ve_g)} {fmt_pct(ve_r)}"
              f"   {fmt_pct_s(ae_v)} {fmt_pct_s(ae_g)} {fmt_pct_s(ae_r)}")

    print(f"{'=' * 78}")
    print()
    print("Notes:")
    print("  .vol + NGSolve: uses curvedElements -> best accuracy")
    print("  .msh + GMSH API: isoparametric HO -> coarse mesh has neg det")
    print("  .msh + ReadGmsh: HO nodes NOT used for curving -> linear accuracy")
    return 0


if __name__ == "__main__":
    sys.exit(main())
