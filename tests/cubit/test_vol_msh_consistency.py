"""
Self-consistency test: .vol vs .msh v2.2 export.

Both export paths start from the same Cubit model and should produce
meshes with identical Volume and Area per label (within tolerance).

This validates the Mesh Exporter by cross-checking two independent
code paths:
  - .vol: Python (cubit_netgen_bridge) + C++ (NetgenCurverPure)
  - .msh: C++ (ExportGmshCommand) + ReadGmsh (Netgen)

Requires: Cubit with radia plugin (.ccm loaded).
"""

import os
import sys
import tempfile

# Add src/radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

# Must import netgen/ngsolve before cubit
import netgen  # noqa: F401
from ngsolve import Mesh, Integrate, CF, BND


MODEL = os.path.join(os.path.dirname(__file__),
                     '../../examples/cubit_panels/inductance/'
                     'inductance_torus.cub5')

THRESHOLD_PCT = 0.1  # .vol and .msh should agree within 0.1%


def test_vol_vs_msh(cub5_path=None, order=2):
    """Compare .vol and .msh v2.2 exports from the same Cubit model."""
    from panels.calc_common import setup_cubit
    model = os.path.abspath(cub5_path or MODEL)
    setup_cubit(cub5_path=model, load_plugins=True)
    import cubit

    tmpdir = tempfile.mkdtemp(prefix="vol_msh_check_")
    vol_path = os.path.join(tmpdir, "model.vol")
    msh_path = os.path.join(tmpdir, "model.msh")

    # --- Export .vol (Python + C++ curving) ---
    from cubit_netgen_bridge import extract_curved_mesh
    build_order = max(order, 2)
    ng_mesh = extract_curved_mesh(cubit, order=build_order)
    ng_mesh.Save(vol_path)
    mesh_vol = Mesh(vol_path)

    # --- Export .msh v2.2 (C++ plugin) ---
    if order >= 2:
        # Set 2nd order element type for blocks
        for bid in cubit.get_block_id_list():
            etype = cubit.get_block_element_type(bid)
            if etype in ("TETRA", "TETRA4", "TET", "TET4"):
                cubit.cmd(f"block {bid} element type tetra10")
            elif etype in ("HEX", "HEX8"):
                cubit.cmd(f"block {bid} element type hex20")
    cubit.cmd(f'export gmsh "{msh_path}" overwrite')

    from netgen.read_gmsh import ReadGmsh
    ng_msh = ReadGmsh(msh_path)
    mesh_msh = Mesh(ng_msh)

    # --- Compare Volume per material ---
    print(f"Model: {os.path.basename(model)}, order={order}")
    print(f"  .vol: {mesh_vol.ne} elements, .msh: {mesh_msh.ne} elements")
    print()

    all_pass = True

    # Volume
    print("  Material            .vol Volume      .msh Volume     Diff")
    print("  " + "-" * 64)
    for mat in mesh_vol.GetMaterials():
        v_vol = Integrate(CF(1), mesh_vol, definedon=mesh_vol.Materials(mat))
        # Find matching material in msh mesh
        v_msh = None
        for mat_msh in mesh_msh.GetMaterials():
            if mat_msh == mat:
                v_msh = Integrate(CF(1), mesh_msh,
                                  definedon=mesh_msh.Materials(mat_msh))
                break
        if v_msh is None:
            # Try by index if names don't match
            mats_msh = mesh_msh.GetMaterials()
            idx = list(mesh_vol.GetMaterials()).index(mat)
            if idx < len(mats_msh):
                mat_msh = mats_msh[idx]
                v_msh = Integrate(CF(1), mesh_msh,
                                  definedon=mesh_msh.Materials(mat_msh))

        if v_msh is not None and v_vol > 0:
            err = (v_msh - v_vol) / v_vol * 100
            flag = " ***" if abs(err) > THRESHOLD_PCT else ""
            if flag:
                all_pass = False
            print(f"  {mat:<20s}{v_vol:>14.6e}  {v_msh:>14.6e}  {err:+.2e}%{flag}")
        else:
            print(f"  {mat:<20s}{v_vol:>14.6e}  {'N/A':>14s}")
    print()

    # Area
    print("  Boundary            .vol Area        .msh Area       Diff")
    print("  " + "-" * 64)
    checked = set()
    for bnd in mesh_vol.GetBoundaries():
        if bnd in checked:
            continue
        checked.add(bnd)
        a_vol = Integrate(CF(1), mesh_vol, BND,
                          definedon=mesh_vol.Boundaries(bnd))
        a_msh = None
        for bnd_msh in mesh_msh.GetBoundaries():
            if bnd_msh == bnd:
                a_msh = Integrate(CF(1), mesh_msh, BND,
                                  definedon=mesh_msh.Boundaries(bnd_msh))
                break

        if a_msh is not None and a_vol > 0:
            err = (a_msh - a_vol) / a_vol * 100
            flag = " ***" if abs(err) > THRESHOLD_PCT else ""
            if flag:
                all_pass = False
            print(f"  {bnd:<20s}{a_vol:>14.6e}  {a_msh:>14.6e}  {err:+.2e}%{flag}")
        else:
            print(f"  {bnd:<20s}{a_vol:>14.6e}  {'N/A':>14s}")
    print()

    if all_pass:
        print("  PASSED: .vol and .msh agree within "
              f"{THRESHOLD_PCT}%")
    else:
        print(f"  FAILED: differences exceed {THRESHOLD_PCT}%")

    return all_pass


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description=".vol vs .msh v2.2 self-consistency test")
    parser.add_argument("--cub5", default=None)
    parser.add_argument("--order", type=int, default=2)
    args = parser.parse_args()
    ok = test_vol_vs_msh(args.cub5, args.order)
    sys.exit(0 if ok else 1)
