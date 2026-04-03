"""Test SetCD2Name for edge labels (BBND) in Cubit -> Netgen export.

Requires: Cubit installed, inductance_torus.cub5 model.
Run with system Python 3.12 (not Cubit's embedded Python).
"""
import sys
import os

# Add src/radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

# Must import netgen/ngsolve before cubit
import netgen  # noqa: F401
import ngsolve

MODEL = os.path.join(os.path.dirname(__file__),
                     '../../examples/cubit_panels/inductance/inductance_torus.cub5')


def test_setcd2name():
    """Test that SetCD2Name works without crash."""
    from panels.calc_common import setup_cubit
    setup_cubit(cub5_path=os.path.abspath(MODEL))
    import cubit

    # List curves for reference
    curve_ids = list(cubit.parse_cubit_list("curve", "all"))
    print(f"Cubit curves ({len(curve_ids)}): {curve_ids}")
    for cid in curve_ids:
        ename = cubit.get_entity_name("curve", cid)
        parent_surfs = list(cubit.parse_cubit_list("surface", f"in curve {cid}"))
        n_edges = len(list(cubit.parse_cubit_list("edge", f"in curve {cid}")))
        print(f"  Curve {cid}: name='{ename}', surfs={parent_surfs}, edges={n_edges}")

    # Add a sideset on a curve to test priority 1
    cubit.cmd('sideset 100 add curve 1')
    cubit.cmd('sideset 100 name "test_edge"')

    # Extract mesh with edge labels
    from cubit_mesh_export import extract_curved_mesh
    print("\nExtracting curved mesh (order=2)...", flush=True)
    sys.stdout.flush()
    ng_mesh = extract_curved_mesh(cubit, order=2)
    print("Mesh extracted OK", flush=True)

    # Create NGSolve mesh
    print("Creating NGSolve mesh...")
    mesh = ngsolve.Mesh(ng_mesh)

    # Check BBoundaries
    bbnd = mesh.GetBBoundaries()
    print(f"\nBBoundaries ({len(bbnd)}): {tuple(bbnd)}")
    print(f"Materials: {tuple(mesh.GetMaterials())}")
    print(f"Boundaries: {tuple(mesh.GetBoundaries())}")

    # Check .vol file contents
    import tempfile
    vol_path = os.path.join(tempfile.gettempdir(), 'test_cd2.vol')
    ng_mesh.Save(vol_path)
    print(f"\nSaved .vol to {vol_path}")

    # Check if cd2names are in the vol file
    with open(vol_path, 'r') as f:
        for line in f:
            if 'cd2name' in line.lower() or 'edgesegment' in line.lower():
                print(f"  VOL: {line.rstrip()}")

    # Also check direct Netgen mesh
    print(f"\nNetgen GetNCD2Names: {ng_mesh.GetNCD2Names()}")
    for i in range(ng_mesh.GetNCD2Names()):
        try:
            print(f"  cd2name[{i}]: '{ng_mesh.GetCD2Name(i)}'")
        except Exception:
            pass

    # Check segments
    segs = list(ng_mesh.Elements1D())
    print(f"\nSegments: {len(segs)}")
    for seg in segs[:5]:
        print(f"  {seg}")

    # Verify sideset name propagated
    if bbnd:
        print("\nSUCCESS: SetCD2Name works!")
        found_test = any('test_edge' in name for name in bbnd)
        print(f"  Sideset 'test_edge' found: {found_test}")
    else:
        print("\nBBoundaries empty - investigating...")


if __name__ == '__main__':
    test_setcd2name()
