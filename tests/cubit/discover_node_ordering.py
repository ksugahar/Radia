"""
Discover how GMSH reads BDF/VTK HEX20 node ordering.

Step 1: Verify linear hex (HEX8) works through BDF
Step 2: Use GMSH to generate HEX20, export to BDF, read back
Step 3: Create manual HEX20 BDF with correct format

Usage:
  python tests/cubit/discover_node_ordering.py
"""

import sys
import os
import numpy as np

out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ho_volume_test')
os.makedirs(out_dir, exist_ok=True)


def gmsh_check(filename):
    """Open file in GMSH, compute volume and check Jacobians."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.open(filename)
    except Exception as e:
        gmsh.finalize()
        return None, -1, str(e)

    etypes, etags, ntags = gmsh.model.mesh.getElements(dim=3)
    if len(etypes) == 0:
        gmsh.finalize()
        return 0.0, 0, "no 3D elements"

    vol = 0.0
    n_neg = 0
    info = []
    for idx, et in enumerate(etypes):
        name, _, order, nn, _, _ = gmsh.model.mesh.getElementProperties(int(et))
        n_elem = len(etags[idx])
        lc, w = gmsh.model.mesh.getIntegrationPoints(int(et), "Gauss4")
        _, det, _ = gmsh.model.mesh.getJacobians(int(et), lc)
        det = np.array(det)
        neg = int((det < 0).sum())
        n_neg += neg
        n_gp = len(w)
        for i in range(n_elem):
            for j in range(n_gp):
                vol += det[i * n_gp + j] * w[j]
        info.append(f"{name}={n_elem}")
    gmsh.finalize()
    return vol, n_neg, ", ".join(info)


def test_linear_hex_bdf():
    """Step 1: Verify GMSH reads linear hex BDF correctly."""
    print("=" * 60)
    print("Step 1: Linear HEX8 BDF (sanity check)")
    print("=" * 60)

    verts = [
        (0, 0, 0), (1, 0, 0), (1, 2, 0), (0, 2, 0),
        (0, 0, 3), (1, 0, 3), (1, 2, 3), (0, 2, 3),
    ]
    V_exact = 6.0

    fname = os.path.join(out_dir, "manual_hex8.bdf")
    with open(fname, 'w') as f:
        f.write("BEGIN BULK\n")
        for i, (x, y, z) in enumerate(verts):
            f.write(f"GRID    {i+1:8d}{0:8d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
        f.write(f"CHEXA   {1:8d}{1:8d}"
                f"{1:8d}{2:8d}{3:8d}{4:8d}{5:8d}{6:8d}+\n")
        f.write(f"+       {7:8d}{8:8d}\n")
        f.write("PSOLID         1       1\n")
        f.write("ENDDATA\n")

    vol, neg, info = gmsh_check(fname)
    print(f"  HEX8 BDF: V={vol}, neg_det={neg}, {info}")
    print(f"  Expected: V={V_exact}")
    assert neg == 0 and abs(vol - V_exact) < 0.01, "Linear hex FAILED!"
    print("  -> PASS\n")


def test_gmsh_generated_hex20():
    """Step 2: Let GMSH generate HEX20, export to BDF, read back."""
    import gmsh

    print("=" * 60)
    print("Step 2: GMSH-generated HEX20 -> BDF round-trip")
    print("=" * 60)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("hex20_gen")
    gmsh.model.occ.addBox(0, 0, 0, 1, 2, 3)
    gmsh.model.occ.synchronize()

    # Transfinite single hex
    for _, tag in gmsh.model.getEntities(1):
        gmsh.model.mesh.setTransfiniteCurve(tag, 2)
    for _, tag in gmsh.model.getEntities(2):
        gmsh.model.mesh.setTransfiniteSurface(tag)
        gmsh.model.mesh.setRecombine(2, tag)
    for _, tag in gmsh.model.getEntities(3):
        gmsh.model.mesh.setTransfiniteVolume(tag)

    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)  # HEX20, not HEX27
    gmsh.model.mesh.generate(3)

    # Check what we got
    etypes, etags, ntags = gmsh.model.mesh.getElements(dim=3)
    for idx, et in enumerate(etypes):
        name, _, order, nn, _, _ = gmsh.model.mesh.getElementProperties(int(et))
        print(f"  Generated: {name} x {len(etags[idx])}, {nn} nodes")

    # Get connectivity and coords
    nt_all, co_all, _ = gmsh.model.mesh.getNodes()
    nc = {}
    for i, tag in enumerate(nt_all):
        nc[int(tag)] = np.array([co_all[3*i], co_all[3*i+1], co_all[3*i+2]])

    # Show GMSH internal connectivity
    for idx, et in enumerate(etypes):
        name, _, _, nn, _, _ = gmsh.model.mesh.getElementProperties(int(et))
        if nn != 20:
            continue
        conn = [int(x) for x in ntags[idx][:nn]]
        print(f"\n  GMSH internal connectivity (node IDs): {conn}")
        vp = [nc[conn[i]] for i in range(8)]
        print(f"  Vertices:")
        for i in range(8):
            print(f"    v{i} = node {conn[i]:3d}: ({vp[i][0]:.1f}, {vp[i][1]:.1f}, {vp[i][2]:.1f})")
        print(f"  Mid-edge nodes:")
        for i in range(8, 20):
            mp = nc[conn[i]]
            for a in range(8):
                for b in range(a+1, 8):
                    expected = 0.5 * (vp[a] + vp[b])
                    if np.linalg.norm(mp - expected) < 0.01:
                        print(f"    gmsh[{i-8:2d}] = mid(v{a},v{b}) node {conn[i]:3d}: "
                              f"({mp[0]:.1f}, {mp[1]:.1f}, {mp[2]:.1f})")
                        break

    # Export BDF
    bdf_file = os.path.join(out_dir, "gmsh_gen_hex20.bdf")
    msh_file = os.path.join(out_dir, "gmsh_gen_hex20.msh")
    gmsh.write(bdf_file)
    gmsh.write(msh_file)
    gmsh.finalize()

    # Show raw BDF content
    print(f"\n  Raw BDF content:")
    with open(bdf_file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(("CHEXA", "+", "GRID")):
                print(f"    {line}")

    # Read back
    for fname, fmt in [(msh_file, "MSH"), (bdf_file, "BDF")]:
        vol, neg, info = gmsh_check(fname)
        if vol is None:
            print(f"\n  {fmt} readback: ERROR: {info}")
        else:
            print(f"\n  {fmt} readback: V={vol:.4f}, neg_det={neg}, {info}")


def test_gmsh_edge_order_all_types():
    """Step 3: Print GMSH reference edge ordering for all quadratic 3D types."""
    import gmsh

    print("\n" + "=" * 60)
    print("Step 3: GMSH Reference Edge Ordering (all types)")
    print("=" * 60)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    types_to_check = [
        (11, "TET10",  4, 6),
        (17, "HEX20",  8, 12),
        (18, "PRISM15", 6, 9),
        (19, "PYRAMID13", 5, 8),
    ]

    nastran_edge_tables = {
        "TET10": [(0,1),(1,2),(0,2),(0,3),(1,3),(2,3)],
        "HEX20": [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)],
        "PRISM15": [(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)],
        "PYRAMID13": [(0,1),(1,2),(2,3),(3,0),(0,4),(1,4),(2,4),(3,4)],
    }

    for gmsh_type, name, nv, ne in types_to_check:
        gname, _, _, nn, ref_coords, _ = gmsh.model.mesh.getElementProperties(gmsh_type)
        ref_pts = np.array(ref_coords).reshape(-1, 3)

        print(f"\n  {name} (GMSH type {gmsh_type}, {gname}, {nn} nodes):")

        ref_verts = ref_pts[:nv]
        gmsh_edges = []
        for i in range(nv, nv + ne):
            mp = ref_pts[i]
            found = False
            for a in range(nv):
                for b in range(a+1, nv):
                    expected = 0.5 * (ref_verts[a] + ref_verts[b])
                    if np.linalg.norm(mp - expected) < 0.01:
                        gmsh_edges.append((a, b))
                        found = True
                        break
                if found:
                    break
            if not found:
                gmsh_edges.append((-1, -1))

        nastran_edges = nastran_edge_tables[name]

        print(f"    {'Pos':>3} {'GMSH':>12} {'Nastran':>12} {'Match':>6}")
        for i in range(ne):
            g = gmsh_edges[i] if i < len(gmsh_edges) else (-1, -1)
            n = nastran_edges[i] if i < len(nastran_edges) else (-1, -1)
            match = "YES" if set(g) == set(n) else "NO"
            print(f"    {i:3d} {str(g):>12} {str(n):>12} {match:>6}")

        # Compute reorder: nastran pos i -> gmsh pos
        reorder = []
        for i, ne_pair in enumerate(nastran_edges):
            ns = set(ne_pair)
            for j, ge in enumerate(gmsh_edges):
                if set(ge) == ns:
                    reorder.append(j)
                    break
            else:
                reorder.append(-1)
        print(f"    Reorder (Nastran->GMSH): {reorder}")

    gmsh.finalize()


def test_cubit_bdf_vs_msh():
    """Step 4: Compare Cubit BDF and MSH output for the cylinder (from earlier test)."""
    print("\n" + "=" * 60)
    print("Step 4: Existing Cubit exports - compare BDF vs MSH")
    print("=" * 60)

    # Check if the cylinder files exist from the earlier test run
    bdf_file = os.path.join(out_dir, "cylinder_Nastran_o2.bdf")
    msh_file = os.path.join(out_dir, "cylinder_GMSH_v4.1_o2.msh")

    if not os.path.exists(bdf_file):
        print(f"  {bdf_file} not found (run test_ho_volume_all_formats.py first)")
        return

    import gmsh

    for fname, fmt in [(msh_file, "MSH"), (bdf_file, "BDF")]:
        if not os.path.exists(fname):
            print(f"  {fname} not found")
            continue
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(fname)

        etypes, etags, ntags = gmsh.model.mesh.getElements(dim=3)
        nt2, co2, _ = gmsh.model.mesh.getNodes()
        nc = {}
        for i, tag in enumerate(nt2):
            nc[int(tag)] = np.array([co2[3*i], co2[3*i+1], co2[3*i+2]])

        for idx, et in enumerate(etypes):
            name, _, order, nn, _, _ = gmsh.model.mesh.getElementProperties(int(et))
            if nn < 20:
                continue
            # Show first element
            conn = [int(x) for x in ntags[idx][:nn]]
            vp = [nc[conn[i]] for i in range(8)]
            print(f"\n  {fmt}: {name} (first element)")
            print(f"    Vertices:")
            for i in range(8):
                print(f"      v{i} = node {conn[i]:4d}: "
                      f"({vp[i][0]:+.5f}, {vp[i][1]:+.5f}, {vp[i][2]:+.5f})")
            print(f"    Mid-edge nodes:")
            for i in range(8, min(20, nn)):
                mp = nc[conn[i]]
                for a in range(8):
                    for b in range(a+1, 8):
                        expected = 0.5 * (vp[a] + vp[b])
                        if np.linalg.norm(mp - expected) < 0.01:
                            print(f"      pos[{i-8:2d}] = mid(v{a},v{b}) node {conn[i]:4d}: "
                                  f"({mp[0]:+.5f}, {mp[1]:+.5f}, {mp[2]:+.5f})")
                            break
                    else:
                        continue
                    break
                else:
                    print(f"      pos[{i-8:2d}] = ??? node {conn[i]:4d}: "
                          f"({mp[0]:+.5f}, {mp[1]:+.5f}, {mp[2]:+.5f})")

        gmsh.finalize()


def main():
    test_linear_hex_bdf()
    test_gmsh_generated_hex20()
    test_gmsh_edge_order_all_types()
    test_cubit_bdf_vs_msh()
    return 0


if __name__ == "__main__":
    sys.exit(main())
