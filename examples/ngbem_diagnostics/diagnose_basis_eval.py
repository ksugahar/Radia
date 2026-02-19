"""
Diagnostic: Evaluate RT0 basis function directly from NGSolve GridFunction
and compare with manual formula to determine correct sign/normalization.

Strategy:
1. Set one DOF to 1 in a GridFunction
2. Evaluate the GF at quadrature points on the supporting triangles
3. Compare with manual RT0: sign * (x - p_opp) / (2*area)
4. Determine the actual sign and any scaling factor
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))


def main():
    from netgen.occ import Box, Pnt, OCCGeometry
    from ngsolve import (Mesh, HDivSurface, GridFunction, BND,
                         BilinearForm, InnerProduct, ds)
    from ngsolve.bem import MaxwellSingleLayerPotentialOperator

    print("=" * 70)
    print("Diagnostic: Direct RT0 basis function evaluation")
    print("=" * 70)

    # Create mesh
    width, height, depth = 0.1, 0.1, 0.002
    maxh = 0.04  # very coarse for inspection

    block = Box(Pnt(-width/2, -height/2, -depth/2),
                Pnt(width/2, height/2, depth/2))
    block.solids.name = "conductor"
    block.faces.name = "surface"
    geo = OCCGeometry(block)
    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)

    fes = HDivSurface(mesh, order=0)
    n_dof = fes.ndof

    print(f"Total DOFs: {n_dof}")

    # Collect active DOFs and their info
    active_dofs = set()
    for el in mesh.Elements(BND):
        for d_raw in fes.GetDofNrs(el):
            d = d_raw if d_raw >= 0 else ~d_raw
            active_dofs.add(d)
    active_dofs = sorted(active_dofs)
    n_active = len(active_dofs)
    print(f"Active boundary DOFs: {n_active}")

    # Build edge info
    edge_vnrs_map = {}
    edge_center_map = {}
    for edge in mesh.edges:
        if edge.nr in set(active_dofs):
            vnrs = [ev.nr for ev in edge.vertices]
            edge_vnrs_map[edge.nr] = set(vnrs)
            p0 = np.array(mesh.vertices[vnrs[0]].point)
            p1 = np.array(mesh.vertices[vnrs[1]].point)
            edge_center_map[edge.nr] = 0.5 * (p0 + p1)

    # For a few DOFs, evaluate the GridFunction at triangle centroids
    # and compare with manual RT0

    gf = GridFunction(fes)
    bnd_elements = list(mesh.Elements(BND))

    print("\n" + "=" * 70)
    print("Test: Set DOF=1 and evaluate GF at triangle centroids")
    print("=" * 70)

    for test_dof in active_dofs[:5]:
        print(f"\n--- DOF {test_dof} ---")

        # Edge vertices
        evnrs = edge_vnrs_map.get(test_dof, set())
        if not evnrs:
            continue
        evnr_list = list(evnrs)
        ep0 = np.array(mesh.vertices[evnr_list[0]].point)
        ep1 = np.array(mesh.vertices[evnr_list[1]].point)
        print(f"  Edge: v{evnr_list[0]} ({ep0}) -> v{evnr_list[1]} ({ep1})")
        print(f"  Edge center: {0.5*(ep0+ep1)}")
        print(f"  Edge length: {np.linalg.norm(ep1-ep0)*1e3:.2f} mm")

        # Set single DOF
        gf.vec[:] = 0
        gf.vec[test_dof] = 1.0

        # Find supporting triangles
        for face_idx, el in enumerate(bnd_elements):
            dofs_raw = fes.GetDofNrs(el)
            dof_nrs = [d if d >= 0 else ~d for d in dofs_raw]
            signs = [+1.0 if d >= 0 else -1.0 for d in dofs_raw]

            if test_dof not in dof_nrs:
                continue

            # Get triangle info
            verts = [mesh.vertices[v.nr] for v in el.vertices]
            coords = [np.array([v.point[0], v.point[1], v.point[2]])
                      for v in verts]
            vnrs = [v.nr for v in el.vertices]

            if len(coords) != 3:
                continue

            v0, v1, v2 = coords
            area = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
            centroid = (v0 + v1 + v2) / 3.0

            # Find p_opp for this DOF in this triangle
            local_idx = dof_nrs.index(test_dof)
            local_sign = signs[local_idx]
            p_opp = None
            for k, vnr in enumerate(vnrs):
                if vnr not in evnrs:
                    p_opp = coords[k]
                    break

            # Manual RT0 at centroid
            if p_opp is not None:
                f_manual = local_sign * (centroid - p_opp) / (2.0 * area)
            else:
                f_manual = np.array([0, 0, 0])

            # Evaluate GridFunction at centroid
            try:
                # Use mesh integration point
                f_gf = gf(mesh(*centroid))
                f_gf = np.array([f_gf[0], f_gf[1], f_gf[2]])
            except Exception as e:
                f_gf = np.array([np.nan, np.nan, np.nan])
                print(f"  [Warning: Could not evaluate GF: {e}]")

            # Compare
            f_manual_mag = np.linalg.norm(f_manual)
            f_gf_mag = np.linalg.norm(f_gf)
            if f_manual_mag > 1e-20 and f_gf_mag > 1e-20:
                ratio_mag = f_gf_mag / f_manual_mag
                cos_angle = np.dot(f_gf, f_manual) / (f_gf_mag * f_manual_mag)
            else:
                ratio_mag = np.nan
                cos_angle = np.nan

            print(f"  Face {face_idx}: verts={vnrs}, area={area:.6e}")
            print(f"    p_opp = v{vnrs[[k for k in range(3) if vnrs[k] not in evnrs][0] if any(vnrs[k] not in evnrs for k in range(3)) else 0]}")
            print(f"    GetDofNrs sign: {local_sign:+.0f}")
            print(f"    f_manual = [{f_manual[0]:.4f}, {f_manual[1]:.4f}, {f_manual[2]:.4f}]")
            print(f"    f_gf     = [{f_gf[0]:.4f}, {f_gf[1]:.4f}, {f_gf[2]:.4f}]")
            print(f"    |f_manual| = {f_manual_mag:.6f}")
            print(f"    |f_gf|     = {f_gf_mag:.6f}")
            print(f"    ratio |gf|/|manual| = {ratio_mag:.4f}")
            print(f"    cos(angle) = {cos_angle:.6f}")

    # Mass matrix check: compare diagonal computed with GF evaluation vs assembly
    print("\n" + "=" * 70)
    print("Mass matrix verification")
    print("=" * 70)

    u, v = fes.TnT()
    M_bf = BilinearForm(fes)
    M_bf += InnerProduct(u.Trace(), v.Trace()) * ds(bonus_intorder=2)
    M_bf.Assemble()

    for test_dof in active_dofs[:3]:
        gf.vec[:] = 0
        gf.vec[test_dof] = 1.0

        # M[i,i] from assembled bilinear form
        xv = M_bf.mat.CreateColVector()
        yv = M_bf.mat.CreateColVector()
        xv[:] = 0
        xv[test_dof] = 1.0
        M_bf.mat.Mult(xv, yv)
        M_diag_assembled = yv[test_dof].real

        # M[i,i] by manual integration
        M_diag_manual = 0.0
        for el in bnd_elements:
            dofs_raw = fes.GetDofNrs(el)
            dof_nrs = [d if d >= 0 else ~d for d in dofs_raw]
            signs = [+1.0 if d >= 0 else -1.0 for d in dofs_raw]

            if test_dof not in dof_nrs:
                continue

            coords = [np.array([mesh.vertices[v.nr].point[0],
                                mesh.vertices[v.nr].point[1],
                                mesh.vertices[v.nr].point[2]])
                      for v in el.vertices]
            vnrs = [v.nr for v in el.vertices]
            if len(coords) != 3:
                continue

            area = 0.5 * np.linalg.norm(np.cross(coords[1]-coords[0], coords[2]-coords[0]))
            local_idx = dof_nrs.index(test_dof)
            local_sign = signs[local_idx]

            evnrs = edge_vnrs_map.get(test_dof, set())
            p_opp = None
            for k, vnr in enumerate(vnrs):
                if vnr not in evnrs:
                    p_opp = coords[k]
                    break

            if p_opp is None:
                continue

            # int_T |f|^2 dS using centroid rule
            centroid = sum(coords) / 3
            f_at_c = local_sign * (centroid - p_opp) / (2 * area)
            M_contribution = np.dot(f_at_c, f_at_c) * area
            M_diag_manual += M_contribution

        print(f"  DOF {test_dof}: M_assembled={M_diag_assembled:.6f}, "
              f"M_manual={M_diag_manual:.6f}, "
              f"ratio={M_diag_assembled/M_diag_manual:.4f}" if abs(M_diag_manual) > 1e-20 else "")

    # Now do the CORRECT V_vec computation using GF evaluation
    # For each pair, set DOF_i = 1, evaluate GF_i at quadrature points,
    # compute ∫∫ G * GF_i · GF_j
    print("\n" + "=" * 70)
    print("V_vec via GF evaluation (using ngsolve's own basis functions)")
    print("=" * 70)

    # This is too expensive for full matrix. Just do one pair.
    # Pick two DOFs that share a triangle (adjacent)
    test_pairs = []
    for el in bnd_elements[:10]:
        dofs_raw = fes.GetDofNrs(el)
        dof_nrs = [d if d >= 0 else ~d for d in dofs_raw]
        dof_nrs_active = [d for d in dof_nrs if d in set(active_dofs)]
        if len(dof_nrs_active) >= 2:
            test_pairs.append((dof_nrs_active[0], dof_nrs_active[1]))
            if len(test_pairs) >= 2:
                break

    # Also add a distant pair
    di_far = active_dofs[0]
    dj_far = active_dofs[n_active//2]
    test_pairs.append((di_far, dj_far))

    # Assemble V_full for comparison
    V_op = MaxwellSingleLayerPotentialOperator(fes, kappa=1.0, intorder=4)
    V_ngmat = V_op.mat
    dof_to_idx = {d: i for i, d in enumerate(active_dofs)}

    V_full = np.zeros((n_active, n_active), dtype=complex)
    xv = V_ngmat.CreateColVector()
    yv = V_ngmat.CreateColVector()
    for i, di in enumerate(active_dofs):
        xv[:] = 0
        xv[di] = 1.0
        V_ngmat.Mult(xv, yv)
        for j, dj in enumerate(active_dofs):
            V_full[j, i] = yv[dj]

    # For the kappa extraction, also get V at kappa=0.01
    V_op2 = MaxwellSingleLayerPotentialOperator(fes, kappa=0.01, intorder=4)
    V_ngmat2 = V_op2.mat
    V_full2 = np.zeros((n_active, n_active), dtype=complex)
    xv2 = V_ngmat2.CreateColVector()
    yv2 = V_ngmat2.CreateColVector()
    for i, di in enumerate(active_dofs):
        xv2[:] = 0
        xv2[di] = 1.0
        V_ngmat2.Mult(xv2, yv2)
        for j, dj in enumerate(active_dofs):
            V_full2[j, i] = yv2[dj]

    # Extract V_vec using: V(k) = -(1/k)*Vv - k*Vd
    # At k1=0.01: V1 = -100*Vv - 0.01*Vd
    # At k2=1.0:  V2 = -Vv - Vd
    # Vv = (V1 - 100*V2) / (-100 + 100) ... hmm det=0 for this form
    # Actually: a1=-100, b1=-0.01, a2=-1, b2=-1
    # det = a1*b2 - a2*b1 = (-100)(-1) - (-1)(-0.01) = 100 - 0.01 = 99.99
    # Vv = (b2*V1 - b1*V2)/det = (-V1 + 0.01*V2)/99.99
    # Vd = (a1*V2 - a2*V1)/det = (-100*V2 + V1)/99.99

    k1, k2 = 0.01, 1.0
    a1, b1 = -1.0/k1, -k1  # -100, -0.01
    a2, b2 = -1.0/k2, -k2  # -1, -1
    det = a1*b2 - a2*b1  # 100 - 0.01 = 99.99
    Vv_ng = (b2 * V_full2 - b1 * V_full) / det
    Vd_ng = (a1 * V_full - a2 * V_full2) / det

    print(f"\nExtracted V_vec from kappa variation:")
    print(f"  Vv_ng diag: real [{Vv_ng.diagonal().real.min():.4f}, {Vv_ng.diagonal().real.max():.4f}]")
    print(f"  Vd_ng diag: real [{Vd_ng.diagonal().real.min():.6f}, {Vd_ng.diagonal().real.max():.6f}]")

    # Check Vv eigenvalues
    Vv_real = Vv_ng.real
    Vv_sym = 0.5 * (Vv_real + Vv_real.T)
    eigvals_vv = np.linalg.eigvalsh(Vv_sym)
    print(f"  Vv eigenvalues: min={eigvals_vv.min():.6f}, max={eigvals_vv.max():.6f}")
    print(f"  Vv positive definite: {eigvals_vv.min() > 0}")

    # The PHYSICAL inductance matrix is L = mu_0 * Vv_ng (in full space)
    # Let's check if this makes sense
    print(f"\n  mu_0 * Vv eigenvalues (H): [{MU_0*eigvals_vv.min():.4e}, {MU_0*eigvals_vv.max():.4e}]")

    # Build loop basis and project
    from scipy.linalg import null_space

    # Build divergence matrix D
    D = np.zeros((len(bnd_elements), n_active))
    for face_idx, el in enumerate(bnd_elements):
        for d_raw in fes.GetDofNrs(el):
            d = d_raw if d_raw >= 0 else ~d_raw
            sign = +1.0 if d_raw >= 0 else -1.0
            if d in dof_to_idx:
                D[face_idx, dof_to_idx[d]] = sign

    T_loop = null_space(D)
    n_loops = T_loop.shape[1]
    print(f"\n  Loop basis: {n_loops} loops from {n_active} edges, {len(bnd_elements)} faces")

    # Since ALL signs are +1, D has all +1 entries. Let me check.
    print(f"  D non-zero values: min={D[D!=0].min():.0f}, max={D[D!=0].max():.0f}")
    print(f"  D entries per row: {[int(np.sum(D[i]!=0)) for i in range(min(5, len(bnd_elements)))]}")

    # Project Vv to loop space
    Vv_LL = T_loop.T @ Vv_real @ T_loop
    eigvals_vv_LL = np.linalg.eigvalsh(0.5*(Vv_LL + Vv_LL.T))
    print(f"\n  Vv_LL eigenvalues: min={eigvals_vv_LL.min():.6f}, max={eigvals_vv_LL.max():.6f}")
    print(f"  mu_0 * Vv_LL eigenvalues (H): [{MU_0*eigvals_vv_LL.min():.4e}, {MU_0*eigvals_vv_LL.max():.4e}]")

    # CHECK: Is the divergence matrix wrong because all signs are +1?
    # On a closed surface, the correct divergence matrix should have
    # opposite signs for the two triangles sharing an edge.
    # If GetDofNrs gives +1 for both, then D*x gives the WRONG divergence.
    #
    # The CORRECT divergence needs the orientation flip from the
    # outward normal convention. Let me build the correct D.

    print("\n" + "=" * 70)
    print("Building CORRECT divergence matrix with proper orientation")
    print("=" * 70)

    # For each triangle T with vertices (v0, v1, v2), the edges are:
    #   e0 = (v1, v2) opposite v0
    #   e1 = (v0, v2) opposite v1
    #   e2 = (v0, v1) opposite v2
    #
    # The RT0 divergence on T for each local edge basis function is +1/area.
    # But the SURFACE divergence of the GLOBAL basis function needs signs
    # that account for edge orientation relative to the triangle.
    #
    # Key: for edge e shared by T1 and T2, the outward edge normal from T1
    # points OPPOSITE to the outward edge normal from T2.
    # So if the global direction is "from T1 to T2", then:
    #   T1: outward = +global → div contribution = +1
    #   T2: outward = -global → div contribution = -1
    #
    # The actual sign can be determined by the cross product of the
    # triangle normal and the edge direction.

    # Build face normals
    face_normals = []
    for el in bnd_elements:
        coords = [np.array([mesh.vertices[v.nr].point[0],
                            mesh.vertices[v.nr].point[1],
                            mesh.vertices[v.nr].point[2]])
                  for v in el.vertices]
        if len(coords) == 3:
            normal = np.cross(coords[1] - coords[0], coords[2] - coords[0])
            normal = normal / np.linalg.norm(normal)
        else:
            normal = np.array([0, 0, 0])
        face_normals.append(normal)

    # For each edge, determine orientation relative to each triangle
    # Edge direction: from vertex with smaller nr to larger nr (convention)
    D_correct = np.zeros((len(bnd_elements), n_active))

    for face_idx, el in enumerate(bnd_elements):
        coords = [np.array([mesh.vertices[v.nr].point[0],
                            mesh.vertices[v.nr].point[1],
                            mesh.vertices[v.nr].point[2]])
                  for v in el.vertices]
        vnrs = [v.nr for v in el.vertices]
        if len(coords) != 3:
            continue

        n_face = face_normals[face_idx]
        dofs_raw = fes.GetDofNrs(el)

        for d_raw in dofs_raw:
            d = d_raw if d_raw >= 0 else ~d_raw
            dof_sign = +1.0 if d_raw >= 0 else -1.0

            if d not in dof_to_idx:
                continue

            idx = dof_to_idx[d]
            evnrs = edge_vnrs_map.get(d, set())

            # Find the two edge vertices on this face
            edge_verts_on_face = [k for k in range(3) if vnrs[k] in evnrs]
            if len(edge_verts_on_face) != 2:
                continue

            k0, k1_idx = edge_verts_on_face
            # Edge direction within this triangle: from v[k0] to v[k1]
            edge_dir = coords[k1_idx] - coords[k0]

            # Outward edge normal (in the surface plane, pointing away from face)
            # = cross(face_normal, edge_direction), then choose sign
            # so it points AWAY from the triangle (toward the opposite vertex)
            n_edge_outward = np.cross(n_face, edge_dir)
            n_edge_outward /= (np.linalg.norm(n_edge_outward) + 1e-30)

            # Check direction: should point toward the opposite vertex
            opp_idx = 3 - k0 - k1_idx  # the remaining vertex
            to_opp = coords[opp_idx] - 0.5 * (coords[k0] + coords[k1_idx])
            # Project to surface plane
            to_opp_surf = to_opp - np.dot(to_opp, n_face) * n_face

            # n_edge_outward should point AWAY from triangle = away from opposite vertex
            if np.dot(n_edge_outward, to_opp_surf) > 0:
                # Points toward opp vertex, need to flip
                n_edge_outward = -n_edge_outward

            # Global edge direction: from smaller to larger vertex nr
            vnr_list = sorted(evnrs)
            p_small = np.array(mesh.vertices[vnr_list[0]].point)
            p_large = np.array(mesh.vertices[vnr_list[1]].point)
            edge_global_dir = p_large - p_small
            edge_global_dir /= (np.linalg.norm(edge_global_dir) + 1e-30)

            # Global edge normal: cross(face_normal, global_edge_dir)
            n_edge_global = np.cross(n_face, edge_global_dir)
            n_edge_global /= (np.linalg.norm(n_edge_global) + 1e-30)

            # Sign: +1 if outward aligns with global normal, -1 otherwise
            orientation_sign = np.sign(np.dot(n_edge_outward, n_edge_global))

            D_correct[face_idx, idx] = orientation_sign

    # Check sign cancellation
    sign_sums_correct = []
    for idx in range(n_active):
        col = D_correct[:, idx]
        nz = col[col != 0]
        sign_sums_correct.append(np.sum(nz))

    print(f"D_correct sign sums: min={min(sign_sums_correct):.0f}, "
          f"max={max(sign_sums_correct):.0f}, "
          f"#zero={sum(1 for s in sign_sums_correct if abs(s) < 0.01)}/{n_active}")

    # Build corrected loop basis
    T_loop_correct = null_space(D_correct)
    n_loops_correct = T_loop_correct.shape[1]
    print(f"Corrected loop basis: {n_loops_correct} loops")

    # Verify div-free
    res = np.linalg.norm(D_correct @ T_loop_correct)
    print(f"||D_correct * T_loop_correct|| = {res:.2e}")

    # Project Vv to corrected loop space
    Vv_LL_correct = T_loop_correct.T @ Vv_real @ T_loop_correct
    eigvals_correct = np.linalg.eigvalsh(0.5*(Vv_LL_correct + Vv_LL_correct.T))
    print(f"\nVv_LL_correct eigenvalues: min={eigvals_correct.min():.6f}, max={eigvals_correct.max():.6f}")
    print(f"mu_0 * Vv_LL_correct eigenvalues (H): "
          f"[{MU_0*eigvals_correct.min():.4e}, {MU_0*eigvals_correct.max():.4e}]")

    # Also project V_full (Maxwell SLP at kappa=1)
    V_LL_correct = T_loop_correct.T @ V_full.real @ T_loop_correct
    eigvals_V_correct = np.linalg.eigvalsh(0.5*(V_LL_correct + V_LL_correct.T))
    print(f"\nV_LL_correct (Maxwell SLP) eigenvalues: "
          f"min={eigvals_V_correct.min():.6f}, max={eigvals_V_correct.max():.6f}")

    # Compare with old (wrong) loop basis
    V_LL_old = T_loop.T @ V_full.real @ T_loop
    eigvals_V_old = np.linalg.eigvalsh(0.5*(V_LL_old + V_LL_old.T))
    print(f"V_LL_old (wrong D) eigenvalues: "
          f"min={eigvals_V_old.min():.6f}, max={eigvals_V_old.max():.6f}")


MU_0 = 4.0 * np.pi * 1e-7

if __name__ == '__main__':
    main()
