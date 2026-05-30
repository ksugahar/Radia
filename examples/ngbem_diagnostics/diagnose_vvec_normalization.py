"""
Diagnostic: Compare ngsolve.bem Maxwell SLP against manually computed V_vec.

Goal: Determine the normalization convention of MaxwellSingleLayerPotentialOperator.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

MU_0 = 4.0 * np.pi * 1e-7


def gauss_triangle_7pt():
    """7-point Gauss quadrature on reference triangle (0,0)-(1,0)-(0,1)."""
    a1 = 0.059715871789770
    b1 = 0.470142064105115
    a2 = 0.797426985353087
    b2 = 0.101286507323456
    pts = np.array([
        [1.0/3, 1.0/3],
        [a1, b1], [b1, a1], [b1, b1],
        [a2, b2], [b2, a2], [b2, b2],
    ])
    w0 = 0.1125
    w1 = 0.06619707639425
    w2 = 0.06296959027241
    wts = np.array([w0, w1, w1, w1, w2, w2, w2]) * 0.5
    return pts, wts


def map_to_physical(pts_ref, v0, v1, v2):
    xi = pts_ref[:, 0:1]
    eta = pts_ref[:, 1:2]
    return np.asarray(v0) + xi * (np.asarray(v1) - v0) + eta * (np.asarray(v2) - v0)


def compute_double_integral(tri_list_i, tri_list_j, mode='vec'):
    """Compute V_vec or V_div entry using 7pt x 7pt quadrature.

    mode='vec':  int int G(r,r') f_i(r) . f_j(r') dS dS'
    mode='div':  int int G(r,r') div_f_i(r) div_f_j(r') dS dS'
    mode='scalar': int int G(r,r') dS dS'  (per triangle pair)
    """
    ref_pts, ref_wts = gauss_triangle_7pt()
    result = 0.0

    for ti in tri_list_i:
        v0i, v1i, v2i = ti['v0'], ti['v1'], ti['v2']
        area_i = ti['area']
        jac_i = 2.0 * area_i
        pts_i = map_to_physical(ref_pts, v0i, v1i, v2i)

        for tj in tri_list_j:
            v0j, v1j, v2j = tj['v0'], tj['v1'], tj['v2']
            area_j = tj['area']
            jac_j = 2.0 * area_j
            pts_j = map_to_physical(ref_pts, v0j, v1j, v2j)

            pair_val = 0.0
            for a in range(len(ref_wts)):
                for b in range(len(ref_wts)):
                    r_ab = np.linalg.norm(pts_i[a] - pts_j[b])
                    if r_ab < 1e-15:
                        continue
                    G = 1.0 / (4.0 * np.pi * r_ab)
                    w = ref_wts[a] * ref_wts[b] * jac_i * jac_j

                    if mode == 'vec':
                        fi = ti['sign'] * (pts_i[a] - ti['p_opp']) / (2.0 * area_i)
                        fj = tj['sign'] * (pts_j[b] - tj['p_opp']) / (2.0 * area_j)
                        pair_val += w * G * np.dot(fi, fj)
                    elif mode == 'div':
                        di = ti['sign'] / area_i
                        dj = tj['sign'] / area_j
                        pair_val += w * G * di * dj
                    elif mode == 'scalar':
                        pair_val += w * G

            if mode in ('div', 'scalar'):
                # Print individual triangle pair contribution for debugging
                pass
            result += pair_val

    return result


def main():
    from netgen.occ import Box, Pnt, OCCGeometry
    from ngsolve import Mesh, HDivSurface, BilinearForm, InnerProduct, ds, BND
    from ngsolve.bem import MaxwellSingleLayerPotentialOperator
    from ngsolve import TaskManager

    print("=" * 70)
    print("Diagnostic: V_vec normalization check (v2)")
    print("=" * 70)

    # Coarse mesh
    width, height, depth = 0.1, 0.1, 0.002
    maxh = 0.03

    block = Box(Pnt(-width/2, -height/2, -depth/2),
                Pnt(width/2, height/2, depth/2))
    block.solids.name = "conductor"
    block.faces.name = "surface"
    geo = OCCGeometry(block)
    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)

    fes = HDivSurface(mesh, order=0)
    u, v = fes.TnT()

    # Collect active DOFs
    active_dof_set = set()
    for el in mesh.Elements(BND):
        for d_raw in fes.GetDofNrs(el):
            d = d_raw if d_raw >= 0 else ~d_raw
            active_dof_set.add(d)
    active_dofs = sorted(active_dof_set)
    n_active = len(active_dofs)
    dof_to_idx = {d: i for i, d in enumerate(active_dofs)}

    print(f"\nActive boundary DOFs (edges): {n_active}")

    # Assemble ngsolve.bem Maxwell SLP at kappa=1
    V_op = MaxwellSingleLayerPotentialOperator(fes, kappa=1.0, intorder=4)
    V_ngmat = V_op.mat

    # Mass matrix
    M_bf = BilinearForm(fes)
    M_bf += InnerProduct(u.Trace(), v.Trace()) * ds(bonus_intorder=2)
    with TaskManager():
        M_bf.Assemble()
        M_ngmat = M_bf.mat

        # Extract to dense numpy
        V_full = np.zeros((n_active, n_active), dtype=complex)
        M_full = np.zeros((n_active, n_active))

        xv = V_ngmat.CreateColVector()
        yv = V_ngmat.CreateColVector()
        ym = M_ngmat.CreateColVector()

        for i, di in enumerate(active_dofs):
            xv[:] = 0
            xv[di] = 1.0
            V_ngmat.Mult(xv, yv)
            M_ngmat.Mult(xv, ym)
            for j, dj in enumerate(active_dofs):
                V_full[j, i] = yv[dj]
                M_full[j, i] = ym[dj].real

        print(f"V_full diag: min={V_full.diagonal().real.min():.4f}, "
              f"max={V_full.diagonal().real.max():.4f}")
        print(f"V_full diag imag: min={V_full.diagonal().imag.min():.4f}, "
              f"max={V_full.diagonal().imag.max():.4f}")

        # Check symmetry
        sym_err = np.linalg.norm(V_full - V_full.T) / np.linalg.norm(V_full)
        print(f"V_full symmetry error: {sym_err:.2e}")

        # Build edge support map
        edge_vnrs_map = {}
        for edge in mesh.edges:
            if edge.nr in dof_to_idx:
                edge_vnrs_map[edge.nr] = {ev.nr for ev in edge.vertices}

        edge_support = {d: [] for d in active_dofs}
        bnd_elements = list(mesh.Elements(BND))

        for face_idx, el in enumerate(bnd_elements):
            verts_pts = [mesh.vertices[v_el.nr].point for v_el in el.vertices]
            coords = [np.array([p[0], p[1], p[2]]) for p in verts_pts]
            vert_nrs = [v_el.nr for v_el in el.vertices]

            if len(coords) != 3:
                continue

            v0, v1, v2 = coords
            area = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))

            dofs_raw = fes.GetDofNrs(el)
            for d_raw in dofs_raw:
                d = d_raw if d_raw >= 0 else ~d_raw
                sign = +1.0 if d_raw >= 0 else -1.0

                if d not in dof_to_idx:
                    continue

                evnrs = edge_vnrs_map.get(d, set())
                p_opp = None
                for k, vnr in enumerate(vert_nrs):
                    if vnr not in evnrs:
                        p_opp = coords[k]
                        break

                if p_opp is not None:
                    edge_support[d].append({
                        'v0': v0, 'v1': v1, 'v2': v2,
                        'p_opp': p_opp, 'area': area, 'sign': sign,
                        'face_idx': face_idx,
                    })

        # Check: each edge should have exactly 2 support triangles with opposite signs
        print("\n--- Edge support check ---")
        sign_issues = 0
        count_issues = 0
        for d in active_dofs[:20]:
            tris = edge_support[d]
            signs = [t['sign'] for t in tris]
            n_tri = len(tris)
            sign_sum = sum(signs)
            if n_tri != 2:
                count_issues += 1
            if abs(sign_sum) > 0.01:
                sign_issues += 1
            if d in list(active_dofs)[:5]:
                print(f"  DOF {d}: {n_tri} triangles, signs={signs}, sum={sign_sum:.0f}")

        all_sign_sums = []
        all_counts = []
        for d in active_dofs:
            tris = edge_support[d]
            all_counts.append(len(tris))
            all_sign_sums.append(sum(t['sign'] for t in tris))

        print(f"\n  Triangle count per edge: min={min(all_counts)}, max={max(all_counts)}, "
              f"mode={max(set(all_counts), key=all_counts.count)}")
        print(f"  Sign sum per edge: min={min(all_sign_sums):.0f}, max={max(all_sign_sums):.0f}, "
              f"#nonzero={sum(1 for s in all_sign_sums if abs(s) > 0.01)}/{n_active}")

        # ===== KEY TEST: Compare V_div for well-separated pair =====
        print("\n" + "=" * 70)
        print("V_div cancellation test for well-separated pair")
        print("=" * 70)

        # Find a specific well-separated pair
        edge_centers = {}
        for d in active_dofs:
            evnrs = edge_vnrs_map.get(d, set())
            if len(evnrs) == 2:
                vnr_list = list(evnrs)
                p0 = np.array(mesh.vertices[vnr_list[0]].point)
                p1 = np.array(mesh.vertices[vnr_list[1]].point)
                edge_centers[d] = 0.5 * (p0 + p1)

        # Find pair at ~80mm distance
        test_i = active_dofs[0]
        best_j = None
        best_dist = 0
        for d in active_dofs[10:]:
            ci = edge_centers.get(test_i)
            cj = edge_centers.get(d)
            if ci is not None and cj is not None:
                dist = np.linalg.norm(ci - cj)
                if 0.07 < dist < 0.09 and dist > best_dist:
                    best_j = d
                    best_dist = dist

        if best_j is not None:
            di, dj = test_i, best_j
            tri_i = edge_support[di]
            tri_j = edge_support[dj]
            dist = best_dist

            print(f"\nEdge pair: DOF {di} and DOF {dj}, distance = {dist*1e3:.1f} mm")
            print(f"  Edge {di}: {len(tri_i)} triangles, signs = {[t['sign'] for t in tri_i]}")
            print(f"  Edge {dj}: {len(tri_j)} triangles, signs = {[t['sign'] for t in tri_j]}")

            # Compute V_div term by term
            print(f"\n  V_div breakdown (per triangle pair):")
            V_div_total = 0.0
            ref_pts, ref_wts = gauss_triangle_7pt()
            for a_idx, ti in enumerate(tri_i):
                for b_idx, tj in enumerate(tri_j):
                    # Compute scalar integral ∫∫ G dS dS'
                    S_ab = compute_double_integral([ti], [tj], mode='scalar')
                    div_i = ti['sign'] / ti['area']
                    div_j = tj['sign'] / tj['area']
                    contribution = div_i * div_j * S_ab

                    # Also compute simple approximation: G(c_i, c_j) * A_i * A_j
                    ci = sum([ti['v0'], ti['v1'], ti['v2']]) / 3
                    cj = sum([tj['v0'], tj['v1'], tj['v2']]) / 3
                    dist_ij = np.linalg.norm(ci - cj)
                    G_approx = 1.0 / (4 * np.pi * dist_ij)
                    S_approx = G_approx * ti['area'] * tj['area']
                    contrib_approx = div_i * div_j * S_approx

                    V_div_total += contribution
                    print(f"    T{a_idx}(sign={ti['sign']:+.0f})-T{b_idx}(sign={tj['sign']:+.0f}): "
                          f"S={S_ab:.6e}, div_i*div_j={div_i*div_j:.2e}, "
                          f"contrib={contribution:.6e}, approx={contrib_approx:.6e}")

            print(f"  V_div total: {V_div_total:.6e}")

            # Compute V_vec
            V_vec_val = compute_double_integral(tri_i, tri_j, mode='vec')
            print(f"  V_vec:       {V_vec_val:.6e}")

            # ngsolve.bem value
            ii = dof_to_idx[di]
            jj = dof_to_idx[dj]
            V_ng = V_full[ii, jj]
            print(f"  V_ngbem:     {V_ng.real:.6e} + {V_ng.imag:.6e}j")

            # Test conventions
            for name, a, b in [
                ("-Vvec - Vdiv", -1, -1),
                ("-Vvec + Vdiv", -1, +1),
                ("+Vvec - Vdiv", +1, -1),
                ("+Vvec + Vdiv", +1, +1),
            ]:
                predicted = a * V_vec_val + b * V_div_total
                print(f"    {name}: predicted={predicted:.6e}, "
                      f"ratio=V_ng/pred={V_ng.real/predicted:.4f}" if abs(predicted)>1e-20 else "")

        # ===== Also try with kappa=0.001 to separate V_vec and V_div =====
        print("\n" + "=" * 70)
        print("Kappa variation test: extract V_vec and V_div from two kappa values")
        print("=" * 70)

        kappas = [0.01, 1.0]
        V_at_kappa = {}
        for kap in kappas:
            V_kap = MaxwellSingleLayerPotentialOperator(fes, kappa=kap, intorder=4)
            mat_kap = V_kap.mat
            V_dense = np.zeros((n_active, n_active), dtype=complex)
            xk = mat_kap.CreateColVector()
            yk = mat_kap.CreateColVector()
            for i, di in enumerate(active_dofs):
                xk[:] = 0
                xk[di] = 1.0
                mat_kap.Mult(xk, yk)
                for j, dj in enumerate(active_dofs):
                    V_dense[j, i] = yk[dj]
            V_at_kappa[kap] = V_dense
            print(f"kappa={kap}: V diag min={V_dense.diagonal().real.min():.4f}, "
                  f"max={V_dense.diagonal().real.max():.4f}")

        k1, k2 = kappas
        V1, V2 = V_at_kappa[k1], V_at_kappa[k2]

        # Test hypothesis: V(kappa) = alpha(kappa)*V_vec + beta(kappa)*V_div
        # At two kappa values: V(k1) = a1*Vv + b1*Vd, V(k2) = a2*Vv + b2*Vd
        # Solve: Vv = (b2*V1 - b1*V2)/(a1*b2 - a2*b1), etc.

        for hyp_name, a_func, b_func in [
            ("V = -(1/k)Vv - k*Vd",  lambda k: -1.0/k,  lambda k: -k),
            ("V = -(1/k)Vv + k*Vd",  lambda k: -1.0/k,  lambda k: +k),
            ("V = +(1/k)Vv - k*Vd",  lambda k: +1.0/k,  lambda k: -k),
            ("V = +(1/k)Vv + k*Vd",  lambda k: +1.0/k,  lambda k: +k),
            ("V = -k*Vv - (1/k)Vd",  lambda k: -k,      lambda k: -1.0/k),
            ("V = -k*Vv + (1/k)Vd",  lambda k: -k,      lambda k: +1.0/k),
            ("V = +k*Vv - (1/k)Vd",  lambda k: +k,      lambda k: -1.0/k),
        ]:
            a1, b1 = a_func(k1), b_func(k1)
            a2, b2 = a_func(k2), b_func(k2)
            det = a1*b2 - a2*b1

            if abs(det) < 1e-20:
                continue

            Vv_extracted = (b2 * V1 - b1 * V2) / det
            Vd_extracted = (a1 * V2 - a2 * V1) / det  # Note: corrected formula

            # Verify: a1*Vv + b1*Vd should = V1
            V1_reconstructed = a1 * Vv_extracted + b1 * Vd_extracted
            err1 = np.linalg.norm(V1 - V1_reconstructed) / np.linalg.norm(V1)

            # Check Vv properties: should be ~real positive semi-def
            Vv_diag = Vv_extracted.diagonal()
            Vd_diag = Vd_extracted.diagonal()

            print(f"\n  {hyp_name}:")
            print(f"    Reconstruction err: {err1:.2e}")
            print(f"    Vv diag: real min={Vv_diag.real.min():.4f}, max={Vv_diag.real.max():.4f}, "
                  f"imag rms={np.sqrt(np.mean(Vv_diag.imag**2)):.4e}")
            print(f"    Vd diag: real min={Vd_diag.real.min():.4f}, max={Vd_diag.real.max():.4f}, "
                  f"imag rms={np.sqrt(np.mean(Vd_diag.imag**2)):.4e}")

            # Compare extracted Vv with manual Vv for the test pair
            if best_j is not None:
                ii = dof_to_idx[test_i]
                jj = dof_to_idx[best_j]
                Vv_ng = Vv_extracted[ii, jj]
                Vd_ng = Vd_extracted[ii, jj]
                print(f"    Test pair: Vv_extracted={Vv_ng.real:.6e}, Vd_extracted={Vd_ng.real:.6e}")
                print(f"    Test pair: Vv_manual   ={V_vec_val:.6e}, Vd_manual   ={V_div_total:.6e}")
                if abs(V_vec_val) > 1e-20:
                    print(f"    Vv ratio (extracted/manual): {Vv_ng.real/V_vec_val:.4f}")
                if abs(V_div_total) > 1e-20:
                    print(f"    Vd ratio (extracted/manual): {Vd_ng.real/V_div_total:.4f}")


if __name__ == '__main__':
    main()
