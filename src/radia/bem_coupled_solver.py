"""
Two-body BEM coupled solver: coil (EFIE) + workpiece (Scalar BIE + SIBC).

Iterative coupling between coil surface current and workpiece eddy current:
  1. EFIE on coil → J_coil (with back-reaction RHS from J_wp)
  2. Biot-Savart J_coil → phi_inc at workpiece
  3. Scalar BIE + SIBC → phi_wp → J_wp
  4. Biot-Savart J_wp → A_back at coil → modify EFIE RHS
  5. Converge → L_total, P, Delta_L

Usage:
    from radia.bem_coupled_solver import CoupledBEMSolver
    solver = CoupledBEMSolver(mesh_coil, mesh_wp, ...)
    result = solver.solve(Z_s, omega)
"""

import math
import time
import numpy as np
from scipy.linalg import solve as scipy_solve, lu_factor, lu_solve

MU_0 = 4e-7 * np.pi


def extract_element_J(mesh, gf_J):
    """Extract per-element current density from HDivSurface GridFunction.

    Returns (centroids, areas, J_vecs) arrays.
    """
    from ngsolve import Integrate, CF, BND

    elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND, element_wise=True)
    elem_Jx = Integrate(gf_J[0], mesh, VOL_or_BND=BND, element_wise=True)
    elem_Jy = Integrate(gf_J[1], mesh, VOL_or_BND=BND, element_wise=True)
    elem_Jz = Integrate(gf_J[2], mesh, VOL_or_BND=BND, element_wise=True)

    centroids, areas, J_vecs = [], [], []
    for el in mesh.Elements(BND):
        area = abs(elem_A[el.nr])
        if area < 1e-30:
            continue
        jvec = np.array([elem_Jx[el.nr], elem_Jy[el.nr],
                         elem_Jz[el.nr]]) / area
        verts = [mesh.vertices[v.nr].point for v in el.vertices]
        c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
        centroids.append(c)
        areas.append(area)
        J_vecs.append(jvec)

    return np.array(centroids), np.array(areas), np.array(J_vecs)


def vector_potential_from_J(obs_points, src_centroids, src_areas, src_J_vecs):
    """Compute vector potential A = (mu0/4pi) * int J/|r-r'| dS' at observation points.

    Args:
        obs_points: (N, 3) observation points
        src_centroids: (M, 3) source element centroids
        src_areas: (M,) source element areas
        src_J_vecs: (M, 3) source current density

    Returns:
        A: (N, 3) vector potential [T*m]
    """
    obs = np.asarray(obs_points)
    src_c = np.asarray(src_centroids)
    src_a = np.asarray(src_areas)
    src_J = np.asarray(src_J_vecs)

    dx = obs[:, None, :] - src_c[None, :, :]       # (N, M, 3)
    r = np.sqrt(np.maximum(np.sum(dx**2, axis=2), 1e-60))  # (N, M)
    # A = (mu0/4pi) * sum_j J_j * area_j / r_j
    weight = (MU_0 / (4 * np.pi)) * src_a[None, :] / r  # (N, M)
    A = np.einsum('nm,mj->nj', weight, src_J)  # (N, 3)
    return A


class CoupledBEMSolver:
    """Two-body BEM: coil (EFIE) + workpiece (Scalar BIE + SIBC).

    Assembles both BEM systems once, then iterates the coupling.
    """

    def __init__(self, mesh_coil, mesh_wp, source_label="source",
                 sink_label="sink", fes_order=0, wp_order=1):
        """Initialize coupled solver.

        Args:
            mesh_coil: NGSolve Mesh for coil (surface mesh with source/sink)
            mesh_wp: NGSolve Mesh for workpiece (surface mesh)
            source_label, sink_label: boundary labels for current injection
            fes_order: HDivSurface order for coil EFIE
            wp_order: H1 order for workpiece BIE
        """
        from ngsolve import (HDivSurface, SurfaceL2, BilinearForm, LinearForm,
                             GridFunction, TaskManager, ds, BND, div)
        from ngsolve.bem import LaplaceSL
        from bem_sibc_solver import ScalarBIESIBCSolver

        self.mesh_coil = mesh_coil
        self.mesh_wp = mesh_wp

        # === Coil EFIE setup ===
        t0 = time.perf_counter()
        fes_J = HDivSurface(mesh_coil, order=fes_order)
        fes_L2 = SurfaceL2(mesh_coil, order=0)
        self.fes_J = fes_J
        self.n_J = fes_J.ndof
        self.n_f = fes_L2.ndof

        # Divergence matrix D
        u_J = fes_J.TrialFunction()
        q = fes_L2.TestFunction()
        bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
        bf_D += div(u_J.Trace()) * q * ds
        bf_D.Assemble()
        from scipy.sparse import coo_matrix
        rows, cols, vals = bf_D.mat.COO()
        self.D = coo_matrix((vals, (rows, cols)),
                            shape=(bf_D.mat.height, bf_D.mat.width)).toarray()

        # LaplaceSL on coil
        jt, jv = fes_J.TnT()
        with TaskManager():
            V_op = LaplaceSL(jt.Trace() * ds, use_fmm=False) * jv.Trace() * ds
        rows, cols, vals = V_op.mat.COO()
        self.SL_coil = coo_matrix((vals, (rows, cols)),
                                  shape=(V_op.mat.height, V_op.mat.width)).toarray()

        # Source/sink RHS
        f_src = LinearForm(fes_L2)
        f_src += q * ds(source_label)
        f_src.Assemble()
        g_src = f_src.vec.FV().NumPy().copy()
        A_src = np.sum(g_src)

        f_snk = LinearForm(fes_L2)
        f_snk += q * ds(sink_label)
        f_snk.Assemble()
        g_snk = f_snk.vec.FV().NumPy().copy()
        A_snk = np.sum(g_snk)

        self.g = g_src / A_src - g_snk / A_snk

        # LU factorize the saddle point matrix (reused each iteration)
        D_red = self.D[:-1, :]
        g_red = self.g[:-1]
        n_J = self.n_J
        n_c = self.n_f - 1
        self.K_saddle = np.block([
            [self.SL_coil, D_red.T],
            [D_red, np.zeros((n_c, n_c))]
        ])
        self.K_lu = lu_factor(self.K_saddle)
        self.g_red = g_red
        self.n_constraint = n_c

        self.t_coil_assembly = time.perf_counter() - t0

        # === Workpiece BIE setup ===
        self.wp_solver = ScalarBIESIBCSolver(mesh_wp, order=wp_order)

        # Workpiece node coordinates
        self.wp_nodes = np.array(
            [[mesh_wp.vertices[i].point[j] for j in range(3)]
             for i in range(mesh_wp.nv)])

        # Coil element data (for Biot-Savart, initially empty)
        self._coil_elems = None

    def solve(self, Z_s, omega, max_iter=10, tol=1e-3, relax=0.5):
        """Solve coupled system iteratively.

        Args:
            Z_s: Surface impedance [Ohm] (complex)
            omega: Angular frequency [rad/s]
            max_iter: Maximum coupling iterations
            tol: Convergence tolerance (relative change in L)
            relax: Under-relaxation for back-reaction

        Returns:
            dict with L_total, L_self, Delta_L, P_total, H_t_rms, iterations
        """
        from ngsolve import GridFunction, Integrate, CF, BND, InnerProduct, grad
        from bem_sibc_solver import (compute_phi_inc_from_surface_J,
                                     biot_savart_surface_current)

        n_J = self.n_J
        n_c = self.n_constraint

        # --- Initial EFIE (no back-reaction) ---
        rhs0 = np.zeros(n_J + n_c)
        rhs0[n_J:] = self.g_red

        x = lu_solve(self.K_lu, rhs0)
        J_coil = x[:n_J]
        L_self = MU_0 * J_coil @ self.SL_coil @ J_coil

        # Set coil GridFunction
        gf_J = GridFunction(self.fes_J)
        gf_J.vec.FV().NumPy()[:] = J_coil

        # Extract coil elements for Biot-Savart
        coil_c, coil_a, coil_J = extract_element_J(self.mesh_coil, gf_J)

        # --- Iteration loop ---
        f_back = np.zeros(n_J)  # back-reaction RHS
        L_prev = L_self

        for iteration in range(max_iter):
            # Forward: coil → workpiece
            phi_inc = compute_phi_inc_from_surface_J(
                self.wp_nodes, coil_c, coil_a, coil_J, n_quad=20)

            # Workpiece BIE-SIBC
            wp_result = self.wp_solver.solve(phi_inc, Z_s=Z_s, omega=omega)
            H_t_rms = wp_result['H_t_rms']
            phi_vec = wp_result['phi_vec']

            # Extract J_wp from phi (tangential gradient)
            # J_wp = -grad_s(phi) per element
            gf_phi_re = GridFunction(self.wp_solver.fes)
            gf_phi_im = GridFunction(self.wp_solver.fes)
            gf_phi_re.vec.FV().NumPy()[:] = phi_vec.real
            gf_phi_im.vec.FV().NumPy()[:] = phi_vec.imag

            elem_A_wp = Integrate(CF(1), self.mesh_wp, BND, element_wise=True)
            wp_c, wp_a, wp_J = [], [], []
            for comp_re, comp_im in [(gf_phi_re, gf_phi_im)]:
                pass  # will compute below

            # Compute per-element J_wp = -grad_s(phi)
            grad_re = [Integrate(grad(gf_phi_re)[i], self.mesh_wp, BND,
                                 element_wise=True) for i in range(3)]
            grad_im = [Integrate(grad(gf_phi_im)[i], self.mesh_wp, BND,
                                 element_wise=True) for i in range(3)]

            wp_c, wp_a, wp_J_re, wp_J_im = [], [], [], []
            for el in self.mesh_wp.Elements(BND):
                area = abs(elem_A_wp[el.nr])
                if area < 1e-30:
                    continue
                # J_wp = -grad_s(phi) (negative tangential gradient)
                j_re = np.array([-grad_re[k][el.nr] / area for k in range(3)])
                j_im = np.array([-grad_im[k][el.nr] / area for k in range(3)])
                verts = [self.mesh_wp.vertices[v.nr].point for v in el.vertices]
                c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
                wp_c.append(c)
                wp_a.append(area)
                wp_J_re.append(j_re)
                wp_J_im.append(j_im)

            wp_c = np.array(wp_c)
            wp_a = np.array(wp_a)
            wp_J_re = np.array(wp_J_re)
            wp_J_im = np.array(wp_J_im)

            # Back-reaction: A_wp at coil element centroids (real part only for DC/MQS)
            # A = mu0 * SL[J_wp], take real part (DC coupling dominates at low freq)
            A_back = vector_potential_from_J(coil_c, wp_c, wp_a, wp_J_re)

            # Project A_back onto coil HDivSurface DOFs via ΔL approach:
            # ΔL = mu0 * integral_S_coil J_coil . A_wp dS
            #    ≈ mu0 * sum_i (J_coil_i . A_wp(c_i)) * area_i
            # This is equivalent to f_back . J_coil = ΔL / mu0
            # For the EFIE RHS modification, we need f_back per DOF.
            # Use: f_back = SL_cross^T @ J_wp_effective
            # But simpler: compute ΔL directly and skip EFIE re-solve.
            #
            # The ΔL approach: L_total = L_self + ΔL
            # ΔL = 2 * mu0 * sum_i (J_coil_i . A_wp_i) * area_coil_i
            Delta_L_elem = MU_0 * np.sum(
                np.sum(coil_J * A_back, axis=1) * coil_a)

            # For f_back projection onto HDivSurface DOFs:
            # Approximate: f_back_dof = sum over elem contributions
            # For each coil element j with area_j and A_back_j:
            #   contribution to DOF i = (A_back_j . n_edge_i) * edge_length * area_j / 3
            # This is complex for general meshes. Instead, use the relation:
            #   f_back . J_coil_vec ≈ sum(J_coil_elem . A_back_elem * area)
            # and distribute proportionally:
            #   f_back ≈ alpha * SL_coil @ J_coil  where alpha = ΔL / (mu0 * L_self)
            # This captures the correct magnitude and distributes via SL structure.
            alpha = Delta_L_elem / max(abs(L_self), 1e-30)
            f_back_new = alpha * (self.SL_coil @ J_coil)

            # Under-relax
            f_back = relax * f_back_new + (1 - relax) * f_back

            # Re-solve EFIE with back-reaction
            rhs = np.zeros(n_J + n_c)
            rhs[:n_J] = -f_back  # negative: SL@J + D^T@p = -f_back
            rhs[n_J:] = self.g_red
            x = lu_solve(self.K_lu, rhs)
            J_coil = x[:n_J]

            # Update coil GridFunction and elements
            gf_J.vec.FV().NumPy()[:] = J_coil
            coil_c, coil_a, coil_J = extract_element_J(self.mesh_coil, gf_J)

            # Compute L
            L_self_new = MU_0 * J_coil @ self.SL_coil @ J_coil
            # Mutual: Delta_L = mu0 * J_coil . A_wp (dot product over coil surface)
            # ΔL = f_back . J_coil (since f_back ≈ ∫ A_wp . v dS)
            Delta_L = float(np.dot(f_back, J_coil))
            L_total = L_self_new + Delta_L

            dL = abs(L_total - L_prev) / max(abs(L_prev), 1e-30)

            P_density = wp_result['P_density']
            P_total = P_density * wp_result['area']

            print(f"  iter {iteration}: L_self={L_self_new*1e9:.2f}nH "
                  f"DeltaL={Delta_L*1e9:.2f}nH L_total={L_total*1e9:.2f}nH "
                  f"H_t={H_t_rms:.2f} P={P_total:.4e} dL={dL:.4e}")

            L_prev = L_total

            if dL < tol and iteration > 0:
                print(f"  Converged at iteration {iteration}")
                break

        return {
            'L_total': float(L_total),
            'L_self': float(L_self_new),
            'Delta_L': float(Delta_L),
            'P_total': float(P_total),
            'P_density': float(P_density),
            'H_t_rms': float(H_t_rms),
            'area_wp': float(wp_result['area']),
            'Z_s': complex(Z_s),
            'iterations': iteration + 1,
            'n_J_coil': self.n_J,
            'n_phi_wp': self.wp_solver.ndof,
        }
