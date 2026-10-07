"""Surface-coil EFIE and scalar-body SIBC with complete mutual reaction.

The coil single-layer matrix uses the unscaled kernel 1/(4 pi r).
Its saddle RHS is -emf_body/(i omega mu0), with terminal constraints unchanged.
The mutual EMF tests the total body surface current and impedance work
against the exact discrete forward incident-field maps. It includes the
handle current during every coupled iteration when loop_dof is enabled.

Coil-current fixed-point and body residuals certify one consistent state.
Terminal work includes complete body reaction; the ideal surface coil has
no ohmic self loss. Body reaction power is checked independently against heat.
"""

import math
import time
import numpy as np
from scipy.linalg import lu_factor, lu_solve

MU_0 = 4e-7 * np.pi


def extract_element_J(mesh, gf_J):
    """Per-element current vectors from an HDivSurface GridFunction.

    Returns ``(centroids, areas, J_vecs)`` arrays. Empty / degenerate
    elements (area < 1e-30) are skipped.
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


class CoupledBEMSolver:
    """Iterative coil EFIE + workpiece scalar BIE + SIBC."""

    def __init__(self, mesh_coil, mesh_wp, source_label="source",
                 sink_label="sink", fes_order=0, wp_order=1,
                 wp_hacapk=False, wp_aca_eps=1e-10, wp_hacapk_leaf=64,
                 wp_hacapk_eta=2.0, wp_gmres_tol=1e-8, wp_gmres_maxiter=500,
                 wp_gmres_restart=80,
                 coil_hacapk=False, coil_aca_eps=1e-8, coil_hacapk_leaf=64,
                 coil_hacapk_eta=2.0):
        from ngsolve import (HDivSurface, SurfaceL2, BilinearForm, LinearForm,
                             ds, BND, div, Compress)
        from ngsolve.bem import LaplaceSL
        from radia.bem_sibc_solver import (ScalarBIESIBCSolver,
                                           SurfacePoissonPhiInc)
        from scipy.sparse import coo_matrix

        self.mesh_coil = mesh_coil
        self.mesh_wp = mesh_wp

        # === Coil EFIE setup (real saddle point factorization) ===
        t0 = time.perf_counter()
        # Compress: the coil mesh is usually the volume .vol whose boundary is
        # the conductor, and HDivSurface on it carries unused interior-edge
        # DOFs (same reasoning as compute_inductance_source_sink).
        fes_J = Compress(HDivSurface(mesh_coil, order=fes_order))
        fes_L2 = SurfaceL2(mesh_coil, order=max(0, fes_order - 1))
        self.fes_J = fes_J
        self.n_J = fes_J.ndof
        self.n_f = fes_L2.ndof

        # Divergence matrix D: n_f x n_J
        u_J = fes_J.TrialFunction()
        q = fes_L2.TestFunction()
        bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
        bf_D += div(u_J.Trace()) * q * ds
        bf_D.Assemble()
        rows, cols, vals = bf_D.mat.COO()
        self.D = coo_matrix((vals, (rows, cols)),
                            shape=(bf_D.mat.height, bf_D.mat.width)).toarray()

        # LaplaceSL on coil
        jt, jv = fes_J.TnT()
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

        if abs(A_src) < 1e-30 or abs(A_snk) < 1e-30:
            raise ValueError(
                f"Source/sink faces empty: A_src={A_src}, A_snk={A_snk}. "
                f"Check that the coil mesh has boundary labels "
                f"'{source_label}' / '{sink_label}'.")

        self.g = g_src / A_src - g_snk / A_snk

        # Coil saddle solve backend.
        #   dense (default): LU-factor [[SL, Dr^T],[Dr, 0]] once and reuse the
        #     factor for every Picard right-hand side (fastest below ~12k n_J).
        #   coil_hacapk: compress SL to an O(N log N) H-matrix and solve the
        #     saddle by loop-COCR (div-free reduction), built ONCE and re-solved
        #     each iteration against the back-reaction rhs_J.  O(N r) storage
        #     instead of the dense LU's O(N^2) -- lifts the coil past the dense
        #     LU memory wall (the dense SL assembly itself still caps at ~12k
        #     n_J; beyond that needs an on-demand H-matrix fill).
        D_red = self.D[:-1, :]
        g_red = self.g[:-1]
        n_c = self.n_f - 1
        self.g_red = g_red
        self.n_constraint = n_c
        self.coil_hacapk = bool(coil_hacapk)
        if self.coil_hacapk:
            from radia.bem.coil_inductance_ngsolve import (
                _LoopReducedSaddle, _dof_cluster_coords)
            # Any fes_order: one cluster point per DOF (edge midpoints for
            # edge DOFs, centroids for face DOFs), attributed through the
            # boundary elements.  Replaces the RT0-only edge-midpoint tree.
            coords = _dof_cluster_coords(mesh_coil, fes_J)
            self._coil_loop = _LoopReducedSaddle(
                self.SL_coil, None, D_red, g_red, 0.0, None, "hacapk",
                coords=coords, hacapk_aca_eps=float(coil_aca_eps),
                hacapk_leaf=int(coil_hacapk_leaf),
                hacapk_eta=float(coil_hacapk_eta))
            # The dense SL is now redundant (the H-matrix holds the compressed
            # operator and serves the energy matvec); free the O(N^2) array.
            self.SL_coil = None
            self.K_lu = None
        else:
            self.K_saddle = np.block([
                [self.SL_coil, D_red.T],
                [D_red, np.zeros((n_c, n_c))]
            ])
            self.K_lu = lu_factor(self.K_saddle)
        self.t_coil_assembly = time.perf_counter() - t0

        # === Workpiece BIE setup ===
        # Default: dense ngsolve.bem scalar BIE (fine for small/moderate wp).
        # wp_hacapk=True: the in-tree Sauter-Schwab Galerkin assembler with an
        # O(N log N) HACApK H-matrix (the weak-path pattern), which scales the
        # workpiece BIE past the ~12k-tri dense-assembly wall (e.g. the 20k-tri
        # reference workpiece).  At order=1 the intree path still creates
        # ``self.fes`` (H1 P1), so the scattered-current extraction in
        # Complete reaction uses the same incident-field maps for either backend.
        self.wp_hacapk = bool(wp_hacapk)
        self._wp_gmres = dict(tol=float(wp_gmres_tol),
                              maxiter=int(wp_gmres_maxiter),
                              restart=int(wp_gmres_restart))
        if self.wp_hacapk:
            self.wp_solver = ScalarBIESIBCSolver(
                mesh_wp, order=wp_order, assemble_dense=True,
                use_intree_bem=True, intree_geom_order=1,
                intree_singular_n_q=6, intree_regular_quad_degree=7,
                use_intree_hacapk=True, hacapk_aca_eps=float(wp_aca_eps),
                hacapk_leaf=int(wp_hacapk_leaf), hacapk_eta=float(wp_hacapk_eta))
        else:
            # In-tree Sauter-Schwab Galerkin dense operators -- the SAME
            # assembler configuration as the weak path's intree-dense
            # backend.  Dense SL/DL also enables the genus-1
            # loop-DOF extension (``loop_dof=True`` in ``solve``).
            self.wp_solver = ScalarBIESIBCSolver(
                mesh_wp, order=wp_order, assemble_dense=True,
                use_intree_bem=True, intree_geom_order=1,
                intree_singular_n_q=6, intree_regular_quad_degree=7)
        self.wp_nodes = np.array(
            [[mesh_wp.vertices[i].point[j] for j in range(3)]
             for i in range(mesh_wp.nv)])

        # Incident-potential reconstruction is basis-determined (same
        # contract as the weak path): P1 -> surface-Poisson psi from the
        # exact vertex H_inc, prepared ONCE (the Picard loop re-derives
        # phi_inc from the updated coil current every iteration, so the
        # cached stiffness factorization is what makes iterations cheap).
        # The former per-iteration path integration (two
        # compute_phi_inc_from_surface_J calls, the dominant per-iteration
        # cost) is no longer part of the production route.
        if int(wp_order) != 1:
            raise ValueError(
                f"CoupledBEMSolver supports wp_order=1 only (got "
                f"{wp_order}): the surface-Poisson phi_inc and the "
                f"scattered-current extraction are P1 vertex-nodal.")
        from ngsolve import BND as _BND
        self.wp_tris = np.array(
            [[v.nr for v in el.vertices] for el in mesh_wp.Elements(_BND)],
            dtype=np.int64)
        self._phi_poisson = SurfacePoissonPhiInc(self.wp_nodes, self.wp_tris)

    def solve(self, Z_s, omega, max_iter=10, tol=1e-3, relax=0.5,
              verbose=False, loop_dof=False):
        """Couple the complete body reaction at each coil-current iterate.

        The coil SL is unscaled: physical f_complete enters as -f_complete/mu0.
        Coil fixed-point and body residuals must pass before returning one state.
        """
        from radia.bem_sibc_solver import H_from_surface_J_complex, A_from_surface_J
        from radia.bem_complete_reaction import (solve_complete_body,
            electric_incident_vertex_load, surface_source_reaction_load,
            surface_current_average_maps, surface_coil_reaction_rhs, _iterate_complete_current)
        from radia.workpiece_surface import _check_sibc_reaction_power
        from ngsolve import GridFunction
        if max_iter < 2 or not np.isfinite(tol) or tol <= 0 or not 0 < relax <= 1:
            raise ValueError('Invalid strong-coupling iteration controls')
        if not np.isfinite(omega) or omega <= 0:
            raise ValueError('Strong coupling needs finite positive frequency')
        z = complex(Z_s)
        n_J, n_c = self.n_J, self.n_constraint
        if self.coil_hacapk:
            def coil_solve(rhs, particular):
                current, _, _ = self._coil_loop.solve(rhs_J=rhs, include_particular=particular)
                return np.real(current)
            def sl_apply(current):
                return self._coil_loop.a11(current)
        else:
            def coil_solve(rhs, particular):
                load = np.zeros(n_J+n_c)
                if rhs is not None:
                    load[:n_J] = rhs
                if particular:
                    load[n_J:] = self.g_red
                return lu_solve(self.K_lu, load)[:n_J]
            def sl_apply(current):
                return self.SL_coil @ current
        air = coil_solve(None, True)
        l_air = MU_0*float(air @ sl_apply(air))
        currents = air.astype(complex)
        maps = surface_current_average_maps(self.fes_J)
        gf = GridFunction(self.fes_J);gf.vec.FV().NumPy()[:] = air
        coil_centers, coil_areas, _ = extract_element_J(self.mesh_coil, gf)
        poisson = self._phi_poisson
        centers = self.wp_nodes[self.wp_tris].mean(axis=1)
        areas = poisson._areas
        def response(currents):
            coil_j = np.stack([m @ currents for m in maps], axis=1)
            incident_h = H_from_surface_J_complex(self.wp_nodes, coil_centers, coil_areas, coil_j)
            phi_inc, _ = poisson(incident_h, max_grad_residual=.10)
            def a_inc(points):
                return A_from_surface_J(points, coil_centers, coil_areas, coil_j)
            body = solve_complete_body(self.wp_solver, poisson, phi_inc, z, omega,
                a_inc, loop_dof=loop_dof, hacapk=self.wp_hacapk, gmres=self._wp_gmres)
            hload = electric_incident_vertex_load(poisson, z, body['field'])
            face_load = surface_source_reaction_load(coil_centers, coil_areas,
                self.wp_nodes, centers, areas, body['current'], hload, omega)
            emf = sum(m.T @ face_load[:, c] for c, m in enumerate(maps))
            return dict(body=body, emf=emf)

        def map_current(emf):
            rhs = surface_coil_reaction_rhs(emf, omega)
            return coil_solve(rhs.real, True)+1j*coil_solve(rhs.imag, False)

        currents, state, residual, iterations = _iterate_complete_current(
            currents, response, map_current, max_iter=max_iter, tol=tol, relax=relax, verbose=verbose)
        body, emf = state['body'], state['emf']
        body_work = np.vdot(currents, emf)
        body_power = float(.5*body_work.real)
        balance = _check_sibc_reaction_power(body['heat'], body_power)
        self_energy = MU_0*float(currents.real @ sl_apply(currents.real)
                                 + currents.imag @ sl_apply(currents.imag))
        port_z = 1j*omega*self_energy+body_work
        return dict(L_air=l_air, R_air=0., L_self=self_energy,
            L_total=float(port_z.imag/omega), R_total=float(port_z.real),
            Delta_L=float(port_z.imag/omega-l_air), Delta_R=float(port_z.real),
            P_total=body['heat'], H_t_rms=body['h_rms'], iterations=iterations,
            converged=True, coupling_residual=residual, body_residual=body['residual'],
            body_reaction_power_W=body_power, body_power_balance_relative_error=balance,
            port_power_W=float(.5*port_z.real), coil_loss_W=0., coil_loss_air_W=0., coil_loss_change_W=0.,
            n_J_coil=n_J, n_phi_wp=self.wp_solver.ndof,
            J_coil_re=currents.real, J_coil_im=currents.imag, body_emf=emf,
            wp_c=centers, wp_a=areas, wp_J_re=body['current'].real,
            wp_J_im=body['current'].imag, wp_q_tri=body['q_tri'], Z_s=z,
            **body['loop_metadata'])
