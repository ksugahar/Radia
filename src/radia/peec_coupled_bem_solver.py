"""PEEC filament coil and scalar-body SIBC with complete mutual reaction.

Each Picard iterate solves the total body field, including a genus-1 handle
current when requested. The complete per-filament body EMF pairs the SAME
unit-source vector potential and surface-Poisson potential as the forward map.
Convergence requires a coil-current fixed point and a checked body solve.

Terminal resistance includes body dissipation AND the change in coil self
loss caused by current redistribution. Body heat is checked against complete
body reaction; it is not used to manufacture the terminal resistance.
"""

import math
import numpy as np

from radia.peec_bundle import solve_loop_bundle


class CoupledPEECBEMSolver:
    """Iterative PEEC filament coil + workpiece scalar BIE + SIBC.

    The coil is given pre-reduced to filament loop form (``R_f``, ``L_f``
    from ``peec_bundle.build_loop_bundle_impedance``) plus the filament
    polylines and the optional per-filament shape-specific internal
    impedance correction ``Zs_fil``.  The workpiece is the same
    ``ScalarBIESIBCSolver`` used by
    ``CoupledBEMSolver`` (dense or intree-HACApK).
    """

    def __init__(self, filament_paths, R_f, L_f, mesh_wp, Zs_fil=None,
                 wp_order=1, wp_hacapk=False, wp_aca_eps=1e-10,
                 wp_hacapk_leaf=64, wp_hacapk_eta=2.0,
                 wp_gmres_tol=1e-8, wp_gmres_maxiter=500,
                 wp_gmres_restart=80):
        from radia.bem_sibc_solver import ScalarBIESIBCSolver

        self.paths = filament_paths
        self.R_f = np.asarray(R_f, dtype=float)
        self.L_f = np.asarray(L_f, dtype=float)
        self.Zs_fil = None if Zs_fil is None else np.asarray(Zs_fil,
                                                             dtype=complex)
        self.n_filaments = len(filament_paths)
        self.mesh_wp = mesh_wp

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
                hacapk_leaf=int(wp_hacapk_leaf),
                hacapk_eta=float(wp_hacapk_eta))
        else:
            # In-tree Sauter-Schwab Galerkin dense operators (same
            # configuration as CoupledBEMSolver / the weak intree-dense
            # backend); dense SL/DL also enables the loop-DOF extension.
            self.wp_solver = ScalarBIESIBCSolver(
                mesh_wp, order=wp_order, assemble_dense=True,
                use_intree_bem=True, intree_geom_order=1,
                intree_singular_n_q=6, intree_regular_quad_degree=7)
        self.wp_nodes = np.array(
            [[mesh_wp.vertices[i].point[j] for j in range(3)]
             for i in range(mesh_wp.nv)])

        # P1 surface-Poisson phi_inc, prepared once (basis-determined;
        # same contract as CoupledBEMSolver -- the Picard loop rebuilds
        # phi_inc from the updated filament currents every iteration).
        if int(wp_order) != 1:
            raise ValueError(
                f"CoupledPEECBEMSolver supports wp_order=1 only (got "
                f"{wp_order}): the surface-Poisson phi_inc and the "
                f"scattered-current extraction are P1 vertex-nodal.")
        from ngsolve import BND as _BND
        from radia.bem_sibc_solver import SurfacePoissonPhiInc
        self.wp_tris = np.array(
            [[v.nr for v in el.vertices] for el in mesh_wp.Elements(_BND)],
            dtype=np.int64)
        self._phi_poisson = SurfacePoissonPhiInc(self.wp_nodes, self.wp_tris)

    def _bundle_solve(self, omega, I_port, emf=None):
        """Solve the loop-bundle at ``omega`` with optional workpiece emf.

        Returns ``(I_f, Z_port)`` with ``Z_port = V_port / I_port``.
        """
        freq = omega / (2.0 * math.pi)
        I_f, V_port = solve_loop_bundle(
            self.R_f, self.L_f, freq, I_port=I_port,
            Zs_fil=self.Zs_fil, emf=emf)
        Z_port = V_port / I_port if I_port else 0.0 + 0.0j
        return I_f, Z_port

    def solve(self, Z_s, omega, I_port=1.0, max_iter=10, tol=1e-3,
              relax=0.5, verbose=False, loop_dof=False):
        """Couple the complete total body field at every current iterate.

        Convergence is the coil-current fixed-point residual, with a 1e-6
        body residual gate. Terminal resistance includes coil-loss redistribution.
        """
        from radia.biot_savart import h_segments_batch
        from radia.bem_loop_extension import A_from_filaments
        from radia.bem_complete_reaction import (electric_incident_vertex_load,
                                                solve_complete_body, _iterate_complete_current)
        from radia.bem_complete_reaction import strong_surface_impedance, strong_impedance_metadata
        from radia.workpiece_surface import _check_sibc_reaction_power
        if max_iter < 2 or not np.isfinite(tol) or tol <= 0 or not 0 < relax <= 1:
            raise ValueError('Invalid strong-coupling iteration controls')
        if not np.isfinite(omega) or omega <= 0 or not np.isfinite(I_port) or I_port == 0:
            raise ValueError('Strong coupling needs finite positive frequency and nonzero current')
        z = strong_surface_impedance(Z_s, len(self.wp_tris))
        poisson = self._phi_poisson
        centers = self.wp_nodes[self.wp_tris].mean(axis=1)
        areas = poisson._areas
        # Exact unit-filament fields use the SAME functions as the weak path.
        unit_h = np.stack([h_segments_batch(path, self.wp_nodes) for path in self.paths])
        unit_a = np.stack([A_from_filaments(centers, [path], [1.]) for path in self.paths])
        currents_air, z_air = self._bundle_solve(omega, I_port)
        currents = currents_air.copy()
        z_fil = self.R_f.astype(complex)+1j*omega*self.L_f
        if self.Zs_fil is not None:
            z_fil += np.diag(self.Zs_fil)
        def response(currents):
            h_inc = np.einsum('k,kvd->vd', currents, unit_h)
            phi_inc, phi_resid = poisson(h_inc, max_grad_residual=.10)
            def a_inc(points):
                return A_from_filaments(points, self.paths, currents)
            body = solve_complete_body(self.wp_solver, poisson, phi_inc, z, omega,
                a_inc, loop_dof=loop_dof, hacapk=self.wp_hacapk, gmres=self._wp_gmres)
            hload = electric_incident_vertex_load(poisson, z, body['field'])
            emf = (1j*omega*np.einsum('ktd,td,t->k', unit_a, body['current'], areas)
                   + np.einsum('kvd,vd->k', unit_h, hload))
            return dict(body=body, emf=emf)

        def map_current(emf):
            return self._bundle_solve(omega, I_port, emf=emf)[0]

        currents, state, residual, iterations = _iterate_complete_current(
            currents, response, map_current, max_iter=max_iter, tol=tol, relax=relax, verbose=verbose)
        body, emf = state['body'], state['emf']
        # All returned body, emf and current quantities are from the SAME state.
        body_reaction_power = float(.5*np.vdot(currents, emf).real)
        balance = _check_sibc_reaction_power(body['heat'], body_reaction_power)
        coil_work = np.vdot(currents, z_fil @ currents)
        z_port = (coil_work+np.vdot(currents, emf))/abs(I_port)**2
        coil_loss = float(.5*coil_work.real)
        coil_loss_air = float(.5*np.vdot(currents_air, z_fil @ currents_air).real)
        return dict(L_air=float(z_air.imag/omega), R_air=float(z_air.real),
            L_total=float(z_port.imag/omega), R_total=float(z_port.real),
            Delta_L=float((z_port-z_air).imag/omega), Delta_R=float((z_port-z_air).real),
            P_total=body['heat'], H_t_rms=body['h_rms'], iterations=iterations,
            converged=True, coupling_residual=residual, body_residual=body['residual'],
            body_reaction_power_W=body_reaction_power, body_power_balance_relative_error=balance,
            coil_loss_W=coil_loss, coil_loss_air_W=coil_loss_air,
            coil_loss_change_W=coil_loss-coil_loss_air,
            port_power_W=float(.5*abs(I_port)**2*z_port.real),
            n_filaments=self.n_filaments, n_phi_wp=self.wp_solver.ndof,
            I_f=currents, body_emf=emf, wp_c=centers, wp_a=areas,
            wp_J_re=body['current'].real, wp_J_im=body['current'].imag,
            wp_q_tri=body['q_tri'], wp_H_t_tri=body['field'],
            **strong_impedance_metadata(z, areas), **body['loop_metadata'])
