"""
PEEC coil <-> workpiece scalar BIE + SIBC strong coupling.

Sibling of ``bem_coupled_solver.CoupledBEMSolver``: SAME iterative
workpiece back-reaction, but the coil is a PEEC filament loop-bundle
instead of an impedance-EFIE surface saddle.

Why this exists
===============

The weak Telegen path (``calc_inductance.py`` ``--coupling-mode weak``)
evaluates the workpiece surface field from the INCIDENT (bare-coil
Biot-Savart) field, so it misses the magnetic workpiece's flux
redistribution / eddy screening and can over-estimate heating under strong
loading.  The BEM-A coil has a self-consistent path via
``CoupledBEMSolver``; this module gives the PEEC coil the corresponding
strong-coupling formulation.

Coupling structure (mirrors CoupledBEMSolver, in impedance space)
=================================================================

The PEEC loop-bundle system is already complex (Z_fil = R_f + jw L_f),
so -- unlike the real EFIE saddle -- there is no Re/Im split: the complex
workpiece back-reaction is injected as one complex per-filament EMF.

    Forward  (coil -> workpiece):
        I_f  ->  phi_inc(obs) = compute_phi_inc_from_filaments(...)
        (identical to the validated weak PEEC forward)

    Workpiece BIE + SIBC (ScalarBIESIBCSolver, shared with CoupledBEMSolver):
        phi_inc  ->  phi_vec, H_t_rms, P_density
        J_wp = n x H_scat = extract_scattered_wp_J(...)

    Back      (workpiece -> coil):
        A_wp(r) = (mu0/4pi) sum_j J_wp[j] area[j] / |r - c_j|
        flux_wp,k = integral over filament k of A_wp . dl        (line integral)
        emf[k]    = jw * flux_wp,k
        (filament analogue of CoupledBEMSolver's f_back = int v . A_wp)

    Re-solve  (coupled bundle):
        [ Z_fil  -1 ] [I_f   ]   [-emf  ]
        [ 1^T     0 ] [V_port] = [I_port]

Picard-iterated (with under-relaxation on emf) until the terminal
impedance Z_port = V_port / I_port converges.

R comes from energy, not the terminal reaction (like CoupledBEMSolver)
=====================================================================

The terminal reflected resistance ``Delta_R = Re(Z_port) - R_air`` from
the coupled loop solve does NOT generally equal
``2 P_wp / |I_port|^2`` (the workpiece SIBC dissipation).  The coupled
solve remains internally reciprocal through the filament back-EMF.  This
is the same
SIBC reaction-vs-power adjoint mismatch that ``CoupledBEMSolver`` faces:
the scalar-potential forward / vector-potential back reaction is a correct
adjoint for the REACTIVE coupling (Delta_L) but not for the DISSIPATIVE
part.  That is exactly why ``CoupledBEMSolver`` reports R from
``2 P_wp / I^2`` (energy) and uses the reaction only for L.  This solver
follows the same rule: the caller (``calc_inductance._assemble_strong_output``)
takes R from ``2 P_wp / I^2``; the coupled solve supplies the coil-current
redistribution + ``Delta_L``.  ``Delta_R`` (from Z_port) is returned only
as a diagnostic.

VALIDATION STATUS (EXPERIMENTAL)
===============================

The strong-coupling heating change is driven by coil-current redistribution
under the back-EMF.  The committed demo is weakly coupled, so it cannot
exercise that effect; it locks only wiring and self-consistency.  A durable
strong-loading reference case still needs to be added to the validation
lane on a compute host.  Until that check passes, treat absolute
strong-PEEC P_wp / Delta_L as UNVERIFIED.

Strong coupling does not supply the harmonic loop current omitted by a
single-valued scalar potential on a flux-linked genus-1 workpiece.  That
limitation is shared by the weak and strong workpiece paths; use the
explicit loop-DOF extension where supported.  The analytic sphere and
shorted-ring goldens lock the genus-0 and loop-extension contracts.

Part of the Radia project.
"""

import math
import time
import numpy as np

from radia.bem_coupled_solver import MU_0, extract_scattered_wp_J
from radia.peec_bundle import solve_loop_bundle


def filament_back_emf(seg_mid, seg_dl, seg_fil, n_fil, wp_c, wp_a,
                      wp_J_complex, omega, chunk_target=4_000_000):
    """Per-filament back-reaction EMF from the workpiece surface current.

    ``emf[k] = jw * flux_wp,k`` with
    ``flux_wp,k = sum_{seg in fil k} A_wp(mid_seg) . dl_seg`` and
    ``A_wp(r) = (mu0/4pi) sum_j wp_J[j] wp_a[j] / |r - wp_c[j]|``.

    Args:
        seg_mid: (S, 3) segment midpoints (flattened over all filaments).
        seg_dl:  (S, 3) segment vectors p2 - p1.
        seg_fil: (S,) int filament index of each segment.
        n_fil:   number of filaments K.
        wp_c:    (M, 3) workpiece panel centroids.
        wp_a:    (M,) workpiece panel areas.
        wp_J_complex: (M, 3) complex scattered surface current density.
        omega:   angular frequency [rad/s].
        chunk_target: cap on ``chunk_rows * M`` to bound the distance
            matrix memory.

    Returns:
        (K,) complex per-filament EMF.
    """
    seg_mid = np.asarray(seg_mid, dtype=float)
    seg_dl = np.asarray(seg_dl, dtype=float)
    wp_c = np.asarray(wp_c, dtype=float)
    wp_a = np.asarray(wp_a, dtype=float)
    wp_J = np.asarray(wp_J_complex, dtype=complex)

    K = int(n_fil)
    flux = np.zeros(K, dtype=complex)
    S = seg_mid.shape[0]
    M = wp_c.shape[0]
    if S == 0 or M == 0:
        return 1j * omega * flux

    inv_4pi_mu0 = MU_0 / (4.0 * np.pi)
    # Source moments J[:,d] * area, one per component.
    JA = wp_J * wp_a[:, None]                       # (M, 3) complex

    chunk = max(1, int(chunk_target // max(M, 1)))
    for s0 in range(0, S, chunk):
        s1 = min(S, s0 + chunk)
        mids = seg_mid[s0:s1]                        # (cs, 3)
        diff = mids[:, None, :] - wp_c[None, :, :]   # (cs, M, 3)
        dist = np.sqrt(np.einsum('smd,smd->sm', diff, diff))  # (cs, M)
        inv = np.where(dist > 1e-30, 1.0 / dist, 0.0)         # (cs, M)
        # A[:, d] = inv_4pi_mu0 * inv @ JA[:, d]
        A = inv_4pi_mu0 * (inv @ JA)                 # (cs, 3) complex
        contrib = np.einsum('sd,sd->s', A, seg_dl[s0:s1])     # (cs,) complex
        np.add.at(flux, seg_fil[s0:s1], contrib)

    return 1j * omega * flux


class CoupledPEECBEMSolver:
    """Iterative PEEC filament coil + workpiece scalar BIE + SIBC.

    The coil is given pre-reduced to filament loop form (``R_f``, ``L_f``
    from ``peec_bundle.build_loop_bundle_impedance``) plus the filament
    polylines and the optional per-filament Bessel skin impedance
    ``Zs_fil``.  The workpiece is the same ``ScalarBIESIBCSolver`` used by
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

        # Flatten filament polylines to segment midpoints / dl / owner index
        # once, for the back-reaction line integral.
        mids, dls, owner = [], [], []
        for k, fil in enumerate(filament_paths):
            for (p1, p2) in fil:
                p1 = np.asarray(p1, dtype=float)
                p2 = np.asarray(p2, dtype=float)
                mids.append(0.5 * (p1 + p2))
                dls.append(p2 - p1)
                owner.append(k)
        self._seg_mid = np.asarray(mids, dtype=float) if mids else \
            np.zeros((0, 3))
        self._seg_dl = np.asarray(dls, dtype=float) if dls else \
            np.zeros((0, 3))
        self._seg_fil = np.asarray(owner, dtype=np.int64)

        # Workpiece BIE (identical construction to CoupledBEMSolver).
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
            self.wp_solver = ScalarBIESIBCSolver(mesh_wp, order=wp_order)
        self.wp_nodes = np.array(
            [[mesh_wp.vertices[i].point[j] for j in range(3)]
             for i in range(mesh_wp.nv)])

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
              relax=0.5, verbose=False):
        """Run the iterative coupled PEEC<->workpiece solve.

        Args:
            Z_s: workpiece surface impedance (complex scalar, or ndarray of
                length ``wp_solver.ndof`` for per-node SIBC).
            omega: angular frequency [rad/s].
            I_port: terminal drive current amplitude [A].
            max_iter, tol, relax: Picard controls (tol on |Z_port|).

        Returns dict with ``L_air``, ``R_air``, ``L_total``, ``R_total``,
        ``Delta_L``, ``Delta_R``, ``P_total``, ``H_t_rms``, ``iterations``,
        ``n_filaments``, ``n_phi_wp`` and the workpiece per-panel viz arrays.
        """
        from radia.bem_sibc_solver import compute_phi_inc_from_filaments

        if max_iter < 2:
            raise ValueError("max_iter must be >= 2 for coupling convergence")
        if tol <= 0.0:
            raise ValueError("tol must be > 0 for coupling convergence")
        if not 0.0 < relax <= 1.0:
            raise ValueError("relax must satisfy 0 < relax <= 1")

        # --- Step 0: uncoupled (air) loop-bundle ---
        _, Z_air = self._bundle_solve(omega, I_port, emf=None)
        L_air = Z_air.imag / omega if omega > 0 else 0.0
        R_air = Z_air.real

        I_f, Z_port = self._bundle_solve(omega, I_port, emf=None)
        emf = np.zeros(self.n_filaments, dtype=complex)
        Z_prev = Z_air
        wp_result = None
        wp_c = wp_a = wp_J_re = wp_J_im = None
        iteration = 0
        dZ_rel = float("inf")
        converged = False

        for iteration in range(max_iter):
            # --- Forward: filament currents -> phi_inc at workpiece ---
            phi_inc = compute_phi_inc_from_filaments(
                self.wp_nodes, self.paths, I_f)

            # --- Workpiece scalar BIE + SIBC ---
            if self.wp_hacapk:
                wp_result = self.wp_solver.solve_hacapk(
                    phi_inc, Z_s=Z_s, omega=omega,
                    tol=self._wp_gmres["tol"],
                    maxiter=self._wp_gmres["maxiter"],
                    restart=self._wp_gmres["restart"])
            else:
                wp_result = self.wp_solver.solve(phi_inc, Z_s=Z_s, omega=omega)
            phi_vec = wp_result['phi_vec']

            # --- Scattered surface current + per-filament back-EMF ---
            wp_c, wp_a, wp_J_re, wp_J_im = extract_scattered_wp_J(
                self.mesh_wp, self.wp_solver.fes, phi_vec, phi_inc)
            emf_new = filament_back_emf(
                self._seg_mid, self._seg_dl, self._seg_fil,
                self.n_filaments, wp_c, wp_a,
                wp_J_re + 1j * wp_J_im, omega)

            if iteration == 0:
                emf = emf_new
            else:
                emf = relax * emf_new + (1.0 - relax) * emf

            # --- Re-solve coupled loop-bundle ---
            I_f, Z_port = self._bundle_solve(omega, I_port, emf=emf)

            dZ_rel = abs(Z_port - Z_prev) / max(abs(Z_prev), 1e-30)
            if verbose:
                L_now = Z_port.imag / omega
                print(f"  iter {iteration}: L_total={L_now * 1e9:.3f}nH "
                      f"R_total={Z_port.real * 1e3:.4f}mOhm "
                      f"dZ={dZ_rel:.3e}")
            Z_prev = Z_port
            if dZ_rel < tol and iteration > 0:
                converged = True
                break

        if not converged:
            raise RuntimeError(
                "PEEC-workpiece strong coupling did not converge: "
                f"relative terminal-impedance change={dZ_rel:.3e}, "
                f"tol={tol:.3e}, max_iter={max_iter}.  Increase "
                "--coupling-max-iter or reduce --coupling-relax.")

        L_total = Z_port.imag / omega if omega > 0 else 0.0
        R_total = Z_port.real
        P_total = wp_result['P_density'] * wp_result['area']
        H_t_rms = wp_result['H_t_rms']

        if isinstance(Z_s, np.ndarray):
            Z_s_out = complex(np.mean(Z_s))
        else:
            Z_s_out = complex(Z_s)

        return {
            'L_air': float(L_air),
            'R_air': float(R_air),
            'L_total': float(L_total),
            'R_total': float(R_total),
            'Delta_L': float(L_total - L_air),
            'Delta_R': float(R_total - R_air),
            'P_total': float(P_total),
            'H_t_rms': float(H_t_rms),
            'iterations': iteration + 1,
            'converged': True,
            'coupling_residual': float(dZ_rel),
            'n_filaments': int(self.n_filaments),
            'n_phi_wp': int(self.wp_solver.ndof),
            'I_f': I_f,
            'wp_c': wp_c,
            'wp_a': wp_a,
            'wp_J_re': wp_J_re,
            'wp_J_im': wp_J_im,
            'Z_s': Z_s_out,
        }
