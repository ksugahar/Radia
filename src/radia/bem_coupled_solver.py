"""
Two-body BEM coupled solver: coil EFIE + workpiece scalar BIE + SIBC.

Iterative coupling with proper per-DOF back-reaction RHS.

Mathematical setup
==================

Coil EFIE saddle point (real, time-harmonic with workpiece reflection):

    [SL_coil  D^T] [J]   [-f_back]
    [D         0 ] [p] = [ g_red ]

solved twice per Picard iteration — once for the real (J_re, f_back_re,
g_red) RHS and once for the imaginary (J_im, f_back_im, 0) RHS — so the
SAME LU factorization can be reused for both.

Workpiece scalar BIE + SIBC (handles complex Z_s, complex phi_inc):

    A_sys phi = M phi_inc                with A_sys = 1/2 M - DL + gamma SL M^-1 K
    -> phi_vec_complex
    -> J_wp = n x H_scat                  where H_scat = -grad_s(phi_total - phi_inc)

Back-reaction
=============

    A_wp(r) = (mu0/4pi) sum_j J_wp[j] * area[j] / |r - centroid[j]|
    f_back[i] = int v_i.Trace() . A_wp dS_coil       (LinearForm assembly)

Real and imaginary parts of A_wp are handled separately. The CoefficientFunction
is built as a sum of M analytic 1/r kernels in (x, y, z), which NGSolve
compiles to a flat expression tree and integrates against the HDivSurface
test functions on the coil surface.

This is the FIX for the v1 bug where f_back was a scalar rescale of
SL @ J_coil. The scalar rescale loses all spatial information about A_wp
at the coil surface and gives the wrong sign of Delta_L.

Inductance from energy
======================

Self magnetic energy (with back-reacted complex coil current):

    W_self = (1/2) mu0 (J_re^T SL_coil J_re + J_im^T SL_coil J_im)

Mutual magnetic energy (coil-workpiece coupling):

    W_mut  = (1/2) (f_back_re^T J_re + f_back_im^T J_im)

Coil terminal effective inductance for unit current amplitude:

    L_total = 2 W_total = mu0 (J_re^T SL J_re + J_im^T SL J_im)
                          + (f_back_re . J_re + f_back_im . J_im)

Reference values:

    L_air     = mu0 J_re_uncoupled^T SL J_re_uncoupled    (no workpiece)
    Delta_L   = L_total - L_air                            (sign-correct)

For non-magnetic conductors with strong screening, J_re drops (current is
pushed out of the original air pattern by the workpiece reflection) AND
the mutual term is negative (Lenz). Result: Delta_L < 0.

For ferromagnetic workpiece, the imaginary part of Z_s gives a strong
in-phase storage in the skin layer, and the mutual term goes positive.
Result: Delta_L > 0.

History
=======

v1 (2026-04 initial): had a SCALAR back-reaction
    f_back_v1 = alpha * SL_coil @ J_coil
which lost all spatial information about A_wp at the coil and produced
wrong-signed Delta_L (verified 2026-04-12: copper at 1 kHz gave +0.72 nH
instead of an expected small negative). See
memory/bem_coupled_solver_existing.md.
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


def assemble_back_reaction_RHS(fes_J, wp_c, wp_a, wp_J):
    """Assemble ``f_back[i] = int v_i.Trace() . A_wp dS_coil``.

    A_wp is the vector potential at the coil surface from the workpiece
    induced surface current::

        A_wp(r) = (mu0/4pi) sum_j wp_J[j] * wp_a[j] / |r - wp_c[j]|

    The CoefficientFunction is built as a sum of M analytic 1/r kernels
    in (x, y, z), then integrated against the HDivSurface test space on
    the coil via a LinearForm assembly.

    Args:
        fes_J: NGSolve HDivSurface space on the coil
        wp_c:  (M, 3) workpiece panel centroids
        wp_a:  (M,)   workpiece panel areas
        wp_J:  (M, 3) workpiece panel REAL surface current density

    Returns:
        np.ndarray of length ``fes_J.ndof``.
    """
    from ngsolve import (CoefficientFunction, LinearForm, InnerProduct,
                         ds, sqrt, x, y, z)

    wp_c = np.asarray(wp_c, dtype=float)
    wp_a = np.asarray(wp_a, dtype=float)
    wp_J = np.asarray(wp_J, dtype=float)
    M = len(wp_c)
    if M == 0:
        return np.zeros(fes_J.ndof)

    A_components = [CoefficientFunction(0.0) for _ in range(3)]
    inv_4pi_mu0 = MU_0 / (4.0 * np.pi)

    for j in range(M):
        cx = float(wp_c[j, 0])
        cy = float(wp_c[j, 1])
        cz = float(wp_c[j, 2])
        r_inv = 1.0 / sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2)
        weight = inv_4pi_mu0 * float(wp_a[j]) * r_inv
        for k in range(3):
            jval = float(wp_J[j, k])
            if jval == 0.0:
                continue
            A_components[k] = A_components[k] + weight * jval

    A_cf = CoefficientFunction(tuple(A_components))

    u_J = fes_J.TestFunction()
    f_form = LinearForm(fes_J)
    f_form += InnerProduct(u_J.Trace(), A_cf) * ds
    f_form.Assemble()
    return f_form.vec.FV().NumPy().copy()


def extract_scattered_wp_J(mesh_wp, wp_fes, phi_vec_complex, phi_inc_complex):
    """Scattered workpiece surface current ``J_wp = n x H_scat`` per BND panel.

    ``H_scat = -grad_s(phi_total - phi_inc)`` computed per element via
    NGSolve element-wise integration of the surface-gradient components.
    Returns ``(centroids, areas, J_re, J_im)`` arrays of shape (M, 3) /
    (M,) over the workpiece boundary elements.

    Shared by ``CoupledBEMSolver`` (EFIE coil) and
    ``CoupledPEECBEMSolver`` (PEEC filament coil): the workpiece side of
    both coupled solvers is identical, so the scattered-current extraction
    lives here as a free function rather than being duplicated.
    """
    from ngsolve import GridFunction, Integrate, CF, BND, grad

    gf_re = GridFunction(wp_fes)
    gf_im = GridFunction(wp_fes)
    gf_re.vec.FV().NumPy()[:] = phi_vec_complex.real - phi_inc_complex.real
    gf_im.vec.FV().NumPy()[:] = phi_vec_complex.imag - phi_inc_complex.imag

    elem_A = Integrate(CF(1), mesh_wp, BND, element_wise=True)
    grad_re = [Integrate(grad(gf_re)[i], mesh_wp, BND,
                         element_wise=True) for i in range(3)]
    grad_im = [Integrate(grad(gf_im)[i], mesh_wp, BND,
                         element_wise=True) for i in range(3)]

    c_list, a_list, jr_list, ji_list = [], [], [], []
    for el in mesh_wp.Elements(BND):
        area = abs(elem_A[el.nr])
        if area < 1e-30:
            continue
        ht_re = np.array([-grad_re[k][el.nr] / area for k in range(3)])
        ht_im = np.array([-grad_im[k][el.nr] / area for k in range(3)])
        verts = [np.array(mesh_wp.vertices[v.nr].point)
                 for v in el.vertices]
        e1 = verts[1] - verts[0]
        e2 = verts[2] - verts[0]
        n_vec = np.cross(e1, e2)
        n_mag = np.linalg.norm(n_vec)
        if n_mag > 1e-30:
            n_vec /= n_mag
        centroid = np.mean(verts, axis=0)
        if np.dot(n_vec, centroid) < 0:
            n_vec = -n_vec
        j_re = np.cross(n_vec, ht_re)
        j_im = np.cross(n_vec, ht_im)
        c_list.append(np.mean([(v[0], v[1], v[2]) for v in verts],
                              axis=0))
        a_list.append(area)
        jr_list.append(j_re)
        ji_list.append(j_im)

    return (np.array(c_list), np.array(a_list),
            np.array(jr_list), np.array(ji_list))


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
                             TaskManager, ds, BND, div)
        from ngsolve.bem import LaplaceSL
        from radia.bem_sibc_solver import (ScalarBIESIBCSolver,
                                           SurfacePoissonPhiInc)
        from scipy.sparse import coo_matrix

        self.mesh_coil = mesh_coil
        self.mesh_wp = mesh_wp

        # === Coil EFIE setup (real saddle point factorization) ===
        t0 = time.perf_counter()
        fes_J = HDivSurface(mesh_coil, order=fes_order)
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
                _LoopReducedSaddle, _edge_midpoint_coords)
            if fes_order != 0:
                raise ValueError(
                    "coil_hacapk requires fes_order==0 (RT0) for the "
                    "edge-midpoint HACApK cluster tree.")
            coords = _edge_midpoint_coords(mesh_coil, self.n_J)
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
        # Takahashi workpiece).  At order=1 the intree path still creates
        # ``self.fes`` (H1 P1), so the scattered-current extraction in
        # ``_extract_wp_J`` is unchanged; only the SL/DL storage + solve differ.
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
            # backend.  (The former default here was the ngsolve.bem
            # column-matvec dense extraction: O(N^3), ~an hour at 3k DOF,
            # with a NaN incident on record -- a superseded route, removed
            # 2026-07-17.)  Dense SL/DL also enables the genus-1
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
        # cost) was removed 2026-07-17 with the weak path's legacy route.
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
        """Run the iterative coupled solve.

        Args:
            Z_s: workpiece surface impedance.
                - **complex scalar**: legacy uniform-Z_s SIBC.
                - **ndarray of length self.wp_solver.ndof (complex)**:
                  per-node Z_s for the per-panel curvature SIBC.
                  A caller that wants per-node SIBC (e.g. the
                  validation-lane
                  ``validation_test/induction_heating/bem_reference/
                  calc_inductance.py::_run_coupled_bem`` local-curvature
                  path) builds this array by computing per-panel local
                  curvature from the workpiece mesh and projecting the
                  resulting per-panel Z_s onto H1 nodes via vertex
                  averaging. The ScalarBIESIBCSolver assembles the
                  Robin term with diag(gamma) so each node sees its
                  own SIBC coefficient.  The production
                  ``panels/calc_inductance.py`` ``--coupling-mode strong``
                  path passes a single global scalar Z_s.
            omega: angular frequency [rad/s]
            max_iter: Picard iteration cap
            tol: relative L_total convergence
            relax: under-relaxation (0..1)

            loop_dof: apply the genus-1 loop-DOF extension
                (``radia.bem_loop_extension.solve_loop_extended``) ONCE on
                the CONVERGED state: the Picard loop itself runs the plain
                scalar BIE (whose L_total / Delta_L convention is the
                validated one), then the shorted-turn current alpha is
                solved against the converged coil current and the reported
                ``P_total`` / ``H_t_rms`` are replaced by the loop-extended
                values (the same dissipation-only convention as the weak
                path's ``--wp-loop-dof``).  Requires the dense wp backend
                (``wp_hacapk=False``) and a genus-1 workpiece; fails loud
                otherwise.  The alpha back-reaction onto the coil current
                is NOT iterated (the coil-current redistribution is a
                percent-level effect where strong coupling applies at all).

        Returns dict with ``L_air``, ``L_total``, ``Delta_L``, ``P_total``,
        ``H_t_rms``, ``iterations``, ``J_coil_re``, ``J_coil_im`` (+
        ``wp_loop_alpha`` / ``wp_loop_theta_jump`` /
        ``wp_loop_cut_n_vertices`` / ``t_loop_dof_s`` when ``loop_dof``).
        """
        from ngsolve import GridFunction

        if loop_dof and self.wp_hacapk:
            raise ValueError(
                "loop_dof=True needs the dense wp backend (wp_hacapk="
                "False): the HACApK backend exposes no dense SL/DL for "
                "the loop column.")

        n_J = self.n_J
        n_c = self.n_constraint

        # Coil saddle solve + magnetic-energy, dispatched to the dense-LU or
        # the loop-COCR (HACApK) backend.  ``include_particular`` picks the
        # constraint RHS: True = terminal drive g_red (the real current);
        # False = zero net current (the imaginary back-reaction current).
        if self.coil_hacapk:
            def _coil_solve(rhs_J, include_particular):
                J, _it, _st = self._coil_loop.solve(
                    rhs_J=rhs_J, include_particular=include_particular)
                return np.real(J)

            def _sl_energy(J):
                return float(np.real(J @ self._coil_loop.a11(J)))
        else:
            def _coil_solve(rhs_J, include_particular):
                rhs = np.zeros(n_J + n_c)
                if rhs_J is not None:
                    rhs[:n_J] = rhs_J
                if include_particular:
                    rhs[n_J:] = self.g_red
                return lu_solve(self.K_lu, rhs)[:n_J]

            def _sl_energy(J):
                return float(J @ self.SL_coil @ J)

        # === Step 0: uncoupled coil solution (air-only L) ===
        J_re_air = _coil_solve(None, True).copy()
        L_air = MU_0 * _sl_energy(J_re_air)

        J_re = J_re_air.copy()
        J_im = np.zeros(n_J)
        f_back_re = np.zeros(n_J)
        f_back_im = np.zeros(n_J)

        gf_J_re = GridFunction(self.fes_J)
        gf_J_im = GridFunction(self.fes_J)
        gf_J_re.vec.FV().NumPy()[:] = J_re
        gf_J_im.vec.FV().NumPy()[:] = J_im

        L_prev = L_air
        Delta_L = 0.0
        wp_result = None
        iteration = 0

        for iteration in range(max_iter):
            # --- Forward: J_coil (complex) -> phi_inc (complex) at workpiece ---
            coil_c, coil_a, coil_J_re_arr = extract_element_J(
                self.mesh_coil, gf_J_re)
            _, _, coil_J_im_arr = extract_element_J(
                self.mesh_coil, gf_J_im)

            # Surface-Poisson psi from the exact vertex H_inc (prepared
            # factorization; 10% grad-consistency gate = the fail-loud
            # contract shared with the weak path).
            from radia.bem_sibc_solver import H_from_surface_J_complex
            H_inc = H_from_surface_J_complex(
                self.wp_nodes, coil_c, coil_a,
                coil_J_re_arr + 1j * coil_J_im_arr)
            phi_inc_cplx, phi_resid = self._phi_poisson(
                H_inc, max_grad_residual=0.10)
            if verbose:
                print(f"  iter {iteration}: phi_inc poisson "
                      f"grad-residual {phi_resid:.1%}")

            # --- Workpiece scalar BIE + SIBC (returns complex phi_vec) ---
            if self.wp_hacapk:
                wp_result = self.wp_solver.solve_hacapk(
                    phi_inc_cplx, Z_s=Z_s, omega=omega,
                    tol=self._wp_gmres["tol"],
                    maxiter=self._wp_gmres["maxiter"],
                    restart=self._wp_gmres["restart"])
            else:
                wp_result = self.wp_solver.solve(
                    phi_inc_cplx, Z_s=Z_s, omega=omega)
            phi_vec = wp_result['phi_vec']

            # --- Extract scattered surface current J_wp ---
            (wp_c, wp_a, wp_J_re_arr,
             wp_J_im_arr) = self._extract_wp_J(phi_vec, phi_inc_cplx)

            # --- Per-DOF back-reaction RHS via LinearForm assembly ---
            f_back_re_new = assemble_back_reaction_RHS(
                self.fes_J, wp_c, wp_a, wp_J_re_arr)
            f_back_im_new = assemble_back_reaction_RHS(
                self.fes_J, wp_c, wp_a, wp_J_im_arr)

            if iteration == 0:
                f_back_re = f_back_re_new
                f_back_im = f_back_im_new
            else:
                f_back_re = (relax * f_back_re_new
                             + (1 - relax) * f_back_re)
                f_back_im = (relax * f_back_im_new
                             + (1 - relax) * f_back_im)

            # --- Re-solve coil EFIE saddle point (real and imag parts) ---
            # SL J_re + D^T p_re = -f_back_re;  D J_re = g_red (terminal drive)
            J_re = _coil_solve(-f_back_re, True)
            # SL J_im + D^T p_im = -f_back_im;  D J_im = 0 (zero net current)
            J_im = _coil_solve(-f_back_im, False)

            gf_J_re.vec.FV().NumPy()[:] = J_re
            gf_J_im.vec.FV().NumPy()[:] = J_im

            # --- Inductance from total magnetic energy ---
            # Self contribution from BOTH real and imag back-reacted current.
            # Mutual contribution from per-DOF back-reaction inner product.
            L_self_now = MU_0 * (_sl_energy(J_re) + _sl_energy(J_im))
            mutual = float(f_back_re @ J_re + f_back_im @ J_im)
            L_total = L_self_now + mutual
            Delta_L = L_total - L_air

            dL_rel = abs(L_total - L_prev) / max(abs(L_prev), 1e-30)
            if verbose:
                print(f"  iter {iteration}: L_self={L_self_now*1e9:.3f}nH "
                      f"mutual={mutual*1e9:+.3f}nH "
                      f"L_total={L_total*1e9:.3f}nH "
                      f"DeltaL={Delta_L*1e9:+.3f}nH dL={dL_rel:.3e}")

            L_prev = L_total
            if dL_rel < tol and iteration > 0:
                break

        P_density = wp_result['P_density']
        P_total = P_density * wp_result['area']
        H_t_rms = wp_result['H_t_rms']

        # Genus-1 loop DOF on the CONVERGED state (see the loop_dof arg
        # docstring): solve the shorted-turn current alpha against the
        # LAST-ITERATION coil current (the same current that produced
        # phi_inc_cplx and wp_result, so the frozen(alpha=0) sub-solve
        # must reproduce wp_result exactly -- the operator cross-check).
        loop_meta = {}
        if loop_dof:
            import time as _time
            from radia.bem_loop_extension import solve_loop_extended
            from radia.bem_sibc_solver import A_from_surface_J

            coil_Jc = coil_J_re_arr + 1j * coil_J_im_arr

            def _A_inc_fn(points):
                A_re = A_from_surface_J(points, coil_c, coil_a,
                                        np.real(coil_Jc))
                A_im = A_from_surface_J(points, coil_c, coil_a,
                                        np.imag(coil_Jc))
                return A_re + 1j * A_im

            t0 = _time.perf_counter()
            loop_out = solve_loop_extended(
                self.wp_solver, phi_inc_cplx, Z_s, omega, _A_inc_fn)
            t_loop = _time.perf_counter() - t0
            frz_rel = (abs(loop_out["P_frozen"] - P_total)
                       / max(abs(P_total), 1e-30))
            if frz_rel > 1e-3:
                raise RuntimeError(
                    f"loop-DOF frozen sub-solve disagrees with the "
                    f"converged plain BIE solve by {frz_rel:.2e} "
                    f"(P {loop_out['P_frozen']:.4e} vs {P_total:.4e} W) "
                    f"-- operator mismatch, refusing to report "
                    f"loop-extended numbers.")
            P_frozen = float(loop_out["P_frozen"])
            P_total = float(loop_out["P_total"])
            H_t_rms = float(loop_out["H_t_rms"])
            loop_meta = {
                'wp_loop_alpha': complex(loop_out["alpha"]),
                'wp_loop_theta_jump': float(loop_out["theta_jump"]),
                'wp_loop_cut_n_vertices': int(loop_out["cut_n_vertices"]),
                # frozen (alpha = 0) = the no-mode physics; the ratio
                # exposes the Lenz-screening effect of the shorted turn
                # (scale-invariant, so current rescaling never touches it)
                'wp_loop_P_frozen': P_frozen,
                'wp_loop_H_t_frozen': float(loop_out["Ht_frozen"]),
                'wp_loop_screening_ratio': P_total / max(P_frozen, 1e-300),
                't_loop_dof_s': float(t_loop),
            }

        # Z_s passthrough: complex scalar in the legacy path, but a
        # ndarray in the per-node path. The caller treats it as opaque.
        if isinstance(Z_s, np.ndarray):
            Z_s_out = complex(np.mean(Z_s))
        else:
            Z_s_out = complex(Z_s)

        return {
            'L_air': float(L_air),
            'L_self': float(L_self_now),
            'L_total': float(L_total),
            'Delta_L': float(Delta_L),
            'P_total': float(P_total),
            'H_t_rms': float(H_t_rms),
            'iterations': iteration + 1,
            'n_J_coil': self.n_J,
            'n_phi_wp': self.wp_solver.ndof,
            'J_coil_re': J_re,
            'J_coil_im': J_im,
            # Workpiece-side per-panel data needed for GMSH
            # visualization (J vector + sigma|J|^2/2 heating density).
            # ``wp_J_re/im_arr`` are arrays of shape (M, 3) of the
            # scattered surface current per BND element on the
            # workpiece mesh, ``wp_c`` are panel centroids and
            # ``wp_a`` are panel areas.
            'wp_c': wp_c,
            'wp_a': wp_a,
            'wp_J_re': wp_J_re_arr,
            'wp_J_im': wp_J_im_arr,
            'Z_s': Z_s_out,
            **loop_meta,
        }

    def _extract_wp_J(self, phi_vec_complex, phi_inc_complex):
        """Extract scattered surface current ``J_wp = n x H_scat`` from the
        SIBC scalar BIE solution.  Thin wrapper over the shared free
        function ``extract_scattered_wp_J`` (behaviour unchanged).
        """
        return extract_scattered_wp_J(
            self.mesh_wp, self.wp_solver.fes,
            phi_vec_complex, phi_inc_complex)
