"""radia.moment_galerkin._solve -- the production symmetric moment-Galerkin MMMM demag solve.

Solves the SPD +N physical demag system  ( (1/chi) M_mass + B^T G B ) m = M_mass H_ext  by the EXISTING
C++ charge-Gram Krylov kernel (the same kernel HDiv-VIM ships): a pure-dipole body uses the mass-Riesz CG
(`solve_linear_material_mass_riesz`, nearly mu_r-flat); a body with ANY quad DOF (wedge/pyramid/hex quad modes)
uses the diagonal-Jacobi CG (`solve_linear_material_auto_prec`, whose C++-computed Jacobi diagonal includes the
demag self-energy, so it stays well-defined even when a quad MASS block is near-singular on a near-symmetric
element).  m = the raw per-element moment amplitudes laid out at the variable per-element DOF offsets dof_off
(3 dipole + (nFace-4) quad per element).  For a curl-free (gradient / uniform) source the demag near-null
(divergence-free circulating) modes are RHS-orthogonal so the +N CG is effectively loop-free; loop-EXCITING
sources (azimuthal / transformer drive) route through HDiv-VIM (radia.vim).

MIXED ELEMENTS: tetrahedron (3 DOF) / pyramid (4) / wedge (4) / hexahedron (5) -- the old yano-type element zoo.
The dipole magnetization M of each element is the first 3 DOF (m[dof_off[e]:dof_off[e]+3]); the quad amplitudes
are the remaining ndof_elem[e]-3.

Validated: cube self-consistency M_z = chi/(1+d chi) H0 (machine precision); demag factor 1/3 reproduced by a
cube meshed as hex / 6 tets / 2 wedges / 6 pyramids (internal charges cancel); N = B^T G B symmetric to 1e-16
on every element type + mixed; the field-from-sigma reconstruction matches rad.Fld to ~1e-13.
"""
import numpy as np
import scipy.sparse as sp
from ._assemble import assemble_moment_system

MU0 = 4e-7 * np.pi


def _chi_of(mu_r, chi):
    if (mu_r is None) == (chi is None):
        raise ValueError("moment_galerkin: pass EXACTLY one of mu_r or chi")
    chi_val = float(chi) if chi is not None else float(mu_r) - 1.0
    if not np.isfinite(chi_val) or chi_val <= 0.0:
        raise ValueError("moment_galerkin: require mu_r > 1 or chi > 0")
    return chi_val


def _pad_field(H_ext, dof_off, ndof_total, n):
    """Per-DOF applied field h: [Hx,Hy,Hz] on each element's 3 dipole DOF, 0 on the quad DOF (a uniform /
    per-element-constant field does not drive the quad rows; the quad RHS INT M_q.H_ext vanishes for a
    within-element-uniform H)."""
    H = np.asarray(H_ext, float)
    per_elem = np.tile(H, (n, 1)) if H.size == 3 else H.reshape(n, 3)
    if per_elem.shape != (n, 3):
        raise ValueError("moment_galerkin: H_ext must be (3,) or (n_elem,3)")
    h = np.zeros(ndof_total)
    for e in range(n):
        h[dof_off[e]:dof_off[e] + 3] = per_elem[e]
    return h


def solve_raw(sys, H_ext, chi, *, tol=1e-9, maxit=2000):
    """Solve for the raw moment amplitudes m (length ndof_total).  Returns (m, iters)."""
    G, B, M_mass = sys["G"], sys["B"], sys["M_mass"]
    n, ndof_total, dof_off, ndof_elem = sys["n_elem"], sys["ndof_total"], sys["dof_off"], sys["ndof_elem"]
    chi = float(chi)
    if not np.isfinite(chi) or chi <= 0.0:
        raise ValueError("moment_galerkin: chi must be positive")
    h = _pad_field(H_ext, dof_off, ndof_total, n)
    rhs = np.asarray(M_mass @ h).ravel()
    inv_chi = 1.0 / chi
    Bc = B.tocsr()
    Mc = sp.coo_matrix(M_mass)
    args = (list(map(int, Bc.indptr)), list(map(int, Bc.indices)), list(map(float, Bc.data)),
            int(B.shape[1]), list(map(int, Mc.row)), list(map(int, Mc.col)), list(map(float, Mc.data)),
            inv_chi, list(map(float, rhs)), tol, int(maxit))
    # Pure dipole -> mass-Riesz CG (mu_r-flat).  Any quad DOF -> diagonal-Jacobi auto_prec (robust to the
    # near-singular quad mass block; its Jacobi diagonal includes the demag self-energy N_diag).
    has_quad = bool(np.any(np.asarray(ndof_elem) > 3))
    res = (G.solve_linear_material_auto_prec(*args) if has_quad
           else G.solve_linear_material_mass_riesz(*args))
    iters = int(res["iters"])
    if iters >= int(maxit):                           # fail-loud (No-Fallbacks): never return a non-converged m
        raise RuntimeError("moment_galerkin: CG did NOT converge in %d iters (n_elem=%d, ndof_total=%d). "
                           "Tighten the ACA (leaf up / eta down) or raise maxit." % (maxit, n, ndof_total))
    return np.asarray(res["m"], float), iters


def solve_assembled(sys, H_ext, chi, *, tol=1e-9, maxit=2000):
    """Solve a pre-assembled `sys`; return (M_dipole (n_elem,3), iters).  M_dipole = the constant (average)
    magnetization per element -- the first 3 (dipole) DOF (the quad amplitudes, if any, stay in solve_raw's m)."""
    m, iters = solve_raw(sys, H_ext, chi, tol=tol, maxit=maxit)
    dof_off, n = sys["dof_off"], sys["n_elem"]
    M = np.array([m[dof_off[e]:dof_off[e] + 3] for e in range(n)])
    return M, iters


def demag_factor(sys, kdir=2):
    """Operator demag factor d = <c, G c>/<m, M_mass m> for a uniform M = e_kdir (cube -> 1/3)."""
    G, B, M_mass, dof_off = sys["G"], sys["B"], sys["M_mass"], sys["dof_off"]
    m = np.zeros(B.shape[1]); m[np.asarray(dof_off) + kdir] = 1.0
    c = np.asarray(B @ m)
    Gc = np.asarray(G.matvec(c.tolist()), float)
    return float(c @ Gc) / float(m @ (M_mass @ m))


def reconstruct_field(sys, m, probes):
    """External magnetic flux density B (Tesla) at probe point(s), reconstructed from the solved moment
    amplitudes via the face surface charges sigma = B m (the SAME charges the demag operator uses):
        B(p) = mu0/(4 pi) * sum_tris sigma_tri * H_unit_charged_triangle(p)
    -- the field-from-sigma reconstruction, anchored to rad.Fld to ~1e-13.  PROBES MUST BE EXTERNAL to the
    iron (inside, add mu0*M); avoid points exactly ON a face plane (the analytic triangle field is singular
    in-plane).  probes : (3,) or (P,3).  Returns (3,) or (P,3)."""
    from radia.vim._field import flat_triangle_charge_field
    B, all_tris = sys["B"], sys["all_tris"]
    sigma = np.asarray(B @ np.asarray(m, float).ravel()).ravel()
    pts = np.atleast_2d(np.asarray(probes, float))
    out = np.zeros((pts.shape[0], 3))
    for i, p in enumerate(pts):
        H = np.zeros(3)
        for t, T in enumerate(all_tris):
            if sigma[t] != 0.0:
                H += sigma[t] * flat_triangle_charge_field(T, p)
        out[i] = MU0 * (1.0 / (4.0 * np.pi)) * H
    return out[0] if np.asarray(probes, float).ndim == 1 else out


def moment_galerkin_demag_solve(elements, *, mu_r=None, chi=None, H_ext=(0.0, 0.0, 0.0), quad=False,
                                tol=1e-9, maxit=2000, eps=1e-9, leaf=40, eta=0.5, near_factor=2.0, far_quad=4):
    """Linear isotropic soft-iron demag on a MIXED polyhedral body via the symmetric moment-Galerkin MMMM.

    elements : list of (nv,3) vertex arrays; nv = 4 (tet) / 5 (pyramid) / 6 (wedge) / 8 (hex), radia corner order.
    quad : False -> 3-DOF dipole/element; True -> + (nFace-4) quad DOF (tet 0, wedge/pyramid 1, hex 2) for
        higher per-element order under skew/gradient loads.
    mu_r OR chi : relative permeability (>1) or susceptibility chi (pass exactly one).
    H_ext : (3,) uniform or (n_elem,3) per-element applied field (A/m).
    near_factor, far_quad : charge-Gram NEAR/FAR fast-build split (defaults 2 / 4 = the precision-preserving fast
        build; pass near_factor=1e30 for the all-analytic Gram).  See assemble_moment_system.

    Returns dict(M=(n_elem,3) dipole magnetization (A/m), m=raw amplitudes, quad_amps ((n_elem,nq) if uniform
                 else a list of per-element quad arrays, or None), iters, demag_factor, n_elem, n_hex, sys).
                 Use reconstruct_field(result['sys'], result['m'], probes) for the external B field."""
    chi = _chi_of(mu_r, chi)
    sys = assemble_moment_system(elements, quad=quad, eps=eps, leaf=leaf, eta=eta,
                                 near_factor=near_factor, far_quad=far_quad)
    m, iters = solve_raw(sys, H_ext, chi, tol=tol, maxit=maxit)
    dof_off, ndof_elem, n = sys["dof_off"], sys["ndof_elem"], sys["n_elem"]
    M = np.array([m[dof_off[e]:dof_off[e] + 3] for e in range(n)])
    quad_amps = None
    if quad:
        quad_amps = [m[dof_off[e] + 3:dof_off[e] + ndof_elem[e]] for e in range(n)]
        if len({len(q) for q in quad_amps}) == 1:     # uniform element type -> rectangular (n_elem, nq)
            quad_amps = np.array(quad_amps)
    return {"M": M, "m": m, "quad_amps": quad_amps, "iters": iters,
            "demag_factor": demag_factor(sys, 2), "n_elem": n, "n_hex": n, "sys": sys}
