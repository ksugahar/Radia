"""radia.moment_galerkin._solve -- the production symmetric moment-Galerkin MMMM demag solve.

Solves the SPD +N physical demag system  ( (1/chi) M_mass + B^T G B ) m = M_mass H_ext  by the EXISTING
C++ charge-Gram Krylov kernel (the same kernel HDiv-VIM ships): the 3-DOF dipole path uses the mass-Riesz CG
(`solve_linear_material_mass_riesz`, nearly mu_r-flat); the 5-DOF (dipole+quad) path uses the diagonal-Jacobi
CG (`solve_linear_material_auto_prec`, whose C++-computed Jacobi diagonal includes the demag self-energy, so
it stays well-defined even when the quad MASS block is near-singular on a cube).  m = per-hex moment
amplitudes (3 or 5/hex).  For a curl-free (gradient / uniform) source the demag near-null (divergence-free
circulating) modes are RHS-orthogonal so the +N CG is effectively loop-free (iterations plateau in mu_r);
loop-EXCITING sources (azimuthal / transformer drive) route through HDiv-VIM (radia.vim).

Validated: cube self-consistency M_z = chi/(1+d chi) H0 (machine precision); cube quad amplitudes ~0
(symmetry); a 3-hex bar 5-DOF stray field is strictly closer to a fine-mesh rad.Solve(yano) truth than
dipole-only at every near probe under an oblique field; the field-from-sigma reconstruction matches rad.Fld
to ~1e-13.
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


def _pad_field(H_ext, n, ndof_per):
    """Per-DOF applied field h: [Hx,Hy,Hz] on the 3 dipole DOF, 0 on the quad DOF (a uniform/per-hex-constant
    field does not drive the quad rows; the quad RHS INT M_q.H_ext vanishes for a within-element-uniform H)."""
    H = np.asarray(H_ext, float)
    per_hex = np.tile(H, (n, 1)) if H.size == 3 else H.reshape(n, 3)
    if per_hex.shape != (n, 3):
        raise ValueError("moment_galerkin: H_ext must be (3,) or (n_hex,3)")
    h = np.zeros(ndof_per * n)
    for e in range(n):
        h[ndof_per * e:ndof_per * e + 3] = per_hex[e]
    return h


def solve_raw(sys, H_ext, chi, *, tol=1e-9, maxit=2000):
    """Solve for the raw moment amplitudes m (length ndof_per*n_hex).  Returns (m, iters)."""
    G, B, M_mass, n, ndof_per = sys["G"], sys["B"], sys["M_mass"], sys["n_hex"], sys["ndof_per"]
    chi = float(chi)
    if not np.isfinite(chi) or chi <= 0.0:
        raise ValueError("moment_galerkin: chi must be positive")
    h = _pad_field(H_ext, n, ndof_per)
    rhs = np.asarray(M_mass @ h).ravel()
    inv_chi = 1.0 / chi
    Bc = B.tocsr()
    Mc = sp.coo_matrix(M_mass)
    args = (list(map(int, Bc.indptr)), list(map(int, Bc.indices)), list(map(float, Bc.data)),
            int(B.shape[1]), list(map(int, Mc.row)), list(map(int, Mc.col)), list(map(float, Mc.data)),
            inv_chi, list(map(float, rhs)), tol, int(maxit))
    # 3-DOF: mass-Riesz CG (mu_r-flat).  5-DOF: diagonal-Jacobi auto_prec (robust to the near-singular quad
    # mass block; its Jacobi diagonal includes the demag self-energy N_diag).
    res = (G.solve_linear_material_auto_prec(*args) if ndof_per == 5
           else G.solve_linear_material_mass_riesz(*args))
    iters = int(res["iters"])
    if iters >= int(maxit):                           # fail-loud (No-Fallbacks): never return a non-converged m
        raise RuntimeError("moment_galerkin: CG did NOT converge in %d iters (n_hex=%d, ndof_per=%d). "
                           "Tighten the ACA (leaf up / eta down) or raise maxit." % (maxit, n, ndof_per))
    return np.asarray(res["m"], float), iters


def solve_assembled(sys, H_ext, chi, *, tol=1e-9, maxit=2000):
    """Solve a pre-assembled `sys`; return (M_dipole (n_hex,3), iters).  M_dipole = the constant (average)
    magnetization per hex -- the dipole amplitudes (the 2 quad amplitudes, if present, are in solve_raw's m)."""
    m, iters = solve_raw(sys, H_ext, chi, tol=tol, maxit=maxit)
    ndof_per, n = sys["ndof_per"], sys["n_hex"]
    M = np.array([m[ndof_per * e:ndof_per * e + 3] for e in range(n)])
    return M, iters


def demag_factor(sys, kdir=2):
    """Operator demag factor d = <c, G c>/<m, M_mass m> for a uniform M = e_kdir (cube -> 1/3)."""
    G, B, M_mass, ndof_per = sys["G"], sys["B"], sys["M_mass"], sys["ndof_per"]
    m = np.zeros(B.shape[1]); m[kdir::ndof_per] = 1.0
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


def moment_galerkin_demag_solve(hexes, *, mu_r=None, chi=None, H_ext=(0.0, 0.0, 0.0), quad=False,
                                tol=1e-9, maxit=2000, eps=1e-9, leaf=40, eta=0.5, near_factor=2.0, far_quad=4):
    """Linear isotropic soft-iron demag on a hexahedral body via the symmetric moment-Galerkin MMMM.

    quad : False -> 3-DOF dipole (constant M/hex); True -> 5-DOF (3 dipole + 2 quad, higher per-element order
        for skew/gradient loads).
    mu_r OR chi : relative permeability (>1) or susceptibility chi (pass exactly one).
    H_ext : (3,) uniform or (n_hex,3) per-hex applied field (A/m).
    near_factor, far_quad : charge-Gram NEAR/FAR fast-build split (defaults 2 / 4 = the precision-preserving fast
        build; pass near_factor=1e30 for the all-analytic Gram).  See assemble_moment_system.

    Returns dict(M=(n_hex,3) dipole magnetization (A/m), m=raw amplitudes, quad_amps=(n_hex,2) or None,
                 iters, demag_factor, n_hex, sys).  Use reconstruct_field(result['sys'], result['m'], probes)
                 for the external B field."""
    chi = _chi_of(mu_r, chi)
    sys = assemble_moment_system(hexes, quad=quad, eps=eps, leaf=leaf, eta=eta,
                                 near_factor=near_factor, far_quad=far_quad)
    m, iters = solve_raw(sys, H_ext, chi, tol=tol, maxit=maxit)
    ndof_per, n = sys["ndof_per"], sys["n_hex"]
    M = np.array([m[ndof_per * e:ndof_per * e + 3] for e in range(n)])
    quad_amps = np.array([m[ndof_per * e + 3:ndof_per * e + 5] for e in range(n)]) if quad else None
    return {"M": M, "m": m, "quad_amps": quad_amps, "iters": iters,
            "demag_factor": demag_factor(sys, 2), "n_hex": n, "sys": sys}
