"""radia.moment_galerkin._solve -- the production symmetric moment-Galerkin MMMM demag solve.

Solves the SPD +N physical demag system

    ( (1/chi) M_mass + B^T G B ) m = M_mass H_ext

by the EXISTING C++ charge-Gram mass-Riesz CG (`_ChargeGramHMatrix.solve_linear_material_mass_riesz`) -- the
same kernel HDiv-VIM ships.  m = per-hex dipole moments (3/hex).  For a curl-free (gradient / uniform) source
the demag near-null (divergence-free circulating) modes are RHS-orthogonal, so this +N CG is effectively
loop-free (the iteration count plateaus in mu_r).  Loop-EXCITING sources (azimuthal / transformer drive)
should be routed through HDiv-VIM (radia.vim), whose RT-flux 'mixed form' quotients the curl space.

Validated (de-risk increment 2): single-cube self-consistency M_z = chi/(1+d*chi)*H0 to machine precision;
a 2-hex bar matches rad.Solve's iron stray field to 0.29% at a far probe (dipole vs full 6-DOF).
"""
import numpy as np
import scipy.sparse as sp
from ._assemble import assemble_moment_system


def _chi_of(mu_r, chi):
    if (mu_r is None) == (chi is None):
        raise ValueError("moment_galerkin: pass EXACTLY one of mu_r or chi")
    chi_val = float(chi) if chi is not None else float(mu_r) - 1.0
    if not np.isfinite(chi_val) or chi_val <= 0.0:
        raise ValueError("moment_galerkin: require mu_r > 1 or chi > 0")
    return chi_val


def solve_assembled(sys, H_ext, chi, *, tol=1e-9, maxit=2000):
    """Solve the +N physical system for a pre-assembled `sys` (dict from assemble_moment_system).

    H_ext : (3,) uniform field OR (n_hex,3) per-hex field (A/m).  Returns (M (n_hex,3), iters)."""
    G, B, M_mass, n = sys["G"], sys["B"], sys["M_mass"], sys["n_hex"]
    chi = float(chi)
    if not np.isfinite(chi) or chi <= 0.0:
        raise ValueError("moment_galerkin: chi must be positive")
    H = np.asarray(H_ext, float)
    h = np.tile(H.ravel(), n) if H.size == 3 else H.ravel()
    if h.size != 3 * n:
        raise ValueError("moment_galerkin: H_ext must be (3,) or (n_hex,3)")
    rhs = np.asarray(M_mass @ h).ravel()
    inv_chi = 1.0 / chi
    Bc = B.tocsr()
    Mc = sp.coo_matrix(M_mass)
    res = G.solve_linear_material_mass_riesz(
        list(map(int, Bc.indptr)), list(map(int, Bc.indices)), list(map(float, Bc.data)),
        int(B.shape[1]), list(map(int, Mc.row)), list(map(int, Mc.col)), list(map(float, Mc.data)),
        inv_chi, list(map(float, rhs)), tol, int(maxit))
    iters = int(res["iters"])
    if iters >= int(maxit):                       # fail-loud (No-Fallbacks): never return a non-converged M
        raise RuntimeError("moment_galerkin_demag_solve: mass-Riesz CG did NOT converge in %d iters "
                           "(n_hex=%d). Tighten the ACA (leaf up / eta down) or raise maxit." % (maxit, n))
    return np.asarray(res["m"], float).reshape(n, 3), iters


def demag_factor(sys, kdir=2):
    """Operator demag factor d = <c, G c>/<m, M_mass m> for a uniform M = e_kdir (cube -> 1/3)."""
    G, B, M_mass = sys["G"], sys["B"], sys["M_mass"]
    m = np.zeros(B.shape[1]); m[kdir::3] = 1.0
    c = np.asarray(B @ m)
    Gc = np.asarray(G.matvec(c.tolist()), float)
    return float(c @ Gc) / float(m @ (M_mass @ m))


def moment_galerkin_demag_solve(hexes, *, mu_r=None, chi=None, H_ext=(0.0, 0.0, 0.0),
                                tol=1e-9, maxit=2000, eps=1e-9, leaf=40, eta=0.5):
    """Linear isotropic soft-iron demag on a hexahedral body via the symmetric moment-Galerkin MMMM.

    Parameters
    ----------
    hexes : sequence of (8,3) hex vertex arrays (meters, Radia/HEX_FACES order).
    mu_r OR chi : relative permeability (>1) or susceptibility chi = mu_r - 1 (pass exactly one).
    H_ext : (3,) uniform or (n_hex,3) per-hex applied field (A/m).
    tol, maxit : CG tolerance / cap.   eps, leaf, eta : charge-Gram ACA (validated defaults).

    Returns
    -------
    dict(M=(n_hex,3) magnetization (A/m), iters=int, demag_factor=float, n_hex=int)
    """
    chi = _chi_of(mu_r, chi)
    sys = assemble_moment_system(hexes, eps=eps, leaf=leaf, eta=eta)
    M, iters = solve_assembled(sys, H_ext, chi, tol=tol, maxit=maxit)
    return {"M": M, "iters": iters, "demag_factor": demag_factor(sys, 2), "n_hex": sys["n_hex"]}
