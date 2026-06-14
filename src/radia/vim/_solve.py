"""radia.vim._solve -- consolidated production demag-solve entry (productionization M1).

`hdiv_demag_solve(mesh, mu_r, H_ext, ...)` is the first shippable HDiv-type VIM Radia API: the LINEAR
soft-iron applied-field demag solve, the candidate replacement for the yano-type MSC hex/wedge
soft-iron demag in `rad.Solve`.

## Formulation (verified-first, 2026-06-15)
The physical linear-material system is the **+N** system

    ((1/chi) M_mass + N) m = M_mass h_ext ,   N = B^T G B ,   chi = mu_r - 1 ,

with `h_ext` the applied field H_ext L2-projected onto RT0 (HDiv order 0).  Verified against the
analytic uniform-sphere magnetization `M = chi/(1 + chi D) H` (D = demag factor ~ 1/3) to <= 1.4e-4
across mu_r 10..1000.  (The **-N** system is the loop-field-null mu_r-independent operator but is
NON-physical for the magnetization -- wrong sign -- so the +N system is the one solved here.)

## Preconditioner -- the FULL HDiv mass inverse (mesh + mu_r robust)
The +N system is solved by GMRES preconditioned with `M_mass^{-1}` (a one-time sparse LU of the local
RT0 HDiv mass).  This is the doc's proven mesh-robust path: `M_mass^{-1}` DEFLATES the de-Rham loops
(ker B), so the iteration count stays bounded as mu_r grows AND as the mesh refines (~tens of iters).
A Jacobi-diagonal preconditioner is NOT mesh-robust (the +N CG count blows up past a few thousand
faces -- measured), which is why the C++ `solve_linear_material` Jacobi kernel is not used here.
GMRES (not CG/MINRES) absorbs the residual ACA asymmetry of the H-matvec demag apply.

KELVIN-LESS: the 1/r charge Gram IS the open boundary (a volume integral method like MMM/MSC); only
the iron is meshed -- no air box / Kelvin needed.

## Scope (M1, this increment)
- LINEAR, single isotropic soft-iron region (scalar `mu_r`).
- Arbitrary applied field `H_ext` (any NGSolve CoefficientFunction: `rad.RadiaField(coil,'h')`,
  analytic, ...).
- Returns per-element M (n_el x 3), the volume-average M, iters, the demag factor.

NONLINEAR (BH table) + per-region materials + the M0 parity gate are the remaining productionization
steps (docs/hdiv_vim/PRODUCTIONIZATION.md).  Until they land, the yano-type MSC stays the `rad.Solve`
default demag backend (`radia.set_demag_backend`); this entry does not touch it.

Per CLAUDE.md "TaskManager Wrap Policy: Caller Wraps, Helper Does NOT" -- this library helper does NOT
open a TaskManager; the caller (calc_*.py / example / test) wraps the call in `with ng.TaskManager():`.
"""
import inspect

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import ngsolve as ng

import radia._radia_pybind as _rp
from . import _core as _tet

# scipy renamed the Krylov tolerance kwarg 'tol' -> 'rtol' (scipy >= 1.12); detect once.
_GMRES_TOL = "rtol" if "rtol" in inspect.signature(spla.gmres).parameters else "tol"


def hdiv_demag_solve(mesh, mu_r, H_ext, *, scalable=True, gram_eps=1e-10, leaf=32, eta=2.0,
                     tol=1e-8, maxit=4000, gmres_restart=400):
    """Linear soft-iron demag solve via the HDiv-type VIM (the +N physical material system).

    Parameters
    ----------
    mesh : ngsolve.Mesh
        Tet mesh of the soft-iron body ONLY (KELVIN-LESS: the 1/r Gram is the open boundary).
    mu_r : float
        Relative permeability (> 1); chi = mu_r - 1, isotropic linear.
    H_ext : ngsolve.CoefficientFunction
        Applied (incident) field H_ext over the mesh (A/m).
    scalable : bool
        True (default): the demag N is applied by the C++ charge-Gram H-matvec (O(N log N), no dense
        N^2).  False: dense exact analytic Gram N (the golden cross-check reference, small meshes).
        Both use the SAME +N system + M_mass^{-1}-preconditioned GMRES.
    gram_eps, leaf, eta : float, int, float
        H-matrix (ACA) parameters for `_ChargeGramHMatrix` (scalable path).  gram_eps=1e-10 keeps the
        H-matvec symmetric enough that GMRES converges in ~tens of iters.
    tol, maxit, gmres_restart : float, int, int
        GMRES relative tolerance, total inner-iteration budget, and restart length.  RAISES if not
        converged (CLAUDE.md No-Fallbacks).

    Returns
    -------
    dict
        M        : (n_el, 3) per-element magnetization (A/m).
        M_avg    : (3,) volume-average magnetization (A/m).
        iters    : GMRES iteration count.
        demag    : demag factor D (Rayleigh quotient; ~1/3 for a sphere).
        ndof, n_el, n_charge : sizes.
    """
    chi = float(mu_r) - 1.0
    if chi <= 0.0:
        raise ValueError("hdiv_demag_solve: mu_r must be > 1 for a soft-iron demag solve (got %r)"
                         % (mu_r,))
    inv_chi = 1.0 / chi

    if scalable:
        d = _tet.build_demag(mesh, skip_dense_gram=True)
        Mm = d["M_mass"]                              # sparse CSR (RT0 HDiv mass is local)
        B = d["B_csr"]
        H = _rp._ChargeGramHMatrix(cell_verts=list(d["cell_verts"]), face_verts=list(d["face_verts"]),
                                   n_el=int(d["n_el"]), eps=gram_eps, leaf=leaf, eta=eta)

        def N_apply(v):                               # demag N v = B^T G (B v) via the H-matvec
            q = B @ np.asarray(v, float)
            return B.T @ np.asarray(H.matvec(q.tolist()), float)
    else:
        d = _tet.build_demag(mesh, analytic_gram=True)
        Mm = d["M_mass"]                              # dense (reference path)
        N = d["N"]

        def N_apply(v):                               # exact dense demag apply
            return N @ np.asarray(v, float)

    n_face = int(d["ndof"]); n_el = int(d["n_el"])

    # source projection: H_ext L2-projected onto RT0 -> h_ext DOFs; physical RHS = M_mass h_ext
    fes = ng.HDiv(mesh, order=0)
    gfH = ng.GridFunction(fes); gfH.Set(H_ext)
    h_ext = gfH.vec.FV().NumPy().copy()
    rhs = np.asarray(Mm @ h_ext).ravel()

    # M_mass^{-1} preconditioner (one-time sparse LU; deflates the de-Rham loops -> mesh/mu_r robust)
    Mfac = spla.splu(sp.csc_matrix(Mm))
    Mprec = spla.LinearOperator((n_face, n_face), matvec=lambda v: Mfac.solve(np.asarray(v, float)))

    def A_apply(v):                                   # the +N system ((1/chi) M_mass + N)
        v = np.asarray(v, float)
        return inv_chi * np.asarray(Mm @ v).ravel() + N_apply(v)

    A = spla.LinearOperator((n_face, n_face), matvec=A_apply)
    it = {"n": 0}
    cycles = max(2, int(np.ceil(maxit / gmres_restart)))
    m, info = spla.gmres(A, rhs, M=Mprec, restart=int(gmres_restart), maxiter=cycles,
                         callback=lambda _xk: it.__setitem__("n", it["n"] + 1),
                         callback_type="pr_norm", **{_GMRES_TOL: tol})
    if info != 0:
        raise RuntimeError("hdiv_demag_solve: +N GMRES did not converge (info=%d, mu_r=%g, n_face=%d, "
                           "iters=%d)" % (info, mu_r, n_face, it["n"]))
    iters = it["n"]

    # per-element M: L2-project the RT0 magnetization onto VectorL2(0).  VectorL2 order-0 DOFs are
    # COMPONENT-major -> reshape(3, n_el).T gives (n_el, 3).
    gfM = ng.GridFunction(fes); gfM.vec.FV().NumPy()[:] = m
    fesM = ng.VectorL2(mesh, order=0); gfMc = ng.GridFunction(fesM); gfMc.Set(gfM)
    M_el = gfMc.vec.FV().NumPy().reshape(3, n_el).T.copy()
    vol = ng.Integrate(ng.CoefficientFunction(1.0), mesh)
    M_avg = np.array([ng.Integrate(gfM[i], mesh) for i in range(3)]) / vol

    # demag factor (Rayleigh quotient of the unit-M projection)
    mu = d["m_unit"]
    D = float((mu @ N_apply(mu)) / float(mu @ np.asarray(Mm @ mu).ravel()))

    return dict(M=M_el, M_avg=M_avg, iters=iters, demag=D,
                ndof=n_face, n_el=n_el, n_charge=int(d["n_charge"]))
