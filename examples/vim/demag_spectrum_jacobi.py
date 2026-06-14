"""
demag_spectrum_jacobi.py -- is the high-mu breakdown the LOOPS (deflatable -> Jacobi)
or a weak-demag CONTINUUM (needs H-ILU)?  Measured on a COMPACT body with Radia's EXACT
ObjRecMag field (no quadrature error), so the spectrum is trustworthy.

System (linear soft magnet, chi = mu_r - 1):  ( (1/chi) I - N ) M = H_ext,
  N = the demag coupling (H_demag at block i from unit M in block j; self term N_ii is
  the exact in-block demag field).  A = (1/chi) I - N.  Eigenvalues lambda = 1/chi - mu_i,
  mu_i = eigenvalues of N (mu=0 = field-null "loops"; mu>0 = charge/star modes).

The Jacobi question (docs/loop_star_breakdown.md "Does it have to be H-ILU?"):
  - if the smallest CHARGE |mu_i| >> 1/chi (a spectral GAP above the loops), then deflating
    the few loop modes bounds cond -> a light (Jacobi) preconditioner suffices for mu_r<=1e4;
  - if mu_i forms a near-continuum down to 0 (Neumann-Poincare accumulation; thin/elongated
    bodies; very fine meshes), no fixed deflation bounds it -> H-ILU / spectral coarse space.

We assemble N exactly, then report (a) the spectrum near 0, (b) cond(A) vs mu_r, (c) cond
after deflating the k near-null modes, (d) Jacobi-preconditioned iteration counts.
"""
import json, os
import numpy as np
import radia as rad

HERE = os.path.dirname(os.path.abspath(__file__))
res = {}

def build_grid_demag(n, a=0.01, distort=0.0, seed=1):
    """n x n x n grid of unit cubes (side a); return N (3n^3 x 3n^3) demag coupling,
    exact via Radia ObjRecMag.  distort>0 jitters block centers (non-uniform geometry)."""
    rng = np.random.default_rng(seed)
    centers = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                c = np.array([i * a, j * a, k * a], float)
                if distort > 0:
                    c = c + distort * a * (rng.random(3) - 0.5)
                centers.append(c)
    centers = np.array(centers)
    nb = len(centers)
    N = np.zeros((3 * nb, 3 * nb))
    dims = [a, a, a]
    for jb in range(nb):
        for b in range(3):                       # source block jb, M = e_b (A/m units = 1)
            rad.UtiDelAll()
            M = [0.0, 0.0, 0.0]; M[b] = 1.0
            src = rad.ObjRecMag(centers[jb].tolist(), dims, M)
            # H field of this single magnetized block at every block center (exact)
            H = np.array(rad.Fld(src, 'h', centers.tolist())).reshape(nb, 3)
            N[:, 3 * jb + b] = H.reshape(-1)     # column: response at all i from (jb,b)
    rad.UtiDelAll()
    return N, nb

def jacobi_pcg_iters(A, b, tol=1e-8, maxit=2000):
    """diagonal (Jacobi) preconditioned GMRES-free test via CG on the normal eqs is
    unstable for non-symmetric A; use a simple right-preconditioned GMRES iteration count."""
    from scipy.sparse.linalg import gmres, LinearOperator
    d = np.diag(A).copy(); d[np.abs(d) < 1e-30] = 1.0
    Minv = LinearOperator(A.shape, matvec=lambda x: x / d)
    it = {"n": 0}
    def cb(xk): it["n"] += 1
    gmres(A, b, M=Minv, rtol=tol, maxiter=maxit, callback=cb, callback_type="pr_norm")
    return it["n"]

for distort, tag in [(0.0, "regular"), (0.3, "distorted")]:
    n = 3
    N, nb = build_grid_demag(n, distort=distort)
    ndof = 3 * nb
    muN = np.linalg.eigvals(N).real
    muN_sorted = np.sort(np.abs(muN))
    # loops = |mu| ~ 0; smallest CHARGE |mu| = first |mu| above a loop tolerance
    loop_tol = 1e-3 * muN_sorted.max()
    n_loop = int(np.sum(muN_sorted < loop_tol))
    smallest_charge_mu = muN_sorted[n_loop] if n_loop < len(muN_sorted) else float("nan")
    print(f"\n=== {tag} {n}x{n}x{n} grid: ndof={ndof}, ||N||~{muN_sorted.max():.3e} ===")
    print(f"  |mu| spectrum (smallest 8): {np.array2string(muN_sorted[:8], precision=3)}")
    print(f"  loops (|mu|<1e-3 max) = {n_loop};  smallest CHARGE |mu| = {smallest_charge_mu:.3e}")

    row = {"tag": tag, "ndof": ndof, "n_loop": n_loop,
           "smallest_charge_mu": float(smallest_charge_mu),
           "muN_max": float(muN_sorted.max())}
    for mur in (1e2, 1e4):
        chi = mur - 1.0
        A = (1.0 / chi) * np.eye(ndof) - N
        sv = np.linalg.svd(A, compute_uv=False)
        cond = sv[0] / sv[-1]
        # deflate the n_loop near-null right-singular modes, recompute cond
        # (cheap proxy for "remove the loops"): zero the smallest singular values
        cond_defl = sv[0] / sv[max(0, len(sv) - 1 - n_loop)]
        b = np.ones(ndof)
        it = jacobi_pcg_iters(A, b)
        gap = smallest_charge_mu * chi      # smallest charge |mu| / (1/chi); >>1 => Jacobi ok
        print(f"  mu_r={mur:.0e}: cond(A)={cond:.2e}, cond(deflate {n_loop})={cond_defl:.2e}, "
              f"Jacobi-GMRES iters={it}, (smallest charge |mu|)/(1/chi)={gap:.1f}")
        row[f"mur_{mur:.0e}"] = {"cond": float(cond), "cond_deflate": float(cond_defl),
                                 "jacobi_iters": int(it), "charge_mu_over_invchi": float(gap)}
    res[tag] = row

print("\n=== READING ===")
print("  If (smallest charge |mu|)/(1/chi) >> 1 AND deflating n_loop modes bounds cond AND")
print("  Jacobi iters stay bounded as mu_r grows -> the breakdown IS the loops -> Jacobi ok.")
print("  If cond stays ~mu_r even after deflation / iters grow -> weak-demag continuum -> H-ILU.")
with open(os.path.join(HERE, "demag_spectrum_jacobi.json"), "w") as f:
    json.dump(res, f, indent=2)
print("saved", os.path.join(HERE, "demag_spectrum_jacobi.json"))
