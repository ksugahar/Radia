"""radia.planar_materials -- shared 2D soft-iron constitutive laws.

ONE source of truth for the material side so a new law is written ONCE:

  * B-H table parse/validate           -> ``hm_arrays`` (H, M, chi0)
  * anhysteretic isotropic scalar law  -> ``law_from_table`` (M_of_h, chi_sec, chi0)
  * per-region (multi-grade) laws      -> ``per_region_chi`` / ``per_region_law`` (element-mats list)
  * anisotropic linear susceptibility  -> ``chi_tensor`` (easy-axis chi_par/chi_perp -> 2x2 X)

The HDiv-VIM planar layer imposes M = X.H (X = chi I for isotropic) via a per-element
(secant) susceptibility.  Dense planar helpers such as anisotropic and hysteresis solvers reuse the
same stateless laws.  Mesh-specific glue (extracting element material names) stays outside this module;
these functions take a plain ``mats`` list of per-element material names.
"""
from __future__ import annotations

import numpy as np

MU0 = 4e-7 * np.pi


# ---- B-H table -> anhysteretic isotropic law ----------------------------------------------------

def hm_arrays(bh_table):
    """[[H,B],...] (A/m, T) -> (H, M, chi0): 0-anchored, strictly-increasing-H validated, M = B/mu0 - H,
    chi0 = M[1]/H[1] the initial susceptibility.  Raises on a non-soft-iron table (B < mu0 H)."""
    tab = np.asarray(bh_table, float)
    if tab.ndim != 2 or tab.shape[1] != 2 or tab.shape[0] < 3:
        raise ValueError("bh_table must be [[H, B], ...] (A/m, T) with >= 3 rows")
    H, Bt = tab[:, 0].copy(), tab[:, 1].copy()
    if H[0] != 0.0:
        H = np.concatenate([[0.0], H]); Bt = np.concatenate([[0.0], Bt])
    if np.any(np.diff(H) <= 0):
        raise ValueError("bh_table H column must be strictly increasing")
    M = Bt / MU0 - H
    if np.any(M < -1e-9):
        raise ValueError("bh_table implies negative magnetization (B < mu0 H) -- not a soft iron")
    return H, M, M[1] / H[1]


def law_from_table(bh_table):
    """[[H,B],...] -> (M_of_h, chi_sec, chi0).  M_of_h(h)=interp with saturation clamp beyond Hmax;
    chi_sec(h)=M(h)/h secant susceptibility (chi0 below the first knee)."""
    H, M, chi0 = hm_arrays(bh_table)

    def M_of_h(h):
        return np.interp(h, H, M)                      # clamps at M[-1] beyond Hmax (saturation)

    def chi_sec(h):
        h = np.asarray(h, float)
        out = np.full_like(h, chi0)
        big = h >= H[1]
        out[big] = np.interp(h[big], H, M) / h[big]
        return out
    return M_of_h, chi_sec, chi0


# ---- per-region (multi-grade) laws over an element-materials list --------------------------------

def region_ids(mats):
    """{material_name: int-array of element indices} preserving first-seen order."""
    d = {}
    for i, m in enumerate(mats):
        d.setdefault(m, []).append(i)
    return {k: np.asarray(v, int) for k, v in d.items()}


def check_regions(mats, provided, what):
    """Raise if any mesh region has no ``what`` entry (fail-loud, lists the missing regions)."""
    missing = set(mats) - set(provided)
    if missing:
        raise ValueError("planar_materials: mesh regions %s have no %s; provided: %s"
                         % (sorted(missing), what, sorted(provided)))


def per_region_chi(mats, mu_r_dict):
    """Per-element chi from a {region: mu_r} dict (multi-grade linear soft iron)."""
    check_regions(mats, mu_r_dict, "mu_r")
    chi = np.empty(len(mats))
    for name, ids in region_ids(mats).items():
        mr = mu_r_dict[name]
        if not mr > 1.0:
            raise ValueError("planar_materials: mu_r[%r] must be > 1 (got %r)" % (name, mr))
        chi[ids] = mr - 1.0
    return chi


def per_region_law(mats, bh_dict):
    """Per-element (M_of_h, chi_sec, chi0_e) from a {region: [[H,B],...]} dict (multi-grade nonlinear).
    M_of_h / chi_sec take a per-element h array and dispatch region-by-region."""
    check_regions(mats, bh_dict, "bh_table")
    rid = region_ids(mats)
    law = {name: hm_arrays(bh_dict[name]) for name in rid}
    chi0_e = np.empty(len(mats))
    for name, ids in rid.items():
        chi0_e[ids] = law[name][2]

    def M_of_h(h):
        out = np.empty(len(mats))
        for name, ids in rid.items():
            H, M, _ = law[name]
            out[ids] = np.interp(h[ids], H, M)
        return out

    def chi_sec(h):
        out = np.empty(len(mats))
        for name, ids in rid.items():
            H, M, c0 = law[name]
            hi = h[ids]; ci = np.full_like(hi, c0)
            big = hi >= H[1]
            ci[big] = np.interp(hi[big], H, M) / hi[big]
            out[ids] = ci
        return out
    return M_of_h, chi_sec, chi0_e


# ---- anisotropic linear susceptibility (GO / oriented laminations) -------------------------------

class PlayHysteresis:
    """Prandtl-Ishlinskii scalar play-hysteresis operator (shared by both planar demag solvers).

    K play operators with thresholds ``eta`` (increasing) and weights ``w``: for input h (the signed
    field along the hysteresis axis) and committed play state p, the trial state is
    p_k' = clip(p_k, h - eta_k, h + eta_k) and M = sum_k w_k p_k'.  The INCREMENTAL susceptibility
    dM/dh = sum_k w_k * 1[|h - p_k| > eta_k] is always >= 0 (even on the descending branch), which is
    exactly what lets a NEWTON demag solve stay well-conditioned where a secant-chi Picard breaks.

    STATE IS EXTERNAL (functional): the solver holds p (n_site, K) and threads it, so the operator is
    reusable/immutable.  eta=0 (all thresholds) reduces to the linear anhysteretic chi = sum w_k.

        play = PlayHysteresis(eta=[0.1,0.3,0.6], w=[3,2,1])
        p = play.fresh_state(n);  M = play.M(H, p);  chi = play.chi_inc(H, p);  p = play.advance(H, p)
    """
    def __init__(self, eta, w):
        self.eta = np.asarray(eta, float)
        self.w = np.asarray(w, float)
        if self.eta.shape != self.w.shape or self.eta.ndim != 1 or len(self.eta) < 1:
            raise ValueError("PlayHysteresis: eta, w must be equal-length 1D arrays")
        if np.any(self.eta < 0) or np.any(np.diff(self.eta) < 0):
            raise ValueError("PlayHysteresis: eta must be >= 0 and non-decreasing")
        if np.any(self.w < 0):
            raise ValueError("PlayHysteresis: weights w must be >= 0 (monotone loop)")
        self.chi0 = float(self.w.sum())               # anhysteretic initial susceptibility

    def fresh_state(self, n):
        return np.zeros((int(n), len(self.eta)))

    def _trial(self, H, p):
        H = np.asarray(H, float)
        return np.minimum(np.maximum(p, H[:, None] - self.eta[None, :]), H[:, None] + self.eta[None, :])

    def M(self, H, p):
        """Magnetisation at signed field H (n,) given committed state p (n,K) -- does NOT advance."""
        return self._trial(H, p) @ self.w

    def chi_inc(self, H, p):
        """Incremental susceptibility dM/dH (n,), >= 0 everywhere."""
        H = np.asarray(H, float)
        return (np.abs(H[:, None] - p) > self.eta[None, :]) @ self.w

    def advance(self, H, p):
        """Return the committed play state after applying field H (n,)."""
        return self._trial(H, p)


def chi_tensor(chi_par, chi_perp, easy_deg=0.0):
    """2x2 susceptibility tensor X for a uniaxially-anisotropic linear material: susceptibility
    ``chi_par`` along the easy axis (angle ``easy_deg`` from +x) and ``chi_perp`` across it, so
    M = X.H.  X = R diag(chi_par, chi_perp) R^T (symmetric positive-definite for chi_par,chi_perp>0).

    The shared spec consumed by dense planar helpers and the HDiv-VIM roadmap; grain-oriented silicon
    steel is the canonical use (chi_par >> chi_perp)."""
    if not (chi_par > 0.0 and chi_perp > 0.0):
        raise ValueError("chi_tensor: chi_par, chi_perp must be > 0 (got %r, %r)" % (chi_par, chi_perp))
    t = np.deg2rad(easy_deg)
    c, s = np.cos(t), np.sin(t)
    R = np.array([[c, -s], [s, c]])
    return R @ np.diag([chi_par, chi_perp]) @ R.T
