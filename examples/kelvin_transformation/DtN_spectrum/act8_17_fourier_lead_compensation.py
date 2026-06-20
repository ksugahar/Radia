# SINGLE-PASS coil, FREE SPACE (no magnetic material): richer wire-path parametrization (Fourier) AND
# lead compensation.  A single circular loop has only 3 DOFs -- too few to hit a higher-order target or
# to buy back the distortion of mandatory feed leads.  Here the working span is a Fourier-shaped arc
#
#     rho(s) = rho0 + sum_k A_k sin(k*pi*s),   z(s) = z0 + sum_k C_k sin(k*pi*s),   s in [0,1]
#
# (the sin(k*pi*s) deviations VANISH at the arc ends s=0,1, so the break/junction geometry -- and hence
# the FIXED feed leads attached there -- is unchanged as the interior is reshaped).  The feed geometry
# (terminals, junctions, rho0, z0, gap) is FIXED = a manufacturing constraint; the optimizer only bends
# the working span via {I, A_k, C_k}.  Free space -> the forward map is pure Biot-Savart (no FEM).
#
# Study A (expressiveness / option 3): a higher-order single-pass design needs higher-order path DOFs --
#   fit a target that requires modes; residual drops as the number of Fourier modes K reaches the target.
# Study B (lead compensation / option 1): leads carry the full current and pass next to the workpiece,
#   distorting the patch field; with only K=0 (a bare arc) the working path cannot undo it, but adding
#   Fourier modes lets the working span CANCEL the lead field at the workpiece -> residual vs K.
import os
import json

os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from scipy.optimize import least_squares

INV4PI = 1.0 / (4.0 * np.pi)
TWO_PI = 2.0 * np.pi
RHO0, Z0, R_TERM = 0.40, 0.12, 0.50
GAP = np.deg2rad(20.0)                       # azimuthal break for the feed
# fixed feed nodes (junctions on the conductor, terminals = busbar)
J_B = np.array([RHO0 * np.cos(GAP / 2), RHO0 * np.sin(GAP / 2), Z0])
J_A = np.array([RHO0 * np.cos(GAP / 2), -RHO0 * np.sin(GAP / 2), Z0])
T_B = np.array([R_TERM * np.cos(GAP / 2), R_TERM * np.sin(GAP / 2), Z0])
T_A = np.array([R_TERM * np.cos(GAP / 2), -R_TERM * np.sin(GAP / 2), Z0])

# workpiece patch right next to the coil (standoff 0.06 below the coil plane) -- the IH target region
_phw = np.linspace(-np.deg2rad(40), np.deg2rad(40), 9)
_rw = np.linspace(0.30, 0.50, 5)
PTS = np.array([[r * np.cos(a), r * np.sin(a), 0.06] for a in _phw for r in _rw])


def _arc(I, A, C, dzJ=0.0, drJ=0.0, M_arc=240):
    """working arc; sin(k*pi*s) deviations vanish at the ends, while drJ/dzJ shift the arc-end
    (junction) baseline so the conductor has shape authority NEAR THE FEED (the leads follow)."""
    s = np.linspace(0.0, 1.0, M_arc + 1)
    phi = GAP / 2 + s * (TWO_PI - GAP)
    dev_r = np.zeros_like(s); dev_z = np.zeros_like(s)
    for k in range(len(A)):
        dev_r = dev_r + A[k] * np.sin((k + 1) * np.pi * s)
    for k in range(len(C)):
        dev_z = dev_z + C[k] * np.sin((k + 1) * np.pi * s)
    rho = (RHO0 + drJ) + dev_r; z = (Z0 + dzJ) + dev_z
    return np.stack([rho * np.cos(phi), rho * np.sin(phi), z], axis=1)


def reality_filaments(I, A, C, dzJ=0.0, drJ=0.0, M_arc=240, ML=16, MC=6):
    """the REAL single-pass coil: Fourier working arc + feed leads (terminal->junction, junction may
    move via drJ/dzJ) + terminal chord, all carry I. Terminals (busbar) stay fixed."""
    arc = _arc(I, A, C, dzJ, drJ, M_arc)
    lead_in = np.linspace(T_B, arc[0], ML + 1)
    lead_out = np.linspace(arc[-1], T_A, ML + 1)
    chord = np.linspace(T_A, T_B, MC + 1)
    pts = np.vstack([lead_in, arc[1:], lead_out[1:], chord[1:]])
    return pts[:-1], pts[1:], I * np.ones(len(pts) - 1)


def ideal_loop_filaments(I, M=240):
    """the lead-free IDEAL: a clean closed circle (rho0, z0), no leads."""
    phi = np.linspace(0.0, TWO_PI, M + 1)
    pts = np.stack([RHO0 * np.cos(phi), RHO0 * np.sin(phi), Z0 * np.ones_like(phi)], axis=1)
    return pts[:-1], pts[1:], I * np.ones(M)


def biot(sp, ep, Iseg, P=PTS):
    Ce = 0.5 * (sp + ep); Kwe = Iseg[:, None] * (ep - sp)
    out = np.zeros((len(P), 3))
    for i, p in enumerate(P):
        d = p - Ce; r3 = (np.einsum('ij,ij->i', d, d)) ** 1.5 + 1e-300
        out[i] = INV4PI * (np.cross(Kwe, d) / r3[:, None]).sum(axis=0)
    return out.reshape(-1)


def unpack(theta, K, free_feed=False):
    """theta = [I, (dzJ, drJ if free_feed), A_1..A_K, C_1..C_K]."""
    I = theta[0]; off = 1
    dzJ, drJ = (theta[1], theta[2]) if free_feed else (0.0, 0.0)
    if free_feed:
        off = 3
    A = list(theta[off:off + K]); C = list(theta[off + K:off + 2 * K])
    return I, A, C, dzJ, drJ


def fit(target, K, free_feed=False):
    theta0 = np.concatenate([[0.8], (np.zeros(2) if free_feed else []), np.zeros(2 * K)])

    def resid(theta):
        I, A, C, dzJ, drJ = unpack(theta, K, free_feed)
        sp, ep, Iseg = reality_filaments(I, A, C, dzJ, drJ)
        return biot(sp, ep, Iseg) - target

    sol = least_squares(resid, theta0, xtol=1e-12, ftol=1e-14, diff_step=1e-4)
    return np.linalg.norm(resid(sol.x)) / np.linalg.norm(target), sol.x


print("SINGLE-PASS coil, FREE SPACE: Fourier path DOFs (option 3) + lead compensation (option 1)")
print("  feed FIXED: rho0=%.2f z0=%.2f gap=%g deg, terminal R=%.2f; workpiece patch=%d pts (standoff 0.06)"
      % (RHO0, Z0, np.rad2deg(GAP), R_TERM, len(PTS)))

# ===== Study A: expressiveness -- a higher-order single-pass design needs higher-order path DOFs =====
A_ref = [0.0, 0.06]; C_ref = [0.04, 0.0]            # reference design: radial mode 2 + axial mode 1
tgt_A = biot(*reality_filaments(1.0, A_ref, C_ref))
print("\n(A) EXPRESSIVENESS: fit a reference single-pass design (needs radial m=2 + axial m=1):")
print("    K modes | nDOF | rel. residual")
resA = {}
for K in (0, 1, 2, 3):
    r, _ = fit(tgt_A, K)
    resA[K] = r
    print("       %d    |  %2d  |   %.3e%s" % (K, 1 + 2 * K, r, "   <- enough DOFs" if K >= 2 else ""))

# ===== Study B: lead compensation -- buy back the mandatory-lead distortion with path shape =====
tgt_B = biot(*ideal_loop_filaments(1.0))            # lead-free ideal field at the workpiece
raw = np.linalg.norm(biot(*reality_filaments(1.0, [], [])) - tgt_B) / np.linalg.norm(tgt_B)
print("\n(B) LEAD COMPENSATION: target = lead-free ideal; reality has FIXED busbar terminals + leads.")
print("    raw distortion (ideal vs bare arc+leads, current-matched only) ~ %.1f%%" % (100 * raw))
print("    K modes | pinned-feed resid | free-feed resid (junction dz,dr free -> shape authority at feed)")
resB, resBf = {}, {}
for K in (0, 1, 2, 3, 4, 6):
    rp, _ = fit(tgt_B, K, free_feed=False)
    rf, _ = fit(tgt_B, K, free_feed=True)
    resB[K], resBf[K] = rp, rf
    print("       %d    |     %.3e     |   %.3e" % (K, rp, rf))

print("\n[RESULT] (A) Fourier path DOFs let a single-pass coil realize a higher-order design"
      " (residual %.1e at K>=2 vs %.1e at K=0).  (B) Active compensation of mandatory feed leads at a"
      " CLOSE workpiece is limited: %.1f%% raw -> %.1f%% (pinned K=6) -> %.1f%% (free-feed K=6); the"
      " lead field is sharply localized under the leads, so geometric cancellation (bifilar/coaxial"
      " leads, smaller gap) is the effective lever, not working-path reshaping."
      % (resA[2], resA[0], 100 * raw, 100 * resB[6], 100 * resBf[6]))

out = {
    "feed": {"rho0": RHO0, "z0": Z0, "gap_deg": np.rad2deg(GAP), "R_term": R_TERM, "patch_pts": len(PTS)},
    "expressiveness": {str(k): v for k, v in resA.items()},
    "compensation": {"raw": raw, "pinned": {str(k): v for k, v in resB.items()},
                     "free_feed": {str(k): v for k, v in resBf.items()}},
}
with open(os.path.join(os.path.dirname(__file__), "demo_sp3_results.json"), "w") as fh:
    json.dump(out, fh, indent=2)
print("  results -> demo_sp3_results.json")
