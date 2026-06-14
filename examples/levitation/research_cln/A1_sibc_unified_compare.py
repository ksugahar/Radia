"""A1 cuboid (17.72×17.72×2 mm Cu) Cauer→SIBC unified Y(jω) comparison.

Compares the unified Cauer + SIBC admittance constructed from
  (A) VIM mpmath HDivDivFreeHex single-element pipeline
      (extract_tau_A1_square_mp.py  --smoke-p3 --save-foster)
  (B) FEM NGSolve Tanimoto-AT alternating Poisson pipeline
      (tanimoto_AT_geom_sweep.py    --geom square)

For each method we:
  1. Reconstruct the Cauer-I rungs (R_2k, L_{2k+1}) up to a chosen stage N
  2. Build the standard J-fraction p_array = [R_0, 1/L_1, R_2, 1/L_3, ...]
  3. Compute the corner frequency ω_N via Matsuo eq (21) at α=1/2
  4. Fix K_SIBC by continuity Y_Cauer(jω_N) = K_SIBC · √(σ/(jω_N μ))
  5. Evaluate Y_unified(jω) over a frequency sweep
  6. Normalize each method's Y by its own Y(0) = 1/R_0 (admittance units differ
     between VIM admittance-per-unit-area and FEM Tanimoto network admittance
     by an overall energy-normalization constant — only the SHAPE of the
     frequency response is convention-independent)

Generates:
  - A1_sibc_unified_compare.json   numerical data
  - A1_sibc_unified_compare.png    log–log |Y(jω)|/|Y(0)| vs f plot
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import List, Tuple

sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)

from mpmath import mp, mpf

# A1 / Cu constants (matched to VIM extract_tau_A1_square_mp.py)
A_M = mpf("17.72e-3")
B_M = mpf("17.72e-3")
C_M = mpf("2e-3")        # slab half-thickness analogue: t = 2 mm, half = 1 mm
C_HALF_M = C_M / mpf(2)  # = 1 mm — used for the asymptotic SIBC reference
SIGMA = mpf("5.8e7")
MU0 = 4 * mp.pi * mpf("1e-7")


# ---------------------------------------------------------------------------
# Cauer-I J-fraction evaluation
# ---------------------------------------------------------------------------

def build_p_array_from_RL(R_list: List[mpf], L_list: List[mpf]) -> List[mpf]:
    """Pack (R_0, R_2, R_4,...) and (L_1, L_3, L_5,...) into the J-fraction
    coefficient sequence used by evaluate_jfraction_Z:
      p = [R_0, 1/L_1, R_2, 1/L_3, R_4, 1/L_5, ...]
    Length is min(2N, 2N+1) depending on whether L_list has the same length
    as R_list or one less.
    """
    n = min(len(R_list), len(L_list))
    p: List[mpf] = []
    for k in range(n):
        p.append(mpf(R_list[k]))
        p.append(mpf(1) / mpf(L_list[k]))
    # If we have one extra R (R_{2N}) without a paired L, append it (terminating
    # the ladder with a final shunt resistor).
    if len(R_list) > n:
        p.append(mpf(R_list[n]))
    return p


def evaluate_jfraction_Z(p_array, s):
    """Z(s) = p[0] + s/(p[1] + s/(p[2] + s/(p[3] + ...))) — Cauer-I J-fraction.

    Matches the convention in extract_tau_A1_square_mp.py.  Z(0) = p[0] = R_0,
    so Y(0) = 1/p[0].
    """
    n = len(p_array)
    if n == 0:
        return mp.mpc(0, 0)
    Z_inner = mpf(p_array[n - 1])
    for i in range(n - 2, -1, -1):
        Z_inner = mpf(p_array[i]) + s / Z_inner
    return Z_inner


def evaluate_Y_truncated(p_array, omega):
    """Y(jω) = 1/Z(jω) of the truncated Cauer ladder."""
    s = mp.mpc(0, omega)
    Z = evaluate_jfraction_Z(p_array, s)
    return mp.mpc(1, 0) / Z


# ---------------------------------------------------------------------------
# Corner frequency (Matsuo SA-26-014 eq 21) and SIBC scaling
# ---------------------------------------------------------------------------

def compute_corner_frequency(R_last: mpf, L_last: mpf, alpha: float = 0.5) -> mpf:
    """ω_N = √(2^(α/2) - 1) · R_last / L_last.  α=1/2 → factor ≈ 0.4350."""
    factor = mp.sqrt(mp.power(mpf(2), mpf(repr(alpha)) / mpf(2)) - mpf(1))
    return factor * mpf(abs(R_last)) / mpf(abs(L_last))


def compute_K_sibc(p_array, omega_N, sigma, mu0_val):
    """Determine SIBC scaling K_sibc by continuity Y_Cauer(jω_N) = K · √(σ/(jω_N μ))."""
    Y_cauer = evaluate_Y_truncated(p_array, omega_N)
    s_N = mp.mpc(0, omega_N)
    Y_sibc_base = mp.sqrt(sigma / (s_N * mu0_val))
    K = Y_cauer / Y_sibc_base
    return K, Y_cauer, Y_sibc_base


def Y_unified_eval(p_array, R_last, L_last, omega_hz, sigma, mu0_val, *,
                   alpha: float = 0.5):
    """Y_unified(j·2πf) — Cauer (ω<ω_N) ⊕ SIBC (ω≥ω_N) basis switch."""
    omega = mpf(repr(2.0)) * mp.pi * mpf(repr(omega_hz))
    omega_N = compute_corner_frequency(R_last, L_last, alpha=alpha)
    K_sibc, _, _ = compute_K_sibc(p_array, omega_N, sigma, mu0_val)
    if abs(omega) < abs(omega_N):
        return evaluate_Y_truncated(p_array, omega), "Cauer", omega_N, K_sibc
    s = mp.mpc(0, omega)
    return (K_sibc * mp.sqrt(sigma / (s * mu0_val)),
            "SIBC", omega_N, K_sibc)


# ---------------------------------------------------------------------------
# Input loading: VIM Foster spectrum & FEM Tanimoto rungs
# ---------------------------------------------------------------------------

def load_vim_rungs_from_foster(
    foster_json: Path, *, n_cauer_stages: int = 6, tau_min_us: float = 0.0,
) -> Tuple[List[dict], List[mpf]]:
    """Load VIM Foster spectrum and re-extract Cauer rungs via QD-Wynn.

    Imports the relevant functions from extract_tau_A1_square_mp.  Returns
    rungs (list of dicts) and the raw J-fraction p_array.
    """
    sys.path.insert(0, str(Path(
        "S:/Radia/01_GitHub/src/ext/radia_vim/scripts"
    ).resolve()))
    from extract_tau_A1_square_mp import (
        cauer_from_spectrum_mp, load_foster_spectrum_json,
    )
    spec = load_foster_spectrum_json(foster_json)
    print(f"VIM Foster spectrum: order={spec['order']}, "
          f"n_dofs={spec['n_dofs']}, dps={spec['dps']}, "
          f"B={spec['B_dir']}")
    print(f"  leading τ = {mp.nstr(spec['leading_tau_us'], 10)} μs")
    rungs, c0_rel, p_array = cauer_from_spectrum_mp(
        spec["tau_us"], spec["g2"],
        n_use=min(80, spec["n_dofs"]),
        n_moments=max(20, 4 * n_cauer_stages),
        max_stages=2 * n_cauer_stages,
        tau_min_us=tau_min_us,
    )
    return rungs, list(p_array)


def load_fem_rungs(fem_json: Path, *, n_stages: int = None,
                    skip_stages: int = 0) -> Tuple[List[dict], List[mpf]]:
    """Load FEM Tanimoto-AT (R_n, L_n) rungs and assemble J-fraction p_array.

    Tanimoto JSON layout (tanimoto_AT_geom_sweep.py):
      Rn = [R_0, R_2, R_4, ...]  (one entry per "stage 0..N" T-equation solve)
      Ln = [L_1, L_3, L_5, ...]  (one entry per A-equation solve)

    skip_stages: drop the first N rungs from Tanimoto output before re-indexing.
    Justification: Tanimoto Stage 0 is the t=0+ skin-layer initial response,
    not the DC-limit Kameari R_0 — empirically τ_pair[0] is ~5–6× the bulk
    leading time constant.  Setting skip_stages=1 treats stages 1, 2, … as
    the Cauer ladder for cross-comparison with VIM Kameari rungs.
    """
    data = json.loads(Path(fem_json).read_text(encoding="utf-8"))
    Rn = [mpf(repr(r)) for r in data["Rn"]]
    Ln = [mpf(repr(l)) for l in data["Ln"]]
    if skip_stages > 0:
        Rn = Rn[skip_stages:]
        Ln = Ln[skip_stages:]
    if n_stages is not None:
        Rn = Rn[: n_stages + 1]
        Ln = Ln[: n_stages]
    rungs = []
    for k in range(min(len(Rn) - 1, len(Ln))):
        R_2k = Rn[k]
        L_k = Ln[k]
        rungs.append({
            "k": k,
            "R_2k": R_2k,
            "L_2k_plus_1": L_k,
            "tau_pair_us": L_k / R_2k * mpf("1e6"),
        })
    p_array = build_p_array_from_RL(Rn[:len(Ln) + 1], Ln)
    return rungs, p_array


# ---------------------------------------------------------------------------
# Frequency sweep + analysis
# ---------------------------------------------------------------------------

def frequency_sweep(p_array, rungs, *, alpha: float = 0.5,
                    f_min_hz: float = 1e0, f_max_hz: float = 1e8,
                    n_omega: int = 80):
    """Evaluate Y_unified(j·2πf) over a log-spaced frequency grid.

    Returns
    -------
    freqs : list of float
    Yabs  : list of float
    Yphase_deg : list of float
    regime : list of str
    omega_N, K_sibc : mpf, mpc
    """
    if not rungs:
        return [], [], [], [], None, None
    R_last = rungs[-1]["R_2k"]
    L_last = rungs[-1]["L_2k_plus_1"]
    omega_N_val = compute_corner_frequency(R_last, L_last, alpha=alpha)
    K_sibc_val, _, _ = compute_K_sibc(p_array, omega_N_val, SIGMA, MU0)
    log_min = math.log10(f_min_hz)
    log_max = math.log10(f_max_hz)
    freqs = [10 ** (log_min + (log_max - log_min) * k / (n_omega - 1))
             for k in range(n_omega)]
    Yabs, Yph, regime = [], [], []
    for f in freqs:
        omega = mpf(repr(2 * math.pi * f))
        if abs(omega) < abs(omega_N_val):
            Y = evaluate_Y_truncated(p_array, omega)
            tag = "Cauer"
        else:
            s = mp.mpc(0, omega)
            Y = K_sibc_val * mp.sqrt(SIGMA / (s * MU0))
            tag = "SIBC"
        Yabs.append(float(abs(Y)))
        Yph.append(float(mp.atan2(mp.im(Y), mp.re(Y))) * 180.0 / math.pi)
        regime.append(tag)
    return freqs, Yabs, Yph, regime, omega_N_val, K_sibc_val


def analytical_sibc_reference(freqs: List[float], *, scale: float = 1.0):
    """|K · √(σ/(jωμ))| at each frequency, with arbitrary scale (set to overlap)."""
    out = []
    for f in freqs:
        omega = mpf(repr(2 * math.pi * f))
        s = mp.mpc(0, omega)
        Y = scale * mp.sqrt(SIGMA / (s * MU0))
        out.append(float(abs(Y)))
    return out


# ---------------------------------------------------------------------------
# Reporting + plotting
# ---------------------------------------------------------------------------

def print_rungs(label, rungs, p_array):
    print()
    print(f"=== {label} ===")
    print(f"  n rungs = {len(rungs)}, p_array length = {len(p_array)}")
    print(f"  R_0 (= 1/Y(0)) = {mp.nstr(rungs[0]['R_2k'], 8)}")
    print(f"  {'k':>3}  {'R_2k':>22}  {'L_(2k+1)':>22}  {'tau_pair (us)':>14}")
    for r in rungs:
        print(f"  {r['k']:>3}  {mp.nstr(r['R_2k'], 10):>22}  "
              f"{mp.nstr(r['L_2k_plus_1'], 10):>22}  "
              f"{mp.nstr(r['tau_pair_us'], 8):>14}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--vim-foster", type=str, required=True,
                    help="VIM Foster spectrum JSON (from --save-foster)")
    ap.add_argument("--fem-results", type=str, required=True,
                    help="FEM Tanimoto-AT results JSON")
    ap.add_argument("--ac-sweep", type=str, default=None,
                    help="Optional NGSolve AC sweep JSON for ground-truth "
                         "overlay (from A1_AC_sweep.py)")
    ap.add_argument("--n-stages", type=int, default=3,
                    help="Number of Cauer stages to use for BOTH methods "
                         "(default 3). Overridden by --vim-n-stages / "
                         "--fem-n-stages.")
    ap.add_argument("--vim-n-stages", type=int, default=None,
                    help="Number of Cauer stages to use for VIM only.")
    ap.add_argument("--fem-n-stages", type=int, default=None,
                    help="Number of Cauer stages to use for FEM only.")
    ap.add_argument("--fem-skip-stages", type=int, default=0,
                    help="Drop first N Tanimoto-AT rungs before forming the "
                         "Cauer ladder (default 0). Use 1 to treat Tanimoto "
                         "Stage 0 as SIBC seed and Stages 1+ as bulk Cauer.")
    ap.add_argument("--tau-min-us", type=float, default=0.0,
                    help="VIM τ-filter (drop spurious K-null modes)")
    ap.add_argument("--dps", type=int, default=50,
                    help="mpmath precision for evaluation (default 50)")
    ap.add_argument("--f-min", type=float, default=1.0)
    ap.add_argument("--f-max", type=float, default=1e8)
    ap.add_argument("--n-omega", type=int, default=120)
    ap.add_argument("--alpha", type=float, default=0.5,
                    help="CPE α for corner frequency (default 0.5 = Warburg/SIBC)")
    ap.add_argument("--out-tag", type=str,
                    default="A1_sibc_unified_compare")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--compact", action="store_true",
                    help="One-panel magnitude-only plot for 1-page digest.")
    args = ap.parse_args()

    mp.dps = args.dps

    print("=" * 78)
    print(" A1 cuboid (17.72×17.72×2 mm Cu) — VIM vs FEM Cauer→SIBC unification")
    print("=" * 78)
    print(f" mpmath dps   = {mp.dps}")
    print(f" stages used  = {args.n_stages}")
    print(f" α (CPE)      = {args.alpha}")
    print(f" σ = {float(SIGMA):.3g} S/m,  μ_0 = {float(MU0):.6g} H/m")
    print(f" c_half = {float(C_HALF_M) * 1000:g} mm "
          f"(reference for SIBC asymptote scale)")
    print(f" f sweep      = [{args.f_min:g}, {args.f_max:g}] Hz, "
          f"n_omega = {args.n_omega}")
    print()

    vim_n = args.vim_n_stages if args.vim_n_stages is not None else args.n_stages
    fem_n = args.fem_n_stages if args.fem_n_stages is not None else args.n_stages

    # --- VIM ---
    print(f"Loading VIM Foster: {args.vim_foster}")
    vim_rungs, vim_p_array = load_vim_rungs_from_foster(
        Path(args.vim_foster),
        n_cauer_stages=max(vim_n, 6),
        tau_min_us=args.tau_min_us,
    )
    vim_rungs = vim_rungs[: vim_n]
    vim_p_array = vim_p_array[: 2 * vim_n]
    print_rungs(f"VIM mpmath HDivDivFreeHex p=3 (n={vim_n})",
                vim_rungs, vim_p_array)

    # --- FEM ---
    print(f"\nLoading FEM Tanimoto-AT: {args.fem_results}")
    fem_rungs, fem_p_array = load_fem_rungs(
        Path(args.fem_results), n_stages=fem_n,
        skip_stages=args.fem_skip_stages,
    )
    fem_label = "FEM NGSolve Tanimoto-AT"
    if args.fem_skip_stages > 0:
        fem_label += f" (skipped first {args.fem_skip_stages} stages)"
    print_rungs(fem_label, fem_rungs, fem_p_array)

    # --- Sweep ---
    print("\nSweeping Y_unified(jω)...")
    vim_freqs, vim_Yabs, vim_Yph, vim_regime, vim_omegaN, vim_K = (
        frequency_sweep(vim_p_array, vim_rungs, alpha=args.alpha,
                        f_min_hz=args.f_min, f_max_hz=args.f_max,
                        n_omega=args.n_omega))
    fem_freqs, fem_Yabs, fem_Yph, fem_regime, fem_omegaN, fem_K = (
        frequency_sweep(fem_p_array, fem_rungs, alpha=args.alpha,
                        f_min_hz=args.f_min, f_max_hz=args.f_max,
                        n_omega=args.n_omega))

    vim_Y0 = float(1.0 / mpf(vim_rungs[0]["R_2k"]))
    fem_Y0 = float(1.0 / mpf(fem_rungs[0]["R_2k"]))

    print(f"\n VIM: ω_N = {mp.nstr(vim_omegaN, 8)} rad/s "
          f"(f_N = {float(vim_omegaN)/(2*math.pi):.3e} Hz),  "
          f"|K_sibc| = {mp.nstr(abs(vim_K), 8)}")
    print(f" FEM: ω_N = {mp.nstr(fem_omegaN, 8)} rad/s "
          f"(f_N = {float(fem_omegaN)/(2*math.pi):.3e} Hz),  "
          f"|K_sibc| = {mp.nstr(abs(fem_K), 8)}")
    print(f"\n VIM Y(0) = 1/R_0 = {vim_Y0:.6e}")
    print(f" FEM Y(0) = 1/R_0 = {fem_Y0:.6e}")
    print(f" ratio Y_FEM(0)/Y_VIM(0) = {fem_Y0 / vim_Y0:.6e}")

    # --- Save JSON ---
    out_tag = args.out_tag
    out_dir = Path(__file__).resolve().parent
    out_json = out_dir / f"{out_tag}.json"
    out = {
        "method_VIM": "mpmath HDivDivFreeHex single-element VIM",
        "method_FEM": "NGSolve HCurl Tanimoto-AT alternating Poisson",
        "geometry_mm": [float(A_M) * 1000, float(B_M) * 1000,
                        float(C_M) * 1000],
        "sigma_S_per_m": float(SIGMA),
        "n_stages_used": args.n_stages,
        "alpha_CPE": args.alpha,
        "vim_rungs": [
            {"k": r["k"],
             "R_2k": mp.nstr(mpf(r["R_2k"]), 20),
             "L_2k_plus_1": mp.nstr(mpf(r["L_2k_plus_1"]), 20),
             "tau_pair_us": mp.nstr(mpf(r["tau_pair_us"]), 12)}
            for r in vim_rungs
        ],
        "fem_rungs": [
            {"k": r["k"],
             "R_2k": mp.nstr(mpf(r["R_2k"]), 20),
             "L_2k_plus_1": mp.nstr(mpf(r["L_2k_plus_1"]), 20),
             "tau_pair_us": mp.nstr(mpf(r["tau_pair_us"]), 12)}
            for r in fem_rungs
        ],
        "vim_omega_N_rad_per_s": mp.nstr(vim_omegaN, 16),
        "vim_f_N_Hz": float(vim_omegaN) / (2 * math.pi),
        "vim_K_sibc_abs": mp.nstr(abs(vim_K), 12),
        "vim_K_sibc_arg_deg": (
            float(mp.atan2(mp.im(vim_K), mp.re(vim_K))) * 180.0 / math.pi),
        "fem_omega_N_rad_per_s": mp.nstr(fem_omegaN, 16),
        "fem_f_N_Hz": float(fem_omegaN) / (2 * math.pi),
        "fem_K_sibc_abs": mp.nstr(abs(fem_K), 12),
        "fem_K_sibc_arg_deg": (
            float(mp.atan2(mp.im(fem_K), mp.re(fem_K))) * 180.0 / math.pi),
        "frequencies_Hz": vim_freqs,
        "Y_VIM_abs": vim_Yabs,
        "Y_VIM_phase_deg": vim_Yph,
        "Y_VIM_regime": vim_regime,
        "Y_FEM_abs": fem_Yabs,
        "Y_FEM_phase_deg": fem_Yph,
        "Y_FEM_regime": fem_regime,
        "vim_Y0": vim_Y0,
        "fem_Y0": fem_Y0,
    }
    out_json.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nSaved data: {out_json}")

    if args.no_plot:
        return

    # --- Plot ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"matplotlib unavailable, skipping plot: {e}")
        return

    if args.compact:
        fig, ax1 = plt.subplots(1, 1, figsize=(5.5, 3.2))
        ax2 = None
    else:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 8.5), sharex=True)

    vim_norm = [y / vim_Y0 for y in vim_Yabs]
    fem_norm = [y / fem_Y0 for y in fem_Yabs]

    # Optional AC sweep overlay — direct NGSolve eddy-current frequency sweep
    # acts as the independent ground truth.
    #
    # For a Foster decomposition with M-orthonormal e_n and g_n = ⟨A_ext, e_n⟩_M,
    # the AC modal amplitudes satisfy c_n(jω) = -jωτ_n g_n / (1+jωτ_n), so
    #   ⟨A_red, A_ext⟩(jω) = Σ c_n g_n = -Σ g_n²·jωτ_n/(1+jωτ_n)
    #                      = -Y_VIM(0) + Y_VIM(jω).
    # Therefore  Y_AC(jω) := |A_ext|²_L2 + ⟨A_red, A_ext⟩(jω)  has the same
    # Foster-mode roll-off as Y_VIM(jω) and starts at |A_ext|²_L2 (DC) →
    # 0 (high-ω full Faraday shielding).  After normalising by their respective
    # DC values, AC ground truth and CLN+SIBC are directly comparable in shape.
    ac_freqs = ac_abs = ac_phase = None
    if args.ac_sweep:
        ac_data = json.loads(Path(args.ac_sweep).read_text(encoding="utf-8"))
        ac_freqs = [r["f_Hz"] for r in ac_data["results"]]
        A_ext_norm2 = ac_data.get("A_ext_norm2_L2", 1.0)
        ac_complex = []
        for r in ac_data["results"]:
            ar = complex(r.get("ared_dot_aext_real", 0.0),
                         r.get("ared_dot_aext_imag", 0.0))
            Y_proxy = A_ext_norm2 + ar
            ac_complex.append(Y_proxy)
        ac_abs = [abs(y) for y in ac_complex]
        ac_phase = [math.degrees(math.atan2(y.imag, y.real))
                    for y in ac_complex]
        ac_Y0 = ac_abs[0]
        ac_norm = [a / ac_Y0 for a in ac_abs]
        ax1.loglog(ac_freqs, ac_norm, "o", color="C2", ms=6, mfc="none",
                   mew=1.6,
                   label=r"NGSolve AC sweep  $|A_{ext}|^2+\langle A_{red},A_{ext}\rangle$")

    ax1.loglog(vim_freqs, vim_norm, "-", color="C0", lw=2,
               label="VIM (mpmath HDivDivFreeHex p=3) + SIBC")
    ax1.loglog(fem_freqs, fem_norm, "--", color="C1", lw=2,
               label="FEM (NGSolve Tanimoto-AT) + SIBC")

    # Ideal SIBC asymptote |1/√(jω)| ∝ 1/√f.  Anchor at each method's ω_N so the
    # asymptote passes through |Y_unified(ω_N)|/|Y(0)| — this visualizes the
    # seam point per method without needing to know absolute units.
    f_N_v = float(vim_omegaN) / (2 * math.pi)
    f_N_f = float(fem_omegaN) / (2 * math.pi)

    # find indexes closest to each ω_N for the anchor points
    iN_v = min(range(len(vim_freqs)), key=lambda i: abs(vim_freqs[i] - f_N_v))
    iN_f = min(range(len(fem_freqs)), key=lambda i: abs(fem_freqs[i] - f_N_f))
    anchor_v = vim_norm[iN_v]
    anchor_f = fem_norm[iN_f]
    vim_sibc = [anchor_v * math.sqrt(f_N_v / f) for f in vim_freqs]
    fem_sibc = [anchor_f * math.sqrt(f_N_f / f) for f in fem_freqs]
    ax1.loglog(vim_freqs, vim_sibc, ":", color="C0", lw=1.0, alpha=0.6,
               label=r"VIM SIBC asymptote $\propto 1/\sqrt{\omega}$")
    ax1.loglog(fem_freqs, fem_sibc, ":", color="C1", lw=1.0, alpha=0.6,
               label=r"FEM SIBC asymptote $\propto 1/\sqrt{\omega}$")

    ax1.axvline(f_N_v, color="C0", ls=":", alpha=0.4,
                label=f"ω_N (VIM) = {f_N_v:.2e} Hz")
    ax1.axvline(f_N_f, color="C1", ls=":", alpha=0.4,
                label=f"ω_N (FEM) = {f_N_f:.2e} Hz")

    ax1.set_ylim(top=2.0)
    ax1.set_ylabel(r"$|Y(j\omega)| / |Y(0)|$  (normalized)")
    if args.vim_n_stages or args.fem_n_stages:
        stage_str = f"VIM {vim_n} + FEM {fem_n} stages"
    else:
        stage_str = f"{args.n_stages} stages"
    ax1.set_title(f"A1 cuboid 17.72×17.72×2 mm Cu — Cauer + SIBC unified Y "
                  f"({stage_str})")
    ax1.legend(fontsize=8, loc="lower left")
    ax1.grid(True, which="both", alpha=0.4)

    if ax2 is None:
        ax1.set_xlabel(r"frequency $f$  [Hz]")
    else:
        if ac_freqs is not None:
            ax2.semilogx(ac_freqs, ac_phase, "o", color="C2", ms=6, mfc="none",
                         mew=1.6, label="NGSolve AC sweep ground truth")
        ax2.semilogx(vim_freqs, vim_Yph, "-", color="C0", lw=2, label="VIM")
        ax2.semilogx(fem_freqs, fem_Yph, "--", color="C1", lw=2, label="FEM")
        ax2.axhline(-45.0, color="gray", ls=":", alpha=0.6,
                    label=r"ideal SIBC limit arg(Y) = -45°")
        ax2.axvline(f_N_v, color="C0", ls=":", alpha=0.4)
        ax2.axvline(f_N_f, color="C1", ls=":", alpha=0.4)
        ax2.set_xlabel(r"frequency $f$  [Hz]")
        ax2.set_ylabel(r"arg($Y$)  [deg]")
        ax2.legend(fontsize=8, loc="lower left")
        ax2.grid(True, which="both", alpha=0.4)

    fig.tight_layout()
    out_png = out_dir / f"{out_tag}.png"
    out_pdf = out_dir / f"{out_tag}.pdf"
    fig.savefig(out_png, dpi=140)
    fig.savefig(out_pdf)
    print(f"Saved plot: {out_png}")
    print(f"Saved plot: {out_pdf}")


if __name__ == "__main__":
    main()
