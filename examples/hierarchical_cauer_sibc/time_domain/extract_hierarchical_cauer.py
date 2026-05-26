"""
Hierarchical Cauer parameter extraction from a Y(omega) frequency sweep.

Pipeline:
  Input:  Y_data = Y(omega_k) sampled on log-spaced omega grid
          + geometry (S = surface area, sigma, mu)
  Stage 1: estimate K_SIBC = S sqrt(sigma/mu) (from geometry alone)
  Stage 2: 1-parameter fit of d (wall-band crossover) by minimising
           |Y_residual(omega) - K_SIBC sqrt(s)/(s+d)|
  Stage 3: bulk Foster fit -- N_CLN modes to capture
           Y_CLN(s) = Y_data(s) - K_SIBC sqrt(s)/(s+d)
           (the bulk part, after removing the Warburg tail)
  Stage 4: diffusive Foster discretization of the Warburg block
           K_SIBC sqrt(s)/(s+d) -> N_xi mode pairs (c_w, xi_w)

Output: SPICE netlist + JSON parameter file ready for
        spice_subcircuit_generator.py.

Usage:
    python extract_hierarchical_cauer.py --input <Y_freq_sweep.npz>
    python extract_hierarchical_cauer.py   (self-test on cylinder Foster)
"""

import argparse
import json
import sys
import numpy as np
from pathlib import Path
from scipy.optimize import minimize_scalar
from scipy import special

# Re-use diffusive quadrature
sys.path.insert(0, str(Path(__file__).parent))
from diffusive_quadrature import fit_diffusive


def Y_warburg(s, K_SIBC, d):
    """Warburg block."""
    return K_SIBC * np.sqrt(s) / (s + d)


def fit_d(omega, Y_data, K_SIBC, d_range=(1e2, 1e8), N_CLN_test=8):
    """
    Stage 2: fit crossover d by minimising the FULL hierarchical-Cauer
    reconstruction error over the entire frequency band.  For each
    candidate d, the bulk Foster modes (Stage 3) are re-fitted from
    the residual; the d that minimises the residual after bulk fit is
    the optimal d (no separate "wall-band only" weighting needed).
    """
    s = 1j * omega

    def err_d(log_d):
        d_try = 10**log_d
        Y_w = Y_warburg(s, K_SIBC, d_try)
        Y_residual = Y_data - Y_w
        xi_CLN_try, c_CLN_try, err_CLN = fit_bulk_foster(
            omega, Y_residual, N_CLN_test, reference_Y=Y_data)
        Y_total = Y_w + np.sum(c_CLN_try[None, :]
                                / (s[:, None] + xi_CLN_try[None, :]), axis=1)
        rel_err = np.abs(Y_total - Y_data) / np.abs(Y_data)
        return np.max(rel_err)

    log_d_lo, log_d_hi = np.log10(d_range)
    res = minimize_scalar(err_d, bounds=(log_d_lo, log_d_hi),
                          method='bounded',
                          options={'xatol': 0.02})
    return 10**res.x, res.fun


def fit_bulk_foster(omega, Y_residual, N_CLN, xi_range_decades=4,
                     reference_Y=None):
    """
    Stage 3: fit N_CLN-mode Foster expansion to Y_residual(s) using
    relative-error weighting (weights = 1/|reference_Y|, defaults to Y_residual).
    Y_CLN(s) = sum_k c_k / (s + xi_k), xi_k log-spaced.
    Returns xi_CLN, c_CLN (real residues by Re/Im stacking).
    """
    omega_min, omega_max = omega.min(), omega.max()
    log_xi_min = np.log10(omega_min) - 1
    log_xi_max = np.log10(omega_max) - 1
    xi_CLN = np.logspace(log_xi_min, log_xi_max, N_CLN)

    s = 1j * omega
    A = 1.0 / (s[:, None] + xi_CLN[None, :])

    # Relative-error weighting
    weights = np.abs(reference_Y if reference_Y is not None else Y_residual)
    weights = np.maximum(weights, weights.max() * 1e-12)
    A_w = A / weights[:, None]
    Y_w = Y_residual / weights

    A_real = np.vstack([A_w.real, A_w.imag])
    b_real = np.hstack([Y_w.real, Y_w.imag])
    c_CLN, _, _, _ = np.linalg.lstsq(A_real, b_real, rcond=None)
    # Evaluate relative fit error
    Y_fit = np.sum(c_CLN[None, :] / (s[:, None] + xi_CLN[None, :]), axis=1)
    rel_err = np.abs(Y_fit - Y_residual) / weights
    err = np.max(rel_err)
    return xi_CLN, c_CLN, err


def emit_spice_netlist(out_path, xi_CLN, c_CLN, xi_w, c_w, K_SIBC, d, comment=""):
    """Generate SPICE subcircuit netlist."""
    lines = [
        "* Hierarchical Cauer SPICE subcircuit",
        f"* {comment}",
        f"* K_SIBC = {K_SIBC:.4e}, d = {d:.3e} rad/s",
        f"* N_CLN bulk modes = {len(xi_CLN)}",
        f"* N_xi Warburg diffusive modes = {len(xi_w)}",
        "",
        ".subckt hierarchical_cauer in out",
        "* Bulk CLN Foster ladder",
    ]
    for k, (xi, c) in enumerate(zip(xi_CLN, c_CLN)):
        lines.append(
            f"G_cln{k:02d} in out LAPLACE {{V(in,out)}} "
            f"{{{c:.6e} / (s + {xi:.6e})}}"
        )
    lines.append("* Warburg diffusive surface ladder")
    for j, (xi_j, c_j) in enumerate(zip(xi_w, c_w)):
        c_scaled = K_SIBC * c_j.real
        lines.append(
            f"G_w{j:02d} in out LAPLACE {{V(in,out)}} "
            f"{{{c_scaled:.6e} / (s + {xi_j:.6e})}}"
        )
    lines.append(".ends hierarchical_cauer")
    lines.append("")
    lines.append(".end")
    out_path.write_text("\n".join(lines))


def self_test_cylinder():
    """Self-test: synthesize cylinder Y(omega), extract Cauer, verify."""
    sigma_cu = 5.96e7
    mu_0 = 4 * np.pi * 1e-7
    a = 5e-3
    K_SIBC_true = 2 * np.pi * a * np.sqrt(sigma_cu / mu_0)

    # Synthesize cylinder Y(omega) via Foster eigenmodes
    j0 = special.jn_zeros(0, 200)
    j1v = special.j1(j0)
    xi_cyl = j0**2 / (a**2 * mu_0 * sigma_cu)
    residue = 4 * np.pi / mu_0
    omega = np.logspace(0, 8, 100) * 2 * np.pi
    s = 1j * omega
    Y_data = np.sum(residue / (s[:, None] + xi_cyl[None, :]), axis=1)

    print(f"Self-test: cylinder Cu, a = {a*1e3:.1f} mm")
    print(f"  K_SIBC (true) = {K_SIBC_true:.4e}")
    print(f"  Y data shape: {Y_data.shape}, omega range "
          f"[{omega.min()/(2*np.pi):.1e}, {omega.max()/(2*np.pi):.1e}] Hz")
    print()

    # --- Stage 1: K_SIBC from geometry ---
    K_SIBC = K_SIBC_true  # in real workflow, computed from mesh
    print(f"Stage 1: K_SIBC = {K_SIBC:.4e}")

    # --- Stage 2: fit d ---
    print(f"Stage 2: fit d...")
    d, err_d = fit_d(omega, Y_data, K_SIBC)
    print(f"  d = {d:.3e} rad/s ({d/(2*np.pi):.3e} Hz)")
    print(f"  residual = {err_d:.3e}")

    # --- Stage 3: Y_CLN = Y_data - K_SIBC sqrt(s)/(s+d) ---
    Y_warb = Y_warburg(s, K_SIBC, d)
    Y_residual = Y_data - Y_warb

    N_CLN = 30   # higher than digest's 10 to allow bulk to dominate
    print(f"Stage 3: fit N_CLN = {N_CLN} bulk Foster modes (relative-error LS)")
    xi_CLN, c_CLN, err_CLN = fit_bulk_foster(omega, Y_residual, N_CLN,
                                              reference_Y=Y_data)
    print(f"  xi_CLN range: [{xi_CLN.min():.2e}, {xi_CLN.max():.2e}] rad/s")
    print(f"  max rel err of bulk fit: {err_CLN:.3e}")

    # --- Stage 4: diffusive Foster of Warburg ---
    print(f"Stage 4: diffusive Foster of K_SIBC sqrt(s)/(s+d), N_xi = 50")
    xi_w, c_w, max_err_w, _, _ = fit_diffusive(d, N_xi=50,
                                                omega_min=omega.min(),
                                                omega_max=omega.max())
    print(f"  Warburg fit max rel err = {max_err_w:.3e}")

    # --- Verify total reconstruction ---
    Y_reconstruct = (np.sum(c_CLN[None, :] / (s[:, None] + xi_CLN[None, :]), axis=1)
                     + K_SIBC * np.sum(c_w[None, :] / (s[:, None] + xi_w[None, :]),
                                       axis=1))
    err_total = np.abs(Y_reconstruct - Y_data) / np.abs(Y_data)
    print()
    print(f"TOTAL reconstruction error:")
    print(f"  max rel err = {np.max(err_total):.3e}")
    print(f"  median      = {np.median(err_total):.3e}")
    print()

    # Emit SPICE
    out_dir = Path(__file__).parent
    emit_spice_netlist(out_dir / "cylinder_extracted.cir",
                        xi_CLN, c_CLN, xi_w, c_w, K_SIBC, d,
                        comment=f"Self-test extraction from cylinder Foster sum")
    print(f"Saved SPICE: cylinder_extracted.cir")

    # JSON params
    params = {
        "K_SIBC": float(K_SIBC),
        "d_rad_s": float(d),
        "d_Hz": float(d / (2 * np.pi)),
        "N_CLN": int(N_CLN),
        "xi_CLN": xi_CLN.tolist(),
        "c_CLN": c_CLN.tolist(),
        "N_xi": int(len(xi_w)),
        "xi_w": xi_w.tolist(),
        "c_w_real": [float(c.real) for c in c_w],
        "max_rel_err_freq": float(np.max(err_total)),
    }
    (out_dir / "cylinder_extracted_params.json").write_text(json.dumps(params, indent=2))
    print(f"Saved JSON params: cylinder_extracted_params.json")
    return params


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="NPZ file with omega, Y_data arrays")
    parser.add_argument("--K_SIBC", type=float,
                       help="Override K_SIBC (default: estimate from geometry)")
    parser.add_argument("--N_CLN", type=int, default=8)
    parser.add_argument("--N_xi", type=int, default=50)
    parser.add_argument("--output_prefix", default="extracted")
    args = parser.parse_args()

    if args.input is None:
        print("(no --input, running self-test on cylinder Foster)")
        print()
        params = self_test_cylinder()
        sys.exit(0)

    # Load data
    print(f"Loading {args.input}...")
    data = np.load(args.input)
    omega = data["omega"]
    Y_data = data["Y_data"] if "Y_data" in data.files else data["Y"]
    K_SIBC = args.K_SIBC or float(data.get("K_SIBC", 0.0))
    if K_SIBC <= 0:
        raise ValueError("K_SIBC must be provided (via --K_SIBC or in input npz).")

    print(f"K_SIBC = {K_SIBC:.4e}")
    print(f"omega range: [{omega.min()/(2*np.pi):.2e}, "
          f"{omega.max()/(2*np.pi):.2e}] Hz")

    # Stage 2 + 3 + 4
    d, _ = fit_d(omega, Y_data, K_SIBC)
    print(f"d = {d:.3e} rad/s")

    s = 1j * omega
    Y_w = Y_warburg(s, K_SIBC, d)
    xi_CLN, c_CLN, err_CLN = fit_bulk_foster(omega, Y_data - Y_w, args.N_CLN)
    print(f"Bulk Foster fit err = {err_CLN:.3e}")

    xi_w, c_w, max_err_w, _, _ = fit_diffusive(d, N_xi=args.N_xi,
                                                omega_min=omega.min(),
                                                omega_max=omega.max())
    print(f"Warburg diffusive fit err = {max_err_w:.3e}")

    # Save
    out_dir = Path(args.input).parent
    emit_spice_netlist(out_dir / f"{args.output_prefix}.cir",
                        xi_CLN, c_CLN, xi_w, c_w, K_SIBC, d,
                        comment=f"Extracted from {args.input}")
    print(f"Saved: {args.output_prefix}.cir")
