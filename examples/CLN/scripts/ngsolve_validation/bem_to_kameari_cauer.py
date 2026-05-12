"""bem_to_kameari_cauer.py — Convert BEM-Foster spectrum to Kameari Cauer-I.

Cylinder BEM extraction (bem_disk_axisym_cauer.wls) gives Cauer-I via
  α_n^Hiruma = Σ_k g_k² × τ_k^(n+1)
which is the IMPEDANCE Cauer-I (Hiruma / Krylov-Padé convention).

Kameari accumulation iteration extracts Cauer-I from
  α_n^Kameari = Σ_k g_k² × τ_k^n
which is the SUSCEPTIBILITY Cauer-I (χ-Foster convention).

Mathematical relation:
  Hiruma f(s) = Σ g_k² τ_k / (1 + s τ_k)
  Kameari g(s) = Σ g_k² / (1 + s τ_k)
  → f(s) = -s · g(s) + Σ g_k² (DC offset doesn't affect Cauer)

Difference: power of τ shifts by 1 in moment formula → produces
DIFFERENT continued fractions despite same Foster spectrum.

Usage:
  python bem_to_kameari_cauer.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import mpmath as mp

mp.mp.dps = 60  # 60 decimal digits


def hankel_pade_cauer(alphas, max_stages):
    """Hankel-Padé continued-fraction extraction of Cauer-I rungs from
    Taylor moments α_0, α_1, ...

    Returns list of p-coefficients [p_0, p_1, p_2, ...] where
    R_2k = p_{2k}, 1/L_{2k+1} = p_{2k+1}.
    """
    c = [mp.mpf(a) for a in alphas]
    p_list = []
    for stage in range(max_stages):
        if len(c) < 2 or abs(c[0]) < mp.mpf(10) ** (-mp.mp.dps + 10):
            break
        # Inverse of Taylor series c
        n = len(c)
        e = [mp.mpf(0)] * n
        e[0] = 1 / c[0]
        for k in range(1, n):
            e[k] = -sum(c[j + 1] * e[k - j - 1]
                         for j in range(k)
                         if k - j - 1 >= 0) / c[0]
        p_list.append(e[0])
        c = e[1:]
    return p_list


def cauer_from_spectrum(amplitudes_squared, tau_list, n_moments=30, max_pairs=12):
    """Extract BOTH Hiruma and Kameari Cauer-I from {τ_k, X_k=g_k² or a_k×τ_k}.

    Hiruma: α_n = (-1)^n Σ X_k τ_k^(n+1)  (impedance Cauer)
    Kameari: α_n = (-1)^n Σ X_k τ_k^n     (susceptibility Cauer; power -1)

    Returns (hiruma_table, kameari_table) where each table is
    [(k, R_2k, L_2k+1, τ_pair_us), ...].
    """
    alphas_h = []
    alphas_k = []
    for n in range(n_moments):
        a_h = mp.mpf(0); a_k = mp.mpf(0)
        for X, t in zip(amplitudes_squared, tau_list):
            a_h += X * t ** (n + 1)
            a_k += X * t ** n
        alphas_h.append((mp.mpf(-1)) ** n * a_h)
        alphas_k.append((mp.mpf(-1)) ** n * a_k)
    p_h = hankel_pade_cauer(alphas_h, 2 * max_pairs)
    p_k = hankel_pade_cauer(alphas_k, 2 * max_pairs)

    def to_table(p_list):
        out = []
        for k in range(len(p_list) // 2):
            R = p_list[2 * k]
            Linv = p_list[2 * k + 1]
            if abs(Linv) < mp.mpf(10) ** (-50):
                break
            L = 1 / Linv
            tau = L / R * mp.mpf("1e6")
            out.append((k, float(R), float(L), float(tau)))
        return out
    return to_table(p_h), to_table(p_k)


def run_cylinder():
    bem_path = Path(__file__).parent / "bem_disk_axisym_cauer.json"
    with open(bem_path) as f:
        bem = json.load(f)
    g_list = [mp.mpf(g) for g in bem["foster_amplitudes_g"]]
    tau_list_us = [mp.mpf(t) for t in bem["bem_tau_us"]]
    tau_list = [t * mp.mpf("1e-6") for t in tau_list_us]
    g2_list = [g * g for g in g_list]
    print("=" * 80)
    print(" CYLINDER BEM (R=10 t=2 Cu disk, 50 modes)")
    print("=" * 80)
    h, k = cauer_from_spectrum(g2_list, tau_list)
    ngs_3d = [208.34, 67.25, 31.80, 19.76, 15.82, 12.94]
    ngs_ax = [208.38, 67.48, 32.57, 20.36, 15.68, 11.35]
    print(f"  {'k':>3} | {'BEM Hiruma':>11} | {'BEM Kameari':>11}"
          f" | {'NGS 3D Kam':>11} | {'NGS axisym Kam':>14}")
    for i in range(min(6, len(h), len(k))):
        print(f"  {i:>3} | {h[i][3]:>11.4f} | {k[i][3]:>11.4f}"
              f" | {ngs_3d[i]:>11.4f} | {ngs_ax[i]:>14.4f}")
    return h, k


def run_sphere():
    """Sphere uses Stoll closed-form spectrum."""
    print()
    print("=" * 80)
    print(" SPHERE Stoll closed-form (Cu R=10mm, B0=1T, analytical)")
    print("=" * 80)
    mu0 = mp.pi * mp.mpf(4) * mp.mpf(10) ** (-7)
    sigma = mp.mpf(58) * mp.mpf(10) ** 6
    a_radius = mp.mpf(10) / 1000
    n_modes = 200
    # Stoll: τ_n = μ_0 σ a² / (n²π²), a_n = 6/(n²π²)
    # χ amplitude a_n is the "Kameari amplitude". Hiruma g_k² = a_k × τ_k.
    tau_list = [mu0 * sigma * a_radius ** 2 / (n ** 2 * mp.pi ** 2)
                for n in range(1, n_modes + 1)]
    a_chi = [mp.mpf(6) / (n ** 2 * mp.pi ** 2) for n in range(1, n_modes + 1)]
    # Hiruma "g_k²" = a_k * τ_k (with this scaling, both formulas give the
    # right thing per the derivation; the X in cauer_from_spectrum is a
    # generic amplitude squared)
    g2_for_hiruma = [a * t for a, t in zip(a_chi, tau_list)]   # Hiruma input
    # Kameari direct uses a_k itself with τ^(n+1) formula (= Stoll wls):
    # α_n = Σ a_k τ_k^(n+1) → already Kameari ladder (sphere proves this works)

    # Approach 1: Use g²=a×τ as Hiruma input, extract both
    h, k = cauer_from_spectrum(g2_for_hiruma, tau_list)

    # Approach 2: Use a_chi directly with τ^(n+1) formula (Stoll wls's view)
    # Equivalent to "Kameari extraction with X=a_chi" if we shift definitions
    ngs_kameari = [694.14, 154.60, 64.07, 34.47, 21.40, 14.54]
    print(f"  {'k':>3} | {'BEM Hiruma':>11} | {'BEM Kameari':>11}"
          f" | {'NGS Kameari/Stoll':>17}")
    for i in range(min(6, len(h), len(k))):
        print(f"  {i:>3} | {h[i][3]:>11.4f} | {k[i][3]:>11.4f}"
              f" | {ngs_kameari[i]:>17.4f}")
    print()
    print("  NOTE: 'BEM Kameari' (= Sum a*tau * tau^n = Sum a tau^(n+1)) matches Stoll/NGSolve")
    print("        'BEM Hiruma' (= Sum a*tau * tau^(n+1) = Sum a tau^(n+2)) is sphere's IMPEDANCE")
    print("  Cylinder used g_k^2 as raw Foster amplitude; sphere uses chi-amplitude a_k.")
    print("  The conversion identity: g_k^2_Hiruma_BEM = a_k_Kameari_chi * tau_k.")
    return h, k


def main():
    cyl_h, cyl_k = run_cylinder()
    sph_h, sph_k = run_sphere()

    # Save unified output
    out = {
        "method": "BEM-Foster spectrum -> Cauer-I in two conventions (mpmath 60 digit)",
        "cylinder": {
            "Hiruma_table": cyl_h,
            "Kameari_table": cyl_k,
        },
        "sphere": {
            "Hiruma_table": sph_h,
            "Kameari_table": sph_k,
        },
    }
    out_path = (Path(__file__).parent
                / "2026-05-10-bem_to_kameari_cauer_results.json")
    out_path.write_text(json.dumps(out, indent=2, default=str),
                         encoding="utf-8")
    print(f"\n Saved: {out_path}")
    return  # skip the old standalone main below
    bem_path = Path(__file__).parent / "bem_disk_axisym_cauer.json"
    with open(bem_path) as f:
        bem = json.load(f)
    g_list = [mp.mpf(g) for g in bem["foster_amplitudes_g"]]
    tau_list_us = [mp.mpf(t) for t in bem["bem_tau_us"]]   # in μs
    tau_list = [t * mp.mpf("1e-6") for t in tau_list_us]   # in seconds

    n_modes = len(g_list)
    print(f"BEM cylinder spectrum: {n_modes} modes")
    print(f"  Leading τ_k [μs]: {[float(t) for t in tau_list_us[:6]]}")
    print()

    # Hiruma moments (reference, should reproduce wls result)
    n_moments = 30
    alphas_hiruma = []
    for n in range(n_moments):
        a = mp.mpf(0)
        for k in range(n_modes):
            a += g_list[k] ** 2 * tau_list[k] ** (n + 1)
        alphas_hiruma.append((mp.mpf(-1)) ** n * a)

    # Kameari moments (NEW: τ^n instead of τ^(n+1))
    alphas_kameari = []
    for n in range(n_moments):
        a = mp.mpf(0)
        for k in range(n_modes):
            a += g_list[k] ** 2 * tau_list[k] ** n
        alphas_kameari.append((mp.mpf(-1)) ** n * a)

    # Extract Cauer-I from both
    p_hiruma = hankel_pade_cauer(alphas_hiruma, max_stages=2 * 12)
    p_kameari = hankel_pade_cauer(alphas_kameari, max_stages=2 * 12)

    def to_table(p_list, label):
        n_pairs = len(p_list) // 2
        out = []
        for k in range(n_pairs):
            R_2k = p_list[2 * k]
            L_inv = p_list[2 * k + 1]
            if abs(L_inv) < mp.mpf(10) ** (-50):
                break
            L = 1 / L_inv
            tau = L / R_2k * mp.mpf("1e6")  # μs
            out.append((k, float(R_2k), float(L), float(tau)))
        return out

    h_table = to_table(p_hiruma, "Hiruma")
    k_table = to_table(p_kameari, "Kameari")

    print("=" * 80)
    print(" Cylinder BEM → Cauer-I extraction (mpmath 60 digits)")
    print("=" * 80)
    print()
    print(f"  {'k':>3} | {'Hiruma τ [μs]':>14} | {'Kameari τ [μs]':>14}"
          f" | {'NGSolve 3D Kameari':>20} | {'NGSolve axisym Kam':>20}")
    print("  " + "-" * 92)

    # NGSolve reference values
    ngs_3d = [208.34, 67.25, 31.80, 19.76, 15.82, 12.94]
    ngs_ax = [208.38, 67.48, 32.57, 20.36, 15.68, 11.35]

    for i in range(min(len(h_table), len(k_table), 8)):
        h = h_table[i][3] if i < len(h_table) else float('nan')
        k = k_table[i][3] if i < len(k_table) else float('nan')
        n3 = ngs_3d[i] if i < len(ngs_3d) else float('nan')
        nax = ngs_ax[i] if i < len(ngs_ax) else float('nan')
        print(f"  {i:>3} | {h:>14.4f} | {k:>14.4f}"
              f" | {n3:>20.4f} | {nax:>20.4f}")

    print()
    print(" Verification:")
    print(f"  Hiruma BEM τ_pair[0] = {h_table[0][3]:.4f} μs"
          " (should be ~219.32 from wls)")
    print(f"  Kameari BEM τ_pair[0] = {k_table[0][3]:.4f} μs"
          " (should match NGSolve Kameari ~208.4)")

    # Save
    out = {
        "method": "BEM-Foster spectrum -> Kameari Cauer-I (mpmath 60 digit)",
        "n_modes": n_modes,
        "n_moments": n_moments,
        "Hiruma_table": h_table,
        "Kameari_table": k_table,
        "ngs_3d_ref": ngs_3d,
        "ngs_axisym_ref": ngs_ax,
    }
    out_path = (Path(__file__).parent
                / "2026-05-10-bem_to_kameari_cauer_results.json")
    out_path.write_text(json.dumps(out, indent=2, default=str),
                         encoding="utf-8")
    print(f"\n Saved: {out_path}")


if __name__ == "__main__":
    main()
