"""Subtract impressed contribution m_imp from NG freq sweep data,
   then Foster fit. Compare with Kameari ladder and ELF Foster.

m_imp(s) = -s sigma V (a^2+b^2) / 48   (analytical, B_0 = 1)

This extracts the genuine eddy current m_ind = m_total - m_imp,
which has clean Foster pole structure.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import csv, json
import numpy as np
from math import pi
from scipy.optimize import least_squares
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az

# Load NG complex freq sweep
ng_csv = Path(__file__).parent / "cuboid_521_vacuum_NG_30f.csv"
data = []
with open(ng_csv) as fp:
    r = csv.reader(fp); next(r)
    for row in r:
        data.append({"freq": float(row[0]),
                     "m_re_NG": float(row[1]), "m_im_NG": float(row[2]),
                     "m_re_ELF": float(row[3]), "m_im_ELF": float(row[4])})

s_arr = np.array([1j * 2*pi*d["freq"] for d in data])
m_NG_total = np.array([d["m_re_NG"] + 1j*d["m_im_NG"] for d in data])
m_ELF = np.array([d["m_re_ELF"] + 1j*d["m_im_ELF"] for d in data])

# Analytical m_imp (B_0 = 1, centered cuboid)
# m_imp = -s * sigma * V * (a^2+b^2) / 48
m_imp = -s_arr * sigma_Cu * V_cond * (ax**2 + ay**2) / 48
print(f"Sample m_imp at f=100 Hz: {m_imp[0]:.4e}")
print(f"Sample m_NG_tot at 100Hz: {m_NG_total[0]:.4e}")
print(f"Predicted match (m_NG ≈ m_imp at low freq): "
      f"NG/imp = {m_NG_total[0] / m_imp[0]:.4f} (should be ~1.0)\n")

# m_ind = m_total - m_imp (genuine eddy current contribution)
m_NG_ind = m_NG_total - m_imp

print(f"After m_imp subtraction:")
for i, f in enumerate([0, 5, 10, 15, 20, 25, 29]):
    d = data[f]
    print(f"  f={d['freq']:.2e} Hz: m_NG_ind = {m_NG_ind[f]:.4e}, |m_NG_ind| = {abs(m_NG_ind[f]):.3e}")

# Convert to chi: chi = m * mu_0 / V (for B_0 = 1)
chi_NG = m_NG_ind * mu0 / V_cond

# Foster fit on M = -chi / (mu_0 s) ... but for chi = m_ind μ_0 / V:
# M_chi = chi / (mu_0 s) = m_ind / (V s)
# Existing fit code uses M = -Y/(mu0 s) where Y is m -> M = -m/(mu_0 s)
# I'll use the same convention: M_NG = -m_NG_ind / (mu_0 s)
M_NG = -m_NG_ind / (mu0 * s_arr)
M_ELF = -m_ELF / (mu0 * s_arr)


def fit_freq(M_arr, N, init_taus=None):
    if init_taus is None:
        init_taus = np.logspace(np.log10(0.5e-6), np.log10(50e-6), N)
    def residual(x):
        taus = np.exp(x[:N])
        a = np.exp(x[N:2*N]) * np.sign(M_arr.real[0])
        fit = np.zeros_like(s_arr, dtype=complex)
        for tk, ak in zip(taus, a):
            fit += ak / (1 + s_arr*tk)
        return np.concatenate([(fit-M_arr).real/np.abs(M_arr),
                               (fit-M_arr).imag/np.abs(M_arr)])
    M0 = abs(M_arr[0])
    x0 = np.concatenate([np.log(init_taus), np.log(np.ones(N)*M0/N)])
    res = least_squares(residual, x0, method='trf', max_nfev=5000)
    taus = np.exp(res.x[:N])
    a = np.exp(res.x[N:2*N]) * np.sign(M_arr.real[0])
    order = np.argsort(-taus)
    M_fit = np.zeros_like(s_arr, dtype=complex)
    for tk, ak in zip(taus[order], a[order]):
        M_fit += ak / (1 + s_arr*tk)
    rms = np.sqrt(np.mean(np.abs(M_fit-M_arr)**2 / np.abs(M_arr)**2)) * 100
    return taus[order], a[order], rms


print("\n" + "=" * 78)
print("FOSTER FIT (after subtracting m_imp from NG)")
print("=" * 78)
for N in [4, 6, 8]:
    taus_ng, a_ng, rms_ng = fit_freq(M_NG, N)
    taus_elf, a_elf, rms_elf = fit_freq(M_ELF, N)
    print(f"\n--- N = {N} ---")
    print(f"  NG (m_ind): tau = {[f'{t*1e6:.3f}' for t in taus_ng]} us, rms = {rms_ng:.3f}%")
    print(f"  ELF:        tau = {[f'{t*1e6:.3f}' for t in taus_elf]} us, rms = {rms_elf:.3f}%")
    print(f"  ratio tau_lead (NG/ELF) = {taus_ng[0]/taus_elf[0]:.3f}")

# 4-stage final
taus_ng4, a_ng4, rms_ng4 = fit_freq(M_NG, 4)
taus_elf4, a_elf4, rms_elf4 = fit_freq(M_ELF, 4)

print("\n" + "=" * 78)
print("CRITICAL COMPARISON (4-stage)")
print("=" * 78)
print(f"  Kameari Stage 0 tau     = 104.37 us (AIR_SCALE=5, accumulated)")
print(f"  Foster-fit NG (m_ind)   = {taus_ng4[0]*1e6:.3f} us (this script, AIR_SCALE=5)")
print(f"  Foster-fit ELF          = {taus_elf4[0]*1e6:.3f} us")
print()
ng_kam = taus_ng4[0]*1e6 / 104.37
ng_elf = taus_ng4[0]*1e6 / taus_elf4[0]/1e6
print(f"  NG / Kameari = {taus_ng4[0]*1e6 / 104.37:.3f}")
print(f"  NG / ELF     = {taus_ng4[0]*1e6 / (taus_elf4[0]*1e6):.3f}")
print()
if abs(taus_ng4[0]*1e6 - 104.37) / 104.37 < 0.10:
    print("  >>> Kameari ladder is FAITHFUL to air-box-truncated NG freq sweep <<<")
    print("  >>> Air-box != BEM is the only issue; Kameari iteration is correct <<<")
elif abs(taus_ng4[0]*1e6 - taus_elf4[0]*1e6) / (taus_elf4[0]*1e6) < 0.10:
    print("  >>> NG freq sweep AGREES with ELF; Kameari iteration produces wrong tau <<<")
    print("  >>> Genuine Kameari iteration breakdown <<<")
else:
    print(f"  >>> THREE-WAY DISAGREEMENT: NG={taus_ng4[0]*1e6:.2f}, "
          f"Kameari=104.37, ELF={taus_elf4[0]*1e6:.2f} us <<<")

out = {
    "method": "NG freq sweep with m_imp subtracted, Foster fit",
    "AIR_SCALE": 5, "N_stages": 4,
    "NG_4stage": {"taus_us": (taus_ng4*1e6).tolist(), "a": a_ng4.tolist(), "rms_pct": float(rms_ng4)},
    "ELF_4stage": {"taus_us": (taus_elf4*1e6).tolist(), "a": a_elf4.tolist(), "rms_pct": float(rms_elf4)},
    "Kameari_stage0_tau_us": 104.37,
    "comparison": {
        "tau_lead_NG_us": float(taus_ng4[0]*1e6),
        "tau_lead_ELF_us": float(taus_elf4[0]*1e6),
        "tau_Kameari_us": 104.37,
        "ratio_NG_over_Kameari": float(taus_ng4[0]*1e6 / 104.37),
        "ratio_NG_over_ELF": float(taus_ng4[0] / taus_elf4[0]),
    },
}
out_path = Path(__file__).parent / "cuboid_521_NG_minus_mimp_foster.json"
out_path.write_text(json.dumps(out, indent=2))
print(f"\nResults: {out_path}")
