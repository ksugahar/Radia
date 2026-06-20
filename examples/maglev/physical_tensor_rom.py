"""physical_tensor_rom.py -- the EXTERIOR-MATCHED physical polarizability tensor
as a passive, stable LTI (the "physical Stoll spectrum -> CLN" route).

Frontier "physical-tensor ROM" of the radia-maglev stack.  The mixed-Galerkin
bulk Foster (bulk_foster_via_eigen / bulk_foster_vector_via_eigen) uses the
INTERIOR-PEC eigenmodes -- a model, not the physical exterior-matched spectrum.
The verified PHYSICAL polarizability tensor is the per-frequency 3D HCurl solve
ellipsoid/ellipsoid_alpha_tensor_3d.py (it carries the air reaction dipole, the
lift / Re[alpha] part).  This script turns that per-frequency solve into a
passive, stable LTI:

    alpha(s) ~ alpha_inf + sum_k g_k/(1 + s tau_k),   g_k >= 0, tau_k > 0,

via the AAA + NNLS recipe in mixed_galerkin/rom_fit.py.  AAA discovers the
DOMINANT real poles (the physical Stoll decay times), a log-spaced filler
captures the tail, and NNLS makes the residues passive.

Why a sample fit and not a Kameari + Kelvin eigen-accumulation: the 3D HCurl
Kameari + Kelvin accumulation structurally BREAKS DOWN on the
isolated-conductor-in-vacuum problem (L_n sign flip at stage 1; see
research_cln/ngsolve_validation/cuboid_521_kameari_kelvin_v15_canonical.py).
This route sidesteps that by building on the verified per-frequency solve.

Two modes:
  (default) analytic sphere Stoll alpha(s) -> ROM.  Pure numpy, ~1 s.  The
            AAA poles reproduce the analytic Stoll tau_n = mu0 sigma a^2/(n pi)^2
            to ~0.00 % and the LTI matches alpha(s) to < 0.2 % over 1 Hz..1 GHz.
  --fem     triaxial ellipsoid: sample the verified 3D HCurl tensor over
            frequency (per axis) -> three ROMs -> a diagonal MIMO LTI.  HEAVY
            (CompactAMS+COCR, ~minutes per (freq, axis); tens of minutes total).

Run:  python physical_tensor_rom.py            (analytic sphere, fast)
      python physical_tensor_rom.py --fem      (triaxial FEM, slow)
"""
from __future__ import annotations

import argparse
import cmath
import json
import math
import os

import numpy as np

from radia.maglev.mixed_galerkin import (
    passive_foster_fit, diagonal_tensor_state_space,
)

HERE = os.path.dirname(os.path.abspath(__file__))
MU0 = 4 * math.pi * 1e-7
SIGMA_CU = 5.8e7
A_SPH = 5e-3


def sphere_alpha_causal(omega, a=A_SPH, sigma=SIGMA_CU):
    """Causal-convention sphere polarizability alpha(j omega) (LHP poles).

    alpha(s)/(4 pi a^3) = -1/2 + sum_n (3/(n pi)^2)/(1 + s tau_n),
    tau_n = mu0 sigma a^2/(n pi)^2  (the analytic Stoll spectrum).
    """
    delta = math.sqrt(2.0 / (omega * MU0 * sigma))
    x = (1 + 1j) * a / delta
    e = cmath.exp(2j * x)
    cot = 1j * (e + 1.0) / (e - 1.0)
    G = -0.5 * (1.0 - 3.0 / x**2 + 3.0 / x * cot)
    return 4 * math.pi * a**3 * G.conjugate()


def run_sphere():
    n = np.arange(1, 6)
    tau_anal = MU0 * SIGMA_CU * A_SPH**2 / (n * math.pi) ** 2
    f = np.logspace(0, 9, 400)
    s = 1j * 2 * math.pi * f
    Y = np.array([sphere_alpha_causal(2 * math.pi * fi) for fi in f])

    rom = passive_foster_fit(s, Y, n_filler=20)
    dom = np.sort(rom.dominant_tau)[::-1]

    print("Analytic sphere Stoll alpha(s) -> passive Foster ROM")
    print(f"  Cu sphere a = {A_SPH*1e3:.1f} mm, sigma = {SIGMA_CU:.2e}")
    print(f"  ROM: {rom.n_states} states, band fit = {rom.band_fit_relerr:.2e}, "
          f"passive = {np.all(rom.g_n >= 0)}")
    print(f"  alpha_inf = {rom.alpha_inf*1e9:+.3f} mm^3  "
          f"(analytic -2 pi a^3 = {-2*math.pi*A_SPH**3*1e9:+.3f})")
    print("  dominant pole (Stoll) vs analytic tau_n:")
    print("    k   tau_ROM[us]   tau_anal[us]   rel.err")
    for k in range(min(len(dom), 4)):
        err = (dom[k] - tau_anal[k]) / tau_anal[k] * 100
        print(f"    {k}   {dom[k]*1e6:9.4f}   {tau_anal[k]*1e6:11.4f}   {err:+6.2f}%")

    # isotropic diagonal LTI (sphere: same ROM on all 3 axes)
    A, B, C, D, ns = diagonal_tensor_state_space([rom, rom, rom])
    print(f"  diagonal MIMO LTI: {ns} states, D_diag = "
          f"{np.round(np.diag(D)*1e9, 1)} mm^3")

    out = {
        "mode": "analytic_sphere_stoll",
        "geometry": "Cu sphere a=5mm", "sigma": SIGMA_CU,
        "rom": {
            "n_states": rom.n_states,
            "band_fit_relerr": rom.band_fit_relerr,
            "alpha_inf_m3": rom.alpha_inf,
            "tau_n_us": (rom.tau_n * 1e6).tolist(),
            "g_n_m3": rom.g_n.tolist(),
            "dominant_tau_us": (dom * 1e6).tolist(),
        },
        "analytic_stoll_tau_us": (tau_anal * 1e6).tolist(),
        "mimo_n_states": int(ns),
        "scope_note": ("exterior-matched physical alpha(s) as a passive LTI; AAA "
                       "poles = physical Stoll decay times.  Kameari+Kelvin "
                       "eigen-accumulation breaks down on this isolated-conductor "
                       "problem, so the LTI is fit from the per-frequency solve."),
    }
    with open(os.path.join(HERE, "physical_tensor_rom_sphere.json"), "w") as fp:
        json.dump(out, fp, indent=2)
    print("  wrote physical_tensor_rom_sphere.json")


def run_fem():
    """Triaxial ellipsoid: verified 3D HCurl tensor over frequency -> diagonal LTI."""
    import sys
    sys.path.insert(0, os.path.join(HERE, "ellipsoid"))
    import ellipsoid_alpha_tensor_3d as ET3
    import ellipsoid_alpha_tensor as ET
    from ngsolve import TaskManager  # noqa: F401  (ET3 wraps the solve internally)

    semi = (5e-3, 3e-3, 1.5e-3)
    freqs = np.logspace(2, 5, 9)
    N = ET.demag_tensor(semi)
    print(f"Triaxial ellipsoid {tuple(x*1e3 for x in semi)} mm -- FEM tensor -> ROM")
    mesh = ET3.build_mesh(semi, box=0.022, maxh=0.012, maxh_cond=5e-4,
                          maxh_shell=7e-4)
    print(f"  mesh ne = {mesh.ne}")

    roms = {}
    rec = {"semi_m": list(semi), "N_demag": N.tolist(), "freqs_Hz": freqs.tolist(),
           "alpha_re_m3": {}, "alpha_im_m3": {}, "rom": {}}
    for d in "xyz":
        s = 1j * 2 * math.pi * freqs
        Y = np.array([ET3.alpha_tensor_component(mesh, 2 * math.pi * f, d)
                      for f in freqs])
        rec["alpha_re_m3"][d] = Y.real.tolist()
        rec["alpha_im_m3"][d] = Y.imag.tolist()
        # passive_foster_fit wants the CAUSAL/passive convention (Im < 0 on the
        # j omega axis: g_k/(1+s tau_k) with g_k>=0).  The 3D HCurl solve returns
        # the physics convention (Im > 0), so conjugate before fitting.
        rom = passive_foster_fit(s, np.conj(Y), n_filler=12)
        roms[d] = rom
        rec["rom"][d] = {"tau_n_us": (rom.tau_n * 1e6).tolist(),
                         "g_n_m3": rom.g_n.tolist(),
                         "alpha_inf_m3": rom.alpha_inf,
                         "dominant_tau_us": (rom.dominant_tau * 1e6).tolist(),
                         "n_states": rom.n_states,
                         "band_fit_relerr": rom.band_fit_relerr}
        print(f"  axis {d}: {rom.n_states} states, band fit {rom.band_fit_relerr:.2e}, "
              f"dominant tau (us) {np.round(rom.dominant_tau[:3]*1e6, 2)}")

    A, B, C, D, ns = diagonal_tensor_state_space([roms["x"], roms["y"], roms["z"]])
    rec["mimo"] = {"n_states": int(ns), "P": 3,
                   "D_diag_m3": np.diag(D).tolist()}
    print(f"  diagonal MIMO LTI: {ns} states, D_diag (mm^3) "
          f"{np.round(np.diag(D)*1e9, 1)}")
    with open(os.path.join(HERE, "physical_tensor_rom_fem.json"), "w") as fp:
        json.dump(rec, fp, indent=2)
    print("  wrote physical_tensor_rom_fem.json")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fem", action="store_true",
                    help="triaxial ellipsoid FEM tensor -> ROM (slow, ~tens of min)")
    args = ap.parse_args()
    if args.fem:
        run_fem()
    else:
        run_sphere()


if __name__ == "__main__":
    main()
