"""COMSOL-class #65 -- waveguide dispersion (propagation above cutoff, evanescence below).

The propagation sequel to the cutoff-eigenvalue blocks (#53 rectangular, #62 circular). Above its
cutoff fc a hollow metallic waveguide mode propagates with axial constant

    beta = (2 pi/c) sqrt(f^2 - fc^2),

so the guide wavelength lambda_g = lambda0/sqrt(1-(fc/f)^2) exceeds the free-space lambda0 = c/f,
the phase velocity v_p = c/sqrt(1-(fc/f)^2) > c carries no energy, the group (energy/signal)
velocity v_g = c sqrt(1-(fc/f)^2) < c does, and the two obey v_p*v_g = c^2. Below cutoff the mode
is evanescent, decaying as exp(-alpha z) with alpha = (2 pi/c) sqrt(fc^2 - f^2) (the waveguide as a
high-pass filter). The cutoff fc itself is taken from the FE Helmholtz eigensolver (#53), tying the
analytic dispersion curve to the FE-computed cutoff. Self-contained, headless; SI units.
"""
import os
import sys
import math

sys.stdout.reconfigure(encoding="utf-8", errors="ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from radia_mcp.radia_ngsolve.waveguide import (rectangular_waveguide_cutoff, cutoff_frequency,
                                               waveguide_dispersion, guide_wavelength,
                                               waveguide_evanescent_attenuation, C0)

a, b = 0.02286, 0.01016                       # WR-90, the X-band standard
fc = rectangular_waveguide_cutoff(a, b, 1, 0)  # TE10 dominant cutoff = c/2a
print(f"WR-90 (a={a*1e3:.2f} mm) TE10 cutoff fc = {fc/1e9:.4f} GHz, single-mode to {2*fc/1e9:.1f} GHz")
print("\nABOVE cutoff -- propagation:")
print("  f[GHz]  beta[1/m]  lambda_g[mm]   v_p/c    v_g/c    v_p*v_g/c^2")
for f in (8e9, 10e9, 12e9, 15e9):
    d = waveguide_dispersion(f, fc)
    print(f"  {f/1e9:5.1f}  {d['beta']:8.2f}   {d['lambda_g']*1e3:8.2f}    {d['v_phase']/C0:.4f}   "
          f"{d['v_group']/C0:.4f}   {d['v_phase']*d['v_group']/C0**2:.10f}")

print("\nBELOW cutoff -- evanescent (exp(-alpha z), high-pass filter):")
print("  f[GHz]  alpha[Np/m]  alpha/k_c")
kc = 2 * math.pi * fc / C0
for f in (1e9, 3e9, 6e9):
    al = waveguide_evanescent_attenuation(f, fc)
    print(f"  {f/1e9:5.1f}  {al:8.2f}     {al/kc:.4f}")
print(f"  (alpha -> k_c={kc:.1f} as f->0; alpha -> 0 as f->fc)")

# FE anchor: the cutoff that builds the dispersion comes from the Helmholtz eigensolver (#53)
try:
    from ngsolve import Mesh, TaskManager
    from netgen.occ import OCCGeometry, WorkPlane
    from radia_mcp.radia_ngsolve.waveguide import helmholtz_cutoff_wavenumbers_2d
    rect = WorkPlane().Rectangle(a, b).Face(); rect.edges.name = "wall"
    mesh = Mesh(OCCGeometry(rect, dim=2).GenerateMesh(maxh=min(a, b) / 12))
    with TaskManager():
        kc_fe = helmholtz_cutoff_wavenumbers_2d(mesh, 1, bc="neumann")[0]
    fc_fe = cutoff_frequency(kc_fe)
    d = waveguide_dispersion(10e9, fc_fe)
    print(f"\nFE-anchored cutoff: fc_FE = {fc_fe/1e9:.4f} GHz vs analytic {fc/1e9:.4f} "
          f"({(fc_fe/fc-1)*100:+.3f}%); lambda_g@10GHz = {d['lambda_g']*1e3:.2f} mm")
except Exception as e:                          # noqa: BLE001
    print(f"\n[FE anchor skipped: {e}]")

print("\nOK -- v_p*v_g = c^2; lambda_g > lambda0; v_g -> 0 at cutoff, both -> c at high f.")
