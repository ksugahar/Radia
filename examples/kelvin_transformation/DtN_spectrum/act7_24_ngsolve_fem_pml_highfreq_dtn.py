# -*- coding: utf-8 -*-
"""
act7_24_ngsolve_fem_pml_highfreq_dtn.py  (Act 7 -- the radiating extension; the PML column made REAL)
=====================================================================================================
The symmetric partner of act7_23_ngsbem_highfreq_dtn: it puts the ACTUAL 3-D NGSolve FEM + PML on the
SAME DtN-spectral yardstick.  In act7_22_dtn_spectrum_consolidated the "PML" column is a radial-FE
PROXY (helm_pml_dtn); here a genuine 3-D NGSolve Helmholtz solve with NGSolve's native radial PML
BUILDS the exterior Dirichlet-to-Neumann operator and we MEASURE its per-multipole defect against
wave_dtn, plus the reflection it implies.

WHY: high-freq IS a studied regime (the act7 radiating thread), and the user's point is that NGSolve
(FEM + PML) AND ngsolve.bem are both working tools there.  act7_23 made the BEM column real; this
makes the FEM+PML column real -- so both halves are on the yardstick, not proxies.

METHOD (exterior acoustic Helmholtz, truncation sphere radius a, wavenumber kappa = ka):
  - mesh ONLY the exterior shell  a <= r <= R_out  (the interior r<a is the un-meshed exterior-DtN
    domain seen from Gamma=r=a);
  - NGSolve native RADIAL PML stretching r>a over the whole shell (mesh.SetPML(pml.Radial(...)));
    outer sphere r=R_out is a Dirichlet wall (the stretched field has decayed);
  - assemble the Helmholtz form  A = int (grad u . grad v - k^2 u v) dx  (complex, PML-stretched);
  - exterior DtN by CONSISTENT FLUX: set the Dirichlet datum u_D = Y_n on Gamma, solve the interior
    (shell) DOFs, and read the weak Neumann reaction  b = A u  on the constrained Gamma DOFs.  For a
    pure spherical-harmonic mode  b = Lambda_n M_Gamma u_D, so
        Lambda_n^FEMPML = - (u_D^T b) / (u_D^T M_Gamma u_D)
    -- the BILINEAR (non-conjugating) Rayleigh quotient (the DtN is complex-SYMMETRIC, not Hermitian,
    so g^T A g is the right form -- a Hermitian g^H A g flips the real part), and the leading MINUS is
    the shell's inner-boundary outward normal (-r at Gamma) vs the exterior normal (+r).
  - exact = radia.open_boundary.wave_dtn(n, ka) = ka h_n^(1)'(ka)/h_n^(1)(ka) (complex / radiating).
  - reflection  R_n = |Lambda_h - Lambda_exact| / |Lambda_h - conj(Lambda_exact)|  (the same mismatch
    as the DtN defect; see act7_22's reflection view).

VERIFIED (asserted): monopole / dipole / quadrupole DtN reldef < 5e-3 vs wave_dtn (coarser than the
~1e-5 BEM of act7_23 because a 3-D FEM+PML carries BOTH volume-discretisation AND PML-truncation
error -- but it is a genuine 3-D solve reproducing the COMPLEX DtN), and low reflection on the
propagating modes.  Needs ngsolve + netgen.occ (6.2.2604+).
"""
import os
import json
import math
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from ngsolve import (Mesh, H1, BilinearForm, GridFunction, TaskManager, CF, z,
                     grad, dx, ds, pml)
from netgen.occ import Sphere, Pnt, OCCGeometry
import radia.open_boundary as ob

N_FAIL = 0


def check(name, cond, detail=""):
    global N_FAIL
    if not cond:
        N_FAIL += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}{('  -- ' + detail) if detail else ''}")


KA = 2.0
A_RAD = 1.0
R_OUT = 3.0
print("=" * 80)
print(f" act7_24_ngsolve_fem_pml_highfreq_dtn : 3-D NGSolve FEM+PML on the DtN yardstick (a=1, ka={KA})")
print("=" * 80)

res = {}
with TaskManager():
    # exterior shell a <= r <= R_out; name inner face Gamma, outer face wall (by area)
    shell = Sphere(Pnt(0, 0, 0), R_OUT) - Sphere(Pnt(0, 0, 0), A_RAD)
    for f in shell.faces:
        af = f.mass
        f.name = "Gamma" if abs(af - 4 * math.pi * A_RAD**2) < abs(af - 4 * math.pi * R_OUT**2) else "wall"
    mesh = Mesh(OCCGeometry(shell).GenerateMesh(maxh=0.4))
    mesh.Curve(3)
    # NGSolve native radial PML: complex-stretch the shell (r > a)
    mesh.SetPML(pml.Radial(origin=(0, 0, 0), rad=A_RAD, alpha=1j), definedon="default")

    fes = H1(mesh, order=3, complex=True, dirichlet="Gamma|wall")
    u, v = fes.TnT()
    Aform = BilinearForm(fes, check_unused=False)
    Aform += (grad(u) * grad(v) - KA * KA * u * v) * dx
    Aform.Assemble()
    mgform = BilinearForm(fes, check_unused=False)          # surface mass on Gamma only
    mgform += u.Trace() * v.Trace() * ds(definedon=mesh.Boundaries("Gamma"))
    mgform.Assemble()
    Ainv = Aform.mat.Inverse(fes.FreeDofs(), inverse="umfpack")   # factor ONCE, reuse per mode

    def npvec(expr):
        w = gfu.vec.CreateVector()
        w.data = expr
        return w.FV().NumPy().copy()

    print(f"  exterior shell [{A_RAD},{R_OUT}] H1(order=3,complex): {fes.ndof} DOF, native radial PML")
    print("  n  mode         FEM+PML DtN         exact wave_dtn      reldef d_n   reflection R_n")
    MODES = [(0, CF(1.0), "monopole"), (1, z, "dipole"),
             (2, 3 * z * z - 1, "quadrupole"), (3, 5 * z * z * z - 3 * z, "octupole")]
    for n, cf, lbl in MODES:
        gfu = GridFunction(fes)
        gfu.Set(cf, definedon=mesh.Boundaries("Gamma"))           # Dirichlet datum on Gamma
        gfu.vec.data += Ainv * (-Aform.mat * gfu.vec)             # solve the shell interior DOFs
        g = gfu.vec.FV().NumPy().copy()
        b = npvec(Aform.mat * gfu.vec)                            # weak Neumann reaction
        mgg = npvec(mgform.mat * gfu.vec)
        lam = complex(-(g @ b) / (g @ mgg))                      # bilinear Rayleigh, inner-normal sign
        exact = complex(ob.wave_dtn(n, KA))
        d = abs(lam - exact) / abs(exact)
        R = abs(lam - exact) / abs(lam - np.conj(exact))
        res[lbl] = {"n": n, "FEM_PML": [lam.real, lam.imag],
                    "exact": [exact.real, exact.imag], "reldef": d, "reflection": R}
        print(f"  {n}  {lbl:10s}  {lam:+.4f}  {exact:+.4f}   {d:.2e}    {R:.2e}")

low = ("monopole", "dipole", "quadrupole")
check("3-D NGSolve FEM+PML reproduces the COMPLEX exterior DtN (monopole/dipole/quadrupole reldef < 5e-3)",
      max(res[l]["reldef"] for l in low) < 5e-3,
      f"max {max(res[l]['reldef'] for l in low):.1e}")
check("3-D NGSolve FEM+PML low-reflection on the propagating modes (R < 5e-2)",
      max(res[l]["reflection"] for l in low) < 5e-2,
      f"max {max(res[l]['reflection'] for l in low):.1e}")
check("octupole (n=3, more radial structure) still resolved (reldef < 2e-2)",
      res["octupole"]["reldef"] < 2e-2, f"{res['octupole']['reldef']:.1e}")

print("\n  => the high-freq 'PML' column is REAL too: a genuine 3-D NGSolve Helmholtz solve with native")
print("     radial PML reproduces the complex exterior DtN (~1e-3; coarser than BEM's ~1e-5 = volume +")
print("     PML-truncation error).  With act7_23 (ngsolve.bem), BOTH NGSolve high-freq tools -- FEM+PML")
print("     and Helmholtz BEM -- are on the DtN yardstick, not proxies.")

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "act7_24_ngsolve_fem_pml_highfreq_dtn.json")
with open(out, "w") as f:
    json.dump({"ka": KA, "a": A_RAD, "R_out": R_OUT, "ndof": fes.ndof, "modes": res}, f, indent=2)
print(f"  wrote {os.path.basename(out)}")

print("\n" + "=" * 80)
print(" ALL CHECKS PASSED" if N_FAIL == 0 else f" {N_FAIL} CHECK(S) FAILED")
print("=" * 80)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
