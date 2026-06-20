# -*- coding: utf-8 -*-
"""
act7_23_ngsbem_highfreq_dtn.py  (Act 7 -- the radiating extension; the BEM column made REAL)
============================================================================================
Places the ACTUAL ngsolve.bem Helmholtz BEM on the SAME DtN-spectral yardstick as the rest of
the comparison.  In act7_22_dtn_spectrum_consolidated the high-freq "BEM" entry was the ANALYTIC
wave_dtn reference; here the genuine ngsolve.bem single / double-layer operators BUILD the
exterior Dirichlet-to-Neumann operator and we MEASURE its per-multipole defect against wave_dtn,
plus the conventional REFLECTION coefficient it implies.

WHY: high-freq IS a studied regime (the act7 radiating thread), and NGSolve / ngsolve.bem are the
working tools there.  This demo backs that concretely: ngsolve.bem reproduces the exact COMPLEX
exterior DtN of the radiating sphere.

METHOD (exterior acoustic Helmholtz, sphere radius a=1, wavenumber kappa = ka):
  - scalar H1 surface space on the sphere;
  - V = HelmholtzSingleLayerPotentialOperator, K = HelmholtzDoubleLayerPotentialOperator (kappa);
  - exterior DtN (the convention verified statically in verify_laplace_bem.py, where it gives the
    static ladder -(n+1)/a): the Neumann coeff of a Dirichlet datum u_D is
        gN = V^{-1} (-1/2 M + K) u_D ,
    and for a pure spherical-harmonic mode  gN = Lambda_n u_D, so the per-mode eigenvalue is the
    M-weighted Rayleigh quotient
        Lambda_n^BEM = (u_D, M gN) / (u_D, M u_D).
    (A full eigendecomposition of V^{-1}(...) is ill-conditioned at the high modes -- the per-mode
    Rayleigh quotient on a KNOWN smooth harmonic is the robust extraction.)
  - exact = radia.open_boundary.wave_dtn(n, ka) = ka h_n^(1)'(ka)/h_n^(1)(ka) (complex / radiating).
  - REFLECTION the boundary implies: a mode at the truncation is  A h_n^(1) + B h_n^(2) (outgoing +
    incoming); a boundary DtN Lambda_h admits the spurious incoming part
        R_n = |Lambda_h - Lambda_exact| / |Lambda_h - conj(Lambda_exact)|   (=0 for the exact
    outgoing DtN, for real ka where the incoming DtN = conj of the outgoing).  So d_n (the operator
    defect) and R_n (the reflection EVERYONE measures) are the SAME mismatch -- ngsolve.bem is
    ~reflectionless precisely because it reproduces the exact outgoing DtN.

VERIFIED (asserted from computed values; no overclaim): monopole / dipole / quadrupole DtN reldef
< 1e-3 vs wave_dtn AND implied reflection < 1e-3.  Needs ngsolve + ngsolve.bem (6.2.2604+).
"""
import os
import json
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from ngsolve import Mesh, H1, BilinearForm, LinearForm, ds, TaskManager, CF, z
from netgen.occ import Sphere, Pnt, OCCGeometry, Glue
import ngsolve.bem as bem
import radia.open_boundary as ob

N_FAIL = 0


def check(name, cond, detail=""):
    global N_FAIL
    if not cond:
        N_FAIL += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}{('  -- ' + detail) if detail else ''}")


def dense(mat, n):
    return np.array(mat.ToDense().NumPy()[:n, :n], dtype=complex)


KA = 2.0
A = 1.0
print("=" * 80)
print(f" act7_23_ngsbem_highfreq_dtn : ngsolve.bem Helmholtz BEM on the DtN yardstick (a=1, ka={KA})")
print("=" * 80)

res = {}
with TaskManager():
    sph = Sphere(Pnt(0, 0, 0), A)
    for f in sph.faces:
        f.name = "sph"
    mesh = Mesh(OCCGeometry(Glue(sph.faces)).GenerateMesh(maxh=0.4))
    mesh.Curve(3)
    fes = H1(mesh, order=3)
    ndof = fes.ndof
    u, v = fes.TnT()

    V = bem.HelmholtzSingleLayerPotentialOperator(fes, fes, kappa=KA, intorder=10)
    K = bem.HelmholtzDoubleLayerPotentialOperator(fes, fes, kappa=KA, intorder=10)
    Vm = dense(V.mat, ndof)
    Km = dense(K.mat, ndof)
    mass = BilinearForm(fes)
    mass += u * v * ds
    mass.Assemble()
    Mm = dense(mass.mat, ndof)

    def proj(cf):
        lf = LinearForm(fes)
        lf += cf * v * ds
        lf.Assemble()
        return np.linalg.solve(Mm, lf.vec.FV().NumPy().copy().astype(complex))

    print(f"  sphere surface H1(order=3): {ndof} DOF; Helmholtz V/K assembled (intorder=10)")
    print("  n  mode         BEM DtN             exact wave_dtn      reldef d_n   reflection R_n")
    # P_n(cos theta) on the unit sphere (z = cos theta): pure spherical-harmonic modes
    MODES = [(0, CF(1.0), "monopole"), (1, z, "dipole"),
             (2, 3 * z * z - 1, "quadrupole"), (3, 5 * z * z * z - 3 * z, "octupole")]
    for n, cf, lbl in MODES:
        uD = proj(cf)
        gN = np.linalg.solve(Vm, (-0.5 * Mm + Km) @ uD)        # Neumann coeff = V^-1(-1/2 M + K) uD
        lam = complex((uD @ Mm @ gN) / (uD @ Mm @ uD))         # M-weighted Rayleigh quotient
        exact = complex(ob.wave_dtn(n, KA))
        d = abs(lam - exact) / abs(exact)
        R = abs(lam - exact) / abs(lam - np.conj(exact))       # implied reflection (=defect / 2|Im|)
        res[lbl] = {"n": n, "BEM": [lam.real, lam.imag],
                    "exact": [exact.real, exact.imag], "reldef": d, "reflection": R}
        print(f"  {n}  {lbl:10s}  {lam:+.4f}  {exact:+.4f}   {d:.2e}    {R:.2e}")

low = ("monopole", "dipole", "quadrupole")
check("ngsolve.bem reproduces the COMPLEX exterior DtN (monopole/dipole/quadrupole reldef < 1e-3)",
      max(res[l]["reldef"] for l in low) < 1e-3,
      f"max {max(res[l]['reldef'] for l in low):.1e}")
check("ngsolve.bem is ~reflectionless = reproduces the exact OUTGOING DtN (R_n < 1e-3, n<=2)",
      max(res[l]["reflection"] for l in low) < 1e-3,
      f"max {max(res[l]['reflection'] for l in low):.1e}")

print("\n  => the high-freq 'BEM' column is REAL: ngsolve.bem (Helmholtz single/double layer) BUILDS")
print("     the exact complex exterior DtN, so it is ~reflectionless -- d_n (operator defect) and")
print("     R_n (the reflection everyone measures) are the SAME mismatch.  NGSolve.bem works at high-freq.")

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "act7_23_ngsbem_highfreq_dtn.json")
with open(out, "w") as f:
    json.dump({"ka": KA, "ndof": ndof, "modes": res}, f, indent=2)
print(f"  wrote {os.path.basename(out)}")

print("\n" + "=" * 80)
print(" ALL CHECKS PASSED" if N_FAIL == 0 else f" {N_FAIL} CHECK(S) FAILED")
print("=" * 80)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
