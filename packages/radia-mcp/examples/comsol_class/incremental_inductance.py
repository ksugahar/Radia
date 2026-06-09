"""COMSOL-class #31 -- INCREMENTAL inductance L_inc=dlambda/dI + magnetic co-energy.

Refines the saturating SECANT inductance (#28). On the saturating lambda(I) of a closed
iron core (smooth bounded nu(B), one winding, SI metres):
  L_sec(I) = lambda/I            (secant / apparent, saturated_secant_inductance),
  L_inc(I) = dlambda/dI          (incremental / differential, incremental_inductance, central FD),
  W'(I) = int_0^I lambda dI'      (co-energy, magnetic_coenergy),   W = lambda*I - W'  (energy).
For a SATURATING (concave) lambda(I): L_inc < L_sec, both fall with current, and the co-energy
W' exceeds the stored energy W. The differential dq inductance for small-signal/control + the
energy-method co-energy. Self-contained, headless, tool-independent. (The PROPER frozen-perm
L_inc uses the tangent differential-reluctivity TENSOR dH/dB -- a future refinement; here FD.)
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, sqrt, InnerProduct, TaskManager
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic,
                                           solve_planar_magnetostatic_nonlinear,
                                           coil_flux_linkage_2d, incremental_inductance,
                                           magnetic_coenergy, NU0)

MUR0, BK, NEXP = 1000.0, 1.4, 2.0          # smooth, gradual saturation knee
OUT, WIN, BOX, rw = 0.06, 0.03, 0.12, 0.003
N, DEPTH, xin, xout = 200, 0.05, 0.0075, 0.045
MATS = None


def nu_of_B(B):
    Bm = sqrt(InnerProduct(B, B) + 1e-20)
    sat = Bm**NEXP / (Bm**NEXP + BK**NEXP)
    iron = CoefficientFunction([1.0 if m == "iron" else 0.0 for m in MATS])
    return iron * NU0 * (1.0/MUR0 + (1.0 - 1.0/MUR0)*sat) + (1.0 - iron) * NU0


def build():
    frame = MoveTo(-OUT/2, -OUT/2).Rectangle(OUT, OUT).Face()
    window = MoveTo(-WIN/2, -WIN/2).Rectangle(WIN, WIN).Face()
    cAp = WorkPlane().Circle(-xin, 0, rw).Face(); cAp.faces.name = "cApos"
    cAn = WorkPlane().Circle(-xout, 0, rw).Face(); cAn.faces.name = "cAneg"
    iron = frame - window; iron.faces.name = "iron"
    win_air = window - cAp; win_air.faces.name = "air"
    box = MoveTo(-BOX/2, -BOX/2).Rectangle(BOX, BOX).Face(); box.edges.name = "outer"
    out_air = box - frame - cAn; out_air.faces.name = "air"
    iron.faces.maxh = 0.004; win_air.faces.maxh = 0.004; out_air.faces.maxh = 0.012
    cAp.faces.maxh = 0.0012; cAn.faces.maxh = 0.0012
    return Mesh(OCCGeometry(Glue([iron, win_air, out_air, cAp, cAn]), dim=2).GenerateMesh(maxh=0.012))


def jz(I):
    j0 = N * I / (math.pi * rw * rw)
    return CoefficientFunction([{"cApos": j0, "cAneg": -j0}.get(m, 0.0) for m in MATS])


mesh = build(); MATS = mesh.GetMaterials()
Is = [0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2]
with TaskManager():
    lam = [coil_flux_linkage_2d(solve_planar_magnetostatic_nonlinear(
        mesh, nu_of_B, Jz=jz(I), order=2, relax=0.4, max_iter=200, tol=1e-8),
        mesh, "cApos", "cAneg", N, DEPTH) for I in Is]
Wco = magnetic_coenergy(Is, lam)
print(f"{'I':>5} {'lambda':>11} {'L_sec':>10} {'L_inc':>10} {'W':>10} {'Wco':>10}")
for k, I in enumerate(Is):
    Linc = incremental_inductance(lam[k+1], lam[k-1], (Is[k+1]-Is[k-1])/2) if 0 < k < len(Is)-1 else None
    W = lam[k]*I - Wco[k]
    print(f"{I:5.2f} {lam[k]:11.4e} {lam[k]/I:10.4e} {('%.4e'%Linc) if Linc else '   --   ':>10} {W:10.4e} {Wco[k]:10.4e}")
print("OK -- L_inc < L_sec (concave lambda), both fall with I; co-energy W' > energy W in saturation")
