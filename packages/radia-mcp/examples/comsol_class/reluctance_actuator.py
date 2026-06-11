"""COMSOL-class #27 -- RELUCTANCE ACTUATOR holding force (electromagnet / relay / solenoid).

The force companion to the iron-yoke gap FIELD (#18): the same window-frame magnetic circuit
(coil NI, iron mu_r, one air gap g) now read for the POLE-FACE HOLDING FORCE via Maxwell stress

    F = B_gap^2 * A_face / (2 mu0),   A_face = pole-face area = b * depth,

with the reluctance model B_gap = mu0 NI/(g + l_fe/mu_r) giving the force law

    F = mu0 (NI)^2 A_face / [2 (g + l_fe/mu_r)^2]   ->   F*(g+l_fe/mu_r)^2 = const  (~ 1/g^2).

Method: current-driven 2D A_z (SI metres -- not scale-invariant), high-mu iron, coil threading
the left leg; the FE holding force (B_gap_FE^2 A_face/2mu0) is matched to the reluctance model
and the 1/g^2 scaling. Self-contained, headless, tool-independent (analytic reluctance ref).
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, grad, Integrate, dx, TaskManager
from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic,
                                           magnetic_circuit_gap_field,
                                           magnetic_circuit_gap_force, NU0, MU0)

W, ww, b = 0.24, 0.12, 0.06                # outer side, window side, pole-face height [m]
mu_r, NI, Rbox, DEPTH = 2000.0, 5000.0, 0.60, 1.0
A_FACE = b * DEPTH


def _rect(x0, y0, w, h):
    return MoveTo(x0, y0).Rectangle(w, h).Face()


def solve_gap(g):
    l_fe = 2.0 * ((W + ww) / 2.0) * 2 - g
    big = _rect(-W/2, -W/2, W, W); window = _rect(-ww/2, -ww/2, ww, ww)
    gap = _rect(-g/2, ww/2, g, b)
    cin = _rect(-ww/2 + 0.005, -0.04, 0.02, 0.08); cout = _rect(-W/2 - 0.025, -0.04, 0.02, 0.08)
    iron = big - window - gap
    gap.faces.name = "gap"; cin.faces.name = "cin"; cout.faces.name = "cout"; iron.faces.name = "iron"
    box = _rect(-Rbox, -Rbox, 2*Rbox, 2*Rbox)
    air = box - iron - gap - cin - cout; air.faces.name = "air"
    for f in (air,):
        f.edges.Max(X).name = "outer"; f.edges.Min(X).name = "outer"
        f.edges.Max(Y).name = "outer"; f.edges.Min(Y).name = "outer"
    mesh = Mesh(OCCGeometry(Glue([air, iron, gap, cin, cout]), dim=2).GenerateMesh(maxh=0.04))
    Ain = Integrate(mesh.MaterialCF({"cin": 1.0}, default=0.0)*dx, mesh)
    Aout = Integrate(mesh.MaterialCF({"cout": 1.0}, default=0.0)*dx, mesh)
    Jz = mesh.MaterialCF({"cin": NI/Ain, "cout": -NI/Aout}, default=0.0)
    nu = mesh.MaterialCF({"iron": NU0/mu_r}, default=NU0)
    with TaskManager():
        gfu = solve_planar_magnetostatic(mesh, nu, Jz=Jz, order=3, dirichlet="outer")
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    Agap = Integrate(mesh.MaterialCF({"gap": 1.0}, default=0.0)*dx, mesh)
    Bx = abs(Integrate(mesh.MaterialCF({"gap": 1.0}, default=0.0)*B[0]*dx, mesh) / Agap)
    F_fe = Bx*Bx*A_FACE/(2*MU0)
    F_model = magnetic_circuit_gap_force(1.0, NI, g, l_fe, mu_r, A_FACE)
    return Bx, F_fe, F_model, g + l_fe/mu_r


print(f"{'g[mm]':>6} {'B_gap[T]':>9} {'F_FE[N]':>10} {'F_model[N]':>11} {'F/Fm':>6} {'F*geff^2':>9}")
rows = []
for g in (0.004, 0.006, 0.009):
    Bx, F_fe, F_model, geff = solve_gap(g); rows.append((F_fe, geff))
    print(f"{g*1e3:6.1f} {Bx:9.4f} {F_fe:10.1f} {F_model:11.1f} {F_fe/F_model:6.3f} {F_fe*geff**2:9.4f}")
vals = [F*ge**2 for F, ge in rows]
spread = (max(vals)-min(vals))/(sum(vals)/len(vals))*100
print(f"reluctance const mu0(NI)^2 A/2 = {MU0*NI**2*A_FACE/2:.4f}; F*geff^2 spread = {spread:.2f}%")
print("OK" if spread < 6 else f"OFF ({spread:.2f}%)")
