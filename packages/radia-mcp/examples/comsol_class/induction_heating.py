"""COMSOL-class multiphysics in NGSolve #1: INDUCTION HEATING (EM -> thermal).

A non-magnetic conducting cylinder in a uniform axial AC field B0. Harmonic
eddy currents (A-Phi) dissipate Joule heat, which drives a steady temperature
field (surface held cool). One-way magneto-thermal coupling -- the canonical
COMSOL induction-heating tutorial, reproduced with the radia_ngsolve helpers.

Validated (this geometry, ~60k tets):
    P/L (central section) vs exact cylinder eddy loss (I0/I1 Bessel) : 8.9 %
    T_centre              vs 1-D radial heat equation                : 9.6 %

Run: python induction_heating.py   (needs ngsolve + radia_mcp on PYTHONPATH)
"""
import numpy as np
from scipy.integrate import quad
from ngsolve import (Mesh, CoefficientFunction, Integrate, IfPos, grad, dx,
                     TaskManager, x, y, z)
from netgen.csg import CSGeometry, Cylinder, OrthoBrick, Plane, Pnt, Vec
from radia_mcp.radia_ngsolve.solve import solve_eddy_current_harmonic_APhi
from radia_mcp.radia_ngsolve.multiphysics import joule_loss_density, solve_heat_steady
from radia_mcp.radia_ngsolve.force import MU0

# -- parameters --
mu0 = MU0; nu0 = 1.0 / mu0
sigma = 1.6e7; f = 500.0; omega = 2 * np.pi * f
B0 = 0.10; R = 0.020; Lh = 0.050; k_th = 50.0; Lb = 0.12
kappa = np.sqrt(1j * omega * sigma * mu0); H0 = B0 / mu0

# -- analytic reference (infinite cylinder, complex power-series Bessel) --
def I0c(zz, n=80):
    s = 1 + 0j; t = 1 + 0j; w = (zz / 2)**2
    for kk in range(1, n): t *= w / kk**2; s += t
    return s
def I1c(zz, n=80):
    s = zz / 2; t = zz / 2; w = (zz / 2)**2
    for kk in range(1, n): t *= w / (kk * (kk + 1)); s += t
    return s
def qabs(r): return 0.5 * abs(-H0 * kappa * I1c(kappa * r) / I0c(kappa * R))**2 / sigma
PL_ref = quad(lambda r: qabs(r) * 2 * np.pi * r, 0, R)[0]
Tc_ref = (1 / k_th) * quad(lambda s: quad(lambda t: qabs(t) * t, 0, s)[0] / s, 1e-6, R)[0]

# -- geometry: finite cylinder (surface named) + air box --
cyl = ((Cylinder(Pnt(0, 0, -1), Pnt(0, 0, 1), R)
        * Plane(Pnt(0, 0, -Lh), Vec(0, 0, -1)) * Plane(Pnt(0, 0, Lh), Vec(0, 0, 1)))
       .mat("cyl").bc("cyl_surf").maxh(0.0025))
box = OrthoBrick(Pnt(-Lb, -Lb, -Lb), Pnt(Lb, Lb, Lb)).bc("outer")
geo = CSGeometry(); geo.Add(cyl); geo.Add((box - cyl).mat("air"))
mesh = Mesh(geo.GenerateMesh(maxh=0.05, grading=0.3)); mesh.Curve(2)

# -- EM: scattered harmonic eddy A-Phi, axial background B0 --
nu_cf = CoefficientFunction([nu0] * len(mesh.GetMaterials()))
sigma_cf = mesh.MaterialCF({"cyl": sigma}, default=0.0)
A0 = CoefficientFunction((-B0 * y / 2, B0 * x / 2, 0.0))   # curl A0 = B0 z
def add_source(lf, test_fn):
    vA, psi = test_fn
    lf += -1j * omega * sigma_cf * A0 * vA * dx
    lf += -1j * omega * sigma_cf * A0 * grad(psi) * dx
with TaskManager():
    gfAPhi = solve_eddy_current_harmonic_APhi(
        mesh, nu_cf, sigma_cf, omega, add_source, order=2, precond="local",
        dirichlet="outer", periodic=False, tol=1e-8, maxsteps=2000, restart=200, reg=1e-10)
gfA, gfPhi = gfAPhi.components

# -- couple: Joule density (NOTE A0 for the scattered solve) -> steady heat --
q = joule_loss_density(gfA, gfPhi, sigma_cf, omega, A0=A0)
gT = solve_heat_steady(mesh, q, k_th, conductor="cyl", dirichlet="cyl_surf")

# -- report --
zc = 0.020
PL_mid = Integrate(q * IfPos(zc * zc - z * z, 1.0, 0.0) * dx, mesh).real / (2 * zc)
T_centre = gT(mesh(1e-4, 1e-4, 0.0))
print(f"P/L (central)  = {PL_mid:8.1f} W/m  | exact {PL_ref:8.1f}  ({abs(PL_mid-PL_ref)/PL_ref*100:.1f} %)")
print(f"T_centre       = {T_centre:8.2f} K    | 1-D   {Tc_ref:8.2f}  ({abs(T_centre-Tc_ref)/Tc_ref*100:.1f} %)")
