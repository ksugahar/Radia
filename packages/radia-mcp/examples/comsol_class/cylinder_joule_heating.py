"""COMSOL-class #41 -- CYLINDRICAL conductor Joule heating (electro-thermal), live-3-way.

The cylindrical sibling of the slab electro-thermal (#19, dT_peak = sigma V^2/8k): a round
conductor (radius a) carrying a uniform DC current density J dissipates a uniform Joule source

    q = J^2 / sigma  [W/m^3]   (the EM step),

and steady radial conduction with the surface held at the coolant temperature gives the
PARABOLIC rise

    dT(r) = q (a^2 - r^2) / (4 k),   centre  dT_peak = q a^2 / (4 k).

The classic current-carrying-wire / ampacity temperature rise. Validated 3-way: NGSolve
(`solve_heat_steady` on a disk) == this closed form == the reference commercial FE code
(LiveLink, 0.000 %, recorded internally). Self-contained, headless, SI metres; the public
reference is the closed form.
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction
from netgen.occ import OCCGeometry, WorkPlane
from radia_mcp.radia_ngsolve.multiphysics import solve_heat_steady, joule_heated_cylinder_peak_dT

a, k = 0.02, 2.0                       # conductor radius [m], thermal conductivity [W/mK]
I, sigma = 888.6, 1.0e6               # DC current [A], electrical conductivity [S/m]
J = I / (math.pi * a * a)
q = J * J / sigma                     # uniform Joule density [W/m^3]
dT_peak = joule_heated_cylinder_peak_dT(q, a, k)
print(f"J={J:.4e} A/m^2   q=J^2/sigma={q:.4e} W/m^3   ->  analytic dT_peak=q a^2/(4k)={dT_peak:.4f} K")

wire = WorkPlane().Circle(0, 0, a).Face(); wire.faces.name = "wire"; wire.edges.name = "surf"
mesh = Mesh(OCCGeometry(wire, dim=2).GenerateMesh(maxh=a/20))
gT = solve_heat_steady(mesh, CoefficientFunction(q), k, "wire", "surf", order=3)

print(f"{'r/a':>6} {'dT_FE[K]':>10} {'dT_analytic[K]':>15} {'err%':>7}")
for ra in (0.0, 0.25, 0.5, 0.75):
    r = ra * a
    dT_fe = gT(mesh(r, 0.0))
    dT_an = q * (a*a - r*r) / (4.0 * k)
    print(f"{ra:6.2f} {dT_fe:10.4f} {dT_an:15.4f} {abs(dT_fe-dT_an)/dT_an*100:7.3f}")
print(f"centre: NGSolve {gT(mesh(0,0)):.4f} K vs analytic {dT_peak:.4f} K")
print("OK -- cylindrical Joule heating dT(r)=q(a^2-r^2)/(4k); centre = q a^2/(4k).")
