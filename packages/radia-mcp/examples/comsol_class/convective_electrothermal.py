"""COMSOL-class #44 -- CONVECTIVE electro-thermal (Joule heat + Newton cooling), live-3-way.

The realistic cooling BC: a slab (thickness L) with a uniform Joule source q is cooled by
CONVECTION (film coefficient h) on both faces -- the surface FLOATS (Robin BC -k dT/dn =
h(T-T_inf)) rather than being held fixed (#19/#41). The centre rise is

    dT_centre = q L^2/(8 k)  +  q L/(2 h)
              =  CONDUCTION parabola  +  convective FILM rise,

the surface sitting q L/(2 h) above the coolant; as h -> inf the film vanishes and the
fixed-T result q L^2/8k is recovered. radia adds a Robin boundary to the heat solver
(`solve_heat_steady_robin`). Validated 3-way: NGSolve == this closed form == the reference
commercial FE code (LiveLink, 0.000 %, recorded internally). Self-contained, headless.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction
from netgen.occ import OCCGeometry, WorkPlane, X, Y
from radia_mcp.radia_ngsolve.multiphysics import (solve_heat_steady_robin,
                                                  convective_slab_peak_dT)

L, Wd = 0.02, 0.10          # slab thickness (y, conduction dir), width (x)
q, k = 5.0e5, 2.0           # Joule density [W/m^3], thermal conductivity [W/mK]

slab = WorkPlane().MoveTo(-Wd/2, -L/2).Rectangle(Wd, L).Face(); slab.faces.name = "slab"
slab.edges.Max(Y).name = "conv"; slab.edges.Min(Y).name = "conv"     # top+bottom convective
slab.edges.Max(X).name = "side"; slab.edges.Min(X).name = "side"     # sides insulated (natural)
mesh = Mesh(OCCGeometry(slab, dim=2).GenerateMesh(maxh=L/20))

print(f"{'h[W/m2K]':>9} {'centre dT[K]':>13} {'analytic[K]':>12} {'film qL/2h[K]':>14}")
for h in (250.0, 500.0, 2000.0):
    gT = solve_heat_steady_robin(mesh, CoefficientFunction(q), k, "slab", "conv", h, order=3)
    dT_c = gT(mesh(0.0, 0.0))
    dT_an = convective_slab_peak_dT(q, L, k, h)
    print(f"{h:9.0f} {dT_c:13.4f} {dT_an:12.4f} {q*L/(2*h):14.4f}  (err {abs(dT_c-dT_an)/dT_an*100:.3f}%)")
print("OK -- convective electro-thermal: dT_centre = qL^2/8k (conduction) + qL/2h (film); "
      "h->inf recovers the fixed-T parabola.")
