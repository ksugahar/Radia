"""COMSOL-class #24 -- TWO-WAY temperature-dependent sigma(T) electro-thermal.

The genuinely 2-way electro-thermal coupling: electrical conductivity falls with temperature,
sigma(T) = sigma0/(1 + alpha T), so the Joule density q = J^2/sigma(T) = (J^2/sigma0)(1+alpha T)
RISES with T and feeds back into the heat -- thermal <-> electrical.  A Joule-heated bar
(length L, cooled ends T=0, sides insulated, uniform current density J) has an EXACT 1-D
closed form (the resulting ODE -k T'' = (J^2/sigma0)(1+alpha T) is linear):

    T(x)  = (1/alpha)[ cos(s x) + ((1-cos(s L))/sin(s L)) sin(s x) - 1 ],   s = sqrt(alpha J^2/(sigma0 k))
    T_max = (1/alpha)(sec(s L/2) - 1)            (-> J^2 L^2/(8 sigma0 k), the parabola, as alpha->0)

radia-ngsolve solves it with `multiphysics.solve_heat_steady_nonlinear` (Picard on the
temperature-dependent source).  Self-contained / headless / tool-independent (analytic ref).
"""
import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction
from netgen.occ import OCCGeometry, MoveTo, X
from radia_mcp.radia_ngsolve.multiphysics import solve_heat_steady_nonlinear

L, hbar = 0.10, 0.01
sigma0, k, alpha, J = 1.0e6, 50.0, 0.005, 8.0e5
s = math.sqrt(alpha * J * J / (sigma0 * k))


def T_closed(x):
    return (1.0/alpha)*(math.cos(s*x) + ((1-math.cos(s*L))/math.sin(s*L))*math.sin(s*x) - 1.0)
Tmax_cf = (1.0/alpha)*(1.0/math.cos(s*L/2.0) - 1.0)
Tmax_parabola = J*J*L*L/(8.0*sigma0*k)            # constant-sigma (no feedback)

face = MoveTo(0, 0).Rectangle(L, hbar).Face(); face.faces.name = "bar"
face.edges.Min(X).name = "ends"; face.edges.Max(X).name = "ends"
mesh = Mesh(OCCGeometry(face, dim=2).GenerateMesh(maxh=L/60))

# q(T) = J^2 / sigma(T) = (J^2/sigma0)(1 + alpha T)  -- the 2-way feedback source
T = solve_heat_steady_nonlinear(mesh, lambda Tcf: (J*J/sigma0)*(1.0 + alpha*Tcf),
                                k, conductor="bar", dirichlet="ends", order=2)
Tmax = T(mesh(L/2, hbar/2))
err = abs(Tmax - Tmax_cf)/Tmax_cf*100
prof = max(abs(T(mesh(x, hbar/2)) - T_closed(x)) for x in
           [0.1*L, 0.25*L, 0.5*L, 0.75*L, 0.9*L])/Tmax_cf*100

print(f"sqrt-arg s*L = {s*L:.3f}")
print(f"T_max: FE={Tmax:.4f} K   closed-form={Tmax_cf:.4f} K   err={err:.3f} %")
print(f"  constant-sigma parabola (no feedback) = {Tmax_parabola:.4f} K  ->  2-way feedback "
      f"raises the peak {100*(Tmax_cf-Tmax_parabola)/Tmax_parabola:.2f} %")
print(f"  profile max err = {prof:.3f} %")
print("OK" if err < 0.5 and prof < 0.5 else f"OFF (Tmax {err:.2f}% / profile {prof:.2f}%)")
