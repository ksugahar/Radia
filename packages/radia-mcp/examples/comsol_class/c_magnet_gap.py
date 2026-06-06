"""COMSOL-class #18 -- iron-yoke window-frame DIPOLE: air-gap field (magnetic circuit).

The accelerator/actuator workhorse: an iron frame (mu_r) closes the flux from a coil and
forces it across a small air gap, where the field does the work. The reluctance model
(Ampere's law around the loop, B continuous) gives

    B_gap = mu0 N I / (gap + iron_path / mu_r)   (magnetic_circuit_gap_field).

Method: 2D planar A_z with high-mu iron (nu = NU0/mu_r in the yoke) and a coil threading
the left leg (+NI inside the window, -NI outside the frame). The gap field is read as the
AVERAGE of B_x over the gap region. The FEM value sits a few % below the reluctance model
-- fringing/leakage widen the effective gap -- which is exactly the physics the lumped
model misses and a COMSOL cross-check would pin down.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
import math
from ngsolve import Mesh, CoefficientFunction, grad, Integrate, dx, TaskManager
from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import solve_planar_magnetostatic, magnetic_circuit_gap_field, NU0, MU0

W, ww, b = 0.24, 0.12, 0.06               # outer side, window side, bar width [m]
g, pw = 0.006, 0.06                        # gap length, pole width (= bar height) [m]
mu_r, NI = 2000.0, 5000.0                  # iron permeability, amp-turns
Rbox = 0.60
l_fe = 2.0 * ((W + ww) / 2.0) * 2 - g      # mean iron path ~ centreline perimeter - gap

def rect(x0, y0, w, h):
    return MoveTo(x0, y0).Rectangle(w, h).Face()

big    = rect(-W/2, -W/2, W, W)
window = rect(-ww/2, -ww/2, ww, ww)
gap    = rect(-g/2, ww/2, g, b)                       # slot in top bar (x in +-g/2, y in [ww/2, W/2])
coil_in  = rect(-ww/2 + 0.005, -0.04, 0.02, 0.08)     # in window, against left leg (+NI)
coil_out = rect(-W/2 - 0.025, -0.04, 0.02, 0.08)      # outside frame, against left leg (-NI)
iron = big - window - gap
gap.faces.name = "gap"; coil_in.faces.name = "cin"; coil_out.faces.name = "cout"
iron.faces.name = "iron"
box = rect(-Rbox, -Rbox, 2*Rbox, 2*Rbox)
air = box - iron - gap - coil_in - coil_out
air.faces.name = "air"
for f in (air,):
    f.edges.Max(X).name = "outer"; f.edges.Min(X).name = "outer"
    f.edges.Max(Y).name = "outer"; f.edges.Min(Y).name = "outer"
geo = OCCGeometry(Glue([air, iron, gap, coil_in, coil_out]), dim=2)
mesh = Mesh(geo.GenerateMesh(maxh=0.04))
print(f"mesh: {mesh.ne} elems, materials={mesh.GetMaterials()}")

# area-normalize the coil so int Jz = +/- NI exactly
Ain  = Integrate(mesh.MaterialCF({"cin": 1.0}, default=0.0) * dx, mesh)
Aout = Integrate(mesh.MaterialCF({"cout": 1.0}, default=0.0) * dx, mesh)
Jz = mesh.MaterialCF({"cin": NI / Ain, "cout": -NI / Aout}, default=0.0)
nu = mesh.MaterialCF({"iron": NU0 / mu_r}, default=NU0)
with TaskManager():
    gfu = solve_planar_magnetostatic(mesh, nu, Jz=Jz, order=3, dirichlet="outer")
B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))      # (Bx, By)

Agap = Integrate(mesh.MaterialCF({"gap": 1.0}, default=0.0) * dx, mesh)
Bx_gap = Integrate(mesh.MaterialCF({"gap": 1.0}, default=0.0) * B[0] * dx, mesh) / Agap
B_ref = magnetic_circuit_gap_field(1.0, NI, g, l_fe, mu_r)  # N*I passed as NI
B_ideal = MU0 * NI / g                                       # mu_r -> inf
err = abs(abs(Bx_gap) - B_ref) / B_ref * 100
print(f"\niron path l_fe ~ {l_fe:.3f} m ; gap g = {g*1e3:.1f} mm ; mu_r = {mu_r:.0f}")
print(f"B_gap FEM      = {abs(Bx_gap):.4f} T  (avg |B_x| over the gap)")
print(f"B_gap reluct.  = {B_ref:.4f} T  = mu0 NI/(g + l_fe/mu_r)")
print(f"B_gap ideal    = {B_ideal:.4f} T  (mu_r -> inf)")
print(f"FEM vs reluctance model: {err:.2f}%  (residual = fringing/leakage)")
print("OK" if err < 8.0 else f"OFF ({err:.2f}%)")
