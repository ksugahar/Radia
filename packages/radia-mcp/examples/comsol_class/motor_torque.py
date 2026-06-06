"""COMSOL-class #2: ROTATING-MACHINE TORQUE (air-gap Maxwell stress).
Core motor physics: a hard PM 'rotor' (Br along x) in a uniform applied field B0
at angle theta. Torque-angle characteristic  tau_z(theta) = m B0 sin(theta),
m = (Br/mu0) * V.  Hard magnet (mu_r=1) -> fields SUPERPOSE, so ONE solve gives
the whole angle sweep (B_total = B_PM + B0(theta)). Torque via eggshell_torque."""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
import numpy as np
from ngsolve import Mesh, CoefficientFunction, curl
from netgen.csg import CSGeometry, Cylinder, OrthoBrick, Plane, Pnt, Vec
from radia_mcp.radia_ngsolve.solve import solve_magnetostatic_Aform, remanent_source
from radia_mcp.radia_ngsolve.force import eggshell_torque, MU0

NU0 = 1.0/MU0
Br = 1.0; B0 = 0.5; R = 0.020; Lh = 0.040           # cylinder axis z, half-length 40mm
V = np.pi*R**2*(2*Lh)
m = (Br/MU0)*V
print(f"magnetic moment m = {m:.3f} A m^2,  tau_max = m B0 = {m*B0:.3f} N m")

t0 = time.time()
cyl = ((Cylinder(Pnt(0,0,-1),Pnt(0,0,1),R)
        * Plane(Pnt(0,0,-Lh),Vec(0,0,-1)) * Plane(Pnt(0,0,Lh),Vec(0,0,1)))
       .mat("magnet").maxh(0.006))
box = OrthoBrick(Pnt(-0.15,-0.15,-0.15), Pnt(0.15,0.15,0.15)).bc("outer")
geo = CSGeometry(); geo.Add(cyl); geo.Add((box-cyl).mat("air"))
mesh = Mesh(geo.GenerateMesh(maxh=0.05, grading=0.3)); mesh.Curve(2)
print(f"mesh: {mesh.ne} elems ({time.time()-t0:.1f}s)")

# one PM solve (Br along x); hard magnet -> B_total = B_PM + B0 superposes
curl_src = remanent_source(mesh, {"magnet": (Br, 0.0, 0.0)})
B_PM = curl(solve_magnetostatic_Aform(mesh, NU0, curl_source=curl_src))
print(f"solve {time.time()-t0:.1f}s")

print("theta   tau_z_fem    tau_z_ref(m B0 sin)   err")
worst = 0.0
for deg in (15, 30, 45, 60, 90):
    th = np.radians(deg)
    B0v = CoefficientFunction((B0*np.cos(th), B0*np.sin(th), 0.0))
    Tx, Ty, Tz = eggshell_torque(B_PM + B0v, mesh, center=(0,0,0),
                                 r_inner=0.055, r_outer=0.085)
    ref = m*B0*np.sin(th)
    err = abs(Tz-ref)/abs(ref)*100
    worst = max(worst, err)
    print(f"{deg:3d}    {Tz:8.3f}     {ref:8.3f}            {err:5.1f}%")
print("OK" if worst < 8 else f"OFF (worst {worst:.1f}%)")
