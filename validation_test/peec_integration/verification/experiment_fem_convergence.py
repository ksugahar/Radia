"""Quick FEM convergence study for L_air.

Usage:
  python test_fem_convergence.py [order] [maxh_mm] [coil_maxh_mm]
  python test_fem_convergence.py 1 5     -> order=1, maxh=5mm (baseline)
  python test_fem_convergence.py 2 8     -> order=2, maxh=8mm
"""
import os, sys, time
import numpy as np

os.environ['NETGENDIR'] = r'S:\NGSolve\01_GitHub\install_ngsolve\bin'
os.environ['PATH'] = (r'S:\NGSolve\01_GitHub\install_ngsolve\bin;'
                       + os.environ.get('PATH', ''))
sys.path.insert(0, r'S:\NGSolve\01_GitHub\install_ngsolve\Lib\site-packages')

from ngsolve import *
from ngsolve import TaskManager
from netgen.occ import *

MU_0 = 4e-7 * np.pi
R_COIL = 20e-3
W_WIRE = 1e-3
H_WIRE = 1e-3
I_COIL = 1.0
AIR_R = 120e-3

# Analytical references
r_equiv = np.sqrt(W_WIRE * H_WIRE / np.pi)
L_no_int = MU_0 * R_COIL * (np.log(8 * R_COIL / r_equiv) - 2)
L_with_int = MU_0 * R_COIL * (np.log(8 * R_COIL / r_equiv) - 7/4)

print(f"Analytical (external only):  {L_no_int * 1e9:.2f} nH")
print(f"Analytical (with internal):  {L_with_int * 1e9:.2f} nH")
print()

# Parameters from command line
order = int(sys.argv[1]) if len(sys.argv) > 1 else 1
maxh = float(sys.argv[2]) * 1e-3 if len(sys.argv) > 2 else 5e-3
coil_maxh = float(sys.argv[3]) * 1e-3 if len(sys.argv) > 3 else W_WIRE

print(f"Parameters: order={order}, maxh={maxh*1e3:.1f}mm, "
      f"coil_maxh={coil_maxh*1e3:.2f}mm, AIR_R={AIR_R*1e3:.0f}mm")
print()

# Build geometry (same as verify_ngsolve_inductance.py)
t0 = time.time()
wp = WorkPlane(Axes((R_COIL, 0, 0), n=Y, h=Z))
rect = wp.Rectangle(W_WIRE, H_WIRE).Face()
rect = rect.Move((-W_WIRE / 2, -H_WIRE / 2, 0))
coil = Revolve(rect, Axis((0, 0, 0), Z), 360)
coil.mat('coil')
coil.faces.maxh = coil_maxh

# Name sphere faces BEFORE boolean ops (critical!)
air_sphere = Sphere(Pnt(0, 0, 0), AIR_R)
air_sphere.faces.name = 'outer'

# Subtract coil from air
air = air_sphere - coil
air.mat('air')

# Glue: air first, then objects
shape = Glue([air, coil])
geo = OCCGeometry(shape)
ngmesh = geo.GenerateMesh(maxh=maxh)
mesh = Mesh(ngmesh)
with TaskManager():
    mesh.Curve(2)

    ne = mesh.ne
    nv = mesh.nv
    materials = mesh.GetMaterials()
    t_mesh = time.time() - t0
    print(f"Mesh: {ne} elements, {nv} vertices ({t_mesh:.1f}s)")
    print(f"Materials: {materials}")
    print()

    # Solve (same formulation as verify_ngsolve_inductance.py)
    t0 = time.time()
    sigma_reg = 1.0
    freq_reg = 0.01
    omega_reg = 2 * np.pi * freq_reg

    # Per-material inv_mu
    inv_mu_vals = [1.0 / MU_0] * len(materials)
    inv_mu = CoefficientFunction(inv_mu_vals)

    fes = HCurl(mesh, order=order, dirichlet='outer', complex=True, nograds=True)
    A_trial, v_test = fes.TnT()

    r_xy = sqrt(x * x + y * y + 1e-30)
    J_mag = I_COIL / (W_WIRE * H_WIRE)
    J_cf = CoefficientFunction((-J_mag * y / r_xy, J_mag * x / r_xy, 0))

    a = BilinearForm(fes)
    a += inv_mu * InnerProduct(curl(A_trial), curl(v_test)) * dx
    a += 1j * omega_reg * sigma_reg * InnerProduct(A_trial, v_test) * dx
    a += 1e-12 * InnerProduct(A_trial, v_test) * dx
    a.Assemble()

    f = LinearForm(fes)
    f += InnerProduct(J_cf, v_test) * dx('coil')
    f.Assemble()

    gfA = GridFunction(fes)
    ndof = fes.ndof
    print(f"DOFs: {ndof}")
    print("Solving...", flush=True)
    gfA.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * f.vec

    # Energy -> inductance
    B = curl(gfA)
    B_conj = Conj(B)
    B_dot_Bconj = B[0] * B_conj[0] + B[1] * B_conj[1] + B[2] * B_conj[2]
    W_stored = 0.25 * Integrate(inv_mu * B_dot_Bconj, mesh).real
    L = 4 * W_stored / I_COIL**2

    t_solve = time.time() - t0
    diff_no_int = (L - L_no_int) / L_no_int * 100
    diff_with_int = (L - L_with_int) / L_with_int * 100

    print(f"W_stored = {W_stored:.6e} J")
    print(f"L = {L * 1e9:.2f} nH ({t_solve:.1f}s)")
    print(f"  vs analytical (no internal):   {diff_no_int:+.1f}%")
    print(f"  vs analytical (with internal): {diff_with_int:+.1f}%")
