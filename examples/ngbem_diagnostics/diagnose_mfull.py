"""Quick diagnostic: check M_full and V_full eigenvalues."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import Mesh, HDivSurface, BilinearForm, InnerProduct, ds, BND
from ngsolve.bem import MaxwellSingleLayerPotentialOperator
from ngsolve import TaskManager
from ngbem_eddy import LoopBasisBuilder

MU_0 = 4e-7 * np.pi

block = Box(Pnt(-0.05, -0.05, -0.001), Pnt(0.05, 0.05, 0.001))
block.solids.name = "conductor"
block.faces.name = "surface"
geo = OCCGeometry(block)
mesh = Mesh(geo.GenerateMesh(maxh=0.025))

fes = HDivSurface(mesh, order=0)
n_dof = fes.ndof
u, v = fes.TnT()

# Build loop basis
loop = LoopBasisBuilder(mesh, fes)
active = loop.active_dofs
n_a = loop.n_active
print(f"DOFs: {n_dof}, active: {n_a}, loops: {loop.n_loops}")

# Mass matrix
M_bf = BilinearForm(fes)
M_bf += InnerProduct(u.Trace(), v.Trace()) * ds(bonus_intorder=2)
with TaskManager():
    M_bf.Assemble()
    M_ngmat = M_bf.mat

    # Extract M using M's own vectors
    x_m = M_ngmat.CreateColVector()
    y_m = M_ngmat.CreateColVector()
    M_full = np.zeros((n_a, n_a))
    for i, di in enumerate(active):
        x_m[:] = 0
        x_m[di] = 1.0
        M_ngmat.Mult(x_m, y_m)
        for j, dj in enumerate(active):
            M_full[j, i] = y_m[dj]

    print(f"\nM_full shape: {M_full.shape}")
    print(f"M_full symmetry: ||M-M^T|| = {np.linalg.norm(M_full - M_full.T):.2e}")
    M_eig = np.linalg.eigvalsh(M_full)
    print(f"M_full eigenvalues: min={M_eig.min():.6e}, max={M_eig.max():.6e}")
    print(f"M_full positive definite: {M_eig.min() > 0}")
    print(f"M_full trace: {np.trace(M_full):.6f}")
    print(f"M_full diagonal range: [{M_full.diagonal().min():.6f}, {M_full.diagonal().max():.6f}]")

    # Project to loop space
    M_LL = loop.T_loop.T @ M_full @ loop.T_loop
    M_LL_eig = np.linalg.eigvalsh(M_LL)
    print(f"\nM_LL eigenvalues: min={M_LL_eig.min():.6e}, max={M_LL_eig.max():.6e}")
    print(f"M_LL positive definite: {M_LL_eig.min() > 0}")

    # Also extract M using V's vectors (the original way)
    V_op = MaxwellSingleLayerPotentialOperator(fes, kappa=1.0, intorder=4)
    V_ngmat = V_op.mat
    x_v = V_ngmat.CreateColVector()
    y_m2 = M_ngmat.CreateColVector()
    M_full2 = np.zeros((n_a, n_a))
    for i, di in enumerate(active):
        x_v[:] = 0
        x_v[di] = 1.0
        M_ngmat.Mult(x_v, y_m2)
        for j, dj in enumerate(active):
            M_full2[j, i] = y_m2[dj]

    M_eig2 = np.linalg.eigvalsh(M_full2)
    diff = np.linalg.norm(M_full - M_full2)
    print(f"\nM_full (V vectors) eigenvalues: min={M_eig2.min():.6e}, max={M_eig2.max():.6e}")
    print(f"||M_proper - M_from_V_vectors|| = {diff:.2e}")

    # V_full
    V_full = np.zeros((n_a, n_a), dtype=complex)
    x_v2 = V_ngmat.CreateColVector()
    y_v = V_ngmat.CreateColVector()
    for i, di in enumerate(active):
        x_v2[:] = 0
        x_v2[di] = 1.0
        V_ngmat.Mult(x_v2, y_v)
        for j, dj in enumerate(active):
            V_full[j, i] = y_v[dj]

    print(f"\nV_full symmetry: ||V-V^T|| = {np.linalg.norm(V_full - V_full.T):.2e}")
    V_real_eig = np.linalg.eigvalsh(0.5*(V_full.real + V_full.real.T))
    print(f"V_full.real eigenvalues: min={V_real_eig.min():.6e}, max={V_real_eig.max():.6e}")

    V_LL = loop.T_loop.T @ V_full @ loop.T_loop
    V_LL_real_sym = 0.5 * (V_LL.real + V_LL.real.T)
    V_LL_eig = np.linalg.eigvalsh(V_LL_real_sym)
    print(f"\nV_LL.real eigenvalues: min={V_LL_eig.min():.6e}, max={V_LL_eig.max():.6e}")
    print(f"V_LL positive definite: {V_LL_eig.min() > 0}")
