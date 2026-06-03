"""team28_cln_force.py -- CLN-reduced TEAM 28 levitation force (coil-driven).

Goal: show that a Cauer-Ladder-Network (CLN) reduced model of the coil-
driven eddy-current disk reproduces the full-FEM Lorentz levitation force.

The axisymmetric mixed phi-B system is, at angular frequency s = j*omega,
    (K + s*N) X = F
with
    K  = the s-INDEPENDENT magnetostatic mixed operator
         nu*r*(B*dB - Bop(u)*dB + Bop(v)*B),
    N  = the conductivity term  v * (sigma * u/r),
    F  = the coil source  v * Jz.
The coil field is the excitation (per the chosen coupling), exactly as in
the full-FEM baseline team28_axisym_fem.py (which reproduces the lab
ground truth to 0.01%).

CLN / Cauer reduction = the Krylov subspace generated from the source by
the magnetostatic-solve / sigma-accumulate recursion
    V_0 = K^{-1} F                       (s=0 coil response, no eddy)
    V_{k+1} = orthonormalize( K^{-1} (N V_k) )
After m stages the reduced model (V^T K V + s V^T N V) y = V^T F gives the
reduced field X_m(s) = V y; we evaluate the Lorentz force from X_m(j*omega)
and watch it converge to the full-FEM force as m grows.

Run:  python team28_cln_force.py
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numpy import pi, linspace
from netgen.occ import MoveTo, Glue, OCCGeometry
from ngsolve import (
    Mesh, H1, L2, VectorL2, GridFunction, BilinearForm, LinearForm,
    Integrate, CF, x, dx, grad, TaskManager, ngsglobals,
)

# ----- TEAM 28 geometry (meters) -----
domain_a, domain_z = 0.120, -0.020
coil_1_w, coil_1_h, coil_1_r = 0.028, 0.052, 0.041
coil_2_w, coil_2_h, coil_2_r = 0.015, 0.052, 0.041 + 0.0465
aluminium_w, aluminium_h, aluminium_z = 0.065, 0.003, 0.0108
N1, I1, N2, I2 = 960, 20.0, 576, 20.0
SIGMA_AL, FREQ, mu0 = 3.4e7, 50.0, 4 * pi * 1e-7


def build_mesh(al_z=aluminium_z):
    semicircle = []
    abc = linspace(1, 1.1, 5)
    for n in range(len(abc)):
        circle = MoveTo(0, domain_z).Circle(abc[n] * domain_a).Face()
        rect = MoveTo(-abc[n] * domain_a, -abc[n] * domain_a + domain_z) \
            .Rectangle(abc[n] * domain_a, 2 * abc[n] * domain_a).Face()
        semicircle.append(circle - rect)
        semicircle[n].bc('Air' if n == 0 else f'u{n}')
        if n == len(abc) - 1:
            semicircle[n].edges.name = 'outer'
    coil_1 = MoveTo(coil_1_r - coil_1_w / 2, -coil_1_h).Rectangle(coil_1_w, coil_1_h).Face(); coil_1.bc('Coil_1')
    coil_2 = MoveTo(coil_2_r - coil_2_w / 2, -coil_2_h).Rectangle(coil_2_w, coil_2_h).Face(); coil_2.bc('Coil_2')
    al = MoveTo(0, al_z).Rectangle(aluminium_w, aluminium_h).Face(); al.bc('Al')
    shape = Glue(semicircle + [coil_1, coil_2, al])
    return Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=0.005))


def Bop(phi):
    g = grad(phi)
    return CF((-1 / x * g[1], 1 / x * g[0]))


def setup(al_z=aluminium_z, p=2):
    ngsglobals.msg_level = 0
    mesh = build_mesh(al_z)
    fesPhi = H1(mesh, order=p, dirichlet="outer", complex=True)
    fesB = VectorL2(mesh, order=p - 1, complex=True)
    fes = fesPhi * fesB
    (u, B), (v, dB) = fes.TnT()
    nu = mesh.MaterialCF(
        {"u1": 1 / (0.285096 * mu0), "u2": 1 / (12.4933 * mu0),
         "u3": 1 / (0.051874 * mu0), "u4": 1 / (147.652 * mu0)}, default=1 / mu0)
    sig = mesh.MaterialCF({"Al": SIGMA_AL}, default=0.0)
    Jz = mesh.MaterialCF(
        {"Coil_1": I1 * N1 / (coil_1_w * coil_1_h),
         "Coil_2": -I2 * N2 / (coil_2_w * coil_2_h)}, default=0.0)

    K = BilinearForm(nu * x * (B * dB - Bop(u) * dB + Bop(v) * B) * dx)
    Nmat = BilinearForm((v * (sig * u / x)) * dx)   # system = K + s*Nmat
    F = LinearForm((v * Jz) * dx)
    with TaskManager():
        K.Assemble(); Nmat.Assemble(); F.Assemble()
    return mesh, fes, fesPhi, fesB, sig, K, Nmat, F


def to_csr(mat, ndof):
    r, c, v = mat.COO()
    return sp.csr_matrix((np.array(v, dtype=complex),
                          (np.array(r), np.array(c))), shape=(ndof, ndof))


def force_from_vec(xfull, fes, fesPhi, fesB, sig, p=2):
    """Lorentz force on the disk from a full-DOF solution vector xfull."""
    mesh = fes.mesh
    gfuB = GridFunction(fes)
    gfuB.vec.FV().NumPy()[:] = xfull
    gfu = GridFunction(fesPhi); gfu.Set(gfuB.components[0])
    gfB1 = GridFunction(fesB);  gfB1.Set(gfuB.components[1])
    s = 2 * pi * FREQ * 1j
    fesL2 = L2(mesh, order=p, dirichlet="outer", complex=True)
    gfJt = GridFunction(fesL2); gfJt.Set(-sig * s * gfu / x)
    fz = Integrate(
        (gfB1[0].real * gfJt.real - gfB1[0].imag * gfJt.imag)
        * (2 * pi * x) * dx(mesh.Materials("Al")), mesh)
    return float(fz.real)


def cln_forces(al_z=aluminium_z, max_stage=10):
    """Return (fz_full, [fz_stage1, ..., fz_stageM]) for the disk at height al_z.

    fz_full = direct (K+sN) solve force; the list is the N-stage CLN/Cauer
    reduced force for N=1..M (M <= max_stage; the Krylov basis stops early
    once it stops growing).
    """
    mesh, fes, fesPhi, fesB, sig, K, Nmat, F = setup(al_z)
    ndof = fes.ndof
    s = 2 * pi * FREQ * 1j
    free = np.array([i for i in range(ndof) if fes.FreeDofs()[i]])

    Kf = to_csr(K.mat, ndof)[free][:, free].tocsc()
    Nf = to_csr(Nmat.mat, ndof)[free][:, free].tocsc()
    Ff = np.array(F.vec.FV().NumPy(), dtype=complex)[free]

    Klu = spla.splu(Kf)                       # K s-independent -> factor once
    xfull = np.zeros(ndof, dtype=complex)
    xfull[free] = spla.spsolve((Kf + s * Nf).tocsc(), Ff)
    fz_full = force_from_vec(xfull, fes, fesPhi, fesB, sig)

    # CLN / Cauer Krylov basis from the coil source
    V = []
    v0 = Klu.solve(Ff); v0 /= np.sqrt(abs(v0 @ v0.conj())); V.append(v0)
    for _ in range(max_stage - 1):
        w = Klu.solve(Nf @ V[-1])
        for vv in V:                          # modified Gram-Schmidt
            w = w - (vv.conj() @ w) * vv
        nrm = np.sqrt(abs(w @ w.conj()))
        if nrm < 1e-12:
            break
        V.append(w / nrm)

    stage_forces = []
    for m in range(1, len(V) + 1):
        Vm = np.array(V[:m]).T
        y = np.linalg.solve(Vm.conj().T @ (Kf @ Vm) + s * (Vm.conj().T @ (Nf @ Vm)),
                            Vm.conj().T @ Ff)
        xr = np.zeros(ndof, dtype=complex); xr[free] = Vm @ y
        stage_forces.append(force_from_vec(xr, fes, fesPhi, fesB, sig))
    return fz_full, stage_forces


def main():
    fz_full, sf = cln_forces()
    print(f"\n direct full-FEM (split K+sN) force = {fz_full:+.4f} N  (ref -2.1928)")
    print(f"\n stages  CLN force [N]   rel.err vs full")
    for m, fz in enumerate(sf, 1):
        print(f"   {m:2d}    {fz:+.5f}     {abs(fz - fz_full)/abs(fz_full)*100:8.3f} %")


if __name__ == "__main__":
    main()
