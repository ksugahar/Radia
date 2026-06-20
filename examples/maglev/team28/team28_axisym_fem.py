"""team28_axisym_fem.py -- repo-clean port of the lab TEAM 28 axisymmetric
full-FEM levitation-force solve (ground-truth baseline for the CLN-reduced
verification).

Ported from the learning material
    W:\\00_CAE\\NGSolve\\01_菅原\\2024_08_TEAM28\\50Hz_可動\\axisymmetric_mixed.ipynb
(axisymmetric MIXED phi-B formulation, anisotropic-nu infinite-element
shell for the open boundary).  This reproduces the full-FEM levitation
force F_z at a single disk height so we have a runnable, self-contained
in-repo baseline before layering the CLN reduction on top.

Reference (from the lab .mat, Fz1 vs dZ): at dZ=0 (disk bottom at 10.8mm)
F_z = -2.1928 N  (sign convention: NEGATIVE F_z = upward lift).  Disk
weight ~1.055 N; stable levitation equilibrium ~dZ=+4mm.

Run:  python team28_axisym_fem.py
"""
from numpy import pi, linspace, array, real, imag
from netgen.occ import MoveTo, Glue, OCCGeometry
from ngsolve import (
    Mesh, H1, L2, VectorL2, GridFunction, BilinearForm,
    LinearForm, CoefficientFunction, Integrate, CF, x, dx, TaskManager,
    ngsglobals,
)

# ----- TEAM 28 geometry (meters), exactly as the lab material -----
domain_a = 0.120          # outer air radius
domain_z = -0.020         # vertical centre offset of the air domain

coil_1_w, coil_1_h, coil_1_r = 0.028, 0.052, 0.041
coil_2_w, coil_2_h, coil_2_r = 0.015, 0.052, 0.041 + 0.0465

aluminium_w = 0.065       # disk radial extent (R = 65 mm)
aluminium_h = 0.003       # disk thickness 3 mm
aluminium_z = 0.0108      # disk bottom height = 10.8 mm  (the 可動 dZ=0 position)

N1, I1 = 960, 20.0        # coil 1: 960 turns, +20 A
N2, I2 = 576, 20.0        # coil 2: 576 turns, -20 A (counter-wound)
SIGMA_AL = 3.4e7          # disk conductivity [S/m]
FREQ = 50.0
mu0 = 4 * pi * 1e-7


def build_mesh():
    semicircle = []
    abc = linspace(1, 1.1, 5)
    for n in range(len(abc)):
        circle = MoveTo(0, domain_z).Circle(abc[n] * domain_a).Face()
        rect = MoveTo(-abc[n] * domain_a, -abc[n] * domain_a + domain_z) \
            .Rectangle(abc[n] * domain_a, 2 * abc[n] * domain_a).Face()
        semicircle.append(circle - rect)
        if n == 0:
            semicircle[n].bc('Air')
        elif n == len(abc) - 1:
            semicircle[n].bc(f'u{n}')
            semicircle[n].edges.name = 'outer'
        else:
            semicircle[n].bc(f'u{n}')
    coil_1 = MoveTo(coil_1_r - coil_1_w / 2, -coil_1_h) \
        .Rectangle(coil_1_w, coil_1_h).Face(); coil_1.bc('Coil_1')
    coil_2 = MoveTo(coil_2_r - coil_2_w / 2, -coil_2_h) \
        .Rectangle(coil_2_w, coil_2_h).Face(); coil_2.bc('Coil_2')
    al = MoveTo(0, aluminium_z).Rectangle(aluminium_w, aluminium_h).Face()
    al.bc('Al')
    shape = Glue(semicircle + [coil_1, coil_2, al])
    return Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=0.005))


def solve_force():
    ngsglobals.msg_level = 0
    mesh = build_mesh()
    print("materials:", mesh.GetMaterials())

    p = 2
    fesPhi = H1(mesh, order=p, dirichlet="outer", complex=True)
    fesB = VectorL2(mesh, order=p - 1, complex=True)
    fes = fesPhi * fesB
    (u, B), (v, dB) = fes.TnT()

    # B from the stream function: B = (-1/r dphi/dz, 1/r dphi/dr)
    from ngsolve import grad

    def Bop(phi):
        g = grad(phi)
        return CF((-1 / x * g[1], 1 / x * g[0]))

    # anisotropic-nu infinite-element shell (open boundary)
    nu = mesh.MaterialCF(
        {"u1": 1 / (0.285096 * mu0), "u2": 1 / (12.4933 * mu0),
         "u3": 1 / (0.051874 * mu0), "u4": 1 / (147.652 * mu0)},
        default=1 / mu0)
    sig = mesh.MaterialCF({"Al": SIGMA_AL}, default=0.0)
    s = 2 * pi * FREQ * 1j

    a = BilinearForm(nu * x * (B * dB - Bop(u) * dB + Bop(v) * B) * dx)
    a += (v * (s * sig * (u / x))) * dx
    Jz = mesh.MaterialCF(
        {"Coil_1": I1 * N1 / (coil_1_w * coil_1_h),
         "Coil_2": -I2 * N2 / (coil_2_w * coil_2_h)}, default=0.0)
    f = LinearForm((v * Jz) * dx)

    with TaskManager():
        a.Assemble()
        f.Assemble()
        gfuB = GridFunction(fes)
        gfuB.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

        gfu = GridFunction(fesPhi); gfu.Set(gfuB.components[0])
        gfB1 = GridFunction(fesB);  gfB1.Set(gfuB.components[1])

        fesL2 = L2(mesh, order=p, dirichlet="outer", complex=True)
        gfJt = GridFunction(fesL2); gfJt.Set(-sig * s * gfu / x)

        # time-average vertical Lorentz force on the disk (lab formula)
        fz = Integrate(
            (real(gfB1[0]) * real(gfJt) - imag(gfB1[0]) * imag(gfJt))
            * (2 * pi * x) * dx(mesh.Materials("Al")), mesh)
    return float(real(fz))


if __name__ == "__main__":
    fz = solve_force()
    ref = -2.1928
    print(f"\nF_z (CLN-baseline full-FEM) = {fz:+.4f} N")
    print(f"reference (lab .mat, dZ=0)  = {ref:+.4f} N")
    print(f"rel. error = {abs(fz - ref) / abs(ref) * 100:.2f} %")
