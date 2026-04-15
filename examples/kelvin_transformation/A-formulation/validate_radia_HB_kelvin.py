"""validate_radia_HB_kelvin.py

Validation of "external source + Kelvin" handling using Radia's
analytical volume-J coil as the source. Both legs are evaluated on
the IDENTICAL Kelvin-transformed mesh (Sugahara 2022 two-sphere)
and the inductance is extracted from the H . B energy integral:

    W = 0.5 * int H . B dV = 0.5 * int nu(r) * |B(r)|^2 dV
    L = 2 W / I^2

  Case A (baseline): full FEM with volume J source on a meshed
                     square 2a x 2a torus (uniform J0 = I/(4a^2)).
                     Solve curl(nu curl A) = J. Energy from curl(A_FEM).

  Case B (external):   no FEM solve. Inject Radia's open-domain A_s
                       (rad.ObjArcCur, same square torus, same J0)
                       onto the same mesh. In the Kelvin exterior
                       domain, evaluate A_s at the Kelvin-mapped
                       physical position and apply the Phase 2
                       1-form pullback. Energy from curl(A_s_comp).

Both cases use IDENTICAL physical current (same uniform-J square
torus) and IDENTICAL mesh + Kelvin convention. They must give the
same L if the Kelvin handling of an external (analytic) source is
correct.
"""

from __future__ import annotations

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src'))
SRC_RADIA = os.path.join(SRC, 'radia')
for p in (SRC, SRC_RADIA):
    if p not in sys.path:
        sys.path.insert(0, p)

import radia as rad   # noqa

from ngsolve import (HCurl, VectorH1, BilinearForm, LinearForm,
                      GridFunction, CoefficientFunction as CF,
                      Integrate, InnerProduct, Conj, curl, grad, dx,
                      x, y, z, sqrt, Periodic, TaskManager, Mesh)
from netgen.occ import (WorkPlane, Axes, Pnt, Dir, Sphere, Axis, Vertex,
                         Glue, OCCGeometry, IdentificationType)


MU_0 = 4e-7 * math.pi
NU_0 = 1.0 / MU_0

R_coil = 0.030
a_coil = 0.003
gap_deg = 5.0
arc_deg = 360.0 - gap_deg
R_K = 0.06
offset = 0.15
maxh = 0.012
maxh_coil = 0.003
order = 1
I_total = 1.0


def build_kelvin_mesh():
    wp = WorkPlane(Axes(p=Pnt(R_coil, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    wp.MoveTo(-a_coil, -a_coil)
    sq = wp.Rectangle(2 * a_coil, 2 * a_coil).Face()
    torus = sq.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), arc_deg)
    torus.name = "coil"; torus.maxh = maxh_coil

    inner = Sphere(Pnt(0, 0, 0), R_K); inner.maxh = maxh
    for f in inner.faces:
        f.name = "kelvin_int"
    inner_air = inner - torus; inner_air.name = "air"

    outer = Sphere(Pnt(offset, 0, 0), R_K); outer.maxh = maxh * 2
    for f in outer.faces:
        f.name = "kelvin_ext"
    outer.name = "kelvin"

    gnd = Vertex(Pnt(offset, 0, 0)); gnd.name = "GND"

    for f in inner_air.faces:
        if f.name == "kelvin_int":
            ki = f; break
    for f in outer.faces:
        if f.name == "kelvin_ext":
            ke = f; break
    ki.Identify(ke, "periodic", IdentificationType.PERIODIC)

    shape = Glue([inner_air, torus, outer, gnd])
    geo = OCCGeometry(shape)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh, grading=0.5)
    mesh = Mesh(ngmesh); mesh.Curve(order + 1)
    return mesh


def kelvin_factor_cf():
    r_prime_sq = (x - offset) ** 2 + y ** 2 + z ** 2 + 1e-24
    return R_K * R_K / r_prime_sq


def nu_cf_kelvin(mesh):
    fac = kelvin_factor_cf()
    nu_dict = {m: (NU_0 * fac if "kelvin" in m.lower() else NU_0)
               for m in mesh.GetMaterials()}
    return mesh.MaterialCF(nu_dict, default=NU_0)


def case_A_fem(mesh):
    """Baseline full-A FEM with volume J on the meshed square torus."""
    nu_cf = nu_cf_kelvin(mesh)
    J0 = I_total / (2 * a_coil) ** 2
    r_xy = sqrt(x * x + y * y)
    r_xy_s = r_xy + 1e-12
    J_src = J0 * CF((-y / r_xy_s, x / r_xy_s, 0))

    fes_base = HCurl(mesh, order=order, dirichlet_bbnd="GND")
    fes = Periodic(fes_base)
    u, v = fes.TnT()

    a_bf = BilinearForm(fes)
    a_bf += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
    a_bf += 1e-8 * NU_0 * u * v * dx
    f_lf = LinearForm(fes)
    f_lf += J_src * v * dx("coil")
    with TaskManager():
        a_bf.Assemble(); f_lf.Assemble()
    gfu = GridFunction(fes)
    with TaskManager():
        gfu.vec.data = a_bf.mat.Inverse(fes.FreeDofs(),
                                         inverse="pardiso") * f_lf.vec
    # H . B = ν * |B|^2  with B = curl(A)
    W_inner = Integrate(
        0.5 * nu_cf * InnerProduct(curl(gfu), Conj(curl(gfu))),
        mesh, definedon=mesh.Materials("air|coil"), order=10).real
    W_kext = Integrate(
        0.5 * nu_cf * InnerProduct(curl(gfu), Conj(curl(gfu))),
        mesh, definedon=mesh.Materials("kelvin"), order=10).real
    L_inner = 2.0 * W_inner / I_total ** 2
    L_kext = 2.0 * W_kext / I_total ** 2
    print(f"  L_inner = {L_inner*1e9:.3f} nH, "
          f"L_kext = {L_kext*1e9:.3f} nH")
    return L_inner + L_kext


def build_radia_coil():
    rad.UtiDelAll()
    J0 = I_total / (2 * a_coil) ** 2
    return rad.ObjArcCur([0, 0, 0],
                          [R_coil - a_coil, R_coil + a_coil],
                          [0, math.radians(arc_deg)],
                          2 * a_coil, 200, 'man', 'z', J0)


def case_B_radia_external(mesh, phase=2):
    """No FEM solve: inject Radia A_s on mesh, compute energy directly.

    In inner ('air' / 'coil'): A_comp(r) = A_radia(r) (Cartesian).
    In Kelvin exterior ('kelvin'): A_comp(r') = pullback of
        A_radia(r_phys = Kelvin_inv(r')).
    """
    coil = build_radia_coil()
    nu_cf = nu_cf_kelvin(mesh)

    # Per-vertex A in COMPUTATIONAL coords:
    nv = mesh.nv
    A_arr = np.zeros((nv, 3))
    # Pre-classify vertices by material via element loop
    from ngsolve import BND, VOL
    vert_in_kelvin = np.zeros(nv, dtype=bool)
    for el in mesh.Elements(VOL):
        if "kelvin" in mesh.GetMaterials()[el.index].lower():
            for v in el.vertices:
                vert_in_kelvin[v.nr] = True

    # GND is the Kelvin-image of physical infinity, where A = 0. We
    # cannot evaluate Kelvin-mapped A there (1/rho' singularity).
    rho_min = R_K / 100.0   # skip vertices closer than this to GND
    skipped = 0

    # Compute A at each vertex
    for vi in range(nv):
        p = mesh.vertices[vi].point
        if vert_in_kelvin[vi]:
            # Kelvin map: r_phys = c + R^2 (r' - c) / |r' - c|^2
            d = np.array([p[0] - offset, p[1], p[2]])
            rho = math.sqrt(np.dot(d, d))
            if rho < rho_min:
                A_arr[vi] = 0.0
                skipped += 1
                continue
            rho_sq = rho * rho
            r_phys = np.array([offset, 0.0, 0.0]) + (R_K * R_K / rho_sq) * d
            A_phys = np.array(rad.Fld(coil, 'a', list(r_phys)))
            # 1-form pullback under Kelvin r -> r' = R^2 r / |r|^2:
            #   A_comp = (rho'/R)^2 * (I - 2 n n^T) * A_phys
            # Note: (rho'/R)^2 (NOT (R/rho')^2; the latter is for scalar
            # potential phi or for B field, not for A 1-form).
            rho_sq = rho * rho   # = |r' - c|^2
            factor = rho_sq / (R_K * R_K)
            if phase == 1:
                A_comp = factor * A_phys
            elif phase == 2:
                n = d / rho
                A_dot_n = float(np.dot(A_phys, n))
                A_refl = A_phys - 2.0 * A_dot_n * n
                A_comp = factor * A_refl
            else:
                A_comp = A_phys
        else:
            A_comp = np.array(rad.Fld(coil, 'a', [p[0], p[1], p[2]]))
        A_arr[vi] = A_comp

    print(f"  skipped {skipped} vertices near GND (rho < R_K/100)")
    print(f"  |A_comp|: min={np.min(np.linalg.norm(A_arr, axis=1)):.3e}  "
          f"max={np.max(np.linalg.norm(A_arr, axis=1)):.3e}  "
          f"mean={np.mean(np.linalg.norm(A_arr, axis=1)):.3e}")
    fes_vec = VectorH1(mesh, order=1)
    gfu_A = GridFunction(fes_vec)
    arr = gfu_A.vec.FV().NumPy()
    arr[:nv] = A_arr[:, 0]
    arr[nv:2*nv] = A_arr[:, 1]
    arr[2*nv:3*nv] = A_arr[:, 2]

    # B = curl(A) via grad(VH1)
    G = grad(gfu_A)
    B = CF((G[2, 1] - G[1, 2],
            G[0, 2] - G[2, 0],
            G[1, 0] - G[0, 1]))
    W_inner = Integrate(0.5 * nu_cf * InnerProduct(B, B),
                         mesh, definedon=mesh.Materials("air|coil"),
                         order=10).real
    W_kext = Integrate(0.5 * nu_cf * InnerProduct(B, B),
                        mesh, definedon=mesh.Materials("kelvin"),
                        order=10).real
    L_inner = 2.0 * W_inner / I_total ** 2
    L_kext = 2.0 * W_kext / I_total ** 2
    print(f"  L_inner = {L_inner*1e9:.3f} nH, "
          f"L_kext = {L_kext*1e9:.3f} nH")
    return L_inner + L_kext


def main():
    print("=" * 60)
    print("Radia A_s + Kelvin FEM (energy = HB integral) vs FEM volume J")
    print("=" * 60)
    print(f"  R_coil={R_coil*1e3:.0f}mm a_coil={a_coil*1e3:.1f}mm "
          f"gap={gap_deg}deg  R_K={R_K*1e3:.0f}mm offset={offset*1e3:.0f}mm")
    print(f"  maxh={maxh*1e3}mm maxh_coil={maxh_coil*1e3}mm order={order}")
    print()

    print("[1/3] Building Kelvin mesh...")
    t0 = time.perf_counter()
    mesh = build_kelvin_mesh()
    print(f"  {mesh.ne} tets, {mesh.nv} verts, mats={mesh.GetMaterials()} "
          f"({time.perf_counter()-t0:.1f}s)")

    print("\n[2/3] Case A: full FEM (volume J meshed coil)")
    t0 = time.perf_counter()
    L_A = case_A_fem(mesh)
    print(f"  L_A = {L_A*1e9:.4f} nH  ({time.perf_counter()-t0:.1f}s)")

    for phase in [1, 2]:
        print(f"\n[3/3] Case B: Radia A_s injected, Phase {phase} pullback")
        t0 = time.perf_counter()
        L_B = case_B_radia_external(mesh, phase=phase)
        print(f"  L_B(phase={phase}) = {L_B*1e9:.4f} nH  "
              f"({time.perf_counter()-t0:.1f}s)")
        print(f"  diff vs Case A = {(L_B - L_A)/L_A*100:+.2f}%")

    print()
    print("=" * 60)


if __name__ == "__main__":
    main()
