"""compare_mesh_vs_peec_same_mesh.py

Direct A/B comparison of two coil representations on the IDENTICAL
mesh, no Kelvin, bounded sphere + Dirichlet A=0 on outer boundary:

  (A) Mesh coil: square-section torus meshed as a volume, uniform
      current density J0 = I / (2a * 2a) on "coil" material.
      Full A-formulation FEM:
        int nu curl(A) . curl(v) dV = int J . v dx("coil")
      Solve for A. Inductance L_A = 2 * int 0.5 nu |curl(A)|^2 dV / I^2.

  (B) Radia analytical: SAME uniform-J square-section torus, but A
      evaluated ANALYTICALLY via rad.ObjArcCur + rad.Fld('a').
      Project A_s onto HCurl, compute
      L_B = 2 * int 0.5 nu |curl(gfu_As)|^2 dV / I^2.

Both cases represent the IDENTICAL physical current distribution
(uniform J in a square-section torus), so with an identical mesh
and identical Dirichlet BC the FEM solve of (A) and the projection
of the analytical A of (B) must agree up to discretization.

Earlier versions used a PEEC filament bundle for (B); that is a
DIFFERENT current distribution (discrete lines vs continuous J)
and unnecessarily muddled the comparison. The filament path is
maintained separately for future PEEC work.
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

import radia  # noqa: F401

from ngsolve import (HCurl, BilinearForm, LinearForm, GridFunction,
                      CoefficientFunction, Integrate, InnerProduct,
                      Conj, curl, dx, x, y, z, sqrt, log,
                      TaskManager, Mesh)
from netgen.occ import (WorkPlane, Axes, Pnt, Dir, Sphere, Axis, Glue,
                         OCCGeometry)


MU_0 = 4e-7 * math.pi
NU_0 = 1.0 / MU_0


def biot_savart_A_cf(filament_paths, currents, xc, yc, zc):
    Ax = CoefficientFunction(0.0)
    Ay = CoefficientFunction(0.0)
    Az = CoefficientFunction(0.0)
    for path, Ik in zip(filament_paths, currents):
        for p1, p2 in path:
            p1 = [float(c) for c in p1]; p2 = [float(c) for c in p2]
            t = [p2[i] - p1[i] for i in range(3)]
            L = math.sqrt(sum(tt * tt for tt in t))
            if L < 1e-14:
                continue
            th = [tt / L for tt in t]
            dx1 = xc - p1[0]; dy1 = yc - p1[1]; dz1 = zc - p1[2]
            dx2 = xc - p2[0]; dy2 = yc - p2[1]; dz2 = zc - p2[2]
            r1 = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1 + 1e-24)
            r2 = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2 + 1e-24)
            a_ = dx1 * th[0] + dy1 * th[1] + dz1 * th[2]
            b_ = dx2 * th[0] + dy2 * th[1] + dz2 * th[2]
            coeff = (MU_0 * float(Ik) / (4 * math.pi)) * \
                log((r1 + a_ + 1e-24) / (r2 + b_ + 1e-24))
            Ax = Ax + coeff * th[0]
            Ay = Ay + coeff * th[1]
            Az = Az + coeff * th[2]
    return CoefficientFunction((Ax, Ay, Az))


def build_shared_mesh(R_coil, a_coil, gap_deg, R_sphere, maxh, maxh_coil):
    """Geometry: SQUARE-cross-section torus (material 'coil') in air sphere.

    Both the volume-J FEM case and the filament-bundle case use a
    2a_coil x 2a_coil square cross-section to eliminate the
    shape-mismatch source of discrepancy. The bundle natively samples
    a rectangular (2a, 2a) bounding box, so matching the torus shape
    to that is a like-for-like comparison.
    """
    wp = WorkPlane(Axes(p=Pnt(R_coil, 0, 0),
                         n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    # Square cross-section centered on the coil centerline:
    wp.MoveTo(-a_coil, -a_coil)
    square = wp.Rectangle(2 * a_coil, 2 * a_coil).Face()
    torus = square.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)),
                             360 - gap_deg)
    torus.name = "coil"
    torus.maxh = maxh_coil

    sphere = Sphere(Pnt(0, 0, 0), R_sphere)
    sphere.maxh = maxh
    for f in sphere.faces:
        f.name = "outer"
    air = sphere - torus
    air.name = "air"

    shape = Glue([air, torus])
    geo = OCCGeometry(shape)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(2)
    return mesh


def main():
    R_coil = 0.030
    a_coil = 0.003
    gap_deg = 5.0
    R_sphere = 0.30
    maxh = 0.015
    maxh_coil = 0.003
    order = 1
    I_total = 1.0
    nw = nh = 3
    n_arc = 30

    L_analytical = MU_0 * R_coil * (math.log(8 * R_coil / a_coil) - 2) \
        * (1.0 - gap_deg / 360.0)

    print("=" * 60)
    print("Mesh-coil (volume J) vs PEEC (bundle filament A_s)")
    print("SAME mesh, no Kelvin, Dirichlet A=0 on outer sphere")
    print("=" * 60)
    print(f"Coil       : R={R_coil*1e3:.0f} mm, a={a_coil*1e3:.1f} mm, "
          f"gap={gap_deg} deg")
    print(f"Air sphere : R={R_sphere*1e3:.0f} mm, maxh={maxh*1e3:.1f} mm")
    print(f"Bundle     : {nw} x {nh} fils, n_arc={n_arc}")
    print(f"L_analytical = {L_analytical*1e9:.2f} nH  (open-domain, torus)")
    print()

    print("[1/5] Building shared mesh...")
    t0 = time.perf_counter()
    mesh = build_shared_mesh(R_coil, a_coil, gap_deg,
                              R_sphere, maxh, maxh_coil)
    print(f"  {mesh.ne} tets, {mesh.nv} vertices, "
          f"mats={mesh.GetMaterials()}, bnds={mesh.GetBoundaries()} "
          f"({time.perf_counter() - t0:.1f}s)")

    # Reluctivity: nu_0 everywhere (no Kelvin)
    nu_cf = CoefficientFunction(NU_0)

    # Shared FE space: HCurl order=1 with Dirichlet on outer sphere
    fes = HCurl(mesh, order=order, dirichlet="outer")
    u, v = fes.TnT()
    print(f"  HCurl ndof = {fes.ndof}")

    # ---------------------------------------------------------------
    # Case A: mesh coil + volume J
    # ---------------------------------------------------------------
    print("\n[2/5] Case A: volume J in coil mesh")
    t0 = time.perf_counter()
    # Square 2a x 2a cross-section (matches the mesh and Radia ObjArcCur)
    J0 = I_total / (2 * a_coil * 2 * a_coil)
    r_xy = sqrt(x * x + y * y)
    J_dir = CoefficientFunction((-y / (r_xy + 1e-12),
                                   x / (r_xy + 1e-12), 0))
    J_source = J0 * J_dir

    a_bf = BilinearForm(fes)
    a_bf += nu_cf * curl(u) * curl(v) * dx
    a_bf += 1e-8 * NU_0 * u * v * dx
    f_A = LinearForm(fes)
    f_A += J_source * v * dx("coil")
    with TaskManager():
        a_bf.Assemble()
        f_A.Assemble()

    gfu_A = GridFunction(fes)
    with TaskManager():
        gfu_A.vec.data = a_bf.mat.Inverse(fes.FreeDofs(),
                                           inverse="pardiso") * f_A.vec

    W_A = Integrate(0.5 * nu_cf
                     * InnerProduct(curl(gfu_A), curl(gfu_A)),
                     mesh, order=8).real
    L_A = 2.0 * W_A / I_total ** 2
    # Mutual-energy definition: L = (1/I^2) * int J . A dV over coil.
    # Independent of BC because J is localized.
    JA_A = Integrate(InnerProduct(J_source, gfu_A),
                      mesh, definedon=mesh.Materials("coil"), order=8).real
    L_A_JA = JA_A / I_total ** 2
    print(f"  W_A = {W_A:.6e} J -> L_A(energy) = {L_A*1e9:.3f} nH")
    print(f"  J.A_coil = {JA_A:.6e} -> L_A(J.A) = {L_A_JA*1e9:.3f} nH")
    print(f"  ({time.perf_counter() - t0:.1f}s)")

    # ---------------------------------------------------------------
    # Case B: same mesh + same uniform-J coil, but A from Radia ObjArcCur
    # ---------------------------------------------------------------
    print("\n[3/5] Case B: Radia ObjArcCur analytical A on SAME mesh")
    t0 = time.perf_counter()
    import radia as rad
    rad.UtiDelAll()
    radia_coil = rad.ObjArcCur(
        [0, 0, 0],                                              # center
        [R_coil - a_coil, R_coil + a_coil],                     # [r_in, r_out]
        [0, np.deg2rad(360 - gap_deg)],                         # phi range
        2 * a_coil,                                             # h
        40,                                                     # nseg
        'man', 'z', J0)
    print(f"  Radia handle = {radia_coil}")

    # Evaluate A at all mesh vertices.
    vert_coords = np.array([[mesh.vertices[i].point[j]
                              for j in range(3)]
                             for i in range(mesh.nv)])
    print(f"  Sampling A at {mesh.nv} vertices...")
    A_vec_at_vertex = np.array([rad.Fld(radia_coil, 'a', list(pt))
                                 for pt in vert_coords])
    # Build a VoxelCoefficient or component CFs from vertex values.
    # Simplest route: wrap via an H1 GridFunction per component, then
    # combine and project to HCurl.
    from ngsolve import VectorH1
    fes_vec = VectorH1(mesh, order=1)
    gfu_A_vec = GridFunction(fes_vec)
    if fes_vec.ndof != mesh.nv * 3:
        raise RuntimeError(
            f"Expected VectorH1 order=1 ndof == 3*nv "
            f"({3*mesh.nv}), got {fes_vec.ndof}")
    # VectorH1 order=1 DOF layout (empirically verified): per-component
    # blocks, i.e. [Ax(v0), Ax(v1), ..., Ax(v_N-1),
    #               Ay(v0), Ay(v1), ..., Ay(v_N-1),
    #               Az(v0), Az(v1), ..., Az(v_N-1)].
    nv = mesh.nv
    arr = gfu_A_vec.vec.FV().NumPy()
    arr[0:nv] = A_vec_at_vertex[:, 0]
    arr[nv:2 * nv] = A_vec_at_vertex[:, 1]
    arr[2 * nv:3 * nv] = A_vec_at_vertex[:, 2]

    # Compute curl directly from VectorH1 gradients (piecewise-constant for
    # order=1 VH1; avoids the HCurl projection loss).
    from ngsolve import grad, CoefficientFunction as CF
    G = grad(gfu_A_vec)   # (3, 3) tensor CF: G[i, j] = d A_i / d x_j
    curl_A = CF((G[2, 1] - G[1, 2],
                 G[0, 2] - G[2, 0],
                 G[1, 0] - G[0, 1]))
    W_B = Integrate(0.5 * nu_cf * InnerProduct(curl_A, curl_A),
                     mesh, order=8).real
    L_B = 2.0 * W_B / I_total ** 2
    # Mutual-energy definition with the SAME J source.
    JA_B = Integrate(InnerProduct(J_source, gfu_A_vec),
                      mesh, definedon=mesh.Materials("coil"), order=8).real
    L_B_JA = JA_B / I_total ** 2
    print(f"  W_B = {W_B:.6e} J -> L_B(energy) = {L_B*1e9:.3f} nH")
    print(f"  J.A_coil = {JA_B:.6e} -> L_B(J.A) = {L_B_JA*1e9:.3f} nH")
    print(f"  ({time.perf_counter() - t0:.1f}s)")

    print()
    print("=" * 60)
    print(f"  {'':25s}  {'L(energy)':>10s}  {'L(J.A)':>10s}")
    print(f"  {'Case A (FEM volume J)':25s}  {L_A*1e9:10.3f}  "
          f"{L_A_JA*1e9:10.3f}")
    print(f"  {'Case B (Radia ObjArcCur)':25s}  {L_B*1e9:10.3f}  "
          f"{L_B_JA*1e9:10.3f}")
    print(f"  L_open (analytical circ) = {L_analytical*1e9:.3f} nH")
    print(f"  L(J.A) B/A ratio = {L_B_JA/L_A_JA:.4f}  "
          f"(this should be ~1 since J is the SAME and BC matters less)")
    print("=" * 60)


if __name__ == "__main__":
    main()
