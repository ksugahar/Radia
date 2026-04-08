"""
Check H1 interpolation accuracy for Biot-Savart A_inc and H_inc.
Compare point-evaluated Biot-Savart with H1 GridFunction at sample points.
"""

import math
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

MU_0 = 4e-7 * np.pi


def main():
    from ngsolve import (Mesh, H1, GridFunction, CF, Integrate, BND, TaskManager,
                         specialcf)
    from netgen.occ import Sphere, Cylinder, Pnt, Dir, OCCGeometry, Glue
    from biot_savart import h_segments_batch, a_segments_batch
    import time

    R_coil = 0.030
    R_wp = 0.010
    H_wp = 0.020
    R_air = 0.10
    I_total = 1.0

    # Generate coil segments
    gap_deg = 5
    gap_rad = math.radians(gap_deg)
    theta_max = 2 * math.pi - gap_rad
    n_seg = 200
    segments = []
    for i in range(n_seg):
        t1 = theta_max * i / n_seg
        t2 = theta_max * (i + 1) / n_seg
        segments.append(((R_coil * math.cos(t1), R_coil * math.sin(t1), 0),
                         (R_coil * math.cos(t2), R_coil * math.sin(t2), 0)))
    seg_arr = np.array(segments)

    # Build mesh (same as fem_scattered_coil.py without Kelvin)
    air_sphere = Sphere(Pnt(0, 0, 0), R_air)
    air_sphere.name = "air"
    wp_cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)
    for f in wp_cyl.faces:
        f.name = "sibc"
        f.maxh = R_wp / 4
    for f in air_sphere.faces:
        f.name = "outer"
    air = air_sphere - wp_cyl
    air.name = "air"
    geo = OCCGeometry(Glue([air]))
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=0.012, grading=0.3)
    mesh = Mesh(ngmesh)
    mesh.Curve(2)

    print(f"Mesh: {mesh.ne} tets, {mesh.nv} verts")

    # H1 interpolation
    nv = mesh.nv
    pts = np.zeros((nv, 3))
    for i in range(nv):
        p = mesh.vertices[i].point
        pts[i] = [p[0], p[1], p[2]]

    H_vals = h_segments_batch(seg_arr, pts, current=I_total)
    A_vals = a_segments_batch(seg_arr, pts, current=I_total)

    fes_h1 = H1(mesh, order=1)
    gf_Hx = GridFunction(fes_h1)
    gf_Hy = GridFunction(fes_h1)
    gf_Hz = GridFunction(fes_h1)
    gf_Ax = GridFunction(fes_h1)
    gf_Ay = GridFunction(fes_h1)
    gf_Az = GridFunction(fes_h1)

    gf_Hx.vec.FV().NumPy()[:nv] = H_vals[:, 0]
    gf_Hy.vec.FV().NumPy()[:nv] = H_vals[:, 1]
    gf_Hz.vec.FV().NumPy()[:nv] = H_vals[:, 2]
    gf_Ax.vec.FV().NumPy()[:nv] = A_vals[:, 0]
    gf_Ay.vec.FV().NumPy()[:nv] = A_vals[:, 1]
    gf_Az.vec.FV().NumPy()[:nv] = A_vals[:, 2]

    # Sample points on workpiece surface
    n_sample = 100
    errs_H = []
    errs_A = []
    for i in range(n_sample):
        # Random point on cylinder lateral surface
        theta = 2 * math.pi * i / n_sample
        z_rand = (np.random.random() - 0.5) * H_wp * 0.9
        x = R_wp * math.cos(theta) * 1.001  # slightly outside
        y = R_wp * math.sin(theta) * 1.001
        z = z_rand

        try:
            mip = mesh(x, y, z)
        except Exception:
            continue

        # Biot-Savart direct
        H_direct = h_segments_batch(seg_arr, np.array([[x, y, z]]), I_total)[0]
        A_direct = a_segments_batch(seg_arr, np.array([[x, y, z]]), I_total)[0]

        # H1 interpolated
        H_interp = np.array([gf_Hx(mip), gf_Hy(mip), gf_Hz(mip)])
        A_interp = np.array([gf_Ax(mip), gf_Ay(mip), gf_Az(mip)])

        err_H = np.linalg.norm(H_interp - H_direct) / max(np.linalg.norm(H_direct), 1e-30)
        err_A = np.linalg.norm(A_interp - A_direct) / max(np.linalg.norm(A_direct), 1e-30)
        errs_H.append(err_H)
        errs_A.append(err_A)

    errs_H = np.array(errs_H)
    errs_A = np.array(errs_A)

    print(f"\nH1 interpolation error on workpiece surface ({len(errs_H)} samples):")
    print(f"  H_inc: mean={np.mean(errs_H)*100:.2f}%, max={np.max(errs_H)*100:.2f}%")
    print(f"  A_inc: mean={np.mean(errs_A)*100:.2f}%, max={np.max(errs_A)*100:.2f}%")

    # Also check: what fraction of |A| is normal component on the cylinder?
    n_mesh = specialcf.normal(3)
    A_inc_cf = CF((gf_Ax, gf_Ay, gf_Az))

    # Sample normal component fraction
    fracs = []
    for i in range(n_sample):
        theta = 2 * math.pi * i / n_sample
        x = R_wp * math.cos(theta) * 1.001
        y = R_wp * math.sin(theta) * 1.001
        z = 0.0  # equator
        try:
            mip = mesh(x, y, z)
        except Exception:
            continue
        A = np.array([gf_Ax(mip), gf_Ay(mip), gf_Az(mip)])
        n = np.array([x, y, 0]) / np.sqrt(x**2 + y**2)  # analytical normal
        An = np.dot(A, n)
        frac = abs(An) / max(np.linalg.norm(A), 1e-30)
        fracs.append(frac)

    fracs = np.array(fracs)
    print(f"\n  |A_n|/|A| on cylinder equator: mean={np.mean(fracs)*100:.1f}%, "
          f"max={np.max(fracs)*100:.1f}%")

    # Direct integration comparison: ∫|A_t|² on sibc
    # Biot-Savart direct
    sibc_region = mesh.Boundaries("sibc")

    # H_inc energy on sibc: ∫ |H|² ds
    H_sq_interp = CF((gf_Hx, gf_Hy, gf_Hz)) * CF((gf_Hx, gf_Hy, gf_Hz))
    H_int = Integrate(H_sq_interp, mesh, BND, definedon=sibc_region).real
    A_wp = Integrate(CF(1), mesh, BND, definedon=sibc_region).real

    # H_inc_rms from interpolated GF
    H_rms_interp = math.sqrt(H_int / A_wp)

    # H_inc_rms from direct Biot-Savart at boundary element centroids
    print(f"\n  A_wp = {A_wp:.6e}")
    print(f"  |H_inc|_rms (H1 interp) = {H_rms_interp:.2f} A/m")


if __name__ == "__main__":
    main()
