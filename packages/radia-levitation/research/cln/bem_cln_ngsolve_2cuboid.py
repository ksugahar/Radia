#!/usr/bin/env python3
"""
bem_cln_ngsolve_2cuboid.py

Phase 3 Task B.2: NGSolve A-formulation full-FEM cross-check of the
3D 2-cuboid BEM-CLN polarizability prediction.

Setup:
    Two Cu cuboids of dimensions 5 x 2 x 1 mm at positions
    (-D/2, 0, 0) and (+D/2, 0, 0) with D = 15 mm separation.
    Surrounded by air + Kelvin sphere (open boundary).
    Applied uniform B_0 along y-axis.

For each test frequency, solve the eddy-current A-formulation
    (1/mu) curl curl A + j omega sigma A = J_source
with the magnetic source applied via the outer Dirichlet condition
on the Kelvin sphere, and extract the induced magnetic moment of
each cuboid via integration over its volume.

Cross-checks the BEM-CLN polarizability prediction
    c_n = alpha(s) B^local_n
    alpha_total = sum_n c_n / B_0
from bem_cln_2cuboid_3D.wls (Mathematica) against direct FEM solve.

Requires:
    - NGSolve 6.2.2603+ (with curvedelements, Kelvin support)
    - Optional: radia.sparsesolv_ngsolve (Compact AMS) for HCurl
      preconditioning

This is a foundation script.  Compute-intensive solves (skin-resolved
mesh + Kelvin sphere + 3-5 frequencies) are typically 10-30 min on
LAB hardware.

Status:
    [ ] script written
    [ ] tested on actual NGSolve
    [ ] cross-check vs Mathematica BEM-CLN
"""

# Note: this is a foundation / template script.  When executing, ensure
# NGSolve >= 6.2.2603 is available and the Kelvin transformation
# tooling is set up.

import math
import os
import sys


def build_geometry_two_cuboids():
    """
    Geometry: 2 Cu cuboids of 5 x 2 x 1 mm, separation D = 15 mm
    centre-to-centre along x.  Surrounded by air + Kelvin sphere.

    Returns the NGSolve OCC geometry object.
    """
    from netgen.occ import Box, Pnt, Sphere, OCCGeometry, Glue
    from ngsolve import Mesh

    # Cuboid dimensions
    a, b, c = 5e-3, 2e-3, 1e-3        # Cu cuboid: 5x2x1 mm
    D = 15e-3                          # centre-to-centre separation

    # Two cuboids centred at (-D/2, 0, 0) and (+D/2, 0, 0)
    cu1 = Box(Pnt(-D/2 - a/2, -b/2, -c/2),
              Pnt(-D/2 + a/2,  b/2,  c/2))
    cu1.mat("cu")
    cu1.maxh = 0.5e-3                  # skin-depth-resolved meshing target

    cu2 = Box(Pnt(D/2 - a/2, -b/2, -c/2),
              Pnt(D/2 + a/2,  b/2,  c/2))
    cu2.mat("cu")
    cu2.maxh = 0.5e-3

    # Air region + Kelvin sphere
    R_inner = 0.05                      # 50 mm sphere
    R_outer = 0.10                      # Kelvin reflection point
    sphere = Sphere(Pnt(0, 0, 0), R_outer)
    air = sphere - cu1 - cu2
    air.mat("air")
    air.maxh = 5e-3

    # Mark outer boundary for Kelvin / Dirichlet
    air.bc("kelvin_outer")

    shape = Glue([cu1, cu2, air])
    return OCCGeometry(shape)


def solve_eddy_2cuboid(geo, freq, sigma_cu=5.8e7, mu0=4*math.pi*1e-7,
                        B0_y=1.0, order=1):
    """
    Solve A-formulation eddy current at frequency freq.
    Apply B_0 along y via Dirichlet on Kelvin outer boundary.
    Extract magnetic moment of each cuboid.

    Returns:
        m1, m2: induced magnetic moments (complex) of cuboid 1 and 2.
        alpha1, alpha2: derived per-cuboid polarizabilities.
    """
    from ngsolve import (Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
                         CoefficientFunction, Integrate, curl,
                         x, y, z, dx, ds,
                         InnerProduct, grad,
                         Preconditioner, CGSolver)

    mesh = Mesh(geo.GenerateMesh(maxh=10e-3))
    print(f"Mesh: {mesh.ne} elements, {mesh.nv} vertices")

    omega = 2 * math.pi * freq
    skin_depth = math.sqrt(2 / (omega * mu0 * sigma_cu))
    print(f"  freq = {freq:.2e} Hz, skin depth = {skin_depth*1e3:.4f} mm")

    # Material coefficients (complex sigma in conductor)
    sigma_cf = CoefficientFunction([sigma_cu if mat == "cu" else 0
                                     for mat in mesh.GetMaterials()])
    nu_cf = CoefficientFunction([1/mu0 for mat in mesh.GetMaterials()])

    # HCurl space with Dirichlet on outer Kelvin boundary
    fes = HCurl(mesh, order=order, dirichlet="kelvin_outer",
                complex=True)
    print(f"  HCurl DOF: {fes.ndof}")

    # Applied A: A_z = -B_0_y * x for uniform B_y (in y direction)
    # Wait: B_y = -dA_z/dx (in 2D thinking) or (curl A)_y in 3D
    # For uniform B_y, choose A = (0, 0, -B_0_y * x) or A = (-B_0_y * z, 0, 0)
    # The simplest is A_ext = (0, 0, -B_0_y * x)
    A_ext = CoefficientFunction((0, 0, -B0_y * x))

    # Bilinear form: (nu curl A) . curl v + j omega sigma A . v
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += nu_cf * InnerProduct(curl(u), curl(v)) * dx
    a += 1j * omega * sigma_cf * InnerProduct(u, v) * dx("cu")

    # Right-hand side: source J = -j omega sigma A_ext on conductor
    f = LinearForm(fes)
    f += -1j * omega * sigma_cf * InnerProduct(A_ext, v) * dx("cu")

    with TaskManager():
        a.Assemble()
        f.Assemble()

        # Solve
        gfu = GridFunction(fes)
        gfu.Set(A_ext, definedon=mesh.Boundaries("kelvin_outer"))   # boundary condition

        # Solve with direct/iterative
        res = f.vec.CreateVector()
        res.data = f.vec - a.mat * gfu.vec

        inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        gfu.vec.data += inv * res

        # Total A field
        A_total = gfu + A_ext

        # Extract induced magnetic moment via volume integration of (r x J)
        # For cuboid k centred at r_k: m_k = (1/2) integral_V (r - r_k) x J dV
        # where J = -j omega sigma A_total in conductor
        J = -1j * omega * sigma_cf * A_total

        # Cuboid centres
        r1_x, r1_y, r1_z = -7.5e-3, 0, 0
        r2_x, r2_y, r2_z =  7.5e-3, 0, 0

        # Magnetic moment of cuboid 1 (only y-component is significant by symmetry)
        # m_y = (1/2) integral [(z - r_z) J_x - (x - r_x) J_z] dV
        m1_y_cf = 0.5 * ((z - r1_z) * J[0] - (x - r1_x) * J[2])
        m2_y_cf = 0.5 * ((z - r2_z) * J[0] - (x - r2_x) * J[2])

        # Hmm — using cuboid-specific markers would be cleaner.  For now,
        # integrate over the "cu" material region (both cuboids combined)
        # and split by spatial position with a CF mask.
        # Cuboid 1 mask: |x - r1_x| < 3 mm (covers cuboid 1, excludes 2)
        from ngsolve import IfPos
        in_cu1 = IfPos(2.5e-3 - (x - r1_x), IfPos(2.5e-3 + (x - r1_x), 1, 0), 0)
        in_cu2 = IfPos(2.5e-3 - (x - r2_x), IfPos(2.5e-3 + (x - r2_x), 1, 0), 0)

        m1_y = Integrate(m1_y_cf * in_cu1, mesh, definedon=mesh.Materials("cu"))
        m2_y = Integrate(m2_y_cf * in_cu2, mesh, definedon=mesh.Materials("cu"))

        # Polarizability: alpha = m_y * mu_0 / B_0_y (volume susceptibility scale)
        alpha1 = m1_y * mu0 / B0_y
        alpha2 = m2_y * mu0 / B0_y

        return m1_y, m2_y, alpha1, alpha2


def main():
    """Run the 2-cuboid eddy-current solve at several frequencies."""
    print("=== bem_cln_ngsolve_2cuboid.py ===")
    print("Phase 3 B.2: NGSolve A-formulation cross-check")
    print()

    try:
        from ngsolve import Mesh
        from ngsolve import TaskManager
        print("NGSolve detected.")
    except ImportError:
        print("ERROR: NGSolve not available.  Install with:")
        print("  pip install ngsolve")
        sys.exit(1)

    print("Building 2-cuboid geometry + Kelvin sphere...")
    geo = build_geometry_two_cuboids()

    # Test frequencies (wall band for 5x2x1 mm Cu)
    freq_list = [1e5, 1e6, 1e7]
    print(f"Test frequencies: {freq_list}")
    print()

    results = []
    for freq in freq_list:
        print(f"=== Solving at f = {freq:.2e} Hz ===")
        try:
            m1, m2, alpha1, alpha2 = solve_eddy_2cuboid(
                geo, freq, order=1)
            alpha_total = alpha1 + alpha2
            print(f"  m1 = {m1:.4e}")
            print(f"  m2 = {m2:.4e}")
            print(f"  alpha_total = {alpha_total:.4e}")
            results.append({
                "freq": freq,
                "alpha1": alpha1,
                "alpha2": alpha2,
                "alpha_total": alpha_total,
            })
        except Exception as e:
            print(f"  ERROR: {e}")
            results.append({"freq": freq, "error": str(e)})

    # Compare with Mathematica BEM-CLN prediction
    # alpha_total prediction from bem_cln_2cuboid_3D.wls, Schur-F basis 3
    # (see data/bem_cln_2cuboid_3D.csv for reference values)
    print()
    print("=== Cross-check vs Mathematica BEM-CLN polarizability ===")
    # TODO: read prediction CSV and tabulate ratio FEM/Mathematica
    # at matching frequencies.

    # Write results CSV
    import csv
    out_csv = os.path.join(os.path.dirname(__file__), "bem_cln_ngsolve_2cuboid.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["freq", "alpha1_re", "alpha1_im", "alpha2_re", "alpha2_im",
                    "alpha_total_re", "alpha_total_im"])
        for r in results:
            if "error" in r:
                continue
            w.writerow([r["freq"],
                        r["alpha1"].real, r["alpha1"].imag,
                        r["alpha2"].real, r["alpha2"].imag,
                        r["alpha_total"].real, r["alpha_total"].imag])

    print(f"Wrote {out_csv}")
    print()
    print("=== DONE ===")


if __name__ == "__main__":
    main()
