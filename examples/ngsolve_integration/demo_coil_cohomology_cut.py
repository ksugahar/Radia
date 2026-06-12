"""Demo: Total scalar potential with a gmsh-free cohomology cut (coil + iron).

This example demonstrates the cohomology-based total scalar potential
method for a coil + iron core magnetostatic problem -- WITHOUT gmsh.

Formulation:
    H = -grad(phi) + NI * h

where h is the cohomology cut basis function (HCurl, curl-free, UNIT
circulation around the coil loop).  The cut is computed by the pure-Python
``radia.cohomology`` engine (combinatorial Hodge Laplacian) -- no gmsh
``computeHomology``, no .msh -> .vol transfer.  No Biot-Savart needed.

Physical model:
    - Cylindrical iron core (mu_r = 1000) in the bore
    - Coil modeled as an annular through-hole in the mesh (NI ampere-turns)
    - Surrounding air box

The air + iron domain with the coil annulus removed is multiply-connected
(b1 = 1); the single cohomology cut carries the NI ampere-turns.

Verification:
    - Ampere's law:  oint H.dl = NI   (exact, by the unit-circulation cut)
    - on-axis B amplified by the iron vs. the same coil with an air core
      (mu_r = 1), via the SAME cohomology solver -- apples-to-apples

Reference:
    Pellikka et al., "Homology and Cohomology Computation in Finite
    Element Modeling," SIAM J. Sci. Comput., 2013.  (The radia.cohomology
    engine realises the same H^1 cut basis via the combinatorial Hodge
    Laplacian -- gmsh-free.)

run:  python demo_coil_cohomology_cut.py
"""

import sys
import os
import math

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src'))

MU_0 = 4 * np.pi * 1e-7


def build_ngsolve_mesh(core_radius=0.008, core_height=0.030,
                       coil_r_inner=0.010, coil_r_outer=0.015,
                       coil_height=0.020, air_size=0.06, maxh=0.006):
    """Air box + iron core, with the coil annulus left as a through-hole.

    The coil region (outer cylinder minus inner cylinder) is SUBTRACTED from
    the air, so the meshed domain (air + iron) is multiply-connected with
    b1 = 1.  No gmsh: the geometry is built with Netgen OCC and meshed by
    Netgen; the cohomology cut comes from ``radia.cohomology``.

    Returns
    -------
    ngsolve.Mesh
        Mesh with materials {'air', 'iron'} and boundary 'outer' on the box.
    """
    import ngsolve as ng
    from netgen.occ import Box, Cylinder, Pnt, Dir, Glue, OCCGeometry

    s = air_size / 2
    ax = Dir(0, 0, 1)

    air_box = Box(Pnt(-s, -s, -s), Pnt(s, s, s))
    iron = Cylinder(Pnt(0, 0, -core_height / 2), ax, r=core_radius, h=core_height)
    coil_outer = Cylinder(Pnt(0, 0, -coil_height / 2), ax, r=coil_r_outer, h=coil_height)
    coil_inner = Cylinder(Pnt(0, 0, -coil_height / 2), ax, r=coil_r_inner, h=coil_height)
    coil = coil_outer - coil_inner            # the annular washer to remove (the hole)

    iron.mat("iron")
    air = (air_box - coil) - iron             # air = box minus coil-hole minus iron
    air.mat("air")

    shape = Glue([air, iron])

    # Name the air-box outer faces "outer" (phi = 0); all other faces natural.
    for f in shape.faces:
        c = f.center
        if max(abs(c.x), abs(c.y), abs(c.z)) > 0.9 * s:
            f.name = "outer"

    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(shape).GenerateMesh(maxh=maxh))
    return mesh


def _ampere_meridian(field, mesh, x0, x1, z0, z1, n=240):
    """oint field.dl around a rectangle in the x-z (y=0) meridian half-plane
    that ENCIRCLES the coil cross-section.  This is the Ampere loop a coil's
    azimuthal current links (an XY circle around z does NOT link it).  For the
    H field this returns the enclosed NI; for the unit cut h0 it returns 1.
    """
    xs = np.linspace(x0, x1, n)
    zs = np.linspace(z0, z1, n)
    total = 0.0
    # bottom (z=z0, +x), top (z=z1, -x): integrate field_x dx
    for k in range(n - 1):
        xm = 0.5 * (xs[k] + xs[k + 1])
        dx_ = xs[k + 1] - xs[k]
        total += field(mesh(xm, 0.0, z0))[0] * dx_
        total -= field(mesh(xm, 0.0, z1))[0] * dx_
    # right (x=x1, +z), left (x=x0, -z): integrate field_z dz
    for k in range(n - 1):
        zm = 0.5 * (zs[k] + zs[k + 1])
        dz_ = zs[k + 1] - zs[k]
        total += field(mesh(x1, 0.0, zm))[2] * dz_
        total -= field(mesh(x0, 0.0, zm))[2] * dz_
    return float(total)


def solve(N=100, I_coil=1.0, core_radius=0.008, core_height=0.030,
          coil_r_inner=0.010, coil_r_outer=0.015, coil_height=0.020,
          air_size=0.06, mu_r=1000.0, maxh=0.006, order=2):
    """Build the mesh, compute the gmsh-free cohomology cut, solve, verify."""
    import ngsolve as ng
    from ngsolve import Integrate, InnerProduct, curl, dx
    from radia.cohomology_cut import CohomologyCutSolver

    NI = N * I_coil
    mesh = build_ngsolve_mesh(core_radius, core_height, coil_r_inner,
                              coil_r_outer, coil_height, air_size, maxh)

    # meridian Ampere rectangle: inner edge in the air gap (iron < x < coil),
    # encircling the coil cross-section, staying inside the air box.
    ax0 = 0.5 * (core_radius + coil_r_inner)         # between iron and coil
    ax1 = coil_r_outer + 0.4 * (air_size / 2 - coil_r_outer)
    az = coil_height / 2 + 0.4 * (air_size / 2 - coil_height / 2)

    with ng.TaskManager():
        solver = CohomologyCutSolver()
        b1 = solver.setup_from_mesh(mesh)            # cohomology cut via radia.cohomology

        # cut basis: curl-free + UNIT circulation around the coil (meridian loop)
        h0 = solver.get_cohomology_basis()[0]
        curl_norm = Integrate(InnerProduct(curl(h0), curl(h0)) * dx, mesh)
        h_circ = abs(_ampere_meridian(h0, mesh, ax0, ax1, -az, az))

        zs = [0.0, 0.010, 0.020]

        # (1) coil + iron core (mu_r):
        solver.solve([NI], {'iron': mu_r}, order=order, dirichlet='outer')
        H = solver.get_H()
        B = solver.get_B()
        amp = abs(_ampere_meridian(H, mesh, ax0, ax1, -az, az))   # = NI carried
        bz_iron = [float(B(mesh(0.0, 0.0, z))[2]) for z in zs]

        # (2) SAME coil, air core (mu_r = 1) -- same solver, apples-to-apples:
        solver.solve([NI], {'iron': 1.0}, order=order, dirichlet='outer')
        B_air = solver.get_B()
        bz_air = [float(B_air(mesh(0.0, 0.0, z))[2]) for z in zs]

    return {
        "b1": int(b1),
        "cut_curl_l2": float(curl_norm),            # ~0 (curl-free)
        "cut_circulation": float(h_circ),           # ~1 (unit circulation)
        "ampere_circulation": float(amp),           # oint H.dl
        "ampere_NI": float(NI),
        "ampere_rel_err": float(abs(amp - NI) / NI),
        "z_points": zs,
        "bz_coil_iron": bz_iron,                    # Tesla, mu_r iron core
        "bz_coil_air": bz_air,                      # Tesla, air core (same solver)
    }


def main():
    print("=" * 70)
    print("  Total scalar potential with a gmsh-free cohomology cut")
    print("  (coil + iron;  H = -grad(phi) + NI * h,  cut via radia.cohomology)")
    print("=" * 70)
    print()

    r = solve()

    print(f"  multiply-connected air+iron (coil hole): b1 = {r['b1']}  (one coil loop)")
    print(f"  cut basis h0:  ||curl h0||^2 = {r['cut_curl_l2']:.2e}  (curl-free), "
          f"oint h0.dl = {r['cut_circulation']:.4f}  (unit)")
    print(f"  Ampere's law:  oint H.dl = {r['ampere_circulation']:.2f}  "
          f"(NI = {r['ampere_NI']:.1f},  rel err {r['ampere_rel_err']:.2e})")
    print()
    print(f"  {'z [mm]':>8s} {'B_z coil+iron':>16s} {'B_z coil+air':>16s} {'ratio':>7s}")
    print("  " + "-" * 50)
    for z, bf, ba in zip(r['z_points'], r['bz_coil_iron'], r['bz_coil_air']):
        ratio = bf / ba if abs(ba) > 1e-12 else float('nan')
        print(f"  {z*1000:8.1f} {bf:16.6f} {ba:16.6f} {ratio:7.1f}")
    print()
    print("  => the gmsh-free cohomology cut carries the NI ampere-turns and makes")
    print("     the magnetic scalar potential single-valued in the coil-linking")
    print("     region; the iron core amplifies the on-axis field (ratio > 1).")


if __name__ == '__main__':
    main()
