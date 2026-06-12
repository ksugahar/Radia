"""T-Omega magnetostatics with a gmsh-free cohomology CUT -- the wire field.

The PRIMARY use of cohomology in the repo is the T-Omega total-scalar-
potential CUT (Kotiuga 1987 / Bossavit / Pellikka 2013), of which the
Clebsch current-linking demo is just one instance.  The total scalar
potential formulation

    H = -grad(phi) + sum_k NI_k h_k

uses cohomology cut basis functions h_k (curl-free, UNIT circulation around
the k-th current loop) to make the otherwise MULTI-VALUED magnetic scalar
potential single-valued in a current-linking (multiply-connected) region.
The cuts are computed by the pure-Python ``radia.cohomology`` engine -- no
gmsh ``computeHomology``, no .msh -> .vol transfer.

Verification (this script): a straight wire carrying current I along z
through an annular air region (b1 = 1).  The cohomology cut carries the I
ampere-turns, and the solved T-Omega field reproduces

    H_phi(r) = I / (2 pi r)      (Ampere's wire field)
    oint H.dl = I                (Ampere's law, exact by the unit-circulation cut)

run:  python tomega_wire.py
"""
import math

import numpy as np


def solve(I=100.0, R_out=0.05, r_w=0.005, L=0.10, maxh=0.008, order=2):
    """Solve the straight-wire T-Omega problem and return the verification dict."""
    import ngsolve as ng
    from netgen.occ import Cylinder, Pnt, Dir, OCCGeometry
    from radia.cohomology_cut import CohomologyCutSolver
    from radia.cohomology import circulation

    ax = Dir(0, 0, 1)
    outer = Cylinder(Pnt(0, 0, -L / 2), ax, r=R_out, h=L)
    inner = Cylinder(Pnt(0, 0, -L / 2), ax, r=r_w, h=L)
    air = (outer - inner)
    air.mat("air")
    for f in air.faces:                      # inner wall = natural; everything else = phi=0
        c = f.center
        rr = math.hypot(c.x, c.y)
        f.name = "wire" if abs(rr - r_w) < 0.3 * r_w else "outer"

    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(air).GenerateMesh(maxh=maxh))
        solver = CohomologyCutSolver()
        b1 = solver.setup_from_mesh(mesh)           # cohomology cut via radia.cohomology
        solver.solve([I], dirichlet="outer", order=order)
        H = solver.get_H()

        # H_phi(r) at mid-height (z=0) vs the analytic wire field I/(2 pi r).
        # At (r,0,0) the azimuthal direction phi_hat = +y, so H_phi = H_y.
        rs = np.linspace(1.6 * r_w, 0.8 * R_out, 8)
        errs = []
        for r in rs:
            Hy = H(mesh(r, 0.0, 0.0))[1]
            ana = I / (2 * math.pi * r)
            errs.append(abs(Hy - ana) / ana)
        wire_field_err = float(np.mean(errs))

        # Ampere's law: oint H.dl around a mid-radius circle (= the NI carried by the cut).
        amp = circulation(H, mesh, 0.0, 0.0, 0.5 * (r_w + R_out))

    return {
        "b1": int(b1),
        "wire_field_err": wire_field_err,           # H_phi vs I/(2 pi r), mean over r
        "ampere_circulation": float(amp),           # oint H.dl
        "ampere_NI": float(I),
        "ampere_rel_err": float(abs(amp - I) / I),
    }


def main():
    r = solve()
    print("T-Omega cohomology cut (gmsh-free) -- straight-wire verification\n")
    print(f"  multiply-connected air (annulus): b1 = {r['b1']}  (one current loop)")
    print(f"  Ampere's law:  oint H.dl = {r['ampere_circulation']:.3f}  "
          f"(NI = {r['ampere_NI']:.1f},  rel err {r['ampere_rel_err']:.2e})")
    print(f"  wire field  :  H_phi vs I/(2 pi r)  mean rel err = "
          f"{r['wire_field_err']:.2e}")
    print("\n  => the cohomology cut carries the NI ampere-turns and makes the")
    print("     magnetic scalar potential single-valued in the current-linking")
    print("     region -- the standard T-Omega use of radia.cohomology.")


if __name__ == "__main__":
    main()
