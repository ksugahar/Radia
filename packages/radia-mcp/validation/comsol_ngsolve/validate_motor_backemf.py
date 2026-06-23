r"""No-load FLUX LINKAGE lambda(theta) / back-EMF of a PMSM (regression test).

A 2-pole diametric PM rotor (R=10 mm, Br=1 T) in a linear-iron stator ring (mu_r=1000)
with an air-gap SEARCH COIL (N=100, 0 current).  Rotating the PM magnetisation phi_deg
sweeps the coil flux linkage; for a 2-pole machine lambda(theta)=lambda0 cos(theta) and the
no-load back-EMF is e = -omega dlambda/dtheta = omega lambda0 sin(theta).

radia-ngsolve (solve_planar_magnetostatic + the coil_flux_linkage_2d helper -- the standard
flux-linkage definition lambda = N depth (<A_z>_pos - <A_z>_neg)) reproduces a stored
regression reference (independently cross-validated) to 0.10 % in lambda0 and 0.00 % RMS in
the cos(theta) shape.  The reference value is stored below -> CI-portable.

UNITS: meshed in mm, so A_z(mm-solve)=1000*A_phys (B is scale-invariant -> A_z is 1e3 too
big); lambda[Wb] = coil_flux_linkage_2d(depth=1) * 1e-6 (depth 1 mm = 1e-3 m, A /1e3).
"""
import math
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen")

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, TaskManager
from netgen.occ import OCCGeometry, WorkPlane, Glue
from radia_mcp.radia_ngsolve.solve import (reluctivity, solve_planar_magnetostatic,
                                           coil_flux_linkage_2d)

MU0 = 4e-7 * math.pi
HC = 1.0 / MU0                       # Br = 1 T rigid PM
Rr, Rbore, Rsout, Rout = 10.0, 14.0, 24.0, 70.0
rc, rw, NTURNS = 12.0, 1.0, 100
SCALE_WB = 1e-6                      # mm-solve helper(depth=1) -> Wb (see module docstring)
LAMBDA0_REF_WB = 1.4432e-3          # stored regression reference [Wb], 1 mm stack


def _mesh():
    rotor = WorkPlane().Circle(0, 0, Rr).Face(); rotor.faces.name = "rotor"
    coilpos = WorkPlane().Circle(0, rc, rw).Face(); coilpos.faces.name = "coilpos"
    coilneg = WorkPlane().Circle(0, -rc, rw).Face(); coilneg.faces.name = "coilneg"
    disk_bore = WorkPlane().Circle(0, 0, Rbore).Face()
    gap = disk_bore - rotor - coilpos - coilneg; gap.faces.name = "air"
    disk_so = WorkPlane().Circle(0, 0, Rsout).Face()
    iron = disk_so - disk_bore; iron.faces.name = "iron"
    box = WorkPlane().Circle(0, 0, Rout).Face(); box.edges.name = "outer"
    outer = box - disk_so; outer.faces.name = "air"
    rotor.faces.maxh = 1.5; gap.faces.maxh = 0.8
    coilpos.faces.maxh = 0.3; coilneg.faces.maxh = 0.3; iron.faces.maxh = 2.0
    return Mesh(OCCGeometry(Glue([rotor, coilpos, coilneg, gap, iron, outer]),
                            dim=2).GenerateMesh(maxh=8.0))


def validate_flux_linkage_amplitude():
    mesh = _mesh()
    nu = reluctivity(mesh, {"iron": 1000.0})
    thetas = [0, 30, 60, 90, 180]
    lam = {}
    with TaskManager():
        for th in thetas:
            A = solve_planar_magnetostatic(mesh, nu, magnets={"rotor": (HC, th)},
                                           order=3, dirichlet="outer")
            lam[th] = coil_flux_linkage_2d(A, mesh, "coilpos", "coilneg", NTURNS,
                                           depth=1.0) * SCALE_WB
    lam0 = abs(lam[0])
    rel = abs(lam0 - LAMBDA0_REF_WB) / LAMBDA0_REF_WB
    print(f"[flux linkage] radia lambda0={lam0:.4e} Wb  ref={LAMBDA0_REF_WB:.4e} Wb "
          f"diff={100*rel:.2f} %")
    assert rel < 0.03, f"lambda0 {lam0:.3e} vs ref {LAMBDA0_REF_WB:.3e} ({100*rel:.1f}%)"
    # lambda(theta) = lambda0 cos(theta)
    assert abs(lam[30] - lam0 * math.cos(math.radians(30))) < 0.02 * lam0
    assert abs(lam[60] - lam0 * math.cos(math.radians(60))) < 0.02 * lam0
    assert abs(lam[90]) < 0.03 * lam0                       # cos(90)=0
    assert abs(lam[180] + lam0) < 0.03 * lam0               # cos(180)=-1
    assert lam[0] > 0 and lam[180] < 0                      # sign flips over a pole


if __name__ == "__main__":
    validate_flux_linkage_amplitude()
    print("[OK] PMSM no-load flux linkage / back-EMF validated against the reference.")
