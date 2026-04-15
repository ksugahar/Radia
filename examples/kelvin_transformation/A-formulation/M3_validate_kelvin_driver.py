"""M3_validate_kelvin_driver.py

Smoke driver for the Layer 4 Kelvin validation harness
(``src/radia/kelvin_validate.compare_against_radia_self_inductance``).

Uses a square-cross-section gapped torus (matching
``validate_radia_vs_fem_kelvin.py``) as the single physical coil: the
FEM leg sees it as a meshed volume-J source; the Radia leg represents
the same current with a ``rad.ObjArcCur``. Both legs read the SAME
Kelvin-transformed mesh; they should agree up to the FEM discretization
error and the Radia reference represents the open-domain truth.

Run:
    python M3_validate_kelvin_driver.py
"""

from __future__ import annotations

import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
for p in (os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src')),
          os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src',
                                        'radia'))):
    if p not in sys.path:
        sys.path.insert(0, p)

import radia as rad

from ngsolve import CoefficientFunction as CF, x, y, sqrt
from netgen.occ import WorkPlane, Axes, Pnt, Dir, Sphere, Axis

from kelvin_validate import compare_against_radia_self_inductance


R_coil, a_coil, gap_deg = 0.030, 0.003, 5.0
arc_deg = 360.0 - gap_deg
R_K, offset = 0.06, (0.15, 0.0, 0.0)
# order=2 (p-refinement) converges much faster than h-refinement on a
# square-section coil: order=1 at this maxh gives -8% vs analytical,
# order=2 gives -0.3% at the SAME mesh (2026-04-15 sweep).
maxh, maxh_coil = 0.012, 0.003
order = 2
I_total = 1.0

MU_0 = 4e-7 * math.pi
J0 = I_total / (2 * a_coil) ** 2   # uniform current density


def build_inner_shapes():
    """Inner domain: sphere(R_K) with a square-section gapped torus in it."""
    wp = WorkPlane(Axes(p=Pnt(R_coil, 0, 0),
                         n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    wp.MoveTo(-a_coil, -a_coil)
    sq = wp.Rectangle(2 * a_coil, 2 * a_coil).Face()
    torus = sq.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), arc_deg)
    torus.name = "coil"
    torus.maxh = maxh_coil

    inner = Sphere(Pnt(0, 0, 0), R_K)
    inner.maxh = maxh
    for f in inner.faces:
        f.name = "kelvin_int"
    air = inner - torus
    air.name = "air"
    return [air, torus]


def build_radia_coil():
    rad.UtiDelAll()
    return rad.ObjArcCur(
        [0, 0, 0],
        [R_coil - a_coil, R_coil + a_coil],
        [0, math.radians(arc_deg)],
        2 * a_coil, 200, 'man', 'z', J0)


def main():
    print("=" * 64)
    print("M3 validation driver: FEM vs Radia on same Kelvin mesh")
    print(f"  square-section gapped torus R={R_coil*1e3:.0f}mm "
          f"a={a_coil*1e3:.1f}mm gap={gap_deg}deg")
    print(f"  Kelvin R_K={R_K*1e3:.0f}mm offset={offset}")
    print("=" * 64)

    r_xy = sqrt(x * x + y * y) + 1e-12
    J_src = J0 * CF((-y / r_xy, x / r_xy, 0))

    radia_handle = build_radia_coil()

    # Analytical reference for a circular-section torus; the square-section
    # coil here has a slightly larger effective inductance, so this is
    # used as a sanity benchmark (expect ~few % deviation).
    L_ref_circ = MU_0 * R_coil \
        * (math.log(8 * R_coil / a_coil) - 2) \
        * (1 - gap_deg / 360)

    t0 = time.perf_counter()
    out = compare_against_radia_self_inductance(
        inner_shapes=build_inner_shapes(),
        radia_handle=radia_handle,
        J_source_cf=J_src,
        kelvin_offset=offset, R_K=R_K, maxh=maxh,
        source_material="coil", order=order, I_total=I_total,
        L_analytical=L_ref_circ,
        verbose=True)
    print(f"  total time: {time.perf_counter()-t0:.1f}s")
    print("=" * 64)
    print("  Note: Radia leg is vertex-projected A (O(h) approx). Use")
    print("        L_analytical or a refinement sweep for tight checks.")


if __name__ == "__main__":
    main()
