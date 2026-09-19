"""Generate reproducible, synthetic beak-fin PEEC conductor STEPs.

Open-ended conductors for testing CAD perimeter recovery, surface-mesh
transitions, and fin-tip current/loss convergence. These are not as-built
induction coils or IH validation results. Coordinates are in mm.

Two variants share one profile: a straight z extrusion, and a revolution
about z that gives the partial-arc conductor a real beak-fin coil wraps
into. The curved one is what exercises the swept-section route end to end
(the straight one takes the z-plane sectioner instead).
"""

from __future__ import annotations

import argparse
from pathlib import Path

from build123d import (
    Axis,
    BuildLine,
    Plane,
    Polyline,
    ThreePointArc,
    export_step,
    extrude,
    make_face,
    revolve,
)

TIP_HALF_HEIGHT = (0.25**2 - 0.05**2) ** 0.5


def _beak_face(plane, radial_offset: float = 0.0):
    """Beak cross-section on ``plane``, shifted along the first axis."""
    def p(u, v):
        return (radial_offset + u, v)

    with BuildLine(plane) as outline:
        Polyline(p(-4.0, -2.0), p(1.6, -2.0), p(1.6, -0.7),
                 p(3.8, -TIP_HALF_HEIGHT))
        ThreePointArc(p(3.8, -TIP_HALF_HEIGHT), p(4.0, 0.0),
                      p(3.8, TIP_HALF_HEIGHT))
        Polyline(p(3.8, TIP_HALF_HEIGHT), p(1.6, 0.7), p(1.6, 2.0),
                 p(-4.0, 2.0), p(-4.0, -2.0))
    return make_face(outline.edges())


def make_beak_fin(length_mm: float = 60.0):
    """Prismatic copper-like conductor with a 0.25 mm rounded beak tip."""
    if length_mm <= 0:
        raise ValueError("length_mm must be positive")
    return extrude(_beak_face(Plane.XY), amount=length_mm)


def make_curved_beak_fin(major_radius_mm: float = 30.0,
                         arc_deg: float = 120.0):
    """The same section revolved about z, beak pointing radially outward.

    Revolution rather than a sweep: OCCT sweeps of a thin profile along a
    tight path degenerate, while a revolve is exact here.
    """
    if major_radius_mm <= 4.0:
        raise ValueError("major radius must clear the section half-width")
    if not 5.0 <= arc_deg <= 355.0:
        raise ValueError("arc_deg must be in [5, 355]")
    return revolve(_beak_face(Plane.XZ, major_radius_mm),
                   axis=Axis.Z, revolution_arc=arc_deg)


def main():
    fixtures = (Path(__file__).resolve().parents[2] / "tests"
                / "coil_from_cad" / "fixtures")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--length-mm", type=float, default=60.0)
    parser.add_argument("--arc-deg", type=float, default=None,
                        help="revolve through this angle instead of extruding")
    parser.add_argument("--major-radius-mm", type=float, default=30.0)
    args = parser.parse_args()

    if args.arc_deg is None:
        shape = make_beak_fin(args.length_mm)
        default = fixtures / "beak_fin_straight.step"
    else:
        shape = make_curved_beak_fin(args.major_radius_mm, args.arc_deg)
        default = fixtures / "beak_fin_curved.step"
    output = args.output or default
    output.parent.mkdir(parents=True, exist_ok=True)
    export_step(shape, str(output))
    print(f"STEP={output} solids={len(shape.solids())} "
          f"volume_mm3={shape.volume:.9g}")


if __name__ == "__main__":
    main()
