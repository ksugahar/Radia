"""Generate reproducible, synthetic beak-fin PEEC conductor STEPs.

Open-ended conductors for testing CAD perimeter recovery, surface-mesh
transitions, and fin-tip current/loss convergence. These are not as-built
induction coils or IH validation results. Coordinates are in mm.

Three variants share one profile. A straight z extrusion; a revolution
about z that gives the partial-arc conductor a real beak-fin coil wraps
into, which is what exercises the swept-section route end to end (the
straight one takes the z-plane sectioner instead); and a z loft whose
beak shortens along the sweep at a fixed tip radius, which is the only
one whose section CHANGES and so the only one that constrains the
station count.
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
    loft,
    make_face,
    revolve,
)

TIP_RADIUS = 0.25
TIP_HALF_HEIGHT = (TIP_RADIUS**2 - 0.05**2) ** 0.5


def _beak_face(plane, radial_offset: float = 0.0, tip_x: float = 4.0):
    """Beak cross-section on ``plane``, shifted along the first axis.

    ``tip_x`` moves how far the beak reaches while the tip circle keeps
    radius ``TIP_RADIUS``, so a loft between two of these varies the fin
    length alone.
    """
    centre = tip_x - TIP_RADIUS

    def p(u, v):
        return (radial_offset + u, v)

    with BuildLine(plane) as outline:
        Polyline(p(-4.0, -2.0), p(1.6, -2.0), p(1.6, -0.7),
                 p(centre + 0.05, -TIP_HALF_HEIGHT))
        ThreePointArc(p(centre + 0.05, -TIP_HALF_HEIGHT), p(tip_x, 0.0),
                      p(centre + 0.05, TIP_HALF_HEIGHT))
        Polyline(p(centre + 0.05, TIP_HALF_HEIGHT), p(1.6, 0.7), p(1.6, 2.0),
                 p(-4.0, 2.0), p(-4.0, -2.0))
    return make_face(outline.edges())


def make_beak_fin(length_mm: float = 60.0, tip_x_mm: float = 4.0):
    """Prismatic copper-like conductor with a 0.25 mm rounded beak tip.

    ``tip_x_mm`` varies how far the beak reaches while the tip circle keeps
    its radius, so a sweep over it changes the fin's length and aspect at a
    fixed local curvature -- which is what separates a geometry effect from
    the curvature effect the surface impedance is sensitive to.
    """
    if length_mm <= 0:
        raise ValueError("length_mm must be positive")
    if not 1.9 < tip_x_mm <= 4.0:
        raise ValueError("the beak tip must stay ahead of its 1.6 mm root")
    return extrude(_beak_face(Plane.XY, tip_x=tip_x_mm), amount=length_mm)


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


def make_tapered_beak_fin(length_mm: float = 60.0,
                          tip_x_end_mm: float = 2.2):
    """Straight sweep whose beak shortens from 4.0 mm to ``tip_x_end_mm``.

    The tip radius is unchanged at both ends, so the only thing varying
    along the sweep is the fin length -- which is what the station-count
    rule keys on.
    """
    if length_mm <= 0:
        raise ValueError("length_mm must be positive")
    if not 1.9 < tip_x_end_mm <= 4.0:
        raise ValueError("the beak tip must stay ahead of its 1.6 mm root")
    return loft([_beak_face(Plane.XY, tip_x=4.0),
                 _beak_face(Plane.XY.offset(length_mm), tip_x=tip_x_end_mm)])


def main():
    fixtures = (Path(__file__).resolve().parents[2] / "tests"
                / "coil_from_cad" / "fixtures")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--length-mm", type=float, default=60.0)
    parser.add_argument("--arc-deg", type=float, default=None,
                        help="revolve through this angle instead of extruding")
    parser.add_argument("--major-radius-mm", type=float, default=30.0)
    parser.add_argument("--tip-x-end-mm", type=float, default=None,
                        help="loft the beak down to this reach instead")
    args = parser.parse_args()

    if args.tip_x_end_mm is not None:
        shape = make_tapered_beak_fin(args.length_mm, args.tip_x_end_mm)
        default = fixtures / "beak_fin_tapered.step"
    elif args.arc_deg is None:
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
