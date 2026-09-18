"""Generate a reproducible, synthetic beak-fin PEEC conductor STEP.

This is a straight, open-ended conductor for testing CAD perimeter recovery,
surface-mesh transitions, and fin-tip current/loss convergence. It is not an
as-built induction coil or an IH validation result. Coordinates are in mm.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from build123d import (
    BuildLine,
    Polyline,
    ThreePointArc,
    export_step,
    extrude,
    make_face,
)


def make_beak_fin(length_mm: float = 60.0):
    """Prismatic copper-like conductor with a 0.25 mm rounded beak tip."""
    if length_mm <= 0:
        raise ValueError("length_mm must be positive")
    tip_lower = (3.8, -(0.25**2 - 0.05**2) ** 0.5)
    tip_upper = (3.8, +(0.25**2 - 0.05**2) ** 0.5)
    with BuildLine() as outline:
        Polyline((-4.0, -2.0), (1.6, -2.0), (1.6, -0.7), tip_lower)
        ThreePointArc(tip_lower, (4.0, 0.0), tip_upper)
        Polyline(tip_upper, (1.6, 0.7), (1.6, 2.0),
                 (-4.0, 2.0), (-4.0, -2.0))
    return extrude(make_face(outline.edges()), amount=length_mm)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=(
        Path(__file__).resolve().parents[2] / "tests" / "coil_from_cad"
        / "fixtures" / "beak_fin_straight.step"))
    parser.add_argument("--length-mm", type=float, default=60.0)
    args = parser.parse_args()
    shape = make_beak_fin(args.length_mm)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    export_step(shape, str(args.output))
    print(f"STEP={args.output} solids={len(shape.solids())} "
          f"volume_mm3={shape.volume:.9g}")


if __name__ == "__main__":
    main()
