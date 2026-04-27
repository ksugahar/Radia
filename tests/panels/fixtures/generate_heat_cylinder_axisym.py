"""Test fixture: 2D axisymmetric workpiece mesh in (r, z).

Half-section of the cylinder R=25mm H=25mm IH workpiece.

Boundaries::

    outer : r=R curve (heating face)
    top   : z=+H/2 curve
    bot   : z=-H/2 curve
    axis  : r=0 curve (symmetry; no BC needed)

Run::

    python tests/panels/fixtures/generate_heat_cylinder_axisym.py

to (re)write::

    tests/panels/fixtures/heat_workpiece_cylinder_R25_H25_axisym.vol
"""

from __future__ import annotations

import os
import sys


def main():
    from netgen.occ import WorkPlane, OCCGeometry

    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here,
                       "heat_workpiece_cylinder_R25_H25_axisym.vol")

    wp = WorkPlane()
    wp.MoveTo(0, -0.0125)
    wp.LineTo(0.025, -0.0125, name="bot")
    wp.LineTo(0.025, +0.0125, name="outer")
    wp.LineTo(0, +0.0125, name="top")
    wp.LineTo(0, -0.0125, name="axis")
    face = wp.Face()
    face.name = "workpiece"

    geo = OCCGeometry(face, dim=2)
    ngm = geo.GenerateMesh(maxh=0.003)
    ngm.Save(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    sys.exit(main() or 0)
