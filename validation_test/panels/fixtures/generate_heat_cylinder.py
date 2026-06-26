"""Test fixture: workpiece-volume cylinder .vol for heat chain tests.

Geometry MUST match the IH panel sample's wp hole exactly so the
cross-mesh qsurf transfer in calc_heat lands on the same physical
surface.  IH samples (ih_peec_bem_coarse.jou,
ih_fem_kelvin_skin_fine.jou) define the workpiece as::

    create cylinder radius 0.025000 height 0.025000

Cubit's ``create cylinder`` is axis-aligned along Z and CENTERED on
the origin (z in [-h/2, +h/2]).  We mirror that exactly here.

Run::

    python validation_test/panels/fixtures/generate_heat_cylinder.py

to (re)write::

    validation_test/panels/fixtures/heat_workpiece_cylinder_R25_H25.vol

The .vol is the target of ``test_heat_chain_golden.py``.
"""
from __future__ import annotations

from ngsolve import TaskManager

import os
import sys


def main():
    from netgen.occ import Cylinder, Z, OCCGeometry

    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "heat_workpiece_cylinder_R25_H25.vol")

    cyl = Cylinder((0, 0, -0.0125), Z, r=0.025, h=0.025)
    cyl.faces.name = "outer"
    cyl.solids.name = "workpiece"

    geo = OCCGeometry(cyl)
    with TaskManager():
        ngm = geo.GenerateMesh(maxh=0.005)
        ngm.Save(out)
        print(f"wrote {out}")


if __name__ == "__main__":
    sys.exit(main() or 0)
