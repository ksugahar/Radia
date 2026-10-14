"""Test fixture generator: a 2-material (coil + workpiece) 3D volume
.vol that the radia-ih thermal step MUST reject.

The thermal step (``calc_heat.py`` / ``calc_heat_axisym.py``) targets the
WORKPIECE SOLID only.  A coil+workpiece mesh must error with a clear
"workpiece only" message rather than silently heating the coil region
as if it were steel (keiko 2026-05-31 incident: a ``coil_wrkpiace_vol``
mesh was passed to the thermal step).  This fixture is the negative-test
input for ``validation_test/panels/test_heat_workpiece_only.py``.

The workpiece solid matches ``generate_heat_cylinder.py`` exactly; a
disjoint "coil" box is added as a SECOND named volume material so the
.vol carries two materials.

Run::

    python validation_test/panels/fixtures/generate_heat_multimaterial.py

to (re)write::

    validation_test/panels/fixtures/heat_multimaterial_coil_wp.vol
"""
from __future__ import annotations

from ngsolve import TaskManager

import os
import sys


def main():
    from netgen.occ import Box, Cylinder, Glue, OCCGeometry, Pnt, Z

    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "heat_multimaterial_coil_wp.vol")

    # Workpiece cylinder -- identical to generate_heat_cylinder.py.
    wp = Cylinder((0, 0, -0.0125), Z, r=0.025, h=0.025)
    wp.solids.name = "workpiece"
    wp.faces.name = "outer"

    # A disjoint "coil" box so the mesh has TWO volume materials --
    # exactly the coil+workpiece case the thermal guard must reject.
    coil = Box(Pnt(0.05, -0.01, -0.01), Pnt(0.07, 0.01, 0.01))
    coil.solids.name = "coil"

    geo = OCCGeometry(Glue([wp, coil]))
    with TaskManager():
        ngm = geo.GenerateMesh(maxh=0.01)
        ngm.Save(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    sys.exit(main() or 0)
