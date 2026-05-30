"""Generate the gapped-torus 2-port fixture via NGSolve OCC.

Geometry matches the PEEC golden fixture
(``src/radia/panels/samples/ih_peec_inductance_coil.step``):
  - Major radius R = 30 mm
  - Minor radius r = 3 mm
  - Gap = 5 deg (so the torus sweep angle is 355 deg)

Surface labels emitted in the .vol:
  - ``body``   -- lateral toroidal surface (Robin SIBC region)
  - ``source`` -- one of the two cap faces at the gap
  - ``sink``   -- the other cap face

Volume material:
  - ``coil``

Output:
  ``gapped_torus_2port.vol`` (next to this script)

Why NGSolve OCC instead of Cubit (2026-05-02):
  v1 prototyping of the cut-surface BEM port solver only needs a
  topologically-correct volumetric .vol with labelled caps.  NGSolve OCC
  via Revolve+Glue handles this in pure Python (no Cubit license).
  A Cubit-meshed .jou with copy-mesh-on-cap (for congruent caps) will
  follow when v1 is validated and we want hex/prism mesh quality.
"""
from __future__ import annotations

import math
import os
import sys

from netgen.occ import (
    OCCGeometry, WorkPlane, Axes, Axis, Pnt, Vec, X, Y, Z, Revolve,
)
from ngsolve import Mesh
from ngsolve import TaskManager


def make_gapped_torus(R_major: float = 0.030,
                       r_minor: float = 0.003,
                       gap_deg: float = 5.0,
                       maxh: float = 0.0015):
    """Build a 355deg gapped torus (volume), labelled for BEM-SIBC port solve.

    Classification of the three boundary faces:
      - body : the only face with mass (= surface area) much larger than
        the cap area pi*r_minor**2 -- this is the lateral toroidal face.
      - source / sink : the two flat caps; we use the y-coordinate of the
        face centre to distinguish them.  By the chosen sweep direction
        (+x -> +y), the START cap (sweep angle 0) sits at y=0 and the END
        cap sits at y < 0 (a 5deg gap goes back into the -y half-space).
    """
    sweep_deg = 360.0 - gap_deg
    # Profile: a disk of radius r_minor in a plane that contains the
    # z-axis -- specifically the y=0 plane, centred at (R, 0, 0).
    wp = WorkPlane(Axes(Pnt(R_major, 0, 0), n=Y, h=X))
    profile = wp.Circle(r_minor).Face()

    # Revolve the profile around the z-axis through `sweep_deg`.
    solid = Revolve(profile, Axis(Pnt(0, 0, 0), Z), sweep_deg)

    cap_area_expected = math.pi * r_minor ** 2
    cap_threshold = 5.0 * cap_area_expected   # body area >> cap area

    body_face = None
    source_face = None
    sink_face = None
    for f in solid.faces:
        if f.mass > cap_threshold:
            body_face = f
        elif f.center[1] >= -1e-9:
            # y >= 0 cap is the start of the sweep -> source
            source_face = f
        else:
            sink_face = f

    if body_face is None or source_face is None or sink_face is None:
        diag = [(f.center, f.mass) for f in solid.faces]
        raise RuntimeError(
            f"Failed to classify torus faces (got body={body_face is not None}, "
            f"source={source_face is not None}, sink={sink_face is not None}). "
            f"Face dump: {diag}")

    body_face.name = "body"
    source_face.name = "source"
    sink_face.name = "sink"
    solid.name = "coil"

    geo = OCCGeometry(solid)
    with TaskManager():
        return Mesh(geo.GenerateMesh(maxh=maxh))


def main():
    print("Generating gapped torus 2-port fixture...")
    mesh = make_gapped_torus(R_major=0.030, r_minor=0.003,
                              gap_deg=5.0, maxh=0.0015)
    materials = list(mesh.GetMaterials())
    bnds = list(mesh.GetBoundaries())
    print(f"  vertices: {mesh.nv}")
    print(f"  volume elements: {mesh.ne}")
    print(f"  surface elements: {mesh.nface}")
    print(f"  materials: {materials}")
    print(f"  BND labels: {bnds}")

    out_dir = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(out_dir, "gapped_torus_2port.vol")
    mesh.ngmesh.Save(out_path)
    print(f"  saved: {out_path}")


if __name__ == "__main__":
    sys.exit(main())
