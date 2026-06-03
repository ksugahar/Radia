#!/usr/bin/env python3
"""Generate the canonical stream-function panel sample meshes (OCC, Cubit-free).

Writes two .vol meshes (gitignored) into the directory given as argv[1]:
  coil_cyl_surf.vol  -- cylinder LATERAL surface (standalone 2D-in-3D),
                        radius 0.15 m, length 0.50 m, maxh 0.04  (the coil)
  eval_dsv.vol       -- DSV sphere VOLUME, radius 0.05 m, maxh 0.025
                        (the evaluation region)

Run inside a SUBPROCESS (the golden test does this) so the NGSolve / Netgen
import stays out of the pytest process (which crashes on the LAB MKL/Qt DLL
shadow).  Deterministic for a fixed geometry + maxh, so the golden bands are
stable.
"""
import os
import sys

from netgen.occ import Cylinder, Sphere, Pnt, Z, OCCGeometry
from ngsolve import TaskManager


def main(outdir):
    os.makedirs(outdir, exist_ok=True)

    with TaskManager():     # caller wraps the mesh-gen (TaskManager-Only policy)
        # coil: cylinder lateral surface -> standalone surface mesh
        a, L = 0.15, 0.50
        cyl = Cylinder(Pnt(0, 0, -L / 2), Z, r=a, h=L)
        lateral = max(cyl.faces, key=lambda f: f.mass)    # 2*pi*a*L > cap area
        coil_ng = OCCGeometry(lateral).GenerateMesh(maxh=0.04)
        coil_path = os.path.join(outdir, "coil_cyl_surf.vol")
        coil_ng.Save(coil_path)

        # eval: DSV sphere volume
        sph = Sphere(Pnt(0, 0, 0), 0.05)
        eval_ng = OCCGeometry(sph).GenerateMesh(maxh=0.025)
        eval_path = os.path.join(outdir, "eval_dsv.vol")
        eval_ng.Save(eval_path)

    print(f"coil_vol={coil_path}")
    print(f"eval_vol={eval_path}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write("usage: make_streamfunction_vol.py <outdir>\n")
        sys.exit(2)
    main(sys.argv[1])
