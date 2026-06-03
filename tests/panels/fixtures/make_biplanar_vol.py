#!/usr/bin/env python3
"""Generate the BIPLANAR stream-function sample meshes (OCC, Cubit-free).

A biplanar coil is the canonical flat-plate gradient / shim geometry: two
parallel square formers above and below the imaging volume, each carrying a
stream function.  This is the multi-surface generality test for the FE-direct
SF designer -- ONE surface mesh holds TWO disconnected plate components, and
``H1(coil, definedon=coil.Boundaries('.*'))`` + the connected-component ``abe``
edge-equipotential BC handle both plates with no special-casing.

Writes two .vol meshes (gitignored) into argv[1]:
  biplanar_coil.vol  -- two 0.30 x 0.30 m square plates at z = +/-0.10 m
                        (one mesh, two components), maxh 0.03  (the coil)
  biplanar_eval.vol  -- DSV sphere VOLUME at the midplane, radius 0.04 m,
                        maxh 0.02                                (eval region)

Run inside a SUBPROCESS (the golden test does) so the NGSolve / Netgen import
stays out of pytest (which crashes on the LAB MKL/Qt DLL shadow).  Deterministic
for fixed geometry + maxh, so the golden bands are stable.
"""
import os
import sys

from netgen.occ import WorkPlane, Axes, Pnt, Z, Sphere, OCCGeometry, Glue
from ngsolve import TaskManager


def _plate(zc, side=0.30):
    return WorkPlane(Axes((0, 0, zc), n=Z)).RectangleC(side, side).Face()


def main(outdir):
    os.makedirs(outdir, exist_ok=True)

    with TaskManager():     # caller wraps the mesh-gen (TaskManager-Only policy)
        # coil: two parallel plates as ONE surface mesh (two components)
        geo = OCCGeometry(Glue([_plate(0.10), _plate(-0.10)]))
        coil_ng = geo.GenerateMesh(maxh=0.03)
        coil_path = os.path.join(outdir, "biplanar_coil.vol")
        coil_ng.Save(coil_path)

        # eval: DSV sphere volume at the midplane between the plates
        sph = Sphere(Pnt(0, 0, 0), 0.04)
        eval_ng = OCCGeometry(sph).GenerateMesh(maxh=0.02)
        eval_path = os.path.join(outdir, "biplanar_eval.vol")
        eval_ng.Save(eval_path)

    print(f"coil_vol={coil_path}")
    print(f"eval_vol={eval_path}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write("usage: make_biplanar_vol.py <outdir>\n")
        sys.exit(2)
    main(sys.argv[1])
