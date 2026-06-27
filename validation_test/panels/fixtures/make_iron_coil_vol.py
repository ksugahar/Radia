#!/usr/bin/env python3
"""Generate the MATERIAL-AWARE stream-function golden meshes (OCC, Cubit-free).

Writes three .vol meshes (gitignored) into argv[1]:
  coil_cyl_surf.vol -- cylinder LATERAL surface, radius 0.15 m, length 0.50 m
                       (the coil; same as make_streamfunction_vol)
  eval_dsv.vol      -- DSV sphere VOLUME, radius 0.05 m (the evaluation region)
  iron_kelvin.vol   -- a spherical IRON shell [0.35, 0.45] enclosing the coil +
                       DSV, with the Kelvin open-boundary ball (material 'kelvin'
                       + periodic kelvin_int<->kelvin_ext + a 'GND' vertex).
                       This is the --iron-vol the material-aware kernel consumes.

The iron Kelvin geometry mirrors the VERIFIED demo_oo /
docs/kelvin/kelvin_dtn_spectrum_archive.ipynb reduced-potential model (the coil at
r~0.29 fits inside the r=0.35 iron inner radius; the DSV r=0.05 sits inside).

Run inside a SUBPROCESS (the golden test does this) so NGSolve/Netgen stays out
of the pytest process.  Deterministic for fixed geometry + maxh.
"""
import os
import sys

import netgen.occ as occ
from netgen.occ import (Cylinder, Sphere, Pnt, Vec, Z, IdentificationType,
                        OCCGeometry)
from ngsolve import TaskManager

B_IN, C_OUT, R_OUT, OFFSET = 0.35, 0.45, 0.70, 2.0    # iron shell, Kelvin ball


def _iron_kelvin_geo():
    sb, sc = Sphere(Pnt(0, 0, 0), B_IN), Sphere(Pnt(0, 0, 0), C_OUT)
    outer = Sphere(Pnt(0, 0, 0), R_OUT)
    for f in outer.faces:
        f.name = "kelvin_int"
    core = sb; core.mat("vac")
    shell = (sc - sb); shell.mat("iron")
    out = (outer - sc); out.mat("vac")
    solids = [core, shell, out]
    kball = Sphere(Pnt(OFFSET, 0, 0), R_OUT)
    for f in kball.faces:
        f.name = "kelvin_ext"
    kball.mat("kelvin")
    gnd = occ.Vertex(Pnt(OFFSET, 0, 0)); gnd.name = "GND"
    fi = [f for s in solids for f in s.faces if f.name == "kelvin_int"][0]
    fe = [f for f in kball.faces if f.name == "kelvin_ext"][0]
    fi.Identify(fe, "kelvin", IdentificationType.PERIODIC,
                occ.gp_Trsf.Translation(Vec(OFFSET, 0, 0)))
    return occ.Glue(solids + [kball, gnd])


def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    with TaskManager():
        a, L = 0.15, 0.50
        cyl = Cylinder(Pnt(0, 0, -L / 2), Z, r=a, h=L)
        lateral = max(cyl.faces, key=lambda f: f.mass)
        coil_path = os.path.join(outdir, "coil_cyl_surf.vol")
        OCCGeometry(lateral).GenerateMesh(maxh=0.05).Save(coil_path)

        eval_path = os.path.join(outdir, "eval_dsv.vol")
        OCCGeometry(Sphere(Pnt(0, 0, 0), 0.05)).GenerateMesh(maxh=0.03).Save(eval_path)

        iron_path = os.path.join(outdir, "iron_kelvin.vol")
        OCCGeometry(_iron_kelvin_geo()).GenerateMesh(maxh=0.16).Save(iron_path)

    print(f"coil_vol={coil_path}")
    print(f"eval_vol={eval_path}")
    print(f"iron_vol={iron_path}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write("usage: make_iron_coil_vol.py <outdir>\n")
        sys.exit(2)
    main(sys.argv[1])
