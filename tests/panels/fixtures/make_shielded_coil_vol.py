#!/usr/bin/env python3
"""Generate mesh files for the active-shielding golden test (OCC, Cubit-free).

Writes four .vol meshes into the directory given as argv[1]:
  primary_cyl.vol   -- PRIMARY coil lateral surface (inner, r=0.15m, L=0.50m)
  shield_cyl.vol    -- SHIELD coil lateral surface  (outer, r=0.20m, L=0.55m)
  dsv.vol           -- DSV sphere volume (target region, r=0.05m)
  external.vol      -- LARGE annular shell OUTSIDE the shield
                       (r=0.235..0.27m, z=+-0.32m -> spans the full coil
                       length and beyond) -- where the stray field must
                       vanish.  CRITICAL: the external null region must
                       COVER the exterior, not just a thin mid-plane slice;
                       a too-small region only nulls the field locally and
                       can make it WORSE elsewhere (measured -- see the
                       active-shielding demo).  With this coverage the shield
                       gives a genuine ~25-140x stray reduction.

Run inside a SUBPROCESS (the golden test does this) so the NGSolve / Netgen
import stays out of the pytest process.  Deterministic for fixed geometry +
maxh, so golden bands are stable.
"""
import os
import sys

from netgen.occ import (Cylinder, Sphere, Pnt, Z, OCCGeometry, WorkPlane,
                        Revolve)
from ngsolve import TaskManager


def main(outdir):
    os.makedirs(outdir, exist_ok=True)

    with TaskManager():
        # Primary coil: inner cylinder lateral surface (2D-in-3D)
        a_p, L_p = 0.15, 0.50
        cyl_p = Cylinder(Pnt(0, 0, -L_p / 2), Z, r=a_p, h=L_p)
        lat_p = max(cyl_p.faces, key=lambda f: f.mass)
        OCCGeometry(lat_p).GenerateMesh(maxh=0.045).Save(
            os.path.join(outdir, "primary_cyl.vol"))

        # Shield coil: outer cylinder lateral surface
        a_s, L_s = 0.20, 0.55
        cyl_s = Cylinder(Pnt(0, 0, -L_s / 2), Z, r=a_s, h=L_s)
        lat_s = max(cyl_s.faces, key=lambda f: f.mass)
        OCCGeometry(lat_s).GenerateMesh(maxh=0.05).Save(
            os.path.join(outdir, "shield_cyl.vol"))

        # DSV: sphere volume inside the primary
        OCCGeometry(Sphere(Pnt(0, 0, 0), 0.05)).GenerateMesh(
            maxh=0.025).Save(os.path.join(outdir, "dsv.vol"))

        # External region: LARGE annular shell outside the shield, spanning
        # the full coil length and beyond (z=+-0.32, longer than shield
        # L=0.55 -> z=+-0.275).  Nulling the field over this whole exterior
        # (not a thin mid-plane slice) is what makes the shielding genuine.
        a_ext_in, a_ext_out, L_ext = 0.235, 0.27, 0.64
        cyl_out = Cylinder(Pnt(0, 0, -L_ext / 2), Z, r=a_ext_out, h=L_ext)
        cyl_in  = Cylinder(Pnt(0, 0, -L_ext / 2), Z, r=a_ext_in,  h=L_ext)
        shell = cyl_out - cyl_in
        OCCGeometry(shell).GenerateMesh(maxh=0.055).Save(
            os.path.join(outdir, "external.vol"))

    print(f"primary_vol={os.path.join(outdir, 'primary_cyl.vol')}")
    print(f"shield_vol={os.path.join(outdir, 'shield_cyl.vol')}")
    print(f"dsv_vol={os.path.join(outdir, 'dsv.vol')}")
    print(f"external_vol={os.path.join(outdir, 'external.vol')}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write("usage: make_shielded_coil_vol.py <outdir>\n")
        sys.exit(2)
    main(sys.argv[1])
