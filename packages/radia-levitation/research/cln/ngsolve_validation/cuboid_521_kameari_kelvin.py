"""Historical field-recurrence experiment retired.

Only the independent geometry, mesh and linear-algebra helpers remain:
build_geo.
There is no reduction experiment entry point or validation result in this file.

Historical method attribution: Kameari, Sugahara.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, r"W:/00_CAE/Radia/01_GitHub/src/radia")

from netgen.occ import Box, Sphere, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
    Integrate, TaskManager, Periodic, ngsglobals,
)
from kelvin_geometry import add_kelvin_exterior_domain
from kelvin_material import make_kelvin_nu_cf, NU_0
from math import pi
import json
import time
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7

# Conductor
ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az

# Kelvin sphere parameters
R_K = 25e-3              # 25 mm sphere = 5x cuboid bounding sphere (~5x scale)
OFFSET = (4 * R_K, 0, 0) # Outer Kelvin sphere: well-separated from inner

H_COND = 0.5e-3          # 0.5 mm in conductor
H_AIR = 4.0e-3           # 4 mm in air (between cuboid and inner sphere boundary)
H_KELVIN = 8.0e-3        # 8 mm in Kelvin exterior
ORDER = 1                # production setting (Sugahara papers)
N_STAGES = 8

R0_anal = 48 / (sigma_Cu * V_cond * (ax**2 + ay**2))


def build_geo():
    """Cuboid + inner air sphere + Kelvin outer sphere (Sugahara two-sphere)."""
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.name = "conductor"
    cuboid.maxh = H_COND

    sphere_inner = Sphere(Pnt(0, 0, 0), R_K)
    for f in sphere_inner.faces:
        f.name = "kelvin_int"
    sphere_inner.name = "air"
    sphere_inner.maxh = H_AIR

    inner_air = sphere_inner - cuboid
    for f in inner_air.faces:
        # Re-tag the spherical face (after subtraction, naming is preserved)
        if f.name not in ("kelvin_int",):
            pass

    geo, info = add_kelvin_exterior_domain(
        [inner_air, cuboid],
        offset=OFFSET,
        R_K=R_K,
        inner_maxh=H_AIR,
        outer_maxh_factor=H_KELVIN/H_AIR,
    )
    return OCCGeometry(geo), info







# CLN experiment entry point retired; independent helpers above are retained.
