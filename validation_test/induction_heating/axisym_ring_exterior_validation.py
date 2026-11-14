"""Validate the axisymmetric ring reference: drive, exterior term, skin.

`radia.eddy_axisym_ring` is the interior-resolved reference that the fin PEEC
is adjudicated against, so its own soundness is checked here against
references that share none of its code:

1. Direct-current inductance of a rectangular-section ring carrying
   ``J ~ 1/r`` (the loop-voltage distribution), by Gauss quadrature of the
   filament-loop mutual inductance (complete elliptic integrals).  This tests
   the stiffness operator, the drive and the exterior term together.
2. Invariance of R and L under the exterior radius.  An exact exterior makes
   the answer independent of where the air is cut.
3. Resistance of a revolved round wire against the exact Bessel internal
   impedance of a straight wire, at a radius where the curvature is a
   fraction of a percent.  This tests the sigma-mass operator in the skin
   regime, which the direct-current check does not reach.

Each row records the measured relative error.  No pass/fail threshold is
applied here: acceptance criteria are owned by the user.
"""
from __future__ import annotations

import argparse
import json
import math
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

MU0 = 4e-7 * math.pi
SIGMA = 5.8e7


def half_disk_mesh(conductor, exterior_radius, maxh_conductor, maxh_air,
                   curve_order=3):
    """Half disk about the axis holding one conductor ``("rect"|"circle", ...)``."""
    import ngsolve as ng
    from netgen.geom2d import SplineGeometry

    geo = SplineGeometry()
    big = exterior_radius
    p = [geo.AppendPoint(*q) for q in [(0, -big), (big, 0), (0, big)]]
    geo.Append(["spline3", p[0], geo.AppendPoint(big, -big), p[1]],
               bc="outer", leftdomain=1, rightdomain=0)
    geo.Append(["spline3", p[1], geo.AppendPoint(big, big), p[2]],
               bc="outer", leftdomain=1, rightdomain=0)
    geo.Append(["line", p[2], p[0]], bc="axis", leftdomain=1, rightdomain=0)
    kind, *args = conductor
    if kind == "rect":
        r0, r1, z0, z1 = args
        geo.AddRectangle((r0, z0), (r1, z1), leftdomain=2, rightdomain=1,
                         bc="surface")
    else:
        center, radius = args
        geo.AddCircle((center, 0.0), radius, leftdomain=2, rightdomain=1,
                      bc="surface")
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "conductor")
    geo.SetDomainMaxH(2, maxh_conductor)
    mesh = ng.Mesh(geo.GenerateMesh(maxh=maxh_air))
    mesh.Curve(curve_order)
    return mesh


def filament_inductance(r0, r1, z0, z1, points):
    """L of a rectangular ring with J ~ 1/r, by tensor Gauss quadrature.

    Source and target rules differ in order by one so no node pair
    coincides; the log singularity is integrable and the sequence converges.
    """
    from scipy.special import ellipe, ellipk

    def rule(n):
        x, w = np.polynomial.legendre.leggauss(n)
        r = 0.5 * (r1 - r0) * x + 0.5 * (r1 + r0)
        z = 0.5 * (z1 - z0) * x + 0.5 * (z1 + z0)
        rr, zz = np.meshgrid(r, z, indexing="ij")
        ww = np.outer(w * 0.5 * (r1 - r0), w * 0.5 * (z1 - z0)) / rr
        return rr.ravel(), zz.ravel(), (ww / ww.sum()).ravel()

    ra, za, wa = rule(points)
    rb, zb, wb = rule(points + 1)
    total = 0.0
    for start in range(0, ra.size, 512):     # bounded memory per block
        sl = slice(start, start + 512)
        a_r, a_z = ra[sl, None], za[sl, None]
        k2 = 4 * a_r * rb / ((a_r + rb) ** 2 + (a_z - zb) ** 2)
        k = np.sqrt(k2)
        mutual = MU0 * np.sqrt(a_r * rb) * ((2 / k - k) * ellipk(k2)
                                            - 2 / k * ellipe(k2))
        total += float(wa[sl] @ mutual @ wb)
    return total


def bessel_resistance_per_metre(frequency_hz, radius):
    """Exact internal resistance per metre of a straight round wire."""
    from scipy.special import jv

    k = np.sqrt(-1j * 2 * math.pi * frequency_hz * MU0 * SIGMA)
    z = k / (2 * math.pi * radius * SIGMA) * jv(0, k * radius) / jv(1, k * radius)
    return float(z.real)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--modes", type=int, default=30)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    from radia.eddy_axisym_ring import solve_axisym_ring

    rect = (0.040, 0.060, -0.010, 0.010)
    quadrature = {n: filament_inductance(*rect, n) for n in (60, 120, 240)}
    l_ref = quadrature[240]

    dc_rows = []
    for ext in (0.1, 0.15, 0.3, 0.6):
        for maxh in (0.0015, 0.00075):
            mesh = half_disk_mesh(("rect", *rect), ext, maxh, ext / 10)
            res = solve_axisym_ring(
                mesh, conductor="conductor", frequency_hz=1.0, sigma=SIGMA,
                exterior_radius=ext, exterior_modes=args.modes)
            dc_rows.append({
                "exterior_radius_m": ext, "maxh_conductor_m": maxh,
                "ndof": res.ndof, "inductance_H": res.inductance_H,
                "inductance_relative_error": res.inductance_H / l_ref - 1.0,
                "resistance_ohm": res.resistance_ohm,
            })
            print(json.dumps({"phase": "dc", **dc_rows[-1]}), flush=True)

    wire = 1e-3
    skin_rows = []
    for ring_radius in (0.1, 0.3):
        for frequency, maxh in ((1e3, wire / 8), (1e5, wire / 16),
                                (1e5, wire / 32)):
            ext = 2 * ring_radius
            mesh = half_disk_mesh(("circle", ring_radius, wire), ext, maxh,
                                  ring_radius / 4)
            res = solve_axisym_ring(
                mesh, conductor="conductor", frequency_hz=frequency,
                sigma=SIGMA, exterior_radius=ext, exterior_modes=args.modes)
            per_metre = res.resistance_ohm / (2 * math.pi * ring_radius)
            exact = bessel_resistance_per_metre(frequency, wire)
            skin_rows.append({
                "ring_radius_m": ring_radius, "wire_radius_m": wire,
                "frequency_hz": frequency, "maxh_conductor_m": maxh,
                "ndof": res.ndof, "resistance_per_metre": per_metre,
                "bessel_per_metre": exact,
                "relative_error": per_metre / exact - 1.0,
                "terminal_vs_energy_resistance": res.terminal_resistance_ohm
                / res.resistance_ohm - 1.0,
            })
            print(json.dumps({"phase": "skin", **skin_rows[-1]}), flush=True)

    fine = [r for r in dc_rows if r["maxh_conductor_m"] == 0.00075]
    spread = (max(r["inductance_H"] for r in fine)
              / min(r["inductance_H"] for r in fine) - 1.0)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    report = {
        "schema": "radia.axisym_ring_exterior_validation.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": head,
        "platform_class": platform.system(),
        "python": sys.version.split()[0],
        "exterior_modes": args.modes,
        "dc_inductance": {
            "section_m": dict(zip(("r0", "r1", "z0", "z1"), rect)),
            "filament_quadrature_H": {str(k): v for k, v in quadrature.items()},
            "levels": dc_rows,
            "exterior_radius_spread_finest": spread,
        },
        "revolved_wire_skin": skin_rows,
        "not_claimed": (
            "the round-wire rows carry the ring's own curvature, of order "
            "the wire radius over the ring radius; they test the skin "
            "operator, not a straight-wire identity"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps({"phase": "complete", "output": str(args.output)}))


if __name__ == "__main__":
    from ngsolve import TaskManager
    with TaskManager():
        main()
