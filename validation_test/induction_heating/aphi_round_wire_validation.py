"""Does the A-V eddy solver reproduce the exact round-wire impedance?

Gate 3 wants an independent reference.  A reference that has not itself been
checked is worth nothing, so before the solver is pointed at the beak fin it
is pointed at the one geometry with a closed-form answer: a solid round wire
carrying a terminal current, whose per-unit-length impedance is the Bessel
solution in `radia.analytical_formulas.conductor_impedance`.

The sweep is over ``a / delta``, which is the only dimensionless group that
matters here.  At small ``a / delta`` the current fills the wire and R tends
to the DC value; at large ``a / delta`` it crowds into the skin and R grows
like ``a / (2 delta)``.  A formulation that gets the gauge or the drive wrong
fails at one end or the other, so both ends are swept.

Only the resistance is compared.  The analytic value is the INTERNAL impedance
per unit length of an infinite wire, while the computed terminal impedance also
contains the external inductance of a finite wire inside its box -- an order
larger and entirely set by the truncation.  The reactance is recorded so the
run is reproducible, and it is not an acceptance quantity.

The resistance is a real test even so.  It is fixed by where the current sits
in the cross-section, which is the whole question a surface-impedance method
assumes the answer to, and at low frequency it must collapse onto
``L / (sigma A)`` of the meshed conductor with no freedom left at all.
"""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))

MU0 = 4e-7 * np.pi


def build_wire_mesh(radius_m, length_m, box_factor, maxh_conductor,
                    maxh_air, layers, growth, curve_order=1):
    """A wire on the z axis inside an air box, with a graded skin layer."""
    from netgen.occ import Box, Cylinder, Glue, OCCGeometry, Pnt, Z

    # The wire spans the box in z, so its two end discs lie in the same planes
    # as the air annuli around them.  A centre-of-face test cannot separate
    # those: an annulus is symmetric, so its centre is the axis too.  Name each
    # solid's faces before gluing instead, where the ownership is unambiguous.
    half = box_factor * radius_m
    tol = 1e-9 * max(1.0, length_m)

    wire = Cylinder(Pnt(0, 0, 0), Z, r=radius_m, h=length_m)
    wire.mat("conductor")
    wire.maxh = maxh_conductor
    for face in wire.faces:
        z = face.center[2]
        face.name = ("source" if abs(z) < tol
                     else "sink" if abs(z - length_m) < tol
                     else "wire_lateral")

    air = Box(Pnt(-half, -half, 0), Pnt(half, half, length_m)) - wire
    air.mat("air")
    air.maxh = maxh_air
    for face in air.faces:
        centre = face.center
        on_box = (abs(abs(centre[0]) - half) < tol
                  or abs(abs(centre[1]) - half) < tol
                  or abs(centre[2]) < tol
                  or abs(centre[2] - length_m) < tol)
        # The end-plane annuli are part of the truncation, not terminals.
        face.name = "outer" if on_box else "wire_lateral"

    shape = Glue([air, wire])
    geo = OCCGeometry(shape)
    ngmesh = geo.GenerateMesh(maxh=maxh_air)
    # Without curving the cylinder is a prism, and the comparison is against a
    # true circle: the error then sits in the geometry, not the solve, and
    # raising the polynomial order does nothing for it.
    if layers > 0:
        ngmesh.BoundaryLayer(boundary="wire_lateral",
                             thickness=[growth ** i for i in range(layers)],
                             material="conductor", domains="conductor",
                             outside=False)
    import ngsolve as ng

    mesh = ng.Mesh(ngmesh)
    if curve_order > 1:
        mesh.Curve(int(curve_order))
    return mesh


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius-mm", type=float, default=1.0)
    parser.add_argument("--length-mm", type=float, default=6.0)
    parser.add_argument("--sigma", type=float, default=5.8e7)
    parser.add_argument("--a-over-delta", type=float, nargs="+",
                        default=[0.5, 1.0, 2.0, 4.0])
    parser.add_argument("--order", type=int, default=3)
    parser.add_argument("--box-factor", type=float, default=6.0)
    parser.add_argument("--maxh-conductor-mm", type=float, default=0.35)
    parser.add_argument("--maxh-air-mm", type=float, default=1.5)
    parser.add_argument("--layers", type=int, default=0)
    parser.add_argument("--growth", type=float, default=0.5)
    parser.add_argument("--curve-order", type=int, default=3,
                        help="curve the faceted cylinder; the geometry is "
                             "the dominant error against a true circle")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    from radia.analytical_formulas.conductor_impedance import (
        cylinder_ac_impedance, cylinder_dc_resistance)
    from radia.eddy_aphi import solve_eddy_aphi

    radius = args.radius_mm * 1e-3
    length = args.length_mm * 1e-3
    mesh = build_wire_mesh(radius, length, args.box_factor,
                           args.maxh_conductor_mm * 1e-3,
                           args.maxh_air_mm * 1e-3, args.layers, args.growth,
                           curve_order=args.curve_order)
    print(json.dumps({"phase": "mesh", "ne": mesh.ne, "nv": mesh.nv,
                      "materials": list(mesh.GetMaterials()),
                      "boundaries": sorted(set(mesh.GetBoundaries()))}),
          flush=True)

    r_dc = cylinder_dc_resistance(radius, args.sigma)
    rows = []
    for ratio in args.a_over_delta:
        delta = radius / ratio
        omega = 2.0 / (delta * delta * MU0 * args.sigma)
        frequency = omega / (2.0 * np.pi)
        t0 = time.perf_counter()
        result = solve_eddy_aphi(mesh, conductor="conductor", source="source",
                                 sink="sink", frequency_hz=frequency,
                                 sigma=args.sigma, order=args.order,
                                 print_dofs=True)
        seconds = time.perf_counter() - t0
        z_num = result.terminal_impedance_ohm / length
        z_exact = cylinder_ac_impedance(radius, args.sigma, omega)
        row = {
            "a_over_delta": ratio,
            "frequency_hz": frequency,
            "skin_depth_mm": delta * 1e3,
            "seconds": seconds,
            "ndof": result.ndof,
            "R_numeric_ohm_per_m": z_num.real,
            "R_exact_ohm_per_m": float(np.real(z_exact)),
            "R_relative_error": float(z_num.real / np.real(z_exact) - 1.0),
            "X_numeric_ohm_per_m": z_num.imag,
            "internal_X_exact_ohm_per_m": float(np.imag(z_exact)),
            "R_over_Rdc_numeric": z_num.real / r_dc,
            "R_over_Rdc_exact": float(np.real(z_exact) / r_dc),
        }
        rows.append(row)
        print(json.dumps({"phase": "level", "a_over_delta": ratio,
                          "R_rel": row["R_relative_error"],
                          "seconds": round(seconds, 1)}), flush=True)

    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    worst_r = max(abs(r["R_relative_error"]) for r in rows)
    report = {
        "schema": "radia.eddy_aphi_round_wire.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": head,
        "host": platform.node(),
        "question": ("before the A-V solver is used as the fin's independent "
                     "reference, does it reproduce the exact Bessel impedance "
                     "of a round wire across the skin-effect range?"),
        "geometry": {"radius_mm": args.radius_mm, "length_mm": args.length_mm,
                     "box_factor": args.box_factor,
                     "maxh_conductor_mm": args.maxh_conductor_mm,
                     "maxh_air_mm": args.maxh_air_mm,
                     "order": args.order, "n_elements": mesh.ne},
        "sigma_S_per_m": args.sigma,
        "R_dc_ohm_per_m": r_dc,
        "levels": rows,
        "worst_R_relative_error": worst_r,
        "not_claimed": ("the reactance is not compared: the computed terminal "
                        "impedance carries the external inductance of a "
                        "finite wire in a box, which the analytic internal "
                        "impedance does not contain"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps({"phase": "complete", "worst_R_relative_error": worst_r,
                                    "output": str(args.output)}))


if __name__ == "__main__":
    main()
