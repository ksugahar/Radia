"""Gate 3: the beak fin against an interior-resolved A-V reference.

BEM-A and the surface PEEC both impose a Leontovich surface impedance.  They
agree with each other to 0.58% on the delivery metrics, but agreement between
two routes that share an assumption says nothing about the assumption.  This
driver solves the same fin with the conductor's interior resolved, so the
current distribution is an output.

The comparison quantities are the ones gate 3 names: terminal impedance, the
beak and tip loss shares, and the field on the workpiece-side probe line.  The
sweep is over frequency, which is also what makes the problem affordable --
the skin depth sets the mesh, and at 150 kHz in copper it is 0.17 mm on a
48 mm fin.  Sweeping downward walks the mesh from tractable to hard while
walking delta / feature from the regime where a surface impedance is a poor
assumption into the one where it is a good one, which is exactly where the two
families should start to disagree.
"""
from __future__ import annotations

import argparse
import json
import math
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
sys.path.insert(0, str(HERE))

MU0 = 4e-7 * math.pi


def build_fin_mesh(step_path, *, box_factor, maxh_conductor_m, maxh_air_m,
                   curve_order):
    """The fin solid inside an air box, terminals on the two end caps."""
    from netgen.occ import Box, Glue, OCCGeometry, Pnt

    import ngsolve as ng

    shape = OCCGeometry(str(step_path)).shape.Scale(Pnt(0, 0, 0), 1e-3)
    box = shape.bounding_box
    span = max(float(box[1][i]) - float(box[0][i]) for i in range(3))
    pad = box_factor * span

    fin = shape
    fin.mat("conductor")
    fin.maxh = maxh_conductor_m
    # The bounding box carries a small pad of its own -- 1e-7 m here -- so it
    # cannot be compared against face centres at any sane tolerance.  The end
    # caps are simply the faces whose centre sits at the extremes in z, with
    # every lateral face at mid-span, so take the extremes from the faces.
    centres_z = [float(face.center[2]) for face in fin.faces]
    z_lo, z_hi = min(centres_z), max(centres_z)
    tol = 1e-6 * (z_hi - z_lo)
    for face in fin.faces:
        z = float(face.center[2])
        face.name = ("source" if abs(z - z_lo) < tol
                     else "sink" if abs(z - z_hi) < tol
                     else "fin_lateral")

    x_lo, x_hi = float(box[0][0]) - pad, float(box[1][0]) + pad
    y_lo, y_hi = float(box[0][1]) - pad, float(box[1][1]) + pad
    air_box = Box(Pnt(x_lo, y_lo, z_lo), Pnt(x_hi, y_hi, z_hi))
    air = air_box - fin
    air.mat("air")
    air.maxh = maxh_air_m
    for face in air.faces:
        centre = face.center
        on_end = (abs(float(centre[2]) - z_lo) < tol
                  or abs(float(centre[2]) - z_hi) < tol)
        on_side = (abs(float(centre[0]) - x_lo) < tol
                   or abs(float(centre[0]) - x_hi) < tol
                   or abs(float(centre[1]) - y_lo) < tol
                   or abs(float(centre[1]) - y_hi) < tol)
        face.name = "outer" if (on_end or on_side) else "fin_lateral"

    geo = OCCGeometry(Glue([air, fin]))
    mesh = ng.Mesh(geo.GenerateMesh(maxh=maxh_air_m))
    if curve_order > 1:
        mesh.Curve(int(curve_order))
    return mesh, (z_lo, z_hi)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", type=Path, default=(
        ROOT / "tests/coil_from_cad/fixtures/beak_fin_48mm.step"))
    parser.add_argument("--frequency-hz", type=float, nargs="+",
                        default=[15_000.0, 150_000.0])
    parser.add_argument("--sigma", type=float, default=5.8e7)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--curve-order", type=int, default=2)
    parser.add_argument("--box-factor", type=float, default=0.35)
    parser.add_argument("--maxh-conductor-mm", type=float, default=0.6)
    parser.add_argument("--maxh-air-mm", type=float, default=4.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    import ngsolve as ng

    from radia.eddy_aphi import solve_eddy_aphi

    mesh, (z_lo, z_hi) = build_fin_mesh(
        args.step, box_factor=args.box_factor,
        maxh_conductor_m=args.maxh_conductor_mm * 1e-3,
        maxh_air_m=args.maxh_air_mm * 1e-3, curve_order=args.curve_order)
    region = mesh.Materials("conductor")
    volume = float(ng.Integrate(ng.CF(1.0), mesh, definedon=region).real)
    length = z_hi - z_lo
    print(json.dumps({"phase": "mesh", "ne": mesh.ne, "nv": mesh.nv,
                      "conductor_volume_mm3": volume * 1e9,
                      "length_mm": length * 1e3,
                      "materials": list(mesh.GetMaterials()),
                      "boundaries": sorted(set(mesh.GetBoundaries()))}),
          flush=True)

    rows = []
    for frequency in args.frequency_hz:
        delta = math.sqrt(2.0 / (2.0 * math.pi * frequency * MU0 * args.sigma))
        t0 = time.perf_counter()
        result = solve_eddy_aphi(mesh, conductor="conductor", source="source",
                                 sink="sink", frequency_hz=frequency,
                                 sigma=args.sigma, order=args.order,
                                 print_dofs=True)
        seconds = time.perf_counter() - t0
        # The loss split that the delivery gate compares, measured in the
        # volume rather than on a surface band.
        j_field = result.current_density()
        loss_density = (j_field * ng.Conj(j_field)).real / (2.0 * args.sigma)
        mid = ng.IfPos((ng.z - (z_lo + 0.35 * length))
                       * ((z_lo + 0.65 * length) - ng.z), 1.0, 0.0)
        total_mid = float(ng.Integrate(loss_density * mid, mesh,
                                       definedon=region).real)
        beak = ng.IfPos(ng.x - 1.6e-3, 1.0, 0.0)
        tip = ng.IfPos(ng.x - 3.5e-3, 1.0, 0.0)
        beak_loss = float(ng.Integrate(loss_density * mid * beak, mesh,
                                       definedon=region).real)
        tip_loss = float(ng.Integrate(loss_density * mid * tip, mesh,
                                      definedon=region).real)
        row = {
            "frequency_hz": frequency,
            "skin_depth_mm": delta * 1e3,
            "seconds": seconds,
            "ndof": result.ndof,
            "R_uohm": result.resistance_ohm * 1e6,
            "L_nH": result.inductance_H * 1e9,
            "mid_span_loss_W": total_mid,
            "beak_loss_fraction": beak_loss / total_mid if total_mid else 0.0,
            "tip_loss_fraction": tip_loss / total_mid if total_mid else 0.0,
        }
        rows.append(row)
        print(json.dumps({"phase": "level", "frequency_hz": frequency,
                          "delta_mm": round(delta * 1e3, 4),
                          "R_uohm": round(row["R_uohm"], 4),
                          "beak_loss_fraction":
                              round(row["beak_loss_fraction"], 5),
                          "seconds": round(seconds, 1)}), flush=True)

    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    report = {
        "schema": "radia.beak_fin_aphi_reference.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": head,
        "host": platform.node(),
        "question": ("BEM-A and the surface PEEC agree with each other, and "
                     "both assume a Leontovich surface impedance.  What does "
                     "an interior-resolved solve say?"),
        "fixture": str(args.step),
        "sigma_S_per_m": args.sigma,
        "mesh": {"n_elements": mesh.ne, "order": args.order,
                 "curve_order": args.curve_order,
                 "maxh_conductor_mm": args.maxh_conductor_mm,
                 "maxh_air_mm": args.maxh_air_mm,
                 "box_factor": args.box_factor,
                 "conductor_volume_mm3": volume * 1e9},
        "reference_validation": ("the solver reproduces the exact Bessel "
                                 "round-wire resistance to 2.3e-5 over "
                                 "a/delta from 0.25 to 4; see "
                                 "results/eddy_aphi_round_wire_20260922.json"),
        "levels": rows,
        "not_claimed": ("the mesh here resolves the skin only where the row "
                        "says so.  A level whose skin depth approaches the "
                        "element size is reported for cost, not for accuracy, "
                        "and the fraction comparison against BEM-A and PEEC is "
                        "only meaningful on the resolved levels"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps({"phase": "complete", "output": str(args.output)}))


if __name__ == "__main__":
    main()
