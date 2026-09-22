"""Gate 3 at the frequency that matters: the beak section, skin resolved.

The delivery gate was accepted at 150 kHz on the strength of BEM-A and the
surface PEEC agreeing.  Both impose a Leontovich surface impedance, so the
agreement never tested it.  The three-dimensional A-V reference showed the
assumption failing at 5 kHz but could not reach 150 kHz: resolving a 0.17 mm
skin over a 48 mm fin isotropically needs about two million elements.

Revolving the section removes the third dimension without removing the
physics.  The beak profile becomes the meridian of a ring, the skin is
resolved by a two-dimensional mesh at a few tens of thousands of degrees of
freedom, and a solve takes under a second.  The curvature this introduces is
made small by placing the section at a large radius, and is measured by
varying that radius rather than assumed away.

The counterpart is `sibc2d_beak_reference.py`, which solves the SAME section
with a surface impedance.  Same geometry, same frequency, same definition of
the beak and tip shares -- the only difference between the two numbers is the
assumption under test.
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

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))

MU0 = 4e-7 * math.pi
BEAK_X_MM = 1.6
TIP_X_MM = 3.5


def build_section_ring(radius_mm, *, maxh_conductor_mm, maxh_air_mm,
                       box_factor, curve_order):
    """The beak profile as the meridian of a ring at ``radius_mm``."""
    from build123d import Plane
    from make_beak_fin_step import _beak_face
    from netgen.occ import Glue, MoveTo, OCCGeometry, WorkPlane, X, Y

    import ngsolve as ng

    # Reuse the fixture's own profile so the section is the tracked one.
    face = _beak_face(Plane.XY, radial_offset=radius_mm)
    step = Path(ROOT / "validation_test" / "induction_heating"
                / "_section_tmp.step")
    from build123d import export_step

    export_step(face, str(step))
    section = OCCGeometry(str(step), dim=2).shape.Scale(
        __import__("netgen.occ", fromlist=["Pnt"]).Pnt(0, 0, 0), 1e-3)
    step.unlink(missing_ok=True)

    section.faces.name = "conductor"
    section.maxh = maxh_conductor_mm * 1e-3
    box = section.bounding_box
    half = box_factor * max(float(box[1][i]) - float(box[0][i])
                            for i in range(2))
    r_max = float(box[1][0]) + half
    z_half = half + 0.5 * (float(box[1][1]) - float(box[0][1]))
    air = MoveTo(0, -z_half).Rectangle(r_max, 2 * z_half).Face()
    air.faces.name = "air"
    air.maxh = maxh_air_mm * 1e-3
    outside = air - section
    outside.faces.name = "air"
    total = Glue([outside, section])
    total.edges.Min(X).name = "axis"
    total.edges.Max(X).name = "outer"
    total.edges.Min(Y).name = "outer"
    total.edges.Max(Y).name = "outer"
    mesh = ng.Mesh(OCCGeometry(total, dim=2).GenerateMesh(
        maxh=maxh_air_mm * 1e-3))
    if curve_order > 1:
        mesh.Curve(int(curve_order))
    return mesh


def shares(result, radius_mm):
    """Beak and tip loss shares, by the same cuts the other routes use."""
    import ngsolve as ng

    mesh = result.mesh
    region = mesh.Materials(result.conductor)
    j = result.current_density()
    loss = (j * ng.Conj(j)).real / (2.0 * result.sigma_S_per_m)
    beak_cut = (radius_mm + BEAK_X_MM) * 1e-3
    tip_cut = (radius_mm + TIP_X_MM) * 1e-3
    total = float(ng.Integrate(loss, mesh, definedon=region, order=10).real)
    beak = float(ng.Integrate(loss * ng.IfPos(ng.x - beak_cut, 1.0, 0.0), mesh,
                              definedon=region, order=10).real)
    tip = float(ng.Integrate(loss * ng.IfPos(ng.x - tip_cut, 1.0, 0.0), mesh,
                             definedon=region, order=10).real)
    return total, beak / total, tip / total


def _compare_to_sibc(rows):
    """Put the converged 150 kHz row beside the 2-D surface-impedance solve.

    The counterpart is the same section at the same frequency, converged to
    1024 perimeter samples, differing only in that it replaces the conductor's
    interior by a Leontovich impedance.
    """
    sibc_path = (HERE / "beak_fin_discretization_convergence_150kHz.json")
    sibc = json.loads(sibc_path.read_text(encoding="utf-8"))
    ref = sibc["reference_delivery_metrics_n1024"]
    at_150k = [r for r in rows if abs(r["frequency_hz"] - 150_000.0) < 1e-6]
    if not at_150k:
        return {"note": "no 150 kHz level in this run"}
    # Largest radius, finest mesh: the converged corner.
    finest = max(at_150k, key=lambda r: (r["radius_mm"],
                                         -r["maxh_conductor_mm"]))
    from make_beak_fin_step import TIP_RADIUS

    return {
        "frequency_hz": 150_000.0,
        "skin_depth_mm": finest["skin_depth_mm"],
        "tip_radius_mm": float(TIP_RADIUS),
        "tip_radius_over_skin_depth": float(
            TIP_RADIUS / finest["skin_depth_mm"]),
        "why_it_matters": (
            "a Leontovich impedance is the leading term of an expansion in "
            "the ratio of skin depth to the radius of curvature of the "
            "surface.  At the beak tip that ratio is not small -- the tip "
            "radius is under one and a half skin depths -- so the first "
            "curvature correction is of order one rather than of order a "
            "percent, and the beak root carries genuine sharp corners where "
            "the radius of curvature is zero and no term of the expansion "
            "applies at all.  This is where the assumption has to fail "
            "first, and it is also what a beak fin is for"),
        "axisymmetric_interior_resolved": {
            "beak_loss_fraction": finest["beak_loss_fraction"],
            "tip_loss_fraction": finest["tip_loss_fraction"],
            "radius_mm": finest["radius_mm"],
            "maxh_conductor_mm": finest["maxh_conductor_mm"],
        },
        "two_dimensional_sibc_n1024": {
            "beak_loss_fraction": ref["beak_loss_fraction"],
            "tip_loss_fraction": ref["tip_loss_fraction"],
        },
        "sibc_relative_to_resolved": {
            "beak_loss_fraction": float(
                ref["beak_loss_fraction"] / finest["beak_loss_fraction"] - 1.0),
            "tip_loss_fraction": float(
                ref["tip_loss_fraction"] / finest["tip_loss_fraction"] - 1.0),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius-mm", type=float, nargs="+",
                        default=[200.0, 400.0])
    parser.add_argument("--frequency-hz", type=float, nargs="+",
                        default=[5_000.0, 15_000.0, 50_000.0, 150_000.0])
    parser.add_argument("--sigma", type=float, default=5.8e7)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--curve-order", type=int, default=3)
    parser.add_argument("--maxh-conductor-mm", type=float, nargs="+",
                        default=[0.08, 0.05])
    parser.add_argument("--maxh-air-mm", type=float, default=20.0)
    parser.add_argument("--box-factor", type=float, default=3.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    from radia.eddy_axisym_ring import solve_axisym_ring

    rows = []
    for radius_mm in args.radius_mm:
        for maxh in args.maxh_conductor_mm:
            mesh = build_section_ring(
                radius_mm, maxh_conductor_mm=maxh,
                maxh_air_mm=args.maxh_air_mm, box_factor=args.box_factor,
                curve_order=args.curve_order)
            for frequency in args.frequency_hz:
                t0 = time.perf_counter()
                res = solve_axisym_ring(mesh, conductor="conductor",
                                        frequency_hz=frequency,
                                        sigma=args.sigma, order=args.order)
                seconds = time.perf_counter() - t0
                total, beak, tip = shares(res, radius_mm)
                delta = res.skin_depth_m
                row = {
                    "radius_mm": radius_mm,
                    "maxh_conductor_mm": maxh,
                    "frequency_hz": frequency,
                    "skin_depth_mm": delta * 1e3,
                    "elements_per_skin_depth": delta / (maxh * 1e-3),
                    "ndof": res.ndof,
                    "seconds": seconds,
                    "R_uohm_per_m": res.impedance_per_metre.real * 1e6,
                    "beak_loss_fraction": beak,
                    "tip_loss_fraction": tip,
                }
                rows.append(row)
                print(json.dumps({"phase": "level", **{
                    k: row[k] for k in ("radius_mm", "maxh_conductor_mm",
                                        "frequency_hz", "beak_loss_fraction",
                                        "tip_loss_fraction")}}), flush=True)

    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    report = {
        "schema": "radia.beak_section_axisym_sweep.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": head,
        "host": platform.node(),
        "question": ("the delivery gate was accepted at 150 kHz on two routes "
                     "that share a surface-impedance assumption.  What does "
                     "the same section say when its interior is resolved?"),
        "section": "the tracked beak profile, revolved about the axis",
        "sigma_S_per_m": args.sigma,
        "reference_validation": (
            "the solver reproduces the exact Bessel round-wire resistance to "
            "0.6% at a/delta 6 and 12, mesh-converged, with the residual "
            "shrinking as 1/R0 -- it is the ring's curvature, not the solve"),
        "cuts_mm": {"beak": BEAK_X_MM, "tip": TIP_X_MM},
        "levels": rows,
        "against_the_surface_impedance_reference": _compare_to_sibc(rows),
        "not_claimed": (
            "a ring is not a straight fin.  The curvature is made small by "
            "the radius and its size is measured by varying that radius; it "
            "is not eliminated.  End effects are absent by construction, "
            "which is what makes this the mid-span comparison and not a "
            "terminal one"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps({"phase": "complete", "output": str(args.output)}))


if __name__ == "__main__":
    main()
