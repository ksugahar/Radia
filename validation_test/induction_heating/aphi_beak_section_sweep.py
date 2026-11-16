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
                       exterior_factor, air_divisor=20,
                       skin_layers_mm=None):
    """The beak profile as the meridian of a ring at ``radius_mm``.

    Returns ``(mesh, exterior_radius_m, exterior_center_z_m)``: the air is a
    half disk about the axis, ``exterior_factor`` times the section's
    farthest reach, closed by the exact exterior term.  The multipole series
    converges as ``1 / exterior_factor`` per mode, so the factor, not a fixed
    margin, keeps the mode count bounded at every ring radius.  Far air is
    meshed no finer than ``1 / air_divisor`` of that radius (default a
    twentieth); a linear-order reference needs a larger divisor.
    ``skin_layers_mm`` inserts that stack of layers inward from the conductor
    surface, so the normal resolution can be refined while the surface nodes,
    the core and the air stay as they are.
    """
    from build123d import Plane
    from make_beak_fin_step import _beak_face
    from netgen.occ import Circle, Glue, MoveTo, OCCGeometry, X

    import ngsolve as ng

    # Reuse the fixture's own profile so the section is the tracked one.
    face = _beak_face(Plane.XY, radial_offset=radius_mm)
    import tempfile

    from build123d import export_step

    # A private directory per call: concurrent jobs must not share (and
    # delete) one STEP file, and the checkout must stay clean during a run.
    with tempfile.TemporaryDirectory(prefix="beak_section_") as tmp:
        step = Path(tmp) / "section.step"
        export_step(face, str(step))
        section = OCCGeometry(str(step), dim=2).shape.Scale(
            __import__("netgen.occ", fromlist=["Pnt"]).Pnt(0, 0, 0), 1e-3)

    section.faces.name = "conductor"
    section.edges.name = "skin"
    section.maxh = maxh_conductor_mm * 1e-3
    box = section.bounding_box
    center_z = 0.5 * (float(box[0][1]) + float(box[1][1]))
    reach = max(math.hypot(float(box[i][0]), float(box[j][1]) - center_z)
                for i in range(2) for j in range(2))
    radius = exterior_factor * reach
    air_h = max(maxh_air_mm * 1e-3, radius / float(air_divisor))
    disk = Circle((0.0, center_z), radius).Face()
    half = MoveTo(0.0, center_z - radius).Rectangle(radius, 2 * radius).Face()
    air = disk * half
    air.edges.name = "outer"
    # The mesh is straight (the Henrotte space refuses curved elements), so
    # the exterior arc is a chord polygon: R / 200 keeps its miss inside the
    # DtN term's 1e-5 geometry tolerance.
    air.edges.maxh = min(air_h, radius / 200.0)
    air.edges.Min(X).name = "axis"
    air.maxh = air_h
    outside = air - section
    outside.faces.name = "air"
    total = Glue([outside, section])
    layers = []
    if skin_layers_mm:
        from netgen.meshing import BoundaryLayerParameters

        layers = [BoundaryLayerParameters(
            boundary="skin", thickness=[t * 1e-3 for t in skin_layers_mm],
            domain="conductor", disable_curving=False)]
    mesh = ng.Mesh(OCCGeometry(total, dim=2).GenerateMesh(
        maxh=air_h, boundary_layers=layers))
    return mesh, radius, center_z


def shares(result, radius_mm):
    """Beak and tip loss shares, by the same cuts the other routes use.

    Dissipation lives in the revolved solid, so it is integrated in the volume
    measure; ``AxisymRingResult`` owns that choice rather than each caller
    rediscovering it.  The weight is a fraction of a percent at these radii,
    but it is a bias, not noise, and it is smallest exactly where the model
    looks most trustworthy.
    """
    total = result.total_loss()
    beak = result.loss_beyond((radius_mm + BEAK_X_MM) * 1e-3)
    tip = result.loss_beyond((radius_mm + TIP_X_MM) * 1e-3)
    return total, beak / total, tip / total


def _compare(rows, planar_path):
    """Put this route beside the planar control, and say what still separates.

    The planar control resolves the same section and imposes the exterior
    field of the same total current on its air boundary, so its truncation
    states the answer instead of obstructing it, and it is anchored on the
    exact Bessel round wire.  This route closes its air with the exact
    exterior term, so the exterior-factor sweep in ``levels`` should show no
    dependence on where the air is cut.
    """
    from make_beak_fin_step import TIP_RADIUS

    at_150k = [r for r in rows if abs(r["frequency_hz"] - 150_000.0) < 1e-6]
    if not at_150k:
        return {"note": "no 150 kHz level in this run"}
    # Largest radius, finest mesh, widest air: the least curved corner.
    finest = max(at_150k, key=lambda r: (r["radius_mm"],
                                         -r["maxh_conductor_mm"],
                                         r["exterior_factor"]))
    out = {
        "frequency_hz": 150_000.0,
        "skin_depth_mm": finest["skin_depth_mm"],
        "tip_radius_mm": float(TIP_RADIUS),
        "tip_radius_over_skin_depth": float(
            TIP_RADIUS / finest["skin_depth_mm"]),
        "why_the_surface_impedance_fails_here": (
            "a Leontovich impedance is the leading term of an expansion in "
            "the ratio of the skin depth to the radius of curvature of the "
            "surface.  At the beak tip that ratio is 0.68 -- the tip radius "
            "is under one and a half skin depths -- so the first curvature "
            "correction is of order one rather than of order a percent, and "
            "the beak root and the far end carry square corners, where the "
            "radius of curvature is zero and no term of the expansion "
            "applies at all"),
        "least_curved_axisymmetric": {
            "beak_loss_fraction": finest["beak_loss_fraction"],
            "tip_loss_fraction": finest["tip_loss_fraction"],
            "radius_mm": finest["radius_mm"],
            "maxh_conductor_mm": finest["maxh_conductor_mm"],
            "exterior_factor": finest["exterior_factor"],
        },
    }
    if not planar_path.exists():
        out["planar_control"] = {"note": f"{planar_path.name} not present"}
        return out
    planar = json.loads(planar_path.read_text(encoding="utf-8"))
    control = planar["control"]["levels"][-1]
    out["planar_control"] = {
        "source": planar_path.name,
        "beak_loss_fraction": control["beak_loss_fraction"],
        "tip_loss_fraction": control["tip_loss_fraction"],
        "worst_error_against_exact_bessel": max(
            abs(row["relative_error"])
            for row in planar["control"]["round_wire"]),
    }
    out["axisymmetric_relative_to_planar"] = {
        "beak_loss_fraction": float(finest["beak_loss_fraction"]
                                    / control["beak_loss_fraction"] - 1.0),
        "tip_loss_fraction": float(finest["tip_loss_fraction"]
                                   / control["tip_loss_fraction"] - 1.0),
        "what_it_is": (
            "the difference between two interior-resolved routes with "
            "independent open boundaries; what remains is the ring's "
            "curvature, measured by the radius sweep"),
    }
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius-mm", type=float, nargs="+",
                        default=[200.0, 400.0, 1600.0, 6400.0, 25600.0])
    parser.add_argument("--frequency-hz", type=float, nargs="+",
                        default=[5_000.0, 15_000.0, 50_000.0, 150_000.0])
    parser.add_argument("--sigma", type=float, default=5.8e7)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--maxh-conductor-mm", type=float, nargs="+",
                        default=[0.08, 0.05])
    parser.add_argument("--maxh-air-mm", type=float, default=20.0)
    parser.add_argument("--exterior-factor", type=float, nargs="+",
                        default=[1.5, 3.0],
                        help="exterior radius over the section's farthest reach "
                             "from the axis centre; the exterior term is "
                             "exact, so results must not depend on it")
    parser.add_argument("--exterior-modes", type=int, default=30)
    parser.add_argument("--planar-control", type=Path, default=(
        HERE / "results" / "sibc_patch_corrected_beak_20260922.json"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    from radia.eddy_axisym_ring import solve_axisym_ring

    rows = []
    for factor in args.exterior_factor:
        for radius_mm in args.radius_mm:
            for maxh in args.maxh_conductor_mm:
                mesh, ext_radius, ext_z = build_section_ring(
                    radius_mm, maxh_conductor_mm=maxh,
                    maxh_air_mm=args.maxh_air_mm, exterior_factor=factor)
                for frequency in args.frequency_hz:
                    t0 = time.perf_counter()
                    res = solve_axisym_ring(
                        mesh, conductor="conductor", frequency_hz=frequency,
                        sigma=args.sigma, order=args.order,
                        exterior_radius=ext_radius,
                        exterior_center_z=ext_z,
                        exterior_modes=args.exterior_modes)
                    seconds = time.perf_counter() - t0
                    total, beak, tip = shares(res, radius_mm)
                    delta = res.skin_depth_m
                    row = {
                        "exterior_factor": factor,
                        "exterior_modes": args.exterior_modes,
                        "radius_mm": radius_mm,
                        "maxh_conductor_mm": maxh,
                        "frequency_hz": frequency,
                        "skin_depth_mm": delta * 1e3,
                        "elements_per_skin_depth": delta / (maxh * 1e-3),
                        "ndof": res.ndof,
                        "seconds": seconds,
                        "R_uohm_per_m": res.impedance_per_metre.real * 1e6,
                        "total_loss_W": total,
                        "beak_loss_fraction": beak,
                        "tip_loss_fraction": tip,
                    }
                    rows.append(row)
                    print(json.dumps({"phase": "level", **{
                        k: row[k] for k in (
                            "exterior_factor", "radius_mm",
                            "maxh_conductor_mm",
                            "frequency_hz", "beak_loss_fraction",
                            "tip_loss_fraction")}}), flush=True)

    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    report = {
        "schema": "radia.beak_section_axisym_sweep.v2",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": head,
        "platform_class": platform.system(),
        "question": ("the delivery gate was accepted at 150 kHz on two routes "
                     "that share a surface-impedance assumption.  What does "
                     "the same section say when its interior is resolved?"),
        "section": "the tracked beak profile, revolved about the axis",
        "sigma_S_per_m": args.sigma,
        "reference_validation": (
            "the A_phi drive and exterior term reproduce the exact Bessel "
            "resistance of a revolved round wire to 1.5e-4 at a/R0 = 1/300 "
            "and 100 kHz, and the direct-current inductance of a "
            "rectangular ring to 4.4e-4 of an independent filament "
            "quadrature over a sixfold range of exterior radius "
            "(axisym_ring_exterior_validation.py)"),
        "exterior": ("exact Dirichlet-to-Neumann term on a half disk about "
                     "the axis, Woodbury-applied"),
        "measure": (
            "dissipation is integrated in the revolved volume measure "
            "2 pi r dr dz, through radia.axisym_measure; the terminal current "
            "is a meridian-measure quantity and is integrated as one"),
        "cuts_mm": {"beak": BEAK_X_MM, "tip": TIP_X_MM},
        "levels": rows,
        "against_the_planar_control": _compare(rows, args.planar_control),
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
    from ngsolve import TaskManager
    with TaskManager():
        main()
