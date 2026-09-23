"""Gate 3's remainder: does the surface impedance fail on delta/R alone?

The 150 kHz comparison said the surface-impedance description of this section
loses 17% of its dissipation and redistributes the rest.  The explanation
offered was geometric -- a Leontovich impedance is the leading term of an
expansion in the skin depth over the local radius of curvature, and at the
beak tip that ratio is 0.68 -- but an explanation that only ever saw one
frequency and one conductivity is a story, not a measurement.

So the skin depth is moved two independent ways.  Frequency and conductivity
enter the physics separately -- one through ``omega`` in the diffusion term,
the other through ``sigma`` in both the diffusion term and the impedance
itself -- and if the explanation is right the error must not care which was
turned, only about the ratio that results.  A conductivity sweep that
reproduced the frequency sweep's error at the same ``delta/R`` would be the
evidence; one that did not would mean the diagnosis is wrong.

The geometry sweep is the control on that control.  ``tip_x`` changes how far
the beak reaches while the tip circle keeps its radius, so it moves the fin's
length and aspect without moving the curvature the diagnosis rests on.  The
error should then move only weakly, and certainly not the way ``delta/R``
moves it.

Everything is scored against the planar interior-resolved control of
``radia.sibc_corner_patch``, which reproduces the exact Bessel round wire to
3.5e-8 and is the reference for this section.
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
SCRATCH = Path(r"C:\temp\beak_validity_map")


def fixture_for(tip_x_mm, length_mm=48.0):
    """A straight prism with the given beak reach, written to scratch.

    Only ``tip_x = 4.0`` is a tracked fixture; the others are generated inputs
    for this sweep and live in scratch, with their parameters recorded in the
    result rather than the file being kept.
    """
    from build123d import export_step

    from make_beak_fin_step import make_beak_fin

    SCRATCH.mkdir(parents=True, exist_ok=True)
    path = SCRATCH / f"beak_fin_length{length_mm:.17g}_tip{tip_x_mm:.17g}mm.step"
    export_step(make_beak_fin(length_mm, tip_x_mm), str(path))
    return path


def section_metres(tip_x_mm):
    """The same section as an OCC face in metres, for the control."""
    from build123d import Plane, export_step
    from netgen.occ import OCCGeometry, Pnt

    from make_beak_fin_step import _beak_face

    SCRATCH.mkdir(parents=True, exist_ok=True)
    tmp = SCRATCH / f"_section_tip{tip_x_mm:.2f}mm.step"
    export_step(_beak_face(Plane.XY, tip_x=tip_x_mm), str(tmp))
    shape = OCCGeometry(str(tmp), dim=2).shape.Scale(Pnt(0, 0, 0), 1e-3)
    tmp.unlink(missing_ok=True)
    return shape


def one_point(tip_x_mm, frequency, sigma, *, n_peri, maxh_over_delta, order,
              air_pad_mm):
    """Surface impedance against the interior-resolved control, at one point."""
    from make_beak_fin_step import TIP_RADIUS
    from radia.sibc_corner_patch import solve_cut_patch
    from sibc2d_beak_reference import solve_beak_sibc_2d

    t0 = time.perf_counter()
    outer = solve_beak_sibc_2d(fixture_for(tip_x_mm), frequency=frequency,
                               sigma=sigma, n_peri=n_peri)
    delta = outer["skin_depth_m"]
    xy = np.asarray(outer["panel_xy_m"], dtype=float)
    ds = np.asarray(outer["panel_ds_m"], dtype=float)
    current = np.asarray(outer["panel_current_A"])
    rs = complex(outer["surface_impedance_ohm"]).real
    sibc_total = float(np.sum(0.5 * rs * np.abs(current) ** 2 / ds))

    # The cuts follow the fin the section analysis actually found, not the
    # fixture constants: a shorter beak puts its tip region somewhere else,
    # and a fixed 3.5 mm cut would fall past the end of it.
    fin = outer["fin"]
    beak_cut = float(fin["root_mid_m"][0])
    tip_cut = float(fin["extremity_m"][0]) - float(fin["tip_span_m"])

    control = solve_cut_patch(
        section_metres(tip_x_mm), x_cut=None, panel_xy=xy, panel_ds=ds,
        panel_current=current, axial_E_field=outer["axial_E_field_V_per_m"],
        surface_impedance=outer["surface_impedance_ohm"],
        omega=outer["omega_rad_per_s"], sigma=sigma,
        maxh_conductor=maxh_over_delta * delta, maxh_air=4.0e-3,
        air_pad=air_pad_mm * 1e-3, order=order)
    scale = 1.0 / abs(control.total_current()) ** 2
    total = control.total_loss() * scale
    beak = control.loss_beyond(beak_cut) * scale / total
    tip = control.loss_beyond(tip_cut) * scale / total
    if min(beak, tip) <= 0.0:
        raise ValueError(
            f"a cut fell outside the conductor: beak at {beak_cut:.4e} m, "
            f"tip at {tip_cut:.4e} m, for a fin reaching "
            f"{float(fin['extremity_m'][0]):.4e} m")

    return {
        "tip_x_mm": tip_x_mm,
        "frequency_hz": frequency,
        "sigma_S_per_m": sigma,
        "skin_depth_mm": delta * 1e3,
        "skin_depth_over_tip_radius": delta * 1e3 / float(TIP_RADIUS),
        "cuts_mm": {"beak": beak_cut * 1e3, "tip": tip_cut * 1e3},
        "fin_length_mm": (float(fin["extremity_m"][0]) - beak_cut) * 1e3,
        "ndof": control.ndof,
        "seconds": time.perf_counter() - t0,
        "exact": {"total_loss_W_per_m": total,
                  "beak_loss_fraction": beak, "tip_loss_fraction": tip},
        "surface_impedance": {
            "total_loss_relative": sibc_total / total - 1.0,
            "beak_relative": outer["beak_loss_fraction"] / beak - 1.0,
            "tip_relative": outer["tip_loss_fraction"] / tip - 1.0,
        },
    }


def _collapse(rows):
    """Do the frequency and conductivity families agree at equal delta/R?

    Pairs are matched on the skin depth itself, which is the quantity the
    diagnosis names.  A pair that agrees says the two ways of moving it are
    interchangeable; a pair that does not says something other than delta/R is
    driving the error.
    """
    by_family = {}
    for row in rows:
        by_family.setdefault(row["family"], []).append(row)
    frequency = by_family.get("frequency", [])
    conductivity = by_family.get("conductivity", [])
    pairs = []
    for a in conductivity:
        match = min(frequency,
                    key=lambda b: abs(b["skin_depth_mm"] / a["skin_depth_mm"]
                                      - 1.0), default=None)
        if match is None:
            continue
        mismatch = abs(match["skin_depth_mm"] / a["skin_depth_mm"] - 1.0)
        if mismatch > 0.02:
            continue
        pairs.append({
            "skin_depth_mm": a["skin_depth_mm"],
            "skin_depth_mismatch": mismatch,
            "by_frequency": {
                "frequency_hz": match["frequency_hz"],
                "sigma_S_per_m": match["sigma_S_per_m"],
                "total_loss_relative":
                    match["surface_impedance"]["total_loss_relative"],
                "beak_relative": match["surface_impedance"]["beak_relative"],
            },
            "by_conductivity": {
                "frequency_hz": a["frequency_hz"],
                "sigma_S_per_m": a["sigma_S_per_m"],
                "total_loss_relative":
                    a["surface_impedance"]["total_loss_relative"],
                "beak_relative": a["surface_impedance"]["beak_relative"],
            },
            "total_loss_error_difference": abs(
                a["surface_impedance"]["total_loss_relative"]
                - match["surface_impedance"]["total_loss_relative"]),
            "beak_error_difference": abs(
                a["surface_impedance"]["beak_relative"]
                - match["surface_impedance"]["beak_relative"]),
        })
    return pairs


def _decay_exponent(rows, key):
    """Fit ``|error| ~ delta^p`` over a family, and say what ``p`` means.

    The exponent is the diagnosis.  A smooth surface whose curvature is merely
    under-resolved gives a first curvature correction of order ``delta / R``,
    so ``p = 1``.  A corner has no radius to divide by; Dauge, Dular,
    Krahenbuhl, Peron, Perrussel and Poignard (2014) put its layer at the
    scale of ``delta`` and its contribution at half-integer order, so a corner
    shows ``p = 1/2``.  Measuring ``p`` therefore says which feature owns the
    error, rather than leaving it to be argued.
    """
    delta = np.array([r["skin_depth_mm"] for r in rows])
    error = np.array([abs(r["surface_impedance"][key]) for r in rows])
    take = error > 0
    if take.sum() < 3:
        return None
    slope, intercept = np.polyfit(np.log(delta[take]), np.log(error[take]), 1)
    fitted = np.exp(intercept) * delta[take] ** slope
    residual = float(np.max(np.abs(fitted / error[take] - 1.0)))
    return {"exponent": float(slope), "worst_relative_residual": residual,
            "reads_as": ("one half is a corner, where there is no radius to "
                         "divide by; one would be a resolved smooth surface "
                         "whose curvature correction was simply dropped")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-peri", type=int, default=1024)
    parser.add_argument("--maxh-over-delta", type=float, default=0.25)
    parser.add_argument("--order", type=int, default=3)
    parser.add_argument("--air-pad-mm", type=float, default=24.0)
    parser.add_argument("--sigma", type=float, default=5.8e7)
    parser.add_argument("--frequency-hz", type=float, nargs="+",
                        default=[15_000.0, 50_000.0, 150_000.0, 450_000.0])
    parser.add_argument("--conductivity-frequency-hz", type=float,
                        default=150_000.0)
    parser.add_argument("--sigma-scale", type=float, nargs="+",
                        default=[10.0, 3.0, 1.0, 1.0 / 3.0])
    parser.add_argument("--tip-x-mm", type=float, nargs="+",
                        default=[4.0, 3.2, 2.4])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    common = dict(n_peri=args.n_peri, maxh_over_delta=args.maxh_over_delta,
                  order=args.order, air_pad_mm=args.air_pad_mm)
    rows = []

    for frequency in args.frequency_hz:
        row = one_point(4.0, frequency, args.sigma, **common)
        row["family"] = "frequency"
        rows.append(row)
        print(json.dumps({"phase": "frequency", **row["surface_impedance"],
                          "skin_depth_mm": row["skin_depth_mm"]}), flush=True)

    for scale in args.sigma_scale:
        row = one_point(4.0, args.conductivity_frequency_hz,
                        args.sigma * scale, **common)
        row["family"] = "conductivity"
        row["sigma_scale"] = scale
        rows.append(row)
        print(json.dumps({"phase": "conductivity", **row["surface_impedance"],
                          "skin_depth_mm": row["skin_depth_mm"]}), flush=True)

    for tip_x_mm in args.tip_x_mm:
        row = one_point(tip_x_mm, 150_000.0, args.sigma, **common)
        row["family"] = "geometry"
        rows.append(row)
        print(json.dumps({"phase": "geometry", "tip_x_mm": tip_x_mm,
                          **row["surface_impedance"]}), flush=True)

    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    geometry = [r for r in rows if r["family"] == "geometry"]
    report = {
        "schema": "radia.beak_sibc_validity_map.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": head,
        "host": platform.node(),
        "question": ("the surface impedance's failure on this section was "
                     "explained by delta over the local radius of curvature.  "
                     "Does the error actually depend on that ratio alone?"),
        "reference": ("the planar interior-resolved control of "
                      "radia.sibc_corner_patch, which reproduces the exact "
                      "Bessel round wire to 3.5e-8"),
        "levels": rows,
        "frequency_against_conductivity": _collapse(rows),
        "how_the_error_decays": {
            "total_loss": _decay_exponent(
                [r for r in rows if r["family"] in ("frequency",
                                                    "conductivity")],
                "total_loss_relative"),
            "over_skin_depths_mm": sorted(
                r["skin_depth_mm"] for r in rows
                if r["family"] in ("frequency", "conductivity")),
        },
        "geometry_spread": {
            "tip_x_mm": [r["tip_x_mm"] for r in geometry],
            "total_loss_relative": [
                r["surface_impedance"]["total_loss_relative"]
                for r in geometry],
            "beak_relative": [r["surface_impedance"]["beak_relative"]
                              for r in geometry],
        },
        "not_claimed": (
            "the geometry sweep moves the beak's reach at a fixed tip radius; "
            "it does not vary the radius itself, so it separates aspect from "
            "curvature rather than mapping curvature.  A shorter beak runs a "
            "steeper flank into the same arc, so the section analysis fits "
            "the tip circle 1.5% large at a 2.4 mm reach and the tip cut, "
            "which is two fitted radii back from the extremity, is that much "
            "wider; the arc itself is exactly TIP_RADIUS by construction at "
            "every reach.  All of it is one section in two dimensions, "
            "mid-span"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps({"phase": "complete", "output": str(args.output)}))


if __name__ == "__main__":
    main()
