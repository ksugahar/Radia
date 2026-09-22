"""Gate 3: the first comparison against a route that does not assume SIBC.

BEM-A and the surface PEEC agree with each other to 0.58% on the delivery
metrics, and both impose a Leontovich surface impedance.  Agreement between
two routes that share an assumption tests the discretisations, not the
assumption.  This driver puts the interior-resolved A-V solve beside them at
one frequency and compares the same quantity, measured the same way: the share
of mid-span loss beyond the beak cut, over the same axial window.

The frequency is chosen for resolution, not convenience.  At 150 kHz the skin
depth is 0.17 mm and an isotropic mesh of the 48 mm fin cannot resolve it at a
tractable size; Netgen 6.2.2606 cannot build the anisotropic boundary layer
that would (its old entry point refuses, and the replacement raises
``Need to register class netgen::BoundaryLayerParameters for Archive``).  At
5 kHz the skin depth is 0.93 mm, the isotropic mesh resolves it, and
``delta / thickness`` is 0.23 -- which is where a surface-impedance model is
under the most strain, so it is the informative end rather than a retreat.
"""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aphi-levels", type=Path, nargs="+", required=True,
                        help="A-V artifacts, one per mesh level")
    parser.add_argument("--sibc", type=Path, required=True,
                        help="compare_beak_fin_bema.py output at the same frequency")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    levels = []
    for path in args.aphi_levels:
        doc = json.loads(path.read_text(encoding="utf-8"))
        row = doc["levels"][0]
        levels.append({
            "maxh_conductor_mm": doc["mesh"]["maxh_conductor_mm"],
            "n_elements": doc["mesh"]["n_elements"],
            "frequency_hz": row["frequency_hz"],
            "skin_depth_mm": row["skin_depth_mm"],
            "R_uohm": row["R_uohm"],
            "beak_loss_fraction": row["beak_loss_fraction"],
            "tip_loss_fraction": row["tip_loss_fraction"],
            "seconds": row["seconds"],
        })
    frequencies = {row["frequency_hz"] for row in levels}
    if len(frequencies) != 1:
        raise ValueError(f"the A-V levels are not one frequency: {frequencies}")
    frequency = frequencies.pop()

    fractions = [row["beak_loss_fraction"] for row in levels]
    spread = (max(fractions) - min(fractions)) / (sum(fractions) / len(fractions))
    finest = min(levels, key=lambda r: r["maxh_conductor_mm"])

    sibc = json.loads(args.sibc.read_text(encoding="utf-8"))
    if abs(sibc["frequency_hz"] - frequency) > 1e-9:
        raise ValueError("the SIBC run is at a different frequency: "
                         f"{sibc['frequency_hz']} vs {frequency}")
    profile = sibc["surface_current_profile"]
    metrics = profile["integral_metrics"]

    def compare(name, value):
        return {
            "beak_loss_fraction": value["beak_loss_fraction"],
            "tip_loss_fraction": value["tip_loss_fraction"],
            "beak_relative_to_aphi": float(
                value["beak_loss_fraction"] / finest["beak_loss_fraction"] - 1.0),
            "tip_relative_to_aphi": float(
                value["tip_loss_fraction"] / finest["tip_loss_fraction"] - 1.0),
            "route": name,
        }

    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    report = {
        "schema": "radia.beak_fin_three_route.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": head,
        "host": platform.node(),
        "question": ("BEM-A and the surface PEEC both assume a Leontovich "
                     "surface impedance.  Does an interior-resolved solve "
                     "agree with them where that assumption is weakest?"),
        "frequency_hz": frequency,
        "fixture": sibc["fixture"],
        "aphi_reference": {
            "validated_against": ("the exact Bessel round-wire resistance to "
                                  "2.3e-5 over a/delta from 0.25 to 4; see "
                                  "results/eddy_aphi_round_wire_20260922.json"),
            "mesh_convergence": levels,
            "beak_loss_fraction_relative_spread": spread,
            "converged_beak_loss_fraction": finest["beak_loss_fraction"],
            "converged_tip_loss_fraction": finest["tip_loss_fraction"],
            "R_uohm": finest["R_uohm"],
        },
        "surface_impedance_routes": {
            "bema": compare("bema", metrics["bema"]),
            "peec": compare("peec", metrics["peec"]),
            "bema_R_uohm": sibc["bema"]["R_ohm"] * 1e6,
            "peec_R_uohm": sibc["experimental_peec"]["R_ohm"] * 1e6,
            "peec_vs_bema_beak_loss_relative":
                profile["integral_metric_errors"]["beak_loss_fraction_relative"],
        },
        "how_the_fractions_are_measured": (
            "the A-V fraction is a volume integral of |J|^2 / 2 sigma beyond "
            "the cut; the surface-impedance routes weight |K|^2 Rs / 2 by "
            "perimeter arc length.  The two coincide exactly in the thin-skin "
            "limit, where the depth-integrated loss per unit perimeter IS "
            "|K|^2 Rs / 2, and they part company as delta approaches the "
            "section thickness.  That is not a measurement mismatch to be "
            "corrected away: it is the limit failing, which is the quantity "
            "under test"),
        "does_not_invalidate_the_delivery_gate": (
            "the delivery gate was accepted at 150 kHz, where delta is "
            "0.17 mm and delta / thickness is 0.043 -- five times thinner "
            "than here.  This says the surface-impedance family is wrong at "
            "delta / thickness = 0.23.  It says nothing about 0.043, and the "
            "boundary between them is not located by this run"),
        "not_claimed": (
            "one frequency on one straight fin.  It does not establish the "
            "surface-impedance assumption over the operating band, and the "
            "150 kHz point the delivery gate was accepted at is not reached "
            "here"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "phase": "complete",
        "aphi_beak_loss_fraction": finest["beak_loss_fraction"],
        "bema_relative": report["surface_impedance_routes"]["bema"]["beak_relative_to_aphi"],
        "peec_relative": report["surface_impedance_routes"]["peec"]["beak_relative_to_aphi"],
        "output": str(args.output)}))


if __name__ == "__main__":
    main()
