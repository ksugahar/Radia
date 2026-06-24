"""Validation-class Carter air-gap and magnetizing-inductance sweep.

This is an example/validation run rather than a pytest test. It links the
standard slotted-air-gap correction to the AC-machine main-flux inductance:

* slot opening raises Carter's coefficient and the effective gap;
* mean permeance and magnetizing inductance fall by the reciprocal factor;
* synchronous dq magnetizing inductance is (m/2) times the per-phase value.

Run:

    python examples/electric_machine/validation_carter_magnetizing_sweep.py
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.solve import (  # noqa: E402
    carter_coefficient,
    effective_air_gap,
    magnetizing_inductance_per_phase,
    slotted_air_gap_permeance_factor,
    synchronous_magnetizing_inductance,
    winding_factor,
)


OUT_JSON = HERE / "validation_carter_magnetizing_sweep_summary.json"

SLOTS = 36
POLES = 4
POLE_PAIRS = POLES // 2
PHASES = 3
GAP_DIAMETER = 0.10
STACK_LENGTH = 0.10
GAP = 0.8e-3
TURNS_PER_PHASE = 120
PITCH_FRACTION = 5.0 / 6.0
SLOT_OPENINGS = (0.0, 0.5e-3, 1.0e-3, 1.5e-3, 2.5e-3, 4.0e-3, 5.5e-3)


def _slot_pitch() -> float:
    return math.pi * GAP_DIAMETER / SLOTS


def _slot_angle_electrical_deg() -> float:
    return 360.0 * POLE_PAIRS / SLOTS


def build_rows() -> list[dict]:
    slot_pitch = _slot_pitch()
    kw1 = winding_factor(
        1,
        slots_per_pole_per_phase=SLOTS / (POLES * PHASES),
        slot_angle_elec_deg=_slot_angle_electrical_deg(),
        pitch_fraction=PITCH_FRACTION,
    )
    smooth_lmu = magnetizing_inductance_per_phase(
        GAP_DIAMETER,
        STACK_LENGTH,
        GAP,
        POLE_PAIRS,
        kw1,
        TURNS_PER_PHASE,
    )
    rows = []
    for opening in SLOT_OPENINGS:
        kc = carter_coefficient(slot_pitch, GAP, opening)
        geff = effective_air_gap(slot_pitch, GAP, opening)
        permeance = slotted_air_gap_permeance_factor(slot_pitch, GAP, opening)
        lmu = magnetizing_inductance_per_phase(
            GAP_DIAMETER,
            STACK_LENGTH,
            geff,
            POLE_PAIRS,
            kw1,
            TURNS_PER_PHASE,
        )
        lm_sync = synchronous_magnetizing_inductance(
            GAP_DIAMETER,
            STACK_LENGTH,
            geff,
            POLE_PAIRS,
            kw1,
            TURNS_PER_PHASE,
            phases=PHASES,
        )
        rows.append({
            "slot_opening_m": opening,
            "slot_opening_over_gap": opening / GAP,
            "carter_coefficient": kc,
            "effective_gap_m": geff,
            "permeance_factor": permeance,
            "L_mu_per_phase_H": lmu,
            "L_m_sync_H": lm_sync,
            "L_mu_over_smooth": lmu / smooth_lmu,
        })
    return rows


def validate(rows: list[dict]) -> dict:
    kc = [row["carter_coefficient"] for row in rows]
    permeance = [row["permeance_factor"] for row in rows]
    lmu = [row["L_mu_per_phase_H"] for row in rows]
    max_lmu_identity_error = max(
        abs(row["L_mu_over_smooth"] - row["permeance_factor"])
        for row in rows
    )
    max_sync_factor_error = max(
        abs(row["L_m_sync_H"] / row["L_mu_per_phase_H"] - PHASES / 2.0)
        for row in rows
    )
    checks = {
        "slot_pitch_m": _slot_pitch(),
        "kw1": winding_factor(
            1,
            slots_per_pole_per_phase=SLOTS / (POLES * PHASES),
            slot_angle_elec_deg=_slot_angle_electrical_deg(),
            pitch_fraction=PITCH_FRACTION,
        ),
        "smooth_L_mu_H": rows[0]["L_mu_per_phase_H"],
        "largest_opening_carter_coefficient": rows[-1]["carter_coefficient"],
        "largest_opening_L_mu_drop_fraction": 1.0 - rows[-1]["L_mu_over_smooth"],
        "max_L_mu_vs_permeance_identity_error": max_lmu_identity_error,
        "max_sync_factor_error": max_sync_factor_error,
        "kc_monotone_increasing": all(a <= b for a, b in zip(kc, kc[1:])),
        "permeance_monotone_decreasing": all(a >= b for a, b in zip(permeance, permeance[1:])),
        "L_mu_monotone_decreasing": all(a >= b for a, b in zip(lmu, lmu[1:])),
    }
    assert math.isclose(rows[0]["carter_coefficient"], 1.0, rel_tol=1e-12)
    assert math.isclose(rows[0]["permeance_factor"], 1.0, rel_tol=1e-12)
    assert checks["kc_monotone_increasing"]
    assert checks["permeance_monotone_decreasing"]
    assert checks["L_mu_monotone_decreasing"]
    assert checks["largest_opening_carter_coefficient"] > 1.3
    assert checks["largest_opening_L_mu_drop_fraction"] > 0.20
    assert max_lmu_identity_error < 1.0e-12
    assert max_sync_factor_error < 1.0e-12
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    rows = build_rows()
    checks = validate(rows)
    summary = {
        "kind": "carter_magnetizing_inductance_validation",
        "validation_class": True,
        "slots": SLOTS,
        "poles": POLES,
        "phases": PHASES,
        "gap_diameter_m": GAP_DIAMETER,
        "stack_length_m": STACK_LENGTH,
        "physical_gap_m": GAP,
        "turns_per_phase": TURNS_PER_PHASE,
        "pitch_fraction": PITCH_FRACTION,
        "rows": rows,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("[Carter air-gap / magnetizing inductance sweep]")
    for row in rows:
        print(
            f"  bo/g={row['slot_opening_over_gap']:5.2f}  "
            f"kC={row['carter_coefficient']:.6f}  "
            f"P/P0={row['permeance_factor']:.6f}  "
            f"Lmu={1e3 * row['L_mu_per_phase_H']:.6f} mH"
        )
    print(
        "[checks] "
        f"kw1={checks['kw1']:.6f}, "
        f"smooth Lmu={1e3 * checks['smooth_L_mu_H']:.6f} mH, "
        f"largest-opening drop={100.0 * checks['largest_opening_L_mu_drop_fraction']:.3f}%"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
