"""Validation-class two-port S-parameter sweep health example.

Run:

    python examples/rf_waveguide/validation_two_port_sparameter_health.py

This example audits a small power-normalized two-port table using three views:
passivity (``|S11|^2 + |S21|^2 <= 1``), reciprocity (``S21 ~= S12``), and the
momentum force implied by one-sided port-1 excitation.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.force import two_port_sparameter_sweep_health_summary  # noqa: E402


OUT_JSON = HERE / "validation_two_port_sparameter_health_summary.json"


def build_summary() -> dict[str, object]:
    power_W = 3.0
    frequencies = [3.0e9, 4.0e9, 5.0e9]
    s11 = [0.05 + 0.0j, 0.20 + 0.0j, 0.45 + 0.0j]
    s21 = [0.95 + 0.0j, 0.80 + 0.0j, 0.55 + 0.0j]
    s12 = [0.95005 + 0.0j, 0.79995 + 0.0j, 0.55001 + 0.0j]
    s22 = [0.05002 + 0.0j, 0.19999 + 0.0j, 0.45003 + 0.0j]
    health = two_port_sparameter_sweep_health_summary(
        frequencies,
        s11,
        s21,
        s12_values=s12,
        s22_values=s22,
        power_incident_W=power_W,
        reciprocity_tolerance=1.0e-3,
        return_symmetry_tolerance=1.0e-3,
    )
    nonreciprocal_probe = two_port_sparameter_sweep_health_summary(
        frequencies,
        s11,
        s21,
        s12_values=[0.95 + 0.0j, 0.72 + 0.0j, 0.55 + 0.0j],
        power_incident_W=power_W,
        reciprocity_tolerance=1.0e-3,
    )
    active_probe = two_port_sparameter_sweep_health_summary(
        [6.0e9],
        [0.8 + 0.0j],
        [0.7 + 0.0j],
        power_incident_W=power_W,
        passivity_tolerance=1.0e-6,
    )

    checks = {
        "n_points": health["n_points"],
        "status": health["status"],
        "passivity_ok": health["passivity_ok"],
        "reciprocity_ok": health["reciprocity_ok"],
        "return_symmetry_ok": health["return_symmetry_ok"],
        "max_s21_s12_abs_error": health["max_s21_s12_abs_error"],
        "max_s11_s22_abs_error": health["max_s11_s22_abs_error"],
        "max_force_frequency_Hz": health["max_force_frequency_Hz"],
        "nonreciprocal_probe_status": nonreciprocal_probe["status"],
        "nonreciprocal_probe_error": nonreciprocal_probe["max_s21_s12_abs_error"],
        "active_probe_status": active_probe["status"],
        "active_probe_excess": active_probe["max_passivity_excess_power_fraction"],
    }

    assert checks["n_points"] == 3
    assert checks["status"] == "ok"
    assert checks["passivity_ok"] is True
    assert checks["reciprocity_ok"] is True
    assert checks["return_symmetry_ok"] is True
    assert abs(float(checks["max_s21_s12_abs_error"]) - 5.0e-5) < 1.0e-12
    assert abs(float(checks["max_s11_s22_abs_error"]) - 3.0e-5) < 1.0e-12
    assert checks["max_force_frequency_Hz"] == 5.0e9
    assert checks["nonreciprocal_probe_status"] == "needs_attention"
    assert abs(float(checks["nonreciprocal_probe_error"]) - 0.08) < 1.0e-12
    assert checks["active_probe_status"] == "needs_attention"
    assert checks["active_probe_excess"] > 0.0

    return {
        "kind": "two_port_sparameter_health_validation",
        "validation_class": True,
        "rf_learning": (
            "two-port S-parameter sweeps should separate passivity, reciprocity, "
            "return symmetry, and momentum-force diagnostics"
        ),
        "checks": checks,
        "health": health,
        "nonreciprocal_probe": nonreciprocal_probe,
        "active_probe": active_probe,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[two-port S-parameter health]")
    print(
        f"  n_points={checks['n_points']} status={checks['status']} "
        f"passivity_ok={checks['passivity_ok']} reciprocity_ok={checks['reciprocity_ok']}"
    )
    print(
        f"  max |S21-S12|={checks['max_s21_s12_abs_error']:.3e} "
        f"max |S11-S22|={checks['max_s11_s22_abs_error']:.3e}"
    )
    print(
        f"  probes: nonreciprocal={checks['nonreciprocal_probe_status']} "
        f"active={checks['active_probe_status']} excess={checks['active_probe_excess']:.3e}"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
