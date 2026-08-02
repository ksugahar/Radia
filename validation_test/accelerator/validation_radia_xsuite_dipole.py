"""Validation-class Radia field solve to Xsuite magnetic tracking bridge.

The case uses a fixed-magnet Radia model and sends its batched ``rad.Fld``
values directly to Xsuite's spatial Boris integrator.  A zero-field control,
step-refinement check, momentum invariant, and tracking-boundary accounting
keep the demonstration independent of any commercial trajectory solver.
"""
from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import radia as rad
import xtrack as xt

from radia.xsuite_bridge import (
    AxisAlignedBox,
    radia_magnetic_fieldmap,
    track_magnetic_fieldmap,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "validation_radia_xsuite_dipole_summary.json"


def _git_head() -> str:
    return subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True
    ).strip()


def _particles():
    return xt.Particles(
        "proton",
        p0c=1.0e9,
        x=[-0.005, 0.0, 0.005],
        y=[0.05, 0.05, 0.05],
    )


def _zero_field(x, y, z):
    return np.zeros_like(x), np.zeros_like(y), np.zeros_like(z)


def _run(fieldmap, n_steps: int):
    particles = _particles()
    result = track_magnetic_fieldmap(
        particles,
        fieldmap,
        s_start_m=-0.12,
        s_end_m=0.12,
        n_steps=n_steps,
        boundary=AxisAlignedBox(
            minimum=(-0.05, 0.0, -0.13),
            maximum=(0.05, 0.10, 0.13),
        ),
    )
    return particles, result


def main() -> int:
    rad.UtiDelAll()
    magnet = rad.magnet_box(
        center=[0.0, 0.0, 0.0],
        dimensions=[0.04, 0.04, 0.08],
        magnetization=[0.0, 9.5e5, 0.0],
    )

    fieldmap = radia_magnetic_fieldmap(magnet)
    off_particles, off = _run(_zero_field, 1200)
    on_coarse_particles, on_coarse = _run(fieldmap, 600)
    on_fine_particles, on_fine = _run(fieldmap, 1200)

    off_xy = np.column_stack((off_particles.x, off_particles.y))
    coarse_xy = np.column_stack((on_coarse_particles.x, on_coarse_particles.y))
    fine_xy = np.column_stack((on_fine_particles.x, on_fine_particles.y))
    endpoint_deflection = np.linalg.norm(fine_xy - off_xy, axis=1)
    refinement_difference = np.linalg.norm(fine_xy - coarse_xy, axis=1)
    delta_before = np.asarray(on_fine["relative_momentum_deviation_before"])
    delta_after = np.asarray(on_fine["relative_momentum_deviation_after"])
    momentum_drift = float(np.max(np.abs(delta_after - delta_before)))
    field_probe = np.asarray(rad.Fld(magnet, "b", [0.0, 0.05, 0.0]), dtype=float)

    checks = {
        "radia_field_nonzero": float(np.linalg.norm(field_probe)) > 0.05,
        "zero_field_is_straight": float(np.max(np.abs(off_xy[:, 0] - [-0.005, 0.0, 0.005])))
        < 1.0e-12,
        "magnetic_field_deflects_every_trajectory": bool(np.all(endpoint_deflection > 1.0e-4)),
        "boris_step_refinement_closes": float(np.max(refinement_difference)) < 5.0e-7,
        "magnetic_tracking_preserves_momentum": momentum_drift <= 2.0e-12,
        "tracking_boundary_has_no_exit": not on_fine["boundary_exit_events"],
        "particle_count_closes": on_fine["particle_count"] == 3,
    }
    summary = {
        "schema": "radia-validation-xsuite-dipole/v1",
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "execution_version": {
            "producer": "validation_radia_xsuite_dipole.py",
            "radia": rad.__version__,
            "radia_git_head": _git_head(),
            "xtrack": xt.__version__,
        },
        "case": "fixed_magnet_external_field_to_xsuite_spatial_boris",
        "units": {"position": "m", "magnetic_flux_density": "T", "momentum": "eV/c"},
        "particle_count": 3,
        "step_counts": {"coarse": 600, "fine": 1200},
        "field_probe_t": field_probe.tolist(),
        "endpoint_deflection_m": endpoint_deflection.tolist(),
        "coarse_fine_endpoint_difference_m": refinement_difference.tolist(),
        "maximum_relative_momentum_drift": momentum_drift,
        "boundary_exit_events": on_fine["boundary_exit_events"],
        "checks": checks,
        "limitations": [
            "magnetic_field_only",
            "electrostatic_acceleration_not_yet_in_this_bridge",
            "particle_matter_interactions_require_xcoll",
        ],
        "pass": all(checks.values()),
    }
    OUTPUT.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    rad.UtiDelAll()
    print(json.dumps(summary, indent=2))
    return 0 if summary["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
