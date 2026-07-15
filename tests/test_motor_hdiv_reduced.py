"""Fast contracts for the planar HDiv reduced-motor operating path."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
from netgen.occ import OCCGeometry, WorkPlane  # noqa: E402

from radia.motor_design import ANALYSIS_HDIV_REDUCED, MotorDesignSpec  # noqa: E402
from radia.motor_hdiv import HDivReducedMotor  # noqa: E402


def test_motor_design_builds_hdiv_reduced_cli_command():
    spec = MotorDesignSpec(
        analysis=ANALYSIS_HDIV_REDUCED,
        rotor_vol="rotor.vol",
        hdiv_mu_r="1200",
        rotor_angle_steps=5,
    )

    command = spec.build_command(python="python", panels_dir="panels")

    assert any("calc_motor_hdiv_reduced.py" in part for part in command)
    assert command[command.index("--mu-r") + 1] == "1200"
    assert command[command.index("--rotor-angle-steps") + 1] == "5"
    assert spec.visible_fields() >= {
        "rotor_vol", "hdiv_h_amplitude", "r_airgap_mid", "energy_delta_deg"
    }
    assert "vol" not in spec.visible_fields()
    defaults = MotorDesignSpec(
        analysis=ANALYSIS_HDIV_REDUCED, rotor_vol="rotor.vol").build_command(
            python="python", panels_dir="panels")
    assert defaults[defaults.index("--mu-r") + 1] == "1000.0"
    assert defaults[defaults.index("--H-amplitude") + 1] == "80000.0"
    with pytest.raises(ValueError, match="Unknown motor analysis"):
        MotorDesignSpec(analysis="unknown").visible_fields()


def test_native_planar_field_cf_matches_explicit_frame_transform():
    mesh = ng.Mesh(
        OCCGeometry(WorkPlane().Ellipse(0.2, 0.1).Face(), dim=2)
        .GenerateMesh(maxh=0.05))
    source_angle = 0.3
    target_angle = -0.2
    point_target = np.array([0.02, 0.01])

    with ng.TaskManager():
        motor = HDivReducedMotor(mesh, 1001.0)
        state = motor.solve_angle(0.0, (8.0e4, 0.0))
        field_cf = motor.body.field_cf(
            state.coefficients,
            source_angle=source_angle,
            target_angle=target_angle,
        )

        delta = source_angle-target_angle
        c, s = math.cos(delta), math.sin(delta)
        rotation = np.array(((c, -s), (s, c)))
        point_source = point_target @ rotation
        field_source = motor.body.H_at([point_source], state.coefficients)[0]
        expected_target = field_source @ rotation.T
        actual_target = np.asarray(
            field_cf(mesh(float(point_target[0]), float(point_target[1]))))

    np.testing.assert_allclose(actual_target, expected_target, rtol=2e-13, atol=2e-10)


def test_hdiv_reduced_cli_round_trip(tmp_path):
    root = Path(__file__).resolve().parents[1]
    vol = tmp_path / "rotor.vol"
    output = tmp_path / "result.json"
    ngmesh = OCCGeometry(
        WorkPlane().Ellipse(0.2, 0.1).Face(), dim=2
    ).GenerateMesh(maxh=0.08)
    ngmesh.Save(str(vol))

    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src") + os.pathsep + env.get("PYTHONPATH", "")
    completed = subprocess.run(
        [
            sys.executable,
            str(root / "src/radia/panels/calc_motor_hdiv_reduced.py"),
            "--vol", str(vol),
            "--mu-r", "1001",
            "--H-amplitude", "80000",
            "--rotor-angle-start-deg", "-20",
            "--rotor-angle-stop-deg", "20",
            "--rotor-angle-steps", "3",
            "--maxwell-radius", "0.28",
            "--energy-delta-deg", "0.1",
            "--circle-points", "360",
            "--output", str(output),
        ],
        cwd=root,
        env=env,
        capture_output=True,
        text=True,
        timeout=120,
        check=True,
    )

    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["analysis"] == "hdiv_reduced_motor"
    assert result["gram_build_count"] == 1
    assert len(result["angles"]) == 3
    assert max(row["torque_spread_relative"] for row in result["angles"]) < 1e-4
    assert "Gram built once" in completed.stderr
