"""Fast contracts for the planar HDiv reduced-motor operating path."""
from __future__ import annotations

import math

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
        vol="rotor.vol",
        mu_r_iron="1200",
        rotor_angle_steps=5,
    )

    command = spec.build_command(python="python", panels_dir="panels")

    assert any("calc_motor_hdiv_reduced.py" in part for part in command)
    assert command[command.index("--mu-r") + 1] == "1200"
    assert command[command.index("--rotor-angle-steps") + 1] == "5"
    assert spec.visible_fields() >= {
        "vol", "h_amplitude", "r_airgap_mid", "energy_delta_deg"
    }
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
