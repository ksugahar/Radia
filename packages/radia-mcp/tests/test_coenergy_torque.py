"""Coenergy-derived virtual-work torque helpers."""

import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    coenergy_torque_from_angle_samples,
    coenergy_torque_summary,
    virtual_work_force_from_displacement_samples,
    virtual_work_force_summary,
    virtual_work_force_sweep_audit_summary,
    virtual_work_symmetric_pair_force_summary,
)


def test_coenergy_torque_periodic_sinusoid_matches_analytic_derivative():
    samples = 1440
    harmonic = 3
    amplitude = 0.75
    angles = [2.0 * math.pi * index / samples for index in range(samples)]
    coenergy = [2.0 - amplitude * math.cos(harmonic * angle) for angle in angles]
    rows = coenergy_torque_from_angle_samples(angles, coenergy, periodic=True)

    exact = [amplitude * harmonic * math.sin(harmonic * angle) for angle in angles]
    max_abs_error = max(abs(row["torque_Nm"] - ref) for row, ref in zip(rows, exact))
    assert rows[0]["stencil"] == "central_periodic"
    assert max_abs_error < 7.0e-5

    summary = coenergy_torque_summary(angles, coenergy, periodic=True)
    assert summary["n_samples"] == samples
    assert summary["torque_peak_abs_Nm"] == pytest.approx(amplitude * harmonic, rel=4.0e-5)
    assert summary["torque_mean_Nm"] == pytest.approx(0.0, abs=1.0e-13)


def test_coenergy_torque_nonperiodic_linear_table_is_exact():
    angles = [0.0, 0.1, 0.25, 0.4]
    coenergy = [1.0 + 4.0 * angle for angle in angles]
    rows = coenergy_torque_from_angle_samples(angles, coenergy)

    assert [row["stencil"] for row in rows] == ["forward", "central", "central", "backward"]
    assert [row["torque_Nm"] for row in rows] == pytest.approx([4.0, 4.0, 4.0, 4.0])


def test_coenergy_torque_rejects_bad_tables():
    with pytest.raises(ValueError):
        coenergy_torque_from_angle_samples([0.0, 1.0], [0.0, 1.0])
    with pytest.raises(ValueError):
        coenergy_torque_from_angle_samples([0.0, 1.0, 2.0], [0.0, 1.0])
    with pytest.raises(ValueError):
        coenergy_torque_from_angle_samples([0.0, 1.0, 1.0], [0.0, 1.0, 2.0])
    with pytest.raises(ValueError):
        coenergy_torque_from_angle_samples([0.0, 1.0, 2.0], [0.0, 1.0, 2.0], periodic=True, period_rad=0.0)


def test_virtual_work_force_linear_coenergy_table_is_exact():
    positions = [-0.002, -0.001, 0.0, 0.001, 0.002]
    expected_force = 12.5
    coenergy = [0.25 + expected_force * x for x in positions]

    rows = virtual_work_force_from_displacement_samples(positions, coenergy)
    summary = virtual_work_force_summary(positions, coenergy, energy_kind="constant_current")

    assert [row["stencil"] for row in rows] == ["forward", "central", "central", "central", "backward"]
    assert [row["energy_kind"] for row in rows] == ["coenergy"] * len(rows)
    assert [row["force_N"] for row in rows] == pytest.approx([expected_force] * len(rows))
    assert summary["virtual_work_identity"] == "F = dW_co/dx at fixed current"
    assert summary["force_mean_N"] == pytest.approx(expected_force)
    assert summary["force_peak_abs_N"] == pytest.approx(expected_force)


def test_virtual_work_force_stored_energy_flips_derivative_sign():
    positions = [-0.002, -0.001, 0.0, 0.001, 0.002]
    expected_force = 8.0
    stored_energy = [0.125 - expected_force * x for x in positions]

    rows = virtual_work_force_from_displacement_samples(
        positions,
        stored_energy,
        energy_kind="field_energy",
    )

    assert [row["energy_kind"] for row in rows] == ["stored_energy"] * len(rows)
    assert [row["denergy_dx_N"] for row in rows] == pytest.approx([-expected_force] * len(rows))
    assert [row["force_N"] for row in rows] == pytest.approx([expected_force] * len(rows))
    assert rows[2]["virtual_work_identity"] == "F = -dW/dx at fixed flux/source-free displacement"


def test_virtual_work_force_quadratic_center_matches_analytic_gradient():
    positions = [-0.02, -0.01, 0.0, 0.01, 0.02]
    stiffness = 300.0
    stored_energy = [0.5 * stiffness * x * x for x in positions]
    rows = virtual_work_force_from_displacement_samples(
        positions,
        stored_energy,
        energy_kind="stored_energy",
    )

    assert rows[1]["force_N"] == pytest.approx(-stiffness * positions[1])
    assert rows[2]["force_N"] == pytest.approx(0.0)
    assert rows[3]["force_N"] == pytest.approx(-stiffness * positions[3])


def test_virtual_work_force_sweep_audit_compares_central_rows():
    positions = [-0.002, -0.001, 0.0, 0.001, 0.002]
    constant_force = 3.0
    stiffness = 2000.0
    coenergy = [
        0.25 + constant_force * x + 0.5 * stiffness * x * x
        for x in positions
    ]
    reference = [constant_force + stiffness * x for x in positions]

    summary = virtual_work_force_sweep_audit_summary(
        positions,
        coenergy,
        reference_force_N=reference,
        force_abs_tolerance_N=1.0e-12,
    )

    assert summary["status"] == "ok"
    assert summary["reference_checked_count"] == 3
    assert summary["comparison_stencils"] == ["central"]
    assert summary["max_reference_force_abs_error_N"] < 1.0e-12
    assert summary["max_abs_force_gradient_N_per_m"] == pytest.approx(stiffness)
    assert summary["rows"][0]["selected_for_reference_check"] is False
    assert summary["rows"][1]["selected_for_reference_check"] is True
    assert summary["rows"][2]["force_gradient_N_per_m"] == pytest.approx(stiffness)

    bad = virtual_work_force_sweep_audit_summary(
        positions,
        coenergy,
        reference_force_N=[value + 0.5 for value in reference],
        force_abs_tolerance_N=1.0e-12,
    )
    assert bad["status"] == "needs_attention"
    assert bad["reference_pass"] is False


def test_virtual_work_symmetric_pair_force_summary_matches_center_gradient():
    h = 1.0e-3
    expected_force = 12.5
    stiffness = 1000.0
    offset = 0.25
    e_minus = offset - expected_force * h + 0.5 * stiffness * h * h
    e_center = offset
    e_plus = offset + expected_force * h + 0.5 * stiffness * h * h

    row = virtual_work_symmetric_pair_force_summary(
        h,
        e_minus,
        e_plus,
        energy_kind="constant_current",
        energy_center_J=e_center,
    )

    assert row["energy_kind"] == "coenergy"
    assert row["virtual_work_identity"] == "F = dW_co/dx at fixed current"
    assert row["denergy_dx_N"] == pytest.approx(expected_force)
    assert row["force_N"] == pytest.approx(expected_force)
    assert row["even_energy_residual_J"] == pytest.approx(0.5 * stiffness * h * h)
    assert row["position_minus_m"] == pytest.approx(-h)
    assert row["position_plus_m"] == pytest.approx(h)

    stored = virtual_work_symmetric_pair_force_summary(
        h,
        offset + expected_force * h,
        offset - expected_force * h,
        energy_kind="stored_energy",
    )
    assert stored["denergy_dx_N"] == pytest.approx(-expected_force)
    assert stored["force_N"] == pytest.approx(expected_force)
    assert stored["even_energy_residual_J"] is None


def test_virtual_work_force_rejects_bad_tables():
    with pytest.raises(ValueError):
        virtual_work_force_from_displacement_samples([0.0, 1.0], [0.0, 1.0])
    with pytest.raises(ValueError):
        virtual_work_force_from_displacement_samples([0.0, 1.0, 2.0], [0.0, 1.0])
    with pytest.raises(ValueError):
        virtual_work_force_from_displacement_samples([0.0, 1.0, 1.0], [0.0, 1.0, 2.0])
    with pytest.raises(ValueError):
        virtual_work_force_from_displacement_samples([0.0, 1.0, 2.0], [0.0, 1.0, 2.0], energy_kind="unknown")
    with pytest.raises(ValueError):
        virtual_work_symmetric_pair_force_summary(0.0, 0.0, 1.0)
    with pytest.raises(ValueError):
        virtual_work_symmetric_pair_force_summary(1.0e-3, 0.0, 1.0, energy_kind="unknown")
