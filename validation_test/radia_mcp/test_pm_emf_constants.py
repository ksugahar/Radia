"""PM flux-linkage, back-EMF, and torque-constant helpers."""

import math
import os
import sys

import pytest

pytest.importorskip("ngsolve")          # solve.py imports ngsolve at module load

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (  # noqa: E402
    dq_torque,
    dq_voltages,
    pm_flux_linkage_constants,
    pm_no_load_back_emf,
)


def test_pm_flux_linkage_constants_peak_and_rms_conversions():
    row = pm_flux_linkage_constants(lambda_m=0.1, pole_pairs=4)
    assert row["back_emf_constant_phase_peak_V_per_rad_per_s_mech"] == pytest.approx(0.4)
    assert row["back_emf_constant_phase_rms_V_per_rad_per_s_mech"] == pytest.approx(0.4 / math.sqrt(2))
    assert row["back_emf_constant_line_line_rms_V_per_rad_per_s_mech"] == pytest.approx(
        math.sqrt(1.5) * 0.4
    )
    assert row["torque_constant_Nm_per_Aq_peak"] == pytest.approx(0.6)
    assert row["torque_constant_Nm_per_Aq_rms"] == pytest.approx(0.6 * math.sqrt(2))
    assert row["Kt_peak_over_phase_peak_Ke"] == pytest.approx(1.5)
    assert row["Kt_rms_over_line_line_rms_Ke"] == pytest.approx(math.sqrt(3.0))


def test_pm_no_load_back_emf_matches_dq_voltage_equation():
    lam = 0.075
    p = 3
    omega_mech = 250.0
    emf = pm_no_load_back_emf(lam, omega_mech, p)
    vd, vq = dq_voltages(R=0.0, Ld=1.0e-3, Lq=2.0e-3,
                         lambda_m=lam, id_=0.0, iq=0.0,
                         omega_e=emf["omega_e_rad_per_s"])
    assert vd == pytest.approx(0.0)
    assert vq == pytest.approx(emf["phase_peak_V"])
    assert emf["phase_rms_V"] == pytest.approx(emf["phase_peak_V"] / math.sqrt(2.0))
    assert emf["line_line_rms_V"] == pytest.approx(math.sqrt(3.0) * emf["phase_rms_V"])
    assert emf["line_line_peak_V"] == pytest.approx(math.sqrt(3.0) * emf["phase_peak_V"])


def test_pm_torque_constant_matches_dq_power_identity():
    lam = 0.08
    p = 5
    iq_peak = 12.0
    omega_mech = 180.0
    constants = pm_flux_linkage_constants(lam, p)
    emf = pm_no_load_back_emf(lam, omega_mech, p)
    torque = dq_torque(lam, Ld=2.0e-3, Lq=2.0e-3, id_=0.0, iq=iq_peak, pole_pairs=p)

    assert torque == pytest.approx(constants["torque_constant_Nm_per_Aq_peak"] * iq_peak)
    mechanical_power = torque * omega_mech
    electrical_airgap_power = 1.5 * emf["phase_peak_V"] * iq_peak
    assert mechanical_power == pytest.approx(electrical_airgap_power)


def test_pm_flux_linkage_constants_reject_bad_pole_pairs():
    with pytest.raises(ValueError):
        pm_flux_linkage_constants(0.1, 0)
    with pytest.raises(ValueError):
        pm_no_load_back_emf(0.1, 100.0, -2)


if __name__ == "__main__":
    test_pm_flux_linkage_constants_peak_and_rms_conversions()
    test_pm_no_load_back_emf_matches_dq_voltage_equation()
    test_pm_torque_constant_matches_dq_power_identity()
    test_pm_flux_linkage_constants_reject_bad_pole_pairs()
    print("[OK] PM Ke/Kt helpers validated against dq voltage, torque, and power identities.")
