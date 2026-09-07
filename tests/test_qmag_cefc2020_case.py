"""Fast contracts of the CEFC 2020 quadrupole validation lane (validation_test/quadrupole_cefc2020).

The heavy HDiv-MMM solves run on hibino; these tests lock what can be checked in seconds: the
racetrack coil set is a clean quadrupole with the documented sign convention, the iron law loads
with its documented shape, and the generated Cubit journal is the conforming imprint/merge +
per-volume z-sweep contract without entity ids.
"""
from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
LANE = REPO / "validation_test" / "quadrupole_cefc2020"


def _load_case():
    spec = importlib.util.spec_from_file_location("qmag_case", LANE / "qmag_case.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["qmag_case"] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def qmag():
    return _load_case()


def test_coil_set_is_a_quadrupole_with_the_documented_sign(qmag):
    rad = pytest.importorskip("radia")
    rad.UtiDelAll()
    coil, manifest = qmag.build_qmag_coils(3.0)
    assert len(manifest["coils"]) == 8
    points = qmag.observation_points()
    field = np.asarray(rad.Fld(coil, "b", points), dtype=float)
    perp = qmag.b_perp(field)
    r = np.arange(-15, 16, dtype=float) * 1.0e-3
    # B_x = G x, B_y = -G y: zero at the centre, odd along the diagonal, B_z = 0 in the midplane.
    # Segmented arcs leave ~1e-7 relative asymmetry between mirrored racetracks.
    assert np.abs(field[15]).max() < 1.0e-6 * np.abs(perp).max()
    assert np.abs(perp + perp[::-1]).max() < 1.0e-6 * np.abs(perp).max()
    assert np.abs(field[:, 2]).max() < 1.0e-6 * np.abs(perp).max()
    assert np.allclose(field[:, 0], -field[:, 1])
    gradient = perp[16:] / r[16:]
    assert gradient.max() < 0.0 and gradient.min() > -1.1 and np.ptp(gradient) / abs(gradient.mean()) < 0.01
    assert perp[-1] < 0.0
    on_axis = np.asarray(rad.Fld(coil, "b", [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0]]), dtype=float)
    assert on_axis[0, 0] > 0.0 and abs(on_axis[0, 1]) < 1.0e-9
    assert on_axis[1, 1] < 0.0 and abs(on_axis[1, 0]) < 1.0e-9
    assert math.isclose(on_axis[0, 0], -on_axis[1, 1], rel_tol=1.0e-9)
    rad.UtiDelAll()


def test_iron_law_is_a_monotone_h_b_table(qmag):
    table = np.asarray(qmag.load_bh_table(), dtype=float)
    assert table.shape == (49, 2)
    assert table[0].tolist() == [0.0, 0.0]
    assert np.all(np.diff(table[:, 0]) > 0.0) and np.all(np.diff(table[:, 1]) > 0.0)
    assert table[-1, 1] == pytest.approx(2.516, abs=1.0e-3) and table[-1, 0] == pytest.approx(2.58e5, rel=1.0e-2)
    mu_r_initial = table[1, 1] / (qmag.MU0 * table[1, 0])
    assert 1.0e4 < mu_r_initial < 1.5e4


def test_three_engine_cube_average_is_exact_for_a_linear_field():
    spec = importlib.util.spec_from_file_location("run_qmag_three_engine", LANE / "run_qmag_three_engine.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["run_qmag_three_engine"] = module
    spec.loader.exec_module(module)
    centres = np.array([[0.001, 0.001, 0.0], [-0.004, -0.004, 0.0]])
    sub = module._cube_points(centres, 2.0e-5)
    assert sub.shape == (16, 3)
    # a quadrupole-like linear field B = (G x, -G y, 0) is averaged exactly by the 2x2x2 Gauss cube
    field = np.column_stack([15.0 * sub[:, 0], -15.0 * sub[:, 1], np.zeros(len(sub))])
    averaged = module._cube_average(field, len(centres))
    expected = np.column_stack([15.0 * centres[:, 0], -15.0 * centres[:, 1], np.zeros(2)])
    assert np.allclose(averaged, expected, atol=1.0e-15)
    assert module._relative_rms(expected, expected) == 0.0
    assert module._relative_rms(expected, 1.01 * expected) == pytest.approx(0.01)


def test_mesh_journal_is_conforming_and_id_free(qmag, tmp_path):
    journal = tmp_path / "qmag.jou"
    qmag.write_mesh_journal(journal, tmp_path / "qmag.vol", size_m=0.01, order=2)
    text = journal.read_text(encoding="ascii")
    assert text.index("import step") < text.index("imprint volume all") < text.index("merge volume all") \
        < text.index("mesh volume all") < text.index("export netgen")
    assert "#{Loop(_n)}" in text and "#{EndLoop}" in text
    sweep = [line for line in text.splitlines() if "scheme sweep" in line]
    assert len(sweep) == 1 and "z_coord < -0.02994" in sweep[0] and "z_coord > 0.02994" in sweep[0]
    assert 'block 1 name "iron"' in text and 'sideset 1 name "outer_boundary"' in text
    assert "order 2 overwrite" in text
    assert not any(line.strip().startswith("volume ") and line.split()[1].isdigit() for line in text.splitlines())
    assert qmag.STEP_PATH.is_file() and qmag.STEP_PATH.stat().st_size > 100_000
    # the TET comparison route shares the import / imprint / merge / labels and only swaps the scheme
    tet = tmp_path / "qmag_tet.jou"
    qmag.write_mesh_journal(tet, tmp_path / "qmag_tet.vol", size_m=0.01, order=2, scheme="tet")
    tet_text = tet.read_text(encoding="ascii")
    assert "volume all scheme tetmesh" in tet_text and "scheme sweep" not in tet_text and "#{Loop" not in tet_text
    assert tet_text.index("merge volume all") < tet_text.index("scheme tetmesh") < tet_text.index("mesh volume all")
    assert 'block 1 name "iron"' in tet_text and "order 2 overwrite" in tet_text
    with pytest.raises(ValueError, match="scheme"):
        qmag.write_mesh_journal(tet, tmp_path / "x.vol", size_m=0.01, order=2, scheme="hex")


def test_three_engine_partial_state_round_trip(tmp_path):
    """A non-converged mixed Omega loop leaves a partial state that only the SAME problem may resume."""
    module = sys.modules.get("run_qmag_three_engine")
    if module is None:
        spec = importlib.util.spec_from_file_location("run_qmag_three_engine", LANE / "run_qmag_three_engine.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules["run_qmag_three_engine"] = module
        spec.loader.exec_module(module)
    identity = {"case": {"name": "J3.0"}, "fem_mesh_sha256": "abc", "fem_order": 2, "nonlinear_tolerance": 2.0e-5}
    stats = {"element_numbers": [3, 1, 2], "mu_r_elements": [10.0, 200.0, 3000.0], "iterations": 80,
             "relative_B_change": 7.8e-4, "contraction_rate_estimate": 0.97, "history": [{"iteration": 1}]}
    path = tmp_path / "three_engine_J3.0.mixed_total_reduced_omega.state.json"
    module._write_state(path, identity, stats, None)
    state = module._read_state(path, identity)
    assert state["converged"] is False and state["schema"] == module.STATE_SCHEMA
    assert state["mu_r_elements"] == stats["mu_r_elements"] and state["element_numbers"] == [3, 1, 2]
    assert state["iterations"] == 80 and state["resumed_iterations"] == 0
    # a second non-converged run accumulates the resumed iterations
    module._write_state(path, identity, stats | {"iterations": 50}, state)
    assert module._read_state(path, identity)["resumed_iterations"] == 80
    assert module._read_state(tmp_path / "missing.json", identity) is None
    with pytest.raises(ValueError, match="different problem"):
        module._read_state(path, identity | {"fem_order": 3})
    path.write_text(path.read_text(encoding="utf-8").replace('"converged": false', '"converged": true'), encoding="utf-8")
    with pytest.raises(ValueError, match="not a partial"):
        module._read_state(path, identity)


def test_multipole_expansion_recovers_a_synthetic_quadrupole_with_b6():
    """A field B_y + i B_x = C2 (z/r) + C6 (z/r)^5 is recovered to round-off by the FFT expansion."""
    spec = importlib.util.spec_from_file_location("run_qmag_multipoles", LANE / "run_qmag_multipoles.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["run_qmag_multipoles"] = module
    spec.loader.exec_module(module)
    points = module.circle_points(0.015, 64)
    z = (points[:, 0] + 1j * points[:, 1]) / 0.015
    complex_field = 0.2 * z + 0.2 * 30.0e-4 * z ** 5 - 0.2 * 2.0e-4 * 1j * z ** 9
    field = np.column_stack([complex_field.imag, complex_field.real, np.zeros(len(z))])
    m = module.multipoles(field)
    assert m["B2_T"] == pytest.approx(0.2)
    assert m["b6_units"] == pytest.approx(30.0, abs=1e-9) and m["a6_units"] == pytest.approx(0.0, abs=1e-9)
    assert m["a10_units"] == pytest.approx(-2.0, abs=1e-9) and m["b10_units"] == pytest.approx(0.0, abs=1e-9)
    assert m["b14_units"] == pytest.approx(0.0, abs=1e-9)
    assert all(m[f"forbidden_b{n}_units"] < 1e-9 for n in (3, 4, 5))
