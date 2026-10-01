"""PRIMA port models: orientation, pole-residue form, LTspice text, LTI exchange."""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pytest

from radia import prima_export as pe

sys.path.insert(0, str(Path(__file__).resolve().parent / "matlab"))
from prima_port_model_python_reference import (  # noqa: E402
    coupled_rl_admittance, rc_two_port_impedance)


def _ltspice_response(text: str, s: complex) -> tuple[str, np.ndarray]:
    """Port matrix that the controlled sources of a generated .subckt implement."""
    orientation = re.search(r", (\w+) realization", text).group(1)
    ports = sorted({int(n) for n in re.findall(r"\bP(\d+)_p\b", text)})
    H = np.zeros((len(ports), len(ports)), dtype=complex)
    pattern = re.compile(r"^[EG]\S+_P(\d+) \S+ \S+ (?:u_P(\d+) 0|P(\d+)_p P\d+_n) (.+)$")
    for line in text.splitlines():
        match = pattern.match(line)
        if not match:
            continue
        k = int(match.group(1)) - 1
        j = int(match.group(2) or match.group(3)) - 1
        value = match.group(4)
        if value.startswith("Laplace="):
            expression = value[len("Laplace="):]
            assert re.fullmatch(r"[0-9eE.+\-*/()s]+", expression)
            H[k, j] += eval(expression, {"__builtins__": {}}, {"s": s})  # noqa: S307
        else:
            H[k, j] += float(value)
    return orientation, H


def test_inverse_is_the_exact_other_orientation():
    model = coupled_rl_admittance()
    inverse = model.inverse()
    assert (model.orientation, inverse.orientation) == ("admittance", "impedance")
    s = 2j * math.pi * 300.0
    assert np.allclose(inverse.response(s) @ model.response(s), np.eye(1), rtol=0, atol=1e-13)
    assert model.oriented("admittance") is model


def test_pole_residue_recovers_the_analytic_coupled_rl_impedance():
    # Z = R1 + sL1 - (sM)^2 / (R2 + sL2) = 0.58 + 9.2e-4 s - 80 / (s + 1000).
    form = pe.pole_residue(coupled_rl_admittance().inverse())
    assert np.allclose(form.poles, [-1000.0], rtol=1e-12)
    assert np.allclose(form.residues.real, [[[-80.0]]], rtol=1e-9)
    assert np.allclose(form.D0, [[0.58]], rtol=1e-12)
    assert np.allclose(form.D1, [[9.2e-4]], rtol=1e-12)


def test_complex_pole_pairs_become_one_real_second_order_section():
    E = np.diag([1e-3, 1e-6])
    A = np.array([[-1.0, -1.0], [1.0, 0.0]])          # series R-L-C admittance
    model = pe.DescriptorPortModel(E, A, np.array([[1.0], [0.0]]), np.array([[1.0, 0.0]]),
                                   np.zeros((1, 1)), "voltage", ("rlc",))
    sections = pe.pole_residue(model).real_sections()
    assert [section["order"] for section in sections] == [2]
    a1, a0 = sections[0]["den"]
    assert a1 == pytest.approx(1e3) and a0 == pytest.approx(1e9)   # s^2 + (R/L) s + 1/(LC)
    assert sections[0]["num"][0][0, 0] == pytest.approx(1e3)       # s / L


def test_unstable_and_malformed_models_fail_loudly():
    with pytest.raises(ValueError, match="unstable poles"):
        pe.pole_residue(pe.DescriptorPortModel(
            np.eye(1), np.eye(1), np.eye(1), np.eye(1), np.zeros((1, 1)), "current", ("p",)))
    with pytest.raises(ValueError, match="input_quantity"):
        pe.DescriptorPortModel(np.eye(1), -np.eye(1), np.eye(1), np.eye(1), np.zeros((1, 1)),
                               "charge", ("p",))
    with pytest.raises(ValueError, match="shape"):
        pe.DescriptorPortModel(np.eye(2), -np.eye(2), np.eye(1), np.eye(1), np.zeros((1, 1)),
                               "current", ("p",))
    with pytest.raises(ValueError, match="real"):
        pe.DescriptorPortModel(np.eye(1), -1j * np.eye(1), np.eye(1), np.eye(1),
                               np.zeros((1, 1)), "current", ("p",))


@pytest.mark.parametrize("model, realization", [
    (coupled_rl_admittance(), "admittance"),
    (coupled_rl_admittance().inverse(), "admittance"),   # improper impedance: realized as Y
    (rc_two_port_impedance(), "impedance"),
])
def test_ltspice_sources_implement_the_port_relation(model, realization):
    text = pe.ltspice_subckt(model, name="DUT")
    assert ".subckt DUT P1_p P1_n" in text and text.rstrip().endswith(".ends DUT")
    assert "+-" not in text
    native = model.oriented(realization)
    for frequency in (1.0, 160.0, 3e3, 1e6):
        s = 2j * math.pi * frequency
        orientation, H = _ltspice_response(text, s)
        assert orientation == realization
        assert np.allclose(H, native.response(s), rtol=1e-9, atol=0)


def test_ltspice_refuses_models_improper_both_ways():
    # Admittance diag(1/(sL), sC): an inductor port and a capacitor port, so
    # Y grows like s at port 2 and Z = diag(sL, 1/(sC)) at port 1.
    inductor = pe.DescriptorPortModel(np.eye(1) * 1e-3, np.zeros((1, 1)), np.eye(1), np.eye(1),
                                      np.zeros((1, 1)), "voltage", ("L",))
    capacitor = pe.DescriptorPortModel(np.eye(1) * 1e-6, np.zeros((1, 1)), np.eye(1), np.eye(1),
                                       np.zeros((1, 1)), "current", ("C",)).inverse()
    blocks = [(inductor.E, capacitor.E), (inductor.A, capacitor.A), (inductor.B, capacitor.B),
              (inductor.C, capacitor.C), (inductor.D, capacitor.D)]
    E, A, B, C, D = (np.block([[x, np.zeros((x.shape[0], y.shape[1]))],
                               [np.zeros((y.shape[0], x.shape[1])), y]]) for x, y in blocks)
    model = pe.DescriptorPortModel(E, A, B, C, D, "voltage", ("L", "C"))
    with pytest.raises(ValueError, match="improper in both orientations"):
        pe.ltspice_subckt(model)


def test_lti_exchange_is_row_major_and_self_checking(tmp_path):
    model = rc_two_port_impedance()
    path = pe.export_prima_lti_json(model, tmp_path / "rc.json", metadata={"note": "x"})
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["schema"] == pe.SCHEMA
    assert (payload["input_unit"], payload["output_unit"]) == ("A", "V")
    G = np.array(payload["arrays"]["A"]["values"]).reshape(payload["arrays"]["A"]["shape"])
    assert np.array_equal(G, model.A)
    frequency = payload["check"]["frequency_hz"][3]
    check = payload["check"]
    response = (np.array(check["response_real"]["values"])
                + 1j * np.array(check["response_imag"]["values"])).reshape(
                    check["response_real"]["shape"])
    assert np.allclose(response[3], model.response(2j * math.pi * frequency), rtol=1e-14)
    assert payload["metadata"] == {"case": "rc_two_port", "note": "x"}


def test_prima_hacapk_descriptor_model_is_stable_and_exported_exactly():
    """The MNA projection keeps every reduced pole stable (the symmetric
    saddle projection left unstable out-of-band poles)."""
    from radia.coil_from_cad import build_peec_from_path, helix_path
    from radia.prima_hacapk import PRIMAHACApKModel

    n_seg = 36
    path = helix_path(radius=0.01, pitch=0.003, n_turns=3, n_points=n_seg + 1)
    topo = build_peec_from_path(path, np.full(n_seg, 1e-3), np.full(n_seg, 1e-3),
                                sigma=5.8e7, nwinc=3, nhinc=3)
    prima = PRIMAHACApKModel.from_topology(topo, q=12)
    assert prima.build_stats["is_descriptor"]
    model = pe.from_prima_hacapk(prima, port_name="helix")
    for frequency in (1e2, 1e4, 1e6):
        assert model.response(2j * math.pi * frequency)[0, 0] == pytest.approx(
            prima.port_impedance(frequency), rel=1e-12)
        assert prima.port_impedance(frequency).real > 0.0
    form = pe.pole_residue(model)
    assert np.all(form.poles.real < 0.0)
    L_dc, R_dc = prima.dc_inductance_resistance()
    z_low = prima.port_impedance(1.0)
    assert R_dc == pytest.approx(z_low.real, rel=1e-6)
    assert L_dc == pytest.approx(z_low.imag / (2.0 * math.pi), rel=1e-6)
    assert "impedance realization" in pe.ltspice_subckt(model, name="HELIX")
