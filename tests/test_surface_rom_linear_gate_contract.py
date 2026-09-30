import importlib.util
import sys
import types

import numpy as np
import pytest


def _load_transient(monkeypatch):
    monkeypatch.setitem(sys.modules, "ngsolve", types.ModuleType("ngsolve"))
    monkeypatch.setitem(sys.modules, "team13_model", types.ModuleType("team13_model"))
    gate_path = SOURCE.parents[2] / "src" / "radia" / "_residual_gate.py"
    gate_spec = importlib.util.spec_from_file_location(
        "radia._residual_gate", gate_path)
    gate = importlib.util.module_from_spec(gate_spec)
    gate_spec.loader.exec_module(gate)
    radia = types.ModuleType("radia")
    radia.__path__ = []
    monkeypatch.setitem(sys.modules, "radia", radia)
    monkeypatch.setitem(sys.modules, "radia._residual_gate", gate)
    spec = importlib.util.spec_from_file_location("surface_rom_transient", SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SOURCE = (Path(__file__).resolve().parents[1] / "validation_test" /
          "surface_rom_team10" / "transient.py")


class _IdentityMatrix:
    @staticmethod
    def CSR():
        return (np.array([1.0]), np.array([0]), np.array([0, 1]))


def test_legacy_surface_rom_fom_rejects_bad_linear_correction(monkeypatch):
    transient = _load_transient(monkeypatch)
    free = np.array([True])
    with pytest.raises(RuntimeError, match="true relative residual"):
        transient._check_linear_correction(
            _IdentityMatrix(), np.array([1.0]), np.array([0.0]),
            np.array([1.0]), free, "deliberately broken solve", 1.0)


def test_legacy_surface_rom_fom_accepts_accurate_linear_correction(monkeypatch):
    transient = _load_transient(monkeypatch)
    relative = transient._check_linear_correction(
        _IdentityMatrix(), np.array([1.0e-8]), np.array([1.0]),
        np.array([1.0]), np.array([True]), "accurate solve", 1.0)
    assert relative == pytest.approx(1.0e-8)
