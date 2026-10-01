"""Write the PRIMA port-model exchanges used by test_prima_port_model.m.

The exporter needs only NumPy/SciPy, so it is loaded from its source file
without importing the native ``radia`` package.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
_SPEC = importlib.util.spec_from_file_location(
    "radia_prima_export", ROOT / "src" / "radia" / "prima_export.py")
prima_export = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = prima_export
_SPEC.loader.exec_module(prima_export)


def coupled_rl_admittance():
    """Two coupled RL loops driven at loop 1: proper admittance, improper impedance."""
    L = np.array([[1e-3, 4e-4], [4e-4, 2e-3]])
    R = np.diag([0.5, 2.0])
    return prima_export.DescriptorPortModel(
        L, -R, np.array([[1.0], [0.0]]), np.array([[1.0, 0.0]]), np.zeros((1, 1)),
        "voltage", ("coil",), {"case": "coupled_rl"})


def rc_two_port_impedance():
    """Two-node RC network with both nodes as ports: proper MIMO impedance."""
    C = np.diag([1e-6, 2e-6])
    G = np.array([[1.5e-3, -1e-3], [-1e-3, 3e-3]])
    return prima_export.DescriptorPortModel(
        C, -G, np.eye(2), np.eye(2), np.zeros((2, 2)), "current",
        ("node1", "node2"), {"case": "rc_two_port"})


def write_all(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    rl = coupled_rl_admittance()
    prima_export.export_prima_lti_json(rl, directory / "rl_admittance.json")
    prima_export.export_prima_lti_json(rl, directory / "rl_impedance.json",
                                       orientation="impedance")
    prima_export.export_prima_lti_json(rc_two_port_impedance(), directory / "rc_impedance.json")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: prima_port_model_python_reference.py <output directory>")
    write_all(Path(sys.argv[1]))
