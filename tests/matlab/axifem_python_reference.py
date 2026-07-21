from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.io import savemat


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import radia.axifem as axifem  # noqa: E402


def evaluate(values: np.ndarray) -> dict[str, np.ndarray]:
    result = axifem.q1_magnetic_element_matrices(*values.tolist())
    return {
        "input": values,
        "stiffness": np.asarray(result["stiffness"], dtype=float),
        "sigma_mass": np.asarray(result["sigma_mass"], dtype=float),
    }


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: axifem_python_reference.py OUTPUT.mat")
    mu0 = 4.0 * np.pi * 1.0e-7
    sigma = 5.8e7
    interior = evaluate(np.array([
        1.0e-3, 2.0e-3, -0.5e-3, 0.5e-3, mu0, sigma], dtype=float))
    axis = evaluate(np.array([
        0.0, 1.0e-3, 0.0, 1.0e-3, mu0, sigma], dtype=float))
    thin = evaluate(np.array([
        2.0e-4, 2.1e-4, -3.0e-3, 3.0e-3, mu0, 1.2e6], dtype=float))
    insulator = evaluate(np.array([
        0.5e-3, 1.7e-3, 1.0e-3, 1.8e-3, 2.5 * mu0, 0.0], dtype=float))
    savemat(sys.argv[1], {
        "interior_input": interior["input"],
        "interior_stiffness": interior["stiffness"],
        "interior_sigma_mass": interior["sigma_mass"],
        "axis_input": axis["input"],
        "axis_stiffness": axis["stiffness"],
        "axis_sigma_mass": axis["sigma_mass"],
        "thin_input": thin["input"],
        "thin_stiffness": thin["stiffness"],
        "thin_sigma_mass": thin["sigma_mass"],
        "insulator_input": insulator["input"],
        "insulator_stiffness": insulator["stiffness"],
        "insulator_sigma_mass": insulator["sigma_mass"],
    })


if __name__ == "__main__":
    main()
