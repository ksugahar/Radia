"""Generate pybind11 reference data for the MATLAB transfer-map MEX test."""

from __future__ import annotations

import sys
import importlib.util
from pathlib import Path

import numpy as np
from scipy.io import savemat


ROOT = Path(__file__).resolve().parents[2]
# Reuse the same-source native provenance gate used by the timing study.
# The helper must not silently compare a current MEX against an installed wheel.
_contract_path = ROOT / "validation_test/accelerator/benchmark_beam_transfer_mex_python.py"
_spec = importlib.util.spec_from_file_location("beam_transfer_native_contract", _contract_path)
_contract = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_contract)
_contract._load_source_backend(ROOT)

from radia.beam import (  # noqa: E402
    canonical_body_hamiltonian_jet,
    propagate_variational_map,
)


def _matlab_segment_axis(values: np.ndarray) -> np.ndarray:
    """Move pybind's leading segment axis to MATLAB's trailing axis."""
    return np.moveaxis(np.asarray(values), 0, -1)


def _nonlinear_case() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    lengths = np.array([0.3, 0.4, 0.2], dtype=float)
    a = np.zeros((3, 6, 6), dtype=float)
    f2 = np.zeros((3, 6, 6, 6), dtype=float)
    f3 = np.zeros((3, 6, 6, 6, 6), dtype=float)
    f2[0, 1, 0, 0] = 2.0
    f2[1, 2, 0, 1] = 1.5
    f2[1, 2, 1, 0] = 1.5
    f3[2, 3, 0, 0, 0] = -0.7
    return lengths, a, f2, f3


def build_reference() -> dict[str, np.ndarray | float]:
    lengths, a, f2, f3 = _nonlinear_case()
    transfer = propagate_variational_map(
        lengths,
        a,
        f2,
        f3,
        ["upstream_sextupole", "downstream_sextupole", "direct_octupole"],
        maximum_order=3,
        maximum_step_m=1.0,
    )
    coefficients = np.array([0.2, 2.4, -0.6, 5.0, 1.5, -7.0, 2.0, 9.0, -3.0])
    canonical = canonical_body_hamiltonian_jet(
        coefficients, 3.0, reference_beta=0.8
    )
    return {
        "lengths_m": lengths,
        "A_per_m": _matlab_segment_axis(a),
        "F2_per_m": _matlab_segment_axis(f2),
        "F3_per_m": _matlab_segment_axis(f3),
        "R": np.asarray(transfer["R"]),
        "T": np.asarray(transfer["T"]),
        "U": np.asarray(transfer["U"]),
        "station_R": _matlab_segment_axis(transfer["station_R"]),
        "region_T": _matlab_segment_axis(transfer["region_T"]),
        "region_U_direct": _matlab_segment_axis(transfer["region_U_direct"]),
        "region_U_local_cascade": _matlab_segment_axis(
            transfer["region_U_local_cascade"]
        ),
        "pair_regions": np.asarray(transfer["pair_regions"], dtype=np.int32) + 1,
        "pair_U_cascade": _matlab_segment_axis(transfer["pair_U_cascade"]),
        "R_composition_error": float(transfer["diagnostics"]["R_composition_error"]),
        "T_reconstruction_error": float(
            transfer["diagnostics"]["T_reconstruction_error"]
        ),
        "U_reconstruction_error": float(
            transfer["diagnostics"]["U_reconstruction_error"]
        ),
        "T_symmetry_defect": float(transfer["diagnostics"]["T_symmetry_defect"]),
        "U_symmetry_defect": float(transfer["diagnostics"]["U_symmetry_defect"]),
        "canonical_coefficients": coefficients,
        "canonical_rigidity_t_m": 3.0,
        "canonical_reference_beta": 0.8,
        "canonical_H2_per_m": np.asarray(canonical["H2_per_m"]),
        "canonical_H3_per_m": np.asarray(canonical["H3_per_m"]),
        "canonical_H4_per_m": np.asarray(canonical["H4_per_m"]),
        "canonical_H5_per_m": np.asarray(canonical["H5_per_m"]),
        "canonical_A_per_m": np.asarray(canonical["A_per_m"]),
        "canonical_F2_per_m": np.asarray(canonical["F2_per_m"]),
        "canonical_F3_per_m": np.asarray(canonical["F3_per_m"]),
        "canonical_F4_per_m": np.asarray(canonical["F4_per_m"]),
    }


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: beam_transfer_python_reference.py OUTPUT.mat")
    savemat(sys.argv[1], build_reference(), do_compression=False, oned_as="column")


if __name__ == "__main__":
    main()
