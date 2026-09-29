"""Acceptance contracts for the isolated IH coarse-solver study."""
import importlib.util
from pathlib import Path

import unittest

PATH = Path(__file__).resolve().parents[1] / "validation_test/induction_heating/bddc_ams_coarse_ih.py"
spec = importlib.util.spec_from_file_location("ih_coarse_study", PATH)
study = importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)


def record(solver="bddc", **changes):
    value = dict(case="tube", mesh="coarse", order=2, solver=solver,
                 process_exit=0, error=None, linear_true_relative_residual=1e-9,
                 P_total=2.0, L=1e-8)
    value.update(changes)
    return value


class TestAcceptance(unittest.TestCase):
    def test_matching_converged_solvers_pass(self):
        self.assertTrue(study.compare_results([record("sparsecholesky"), record()])["all_passed"])

    def test_invalid_result_cannot_pass(self):
        for changes in (
            {"process_exit": 1}, {"timeout_s": 10}, {"error": "failed"},
            {"linear_true_relative_residual": 2e-7},
            {"linear_true_relative_residual": float("nan")},
            {"P_total": float("nan")}, {"P_total": 2.01}, {"L": 1.1e-8},
        ):
            with self.subTest(changes=changes):
                self.assertFalse(study.compare_results([record("sparsecholesky"), record(**changes)])["all_passed"])

    def test_missing_same_mesh_direct_reference_cannot_pass(self):
        self.assertFalse(study.compare_results([record()])["all_passed"])
        self.assertFalse(study.compare_results([record("sparsecholesky", mesh="fine"), record()])["all_passed"])
        self.assertFalse(study.compare_results([])["all_passed"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
