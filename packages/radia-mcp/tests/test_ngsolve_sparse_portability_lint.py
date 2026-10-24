"""Guards: compiler-generated sparse names cannot enter portable NGSolve scripts.
Assumes: static imports/attributes in parseable Python; no numerical imports.
"""
import pytest
from radia_mcp.radia_ngsolve.rules import (
    ALL_RULES, check_ngsolve_sparse_type_portability,
)


@pytest.mark.parametrize("prefix", ["SparseMatrix", "SparseMatrixSymmetric"])
@pytest.mark.parametrize("suffix", ["double", "d"])
def test_generated_sparse_names_flagged_with_import_alias(prefix, suffix):
    code = "import ngsolve as fem\nfem.la." + prefix + suffix + ".CreateFromCOO(i,j,v)"
    findings = check_ngsolve_sparse_type_portability("sample.py", code.splitlines())
    assert [(row["line"], row["severity"], row["rule"]) for row in findings] == [
        (2, "HIGH", "ngsolve-sparse-type-portability")]


@pytest.mark.parametrize("code", [
    "from ngsolve import la as sparse\nsparse.{name}.CreateFromCOO(i,j,v)",
    "import ngsolve.la as sparse\nsparse.{name}.CreateFromCOO(i,j,v)",
    "from ngsolve.la import {name} as RealMatrix",
])
def test_import_variants_are_checked(code):
    code = code.format(name="SparseMatrix" + "double")
    findings = check_ngsolve_sparse_type_portability("sample.py", code.splitlines())
    assert len(findings) == 1, findings
    assert "real_sparse_matrix_type" in findings[0]["message"]


def test_comments_strings_and_unrelated_classes_are_not_flagged():
    name = "SparseMatrix" + "double"
    code = ("# ng.la." + name + "\nexplanation = 'ng.la." + name
            + "'\nother_library." + name + "()\n"
            + "from radia.ngsolve_sparse import real_sparse_matrix_type\n"
            + "real_sparse_matrix_type(symmetric=True).CreateFromCOO(i,j,v)")
    assert check_ngsolve_sparse_type_portability("sample.py", code.splitlines()) == []


def test_rule_runs_through_registered_script_lint_rules():
    code = "import ngsolve as ng\nng.la." + "SparseMatrix" + "double"
    findings = [row for rule in ALL_RULES for row in rule("sample.py", code.splitlines())]
    assert any(row["rule"] == "ngsolve-sparse-type-portability" for row in findings)
