"""CPU distribution dependency ownership.

Guards: Accidental coupling of CPU Radia requirements to the independently
released MATLAB distribution, and drift of the upstream Optuna extra.
Assumes: Root pyproject.toml is the CPU dependency source of truth. This source
contract reads only that metadata and does not import either runtime package.
"""

from pathlib import Path
import tomllib

import pytest
from packaging.requirements import Requirement
from packaging.utils import canonicalize_name


def _assert_cpu_dependency_boundary(project):
    extras = project["optional-dependencies"]
    assert extras["optuna"] == ["optuna==5.0.0"], extras["optuna"]
    assert "optuna-upstream" not in extras, sorted(extras)
    requirements = project["dependencies"] + [
        requirement for group in extras.values() for requirement in group
    ]
    names = [canonicalize_name(Requirement(req).name) for req in requirements]
    assert "radia-optuna" not in names, (
        f"CPU radia must not depend on radia-optuna: {requirements}"
    )


def test_cpu_optuna_extra_uses_only_the_upstream_distribution():
    project = tomllib.loads(
        (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text(
            encoding="utf-8"
        )
    )["project"]
    _assert_cpu_dependency_boundary(project)


@pytest.mark.parametrize("requirement", [
    "radia-optuna@https://x/y.whl",
    " radia.optuna(>=1)",
    "radia--optuna>=1",
])
def test_cpu_boundary_rejects_pep508_distribution_spellings(requirement):
    project = {"dependencies": [requirement],
               "optional-dependencies": {"optuna": ["optuna==5.0.0"]}}
    with pytest.raises(AssertionError, match="CPU radia must not depend"):
        _assert_cpu_dependency_boundary(project)


@pytest.mark.parametrize("extra", ["optuna", "optuna-upstream"])
def test_cpu_boundary_rejects_invalid_optuna_extras(extra):
    project = {"dependencies": [],
               "optional-dependencies": {"optuna": ["optuna==5.0.0"]}}
    project["optional-dependencies"][extra] = ["radia-optuna"]
    with pytest.raises(AssertionError):
        _assert_cpu_dependency_boundary(project)
