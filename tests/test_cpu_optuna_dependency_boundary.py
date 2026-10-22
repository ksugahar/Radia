"""CPU distribution dependency ownership.

Guards: Accidental coupling of CPU Radia requirements to the independently
released MATLAB distribution, and drift of the upstream Optuna extra.
Assumes: Root pyproject.toml is the CPU dependency source of truth. This source
contract reads only that metadata and does not import either runtime package.
"""

from pathlib import Path
import re
import tomllib


def test_cpu_optuna_extra_uses_only_the_upstream_distribution():
    project = tomllib.loads(
        (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text(
            encoding="utf-8"
        )
    )["project"]
    extras = project["optional-dependencies"]
    assert extras["optuna"] == ["optuna==5.0.0"], extras["optuna"]
    assert "optuna-upstream" not in extras, sorted(extras)
    requirements = project["dependencies"] + [
        requirement for group in extras.values() for requirement in group
    ]
    names = [re.split(r"[\s\[<>=!~;]", req, maxsplit=1)[0].lower().replace("_", "-")
             for req in requirements]
    assert "radia-optuna" not in names, requirements
