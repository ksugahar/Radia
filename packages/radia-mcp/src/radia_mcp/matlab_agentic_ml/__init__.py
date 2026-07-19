"""Public MATLAB ML/RL workflow guidance and result-contract checks."""

from .artifact_gate import (
    validate_matlab_ml_rl_artifact,
    validate_matlab_ml_rl_v44_identity,
    validate_matlab_ml_rl_v45_identity,
    validate_matlab_ml_rl_v46_identity,
    validate_matlab_ml_rl_v47_identity,
    validate_matlab_ml_rl_v48_identity,
    validate_matlab_ml_rl_v49_identity,
)
from .v50_identity import validate_matlab_ml_rl_v50_identity


def matlab_agentic_ml_guide() -> str:
    """Return the public MATLAB machine-learning/RL workflow guide."""
    from pathlib import Path

    return (Path(__file__).with_name("skill.md")).read_text(encoding="utf-8")


__all__ = [
    "matlab_agentic_ml_guide",
    "validate_matlab_ml_rl_artifact",
    "validate_matlab_ml_rl_v44_identity",
    "validate_matlab_ml_rl_v45_identity",
    "validate_matlab_ml_rl_v46_identity",
    "validate_matlab_ml_rl_v47_identity",
    "validate_matlab_ml_rl_v48_identity",
    "validate_matlab_ml_rl_v49_identity",
    "validate_matlab_ml_rl_v50_identity",
]
