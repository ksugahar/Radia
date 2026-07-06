"""MATLAB acoustic FEM-BEM agent guidance."""

from __future__ import annotations

from pathlib import Path


_HERE = Path(__file__).resolve().parent


def matlab_acoustic_fembem_agent_guide() -> str:
    """Return the MATLAB acoustic FEM-BEM agent skill document."""
    return (_HERE / "skill.md").read_text(encoding="utf-8")


__all__ = ["matlab_acoustic_fembem_agent_guide"]
