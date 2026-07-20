"""MCP integration for the readable MATLAB acoustic FEM-BEM solver."""

from pathlib import Path

from .integration import (
    acoustic_fembem_extension_contract,
    acoustic_fembem_extension_path,
    acoustic_fembem_server_config,
)


_HERE = Path(__file__).resolve().parent


def acoustic_fembem_agent_guide() -> str:
    """Return the educational acoustic FEM-BEM agent guide."""
    return (_HERE / "skill.md").read_text(encoding="utf-8")


__all__ = [
    "acoustic_fembem_agent_guide",
    "acoustic_fembem_extension_contract",
    "acoustic_fembem_extension_path",
    "acoustic_fembem_server_config",
]
