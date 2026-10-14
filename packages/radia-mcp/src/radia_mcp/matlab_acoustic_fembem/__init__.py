"""Compatibility imports for the former combined MATLAB/FEM-BEM namespace."""

from ..acoustic_fembem import (
    acoustic_fembem_agent_guide as matlab_acoustic_fembem_agent_guide,
    acoustic_fembem_extension_contract as matlab_acoustic_fembem_extension_contract,
    acoustic_fembem_extension_path as matlab_acoustic_fembem_extension_path,
    acoustic_fembem_server_config as matlab_acoustic_fembem_server_config,
)


__all__ = [
    "matlab_acoustic_fembem_agent_guide",
    "matlab_acoustic_fembem_extension_contract",
    "matlab_acoustic_fembem_extension_path",
    "matlab_acoustic_fembem_server_config",
]
