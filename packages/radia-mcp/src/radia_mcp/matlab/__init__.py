from pathlib import Path
from .runtime import (
    matlab_extension_contract,
    matlab_extension_path,
    matlab_official_server_config,
    matlab_radia_acoustic_interface_contract,
    matlab_radia_mex_contract,
    matlab_optuna_simulink_contract,
    matlab_simulink_library_contract,
)
from .optimize import matlab_cad_topology_build, matlab_optimize_build, matlab_optimize_resume, matlab_sheet_metal_topology_build
def matlab_agent_guide(): return (Path(__file__).parent/"skill.md").read_text(encoding="utf-8")
__all__=["matlab_agent_guide","matlab_extension_contract","matlab_extension_path","matlab_official_server_config","matlab_radia_acoustic_interface_contract","matlab_radia_mex_contract","matlab_optuna_simulink_contract","matlab_simulink_library_contract","matlab_optimize_build","matlab_optimize_resume","matlab_cad_topology_build","matlab_sheet_metal_topology_build"]
