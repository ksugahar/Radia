from pathlib import Path
from .runtime import matlab_extension_contract, matlab_extension_path, matlab_official_server_config, matlab_radia_acoustic_interface_contract
def matlab_agent_guide(): return (Path(__file__).parent/"skill.md").read_text(encoding="utf-8")
__all__=["matlab_agent_guide","matlab_extension_contract","matlab_extension_path","matlab_official_server_config","matlab_radia_acoustic_interface_contract"]
