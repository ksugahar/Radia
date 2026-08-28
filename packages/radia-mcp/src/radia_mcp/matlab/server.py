import json,sys
from mcp.server.fastmcp import FastMCP
from radia_mcp.common import register_status_tool
from radia_mcp.matlab_agentic_ml import (
    matlab_agentic_ml_guide as _ml_guide,
    validate_matlab_ml_rl_artifact as _validate_ml_rl_artifact,
    validate_aicia_catalog as _validate_aicia_catalog,
)
from . import matlab_agent_guide as _guide
from .optuna_boundary import matlab_optuna_mcp_route as _optuna_mcp_route
from .optuna_quality import (
    matlab_optuna_benchmark_plan as _optuna_benchmark_plan,
    matlab_optuna_health as _optuna_health,
    matlab_optuna_oracle_plan as _optuna_oracle_plan,
    matlab_optuna_release_gate as _optuna_release_gate,
)
from .optimize import matlab_cad_topology_build as _cad_topology_build, matlab_optimize_build as _optimize_build, matlab_optimize_resume as _optimize_resume, matlab_sheet_metal_topology_build as _sheet_metal_topology_build
from .optuna_oracle import (
    matlab_optuna_compatibility_contract as _optuna_compatibility,
    matlab_optuna_oracle_audit as _optuna_audit,
)
from .runtime import (
    matlab_extension_contract as _contract,
    matlab_official_server_config as _config,
    matlab_radia_acoustic_interface_contract as _boundary,
    matlab_radia_mex_contract as _mex_contract,
    matlab_optuna_simulink_contract as _optuna_contract,
    matlab_simulink_library_contract as _simulink_library_contract,
)
mcp=FastMCP("mcp-server-radia-matlab")
@mcp.tool()
def matlab_agent_guide()->str: return _guide()
@mcp.tool()
def matlab_agentic_ml_guide()->str: return _ml_guide()
@mcp.tool()
def matlab_ml_rl_artifact_gate(artifact_json:str)->str:
    artifact=json.loads(artifact_json)
    return json.dumps(_validate_ml_rl_artifact(artifact),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_aicia_catalog_gate(catalog_json:str)->str:
    """Validate full-channel metadata scope and solver-gated CAE promotion."""
    return json.dumps(_validate_aicia_catalog(json.loads(catalog_json)),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_extension_contract()->str: return json.dumps(_contract(),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_official_server_config(profile:str="existing",include_generic_extension:bool=False)->str: return json.dumps(_config(profile,include_generic_extension=include_generic_extension),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_radia_acoustic_interface_contract()->str: return json.dumps(_boundary(),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_radia_mex_contract(topic:str="all")->str:
    """Expose the shared Radia/NGSolve Python-to-MATLAB MEX capability contract."""
    return json.dumps(_mex_contract(topic),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_optuna_simulink_contract()->str:
    """Describe the table-backed MATLAB Optuna and Simulink workflow."""
    return json.dumps(_optuna_contract(),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_optuna_mcp_route(topic:str="overview")->str:
    """Route shared Optuna tools upstream and MATLAB differences to Radia."""
    return json.dumps(_optuna_mcp_route(topic),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_optuna_health(distribution_path:str="")->str:
    """Check manifest, full upstream API mapping, MEX, Simulink, and notices."""
    return json.dumps(_optuna_health(distribution_path),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_optuna_oracle_plan(scope:str="all",repository_path:str="",output_path:str="")->str:
    """Build an official-MATLAB-MCP-ready Optuna differential test plan."""
    return json.dumps(_optuna_oracle_plan(scope,repository_path,output_path),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_optuna_benchmark_plan(repository_path:str="",output_directory:str=r"C:\temp\radia-optuna-benchmark")->str:
    """Build the same-host cold/warmed MATLAB-versus-Optuna benchmark plan."""
    return json.dumps(_optuna_benchmark_plan(repository_path,output_directory),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_optuna_release_gate(evidence_json:str,repository_path:str="",max_warmed_time_ratio:float=1.0)->str:
    """Validate installed-wheel, Simulink, resume, oracle, and speed evidence."""
    return json.dumps(_optuna_release_gate(evidence_json,repository_path,max_warmed_time_ratio),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_optuna_compatibility_contract(repo_root:str="")->str:
    """Report verified Optuna 4.9.0 parity, MATLAB-only behavior, and gaps."""
    return json.dumps(_optuna_compatibility(repo_root or None),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_optuna_oracle_audit(repo_root:str="")->str:
    """Audit Optuna fixtures, versions, hashes, policy sync, and test coverage."""
    return json.dumps(_optuna_audit(repo_root or None),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_simulink_library_contract()->str:
    """Describe Radia application blocks, Library Browser registration, and LTspice compatibility."""
    return json.dumps(_simulink_library_contract(),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_optimize_build(spec_json:str)->str:
    """Build validated MATLAB code for objective, Simulink, or LTspice optimization."""
    return json.dumps(_optimize_build(spec_json),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_optimize_resume(storage_path:str,n_trials:int,parallel:bool=False)->str:
    """Build official-MATLAB-MCP-ready code to resume a persisted Study."""
    return json.dumps(_optimize_resume(storage_path,n_trials,parallel=parallel),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_cad_topology_build(spec_json:str)->str:
    """Build a Cubit + Radia-VIM linearization + LP topology workflow."""
    return json.dumps(_cad_topology_build(spec_json),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_sheet_metal_topology_build(spec_json:str)->str:
    """Build a Radia-VIM + LP + adaptive NGSolve/Cubit sheet-metal workflow."""
    return json.dumps(_sheet_metal_topology_build(spec_json),ensure_ascii=False,indent=2)
register_status_tool(mcp,server_name="mcp-server-radia-matlab",description="Official MATLAB MCP composition, Radia/NGSolve MEX bridge, table-backed optimization, Simulink, and generic ML/RL gates",subpackage="radia_mcp.matlab",related_servers=["radia-ngsolve"])
def main():
    if "--selftest" in sys.argv: c=_contract(); assert c["ok"] and c["tool_count"]==43; print("Radia MATLAB MCP server self-test: OK"); return
    mcp.run()
