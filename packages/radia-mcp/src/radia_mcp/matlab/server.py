import json,sys
from mcp.server.fastmcp import FastMCP
from radia_mcp.common import register_status_tool
from . import matlab_agent_guide as _guide
from .runtime import matlab_extension_contract as _contract,matlab_official_server_config as _config,matlab_radia_acoustic_interface_contract as _boundary
mcp=FastMCP("mcp-server-radia-matlab")
@mcp.tool()
def matlab_agent_guide()->str: return _guide()
@mcp.tool()
def matlab_extension_contract()->str: return json.dumps(_contract(),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_official_server_config(profile:str="existing",include_generic_extension:bool=False)->str: return json.dumps(_config(profile,include_generic_extension=include_generic_extension),ensure_ascii=False,indent=2)
@mcp.tool()
def matlab_radia_acoustic_interface_contract()->str: return json.dumps(_boundary(),ensure_ascii=False,indent=2)
register_status_tool(mcp,server_name="mcp-server-radia-matlab",description="Official MATLAB MCP composition and generic ML/RL gates",subpackage="radia_mcp.matlab",related_servers=["radia-ngsolve"])
def main():
    if "--selftest" in sys.argv: c=_contract(); assert c["ok"] and c["tool_count"]==43; print("Radia MATLAB MCP server self-test: OK"); return
    mcp.run()
