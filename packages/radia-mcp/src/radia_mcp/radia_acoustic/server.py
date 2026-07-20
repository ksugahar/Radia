"""Production Radia acoustic MCP server (NGSolve/ngsolve.bem first)."""

from __future__ import annotations
import json, sys
from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool
from . import acoustic_capabilities as _capabilities, acoustic_usage as _usage, cq_grid_gate as _cq, fsi_preflight_gate as _fsi

mcp = FastMCP("mcp-server-radia-acoustic")

@mcp.tool()
def radia_acoustic_usage(topic: str = "overview") -> str:
    """Get NGSolve BEM, FSI, CQ, or validation guidance."""
    try: return _usage(topic)
    except ValueError as exc: return f"ERROR: {exc}"

@mcp.tool()
def radia_acoustic_capabilities() -> str:
    """List production APIs and numerical-backend ownership."""
    return json.dumps(_capabilities(), ensure_ascii=False, indent=2)

@mcp.tool()
def radia_acoustic_fsi_preflight(wavenumber: float, c_longitudinal: float = 2.0, c_transverse: float = 1.0, solid_density: float = 1.5, fluid_density: float = 1.0, order: int = 1, boundary: str = "gamma", radius_deviation: float = 0.0, radius_tolerance: float = 0.03) -> str:
    """Validate FSI/DtN inputs before an expensive NGSolve solve."""
    return json.dumps(_fsi(wavenumber=wavenumber,c_longitudinal=c_longitudinal,c_transverse=c_transverse,solid_density=solid_density,fluid_density=fluid_density,order=order,boundary=boundary,radius_deviation=radius_deviation,radius_tolerance=radius_tolerance), indent=2)

@mcp.tool()
def radia_acoustic_cq_grid(num_time: int, time_step: float, sound_speed: float = 1.0, method: str = "BDF2") -> str:
    """Build and validate Lubich CQ Laplace/complex-wavenumber grids."""
    return json.dumps(_cq(num_time=num_time,time_step=time_step,sound_speed=sound_speed,method=method), indent=2)

register_status_tool(mcp,server_name="mcp-server-radia-acoustic",description="Production NGSolve/ngsolve.bem acoustics, FSI, and convolution quadrature",subpackage="radia_mcp.radia_acoustic",related_servers=["radia-ngsolve","bem"])
def main():
    if "--selftest" in sys.argv:
        assert _cq(num_time=16,time_step=.1)["ok"] and _fsi(wavenumber=2.0)["ok"]
        print("Radia acoustic MCP server self-test: OK"); return
    mcp.run()
if __name__ == "__main__": main()
