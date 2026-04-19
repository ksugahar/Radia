"""CadQuery → STEP → Cubit hex mesh end-to-end via cadquery_to_cubit_hex.

CadQuery interop demo: AI authors a bracket in cadquery syntax, our
MCP server routes through execute_cadquery → STEP → cubit_mesh_auto →
live GUI. One call covers the full pipeline.
"""
from __future__ import annotations
import subprocess, sys, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
subprocess.run(["taskkill", "/F", "/IM", "coreform_cubit.exe"],
               capture_output=True, check=False)

# Simple L-shaped bracket: 2D profile extruded → clean prism → hex sweeps trivially.
CQ_SCRIPT = """
# L-bracket: 2D polyline profile, extruded
pts = [(0, 0), (60, 0), (60, 10), (10, 10), (10, 40), (0, 40), (0, 0)]
result = (
    Workplane("XY")
    .polyline(pts).close()
    .extrude(15)
)
"""

from radia_mcp.build123d import server as bs

print("=== cadquery_to_cubit_hex: L-bracket ===")
r = json.loads(bs.cadquery_to_cubit_hex(
    script=CQ_SCRIPT,
    target_size=2.0,
    prefer="hex",
    commit_to_gui=True,
    timeout_s=300,
))

print(f"\n-- status: {r.get('status')}")
cq = r.get("cadquery", {})
print(f"-- cadquery:")
print(f"   valid={cq.get('is_valid')} volume={cq.get('volume')} faces={cq.get('face_count')}")
print(f"   bbox={cq.get('bounding_box')}")

print(f"-- step: {r.get('step_path')}")

mesh = r.get("cubit_mesh_auto", {})
print(f"\n-- ladder attempts:")
for a in mesh.get("attempts", []):
    flag = "WIN" if a.get("hexes",0) > 0 else "..."
    print(f"   {flag} {a['scheme']:10s} hex={a['hexes']:>5d} tet={a['tets']:>5d} nodes={a['nodes']:>6d}")
winner = mesh.get("winner")
if winner:
    print(f"\n-- winner: scheme={winner['scheme']} hex={winner['hexes']} tet={winner['tets']}")
gui = mesh.get("gui", {})
print(f"-- GUI replay: summary={gui.get('summary')}")
