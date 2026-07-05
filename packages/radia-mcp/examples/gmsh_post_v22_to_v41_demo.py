"""gmsh_post end-to-end smoke demo for a user-provided Gmsh v2.2 mesh."""
from __future__ import annotations
import os
import sys, json, io
from pathlib import Path
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

from radia_mcp.gmsh_post import server as g

V22 = os.environ.get("RADIA_MCP_GMSH_POST_V22_MSH", "")
OUT_DIR = Path(os.environ.get("RADIA_MCP_GMSH_POST_OUT", r"C:\temp"))
V41_OUT = str(OUT_DIR / "roundtrip_v41.msh")
V41_WITH_DATA = str(OUT_DIR / "roundtrip_v41_with_data.msh")

if not V22:
    raise SystemExit(
        "Set RADIA_MCP_GMSH_POST_V22_MSH to a local Gmsh v2.2 .msh file "
        "before running this smoke demo."
    )

# 1. inspect v2.2 source
print("=== 1. inspect v2.2 Cubit export ===")
r = json.loads(g.gmsh_post_inspect(V22))
print(f"  version: {r['format_version']}  binary: {r['binary']}")
print(f"  sections: {r['sections']}")
print(f"  nodes: {r.get('node_count')}, elements: {r.get('total_elements')}")
print(f"  element_counts: {r.get('element_counts')}")
print(f"  physical_groups: {len(r.get('physical_groups', []))}")
print(f"  entities by dim: {r.get('entities',{}).get('by_dim')}")

# 2. validate v2.2 (should fail: wrong version)
print("\n=== 2. validate v2.2 (should ERROR on version) ===")
r = json.loads(g.gmsh_post_validate(V22))
print(f"  passed_all: {r['passed_all']}  errors: {r['errors']}  warns: {r['warnings']}")
for f in r['findings']:
    print(f"  [{f['severity']}] {f['rule']}: {f['message']}")

# 3. convert v2.2 → v4.1
print("\n=== 3. convert v2.2 → v4.1 ===")
r = json.loads(g.gmsh_post_convert(V22, V41_OUT))
print(f"  status={r['status']} output_version={r.get('output_version')}")
print(f"  {r.get('input_size')} bytes → {r.get('output_size')} bytes")

# 4. validate v4.1
print("\n=== 4. validate v4.1 (should PASS) ===")
r = json.loads(g.gmsh_post_validate(V41_OUT))
print(f"  passed_all: {r['passed_all']}  errors: {r['errors']}  warns: {r['warnings']}")
for f in r['findings']:
    print(f"  [{f['severity']}] {f['rule']}: {f['message']}")

# 5. write $NodeData — scalar field (r² distance from origin) on each node
print("\n=== 5. write $NodeData (r² at nodes) ===")
r = json.loads(g.gmsh_post_inspect(V41_OUT))
n_count = r['node_count']
# synthesize per-node values using gmsh directly
import gmsh
gmsh.initialize([], False)
gmsh.option.setNumber('General.Terminal', 0)
gmsh.open(V41_OUT)
tags, coords, _ = gmsh.model.mesh.getNodes()
values = {}
for i, t in enumerate(tags):
    x, y, z = coords[3*i], coords[3*i+1], coords[3*i+2]
    values[int(t)] = [x*x + y*y + z*z]
gmsh.finalize()
r = json.loads(g.gmsh_post_write_node_data(
    msh_path=V41_OUT, out_path=V41_WITH_DATA,
    view_name='r_squared',
    values=values,
))
print(f"  status={r['status']} values_written={r.get('values_written')}")

# 6. re-inspect the result
print("\n=== 6. inspect output with $NodeData ===")
r = json.loads(g.gmsh_post_inspect(V41_WITH_DATA))
print(f"  version: {r['format_version']}  sections: {r['sections']}")
print(f"  post_views: {r.get('post_views')}")
