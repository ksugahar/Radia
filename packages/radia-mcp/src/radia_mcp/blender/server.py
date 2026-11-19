"""MCP Server: radia_mcp.blender  (3-D model animation through headless Blender)

Gmsh owns Radia mesh/field inspection and sections; this domain owns 3-D model
animation.  One typed operation, ``blender_field_animation``, maps a solver
dataset (``radia.field-animation.v1`` HDF5: one solved vector field, J in A/m^2
or B in T, optionally with the magnet pose of the same solver steps) to a ``.blend``, preview frames and an MP4,
then reopens the ``.blend`` in a separate process and compares it with the
dataset.  The Blender worker is optional, runs as an explicit Python 3.13 +
bpy interpreter with ``-I -B`` and is not shipped in the wheel.  There is no
generic script execution and no socket listener.

Usage:
    radia-mcp serve radia-mcp-postprocess --profile blender
    python -m radia_mcp.blender.server --selftest   # no bpy needed
"""
import json
import sys

from mcp.server.fastmcp import FastMCP

from radia_mcp.common import register_status_tool
from .knowledge import GUIDE

mcp = FastMCP("blender-mcp")

@mcp.tool(meta={"caeai.maturity": "artifact_gate"})
def blender_field_animation(dataset_path: str, out_dir: str, spec_json: str = "{}") -> str:
    """Animate a solved vector field (J or B, optional moving magnet) in Blender, then verify it.

    Input: a `radia.field-animation.v1` HDF5 file (SI; /time_s, /magnet/center_m,
    /samples/id, /samples/point_m, exactly one of /field/J_A_per_m2 or /field/B_T,
    optional /magnet/center_m and /magnet/rotation), for
    example from docs/maglev/demos/rotating_magnet_eddy_animation_export.py.
    Writes <stem>.blend, <stem>_stepNNNN.png previews and <stem>.mp4 into
    out_dir, then reopens the .blend in a separate headless process and
    compares frame count, fps, physical time per frame, the fixed |J| scale,
    magnet pose and the signed J vector of every sample identifier.

    spec_json (optional keys): fps, samples, threads, resolution [w, h], title,
    arrow_length_m, arrow_width_m, glyph_stride, preview_steps [solver steps],
    video (bool), field_scale ("global" or a positive value in the field unit).

    Requires RADIA_BLENDER_PYTHON (Python 3.13 with bpy) and, outside a source
    checkout, RADIA_BLENDER_WORKER.  Returns JSON with artifact paths, SHA-256,
    dataset provenance, timing and readback errors.  status "mismatch" means
    the reopened .blend did not reproduce the dataset.  Rendering quality is
    not numerical validation.
    """
    # Imported on use: the module needs numpy (optional), and tool discovery must
    # work under the minimum SDK without optional numeric dependencies.
    from .field_animation import BlenderWorkerError, field_animation

    try:
        spec = json.loads(spec_json or "{}")
        if not isinstance(spec, dict):
            raise ValueError("spec_json must be a JSON object")
        return json.dumps(field_animation(dataset_path, out_dir, spec), indent=1)
    except (BlenderWorkerError, ValueError, OSError, KeyError) as exc:
        return json.dumps({"status": "error", "kind": type(exc).__name__, "message": str(exc)})


@mcp.tool()
def blender_usage() -> str:
    """Usage, input contract, fidelity rules and headless verification for Blender animation."""
    return GUIDE


# Registered last: it attaches titles/annotations to the tools above.
register_status_tool(
    mcp,
    server_name="blender-mcp", tool_name="blender_status",
    description="Solver-faithful 3-D model animation through an optional headless Blender worker",
    subpackage="radia_mcp.blender",
    related_servers=["gmsh-mcp"],
    optional_deps=["h5py"],
)


def main():
    if "--selftest" in sys.argv:
        print("blender MCP server self-test:")
        print(f"  guide: {len(GUIDE)} chars")
        from .field_animation import BlenderWorkerError, resolve_worker
        try:
            py, worker, version = resolve_worker()
            print(f"  worker: {worker} (bpy {version})")
        except BlenderWorkerError as exc:
            print(f"  worker: unavailable ({exc})")
        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
