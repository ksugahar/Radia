"""
Export and Workspace Tools for Radia MCP Server

Tools for exporting results and managing shared workspace:
- VTK export for geometry and fields
- Shared workspace export for NGSolve coupling
- Session management
"""

import sys
import json
import tempfile
from pathlib import Path
from typing import Any, Dict, List

# Add Radia to path
radia_path = Path(__file__).parent.parent.parent / "build" / "Release"
if radia_path.exists() and str(radia_path) not in sys.path:
    sys.path.insert(0, str(radia_path))

shared_path = Path(__file__).parent.parent.parent / "mcp_shared"
if shared_path.exists() and str(shared_path) not in sys.path:
    sys.path.insert(0, str(shared_path))

try:
    import radia as rad
except ImportError:
    rad = None

try:
    from workspace import SharedWorkspace
    WORKSPACE_AVAILABLE = True
except ImportError:
    WORKSPACE_AVAILABLE = False

from mcp.types import Tool


def get_tools() -> List[Tool]:
    """Get list of export and workspace tools."""
    tools = [
        Tool(
            name="radia_export_vtk_legacy",
            description="Export geometry using Radia's built-in VTK export (ObjDrwVTK).",
            inputSchema={
                "type": "object",
                "properties": {
                    "object_name": {
                        "type": "string",
                        "description": "Name of the object to export"
                    },
                    "output_file": {
                        "type": "string",
                        "description": "Output VTK file path"
                    }
                },
                "required": ["object_name", "output_file"]
            }
        ),
        Tool(
            name="radia_export_field_data_json",
            description="Export computed field data to JSON file.",
            inputSchema={
                "type": "object",
                "properties": {
                    "field_data": {
                        "type": "object",
                        "description": "Field data to export (from field computation results)"
                    },
                    "output_file": {
                        "type": "string",
                        "description": "Output JSON file path"
                    }
                },
                "required": ["field_data", "output_file"]
            }
        ),
    ]

    # Add workspace tools if available
    if WORKSPACE_AVAILABLE:
        tools.extend([
            Tool(
                name="radia_workspace_create_session",
                description="Create a new shared workspace session for Radia-NGSolve coupling.",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "session_name": {
                            "type": "string",
                            "description": "Human-readable session name (optional)"
                        }
                    },
                    "required": []
                }
            ),
            Tool(
                name="radia_workspace_export_object",
                description="Export Radia object to shared workspace for NGSolve server access.",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "object_name": {
                            "type": "string",
                            "description": "Name of Radia object to export"
                        },
                        "session_id": {
                            "type": "string",
                            "description": "Session ID (if not provided, creates new session)"
                        },
                        "export_geometry": {
                            "type": "boolean",
                            "description": "Export geometry to VTK",
                            "default": True
                        },
                        "export_fields": {
                            "type": "boolean",
                            "description": "Pre-compute and export field data on a grid",
                            "default": False
                        },
                        "field_grid_size": {
                            "type": "integer",
                            "description": "Grid points per dimension for field export",
                            "default": 20
                        }
                    },
                    "required": ["object_name"]
                }
            ),
            Tool(
                name="radia_workspace_list_sessions",
                description="List all available workspace sessions.",
                inputSchema={
                    "type": "object",
                    "properties": {},
                    "required": []
                }
            ),
            Tool(
                name="radia_workspace_get_session_info",
                description="Get detailed information about a workspace session.",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "session_id": {
                            "type": "string",
                            "description": "Session ID"
                        }
                    },
                    "required": ["session_id"]
                }
            ),
        ])

    return tools


async def execute(name: str, arguments: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Execute an export or workspace tool."""
    if rad is None:
        return {"error": "Radia module not available. Please build Radia first."}

    try:
        if name == "radia_export_vtk_legacy":
            return _export_vtk_legacy(arguments, state)
        elif name == "radia_export_field_data_json":
            return _export_field_json(arguments, state)
        elif name == "radia_workspace_create_session":
            return _workspace_create_session(arguments, state)
        elif name == "radia_workspace_export_object":
            return _workspace_export_object(arguments, state)
        elif name == "radia_workspace_list_sessions":
            return _workspace_list_sessions(arguments, state)
        elif name == "radia_workspace_get_session_info":
            return _workspace_get_session_info(arguments, state)
        else:
            return {"error": f"Unknown export tool: {name}"}
    except Exception as e:
        return {"error": str(e), "tool": name}


def _export_vtk_legacy(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Export geometry using built-in ObjDrwVTK."""
    obj_name = args["object_name"]
    if obj_name not in state:
        return {"error": f"Object '{obj_name}' not found in state"}

    obj = state[obj_name]
    output_file = args["output_file"]

    rad.ObjDrwVTK(obj, output_file)

    return {
        "success": True,
        "object_name": obj_name,
        "output_file": output_file,
        "method": "ObjDrwVTK",
        "message": f"Geometry exported to {output_file}"
    }


def _export_field_json(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Export field data to JSON."""
    field_data = args["field_data"]
    output_file = args["output_file"]

    with open(output_file, 'w') as f:
        json.dump(field_data, f, indent=2)

    return {
        "success": True,
        "output_file": output_file,
        "num_records": len(field_data.get("field_values", [])),
        "message": f"Field data exported to {output_file}"
    }


def _workspace_create_session(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create a new workspace session."""
    if not WORKSPACE_AVAILABLE:
        return {"error": "Shared workspace not available"}

    workspace = SharedWorkspace()
    session_name = args.get("session_name")
    session_id = workspace.create_session(session_name)

    # Store session ID in server state
    if "workspace_sessions" not in state:
        state["workspace_sessions"] = []
    state["workspace_sessions"].append(session_id)
    state["current_session"] = session_id

    return {
        "success": True,
        "session_id": session_id,
        "session_name": session_name or f"Session {session_id[:8]}",
        "workspace_location": str(workspace.base_path),
        "message": f"Created workspace session: {session_id}"
    }


def _workspace_export_object(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Export Radia object to shared workspace."""
    if not WORKSPACE_AVAILABLE:
        return {"error": "Shared workspace not available"}

    workspace = SharedWorkspace()
    obj_name = args["object_name"]

    if obj_name not in state:
        return {"error": f"Object '{obj_name}' not found in state"}

    # Get or create session
    session_id = args.get("session_id")
    if not session_id:
        session_id = workspace.create_session()
        state["current_session"] = session_id

    radia_obj = state[obj_name]

    # Export geometry to VTK
    vtk_file = None
    if args.get("export_geometry", True):
        with tempfile.NamedTemporaryFile(suffix=".vtk", delete=False) as f:
            vtk_file = f.name
        rad.ObjDrwVTK(radia_obj, vtk_file)

    # Export object metadata
    obj_file = workspace.export_radia_object(
        session_id,
        obj_name,
        radia_obj,
        vtk_file
    )

    # Pre-compute fields if requested
    field_file = None
    if args.get("export_fields", False):
        import numpy as np
        grid_size = args.get("field_grid_size", 20)

        # Compute field on a grid
        x = np.linspace(-0.2, 0.2, grid_size)
        y = np.linspace(-0.2, 0.2, grid_size)
        z = np.linspace(-0.2, 0.2, grid_size)

        points = []
        field_values = []

        for zi in z:
            for yi in y:
                for xi in x:
                    pt = [float(xi), float(yi), float(zi)]
                    points.append(pt)
                    try:
                        field = rad.Fld(radia_obj, 'b', pt)
                        field_values.append([float(f) for f in field])
                    except:
                        field_values.append([0.0, 0.0, 0.0])

        field_file = workspace.export_field_data(
            session_id,
            obj_name,
            points,
            field_values,
            field_type='b'
        )

    return {
        "success": True,
        "session_id": session_id,
        "object_name": obj_name,
        "exported_files": {
            "metadata": obj_file,
            "geometry": vtk_file if vtk_file else None,
            "fields": field_file if field_file else None
        },
        "message": f"Object '{obj_name}' exported to workspace session {session_id}"
    }


def _workspace_list_sessions(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """List all workspace sessions."""
    if not WORKSPACE_AVAILABLE:
        return {"error": "Shared workspace not available"}

    workspace = SharedWorkspace()
    sessions = workspace.list_sessions()

    return {
        "success": True,
        "num_sessions": len(sessions),
        "sessions": sessions,
        "workspace_location": str(workspace.base_path)
    }


def _workspace_get_session_info(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Get session information."""
    if not WORKSPACE_AVAILABLE:
        return {"error": "Shared workspace not available"}

    workspace = SharedWorkspace()
    session_id = args["session_id"]

    try:
        info = workspace.get_session_info(session_id)
        return {
            "success": True,
            "session_info": info
        }
    except Exception as e:
        return {"error": f"Session '{session_id}' not found: {str(e)}"}
