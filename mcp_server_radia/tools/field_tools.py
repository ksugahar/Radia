"""
Field Calculation Tools for Radia MCP Server

Tools for computing magnetic fields and potentials:
- B (magnetic flux density)
- H (magnetic field intensity)
- A (vector potential)
- Phi (scalar potential)
"""

import sys
from pathlib import Path
from typing import Any, Dict, List

# Add Radia to path
radia_path = Path(__file__).parent.parent.parent / "build" / "Release"
if radia_path.exists() and str(radia_path) not in sys.path:
    sys.path.insert(0, str(radia_path))

try:
    import radia as rad
except ImportError:
    rad = None

from mcp.types import Tool


def get_tools() -> List[Tool]:
    """Get list of field calculation tools."""
    return [
        Tool(
            name="radia_field_compute",
            description="Compute magnetic field or potential at a point using rad.Fld. "
                       "Field types: 'b' (flux density in T), 'h' (field intensity in A/m), "
                       "'a' (vector potential in T*m), 'phi' (scalar potential).",
            inputSchema={
                "type": "object",
                "properties": {
                    "object_name": {
                        "type": "string",
                        "description": "Name of the object to compute field from"
                    },
                    "field_type": {
                        "type": "string",
                        "enum": ["b", "h", "a", "phi"],
                        "description": "Type of field to compute"
                    },
                    "point": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Evaluation point [x, y, z] in current units"
                    }
                },
                "required": ["object_name", "field_type", "point"]
            }
        ),
        Tool(
            name="radia_field_compute_line",
            description="Compute field along a line with multiple points.",
            inputSchema={
                "type": "object",
                "properties": {
                    "object_name": {
                        "type": "string",
                        "description": "Name of the object to compute field from"
                    },
                    "field_type": {
                        "type": "string",
                        "enum": ["b", "h", "a", "phi"],
                        "description": "Type of field to compute"
                    },
                    "points": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": {"type": "number"},
                            "minItems": 3,
                            "maxItems": 3
                        },
                        "description": "List of evaluation points [[x,y,z], ...]"
                    }
                },
                "required": ["object_name", "field_type", "points"]
            }
        ),
        Tool(
            name="radia_field_compute_grid",
            description="Compute field on a regular 3D grid.",
            inputSchema={
                "type": "object",
                "properties": {
                    "object_name": {
                        "type": "string",
                        "description": "Name of the object to compute field from"
                    },
                    "field_type": {
                        "type": "string",
                        "enum": ["b", "h", "a", "phi"],
                        "description": "Type of field to compute"
                    },
                    "x_range": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "[x_min, x_max, num_points]"
                    },
                    "y_range": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "[y_min, y_max, num_points]"
                    },
                    "z_range": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "[z_min, z_max, num_points]"
                    }
                },
                "required": ["object_name", "field_type", "x_range", "y_range", "z_range"]
            }
        ),
        Tool(
            name="radia_field_force",
            description="Compute force on an object using rad.FldFrc.",
            inputSchema={
                "type": "object",
                "properties": {
                    "object_name": {
                        "type": "string",
                        "description": "Name of the object to compute force on"
                    },
                    "source_name": {
                        "type": "string",
                        "description": "Name of the source object creating the field (optional)",
                        "default": None
                    }
                },
                "required": ["object_name"]
            }
        ),
        Tool(
            name="radia_field_torque",
            description="Compute torque on an object using rad.FldTrq.",
            inputSchema={
                "type": "object",
                "properties": {
                    "object_name": {
                        "type": "string",
                        "description": "Name of the object to compute torque on"
                    },
                    "pivot_point": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Pivot point [x, y, z] for torque calculation"
                    },
                    "source_name": {
                        "type": "string",
                        "description": "Name of the source object (optional)",
                        "default": None
                    }
                },
                "required": ["object_name", "pivot_point"]
            }
        ),
    ]


async def execute(name: str, arguments: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Execute a field calculation tool."""
    if rad is None:
        return {"error": "Radia module not available. Please build Radia first."}

    try:
        if name == "radia_field_compute":
            return _compute_field(arguments, state)
        elif name == "radia_field_compute_line":
            return _compute_field_line(arguments, state)
        elif name == "radia_field_compute_grid":
            return _compute_field_grid(arguments, state)
        elif name == "radia_field_force":
            return _compute_force(arguments, state)
        elif name == "radia_field_torque":
            return _compute_torque(arguments, state)
        else:
            return {"error": f"Unknown field tool: {name}"}
    except Exception as e:
        return {"error": str(e), "tool": name}


def _compute_field(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Compute field at a single point."""
    obj_name = args["object_name"]
    if obj_name not in state:
        return {"error": f"Object '{obj_name}' not found in state"}

    obj = state[obj_name]
    field_type = args["field_type"]
    point = args["point"]

    field = rad.Fld(obj, field_type, point)

    return {
        "success": True,
        "object_name": obj_name,
        "field_type": field_type,
        "point": point,
        "field_value": field,
        "units": {
            "b": "Tesla",
            "h": "A/m",
            "a": "T*m",
            "phi": "dimensionless"
        }[field_type]
    }


def _compute_field_line(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Compute field along a line."""
    obj_name = args["object_name"]
    if obj_name not in state:
        return {"error": f"Object '{obj_name}' not found in state"}

    obj = state[obj_name]
    field_type = args["field_type"]
    points = args["points"]

    fields = []
    for point in points:
        field = rad.Fld(obj, field_type, point)
        fields.append(field)

    return {
        "success": True,
        "object_name": obj_name,
        "field_type": field_type,
        "num_points": len(points),
        "points": points,
        "field_values": fields
    }


def _compute_field_grid(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Compute field on a regular grid."""
    import numpy as np

    obj_name = args["object_name"]
    if obj_name not in state:
        return {"error": f"Object '{obj_name}' not found in state"}

    obj = state[obj_name]
    field_type = args["field_type"]

    x_range = args["x_range"]
    y_range = args["y_range"]
    z_range = args["z_range"]

    x = np.linspace(x_range[0], x_range[1], int(x_range[2]))
    y = np.linspace(y_range[0], y_range[1], int(y_range[2]))
    z = np.linspace(z_range[0], z_range[1], int(z_range[2]))

    fields = []
    grid_points = []

    for zi in z:
        for yi in y:
            for xi in x:
                point = [xi, yi, zi]
                grid_points.append(point)
                field = rad.Fld(obj, field_type, point)
                fields.append(field)

    return {
        "success": True,
        "object_name": obj_name,
        "field_type": field_type,
        "grid_shape": [int(x_range[2]), int(y_range[2]), int(z_range[2])],
        "num_points": len(fields),
        "field_values": fields,
        "note": "Field values are in flattened (Z-Y-X) order"
    }


def _compute_force(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Compute force on an object."""
    obj_name = args["object_name"]
    if obj_name not in state:
        return {"error": f"Object '{obj_name}' not found in state"}

    obj = state[obj_name]
    source_name = args.get("source_name")

    if source_name:
        if source_name not in state:
            return {"error": f"Source object '{source_name}' not found in state"}
        source = state[source_name]
        force = rad.FldFrc(obj, source)
    else:
        force = rad.FldFrc(obj)

    return {
        "success": True,
        "object_name": obj_name,
        "source_name": source_name,
        "force": force,
        "units": "Newton"
    }


def _compute_torque(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Compute torque on an object."""
    obj_name = args["object_name"]
    if obj_name not in state:
        return {"error": f"Object '{obj_name}' not found in state"}

    obj = state[obj_name]
    pivot = args["pivot_point"]
    source_name = args.get("source_name")

    if source_name:
        if source_name not in state:
            return {"error": f"Source object '{source_name}' not found in state"}
        source = state[source_name]
        torque = rad.FldTrq(obj, source, pivot)
    else:
        torque = rad.FldTrq(obj, pivot)

    return {
        "success": True,
        "object_name": obj_name,
        "pivot_point": pivot,
        "source_name": source_name,
        "torque": torque,
        "units": "N*m"
    }
