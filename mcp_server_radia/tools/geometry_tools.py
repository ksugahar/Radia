"""
Geometry Tools for Radia MCP Server

Tools for creating magnetic geometry objects:
- Rectangular magnets
- Hexahedral elements
- Tetrahedral elements
- Coils and conductors
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
    """Get list of geometry creation tools."""
    return [
        Tool(
            name="radia_geometry_create_recmag",
            description="Create a rectangular permanent magnet block using ObjRecMag. "
                       "Supports fixed magnetization for permanent magnets.",
            inputSchema={
                "type": "object",
                "properties": {
                    "center": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Center position [x, y, z] in current units"
                    },
                    "dimensions": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Dimensions [width, height, depth] in current units"
                    },
                    "magnetization": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Magnetization vector [Mx, My, Mz] in A/m"
                    },
                    "name": {
                        "type": "string",
                        "description": "Object name for reference in state",
                        "default": "recmag"
                    }
                },
                "required": ["center", "dimensions", "magnetization"]
            }
        ),
        Tool(
            name="radia_geometry_create_hexahedron",
            description="Create a hexahedral element using ObjHexahedron. "
                       "Requires 8 vertices defining the hexahedron corners.",
            inputSchema={
                "type": "object",
                "properties": {
                    "vertices": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": {"type": "number"},
                            "minItems": 3,
                            "maxItems": 3
                        },
                        "minItems": 8,
                        "maxItems": 8,
                        "description": "8 vertices [[x,y,z], ...] defining the hexahedron"
                    },
                    "magnetization": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Magnetization vector [Mx, My, Mz] in A/m"
                    },
                    "name": {
                        "type": "string",
                        "description": "Object name for reference in state",
                        "default": "hexahedron"
                    }
                },
                "required": ["vertices", "magnetization"]
            }
        ),
        Tool(
            name="radia_geometry_create_tetrahedron",
            description="Create a tetrahedral element using ObjTetrahedron. "
                       "Requires 4 vertices defining the tetrahedron corners.",
            inputSchema={
                "type": "object",
                "properties": {
                    "vertices": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": {"type": "number"},
                            "minItems": 3,
                            "maxItems": 3
                        },
                        "minItems": 4,
                        "maxItems": 4,
                        "description": "4 vertices [[x,y,z], ...] defining the tetrahedron"
                    },
                    "magnetization": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Magnetization vector [Mx, My, Mz] in A/m"
                    },
                    "name": {
                        "type": "string",
                        "description": "Object name for reference in state",
                        "default": "tetrahedron"
                    }
                },
                "required": ["vertices", "magnetization"]
            }
        ),
        Tool(
            name="radia_geometry_create_racetrack_coil",
            description="Create a race-track coil using ObjRaceTrk.",
            inputSchema={
                "type": "object",
                "properties": {
                    "center": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Center position [x, y, z]"
                    },
                    "radii": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 2,
                        "maxItems": 2,
                        "description": "Inner radii [R_min, R_max]"
                    },
                    "lengths": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 2,
                        "maxItems": 2,
                        "description": "Straight section lengths [Lx, Ly]"
                    },
                    "height": {
                        "type": "number",
                        "description": "Coil height"
                    },
                    "current": {
                        "type": "number",
                        "description": "Current in Amperes"
                    },
                    "name": {
                        "type": "string",
                        "description": "Object name for reference in state",
                        "default": "racetrack_coil"
                    }
                },
                "required": ["center", "radii", "lengths", "height", "current"]
            }
        ),
        Tool(
            name="radia_geometry_create_container",
            description="Create a container (ObjCnt) to group multiple objects together.",
            inputSchema={
                "type": "object",
                "properties": {
                    "object_names": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Names of objects to include in container"
                    },
                    "name": {
                        "type": "string",
                        "description": "Container name for reference in state",
                        "default": "container"
                    }
                },
                "required": ["object_names"]
            }
        ),
        Tool(
            name="radia_geometry_set_units",
            description="Set field units for Radia. Use 'm' for NGSolve integration.",
            inputSchema={
                "type": "object",
                "properties": {
                    "units": {
                        "type": "string",
                        "enum": ["m", "mm"],
                        "description": "Length units: 'm' (meters) or 'mm' (millimeters)"
                    }
                },
                "required": ["units"]
            }
        ),
    ]


async def execute(name: str, arguments: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Execute a geometry tool."""
    if rad is None:
        return {"error": "Radia module not available. Please build Radia first."}

    try:
        if name == "radia_geometry_create_recmag":
            return _create_recmag(arguments, state)
        elif name == "radia_geometry_create_hexahedron":
            return _create_hexahedron(arguments, state)
        elif name == "radia_geometry_create_tetrahedron":
            return _create_tetrahedron(arguments, state)
        elif name == "radia_geometry_create_racetrack_coil":
            return _create_racetrack_coil(arguments, state)
        elif name == "radia_geometry_create_container":
            return _create_container(arguments, state)
        elif name == "radia_geometry_set_units":
            return _set_units(arguments, state)
        else:
            return {"error": f"Unknown geometry tool: {name}"}
    except Exception as e:
        return {"error": str(e), "tool": name}


def _create_recmag(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create a rectangular magnet."""
    obj = rad.ObjRecMag(args["center"], args["dimensions"], args["magnetization"])
    obj_name = args.get("name", "recmag")
    state[obj_name] = obj

    return {
        "success": True,
        "object_name": obj_name,
        "object_id": obj,
        "type": "rectangular_magnet",
        "center": args["center"],
        "dimensions": args["dimensions"],
        "magnetization": args["magnetization"]
    }


def _create_hexahedron(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create a hexahedral element."""
    obj = rad.ObjHexahedron(args["vertices"], args["magnetization"])
    obj_name = args.get("name", "hexahedron")
    state[obj_name] = obj

    return {
        "success": True,
        "object_name": obj_name,
        "object_id": obj,
        "type": "hexahedron",
        "num_vertices": len(args["vertices"]),
        "magnetization": args["magnetization"]
    }


def _create_tetrahedron(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create a tetrahedral element."""
    obj = rad.ObjTetrahedron(args["vertices"], args["magnetization"])
    obj_name = args.get("name", "tetrahedron")
    state[obj_name] = obj

    return {
        "success": True,
        "object_name": obj_name,
        "object_id": obj,
        "type": "tetrahedron",
        "num_vertices": len(args["vertices"]),
        "magnetization": args["magnetization"]
    }


def _create_racetrack_coil(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create a race-track coil."""
    obj = rad.ObjRaceTrk(
        args["center"],
        args["radii"],
        args["lengths"],
        args["height"],
        3.0,  # Default curvature radius
        args["current"],
        'man'  # Manual mode
    )
    obj_name = args.get("name", "racetrack_coil")
    state[obj_name] = obj

    return {
        "success": True,
        "object_name": obj_name,
        "object_id": obj,
        "type": "racetrack_coil",
        "center": args["center"],
        "current": args["current"]
    }


def _create_container(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create a container from multiple objects."""
    obj_names = args["object_names"]
    objects = []

    for obj_name in obj_names:
        if obj_name not in state:
            return {"error": f"Object '{obj_name}' not found in state"}
        objects.append(state[obj_name])

    container = rad.ObjCnt(objects)
    container_name = args.get("name", "container")
    state[container_name] = container

    return {
        "success": True,
        "object_name": container_name,
        "object_id": container,
        "type": "container",
        "num_objects": len(objects),
        "contained_objects": obj_names
    }


def _set_units(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Set Radia field units."""
    units = args["units"]
    result = rad.FldUnits(units)

    return {
        "success": True,
        "units": units,
        "radia_response": result
    }
