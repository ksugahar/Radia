"""
NGSolve Integration Tools for Radia MCP Server

Tools for creating NGSolve-compatible field representations from Radia objects:
- RadiaField CoefficientFunction creation
- Batch field evaluation for efficiency
- H-matrix acceleration support
- Field data export for NGSolve coupling
"""

import sys
from typing import Any, Dict, List

try:
    import radia as rad
    import numpy as np
    RADIA_AVAILABLE = True
except ImportError:
    RADIA_AVAILABLE = False

from mcp.types import Tool


def get_tools() -> List[Tool]:
    """Get list of NGSolve integration tools."""
    if not RADIA_AVAILABLE:
        return []

    return [
        Tool(
            name="radia_ngsolve_create_field",
            description="Create field data for NGSolve CoefficientFunction integration. Returns field values on a grid.",
            inputSchema={
                "type": "object",
                "properties": {
                    "radia_object_name": {
                        "type": "string",
                        "description": "Name of the Radia object"
                    },
                    "field_type": {
                        "type": "string",
                        "enum": ["b", "h", "a", "m"],
                        "description": "Field type: 'b' (B field), 'h' (H field), 'a' (vector potential), 'm' (magnetization)",
                        "default": "b"
                    },
                    "grid_points": {
                        "type": "array",
                        "description": "List of [x,y,z] points where to evaluate field",
                        "items": {
                            "type": "array",
                            "items": {"type": "number"},
                            "minItems": 3,
                            "maxItems": 3
                        }
                    },
                    "use_hmatrix": {
                        "type": "boolean",
                        "description": "Use H-matrix acceleration for large problems",
                        "default": False
                    }
                },
                "required": ["radia_object_name", "grid_points"]
            }
        ),
        Tool(
            name="radia_ngsolve_batch_evaluate",
            description="Batch field evaluation using H-matrix acceleration for efficient NGSolve integration.",
            inputSchema={
                "type": "object",
                "properties": {
                    "radia_object_name": {
                        "type": "string",
                        "description": "Name of the Radia object"
                    },
                    "field_type": {
                        "type": "string",
                        "enum": ["b", "h", "a", "m"],
                        "description": "Field type to evaluate",
                        "default": "b"
                    },
                    "points": {
                        "type": "array",
                        "description": "Array of [x,y,z] evaluation points",
                        "items": {
                            "type": "array",
                            "items": {"type": "number"},
                            "minItems": 3,
                            "maxItems": 3
                        }
                    },
                    "use_hmatrix": {
                        "type": "boolean",
                        "description": "Enable H-matrix acceleration (recommended for >100 points)",
                        "default": True
                    },
                    "hmatrix_precision": {
                        "type": "number",
                        "description": "H-matrix precision (typically 1e-6)",
                        "default": 1e-6
                    }
                },
                "required": ["radia_object_name", "points"]
            }
        ),
        Tool(
            name="radia_ngsolve_enable_hmatrix",
            description="Enable H-matrix acceleration for field evaluation (speeds up batch evaluation).",
            inputSchema={
                "type": "object",
                "properties": {
                    "enable": {
                        "type": "boolean",
                        "description": "Enable (true) or disable (false) H-matrix",
                        "default": True
                    },
                    "precision": {
                        "type": "number",
                        "description": "Precision for H-matrix approximation",
                        "default": 1e-6
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="radia_ngsolve_export_field_grid",
            description="Export field values on a structured grid for NGSolve import.",
            inputSchema={
                "type": "object",
                "properties": {
                    "radia_object_name": {
                        "type": "string",
                        "description": "Name of the Radia object"
                    },
                    "field_type": {
                        "type": "string",
                        "enum": ["b", "h", "a"],
                        "description": "Field type",
                        "default": "b"
                    },
                    "bbox_min": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Minimum corner [x,y,z] of evaluation box"
                    },
                    "bbox_max": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Maximum corner [x,y,z] of evaluation box"
                    },
                    "grid_size": {
                        "type": "array",
                        "items": {"type": "integer"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Number of points in [x,y,z] directions",
                        "default": [10, 10, 10]
                    },
                    "output_file": {
                        "type": "string",
                        "description": "Output NPZ file path for field data"
                    }
                },
                "required": ["radia_object_name", "bbox_min", "bbox_max", "output_file"]
            }
        ),
        Tool(
            name="radia_ngsolve_mesh_field_sample",
            description="Sample Radia field values at NGSolve mesh vertices (requires mesh vertex coordinates).",
            inputSchema={
                "type": "object",
                "properties": {
                    "radia_object_name": {
                        "type": "string",
                        "description": "Name of the Radia object"
                    },
                    "field_type": {
                        "type": "string",
                        "enum": ["b", "h", "a"],
                        "description": "Field type",
                        "default": "b"
                    },
                    "mesh_vertices": {
                        "type": "array",
                        "description": "Array of mesh vertex coordinates [[x,y,z], ...]",
                        "items": {
                            "type": "array",
                            "items": {"type": "number"},
                            "minItems": 3,
                            "maxItems": 3
                        }
                    }
                },
                "required": ["radia_object_name", "mesh_vertices"]
            }
        ),
    ]


async def execute(name: str, arguments: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Execute an NGSolve integration tool."""
    if not RADIA_AVAILABLE:
        return {
            "error": "Radia not available",
            "message": "Build Radia first using BuildMSVC.ps1"
        }

    try:
        if name == "radia_ngsolve_create_field":
            return _create_field(arguments, state)
        elif name == "radia_ngsolve_batch_evaluate":
            return _batch_evaluate(arguments, state)
        elif name == "radia_ngsolve_enable_hmatrix":
            return _enable_hmatrix(arguments, state)
        elif name == "radia_ngsolve_export_field_grid":
            return _export_field_grid(arguments, state)
        elif name == "radia_ngsolve_mesh_field_sample":
            return _mesh_field_sample(arguments, state)
        else:
            return {"error": f"Unknown NGSolve integration tool: {name}"}
    except Exception as e:
        return {"error": str(e), "tool": name, "traceback": __import__('traceback').format_exc()}


def _create_field(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create field data for NGSolve integration."""
    obj_name = args["radia_object_name"]
    field_type = args.get("field_type", "b")
    grid_points = args["grid_points"]
    use_hmatrix = args.get("use_hmatrix", False)

    if obj_name not in state:
        return {"error": f"Radia object '{obj_name}' not found"}

    radia_obj = state[obj_name]

    # Evaluate field at grid points
    if use_hmatrix and len(grid_points) > 100:
        # Use batch evaluation with H-matrix
        field_values = rad.FldBatch(radia_obj, field_type, grid_points, 1)
    else:
        # Standard evaluation
        field_values = [rad.Fld(radia_obj, field_type, pt) for pt in grid_points]

    field_values = np.array(field_values)

    return {
        "success": True,
        "field_type": field_type,
        "num_points": len(grid_points),
        "field_values": field_values.tolist(),
        "statistics": {
            "min": field_values.min(axis=0).tolist(),
            "max": field_values.max(axis=0).tolist(),
            "mean": field_values.mean(axis=0).tolist()
        },
        "message": f"Field evaluated at {len(grid_points)} points"
    }


def _batch_evaluate(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Batch field evaluation with H-matrix acceleration."""
    import time

    obj_name = args["radia_object_name"]
    field_type = args.get("field_type", "b")
    points = args["points"]
    use_hmatrix = args.get("use_hmatrix", True)
    hmatrix_precision = args.get("hmatrix_precision", 1e-6)

    if obj_name not in state:
        return {"error": f"Radia object '{obj_name}' not found"}

    radia_obj = state[obj_name]

    # Enable H-matrix if requested
    if use_hmatrix:
        rad.SetHMatrixFieldEval(1, hmatrix_precision)

    # Batch evaluation
    t0 = time.time()
    field_values = rad.FldBatch(radia_obj, field_type, points, int(use_hmatrix))
    eval_time = time.time() - t0

    field_values = np.array(field_values)

    return {
        "success": True,
        "field_type": field_type,
        "num_points": len(points),
        "evaluation_time_ms": eval_time * 1000,
        "points_per_second": len(points) / eval_time if eval_time > 0 else 0,
        "hmatrix_enabled": use_hmatrix,
        "field_values": field_values.tolist(),
        "statistics": {
            "min": field_values.min(axis=0).tolist(),
            "max": field_values.max(axis=0).tolist(),
            "mean": field_values.mean(axis=0).tolist(),
            "std": field_values.std(axis=0).tolist()
        }
    }


def _enable_hmatrix(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Enable or disable H-matrix acceleration."""
    enable = args.get("enable", True)
    precision = args.get("precision", 1e-6)

    if enable:
        rad.SetHMatrixFieldEval(1, precision)
        return {
            "success": True,
            "hmatrix_enabled": True,
            "precision": precision,
            "message": f"H-matrix enabled with precision {precision}"
        }
    else:
        rad.SetHMatrixFieldEval(0, 0)
        return {
            "success": True,
            "hmatrix_enabled": False,
            "message": "H-matrix disabled"
        }


def _export_field_grid(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Export field on structured grid to NPZ file."""
    obj_name = args["radia_object_name"]
    field_type = args.get("field_type", "b")
    bbox_min = np.array(args["bbox_min"])
    bbox_max = np.array(args["bbox_max"])
    grid_size = args.get("grid_size", [10, 10, 10])
    output_file = args["output_file"]

    if obj_name not in state:
        return {"error": f"Radia object '{obj_name}' not found"}

    radia_obj = state[obj_name]

    # Create structured grid
    x = np.linspace(bbox_min[0], bbox_max[0], grid_size[0])
    y = np.linspace(bbox_min[1], bbox_max[1], grid_size[1])
    z = np.linspace(bbox_min[2], bbox_max[2], grid_size[2])

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()]).tolist()

    # Batch evaluate
    rad.SetHMatrixFieldEval(1, 1e-6)
    field_values = rad.FldBatch(radia_obj, field_type, points, 1)
    field_values = np.array(field_values)

    # Reshape to grid
    field_grid = field_values.reshape(grid_size[0], grid_size[1], grid_size[2], 3)

    # Save to NPZ
    np.savez(output_file,
             x=x, y=y, z=z,
             field_type=field_type,
             field_values=field_grid,
             bbox_min=bbox_min,
             bbox_max=bbox_max)

    return {
        "success": True,
        "output_file": output_file,
        "field_type": field_type,
        "grid_size": grid_size,
        "total_points": len(points),
        "bbox": {
            "min": bbox_min.tolist(),
            "max": bbox_max.tolist()
        },
        "statistics": {
            "min": field_values.min(axis=0).tolist(),
            "max": field_values.max(axis=0).tolist(),
            "mean": field_values.mean(axis=0).tolist()
        },
        "message": f"Field grid exported to {output_file}"
    }


def _mesh_field_sample(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Sample field at NGSolve mesh vertices."""
    obj_name = args["radia_object_name"]
    field_type = args.get("field_type", "b")
    mesh_vertices = args["mesh_vertices"]

    if obj_name not in state:
        return {"error": f"Radia object '{obj_name}' not found"}

    radia_obj = state[obj_name]

    # Use H-matrix for efficiency
    rad.SetHMatrixFieldEval(1, 1e-6)
    field_values = rad.FldBatch(radia_obj, field_type, mesh_vertices, 1)
    field_values = np.array(field_values)

    return {
        "success": True,
        "field_type": field_type,
        "num_vertices": len(mesh_vertices),
        "field_values": field_values.tolist(),
        "statistics": {
            "min": field_values.min(axis=0).tolist(),
            "max": field_values.max(axis=0).tolist(),
            "mean": field_values.mean(axis=0).tolist(),
            "std": field_values.std(axis=0).tolist()
        },
        "message": f"Field sampled at {len(mesh_vertices)} mesh vertices"
    }
