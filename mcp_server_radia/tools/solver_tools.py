"""
Solver Tools for Radia MCP Server

Tools for solving magnetostatic problems:
- Solve with various methods (LU, BiCGSTAB, HACApK)
- Solver parameter configuration
- Background fields
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
    """Get list of solver tools."""
    return [
        Tool(
            name="radia_solver_solve",
            description="Solve the magnetostatic problem using rad.Solve. "
                       "Method 0=LU (direct), 1=BiCGSTAB (iterative), 2=HACApK (H-matrix accelerated).",
            inputSchema={
                "type": "object",
                "properties": {
                    "object_name": {
                        "type": "string",
                        "description": "Name of the object or container to solve"
                    },
                    "precision": {
                        "type": "number",
                        "description": "Convergence tolerance",
                        "default": 0.0001
                    },
                    "max_iterations": {
                        "type": "integer",
                        "description": "Maximum number of iterations",
                        "default": 1000
                    },
                    "method": {
                        "type": "integer",
                        "enum": [0, 1, 2],
                        "description": "Solver method: 0=LU, 1=BiCGSTAB, 2=HACApK",
                        "default": 1
                    }
                },
                "required": ["object_name"]
            }
        ),
        Tool(
            name="radia_solver_set_hacapk_params",
            description="Set H-matrix (HACApK) solver parameters for Method 2.",
            inputSchema={
                "type": "object",
                "properties": {
                    "eps": {
                        "type": "number",
                        "description": "ACA tolerance (1e-6 to 1e-2, default 1e-4)",
                        "default": 0.0001
                    },
                    "leaf_size": {
                        "type": "integer",
                        "description": "Minimum cluster size (5-50, default 10)",
                        "default": 10
                    },
                    "eta": {
                        "type": "number",
                        "description": "Admissibility parameter (1.0-3.0, default 2.0)",
                        "default": 2.0
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="radia_solver_set_relax_param",
            description="Set under-relaxation parameter for difficult nonlinear problems (0.0-1.0).",
            inputSchema={
                "type": "object",
                "properties": {
                    "relax_param": {
                        "type": "number",
                        "description": "Relaxation parameter (0.0 = full step, higher = more damping)",
                        "minimum": 0.0,
                        "maximum": 1.0,
                        "default": 0.0
                    }
                },
                "required": ["relax_param"]
            }
        ),
        Tool(
            name="radia_solver_create_background_field",
            description="Create a background field using ObjBckg with a Python callback.",
            inputSchema={
                "type": "object",
                "properties": {
                    "field_type": {
                        "type": "string",
                        "enum": ["uniform", "gradient", "quadrupole"],
                        "description": "Type of background field"
                    },
                    "parameters": {
                        "type": "object",
                        "description": "Field-specific parameters",
                        "properties": {
                            "bx": {"type": "number", "description": "X component (Tesla) for uniform"},
                            "by": {"type": "number", "description": "Y component (Tesla) for uniform"},
                            "bz": {"type": "number", "description": "Z component (Tesla) for uniform"},
                            "gradient": {"type": "number", "description": "Gradient (T/m) for quadrupole"}
                        }
                    },
                    "name": {
                        "type": "string",
                        "description": "Background field name",
                        "default": "background_field"
                    }
                },
                "required": ["field_type", "parameters"]
            }
        ),
        Tool(
            name="radia_solver_get_hacapk_stats",
            description="Get H-matrix statistics after solving with Method 2.",
            inputSchema={
                "type": "object",
                "properties": {},
                "required": []
            }
        ),
    ]


async def execute(name: str, arguments: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Execute a solver tool."""
    if rad is None:
        return {"error": "Radia module not available. Please build Radia first."}

    try:
        if name == "radia_solver_solve":
            return _solve(arguments, state)
        elif name == "radia_solver_set_hacapk_params":
            return _set_hacapk_params(arguments, state)
        elif name == "radia_solver_set_relax_param":
            return _set_relax_param(arguments, state)
        elif name == "radia_solver_create_background_field":
            return _create_background_field(arguments, state)
        elif name == "radia_solver_get_hacapk_stats":
            return _get_hacapk_stats(arguments, state)
        else:
            return {"error": f"Unknown solver tool: {name}"}
    except Exception as e:
        return {"error": str(e), "tool": name}


def _solve(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Solve the magnetostatic problem."""
    obj_name = args["object_name"]
    if obj_name not in state:
        return {"error": f"Object '{obj_name}' not found in state"}

    obj = state[obj_name]
    precision = args.get("precision", 0.0001)
    max_iter = args.get("max_iterations", 1000)
    method = args.get("method", 1)

    result = rad.Solve(obj, precision, max_iter, method)

    return {
        "success": True,
        "object_name": obj_name,
        "solver_result": result,
        "precision": precision,
        "max_iterations": max_iter,
        "method": method,
        "method_name": {0: "LU", 1: "BiCGSTAB", 2: "HACApK"}[method]
    }


def _set_hacapk_params(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Set HACApK solver parameters."""
    eps = args.get("eps", 0.0001)
    leaf_size = args.get("leaf_size", 10)
    eta = args.get("eta", 2.0)

    rad.SetHACApKParams(eps, leaf_size, eta)

    return {
        "success": True,
        "eps": eps,
        "leaf_size": leaf_size,
        "eta": eta,
        "message": "HACApK parameters set successfully"
    }


def _set_relax_param(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Set relaxation parameter."""
    relax = args["relax_param"]
    rad.SetRelaxParam(relax)

    return {
        "success": True,
        "relax_param": relax,
        "message": f"Relaxation parameter set to {relax}"
    }


def _create_background_field(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create a background field."""
    field_type = args["field_type"]
    params = args["parameters"]

    if field_type == "uniform":
        bx = params.get("bx", 0.0)
        by = params.get("by", 0.0)
        bz = params.get("bz", 0.0)
        bkg = rad.ObjBckg(lambda p: [bx, by, bz])
        field_desc = f"Uniform field [{bx}, {by}, {bz}] T"

    elif field_type == "quadrupole":
        G = params.get("gradient", 10.0)  # T/m
        def quadrupole_field(point):
            x, y, z = point
            return [G * y, G * x, 0]
        bkg = rad.ObjBckg(quadrupole_field)
        field_desc = f"Quadrupole field with gradient {G} T/m"

    else:
        return {"error": f"Unknown field type: {field_type}"}

    bkg_name = args.get("name", "background_field")
    state[bkg_name] = bkg

    return {
        "success": True,
        "background_field_name": bkg_name,
        "object_id": bkg,
        "field_type": field_type,
        "description": field_desc
    }


def _get_hacapk_stats(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Get HACApK statistics."""
    try:
        stats = rad.GetHACApKStats()
        return {
            "success": True,
            "stats": stats
        }
    except AttributeError:
        return {
            "success": False,
            "message": "GetHACApKStats not available. Make sure to solve with Method 2 first."
        }
