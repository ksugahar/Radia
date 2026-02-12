"""
Material Tools for Radia MCP Server

Tools for defining and applying material properties:
- Linear materials (soft magnetic materials)
- Nonlinear materials (B-H curves)
- Permanent magnet materials
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
    """Get list of material definition tools."""
    return [
        Tool(
            name="radia_material_create_linear",
            description="Create a linear isotropic magnetic material using MatLin. "
                       "For soft magnetic materials (iron, steel, mu-metal).",
            inputSchema={
                "type": "object",
                "properties": {
                    "mu_r": {
                        "type": "number",
                        "description": "Relative permeability (mu_r >= 1)",
                        "minimum": 1.0
                    },
                    "name": {
                        "type": "string",
                        "description": "Material name for reference in state",
                        "default": "linear_material"
                    }
                },
                "required": ["mu_r"]
            }
        ),
        Tool(
            name="radia_material_create_linear_anisotropic",
            description="Create an anisotropic linear material with easy axis.",
            inputSchema={
                "type": "object",
                "properties": {
                    "mu_parallel": {
                        "type": "number",
                        "description": "Parallel relative permeability",
                        "minimum": 1.0
                    },
                    "mu_perpendicular": {
                        "type": "number",
                        "description": "Perpendicular relative permeability",
                        "minimum": 1.0
                    },
                    "easy_axis": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Easy axis direction [ex, ey, ez]"
                    },
                    "name": {
                        "type": "string",
                        "description": "Material name for reference in state",
                        "default": "anisotropic_material"
                    }
                },
                "required": ["mu_parallel", "mu_perpendicular", "easy_axis"]
            }
        ),
        Tool(
            name="radia_material_create_nonlinear",
            description="Create a nonlinear isotropic material using MatSatIsoTab with B-H curve.",
            inputSchema={
                "type": "object",
                "properties": {
                    "bh_curve": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": {"type": "number"},
                            "minItems": 2,
                            "maxItems": 2
                        },
                        "description": "B-H curve data [[H1, B1], [H2, B2], ...] where H is in A/m and B is in Tesla"
                    },
                    "name": {
                        "type": "string",
                        "description": "Material name for reference in state",
                        "default": "nonlinear_material"
                    }
                },
                "required": ["bh_curve"]
            }
        ),
        Tool(
            name="radia_material_create_permanent_magnet_fixed",
            description="Create a permanent magnet material with fixed magnetization (MatMagFixed).",
            inputSchema={
                "type": "object",
                "properties": {
                    "magnetization": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Magnetization vector [Mx, My, Mz] in A/m"
                    },
                    "name": {
                        "type": "string",
                        "description": "Material name for reference in state",
                        "default": "pm_fixed"
                    }
                },
                "required": ["magnetization"]
            }
        ),
        Tool(
            name="radia_material_create_permanent_magnet_linear",
            description="Create a permanent magnet material with linear demagnetization (MatMagLinear).",
            inputSchema={
                "type": "object",
                "properties": {
                    "br": {
                        "type": "number",
                        "description": "Remanence flux density Br in Tesla"
                    },
                    "hc": {
                        "type": "number",
                        "description": "Coercivity Hc in A/m"
                    },
                    "easy_axis": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Easy axis direction [ex, ey, ez]"
                    },
                    "name": {
                        "type": "string",
                        "description": "Material name for reference in state",
                        "default": "pm_linear"
                    }
                },
                "required": ["br", "hc", "easy_axis"]
            }
        ),
        Tool(
            name="radia_material_apply",
            description="Apply a material to a Radia object using MatApl.",
            inputSchema={
                "type": "object",
                "properties": {
                    "object_name": {
                        "type": "string",
                        "description": "Name of the object to apply material to"
                    },
                    "material_name": {
                        "type": "string",
                        "description": "Name of the material to apply"
                    }
                },
                "required": ["object_name", "material_name"]
            }
        ),
    ]


async def execute(name: str, arguments: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Execute a material tool."""
    if rad is None:
        return {"error": "Radia module not available. Please build Radia first."}

    try:
        if name == "radia_material_create_linear":
            return _create_linear(arguments, state)
        elif name == "radia_material_create_linear_anisotropic":
            return _create_linear_anisotropic(arguments, state)
        elif name == "radia_material_create_nonlinear":
            return _create_nonlinear(arguments, state)
        elif name == "radia_material_create_permanent_magnet_fixed":
            return _create_pm_fixed(arguments, state)
        elif name == "radia_material_create_permanent_magnet_linear":
            return _create_pm_linear(arguments, state)
        elif name == "radia_material_apply":
            return _apply_material(arguments, state)
        else:
            return {"error": f"Unknown material tool: {name}"}
    except Exception as e:
        return {"error": str(e), "tool": name}


def _create_linear(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create a linear isotropic material."""
    mat = rad.MatLin(args["mu_r"])
    mat_name = args.get("name", "linear_material")
    state[mat_name] = mat

    return {
        "success": True,
        "material_name": mat_name,
        "material_id": mat,
        "type": "linear_isotropic",
        "mu_r": args["mu_r"]
    }


def _create_linear_anisotropic(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create an anisotropic linear material."""
    mat = rad.MatLin(
        [args["mu_parallel"], args["mu_perpendicular"]],
        args["easy_axis"]
    )
    mat_name = args.get("name", "anisotropic_material")
    state[mat_name] = mat

    return {
        "success": True,
        "material_name": mat_name,
        "material_id": mat,
        "type": "linear_anisotropic",
        "mu_parallel": args["mu_parallel"],
        "mu_perpendicular": args["mu_perpendicular"],
        "easy_axis": args["easy_axis"]
    }


def _create_nonlinear(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create a nonlinear material with B-H curve."""
    mat = rad.MatSatIsoTab(args["bh_curve"])
    mat_name = args.get("name", "nonlinear_material")
    state[mat_name] = mat

    return {
        "success": True,
        "material_name": mat_name,
        "material_id": mat,
        "type": "nonlinear_isotropic",
        "num_bh_points": len(args["bh_curve"])
    }


def _create_pm_fixed(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create a permanent magnet material with fixed magnetization."""
    mat = rad.MatMagFixed(args["magnetization"])
    mat_name = args.get("name", "pm_fixed")
    state[mat_name] = mat

    return {
        "success": True,
        "material_name": mat_name,
        "material_id": mat,
        "type": "permanent_magnet_fixed",
        "magnetization": args["magnetization"]
    }


def _create_pm_linear(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create a permanent magnet material with linear demagnetization."""
    mat = rad.MatMagLinear(args["br"], args["hc"], args["easy_axis"])
    mat_name = args.get("name", "pm_linear")
    state[mat_name] = mat

    return {
        "success": True,
        "material_name": mat_name,
        "material_id": mat,
        "type": "permanent_magnet_linear",
        "br": args["br"],
        "hc": args["hc"],
        "easy_axis": args["easy_axis"]
    }


def _apply_material(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Apply material to an object."""
    obj_name = args["object_name"]
    mat_name = args["material_name"]

    if obj_name not in state:
        return {"error": f"Object '{obj_name}' not found in state"}
    if mat_name not in state:
        return {"error": f"Material '{mat_name}' not found in state"}

    obj = state[obj_name]
    mat = state[mat_name]

    rad.MatApl(obj, mat)

    return {
        "success": True,
        "object_name": obj_name,
        "material_name": mat_name,
        "message": f"Material '{mat_name}' applied to object '{obj_name}'"
    }
