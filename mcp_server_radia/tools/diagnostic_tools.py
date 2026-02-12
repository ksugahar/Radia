"""
Diagnostic and Debugging Tools for Radia MCP Server

Tools for server health check, state inspection, and debugging.
"""

from typing import Any, Dict, List
from mcp.types import Tool

__version__ = "1.0.0"


def get_tools() -> List[Tool]:
    """Get list of diagnostic tools."""
    return [
        Tool(
            name="radia_server_info",
            description="Get server version and status information",
            inputSchema={
                "type": "object",
                "properties": {},
                "required": []
            }
        ),
        Tool(
            name="radia_list_objects",
            description="List all Radia objects currently in server state",
            inputSchema={
                "type": "object",
                "properties": {},
                "required": []
            }
        ),
        Tool(
            name="radia_get_object_info",
            description="Get detailed information about a specific Radia object",
            inputSchema={
                "type": "object",
                "properties": {
                    "object_name": {
                        "type": "string",
                        "description": "Name of the Radia object"
                    }
                },
                "required": ["object_name"]
            }
        ),
        Tool(
            name="radia_clear_state",
            description="Clear all objects from server state (reset)",
            inputSchema={
                "type": "object",
                "properties": {
                    "confirm": {
                        "type": "boolean",
                        "description": "Must be true to confirm clearing all state",
                        "default": False
                    }
                },
                "required": ["confirm"]
            }
        ),
    ]


async def execute(name: str, arguments: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Execute a diagnostic tool."""
    try:
        if name == "radia_server_info":
            return _server_info(state)
        elif name == "radia_list_objects":
            return _list_objects(state)
        elif name == "radia_get_object_info":
            return _get_object_info(arguments, state)
        elif name == "radia_clear_state":
            return _clear_state(arguments, state)
        else:
            return {"error": f"Unknown diagnostic tool: {name}"}
    except Exception as e:
        return {"error": str(e), "tool": name, "traceback": __import__('traceback').format_exc()}


def _server_info(state: Dict[str, Any]) -> Dict[str, Any]:
    """Get server information."""
    try:
        import radia as rad
        radia_available = True
        radia_version = getattr(rad, '__version__', 'unknown')
    except ImportError:
        radia_available = False
        radia_version = None

    return {
        "success": True,
        "server": "Radia MCP Server",
        "version": __version__,
        "radia_available": radia_available,
        "radia_version": radia_version,
        "state_objects": len(state),
        "object_names": list(state.keys())
    }


def _list_objects(state: Dict[str, Any]) -> Dict[str, Any]:
    """List all objects in state."""
    objects = []

    for name, obj in state.items():
        obj_info = {
            "name": name,
            "type": type(obj).__name__
        }

        # Add Radia-specific info if it's a Radia object
        if isinstance(obj, int):
            obj_info["radia_id"] = obj

        objects.append(obj_info)

    return {
        "success": True,
        "total_objects": len(objects),
        "objects": objects
    }


def _get_object_info(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Get detailed information about a specific object."""
    obj_name = args["object_name"]

    if obj_name not in state:
        return {"error": f"Object '{obj_name}' not found in state"}

    obj = state[obj_name]

    info = {
        "success": True,
        "name": obj_name,
        "type": type(obj).__name__,
        "value": str(obj)
    }

    # Add Radia-specific information
    if isinstance(obj, int):
        try:
            import radia as rad
            # Try to get object properties
            info["radia_id"] = obj
            # Note: Radia doesn't have a direct way to query object properties
            # but we can note it's a valid Radia object ID
        except ImportError:
            pass

    return info


def _clear_state(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Clear all objects from state."""
    if not args.get("confirm", False):
        return {
            "error": "Must set confirm=true to clear state",
            "warning": "This will delete all objects from server memory"
        }

    # Store count before clearing
    object_count = len(state)
    object_names = list(state.keys())

    # Clear the state dictionary
    state.clear()

    # Also call Radia's cleanup
    try:
        import radia as rad
        rad.UtiDelAll()
    except ImportError:
        pass

    return {
        "success": True,
        "message": "State cleared successfully",
        "objects_removed": object_count,
        "removed_names": object_names
    }
