"""Shared loader for calculation regression and independent analytic validation."""
from contextlib import contextmanager
import importlib
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = Path(__file__).resolve().parent
CASES = ('validation_p1_surface_triangle_element_summary', 'validation_vol_boundary_condition_assignment', 'validation_vol_boundary_edge_inventory', 'validation_vol_boundary_inventory', 'validation_vol_boundary_normal_vectors', 'validation_vol_boundary_oriented_edges', 'validation_vol_boundary_pressure_force', 'validation_vol_boundary_pressure_moment', 'validation_vol_boundary_pressure_resultant', 'validation_vol_boundary_tet_face_incidence', 'validation_vol_boundary_traction_moment', 'validation_vol_fem_bem_topology', 'validation_vol_material_interface', 'validation_vol_mesh_health', 'validation_vol_p1_trace_matrix', 'validation_vol_surface_closure', 'validation_vol_surface_triangle_quality', 'validation_vol_tet_quality')

@contextmanager
def example_module(name):
    if name not in CASES:
        raise ValueError(name)
    import radia_mcp
    expected = ROOT / "packages/radia-mcp/src/radia_mcp/__init__.py"
    if Path(radia_mcp.__file__).resolve() != expected:
        raise RuntimeError(f"Unexpected radia_mcp source: {radia_mcp.__file__}; expected {expected}")
    before_modules = set(sys.modules)
    before_path = sys.path[:]
    before_bytecode = sys.dont_write_bytecode
    try:
        sys.dont_write_bytecode = True
        sys.path.insert(0, str(EXAMPLES))
        module = importlib.import_module(name)
        if Path(module.__file__).resolve() != EXAMPLES / (name + ".py"):
            raise RuntimeError(f"Unexpected example source: {module.__file__}")
        yield module
    finally:
        sys.path[:] = before_path
        sys.dont_write_bytecode = before_bytecode
        for key in set(sys.modules) - before_modules:
            if key in CASES:
                del sys.modules[key]
