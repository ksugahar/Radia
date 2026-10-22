"""Collect the documentation examples' existing analytic guards without CLI writes.

These eighteen dependency-free examples compare production mesh bookkeeping
with explicit tetrahedron/box geometry, topology and force identities. Their
build_summary functions own the unchanged analytic limits; no saved JSON is
used as an oracle. The examples remain documentation/result generators, while
this parametrized module is their single collected correctness gate.
"""

import importlib
import json
from pathlib import Path
import sys

import pytest


EXAMPLES = Path(__file__).resolve().parents[2] / "validation_test" / "cubit_mesh_export"
CASES = (
    "validation_p1_surface_triangle_element_summary",
    "validation_vol_boundary_condition_assignment",
    "validation_vol_boundary_edge_inventory",
    "validation_vol_boundary_inventory",
    "validation_vol_boundary_normal_vectors",
    "validation_vol_boundary_oriented_edges",
    "validation_vol_boundary_pressure_force",
    "validation_vol_boundary_pressure_moment",
    "validation_vol_boundary_pressure_resultant",
    "validation_vol_boundary_tet_face_incidence",
    "validation_vol_boundary_traction_moment",
    "validation_vol_fem_bem_topology",
    "validation_vol_material_interface",
    "validation_vol_mesh_health",
    "validation_vol_p1_trace_matrix",
    "validation_vol_surface_closure",
    "validation_vol_surface_triangle_quality",
    "validation_vol_tet_quality",
)


@pytest.mark.parametrize("module_name", CASES, ids=CASES)
def test_analytic_summary(module_name, tmp_path, monkeypatch):
    """Run the real independent guards; persist only this run's temporary output."""
    monkeypatch.setattr(sys, "dont_write_bytecode", True)
    monkeypatch.syspath_prepend(str(EXAMPLES))
    before = set(sys.modules)
    try:
        module = importlib.import_module(module_name)
        assert Path(module.__file__).resolve() == EXAMPLES / (module_name + ".py")
        summary = module.build_summary()
        output = tmp_path / (module_name + ".json")
        output.write_text(json.dumps(summary, sort_keys=True), encoding="utf-8")
    finally:
        # Examples import other examples by bare name. Do not leak those modules
        # into unrelated tests or a later run from another source checkout.
        for name in set(sys.modules) - before:
            if name in CASES:
                del sys.modules[name]
