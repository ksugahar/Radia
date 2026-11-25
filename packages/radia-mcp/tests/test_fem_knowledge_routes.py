"""Every topic text of the retired ngsolve-mcp servers stays reachable from radia-mcp.

The fixture holds the substantive FEM topics from 5813761fe.
Each text must be returned by the named Radia MCP route.
"""

import json
import pathlib

from radia_mcp.fem.gauge_open_boundary_knowledge import get_gauge_open_boundary_knowledge
from radia_mcp.motor.age_quality_knowledge import get_age_quality_report
from radia_mcp.radia_ngsolve.knowledge.dtn_coarse_mesh import get_dtn_coarse_mesh_documentation
from radia_mcp.radia_ngsolve.knowledge.kelvin import get_kelvin_documentation
from radia_mcp.radia_ngsolve.knowledge.ngsolve import get_ngsolve_documentation
from radia_mcp.radia_ngsolve.knowledge.sparsesolv import get_sparsesolv_documentation

FIXTURE = pathlib.Path(__file__).parent / "data" / "fem_migrated_topics.json"

# server -> (radia-mcp dispatcher, {old topic: new topic})
ROUTES = {
    "core": (get_ngsolve_documentation, {}),
    "weak_form": (get_ngsolve_documentation, {}),
    "kelvin": (get_kelvin_documentation, {}),
    "dtn": (get_dtn_coarse_mesh_documentation, {}),
    "gauging": (get_gauge_open_boundary_knowledge, {"open_boundary": "gauge_open_boundary"}),
    "age": (get_age_quality_report, {}),
}
def test_every_retired_topic_has_a_radia_mcp_route():
    tables = json.loads(FIXTURE.read_text(encoding="utf-8"))
    assert set(tables) == set(ROUTES)
    missing = []
    for server, topics in tables.items():
        dispatch, renamed = ROUTES[server]
        for topic, text in topics.items():
            if text not in dispatch(renamed.get(topic, topic)):
                missing.append(f"{server}:{topic}")
    assert not missing, missing


def test_unlanded_ngsolve_mcp_know_how_is_in_radia_mcp():
    """AMS and Newton notes that were only on the unpushed ngsolve-mcp branch."""
    ams = get_sparsesolv_documentation("ams")
    for key in ("auto_line_smoothing", "stable_pmis_rank", "line_max_length", 'Backend="gpu"'):
        assert key in ams, key
    nonlinear = get_ngsolve_documentation("nonlinear")
    assert "d nu / d(B^2) >= 0" in nonlinear
    assert "ramp the excitation" in nonlinear
