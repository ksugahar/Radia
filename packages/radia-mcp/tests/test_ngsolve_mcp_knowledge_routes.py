"""Every topic text of the retired ngsolve-mcp servers stays reachable from radia-mcp.

The fixture holds the ngsolve-mcp topic tables as of 5813761fe. Index and
pointer entries carry no knowledge and are listed explicitly; every other
text must be returned by the named radia-mcp route.
"""

import json
import pathlib

from radia_mcp.fem.gauge_open_boundary_knowledge import get_gauge_open_boundary_knowledge
from radia_mcp.motor.age_quality_knowledge import get_age_quality_report
from radia_mcp.radia_ngsolve.knowledge.dtn_coarse_mesh import get_dtn_coarse_mesh_documentation
from radia_mcp.radia_ngsolve.knowledge.kelvin import get_kelvin_documentation
from radia_mcp.radia_ngsolve.knowledge.ngsolve import get_ngsolve_documentation

FIXTURE = pathlib.Path(__file__).parent / "data" / "ngsolve_mcp_retired_topics.json"

# server -> (radia-mcp dispatcher, {old topic: new topic})
ROUTES = {
    "core": (get_ngsolve_documentation, {}),
    "weak_form": (get_ngsolve_documentation, {}),
    "kelvin": (get_kelvin_documentation, {}),
    "dtn": (get_dtn_coarse_mesh_documentation, {}),
    "gauging": (get_gauge_open_boundary_knowledge, {"open_boundary": "gauge_open_boundary"}),
    "age": (get_age_quality_report, {}),
}
NO_KNOWLEDGE = {
    ("core", "index"), ("core", "eddy_current"), ("weak_form", "overview"),
    ("kelvin", "overview"), ("dtn", "overview"), ("gauging", "overview"), ("age", "overview"),
}


def test_every_retired_topic_has_a_radia_mcp_route():
    tables = json.loads(FIXTURE.read_text(encoding="utf-8"))
    assert set(tables) == set(ROUTES)
    missing = []
    for server, topics in tables.items():
        dispatch, renamed = ROUTES[server]
        for topic, text in topics.items():
            if (server, topic) in NO_KNOWLEDGE:
                continue
            if text not in dispatch(renamed.get(topic, topic)):
                missing.append(f"{server}:{topic}")
    assert not missing, missing
