"""Guards: shipped lookup scope/persistence warnings. Assumes: no solver accuracy claim."""
import ast
from pathlib import Path

import pytest

from radia_mcp.radia_ngsolve.knowledge.kelvin import get_kelvin_documentation
from radia_mcp.radia_ngsolve.knowledge.kelvin_identify_post_hoc import (
    TOPICS, get_post_hoc_documentation,
)


@pytest.mark.parametrize('lookup,topic', [
    (get_kelvin_documentation, 'curved_geometry'),
    (get_kelvin_documentation, 'identify'),
    (get_kelvin_documentation, 'mesh_control'),
    (get_kelvin_documentation, 'all'),
    *[(get_post_hoc_documentation, t) for t in TOPICS],
])
def test_existing_lookups_carry_curved_route_limitations(lookup, topic):
    """Guards: all entry recipes warn before reuse. Assumes: neutral artifact checks."""
    text = lookup(topic)
    for required in ('Keep the pinned Netgen/NGSolve runtime unchanged',
                     'version/route-specific', 'Load/Save',
                     'NOT a general repair', 'caller-owned NGSolve',
                     'not field-accuracy adoption'):
        assert required in text, (lookup.__name__, topic, required)


def test_both_tools_share_the_same_geometry_contract():
    """Guards: independently drifting explanations. Assumes: existing public tool names."""
    text = get_kelvin_documentation('curved_geometry')
    assert text == get_post_hoc_documentation('curved_geometry')
    for required in ('polynomial q1/q2/q3 tetrahedra', 'linear triangular/tetrahedral',
                     'name `kt`, group 1', 'Source/destination overlap',
                     'q3 odd mode', 'non-target archive', 'tangential traces',
                     'WITH gradients', '|x_h| != R'):
        assert required in text, required


def test_post_hoc_does_not_approve_partial_or_conflicting_pairs():
    """Guards: unsafe fallback recipe. Assumes: helper code is not changed here."""
    text = get_post_hoc_documentation('caveats')
    assert 'Nonzero unmatched vertices are incomplete coverage' in text
    assert 'appending is not repair' in text
    assert 'tolerates a few unmatched vertices' not in text
    assert 'the FES result is still correct' not in text


def test_mapped_point_example_retains_its_transformation():
    """Guards: lifetime guidance removal/reversal. Assumes: no native execution in lookup."""
    text = get_kelvin_documentation('curved_geometry')
    start = text.index('    trafo = mesh.GetTrafo(element)')
    stop = text.index('\n\nKeep mesh/geometry', start)
    snippet = '\n'.join(line[4:] for line in text[start:stop].splitlines())
    tree = ast.parse(snippet)
    assert isinstance(tree.body[0], ast.Assign)
    assert tree.body[0].targets[0].id == 'trafo'
    assert tree.body[1].value.func.id == 'trafo'
    assert tree.body[-2].targets[0].id == 'mapped_point'
    assert tree.body[-1].targets[0].id == 'trafo'


def test_existing_server_dispatch_reaches_the_owned_lookup():
    """Guards: knowledge invisible through existing dispatch. Assumes: transport unchanged."""
    server = Path(__file__).resolve().parents[1] / 'src/radia_mcp/radia_ngsolve/server.py'
    tree = ast.parse(server.read_text(encoding='utf-8'))
    functions = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
    for name, getter in [('kelvin_transformation', 'get_kelvin_documentation'),
                         ('kelvin_identify_post_hoc', '_get_kelvin_identify_post_hoc_doc')]:
        node = functions[name]
        assert 'curved_geometry' in ast.get_docstring(node)
        calls = [n for n in ast.walk(node) if isinstance(n, ast.Call)
                 and isinstance(n.func, ast.Name) and n.func.id == getter]
        assert any(n.args and isinstance(n.args[0], ast.Name)
                   and n.args[0].id == 'topic' for n in calls), name
