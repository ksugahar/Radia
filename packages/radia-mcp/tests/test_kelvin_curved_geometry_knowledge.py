"""Guards: shipped lookup scope/persistence warnings. Assumes: no solver accuracy claim."""
import ast
import asyncio
import hashlib
import json
import os
from pathlib import Path
import sys

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


def test_documentation_is_linked_to_existing_lookup_and_publication_limits():
    """Guards: a private draft mistaken for delivered docs. Assumes: package-owned guide."""
    package = Path(__file__).resolve().parents[1]
    guide = package / 'docs/kelvin-periodic-geometry.md'
    text = guide.read_text(encoding='utf-8')
    assert 'docs/kelvin-periodic-geometry.md' in (package / 'README.md').read_text(encoding='utf-8')
    for token in ('mcp-server-cubit', 'mcp-server-radia-ngsolve',
                  'kelvin_transformation(topic="curved_geometry")',
                  'kelvin_identify_post_hoc(topic="curved_geometry")',
                  'model_capture.md', 'actual client', 'Static source tests',
                  'Private engineering models'):
        assert token in text, token
    knowledge = get_kelvin_documentation('curved_geometry')
    assert 'before publishing the final `.vol`' in knowledge
    assert 'Producer build success is not upload, CI, main integration' in knowledge


async def _kelvin_stdio_probe():
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    package = Path(__file__).resolve().parents[1]
    server = package / 'src/radia_mcp/radia_ngsolve/server.py'
    env = os.environ.copy()
    env['PYTHONPATH'] = os.pathsep.join([str(package / 'src'), env.get('PYTHONPATH', '')])
    params = StdioServerParameters(command=sys.executable,
                                   args=['-B', '-m', 'radia_mcp.radia_ngsolve.server'],
                                   cwd=str(package), env=env)
    async with stdio_client(params) as streams:
        async with ClientSession(*streams) as session:
            await session.initialize()
            listed = {t.name: t for t in (await session.list_tools()).tools}
            for name in ('radia_ngsolve_status', 'kelvin_transformation', 'kelvin_identify_post_hoc'):
                assert name in listed, name
            for name in ('radia_ngsolve_status', 'kelvin_transformation'):
                assert listed[name].annotations.readOnlyHint is True, name
            # The same existing tool also has topic='verify', which adds in-memory
            # identifications. Its tool-wide annotation must not promise read-only.
            assert listed['kelvin_identify_post_hoc'].annotations.readOnlyHint is False
            status = await session.call_tool('radia_ngsolve_status', {})
            assert not status.isError
            payload = status.structuredContent or json.loads(status.content[0].text)
            provenance = payload['runtime_provenance']
            assert Path(provenance['module_file']).resolve() == server.resolve()
            assert provenance['module_file_sha256'] == hashlib.sha256(server.read_bytes()).hexdigest()
            assert Path(provenance['python_executable']).resolve() == Path(sys.executable).resolve()
            for name in ('kelvin_transformation', 'kelvin_identify_post_hoc'):
                called = await session.call_tool(name, {'topic': 'curved_geometry'})
                assert not called.isError, name
                text = '\n'.join(c.text for c in called.content if c.type == 'text')
                assert 'Curved periodic geometry: producer and consumer contract' in text
                assert 'before publishing the final `.vol`' in text
                assert 'Producer build success is not upload' in text


def test_candidate_server_lists_and_serves_kelvin_knowledge_over_stdio():
    """Guards: client-visible lookup and identity. Assumes: candidate server, not deployment."""
    asyncio.run(asyncio.wait_for(_kelvin_stdio_probe(), timeout=60))
