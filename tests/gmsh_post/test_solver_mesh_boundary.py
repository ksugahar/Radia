"""Guards: PEEC accepts checked solver meshes, never Gmsh file inputs.
Assumes: native Netgen fixtures preserve the adopted analytic discretizations.
"""
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]


def source_module(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


peec = source_module('peec_solver_input_boundary', 'src/radia/peec_mesh_import.py')
parser = source_module('fasthenry_solver_input_boundary', 'src/radia/fasthenry_parser.py')
hiruma = source_module('hiruma_solver_input_boundary', 'validation_test/sparsesolv/hiruma/prepare_hiruma_vol.py')


@pytest.mark.parametrize('value', ['display.msh', Path('display.msh'), Path('solver.vol')])
def test_peec_requires_an_already_checked_mesh_object(value, monkeypatch):
    monkeypatch.setitem(sys.modules, 'gmsh', None)
    with pytest.raises(TypeError, match='already checked'):
        peec.surface_mesh_to_peec(value)


@pytest.mark.parametrize('directive', ['panel', 'magnetic'])
def test_fasthenry_mesh_syntax_is_removed_before_native_access(monkeypatch, directive):
    monkeypatch.setitem(sys.modules, 'gmsh', None)
    with pytest.raises(ValueError, match=f'Unknown .{directive} type: mesh'):
        parser.FastHenryParser().parse_string(f'.{directive}\ntype=mesh\nfile=display.msh\n.end{directive}\n')


def test_fasthenry_explicit_panel_geometry_is_unchanged():
    model = parser.FastHenryParser()
    model.parse_string('.panel\ntype=triangle\nvertices=0,0,0,1,0,0,0,1,0\n.endpanel\n')
    assert model.panel_blocks[0]['vertices'] == [[0.,0.,0.], [1.,0.,0.], [0.,1.,0.]]


def test_hiruma_normalizer_rejects_gmsh_without_loading_mesh(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, 'netgen.meshing', None)
    with pytest.raises(ValueError, match='expected .vol input'):
        hiruma._load(tmp_path/'display.msh')


def test_checked_vol_surface_extraction_uses_ngsolve_boundary_semantics(tmp_path):
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    mesh = MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1)
    path = tmp_path/'cube.vol'
    mesh.ngmesh.Save(str(path))
    mesh = ng.Mesh(str(path))
    nodes, faces = peec._extract_ngsolve_surface(mesh)
    expected = [tuple(v.nr for v in element.vertices) for element in mesh.Elements(ng.BND)]
    assert faces == expected
    assert len(faces) == mesh.GetNE(ng.BND) == 6
    np.testing.assert_array_equal(nodes, [v.point for v in mesh.vertices])


@pytest.fixture
def recording_builder(monkeypatch):
    class Builder:
        def __init__(self):
            self.nodes, self.segments, self.panels, self.ports = [], [], [], []

        def add_node_at(self, *point):
            self.nodes.append(point)
            return len(self.nodes) - 1

        def add_connected_segment(self, *segment):
            self.segments.append(segment)

        def add_panel(self, vertices):
            self.panels.append(vertices)

        def add_port(self, *port):
            self.ports.append(port)

    monkeypatch.setitem(sys.modules, 'radia.peec_matrices',
                        SimpleNamespace(PEECBuilder=Builder))
    return Builder


def test_ngsolve_2d_mesh_builds_volume_cells_as_surface_panels(recording_builder):
    from ngsolve.meshes import MakeStructured2DMesh
    mesh = MakeStructured2DMesh(quads=True, nx=1, ny=1)
    builder = peec.surface_mesh_to_peec(mesh, ports=[])
    assert set(builder.nodes) == {(0., 0., 0.), (1., 0., 0.),
                                  (0., 1., 0.), (1., 1., 0.)}
    assert len(builder.panels) == 1
    assert len(builder.panels[0]) == 4
    assert len(builder.segments) == 4


def test_explicit_dictionary_preserves_geometry_scaling_and_ports(recording_builder):
    builder = peec.surface_mesh_to_peec(
        {'nodes': [[0, 0, 0], [1000, 0, 0], [0, 1000, 0]],
         'faces': [[0, 1, 2]]}, unit_scale=1e-3, ports=[(0, 2)])
    assert builder.nodes == [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.)]
    assert len(builder.segments) == 3
    assert builder.panels == [[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]]]
    assert builder.ports == [(0, 2)]


@pytest.mark.parametrize('closed,edges', [(False, [[1, 2], [2, 3]]),
                                       (True, [[1, 2], [2, 3], [3, 1]])])
def test_centerline_closed_contract_matches_connectivity(closed, edges):
    helper = source_module('centerline_contract_boundary',
                           'validation_test/peec_integration/verification/centerline_vol.py')
    nodes = {1: [0, 0, 0], 2: [1, 0, 0], 3: [0, 1, 0]}
    contract = {'units': 'm', 'nodes': 3, 'segments': len(edges), 'closed': closed}
    helper._validate_centerline_contract(nodes, edges, contract)
    with pytest.raises(ValueError, match='contract closed'):
        helper._validate_centerline_contract(nodes, edges, {**contract, 'closed': not closed})


@pytest.mark.parametrize('units', ['mm', None])
def test_demo_rejects_non_si_contract_before_native_loading(tmp_path, monkeypatch, units):
    import json
    helper = source_module('demo_si_input_boundary',
                           'docs/peec_integration/demos/gmsh_models/solver_vol.py')
    contract = tmp_path / 'labels.json'
    contract.write_text(json.dumps({'units': units}), encoding='utf-8')
    monkeypatch.setattr(sys, 'argv', ['demo', '--mesh', str(tmp_path / 'mesh.vol'),
                                     '--contract', str(contract)])
    monkeypatch.setitem(sys.modules, 'ngsolve', None)
    with pytest.raises(ValueError, match='SI units'):
        helper.load_solver_mesh()


@pytest.mark.parametrize('kind,node_count,edge_count', [('wire',11,10),('coil',36,36)])
def test_netgen_centerlines_preserve_existing_analytic_discretizations(kind, node_count, edge_count):
    import ngsolve as ng
    path = ROOT/'validation_test/peec_integration/verification/netgen_mesh_generation'/f'circular_{kind}_centerline.vol'
    # Binary Git storage must preserve portable LF, including boundary labels.
    assert b'\r' not in path.read_bytes(), path
    import subprocess
    relative = path.relative_to(ROOT).as_posix()
    committed = subprocess.check_output(['git', '-c', f'safe.directory={ROOT.as_posix()}', 'show', f'HEAD:{relative}'], cwd=ROOT)
    assert b'\r' not in committed, relative
    mesh = ng.Mesh(str(path))
    points = {i: np.asarray(point.p) for i, point in enumerate(mesh.ngmesh.Points(), 1)}
    edges = [tuple(vertex.nr for vertex in element.vertices) for element in mesh.ngmesh.Elements1D()]
    length = sum(np.linalg.norm(points[b]-points[a]) for a, b in edges)
    assert len(points) == node_count and len(edges) == edge_count
    assert mesh.GetMaterials() == (f'{kind}_centerline',)
    assert mesh.GetBoundaries() == (f'{kind}_centerline',)
    if kind == 'wire':
        assert length == pytest.approx(1.0)
    else:
        # Existing geometric tolerance, copied from the retired reader test.
        assert length == pytest.approx(2.0*np.pi*0.05, rel=2e-3)
