# Cubit-owned implementation; intentionally maintained independently of Radia MCP.
# Derived support retains BSD-3-Clause terms: see LICENSE-BSD-3-Clause.txt.
"""Mesh-only Gmsh 4.1 handoff for the Netgen comparison route (TET, orders 1-3).

NGSolve owns curved geometry evaluation; Gmsh owns its reference-node ordering.
This private adapter deliberately excludes Radia's field/postprocessing APIs.
The artifact is only for the Gmsh quality referee, not a general interchange mesh.
"""
from functools import lru_cache
import math
from pathlib import Path
from ._gmsh_subprocess import run_gmsh_json_subprocess

_REFERENCE_SCRIPT = '''
import json, sys
import gmsh
gmsh.initialize(['-noconfig'])
try:
    element_type = gmsh.model.mesh.getElementType('Tetrahedron', int(sys.argv[1]))
    _, dim, order, count, coordinates, _ = gmsh.model.mesh.getElementProperties(element_type)
    with open(sys.argv[2], 'w', encoding='utf-8') as f:
        json.dump(dict(ok=True, type=element_type, dim=dim, order=order,
                       points=coordinates.reshape(count, dim).tolist()), f)
finally:
    gmsh.finalize()
'''


@lru_cache(maxsize=3)
def _reference_nodes(order: int) -> dict:
    if order not in (1, 2, 3):
        raise ValueError('Comparison supports tetrahedral orders 1-3')
    reference = run_gmsh_json_subprocess(_REFERENCE_SCRIPT, [str(order)], timeout_s=60,
                                         prefix='cubit_mcp_gmsh_reference_')
    if not reference.get('ok'):
        raise RuntimeError('Gmsh reference nodes unavailable: ' + str(reference))
    return reference


def _write_tet_msh(mesh, out_msh: Path, order: int = 1) -> dict:
    """Write and validate the private quality reference, never a solver mesh."""
    reference = _reference_nodes(order)
    from ngsolve import TaskManager, IntegrationRule, VOL
    points = reference['points']
    rule = IntegrationRule(points, [1.0] * len(points))
    nodes, elements, tags = [], [], {}
    with TaskManager():
        if order > 1:
            mesh.Curve(order)
        for element in mesh.Elements(VOL):
            if len(element.vertices) != 4:
                raise ValueError('Netgen comparison requires tetrahedra')
            trafo = mesh.GetTrafo(element)
            conn = []
            for index, (x, y, z) in enumerate(points):
                # NGSolve TET vertex coordinates are e1, e2, e3, origin.
                # Orders 1-3 have barycentric nodes on the 1/order lattice;
                # extending the supported orders requires revalidating this key.
                weights = (x, y, z, 1-x-y-z)
                key = tuple(sorted((vertex.nr, round(weight * order))
                                   for vertex, weight in zip(element.vertices, weights)
                                   if abs(weight) > 1e-10))
                if key not in tags:
                    tags[key] = len(nodes) + 1
                    nodes.append(tuple(trafo(rule[index]).point))
                conn.append(tags[key])
            elements.append(conn)
        summary = {'n_elements': mesh.ne, 'n_vertices': mesh.nv}
    if not nodes or not elements:
        raise ValueError('Netgen comparison produced an empty volume mesh')
    # One discrete volume suffices for this mesh-quality-only referee.
    bounds = [min(p[i] for p in nodes) for i in range(3)] + [max(p[i] for p in nodes) for i in range(3)]
    lines = ['$MeshFormat', '4.1 0 8', '$EndMeshFormat', '$Entities', '0 0 0 1',
             '1 ' + ' '.join(format(v, '.17g') for v in bounds) + ' 0 0',
             '$EndEntities', '$Nodes',
             f'1 {len(nodes)} 1 {len(nodes)}', f'3 1 0 {len(nodes)}']
    lines.extend(str(i) for i in range(1, len(nodes)+1))
    lines.extend(' '.join(format(value, '.17g') for value in point) for point in nodes)
    lines.extend(['$EndNodes', '$Elements', f'1 {len(elements)} 1 {len(elements)}',
                  f"3 1 {reference['type']} {len(elements)}"])
    lines.extend(f"{i} " + ' '.join(map(str, conn)) for i, conn in enumerate(elements, 1))
    lines.append('$EndElements')
    Path(out_msh).write_text('\n'.join(lines) + '\n', encoding='ascii')
    from .mesh_quality import mesh_total_volume
    quality = mesh_total_volume(out_msh, quadrature='Gauss10')
    minimum = quality.get('min_jacobian_det')
    if not quality.get('ok') or minimum is None or not math.isfinite(minimum) or minimum <= 0:
        raise ValueError('Invalid Netgen reference Jacobian: ' + str(quality))
    summary['min_jacobian_det'] = minimum
    return summary


def step_to_msh(step_path: Path, maxh: float, out_msh: Path, order: int = 1) -> dict:
    """Create the Netgen tetrahedral reference mesh without a Radia installation."""
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh, TaskManager
    with TaskManager():
        mesh = Mesh(OCCGeometry(str(step_path)).GenerateMesh(maxh=maxh))
    return _write_tet_msh(mesh, out_msh, order)
