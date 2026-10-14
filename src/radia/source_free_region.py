"""Geometric and real-cohomology prerequisites for a source-free total field."""
from itertools import combinations
import hashlib
import numpy as np


def _box_separated(tetrahedron, lower, upper, tolerance):
    """Separating-axis certificate for a tetrahedron and an axis-aligned box."""
    center = (lower + upper) / 2
    half = (upper - lower) / 2
    edges = [tetrahedron[j] - tetrahedron[i] for i, j in combinations(range(4), 2)]
    axes = list(np.eye(3))
    axes += [np.cross(tetrahedron[j] - tetrahedron[i], tetrahedron[k] - tetrahedron[i])
             for i, j, k in combinations(range(4), 3)]
    axes += [np.cross(edge, axis) for edge in edges for axis in np.eye(3)]
    for axis in axes:
        length = np.linalg.norm(axis)
        if length == 0:
            continue
        axis = axis / length
        projected = (tetrahedron - center) @ axis
        radius = half @ np.abs(axis)
        if projected.min() > radius + tolerance or projected.max() < -radius - tolerance:
            return True
    return False


def certify_current_free_total_region(mesh, materials, source_support_boxes, *, tolerance=1e-10):
    """Fail unless declared current supports are disjoint and real b1 is zero.

    ``source_support_boxes`` is an SI array (N,2,3) of conservative enclosures
    of **all** source-current support. Its completeness is caller-owned physical
    input, like the excitation itself; this function cannot infer it from an
    arbitrary field CoefficientFunction. Touching or uncertain separation fails.
    Requires a valid conforming straight tetrahedral volume mesh. GF(2) b1=0
    implies real b1=0, sufficient for a closed real H one-form to be exact;
    this does not assert simple connectivity or trivial integral homology.
    """
    import ngsolve as ng
    boxes = np.asarray(source_support_boxes, dtype=float)
    names = tuple(materials)
    if not names or not set(names) <= set(mesh.GetMaterials()):
        raise ValueError("current-free total materials must exist")
    if (boxes.ndim != 3 or boxes.shape[1:] != (2, 3) or len(boxes) == 0
            or not np.isfinite(boxes).all() or np.any(boxes[:, 1] <= boxes[:, 0])):
        raise ValueError("declare finite positive-volume source support boxes (N,2,3)")
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("source exclusion tolerance must be positive and finite")
    if mesh.GetCurveOrder() != 1:
        raise ValueError("current-free certificate supports straight tetrahedra only")
    cells = [e for e in mesh.Elements(ng.VOL) if e.mat in names]
    if not cells or any(e.type != ng.ET.TET for e in cells):
        raise ValueError("current-free total region requires nonempty tetrahedra")
    vertices = set(); edges = set(); faces = {}; parent = {}; links = {}
    def find(v):
        parent.setdefault(v, v)
        while parent[v] != v:
            parent[v] = parent[parent[v]]; v = parent[v]
        return v
    for element in cells:
        cell = sorted(v.nr for v in element.vertices)
        xyz = np.asarray([mesh[v].point for v in element.vertices], dtype=float)
        if np.linalg.matrix_rank(xyz[1:] - xyz[0]) != 3:
            raise ValueError("degenerate total tetrahedron")
        candidates = np.all(boxes[:, 1] >= xyz.min(axis=0)-tolerance, axis=1) & np.all(
            boxes[:, 0] <= xyz.max(axis=0)+tolerance, axis=1)
        for lower, upper in boxes[candidates]:
            if not _box_separated(xyz, lower, upper, tolerance):
                raise ValueError("source support intersects or touches total region; current-free mode rejected")
        vertices.update(cell)
        for vertex in cell:
            links.setdefault(vertex, []).append(tuple(v for v in cell if v != vertex))
        for edge in combinations(cell, 2):
            edges.add(edge); parent[find(edge[0])] = find(edge[1])
        for face in combinations(cell, 3):
            faces[face] = faces.get(face, 0) + 1
    if any(count > 2 for count in faces.values()):
        raise ValueError("nonmanifold total-region face")
    # A valid volume vertex has a connected sphere or disk as its link.
    # This rejects components touching only at an edge or a vertex.
    for triangles in links.values():
        link_edges = {}; neighbors = {}
        for triangle in triangles:
            for a, b in combinations(triangle, 2):
                edge = (a, b)
                link_edges[edge] = link_edges.get(edge, 0) + 1
                neighbors.setdefault(a, set()).add(b); neighbors.setdefault(b, set()).add(a)
        visited = set(); stack = [next(iter(neighbors))]
        while stack:
            vertex = stack.pop()
            if vertex not in visited:
                visited.add(vertex); stack.extend(neighbors[vertex] - visited)
        boundary_degree = {}
        for (a, b), count in link_edges.items():
            if count == 1:
                boundary_degree[a] = boundary_degree.get(a, 0) + 1
                boundary_degree[b] = boundary_degree.get(b, 0) + 1
        euler = len(neighbors) - len(link_edges) + len(triangles)
        if (len(visited) != len(neighbors) or any(n > 2 for n in link_edges.values())
                or any(n != 2 for n in boundary_degree.values())
                or euler != (1 if boundary_degree else 2)):
            raise ValueError("nonmanifold total-region vertex link")
    edge_ids = {edge: i for i, edge in enumerate(sorted(edges))}; pivots = {}
    for face in faces:
        row = 0
        for edge in combinations(face, 2):
            row ^= 1 << edge_ids[edge]
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row; break
            row ^= pivots[pivot]
    components = len({find(v) for v in vertices})
    b1 = len(edges) - len(vertices) + components - len(pivots)
    # Integer boundary matrices: a nonzero GF(2) minor is an odd integer
    # determinant, hence nonzero over R. Thus real boundary ranks cannot be
    # smaller and real b1 cannot exceed b1_GF2. Zero is a sufficient gate.
    if b1:
        raise ValueError("nonzero first homology; current-free total mode requires b1_GF2=0")
    return dict(passed=True, scope="declared current supports disjoint; real first cohomology vanishes",
                betti_1_GF2=b1, connected_components=components,
                rank_boundary_1_GF2=len(vertices)-components, rank_boundary_2_GF2=len(pivots),
                source_support_box_count=len(boxes), tolerance_m=float(tolerance),
                source_support_sha256=hashlib.sha256(boxes.astype('<f8').tobytes()).hexdigest(),
                source_support_completeness="caller-owned conservative cover of all currents")
