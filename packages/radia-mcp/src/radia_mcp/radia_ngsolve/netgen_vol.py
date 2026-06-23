"""Lightweight parser for tri/tet Netgen ``.vol`` meshes.

The helper reads the first-order topology sections that are useful for
license-free mesh validation and FEM/BEM trace bookkeeping.  It intentionally
accepts only triangular surface elements and tetrahedral volume elements:
quad/hex/pyramid/wedge meshes should be fixed at the mesh-generation step,
not silently split by a downstream solver adapter.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path


@dataclass(frozen=True)
class NetgenSurfaceTriangle:
    """Boundary triangle record from a Netgen ``surfaceelements`` section."""

    surfnr: int
    bcnr: int
    domin: int
    domout: int
    nodes: tuple[int, int, int]


@dataclass(frozen=True)
class NetgenTetrahedron:
    """Tetrahedron record from a Netgen ``volumeelements`` section."""

    matnr: int
    nodes: tuple[int, int, int, int]


@dataclass(frozen=True)
class NetgenTriTetVolMesh:
    """Tri/tet-only view of a Netgen ``.vol`` mesh."""

    points: tuple[tuple[float, float, float], ...]
    surface_triangles: tuple[NetgenSurfaceTriangle, ...]
    tetrahedra: tuple[NetgenTetrahedron, ...]
    materials: dict[int, str]
    boundary_names: dict[int, str]
    face_descriptors: tuple[tuple[int, ...], ...] = ()
    dimension: int = 3
    geomtype: int | None = None
    source: str | None = None

    def summary(self) -> dict[str, int]:
        """Return a compact mesh inventory."""

        return {
            "points": len(self.points),
            "surface_triangles": len(self.surface_triangles),
            "tetrahedra": len(self.tetrahedra),
            "materials": len(self.materials),
            "boundary_names": len(self.boundary_names),
        }

    def trace_node_ids(self) -> tuple[int, ...]:
        """Return sorted one-based node ids used by boundary triangles."""

        return tuple(sorted({node for tri in self.surface_triangles for node in tri.nodes}))

    def bounding_box(self) -> dict[str, tuple[float, float]]:
        """Return axis-aligned coordinate bounds."""

        if not self.points:
            raise ValueError("cannot compute a bounding box for a mesh with no points")
        xs = [p[0] for p in self.points]
        ys = [p[1] for p in self.points]
        zs = [p[2] for p in self.points]
        return {"x": (min(xs), max(xs)), "y": (min(ys), max(ys)), "z": (min(zs), max(zs))}

    def tetrahedron_signed_volumes(self) -> tuple[float, ...]:
        """Return signed tetrahedron volumes using the stored node orientation."""

        volumes: list[float] = []
        for tet in self.tetrahedra:
            a, b, c, d = (self.points[node - 1] for node in tet.nodes)
            ab = _sub(b, a)
            ac = _sub(c, a)
            ad = _sub(d, a)
            volumes.append(_dot(ab, _cross(ac, ad)) / 6.0)
        return tuple(volumes)

    def tetrahedron_volumes(self) -> tuple[float, ...]:
        """Return positive tetrahedron volumes."""

        return tuple(abs(v) for v in self.tetrahedron_signed_volumes())

    def total_volume(self) -> float:
        """Return the sum of positive tetrahedron volumes."""

        return sum(self.tetrahedron_volumes())

    def surface_triangle_areas(self) -> tuple[float, ...]:
        """Return positive boundary triangle areas."""

        areas: list[float] = []
        for tri in self.surface_triangles:
            a, b, c = (self.points[node - 1] for node in tri.nodes)
            areas.append(0.5 * _norm(_cross(_sub(b, a), _sub(c, a))))
        return tuple(areas)

    def total_surface_area(self) -> float:
        """Return the sum of boundary triangle areas."""

        return sum(self.surface_triangle_areas())

    def surface_area_by_boundary_number(self) -> dict[int, float]:
        """Return boundary surface area grouped by Netgen boundary number."""

        areas: dict[int, float] = {}
        for tri, area in zip(self.surface_triangles, self.surface_triangle_areas()):
            areas[tri.bcnr] = areas.get(tri.bcnr, 0.0) + area
        return areas

    def fem_bem_trace_view(self) -> dict[str, object]:
        """Return shared-node volume/surface connectivity for FEM/BEM coupling."""

        trace_nodes = list(self.trace_node_ids())
        return {
            "points": [list(point) for point in self.points],
            "tetrahedra": [list(tet.nodes) for tet in self.tetrahedra],
            "tetrahedron_material_numbers": [tet.matnr for tet in self.tetrahedra],
            "surface_triangles": [list(tri.nodes) for tri in self.surface_triangles],
            "surface_boundary_numbers": [tri.bcnr for tri in self.surface_triangles],
            "trace_node_ids": trace_nodes,
            "boundary_names": dict(self.boundary_names),
            "materials": dict(self.materials),
            "total_volume": self.total_volume(),
            "total_surface_area": self.total_surface_area(),
            "policy": "netgen_vol_tri_tet_only_shared_one_based_nodes",
        }

    def first_order_fem_bem_topology(self) -> dict[str, object]:
        """Return first-order H1/HCurl/P1/RWG topology with one-based ids.

        This mirrors the small MATLAB prototype API: H1 uses volume nodes,
        HCurl uses first-order tetrahedron edges, scalar BEM uses compacted
        boundary nodes, and RWG uses closed-manifold boundary edges.
        """

        volume_edges, tet_edges, tet_edge_signs = _build_tet_edges(
            [tet.nodes for tet in self.tetrahedra]
        )
        trace_nodes = self.trace_node_ids()
        trace_node_to_local = {node: i for i, node in enumerate(trace_nodes, start=1)}
        boundary_triangles = tuple(
            tuple(trace_node_to_local[node] for node in tri.nodes) for tri in self.surface_triangles
        )
        surface_edges_local, tri_edges, tri_edge_signs, edge_triangles, opposite_vertices = _build_tri_edges(
            boundary_triangles
        )
        surface_edges_global = tuple(
            tuple(trace_nodes[node - 1] for node in edge) for edge in surface_edges_local
        )
        volume_edge_to_id = {edge: i for i, edge in enumerate(volume_edges, start=1)}
        try:
            surface_to_hcurl = tuple(volume_edge_to_id[edge] for edge in surface_edges_global)
        except KeyError as exc:
            raise ValueError("a boundary RWG edge is not present in the volume HCurl edge set") from exc

        rwg_dof_edge_ids = tuple(i for i, adjacent in enumerate(edge_triangles, start=1) if all(adjacent))
        rwg_to_hcurl_edge_ids = tuple(surface_to_hcurl[i - 1] for i in rwg_dof_edge_ids)

        return {
            "h1": {
                "node_ids": list(range(1, len(self.points) + 1)),
                "trace_node_ids": list(trace_nodes),
            },
            "hcurl": {
                "edges": [list(edge) for edge in volume_edges],
                "tet_edges": [list(row) for row in tet_edges],
                "tet_edge_signs": [list(row) for row in tet_edge_signs],
            },
            "scalar_bem": {
                "node_ids": list(range(1, len(trace_nodes) + 1)),
                "global_node_ids": list(trace_nodes),
                "triangles": [list(tri) for tri in boundary_triangles],
            },
            "rwg": {
                "edges_local": [list(edge) for edge in surface_edges_local],
                "edges_global": [list(edge) for edge in surface_edges_global],
                "tri_edges": [list(row) for row in tri_edges],
                "tri_edge_signs": [list(row) for row in tri_edge_signs],
                "edge_triangles": [list(row) for row in edge_triangles],
                "opposite_vertices_local": [list(row) for row in opposite_vertices],
                "dof_edge_ids": list(rwg_dof_edge_ids),
                "dof_edges_local": [list(surface_edges_local[i - 1]) for i in rwg_dof_edge_ids],
                "dof_edges_global": [list(surface_edges_global[i - 1]) for i in rwg_dof_edge_ids],
                "hcurl_edge_ids": list(rwg_to_hcurl_edge_ids),
            },
            "trace": {
                "h1_to_scalar_bem_rows": list(range(1, len(trace_nodes) + 1)),
                "h1_to_scalar_bem_cols": list(trace_nodes),
                "rwg_to_hcurl_edge_ids": list(rwg_to_hcurl_edge_ids),
            },
            "policy": "first_order_h1_p1_hcurl_nedelec0_bem_p1_rwg_only",
        }


def read_netgen_tri_tet_vol(path: str | Path) -> NetgenTriTetVolMesh:
    """Read a tri/tet Netgen ``.vol`` file from disk."""

    p = Path(path)
    return parse_netgen_tri_tet_vol(p.read_text(encoding="utf-8"), source=str(p))


def parse_netgen_tri_tet_vol(text: str, source: str | None = None) -> NetgenTriTetVolMesh:
    """Parse a Netgen ``.vol`` mesh, rejecting non-tri/tet element records."""

    lines = [line.strip() for line in text.splitlines()]
    index = 0
    dimension = 3
    geomtype: int | None = None
    face_descriptors: list[tuple[int, ...]] = []
    surface_triangles: list[NetgenSurfaceTriangle] = []
    tetrahedra: list[NetgenTetrahedron] = []
    points: list[tuple[float, float, float]] = []
    materials: dict[int, str] = {}
    boundary_names: dict[int, str] = {}

    def next_data_line() -> str | None:
        nonlocal index
        while index < len(lines):
            line = lines[index]
            index += 1
            if not line or line.startswith("#"):
                continue
            return line
        return None

    def next_required_line(context: str) -> str:
        line = next_data_line()
        if line is None:
            raise ValueError(f"unexpected end of .vol while reading {context}")
        return line

    def next_count(context: str) -> int:
        return int(next_required_line(context).split()[0])

    def skip_counted_section(context: str) -> None:
        count = next_count(context)
        for _ in range(count):
            next_required_line(context)

    while True:
        key = next_data_line()
        if key is None:
            break
        key_lower = key.lower()
        if key_lower == "endmesh":
            break
        if key_lower == "mesh3d":
            continue
        if key_lower == "dimension":
            dimension = int(next_required_line("dimension").split()[0])
        elif key_lower == "geomtype":
            geomtype = int(next_required_line("geomtype").split()[0])
        elif key_lower == "facedescriptors":
            count = next_count("facedescriptors")
            for _ in range(count):
                face_descriptors.append(tuple(int(v) for v in next_required_line("facedescriptors").split()))
        elif key_lower == "surfaceelements":
            count = next_count("surfaceelements")
            for i in range(1, count + 1):
                values = [int(v) for v in next_required_line("surfaceelements").split()]
                if len(values) < 6:
                    raise ValueError("surface element line is too short")
                node_count = values[4]
                nodes = tuple(values[5 : 5 + node_count])
                if len(nodes) != node_count:
                    raise ValueError("surface element node count does not match record")
                if node_count != 3:
                    raise ValueError(
                        f"Netgen .vol tri/tet-only policy rejected surface element {i}: "
                        f"expected 3 nodes, got {node_count}"
                    )
                surface_triangles.append(NetgenSurfaceTriangle(values[0], values[1], values[2], values[3], nodes))
        elif key_lower == "volumeelements":
            count = next_count("volumeelements")
            for i in range(1, count + 1):
                values = [int(v) for v in next_required_line("volumeelements").split()]
                if len(values) < 3:
                    raise ValueError("volume element line is too short")
                node_count = values[1]
                nodes = tuple(values[2 : 2 + node_count])
                if len(nodes) != node_count:
                    raise ValueError("volume element node count does not match record")
                if node_count != 4:
                    raise ValueError(
                        f"Netgen .vol tri/tet-only policy rejected volume element {i}: "
                        f"expected 4 nodes, got {node_count}"
                    )
                tetrahedra.append(NetgenTetrahedron(values[0], nodes))
        elif key_lower == "points":
            count = next_count("points")
            for _ in range(count):
                values = [float(v) for v in next_required_line("points").split()]
                if len(values) < 3:
                    raise ValueError("point line must contain x y z")
                points.append((values[0], values[1], values[2]))
        elif key_lower == "materials":
            count = next_count("materials")
            for _ in range(count):
                record = next_required_line("materials").split(maxsplit=1)
                materials[int(record[0])] = record[1] if len(record) > 1 else f"material_{record[0]}"
        elif key_lower == "bcnames":
            count = next_count("bcnames")
            for _ in range(count):
                record = next_required_line("bcnames").split(maxsplit=1)
                boundary_names[int(record[0])] = record[1] if len(record) > 1 else f"boundary_{record[0]}"
        elif key_lower in {"edgesegmentsgi2", "pointelements", "face_colours", "face_transparencies"}:
            skip_counted_section(key_lower)
        elif key_lower == "curvedelements":
            break
        else:
            raise ValueError(f"unsupported or unexpected .vol section: {key}")

    mesh = NetgenTriTetVolMesh(
        points=tuple(points),
        surface_triangles=tuple(surface_triangles),
        tetrahedra=tuple(tetrahedra),
        materials=materials,
        boundary_names=boundary_names,
        face_descriptors=tuple(face_descriptors),
        dimension=dimension,
        geomtype=geomtype,
        source=source,
    )
    _validate_node_references(mesh)
    return mesh


def _validate_node_references(mesh: NetgenTriTetVolMesh) -> None:
    npoints = len(mesh.points)
    for section_name, elements in (
        ("surface element", mesh.surface_triangles),
        ("volume element", mesh.tetrahedra),
    ):
        for i, elem in enumerate(elements, start=1):
            for node in elem.nodes:
                if node < 1 or node > npoints:
                    raise ValueError(f"{section_name} {i} references node {node}, but point count is {npoints}")


def _build_tet_edges(
    tetrahedra: list[tuple[int, int, int, int]],
) -> tuple[
    tuple[tuple[int, int], ...],
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, ...], ...],
]:
    local_pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    raw_edges: list[tuple[int, int]] = []
    raw_signs: list[int] = []
    for tet in tetrahedra:
        for a, b in local_pairs:
            pair = (tet[a], tet[b])
            sorted_pair = tuple(sorted(pair))
            raw_edges.append(sorted_pair)
            raw_signs.append(1 if pair == sorted_pair else -1)

    edges = tuple(sorted(set(raw_edges)))
    edge_to_id = {edge: i for i, edge in enumerate(edges, start=1)}
    tet_edges: list[tuple[int, ...]] = []
    tet_edge_signs: list[tuple[int, ...]] = []
    index = 0
    for _ in tetrahedra:
        tet_edges.append(
            tuple(edge_to_id[raw_edges[index + k]] for k in range(len(local_pairs)))
        )
        tet_edge_signs.append(
            tuple(raw_signs[index + k] for k in range(len(local_pairs)))
        )
        index += len(local_pairs)

    return edges, tuple(tet_edges), tuple(tet_edge_signs)


def _build_tri_edges(
    triangles: tuple[tuple[int, int, int], ...],
) -> tuple[
    tuple[tuple[int, int], ...],
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, int], ...],
    tuple[tuple[int, int], ...],
]:
    local_pairs = ((0, 1), (1, 2), (2, 0))
    raw_edges: list[tuple[int, int]] = []
    raw_signs: list[int] = []
    raw_opposites: list[int] = []
    for tri in triangles:
        for a, b in local_pairs:
            pair = (tri[a], tri[b])
            sorted_pair = tuple(sorted(pair))
            raw_edges.append(sorted_pair)
            raw_signs.append(1 if pair == sorted_pair else -1)
            raw_opposites.append(next(tri[i] for i in range(3) if i not in (a, b)))

    edges = tuple(sorted(set(raw_edges)))
    edge_to_id = {edge: i for i, edge in enumerate(edges, start=1)}
    tri_edges: list[tuple[int, ...]] = []
    tri_edge_signs: list[tuple[int, ...]] = []
    edge_triangles = [[0, 0] for _ in edges]
    opposite_vertices = [[0, 0] for _ in edges]

    index = 0
    for tri_id, _tri in enumerate(triangles, start=1):
        tri_edges.append(
            tuple(edge_to_id[raw_edges[index + k]] for k in range(len(local_pairs)))
        )
        tri_edge_signs.append(
            tuple(raw_signs[index + k] for k in range(len(local_pairs)))
        )
        for k in range(len(local_pairs)):
            edge_id = edge_to_id[raw_edges[index + k]]
            slots = edge_triangles[edge_id - 1]
            try:
                slot = slots.index(0)
            except ValueError as exc:
                raise ValueError("surface edge belongs to more than two triangles") from exc
            slots[slot] = tri_id
            opposite_vertices[edge_id - 1][slot] = raw_opposites[index + k]
        index += len(local_pairs)

    return (
        edges,
        tuple(tri_edges),
        tuple(tri_edge_signs),
        tuple(tuple(row) for row in edge_triangles),
        tuple(tuple(row) for row in opposite_vertices),
    )


def _sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a: tuple[float, float, float]) -> float:
    return math.sqrt(_dot(a, a))
