"""Lightweight parser for tri/tet Netgen ``.vol`` meshes.

The helper reads the first-order topology sections that are useful for
license-free mesh validation and FEM/BEM trace bookkeeping.  It intentionally
accepts only triangular surface elements and tetrahedral volume elements:
quad/hex/pyramid/wedge meshes should be fixed at the mesh-generation step,
not silently split by a downstream solver adapter.
"""

from __future__ import annotations

from dataclasses import dataclass
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
            "policy": "netgen_vol_tri_tet_only_shared_one_based_nodes",
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
