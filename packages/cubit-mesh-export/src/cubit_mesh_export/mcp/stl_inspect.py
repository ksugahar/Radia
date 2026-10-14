"""Dependency-free STL closure, orientation, volume, and bounds inspection."""

from __future__ import annotations

import math
import struct
from collections import Counter, defaultdict
from collections.abc import Iterable
from pathlib import Path
from typing import Any

# The inspector materializes vertices and edge incidence in Python.  Keep a
# deterministic ceiling that bounds memory before Cubit sees the artifact.
_MAX_STL_BYTES = 32 * 1024 * 1024
_MAX_TRIANGLES = 250_000


def _binary_triangles(data: bytes) -> list[tuple[tuple[float, ...], ...]] | None:
    if len(data) < 84:
        return None
    count = struct.unpack_from("<I", data, 80)[0]
    if 84 + 50 * count != len(data):
        return None
    if count > _MAX_TRIANGLES:
        raise ValueError(f"STL exceeds the {_MAX_TRIANGLES} triangle limit")
    triangles = []
    for index in range(count):
        values = struct.unpack_from("<12fH", data, 84 + 50 * index)
        triangles.append((tuple(values[3:6]), tuple(values[6:9]),
                          tuple(values[9:12])))
    return triangles


def _ascii_triangles(data: bytes) -> list[tuple[tuple[float, ...], ...]]:
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError as exc:
        raise ValueError("STL is neither a valid binary nor ASCII file") from exc
    vertices = []
    for line in text.splitlines():
        fields = line.strip().split()
        if fields and fields[0].lower() == "vertex":
            if len(fields) != 4:
                raise ValueError("malformed ASCII STL vertex row")
            try:
                vertices.append(tuple(float(value) for value in fields[1:]))
            except ValueError as exc:
                raise ValueError("non-numeric ASCII STL vertex") from exc
    if not vertices or len(vertices) % 3:
        raise ValueError("ASCII STL does not contain complete triangle facets")
    count = len(vertices) // 3
    if count > _MAX_TRIANGLES:
        raise ValueError(f"STL exceeds the {_MAX_TRIANGLES} triangle limit")
    return [tuple(vertices[index:index + 3])
            for index in range(0, len(vertices), 3)]


def _cross(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, ...]:
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _sub(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, ...]:
    return tuple(x - y for x, y in zip(a, b))


def _dot(a: tuple[float, ...], b: tuple[float, ...]) -> float:
    return sum(x * y for x, y in zip(a, b))


def _point_inside_component(point, component_faces, *, scale: float) -> bool:
    """Odd-even ray test with duplicate coplanar triangle hits collapsed."""
    direction = (1.0, 0.347, 0.719)
    hits = []
    bary_tol = 1.0e-12
    for face in component_faces:
        edge1 = _sub(face[1], face[0])
        edge2 = _sub(face[2], face[0])
        h = _cross(direction, edge2)
        det = _dot(edge1, h)
        if abs(det) <= 1.0e-15:
            continue
        inv_det = 1.0 / det
        rel = _sub(point, face[0])
        u = inv_det * _dot(rel, h)
        if u < -bary_tol or u > 1.0 + bary_tol:
            continue
        q = _cross(rel, edge1)
        v = inv_det * _dot(direction, q)
        if v < -bary_tol or u + v > 1.0 + bary_tol:
            continue
        distance = inv_det * _dot(edge2, q)
        if distance > max(scale * 1.0e-12, 1.0e-15):
            hits.append(distance)
    hits.sort()
    unique_hits = []
    merge_tol = max(scale * 1.0e-10, 1.0e-13)
    for distance in hits:
        if not unique_hits or abs(distance - unique_hits[-1]) > merge_tol:
            unique_hits.append(distance)
    return len(unique_hits) % 2 == 1


def _surface_metrics(
        triangles: Iterable[tuple[tuple[float, ...], ...]]) -> dict[str, Any]:
    faces = list(triangles)
    points = [point for face in faces for point in face]
    if any(len(point) != 3 or not all(math.isfinite(value) for value in point)
           for point in points):
        raise ValueError("STL contains non-finite or non-3D vertices")
    lo = [min(point[axis] for point in points) for axis in range(3)]
    hi = [max(point[axis] for point in points) for axis in range(3)]
    diagonal = math.sqrt(sum((b - a) ** 2 for a, b in zip(lo, hi)))
    if diagonal <= 0.0:
        raise ValueError("STL bounding box is degenerate")

    vertex_ids: dict[tuple[float, ...], int] = {}
    indexed = []
    edge_counts: Counter[tuple[int, int]] = Counter()
    edge_orientation: Counter[tuple[int, int]] = Counter()
    edge_owner: dict[tuple[int, int], int] = {}
    parent = list(range(len(faces)))
    component_size = [1] * len(faces)
    signed_face_volumes: list[float] = []

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root == right_root:
            return
        if component_size[left_root] < component_size[right_root]:
            left_root, right_root = right_root, left_root
        parent[right_root] = left_root
        component_size[left_root] += component_size[right_root]

    area_floor = max(diagonal * diagonal * 1.0e-24, 1.0e-30)
    for face_index, face in enumerate(faces):
        ids = []
        for point in face:
            if point not in vertex_ids:
                vertex_ids[point] = len(vertex_ids)
            ids.append(vertex_ids[point])
        if len(set(ids)) != 3:
            raise ValueError("STL contains a triangle with repeated vertices")
        ab, ac = _sub(face[1], face[0]), _sub(face[2], face[0])
        cross = _cross(ab, ac)
        if _dot(cross, cross) <= area_floor:
            raise ValueError("STL contains a degenerate triangle")
        indexed.append(tuple(ids))
        signed_face_volumes.append(
            _dot(face[0], _cross(face[1], face[2])) / 6.0)
        for start, end in ((ids[0], ids[1]), (ids[1], ids[2]),
                           (ids[2], ids[0])):
            edge = (min(start, end), max(start, end))
            edge_counts[edge] += 1
            edge_orientation[edge] += 1 if start < end else -1
            if edge in edge_owner:
                union(face_index, edge_owner[edge])
            else:
                edge_owner[edge] = face_index

    open_edges = sum(count == 1 for count in edge_counts.values())
    nonmanifold_edges = sum(count > 2 for count in edge_counts.values())
    boundary_defects = open_edges + nonmanifold_edges
    orientation_errors = sum(
        edge_counts[edge] == 2 and edge_orientation[edge] != 0
        for edge in edge_counts)
    component_volumes: dict[int, float] = defaultdict(float)
    component_face_indices: dict[int, list[int]] = defaultdict(list)
    for face_index, signed_volume in enumerate(signed_face_volumes):
        root = find(face_index)
        component_volumes[root] += signed_volume
        component_face_indices[root].append(face_index)

    component_bounds = {}
    for root, indices in component_face_indices.items():
        component_points = [point for index in indices for point in faces[index]]
        component_bounds[root] = (
            tuple(min(point[axis] for point in component_points)
                  for axis in range(3)),
            tuple(max(point[axis] for point in component_points)
                  for axis in range(3)),
        )
    component_depths = {}
    for root, indices in component_face_indices.items():
        representative = faces[indices[0]][0]
        depth = 0
        for container, container_indices in component_face_indices.items():
            if container == root:
                continue
            c_lo, c_hi = component_bounds[container]
            if not all(c_lo[axis] < representative[axis] < c_hi[axis]
                       for axis in range(3)):
                continue
            container_faces = [faces[index] for index in container_indices]
            if _point_inside_component(
                    representative, container_faces, scale=diagonal):
                depth += 1
        component_depths[root] = depth
    volume = abs(sum(
        (-1.0 if component_depths[root] % 2 else 1.0) * abs(value)
        for root, value in component_volumes.items()
    ))
    return {
        "ok": True,
        "watertight": boundary_defects == 0,
        "winding_consistent": boundary_defects == 0 and orientation_errors == 0,
        "volume": float(volume),
        "bounds": [lo, hi],
        "triangle_count": len(indexed),
        "faces": len(indexed),
        "vertices": len(vertex_ids),
        "open_edge_count": int(open_edges),
        "nonmanifold_edge_count": int(nonmanifold_edges),
        "nonmanifold_or_open_edges": int(boundary_defects),
        "orientation_error_edges": int(orientation_errors),
        "connected_components": len(component_volumes),
        "component_nesting_depths": sorted(component_depths.values()),
    }


def inspect_stl(path: str | Path, *, timeout_s: float = 120.0) -> dict[str, Any]:
    """Inspect binary or ASCII STL without importing NumPy in the MCP process.

    ``timeout_s`` remains in the public helper signature for compatibility;
    deterministic size and triangle limits bound the parser instead.
    """
    del timeout_s
    p = Path(path)
    try:
        size = p.stat().st_size
        if size <= 0:
            raise ValueError("STL is empty")
        if size > _MAX_STL_BYTES:
            raise ValueError(
                f"STL exceeds the {_MAX_STL_BYTES} byte inspection limit")
        data = p.read_bytes()
        triangles = _binary_triangles(data)
        if triangles is None:
            triangles = _ascii_triangles(data)
        return _surface_metrics(triangles)
    except (OSError, ValueError, struct.error) as exc:
        return {
            "ok": False, "kind": "input",
            "error": f"{type(exc).__name__}: {exc}",
        }
