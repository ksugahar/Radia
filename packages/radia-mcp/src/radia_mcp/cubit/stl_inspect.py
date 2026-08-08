"""Dependency-free STL closure, orientation, volume, and bounds inspection."""

from __future__ import annotations

from collections import Counter
import math
from pathlib import Path
import struct
from typing import Any, Iterable


_MAX_STL_BYTES = 256 * 1024 * 1024
_MAX_TRIANGLES = 5_000_000


def _binary_triangles(data: bytes) -> list[tuple[tuple[float, ...], ...]] | None:
    if len(data) < 84:
        return None
    count = struct.unpack_from("<I", data, 80)[0]
    if count > _MAX_TRIANGLES or 84 + 50 * count != len(data):
        return None
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
    signed_volume = 0.0
    area_floor = max(diagonal * diagonal * 1.0e-24, 1.0e-30)
    for face in faces:
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
        signed_volume += _dot(face[0], _cross(face[1], face[2])) / 6.0
        for start, end in ((ids[0], ids[1]), (ids[1], ids[2]),
                           (ids[2], ids[0])):
            edge = (min(start, end), max(start, end))
            edge_counts[edge] += 1
            edge_orientation[edge] += 1 if start < end else -1

    open_edges = sum(count == 1 for count in edge_counts.values())
    nonmanifold_edges = sum(count > 2 for count in edge_counts.values())
    boundary_defects = open_edges + nonmanifold_edges
    orientation_errors = sum(
        edge_counts[edge] == 2 and edge_orientation[edge] != 0
        for edge in edge_counts)
    return {
        "ok": True,
        "watertight": boundary_defects == 0,
        "winding_consistent": boundary_defects == 0 and orientation_errors == 0,
        "volume": abs(float(signed_volume)),
        "bounds": [lo, hi],
        "triangle_count": len(indexed),
        "faces": len(indexed),
        "vertices": len(vertex_ids),
        "open_edge_count": int(open_edges),
        "nonmanifold_edge_count": int(nonmanifold_edges),
        "nonmanifold_or_open_edges": int(boundary_defects),
        "orientation_error_edges": int(orientation_errors),
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
