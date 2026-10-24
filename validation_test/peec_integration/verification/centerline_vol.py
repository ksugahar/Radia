"""Checked native Netgen 1D .vol inputs for the existing PEEC cases."""
import json
from pathlib import Path
import numpy as np


def _validate_centerline_contract(nodes, edges, expected):
    if expected["units"] != "m" or len(nodes) != expected["nodes"] or len(edges) != expected["segments"]:
        raise ValueError("Centerline .vol topology/SI contract mismatch")
    if not isinstance(expected["closed"], bool):
        raise ValueError("Centerline contract closed must be boolean")
    degree = {tag: 0 for tag in nodes}
    adjacent = {tag: set() for tag in nodes}
    for a, b in edges:
        if a not in nodes or b not in nodes or a == b:
            raise ValueError("Centerline has invalid segment endpoints")
        degree[a] += 1
        degree[b] += 1
        adjacent[a].add(b)
        adjacent[b].add(a)
    visited = set()
    pending = [next(iter(nodes))] if nodes else []
    while pending:
        node = pending.pop()
        if node not in visited:
            visited.add(node)
            pending.extend(adjacent[node] - visited)
    endpoints = sum(value == 1 for value in degree.values())
    expected_endpoints = 0 if expected["closed"] else 2
    if (not nodes or len(visited) != len(nodes) or endpoints != expected_endpoints
            or any(value not in (1, 2) for value in degree.values())):
        raise ValueError("Centerline connectivity does not match contract closed")


def read_centerline(path):
    from netgen.meshing import Mesh
    from cubit_mesh_export.check import check_label_contract
    path = Path(path)
    if path.suffix.lower() != ".vol":
        raise ValueError("PEEC solver centerlines require Netgen .vol input")
    contract = json.loads(path.with_suffix(".labels.json").read_text(encoding="utf-8"))
    report = check_label_contract(path, contract=contract)
    if not report["passed"]:
        raise ValueError(f"Centerline .vol label contract failed: {report}")
    mesh = Mesh(dim=1)
    mesh.Load(str(path))
    if len(list(mesh.Elements2D())) or len(list(mesh.Elements3D())):
        raise ValueError("Expected a native one-dimensional centerline .vol")
    nodes = {i: list(point.p) for i, point in enumerate(mesh.Points(), 1)}
    edges = [[point.nr for point in element.vertices] for element in mesh.Elements1D()]
    expected = contract["centerline"]
    _validate_centerline_contract(nodes, edges, expected)
    segments = [[nodes[a], nodes[b]] for a, b in edges]
    length = float(sum(np.linalg.norm(np.asarray(a)-np.asarray(b)) for a,b in segments))
    return nodes, edges, segments, length
