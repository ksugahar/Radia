"""Thin Netgen input helpers for the PEEC demonstrations; SI .vol only."""
import argparse
import json
from pathlib import Path
import numpy as np


def load_solver_mesh():
    parser = argparse.ArgumentParser(description="PEEC from a checked Cubit/Netgen .vol")
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    args = parser.parse_args()
    if args.mesh.suffix.lower() != ".vol":
        parser.error("solver input must be a checked Netgen .vol; Gmsh is display-only")
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    units = contract.get("centerline", {}).get("units", contract.get("units"))
    if units != "m" or ("units" in contract and contract["units"] != "m"):
        raise ValueError("PEEC solver .vol contract must declare SI units 'm'")
    from ngsolve import Mesh
    from cubit_mesh_export.check import check_label_contract
    mesh = Mesh(str(args.mesh))
    report = check_label_contract(mesh, contract=contract)
    if not report["passed"]:
        raise ValueError(f".vol label contract failed: {report}")
    return args.mesh, mesh


def vertex_data(mesh):
    points = np.asarray([point.p for point in mesh.ngmesh.Points()], dtype=float)
    return np.arange(1, len(points)+1, dtype=np.int64), points


def element_groups(mesh):
    """Group native Netgen edges/linear panels by neutral cell family."""
    grouped = {}
    for dimension, elements in [(1, mesh.ngmesh.Elements1D()), (2, mesh.ngmesh.Elements2D())]:
        for number, element in enumerate(elements, 1):
            nodes = [vertex.nr for vertex in element.vertices]
            key = {(1,2): "line", (2,3): "triangle", (2,4): "quad"}.get((dimension,len(nodes)))
            if key is None or len(element.points) != len(nodes):
                raise ValueError("PEEC demo requires explicit linear line/triangle/quad cells")
            tags, connectivity = grouped.setdefault(key, ([], []))
            tags.append(number)
            connectivity.extend(nodes)
    families = list(grouped)
    return families, [np.asarray(grouped[k][0]) for k in families], [np.asarray(grouped[k][1]) for k in families]
