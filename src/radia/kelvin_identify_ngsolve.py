"""Read-only verification of producer-owned Kelvin periodic point contracts.

Mesh-only verification proves native structure, not a persisted identification
name. Use load_verified_kelvin_mesh for named proof at an immutable file boundary.
No matching, tolerance, identification repair or mesh serialization is performed.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Tuple

import numpy as np

__all__ = ["verify_kelvin_identification", "load_verified_kelvin_mesh", "detect_kelvin_offset"]

def detect_kelvin_offset(mesh) -> Tuple[float, float, float]:
    """Observe mean boundary-coordinate displacement; never match or repair IDs.

    The returned value is not evidence of complete correspondence, geometry
    accuracy or an accepted translation residual. Missing labels return zero for
    non-Kelvin callers. Producer-contract verification is a separate operation.
    """
    from ngsolve import BND

    boundaries = mesh.GetBoundaries()
    if "kelvin_int" not in boundaries or "kelvin_ext" not in boundaries:
        return (0.0, 0.0, 0.0)

    int_pts: list = []
    ext_pts: list = []
    for el in mesh.Elements(BND):
        bc = boundaries[el.index]
        if bc == "kelvin_int":
            int_pts.extend(mesh.vertices[v.nr].point for v in el.vertices)
        elif bc == "kelvin_ext":
            ext_pts.extend(mesh.vertices[v.nr].point for v in el.vertices)

    if not int_pts or not ext_pts:
        return (0.0, 0.0, 0.0)

    int_mean = np.mean(int_pts, axis=0)
    ext_mean = np.mean(ext_pts, axis=0)
    return (float(ext_mean[0] - int_mean[0]),
            float(ext_mean[1] - int_mean[1]),
            float(ext_mean[2] - int_mean[2]))


def _point_nr(value):
    """Pinned Netgen PointId numbers are one-based (not NGSolve VertexId)."""
    nr = value.nr
    if type(nr) is not int or nr < 1:
        raise ValueError(f"Invalid one-based PointId: {nr!r}")
    return nr


def _native_rows(ngmesh):
    rows = []
    for row in ngmesh.GetIdentifications():
        if len(row) != 3:
            raise ValueError(f"Expected pinned identification triple, actual={row!r}")
        group = row[2]
        if type(group) is not int or group < 1:
            raise ValueError(f"Invalid identification group: {group!r}")
        rows.append((_point_nr(row[0]), _point_nr(row[1]), group))
    if len(rows) != len(set(rows)):
        raise ValueError("Duplicate native identification triple")
    return rows


def _boundary_points(mesh, master, slave):
    if not master or not slave or master == slave:
        raise ValueError("Distinct explicit master/slave boundary labels required")
    ngmesh = mesh.ngmesh
    descriptors = list(ngmesh.FaceDescriptors())
    labels = list(mesh.GetBoundaries())
    bc_labels = {}
    names_by_bc = {}
    fd_table = []
    points = {master: set(), slave: set()}
    for index, descriptor in enumerate(descriptors, 1):
        bc = descriptor.bc
        if type(bc) is not int or bc < 1:
            raise ValueError(f"FD {index}: invalid BC property {bc!r}")
        # The descriptor owns its label; GetBoundaries must not be indexed by
        # FD or BC under an assumed equality of those separate namespaces.
        label = descriptor.bcname
        if not isinstance(label, str) or not label or label not in labels:
            raise ValueError(f"FD {index}/BC {bc}: inconsistent boundary label")
        names_by_bc.setdefault(bc, set()).add(label)
        if len(names_by_bc[bc]) != 1:
            raise ValueError(f"BC {bc}: conflicting descriptor labels")
        bc_labels.setdefault(label, set()).add(bc)
        fd_table.append({"fd": index, "bc": bc, "label": label})
    for label in points:
        # A named boundary may intentionally span several BC patches. Union
        # those exact labels; reject distinct case-fold aliases, not split CAD.
        aliases = {actual for actual in bc_labels if actual.casefold() == label.casefold()}
        if not bc_labels.get(label) or aliases != {label}:
            raise ValueError(f"Missing or ambiguous boundary label {label!r}")
    for element in ngmesh.Elements2D():
        fd = element.index
        if not 1 <= fd <= len(fd_table):
            raise ValueError(f"Unknown surface FD {fd}")
        label = fd_table[fd - 1]["label"]
        if label in points:
            points[label].update(_point_nr(v) for v in element.vertices)
    if not points[master] or not points[slave] or points[master] & points[slave]:
        raise ValueError("Empty or overlapping master/slave point sets")
    return points[master], points[slave], fd_table


def verify_kelvin_identification(mesh, *, master, slave, group=1):
    """Require directed complete native periodic coverage without modifying mesh.

    The name cannot be read through the pinned Netgen Python API. This function
    therefore never claims persisted-name proof, including for OCC input meshes.
    Other groups may exist; only the explicitly selected group is verified here.
    Case-specific sole-group, pair-oracle and curved-map checks remain separate.
    """
    from netgen.meshing import IdentificationType

    if type(group) is not int or group < 1:
        raise ValueError("Positive explicit identification group required")
    ngmesh = mesh.ngmesh
    if ngmesh.GetNrIdentifications() < group:
        raise ValueError(f"Missing persisted group {group}")
    if ngmesh.GetIdentificationType(group) != IdentificationType.PERIODIC:
        raise ValueError(f"Group {group}: expected PERIODIC")
    master_points, slave_points, table = _boundary_points(mesh, master, slave)
    rows = _native_rows(ngmesh)
    selected = [(a, b) for a, b, nr in rows if nr == group]
    masters = [a for a, b in selected]
    slaves = [b for a, b in selected]
    if any(a == b for a, b in selected):
        raise ValueError("Self identification is unsupported")
    if len(masters) != len(set(masters)) or len(slaves) != len(set(slaves)):
        raise ValueError("Identification is not a directed bijection")
    if set(masters) != master_points or set(slaves) != slave_points:
        raise ValueError(
            f"Incomplete/foreign/reversed periodic coverage: actual master/slave "
            f"{len(set(masters))}/{len(set(slaves))}, expected "
            f"{len(master_points)}/{len(slave_points)}; "
            f"missing={sorted(master_points-set(masters))[:8]}/"
            f"{sorted(slave_points-set(slaves))[:8]}, "
            f"foreign={sorted(set(masters)-master_points)[:8]}/"
            f"{sorted(set(slaves)-slave_points)[:8]}")
    return {"verified": True, "group": group, "type": "PERIODIC", "n_pairs": len(selected),
            "master": master, "slave": slave, "master_points": len(master_points),
            "slave_points": len(slave_points), "fd_bc_table": table,
            "point_index_convention": "Netgen PointId.nr one-based",
            "name_status": "NOT_VERIFIABLE_FROM_THIS_API", "mesh_modified": False}


def _verify_kelvin_structure(mesh, *, boundaries=("kelvin_int", "kelvin_ext"), group=1):
    """Internal in-memory/OCC check: observe native direction, not provenance.

    Neither orientation is presumed by an arbitrary in-memory mesh. Require one
    consistent complete directed relation, report its actual orientation, and
    leave expected producer orientation/name proof to the strict file loader.
    """
    first, second = boundaries
    left, right, _ = _boundary_points(mesh, first, second)
    selected = [(a, b) for a, b, nr in _native_rows(mesh.ngmesh) if nr == group]
    starts, ends = {a for a, b in selected}, {b for a, b in selected}
    if starts == left and ends == right:
        master, slave = first, second
    elif starts == right and ends == left:
        master, slave = second, first
    else:
        raise ValueError("Incomplete/foreign/inconsistently directed native Kelvin structure")
    evidence = verify_kelvin_identification(mesh, master=master, slave=slave, group=group)
    evidence["direction_status"] = "OBSERVED_NATIVE_NOT_PRODUCER_ORIENTATION_PROOF"
    return evidence


def _persisted_sections(raw):
    lines = raw.decode("utf-8", errors="strict").splitlines()
    sections = {}
    keys = {"identifications", "identificationtypes", "identificationnames"}
    i = 0
    while i < len(lines):
        key = lines[i]
        i += 1
        if key not in keys:
            continue
        if key in sections or i >= len(lines) or not lines[i].strip().isdigit():
            raise ValueError("Duplicate/malformed periodic section")
        count = int(lines[i]); i += 1
        nrows = count if key == "identifications" else 1
        if count < 1 or i + nrows > len(lines):
            raise ValueError("Empty/truncated periodic section")
        block = [line.split() for line in lines[i:i+nrows]]; i += nrows
        if key == "identifications":
            if any(len(row) != 3 or any(not x.isdigit() or int(x)<1 for x in row) for row in block):
                raise ValueError("Malformed persisted identification triple")
            value = [tuple(map(int, row)) for row in block]
            if len(value) != len(set(value)):
                raise ValueError("Duplicate persisted identification triple")
        else:
            if len(block[0]) != count:
                raise ValueError("Malformed persisted type/name count")
            value = block[0]
        sections[key] = value
    if set(sections) != keys:
        raise ValueError("Missing persisted identification/name/type section")
    return sections


def load_verified_kelvin_mesh(path, master, slave, group, name):
    """Load exactly one immutable .vol and prove its named producer contract.

    Returns (loaded_mesh, evidence). The evidence is never accepted as proof for
    some unrelated in-memory mesh. No Save, Curve or identification repair occurs.
    """
    from ngsolve import Mesh

    if type(group) is not int or group < 1 or not name or len(name.split()) != 1:
        raise ValueError("Positive group and explicit single-token name required")
    path = Path(path).resolve(strict=True)
    if path.suffix != ".vol":
        raise ValueError("Named periodic verification requires an uncompressed .vol")
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    sections = _persisted_sections(raw)
    names, types = sections["identificationnames"], sections["identificationtypes"]
    if len(names) != len(types) or group > len(names) or names[group-1] != name or types[group-1] != "2":
        raise ValueError(f"Expected persisted PERIODIC group {group} name {name!r}")
    if names.count(name) != 1:
        raise ValueError("Ambiguous persisted identification name")
    rows = sections["identifications"]
    if any(nr > len(names) for a, b, nr in rows):
        raise ValueError("Persisted row refers to unknown group")
    mesh = Mesh(str(path))
    if path.read_bytes() != raw:
        raise ValueError("Input changed while loading named periodic mesh")
    evidence = verify_kelvin_identification(mesh, master=master, slave=slave, group=group)
    if sorted(_native_rows(mesh.ngmesh)) != sorted(rows) or mesh.ngmesh.GetNrIdentifications() != len(names):
        raise ValueError("Native periodic rows differ from persisted file")
    if path.read_bytes() != raw:
        raise ValueError("Input changed during periodic verification")
    evidence.update(path=str(path), input_sha256=digest, name=name,
                    name_status="VERIFIED_PERSISTED_FILE_AND_NATIVE_ROWS", input_unchanged=True)
    return mesh, evidence
