"""Small GMSH MSH v4.1 ASCII reader for PEEC examples.

The PEEC demos only need fixture data: nodes, point/line/triangle/quad
elements, and optional PhysicalNames.  Reading that subset directly keeps
public examples off the pip ``gmsh`` runtime while preserving the existing
``.msh`` fixtures.

Scope: MSH **v4.1 ASCII** only, matching the repository-wide format policy
(GMSH .msh Format Version Policy, 2026-04-15).  v2.2 and binary files are
rejected loudly -- re-export the mesh instead of extending this reader.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


_DIM_BY_ELEMENT_TYPE = {
    15: 0,  # 1-node point
    1: 1,  # 2-node line
    2: 2,  # 3-node triangle
    3: 2,  # 4-node quadrangle
    4: 3,  # 4-node tetrahedron
    5: 3,  # 8-node hexahedron
}

_NODES_PER_ELEMENT_TYPE = {15: 1, 1: 2, 2: 3, 3: 4, 4: 4, 5: 8}


@dataclass(frozen=True)
class GmshElement:
    tag: int
    element_type: int
    tags: tuple[int, ...]  # (physical_tag, elementary/entity_tag)
    nodes: tuple[int, ...]

    @property
    def dim(self) -> int:
        return _DIM_BY_ELEMENT_TYPE.get(self.element_type, -1)

    @property
    def physical_tag(self) -> int | None:
        return self.tags[0] if self.tags else None

    @property
    def elementary_tag(self) -> int | None:
        return self.tags[1] if len(self.tags) > 1 else None


class GmshAsciiMesh:
    def __init__(
        self,
        nodes: list[tuple[int, float, float, float]],
        elements: list[GmshElement],
        physical_names: dict[tuple[int, int], str],
    ) -> None:
        self.node_tags = np.asarray([row[0] for row in nodes], dtype=np.int64)
        self.coords = np.asarray([row[1:] for row in nodes], dtype=float)
        self.elements = elements
        self.physical_names = physical_names
        self._node_index = {int(tag): i for i, tag in enumerate(self.node_tags)}
        self._element_by_tag = {elem.tag: elem for elem in elements}

    @property
    def node_coords_flat(self) -> np.ndarray:
        return self.coords.reshape(-1)

    def coord(self, node_tag: int) -> np.ndarray:
        return self.coords[self._node_index[int(node_tag)]]

    def grouped_elements(self) -> tuple[np.ndarray, list[np.ndarray], list[np.ndarray]]:
        """Return arrays shaped like ``gmsh.model.mesh.getElements()``."""
        by_type: dict[int, list[GmshElement]] = {}
        for elem in self.elements:
            by_type.setdefault(elem.element_type, []).append(elem)

        elem_types = np.asarray(sorted(by_type), dtype=np.int32)
        elem_tags = []
        elem_node_tags = []
        for elem_type in elem_types:
            group = by_type[int(elem_type)]
            elem_tags.append(np.asarray([elem.tag for elem in group], dtype=np.int64))
            flat_nodes = [node for elem in group for node in elem.nodes]
            elem_node_tags.append(np.asarray(flat_nodes, dtype=np.int64))
        return elem_types, elem_tags, elem_node_tags

    def elements_of_type(self, element_type: int) -> list[GmshElement]:
        return [elem for elem in self.elements if elem.element_type == element_type]

    def physical_groups(self) -> list[tuple[int, int]]:
        groups = set(self.physical_names)
        for elem in self.elements:
            if elem.physical_tag and elem.dim >= 0:
                groups.add((elem.dim, int(elem.physical_tag)))
        return sorted(groups)

    def physical_name(self, dim: int, tag: int) -> str:
        return self.physical_names.get((int(dim), int(tag)), f"physical_{dim}_{tag}")

    def elements_for_physical_group(self, dim: int, tag: int) -> list[GmshElement]:
        return [
            elem
            for elem in self.elements
            if elem.dim == int(dim) and elem.physical_tag == int(tag)
        ]

    def element(self, element_tag: int) -> GmshElement:
        return self._element_by_tag[int(element_tag)]


def _section_body(lines: list[str], name: str) -> list[str] | None:
    """Return the stripped lines between ``$name`` and ``$Endname``."""
    try:
        start = lines.index(f"${name}") + 1
    except ValueError:
        return None
    end = lines.index(f"$End{name}", start)
    return lines[start:end]


def _parse_physical_names(lines: list[str]) -> dict[tuple[int, int], str]:
    body = _section_body(lines, "PhysicalNames")
    if body is None:
        return {}
    count = int(body[0])
    names: dict[tuple[int, int], str] = {}
    for row in body[1:1 + count]:
        dim_text, tag_text, name_text = row.split(maxsplit=2)
        names[(int(dim_text), int(tag_text))] = name_text.strip().strip('"')
    return names


def _parse_entity_physicals(lines: list[str]) -> dict[tuple[int, int], int]:
    """Map (dim, entity_tag) -> physical tag (0 when untagged).

    The compact example API represents one physical tag per element, so an
    entity assigned to multiple groups must fail loudly instead of silently
    discarding all but the first tag.
    """
    body = _section_body(lines, "Entities")
    if body is None:
        return {}
    toks = iter(" ".join(body).split())
    counts = [int(next(toks)) for _ in range(4)]
    physicals: dict[tuple[int, int], int] = {}
    # points: tag x y z numPhys phys...
    for _ in range(counts[0]):
        tag = int(next(toks))
        for _ in range(3):
            next(toks)
        n_phys = int(next(toks))
        phys = [int(next(toks)) for _ in range(n_phys)]
        if len(phys) > 1:
            raise ValueError(
                f"point entity {tag} has multiple physical tags {phys}; "
                "this example reader supports exactly one")
        physicals[(0, tag)] = phys[0] if phys else 0
    # curves / surfaces / volumes: tag bbox(6) numPhys phys... numBnd bnd...
    for dim in (1, 2, 3):
        for _ in range(counts[dim]):
            tag = int(next(toks))
            for _ in range(6):
                next(toks)
            n_phys = int(next(toks))
            phys = [int(next(toks)) for _ in range(n_phys)]
            if len(phys) > 1:
                raise ValueError(
                    f"dim{dim} entity {tag} has multiple physical tags "
                    f"{phys}; this example reader supports exactly one")
            physicals[(dim, tag)] = phys[0] if phys else 0
            n_bnd = int(next(toks))
            for _ in range(n_bnd):
                next(toks)
    return physicals


def _parse_nodes(lines: list[str]) -> list[tuple[int, float, float, float]]:
    body = _section_body(lines, "Nodes")
    if body is None:
        raise ValueError("$Nodes section missing")
    n_blocks, n_nodes = (int(v) for v in body[0].split()[:2])
    nodes: list[tuple[int, float, float, float]] = []
    row = 1
    for _ in range(n_blocks):
        _dim, _etag, parametric, n_in_block = (int(v) for v in body[row].split())
        if parametric:
            raise ValueError("parametric $Nodes blocks are not supported")
        row += 1
        tags = [int(body[row + i]) for i in range(n_in_block)]
        row += n_in_block
        for i in range(n_in_block):
            x, y, z = (float(v) for v in body[row + i].split()[:3])
            nodes.append((tags[i], x, y, z))
        row += n_in_block
    if len(nodes) != n_nodes:
        raise ValueError(
            f"$Nodes declares {n_nodes} nodes but {len(nodes)} were parsed")
    return nodes


def _parse_elements(
    lines: list[str],
    entity_physicals: dict[tuple[int, int], int],
) -> list[GmshElement]:
    body = _section_body(lines, "Elements")
    if body is None:
        raise ValueError("$Elements section missing")
    n_blocks, n_elems = (int(v) for v in body[0].split()[:2])
    elements: list[GmshElement] = []
    row = 1
    for _ in range(n_blocks):
        dim, entity_tag, elem_type, n_in_block = (int(v) for v in body[row].split())
        row += 1
        if elem_type not in _NODES_PER_ELEMENT_TYPE:
            raise ValueError(
                f"unsupported element type {elem_type}; extend "
                f"_NODES_PER_ELEMENT_TYPE if this type is legitimate")
        n_per = _NODES_PER_ELEMENT_TYPE[elem_type]
        physical_tag = entity_physicals.get((dim, entity_tag), 0)
        for i in range(n_in_block):
            parts = body[row + i].split()
            elem_tag = int(parts[0])
            elem_nodes = tuple(int(v) for v in parts[1:1 + n_per])
            elements.append(
                GmshElement(elem_tag, elem_type, (physical_tag, entity_tag),
                            elem_nodes))
        row += n_in_block
    if len(elements) != n_elems:
        raise ValueError(
            f"$Elements declares {n_elems} elements but {len(elements)} "
            f"were parsed")
    return elements


def read_gmsh_ascii(path: str | Path) -> GmshAsciiMesh:
    """Read a GMSH MSH v4.1 ASCII mesh (nodes, elements, physical names)."""
    lines = [ln.strip() for ln in
             Path(path).read_text(errors="replace").splitlines()]

    fmt = _section_body(lines, "MeshFormat")
    if fmt is None:
        raise ValueError(f"$MeshFormat section missing: {path}")
    version, file_type = fmt[0].split()[:2]
    if version != "4.1" or file_type != "0":
        raise ValueError(
            f"unsupported MSH format {version!r} (file_type={file_type}) in "
            f"{path}: this reader accepts ASCII v4.1 only (repo format "
            f"policy); re-export the mesh as v4.1")

    physical_names = _parse_physical_names(lines)
    entity_physicals = _parse_entity_physicals(lines)
    nodes = _parse_nodes(lines)
    elements = _parse_elements(lines, entity_physicals)
    return GmshAsciiMesh(nodes, elements, physical_names)
