"""Small GMSH ASCII reader for PEEC examples.

The PEEC demos only need fixture data: nodes, line/triangle elements, and
optional PhysicalNames.  Reading that subset directly keeps public examples off
the pip ``gmsh`` runtime while preserving the existing ``.msh`` fixtures.
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


@dataclass(frozen=True)
class GmshElement:
    tag: int
    element_type: int
    tags: tuple[int, ...]
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
            if elem.physical_tag is not None and elem.dim >= 0:
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


def _parse_physical_names(lines: list[str]) -> dict[tuple[int, int], str]:
    if "$PhysicalNames" not in lines:
        return {}
    start = lines.index("$PhysicalNames") + 1
    count = int(lines[start])
    names: dict[tuple[int, int], str] = {}
    for row in lines[start + 1:start + 1 + count]:
        dim_text, tag_text, name_text = row.split(maxsplit=2)
        names[(int(dim_text), int(tag_text))] = name_text.strip().strip('"')
    return names


def read_gmsh_ascii(path: str | Path) -> GmshAsciiMesh:
    """Read GMSH v2.2 or legacy Netgen ``$NOD/$ELM`` ASCII meshes."""
    lines = Path(path).read_text(errors="replace").splitlines()
    physical_names = _parse_physical_names(lines)

    if "$Nodes" in lines:
        node_start = lines.index("$Nodes") + 1
        n_nodes = int(lines[node_start])
        nodes = []
        for row in lines[node_start + 1:node_start + 1 + n_nodes]:
            parts = row.split()
            nodes.append((int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])))

        elem_start = lines.index("$Elements") + 1
        n_elems = int(lines[elem_start])
        elements = []
        for row in lines[elem_start + 1:elem_start + 1 + n_elems]:
            parts = row.split()
            elem_tag = int(parts[0])
            elem_type = int(parts[1])
            n_tags = int(parts[2])
            tags = tuple(int(v) for v in parts[3:3 + n_tags])
            elem_nodes = tuple(int(v) for v in parts[3 + n_tags:])
            elements.append(GmshElement(elem_tag, elem_type, tags, elem_nodes))
        return GmshAsciiMesh(nodes, elements, physical_names)

    if "$NOD" in lines:
        node_start = lines.index("$NOD") + 1
        n_nodes = int(lines[node_start])
        nodes = []
        for row in lines[node_start + 1:node_start + 1 + n_nodes]:
            parts = row.split()
            nodes.append((int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])))

        elem_start = lines.index("$ELM") + 1
        n_elems = int(lines[elem_start])
        elements = []
        for row in lines[elem_start + 1:elem_start + 1 + n_elems]:
            parts = row.split()
            elem_tag = int(parts[0])
            elem_type = int(parts[1])
            physical_tag = int(parts[2])
            elementary_tag = int(parts[3])
            n_elem_nodes = int(parts[4])
            elem_nodes = tuple(int(v) for v in parts[5:5 + n_elem_nodes])
            elements.append(
                GmshElement(elem_tag, elem_type, (physical_tag, elementary_tag), elem_nodes)
            )
        return GmshAsciiMesh(nodes, elements, physical_names)

    raise ValueError(f"Unsupported GMSH ASCII mesh format: {path}")
