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

    def trace_node_ids_by_boundary_number(self) -> dict[int, tuple[int, ...]]:
        """Return one-based trace node ids grouped by Netgen boundary number."""

        groups: dict[int, set[int]] = {}
        for tri in self.surface_triangles:
            groups.setdefault(tri.bcnr, set()).update(tri.nodes)
        return {bcnr: tuple(sorted(nodes)) for bcnr, nodes in sorted(groups.items())}

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

    def tetrahedron_edge_lengths(self) -> tuple[tuple[float, ...], ...]:
        """Return the six edge lengths of each tetrahedron.

        The local edge order is ``(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)``
        in the tetrahedron's stored one-based node order.
        """

        local_pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
        rows: list[tuple[float, ...]] = []
        for tet in self.tetrahedra:
            coords = [self.points[node - 1] for node in tet.nodes]
            rows.append(tuple(_norm(_sub(coords[j], coords[i])) for i, j in local_pairs))
        return tuple(rows)

    def tetrahedron_edge_length_ratios(self) -> tuple[float, ...]:
        """Return ``max(edge_length) / min(edge_length)`` for each tetrahedron."""

        ratios: list[float] = []
        for lengths in self.tetrahedron_edge_lengths():
            min_len = min(lengths)
            if min_len <= 0.0:
                raise ValueError("tetrahedron contains a zero-length edge")
            ratios.append(max(lengths) / min_len)
        return tuple(ratios)

    def tetrahedron_edge_length_summary(self) -> dict[str, float | int | None]:
        """Return compact edge-length and edge-ratio summary statistics."""

        rows = self.tetrahedron_edge_lengths()
        if not rows:
            return {
                "tetrahedra": 0,
                "min_edge": None,
                "max_edge": None,
                "mean_edge": None,
                "max_edge_ratio": None,
            }
        lengths = [value for row in rows for value in row]
        ratios = self.tetrahedron_edge_length_ratios()
        return {
            "tetrahedra": len(rows),
            "min_edge": min(lengths),
            "max_edge": max(lengths),
            "mean_edge": sum(lengths) / len(lengths),
            "max_edge_ratio": max(ratios),
        }

    def tetrahedron_quality_rows(self) -> tuple[dict[str, float | int], ...]:
        """Return per-tetrahedron volume/radius quality records.

        ``radius_ratio_quality = 3 * inradius / circumradius`` is one for an
        equilateral tetrahedron and tends to zero for flat/sliver elements.
        ``corner_normalized_jacobians`` are the absolute triple product of the
        three incident edges at each corner divided by the product of their
        lengths.  They are one for a locally orthogonal corner and approach zero
        as the corner becomes flat.  Together these form a compact, public-safe
        companion to Cubit/Coreform tet quality checks.
        """

        rows: list[dict[str, float | int]] = []
        edge_rows = self.tetrahedron_edge_lengths()
        signed_volumes = self.tetrahedron_signed_volumes()
        for index, (tet, edge_lengths, signed_volume) in enumerate(
            zip(self.tetrahedra, edge_rows, signed_volumes), start=1
        ):
            coords = tuple(self.points[node - 1] for node in tet.nodes)
            volume = abs(signed_volume)
            if volume <= 0.0:
                raise ValueError(f"tetrahedron {index} has zero volume")
            surface_area = _tetrahedron_surface_area(coords)
            inradius = 3.0 * volume / surface_area
            circumradius = _tetrahedron_circumradius(coords)
            if circumradius <= 0.0:
                raise ValueError(f"tetrahedron {index} has zero circumradius")
            min_edge = min(edge_lengths)
            max_edge = max(edge_lengths)
            corner_jacobians = _tetrahedron_corner_normalized_jacobians(coords)
            rows.append({
                "tetrahedron": index,
                "volume": volume,
                "surface_area": surface_area,
                "inradius": inradius,
                "circumradius": circumradius,
                "radius_ratio_quality": 3.0 * inradius / circumradius,
                "corner_normalized_jacobians": list(corner_jacobians),
                "min_normalized_corner_jacobian": min(corner_jacobians),
                "max_normalized_corner_jacobian": max(corner_jacobians),
                "min_edge": min_edge,
                "max_edge": max_edge,
                "edge_ratio": max_edge / min_edge,
            })
        return tuple(rows)

    def tetrahedron_quality_summary(self) -> dict[str, float | int | None]:
        """Return compact radius-ratio and edge-ratio quality statistics."""

        rows = self.tetrahedron_quality_rows()
        if not rows:
            return {
                "tetrahedra": 0,
                "min_radius_ratio_quality": None,
                "max_radius_ratio_quality": None,
                "mean_radius_ratio_quality": None,
                "min_inradius": None,
                "max_circumradius": None,
                "min_normalized_corner_jacobian": None,
                "mean_min_normalized_corner_jacobian": None,
                "max_edge_ratio": None,
            }
        qualities = [float(row["radius_ratio_quality"]) for row in rows]
        corner_minima = [float(row["min_normalized_corner_jacobian"]) for row in rows]
        inradii = [float(row["inradius"]) for row in rows]
        circumradii = [float(row["circumradius"]) for row in rows]
        edge_ratios = [float(row["edge_ratio"]) for row in rows]
        return {
            "tetrahedra": len(rows),
            "min_radius_ratio_quality": min(qualities),
            "max_radius_ratio_quality": max(qualities),
            "mean_radius_ratio_quality": sum(qualities) / len(qualities),
            "min_inradius": min(inradii),
            "max_circumradius": max(circumradii),
            "min_normalized_corner_jacobian": min(corner_minima),
            "mean_min_normalized_corner_jacobian": sum(corner_minima) / len(corner_minima),
            "max_edge_ratio": max(edge_ratios),
        }

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

    def surface_triangle_edge_lengths(self) -> tuple[tuple[float, float, float], ...]:
        """Return the three edge lengths of each boundary triangle."""

        rows: list[tuple[float, float, float]] = []
        for tri in self.surface_triangles:
            a, b, c = (self.points[node - 1] for node in tri.nodes)
            rows.append((_norm(_sub(b, a)), _norm(_sub(c, b)), _norm(_sub(a, c))))
        return tuple(rows)

    def surface_triangle_quality_rows(self) -> tuple[dict[str, float | int], ...]:
        """Return per-boundary-triangle quality records.

        ``radius_ratio_quality = 2 * inradius / circumradius`` is one for an
        equilateral triangle and tends to zero for sliver triangles.  This is a
        compact, solver-neutral companion to Cubit/Coreform surface-mesh quality
        checks before a boundary mesh is used as a scalar-BEM or RWG trace.
        """

        rows: list[dict[str, float | int]] = []
        for index, (tri, edge_lengths, area) in enumerate(
            zip(self.surface_triangles, self.surface_triangle_edge_lengths(), self.surface_triangle_areas()),
            start=1,
        ):
            if area <= 0.0:
                raise ValueError(f"surface triangle {index} has zero area")
            min_edge = min(edge_lengths)
            max_edge = max(edge_lengths)
            if min_edge <= 0.0:
                raise ValueError(f"surface triangle {index} contains a zero-length edge")
            semiperimeter = 0.5 * sum(edge_lengths)
            inradius = area / semiperimeter
            circumradius = edge_lengths[0] * edge_lengths[1] * edge_lengths[2] / (4.0 * area)
            angles = _triangle_angles_degrees(edge_lengths)
            rows.append({
                "surface_triangle": index,
                "boundary_number": tri.bcnr,
                "surfnr": tri.surfnr,
                "area": area,
                "inradius": inradius,
                "circumradius": circumradius,
                "radius_ratio_quality": 2.0 * inradius / circumradius,
                "min_edge": min_edge,
                "max_edge": max_edge,
                "edge_ratio": max_edge / min_edge,
                "min_angle_deg": min(angles),
                "max_angle_deg": max(angles),
            })
        return tuple(rows)

    def surface_triangle_quality_summary(self) -> dict[str, float | int | None]:
        """Return compact boundary-triangle quality statistics."""

        rows = self.surface_triangle_quality_rows()
        if not rows:
            return {
                "surface_triangles": 0,
                "min_radius_ratio_quality": None,
                "max_radius_ratio_quality": None,
                "mean_radius_ratio_quality": None,
                "min_area": None,
                "max_edge_ratio": None,
                "min_angle_deg": None,
                "max_angle_deg": None,
            }
        qualities = [float(row["radius_ratio_quality"]) for row in rows]
        areas = [float(row["area"]) for row in rows]
        edge_ratios = [float(row["edge_ratio"]) for row in rows]
        min_angles = [float(row["min_angle_deg"]) for row in rows]
        max_angles = [float(row["max_angle_deg"]) for row in rows]
        return {
            "surface_triangles": len(rows),
            "min_radius_ratio_quality": min(qualities),
            "max_radius_ratio_quality": max(qualities),
            "mean_radius_ratio_quality": sum(qualities) / len(qualities),
            "min_area": min(areas),
            "max_edge_ratio": max(edge_ratios),
            "min_angle_deg": min(min_angles),
            "max_angle_deg": max(max_angles),
        }

    def total_surface_area(self) -> float:
        """Return the sum of boundary triangle areas."""

        return sum(self.surface_triangle_areas())

    def surface_area_by_boundary_number(self) -> dict[int, float]:
        """Return boundary surface area grouped by Netgen boundary number."""

        areas: dict[int, float] = {}
        for tri, area in zip(self.surface_triangles, self.surface_triangle_areas()):
            areas[tri.bcnr] = areas.get(tri.bcnr, 0.0) + area
        return areas

    def boundary_summary_rows(self) -> tuple[dict[str, object], ...]:
        """Return named boundary inventory rows for FEM/BEM conditions.

        Coreform/Cubit sidesets and Netgen boundary names arrive as integer
        boundary numbers plus optional names.  This compact table keeps the
        per-boundary triangle count, area, and one-based trace nodes together
        so teaching scripts can build readable boundary-condition maps without
        re-walking the raw triangle records.
        """

        areas = self.surface_area_by_boundary_number()
        nodes = self.trace_node_ids_by_boundary_number()
        triangle_counts: dict[int, int] = {}
        for tri in self.surface_triangles:
            triangle_counts[tri.bcnr] = triangle_counts.get(tri.bcnr, 0) + 1
        boundary_numbers = sorted(
            set(self.boundary_names) | set(areas) | set(nodes) | set(triangle_counts)
        )
        rows: list[dict[str, object]] = []
        for bcnr in boundary_numbers:
            trace_nodes = nodes.get(bcnr, ())
            rows.append({
                "boundary_number": bcnr,
                "name": self.boundary_names.get(bcnr, f"boundary_{bcnr}"),
                "surface_triangles": triangle_counts.get(bcnr, 0),
                "surface_area": areas.get(bcnr, 0.0),
                "trace_node_count": len(trace_nodes),
                "trace_node_ids": list(trace_nodes),
            })
        return tuple(rows)

    def domain_boundary_incidence_rows(self) -> tuple[dict[str, object], ...]:
        """Return boundary rows grouped by adjacent volume domains.

        Netgen ``surfaceelements`` carry ``domin`` and ``domout`` integers.  For
        Coreform/Cubit exports these are the small but important bridge from a
        named sideset to the volume material(s) it touches: ``domout == 0`` is an
        exterior boundary, while two nonzero domains mark an interface.  Keeping
        this as a readable table makes FEM material assignment and BEM/interface
        operator setup share the same one-based node ids.
        """

        groups: dict[tuple[int, int, int], list[tuple[NetgenSurfaceTriangle, float]]] = {}
        for tri, area in zip(self.surface_triangles, self.surface_triangle_areas()):
            groups.setdefault((tri.bcnr, tri.domin, tri.domout), []).append((tri, area))

        rows: list[dict[str, object]] = []
        for (bcnr, domin, domout), entries in sorted(groups.items()):
            node_ids = sorted({node for tri, _area in entries for node in tri.nodes})
            rows.append({
                "boundary_number": bcnr,
                "name": self.boundary_names.get(bcnr, f"boundary_{bcnr}"),
                "domin": domin,
                "domout": domout,
                "domin_material": self.materials.get(domin) if domin else None,
                "domout_material": self.materials.get(domout) if domout else None,
                "kind": "exterior" if 0 in (domin, domout) else "interface",
                "surface_triangles": len(entries),
                "surface_area": sum(area for _tri, area in entries),
                "trace_node_count": len(node_ids),
                "trace_node_ids": node_ids,
            })
        return tuple(rows)

    def material_summary_rows(self) -> tuple[dict[str, object], ...]:
        """Return volume/material inventory with boundary incidence.

        The rows are keyed by ``volumeelements`` material number and include
        tetrahedron count, volume, node ids, touching boundary names, exterior
        area, and interface area.  This is the volume-side companion to
        :meth:`boundary_summary_rows`.
        """

        signed_volumes = self.tetrahedron_signed_volumes()
        material_numbers = sorted(
            set(self.materials) | {tet.matnr for tet in self.tetrahedra}
        )
        total = self.total_volume()
        rows: list[dict[str, object]] = []
        for matnr in material_numbers:
            tet_ids = [
                index
                for index, tet in enumerate(self.tetrahedra, start=1)
                if tet.matnr == matnr
            ]
            node_ids = sorted({
                node
                for tet_id in tet_ids
                for node in self.tetrahedra[tet_id - 1].nodes
            })
            volume = sum(abs(signed_volumes[tet_id - 1]) for tet_id in tet_ids)
            boundary_numbers: set[int] = set()
            exterior_boundary_numbers: set[int] = set()
            interface_boundary_numbers: set[int] = set()
            neighboring_material_numbers: set[int] = set()
            exterior_area = 0.0
            interface_area = 0.0
            for tri, area in zip(self.surface_triangles, self.surface_triangle_areas()):
                if matnr not in (tri.domin, tri.domout):
                    continue
                boundary_numbers.add(tri.bcnr)
                other = tri.domout if tri.domin == matnr else tri.domin
                if other == 0:
                    exterior_boundary_numbers.add(tri.bcnr)
                    exterior_area += area
                else:
                    interface_boundary_numbers.add(tri.bcnr)
                    neighboring_material_numbers.add(other)
                    interface_area += area

            rows.append({
                "material_number": matnr,
                "name": self.materials.get(matnr, f"material_{matnr}"),
                "tetrahedra": len(tet_ids),
                "tetrahedron_ids": tet_ids,
                "volume": volume,
                "volume_fraction": volume / total if total > 0.0 else None,
                "node_count": len(node_ids),
                "node_ids": node_ids,
                "boundary_numbers": sorted(boundary_numbers),
                "boundary_names": [
                    self.boundary_names.get(bcnr, f"boundary_{bcnr}")
                    for bcnr in sorted(boundary_numbers)
                ],
                "exterior_boundary_numbers": sorted(exterior_boundary_numbers),
                "interface_boundary_numbers": sorted(interface_boundary_numbers),
                "neighboring_material_numbers": sorted(neighboring_material_numbers),
                "exterior_surface_area": exterior_area,
                "interface_surface_area": interface_area,
            })
        return tuple(rows)

    def surface_triangle_area_vectors(self) -> tuple[tuple[float, float, float], ...]:
        """Return oriented boundary triangle area vectors.

        Each vector is ``0.5 * cross(b-a, c-a)`` in the stored node order. For a
        consistently oriented closed surface, the sum of these vectors should be
        close to zero. This is a useful early gate for BEM normal conventions.
        """

        vectors: list[tuple[float, float, float]] = []
        for tri in self.surface_triangles:
            a, b, c = (self.points[node - 1] for node in tri.nodes)
            vectors.append(_scale(_cross(_sub(b, a), _sub(c, a)), 0.5))
        return tuple(vectors)

    def surface_vector_area(self) -> tuple[float, float, float]:
        """Return the vector sum of oriented boundary triangle areas."""

        total = (0.0, 0.0, 0.0)
        for vec in self.surface_triangle_area_vectors():
            total = _add(total, vec)
        return total

    def surface_signed_volume_from_triangles(self) -> float:
        """Return oriented volume enclosed by boundary triangles.

        The divergence-theorem formula ``sum(dot(a, cross(b, c))/6)`` is
        positive for outward-oriented closed surfaces and negative for inward
        orientation. Its absolute value should match ``total_volume()`` when
        the surface triangles close the tetrahedral volume.
        """

        volume = 0.0
        for tri in self.surface_triangles:
            a, b, c = (self.points[node - 1] for node in tri.nodes)
            volume += _dot(a, _cross(b, c)) / 6.0
        return volume

    def surface_closure_summary(self) -> dict[str, object]:
        """Return closure/orientation checks for the boundary triangle mesh."""

        total_area = self.total_surface_area()
        vector_area = self.surface_vector_area()
        vector_norm = _norm(vector_area)
        signed_surface_volume = self.surface_signed_volume_from_triangles()
        tet_volume = self.total_volume()
        abs_volume_error = abs(abs(signed_surface_volume) - tet_volume)
        rel_volume_error = abs_volume_error / tet_volume if tet_volume > 0.0 else None
        if signed_surface_volume > 0.0:
            orientation = "outward"
        elif signed_surface_volume < 0.0:
            orientation = "inward"
        else:
            orientation = "zero_or_open"
        return {
            "surface_triangles": len(self.surface_triangles),
            "tetrahedra": len(self.tetrahedra),
            "total_surface_area": total_area,
            "surface_vector_area": vector_area,
            "surface_vector_area_norm": vector_norm,
            "surface_vector_area_norm_over_area": (
                vector_norm / total_area if total_area > 0.0 else None
            ),
            "surface_signed_volume": signed_surface_volume,
            "surface_abs_volume": abs(signed_surface_volume),
            "tetrahedron_total_volume": tet_volume,
            "surface_abs_volume_error": abs_volume_error,
            "surface_abs_volume_rel_error": rel_volume_error,
            "boundary_orientation": orientation,
        }

    def surface_edge_manifold_summary(self) -> dict[str, object]:
        """Return boundary-edge adjacency and Euler-characteristic checks.

        A closed first-order RWG surface should have every boundary edge shared
        by exactly two triangles.  The Euler characteristic is a compact sanity
        check for simple connected surfaces: a triangulated sphere-like boundary
        has ``V - E + F = 2``.
        """

        trace_nodes = self.trace_node_ids()
        if not self.surface_triangles:
            return {
                "trace_nodes": len(trace_nodes),
                "surface_edges": 0,
                "surface_triangles": 0,
                "closed_edges": 0,
                "open_edges": 0,
                "is_closed_manifold": False,
                "euler_characteristic": len(trace_nodes),
            }
        trace_node_to_local = {node: i for i, node in enumerate(trace_nodes, start=1)}
        boundary_triangles = tuple(
            tuple(trace_node_to_local[node] for node in tri.nodes) for tri in self.surface_triangles
        )
        surface_edges, _tri_edges, _tri_edge_signs, edge_triangles, _opposites = _build_tri_edges(
            boundary_triangles
        )
        adjacency_counts = tuple(sum(1 for tri_id in adjacent if tri_id) for adjacent in edge_triangles)
        closed_edges = sum(1 for count in adjacency_counts if count == 2)
        open_edges = sum(1 for count in adjacency_counts if count == 1)
        return {
            "trace_nodes": len(trace_nodes),
            "surface_edges": len(surface_edges),
            "surface_triangles": len(self.surface_triangles),
            "closed_edges": closed_edges,
            "open_edges": open_edges,
            "is_closed_manifold": open_edges == 0 and closed_edges == len(surface_edges),
            "euler_characteristic": len(trace_nodes) - len(surface_edges) + len(self.surface_triangles),
        }

    def surface_connected_components(self) -> tuple[dict[str, object], ...]:
        """Return connected boundary-triangle components for BEM block setup.

        Components are built by triangle edge adjacency on the trace surface.
        This is the small-readable counterpart to a BEM mesh partition: one
        component may be an exterior truncation surface, another a disconnected
        body, and each can carry its own boundary condition/operator block.
        """

        if not self.surface_triangles:
            return ()

        trace_nodes = self.trace_node_ids()
        trace_node_to_local = {node: i for i, node in enumerate(trace_nodes, start=1)}
        boundary_triangles = tuple(
            tuple(trace_node_to_local[node] for node in tri.nodes) for tri in self.surface_triangles
        )
        surface_edges, tri_edges, _tri_edge_signs, edge_triangles, _opposites = _build_tri_edges(
            boundary_triangles
        )

        adjacency: list[set[int]] = [set() for _ in self.surface_triangles]
        for adjacent in edge_triangles:
            tri_ids = [tri_id for tri_id in adjacent if tri_id]
            for tri_id in tri_ids:
                adjacency[tri_id - 1].update(other for other in tri_ids if other != tri_id)

        seen: set[int] = set()
        components: list[list[int]] = []
        for tri_id in range(1, len(self.surface_triangles) + 1):
            if tri_id in seen:
                continue
            stack = [tri_id]
            seen.add(tri_id)
            component: list[int] = []
            while stack:
                current = stack.pop()
                component.append(current)
                for other in sorted(adjacency[current - 1]):
                    if other not in seen:
                        seen.add(other)
                        stack.append(other)
            components.append(sorted(component))

        areas = self.surface_triangle_areas()
        area_vectors = self.surface_triangle_area_vectors()
        rows: list[dict[str, object]] = []
        for comp_index, tri_ids in enumerate(sorted(components, key=lambda row: row[0]), start=1):
            tri_set = set(tri_ids)
            node_ids = sorted({
                node
                for tri_id in tri_ids
                for node in self.surface_triangles[tri_id - 1].nodes
            })
            boundary_numbers = sorted({
                self.surface_triangles[tri_id - 1].bcnr for tri_id in tri_ids
            })
            edge_ids = sorted({
                edge_id
                for tri_id in tri_ids
                for edge_id in tri_edges[tri_id - 1]
            })
            adjacency_counts = [
                sum(1 for adjacent_tri_id in edge_triangles[edge_id - 1] if adjacent_tri_id in tri_set)
                for edge_id in edge_ids
            ]
            closed_edges = sum(1 for count in adjacency_counts if count == 2)
            open_edges = sum(1 for count in adjacency_counts if count == 1)
            vector_area = (0.0, 0.0, 0.0)
            signed_volume = 0.0
            surface_area = 0.0
            for tri_id in tri_ids:
                tri = self.surface_triangles[tri_id - 1]
                a, b, c = (self.points[node - 1] for node in tri.nodes)
                vector_area = _add(vector_area, area_vectors[tri_id - 1])
                signed_volume += _dot(a, _cross(b, c)) / 6.0
                surface_area += areas[tri_id - 1]

            rows.append({
                "component": comp_index,
                "surface_triangles": len(tri_ids),
                "triangle_ids": tri_ids,
                "surface_edges": len(edge_ids),
                "trace_node_count": len(node_ids),
                "trace_node_ids": node_ids,
                "boundary_numbers": boundary_numbers,
                "boundary_names": [
                    self.boundary_names.get(bcnr, f"boundary_{bcnr}")
                    for bcnr in boundary_numbers
                ],
                "surface_area": surface_area,
                "surface_vector_area": vector_area,
                "surface_vector_area_norm": _norm(vector_area),
                "surface_signed_volume": signed_volume,
                "surface_abs_volume": abs(signed_volume),
                "closed_edges": closed_edges,
                "open_edges": open_edges,
                "is_closed_manifold": open_edges == 0 and closed_edges == len(edge_ids),
                "euler_characteristic": len(node_ids) - len(edge_ids) + len(tri_ids),
            })
        return tuple(rows)

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
            "trace_node_ids_by_boundary_number": {
                bcnr: list(nodes) for bcnr, nodes in self.trace_node_ids_by_boundary_number().items()
            },
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
        elif key_lower in {"surfaceelements", "surfaceelementsuv"}:
            count = next_count("surfaceelements")
            for i in range(1, count + 1):
                record = next_required_line("surfaceelements").split()
                if len(record) < 6:
                    raise ValueError("surface element line is too short")
                head = [int(v) for v in record[:5]]
                node_count = head[4]
                nodes = tuple(int(v) for v in record[5 : 5 + node_count])
                if len(nodes) != node_count:
                    raise ValueError("surface element node count does not match record")
                if node_count != 3:
                    raise ValueError(
                        f"Netgen .vol tri/tet-only policy rejected surface element {i}: "
                        f"expected 3 nodes, got {node_count}"
                    )
                surface_triangles.append(NetgenSurfaceTriangle(head[0], head[1], head[2], head[3], nodes))
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


def _add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _scale(a: tuple[float, float, float], value: float) -> tuple[float, float, float]:
    return (value * a[0], value * a[1], value * a[2])


def _norm(a: tuple[float, float, float]) -> float:
    return math.sqrt(_dot(a, a))


def _triangle_area(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
) -> float:
    return 0.5 * _norm(_cross(_sub(b, a), _sub(c, a)))


def _triangle_angles_degrees(edge_lengths: tuple[float, float, float]) -> tuple[float, float, float]:
    angles: list[float] = []
    for i, opposite in enumerate(edge_lengths):
        adjacent = [edge_lengths[j] for j in range(3) if j != i]
        denom = 2.0 * adjacent[0] * adjacent[1]
        if denom <= 0.0:
            raise ValueError("triangle contains a zero-length edge")
        value = (adjacent[0] ** 2 + adjacent[1] ** 2 - opposite ** 2) / denom
        value = max(-1.0, min(1.0, value))
        angles.append(math.degrees(math.acos(value)))
    return (angles[0], angles[1], angles[2])


def _tetrahedron_surface_area(points: tuple[tuple[float, float, float], ...]) -> float:
    a, b, c, d = points
    return (
        _triangle_area(a, b, c)
        + _triangle_area(a, b, d)
        + _triangle_area(a, c, d)
        + _triangle_area(b, c, d)
    )


def _det3(rows: tuple[tuple[float, float, float], ...]) -> float:
    return (
        rows[0][0] * (rows[1][1] * rows[2][2] - rows[1][2] * rows[2][1])
        - rows[0][1] * (rows[1][0] * rows[2][2] - rows[1][2] * rows[2][0])
        + rows[0][2] * (rows[1][0] * rows[2][1] - rows[1][1] * rows[2][0])
    )


def _solve3_rows(
    rows: tuple[tuple[float, float, float], ...],
    rhs: tuple[float, float, float],
) -> tuple[float, float, float]:
    det = _det3(rows)
    if det == 0.0:
        raise ValueError("singular 3x3 system")
    cols = []
    for col in range(3):
        replaced = tuple(
            tuple(rhs[i] if j == col else rows[i][j] for j in range(3))
            for i in range(3)
        )
        cols.append(_det3(replaced) / det)
    return (cols[0], cols[1], cols[2])


def _tetrahedron_circumradius(points: tuple[tuple[float, float, float], ...]) -> float:
    a, b, c, d = points
    edges = (_sub(b, a), _sub(c, a), _sub(d, a))
    rows = tuple(tuple(2.0 * value for value in edge) for edge in edges)
    rhs = tuple(_dot(edge, edge) for edge in edges)
    center_rel = _solve3_rows(rows, rhs)
    return _norm(center_rel)


def _tetrahedron_corner_normalized_jacobians(
    points: tuple[tuple[float, float, float], ...],
) -> tuple[float, float, float, float]:
    rows: list[float] = []
    for i, origin in enumerate(points):
        edges = tuple(_sub(points[j], origin) for j in range(4) if j != i)
        denom = 1.0
        for edge in edges:
            length = _norm(edge)
            if length <= 0.0:
                raise ValueError("tetrahedron contains a zero-length edge")
            denom *= length
        rows.append(abs(_det3(edges)) / denom)
    return (rows[0], rows[1], rows[2], rows[3])
