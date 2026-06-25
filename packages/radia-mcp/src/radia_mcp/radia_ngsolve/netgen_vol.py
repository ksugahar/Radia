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

    def worst_tetrahedra_by_quality(self, limit: int = 5) -> tuple[dict[str, object], ...]:
        """Return the lowest-quality tetrahedra first.

        The ordering is intentionally simple and teachable: radius-ratio
        quality first, corner-Jacobian quality second, then edge ratio as a
        tie-breaker.  These rows make a large exported mesh debuggable without
        dumping every element.
        """

        if limit < 0:
            raise ValueError("limit must be non-negative")
        rows = sorted(
            self.tetrahedron_quality_rows(),
            key=lambda row: (
                float(row["radius_ratio_quality"]),
                float(row["min_normalized_corner_jacobian"]),
                -float(row["edge_ratio"]),
                int(row["tetrahedron"]),
            ),
        )
        return tuple(rows[:limit])

    def worst_surface_triangles_by_quality(self, limit: int = 5) -> tuple[dict[str, object], ...]:
        """Return the lowest-quality boundary triangles first."""

        if limit < 0:
            raise ValueError("limit must be non-negative")
        rows = sorted(
            self.surface_triangle_quality_rows(),
            key=lambda row: (
                float(row["radius_ratio_quality"]),
                float(row["min_angle_deg"]),
                -float(row["edge_ratio"]),
                int(row["surface_triangle"]),
            ),
        )
        return tuple(rows[:limit])

    def total_surface_area(self) -> float:
        """Return the sum of boundary triangle areas."""

        return sum(self.surface_triangle_areas())

    def surface_area_by_boundary_number(self) -> dict[int, float]:
        """Return boundary surface area grouped by Netgen boundary number."""

        areas: dict[int, float] = {}
        for tri, area in zip(self.surface_triangles, self.surface_triangle_areas()):
            areas[tri.bcnr] = areas.get(tri.bcnr, 0.0) + area
        return areas

    def surface_vector_area_by_boundary_number(self) -> dict[int, tuple[float, float, float]]:
        """Return oriented boundary area vectors grouped by boundary number.

        For a planar sideset with consistently oriented triangles, the vector
        norm equals the scalar area and the unit vector is the sideset normal.
        This is the compact `.vol`-side gate needed before Maxwell-stress force
        integration over named Cubit/Coreform sidesets.
        """

        vectors: dict[int, tuple[float, float, float]] = {}
        for tri, vector in zip(self.surface_triangles, self.surface_triangle_area_vectors()):
            vectors[tri.bcnr] = _add(vectors.get(tri.bcnr, (0.0, 0.0, 0.0)), vector)
        return vectors

    def boundary_normal_summary_rows(self) -> tuple[dict[str, object], ...]:
        """Return per-boundary oriented area vector and normal summaries."""

        areas = self.surface_area_by_boundary_number()
        vectors = self.surface_vector_area_by_boundary_number()
        triangle_counts: dict[int, int] = {}
        for tri in self.surface_triangles:
            triangle_counts[tri.bcnr] = triangle_counts.get(tri.bcnr, 0) + 1
        boundary_numbers = sorted(set(self.boundary_names) | set(areas) | set(vectors) | set(triangle_counts))
        rows: list[dict[str, object]] = []
        for bcnr in boundary_numbers:
            vector = vectors.get(bcnr, (0.0, 0.0, 0.0))
            vector_norm = _norm(vector)
            scalar_area = areas.get(bcnr, 0.0)
            unit_normal = (
                tuple(component / vector_norm for component in vector)
                if vector_norm > 0.0
                else None
            )
            rows.append({
                "boundary_number": bcnr,
                "name": self.boundary_names.get(bcnr, f"boundary_{bcnr}"),
                "surface_triangles": triangle_counts.get(bcnr, 0),
                "surface_area": scalar_area,
                "vector_area": vector,
                "vector_area_norm": vector_norm,
                "vector_area_norm_over_area": (
                    vector_norm / scalar_area if scalar_area > 0.0 else None
                ),
                "unit_normal": unit_normal,
            })
        return tuple(rows)

    def _boundary_pressure_value(
        self,
        normal_row: dict[str, object],
        pressure_by_boundary: dict[int | str, float],
        default_pressure: float | None,
    ) -> tuple[float, str]:
        bcnr = int(normal_row["boundary_number"])
        name = str(normal_row["name"])
        if bcnr in pressure_by_boundary:
            return float(pressure_by_boundary[bcnr]), "boundary_number"
        if name in pressure_by_boundary:
            return float(pressure_by_boundary[name]), "name"
        if default_pressure is not None:
            return float(default_pressure), "default"
        raise KeyError(f"missing pressure for boundary {bcnr} ({name})")

    def _boundary_vector_value(
        self,
        normal_row: dict[str, object],
        vector_by_boundary: dict[int | str, tuple[float, float, float]],
        default_vector: tuple[float, float, float] | None,
        value_name: str,
    ) -> tuple[tuple[float, float, float], str]:
        bcnr = int(normal_row["boundary_number"])
        name = str(normal_row["name"])
        if bcnr in vector_by_boundary:
            value = vector_by_boundary[bcnr]
            source = "boundary_number"
        elif name in vector_by_boundary:
            value = vector_by_boundary[name]
            source = "name"
        elif default_vector is not None:
            value = default_vector
            source = "default"
        else:
            raise KeyError(f"missing {value_name} for boundary {bcnr} ({name})")
        vector = tuple(float(component) for component in value)
        if len(vector) != 3:
            raise ValueError(f"{value_name} values must have three components")
        return vector, source

    def boundary_pressure_force_rows(
        self,
        pressure_by_boundary: dict[int | str, float],
        default_pressure: float | None = 0.0,
    ) -> tuple[dict[str, object], ...]:
        """Return per-boundary force rows from scalar pressure and vector area.

        Pressure is positive along the boundary's oriented outward normal:
        ``force = pressure * vector_area``.  Keys may be Netgen boundary numbers
        or boundary names.  Set ``default_pressure=None`` to require every
        boundary to be explicitly present.
        """

        rows: list[dict[str, object]] = []
        for normal_row in self.boundary_normal_summary_rows():
            bcnr = int(normal_row["boundary_number"])
            name = str(normal_row["name"])
            pressure, source = self._boundary_pressure_value(
                normal_row,
                pressure_by_boundary,
                default_pressure,
            )
            vector = tuple(float(value) for value in normal_row["vector_area"])
            force = tuple(pressure * value for value in vector)
            rows.append({
                "boundary_number": bcnr,
                "name": name,
                "pressure_Pa": pressure,
                "pressure_source": source,
                "surface_area": normal_row["surface_area"],
                "vector_area": vector,
                "unit_normal": normal_row["unit_normal"],
                "force_N": force,
                "force_magnitude_N": _norm(force),
            })
        return tuple(rows)

    def boundary_traction_force_moment_rows(
        self,
        traction_by_boundary: dict[int | str, tuple[float, float, float]],
        default_traction: tuple[float, float, float] | None = (0.0, 0.0, 0.0),
        pivot_m: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> tuple[dict[str, object], ...]:
        """Return per-boundary vector-traction resultant force and moment rows.

        The traction vector is a constant global vector [N/m2] on each named
        boundary.  Each triangle contributes ``area * traction`` and
        ``(centroid - pivot) x dF``.  This is the vector-load companion to
        scalar pressure rows and is useful for readable FEM/BEM examples,
        imported sidesets, and field-postprocessing checks where the traction
        direction is already known.
        """

        pivot = tuple(float(value) for value in pivot_m)
        if len(pivot) != 3:
            raise ValueError("pivot_m must have three components")

        normal_rows = {
            int(row["boundary_number"]): row
            for row in self.boundary_normal_summary_rows()
        }
        accum: dict[int, dict[str, object]] = {}
        for bcnr, normal_row in normal_rows.items():
            traction, source = self._boundary_vector_value(
                normal_row,
                traction_by_boundary,
                default_traction,
                "traction",
            )
            accum[bcnr] = {
                "boundary_number": bcnr,
                "name": str(normal_row["name"]),
                "traction_N_per_m2": traction,
                "traction_source": source,
                "surface_area": float(normal_row["surface_area"]),
                "vector_area": tuple(float(value) for value in normal_row["vector_area"]),
                "unit_normal": normal_row["unit_normal"],
                "_force": (0.0, 0.0, 0.0),
                "_moment": (0.0, 0.0, 0.0),
                "_centroid_weight": (0.0, 0.0, 0.0),
            }

        for tri, area_vector in zip(self.surface_triangles, self.surface_triangle_area_vectors()):
            row = accum[tri.bcnr]
            traction = tuple(float(value) for value in row["traction_N_per_m2"])
            area = _norm(area_vector)
            force = _scale(traction, area)
            a, b, c = (self.points[node - 1] for node in tri.nodes)
            centroid = (
                (a[0] + b[0] + c[0]) / 3.0,
                (a[1] + b[1] + c[1]) / 3.0,
                (a[2] + b[2] + c[2]) / 3.0,
            )
            moment = _cross(_sub(centroid, pivot), force)
            row["_force"] = _add(row["_force"], force)
            row["_moment"] = _add(row["_moment"], moment)
            row["_centroid_weight"] = _add(row["_centroid_weight"], _scale(centroid, area))

        rows: list[dict[str, object]] = []
        for bcnr in sorted(accum):
            row = accum[bcnr]
            force = tuple(float(value) for value in row["_force"])
            moment = tuple(float(value) for value in row["_moment"])
            area = float(row["surface_area"])
            centroid = (
                _scale(row["_centroid_weight"], 1.0 / area)
                if area > 0.0
                else None
            )
            rows.append({
                "boundary_number": row["boundary_number"],
                "name": row["name"],
                "traction_N_per_m2": row["traction_N_per_m2"],
                "traction_source": row["traction_source"],
                "surface_area": area,
                "vector_area": row["vector_area"],
                "unit_normal": row["unit_normal"],
                "centroid_m": centroid,
                "force_N": force,
                "force_magnitude_N": _norm(force),
                "pivot_m": pivot,
                "moment_about_pivot_Nm": moment,
                "moment_magnitude_Nm": _norm(moment),
            })
        return tuple(rows)

    def boundary_pressure_force_moment_rows(
        self,
        pressure_by_boundary: dict[int | str, float],
        default_pressure: float | None = 0.0,
        pivot_m: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> tuple[dict[str, object], ...]:
        """Return per-boundary pressure resultant force and moment rows.

        Pressure is constant on each named boundary and positive along the
        oriented outward normal.  Force is integrated from each triangle's
        area vector; moment is integrated as ``sum((centroid - pivot) x dF)``.
        This makes the rows suitable for Maxwell-pressure, electrostatic
        pressure, radiation-pressure, and acoustic-pressure teaching examples.
        """

        pivot = tuple(float(value) for value in pivot_m)
        if len(pivot) != 3:
            raise ValueError("pivot_m must have three components")

        normal_rows = {
            int(row["boundary_number"]): row
            for row in self.boundary_normal_summary_rows()
        }
        accum: dict[int, dict[str, object]] = {}
        for bcnr, normal_row in normal_rows.items():
            pressure, source = self._boundary_pressure_value(
                normal_row,
                pressure_by_boundary,
                default_pressure,
            )
            accum[bcnr] = {
                "boundary_number": bcnr,
                "name": str(normal_row["name"]),
                "pressure_Pa": pressure,
                "pressure_source": source,
                "surface_area": float(normal_row["surface_area"]),
                "vector_area": tuple(float(value) for value in normal_row["vector_area"]),
                "unit_normal": normal_row["unit_normal"],
                "_force": (0.0, 0.0, 0.0),
                "_moment": (0.0, 0.0, 0.0),
                "_centroid_weight": (0.0, 0.0, 0.0),
            }

        for tri, area_vector in zip(self.surface_triangles, self.surface_triangle_area_vectors()):
            row = accum[tri.bcnr]
            pressure = float(row["pressure_Pa"])
            force = _scale(area_vector, pressure)
            a, b, c = (self.points[node - 1] for node in tri.nodes)
            centroid = (
                (a[0] + b[0] + c[0]) / 3.0,
                (a[1] + b[1] + c[1]) / 3.0,
                (a[2] + b[2] + c[2]) / 3.0,
            )
            area = _norm(area_vector)
            moment = _cross(_sub(centroid, pivot), force)
            row["_force"] = _add(row["_force"], force)
            row["_moment"] = _add(row["_moment"], moment)
            row["_centroid_weight"] = _add(row["_centroid_weight"], _scale(centroid, area))

        rows: list[dict[str, object]] = []
        for bcnr in sorted(accum):
            row = accum[bcnr]
            force = tuple(float(value) for value in row["_force"])
            moment = tuple(float(value) for value in row["_moment"])
            area = float(row["surface_area"])
            centroid = (
                _scale(row["_centroid_weight"], 1.0 / area)
                if area > 0.0
                else None
            )
            rows.append({
                "boundary_number": row["boundary_number"],
                "name": row["name"],
                "pressure_Pa": row["pressure_Pa"],
                "pressure_source": row["pressure_source"],
                "surface_area": area,
                "vector_area": row["vector_area"],
                "unit_normal": row["unit_normal"],
                "centroid_m": centroid,
                "force_N": force,
                "force_magnitude_N": _norm(force),
                "pivot_m": pivot,
                "moment_about_pivot_Nm": moment,
                "moment_magnitude_Nm": _norm(moment),
            })
        return tuple(rows)

    def boundary_pressure_resultant_summary(
        self,
        pressure_by_boundary: dict[int | str, float],
        default_pressure: float | None = 0.0,
        pivot_m: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> dict[str, object]:
        """Return pressure rows plus net force/moment balance metrics.

        This is the compact final reduction after importing named
        Cubit/Coreform sidesets from a tri/tet Netgen ``.vol`` mesh.  A uniform
        pressure over a closed, consistently oriented surface should have
        near-zero total force and moment; a one-sided pressure load should
        reduce to ``pressure * oriented vector area`` and the corresponding
        pivot moment.
        """

        rows = self.boundary_pressure_force_moment_rows(
            pressure_by_boundary,
            default_pressure=default_pressure,
            pivot_m=pivot_m,
        )
        total_force = tuple(
            sum(float(row["force_N"][axis]) for row in rows)
            for axis in range(3)
        )
        total_moment = tuple(
            sum(float(row["moment_about_pivot_Nm"][axis]) for row in rows)
            for axis in range(3)
        )
        total_force_norm = _norm(total_force)
        total_moment_norm = _norm(total_moment)
        absolute_force_sum = sum(float(row["force_magnitude_N"]) for row in rows)
        absolute_moment_sum = sum(float(row["moment_magnitude_Nm"]) for row in rows)
        surface_vector_area = self.surface_vector_area()
        total_area = self.total_surface_area()

        return {
            "boundary_count": len(rows),
            "rows": rows,
            "pivot_m": tuple(float(value) for value in pivot_m),
            "total_force_N": total_force,
            "total_force_magnitude_N": total_force_norm,
            "total_moment_about_pivot_Nm": total_moment,
            "total_moment_magnitude_Nm": total_moment_norm,
            "absolute_force_sum_N": absolute_force_sum,
            "absolute_moment_sum_Nm": absolute_moment_sum,
            "force_balance_ratio": (
                total_force_norm / absolute_force_sum
                if absolute_force_sum > 0.0
                else 0.0
            ),
            "moment_balance_ratio": (
                total_moment_norm / absolute_moment_sum
                if absolute_moment_sum > 0.0
                else 0.0
            ),
            "surface_vector_area": surface_vector_area,
            "surface_vector_area_norm": _norm(surface_vector_area),
            "surface_vector_area_norm_over_area": (
                _norm(surface_vector_area) / total_area if total_area > 0.0 else None
            ),
        }

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

    def boundary_edge_inventory_rows(self) -> tuple[dict[str, object], ...]:
        """Return boundary-local edge inventory rows.

        A named boundary made from triangles has perimeter edges used by one
        triangle in that boundary and split/diagonal edges used by two triangles
        in that same boundary.  Keeping these separate is useful when auditing
        sidesets exported from a CAD/mesh tool: the perimeter is the physical
        boundary curve, while the shared edges are triangulation details.
        """

        groups: dict[int, list[NetgenSurfaceTriangle]] = {}
        area_by_boundary: dict[int, float] = {}
        for tri, area in zip(self.surface_triangles, self.surface_triangle_areas()):
            groups.setdefault(tri.bcnr, []).append(tri)
            area_by_boundary[tri.bcnr] = area_by_boundary.get(tri.bcnr, 0.0) + area

        rows: list[dict[str, object]] = []
        for bcnr in sorted(set(self.boundary_names) | set(groups)):
            edge_counts: dict[tuple[int, int], int] = {}
            for tri in groups.get(bcnr, []):
                n0, n1, n2 = tri.nodes
                for edge in ((n0, n1), (n1, n2), (n2, n0)):
                    key = tuple(sorted(edge))
                    edge_counts[key] = edge_counts.get(key, 0) + 1

            edge_lengths = {
                edge: _norm(_sub(self.points[edge[1] - 1], self.points[edge[0] - 1]))
                for edge in edge_counts
            }
            perimeter_edges = [edge for edge, count in edge_counts.items() if count == 1]
            shared_edges = [edge for edge, count in edge_counts.items() if count == 2]
            overused_edges = [edge for edge, count in edge_counts.items() if count > 2]
            rows.append({
                "boundary_number": bcnr,
                "name": self.boundary_names.get(bcnr, f"boundary_{bcnr}"),
                "surface_triangles": len(groups.get(bcnr, [])),
                "surface_area": area_by_boundary.get(bcnr, 0.0),
                "unique_boundary_edges": len(edge_counts),
                "perimeter_edges": len(perimeter_edges),
                "shared_diagonal_edges": len(shared_edges),
                "overused_edges": len(overused_edges),
                "perimeter_edge_length_sum_m": sum(edge_lengths[edge] for edge in perimeter_edges),
                "shared_diagonal_edge_length_sum_m": sum(edge_lengths[edge] for edge in shared_edges),
                "total_unique_edge_length_sum_m": sum(edge_lengths.values()),
                "min_edge_length_m": min(edge_lengths.values()) if edge_lengths else None,
                "max_edge_length_m": max(edge_lengths.values()) if edge_lengths else None,
                "perimeter_edge_nodes": [list(edge) for edge in sorted(perimeter_edges)],
                "shared_diagonal_edge_nodes": [list(edge) for edge in sorted(shared_edges)],
                "overused_edge_nodes": [list(edge) for edge in sorted(overused_edges)],
            })
        return tuple(rows)

    def boundary_edge_inventory_summary(self) -> dict[str, object]:
        """Return compact boundary-local edge inventory totals."""

        rows = self.boundary_edge_inventory_rows()
        overused = [row for row in rows if int(row["overused_edges"]) > 0]
        return {
            "policy": "netgen_vol_boundary_local_edge_inventory",
            "boundary_count": len(rows),
            "surface_triangles": len(self.surface_triangles),
            "unique_boundary_edges_total": sum(int(row["unique_boundary_edges"]) for row in rows),
            "perimeter_edges_total": sum(int(row["perimeter_edges"]) for row in rows),
            "shared_diagonal_edges_total": sum(int(row["shared_diagonal_edges"]) for row in rows),
            "overused_edges_total": sum(int(row["overused_edges"]) for row in rows),
            "has_overused_boundary_edges": bool(overused),
            "total_perimeter_edge_length_m": sum(float(row["perimeter_edge_length_sum_m"]) for row in rows),
            "total_shared_diagonal_edge_length_m": sum(float(row["shared_diagonal_edge_length_sum_m"]) for row in rows),
            "rows": rows,
        }

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

    def boundary_condition_assignment_summary(
        self,
        condition_by_boundary: dict[int | str, str],
        default_condition: str | None = None,
    ) -> dict[str, object]:
        """Audit boundary-condition labels against named `.vol` boundaries.

        Keys in ``condition_by_boundary`` may be Netgen boundary numbers or
        boundary names.  Values are intentionally plain labels, not solver
        objects.  The summary catches missing boundary assignments and unknown
        keys before a readable FEM/BEM script turns the labels into Dirichlet,
        Neumann, impedance, or coupling operators.
        """

        known_numbers = {int(row["boundary_number"]) for row in self.boundary_summary_rows()}
        known_names = {str(row["name"]) for row in self.boundary_summary_rows()}
        unknown_keys: list[int | str] = []
        for key in condition_by_boundary:
            if isinstance(key, int):
                if key not in known_numbers:
                    unknown_keys.append(key)
            elif isinstance(key, str):
                if key not in known_names:
                    unknown_keys.append(key)
            else:
                raise TypeError("boundary condition keys must be boundary numbers or names")

        incidence_by_bcnr: dict[int, list[dict[str, object]]] = {}
        for row in self.domain_boundary_incidence_rows():
            incidence_by_bcnr.setdefault(int(row["boundary_number"]), []).append(row)

        rows: list[dict[str, object]] = []
        condition_counts: dict[str, int] = {}
        missing = 0
        for boundary in self.boundary_summary_rows():
            bcnr = int(boundary["boundary_number"])
            name = str(boundary["name"])
            if bcnr in condition_by_boundary:
                condition = str(condition_by_boundary[bcnr])
                source = "boundary_number"
            elif name in condition_by_boundary:
                condition = str(condition_by_boundary[name])
                source = "boundary_name"
            elif default_condition is not None:
                condition = str(default_condition)
                source = "default"
            else:
                condition = None
                source = "missing"
                missing += 1
            if condition is not None:
                condition_counts[condition] = condition_counts.get(condition, 0) + 1

            incidence_rows = incidence_by_bcnr.get(bcnr, [])
            rows.append({
                "boundary_number": bcnr,
                "name": name,
                "condition": condition,
                "condition_source": source,
                "surface_triangles": boundary["surface_triangles"],
                "surface_area": boundary["surface_area"],
                "trace_node_count": boundary["trace_node_count"],
                "trace_node_ids": boundary["trace_node_ids"],
                "incidence_kinds": sorted({str(row["kind"]) for row in incidence_rows}),
                "adjacent_material_numbers": sorted({
                    int(domain)
                    for row in incidence_rows
                    for domain in (int(row["domin"]), int(row["domout"]))
                    if domain != 0
                }),
                "adjacent_material_names": sorted({
                    str(row["domin_material"])
                    for row in incidence_rows
                    if row["domin_material"] is not None
                } | {
                    str(row["domout_material"])
                    for row in incidence_rows
                    if row["domout_material"] is not None
                }),
            })

        return {
            "policy": "netgen_vol_boundary_conditions_are_assigned_by_number_or_name",
            "boundary_count": len(rows),
            "assigned_boundary_count": len(rows) - missing,
            "missing_boundary_count": missing,
            "unknown_condition_keys": unknown_keys,
            "unknown_condition_key_count": len(unknown_keys),
            "condition_counts": dict(sorted(condition_counts.items())),
            "ok": missing == 0 and not unknown_keys,
            "rows": tuple(rows),
        }

    def boundary_tet_face_incidence_rows(self) -> tuple[dict[str, object], ...]:
        """Match each boundary triangle to adjacent tetrahedron faces.

        This is a `.vol` export sanity gate for FEM/BEM coupling: every
        exterior surface triangle should match exactly one volume tetrahedron
        face, while a material interface triangle should match two.  The rows
        also compare the matched tetrahedron material numbers with the
        ``domin/domout`` domain ids stored by Netgen/Coreform export.
        """

        tet_faces: dict[tuple[int, int, int], list[dict[str, object]]] = {}
        local_faces = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))
        for tet_index, tet in enumerate(self.tetrahedra, start=1):
            for opposite, local in enumerate(local_faces):
                nodes = tuple(tet.nodes[i] for i in local)
                key = tuple(sorted(nodes))
                tet_faces.setdefault(key, []).append({
                    "tetrahedron": tet_index,
                    "material_number": tet.matnr,
                    "material_name": self.materials.get(tet.matnr, f"material_{tet.matnr}"),
                    "opposite_local_node": opposite,
                    "face_local_nodes": list(local),
                    "face_nodes": list(nodes),
                })

        rows: list[dict[str, object]] = []
        for tri_index, tri in enumerate(self.surface_triangles, start=1):
            key = tuple(sorted(tri.nodes))
            adjacent = tet_faces.get(key, [])
            adjacent_materials = sorted({int(item["material_number"]) for item in adjacent})
            declared_domains = sorted(domain for domain in (tri.domin, tri.domout) if domain != 0)
            if len(adjacent) == 0:
                kind = "orphan"
            elif len(adjacent) == 1:
                kind = "exterior"
            elif len(adjacent) == 2:
                kind = "interface"
            else:
                kind = "overconnected"
            rows.append({
                "surface_triangle": tri_index,
                "boundary_number": tri.bcnr,
                "name": self.boundary_names.get(tri.bcnr, f"boundary_{tri.bcnr}"),
                "nodes": list(tri.nodes),
                "domin": tri.domin,
                "domout": tri.domout,
                "declared_domain_numbers": declared_domains,
                "declared_domain_names": [
                    self.materials.get(domain, f"material_{domain}")
                    for domain in declared_domains
                ],
                "adjacent_tetrahedron_faces": adjacent,
                "adjacent_tetrahedron_count": len(adjacent),
                "adjacent_material_numbers": adjacent_materials,
                "adjacent_material_names": [
                    self.materials.get(matnr, f"material_{matnr}")
                    for matnr in adjacent_materials
                ],
                "domain_material_match": declared_domains == adjacent_materials,
                "kind": kind,
            })
        return tuple(rows)

    def boundary_tet_face_incidence_summary(self) -> dict[str, object]:
        """Return compact counts for boundary-triangle to tetra-face incidence."""

        rows = self.boundary_tet_face_incidence_rows()
        kind_counts: dict[str, int] = {}
        for row in rows:
            kind = str(row["kind"])
            kind_counts[kind] = kind_counts.get(kind, 0) + 1
        mismatches = [row for row in rows if not bool(row["domain_material_match"])]
        return {
            "surface_triangles": len(rows),
            "tetrahedra": len(self.tetrahedra),
            "exterior_surface_triangles": kind_counts.get("exterior", 0),
            "interface_surface_triangles": kind_counts.get("interface", 0),
            "orphan_surface_triangles": kind_counts.get("orphan", 0),
            "overconnected_surface_triangles": kind_counts.get("overconnected", 0),
            "domain_material_mismatch_count": len(mismatches),
            "max_adjacent_tetrahedra": max(
                (int(row["adjacent_tetrahedron_count"]) for row in rows),
                default=0,
            ),
            "is_volume_boundary_consistent": (
                kind_counts.get("orphan", 0) == 0
                and kind_counts.get("overconnected", 0) == 0
                and not mismatches
            ),
            "rows": rows,
            "policy": "surface_triangles_match_volume_tet_faces_for_tri_tet_vol",
        }

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

    def mesh_health_summary(
        self,
        *,
        min_surface_triangle_quality: float = 1.0e-8,
        min_tetrahedron_quality: float = 1.0e-8,
        closure_relative_tolerance: float = 1.0e-9,
        worst_limit: int = 5,
    ) -> dict[str, object]:
        """Return a readable go/no-go report for first-order FEM/BEM use.

        This bundles the checks that matter before a Cubit/Coreform Netgen
        ``.vol`` export is used as a shared FEM/BEM trace: volume tets exist,
        boundary triangles exist, the surface closes, boundary triangles match
        tetrahedron faces, and the worst shape-quality rows are easy to find.
        """

        if min_surface_triangle_quality < 0.0:
            raise ValueError("min_surface_triangle_quality must be non-negative")
        if min_tetrahedron_quality < 0.0:
            raise ValueError("min_tetrahedron_quality must be non-negative")
        if closure_relative_tolerance < 0.0:
            raise ValueError("closure_relative_tolerance must be non-negative")
        if worst_limit < 0:
            raise ValueError("worst_limit must be non-negative")

        inventory = self.summary()
        surface_quality = self.surface_triangle_quality_summary()
        tet_quality = self.tetrahedron_quality_summary()
        closure = self.surface_closure_summary()
        manifold = self.surface_edge_manifold_summary()
        incidence_full = self.boundary_tet_face_incidence_summary()
        incidence = {
            key: value
            for key, value in incidence_full.items()
            if key != "rows"
        }

        surface_min_quality = surface_quality["min_radius_ratio_quality"]
        tet_min_quality = tet_quality["min_radius_ratio_quality"]
        closure_rel_error = closure["surface_abs_volume_rel_error"]
        vector_area_rel = closure["surface_vector_area_norm_over_area"]

        checks: dict[str, bool] = {}
        issues: list[str] = []

        def add_check(name: str, ok: bool, issue: str) -> None:
            checks[name] = ok
            if not ok:
                issues.append(issue)

        add_check("has_points", inventory["points"] > 0, "mesh has no points")
        add_check(
            "has_surface_triangles",
            inventory["surface_triangles"] > 0,
            "mesh has no boundary triangles",
        )
        add_check(
            "has_tetrahedra",
            inventory["tetrahedra"] > 0,
            "mesh has no volume tetrahedra",
        )
        add_check(
            "surface_is_closed_manifold",
            bool(manifold["is_closed_manifold"]),
            "boundary triangles are not a closed edge manifold",
        )
        add_check(
            "surface_vector_area_closes",
            vector_area_rel is not None and float(vector_area_rel) <= closure_relative_tolerance,
            "oriented boundary vector area does not close",
        )
        add_check(
            "surface_volume_matches_tetrahedra",
            closure_rel_error is not None and float(closure_rel_error) <= closure_relative_tolerance,
            "boundary-triangle enclosed volume does not match tetrahedron volume",
        )
        add_check(
            "boundary_faces_match_tetrahedra",
            bool(incidence["is_volume_boundary_consistent"]),
            "at least one boundary triangle is orphaned, overconnected, or material-mismatched",
        )
        add_check(
            "surface_triangle_quality_above_threshold",
            surface_min_quality is not None and float(surface_min_quality) >= min_surface_triangle_quality,
            "at least one boundary triangle is below the surface quality threshold",
        )
        add_check(
            "tetrahedron_quality_above_threshold",
            tet_min_quality is not None and float(tet_min_quality) >= min_tetrahedron_quality,
            "at least one tetrahedron is below the volume quality threshold",
        )

        ok = all(checks.values())
        return {
            "policy": "netgen_vol_tri_tet_first_order_fem_bem_mesh_health",
            "status": "ok" if ok else "needs_attention",
            "ok_for_first_order_fem_bem": ok,
            "inventory": inventory,
            "thresholds": {
                "min_surface_triangle_quality": min_surface_triangle_quality,
                "min_tetrahedron_quality": min_tetrahedron_quality,
                "closure_relative_tolerance": closure_relative_tolerance,
            },
            "checks": checks,
            "issues": issues,
            "surface_triangle_quality": surface_quality,
            "tetrahedron_quality": tet_quality,
            "closure": closure,
            "manifold": manifold,
            "boundary_tet_face_incidence": incidence,
            "worst_surface_triangles": self.worst_surface_triangles_by_quality(worst_limit),
            "worst_tetrahedra": self.worst_tetrahedra_by_quality(worst_limit),
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

    def p1_fem_bem_trace_matrix_summary(self) -> dict[str, object]:
        """Return a one-based sparse P1 H1-to-scalar-BEM trace matrix.

        For first-order nodal FEM and first-order scalar BEM on the same
        tri/tet ``.vol`` mesh, the trace operator is just a boolean gather:
        each compact boundary node row points to exactly one volume H1 node
        column.  The returned one-based COO arrays can be used directly as
        ``sparse(rows, cols, values, nTraceNodes, nVolumeNodes)`` in numerical
        environments that use one-based sparse indexing.
        """

        trace_nodes = list(self.trace_node_ids())
        volume_nodes = list(range(1, len(self.points) + 1))
        trace_node_set = set(trace_nodes)
        interior_nodes = [node for node in volume_nodes if node not in trace_node_set]
        rows = list(range(1, len(trace_nodes) + 1))
        cols = list(trace_nodes)
        values = [1.0 for _node in trace_nodes]
        return {
            "policy": "p1_h1_to_scalar_bem_trace_is_boolean_gather",
            "sparse_coo_call": "T = sparse(rows, cols, values, nTraceNodes, nVolumeNodes)",
            "n_volume_nodes": len(volume_nodes),
            "n_trace_nodes": len(trace_nodes),
            "matrix_shape": [len(trace_nodes), len(volume_nodes)],
            "nnz": len(trace_nodes),
            "rows": rows,
            "cols": cols,
            "values": values,
            "trace_node_ids": trace_nodes,
            "interior_node_ids": interior_nodes,
            "is_boolean_gather": True,
            "row_nnz_min": 1 if trace_nodes else 0,
            "row_nnz_max": 1 if trace_nodes else 0,
            "boundary_column_nnz_max": 1 if trace_nodes else 0,
            "interior_column_nnz_max": 0,
            "surface_triangles_local": [
                [trace_nodes.index(node) + 1 for node in tri.nodes]
                for tri in self.surface_triangles
            ],
            "surface_triangles_global": [list(tri.nodes) for tri in self.surface_triangles],
        }

    def first_order_fem_bem_topology(self) -> dict[str, object]:
        """Return first-order H1/HCurl/P1/RWG topology with one-based ids.

        This mirrors a small readable prototype API: H1 uses volume nodes,
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

    def boundary_oriented_edge_summary(self) -> dict[str, object]:
        """Return one row per oriented boundary-triangle edge.

        This is a teaching-friendly expansion of the RWG part of
        :meth:`first_order_fem_bem_topology`: each surface triangle contributes
        three local oriented edges, each linked to the compact boundary edge id,
        the global volume-node edge, the sign relative to the sorted global edge,
        and the matching HCurl edge id when the boundary edge is a closed RWG
        degree of freedom.
        """

        topology = self.first_order_fem_bem_topology()
        trace_nodes = tuple(topology["scalar_bem"]["global_node_ids"])
        rwg = topology["rwg"]
        surface_edges_global = tuple(tuple(edge) for edge in rwg["edges_global"])
        hcurl_by_rwg_edge = {
            int(rwg_edge_id): int(hcurl_edge_id)
            for rwg_edge_id, hcurl_edge_id in zip(rwg["dof_edge_ids"], rwg["hcurl_edge_ids"])
        }
        dof_edge_ids = set(hcurl_by_rwg_edge)

        rows: list[dict[str, object]] = []
        sign_counts = {-1: 0, 1: 0}
        for tri_index, tri in enumerate(self.surface_triangles, start=1):
            for local_edge_index, (edge_id, sign) in enumerate(
                zip(rwg["tri_edges"][tri_index - 1], rwg["tri_edge_signs"][tri_index - 1]),
                start=1,
            ):
                edge_id = int(edge_id)
                sign = int(sign)
                sign_counts[sign] = sign_counts.get(sign, 0) + 1
                edge_global = surface_edges_global[edge_id - 1]
                oriented_edge = edge_global if sign > 0 else tuple(reversed(edge_global))
                adjacent = [int(value) for value in rwg["edge_triangles"][edge_id - 1] if int(value) != 0]
                opposites_local = [
                    int(value)
                    for value in rwg["opposite_vertices_local"][edge_id - 1]
                    if int(value) != 0
                ]
                opposites_global = [int(trace_nodes[value - 1]) for value in opposites_local]
                p0, p1 = (self.points[node - 1] for node in edge_global)
                rows.append({
                    "surface_triangle": tri_index,
                    "local_edge": local_edge_index,
                    "boundary_number": tri.bcnr,
                    "name": self.boundary_names.get(tri.bcnr, f"boundary_{tri.bcnr}"),
                    "triangle_nodes_global": list(tri.nodes),
                    "rwg_edge_id": edge_id,
                    "is_rwg_dof": edge_id in dof_edge_ids,
                    "hcurl_edge_id": hcurl_by_rwg_edge.get(edge_id),
                    "edge_nodes_global": list(edge_global),
                    "oriented_edge_nodes_global": list(oriented_edge),
                    "orientation_sign": sign,
                    "edge_length_m": _norm(_sub(p1, p0)),
                    "adjacent_surface_triangles": adjacent,
                    "adjacent_surface_triangle_count": len(adjacent),
                    "opposite_vertices_global": opposites_global,
                })

        manifold = self.surface_edge_manifold_summary()
        return {
            "policy": "boundary_triangle_oriented_edges_for_first_order_rwg_trace",
            "surface_triangles": len(self.surface_triangles),
            "oriented_edge_rows": len(rows),
            "surface_edges": len(surface_edges_global),
            "rwg_dof_edges": len(dof_edge_ids),
            "hcurl_trace_edges": len(hcurl_by_rwg_edge),
            "open_edges": manifold["open_edges"],
            "closed_edges": manifold["closed_edges"],
            "is_closed_manifold": manifold["is_closed_manifold"],
            "orientation_sign_counts": {str(key): value for key, value in sorted(sign_counts.items())},
            "max_edge_length_m": max((float(row["edge_length_m"]) for row in rows), default=0.0),
            "min_edge_length_m": min((float(row["edge_length_m"]) for row in rows), default=0.0),
            "rows": rows,
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
