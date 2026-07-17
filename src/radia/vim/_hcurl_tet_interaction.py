"""Analytic affine-tetrahedron interaction for reduced HCurl currents.

The sampled Laplace backend remains useful for fast exploration.  This module
is the acceptance/production reference for volume-only HCurl-VIM: reconstruct
``curl(T)`` as element-local reference polynomials, evaluate every source
tetrahedron with Radia's analytic Newton-potential moments, and quadrature only
the smooth outer integral.  No diagonal epsilon is used.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

import radia._radia_pybind as _rp

from ._eddy_hybrid import (
    MU0,
    NgsolveHCurlCellFamilies,
    SampledCurrentBasis,
)
from ._vim import _f64_buffer, _i32_buffer, _monos_vol, _outer_tet, _tet_ref
from ._vim import (
    _IR_TET_NODES,
    _g01,
    _outer_tri,
    _trafo_lattice_nodes,
    _tri_ref,
)


def _labels(value) -> set[str] | None:
    if value is None:
        return None
    if isinstance(value, str):
        return {value}
    return {str(item) for item in value}


def _mode_gridfunctions(fes, vectors):
    import ngsolve as ng

    values = np.asarray(vectors)
    if values.ndim == 1:
        values = values[:, np.newaxis]
    if values.ndim != 2 or values.shape[0] != fes.ndof:
        raise ValueError(f"vectors must have shape ({fes.ndof}, n_modes)")
    gridfunctions = []
    currents = []
    for column in range(values.shape[1]):
        gf = ng.GridFunction(fes)
        target = gf.vec.FV().NumPy()
        if np.iscomplexobj(values) and not np.iscomplexobj(target):
            imag = float(np.max(np.abs(np.imag(values[:, column]))))
            if imag > 1.0e-13:
                raise ValueError("the affine tet interaction currently requires real basis vectors")
        target[:] = np.real(values[:, column])
        gridfunctions.append(gf)
        currents.append(ng.curl(gf))
    return gridfunctions, currents


_REFERENCE_VERTICES = {
    "tet": np.asarray(
        ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, 0.0))
    ),
    "hex": np.asarray(
        (
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 1.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 1.0),
            (1.0, 1.0, 1.0),
            (0.0, 1.0, 1.0),
        )
    ),
    "wedge": np.asarray(
        (
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 1.0),
            (0.0, 1.0, 1.0),
            (0.0, 0.0, 1.0),
        )
    ),
    "pyramid": np.asarray(
        (
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 1.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        )
    ),
}

_CELL_SUBTETS = {
    "tet": ((0, 1, 2, 3),),
    "hex": (
        (0, 1, 2, 6),
        (0, 2, 3, 6),
        (0, 3, 7, 6),
        (0, 7, 4, 6),
        (0, 4, 5, 6),
        (0, 5, 1, 6),
    ),
    "wedge": ((0, 1, 2, 5), (0, 1, 5, 4), (0, 4, 5, 3)),
    "pyramid": ((0, 1, 2, 4), (0, 2, 3, 4)),
}

_ELEMENT_TYPE_TO_FAMILY = {
    "ET.TET": "tet",
    "ET.HEX": "hex",
    "ET.PRISM": "wedge",
    "ET.PYRAMID": "pyramid",
}


def _positive_reference_tet(vertices: np.ndarray) -> np.ndarray:
    tet = np.array(vertices, dtype=float, copy=True)
    determinant = float(
        np.linalg.det(
            np.column_stack(
                (tet[1] - tet[0], tet[2] - tet[0], tet[3] - tet[0])
            )
        )
    )
    if abs(determinant) <= 1.0e-18:
        raise ValueError("degenerate reference sub-tetrahedron")
    if determinant < 0.0:
        tet[[1, 2]] = tet[[2, 1]]
    return tet


def _red_refine_reference_tet(vertices: np.ndarray) -> tuple[np.ndarray, ...]:
    """Split one tetrahedron into eight midpoint children."""

    v0, v1, v2, v3 = np.asarray(vertices, dtype=float)
    m01 = 0.5 * (v0 + v1)
    m02 = 0.5 * (v0 + v2)
    m03 = 0.5 * (v0 + v3)
    m12 = 0.5 * (v1 + v2)
    m13 = 0.5 * (v1 + v3)
    m23 = 0.5 * (v2 + v3)
    children = (
        (v0, m01, m02, m03),
        (m01, v1, m12, m13),
        (m02, m12, v2, m23),
        (m03, m13, m23, v3),
        (m01, m02, m03, m23),
        (m01, m02, m12, m23),
        (m01, m03, m13, m23),
        (m01, m12, m13, m23),
    )
    return tuple(_positive_reference_tet(child) for child in children)


def _reference_subtets(
    family: str,
    subdivision_level: int,
    refinement_strategy: str = "pyramid-apex",
) -> tuple[np.ndarray, ...]:
    reference_vertices = _REFERENCE_VERTICES[family]
    tets = tuple(
        _positive_reference_tet(reference_vertices[np.asarray(indices)])
        for indices in _CELL_SUBTETS[family]
    )
    if refinement_strategy not in {"pyramid-apex", "uniform"}:
        raise ValueError(
            "refinement_strategy must be 'pyramid-apex' or 'uniform'"
        )
    apex = np.asarray((0.0, 0.0, 1.0))
    for _ in range(int(subdivision_level)):
        refined: list[np.ndarray] = []
        for parent in tets:
            refine = refinement_strategy == "uniform" or (
                family == "pyramid"
                and np.any(np.linalg.norm(parent - apex, axis=1) <= 1.0e-14)
            )
            if refine:
                refined.extend(_red_refine_reference_tet(parent))
            else:
                refined.append(parent)
        tets = tuple(refined)
    return tets


def _tet_bernstein_exponents(degree: int) -> list[tuple[int, int, int, int]]:
    return [
        (i0, i1, i2, degree - i0 - i1 - i2)
        for i0 in range(degree + 1)
        for i1 in range(degree + 1 - i0)
        for i2 in range(degree + 1 - i0 - i1)
    ]


def _tet_bernstein_vandermonde(
    degree: int,
    points: np.ndarray,
) -> np.ndarray:
    bernstein_exponents = _tet_bernstein_exponents(degree)
    barycentric = np.column_stack(
        (1.0 - np.sum(points, axis=1), points[:, 0], points[:, 1], points[:, 2])
    )
    vandermonde = np.empty((len(points), len(bernstein_exponents)), dtype=float)
    factorial_degree = math.factorial(degree)
    for column, exponent in enumerate(bernstein_exponents):
        multinomial = factorial_degree
        for value in exponent:
            multinomial /= math.factorial(value)
        result = np.full(len(points), multinomial, dtype=float)
        for coordinate, power in enumerate(exponent):
            if power:
                result *= barycentric[:, coordinate] ** power
        vandermonde[:, column] = result
    return vandermonde


def _tet_bernstein_projection(
    degree: int,
    validation_points: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """Return stable lattice fit data and an exact monomial conversion."""

    degree = int(degree)
    monomials = _monos_vol(degree)
    monomial_index = {exponent: index for index, exponent in enumerate(monomials)}
    bernstein_exponents = _tet_bernstein_exponents(degree)
    denominator = degree + 2
    fit_points = np.asarray(
        [
            (i / denominator, j / denominator, k / denominator)
            for i in range(denominator + 1)
            for j in range(denominator + 1 - i)
            for k in range(denominator + 1 - i - j)
        ],
        dtype=float,
    )
    fit_vandermonde = _tet_bernstein_vandermonde(degree, fit_points)
    validation_vandermonde = _tet_bernstein_vandermonde(
        degree,
        validation_points,
    )
    expansion = np.zeros(
        (len(bernstein_exponents), len(monomials)),
        dtype=float,
    )
    factorial_degree = math.factorial(degree)
    for column, exponent in enumerate(bernstein_exponents):
        multinomial = factorial_degree
        for value in exponent:
            multinomial /= math.factorial(value)
        i0, i1, i2, i3 = exponent
        for r in range(i0 + 1):
            for s in range(i0 + 1 - r):
                for t in range(i0 + 1 - r - s):
                    remainder = i0 - r - s - t
                    coefficient = multinomial * math.factorial(i0)
                    coefficient /= (
                        math.factorial(r)
                        * math.factorial(s)
                        * math.factorial(t)
                        * math.factorial(remainder)
                    )
                    coefficient *= (-1.0) ** (r + s + t)
                    target = (i1 + r, i2 + s, i3 + t)
                    expansion[column, monomial_index[target]] += coefficient

    left, singular_values, right_h = np.linalg.svd(
        fit_vandermonde,
        full_matrices=False,
    )
    if singular_values.size == 0 or singular_values[-1] <= 0.0:
        condition = float("inf")
        fit_operator = np.empty((fit_vandermonde.shape[1], len(fit_points)))
    else:
        condition = float(singular_values[0] / singular_values[-1])
        cutoff = 1.0e-13 * singular_values[0]
        inverse = np.where(singular_values > cutoff, 1.0 / singular_values, 0.0)
        fit_operator = (right_h.T * inverse) @ left.T
    return (
        fit_points,
        fit_operator,
        validation_vandermonde,
        expansion,
        condition,
    )


def _project_reference_currents(
    mesh,
    fes,
    vectors,
    *,
    degree: int,
    projection_quad: int,
    materials,
):
    import ngsolve as ng

    selected = _labels(materials)
    elements = [
        element
        for element in mesh.Elements(ng.VOL)
        if selected is None or str(element.mat) in selected
    ]
    if not elements:
        raise ValueError("the selected conductor region contains no volume elements")
    if any(len(element.vertices) != 4 for element in elements):
        raise NotImplementedError("HCurl analytic volume interaction currently requires tetrahedra")

    gridfunctions, currents = _mode_gridfunctions(fes, vectors)
    n_modes = len(currents)
    monomials = _monos_vol(degree)
    ref_points, ref_weights = _tet_ref(projection_quad)
    integration_rule = ng.IntegrationRule(
        [tuple(point) for point in ref_points],
        list(ref_weights),
    )
    coordinate = ng.CF((ng.x, ng.y, ng.z))
    coefficients = np.zeros(
        (n_modes, len(elements), len(monomials), 3),
        dtype=float,
    )
    cell_verts = np.empty((len(elements), 4, 3), dtype=float)
    residual_sq = 0.0
    field_sq = 0.0
    max_vandermonde_condition = 0.0

    for cell, element in enumerate(elements):
        vertices = np.asarray([mesh[v].point for v in element.vertices], dtype=float)
        cell_verts[cell] = vertices
        jacobian = np.column_stack(
            (vertices[1] - vertices[0], vertices[2] - vertices[0], vertices[3] - vertices[0])
        )
        det = float(np.linalg.det(jacobian))
        if not np.isfinite(det) or abs(det) <= 1.0e-18:
            raise ValueError(f"degenerate tetrahedron at selected cell {cell}")

        mapped_rule = mesh.GetTrafo(element)(integration_rule)
        physical_points = np.asarray(coordinate(mapped_rule), dtype=float).reshape(-1, 3)
        local_points = np.linalg.solve(
            jacobian,
            (physical_points - vertices[0]).T,
        ).T
        vandermonde = np.asarray(
            [
                [
                    point[0] ** i * point[1] ** j * point[2] ** k
                    for i, j, k in monomials
                ]
                for point in local_points
            ],
            dtype=float,
        )
        condition = float(np.linalg.cond(vandermonde))
        max_vandermonde_condition = max(max_vandermonde_condition, condition)
        if not np.isfinite(condition):
            raise ValueError(f"singular current projection at selected cell {cell}")

        sampled = np.stack(
            [np.asarray(current(mapped_rule), dtype=float).reshape(-1, 3) for current in currents],
            axis=1,
        )
        flat_sampled = sampled.reshape(sampled.shape[0], -1)
        local_coefficients, *_ = np.linalg.lstsq(
            vandermonde,
            flat_sampled,
            rcond=1.0e-13,
        )
        reconstructed = vandermonde @ local_coefficients
        residual_sq += float(np.sum((reconstructed - flat_sampled) ** 2))
        field_sq += float(np.sum(flat_sampled**2))
        coefficients[:, cell, :, :] = local_coefficients.reshape(
            len(monomials), n_modes, 3
        ).transpose(1, 0, 2)

    relative_residual = float(np.sqrt(residual_sq / max(field_sq, np.finfo(float).tiny)))
    del gridfunctions
    return {
        "cell_verts": cell_verts,
        "coefficients": coefficients,
        "exponents": np.asarray(monomials, dtype=np.int32),
        "relative_residual": relative_residual,
        "max_vandermonde_condition": max_vandermonde_condition,
        "cell_count": len(elements),
        "mode_count": n_modes,
    }


def _project_cell_currents_to_subtets(
    mesh,
    fes,
    vectors,
    *,
    degree: int,
    projection_quad: int,
    subdivision_level: int,
    refinement_strategy: str = "pyramid-apex",
    max_subtets: int = 512,
    max_dense_moment_pairs: int = 20_000_000,
    materials,
):
    """Project reduced HCurl currents onto a tetrahedralization of 3-D cells."""

    import ngsolve as ng

    selected = _labels(materials)
    elements = [
        element
        for element in mesh.Elements(ng.VOL)
        if selected is None or str(element.mat) in selected
    ]
    if not elements:
        raise ValueError("the selected conductor region contains no volume elements")
    families = []
    for element in elements:
        family = _ELEMENT_TYPE_TO_FAMILY.get(str(element.type))
        if family is None:
            raise NotImplementedError(
                f"HCurl cell interaction does not support {element.type}"
            )
        families.append(family)

    gridfunctions, currents = _mode_gridfunctions(fes, vectors)
    n_modes = len(currents)
    monomials = _monos_vol(degree)
    local_points, local_weights = _tet_ref(projection_quad)
    (
        fit_points,
        fit_operator,
        validation_vandermonde,
        bernstein_to_monomial,
        condition,
    ) = _tet_bernstein_projection(
        degree,
        local_points,
    )
    if not np.isfinite(condition):
        raise ValueError("singular sub-tetrahedron current projection")

    subdivision_level = int(subdivision_level)
    if subdivision_level < 0:
        raise ValueError("subdivision_level must be non-negative")
    family_reference_tets = {
        family: _reference_subtets(
            family,
            subdivision_level,
            refinement_strategy,
        )
        for family in set(families)
    }
    subtet_count = sum(len(family_reference_tets[family]) for family in families)
    max_subtets = int(max_subtets)
    if max_subtets <= 0:
        raise ValueError("max_subtets must be positive")
    if subtet_count > max_subtets:
        raise ValueError(
            "adaptive HCurl analytic interaction would require "
            f"{subtet_count} leaf tetrahedra, above max_subtets={max_subtets}"
        )
    max_dense_moment_pairs = int(max_dense_moment_pairs)
    if max_dense_moment_pairs <= 0:
        raise ValueError("max_dense_moment_pairs must be positive")
    dense_moment_pairs = subtet_count * subtet_count * len(monomials)
    if dense_moment_pairs > max_dense_moment_pairs:
        raise ValueError(
            "adaptive HCurl analytic interaction would require "
            f"{dense_moment_pairs} leaf-pair polynomial moments, above "
            f"max_dense_moment_pairs={max_dense_moment_pairs}"
        )
    coefficients = np.zeros(
        (n_modes, subtet_count, len(monomials), 3),
        dtype=float,
    )
    cell_verts = np.empty((subtet_count, 4, 3), dtype=float)
    coordinate = ng.CF((ng.x, ng.y, ng.z))
    residual_sq = 0.0
    field_sq = 0.0
    geometry_error = 0.0
    geometry_scale = 0.0
    family_subtet_counts: dict[str, int] = {}
    subtet = 0

    for element, family in zip(elements, families):
        physical_vertices = np.asarray(
            [mesh[vertex].point for vertex in element.vertices],
            dtype=float,
        )
        reference_vertices = _REFERENCE_VERTICES[family]
        if physical_vertices.shape != reference_vertices.shape:
            raise RuntimeError(
                f"{family} vertex contract changed: {physical_vertices.shape}"
            )
        element_scale = float(
            np.max(
                np.linalg.norm(
                    physical_vertices[:, None, :] - physical_vertices[None, :, :],
                    axis=2,
                )
            )
        )
        geometry_scale = max(geometry_scale, element_scale)
        trafo = mesh.GetTrafo(element)
        for ref_tet in family_reference_tets[family]:
            vertex_rule = ng.IntegrationRule(
                [tuple(float(value) for value in point) for point in ref_tet],
                [1.0] * 4,
            )
            physical_tet = np.array(
                coordinate(trafo(vertex_rule)),
                dtype=float,
                copy=True,
            ).reshape(4, 3)
            cell_verts[subtet] = physical_tet
            jacobian = np.column_stack(
                (
                    physical_tet[1] - physical_tet[0],
                    physical_tet[2] - physical_tet[0],
                    physical_tet[3] - physical_tet[0],
                )
            )
            abs_determinant = abs(float(np.linalg.det(jacobian)))
            if abs_determinant <= 1.0e-18:
                raise ValueError(f"degenerate {family} sub-tetrahedron {subtet}")

            fit_element_points = (
                ref_tet[0]
                + fit_points[:, 0, None] * (ref_tet[1] - ref_tet[0])
                + fit_points[:, 1, None] * (ref_tet[2] - ref_tet[0])
                + fit_points[:, 2, None] * (ref_tet[3] - ref_tet[0])
            )
            fit_rule = ng.IntegrationRule(
                [
                    tuple(float(value) for value in point)
                    for point in fit_element_points
                ],
                [1.0] * len(fit_element_points),
            )
            mapped_fit_rule = trafo(fit_rule)
            sampled_fit = np.stack(
                [
                    np.asarray(current(mapped_fit_rule), dtype=float).reshape(-1, 3)
                    for current in currents
                ],
                axis=1,
            ).reshape(len(fit_points), -1)
            bernstein_coefficients = fit_operator @ sampled_fit

            element_points = (
                ref_tet[0]
                + local_points[:, 0, None] * (ref_tet[1] - ref_tet[0])
                + local_points[:, 1, None] * (ref_tet[2] - ref_tet[0])
                + local_points[:, 2, None] * (ref_tet[3] - ref_tet[0])
            )
            rule = ng.IntegrationRule(
                [tuple(float(value) for value in point) for point in element_points],
                [float(weight) for weight in local_weights],
            )
            mapped_rule = trafo(rule)
            mapped_points = np.asarray(coordinate(mapped_rule), dtype=float).reshape(-1, 3)
            affine_points = (
                physical_tet[0]
                + local_points[:, 0, None] * (physical_tet[1] - physical_tet[0])
                + local_points[:, 1, None] * (physical_tet[2] - physical_tet[0])
                + local_points[:, 2, None] * (physical_tet[3] - physical_tet[0])
            )
            geometry_error = max(
                geometry_error,
                float(np.max(np.linalg.norm(mapped_points - affine_points, axis=1))),
            )

            sampled = np.stack(
                [
                    np.asarray(current(mapped_rule), dtype=float).reshape(-1, 3)
                    for current in currents
                ],
                axis=1,
            )
            flat_sampled = sampled.reshape(sampled.shape[0], -1)
            reconstructed = validation_vandermonde @ bernstein_coefficients
            residual_sq += abs_determinant * float(
                np.sum(
                    local_weights[:, None]
                    * (reconstructed - flat_sampled) ** 2
                )
            )
            field_sq += abs_determinant * float(
                np.sum(local_weights[:, None] * flat_sampled**2)
            )
            local_coefficients = (
                bernstein_to_monomial.T @ bernstein_coefficients
            )
            coefficients[:, subtet, :, :] = local_coefficients.reshape(
                len(monomials), n_modes, 3
            ).transpose(1, 0, 2)
            family_subtet_counts[family] = family_subtet_counts.get(family, 0) + 1
            subtet += 1

    relative_residual = float(
        np.sqrt(residual_sq / max(field_sq, np.finfo(float).tiny))
    )
    relative_geometry_error = float(
        geometry_error / max(geometry_scale, np.finfo(float).tiny)
    )
    del gridfunctions
    return {
        "cell_verts": cell_verts,
        "coefficients": coefficients,
        "exponents": np.asarray(monomials, dtype=np.int32),
        "relative_residual": relative_residual,
        "relative_geometry_error": relative_geometry_error,
        "max_vandermonde_condition": condition,
        "cell_count": len(elements),
        "subtet_count": subtet_count,
        "mode_count": n_modes,
        "family_subtet_counts": family_subtet_counts,
        "subdivision_level": subdivision_level,
        "refinement_strategy": refinement_strategy,
        "dense_moment_pairs": dense_moment_pairs,
    }


def _p2_tet_map(nodes: np.ndarray, points: np.ndarray) -> np.ndarray:
    """Evaluate the C++ curved-tet P2 map convention in NumPy."""

    xi = np.asarray(points, dtype=float)
    barycentric = np.column_stack(
        (1.0 - np.sum(xi, axis=1), xi[:, 0], xi[:, 1], xi[:, 2])
    )
    shapes = np.column_stack(
        (
            barycentric[:, 0] * (2.0 * barycentric[:, 0] - 1.0),
            barycentric[:, 1] * (2.0 * barycentric[:, 1] - 1.0),
            barycentric[:, 2] * (2.0 * barycentric[:, 2] - 1.0),
            barycentric[:, 3] * (2.0 * barycentric[:, 3] - 1.0),
            4.0 * barycentric[:, 0] * barycentric[:, 1],
            4.0 * barycentric[:, 1] * barycentric[:, 2],
            4.0 * barycentric[:, 2] * barycentric[:, 0],
            4.0 * barycentric[:, 0] * barycentric[:, 3],
            4.0 * barycentric[:, 1] * barycentric[:, 3],
            4.0 * barycentric[:, 2] * barycentric[:, 3],
        )
    )
    return shapes @ np.asarray(nodes, dtype=float)


def _project_curved_tet_reference_currents(
    mesh,
    fes,
    vectors,
    *,
    degree: int,
    projection_quad: int,
    materials,
):
    """Project K=J|det(dX/dxi)|, the curl-Piola reference density."""

    import ngsolve as ng

    selected = _labels(materials)
    elements = [
        element
        for element in mesh.Elements(ng.VOL)
        if selected is None or str(element.mat) in selected
    ]
    if not elements:
        raise ValueError("the selected conductor region contains no volume elements")
    if any(str(element.type) != "ET.TET" for element in elements):
        raise NotImplementedError(
            "the exact curved reference-density path currently requires P2 tetrahedra"
        )

    gridfunctions, currents = _mode_gridfunctions(fes, vectors)
    n_modes = len(currents)
    monomials = _monos_vol(degree)
    validation_points, validation_weights = _tet_ref(projection_quad)
    (
        fit_points,
        fit_operator,
        validation_vandermonde,
        bernstein_to_monomial,
        condition,
    ) = _tet_bernstein_projection(degree, validation_points)
    if not np.isfinite(condition):
        raise ValueError("singular curved reference-current projection")

    fit_rule = ng.IntegrationRule(
        [tuple(float(value) for value in point) for point in fit_points],
        [1.0] * len(fit_points),
    )
    validation_rule = ng.IntegrationRule(
        [tuple(float(value) for value in point) for point in validation_points],
        [float(weight) for weight in validation_weights],
    )
    coordinate = ng.CF((ng.x, ng.y, ng.z))
    jacobian_determinant = ng.Det(ng.specialcf.JacobianMatrix(3))
    coefficients = np.zeros(
        (n_modes, len(elements), len(monomials), 3),
        dtype=float,
    )
    cell_nodes = np.empty((len(elements), 10, 3), dtype=float)
    cell_vertices: list[int] = []
    residual_sq = 0.0
    field_sq = 0.0
    geometry_error = 0.0
    geometry_scale = 0.0

    for cell, element in enumerate(elements):
        trafo = mesh.GetTrafo(element)
        nodes = np.asarray(
            _trafo_lattice_nodes(mesh, element, _IR_TET_NODES),
            dtype=float,
        ).reshape(10, 3)
        cell_nodes[cell] = nodes
        cell_vertices.extend(int(vertex.nr) for vertex in element.vertices)
        geometry_scale = max(
            geometry_scale,
            float(
                np.max(
                    np.linalg.norm(
                        nodes[:, None, :] - nodes[None, :, :],
                        axis=2,
                    )
                )
            ),
        )

        mapped_fit = trafo(fit_rule)
        fit_measure = np.abs(
            np.asarray(jacobian_determinant(mapped_fit), dtype=float).reshape(-1)
        )
        sampled_fit = np.stack(
            [
                np.asarray(current(mapped_fit), dtype=float).reshape(-1, 3)
                for current in currents
            ],
            axis=1,
        )
        reference_fit = sampled_fit * fit_measure[:, None, None]
        bernstein_coefficients = fit_operator @ reference_fit.reshape(
            len(fit_points), -1
        )

        mapped_validation = trafo(validation_rule)
        validation_measure = np.abs(
            np.asarray(
                jacobian_determinant(mapped_validation),
                dtype=float,
            ).reshape(-1)
        )
        sampled_validation = np.stack(
            [
                np.asarray(current(mapped_validation), dtype=float).reshape(-1, 3)
                for current in currents
            ],
            axis=1,
        )
        reference_validation = sampled_validation * validation_measure[:, None, None]
        flat_validation = reference_validation.reshape(len(validation_points), -1)
        reconstructed = validation_vandermonde @ bernstein_coefficients
        residual_sq += float(
            np.sum(
                validation_weights[:, None]
                * (reconstructed - flat_validation) ** 2
            )
        )
        field_sq += float(
            np.sum(validation_weights[:, None] * flat_validation**2)
        )
        local_coefficients = bernstein_to_monomial.T @ bernstein_coefficients
        coefficients[:, cell, :, :] = local_coefficients.reshape(
            len(monomials), n_modes, 3
        ).transpose(1, 0, 2)

        mapped_points = np.asarray(
            coordinate(mapped_validation), dtype=float
        ).reshape(-1, 3)
        geometry_error = max(
            geometry_error,
            float(
                np.max(
                    np.linalg.norm(
                        mapped_points - _p2_tet_map(nodes, validation_points),
                        axis=1,
                    )
                )
            ),
        )

    del gridfunctions
    return {
        "cell_nodes": cell_nodes,
        "cell_vertices": np.asarray(cell_vertices, dtype=np.int32),
        "coefficients": coefficients,
        "exponents": np.asarray(monomials, dtype=np.int32),
        "relative_residual": float(
            np.sqrt(residual_sq / max(field_sq, np.finfo(float).tiny))
        ),
        "relative_geometry_error": float(
            geometry_error / max(geometry_scale, np.finfo(float).tiny)
        ),
        "max_vandermonde_condition": float(condition),
        "cell_count": len(elements),
        "mode_count": n_modes,
        "charge_count": len(elements) * len(monomials),
    }


def _curved_reference_density_matrix(
    projected,
    *,
    outer_quad: int,
    curve_gauss: int,
    far_quad: int,
    ho_far_factor: float,
    eps: float,
    leafsize: int,
    eta: float,
) -> np.ndarray:
    """Contract three scalar curved reference-density Grams."""

    outer_quad = int(outer_quad)
    curve_gauss = int(curve_gauss)
    far_quad = int(far_quad)
    leafsize = int(leafsize)
    eps = float(eps)
    eta = float(eta)
    ho_far_factor = float(ho_far_factor)
    if outer_quad < 1 or curve_gauss < 1 or far_quad < 1:
        raise ValueError("curved HCurl quadrature orders must be positive")
    if leafsize < 1:
        raise ValueError("hmatrix_leafsize must be positive")
    if not np.isfinite(eps) or eps <= 0.0:
        raise ValueError("hmatrix_eps must be positive")
    if not np.isfinite(eta) or eta <= 0.0:
        raise ValueError("hmatrix_eta must be positive")
    if np.isnan(ho_far_factor) or ho_far_factor <= 0.0:
        raise ValueError("ho_far_factor must be positive")

    outer_tet_points, outer_tet_weights = _outer_tet(outer_quad)
    outer_tri_points, outer_tri_weights = _outer_tri(outer_quad)
    far_tet_points, far_tet_weights = _tet_ref(far_quad)
    far_tri_points, far_tri_weights = _tri_ref(far_quad)
    curve_gl, curve_gw = _g01(curve_gauss)
    n_cells = int(projected["cell_count"])
    exponents = np.asarray(projected["exponents"], dtype=np.int32)
    n_monomials = len(exponents)
    charge_host = np.repeat(np.arange(n_cells, dtype=np.int32), n_monomials)
    charge_kind = np.zeros(len(charge_host), dtype=np.int32)
    charge_exponents = np.tile(exponents, (n_cells, 1))
    empty_f64 = np.empty(0, dtype=float)
    empty_i32 = np.empty(0, dtype=np.int32)

    gram = _rp._ChargeGramHMatrix(
        cell_nodes=_f64_buffer(projected["cell_nodes"]),
        face_nodes=empty_f64,
        cell_vertices=_i32_buffer(projected["cell_vertices"]),
        face_vertices=empty_i32,
        n_el=n_cells,
        curve_order=2,
        charge_host=_i32_buffer(charge_host),
        charge_kind=_i32_buffer(charge_kind),
        charge_expo=_i32_buffer(charge_exponents),
        ref_tet_pts=_f64_buffer(outer_tet_points),
        ref_tet_w=_f64_buffer(outer_tet_weights),
        ref_tri_pts=_f64_buffer(outer_tri_points),
        ref_tri_w=_f64_buffer(outer_tri_weights),
        curve_gl=_f64_buffer(curve_gl),
        curve_gw=_f64_buffer(curve_gw),
        ref_tet_pts_lo=_f64_buffer(far_tet_points),
        ref_tet_w_lo=_f64_buffer(far_tet_weights),
        ref_tri_pts_lo=_f64_buffer(far_tri_points),
        ref_tri_w_lo=_f64_buffer(far_tri_weights),
        ho_far_factor=ho_far_factor,
        reference_density=True,
        eps=eps,
        leaf=leafsize,
        eta=eta,
        build=True,
    )
    coefficients = np.asarray(projected["coefficients"], dtype=float)
    n_modes = int(projected["mode_count"])
    matrix = np.zeros((n_modes, n_modes), dtype=float)
    for component in range(3):
        charge_map = coefficients[:, :, :, component].reshape(
            n_modes, -1
        ).T
        applied = np.column_stack(
            [gram.matvec_sym(charge_map[:, mode]) for mode in range(n_modes)]
        )
        matrix += charge_map.T @ applied
    return 0.5 * (matrix + matrix.T)


@dataclass(frozen=True)
class HCurlTetVolumeInteraction:
    """Precomputed reduced inductance block with analytic tet self terms."""

    basis: SampledCurrentBasis
    matrix: np.ndarray
    polynomial_degree: int
    projection_relative_residual: float
    projection_tolerance: float
    projection_quadrature_points: int
    outer_quadrature_points: int
    cell_count: int
    max_vandermonde_condition: float
    mu: float = MU0
    name: str = "analytic-affine-tet-hcurl"

    def __post_init__(self) -> None:
        if not isinstance(self.basis, SampledCurrentBasis) or self.basis.kind != "volume":
            raise TypeError("basis must be a volume SampledCurrentBasis")
        matrix = np.asarray(self.matrix)
        if matrix.shape != (self.basis.n_modes, self.basis.n_modes):
            raise ValueError("matrix shape must match the registered basis")
        if not np.all(np.isfinite(matrix)):
            raise ValueError("matrix contains non-finite values")
        object.__setattr__(self, "matrix", np.array(matrix, copy=True))

    def __call__(self, left: SampledCurrentBasis, right: SampledCurrentBasis) -> np.ndarray:
        if left is not self.basis or right is not self.basis:
            raise ValueError("analytic tet interaction is registered for one volume basis")
        return np.array(self.matrix, copy=True)

    def diagnostics(self) -> dict[str, object]:
        eigenvalues = np.linalg.eigvalsh(np.real_if_close(self.matrix))
        return {
            "kind": type(self).__name__,
            "backend": self.name,
            "singular_self_treatment": "analytic-reference-moments-through-degree-6",
            "cell_count": int(self.cell_count),
            "mode_count": int(self.basis.n_modes),
            "polynomial_degree": int(self.polynomial_degree),
            "projection_relative_residual": float(self.projection_relative_residual),
            "projection_tolerance": float(self.projection_tolerance),
            "projection_quadrature_points": int(self.projection_quadrature_points),
            "outer_quadrature_points": int(self.outer_quadrature_points),
            "max_vandermonde_condition": float(self.max_vandermonde_condition),
            "mu_H_per_m": float(self.mu),
            "minimum_eigenvalue_H": float(eigenvalues[0]),
            "maximum_eigenvalue_H": float(eigenvalues[-1]),
            "kernel_epsilon_m": None,
            "geometry": "affine-tetrahedron",
        }


@dataclass(frozen=True)
class HCurlCellVolumeInteraction:
    """Reduced 3-D HCurl interaction assembled on analytic sub-tetrahedra."""

    basis: SampledCurrentBasis
    matrix: np.ndarray
    polynomial_degree: int
    projection_relative_residual: float
    projection_tolerance: float
    geometry_relative_residual: float
    geometry_tolerance: float
    projection_quadrature_points: int
    outer_quadrature_points: int
    cell_count: int
    subtet_count: int
    family_counts: dict[str, int]
    family_subtet_counts: dict[str, int]
    max_vandermonde_condition: float
    projection_residual_history: tuple[float, ...] = ()
    geometry_residual_history: tuple[float, ...] = ()
    subdivision_level: int = 0
    subdivision_strategy: str = "pyramid-apex-midpoint"
    required_polynomial_degree: int = 0
    degree_capped: bool = False
    geometry_order: int = 1
    max_subtets: int = 512
    max_dense_moment_pairs: int = 20_000_000
    dense_moment_pairs: int = 0
    charge_count: int = 0
    geometry_backend: str = "piecewise-affine-subtet"
    mu: float = MU0
    name: str = "analytic-high-order-cell-subtet-hcurl"

    def __post_init__(self) -> None:
        if not isinstance(self.basis, SampledCurrentBasis) or self.basis.kind != "volume":
            raise TypeError("basis must be a volume SampledCurrentBasis")
        matrix = np.asarray(self.matrix)
        if matrix.shape != (self.basis.n_modes, self.basis.n_modes):
            raise ValueError("matrix shape must match the registered basis")
        if not np.all(np.isfinite(matrix)):
            raise ValueError("matrix contains non-finite values")
        object.__setattr__(self, "matrix", np.array(matrix, copy=True))
        object.__setattr__(
            self,
            "family_counts",
            {str(name): int(count) for name, count in self.family_counts.items()},
        )
        object.__setattr__(
            self,
            "family_subtet_counts",
            {
                str(name): int(count)
                for name, count in self.family_subtet_counts.items()
            },
        )
        object.__setattr__(
            self,
            "projection_residual_history",
            tuple(float(value) for value in self.projection_residual_history),
        )
        object.__setattr__(
            self,
            "geometry_residual_history",
            tuple(float(value) for value in self.geometry_residual_history),
        )

    def __call__(self, left: SampledCurrentBasis, right: SampledCurrentBasis) -> np.ndarray:
        if left is not self.basis or right is not self.basis:
            raise ValueError("analytic cell interaction is registered for one volume basis")
        return np.array(self.matrix, copy=True)

    def diagnostics(self) -> dict[str, object]:
        eigenvalues = np.linalg.eigvalsh(np.real_if_close(self.matrix))
        return {
            "kind": type(self).__name__,
            "backend": self.name,
            "singular_self_treatment": (
                "curved-p2-reference-density-duffy"
                if self.geometry_backend == "curved-p2-reference-density"
                else "analytic-subtet-reference-moments-through-degree-18"
            ),
            "cell_count": int(self.cell_count),
            "subtet_count": int(self.subtet_count),
            "family_counts": dict(self.family_counts),
            "family_subtet_counts": dict(self.family_subtet_counts),
            "mode_count": int(self.basis.n_modes),
            "polynomial_degree": int(self.polynomial_degree),
            "required_polynomial_degree": int(self.required_polynomial_degree),
            "degree_capped": bool(self.degree_capped),
            "projection_relative_residual": float(self.projection_relative_residual),
            "projection_residual_history": self.projection_residual_history,
            "projection_tolerance": float(self.projection_tolerance),
            "geometry_relative_residual": float(self.geometry_relative_residual),
            "geometry_residual_history": self.geometry_residual_history,
            "geometry_tolerance": float(self.geometry_tolerance),
            "projection_quadrature_points_per_subtet": int(
                self.projection_quadrature_points
            ),
            "outer_quadrature_points": int(self.outer_quadrature_points),
            "max_vandermonde_condition": float(self.max_vandermonde_condition),
            "subdivision_level": int(self.subdivision_level),
            "subdivision_strategy": self.subdivision_strategy,
            "geometry_order": int(self.geometry_order),
            "max_subtets": int(self.max_subtets),
            "max_dense_moment_pairs": int(self.max_dense_moment_pairs),
            "dense_moment_pairs": int(self.dense_moment_pairs),
            "charge_count": int(self.charge_count),
            "mu_H_per_m": float(self.mu),
            "minimum_eigenvalue_H": float(eigenvalues[0]),
            "maximum_eigenvalue_H": float(eigenvalues[-1]),
            "kernel_epsilon_m": None,
            "geometry": self.geometry_backend,
        }


def NgsolveHCurlTetVolumeInteraction(
    mesh,
    fes,
    vectors,
    basis: SampledCurrentBasis,
    *,
    degree: int | None = None,
    projection_quad: int | None = None,
    outer_quad: int | None = None,
    projection_tolerance: float = 1.0e-9,
    materials=None,
    mu: float = MU0,
) -> HCurlTetVolumeInteraction:
    """Build the epsilon-free reduced volume interaction for an HCurl basis.

    The default degree is ``fes.globalorder - 1``, which exactly represents
    ``curl(T)`` on affine tetrahedra.  Supplying a lower degree is an explicit
    model reduction and is accepted only when its measured projection residual
    is below ``projection_tolerance``.
    """

    if int(mesh.dim) != 3:
        raise ValueError("HCurl tet volume interaction requires a 3-D mesh")
    curve_order = int(mesh.GetCurveOrder())
    if curve_order > 1:
        raise NotImplementedError(
            "analytic HCurl tet volume interaction currently requires affine geometry; "
            "use the sampled backend for exploration until the curved P2 moment path is selected"
        )
    if not isinstance(basis, SampledCurrentBasis) or basis.kind != "volume":
        raise TypeError("basis must be a volume SampledCurrentBasis")
    parent_order = int(getattr(fes, "globalorder", 0))
    if degree is None:
        degree = max(parent_order - 1, 0)
    degree = int(degree)
    if degree < 0 or degree > 6:
        raise ValueError("degree must be in [0, 6]")
    if projection_quad is None:
        projection_quad = max(degree + 2, 4)
    if outer_quad is None:
        outer_quad = max(degree + 1, 4)
    projection_quad = int(projection_quad)
    outer_quad = int(outer_quad)
    if projection_quad < 1 or outer_quad < 1:
        raise ValueError("projection_quad and outer_quad must be positive")
    projection_tolerance = float(projection_tolerance)
    if not np.isfinite(projection_tolerance) or projection_tolerance <= 0.0:
        raise ValueError("projection_tolerance must be positive")
    mu = float(mu)
    if not np.isfinite(mu) or mu <= 0.0:
        raise ValueError("mu must be positive")

    projected = _project_reference_currents(
        mesh,
        fes,
        vectors,
        degree=degree,
        projection_quad=projection_quad,
        materials=materials,
    )
    if projected["mode_count"] != basis.n_modes:
        raise ValueError("vectors and sampled basis have different mode counts")
    residual = float(projected["relative_residual"])
    if residual > projection_tolerance:
        raise ValueError(
            "HCurl current projection residual exceeds tolerance: "
            f"{residual:.6e} > {projection_tolerance:.6e} at degree {degree}"
        )

    outer_points, outer_weights = _outer_tet(outer_quad)
    matrix = mu * np.asarray(
        _rp._TetHCurlReducedGram(
            _f64_buffer(projected["cell_verts"]),
            _i32_buffer(projected["exponents"]),
            _f64_buffer(projected["coefficients"]),
            int(basis.n_modes),
            _f64_buffer(outer_points),
            _f64_buffer(outer_weights),
        ),
        dtype=float,
    )
    matrix = 0.5 * (matrix + matrix.T)
    return HCurlTetVolumeInteraction(
        basis=basis,
        matrix=matrix,
        polynomial_degree=degree,
        projection_relative_residual=residual,
        projection_tolerance=projection_tolerance,
        projection_quadrature_points=int(len(_tet_ref(projection_quad)[1])),
        outer_quadrature_points=int(len(outer_weights)),
        cell_count=int(projected["cell_count"]),
        max_vandermonde_condition=float(projected["max_vandermonde_condition"]),
        mu=mu,
    )


def NgsolveHCurlCellVolumeInteraction(
    mesh,
    fes,
    vectors,
    basis: SampledCurrentBasis,
    *,
    degree: int | None = None,
    projection_quad: int | None = None,
    outer_quad: int | None = None,
    projection_tolerance: float = 1.0e-4,
    geometry_tolerance: float = 1.0e-3,
    max_subdivision_levels: int = 8,
    max_subtets: int = 512,
    max_dense_moment_pairs: int = 20_000_000,
    max_charges: int = 4096,
    curve_gauss: int = 8,
    far_quad: int = 3,
    ho_far_factor: float = 2.0,
    hmatrix_eps: float = 1.0e-8,
    hmatrix_leafsize: int = 32,
    hmatrix_eta: float = 2.0,
    materials=None,
    mu: float = MU0,
) -> HCurlCellVolumeInteraction:
    """Build an epsilon-free reduced interaction for TET/HEX/WEDGE/PYRAMID.

    Each NGSolve cell is partitioned into affine tetrahedra.  A stable
    tetrahedral Bernstein fit is converted exactly to reference monomials and
    passed to analytic Newton-potential moments through total degree 18.  For
    a p-order parent, the automatic required degrees are p-1 (TET), 2p
    (WEDGE), and 3p (HEX/PYRAMID).  Required degrees above 18 use an hp path:
    analytic degree-18 leaf moments plus uniform h refinement.  PYRAMID's
    rational p<=6 apex modes use the cheaper apex-only refinement.

    The default projection tolerance 1e-4 accepts the p=6 pyramid at degree 18
    without refinement.  Tighter tolerances may invoke up to
    ``max_subdivision_levels`` apex refinements and should be paired with an
    outer-quadrature convergence check.  P2 tetrahedra instead project the
    curl-Piola reference density ``K=J*abs(det(dX/dxi))`` and use the exact P2
    geometry in Radia's curved Duffy/H-matrix Gram.  Other warped maps use the
    residual-controlled piecewise-affine loop.  ``max_subtets`` and
    ``max_charges`` prevent accidental Gram explosions.  No diagonal kernel
    epsilon is used.
    """

    if int(mesh.dim) != 3:
        raise ValueError("HCurl cell volume interaction requires a 3-D mesh")
    if not isinstance(basis, SampledCurrentBasis) or basis.kind != "volume":
        raise TypeError("basis must be a volume SampledCurrentBasis")
    inventory = NgsolveHCurlCellFamilies(mesh, materials=materials)
    parent_order = int(getattr(fes, "globalorder", 0))
    curve_order = int(mesh.GetCurveOrder())
    if curve_order == 2 and inventory.families == ("tet",):
        curved_degree = parent_order if degree is None else int(degree)
        if curved_degree < 0 or curved_degree > 18:
            raise ValueError("degree must be in [0, 18]")
        curved_projection_quad = (
            max((curved_degree + 3) // 2, 4)
            if projection_quad is None
            else int(projection_quad)
        )
        curved_outer_quad = 5 if outer_quad is None else int(outer_quad)
        if curved_projection_quad < 1 or curved_outer_quad < 1:
            raise ValueError("projection_quad and outer_quad must be positive")
        projection_tolerance = float(projection_tolerance)
        geometry_tolerance = float(geometry_tolerance)
        mu = float(mu)
        if not np.isfinite(projection_tolerance) or projection_tolerance <= 0.0:
            raise ValueError("projection_tolerance must be positive")
        if not np.isfinite(geometry_tolerance) or geometry_tolerance <= 0.0:
            raise ValueError("geometry_tolerance must be positive")
        if not np.isfinite(mu) or mu <= 0.0:
            raise ValueError("mu must be positive")
        projected = _project_curved_tet_reference_currents(
            mesh,
            fes,
            vectors,
            degree=curved_degree,
            projection_quad=curved_projection_quad,
            materials=materials,
        )
        if projected["mode_count"] != basis.n_modes:
            raise ValueError("vectors and sampled basis have different mode counts")
        residual = float(projected["relative_residual"])
        geometry_residual = float(projected["relative_geometry_error"])
        if residual > projection_tolerance:
            raise ValueError(
                "curved HCurl reference-density projection residual "
                f"{residual:.6e} exceeds tolerance {projection_tolerance:.6e} "
                f"at degree {curved_degree}"
            )
        if geometry_residual > geometry_tolerance:
            raise ValueError(
                "NGSolve and curved P2 Gram geometry disagree: residual "
                f"{geometry_residual:.6e} exceeds tolerance {geometry_tolerance:.6e}"
            )
        max_charges = int(max_charges)
        if max_charges <= 0:
            raise ValueError("max_charges must be positive")
        if int(projected["charge_count"]) > max_charges:
            raise ValueError(
                "curved HCurl reference-density Gram would require "
                f"{projected['charge_count']} scalar charges, "
                f"above max_charges={max_charges}"
            )
        matrix = mu * _curved_reference_density_matrix(
            projected,
            outer_quad=curved_outer_quad,
            curve_gauss=int(curve_gauss),
            far_quad=int(far_quad),
            ho_far_factor=float(ho_far_factor),
            eps=float(hmatrix_eps),
            leafsize=int(hmatrix_leafsize),
            eta=float(hmatrix_eta),
        )
        outer_weights = _outer_tet(curved_outer_quad)[1]
        return HCurlCellVolumeInteraction(
            basis=basis,
            matrix=matrix,
            polynomial_degree=curved_degree,
            projection_relative_residual=residual,
            projection_tolerance=projection_tolerance,
            geometry_relative_residual=geometry_residual,
            geometry_tolerance=geometry_tolerance,
            projection_quadrature_points=int(
                len(_tet_ref(curved_projection_quad)[1])
            ),
            outer_quadrature_points=int(len(outer_weights)),
            cell_count=int(projected["cell_count"]),
            subtet_count=int(projected["cell_count"]),
            family_counts=dict(inventory.counts),
            family_subtet_counts={"tet": int(projected["cell_count"])},
            max_vandermonde_condition=float(
                projected["max_vandermonde_condition"]
            ),
            projection_residual_history=(residual,),
            geometry_residual_history=(geometry_residual,),
            subdivision_level=0,
            subdivision_strategy="curved-p2-reference-density",
            required_polynomial_degree=parent_order,
            degree_capped=curved_degree < parent_order,
            geometry_order=curve_order,
            max_subtets=int(max_subtets),
            max_dense_moment_pairs=int(max_dense_moment_pairs),
            dense_moment_pairs=0,
            charge_count=int(projected["charge_count"]),
            geometry_backend="curved-p2-reference-density",
            mu=mu,
            name="curved-p2-reference-density-hcurl-hmatrix",
        )
    if curve_order > 2 and inventory.families == ("tet",):
        raise NotImplementedError(
            "exact curved HCurl reference-density integration currently supports geometry order 2"
        )
    required = {
        "tet": max(parent_order - 1, 0),
        "wedge": 2 * parent_order,
        "hex": 3 * parent_order,
        "pyramid": 3 * parent_order,
    }
    required_degree = max(required[family] for family in inventory.families)
    if degree is None:
        degree = min(required_degree, 18)
    degree = int(degree)
    if degree < 0 or degree > 18:
        raise ValueError("degree must be in [0, 18]")
    degree_capped = degree < required_degree
    if projection_quad is None:
        projection_quad = max((degree + 3) // 2, 4)
    if outer_quad is None:
        outer_quad = 4
    projection_quad = int(projection_quad)
    outer_quad = int(outer_quad)
    if projection_quad < 1 or outer_quad < 1:
        raise ValueError("projection_quad and outer_quad must be positive")
    max_subdivision_levels = int(max_subdivision_levels)
    if max_subdivision_levels < 0:
        raise ValueError("max_subdivision_levels must be non-negative")
    max_subtets = int(max_subtets)
    if max_subtets <= 0:
        raise ValueError("max_subtets must be positive")
    projection_tolerance = float(projection_tolerance)
    geometry_tolerance = float(geometry_tolerance)
    if not np.isfinite(projection_tolerance) or projection_tolerance <= 0.0:
        raise ValueError("projection_tolerance must be positive")
    if not np.isfinite(geometry_tolerance) or geometry_tolerance <= 0.0:
        raise ValueError("geometry_tolerance must be positive")
    mu = float(mu)
    if not np.isfinite(mu) or mu <= 0.0:
        raise ValueError("mu must be positive")

    residual_history: list[float] = []
    projected = None
    geometry_history: list[float] = []
    refinement_strategy = (
        "uniform"
        if curve_order > 1 or degree_capped
        else "pyramid-apex"
    )
    for subdivision_level in range(max_subdivision_levels + 1):
        projected = _project_cell_currents_to_subtets(
            mesh,
            fes,
            vectors,
            degree=degree,
            projection_quad=projection_quad,
            subdivision_level=subdivision_level,
            refinement_strategy=refinement_strategy,
            max_subtets=max_subtets,
            max_dense_moment_pairs=max_dense_moment_pairs,
            materials=materials,
        )
        if projected["mode_count"] != basis.n_modes:
            raise ValueError("vectors and sampled basis have different mode counts")
        geometry_residual = float(projected["relative_geometry_error"])
        residual = float(projected["relative_residual"])
        residual_history.append(residual)
        geometry_history.append(geometry_residual)
        if (
            residual <= projection_tolerance
            and geometry_residual <= geometry_tolerance
        ):
            break
        if subdivision_level == 0 and geometry_residual > geometry_tolerance:
            refinement_strategy = "uniform"
    assert projected is not None
    if residual > projection_tolerance or geometry_residual > geometry_tolerance:
        raise ValueError(
            "HCurl hp cell projection did not reach both tolerances after "
            f"{projected['subdivision_level']} subdivision levels: "
            f"current {residual:.6e} (tol {projection_tolerance:.6e}), "
            f"geometry {geometry_residual:.6e} (tol {geometry_tolerance:.6e}), "
            f"degree {degree}, strategy {refinement_strategy}"
        )

    outer_points, outer_weights = _outer_tet(outer_quad)
    matrix = mu * np.asarray(
        _rp._TetHCurlReducedGram(
            _f64_buffer(projected["cell_verts"]),
            _i32_buffer(projected["exponents"]),
            _f64_buffer(projected["coefficients"]),
            int(basis.n_modes),
            _f64_buffer(outer_points),
            _f64_buffer(outer_weights),
        ),
        dtype=float,
    )
    matrix = 0.5 * (matrix + matrix.T)
    return HCurlCellVolumeInteraction(
        basis=basis,
        matrix=matrix,
        polynomial_degree=degree,
        projection_relative_residual=residual,
        projection_tolerance=projection_tolerance,
        geometry_relative_residual=geometry_residual,
        geometry_tolerance=geometry_tolerance,
        projection_quadrature_points=int(len(_tet_ref(projection_quad)[1])),
        outer_quadrature_points=int(len(outer_weights)),
        cell_count=int(projected["cell_count"]),
        subtet_count=int(projected["subtet_count"]),
        family_counts=dict(inventory.counts),
        family_subtet_counts=dict(projected["family_subtet_counts"]),
        max_vandermonde_condition=float(projected["max_vandermonde_condition"]),
        projection_residual_history=tuple(residual_history),
        geometry_residual_history=tuple(geometry_history),
        subdivision_level=int(projected["subdivision_level"]),
        subdivision_strategy=str(projected["refinement_strategy"]),
        required_polynomial_degree=required_degree,
        degree_capped=degree_capped,
        geometry_order=curve_order,
        max_subtets=max_subtets,
        max_dense_moment_pairs=max_dense_moment_pairs,
        dense_moment_pairs=int(projected["dense_moment_pairs"]),
        mu=mu,
    )


__all__ = [
    "HCurlTetVolumeInteraction",
    "HCurlCellVolumeInteraction",
    "NgsolveHCurlTetVolumeInteraction",
    "NgsolveHCurlCellVolumeInteraction",
]
