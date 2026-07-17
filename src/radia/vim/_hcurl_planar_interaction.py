"""Planar log-kernel interaction for reduced TRIG/QUAD HCurl currents."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

import radia._radia_pybind as _rp

from ._eddy_hybrid import MU0, SampledCurrentBasis
from ._vim import (
    _EMPTY_F64,
    _EMPTY_I32,
    _SYM5_TRI,
    _f64_buffer,
    _fit_geometry_map_2d,
    _g01,
    _i32_buffer,
    _prod_tri01,
)


def _labels(value) -> set[str] | None:
    if value is None:
        return None
    if isinstance(value, str):
        return {value}
    return {str(item) for item in value}


def _monomials(family: str, degree: int) -> tuple[tuple[int, int], ...]:
    if family == "trig":
        return tuple(
            (i, j)
            for i in range(degree + 1)
            for j in range(degree + 1 - i)
        )
    return tuple(
        (i, j)
        for i in range(degree + 1)
        for j in range(degree + 1)
    )


def _gridfunction_currents(fes, vectors):
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
                raise ValueError("planar HCurl interaction currently requires real basis vectors")
        target[:] = np.real(values[:, column])
        gridfunctions.append(gf)
        currents.append(ng.curl(gf))
    return gridfunctions, currents


def _project_reference_current_density(
    mesh,
    fes,
    vectors,
    *,
    degree: int,
    projection_order: int,
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
        raise ValueError("the selected conductor region contains no planar cells")
    gridfunctions, currents = _gridfunction_currents(fes, vectors)
    n_modes = len(currents)
    geometry_order = max(1, int(mesh.GetCurveOrder()))
    geometry_slots = (geometry_order + 1) ** 2
    cell_maps = []
    cell_types = []
    charge_host = []
    charge_kind = []
    charge_exponents = []
    coefficient_rows = []
    family_counts: dict[str, int] = {}
    residual_sq = 0.0
    field_sq = 0.0
    max_condition = 0.0

    for cell, element in enumerate(elements):
        vertex_count = len(element.vertices)
        if vertex_count == 3:
            family = "trig"
            cell_type = 0
        elif vertex_count == 4:
            family = "quad"
            cell_type = 1
        else:
            raise NotImplementedError(
                f"planar HCurl interaction does not support {element.type}"
            )
        family_counts[family] = family_counts.get(family, 0) + 1
        cell_types.append(cell_type)
        element_id = ng.ElementId(ng.VOL, int(element.nr))
        fitted_map = _fit_geometry_map_2d(
            mesh,
            element_id,
            cell_type,
            geometry_order,
        )
        map_slot = np.zeros((geometry_slots, 2), dtype=float)
        map_slot[: fitted_map.shape[0]] = fitted_map
        cell_maps.append(map_slot)

        exponents = _monomials(family, degree)
        rule = ng.IntegrationRule(element.type, projection_order)
        reference_points = np.asarray(
            [(float(ip.point[0]), float(ip.point[1])) for ip in rule],
            dtype=float,
        )
        reference_weights = np.asarray([float(ip.weight) for ip in rule])
        vandermonde = np.asarray(
            [
                [point[0] ** i * point[1] ** j for i, j in exponents]
                for point in reference_points
            ],
            dtype=float,
        )
        weighted_vandermonde = vandermonde * np.sqrt(reference_weights)[:, None]
        condition = float(np.linalg.cond(weighted_vandermonde))
        max_condition = max(max_condition, condition)
        if not np.isfinite(condition):
            raise ValueError(f"singular planar current projection at cell {cell}")

        trafo = mesh.GetTrafo(element)
        sampled = np.empty((len(rule), n_modes), dtype=float)
        for point_index, ip in enumerate(rule):
            mapped = trafo(ip)
            jacobian_measure = float(mapped.measure)
            for mode, current in enumerate(currents):
                value = np.asarray(current(mapped), dtype=float).reshape(-1)
                if value.shape != (1,):
                    raise ValueError(
                        "A planar HCurl curl basis must evaluate to one scalar Jz value."
                    )
                sampled[point_index, mode] = float(value[0]) * jacobian_measure
        weighted_sampled = sampled * np.sqrt(reference_weights)[:, None]
        coefficients, *_ = np.linalg.lstsq(
            weighted_vandermonde,
            weighted_sampled,
            rcond=1.0e-13,
        )
        reconstructed = vandermonde @ coefficients
        residual_sq += float(
            np.sum(reference_weights[:, None] * (reconstructed - sampled) ** 2)
        )
        field_sq += float(np.sum(reference_weights[:, None] * sampled**2))
        for exponent, row in zip(exponents, coefficients):
            charge_host.append(cell)
            charge_kind.append(0)
            charge_exponents.extend((exponent[0], exponent[1], 0))
            coefficient_rows.append(row)

    del gridfunctions
    return {
        "cell_map": _f64_buffer(np.concatenate([item.ravel() for item in cell_maps])),
        "cell_type": _i32_buffer(cell_types),
        "charge_host": _i32_buffer(charge_host),
        "charge_kind": _i32_buffer(charge_kind),
        "charge_exponents": _i32_buffer(charge_exponents),
        "coefficients": np.asarray(coefficient_rows, dtype=float),
        "n_cells": len(elements),
        "n_modes": n_modes,
        "geometry_order": geometry_order,
        "family_counts": family_counts,
        "projection_relative_residual": float(
            np.sqrt(residual_sq / max(field_sq, np.finfo(float).tiny))
        ),
        "max_vandermonde_condition": max_condition,
    }


@dataclass(frozen=True)
class HCurlPlanarVolumeInteraction:
    """Reduced TRIG/QUAD inductance using the planar log Green kernel."""

    basis: SampledCurrentBasis
    matrix: np.ndarray
    polynomial_degree: int
    projection_relative_residual: float
    projection_tolerance: float
    projection_order: int
    geometry_order: int
    cell_count: int
    charge_mode_count: int
    family_counts: dict[str, int]
    max_vandermonde_condition: float
    net_current_per_mode: np.ndarray
    mu: float = MU0
    name: str = "analytic-planar-log-hcurl"

    def __post_init__(self) -> None:
        matrix = np.asarray(self.matrix)
        if matrix.shape != (self.basis.n_modes, self.basis.n_modes):
            raise ValueError("matrix shape must match the registered basis")
        if not np.all(np.isfinite(matrix)):
            raise ValueError("matrix contains non-finite values")
        object.__setattr__(self, "matrix", np.array(matrix, copy=True))
        object.__setattr__(
            self,
            "net_current_per_mode",
            np.asarray(self.net_current_per_mode, dtype=float).copy(),
        )

    def __call__(self, left: SampledCurrentBasis, right: SampledCurrentBasis) -> np.ndarray:
        if left is not self.basis or right is not self.basis:
            raise ValueError("planar interaction is registered for one volume basis")
        return np.array(self.matrix, copy=True)

    def diagnostics(self) -> dict[str, object]:
        eigenvalues = np.linalg.eigvalsh(np.real_if_close(self.matrix))
        return {
            "kind": type(self).__name__,
            "backend": self.name,
            "kernel": "-log(r / 1m)/(2*pi)",
            "cell_count": int(self.cell_count),
            "charge_mode_count": int(self.charge_mode_count),
            "family_counts": dict(self.family_counts),
            "mode_count": int(self.basis.n_modes),
            "polynomial_degree": int(self.polynomial_degree),
            "projection_relative_residual": float(self.projection_relative_residual),
            "projection_tolerance": float(self.projection_tolerance),
            "projection_order": int(self.projection_order),
            "geometry_order": int(self.geometry_order),
            "max_vandermonde_condition": float(self.max_vandermonde_condition),
            "net_current_per_mode_A_per_m": self.net_current_per_mode.tolist(),
            "mu_H_per_m": float(self.mu),
            "minimum_eigenvalue_H_per_m": float(eigenvalues[0]),
            "maximum_eigenvalue_H_per_m": float(eigenvalues[-1]),
            "kernel_epsilon_m": None,
            "geometry": "planar-trig-quad",
        }


def NgsolveHCurlPlanarVolumeInteraction(
    mesh,
    fes,
    vectors,
    basis: SampledCurrentBasis,
    *,
    degree: int | None = None,
    projection_order: int | None = None,
    projection_tolerance: float = 1.0e-9,
    materials=None,
    mu: float = MU0,
    eps: float = 1.0e-12,
    leafsize: int = 64,
    eta: float = 2.0,
) -> HCurlPlanarVolumeInteraction:
    """Build the epsilon-free planar interaction for TRIG/QUAD HCurl modes."""

    if int(mesh.dim) != 2:
        raise ValueError("planar HCurl interaction requires a 2-D mesh")
    if not isinstance(basis, SampledCurrentBasis) or basis.kind != "volume":
        raise TypeError("basis must be a volume SampledCurrentBasis")
    parent_order = int(getattr(fes, "globalorder", 0))
    if degree is None:
        degree = max(parent_order, 1)
    degree = int(degree)
    if degree < 0:
        raise ValueError("degree must be non-negative")
    if projection_order is None:
        projection_order = max(2 * degree + 4, 6)
    projection_order = int(projection_order)
    if projection_order < 1:
        raise ValueError("projection_order must be positive")
    projection_tolerance = float(projection_tolerance)
    if projection_tolerance <= 0.0 or not np.isfinite(projection_tolerance):
        raise ValueError("projection_tolerance must be positive")
    mu = float(mu)
    if mu <= 0.0 or not np.isfinite(mu):
        raise ValueError("mu must be positive")

    projected = _project_reference_current_density(
        mesh,
        fes,
        vectors,
        degree=degree,
        projection_order=projection_order,
        materials=materials,
    )
    if projected["n_modes"] != basis.n_modes:
        raise ValueError("vectors and sampled basis have different mode counts")
    residual = float(projected["projection_relative_residual"])
    if residual > projection_tolerance:
        raise ValueError(
            "planar HCurl current projection residual exceeds tolerance: "
            f"{residual:.6e} > {projection_tolerance:.6e} at degree {degree}"
        )

    outer_points, outer_weights = _prod_tri01(4 if degree <= 1 else 6)
    quad_points, quad_weights = _g01(4 if degree <= 1 else 6)
    edge_points, edge_weights = _g01(12 if degree <= 1 else 16)
    inner_points, inner_weights = _g01(8 if degree <= 1 else 10)
    gram = _rp._ChargeGramHMatrix(
        dim2=2,
        geometry_order=int(projected["geometry_order"]),
        cell_map=projected["cell_map"],
        cell_type=projected["cell_type"],
        edge_map=_EMPTY_F64,
        n_el=int(projected["n_cells"]),
        n_be=0,
        charge_host=projected["charge_host"],
        charge_kind=projected["charge_kind"],
        charge_expo=projected["charge_exponents"],
        sym_tri_pts=_f64_buffer(outer_points),
        sym_tri_w=_f64_buffer(outer_weights),
        gl_quad=_f64_buffer(quad_points),
        gw_quad=_f64_buffer(quad_weights),
        gl_edge=_f64_buffer(edge_points),
        gw_edge=_f64_buffer(edge_weights),
        gl_in=_f64_buffer(inner_points),
        gw_in=_f64_buffer(inner_weights),
        far_tri_pts=_f64_buffer(_SYM5_TRI[0]),
        far_tri_w=_f64_buffer(_SYM5_TRI[1]),
        image_masks=_EMPTY_I32,
        image_signs=_EMPTY_F64,
        eps=float(eps),
        leaf=max(64, int(leafsize)),
        eta=float(eta),
        build=True,
    )
    coefficients = projected["coefficients"]
    applied = np.column_stack(
        [np.asarray(gram.matvec_sym(coefficients[:, mode])) for mode in range(basis.n_modes)]
    )
    matrix = mu * (coefficients.T @ applied)
    matrix = 0.5 * (matrix + matrix.T)
    net_current = np.einsum("ai,i->a", basis.modes[:, :, 2], basis.weights)
    return HCurlPlanarVolumeInteraction(
        basis=basis,
        matrix=matrix,
        polynomial_degree=degree,
        projection_relative_residual=residual,
        projection_tolerance=projection_tolerance,
        projection_order=projection_order,
        geometry_order=int(projected["geometry_order"]),
        cell_count=int(projected["n_cells"]),
        charge_mode_count=int(coefficients.shape[0]),
        family_counts=dict(projected["family_counts"]),
        max_vandermonde_condition=float(projected["max_vandermonde_condition"]),
        net_current_per_mode=net_current,
        mu=mu,
    )


__all__ = ["HCurlPlanarVolumeInteraction", "NgsolveHCurlPlanarVolumeInteraction"]
