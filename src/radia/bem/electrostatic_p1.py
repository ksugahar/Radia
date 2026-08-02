"""Readable P1 Galerkin BEM helpers for electrostatic conductors."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .sibc_hacapk import assemble_SL_dense


EPSILON_0 = 8.8541878128e-12


@dataclass(frozen=True)
class ElectrostaticP1Result:
    """Surface-charge solution for prescribed conductor potentials."""

    charge_density_coefficients: np.ndarray
    nodal_area_weights: np.ndarray
    potential_rhs: np.ndarray
    single_layer_matrix: np.ndarray

    @property
    def total_charge_c(self) -> float:
        return float(self.nodal_area_weights @ self.charge_density_coefficients)

    def charge_on_vertices(self, vertex_ids: np.ndarray) -> float:
        ids = np.asarray(vertex_ids, dtype=np.int64)
        return float(
            self.nodal_area_weights[ids]
            @ self.charge_density_coefficients[ids]
        )


def _surface_mass_matrix(vertices: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    mass = np.zeros((len(vertices), len(vertices)), dtype=float)
    local = np.array(
        [[2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]],
        dtype=float,
    )
    for tri in triangles:
        points = vertices[tri]
        area = 0.5 * np.linalg.norm(
            np.cross(points[1] - points[0], points[2] - points[0])
        )
        if not np.isfinite(area) or area <= 0.0:
            raise ValueError("triangles must have positive finite area")
        mass[np.ix_(tri, tri)] += (area / 12.0) * local
    return mass


def solve_prescribed_potential_p1(
    vertices: np.ndarray,
    triangles: np.ndarray,
    nodal_potential_v: np.ndarray,
    *,
    permittivity_f_per_m: float = EPSILON_0,
    regular_quad_degree: int = 4,
    singular_n_q: int = 5,
) -> ElectrostaticP1Result:
    """Solve ``V = S sigma / epsilon`` on a triangulated conductor surface.

    The unknown charge density uses continuous P1 basis functions.  Each
    disconnected conductor can carry an independently prescribed constant
    potential by assigning that value to its surface vertices.
    """
    verts = np.asarray(vertices, dtype=float)
    triangle_indices = np.asarray(triangles)
    potential = np.asarray(nodal_potential_v, dtype=float)
    if verts.ndim != 2 or verts.shape[1] != 3:
        raise ValueError("vertices must have shape (n, 3)")
    if (
        triangle_indices.ndim != 2
        or triangle_indices.shape[1] != 3
        or len(triangle_indices) == 0
    ):
        raise ValueError("triangles must have shape (m, 3) with m > 0")
    if np.issubdtype(triangle_indices.dtype, np.bool_) or not np.issubdtype(
        triangle_indices.dtype, np.integer
    ):
        raise ValueError("triangle vertex indices must be integers")
    if potential.shape != (len(verts),):
        raise ValueError("nodal_potential_v must have one value per vertex")
    if not np.all(np.isfinite(verts)) or not np.all(np.isfinite(potential)):
        raise ValueError("vertices and potentials must be finite")
    if np.min(triangle_indices) < 0 or np.max(triangle_indices) >= len(verts):
        raise ValueError("triangle vertex index is out of range")
    tris = triangle_indices.astype(np.int64, copy=False)
    epsilon = float(permittivity_f_per_m)
    if not np.isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError("permittivity_f_per_m must be positive and finite")

    mass = _surface_mass_matrix(verts, tris)
    single_layer = assemble_SL_dense(
        verts,
        tris,
        regular_quad_degree=regular_quad_degree,
        include_singular=True,
        singular_n_q=singular_n_q,
    )
    rhs = epsilon * (mass @ potential)
    coefficients = np.linalg.solve(single_layer, rhs)
    nodal_area_weights = mass @ np.ones(len(verts), dtype=float)
    return ElectrostaticP1Result(
        charge_density_coefficients=coefficients,
        nodal_area_weights=nodal_area_weights,
        potential_rhs=rhs,
        single_layer_matrix=single_layer,
    )
