"""Hybrid eddy-current VIM primitives.

This module is the small, composable core for the motor-oriented eddy-current
idea: build a reduced current basis from a high-order HCurl/T space, enrich it
with a surface-Omega/SIBC basis, then assemble the VIM interaction matrices on
that reduced basis.

The code intentionally works on sampled basis functions.  NGSolve owns the
high-order FE construction and projection; Radia owns the sampled VIM operator
and the reduced circuit-facing matrices.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import hashlib
import json
import os
from pathlib import Path
from typing import Iterable

import numpy as np


MU0 = 4.0e-7 * np.pi


def SkinDepth(frequency_hz: float, sigma: float, mu: float = MU0) -> float:
    """Return the plane-conductor skin depth in metres.

    ``delta = sqrt(2 / (omega * mu * sigma))``.  This helper is deliberately
    separate from :func:`SkinImpedance`: adjacency to air makes a face an SIBC
    *candidate*, while the thickness-to-skin-depth ratio decides whether the
    half-space approximation is admissible for the body at this frequency.
    """

    frequency = float(frequency_hz)
    conductivity = float(sigma)
    permeability = float(mu)
    if not np.isfinite(frequency) or frequency <= 0.0:
        raise ValueError("frequency_hz must be positive")
    if not np.isfinite(conductivity) or conductivity <= 0.0:
        raise ValueError("sigma must be positive")
    if not np.isfinite(permeability) or permeability <= 0.0:
        raise ValueError("mu must be positive")
    omega = 2.0 * np.pi * frequency
    return float(np.sqrt(2.0 / (omega * permeability * conductivity)))


@dataclass(frozen=True)
class EddySIBCApplicability:
    """Body-level gate between volumetric VIM and half-space SIBC.

    Face adjacency remains the topology gate: only conductor-air/exterior
    faces can ever carry SIBC modes.  This second, physical gate prevents a
    thin conductor whose opposite skin layers overlap from being converted to
    a surface-only model merely because all its faces touch air.

    ``min_thickness_to_skin_depth`` is an explicit engineering tolerance, not
    a claim of a universal transition.  The default value three is
    conservative enough for routing; convergence against the volumetric path
    remains the acceptance criterion near the transition.
    """

    frequency_hz: float
    sigma: float
    characteristic_thickness_m: float
    mu: float = MU0
    min_thickness_to_skin_depth: float = 3.0
    characteristic_curvature_radius_m: float | None = None
    min_curvature_radius_to_skin_depth: float = 3.0

    def __post_init__(self) -> None:
        values = {
            "frequency_hz": self.frequency_hz,
            "sigma": self.sigma,
            "characteristic_thickness_m": self.characteristic_thickness_m,
            "mu": self.mu,
            "min_thickness_to_skin_depth": self.min_thickness_to_skin_depth,
            "min_curvature_radius_to_skin_depth": self.min_curvature_radius_to_skin_depth,
        }
        for name, value in values.items():
            number = float(value)
            if not np.isfinite(number) or number <= 0.0:
                raise ValueError(f"{name} must be positive")
            object.__setattr__(self, name, number)
        if self.characteristic_curvature_radius_m is not None:
            radius = float(self.characteristic_curvature_radius_m)
            if not np.isfinite(radius) or radius <= 0.0:
                raise ValueError("characteristic_curvature_radius_m must be positive")
            object.__setattr__(self, "characteristic_curvature_radius_m", radius)

    @property
    def skin_depth_m(self) -> float:
        return SkinDepth(self.frequency_hz, self.sigma, self.mu)

    @property
    def thickness_to_skin_depth(self) -> float:
        return float(self.characteristic_thickness_m / self.skin_depth_m)

    @property
    def curvature_radius_to_skin_depth(self) -> float | None:
        if self.characteristic_curvature_radius_m is None:
            return None
        return float(self.characteristic_curvature_radius_m / self.skin_depth_m)

    @property
    def thickness_separated(self) -> bool:
        return self.thickness_to_skin_depth >= self.min_thickness_to_skin_depth

    @property
    def curvature_resolved(self) -> bool:
        ratio = self.curvature_radius_to_skin_depth
        return ratio is None or ratio >= self.min_curvature_radius_to_skin_depth

    @property
    def sibc_applicable(self) -> bool:
        return self.thickness_separated and self.curvature_resolved

    @property
    def selected_model(self) -> str:
        return "sibc" if self.sibc_applicable else "volumetric"

    def diagnostics(self) -> dict[str, float | bool | str | None]:
        return {
            "frequency_hz": self.frequency_hz,
            "sigma_S_per_m": self.sigma,
            "mu_H_per_m": self.mu,
            "skin_depth_m": self.skin_depth_m,
            "characteristic_thickness_m": self.characteristic_thickness_m,
            "thickness_to_skin_depth": self.thickness_to_skin_depth,
            "min_thickness_to_skin_depth": self.min_thickness_to_skin_depth,
            "characteristic_curvature_radius_m": self.characteristic_curvature_radius_m,
            "curvature_radius_to_skin_depth": self.curvature_radius_to_skin_depth,
            "min_curvature_radius_to_skin_depth": self.min_curvature_radius_to_skin_depth,
            "thickness_separated": self.thickness_separated,
            "curvature_resolved": self.curvature_resolved,
            "sibc_applicable": self.sibc_applicable,
            "selected_model": self.selected_model,
        }


def _nonnegative_int(value, name: str) -> int:
    out = int(value)
    if out < 0:
        raise ValueError(f"{name} must be non-negative")
    return out


def EddyTracePolynomialDim(degree: int, *, face_family: str = "simplex") -> int:
    """Return the number of scalar trace polynomials on one face.

    ``face_family="simplex"`` uses ``dim P_r(face) = (r+1)(r+2)/2`` for a
    triangular face.  ``face_family="tensor"`` uses ``dim Q_r(face)=(r+1)^2``
    for a tensor-product quadrilateral face.  The helper is intentionally
    separate from the reduced basis builder: it is an order ledger, not a
    numerical claim that a particular p is optimal.
    """

    r = _nonnegative_int(degree, "degree")
    if face_family == "simplex":
        return (r + 1) * (r + 2) // 2
    if face_family == "tensor":
        return (r + 1) * (r + 1)
    raise ValueError("face_family must be 'simplex' or 'tensor'")


@dataclass(frozen=True)
class EddyParentOrderLedger:
    """Symbolic parent-order ledger for topology-aware eddy bubbling.

    ``p`` is treated as a parent-space admissibility order.  It is not the CLN
    stage count and it is not automatically optimal.  The minimum admissible
    parent order is

    ``max(bulk_degree, bridge_trace_degree, surface_current_degree)``.

    The surface-Omega scalar potential degree is one higher than the retained
    surface-current degree because ``K = n x grad_Gamma Omega``.
    """

    bulk_degree: int
    bridge_trace_degree: int = 0
    surface_current_degree: int = 0
    face_family: str = "simplex"

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "bulk_degree",
            _nonnegative_int(self.bulk_degree, "bulk_degree"),
        )
        object.__setattr__(
            self,
            "bridge_trace_degree",
            _nonnegative_int(self.bridge_trace_degree, "bridge_trace_degree"),
        )
        object.__setattr__(
            self,
            "surface_current_degree",
            _nonnegative_int(self.surface_current_degree, "surface_current_degree"),
        )
        EddyTracePolynomialDim(self.bridge_trace_degree, face_family=self.face_family)

    @property
    def required_parent_order(self) -> int:
        """Smallest parent order admitted by the symbolic ledger."""

        return max(
            self.bulk_degree,
            self.bridge_trace_degree,
            self.surface_current_degree,
        )

    @property
    def surface_omega_degree(self) -> int:
        """Scalar surface-Omega degree required for the retained current trace."""

        return self.surface_current_degree + 1

    @property
    def bridge_trace_dim(self) -> int:
        """Number of bridge trace modes per conductor-graph cycle."""

        return EddyTracePolynomialDim(
            self.bridge_trace_degree,
            face_family=self.face_family,
        )

    def is_parent_order_admissible(self, parent_order: int) -> bool:
        """Return whether ``parent_order`` satisfies the ledger."""

        p = _nonnegative_int(parent_order, "parent_order")
        return p >= self.required_parent_order

    def bridge_modes(self, cycle_rank: int) -> int:
        """Return bridge modes implied by the cycle rank and trace degree."""

        return _nonnegative_int(cycle_rank, "cycle_rank") * self.bridge_trace_dim

    def estimated_reduced_modes(
        self,
        *,
        evrs_rank: int,
        cycle_rank: int,
        surface_modes: int,
        non_sibc_trace_modes: int = 0,
    ) -> int:
        """Return the ledger mode count for a class-wise reduced basis."""

        return (
            _nonnegative_int(evrs_rank, "evrs_rank")
            + self.bridge_modes(cycle_rank)
            + _nonnegative_int(surface_modes, "surface_modes")
            + _nonnegative_int(non_sibc_trace_modes, "non_sibc_trace_modes")
        )

    def diagnostics(
        self,
        *,
        parent_order: int | None = None,
        evrs_rank: int | None = None,
        cycle_rank: int | None = None,
        surface_modes: int | None = None,
        non_sibc_trace_modes: int = 0,
    ) -> dict[str, int | bool | str | None]:
        """Return order and optional mode-count diagnostics."""

        info: dict[str, int | bool | str | None] = {
            "bulk_degree": self.bulk_degree,
            "bridge_trace_degree": self.bridge_trace_degree,
            "surface_current_degree": self.surface_current_degree,
            "surface_omega_degree": self.surface_omega_degree,
            "face_family": self.face_family,
            "bridge_trace_dim": self.bridge_trace_dim,
            "required_parent_order": self.required_parent_order,
            "parent_order": None,
            "parent_order_admissible": None,
            "parent_order_excess": None,
        }
        if parent_order is not None:
            p = _nonnegative_int(parent_order, "parent_order")
            info["parent_order"] = p
            info["parent_order_admissible"] = p >= self.required_parent_order
            info["parent_order_excess"] = p - self.required_parent_order
        if evrs_rank is not None and cycle_rank is not None and surface_modes is not None:
            info["evrs_rank"] = _nonnegative_int(evrs_rank, "evrs_rank")
            info["cycle_rank"] = _nonnegative_int(cycle_rank, "cycle_rank")
            info["surface_modes"] = _nonnegative_int(surface_modes, "surface_modes")
            info["non_sibc_trace_modes"] = _nonnegative_int(
                non_sibc_trace_modes,
                "non_sibc_trace_modes",
            )
            info["bridge_modes"] = self.bridge_modes(info["cycle_rank"])
            info["estimated_reduced_modes"] = self.estimated_reduced_modes(
                evrs_rank=info["evrs_rank"],
                cycle_rank=info["cycle_rank"],
                surface_modes=info["surface_modes"],
                non_sibc_trace_modes=info["non_sibc_trace_modes"],
            )
        return info


def _radia_cpp_kernel(name: str):
    try:
        from radia import _radia_pybind as _radia_cpp
    except Exception:
        return None
    return getattr(_radia_cpp, name, None)


def _as_points(points, name: str) -> np.ndarray:
    arr = np.asarray(points, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n, 3)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def _as_weights(weights, n: int, name: str) -> np.ndarray:
    arr = np.asarray(weights, dtype=float)
    if arr.ndim != 1 or arr.shape[0] != n:
        raise ValueError(f"{name} must have shape ({n},)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be positive")
    return arr


def _as_modes(modes, n: int, name: str) -> np.ndarray:
    arr = np.asarray(modes)
    if arr.ndim == 2 and arr.shape == (n, 3):
        arr = arr[np.newaxis, :, :]
    if arr.ndim != 3 or arr.shape[1:] != (n, 3):
        raise ValueError(f"{name} must have shape (m, {n}, 3)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def _mode_names(names: Iterable[str] | None, count: int) -> tuple[str, ...]:
    if names is None:
        return tuple(f"mode{i}" for i in range(count))
    out = tuple(str(name) for name in names)
    if len(out) != count:
        raise ValueError(f"names must contain {count} entries")
    return out


@dataclass(frozen=True)
class SampledCurrentBasis:
    """Current basis sampled on a volume or surface quadrature rule.

    Parameters
    ----------
    points:
        Quadrature points, shape ``(n, 3)``.
    weights:
        Positive volume or surface weights, shape ``(n,)``.
    modes:
        Current-density samples, shape ``(m, n, 3)``.  For ``kind="volume"``
        the modes are volume current densities J [A/m^2].  For
        ``kind="surface"`` the modes are surface currents K [A/m].
    kind:
        Either ``"volume"`` for T-method bulk currents or ``"surface"`` for
        surface-Omega/SIBC currents.
    """

    points: np.ndarray
    weights: np.ndarray
    modes: np.ndarray
    kind: str
    names: tuple[str, ...]

    def __post_init__(self) -> None:
        points = _as_points(self.points, "points")
        weights = _as_weights(self.weights, points.shape[0], "weights")
        modes = _as_modes(self.modes, points.shape[0], "modes")
        if self.kind not in {"volume", "surface"}:
            raise ValueError("kind must be 'volume' or 'surface'")
        names = _mode_names(self.names, modes.shape[0])
        object.__setattr__(self, "points", points)
        object.__setattr__(self, "weights", weights)
        object.__setattr__(self, "modes", modes)
        object.__setattr__(self, "names", names)

    @property
    def n_modes(self) -> int:
        return int(self.modes.shape[0])

    @property
    def n_samples(self) -> int:
        return int(self.points.shape[0])

    def mass_matrix(self) -> np.ndarray:
        """Return ``int mode_i dot mode_j dV`` or ``dS`` on this quadrature."""

        return np.einsum(
            "aik,bik,i->ab", self.modes.conj(), self.modes, self.weights
        )


@dataclass(frozen=True)
class SampledMagnetizationBasis:
    """Magnetization basis sampled on a volume quadrature rule.

    This is the HDiv-VIM coexistence hook.  The modes represent magnetization
    ``M`` [A/m], typically sampled from an NGSolve ``HDiv`` GridFunction or from
    an already reduced Radia VIM magnetic basis.  Coupling to the eddy branch is
    formed through ``int M_i dot B[J_j] dV``.
    """

    points: np.ndarray
    weights: np.ndarray
    modes: np.ndarray
    names: tuple[str, ...]

    def __post_init__(self) -> None:
        points = _as_points(self.points, "points")
        weights = _as_weights(self.weights, points.shape[0], "weights")
        modes = _as_modes(self.modes, points.shape[0], "modes")
        names = _mode_names(self.names, modes.shape[0])
        object.__setattr__(self, "points", points)
        object.__setattr__(self, "weights", weights)
        object.__setattr__(self, "modes", modes)
        object.__setattr__(self, "names", names)

    @property
    def n_modes(self) -> int:
        return int(self.modes.shape[0])

    @property
    def n_samples(self) -> int:
        return int(self.points.shape[0])

    def mass_matrix(self) -> np.ndarray:
        """Return ``int M_i dot M_j dV`` on this quadrature."""

        return np.einsum(
            "aik,bik,i->ab", self.modes.conj(), self.modes, self.weights
        )


def VolumeCurrentBasis(points, weights, current_modes, names=None) -> SampledCurrentBasis:
    """Create sampled T-method bulk-current modes."""

    points_arr = _as_points(points, "points")
    modes = _as_modes(current_modes, points_arr.shape[0], "current_modes")
    return SampledCurrentBasis(
        points=points_arr,
        weights=_as_weights(weights, points_arr.shape[0], "weights"),
        modes=modes,
        kind="volume",
        names=_mode_names(names, modes.shape[0]),
    )


def MagnetizationBasis(points, weights, magnetization_modes, names=None) -> SampledMagnetizationBasis:
    """Create sampled HDiv-compatible magnetization modes."""

    points_arr = _as_points(points, "points")
    modes = _as_modes(magnetization_modes, points_arr.shape[0], "magnetization_modes")
    return SampledMagnetizationBasis(
        points=points_arr,
        weights=_as_weights(weights, points_arr.shape[0], "weights"),
        modes=modes,
        names=_mode_names(names, modes.shape[0]),
    )


def SurfaceOmegaBasis(points, weights, normals, grad_omega_modes, names=None) -> SampledCurrentBasis:
    """Create surface-Omega/SIBC current modes ``K = n x grad_Gamma Omega``.

    ``grad_omega_modes`` may include a normal component; it is projected onto
    the tangent plane before the cross product is formed.
    """

    points_arr = _as_points(points, "points")
    normals_arr = _as_points(normals, "normals")
    if normals_arr.shape[0] != points_arr.shape[0]:
        raise ValueError("normals must have the same sample count as points")
    normal_norm = np.linalg.norm(normals_arr, axis=1)
    if np.any(normal_norm <= 0.0):
        raise ValueError("normals must be non-zero")
    normals_arr = normals_arr / normal_norm[:, np.newaxis]
    grads = _as_modes(grad_omega_modes, points_arr.shape[0], "grad_omega_modes")
    normal_part = np.einsum("aik,ik->ai", grads, normals_arr)
    tangent_grads = grads - normal_part[:, :, np.newaxis] * normals_arr[np.newaxis, :, :]
    currents = np.cross(normals_arr[np.newaxis, :, :], tangent_grads)
    return SampledCurrentBasis(
        points=points_arr,
        weights=_as_weights(weights, points_arr.shape[0], "weights"),
        modes=currents,
        kind="surface",
        names=_mode_names(names, currents.shape[0]),
    )


def _ngsolve_vb(vb):
    import ngsolve as ng

    if vb in (ng.VOL, "VOL", "vol", None):
        return ng.VOL
    if vb in (ng.BND, "BND", "bnd"):
        return ng.BND
    raise ValueError("vb must be 'VOL' or 'BND'")


def _label_filter(labels):
    if labels is None:
        return None
    if isinstance(labels, str):
        return {labels}
    return {str(label) for label in labels}


def _vector_cf_value(value, name: str) -> np.ndarray:
    arr = np.asarray(value).reshape(-1)
    if arr.shape == (2,):
        arr = np.array((arr[0], arr[1], 0.0), dtype=arr.dtype)
    if arr.shape != (3,):
        raise ValueError(f"{name} must evaluate to a 2- or 3-vector")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} evaluated to a non-finite value")
    return arr


def _physical_point3(value, name: str) -> np.ndarray:
    # NGSolve may expose mapped coordinates through a reused temporary buffer.
    # Own the values before retaining them across quadrature iterations.
    point = np.array(value, dtype=float, copy=True).reshape(-1)
    if point.shape == (2,):
        point = np.array((point[0], point[1], 0.0), dtype=float)
    if point.shape != (3,):
        raise ValueError(f"{name} must evaluate to a 2-D or 3-D point")
    if not np.all(np.isfinite(point)):
        raise ValueError(f"{name} evaluated to a non-finite point")
    return point


_HCURL_CELL_FAMILY_NAMES = {
    "ET.TET": "tet",
    "ET.HEX": "hex",
    "ET.PRISM": "wedge",
    "ET.PYRAMID": "pyramid",
    "ET.TRIG": "trig",
    "ET.QUAD": "quad",
}


@dataclass(frozen=True)
class HCurlCellFamilyInventory:
    """NGSolve cell families carried by one HCurl Eddy Bubble parent mesh."""

    dimension: int
    counts: dict[str, int]

    def __post_init__(self) -> None:
        dimension = int(self.dimension)
        if dimension not in (2, 3):
            raise ValueError("HCurl Eddy Bubble supports mesh dimension 2 or 3")
        counts = {str(name): int(count) for name, count in self.counts.items()}
        if not counts or any(count <= 0 for count in counts.values()):
            raise ValueError("cell-family counts must be positive")
        allowed = {"trig", "quad"} if dimension == 2 else {
            "tet", "hex", "pyramid", "wedge"
        }
        unsupported = set(counts) - allowed
        if unsupported:
            raise ValueError(f"unsupported HCurl cell families: {sorted(unsupported)}")
        object.__setattr__(self, "dimension", dimension)
        object.__setattr__(self, "counts", counts)

    @property
    def families(self) -> tuple[str, ...]:
        return tuple(sorted(self.counts))

    @property
    def cell_count(self) -> int:
        return int(sum(self.counts.values()))

    @property
    def mixed(self) -> bool:
        return len(self.counts) > 1

    def diagnostics(self) -> dict[str, object]:
        return {
            "dimension": self.dimension,
            "families": self.families,
            "family_counts": dict(self.counts),
            "cell_count": self.cell_count,
            "mixed_cell_families": self.mixed,
            "eddy_bubble_supported": True,
            "current_representation": (
                "three-component-vector"
                if self.dimension == 3
                else "out-of-plane-Jz"
            ),
            "vim_interaction_family": (
                "analytic-tet-decomposition"
                if self.dimension == 3
                else "planar-log-kernel"
            ),
        }


def NgsolveHCurlCellFamilies(mesh, *, materials=None) -> HCurlCellFamilyInventory:
    """Inventory supported T/HCurl parent cells in an NGSolve mesh."""

    import ngsolve as ng

    dimension = int(mesh.dim)
    selected = _label_filter(materials)
    counts: dict[str, int] = {}
    for element in mesh.Elements(ng.VOL):
        if selected is not None and str(element.mat) not in selected:
            continue
        element_type = str(element.type)
        family = _HCURL_CELL_FAMILY_NAMES.get(element_type)
        if family is None:
            raise NotImplementedError(
                f"HCurl Eddy Bubble does not support NGSolve cell {element_type}"
            )
        counts[family] = counts.get(family, 0) + 1
    if not counts:
        raise ValueError("the selected conductor region contains no volume cells")
    return HCurlCellFamilyInventory(dimension=dimension, counts=counts)


@dataclass(frozen=True)
class EddyFaceTopology:
    """Topology classification of one mesh face for eddy-current reduction.

    The classification separates air/exterior conductor faces, where SIBC is
    allowed, from conductor-conductor faces, where loop currents must be
    allowed to cross the face.  Conductor-insulator faces may still need a
    boundary trace basis, but they are not SIBC half-space faces.  This
    prevents a local eddy-bubble elimination from accidentally cutting an
    eddy-current loop.
    """

    face_nr: int
    role: str
    volume_elements: tuple[int, ...]
    volume_materials: tuple[str, ...]
    boundary_elements: tuple[int, ...] = ()
    boundary_labels: tuple[str, ...] = ()

    @property
    def has_conductor(self) -> bool:
        return self.role != "nonconductive"

    @property
    def requires_surface_basis(self) -> bool:
        """Return whether this conductor face may need a boundary trace basis.

        Only ``is_sibc_face`` may receive SIBC termination.  Insulating
        neighbors are trace-boundary faces, not air half-spaces.
        """

        return self.role in {
            "conductor-exterior",
            "conductor-air",
            "conductor-insulator",
        }

    @property
    def is_sibc_face(self) -> bool:
        return self.role in {
            "conductor-exterior",
            "conductor-air",
        }

    @property
    def requires_loop_bridge(self) -> bool:
        return self.role in {
            "conductor-conductor",
            "conductive-interface",
        }

    @property
    def can_sibc_terminate(self) -> bool:
        return self.is_sibc_face


@dataclass(frozen=True)
class EddyConductorGraphEdge:
    """One conductor-conductor face as an edge in the conductive graph."""

    face_nr: int
    left_element: int
    right_element: int
    left_material: str
    right_material: str

    @property
    def endpoints(self) -> tuple[int, int]:
        return (self.left_element, self.right_element)

    @property
    def is_material_interface(self) -> bool:
        return self.left_material != self.right_material


@dataclass(frozen=True)
class EddyConductorCycle:
    """A fundamental cycle in the conductor adjacency graph."""

    edge_indices: tuple[int, ...]
    face_nrs: tuple[int, ...]
    elements: tuple[int, ...]
    edge_signs: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        if self.edge_signs:
            if len(self.edge_signs) != len(self.edge_indices):
                raise ValueError("edge_signs must match edge_indices length")
            if any(sign not in (-1, 1) for sign in self.edge_signs):
                raise ValueError("edge_signs must contain only -1 or 1")


@dataclass(frozen=True)
class EddyConductorGraph:
    """Conductor adjacency graph used to reduce loop-bridge DoFs."""

    nodes: tuple[int, ...]
    edges: tuple[EddyConductorGraphEdge, ...]

    @property
    def node_count(self) -> int:
        return len(self.nodes)

    @property
    def edge_count(self) -> int:
        return len(self.edges)

    def _component_roots(self) -> set[int]:
        parent = {node: node for node in self.nodes}

        def find(node: int) -> int:
            while parent[node] != node:
                parent[node] = parent[parent[node]]
                node = parent[node]
            return node

        def union(left: int, right: int) -> None:
            root_left = find(left)
            root_right = find(right)
            if root_left != root_right:
                parent[root_right] = root_left

        for edge in self.edges:
            union(edge.left_element, edge.right_element)
        return {find(node) for node in self.nodes}

    @property
    def component_count(self) -> int:
        return len(self._component_roots())

    @property
    def cycle_rank(self) -> int:
        return max(0, self.edge_count - self.node_count + self.component_count)

    def cycle_basis(self) -> tuple[EddyConductorCycle, ...]:
        """Return a fundamental cycle basis from a spanning forest.

        The cycles are graph-topological.  They are the production target for
        replacing a large set of conductor-conductor bridge DoFs with a much
        smaller loop basis that still allows eddy currents to circulate.
        """

        node_set = set(self.nodes)
        parent = {node: node for node in self.nodes}

        def find(node: int) -> int:
            while parent[node] != node:
                parent[node] = parent[parent[node]]
                node = parent[node]
            return node

        def union(left: int, right: int) -> bool:
            root_left = find(left)
            root_right = find(right)
            if root_left == root_right:
                return False
            parent[root_right] = root_left
            return True

        tree_adjacency: dict[int, list[tuple[int, int]]] = {node: [] for node in node_set}
        cycles: list[EddyConductorCycle] = []

        def tree_path(start: int, goal: int) -> tuple[list[int], list[int], list[int]]:
            queue = [start]
            seen = {start}
            prev: dict[int, tuple[int, int]] = {}
            for node in queue:
                if node == goal:
                    break
                for nxt, edge_idx in tree_adjacency[node]:
                    if nxt in seen:
                        continue
                    seen.add(nxt)
                    prev[nxt] = (node, edge_idx)
                    queue.append(nxt)
            if goal not in seen:
                raise RuntimeError("spanning-forest path was not found")
            elements = [goal]
            edge_indices: list[int] = []
            edge_signs: list[int] = []
            node = goal
            while node != start:
                node_prev, edge_idx = prev[node]
                edge = self.edges[edge_idx]
                edge_indices.append(edge_idx)
                edge_signs.append(1 if edge.endpoints == (node_prev, node) else -1)
                elements.append(node_prev)
                node = node_prev
            elements.reverse()
            edge_indices.reverse()
            edge_signs.reverse()
            return elements, edge_indices, edge_signs

        for edge_idx, edge in enumerate(self.edges):
            left, right = edge.endpoints
            if union(left, right):
                tree_adjacency[left].append((right, edge_idx))
                tree_adjacency[right].append((left, edge_idx))
                continue
            elements, path_edges, path_signs = tree_path(left, right)
            cycle_edges = tuple(path_edges + [edge_idx])
            cycle_signs = tuple(path_signs + [-1])
            cycles.append(
                EddyConductorCycle(
                    edge_indices=cycle_edges,
                    face_nrs=tuple(self.edges[i].face_nr for i in cycle_edges),
                    elements=tuple(elements),
                    edge_signs=cycle_signs,
                )
            )
        return tuple(cycles)

    def cycle_edge_matrix(self, *, dtype=float) -> np.ndarray:
        """Return signed cycle-edge incidence, shape ``(n_cycles, n_edges)``."""

        cycles = self.cycle_basis()
        matrix = np.zeros((len(cycles), self.edge_count), dtype=dtype)
        for i, cycle in enumerate(cycles):
            for edge_idx, sign in zip(cycle.edge_indices, cycle.edge_signs):
                matrix[i, edge_idx] = sign
        return matrix

    def diagnostics(self) -> dict[str, int]:
        cycles = self.cycle_basis()
        return {
            "node_count": self.node_count,
            "edge_count": self.edge_count,
            "component_count": self.component_count,
            "cycle_rank": self.cycle_rank,
            "fundamental_cycle_count": len(cycles),
            "material_interface_edge_count": sum(
                1 for edge in self.edges if edge.is_material_interface
            ),
        }


@dataclass(frozen=True)
class EddyMeshTopology:
    """Conductive adjacency summary used by topology-aware EVRS pruning."""

    faces: tuple[EddyFaceTopology, ...]
    conductive_materials: tuple[str, ...]
    air_materials: tuple[str, ...]

    def faces_by_role(self, role: str) -> tuple[EddyFaceTopology, ...]:
        """Return all faces with the requested role."""

        return tuple(face for face in self.faces if face.role == role)

    @property
    def surface_faces(self) -> tuple[EddyFaceTopology, ...]:
        """Faces where surface-Omega/SIBC basis functions may be required."""

        return tuple(face for face in self.faces if face.requires_surface_basis)

    @property
    def sibc_faces(self) -> tuple[EddyFaceTopology, ...]:
        """Air/exterior conductor faces where SIBC termination is allowed."""

        return tuple(face for face in self.faces if face.is_sibc_face)

    @property
    def non_sibc_trace_faces(self) -> tuple[EddyFaceTopology, ...]:
        """Boundary-trace conductor faces that are not SIBC half-space faces."""

        return tuple(
            face for face in self.surface_faces if not face.is_sibc_face
        )

    @property
    def loop_bridge_faces(self) -> tuple[EddyFaceTopology, ...]:
        """Faces where conductor-conductor loop closure must be preserved."""

        return tuple(face for face in self.faces if face.requires_loop_bridge)

    def conductor_graph(self) -> EddyConductorGraph:
        """Return the conductive adjacency graph induced by loop-bridge faces."""

        nodes = {
            element
            for face in self.faces
            for element, material in zip(face.volume_elements, face.volume_materials)
            if material in self.conductive_materials
        }
        edges: list[EddyConductorGraphEdge] = []
        for face in self.loop_bridge_faces:
            conductive = [
                (element, material)
                for element, material in zip(face.volume_elements, face.volume_materials)
                if material in self.conductive_materials
            ]
            for i, (left, left_material) in enumerate(conductive):
                for right, right_material in conductive[i + 1 :]:
                    if left == right:
                        continue
                    if right < left:
                        left, right = right, left
                        left_material, right_material = right_material, left_material
                    edges.append(
                        EddyConductorGraphEdge(
                            face_nr=face.face_nr,
                            left_element=left,
                            right_element=right,
                            left_material=left_material,
                            right_material=right_material,
                        )
                    )
        edges.sort(key=lambda edge: (edge.left_element, edge.right_element, edge.face_nr))
        return EddyConductorGraph(nodes=tuple(sorted(nodes)), edges=tuple(edges))

    def diagnostics(self) -> dict[str, int | tuple[str, ...]]:
        """Return compact topology counts for validation JSON."""

        roles = sorted({face.role for face in self.faces})
        graph = self.conductor_graph()
        return {
            "face_count": len(self.faces),
            "roles": tuple(roles),
            "conductive_materials": self.conductive_materials,
            "air_materials": self.air_materials,
            "conductive_element_count": graph.node_count,
            "conductive_component_count": graph.component_count,
            "conductive_graph_edge_count": graph.edge_count,
            "conductive_graph_cycle_rank": graph.cycle_rank,
            "surface_face_count": len(self.surface_faces),
            "sibc_face_count": len(self.sibc_faces),
            "non_sibc_trace_face_count": len(self.non_sibc_trace_faces),
            "loop_bridge_face_count": len(self.loop_bridge_faces),
            "conductor_exterior_face_count": len(self.faces_by_role("conductor-exterior")),
            "conductor_air_face_count": len(self.faces_by_role("conductor-air")),
            "conductor_insulator_face_count": len(self.faces_by_role("conductor-insulator")),
            "conductor_conductor_face_count": len(self.faces_by_role("conductor-conductor")),
            "conductive_interface_face_count": len(self.faces_by_role("conductive-interface")),
        }


def _node_nr(node) -> int:
    nr = getattr(node, "nr", None)
    if nr is None:
        return int(node)
    return int(nr)


def _element_nr(element) -> int:
    nr = getattr(element, "nr", None)
    if nr is None:
        return int(getattr(element, "index"))
    return int(nr)


def _material_name(element) -> str:
    return str(getattr(element, "mat"))


def _classify_eddy_face(
    volume_records: tuple[tuple[int, str], ...],
    boundary_records: tuple[tuple[int, str], ...],
    conductive_materials: set[str],
    air_materials: set[str],
) -> str:
    conductive = [record for record in volume_records if record[1] in conductive_materials]
    if not conductive:
        return "nonconductive"
    if len(conductive) >= 2:
        materials = {material for _, material in conductive}
        return "conductor-conductor" if len(materials) == 1 else "conductive-interface"

    nonconductive_neighbors = [
        material for _, material in volume_records if material not in conductive_materials
    ]
    if any(material in air_materials for material in nonconductive_neighbors):
        return "conductor-air"
    if nonconductive_neighbors:
        return "conductor-insulator"
    return "conductor-exterior"


def ClassifyNgsolveEddyTopology(
    mesh,
    conductive_materials,
    *,
    air_materials=("air", "vacuum"),
) -> EddyMeshTopology:
    """Classify conductive face adjacency in an NGSolve mesh.

    The resulting topology is intentionally independent of a particular basis.
    It tells the reducer where surface/SIBC modes are allowed and where
    conductor-conductor loop bridges must be kept visible during Schur or EVRS
    elimination.
    """

    import ngsolve as ng

    conductive = _label_filter(conductive_materials)
    if not conductive:
        raise ValueError("conductive_materials must not be empty")
    air = _label_filter(air_materials) or set()

    dimension = int(mesh.dim)
    if dimension not in (2, 3):
        raise ValueError("eddy topology requires a 2-D or 3-D mesh")

    def codimension_one_entities(element):
        return element.faces if dimension == 3 else element.edges

    face_to_volume: dict[int, list[tuple[int, str]]] = {}
    for element in mesh.Elements(ng.VOL):
        record = (_element_nr(element), _material_name(element))
        for face in codimension_one_entities(element):
            face_to_volume.setdefault(_node_nr(face), []).append(record)

    face_to_boundary: dict[int, list[tuple[int, str]]] = {}
    for element in mesh.Elements(ng.BND):
        record = (_element_nr(element), _material_name(element))
        for face in codimension_one_entities(element):
            face_to_boundary.setdefault(_node_nr(face), []).append(record)

    faces: list[EddyFaceTopology] = []
    for face_nr in sorted(set(face_to_volume) | set(face_to_boundary)):
        volume_records = tuple(face_to_volume.get(face_nr, ()))
        boundary_records = tuple(face_to_boundary.get(face_nr, ()))
        role = _classify_eddy_face(
            volume_records,
            boundary_records,
            conductive,
            air,
        )
        faces.append(
            EddyFaceTopology(
                face_nr=face_nr,
                role=role,
                volume_elements=tuple(record[0] for record in volume_records),
                volume_materials=tuple(record[1] for record in volume_records),
                boundary_elements=tuple(record[0] for record in boundary_records),
                boundary_labels=tuple(record[1] for record in boundary_records),
            )
        )

    return EddyMeshTopology(
        faces=tuple(faces),
        conductive_materials=tuple(sorted(conductive)),
        air_materials=tuple(sorted(air)),
    )


def _mark_dofs(mask: np.ndarray, dofs) -> None:
    for dof in dofs:
        idx = int(dof)
        if idx >= 0 and idx < mask.shape[0]:
            mask[idx] = True


@dataclass(frozen=True)
class EddyDofPolicy:
    """Topology-aware HCurl DoF masks for eddy-current reduction.

    The masks are not a final basis by themselves.  They are the policy layer
    that tells the EVRS reducer which DoFs are surface/SIBC candidates, which
    DoFs protect conductor-conductor loop closure, and which DoFs can be
    aggressively response-compressed as ordinary interior eddy bubbles.
    """

    free: np.ndarray
    sibc_surface: np.ndarray
    surface_candidate: np.ndarray
    loop_bridge: np.ndarray
    local_bubble: np.ndarray
    interface: np.ndarray
    wirebasket: np.ndarray

    def __post_init__(self) -> None:
        masks = {
            "free": self.free,
            "sibc_surface": self.sibc_surface,
            "surface_candidate": self.surface_candidate,
            "loop_bridge": self.loop_bridge,
            "local_bubble": self.local_bubble,
            "interface": self.interface,
            "wirebasket": self.wirebasket,
        }
        n = None
        for name, value in masks.items():
            arr = np.asarray(value, dtype=bool)
            if arr.ndim != 1:
                raise ValueError(f"{name} mask must be 1-dimensional")
            if n is None:
                n = arr.shape[0]
            elif arr.shape[0] != n:
                raise ValueError("all masks must have the same length")
            object.__setattr__(self, name, arr)

    @property
    def ndof(self) -> int:
        return int(self.free.shape[0])

    @property
    def topology_protected(self) -> np.ndarray:
        """Free DoFs that must preserve conductor-conductor loop closure."""

        return self.free & self.loop_bridge

    @property
    def sibc_surface_only(self) -> np.ndarray:
        """Free SIBC-surface DoFs that are not also loop-bridge DoFs."""

        return self.free & self.sibc_surface & ~self.loop_bridge

    @property
    def non_sibc_trace(self) -> np.ndarray:
        """Free boundary-trace DoFs that are not SIBC and not loop bridges."""

        return self.free & self.surface_candidate & ~self.sibc_surface & ~self.loop_bridge

    @property
    def exterior_surface_only(self) -> np.ndarray:
        """Free exterior surface DoFs excluding loop-bridge DoFs."""

        return self.free & self.surface_candidate & ~self.loop_bridge

    @property
    def ordinary_evrs_candidate(self) -> np.ndarray:
        """Free bulk DoFs that are not boundary traces or loop bridges."""

        return self.free & ~self.surface_candidate & ~self.loop_bridge

    @property
    def interior_evrs_candidate(self) -> np.ndarray:
        """Alias for ordinary EVRS candidates used by validation JSON."""

        return self.ordinary_evrs_candidate

    def diagnostics(self) -> dict[str, int | float]:
        """Return count diagnostics for topology-aware DoF reduction."""

        free_count = int(np.count_nonzero(self.free))
        sibc_count = int(np.count_nonzero(self.free & self.sibc_surface))
        surface_count = int(np.count_nonzero(self.free & self.surface_candidate))
        bridge_count = int(np.count_nonzero(self.topology_protected))
        overlap_count = int(np.count_nonzero(self.free & self.sibc_surface & self.loop_bridge))
        surface_only_count = int(np.count_nonzero(self.sibc_surface_only))
        non_sibc_trace_count = int(np.count_nonzero(self.non_sibc_trace))
        ordinary_count = int(np.count_nonzero(self.ordinary_evrs_candidate))
        partition_count = bridge_count + surface_only_count + non_sibc_trace_count + ordinary_count
        return {
            "ndof": self.ndof,
            "free_dofs": free_count,
            "sibc_surface_dofs": sibc_count,
            "surface_candidate_dofs": surface_count,
            "non_sibc_trace_dofs": non_sibc_trace_count,
            "loop_bridge_dofs": bridge_count,
            "sibc_loop_bridge_overlap_dofs": overlap_count,
            "sibc_surface_only_dofs": surface_only_count,
            "ordinary_evrs_candidate_dofs": ordinary_count,
            "partitioned_free_dofs": partition_count,
            "local_bubble_dofs": int(np.count_nonzero(self.free & self.local_bubble)),
            "interface_dofs": int(np.count_nonzero(self.free & self.interface)),
            "wirebasket_dofs": int(np.count_nonzero(self.free & self.wirebasket)),
            "loop_bridge_fraction": float(bridge_count / free_count) if free_count else 0.0,
            "sibc_surface_fraction": float(sibc_count / free_count) if free_count else 0.0,
            "non_sibc_trace_fraction": (
                float(non_sibc_trace_count / free_count) if free_count else 0.0
            ),
            "ordinary_evrs_fraction": float(ordinary_count / free_count) if free_count else 0.0,
        }

    def reduction_plan(
        self,
        *,
        evrs_rank: int | None = None,
        surface_modes: int | None = None,
        non_sibc_trace_modes: int | None = None,
        loop_bridge_modes: int | None = None,
        bridge_strategy: str | None = None,
    ) -> "EddyReductionPlan":
        """Return the conservative production reduction plan.

        The plan keeps conductor-conductor loop-bridge DoFs structurally
        protected, routes air-touching conductor-face DoFs to the SIBC/surface
        branch, keeps non-air boundary traces separate from SIBC, and leaves
        the remaining interior DoFs as ordinary EVRS candidates.
        """

        return EddyReductionPlan(
            free=self.free,
            loop_bridge_keep=self.topology_protected,
            sibc_surface_trace=self.free & self.sibc_surface,
            sibc_surface_only=self.sibc_surface_only,
            non_sibc_trace=self.non_sibc_trace,
            ordinary_evrs_candidate=self.ordinary_evrs_candidate,
            evrs_rank=evrs_rank,
            surface_modes=surface_modes,
            non_sibc_trace_modes=non_sibc_trace_modes,
            loop_bridge_modes=loop_bridge_modes,
            bridge_strategy=bridge_strategy,
        )


@dataclass(frozen=True)
class EddyReductionPlan:
    """Conservative production plan for topology-aware eddy DoF reduction."""

    free: np.ndarray
    loop_bridge_keep: np.ndarray
    sibc_surface_trace: np.ndarray
    sibc_surface_only: np.ndarray
    non_sibc_trace: np.ndarray
    ordinary_evrs_candidate: np.ndarray
    evrs_rank: int | None = None
    surface_modes: int | None = None
    non_sibc_trace_modes: int | None = None
    loop_bridge_modes: int | None = None
    bridge_strategy: str | None = None

    def __post_init__(self) -> None:
        masks = {
            "free": self.free,
            "loop_bridge_keep": self.loop_bridge_keep,
            "sibc_surface_trace": self.sibc_surface_trace,
            "sibc_surface_only": self.sibc_surface_only,
            "non_sibc_trace": self.non_sibc_trace,
            "ordinary_evrs_candidate": self.ordinary_evrs_candidate,
        }
        n = None
        for name, value in masks.items():
            arr = np.asarray(value, dtype=bool)
            if arr.ndim != 1:
                raise ValueError(f"{name} mask must be 1-dimensional")
            if n is None:
                n = arr.shape[0]
            elif arr.shape[0] != n:
                raise ValueError("all masks must have the same length")
            object.__setattr__(self, name, arr)
        if self.evrs_rank is not None and self.evrs_rank < 0:
            raise ValueError("evrs_rank must be non-negative")
        if self.surface_modes is not None and self.surface_modes < 0:
            raise ValueError("surface_modes must be non-negative")
        if self.non_sibc_trace_modes is not None and self.non_sibc_trace_modes < 0:
            raise ValueError("non_sibc_trace_modes must be non-negative")
        if self.loop_bridge_modes is not None and self.loop_bridge_modes < 0:
            raise ValueError("loop_bridge_modes must be non-negative")

    @property
    def ndof(self) -> int:
        return int(self.free.shape[0])

    @property
    def conservative_structural_keep(self) -> np.ndarray:
        """DoFs not allowed to disappear through local eddy-bubble pruning."""

        return self.free & (
            self.loop_bridge_keep | self.sibc_surface_only | self.non_sibc_trace
        )

    def diagnostics(self) -> dict[str, int | float | str | None]:
        """Return plan counts and optional reduced-mode estimate."""

        free_count = int(np.count_nonzero(self.free))
        bridge_count = int(np.count_nonzero(self.free & self.loop_bridge_keep))
        surface_trace_count = int(np.count_nonzero(self.free & self.sibc_surface_trace))
        surface_only_count = int(np.count_nonzero(self.free & self.sibc_surface_only))
        non_sibc_trace_count = int(np.count_nonzero(self.free & self.non_sibc_trace))
        ordinary_count = int(np.count_nonzero(self.free & self.ordinary_evrs_candidate))
        protected_count = int(np.count_nonzero(self.conservative_structural_keep))
        overlap_count = int(
            np.count_nonzero(self.free & self.loop_bridge_keep & self.sibc_surface_trace)
        )
        bridge_modes = (
            bridge_count if self.loop_bridge_modes is None else int(self.loop_bridge_modes)
        )
        non_sibc_trace_modes = (
            non_sibc_trace_count
            if self.non_sibc_trace_modes is None
            else int(self.non_sibc_trace_modes)
        )
        bridge_strategy = self.bridge_strategy
        if bridge_strategy is None:
            bridge_strategy = (
                "conservative-dof-keep"
                if self.loop_bridge_modes is None
                else "cycle-basis"
            )
        conservative_estimated = None
        estimated = None
        estimated_ratio = None
        if self.evrs_rank is not None and self.surface_modes is not None:
            conservative_estimated = (
                bridge_count
                + non_sibc_trace_count
                + int(self.evrs_rank)
                + int(self.surface_modes)
            )
            estimated = (
                bridge_modes
                + non_sibc_trace_modes
                + int(self.evrs_rank)
                + int(self.surface_modes)
            )
            estimated_ratio = float(estimated / free_count) if free_count else 0.0
        return {
            "rule": (
                "conductor-air/exterior-to-sibc; conductor-insulator-to-non-sibc-trace; "
                "conductor-conductor-to-loop-bridge; interior-to-evrs"
            ),
            "free_dofs": free_count,
            "loop_bridge_keep_dofs": bridge_count,
            "loop_bridge_reduced_modes": bridge_modes,
            "loop_bridge_reduction_strategy": bridge_strategy,
            "loop_bridge_reduction_savings": max(0, bridge_count - bridge_modes),
            "sibc_surface_trace_dofs": surface_trace_count,
            "sibc_surface_only_dofs": surface_only_count,
            "non_sibc_trace_dofs": non_sibc_trace_count,
            "non_sibc_trace_modes": non_sibc_trace_modes,
            "sibc_loop_bridge_overlap_dofs": overlap_count,
            "ordinary_evrs_candidate_dofs": ordinary_count,
            "conservative_structural_keep_dofs": protected_count,
            "evrs_rank": self.evrs_rank,
            "surface_modes": self.surface_modes,
            "conservative_estimated_reduced_modes": conservative_estimated,
            "estimated_reduced_modes": estimated,
            "estimated_reduction_ratio": estimated_ratio,
            "ordinary_evrs_candidate_fraction": (
                float(ordinary_count / free_count) if free_count else 0.0
            ),
            "conservative_structural_keep_fraction": (
                float(protected_count / free_count) if free_count else 0.0
            ),
        }


@dataclass(frozen=True)
class EddyBubbleDecomposition:
    """Named decomposition for topology-aware eddy bubbling.

    Eddy bubbling is the production reduction rule that separates HCurl parent
    DoFs into structural classes before any EVRS compression:

    * air/exterior conductor traces -> SIBC/surface-Omega branch,
    * conductor-insulator traces -> non-SIBC boundary trace branch,
    * conductor-conductor traces -> topology-protected bridge/cycle branch,
    * remaining bulk DoFs -> EVRS / eddy-bubble candidate branch.

    The object does not assemble a numerical basis by itself.  It records the
    class split and the conservative reduction plan that downstream VIM
    builders use.
    """

    policy: EddyDofPolicy
    plan: EddyReductionPlan
    topology: EddyMeshTopology | None = None
    conductor_graph: EddyConductorGraph | None = None
    parent_order: int | None = None
    parent_order_ledger: EddyParentOrderLedger | None = None

    @property
    def structural_keep(self) -> np.ndarray:
        """DoFs protected from ordinary eddy-bubble elimination."""

        return self.plan.conservative_structural_keep

    @property
    def eddy_bubble_candidate(self) -> np.ndarray:
        """Bulk DoFs eligible for EVRS/eddy-bubble compression."""

        return self.plan.ordinary_evrs_candidate

    @property
    def removable_candidate(self) -> np.ndarray:
        """Alias for the ordinary bulk eddy-bubble candidate mask."""

        return self.eddy_bubble_candidate

    def diagnostics(self) -> dict[str, object]:
        """Return class counts and reduction estimates for eddy bubbling."""

        plan_info = self.plan.diagnostics()
        policy_info = self.policy.diagnostics()
        free_count = int(plan_info["free_dofs"])
        keep_count = int(plan_info["conservative_structural_keep_dofs"])
        bubble_count = int(plan_info["ordinary_evrs_candidate_dofs"])
        info: dict[str, object] = {
            "rule": "topology-aware-eddy-bubbling",
            "free_dofs": free_count,
            "structural_keep_dofs": keep_count,
            "eddy_bubble_candidate_dofs": bubble_count,
            "eddy_bubble_candidate_fraction": (
                float(bubble_count / free_count) if free_count else 0.0
            ),
            "classes": {
                "sibc_surface": int(plan_info["sibc_surface_only_dofs"]),
                "non_sibc_trace": int(plan_info["non_sibc_trace_dofs"]),
                "loop_bridge": int(plan_info["loop_bridge_keep_dofs"]),
                "ordinary_bulk_eddy_bubble": bubble_count,
            },
            "policy": policy_info,
            "plan": plan_info,
        }
        if self.topology is not None:
            info["topology"] = self.topology.diagnostics()
        if self.conductor_graph is not None:
            info["conductor_graph"] = self.conductor_graph.diagnostics()
        if self.parent_order is not None:
            info["parent_order"] = _nonnegative_int(self.parent_order, "parent_order")
        if self.parent_order_ledger is not None:
            cycle_rank = None
            if self.conductor_graph is not None:
                cycle_rank = self.conductor_graph.cycle_rank
            info["parent_order_ledger"] = self.parent_order_ledger.diagnostics(
                parent_order=self.parent_order,
                evrs_rank=plan_info["evrs_rank"],
                cycle_rank=cycle_rank,
                surface_modes=plan_info["surface_modes"],
                non_sibc_trace_modes=plan_info["non_sibc_trace_modes"],
            )
        return info


def EddyBubbleReduction(
    policy: EddyDofPolicy,
    *,
    topology: EddyMeshTopology | None = None,
    conductor_graph: EddyConductorGraph | None = None,
    evrs_rank: int | None = None,
    surface_modes: int | None = None,
    non_sibc_trace_modes: int | None = None,
    loop_bridge_modes: int | None = None,
    bridge_strategy: str | None = None,
    parent_order: int | None = None,
    parent_order_ledger: EddyParentOrderLedger | None = None,
) -> EddyBubbleDecomposition:
    """Create an explicit eddy-bubbling decomposition from a DoF policy."""

    if conductor_graph is None and topology is not None:
        conductor_graph = topology.conductor_graph()
    plan = policy.reduction_plan(
        evrs_rank=evrs_rank,
        surface_modes=surface_modes,
        non_sibc_trace_modes=non_sibc_trace_modes,
        loop_bridge_modes=loop_bridge_modes,
        bridge_strategy=bridge_strategy,
    )
    return EddyBubbleDecomposition(
        policy=policy,
        plan=plan,
        topology=topology,
        conductor_graph=conductor_graph,
        parent_order=parent_order,
        parent_order_ledger=parent_order_ledger,
    )


def _ngsolve_fes_order(fes) -> int | None:
    """Best-effort extraction of an NGSolve finite-element parent order."""

    for name in ("globalorder", "order"):
        value = getattr(fes, name, None)
        if callable(value):
            try:
                value = value()
            except TypeError:
                continue
        if value is None:
            continue
        try:
            return int(value)
        except (TypeError, ValueError):
            continue
    for name in ("GetOrder",):
        func = getattr(fes, name, None)
        if not callable(func):
            continue
        try:
            return int(func())
        except (TypeError, ValueError):
            continue
    return None


def NgsolveEddyDofPolicy(
    mesh,
    fes,
    topology: EddyMeshTopology | None = None,
    *,
    conductive_materials=None,
    air_materials=("air", "vacuum"),
    free_dofs=None,
    include_edge_dofs: bool = True,
) -> EddyDofPolicy:
    """Build topology-aware HCurl DoF masks from face adjacency.

    Air/exterior conductor faces become SIBC-surface candidates.  Conductive
    interior faces become loop-bridge DoFs that should not be blindly removed
    by local high-order pruning, because they can close an eddy-current loop.
    """

    if topology is None:
        if conductive_materials is None:
            raise ValueError("conductive_materials is required when topology is not supplied")
        topology = ClassifyNgsolveEddyTopology(
            mesh,
            conductive_materials,
            air_materials=air_materials,
        )
    if int(fes.ndof) <= 0:
        raise ValueError("fes must have positive ndof")

    coupling = NgsolveCouplingDofMasks(fes, free_dofs=free_dofs)
    n = int(fes.ndof)
    sibc = np.zeros(n, dtype=bool)
    surface = np.zeros(n, dtype=bool)
    bridge = np.zeros(n, dtype=bool)

    dimension = int(mesh.dim)

    def mark_face(mask: np.ndarray, face_nr: int) -> None:
        entity = (
            mesh.faces[int(face_nr)]
            if dimension == 3
            else mesh.edges[int(face_nr)]
        )
        _mark_dofs(mask, fes.GetDofNrs(entity))
        if dimension == 3 and include_edge_dofs:
            for edge in entity.edges:
                _mark_dofs(mask, fes.GetDofNrs(edge))

    for face in topology.surface_faces:
        mark_face(surface, face.face_nr)
    for face in topology.sibc_faces:
        mark_face(sibc, face.face_nr)
    for face in topology.loop_bridge_faces:
        mark_face(bridge, face.face_nr)

    return EddyDofPolicy(
        free=coupling["free"],
        sibc_surface=sibc,
        surface_candidate=surface,
        loop_bridge=bridge,
        local_bubble=coupling["local_bubble"],
        interface=coupling["interface"],
        wirebasket=coupling["wirebasket"],
    )


def NgsolveEddyBubbleReduction(
    mesh,
    fes,
    topology: EddyMeshTopology | None = None,
    *,
    conductive_materials=None,
    air_materials=("air", "vacuum"),
    free_dofs=None,
    include_edge_dofs: bool = True,
    evrs_rank: int | None = None,
    surface_modes: int | None = None,
    non_sibc_trace_modes: int | None = None,
    loop_bridge_modes: int | None = None,
    bridge_strategy: str | None = None,
    parent_order: int | None = None,
    parent_order_ledger: EddyParentOrderLedger | None = None,
) -> EddyBubbleDecomposition:
    """Build the topology-aware eddy-bubbling decomposition for NGSolve HCurl."""

    if topology is None:
        if conductive_materials is None:
            raise ValueError("conductive_materials is required when topology is not supplied")
        topology = ClassifyNgsolveEddyTopology(
            mesh,
            conductive_materials,
            air_materials=air_materials,
        )
    policy = NgsolveEddyDofPolicy(
        mesh,
        fes,
        topology,
        free_dofs=free_dofs,
        include_edge_dofs=include_edge_dofs,
    )
    graph = topology.conductor_graph()
    return EddyBubbleReduction(
        policy,
        topology=topology,
        conductor_graph=graph,
        evrs_rank=evrs_rank,
        surface_modes=surface_modes,
        non_sibc_trace_modes=non_sibc_trace_modes,
        loop_bridge_modes=loop_bridge_modes,
        bridge_strategy=bridge_strategy,
        parent_order=(
            _ngsolve_fes_order(fes) if parent_order is None else parent_order
        ),
        parent_order_ledger=parent_order_ledger,
    )


def _mesh_node_point(mesh, node) -> np.ndarray:
    point = mesh.vertices[_node_nr(node)].point
    values = [float(value) for value in point]
    if len(values) == 2:
        values.append(0.0)
    if len(values) != 3:
        raise ValueError("mesh vertex must have two or three coordinates")
    return np.asarray(values, dtype=float)


def _mesh_face_points(mesh, face_nr: int) -> np.ndarray:
    face = mesh.faces[int(face_nr)]
    return np.vstack([_mesh_node_point(mesh, vertex) for vertex in face.vertices])


def _mesh_edge_center_length(mesh, edge_nr: int) -> tuple[np.ndarray, float]:
    edge = mesh.edges[int(edge_nr)]
    points = np.vstack([_mesh_node_point(mesh, vertex) for vertex in edge.vertices])
    if points.shape[0] != 2:
        raise ValueError("edge must have exactly two vertices")
    length = float(np.linalg.norm(points[1] - points[0]))
    if length <= 0.0:
        raise ValueError(f"edge {edge_nr} has zero length")
    return np.mean(points, axis=0), length


def _mesh_face_center_area(mesh, face_nr: int) -> tuple[np.ndarray, float]:
    points = _mesh_face_points(mesh, face_nr)
    center = np.mean(points, axis=0)
    if points.shape[0] < 3:
        raise ValueError("face must have at least three vertices")
    area_vector = np.zeros(3, dtype=float)
    origin = points[0]
    for i in range(1, points.shape[0] - 1):
        area_vector += 0.5 * np.cross(points[i] - origin, points[i + 1] - origin)
    area = float(np.linalg.norm(area_vector))
    if area <= 0.0:
        raise ValueError(f"face {face_nr} has zero area")
    return center, area


def _ngsolve_curve_order(mesh) -> int:
    getter = getattr(mesh, "GetCurveOrder", None)
    if getter is None:
        return 1
    try:
        return max(1, int(getter()))
    except (TypeError, ValueError, RuntimeError):
        return 1


def _ngsolve_geometry_intorder(mesh, intorder: int | None) -> int:
    if intorder is not None:
        if intorder < 0:
            raise ValueError("geometry_intorder must be non-negative")
        return int(intorder)
    return max(2, 2 * _ngsolve_curve_order(mesh) + 2)


def _ngsolve_element_centroids(mesh, intorder: int) -> dict[int, np.ndarray]:
    import ngsolve as ng

    centroids: dict[int, np.ndarray] = {}
    for element in mesh.Elements(ng.VOL):
        trafo = mesh.GetTrafo(element)
        volume = 0.0
        first_moment = np.zeros(3, dtype=float)
        for ip in ng.IntegrationRule(element.type, intorder):
            mip = trafo(ip)
            weight = float(ip.weight * mip.measure)
            volume += weight
            first_moment += weight * _physical_point3(
                mip.point,
                "element quadrature point",
            )
        if volume <= 0.0:
            raise ValueError(f"element {_element_nr(element)} has zero volume")
        centroids[_element_nr(element)] = first_moment / volume
    return centroids


def _ngsolve_tet_reference_vertices(mesh, element) -> dict[int, np.ndarray]:
    import ngsolve as ng

    reference = np.array(
        ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
        dtype=float,
    )
    vertices = tuple(element.vertices)
    if len(vertices) != 4:
        raise ValueError("curved bridge integration currently requires tetrahedra")
    physical = np.vstack([_mesh_node_point(mesh, vertex) for vertex in vertices])
    rule = ng.IntegrationRule(
        [tuple(float(value) for value in point) for point in reference],
        [1.0] * len(reference),
    )
    coordinate = ng.CF((ng.x, ng.y, ng.z))
    previous = None
    mapped = None
    for _ in range(4):
        current = np.asarray(
            coordinate(mesh.GetTrafo(element)(rule)),
            dtype=float,
        ).reshape(-1, 3)
        if (
            previous is not None
            and np.array_equal(previous, current)
            and np.all(np.isfinite(current))
        ):
            mapped = current
            break
        previous = current
    if mapped is None:
        raise RuntimeError("tetrahedron reference-vertex mapping remained unstable")

    result: dict[int, np.ndarray] = {}
    used: set[int] = set()
    scale = max(1.0, float(np.max(np.linalg.norm(physical, axis=1))))
    for local_index, point in enumerate(mapped):
        distances = np.linalg.norm(physical - point, axis=1)
        vertex_index = int(np.argmin(distances))
        if vertex_index in used or distances[vertex_index] > 1.0e-9 * scale:
            raise ValueError("could not match tetrahedron reference and mesh vertices")
        used.add(vertex_index)
        result[_node_nr(vertices[vertex_index])] = reference[local_index]
    return result


def _ngsolve_curved_tet_face_center_area(
    mesh,
    face_nr: int,
    element,
    intorder: int,
) -> tuple[np.ndarray, float]:
    import ngsolve as ng

    face_vertices = tuple(mesh.faces[int(face_nr)].vertices)
    if len(face_vertices) != 3:
        raise ValueError("curved tetrahedron face must have three vertices")
    reference_vertices = _ngsolve_tet_reference_vertices(mesh, element)
    try:
        r0, r1, r2 = (
            reference_vertices[_node_nr(vertex)] for vertex in face_vertices
        )
    except KeyError as exc:
        raise ValueError(f"face {face_nr} does not belong to the supplied element") from exc

    tangent0 = r1 - r0
    tangent1 = r2 - r0
    triangle_rule = ng.IntegrationRule(ng.ET.TRIG, intorder)
    reference_points = [
        tuple(
            float(value)
            for value in r0 + ip.point[0] * tangent0 + ip.point[1] * tangent1
        )
        for ip in triangle_rule
    ]
    mapped_rule = ng.IntegrationRule(
        reference_points,
        [float(ip.weight) for ip in triangle_rule],
    )
    trafo = mesh.GetTrafo(element)
    area = 0.0
    first_moment = np.zeros(3, dtype=float)
    for ip in mapped_rule:
        mip = trafo(ip)
        jacobian = np.asarray(mip.jacobi, dtype=float)
        surface_measure = float(
            np.linalg.norm(np.cross(jacobian @ tangent0, jacobian @ tangent1))
        )
        weight = float(ip.weight) * surface_measure
        area += weight
        first_moment += weight * np.asarray(mip.point, dtype=float)
    if area <= 0.0:
        raise ValueError(f"face {face_nr} has zero curved area")
    return first_moment / area, area


def _ngsolve_face_center_area(
    mesh,
    face_nr: int,
    element,
    intorder: int,
) -> tuple[np.ndarray, float]:
    if int(mesh.dim) == 2:
        return _mesh_edge_center_length(mesh, face_nr)
    if len(tuple(element.vertices)) == 4 and len(tuple(mesh.faces[int(face_nr)].vertices)) == 3:
        return _ngsolve_curved_tet_face_center_area(
            mesh,
            face_nr,
            element,
            intorder,
        )
    return _mesh_face_center_area(mesh, face_nr)


def NgsolveBridgeCycleCurrentBasis(
    mesh,
    topology: EddyMeshTopology | None = None,
    *,
    conductive_materials=None,
    air_materials=("air", "vacuum"),
    current_scale: float = 1.0,
    geometry_intorder: int | None = None,
    names=None,
) -> SampledCurrentBasis:
    """Build a coarse bridge-current basis from conductor graph cycles.

    Each conductor-conductor face is sampled at its area centroid.  A cycle mode
    places oriented current density along the segment between adjacent element
    volume centroids, with a dual-volume weight
    ``area(face) * centroid_distance``.  Curved tetrahedral maps are integrated
    with their pointwise Jacobian; ``geometry_intorder=None`` selects an order
    from ``mesh.GetCurveOrder()``.
    """

    if current_scale == 0.0 or not np.isfinite(current_scale):
        raise ValueError("current_scale must be finite and non-zero")
    if topology is None:
        if conductive_materials is None:
            raise ValueError("conductive_materials is required when topology is not supplied")
        topology = ClassifyNgsolveEddyTopology(
            mesh,
            conductive_materials,
            air_materials=air_materials,
        )

    graph = topology.conductor_graph()
    cycles = graph.cycle_basis()
    if graph.edge_count == 0:
        return VolumeCurrentBasis(
            np.zeros((0, 3)),
            np.zeros(0),
            np.zeros((len(cycles), 0, 3)),
            names=_mode_names(names, len(cycles)),
        )

    import ngsolve as ng

    geometry_intorder = _ngsolve_geometry_intorder(mesh, geometry_intorder)
    elements = {_element_nr(element): element for element in mesh.Elements(ng.VOL)}
    centroids = _ngsolve_element_centroids(mesh, geometry_intorder)
    points = np.zeros((graph.edge_count, 3), dtype=float)
    weights = np.zeros(graph.edge_count, dtype=float)
    directions = np.zeros((graph.edge_count, 3), dtype=float)
    for i, edge in enumerate(graph.edges):
        center, area = _ngsolve_face_center_area(
            mesh,
            edge.face_nr,
            elements[edge.left_element],
            geometry_intorder,
        )
        left_center = centroids[edge.left_element]
        right_center = centroids[edge.right_element]
        direction = right_center - left_center
        length = float(np.linalg.norm(direction))
        if length <= 0.0:
            raise ValueError("adjacent element centroids must be distinct")
        points[i] = center
        weights[i] = area * length
        directions[i] = direction / length

    modes = np.zeros((len(cycles), graph.edge_count, 3), dtype=float)
    for cycle_index, cycle in enumerate(cycles):
        for edge_idx, sign in zip(cycle.edge_indices, cycle.edge_signs):
            modes[cycle_index, edge_idx, :] = float(sign) * current_scale * directions[edge_idx]

    return VolumeCurrentBasis(
        points,
        weights,
        modes,
        names=_mode_names(names, len(cycles)),
    )


def SampleNgsolveVectorCFs(
    mesh,
    coefficient_functions,
    *,
    vb="VOL",
    intorder: int = 2,
    materials=None,
    boundaries=None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sample NGSolve vector coefficient functions on volume/boundary quadrature.

    Returns ``(points, weights, modes)`` where ``modes`` has shape
    ``(n_modes, n_samples, 3)``.  This is the bridge from high-order NGSolve
    fields such as ``curl(T)`` or ``grad(Omega).Trace()`` into the sampled
    reduced VIM basis.
    """

    import ngsolve as ng

    if intorder < 0:
        raise ValueError("intorder must be non-negative")
    vb_ng = _ngsolve_vb(vb)
    cfs = tuple(coefficient_functions)
    if not cfs:
        raise ValueError("coefficient_functions must not be empty")
    labels = _label_filter(materials if vb_ng == ng.VOL else boundaries)

    points: list[np.ndarray] = []
    weights: list[float] = []
    values: list[list[np.ndarray]] = [[] for _ in cfs]
    for el in mesh.Elements(vb_ng):
        label = str(el.mat)
        if labels is not None and label not in labels:
            continue
        ir = ng.IntegrationRule(el.type, intorder)
        trafo = mesh.GetTrafo(el)
        for ip in ir:
            mip = trafo(ip)
            points.append(_physical_point3(mip.point, "mapped quadrature point"))
            weights.append(float(ip.weight * mip.measure))
            for i, cf in enumerate(cfs):
                values[i].append(
                    np.array(
                        _vector_cf_value(cf(mip), f"coefficient_functions[{i}]"),
                        copy=True,
                    )
                )

    if not points:
        raise ValueError("sampling region produced no quadrature points")
    return (
        np.asarray(points, dtype=float),
        np.asarray(weights, dtype=float),
        np.asarray(values),
    )


def NgsolveVolumeCurrentBasis(
    mesh,
    current_modes,
    *,
    intorder: int = 2,
    materials=None,
    names=None,
) -> SampledCurrentBasis:
    """Sample NGSolve T-method volume-current modes into a VIM basis."""

    points, weights, modes = SampleNgsolveVectorCFs(
        mesh,
        current_modes,
        vb="VOL",
        intorder=intorder,
        materials=materials,
    )
    return VolumeCurrentBasis(points, weights, modes, names=names)


def NgsolveMagnetizationBasis(
    mesh,
    magnetization_modes,
    *,
    intorder: int = 2,
    materials=None,
    names=None,
) -> SampledMagnetizationBasis:
    """Sample NGSolve vector coefficient functions as magnetization modes."""

    points, weights, modes = SampleNgsolveVectorCFs(
        mesh,
        magnetization_modes,
        vb="VOL",
        intorder=intorder,
        materials=materials,
    )
    return MagnetizationBasis(points, weights, modes, names=names)


def NgsolveHDivMagnetizationBasis(
    mesh,
    fes,
    vectors,
    *,
    intorder: int = 2,
    materials=None,
    names=None,
) -> SampledMagnetizationBasis:
    """Convert HDiv coefficient vectors into sampled magnetization modes.

    ``vectors`` is shaped ``(fes.ndof, n_modes)``.  This is the light bridge
    from Radia/NGSolve HDiv magnetization coordinates into the eddy-current
    coupling helpers in this module.
    """

    import ngsolve as ng

    arr = np.asarray(vectors)
    if arr.ndim == 1:
        arr = arr[:, np.newaxis]
    if arr.ndim != 2 or arr.shape[0] != fes.ndof:
        raise ValueError(f"vectors must have shape ({fes.ndof}, n_modes)")
    gridfunctions = []
    modes = []
    for i in range(arr.shape[1]):
        gf = ng.GridFunction(fes)
        coeffs = gf.vec.FV().NumPy()
        coeffs[:] = arr[:, i]
        gridfunctions.append(gf)
        modes.append(gf)
    return NgsolveMagnetizationBasis(
        mesh,
        modes,
        intorder=intorder,
        materials=materials,
        names=names,
    )


def NgsolveHCurlCurlBasis(
    mesh,
    fes,
    vectors,
    *,
    intorder: int = 2,
    materials=None,
    names=None,
) -> SampledCurrentBasis:
    """Convert HCurl coefficient vectors into sampled ``curl(T)`` modes.

    ``vectors`` is shaped ``(fes.ndof, n_modes)``.  This is the direct handoff
    from :func:`BlockKrylovBasis` or an NGSolve eigensolve to the T-method
    volume-current basis used by the hybrid VIM.
    """

    import ngsolve as ng

    arr = np.asarray(vectors)
    if arr.ndim == 1:
        arr = arr[:, np.newaxis]
    if arr.ndim != 2 or arr.shape[0] != fes.ndof:
        raise ValueError(f"vectors must have shape ({fes.ndof}, n_modes)")
    gridfunctions = []
    curl_modes = []
    for i in range(arr.shape[1]):
        gf = ng.GridFunction(fes)
        coeffs = gf.vec.FV().NumPy()
        coeffs[:] = arr[:, i]
        gridfunctions.append(gf)
        current = ng.curl(gf)
        if int(mesh.dim) == 2:
            current = ng.CF((0.0, 0.0, current))
        curl_modes.append(current)
    return NgsolveVolumeCurrentBasis(
        mesh,
        curl_modes,
        intorder=intorder,
        materials=materials,
        names=names,
    )


def NgsolveSurfaceOmegaBasis(
    mesh,
    grad_omega_modes,
    *,
    intorder: int = 2,
    boundaries=None,
    normal_cf=None,
    names=None,
) -> SampledCurrentBasis:
    """Sample surface-Omega gradients into SIBC surface-current modes."""

    import ngsolve as ng

    grad_omega_modes = tuple(grad_omega_modes)
    if not grad_omega_modes:
        return SampledCurrentBasis(
            points=np.zeros((0, 3), dtype=float),
            weights=np.zeros(0, dtype=float),
            modes=np.zeros((0, 0, 3), dtype=float),
            kind="surface",
            names=_mode_names(names, 0),
        )
    points, weights, grads = SampleNgsolveVectorCFs(
        mesh,
        grad_omega_modes,
        vb="BND",
        intorder=intorder,
        boundaries=boundaries,
    )
    if normal_cf is None:
        normal_cf = ng.specialcf.normal(int(mesh.dim))
    normal_points, normal_weights, normals = SampleNgsolveVectorCFs(
        mesh,
        (normal_cf,),
        vb="BND",
        intorder=intorder,
        boundaries=boundaries,
    )
    if (
        normal_points.shape != points.shape
        or not np.allclose(normal_points, points)
        or not np.allclose(normal_weights, weights)
    ):
        raise RuntimeError("normal sampling did not match gradient sampling")
    return SurfaceOmegaBasis(points, weights, normals[0], grads, names=names)


def _call_or_value(obj, name: str):
    value = getattr(obj, name, None)
    return value() if callable(value) else value


def _matrix_size_from_coo(matrix, rows: np.ndarray, cols: np.ndarray) -> tuple[int, int]:
    height = _call_or_value(matrix, "height")
    width = _call_or_value(matrix, "width")
    if height is None:
        height = int(rows.max()) + 1 if rows.size else 0
    if width is None:
        width = int(cols.max()) + 1 if cols.size else 0
    return int(height), int(width)


def NgsolveMatrixToDense(matrix, *, dtype=None) -> np.ndarray:
    """Return an NGSolve matrix or bilinear form as a dense NumPy array.

    This is intentionally a small-system bridge for basis construction and
    validation.  Production large systems should use NGSolve's sparse solvers;
    the reduced VIM basis itself is small after the response compression.
    """

    mat = getattr(matrix, "mat", matrix)
    if hasattr(mat, "COO"):
        rows, cols, vals = mat.COO()
        rows = np.asarray(rows, dtype=int)
        cols = np.asarray(cols, dtype=int)
        vals = np.asarray(vals)
        if rows.shape != cols.shape or rows.shape != vals.shape:
            raise ValueError("NGSolve COO arrays must have matching shapes")
        height, width = _matrix_size_from_coo(mat, rows, cols)
        out_dtype = np.dtype(dtype) if dtype is not None else np.result_type(vals, float)
        dense = np.zeros((height, width), dtype=out_dtype)
        np.add.at(dense, (rows, cols), vals.astype(out_dtype, copy=False))
        return dense

    arr = np.asarray(matrix)
    if arr.ndim != 2:
        raise ValueError("matrix must be 2-dimensional or expose COO()")
    return np.array(arr, dtype=dtype, copy=True)


def _is_ngsolve_vector_like(obj) -> bool:
    return hasattr(obj, "vec") or hasattr(obj, "FV") or hasattr(obj, "NumPy")


def NgsolveVectorToArray(vector, *, dtype=None) -> np.ndarray:
    """Return an NGSolve vector, grid-function vector, or linear form as 1-D NumPy."""

    vec = getattr(vector, "vec", vector)
    if hasattr(vec, "FV"):
        arr = vec.FV().NumPy()
    elif hasattr(vec, "NumPy"):
        arr = vec.NumPy()
    else:
        arr = np.asarray(vec)
    arr = np.asarray(arr)
    if arr.ndim != 1:
        raise ValueError("vector must be 1-dimensional")
    return np.array(arr, dtype=dtype, copy=True)


def _ports_to_matrix(ports, *, dtype=None) -> np.ndarray:
    if _is_ngsolve_vector_like(ports):
        return NgsolveVectorToArray(ports, dtype=dtype)[:, np.newaxis]
    if isinstance(ports, (list, tuple)) and ports and any(
        _is_ngsolve_vector_like(port) for port in ports
    ):
        return np.column_stack(
            [NgsolveVectorToArray(port, dtype=dtype) for port in ports]
        )
    arr = np.asarray(ports)
    if dtype is not None:
        arr = arr.astype(dtype, copy=False)
    return arr


def NgsolveCouplingDofMasks(fes, free_dofs=None) -> dict[str, np.ndarray]:
    """Return NGSolve coupling-type masks for local static condensation.

    ``local_bubble`` corresponds to NGSolve ``LOCAL_DOF`` entries, i.e. the
    usual element-local high-order bubbles that can be removed by NGSolve
    static condensation.  This is distinct from Radia's response-compression
    terminology, where the removed high-order response components are called
    Eddy-Invisible DoFs (EIDs), or eddy bubbles for short.
    """

    n = int(fes.ndof)
    if free_dofs is None:
        free = np.asarray(fes.FreeDofs(False), dtype=bool)
    else:
        free = np.asarray(free_dofs, dtype=bool)
        if free.shape != (n,):
            raise ValueError(f"free_dofs boolean mask must have shape ({n},)")
    local = np.zeros(n, dtype=bool)
    interface = np.zeros(n, dtype=bool)
    wirebasket = np.zeros(n, dtype=bool)
    hidden = np.zeros(n, dtype=bool)
    for i in range(n):
        coupling = str(fes.CouplingType(i))
        if coupling.endswith("LOCAL_DOF"):
            local[i] = True
        elif coupling.endswith("INTERFACE_DOF"):
            interface[i] = True
        elif coupling.endswith("WIREBASKET_DOF"):
            wirebasket[i] = True
        else:
            hidden[i] = True
    local_bubble = free & local
    keep = free & ~local
    return {
        "free": free,
        "local": local,
        "interface": interface,
        "wirebasket": wirebasket,
        "local_bubble": local_bubble,
        "keep": keep,
        "other": hidden,
    }


@dataclass(frozen=True)
class ResponseBasis:
    """Dense response basis returned by :func:`BlockKrylovBasis`."""

    vectors: np.ndarray
    active_dofs: np.ndarray

    @property
    def rank(self) -> int:
        return int(self.vectors.shape[1])

    @property
    def ndof(self) -> int:
        return int(self.vectors.shape[0])

    @property
    def active_count(self) -> int:
        return int(self.active_dofs.shape[0])

    @property
    def compression_ratio(self) -> float:
        return float(self.rank / self.ndof) if self.ndof else 0.0

    @property
    def inactive_dofs(self) -> int:
        return self.ndof - self.active_count

    @property
    def port_visible_dofs(self) -> int:
        return self.rank

    @property
    def port_invisible_dofs(self) -> int:
        return self.active_count - self.rank

    @property
    def eddy_visible_dofs(self) -> int:
        """Number of retained Eddy-Visible Response Space coordinates."""

        return self.rank

    @property
    def eddy_invisible_dofs(self) -> int:
        """Number of eliminated Eddy-Invisible DoFs (eddy bubbles)."""

        return self.active_count - self.rank

    def diagnostics(self) -> dict[str, int | float]:
        """Return compact DoF-reduction diagnostics."""

        return {
            "ndof": self.ndof,
            "active_dofs": self.active_count,
            "rank": self.rank,
            "port_visible_dofs": self.port_visible_dofs,
            "eddy_visible_dofs": self.eddy_visible_dofs,
            "compression_ratio": self.compression_ratio,
            "inactive_dofs": self.inactive_dofs,
            "port_invisible_dofs": self.port_invisible_dofs,
            "eddy_invisible_dofs": self.eddy_invisible_dofs,
            "eliminated_dofs": self.ndof - self.rank,
        }


@dataclass(frozen=True)
class EVRSBasis(ResponseBasis):
    """Eddy-Visible Response Space basis embedded in a high-order parent space.

    ``vectors`` contains full parent-space coefficient columns, while
    ``active_dofs`` records the parent DoFs that participated in the response
    solve.  The retained columns are the EVRS coordinates; active DoFs that do
    not enter this span are Eddy-Invisible DoFs (EIDs, informal: eddy bubbles).
    """

    port_count: int | None = None
    krylov_steps: int | None = None
    construction: str = "block-krylov"
    parent_space: str = "HCurl"
    pre_current_gram_rank: int | None = None
    current_gram_rtol: float | None = None
    current_gram_relative_eigenvalues: tuple[float, ...] | None = None

    def diagnostics(self) -> dict[str, object]:
        """Return DoF-reduction diagnostics plus EVRS construction metadata."""

        info = dict(super().diagnostics())
        info.update(
            {
                "port_count": self.port_count,
                "krylov_steps": self.krylov_steps,
                "construction": self.construction,
                "parent_space": self.parent_space,
            }
        )
        if self.pre_current_gram_rank is not None:
            info.update(
                {
                    "pre_current_gram_rank": self.pre_current_gram_rank,
                    "current_gram_rank": self.rank,
                    "current_gram_rtol": self.current_gram_rtol,
                    "current_gram_relative_eigenvalues": list(
                        self.current_gram_relative_eigenvalues or ()
                    ),
                }
            )
        return info


def CompressHCurlResponseInCurrentGram(
    response_basis: ResponseBasis,
    current_basis: SampledCurrentBasis,
    *,
    rtol: float = 1.0e-10,
) -> tuple[ResponseBasis, SampledCurrentBasis]:
    """Remove response directions that are dependent after ``curl(T)``.

    The parent response vectors may be independent in an HCurl mass metric but
    map to zero or nearly dependent current fields.  This helper eigendecomposes
    the sampled current Gram matrix, drops its relative null space, and applies
    the same whitening transform to both parent ``T`` vectors and current
    samples.  The retained current basis therefore has an identity Gram matrix.
    """

    if not isinstance(response_basis, ResponseBasis):
        raise TypeError("response_basis must be a ResponseBasis")
    if not isinstance(current_basis, SampledCurrentBasis):
        raise TypeError("current_basis must be a SampledCurrentBasis")
    rtol = float(rtol)
    if not np.isfinite(rtol) or rtol <= 0.0:
        raise ValueError("rtol must be positive")
    if response_basis.rank != current_basis.n_modes:
        raise ValueError("response and current basis mode counts must match")

    gram = current_basis.mass_matrix()
    gram = 0.5 * (gram + gram.conj().T)
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    scale = max(float(np.max(eigenvalues)), 0.0)
    if scale <= 0.0:
        raise ValueError("current response space has zero Gram energy")
    negative_tolerance = (
        100.0 * max(gram.shape[0], 1) * np.finfo(float).eps * scale
    )
    if float(np.min(eigenvalues)) < -negative_tolerance:
        raise ValueError("current Gram matrix is not positive semidefinite")
    keep = eigenvalues > rtol * scale
    if not np.any(keep):
        raise ValueError("current Gram rank tolerance eliminated every response")

    retained_values = eigenvalues[keep]
    transform = eigenvectors[:, keep] / np.sqrt(retained_values)[np.newaxis, :]
    response_vectors = response_basis.vectors @ transform
    current_modes = np.einsum("mk,mpc->kpc", transform, current_basis.modes)
    current = SampledCurrentBasis(
        points=current_basis.points,
        weights=current_basis.weights,
        modes=current_modes,
        kind=current_basis.kind,
        names=tuple(
            f"{current_basis.kind}_current_gram_{index}"
            for index in range(retained_values.size)
        ),
    )
    relative_eigenvalues = tuple(
        float(value / scale) for value in eigenvalues[::-1]
    )
    if isinstance(response_basis, EVRSBasis):
        response = EVRSBasis(
            vectors=response_vectors,
            active_dofs=response_basis.active_dofs,
            port_count=response_basis.port_count,
            krylov_steps=response_basis.krylov_steps,
            construction=response_basis.construction + "+current-gram",
            parent_space=response_basis.parent_space,
            pre_current_gram_rank=response_basis.rank,
            current_gram_rtol=rtol,
            current_gram_relative_eigenvalues=relative_eigenvalues,
        )
    else:
        response = ResponseBasis(
            vectors=response_vectors,
            active_dofs=response_basis.active_dofs,
        )
    return response, current


def _active_indices(free_dofs, n: int) -> np.ndarray:
    if free_dofs is None:
        return np.arange(n, dtype=int)
    arr = np.asarray(free_dofs)
    if arr.dtype == bool:
        if arr.ndim != 1 or arr.shape[0] != n:
            raise ValueError(f"free_dofs boolean mask must have shape ({n},)")
        return np.flatnonzero(arr)
    idx = arr.astype(int, copy=False).ravel()
    if np.any(idx < 0) or np.any(idx >= n):
        raise ValueError("free_dofs contains an out-of-range index")
    return idx


def _metric_inner(metric: np.ndarray, a: np.ndarray, b: np.ndarray):
    return np.vdot(a, metric @ b)


def _append_metric_orthonormal(
    columns: list[np.ndarray],
    candidate: np.ndarray,
    metric: np.ndarray,
    *,
    rtol: float,
) -> None:
    v = np.array(candidate, copy=True)
    candidate_norm2 = _metric_inner(metric, v, v)
    candidate_norm = float(np.sqrt(max(np.real(candidate_norm2), 0.0)))
    if candidate_norm == 0.0:
        return
    for _ in range(2):
        for q in columns:
            v -= q * _metric_inner(metric, q, v)
    norm2 = _metric_inner(metric, v, v)
    norm = float(np.sqrt(max(np.real(norm2), 0.0)))
    relative_floor = max(
        rtol,
        10.0 * candidate.size * np.finfo(float).eps,
    )
    if norm > relative_floor * candidate_norm:
        columns.append(v / norm)


def BlockKrylovBasis(
    stiffness,
    mass,
    ports,
    steps: int,
    *,
    free_dofs=None,
    metric=None,
    rtol: float = 1.0e-12,
) -> EVRSBasis:
    """Generate a reduced response basis from ``K^{-1} B`` block Krylov vectors.

    This is the algebraic hook for the NGSolve side.  Build high-order HCurl
    matrices in NGSolve, apply tree-cotree/``nograds`` through ``free_dofs``,
    and pass the reduced matrices plus low-order external-field ports here.
    The returned columns are orthonormal in the supplied metric, by default the
    conductive mass matrix.
    """

    if steps < 1:
        raise ValueError("steps must be >= 1")
    if rtol <= 0.0:
        raise ValueError("rtol must be positive")
    k = np.asarray(stiffness)
    m = np.asarray(mass)
    if k.ndim != 2 or k.shape[0] != k.shape[1]:
        raise ValueError("stiffness must be square")
    if m.shape != k.shape:
        raise ValueError("mass must have the same shape as stiffness")
    n = k.shape[0]
    b = np.asarray(ports)
    if b.ndim == 1:
        b = b[:, np.newaxis]
    if b.ndim != 2 or b.shape[0] != n:
        raise ValueError(f"ports must have shape ({n}, p)")

    active = _active_indices(free_dofs, n)
    if active.size == 0:
        raise ValueError("no active dofs")
    idx = np.ix_(active, active)
    kf = k[idx]
    mf = m[idx]
    gf = mf if metric is None else np.asarray(metric)[idx]
    current = b[active, :]
    columns: list[np.ndarray] = []
    for _ in range(steps):
        solved = np.linalg.solve(kf, current)
        for col in range(solved.shape[1]):
            _append_metric_orthonormal(
                columns, solved[:, col], gf, rtol=rtol
            )
        current = mf @ solved

    if columns:
        restricted = np.column_stack(columns)
    else:
        restricted = np.zeros((active.size, 0), dtype=np.result_type(k, m, b))
    full = np.zeros((n, restricted.shape[1]), dtype=restricted.dtype)
    full[active, :] = restricted
    return EVRSBasis(
        vectors=full,
        active_dofs=active,
        port_count=int(b.shape[1]),
        krylov_steps=steps,
        construction="block-krylov",
        parent_space="HCurl",
    )


def NgsolveBlockKrylovBasis(
    stiffness,
    mass,
    ports,
    steps: int,
    *,
    free_dofs=None,
    metric=None,
    rtol: float = 1.0e-12,
) -> EVRSBasis:
    """Generate a block-Krylov response basis from NGSolve matrices/vectors.

    ``stiffness`` and ``mass`` may be assembled NGSolve bilinear forms,
    NGSolve matrices, or NumPy arrays.  ``ports`` may be a single NGSolve
    vector/linear form, a sequence of such vectors, or a NumPy array with shape
    ``(ndof, n_ports)``.  The output coefficient columns can be passed directly
    to :func:`NgsolveHCurlCurlBasis`.
    """

    k = NgsolveMatrixToDense(stiffness)
    m = NgsolveMatrixToDense(mass)
    dtype = np.result_type(k, m)
    b = _ports_to_matrix(ports, dtype=dtype)
    g = None if metric is None else NgsolveMatrixToDense(metric, dtype=dtype)
    basis = BlockKrylovBasis(
        k,
        m,
        b,
        steps,
        free_dofs=free_dofs,
        metric=g,
        rtol=rtol,
    )
    return EVRSBasis(
        vectors=basis.vectors,
        active_dofs=basis.active_dofs,
        port_count=basis.port_count,
        krylov_steps=basis.krylov_steps,
        construction="ngsolve-dense-block-krylov",
        parent_space="HCurl",
    )


def _ngsolve_matrix_object(matrix):
    return getattr(matrix, "mat", matrix)


def _matrix_height(matrix) -> int:
    height = _call_or_value(matrix, "height")
    if height is None:
        raise ValueError("NGSolve matrix does not expose height")
    return int(height)


def _ngsolve_array_to_vector(matrix, values: np.ndarray):
    vec = matrix.CreateColVector()
    data = vec.FV().NumPy()
    if data.shape[0] != values.shape[0]:
        raise ValueError(
            f"array length {values.shape[0]} does not match vector length {data.shape[0]}"
        )
    data[:] = values
    return vec


def _ngsolve_apply_to_array(operator, values: np.ndarray, prototype_matrix) -> np.ndarray:
    vec = _ngsolve_array_to_vector(prototype_matrix, values)
    out = vec.CreateVector()
    out.data = operator * vec
    return np.array(out.FV().NumPy(), copy=True)


def _hdiv_reduction_vectors(fes, vectors) -> np.ndarray:
    arr = np.asarray(vectors)
    if arr.ndim == 1:
        arr = arr[:, np.newaxis]
    if arr.ndim != 2 or arr.shape[0] != fes.ndof:
        raise ValueError(f"vectors must have shape ({fes.ndof}, n_modes)")
    if arr.shape[1] == 0:
        raise ValueError("vectors must contain at least one reduced mode")
    if not np.all(np.isfinite(arr)):
        raise ValueError("vectors contain non-finite values")
    return np.array(arr, copy=True)


def _ngsolve_material_dx(mesh, materials):
    import ngsolve as ng

    if materials is None:
        return ng.dx
    labels = _label_tuple(materials, "materials")
    if not labels:
        return ng.dx
    return ng.dx(definedon=mesh.Materials("|".join(labels)))


def _ngsolve_reduced_operator(operator, vectors: np.ndarray, prototype_matrix) -> np.ndarray:
    applied = np.column_stack(
        [
            _ngsolve_apply_to_array(operator, vectors[:, i], prototype_matrix)
            for i in range(vectors.shape[1])
        ]
    )
    reduced = vectors.conj().T @ applied
    return 0.5 * (reduced + reduced.conj().T)


@dataclass(frozen=True)
class HDivMultipolePortSet:
    """Regular-solid-harmonic applied-H ports for motor response training."""

    fields: tuple[object, ...]
    names: tuple[str, ...]
    degrees: tuple[int, ...]
    center: tuple[float, float, float]
    radius: float

    def __post_init__(self) -> None:
        fields = tuple(self.fields)
        names = tuple(str(name) for name in self.names)
        degrees = tuple(int(degree) for degree in self.degrees)
        if not fields or len(names) != len(fields) or len(degrees) != len(fields):
            raise ValueError("fields, names, and degrees must have the same non-zero length")
        if any(degree < 1 or degree > 3 for degree in degrees):
            raise ValueError("regular-solid-harmonic degrees must lie in [1, 3]")
        center = tuple(float(value) for value in self.center)
        if len(center) != 3 or not np.all(np.isfinite(center)):
            raise ValueError("center must contain three finite coordinates")
        radius = float(self.radius)
        if not np.isfinite(radius) or radius <= 0.0:
            raise ValueError("radius must be positive")
        object.__setattr__(self, "fields", fields)
        object.__setattr__(self, "names", names)
        object.__setattr__(self, "degrees", degrees)
        object.__setattr__(self, "center", center)
        object.__setattr__(self, "radius", radius)

    @property
    def count(self) -> int:
        return len(self.fields)

    def diagnostics(self) -> dict[str, object]:
        counts = {
            str(degree): self.degrees.count(degree)
            for degree in sorted(set(self.degrees))
        }
        return {
            "family": "regular-solid-harmonic-gradient",
            "count": self.count,
            "degrees": list(self.degrees),
            "degree_counts": counts,
            "names": list(self.names),
            "center": list(self.center),
            "radius": float(self.radius),
        }


@dataclass(frozen=True)
class PlanarHarmonicPortSet:
    """Real 2-D harmonic-gradient applied-H ports."""

    fields: tuple[object, ...]
    names: tuple[str, ...]
    degrees: tuple[int, ...]
    center: tuple[float, float]
    radius: float

    def __post_init__(self) -> None:
        fields = tuple(self.fields)
        names = tuple(str(name) for name in self.names)
        degrees = tuple(int(degree) for degree in self.degrees)
        if not fields or len(names) != len(fields) or len(degrees) != len(fields):
            raise ValueError("fields, names, and degrees must have the same non-zero length")
        if any(degree < 1 for degree in degrees):
            raise ValueError("planar harmonic degrees must be positive")
        center = tuple(float(value) for value in self.center)
        if len(center) != 2 or not np.all(np.isfinite(center)):
            raise ValueError("center must contain two finite coordinates")
        radius = float(self.radius)
        if not np.isfinite(radius) or radius <= 0.0:
            raise ValueError("radius must be positive")
        object.__setattr__(self, "fields", fields)
        object.__setattr__(self, "names", names)
        object.__setattr__(self, "degrees", degrees)
        object.__setattr__(self, "center", center)
        object.__setattr__(self, "radius", radius)

    @property
    def count(self) -> int:
        return len(self.fields)

    @property
    def max_degree(self) -> int:
        return max(self.degrees)

    def diagnostics(self) -> dict[str, object]:
        return {
            "family": "planar-harmonic-gradient",
            "count": self.count,
            "max_degree": self.max_degree,
            "degrees": list(self.degrees),
            "degree_counts": {
                str(degree): self.degrees.count(degree)
                for degree in sorted(set(self.degrees))
            },
            "names": list(self.names),
            "center": list(self.center),
            "radius": float(self.radius),
        }


def NgsolvePlanarHarmonicPorts(
    mesh,
    *,
    max_degree: int = 3,
    center=None,
    radius=None,
) -> PlanarHarmonicPortSet:
    """Return gradients of ``Re/Im[((x+iy)/radius)^l]`` for ``l=1..L``."""

    import ngsolve as ng

    if mesh.dim != 2:
        raise ValueError("planar harmonic ports require a 2-D mesh")
    max_degree = int(max_degree)
    if max_degree < 1:
        raise ValueError("max_degree must be positive")
    coordinates = np.asarray(
        [tuple(vertex.point)[:2] for vertex in mesh.vertices],
        dtype=float,
    )
    if coordinates.ndim != 2 or coordinates.shape[1] != 2 or not coordinates.size:
        raise ValueError("mesh must expose two-dimensional vertices")
    if center is None:
        center_values = 0.5 * (
            np.min(coordinates, axis=0) + np.max(coordinates, axis=0)
        )
    else:
        center_values = np.asarray(center, dtype=float).reshape(-1)
    if center_values.size != 2 or not np.all(np.isfinite(center_values)):
        raise ValueError("center must contain two finite coordinates")
    if radius is None:
        radius_value = float(
            np.max(np.linalg.norm(coordinates - center_values, axis=1))
        )
    else:
        radius_value = float(radius)
    if not np.isfinite(radius_value) or radius_value <= 0.0:
        raise ValueError("radius must be positive")

    x = (ng.x - float(center_values[0])) / radius_value
    y = (ng.y - float(center_values[1])) / radius_value
    real_power = 1.0
    imag_power = 0.0
    fields = []
    names = []
    degrees = []
    for degree in range(1, max_degree + 1):
        fields.extend(
            (
                ng.CoefficientFunction(
                    (degree * real_power, -degree * imag_power)
                ),
                ng.CoefficientFunction(
                    (degree * imag_power, degree * real_power)
                ),
            )
        )
        names.extend((f"ph_l{degree}_cos", f"ph_l{degree}_sin"))
        degrees.extend((degree, degree))
        real_power, imag_power = (
            real_power * x - imag_power * y,
            real_power * y + imag_power * x,
        )
    return PlanarHarmonicPortSet(
        fields=tuple(fields),
        names=tuple(names),
        degrees=tuple(degrees),
        center=tuple(center_values),
        radius=radius_value,
    )


def NgsolveHDivRegularSolidHarmonicPorts(
    mesh,
    *,
    max_degree: int = 3,
    center=None,
    radius=None,
) -> HDivMultipolePortSet:
    """Return real 3-D harmonic-gradient ports through degree three.

    The scalar potentials are regular solid harmonics in normalized Cartesian
    coordinates.  Their gradients span uniform, linear-gradient, and quadratic
    applied-H content without introducing fictitious sources inside the body.
    """

    import ngsolve as ng

    if mesh.dim != 3:
        raise ValueError("regular-solid-harmonic HDiv ports require a 3-D mesh")
    max_degree = int(max_degree)
    if max_degree < 1 or max_degree > 3:
        raise ValueError("max_degree must lie in [1, 3]")
    coordinates = np.asarray(
        [tuple(vertex.point) for vertex in mesh.vertices],
        dtype=float,
    )
    if coordinates.ndim != 2 or coordinates.shape[1] != 3 or not coordinates.size:
        raise ValueError("mesh must expose three-dimensional vertices")
    if center is None:
        center_values = 0.5 * (
            np.min(coordinates, axis=0) + np.max(coordinates, axis=0)
        )
    else:
        center_values = np.asarray(center, dtype=float).reshape(-1)
    if center_values.size != 3 or not np.all(np.isfinite(center_values)):
        raise ValueError("center must contain three finite coordinates")
    if radius is None:
        radius_value = float(
            np.max(np.linalg.norm(coordinates - center_values, axis=1))
        )
    else:
        radius_value = float(radius)
    if not np.isfinite(radius_value) or radius_value <= 0.0:
        raise ValueError("radius must be positive")

    x = (ng.x - float(center_values[0])) / radius_value
    y = (ng.y - float(center_values[1])) / radius_value
    z = (ng.z - float(center_values[2])) / radius_value
    fields = [
        ng.CoefficientFunction((1.0, 0.0, 0.0)),
        ng.CoefficientFunction((0.0, 1.0, 0.0)),
        ng.CoefficientFunction((0.0, 0.0, 1.0)),
    ]
    names = ["rsh_l1_x", "rsh_l1_y", "rsh_l1_z"]
    degrees = [1, 1, 1]
    if max_degree >= 2:
        fields.extend(
            (
                ng.CoefficientFunction((y, x, 0.0)),
                ng.CoefficientFunction((z, 0.0, x)),
                ng.CoefficientFunction((0.0, z, y)),
                ng.CoefficientFunction((2.0 * x, -2.0 * y, 0.0)),
                ng.CoefficientFunction((-2.0 * x, -2.0 * y, 4.0 * z)),
            )
        )
        names.extend(
            ("rsh_l2_xy", "rsh_l2_xz", "rsh_l2_yz", "rsh_l2_x2_y2", "rsh_l2_z2")
        )
        degrees.extend((2,) * 5)
    if max_degree >= 3:
        fields.extend(
            (
                ng.CoefficientFunction(
                    (-6.0 * x * z, -6.0 * y * z, 6.0 * z * z - 3.0 * x * x - 3.0 * y * y)
                ),
                ng.CoefficientFunction(
                    (4.0 * z * z - 3.0 * x * x - y * y, -2.0 * x * y, 8.0 * x * z)
                ),
                ng.CoefficientFunction(
                    (-2.0 * x * y, 4.0 * z * z - x * x - 3.0 * y * y, 8.0 * y * z)
                ),
                ng.CoefficientFunction((2.0 * x * z, -2.0 * y * z, x * x - y * y)),
                ng.CoefficientFunction((2.0 * y * z, 2.0 * x * z, 2.0 * x * y)),
                ng.CoefficientFunction((3.0 * x * x - 3.0 * y * y, -6.0 * x * y, 0.0)),
                ng.CoefficientFunction((6.0 * x * y, 3.0 * x * x - 3.0 * y * y, 0.0)),
            )
        )
        names.extend(
            (
                "rsh_l3_z3",
                "rsh_l3_xz2",
                "rsh_l3_yz2",
                "rsh_l3_z_x2_y2",
                "rsh_l3_xyz",
                "rsh_l3_x3_3xy2",
                "rsh_l3_3x2y_y3",
            )
        )
        degrees.extend((3,) * 7)
    return HDivMultipolePortSet(
        fields=tuple(fields),
        names=tuple(names),
        degrees=tuple(degrees),
        center=tuple(center_values),
        radius=radius_value,
    )


def _ngsolve_field_tuple(fields):
    import ngsolve as ng

    if fields is None:
        return ()
    if isinstance(fields, HDivMultipolePortSet):
        return fields.fields
    if isinstance(fields, (tuple, list)):
        if len(fields) in (2, 3) and all(np.isscalar(value) for value in fields):
            return (ng.CoefficientFunction(tuple(fields)),)
        return tuple(fields)
    return (fields,)


def _ngsolve_parent_external_field_rhs(mesh, fes, external_fields, *, materials=None):
    import ngsolve as ng

    fields = _ngsolve_field_tuple(external_fields)
    if not fields:
        raise ValueError("external_fields must not be empty")
    test = fes.TestFunction()
    measure = _ngsolve_material_dx(mesh, materials)
    columns = []
    for field in fields:
        linear = ng.LinearForm(fes)
        linear += field * test * measure
        linear.Assemble()
        columns.append(NgsolveVectorToArray(linear))
    return np.column_stack(columns)


def NgsolveHDivExternalFieldRHS(
    mesh,
    fes,
    vectors,
    external_fields,
    *,
    materials=None,
) -> np.ndarray:
    """Project one or more applied H fields onto reduced HDiv coordinates.

    The result is ``Q^* M H_ext`` and therefore matches the right-hand side of
    the production linear HDiv-MMM system ``((1/chi) M + N) m = M H_ext``.
    This function can be called repeatedly for rotor-angle or stator-current
    sweeps without rebuilding the demagnetizing operator.
    """

    q = _hdiv_reduction_vectors(fes, vectors)
    parent_rhs = _ngsolve_parent_external_field_rhs(
        mesh,
        fes,
        external_fields,
        materials=materials,
    )
    return q.conj().T @ parent_rhs


@dataclass(frozen=True)
class HDivMMMReducedModel:
    """Projected production HDiv-MMM material/demagnetizing model."""

    mesh: object
    fes: object
    parent_vectors: np.ndarray
    magnetization_basis: SampledMagnetizationBasis
    mass: np.ndarray
    demag: np.ndarray
    magnetic_operator: np.ndarray
    magnetic_rhs: np.ndarray | None
    mu_r: float
    materials: tuple[str, ...]
    demag_backend: object
    basis_generation: dict[str, object] | None = None
    parent_family: str = "unspecified"
    parent_order: int | None = None

    def __post_init__(self) -> None:
        q = _hdiv_reduction_vectors(self.fes, self.parent_vectors)
        modes = self.magnetization_basis.n_modes
        if q.shape[1] != modes:
            raise ValueError("parent_vectors columns must match magnetization basis modes")
        mass = _square_matrix(self.mass, modes, "mass")
        demag = _square_matrix(self.demag, modes, "demag")
        operator = _square_matrix(
            self.magnetic_operator,
            modes,
            "magnetic_operator",
        )
        rhs = self.magnetic_rhs
        if rhs is not None:
            rhs = _port_rhs_matrix(rhs, modes)
        if not np.isfinite(self.mu_r) or self.mu_r <= 1.0:
            raise ValueError("mu_r must be greater than 1")
        parent_family = str(self.parent_family).strip().upper()
        if parent_family not in ("BDM", "RT", "UNSPECIFIED"):
            raise ValueError("parent_family must be 'BDM', 'RT', or 'unspecified'")
        parent_order = self.parent_order
        if parent_order is not None:
            parent_order = int(parent_order)
            if parent_order < 1:
                raise ValueError("parent_order must be positive or None")
        object.__setattr__(self, "parent_vectors", q)
        object.__setattr__(self, "mass", mass)
        object.__setattr__(self, "demag", demag)
        object.__setattr__(self, "magnetic_operator", operator)
        object.__setattr__(self, "magnetic_rhs", rhs)
        object.__setattr__(self, "materials", tuple(self.materials))
        object.__setattr__(self, "parent_family", parent_family)
        object.__setattr__(self, "parent_order", parent_order)
        object.__setattr__(
            self,
            "basis_generation",
            None if self.basis_generation is None else dict(self.basis_generation),
        )

    @property
    def n_modes(self) -> int:
        return self.magnetization_basis.n_modes

    @property
    def parent_ndof(self) -> int:
        return int(self.parent_vectors.shape[0])

    @property
    def chi(self) -> float:
        return self.mu_r - 1.0

    @property
    def basis(self) -> SampledMagnetizationBasis:
        return self.magnetization_basis

    def external_field_rhs(self, external_fields, *, materials=None) -> np.ndarray:
        """Project new rotor/stator applied-H columns without rebuilding N."""

        if materials is None:
            materials = self.materials
        return NgsolveHDivExternalFieldRHS(
            self.mesh,
            self.fes,
            self.parent_vectors,
            external_fields,
            materials=materials,
        )

    def reconstruct_parent(self, reduced_coefficients) -> np.ndarray:
        """Lift reduced magnetization coordinates to the parent HDiv space."""

        coefficients = np.asarray(reduced_coefficients)
        if coefficients.ndim == 1:
            coefficients = coefficients[:, np.newaxis]
        if coefficients.ndim != 2 or coefficients.shape[0] != self.n_modes:
            raise ValueError(
                f"reduced_coefficients must have shape ({self.n_modes}, n_rhs)"
            )
        return self.parent_vectors @ coefficients

    def diagnostics(self) -> dict[str, object]:
        native_gram = getattr(self.demag_backend, "_G", None)
        return {
            "parent_space": "HDiv",
            "parent_family": self.parent_family,
            "parent_order": self.parent_order,
            "parent_ndof": self.parent_ndof,
            "reduced_modes": self.n_modes,
            "compression_ratio": float(self.n_modes / self.parent_ndof),
            "mu_r": float(self.mu_r),
            "chi": float(self.chi),
            "materials": list(self.materials),
            "mass_hermitian_error": _relative_hermitian_error(self.mass),
            "demag_hermitian_error": _relative_hermitian_error(self.demag),
            "operator_hermitian_error": _relative_hermitian_error(
                self.magnetic_operator
            ),
            "min_mass_eigenvalue": _min_hermitian_eigenvalue(self.mass),
            "min_demag_eigenvalue": _min_hermitian_eigenvalue(self.demag),
            "min_operator_eigenvalue": _min_hermitian_eigenvalue(
                self.magnetic_operator
            ),
            "demag_frobenius_norm": float(np.linalg.norm(self.demag)),
            "has_rhs": self.magnetic_rhs is not None,
            "rhs_columns": (
                0 if self.magnetic_rhs is None else int(self.magnetic_rhs.shape[1])
            ),
            "demag_backend": self.demag_backend.__class__.__name__,
            "demag_hmatrix_backend": (
                None if native_gram is None else native_gram.__class__.__name__
            ),
            "demag_hmatrix_active": native_gram is not None,
            "basis_generation": self.basis_generation,
        }


def NgsolveHDivMMMReduction(
    mesh,
    fes,
    vectors,
    *,
    mu_r: float,
    external_fields=None,
    intorder: int = 2,
    materials=None,
    names=None,
    mass=None,
    demag_operator=None,
    demag_intorder=None,
    demag_eps: float = 1.0e-7,
    demag_leafsize: int = 16,
    demag_eta: float = 2.0,
    demag_far_quad: int = 3,
    demag_ho_far_factor: float = 2.0,
    basis_generation=None,
    parent_family: str = "unspecified",
    parent_order: int | None = None,
) -> HDivMMMReducedModel:
    """Project Radia's production HDiv mass and demag operator onto ``vectors``.

    No full dense parent matrix is formed.  ``M`` and ``N=B^T G B`` are
    applied column-by-column through NGSolve/Radia native matrices, then
    Galerkin-projected as ``Q^* M Q`` and ``Q^* N Q``.
    """

    import ngsolve as ng

    if mesh.dim != 3:
        raise ValueError("NgsolveHDivMMMReduction currently requires a 3-D mesh")
    mu_r = float(mu_r)
    if not np.isfinite(mu_r) or mu_r <= 1.0:
        raise ValueError("mu_r must be greater than 1")
    q = _hdiv_reduction_vectors(fes, vectors)
    labels = _label_tuple(materials, "materials")
    basis = NgsolveHDivMagnetizationBasis(
        mesh,
        fes,
        q,
        intorder=intorder,
        materials=materials,
        names=names,
    )

    if mass is None:
        trial, test = fes.TnT()
        mass = ng.BilinearForm(fes)
        mass += trial * test * _ngsolve_material_dx(mesh, materials)
        mass.Assemble()
    mass_matrix = _ngsolve_matrix_object(mass)

    if demag_operator is None:
        from ._vim import DemagOperator

        demag_operator = DemagOperator(
            fes,
            intorder=demag_intorder,
            eps=demag_eps,
            leafsize=demag_leafsize,
            eta=demag_eta,
            far_quad=demag_far_quad,
            ho_far_factor=demag_ho_far_factor,
        )
    demag_matrix = _ngsolve_matrix_object(demag_operator)
    reduced_mass = _ngsolve_reduced_operator(mass_matrix, q, mass_matrix)
    native_gram = getattr(demag_operator, "_G", None)
    if native_gram is not None and hasattr(native_gram, "apply_configured_demag"):
        applied_demag = np.column_stack(
            [
                np.asarray(native_gram.apply_configured_demag(q[:, i], True))
                for i in range(q.shape[1])
            ]
        )
        reduced_demag = q.conj().T @ applied_demag
        reduced_demag = 0.5 * (reduced_demag + reduced_demag.conj().T)
    else:
        reduced_demag = _ngsolve_reduced_operator(demag_matrix, q, mass_matrix)
    reduced_operator = reduced_mass / (mu_r - 1.0) + reduced_demag
    rhs = None
    if external_fields is not None:
        rhs = NgsolveHDivExternalFieldRHS(
            mesh,
            fes,
            q,
            external_fields,
            materials=materials,
        )
    return HDivMMMReducedModel(
        mesh=mesh,
        fes=fes,
        parent_vectors=q,
        magnetization_basis=basis,
        mass=reduced_mass,
        demag=reduced_demag,
        magnetic_operator=reduced_operator,
        magnetic_rhs=rhs,
        mu_r=mu_r,
        materials=labels,
        demag_backend=demag_operator,
        basis_generation=basis_generation,
        parent_family=parent_family,
        parent_order=parent_order,
    )


def _hdiv_response_field_names(fields, supplied_names, prefix: str) -> tuple[str, ...]:
    values = _ngsolve_field_tuple(fields)
    if supplied_names is None:
        if isinstance(fields, HDivMultipolePortSet):
            return fields.names
        return tuple(f"{prefix}_{index}" for index in range(len(values)))
    names = tuple(str(name) for name in supplied_names)
    if len(names) != len(values):
        raise ValueError(f"{prefix}_names must match the number of fields")
    return names


def _hdiv_response_apply_columns(apply, vectors: np.ndarray) -> np.ndarray:
    return np.column_stack([apply(vectors[:, index]) for index in range(vectors.shape[1])])


def _protected_magnetic_response_pod(
    snapshot_matrix,
    parent_rhs,
    apply_operator,
    *,
    physical_count: int,
    port_names,
    training_names,
    port_weights=None,
    normalize_ports: bool = True,
    pod_rtol: float = 1.0e-10,
    max_modes=None,
):
    """Protect physical responses, then energy-POD their training complement."""

    snapshots = np.asarray(snapshot_matrix)
    rhs = np.asarray(parent_rhs)
    if snapshots.ndim != 2 or rhs.shape != snapshots.shape:
        raise ValueError("snapshot_matrix and parent_rhs must have the same 2-D shape")
    count = snapshots.shape[1]
    physical_count = int(physical_count)
    if count < 1 or physical_count < 0 or physical_count > count:
        raise ValueError("physical_count must lie between zero and the port count")
    port_names = tuple(str(name) for name in port_names)
    training_names = tuple(str(name) for name in training_names)
    if len(port_names) != count or len(training_names) != count - physical_count:
        raise ValueError("port names do not match the snapshot partition")
    pod_rtol = float(pod_rtol)
    if pod_rtol <= 0.0:
        raise ValueError("pod_rtol must be positive")
    if max_modes is not None and int(max_modes) < 1:
        raise ValueError("max_modes must be positive or None")
    if port_weights is None:
        weights = np.ones(count)
    else:
        weights = np.asarray(port_weights, dtype=float).reshape(-1)
        if (
            weights.size != count
            or not np.all(np.isfinite(weights))
            or np.any(weights <= 0.0)
        ):
            raise ValueError("port_weights must contain one positive finite value per field")

    energies = np.real(np.sum(snapshots.conj() * rhs, axis=0))
    energy_floor = np.finfo(float).eps * max(float(np.max(np.abs(energies))), 1.0)
    if np.any(energies <= energy_floor):
        bad = [port_names[index] for index in np.flatnonzero(energies <= energy_floor)]
        raise ValueError(f"response ports have zero or non-positive magnetic energy: {bad}")
    snapshot_actions = _hdiv_response_apply_columns(apply_operator, snapshots)

    def energy_pod(values, actions, scales):
        if values.shape[1] == 0:
            empty = np.zeros((snapshots.shape[0], 0), dtype=values.dtype)
            return empty, empty.copy(), np.zeros(0), 0
        weighted_values = values * scales[np.newaxis, :]
        weighted_actions = actions * scales[np.newaxis, :]
        correlation = weighted_values.conj().T @ weighted_actions
        correlation = 0.5 * (correlation + correlation.conj().T)
        eigenvalues, eigenvectors = np.linalg.eigh(correlation)
        order = np.argsort(eigenvalues)[::-1]
        eigenvalues = np.real(eigenvalues[order])
        eigenvectors = eigenvectors[:, order]
        threshold = pod_rtol * max(float(eigenvalues[0]), np.finfo(float).tiny)
        rank = int(np.count_nonzero(eigenvalues > threshold))
        if rank == 0:
            empty = np.zeros((snapshots.shape[0], 0), dtype=values.dtype)
            return empty, empty.copy(), eigenvalues, 0
        transform = eigenvectors[:, :rank] @ np.diag(
            1.0 / np.sqrt(eigenvalues[:rank])
        )
        basis = weighted_values @ transform
        basis_actions = weighted_actions @ transform
        basis_gram = basis.conj().T @ basis_actions
        basis_gram = 0.5 * (basis_gram + basis_gram.conj().T)
        try:
            cholesky = np.linalg.cholesky(basis_gram)
        except np.linalg.LinAlgError as exc:
            raise RuntimeError("response basis lost magnetic-energy independence") from exc
        correction = np.linalg.inv(cholesky.conj().T)
        return basis @ correction, basis_actions @ correction, eigenvalues, rank

    if physical_count:
        physical_scales = weights[:physical_count] / (
            np.sqrt(energies[:physical_count])
            if normalize_ports
            else np.ones(physical_count)
        )
        protected_vectors, protected_actions, protected_eigenvalues, protected_rank = (
            energy_pod(
                snapshots[:, :physical_count],
                snapshot_actions[:, :physical_count],
                physical_scales,
            )
        )
    else:
        protected_vectors = np.zeros((snapshots.shape[0], 0))
        protected_actions = np.zeros((snapshots.shape[0], 0))
        protected_eigenvalues = np.zeros(0)
        protected_rank = 0

    training_values = snapshots[:, physical_count:]
    training_actions = snapshot_actions[:, physical_count:]
    if protected_rank and training_values.shape[1]:
        protected_gram = protected_vectors.conj().T @ protected_actions
        projection = np.linalg.solve(
            protected_gram,
            protected_vectors.conj().T @ training_actions,
        )
        training_values = training_values - protected_vectors @ projection
        training_actions = training_actions - protected_actions @ projection
    training_energies = np.real(np.sum(training_values.conj() * training_actions, axis=0))
    if training_energies.size:
        training_floor = np.finfo(float).eps * max(
            float(np.max(np.abs(energies[physical_count:]))),
            1.0,
        )
        active_training = training_energies > training_floor
    else:
        active_training = np.zeros(0, dtype=bool)
    dependent_training_names = [
        training_names[index] for index in np.flatnonzero(~active_training)
    ]
    active_training_values = training_values[:, active_training]
    active_training_actions = training_actions[:, active_training]
    active_training_energies = training_energies[active_training]
    active_training_weights = weights[physical_count:][active_training]
    training_scales = active_training_weights / (
        np.sqrt(active_training_energies)
        if normalize_ports
        else np.ones(active_training_energies.size)
    )
    training_vectors, training_basis_actions, training_eigenvalues, training_rank = (
        energy_pod(
            active_training_values,
            active_training_actions,
            training_scales,
        )
    )

    available_vectors = np.column_stack((protected_vectors, training_vectors))
    available_actions = np.column_stack((protected_actions, training_basis_actions))
    available = int(available_vectors.shape[1])
    if available < 1:
        raise RuntimeError("response POD retained no magnetic modes")
    available_gram = available_vectors.conj().T @ available_actions
    available_gram = 0.5 * (available_gram + available_gram.conj().T)
    try:
        gram_cholesky = np.linalg.cholesky(available_gram)
    except np.linalg.LinAlgError as exc:
        raise RuntimeError("response basis lost magnetic-energy independence") from exc
    correction = np.linalg.inv(gram_cholesky.conj().T)
    available_vectors = available_vectors @ correction
    available_actions = available_actions @ correction
    available_gram = available_vectors.conj().T @ available_actions
    available_gram = 0.5 * (available_gram + available_gram.conj().T)
    retained = available if max_modes is None else min(available, int(max_modes))
    if retained < protected_rank:
        raise ValueError(
            f"max_modes={retained} cannot discard {protected_rank} protected physical-response modes"
        )
    vectors = available_vectors[:, :retained]
    actions = available_actions[:, :retained]
    gram = vectors.conj().T @ actions
    energy_orthonormality_error = float(
        np.linalg.norm(gram - np.eye(retained)) / max(np.sqrt(retained), 1.0)
    )

    protected_response_error = None
    if physical_count and protected_rank:
        protected_basis = available_vectors[:, :protected_rank]
        protected_basis_gram = available_gram[:protected_rank, :protected_rank]
        protected_errors = []
        for index in range(physical_count):
            coefficients = np.linalg.solve(
                protected_basis_gram,
                protected_basis.conj().T @ rhs[:, index],
            )
            error = snapshots[:, index] - protected_basis @ coefficients
            error_energy = max(float(np.real(np.vdot(error, apply_operator(error)))), 0.0)
            protected_errors.append(float(np.sqrt(error_energy / energies[index])))
        protected_response_error = float(max(protected_errors))

    truncation_curve = []
    for rank in range(1, available + 1):
        rank_vectors = available_vectors[:, :rank]
        rank_gram = available_gram[:rank, :rank]
        reduced_rhs = rank_vectors.conj().T @ rhs
        reduced_coefficients = np.linalg.solve(rank_gram, reduced_rhs)
        captured = np.real(np.sum(reduced_rhs.conj() * reduced_coefficients, axis=0))
        errors = np.sqrt(np.maximum(energies - captured, 0.0) / energies)
        truncation_curve.append(
            {
                "rank": rank,
                "max_all_response_relative_energy_error": float(np.max(errors)),
                "rms_all_response_relative_energy_error": float(
                    np.sqrt(np.mean(errors * errors))
                ),
                "max_physical_response_relative_energy_error": (
                    None
                    if not physical_count
                    else (
                        protected_response_error
                        if rank >= protected_rank
                        else float(np.max(errors[:physical_count]))
                    )
                ),
            }
        )
    response_errors = []
    for index in range(count):
        coefficients = np.linalg.solve(gram, vectors.conj().T @ rhs[:, index])
        error = snapshots[:, index] - vectors @ coefficients
        error_energy = max(float(np.real(np.vdot(error, apply_operator(error)))), 0.0)
        response_errors.append(float(np.sqrt(error_energy / energies[index])))
    retained_curve = truncation_curve[retained - 1]
    retained_curve["max_all_response_relative_energy_error"] = float(max(response_errors))
    retained_curve["rms_all_response_relative_energy_error"] = float(
        np.sqrt(np.mean(np.square(response_errors)))
    )
    retained_curve["max_physical_response_relative_energy_error"] = (
        None
        if not physical_count
        else float(max(response_errors[:physical_count]))
    )
    return vectors, {
        "normalize_ports": bool(normalize_ports),
        "port_weights": weights.tolist(),
        "pod_rtol": pod_rtol,
        "protected_physical_modes": protected_rank,
        "protected_eigenvalues": protected_eigenvalues.tolist(),
        "training_response_modes": training_rank,
        "pod_eigenvalues": training_eigenvalues.tolist(),
        "available_modes": available,
        "retained_modes": retained,
        "discarded_modes": int(count - retained),
        "dependent_port_directions": int(count - available),
        "dependent_training_ports": dependent_training_names,
        "pod_truncated_modes": int(available - retained),
        "pod_truncation_curve": truncation_curve,
        "energy_orthonormality_error": energy_orthonormality_error,
        "response_relative_energy_errors": response_errors,
        "max_response_relative_energy_error": float(max(response_errors)),
    }


def NgsolveHDivMMMResponseReduction(
    mesh,
    fes,
    *,
    mu_r: float,
    external_fields=None,
    training_fields=None,
    external_names=None,
    training_names=None,
    port_weights=None,
    normalize_ports: bool = True,
    pod_rtol: float = 1.0e-10,
    max_modes=None,
    solve_tol: float = 1.0e-10,
    solve_maxit: int = 5000,
    inverse: str = "sparsecholesky",
    intorder: int = 2,
    materials=None,
    mass=None,
    demag_operator=None,
    demag_intorder=None,
    demag_eps: float = 1.0e-7,
    demag_leafsize: int = 16,
    demag_eta: float = 2.0,
    demag_far_quad: int = 3,
    demag_ho_far_factor: float = 2.0,
    parent_family: str = "unspecified",
    parent_order: int | None = None,
) -> HDivMMMReducedModel:
    """Build an operator- and excitation-adapted HDiv-MMM response basis.

    Full parent snapshots solve ``(M/chi + N) m_j = M H_j``.  The snapshots
    are normalized per port by default, compressed in the same magnetic-energy
    inner product, and returned through the ordinary :class:`HDivMMMReducedModel`.
    ``training_fields`` enrich the basis but do not add physical RHS columns;
    this is where rotor-angle samples and regular-solid-harmonic ports belong.
    """

    import ngsolve as ng

    if mesh.dim != 3:
        raise ValueError("NgsolveHDivMMMResponseReduction currently requires a 3-D mesh")
    mu_r = float(mu_r)
    if not np.isfinite(mu_r) or mu_r <= 1.0:
        raise ValueError("mu_r must be greater than 1")
    pod_rtol = float(pod_rtol)
    solve_tol = float(solve_tol)
    solve_maxit = int(solve_maxit)
    if pod_rtol <= 0.0 or solve_tol <= 0.0 or solve_maxit < 1:
        raise ValueError("pod_rtol, solve_tol, and solve_maxit must be positive")
    if max_modes is not None and int(max_modes) < 1:
        raise ValueError("max_modes must be positive or None")

    physical_fields = _ngsolve_field_tuple(external_fields)
    enrichment_fields = _ngsolve_field_tuple(training_fields)
    if not physical_fields and not enrichment_fields:
        raise ValueError("at least one external or training field is required")
    physical_names = _hdiv_response_field_names(
        external_fields,
        external_names,
        "external",
    )
    enrichment_names = _hdiv_response_field_names(
        training_fields,
        training_names,
        "training",
    )
    fields = physical_fields + enrichment_fields
    port_names = physical_names + enrichment_names
    count = len(fields)
    if port_weights is None:
        weights = np.ones(count)
    else:
        weights = np.asarray(port_weights, dtype=float).reshape(-1)
        if weights.size != count or not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
            raise ValueError("port_weights must contain one positive finite value per field")

    custom_mass = mass is not None
    if mass is None:
        trial, test = fes.TnT()
        mass = ng.BilinearForm(fes)
        mass += trial * test * _ngsolve_material_dx(mesh, materials)
        mass.Assemble()
    mass_matrix = _ngsolve_matrix_object(mass)
    if demag_operator is None:
        from ._vim import DemagOperator

        demag_operator = DemagOperator(
            fes,
            intorder=demag_intorder,
            eps=demag_eps,
            leafsize=demag_leafsize,
            eta=demag_eta,
            far_quad=demag_far_quad,
            ho_far_factor=demag_ho_far_factor,
        )
    demag_matrix = _ngsolve_matrix_object(demag_operator)
    parent_rhs = _ngsolve_parent_external_field_rhs(
        mesh,
        fes,
        fields,
        materials=materials,
    )
    if np.iscomplexobj(parent_rhs) and np.max(np.abs(parent_rhs.imag)) > 1.0e-13:
        raise ValueError("response training currently requires real applied-H fields")
    parent_rhs = np.asarray(parent_rhs.real, dtype=float)
    inv_chi = 1.0 / (mu_r - 1.0)

    labels = _label_tuple(materials, "materials")
    mesh_labels = tuple(str(label) for label in mesh.GetMaterials())
    whole_mesh_mass = not labels or set(labels) == set(mesh_labels)
    native_gram = getattr(demag_operator, "_G", None)
    use_native = (
        not custom_mass
        and whole_mesh_mass
        and native_gram is not None
        and hasattr(native_gram, "solve_configured_linear_material_mass_riesz")
    )
    snapshots = []
    iterations = []
    residuals = []
    if use_native:
        backend = "radia-cpp-mass-riesz-cg"

        def apply_operator(values):
            return np.asarray(native_gram.apply_configured_demag(values, True)) + (
                inv_chi
                * np.asarray(native_gram.apply_configured_geometry_mass(values))
            )

        for index in range(count):
            result = native_gram.solve_configured_linear_material_mass_riesz(
                inv_chi,
                parent_rhs[:, index],
                solve_tol,
                solve_maxit,
                True,
            )
            snapshot = np.asarray(result["m"], dtype=float)
            snapshots.append(snapshot)
            iterations.append(int(result["iters"]))
            residuals.append(
                float(
                    np.linalg.norm(apply_operator(snapshot) - parent_rhs[:, index])
                    / max(np.linalg.norm(parent_rhs[:, index]), np.finfo(float).tiny)
                )
            )
    else:
        backend = "ngsolve-mass-preconditioned-cg"
        system_matrix = demag_matrix + inv_chi * mass_matrix
        preconditioner = mass_matrix.Inverse(inverse=inverse)
        solver = ng.CGSolver(
            system_matrix,
            preconditioner,
            printrates=False,
            precision=solve_tol,
            maxsteps=solve_maxit,
        )

        def apply_operator(values):
            return _ngsolve_apply_to_array(system_matrix, values, mass_matrix)

        for index in range(count):
            rhs_vector = _ngsolve_array_to_vector(mass_matrix, parent_rhs[:, index])
            solved = solver * rhs_vector
            snapshot = np.array(solved.FV().NumPy(), copy=True)
            snapshots.append(snapshot)
            iterations.append(int(solver.GetSteps()))
            residuals.append(
                float(
                    np.linalg.norm(apply_operator(snapshot) - parent_rhs[:, index])
                    / max(np.linalg.norm(parent_rhs[:, index]), np.finfo(float).tiny)
                )
            )
    snapshot_matrix = np.column_stack(snapshots)
    if max(residuals) > max(20.0 * solve_tol, 1.0e-8):
        raise RuntimeError(
            "HDiv-MMM response snapshot solve did not reach the requested residual; "
            f"max relative residual={max(residuals):.3e}"
        )

    energies = np.real(np.sum(snapshot_matrix.conj() * parent_rhs, axis=0))
    energy_floor = np.finfo(float).eps * max(float(np.max(np.abs(energies))), 1.0)
    if np.any(energies <= energy_floor):
        bad = [port_names[index] for index in np.flatnonzero(energies <= energy_floor)]
        raise ValueError(f"response ports have zero or non-positive magnetic energy: {bad}")
    snapshot_actions = _hdiv_response_apply_columns(apply_operator, snapshot_matrix)

    def energy_pod(values, actions, scales):
        if values.shape[1] == 0:
            empty = np.zeros((fes.ndof, 0), dtype=values.dtype)
            return empty, empty.copy(), np.zeros(0), 0
        weighted_values = values * scales[np.newaxis, :]
        weighted_actions = actions * scales[np.newaxis, :]
        correlation = weighted_values.conj().T @ weighted_actions
        correlation = 0.5 * (correlation + correlation.conj().T)
        eigenvalues, eigenvectors = np.linalg.eigh(correlation)
        order = np.argsort(eigenvalues)[::-1]
        eigenvalues = np.real(eigenvalues[order])
        eigenvectors = eigenvectors[:, order]
        threshold = pod_rtol * max(float(eigenvalues[0]), np.finfo(float).tiny)
        rank = int(np.count_nonzero(eigenvalues > threshold))
        if rank == 0:
            empty = np.zeros((fes.ndof, 0), dtype=values.dtype)
            return empty, empty.copy(), eigenvalues, 0
        transform = eigenvectors[:, :rank] @ np.diag(
            1.0 / np.sqrt(eigenvalues[:rank])
        )
        basis = weighted_values @ transform
        basis_actions = weighted_actions @ transform
        basis_gram = basis.conj().T @ basis_actions
        basis_gram = 0.5 * (basis_gram + basis_gram.conj().T)
        try:
            cholesky = np.linalg.cholesky(basis_gram)
        except np.linalg.LinAlgError as exc:
            raise RuntimeError("response basis lost magnetic-energy independence") from exc
        correction = np.linalg.inv(cholesky.conj().T)
        return basis @ correction, basis_actions @ correction, eigenvalues, rank

    physical_count = len(physical_fields)
    if physical_count:
        physical_scales = weights[:physical_count] / (
            np.sqrt(energies[:physical_count])
            if normalize_ports
            else np.ones(physical_count)
        )
        protected_vectors, protected_actions, protected_eigenvalues, protected_rank = (
            energy_pod(
                snapshot_matrix[:, :physical_count],
                snapshot_actions[:, :physical_count],
                physical_scales,
            )
        )
    else:
        protected_vectors = np.zeros((fes.ndof, 0))
        protected_actions = np.zeros((fes.ndof, 0))
        protected_eigenvalues = np.zeros(0)
        protected_rank = 0

    training_values = snapshot_matrix[:, physical_count:]
    training_actions = snapshot_actions[:, physical_count:]
    if protected_rank and training_values.shape[1]:
        protected_gram = protected_vectors.conj().T @ protected_actions
        projection = np.linalg.solve(
            protected_gram,
            protected_vectors.conj().T @ training_actions,
        )
        training_values = training_values - protected_vectors @ projection
        training_actions = training_actions - protected_actions @ projection
    training_energies = np.real(
        np.sum(training_values.conj() * training_actions, axis=0)
    )
    if training_energies.size:
        training_floor = np.finfo(float).eps * max(
            float(np.max(np.abs(energies[physical_count:]))),
            1.0,
        )
        active_training = training_energies > training_floor
    else:
        active_training = np.zeros(0, dtype=bool)
    dependent_training_names = [
        enrichment_names[index]
        for index in np.flatnonzero(~active_training)
    ]
    active_training_values = training_values[:, active_training]
    active_training_actions = training_actions[:, active_training]
    active_training_energies = training_energies[active_training]
    active_training_weights = weights[physical_count:][active_training]
    training_scales = active_training_weights / (
        np.sqrt(active_training_energies)
        if normalize_ports
        else np.ones(active_training_energies.size)
    )
    training_vectors, training_basis_actions, training_eigenvalues, training_rank = (
        energy_pod(
            active_training_values,
            active_training_actions,
            training_scales,
        )
    )

    available_vectors = np.column_stack((protected_vectors, training_vectors))
    available_actions = np.column_stack((protected_actions, training_basis_actions))
    available = int(available_vectors.shape[1])
    if available < 1:
        raise RuntimeError("response POD retained no HDiv modes")
    available_gram = available_vectors.conj().T @ available_actions
    available_gram = 0.5 * (available_gram + available_gram.conj().T)
    try:
        gram_cholesky = np.linalg.cholesky(available_gram)
    except np.linalg.LinAlgError as exc:
        raise RuntimeError("response basis lost magnetic-energy independence") from exc
    correction = np.linalg.inv(gram_cholesky.conj().T)
    available_vectors = available_vectors @ correction
    available_actions = available_actions @ correction
    available_gram = available_vectors.conj().T @ available_actions
    available_gram = 0.5 * (available_gram + available_gram.conj().T)
    retained = available if max_modes is None else min(available, int(max_modes))
    if retained < protected_rank:
        raise ValueError(
            f"max_modes={retained} cannot discard {protected_rank} protected physical-response modes"
        )
    vectors = available_vectors[:, :retained]
    actions = available_actions[:, :retained]
    gram = vectors.conj().T @ actions
    energy_orthonormality_error = float(
        np.linalg.norm(gram - np.eye(retained)) / max(np.sqrt(retained), 1.0)
    )

    protected_response_error = None
    if physical_count and protected_rank:
        protected_basis = available_vectors[:, :protected_rank]
        protected_basis_gram = available_gram[:protected_rank, :protected_rank]
        protected_errors = []
        for index in range(physical_count):
            coefficients = np.linalg.solve(
                protected_basis_gram,
                protected_basis.conj().T @ parent_rhs[:, index],
            )
            error = snapshot_matrix[:, index] - protected_basis @ coefficients
            error_energy = max(
                float(np.real(np.vdot(error, apply_operator(error)))),
                0.0,
            )
            protected_errors.append(float(np.sqrt(error_energy / energies[index])))
        protected_response_error = float(max(protected_errors))

    truncation_curve = []
    for rank in range(1, available + 1):
        rank_vectors = available_vectors[:, :rank]
        rank_gram = available_gram[:rank, :rank]
        reduced_rhs = rank_vectors.conj().T @ parent_rhs
        reduced_coefficients = np.linalg.solve(rank_gram, reduced_rhs)
        captured = np.real(
            np.sum(reduced_rhs.conj() * reduced_coefficients, axis=0)
        )
        errors = np.sqrt(np.maximum(energies - captured, 0.0) / energies)
        row = {
            "rank": rank,
            "max_all_response_relative_energy_error": float(np.max(errors)),
            "rms_all_response_relative_energy_error": float(
                np.sqrt(np.mean(errors * errors))
            ),
            "max_physical_response_relative_energy_error": (
                None
                if not physical_fields
                else (
                    protected_response_error
                    if rank >= protected_rank
                    else float(np.max(errors[:physical_count]))
                )
            ),
        }
        truncation_curve.append(row)
    response_errors = []
    for index in range(count):
        coefficients = np.linalg.solve(
            gram,
            vectors.conj().T @ parent_rhs[:, index],
        )
        error = snapshot_matrix[:, index] - vectors @ coefficients
        error_energy = max(float(np.real(np.vdot(error, apply_operator(error)))), 0.0)
        response_errors.append(float(np.sqrt(error_energy / energies[index])))
    retained_curve = truncation_curve[retained - 1]
    retained_curve["max_all_response_relative_energy_error"] = float(
        max(response_errors)
    )
    retained_curve["rms_all_response_relative_energy_error"] = float(
        np.sqrt(np.mean(np.square(response_errors)))
    )
    retained_curve["max_physical_response_relative_energy_error"] = (
        None
        if not physical_fields
        else float(max(response_errors[: len(physical_fields)]))
    )
    generation = {
        "construction": "hdiv-mmm-response-energy-pod",
        "snapshot_backend": backend,
        "snapshot_port_count": count,
        "physical_rhs_columns": len(physical_fields),
        "training_port_count": len(enrichment_fields),
        "port_names": list(port_names),
        "normalize_ports": bool(normalize_ports),
        "port_weights": weights.tolist(),
        "snapshot_iterations": iterations,
        "snapshot_relative_residuals": residuals,
        "max_snapshot_relative_residual": float(max(residuals)),
        "pod_rtol": pod_rtol,
        "protected_physical_modes": protected_rank,
        "protected_eigenvalues": protected_eigenvalues.tolist(),
        "training_response_modes": training_rank,
        "pod_eigenvalues": training_eigenvalues.tolist(),
        "available_modes": available,
        "retained_modes": retained,
        "discarded_modes": int(count - retained),
        "dependent_port_directions": int(count - available),
        "dependent_training_ports": dependent_training_names,
        "pod_truncated_modes": int(available - retained),
        "pod_truncation_curve": truncation_curve,
        "energy_orthonormality_error": energy_orthonormality_error,
        "response_relative_energy_errors": response_errors,
        "max_response_relative_energy_error": float(max(response_errors)),
        "multipole_training": (
            training_fields.diagnostics()
            if isinstance(training_fields, HDivMultipolePortSet)
            else None
        ),
    }
    return NgsolveHDivMMMReduction(
        mesh,
        fes,
        vectors,
        mu_r=mu_r,
        external_fields=physical_fields or None,
        intorder=intorder,
        materials=materials,
        names=tuple(f"M_response_{index}" for index in range(retained)),
        mass=mass,
        demag_operator=demag_operator,
        basis_generation=generation,
        parent_family=parent_family,
        parent_order=parent_order,
    )


def NgsolveBDMHDivMMMResponseReduction(
    mesh,
    *,
    order: int = 1,
    mu_r: float,
    external_fields=None,
    training_fields=None,
    external_names=None,
    training_names=None,
    port_weights=None,
    normalize_ports: bool = True,
    pod_rtol: float = 1.0e-10,
    max_modes=None,
    solve_tol: float = 1.0e-10,
    solve_maxit: int = 5000,
    inverse: str = "sparsecholesky",
    intorder: int = 2,
    materials=None,
    mass=None,
    demag_operator=None,
    demag_intorder=None,
    demag_eps: float = 1.0e-7,
    demag_leafsize: int = 16,
    demag_eta: float = 2.0,
    demag_far_quad: int = 3,
    demag_ho_far_factor: float = 2.0,
) -> HDivMMMReducedModel:
    """Build the production BDM-MMM response reduction on ``mesh``.

    The parent space is deliberately constructed as bare
    ``ngsolve.HDiv(mesh, order=order)``.  In NGSolve this selects BDM on
    simplex cells; no ``RT=True`` flag is forwarded.  The returned diagnostics
    lock both the family and parent order so a later mixed solve cannot silently
    relabel an RT comparison as production BDM.
    """

    import ngsolve as ng

    order = int(order)
    if order < 1:
        raise ValueError("order must be positive")
    if mesh.dim != 3:
        raise ValueError("NgsolveBDMHDivMMMResponseReduction requires a 3-D mesh")
    fes = ng.HDiv(mesh, order=order)
    return NgsolveHDivMMMResponseReduction(
        mesh,
        fes,
        mu_r=mu_r,
        external_fields=external_fields,
        training_fields=training_fields,
        external_names=external_names,
        training_names=training_names,
        port_weights=port_weights,
        normalize_ports=normalize_ports,
        pod_rtol=pod_rtol,
        max_modes=max_modes,
        solve_tol=solve_tol,
        solve_maxit=solve_maxit,
        inverse=inverse,
        intorder=intorder,
        materials=materials,
        mass=mass,
        demag_operator=demag_operator,
        demag_intorder=demag_intorder,
        demag_eps=demag_eps,
        demag_leafsize=demag_leafsize,
        demag_eta=demag_eta,
        demag_far_quad=demag_far_quad,
        demag_ho_far_factor=demag_ho_far_factor,
        parent_family="BDM",
        parent_order=order,
    )


def _ngsolve_planar_field_tuple(fields):
    import ngsolve as ng

    if fields is None:
        return ()
    if isinstance(fields, PlanarHarmonicPortSet):
        return fields.fields
    if isinstance(fields, (tuple, list)):
        if len(fields) == 2 and all(np.isscalar(value) for value in fields):
            return (ng.CoefficientFunction(tuple(fields)),)
        return tuple(fields)
    return (fields,)


def _planar_response_field_names(fields, supplied_names, prefix: str) -> tuple[str, ...]:
    values = _ngsolve_planar_field_tuple(fields)
    if supplied_names is None:
        if isinstance(fields, PlanarHarmonicPortSet):
            return fields.names
        return tuple(f"{prefix}_{index}" for index in range(len(values)))
    names = tuple(str(name) for name in supplied_names)
    if len(names) != len(values):
        raise ValueError(f"{prefix}_names must match the number of fields")
    return names


@dataclass(frozen=True)
class PlanarHDivMMMReducedSolution:
    """Reduced and parent-space magnetization coefficients for a 2-D solve."""

    reduced_coefficients: np.ndarray
    parent_coefficients: np.ndarray
    residual_relative_norm: float

    @property
    def n_excitations(self) -> int:
        return int(self.reduced_coefficients.shape[1])

    def diagnostics(self) -> dict[str, object]:
        return {
            "n_excitations": self.n_excitations,
            "reduced_modes": int(self.reduced_coefficients.shape[0]),
            "parent_dofs": int(self.parent_coefficients.shape[0]),
            "residual_relative_norm": float(self.residual_relative_norm),
        }


@dataclass(frozen=True)
class PlanarHDivMMMReducedModel:
    """Protected response-POD reduction of a production PlanarDemagBody."""

    mesh: object
    body: object
    parent_vectors: np.ndarray
    mass: np.ndarray
    demag: np.ndarray
    magnetic_operator: np.ndarray
    magnetic_rhs: np.ndarray | None
    mu_r: float
    basis_generation: dict[str, object]

    def __post_init__(self) -> None:
        vectors = _hdiv_reduction_vectors(self.body.fes, self.parent_vectors)
        modes = vectors.shape[1]
        mass = _square_matrix(self.mass, modes, "mass")
        demag = _square_matrix(self.demag, modes, "demag")
        operator = _square_matrix(self.magnetic_operator, modes, "magnetic_operator")
        rhs = self.magnetic_rhs
        if rhs is not None:
            rhs = _port_rhs_matrix(rhs, modes)
        if not np.isfinite(self.mu_r) or self.mu_r <= 1.0:
            raise ValueError("mu_r must be greater than 1")
        object.__setattr__(self, "parent_vectors", vectors)
        object.__setattr__(self, "mass", mass)
        object.__setattr__(self, "demag", demag)
        object.__setattr__(self, "magnetic_operator", operator)
        object.__setattr__(self, "magnetic_rhs", rhs)
        object.__setattr__(self, "basis_generation", dict(self.basis_generation))

    @property
    def n_modes(self) -> int:
        return int(self.parent_vectors.shape[1])

    @property
    def parent_ndof(self) -> int:
        return int(self.parent_vectors.shape[0])

    @property
    def chi(self) -> float:
        return float(self.mu_r - 1.0)

    def reconstruct_parent(self, reduced_coefficients) -> np.ndarray:
        coefficients = np.asarray(reduced_coefficients)
        if coefficients.ndim == 1:
            coefficients = coefficients[:, np.newaxis]
        if coefficients.ndim != 2 or coefficients.shape[0] != self.n_modes:
            raise ValueError(
                f"reduced_coefficients must have shape ({self.n_modes}, n_rhs)"
            )
        return self.parent_vectors @ coefficients

    def external_field_rhs(self, external_fields) -> np.ndarray:
        fields = _ngsolve_planar_field_tuple(external_fields)
        if not fields:
            raise ValueError("external_fields must not be empty")
        projected = np.column_stack([self.body.project(field) for field in fields])
        parent_rhs = np.column_stack(
            [np.asarray(self.body.Mm @ projected[:, index]).reshape(-1)
             for index in range(projected.shape[1])]
        )
        return self.parent_vectors.conj().T @ parent_rhs

    def solve(self, external_fields=None) -> PlanarHDivMMMReducedSolution:
        rhs = (
            self.magnetic_rhs
            if external_fields is None
            else self.external_field_rhs(external_fields)
        )
        if rhs is None:
            raise ValueError("no stored magnetic RHS; pass external_fields")
        coefficients = np.linalg.solve(self.magnetic_operator, rhs)
        residual = self.magnetic_operator @ coefficients - rhs
        relative = float(
            np.linalg.norm(residual)
            / max(np.linalg.norm(rhs), np.finfo(float).tiny)
        )
        return PlanarHDivMMMReducedSolution(
            reduced_coefficients=coefficients,
            parent_coefficients=self.reconstruct_parent(coefficients),
            residual_relative_norm=relative,
        )

    def diagnostics(self) -> dict[str, object]:
        return {
            "parent_space": "PlanarHDiv",
            "parent_family": (
                "Raviart-Thomas" if getattr(self.body, "rt", False) else "BDM"
            ),
            "parent_order": int(self.body.order),
            "parent_ndof": self.parent_ndof,
            "reduced_modes": self.n_modes,
            "compression_ratio": float(self.n_modes / self.parent_ndof),
            "mu_r": float(self.mu_r),
            "chi": self.chi,
            "mass_hermitian_error": _relative_hermitian_error(self.mass),
            "demag_hermitian_error": _relative_hermitian_error(self.demag),
            "operator_hermitian_error": _relative_hermitian_error(
                self.magnetic_operator
            ),
            "min_operator_eigenvalue": _min_hermitian_eigenvalue(
                self.magnetic_operator
            ),
            "has_rhs": self.magnetic_rhs is not None,
            "rhs_columns": (
                0 if self.magnetic_rhs is None else int(self.magnetic_rhs.shape[1])
            ),
            "basis_generation": self.basis_generation,
        }


def NgsolvePlanarHDivMMMResponseReduction(
    mesh,
    *,
    mu_r: float,
    external_fields=None,
    training_fields=None,
    external_names=None,
    training_names=None,
    port_weights=None,
    normalize_ports: bool = True,
    pod_rtol: float = 1.0e-10,
    max_modes=None,
    body=None,
    order: int = 1,
    rt: bool = False,
    eta: float = 2.0,
    cg_tol: float = 1.0e-10,
    cg_maxit: int = 5000,
) -> PlanarHDivMMMReducedModel:
    """Build the 2-D counterpart of the protected 3-D HDiv response basis."""

    if mesh.dim != 2:
        raise ValueError("NgsolvePlanarHDivMMMResponseReduction requires a 2-D mesh")
    mu_r = float(mu_r)
    if not np.isfinite(mu_r) or mu_r <= 1.0:
        raise ValueError("mu_r must be greater than 1")
    physical_fields = _ngsolve_planar_field_tuple(external_fields)
    enrichment_fields = _ngsolve_planar_field_tuple(training_fields)
    if not physical_fields and not enrichment_fields:
        raise ValueError("at least one external or training field is required")
    physical_names = _planar_response_field_names(
        external_fields,
        external_names,
        "external",
    )
    enrichment_names = _planar_response_field_names(
        training_fields,
        training_names,
        "training",
    )
    if body is None:
        from ._vim2d import PlanarDemagBody

        body = PlanarDemagBody(
            mesh,
            order=order,
            rt=rt,
            eta=eta,
            cg_tol=cg_tol,
            cg_maxit=cg_maxit,
        )
    elif body.mesh is not mesh:
        raise ValueError("body must have been built from mesh")
    if (
        isinstance(training_fields, PlanarHarmonicPortSet)
        and body.order < training_fields.max_degree - 1
    ):
        raise ValueError(
            "Planar HDiv order must be at least harmonic max_degree - 1 "
            "for polynomially admissible training"
        )

    fields = physical_fields + enrichment_fields
    port_names = physical_names + enrichment_names
    projected = np.column_stack([body.project(field) for field in fields])
    parent_rhs = np.column_stack(
        [np.asarray(body.Mm @ projected[:, index]).reshape(-1)
         for index in range(projected.shape[1])]
    )
    snapshots = []
    iterations = []
    residuals = []
    inv_chi = 1.0 / (mu_r - 1.0)

    def apply_operator(values):
        return np.asarray(body.apply_demag(values)) + inv_chi * np.asarray(
            body.Mm @ values
        ).reshape(-1)

    for index in range(projected.shape[1]):
        snapshot = np.asarray(body.solve_linear(mu_r - 1.0, projected[:, index]))
        snapshots.append(snapshot)
        iterations.append(int(body.last_linear_iterations))
        residuals.append(
            float(
                np.linalg.norm(apply_operator(snapshot) - parent_rhs[:, index])
                / max(np.linalg.norm(parent_rhs[:, index]), np.finfo(float).tiny)
            )
        )
    snapshot_matrix = np.column_stack(snapshots)
    if max(residuals) > max(20.0 * float(body.cg_tol), 1.0e-8):
        raise RuntimeError(
            "Planar HDiv-MMM response snapshot solve did not converge; "
            f"max relative residual={max(residuals):.3e}"
        )
    vectors, pod_diagnostics = _protected_magnetic_response_pod(
        snapshot_matrix,
        parent_rhs,
        apply_operator,
        physical_count=len(physical_fields),
        port_names=port_names,
        training_names=enrichment_names,
        port_weights=port_weights,
        normalize_ports=normalize_ports,
        pod_rtol=pod_rtol,
        max_modes=max_modes,
    )
    mass_actions = np.asarray(body.Mm @ vectors)
    demag_actions = np.column_stack(
        [np.asarray(body.apply_demag(vectors[:, index]))
         for index in range(vectors.shape[1])]
    )
    reduced_mass = vectors.conj().T @ mass_actions
    reduced_demag = vectors.conj().T @ demag_actions
    reduced_mass = 0.5 * (reduced_mass + reduced_mass.conj().T)
    reduced_demag = 0.5 * (reduced_demag + reduced_demag.conj().T)
    reduced_operator = reduced_mass * inv_chi + reduced_demag
    reduced_rhs = (
        None
        if not physical_fields
        else vectors.conj().T @ parent_rhs[:, : len(physical_fields)]
    )
    generation = {
        "construction": "planar-hdiv-mmm-response-energy-pod",
        "snapshot_backend": "radia-cpp-mass-riesz-cg-2d",
        "snapshot_port_count": len(fields),
        "physical_rhs_columns": len(physical_fields),
        "training_port_count": len(enrichment_fields),
        "port_names": list(port_names),
        "snapshot_iterations": iterations,
        "snapshot_relative_residuals": residuals,
        "max_snapshot_relative_residual": float(max(residuals)),
        "planar_harmonic_training": (
            training_fields.diagnostics()
            if isinstance(training_fields, PlanarHarmonicPortSet)
            else None
        ),
        **pod_diagnostics,
    }
    return PlanarHDivMMMReducedModel(
        mesh=mesh,
        body=body,
        parent_vectors=vectors,
        mass=reduced_mass,
        demag=reduced_demag,
        magnetic_operator=reduced_operator,
        magnetic_rhs=reduced_rhs,
        mu_r=mu_r,
        basis_generation=generation,
    )


def _ngsolve_static_condensed_solve(form, inverse_operator, rhs_values: np.ndarray) -> np.ndarray:
    kmat = _ngsolve_matrix_object(form)
    rhs = _ngsolve_array_to_vector(kmat, rhs_values)
    rhs.data += form.harmonic_extension_trans * rhs
    solution = rhs.CreateVector()
    solution.data = inverse_operator * rhs
    solution.data += form.harmonic_extension * solution
    solution.data += form.inner_solve * rhs
    return np.array(solution.FV().NumPy(), copy=True)


def _append_inner_orthonormal(
    columns: list[np.ndarray],
    candidate: np.ndarray,
    inner,
    *,
    rtol: float,
) -> None:
    v = np.array(candidate, copy=True)
    candidate_norm2 = inner(v, v)
    candidate_norm = float(np.sqrt(max(np.real(candidate_norm2), 0.0)))
    if candidate_norm == 0.0:
        return
    for _ in range(2):
        for q in columns:
            v -= q * inner(q, v)
    norm2 = inner(v, v)
    norm = float(np.sqrt(max(np.real(norm2), 0.0)))
    relative_floor = max(
        rtol,
        10.0 * candidate.size * np.finfo(float).eps,
    )
    if norm > relative_floor * candidate_norm:
        columns.append(v / norm)


def NgsolveOperatorBlockKrylovBasis(
    stiffness,
    mass,
    ports,
    steps: int,
    *,
    free_dofs=None,
    inverse: str = "sparsecholesky",
    rtol: float = 1.0e-12,
) -> EVRSBasis:
    """Operator-based NGSolve block-Krylov basis for higher-order spaces.

    Unlike :func:`NgsolveBlockKrylovBasis`, this helper does not convert the
    NGSolve matrices to dense NumPy arrays.  It uses ``stiffness.mat.Inverse``
    and repeated ``mass.mat`` matvecs, while still returning ordinary dense
    coefficient columns for the reduced basis.  Use it when the HCurl parent
    space is high-order and only the response basis is meant to stay small.
    """

    if steps < 1:
        raise ValueError("steps must be >= 1")
    if rtol <= 0.0:
        raise ValueError("rtol must be positive")
    kmat = _ngsolve_matrix_object(stiffness)
    mmat = _ngsolve_matrix_object(mass)
    n = _matrix_height(kmat)
    b = _ports_to_matrix(ports)
    if b.ndim == 1:
        b = b[:, np.newaxis]
    if b.ndim != 2 or b.shape[0] != n:
        raise ValueError(f"ports must have shape ({n}, p)")

    active = _active_indices(free_dofs, n)
    if active.size == 0:
        raise ValueError("no active dofs")
    active_mask = np.zeros(n, dtype=bool)
    active_mask[active] = True

    def constrained(values):
        out = np.asarray(values).copy()
        out[~active_mask] = 0.0
        return out

    def mass_apply(values):
        return constrained(_ngsolve_apply_to_array(mmat, constrained(values), kmat))

    def metric_inner(a, b_):
        return np.vdot(a, mass_apply(b_))

    if free_dofs is None:
        inv = kmat.Inverse(inverse=inverse)
    else:
        inv = kmat.Inverse(freedofs=free_dofs, inverse=inverse)

    current = np.column_stack([constrained(b[:, i]) for i in range(b.shape[1])])
    columns: list[np.ndarray] = []
    for _ in range(steps):
        solved = np.column_stack(
            [
                constrained(_ngsolve_apply_to_array(inv, current[:, i], kmat))
                for i in range(current.shape[1])
            ]
        )
        for col in range(solved.shape[1]):
            _append_inner_orthonormal(
                columns, solved[:, col], metric_inner, rtol=rtol
            )
        current = np.column_stack([mass_apply(solved[:, i]) for i in range(solved.shape[1])])

    if columns:
        vectors = np.column_stack(columns)
    else:
        vectors = np.zeros((n, 0), dtype=np.result_type(b, float))
    return EVRSBasis(
        vectors=vectors,
        active_dofs=active,
        port_count=int(b.shape[1]),
        krylov_steps=steps,
        construction="ngsolve-operator-block-krylov",
        parent_space="HCurl",
    )


def NgsolveStaticCondensedBlockKrylovBasis(
    stiffness,
    mass,
    ports,
    steps: int,
    *,
    free_dofs=None,
    inverse: str = "sparsecholesky",
    rtol: float = 1.0e-12,
) -> EVRSBasis:
    """NGSolve Schur/static-condensed block-Krylov basis.

    ``stiffness`` must be an assembled ``BilinearForm(..., condense=True)``.
    Its NGSolve static-condensation operators eliminate NGSolve ``LOCAL_DOF``
    element bubbles exactly in the stiffness solve, then reconstruct full HCurl
    coefficients for sampling ``curl(T)``.  Runtime unknowns are the non-local
    ``free_dofs`` (usually ``fes.FreeDofs(coupling=True)``), while the returned
    vectors remain full length so they can be evaluated by NGSolve.
    """

    if not all(
        hasattr(stiffness, name)
        for name in ("harmonic_extension", "harmonic_extension_trans", "inner_solve")
    ):
        raise ValueError("stiffness must be a condensed NGSolve BilinearForm")
    if steps < 1:
        raise ValueError("steps must be >= 1")
    if rtol <= 0.0:
        raise ValueError("rtol must be positive")
    kmat = _ngsolve_matrix_object(stiffness)
    mmat = _ngsolve_matrix_object(mass)
    n = _matrix_height(kmat)
    b = _ports_to_matrix(ports)
    if b.ndim == 1:
        b = b[:, np.newaxis]
    if b.ndim != 2 or b.shape[0] != n:
        raise ValueError(f"ports must have shape ({n}, p)")

    if free_dofs is None:
        space = getattr(stiffness, "space", None)
        if space is None:
            raise ValueError("free_dofs is required when stiffness has no space")
        free_dofs = space.FreeDofs(True)
    active = _active_indices(free_dofs, n)
    if active.size == 0:
        raise ValueError("no active dofs")
    active_mask = np.zeros(n, dtype=bool)
    active_mask[active] = True

    def constrained(values):
        out = np.asarray(values).copy()
        out[~active_mask] = 0.0
        return out

    def mass_apply(values):
        return _ngsolve_apply_to_array(mmat, np.asarray(values), kmat)

    def metric_inner(a, b_):
        return np.vdot(a, mass_apply(b_))

    inv = kmat.Inverse(freedofs=free_dofs, inverse=inverse)
    current = np.column_stack([b[:, i] for i in range(b.shape[1])])
    columns: list[np.ndarray] = []
    for _ in range(steps):
        solved = np.column_stack(
            [
                _ngsolve_static_condensed_solve(stiffness, inv, current[:, i])
                for i in range(current.shape[1])
            ]
        )
        for col in range(solved.shape[1]):
            _append_inner_orthonormal(
                columns, solved[:, col], metric_inner, rtol=rtol
            )
        current = np.column_stack([mass_apply(solved[:, i]) for i in range(solved.shape[1])])

    if columns:
        vectors = np.column_stack(columns)
    else:
        vectors = np.zeros((n, 0), dtype=np.result_type(b, float))
    return EVRSBasis(
        vectors=vectors,
        active_dofs=active,
        port_count=int(b.shape[1]),
        krylov_steps=steps,
        construction="ngsolve-static-condensed-block-krylov",
        parent_space="HCurl",
    )


@dataclass(frozen=True)
class EddyBubbleHCurlBasis:
    """Production HCurl eddy-current basis after eddy-bubble elimination.

    The object is the reusable basis-function artifact behind the current
    research path:

    ``high-order HCurl parent T -> EVRS coefficients -> J = curl(T) samples``.

    ``eddy_bubbling`` records which parent DoFs were structural keeps
    (SIBC/non-SIBC trace/bridge) and which were ordinary bulk eddy-bubble
    candidates.  ``current_basis`` is the small VIM/MMM-facing current basis
    that survives the reduction.
    """

    response_basis: EVRSBasis
    current_basis: SampledCurrentBasis
    eddy_bubbling: EddyBubbleDecomposition
    material_model: object | None = None
    cell_families: HCurlCellFamilyInventory | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.response_basis, EVRSBasis):
            raise TypeError("response_basis must be an EVRSBasis")
        if not isinstance(self.current_basis, SampledCurrentBasis):
            raise TypeError("current_basis must be a SampledCurrentBasis")
        if not isinstance(self.eddy_bubbling, EddyBubbleDecomposition):
            raise TypeError("eddy_bubbling must be an EddyBubbleDecomposition")
        if self.cell_families is not None and not isinstance(
            self.cell_families, HCurlCellFamilyInventory
        ):
            raise TypeError("cell_families must be an HCurlCellFamilyInventory")
        if self.current_basis.n_modes != self.response_basis.rank:
            raise ValueError("current_basis modes must match response_basis rank")

    @property
    def n_modes(self) -> int:
        return self.current_basis.n_modes

    @property
    def rank(self) -> int:
        return self.response_basis.rank

    @property
    def active_dofs(self) -> int:
        return self.response_basis.active_count

    @property
    def eddy_basis(self) -> SampledCurrentBasis:
        """Alias used by HDiv-MMM coupling code."""

        return self.current_basis

    @property
    def volume_basis(self) -> SampledCurrentBasis:
        """Alias used by VIM assembly code."""

        return self.current_basis

    def assemble_vim(
        self,
        *,
        sigma: float,
        kernel_epsilon: float | None = None,
        interaction=None,
    ) -> "HybridVIMSystem":
        """Assemble the reduced HCurl-VIM eddy impedance on this basis."""

        if (
            self.cell_families is not None
            and self.cell_families.dimension == 2
            and interaction is None
        ):
            raise ValueError(
                "2-D HCurl Eddy Bubble requires an explicit planar-log interaction; "
                "the default sampled 3-D Laplace kernel is not physically admissible"
            )

        return AssembleHybridVIM(
            self.current_basis,
            sigma=sigma,
            kernel_epsilon=kernel_epsilon,
            interaction=interaction,
        )

    def tet_volume_interaction(self, mesh, fes, **kwargs):
        """Build the epsilon-free affine-tet inductance backend for this basis.

        The response vectors are the same compressed vectors that generated
        ``current_basis``, so the analytic Gram and the sampled resistance/port
        blocks share one reduced coordinate system.
        """

        from ._hcurl_tet_interaction import NgsolveHCurlTetVolumeInteraction

        return NgsolveHCurlTetVolumeInteraction(
            mesh,
            fes,
            self.response_basis.vectors,
            self.current_basis,
            **kwargs,
        )

    def cell_volume_interaction(self, mesh, fes, **kwargs):
        """Build the epsilon-free 3-D affine, curved, or hp cell interaction."""

        from ._hcurl_tet_interaction import NgsolveHCurlCellVolumeInteraction

        return NgsolveHCurlCellVolumeInteraction(
            mesh,
            fes,
            self.response_basis.vectors,
            self.current_basis,
            **kwargs,
        )

    def planar_volume_interaction(self, mesh, fes, **kwargs):
        """Build the epsilon-free planar log interaction for TRIG/QUAD."""

        from ._hcurl_planar_interaction import NgsolveHCurlPlanarVolumeInteraction

        return NgsolveHCurlPlanarVolumeInteraction(
            mesh,
            fes,
            self.response_basis.vectors,
            self.current_basis,
            **kwargs,
        )

    def couple_hdiv_mmm(
        self,
        magnetization_basis: SampledMagnetizationBasis,
        *,
        eddy_system: "HybridVIMSystem | None" = None,
        sigma: float | None = None,
        material_model: object | None = None,
        mu: float = MU0,
        kernel_epsilon: float = 0.0,
    ) -> "HCurlVIMHDivMMMSystem":
        """Couple this reduced eddy basis to an HDiv-MMM magnetization basis."""

        if eddy_system is None and sigma is not None:
            eddy_system = self.assemble_vim(sigma=sigma, kernel_epsilon=kernel_epsilon)
        if material_model is None:
            material_model = self.material_model
        return CoupleHCurlVIMWithHDivMMM(
            magnetization_basis,
            self.current_basis,
            eddy_system=eddy_system,
            material_model=material_model,
            mu=mu,
            kernel_epsilon=kernel_epsilon,
        )

    def diagnostics(self) -> dict[str, object]:
        """Return basis, reduction, and material diagnostics for production logs."""

        info: dict[str, object] = {
            "kind": "EddyBubbleHCurlBasis",
            "modes": self.n_modes,
            "rank": self.rank,
            "active_dofs": self.active_dofs,
            "response_basis": self.response_basis.diagnostics(),
            "current_basis": _sampled_current_basis_diagnostics(self.current_basis),
            "eddy_bubbling": self.eddy_bubbling.diagnostics(),
            "cell_families": (
                self.cell_families.diagnostics()
                if self.cell_families is not None
                else None
            ),
            "has_material_model": self.material_model is not None,
            "has_shared_mesh_material_model": isinstance(
                self.material_model,
                SharedMeshMaterialModel,
            ),
        }
        if isinstance(self.material_model, SharedMeshMaterialModel):
            info["material_model"] = self.material_model.diagnostics()
        return info


def _ngsolve_response_basis_for_eddy_bubbling(
    stiffness,
    mass,
    ports,
    *,
    steps: int,
    free_dofs,
    condense: bool,
    response_backend: str,
    inverse: str,
    rtol: float,
) -> EVRSBasis:
    backend = response_backend
    if backend == "auto":
        backend = "static-condensed" if condense else "operator"
    if backend == "static-condensed":
        return NgsolveStaticCondensedBlockKrylovBasis(
            stiffness,
            mass,
            ports,
            steps=steps,
            free_dofs=free_dofs,
            inverse=inverse,
            rtol=rtol,
        )
    if backend == "operator":
        return NgsolveOperatorBlockKrylovBasis(
            stiffness,
            mass,
            ports,
            steps=steps,
            free_dofs=free_dofs,
            inverse=inverse,
            rtol=rtol,
        )
    if backend == "dense":
        return NgsolveBlockKrylovBasis(
            stiffness,
            mass,
            ports,
            steps=steps,
            free_dofs=free_dofs,
            rtol=rtol,
        )
    raise ValueError("response_backend must be 'auto', 'operator', 'static-condensed', or 'dense'")


def NgsolveEddyBubbleHCurlBasis(
    mesh,
    fes,
    stiffness,
    mass,
    ports,
    *,
    steps: int,
    conductive_materials,
    air_materials=("air", "vacuum"),
    volume_materials=None,
    topology: EddyMeshTopology | None = None,
    intorder: int = 2,
    free_dofs=None,
    condense: bool = False,
    response_backend: str = "auto",
    inverse: str = "sparsecholesky",
    rtol: float = 1.0e-10,
    current_gram_rtol: float = 1.0e-10,
    names=None,
    include_edge_dofs: bool = True,
    parent_order: int | None = None,
    parent_order_ledger: EddyParentOrderLedger | None = None,
    surface_modes: int = 0,
    non_sibc_trace_modes: int | None = None,
    loop_bridge_modes: int | None = None,
    bridge_strategy: str | None = None,
    material_model: object | None = None,
) -> EddyBubbleHCurlBasis:
    """Build the production eddy-bubbled HCurl basis from an NGSolve parent.

    This is the high-level basis constructor used before assembling VIM or
    coupling to HDiv-MMM.  It generates an EVRS response basis in the high-order
    parent ``fes``, samples the physical current basis ``J = curl(T)``, and
    records the topology-aware eddy-bubbling split.
    """

    if steps < 1:
        raise ValueError("steps must be >= 1")
    if topology is None:
        topology = ClassifyNgsolveEddyTopology(
            mesh,
            conductive_materials,
            air_materials=air_materials,
        )
    if free_dofs is None:
        free_dofs = fes.FreeDofs(bool(condense))
    if volume_materials is None:
        volume_materials = conductive_materials
    cell_families = NgsolveHCurlCellFamilies(
        mesh,
        materials=volume_materials,
    )

    response = _ngsolve_response_basis_for_eddy_bubbling(
        stiffness,
        mass,
        ports,
        steps=steps,
        free_dofs=free_dofs,
        condense=condense,
        response_backend=response_backend,
        inverse=inverse,
        rtol=rtol,
    )
    if names is None:
        p = _ngsolve_fes_order(fes) if parent_order is None else parent_order
        prefix = f"eddy_p{p}_n{steps}" if p is not None else f"eddy_n{steps}"
        names = [f"{prefix}_{i}" for i in range(response.rank)]
    current = NgsolveHCurlCurlBasis(
        mesh,
        fes,
        response.vectors,
        intorder=intorder,
        materials=volume_materials,
        names=names,
    )
    response, current = CompressHCurlResponseInCurrentGram(
        response,
        current,
        rtol=current_gram_rtol,
    )
    graph = topology.conductor_graph()
    if loop_bridge_modes is None:
        loop_bridge_modes = graph.cycle_rank
    if bridge_strategy is None:
        bridge_strategy = "cycle-basis"
    bubbling = NgsolveEddyBubbleReduction(
        mesh,
        fes,
        topology,
        free_dofs=free_dofs,
        include_edge_dofs=include_edge_dofs,
        evrs_rank=response.rank,
        surface_modes=surface_modes,
        non_sibc_trace_modes=non_sibc_trace_modes,
        loop_bridge_modes=loop_bridge_modes,
        bridge_strategy=bridge_strategy,
        parent_order=parent_order,
        parent_order_ledger=parent_order_ledger,
    )
    return EddyBubbleHCurlBasis(
        response_basis=response,
        current_basis=current,
        eddy_bubbling=bubbling,
        material_model=material_model,
        cell_families=cell_families,
    )


def _default_kernel_epsilon(bases: tuple[SampledCurrentBasis, ...]) -> float:
    points = np.vstack([basis.points for basis in bases])
    if points.shape[0] < 2:
        return 1.0
    diff = points[:, np.newaxis, :] - points[np.newaxis, :, :]
    dist = np.linalg.norm(diff, axis=2)
    positive = dist[dist > 0.0]
    if positive.size:
        return 0.25 * float(np.min(positive))
    return 1.0


def _interaction_block(
    left: SampledCurrentBasis,
    right: SampledCurrentBasis,
    *,
    mu: float,
    kernel_epsilon: float,
) -> np.ndarray:
    diff = left.points[:, np.newaxis, :] - right.points[np.newaxis, :, :]
    dist = np.sqrt(np.einsum("ijk,ijk->ij", diff, diff) + kernel_epsilon**2)
    kernel = (mu / (4.0 * np.pi)) * (
        left.weights[:, np.newaxis] * right.weights[np.newaxis, :]
    ) / dist
    return np.einsum(
        "aik,bjk,ij->ab", left.modes.conj(), right.modes, kernel
    )


def _planar_interaction_block(
    left: SampledCurrentBasis,
    right: SampledCurrentBasis,
    *,
    mu: float,
    kernel_epsilon: float,
    reference_length: float,
) -> np.ndarray:
    diff = left.points[:, np.newaxis, :2] - right.points[np.newaxis, :, :2]
    distance = np.sqrt(np.einsum("ijk,ijk->ij", diff, diff) + kernel_epsilon**2)
    kernel = -(mu / (2.0 * np.pi)) * (
        left.weights[:, np.newaxis] * right.weights[np.newaxis, :]
    ) * np.log(distance / reference_length)
    return np.einsum(
        "aik,bjk,ij->ab", left.modes.conj(), right.modes, kernel
    )


@dataclass(frozen=True)
class SampledLaplaceInteraction:
    """Callable sampled Laplace single-layer interaction backend.

    This is the dense quadrature reference backend made explicit.  Production
    high-level builders select a HACApK-compressed interaction unless the
    caller deliberately supplies this backend for a small-ROM comparison.
    """

    mu: float = MU0
    kernel_epsilon: float = 1.0
    name: str = "sampled-laplace"

    def __post_init__(self) -> None:
        if self.mu <= 0.0:
            raise ValueError("mu must be positive")
        if self.kernel_epsilon <= 0.0:
            raise ValueError("kernel_epsilon must be positive")
        if not self.name:
            raise ValueError("name must not be empty")

    def __call__(self, left: SampledCurrentBasis, right: SampledCurrentBasis) -> np.ndarray:
        return _interaction_block(
            left,
            right,
            mu=self.mu,
            kernel_epsilon=self.kernel_epsilon,
        )


@dataclass(frozen=True)
class HACApKSampledLaplaceInteraction:
    """Build one projected HACApK operator for hybrid current bases.

    The scalar H-matrix stores the stable full sampled Laplace Gram.  In the
    default ``cross_only`` composition mode, its small reduced diagonal blocks
    are subtracted and replaced by ``diagonal_interaction``.  Thus every cross
    block remains in one HACApK apply without imposing a partition-zero pattern
    on ACA, while exact high-order HCurl diagonal operators are not double
    counted.  Set ``cross_only=False`` to use the sampled HACApK diagonal too.

    ``diagonal_interaction`` may be a callable backend for the diagonal basis
    pairs.  It is useful for retaining an analytic HCurl volume operator while
    the volume/bridge/SIBC coupling is supplied by this class.  If omitted,
    diagonal blocks use the sampled Laplace reference quadrature.
    """

    mu: float = MU0
    kernel_epsilon: float = 1.0
    aca_eps: float = 1.0e-11
    leaf_size: int = 64
    eta: float = 2.0
    cross_only: bool = True
    diagonal_interaction: object | None = None
    diagonal_bases: tuple[SampledCurrentBasis, ...] = ()
    name: str = "hacapk-sampled-laplace"

    def __post_init__(self) -> None:
        positive = {
            "mu": self.mu,
            "kernel_epsilon": self.kernel_epsilon,
            "aca_eps": self.aca_eps,
            "eta": self.eta,
        }
        for key, value in positive.items():
            number = float(value)
            if not np.isfinite(number) or number <= 0.0:
                raise ValueError(f"{key} must be positive")
            object.__setattr__(self, key, number)
        leaf_size = int(self.leaf_size)
        if leaf_size < 1:
            raise ValueError("leaf_size must be >= 1")
        object.__setattr__(self, "leaf_size", leaf_size)
        if self.diagonal_interaction is not None and not callable(self.diagonal_interaction):
            raise TypeError("diagonal_interaction must be callable")
        diagonal_bases = tuple(self.diagonal_bases)
        if any(not isinstance(basis, SampledCurrentBasis) for basis in diagonal_bases):
            raise TypeError("diagonal_bases must contain SampledCurrentBasis objects")
        if diagonal_bases and self.diagonal_interaction is None:
            raise ValueError("diagonal_bases requires diagonal_interaction")
        object.__setattr__(self, "diagonal_bases", diagonal_bases)
        if not self.name:
            raise ValueError("name must not be empty")

    @property
    def operator_scope(self) -> str:
        return "cross" if self.cross_only else "full"

    def build_operator(self, bases: tuple[SampledCurrentBasis, ...]):
        return _build_sampled_hacapk_operator(
            bases,
            mu=self.mu,
            kernel_epsilon=self.kernel_epsilon,
            aca_eps=self.aca_eps,
            leaf_size=self.leaf_size,
            eta=self.eta,
            cross_only=self.cross_only,
        )

    def __call__(self, left: SampledCurrentBasis, right: SampledCurrentBasis):
        sampled = _interaction_block(
            left,
            right,
            mu=self.mu,
            kernel_epsilon=self.kernel_epsilon,
        )
        if not self.cross_only or left is not right:
            return sampled
        use_override = self.diagonal_interaction is not None and (
            not self.diagonal_bases
            or any(left is basis for basis in self.diagonal_bases)
        )
        if not use_override:
            return np.zeros_like(sampled)
        desired = self.diagonal_interaction(left, right)
        if _is_matrix_free_operator(desired):
            return ReducedBlockHMatrixOperator(
                -sampled,
                [(0, left.n_modes, desired)],
            )
        return np.asarray(desired) - sampled


@dataclass(frozen=True)
class HACApKSampledPlanarLogInteraction:
    """Projected HACApK cross interaction for a 2-D ``-log(r/L)`` kernel."""

    mu: float = MU0
    kernel_epsilon: float = 1.0
    reference_length: float = 1.0
    aca_eps: float = 1.0e-11
    leaf_size: int = 64
    eta: float = 2.0
    cross_only: bool = True
    diagonal_interaction: object | None = None
    diagonal_bases: tuple[SampledCurrentBasis, ...] = ()
    name: str = "hacapk-sampled-planar-log"

    def __post_init__(self) -> None:
        for key in (
            "mu", "kernel_epsilon", "reference_length", "aca_eps", "eta"
        ):
            number = float(getattr(self, key))
            if not np.isfinite(number) or number <= 0.0:
                raise ValueError(f"{key} must be positive")
            object.__setattr__(self, key, number)
        leaf_size = int(self.leaf_size)
        if leaf_size < 1:
            raise ValueError("leaf_size must be >= 1")
        object.__setattr__(self, "leaf_size", leaf_size)
        if self.diagonal_interaction is not None and not callable(
            self.diagonal_interaction
        ):
            raise TypeError("diagonal_interaction must be callable")
        diagonal_bases = tuple(self.diagonal_bases)
        if any(
            not isinstance(basis, SampledCurrentBasis)
            for basis in diagonal_bases
        ):
            raise TypeError(
                "diagonal_bases must contain SampledCurrentBasis objects"
            )
        if diagonal_bases and self.diagonal_interaction is None:
            raise ValueError("diagonal_bases requires diagonal_interaction")
        object.__setattr__(self, "diagonal_bases", diagonal_bases)
        if not self.name:
            raise ValueError("name must not be empty")

    @property
    def operator_scope(self) -> str:
        return "cross" if self.cross_only else "full"

    def build_operator(self, bases: tuple[SampledCurrentBasis, ...]):
        return _build_sampled_hacapk_operator(
            bases,
            mu=self.mu,
            kernel_epsilon=self.kernel_epsilon,
            aca_eps=self.aca_eps,
            leaf_size=self.leaf_size,
            eta=self.eta,
            cross_only=self.cross_only,
            kernel="planar-log",
            reference_length=self.reference_length,
        )

    def __call__(self, left: SampledCurrentBasis, right: SampledCurrentBasis):
        sampled = _planar_interaction_block(
            left,
            right,
            mu=self.mu,
            kernel_epsilon=self.kernel_epsilon,
            reference_length=self.reference_length,
        )
        if not self.cross_only or left is not right:
            return sampled
        use_override = self.diagonal_interaction is not None and (
            not self.diagonal_bases
            or any(left is basis for basis in self.diagonal_bases)
        )
        if not use_override:
            return np.zeros_like(sampled)
        desired = self.diagonal_interaction(left, right)
        if _is_matrix_free_operator(desired):
            return ReducedBlockHMatrixOperator(
                -sampled,
                [(0, left.n_modes, desired)],
            )
        return np.asarray(desired) - sampled


def _basis_offsets(bases: tuple[SampledCurrentBasis, ...]) -> list[tuple[int, int]]:
    offsets: list[tuple[int, int]] = []
    start = 0
    for basis in bases:
        stop = start + basis.n_modes
        offsets.append((start, stop))
        start = stop
    return offsets


@dataclass(frozen=True)
class ReducedInteractionMatrix:
    """Callable backend backed by a precomputed reduced interaction matrix.

    Use this when a BEM/H-matrix implementation has already projected the
    Laplace single-layer operator onto the reduced basis.  The object slices the
    full reduced matrix by basis identity, so :func:`AssembleHybridVIM` can
    assemble resistance and SIBC surface terms while reusing the supplied BEM
    inductance blocks.
    """

    bases: tuple[SampledCurrentBasis, ...]
    matrix: np.ndarray
    name: str = "reduced-interaction-matrix"

    def __post_init__(self) -> None:
        bases = tuple(self.bases)
        if not bases:
            raise ValueError("bases must not be empty")
        for basis in bases:
            if not isinstance(basis, SampledCurrentBasis):
                raise TypeError("bases must contain SampledCurrentBasis objects")
        total = sum(basis.n_modes for basis in bases)
        matrix = np.asarray(self.matrix)
        if matrix.shape != (total, total):
            raise ValueError(f"matrix must have shape ({total}, {total})")
        if not np.all(np.isfinite(matrix)):
            raise ValueError("matrix contains non-finite values")
        if not self.name:
            raise ValueError("name must not be empty")
        object.__setattr__(self, "bases", bases)
        object.__setattr__(self, "matrix", np.array(matrix, copy=True))

    def _index(self, basis: SampledCurrentBasis) -> int:
        for i, candidate in enumerate(self.bases):
            if candidate is basis:
                return i
        raise ValueError("basis is not registered in this interaction matrix")

    def __call__(self, left: SampledCurrentBasis, right: SampledCurrentBasis) -> np.ndarray:
        offsets = _basis_offsets(self.bases)
        a0, a1 = offsets[self._index(left)]
        b0, b1 = offsets[self._index(right)]
        return self.matrix[a0:a1, b0:b1]


class NGSolveProjectedOperator:
    """Matrix-free reduced operator backed by an NGSolve ``BaseMatrix``."""

    def __init__(self, matrix, *, scale=1.0, name="ngsolve-projected"):
        self._native_matrix = matrix
        self._scale = complex(scale)
        self._name = str(name)
        self.shape = (int(matrix.height), int(matrix.width))
        if self.shape[0] != self.shape[1]:
            raise ValueError("projected NGSolve operator must be square")
        self.dtype = np.dtype(complex)

    @property
    def T(self):
        return self

    @property
    def H(self):
        return self

    @property
    def matrix_free(self) -> bool:
        return True

    @property
    def is_hermitian(self) -> bool:
        return np.isreal(self._scale)

    def matvec(self, vector):
        values = np.asarray(vector, dtype=np.complex128)
        if values.shape != (self.shape[1],):
            raise ValueError(f"vector must have shape ({self.shape[1]},)")
        source = self._native_matrix.CreateColVector()
        target = self._native_matrix.CreateRowVector()
        source.FV().NumPy()[:] = values
        target.data = self._native_matrix * source
        return self._scale * target.FV().NumPy().copy()

    def matmat(self, matrix):
        values = np.asarray(matrix)
        if values.ndim == 1:
            return self.matvec(values)
        if values.ndim != 2 or values.shape[0] != self.shape[1]:
            raise ValueError(f"matrix must have shape ({self.shape[1]}, n_rhs)")
        return np.column_stack(
            [self.matvec(values[:, column]) for column in range(values.shape[1])]
        )

    def __matmul__(self, other):
        values = np.asarray(other)
        return self.matvec(values) if values.ndim == 1 else self.matmat(values)

    def to_dense(self):
        return self.matmat(np.eye(self.shape[1], dtype=complex))

    def stats(self):
        return {
            "operator": self._name,
            "matrix_free": True,
            "mode_count": int(self.shape[0]),
            "backend": "ngsolve-base-matrix",
        }


@dataclass(frozen=True)
class NGSolveProjectedInteraction:
    """Galerkin-project an actual NGSolve BEM/DtN ``BaseMatrix``.

    ``projections[i]`` maps reduced coefficients of ``bases[i]`` to the parent
    NGSolve space.  The C++ adapter applies ``P.H @ A @ P`` without exporting
    the parent BEM matrix or materializing the reduced interaction.
    """

    bases: tuple[SampledCurrentBasis, ...]
    parent_matrix: object
    projections: tuple[np.ndarray, ...]
    scale: complex = 1.0
    name: str = "ngsolve-bem-dtn-projected"

    def __post_init__(self) -> None:
        bases = tuple(self.bases)
        projections = tuple(np.asarray(item) for item in self.projections)
        if not bases or len(bases) != len(projections):
            raise ValueError("bases and projections must have equal non-zero length")
        parent_rows = int(self.parent_matrix.height)
        if parent_rows != int(self.parent_matrix.width):
            raise ValueError("parent_matrix must be square")
        for basis, projection in zip(bases, projections):
            if not isinstance(basis, SampledCurrentBasis):
                raise TypeError("bases must contain SampledCurrentBasis objects")
            if projection.shape != (parent_rows, basis.n_modes):
                raise ValueError(
                    "each projection must have shape "
                    f"({parent_rows}, basis.n_modes)"
                )
            if not np.all(np.isfinite(projection)):
                raise ValueError("projections must contain only finite values")
        object.__setattr__(self, "bases", bases)
        object.__setattr__(self, "projections", projections)
        if not self.name:
            raise ValueError("name must not be empty")

    @classmethod
    def from_laplace_sl(
        cls,
        fes,
        bases,
        projections,
        *,
        use_fmm: bool = False,
        scale=1.0,
    ):
        """Assemble ``ngsolve.bem.LaplaceSL`` and retain its BaseMatrix."""

        from ngsolve import ds
        from ngsolve.bem import LaplaceSL

        test, trial = fes.TnT()
        parent = (
            LaplaceSL(test.Trace() * ds, use_fmm=bool(use_fmm))
            * trial.Trace()
            * ds
        ).mat
        return cls(
            tuple(bases),
            parent,
            tuple(projections),
            scale=scale,
            name="ngsolve-bem-laplace-sl-projected",
        )

    @property
    def operator_scope(self) -> str:
        return "full"

    def build_operator(self, bases):
        bases = tuple(bases)
        if len(bases) != len(self.bases) or any(
            left is not right for left, right in zip(bases, self.bases)
        ):
            raise ValueError("assembled bases do not match projected BEM bases")
        projection = np.column_stack(self.projections)
        matrix_type = _radia_cpp_kernel("_ProjectedBaseMatrix")
        if matrix_type is None:
            raise RuntimeError(
                "Radia C++ extension lacks the NGSolve projected BaseMatrix adapter"
            )
        native = matrix_type(
            self.parent_matrix,
            np.ascontiguousarray(projection, dtype=np.complex128),
        )
        return NGSolveProjectedOperator(
            native,
            scale=self.scale,
            name=self.name,
        )

    def __call__(self, left, right):
        raise RuntimeError(
            "NGSolveProjectedInteraction is a full matrix-free operator; "
            "use build_operator(bases)"
        )


def _is_matrix_free_operator(value) -> bool:
    return (
        hasattr(value, "shape")
        and callable(getattr(value, "matvec", None))
        and not isinstance(value, np.ndarray)
    )


class SampledHACApKOperator:
    """Projected ``mu * sum_c B_c.T G B_c`` sampled-current operator."""

    def __init__(self, gram, charge_map, *, shape, mu: float, metadata=None):
        indptr, indices, data = charge_map
        rows, columns = (int(shape[0]), int(shape[1]))
        indptr = np.asarray(indptr, dtype=np.int32)
        indices = np.asarray(indices, dtype=np.int32)
        data = np.asarray(data, dtype=float)
        if rows != 3 * int(gram.ndof()):
            raise ValueError("charge-map rows must equal three times the HACApK point count")
        if columns < 1:
            raise ValueError("charge-map must contain at least one reduced mode")
        if indptr.shape != (rows + 1,) or indices.shape != data.shape:
            raise ValueError("invalid charge-map CSR arrays")
        gram.configure_vector_charge_map(
            indptr,
            indices,
            data,
            columns,
            3,
        )
        self._gram = gram
        self._mu = float(mu)
        self._metadata = {} if metadata is None else dict(metadata)
        self._point_count = int(gram.ndof())
        self._map_entries = int(data.size)
        self._map_bytes = int(indptr.nbytes + indices.nbytes + data.nbytes)
        self.shape = (columns, columns)
        self.dtype = np.dtype(float)

    @property
    def T(self):
        return self

    @property
    def H(self):
        return self

    @property
    def matrix_free(self) -> bool:
        return True

    @property
    def is_hermitian(self) -> bool:
        return True

    def _apply_real(self, vector) -> np.ndarray:
        return np.asarray(
            self._gram.apply_configured_demag(np.asarray(vector, dtype=float), True)
        )

    def matvec(self, vector):
        x = np.asarray(vector)
        if x.ndim != 1 or x.shape[0] != self.shape[1]:
            raise ValueError(f"vector must have shape ({self.shape[1]},)")
        if np.iscomplexobj(x):
            value = self._apply_real(x.real) + 1j * self._apply_real(x.imag)
        else:
            value = self._apply_real(x)
        return self._mu * value

    def matmat(self, matrix):
        values = np.asarray(matrix)
        if values.ndim == 1:
            return self.matvec(values)
        if values.ndim != 2 or values.shape[0] != self.shape[1]:
            raise ValueError(f"matrix must have shape ({self.shape[1]}, n_rhs)")
        return np.column_stack(
            [self.matvec(values[:, column]) for column in range(values.shape[1])]
        )

    def __matmul__(self, other):
        values = np.asarray(other)
        return self.matvec(values) if values.ndim == 1 else self.matmat(values)

    def to_dense(self):
        dense = self.matmat(np.eye(self.shape[1]))
        return 0.5 * (dense + dense.conj().T)

    def __array__(self, dtype=None, copy=None):
        dense = self.to_dense()
        if dtype is not None:
            dense = dense.astype(dtype, copy=False)
        if copy is False:
            return dense
        return np.array(dense, copy=True)

    def __getitem__(self, key):
        return self.to_dense()[key]

    def stats(self):
        values = dict(self._metadata)
        values.update(
            {
                "operator": "sampled-sum-component-BtGB",
                "matrix_free": True,
                "mode_count": int(self.shape[0]),
                "sample_point_count": self._point_count,
                "charge_map_entries": self._map_entries,
                "charge_map_bytes": self._map_bytes,
                "charge_map_backend": "native-cpp-component-csr",
                "hmatrix": dict(self._gram.stats()),
            }
        )
        return values


def _build_sampled_hacapk_operator(
    bases,
    *,
    mu: float,
    kernel_epsilon: float,
    aca_eps: float,
    leaf_size: int,
    eta: float,
    cross_only: bool,
    kernel: str = "laplace",
    reference_length: float = 1.0,
):
    bases = tuple(bases)
    active = [
        (index, basis)
        for index, basis in enumerate(bases)
        if basis.n_modes > 0 and basis.n_samples > 0
    ]
    total_modes = sum(basis.n_modes for basis in bases)
    if not active or total_modes < 1:
        return None
    for _, basis in active:
        if np.iscomplexobj(basis.modes) and np.any(np.imag(basis.modes) != 0.0):
            raise ValueError("sampled HACApK interaction currently requires real current modes")

    points = np.vstack([basis.points for _, basis in active])
    weights = np.concatenate([basis.weights for _, basis in active])
    point_offsets = {}
    cursor = 0
    for index, basis in active:
        point_offsets[index] = cursor
        cursor += basis.n_samples
    mode_offsets = _basis_offsets(bases)

    indptr = [0]
    indices = []
    data = []
    point_count = points.shape[0]
    for component in range(3):
        for index, basis in active:
            values = np.asarray(basis.modes[:, :, component].T.real, dtype=float)
            for row in values:
                nonzero = np.flatnonzero(row)
                indices.extend(mode_offsets[index][0] + nonzero)
                data.extend(row[nonzero])
                indptr.append(len(indices))

    gram_type = _radia_cpp_kernel("_ChargeGramHMatrix")
    factory_name = (
        "from_sampled_laplace"
        if kernel == "laplace"
        else "from_sampled_planar_log"
    )
    if kernel not in {"laplace", "planar-log"}:
        raise ValueError("kernel must be 'laplace' or 'planar-log'")
    if gram_type is None or not hasattr(gram_type, factory_name):
        raise RuntimeError(
            f"Radia C++ extension lacks sampled {kernel} HACApK support; rebuild _radia_pybind"
        )
    factory = getattr(gram_type, factory_name)
    arguments = [
        np.ascontiguousarray(points.reshape(-1), dtype=float),
        np.ascontiguousarray(weights, dtype=float),
        float(kernel_epsilon),
    ]
    if kernel == "planar-log":
        arguments.append(float(reference_length))
    gram = factory(
        *arguments,
        float(aca_eps),
        int(leaf_size),
        float(eta),
        True,
    )
    return SampledHACApKOperator(
        gram,
        (
            np.asarray(indptr, dtype=np.int32),
            np.asarray(indices, dtype=np.int32),
            np.asarray(data, dtype=float),
        ),
        shape=(3 * point_count, total_modes),
        mu=mu,
        metadata={
            "backend": f"hacapk-sampled-{kernel}",
            "cross_only": bool(cross_only),
            "kernel_scope": "full-stable-gram",
            "diagonal_policy": (
                "reduced-correction" if cross_only else "sampled-hacapk"
            ),
            "basis_count": len(bases),
            "active_basis_count": len(active),
            "kernel_epsilon": float(kernel_epsilon),
            "reference_length": (
                float(reference_length) if kernel == "planar-log" else None
            ),
        },
    )


class ReducedBlockHMatrixOperator:
    """Dense remainder plus additive symmetric matrix-free block operators."""

    def __init__(self, dense, operator_blocks):
        base = np.asarray(dense)
        if base.ndim != 2 or base.shape[0] != base.shape[1]:
            raise ValueError("dense reduced block matrix must be square")
        self._dense = np.array(base, copy=True)
        self._operator_blocks = tuple(operator_blocks)
        self.shape = self._dense.shape
        self.dtype = np.result_type(self._dense.dtype, float)

    @property
    def T(self):
        return self

    @property
    def H(self):
        return self

    @property
    def matrix_free(self) -> bool:
        return True

    @property
    def is_hermitian(self) -> bool:
        return True

    def matvec(self, vector):
        x = np.asarray(vector)
        if x.ndim != 1 or x.shape[0] != self.shape[1]:
            raise ValueError(f"vector must have shape ({self.shape[1]},)")
        result = self._dense @ x
        for start, stop, operator in self._operator_blocks:
            result[start:stop] += operator.matvec(x[start:stop])
        return result

    def matmat(self, matrix):
        values = np.asarray(matrix)
        if values.ndim == 1:
            return self.matvec(values)
        if values.ndim != 2 or values.shape[0] != self.shape[1]:
            raise ValueError(f"matrix must have shape ({self.shape[1]}, n_rhs)")
        result = self._dense @ values
        for start, stop, operator in self._operator_blocks:
            result[start:stop, :] += operator.matmat(values[start:stop, :])
        return result

    def __matmul__(self, other):
        return self.matmat(other)

    def to_dense(self):
        dense = np.array(self._dense, copy=True)
        for start, stop, operator in self._operator_blocks:
            dense[start:stop, start:stop] += operator.to_dense()
        return 0.5 * (dense + dense.conj().T)

    def __array__(self, dtype=None, copy=None):
        dense = self.to_dense()
        if dtype is not None:
            dense = dense.astype(dtype, copy=False)
        if copy is False:
            return dense
        return np.array(dense, copy=True)

    def __getitem__(self, key):
        return self.to_dense()[key]

    def stats(self):
        return {
            "matrix_free": True,
            "mode_count": int(self.shape[0]),
            "operator_block_count": len(self._operator_blocks),
            "dense_remainder_bytes": int(self._dense.nbytes),
            "dense_cross_bytes": int(self._dense.nbytes),
            "operator_blocks": [
                {
                    "start": int(start),
                    "stop": int(stop),
                    "stats": operator.stats(),
                }
                for start, stop, operator in self._operator_blocks
            ],
        }


class ShiftedReducedOperator:
    """Matrix-free ``dense + scale * inductance`` reduced impedance."""

    def __init__(self, dense, inductance, scale):
        self._dense = np.asarray(dense)
        self._inductance = inductance
        self._scale = scale
        self.shape = self._dense.shape
        self.dtype = np.result_type(self._dense.dtype, scale, inductance.dtype)

    def matvec(self, vector):
        x = np.asarray(vector)
        return self._dense @ x + self._scale * self._inductance.matvec(x)

    def matmat(self, matrix):
        values = np.asarray(matrix)
        return self._dense @ values + self._scale * self._inductance.matmat(values)

    def __matmul__(self, other):
        values = np.asarray(other)
        return self.matvec(values) if values.ndim == 1 else self.matmat(values)

    def to_dense(self):
        return self._dense + self._scale * self._inductance.to_dense()

    def __array__(self, dtype=None, copy=None):
        dense = self.to_dense()
        if dtype is not None:
            dense = dense.astype(dtype, copy=False)
        if copy is False:
            return dense
        return np.array(dense, copy=True)

    def __getitem__(self, key):
        return self.to_dense()[key]


class CoupledReducedOperator:
    """Native-ready HDiv/HCurl block operator with an HCurl sub-operator."""

    def __init__(self, magnetic, upper, lower, eddy):
        upper = np.asarray(upper)
        lower = np.asarray(lower)
        magnetic_shape = tuple(magnetic.shape)
        if len(magnetic_shape) != 2 or magnetic_shape[0] != magnetic_shape[1]:
            raise ValueError("magnetic block must be square")
        if upper.shape != (magnetic_shape[0], eddy.shape[1]):
            raise ValueError("upper coupling block shape mismatch")
        if lower.shape != (eddy.shape[0], magnetic_shape[1]):
            raise ValueError("lower coupling block shape mismatch")
        self._magnetic = (
            magnetic if _is_matrix_free_operator(magnetic)
            else np.array(magnetic, copy=True)
        )
        self._upper = np.array(upper, copy=True)
        self._lower = np.array(lower, copy=True)
        self._eddy = eddy
        self._magnetic_size = int(magnetic_shape[0])
        total = self._magnetic_size + int(eddy.shape[0])
        self.shape = (total, total)
        magnetic_dtype = (
            magnetic.dtype
            if hasattr(magnetic, "dtype")
            else np.asarray(magnetic).dtype
        )
        self.dtype = np.result_type(magnetic_dtype, upper, lower, eddy.dtype)

    @property
    def matrix_free(self) -> bool:
        return True

    def matvec(self, vector):
        values = np.asarray(vector)
        if values.shape != (self.shape[1],):
            raise ValueError(f"vector must have shape ({self.shape[1]},)")
        mh = self._magnetic_size
        magnetic_action = (
            self._magnetic.matvec(values[:mh])
            if _is_matrix_free_operator(self._magnetic)
            else self._magnetic @ values[:mh]
        )
        magnetic = magnetic_action + self._upper @ values[mh:]
        eddy = self._lower @ values[:mh] + self._eddy.matvec(values[mh:])
        return np.concatenate((magnetic, eddy))

    def matmat(self, matrix):
        values = np.asarray(matrix)
        if values.ndim == 1:
            return self.matvec(values)
        if values.ndim != 2 or values.shape[0] != self.shape[1]:
            raise ValueError(f"matrix must have shape ({self.shape[1]}, n_rhs)")
        return np.column_stack(
            [self.matvec(values[:, column]) for column in range(values.shape[1])]
        )

    def __matmul__(self, other):
        values = np.asarray(other)
        return self.matvec(values) if values.ndim == 1 else self.matmat(values)

    def to_dense(self):
        mh = self._magnetic_size
        out = np.empty(self.shape, dtype=self.dtype)
        out[:mh, :mh] = (
            self._magnetic.to_dense()
            if _is_matrix_free_operator(self._magnetic)
            else self._magnetic
        )
        out[:mh, mh:] = self._upper
        out[mh:, :mh] = self._lower
        out[mh:, mh:] = self._eddy.to_dense()
        return out


def _flatten_native_operator(operator, *, offset=0, scale=1.0):
    """Return ``(dense, BaseMatrix terms)`` without materializing H-matrices."""

    n = int(operator.shape[0])
    if tuple(operator.shape) != (n, n):
        raise ValueError("native reduced operator must be square")
    if isinstance(operator, NGSolveProjectedOperator):
        return (
            np.zeros((n, n), dtype=complex),
            [(operator._native_matrix, offset, offset + n, scale * operator._scale)],
        )
    if isinstance(operator, ShiftedReducedOperator):
        dense, terms = _flatten_native_operator(
            operator._inductance,
            offset=offset,
            scale=scale * operator._scale,
        )
        dense += scale * np.asarray(operator._dense)
        return dense, terms
    if isinstance(operator, ReducedBlockHMatrixOperator):
        dense = scale * np.asarray(operator._dense, dtype=complex)
        terms = []
        for start, stop, block in operator._operator_blocks:
            block_dense, block_terms = _flatten_native_operator(
                block,
                offset=offset + start,
                scale=scale,
            )
            dense[start:stop, start:stop] += block_dense
            terms.extend(block_terms)
        return dense, terms
    if isinstance(operator, CoupledReducedOperator):
        mh = operator._magnetic_size
        dense = np.zeros(operator.shape, dtype=complex)
        dense[:mh, mh:] = scale * operator._upper
        dense[mh:, :mh] = scale * operator._lower
        if _is_matrix_free_operator(operator._magnetic):
            magnetic_dense, magnetic_terms = _flatten_native_operator(
                operator._magnetic,
                offset=offset,
                scale=scale,
            )
            dense[:mh, :mh] = magnetic_dense
        else:
            dense[:mh, :mh] = scale * operator._magnetic
            magnetic_terms = []
        eddy_dense, eddy_terms = _flatten_native_operator(
            operator._eddy,
            offset=offset + mh,
            scale=scale,
        )
        dense[mh:, mh:] = eddy_dense
        return dense, magnetic_terms + eddy_terms
    gram = getattr(operator, "_gram", None)
    local_scale = getattr(operator, "_mu", None)
    if gram is not None and local_scale is not None:
        return (
            np.zeros((n, n), dtype=complex),
            [(gram.demag_matrix(), offset, offset + n, scale * local_scale)],
        )
    raise TypeError(
        f"matrix-free operator {type(operator).__name__} has no native BaseMatrix path"
    )


def _native_reduced_base_matrix(operator):
    """Build one C++/NGSolve BaseMatrix from independent reduced terms."""

    dense, terms = _flatten_native_operator(operator)
    matrix_type = _radia_cpp_kernel("_ReducedBlockMatrix")
    if matrix_type is None:
        raise RuntimeError(
            "Radia C++ extension lacks the reduced NGSolve BaseMatrix adapter"
        )
    return matrix_type(
        np.ascontiguousarray(dense, dtype=np.complex128),
        [term[0] for term in terms],
        np.asarray([term[1] for term in terms], dtype=np.int32),
        np.asarray([term[2] for term in terms], dtype=np.int32),
        np.asarray([term[3] for term in terms], dtype=np.complex128),
    )


def _solve_native_reduced_operator(
    operator,
    rhs,
    *,
    solver: str = "gmres",
    tolerance: float = 1.0e-10,
    max_iterations: int | None = None,
    restart: int | None = None,
):
    """Solve a matrix-free reduced system in NGSolve's BaseMatrix style."""

    import ngsolve as ng
    from radia import sparsesolv_ngsolve

    solver = str(solver).lower()
    if solver not in {"gmres", "cocr"}:
        raise ValueError("solver must be 'gmres' or 'cocr'")
    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance must be positive")
    native = _native_reduced_base_matrix(operator)
    preconditioner = native.diagonal_preconditioner()
    n = int(native.height)
    values = np.asarray(rhs)
    vector_result = values.ndim == 1
    if vector_result:
        values = values[:, np.newaxis]
    if values.ndim != 2 or values.shape[0] != n:
        raise ValueError(f"rhs must have shape ({n}, n_rhs)")
    max_iterations = max(100, 10 * n) if max_iterations is None else int(max_iterations)
    restart = min(max_iterations, 50) if restart is None else int(restart)
    if max_iterations <= 0 or restart <= 0:
        raise ValueError("max_iterations and restart must be positive")
    if solver == "gmres":
        inverse = sparsesolv_ngsolve.GMRESSolver(
            native,
            preconditioner,
            maxiter=max_iterations,
            tol=float(tolerance),
            restart=restart,
        )
    else:
        inverse = sparsesolv_ngsolve.COCRSolver(
            native,
            preconditioner,
            maxiter=max_iterations,
            tol=float(tolerance),
        )
    solutions = []
    residuals = []
    iteration_counts = []
    refinement_counts = []
    with ng.TaskManager():
        for column in range(values.shape[1]):
            source = native.CreateColVector()
            source.FV().NumPy()[:] = np.asarray(
                values[:, column], dtype=np.complex128
            )
            solution = native.CreateColVector()
            solution.FV().NumPy()[:] = 0.0
            residual = native.CreateRowVector()
            denominator = max(float(np.linalg.norm(source.FV().NumPy())), 1.0)
            total_iterations = 0
            relative_residual = np.inf
            refinement_count = 0
            correction_rhs = source
            for refinement_count in range(8):
                correction = native.CreateColVector()
                correction.data = inverse * correction_rhs
                solution.data += correction
                solver_iterations = int(getattr(inverse, "iterations", -1))
                if solver_iterations >= 0:
                    total_iterations += solver_iterations
                residual.data = source - native * solution
                relative_residual = float(
                    np.linalg.norm(residual.FV().NumPy()) / denominator
                )
                if relative_residual <= max(10.0 * tolerance, 1.0e-12):
                    break
                correction_rhs = residual
            if relative_residual > max(10.0 * tolerance, 1.0e-12):
                raise RuntimeError(
                    f"native reduced {solver.upper()} residual {relative_residual:.3e} "
                    f"exceeds tolerance {tolerance:.3e}"
                )
            solutions.append(solution.FV().NumPy().copy())
            residuals.append(relative_residual)
            iteration_counts.append(total_iterations)
            refinement_counts.append(refinement_count)
    result = np.column_stack(solutions)
    diagnostics = {
        "backend": f"ngsolve-base-matrix-{solver}",
        "iterations": max(iteration_counts, default=0),
        "defect_corrections": max(refinement_counts, default=0),
        "relative_residual_max": max(residuals, default=0.0),
        "native_term_count": int(native.term_count),
        "preconditioner": "native-reduced-block-jacobi",
    }
    return (result[:, 0] if vector_result else result), diagnostics


def _validate_interaction_block(block, left: SampledCurrentBasis, right: SampledCurrentBasis):
    if _is_matrix_free_operator(block):
        if tuple(block.shape) != (left.n_modes, right.n_modes):
            raise ValueError(
                "interaction backend returned an operator with shape "
                f"{block.shape}, expected {(left.n_modes, right.n_modes)}"
            )
        return block
    arr = np.asarray(block)
    if arr.shape != (left.n_modes, right.n_modes):
        raise ValueError(
            "interaction backend returned a block with shape "
            f"{arr.shape}, expected {(left.n_modes, right.n_modes)}"
        )
    if not np.all(np.isfinite(arr)):
        raise ValueError("interaction backend returned non-finite values")
    return arr


def CurrentMagneticFluxDensitySamples(
    current_basis: SampledCurrentBasis,
    target_points,
    *,
    mu: float = MU0,
    kernel_epsilon: float = 0.0,
) -> np.ndarray:
    """Evaluate ``B`` from sampled volume/surface currents at target points.

    Returns an array with shape ``(n_current_modes, n_targets, 3)`` using the
    quasi-static Biot-Savart kernel

    ``B(x) = mu/(4*pi) int J(y) x (x-y) / |x-y|^3 dy``.

    The same routine accepts surface-Omega currents because their quadrature
    weights are surface weights and their modes are surface current densities.
    """

    if not isinstance(current_basis, SampledCurrentBasis):
        raise TypeError("current_basis must be a SampledCurrentBasis")
    if mu <= 0.0:
        raise ValueError("mu must be positive")
    if kernel_epsilon < 0.0:
        raise ValueError("kernel_epsilon must be non-negative")
    targets = _as_points(target_points, "target_points")
    diff = targets[:, np.newaxis, :] - current_basis.points[np.newaxis, :, :]
    radius2 = np.einsum("ijk,ijk->ij", diff, diff) + kernel_epsilon**2
    if np.any(radius2 <= 0.0):
        raise ValueError(
            "target_points contain a source point; use a positive kernel_epsilon "
            "or a proper singular quadrature backend"
        )
    kernel = current_basis.weights[np.newaxis, :] / (radius2 ** 1.5)
    cross = np.cross(
        current_basis.modes[:, np.newaxis, :, :],
        diff[np.newaxis, :, :, :],
    )
    return (mu / (4.0 * np.pi)) * np.einsum("misk,is->mik", cross, kernel)


def MagnetizationCurrentCoupling(
    magnetization_basis: SampledMagnetizationBasis,
    current_basis: SampledCurrentBasis,
    *,
    mu: float = MU0,
    kernel_epsilon: float = 0.0,
) -> np.ndarray:
    """Return the rectangular coupling ``C_ij = int M_i dot B[J_j] dV``.

    Rows correspond to magnetization/HDiv modes and columns to eddy-current
    modes.  This is the first thin bridge between Radia's HDiv-VIM magnetic
    branch and the high-order HCurl response-compressed eddy branch.
    """

    if not isinstance(magnetization_basis, SampledMagnetizationBasis):
        raise TypeError("magnetization_basis must be a SampledMagnetizationBasis")
    fields = CurrentMagneticFluxDensitySamples(
        current_basis,
        magnetization_basis.points,
        mu=mu,
        kernel_epsilon=kernel_epsilon,
    )
    return np.einsum(
        "aik,bik,i->ab",
        magnetization_basis.modes.conj(),
        fields,
        magnetization_basis.weights,
    )


def _port_rhs_matrix(rhs, n: int) -> np.ndarray:
    matrix = np.asarray(rhs)
    if matrix.ndim == 1:
        matrix = matrix[:, np.newaxis]
    if matrix.ndim != 2 or matrix.shape[0] != n:
        raise ValueError(f"rhs must have shape ({n}, n_ports)")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("rhs contains non-finite values")
    return matrix


def _as_finite_matrix(matrix, name: str) -> np.ndarray:
    arr = np.asarray(matrix)
    if arr.ndim != 2:
        raise ValueError(f"{name} must be a 2D matrix")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def _square_matrix(matrix, n: int, name: str) -> np.ndarray:
    arr = _as_finite_matrix(matrix, name)
    if arr.shape != (n, n):
        raise ValueError(f"{name} must have shape ({n}, {n})")
    return arr


def _solve_reduced_linear(operator, rhs) -> np.ndarray:
    """Solve a small reduced system through the native kernel when available."""

    matrix = np.asarray(operator)
    values = np.asarray(rhs)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("operator must be a square 2-D matrix")
    vector_rhs = values.ndim == 1
    if vector_rhs:
        values = values[:, np.newaxis]
    if values.ndim != 2 or values.shape[0] != matrix.shape[0]:
        raise ValueError("rhs rows must match operator rows")
    if not np.all(np.isfinite(matrix)) or not np.all(np.isfinite(values)):
        raise ValueError("operator and rhs must contain only finite values")

    func = _radia_cpp_kernel("_HybridVIMSolve")
    if func is not None and (np.iscomplexobj(matrix) or np.iscomplexobj(values)):
        solution = func(
            np.ascontiguousarray(matrix, dtype=np.complex128),
            np.ascontiguousarray(values, dtype=np.complex128),
        )
    else:
        solution = np.linalg.solve(matrix, values)
    solution = np.asarray(solution)
    return solution[:, 0] if vector_rhs else solution


def _label_tuple(labels, name: str) -> tuple[str, ...]:
    if labels is None:
        return ()
    if isinstance(labels, str):
        return (labels,)
    try:
        return tuple(str(label) for label in labels)
    except TypeError as exc:
        raise ValueError(f"{name} must be a string or iterable of strings") from exc


def _check_positive_material_coeff(value, name: str) -> None:
    if value is None:
        return
    if isinstance(value, dict):
        for key, item in value.items():
            _check_positive_material_coeff(item, f"{name}[{key!r}]")
        return
    if np.isscalar(value):
        if not np.isfinite(value) or float(value) <= 0.0:
            raise ValueError(f"{name} must be positive")


def _as_real_matrix(matrix, name: str) -> np.ndarray:
    arr = np.asarray(matrix, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"{name} must be a 2D matrix")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return np.ascontiguousarray(arr)


def _fro_norm(matrix: np.ndarray) -> float:
    return float(np.linalg.norm(matrix))


def _evrs_tmethod_algebra_numpy(
    curl_map,
    div_map,
    grad_map,
    evrs_map,
    resistance_current,
    inductance_current,
    port_current,
) -> dict[str, object]:
    c = _as_real_matrix(curl_map, "curl_map")
    d = _as_real_matrix(div_map, "div_map")
    g = _as_real_matrix(grad_map, "grad_map")
    q = _as_real_matrix(evrs_map, "evrs_map")
    mr = _as_real_matrix(resistance_current, "resistance_current")
    ml = _as_real_matrix(inductance_current, "inductance_current")
    p = _as_real_matrix(port_current, "port_current")
    n_current, n_t = c.shape
    if d.shape[1] != n_current:
        raise ValueError("div_map columns must match curl_map rows")
    if g.shape[0] != n_t:
        raise ValueError("grad_map rows must match curl_map columns")
    if q.shape[0] != n_t:
        raise ValueError("evrs_map rows must match curl_map columns")
    if mr.shape != (n_current, n_current):
        raise ValueError("resistance_current must be square in current space")
    if ml.shape != (n_current, n_current):
        raise ValueError("inductance_current must be square in current space")
    if p.shape[0] != n_current:
        raise ValueError("port_current rows must match current-space dimension")

    current_evrs = c @ q
    resistance_t = c.T @ mr @ c
    inductance_t = c.T @ ml @ c
    resistance_evrs = q.T @ resistance_t @ q
    inductance_evrs = q.T @ inductance_t @ q
    port_t = c.T @ p
    port_evrs = q.T @ port_t
    resistance_current_evrs = current_evrs.T @ mr @ current_evrs
    inductance_current_evrs = current_evrs.T @ ml @ current_evrs
    diagnostics = {
        "n_current": n_current,
        "n_t": n_t,
        "n_phi": int(g.shape[1]),
        "n_evrs": int(q.shape[1]),
        "n_ports": int(p.shape[1]),
        "n_rho": int(d.shape[0]),
        "div_curl_norm": _fro_norm(d @ c),
        "div_evrs_norm": _fro_norm(d @ current_evrs),
        "resistance_gauge_norm": _fro_norm(resistance_t @ g),
        "inductance_gauge_norm": _fro_norm(inductance_t @ g),
        "port_gauge_norm": _fro_norm(g.T @ port_t),
        "resistance_symmetry_norm": _fro_norm(resistance_t - resistance_t.T),
        "inductance_symmetry_norm": _fro_norm(inductance_t - inductance_t.T),
        "evrs_resistance_symmetry_norm": _fro_norm(resistance_evrs - resistance_evrs.T),
        "evrs_inductance_symmetry_norm": _fro_norm(inductance_evrs - inductance_evrs.T),
        "evrs_resistance_galerkin_residual": _fro_norm(
            resistance_evrs - resistance_current_evrs
        ),
        "evrs_inductance_galerkin_residual": _fro_norm(
            inductance_evrs - inductance_current_evrs
        ),
    }
    return {
        "current_evrs": current_evrs,
        "resistance_t": resistance_t,
        "inductance_t": inductance_t,
        "resistance_evrs": resistance_evrs,
        "inductance_evrs": inductance_evrs,
        "port_t": port_t,
        "port_evrs": port_evrs,
        "diagnostics": diagnostics,
    }


def EVRSTMethodAlgebra(
    curl_map,
    div_map,
    grad_map,
    evrs_map,
    resistance_current,
    inductance_current,
    port_current,
    *,
    backend: str = "auto",
) -> dict[str, object]:
    """Evaluate the EVRS/T-method projection algebra.

    ``T`` is the HCurl variable, ``J = C T`` is the current-space variable, and
    ``Q`` maps retained EVRS coordinates into the parent T space.  The returned
    matrices include ``C Q``, ``C^T M C``, ``Q^T C^T M C Q`` and the reduced
    port RHS ``Q^T C^T P``.
    """

    backend = backend.lower()
    if backend not in {"auto", "cpp", "python"}:
        raise ValueError("backend must be 'auto', 'cpp', or 'python'")
    args = (
        _as_real_matrix(curl_map, "curl_map"),
        _as_real_matrix(div_map, "div_map"),
        _as_real_matrix(grad_map, "grad_map"),
        _as_real_matrix(evrs_map, "evrs_map"),
        _as_real_matrix(resistance_current, "resistance_current"),
        _as_real_matrix(inductance_current, "inductance_current"),
        _as_real_matrix(port_current, "port_current"),
    )
    if backend in {"auto", "cpp"}:
        try:
            from radia import _radia_pybind as _radia_cpp
        except Exception as exc:
            if backend == "cpp":
                raise RuntimeError("Radia C++ extension is not importable") from exc
        else:
            func = getattr(_radia_cpp, "_EVRSTMethodAlgebra", None)
            if func is not None:
                return func(*args)
            if backend == "cpp":
                raise RuntimeError("Radia C++ extension lacks _EVRSTMethodAlgebra; rebuild _radia_pybind")
    return _evrs_tmethod_algebra_numpy(*args)


def _relative_hermitian_error(matrix) -> float:
    if _is_matrix_free_operator(matrix):
        if getattr(matrix, "is_hermitian", False):
            return 0.0
        if callable(getattr(matrix, "to_dense", None)):
            return _relative_hermitian_error(matrix.to_dense())
        raise TypeError("matrix-free operator must declare is_hermitian or provide to_dense")
    norm = float(np.linalg.norm(matrix))
    if norm == 0.0:
        return 0.0
    return float(np.linalg.norm(matrix - matrix.conj().T) / norm)


def _min_hermitian_eigenvalue(matrix) -> float:
    if _is_matrix_free_operator(matrix):
        if not getattr(matrix, "is_hermitian", False):
            raise TypeError("matrix-free inductance operator must declare is_hermitian")
        n = int(matrix.shape[0])
        if n == 0:
            return 0.0
        # This is an explicit reduced-order diagnostic, not the production
        # solve path.  Keep SciPy out of the HCurl solver and inspect the small
        # retained operator with NumPy only when diagnostics are requested.
        dense = matrix.matmat(np.eye(n))
        hermitian = 0.5 * (dense + dense.conj().T)
        return float(np.min(np.linalg.eigvalsh(hermitian).real))
    if matrix.size == 0:
        return 0.0
    hermitian = 0.5 * (matrix + matrix.conj().T)
    return float(np.min(np.linalg.eigvalsh(hermitian).real))


def _interaction_backend_name(interaction) -> str:
    if interaction is None:
        return "sampled-laplace"
    name = getattr(interaction, "name", None)
    if name:
        return str(name)
    name = getattr(interaction, "__name__", None)
    if name:
        return str(name)
    return interaction.__class__.__name__


@dataclass(frozen=True)
class MixedGalerkinOrthogonalization:
    """Block-biorthogonal trial/test bases for an IGTE mixed reduction."""

    trial_transform: np.ndarray
    test_transform: np.ndarray
    reduced_operator: np.ndarray
    keep_indices: np.ndarray
    eliminate_indices: np.ndarray
    diagnostics_data: dict[str, float | int]

    @property
    def rank(self) -> int:
        return int(self.trial_transform.shape[1])

    def diagnostics(self) -> dict[str, float | int]:
        return dict(self.diagnostics_data)

    def solve(self, operator, rhs) -> np.ndarray:
        """Solve with the Schur basis and restore the eliminated-block lift.

        ``T c`` is the homogeneous harmonic extension of the retained
        coefficients.  When the eliminated block is directly excited, the
        complete solution also contains the particular lift
        ``Z_ee^-1 b_e``.  Keeping this term explicit prevents a correct Schur
        operator from being paired with an incomplete reconstructed field.
        """

        matrix = np.asarray(operator)
        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
            raise ValueError("operator must be square")
        if matrix.shape[0] != self.trial_transform.shape[0]:
            raise ValueError("operator size does not match the mixed basis")
        values = np.asarray(rhs)
        vector_rhs = values.ndim == 1
        if vector_rhs:
            values = values[:, np.newaxis]
        if values.ndim != 2 or values.shape[0] != matrix.shape[0]:
            raise ValueError("rhs must have shape (n_modes,) or (n_modes, n_rhs)")

        scale = max(float(np.linalg.norm(matrix)), np.finfo(float).tiny)
        trial_defect = float(
            np.linalg.norm(matrix[self.eliminate_indices, :] @ self.trial_transform)
            / scale
        )
        test_defect = float(
            np.linalg.norm(
                self.test_transform.T @ matrix[:, self.eliminate_indices]
            )
            / scale
        )
        projected = self.test_transform.T @ matrix @ self.trial_transform
        projected_scale = max(
            float(np.linalg.norm(self.reduced_operator)),
            np.finfo(float).tiny,
        )
        projected_error = float(
            np.linalg.norm(projected - self.reduced_operator) / projected_scale
        )
        if max(trial_defect, test_defect, projected_error) > 1.0e-10:
            raise ValueError(
                "operator does not match the matrix used to construct this "
                "mixed-Galerkin orthogonalization"
            )

        reduced_rhs = self.test_transform.T @ values
        reduced_coefficients = np.linalg.solve(self.reduced_operator, reduced_rhs)
        solution = self.trial_transform @ reduced_coefficients
        eliminated = self.eliminate_indices
        z_ee = matrix[np.ix_(eliminated, eliminated)]
        solution[eliminated, :] += np.linalg.solve(z_ee, values[eliminated, :])
        return solution[:, 0] if vector_rhs else solution


@dataclass(frozen=True)
class SurfaceImpedanceGram:
    """Preassembled local SIBC term ``int_Gamma Zs(x) jt.jv dS``.

    A scalar surface impedance can be multiplied by ``surface_mass``.  ESIM,
    curvature corrections, and material patches instead produce a value per
    surface quadrature sample.  This object keeps that weighted Gram explicit
    so local impedances are not silently replaced by an arithmetic mean.
    """

    matrix: np.ndarray
    sample_values: np.ndarray
    surface_mass_reference: np.ndarray
    label: str = "local-sibc"

    def __post_init__(self) -> None:
        matrix = np.asarray(self.matrix, dtype=complex)
        values = np.asarray(self.sample_values, dtype=complex).reshape(-1)
        surface_mass = np.asarray(self.surface_mass_reference, dtype=complex)
        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
            raise ValueError("surface impedance Gram matrix must be square")
        if surface_mass.shape != matrix.shape:
            raise ValueError("surface mass reference must match the Gram matrix")
        if values.size == 0:
            raise ValueError("surface impedance samples must not be empty")
        if not (
            np.all(np.isfinite(matrix))
            and np.all(np.isfinite(values))
            and np.all(np.isfinite(surface_mass))
        ):
            raise ValueError("surface impedance Gram contains non-finite values")
        object.__setattr__(self, "matrix", np.array(matrix, copy=True))
        object.__setattr__(self, "sample_values", np.array(values, copy=True))
        object.__setattr__(
            self,
            "surface_mass_reference",
            np.array(surface_mass, copy=True),
        )

    @property
    def dissipative_matrix(self) -> np.ndarray:
        """Return the Hermitian dissipative part used for average loss."""

        return 0.5 * (self.matrix + self.matrix.conj().T)

    def diagnostics(self, *, passive_tol: float = 1.0e-10) -> dict[str, object]:
        """Return JSON-ready local-impedance and passivity diagnostics."""

        min_dissipative = _min_hermitian_eigenvalue(self.dissipative_matrix)
        return {
            "label": self.label,
            "modes": int(self.matrix.shape[0]),
            "samples": int(self.sample_values.size),
            "real_min": float(np.min(self.sample_values.real)),
            "real_max": float(np.max(self.sample_values.real)),
            "imag_min": float(np.min(self.sample_values.imag)),
            "imag_max": float(np.max(self.sample_values.imag)),
            "minimum_dissipative_eigenvalue": min_dissipative,
            "passive": bool(min_dissipative >= -passive_tol),
        }


def _surface_impedance_term(surface_mass, surface_impedance) -> np.ndarray:
    """Resolve a scalar or preassembled local SIBC contribution."""

    if isinstance(surface_impedance, SurfaceImpedanceGram):
        if surface_impedance.matrix.shape != surface_mass.shape:
            raise ValueError(
                "surface impedance Gram shape does not match the hybrid system"
            )
        scale = max(float(np.linalg.norm(surface_mass)), np.finfo(float).tiny)
        if (
            np.linalg.norm(surface_impedance.surface_mass_reference - surface_mass)
            / scale
            > 1.0e-10
        ):
            raise ValueError(
                "surface impedance Gram was assembled for a different hybrid system"
            )
        return surface_impedance.matrix
    if np.isscalar(surface_impedance):
        return complex(surface_impedance) * surface_mass
    raise TypeError(
        "surface_impedance must be a scalar or SurfaceImpedanceGram; "
        "assemble spatial values with AssembleSurfaceImpedanceGram"
    )


def _surface_impedance_dissipative_term(surface_mass, surface_impedance) -> np.ndarray:
    """Return the Hermitian loss contribution of a SIBC term."""

    term = _surface_impedance_term(surface_mass, surface_impedance)
    return 0.5 * (term + term.conj().T)


def _mixed_galerkin_orthogonalize_operator(
    operator,
    keep_indices,
    eliminate_indices,
) -> MixedGalerkinOrthogonalization:
    """Return exact block-biorthogonal bases for a square operator."""

    matrix = np.asarray(operator)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("operator must be square")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("operator contains non-finite values")
    keep = np.asarray(keep_indices, dtype=int).reshape(-1)
    eliminate = np.asarray(eliminate_indices, dtype=int).reshape(-1)
    if keep.size == 0:
        raise ValueError("keep_indices must not be empty")
    if eliminate.size == 0:
        raise ValueError("eliminate_indices must not be empty")
    combined = np.concatenate((keep, eliminate))
    if np.any(combined < 0) or np.any(combined >= matrix.shape[0]):
        raise ValueError("mixed Galerkin indices are out of range")
    if np.unique(combined).size != combined.size:
        raise ValueError("keep_indices and eliminate_indices must not overlap or repeat")
    if combined.size != matrix.shape[0]:
        raise ValueError("keep_indices and eliminate_indices must cover all modes")

    z_kk = matrix[np.ix_(keep, keep)]
    z_ke = matrix[np.ix_(keep, eliminate)]
    z_ek = matrix[np.ix_(eliminate, keep)]
    z_ee = matrix[np.ix_(eliminate, eliminate)]
    trial = np.zeros((matrix.shape[0], keep.size), dtype=matrix.dtype)
    test = np.zeros_like(trial)
    trial[keep, :] = np.eye(keep.size, dtype=matrix.dtype)
    test[keep, :] = np.eye(keep.size, dtype=matrix.dtype)
    trial[eliminate, :] = -np.linalg.solve(z_ee, z_ek)
    test[eliminate, :] = -np.linalg.solve(z_ee.T, z_ke.T)

    reduced = test.T @ matrix @ trial
    schur = z_kk - z_ke @ np.linalg.solve(z_ee, z_ek)
    scale = max(float(np.linalg.norm(matrix)), np.finfo(float).tiny)
    reduced_scale = max(float(np.linalg.norm(schur)), np.finfo(float).tiny)
    diagnostics = {
        "full_modes": int(matrix.shape[0]),
        "retained_modes": int(keep.size),
        "eliminated_modes": int(eliminate.size),
        "trial_orthogonality_relative_defect": float(
            np.linalg.norm(matrix[np.ix_(eliminate, np.arange(matrix.shape[0]))] @ trial)
            / scale
        ),
        "test_orthogonality_relative_defect": float(
            np.linalg.norm(test.T @ matrix[np.ix_(np.arange(matrix.shape[0]), eliminate)])
            / scale
        ),
        "schur_relative_error": float(
            np.linalg.norm(reduced - schur) / reduced_scale
        ),
        "trial_test_relative_difference": float(
            np.linalg.norm(trial - test)
            / max(float(np.linalg.norm(trial)), np.finfo(float).tiny)
        ),
        "full_operator_condition": float(np.linalg.cond(matrix)),
        "reduced_operator_condition": float(np.linalg.cond(reduced)),
    }
    return MixedGalerkinOrthogonalization(
        trial_transform=trial,
        test_transform=test,
        reduced_operator=reduced,
        keep_indices=keep,
        eliminate_indices=eliminate,
        diagnostics_data=diagnostics,
    )


@dataclass(frozen=True)
class HybridVIMSystem:
    """Reduced eddy-current VIM matrices on bulk-T plus surface-Omega bases."""

    resistance: np.ndarray
    inductance: object
    surface_mass: np.ndarray
    basis_names: tuple[str, ...]
    blocks: dict[str, tuple[int, int]]
    interaction_backend: str = "unknown"

    @property
    def n_modes(self) -> int:
        return int(self.resistance.shape[0])

    def impedance(self, s, *, surface_impedance=0.0) -> np.ndarray:
        """Return ``R + s L + Zs(s) Ms``."""

        inductance = (
            self.inductance.to_dense()
            if _is_matrix_free_operator(self.inductance)
            else self.inductance
        )
        return (
            self.resistance
            + s * inductance
            + _surface_impedance_term(self.surface_mass, surface_impedance)
        )

    def impedance_operator(self, s, *, surface_impedance=0.0):
        """Return the impedance without materializing HACApK-backed inductance."""

        dense = self.resistance + _surface_impedance_term(
            self.surface_mass,
            surface_impedance,
        )
        if _is_matrix_free_operator(self.inductance):
            return ShiftedReducedOperator(dense, self.inductance, s)
        return dense + s * self.inductance

    def block_slice(self, name: str) -> slice:
        """Return the coefficient slice for a named reduced basis block."""

        if name not in self.blocks:
            known = ", ".join(sorted(self.blocks))
            raise KeyError(f"unknown block {name!r}; known blocks: {known}")
        start, stop = self.blocks[name]
        return slice(start, stop)

    def block_indices(self, blocks) -> np.ndarray:
        """Return concatenated coefficient indices for one or more blocks."""

        if isinstance(blocks, str):
            names = (blocks,)
        else:
            names = tuple(blocks)
        if not names:
            raise ValueError("at least one block name is required")
        if len(set(names)) != len(names):
            raise ValueError("block names must be unique")
        ranges = []
        for name in names:
            block = self.block_slice(name)
            ranges.append(np.arange(block.start, block.stop, dtype=int))
        return np.concatenate(ranges)

    def block_matrix(
        self,
        row_block: str,
        col_block: str,
        s,
        *,
        surface_impedance=0.0,
    ) -> np.ndarray:
        """Return one named block of ``R + s L + Zs M_Gamma``."""

        z = self.impedance(s, surface_impedance=surface_impedance)
        rows = self.block_slice(row_block)
        cols = self.block_slice(col_block)
        return z[rows, cols]

    def block_matrix_blocks(
        self,
        row_blocks,
        col_blocks,
        s,
        *,
        surface_impedance=0.0,
    ) -> np.ndarray:
        """Return a matrix block assembled from multiple named blocks."""

        z = self.impedance(s, surface_impedance=surface_impedance)
        rows = self.block_indices(row_blocks)
        cols = self.block_indices(col_blocks)
        return z[np.ix_(rows, cols)]

    def schur_complement(
        self,
        keep_block: str,
        eliminate_block: str,
        s,
        *,
        surface_impedance=0.0,
    ) -> np.ndarray:
        """Return the named-block Schur complement.

        This is the production hook for the IGTE mixed Galerkin reduction:
        keep the surface/SIBC block and eliminate the finite bulk CLN block to
        obtain the surface-port DtN/admittance block, or do the reverse for a
        bulk-only effective operator.
        """

        z = self.impedance(s, surface_impedance=surface_impedance)
        keep = self.block_slice(keep_block)
        elim = self.block_slice(eliminate_block)
        z_kk = z[keep, keep]
        z_ke = z[keep, elim]
        z_ek = z[elim, keep]
        z_ee = z[elim, elim]
        func = _radia_cpp_kernel("_HybridVIMSchurComplement")
        if func is not None:
            return func(
                np.ascontiguousarray(z_kk, dtype=np.complex128),
                np.ascontiguousarray(z_ke, dtype=np.complex128),
                np.ascontiguousarray(z_ek, dtype=np.complex128),
                np.ascontiguousarray(z_ee, dtype=np.complex128),
            )
        return z_kk - z_ke @ np.linalg.solve(z_ee, z_ek)

    def schur_complement_blocks(
        self,
        keep_blocks,
        eliminate_blocks,
        s,
        *,
        surface_impedance=0.0,
    ) -> np.ndarray:
        """Return a Schur complement over multiple named block groups."""

        z = self.impedance(s, surface_impedance=surface_impedance)
        keep = self.block_indices(keep_blocks)
        elim = self.block_indices(eliminate_blocks)
        if np.intersect1d(keep, elim).size:
            raise ValueError("keep_blocks and eliminate_blocks must not overlap")
        z_kk = z[np.ix_(keep, keep)]
        z_ke = z[np.ix_(keep, elim)]
        z_ek = z[np.ix_(elim, keep)]
        z_ee = z[np.ix_(elim, elim)]
        func = _radia_cpp_kernel("_HybridVIMSchurComplement")
        if func is not None:
            return func(
                np.ascontiguousarray(z_kk, dtype=np.complex128),
                np.ascontiguousarray(z_ke, dtype=np.complex128),
                np.ascontiguousarray(z_ek, dtype=np.complex128),
                np.ascontiguousarray(z_ee, dtype=np.complex128),
            )
        return z_kk - z_ke @ np.linalg.solve(z_ee, z_ek)

    def mixed_galerkin_orthogonalization(
        self,
        keep_blocks,
        eliminate_blocks,
        s,
        *,
        surface_impedance=0.0,
    ) -> MixedGalerkinOrthogonalization:
        """Biorthogonalize kept trial/test bases against eliminated blocks.

        The trial transform enforces ``Z_elim @ T = 0`` and the test transform
        enforces ``W.T @ Z[:, elim] = 0``.  Their mixed Galerkin projection
        ``W.T @ Z @ T`` is exactly the block Schur complement.  For reciprocal
        complex-symmetric VIM operators, ``W == T`` and this becomes ordinary
        block Gram orthogonalization.
        """

        z = self.impedance(s, surface_impedance=surface_impedance)
        keep = self.block_indices(keep_blocks)
        eliminate = self.block_indices(eliminate_blocks)
        if np.intersect1d(keep, eliminate).size:
            raise ValueError("keep_blocks and eliminate_blocks must not overlap")
        if keep.size + eliminate.size != self.n_modes:
            raise ValueError("keep_blocks and eliminate_blocks must cover all modes")

        return _mixed_galerkin_orthogonalize_operator(
            z,
            keep,
            eliminate,
        )

    def block_rhs(self, **block_rhs) -> np.ndarray:
        """Assemble a full reduced right-hand side from named block pieces.

        Each value may be either a vector with length equal to that block size
        or a matrix with shape ``(block_size, n_ports)``.  Missing blocks are
        filled with zeros.  This is the practical glue for mixed bulk/surface
        systems where the volume EVRS ports and the surface-Omega/SIBC ports
        are projected on different quadrature rules.
        """

        if not block_rhs:
            raise ValueError("at least one block RHS must be supplied")

        pieces: dict[str, np.ndarray] = {}
        n_ports: int | None = None
        all_vectors = True
        for name, value in block_rhs.items():
            block = self.block_slice(name)
            block_size = block.stop - block.start
            arr = np.asarray(value)
            if arr.ndim == 1:
                arr = arr[:, np.newaxis]
            elif arr.ndim == 2:
                all_vectors = False
            else:
                raise ValueError(f"RHS for block {name!r} must be 1-D or 2-D")
            if arr.shape[0] != block_size:
                raise ValueError(
                    f"RHS for block {name!r} has {arr.shape[0]} rows; "
                    f"expected {block_size}"
                )
            if n_ports is None:
                n_ports = int(arr.shape[1])
            elif arr.shape[1] != n_ports:
                raise ValueError("all block RHS matrices must have the same column count")
            pieces[name] = arr

        assert n_ports is not None
        dtype = np.result_type(*(piece.dtype for piece in pieces.values()))
        rhs = np.zeros((self.n_modes, n_ports), dtype=dtype)
        for name, piece in pieces.items():
            rhs[self.block_slice(name), :] = piece
        if n_ports == 1 and all_vectors:
            return rhs[:, 0]
        return rhs

    def solve(
        self,
        s,
        rhs,
        *,
        surface_impedance=0.0,
        solver: str = "cocr",
        tolerance: float = 1.0e-10,
        max_iterations: int | None = None,
        restart: int | None = None,
        return_diagnostics: bool = False,
    ) -> np.ndarray:
        """Solve through dense C++ or native NGSolve COCR/GMRES."""

        z = self.impedance_operator(s, surface_impedance=surface_impedance)
        values = np.asarray(rhs)
        if not _is_matrix_free_operator(z):
            result = _solve_reduced_linear(z, values)
            diagnostics = {
                "backend": "native-dense-reduced-lu",
                "iterations": 0,
                "relative_residual_max": 0.0,
            }
        else:
            result, diagnostics = _solve_native_reduced_operator(
                z,
                values,
                solver=solver,
                tolerance=tolerance,
                max_iterations=max_iterations,
                restart=restart,
            )
        return (result, diagnostics) if return_diagnostics else result

    def port_admittance(self, s, rhs, *, surface_impedance=0.0) -> np.ndarray:
        """Return the reduced port admittance ``B^* Z(s)^-1 B``."""

        b = _port_rhs_matrix(rhs, self.n_modes)
        x = self.solve(s, b, surface_impedance=surface_impedance)
        return b.conj().T @ x

    def port_impedance(self, s, rhs, *, surface_impedance=0.0) -> np.ndarray:
        """Return the reduced port impedance, the inverse of port admittance."""

        return np.linalg.inv(
            self.port_admittance(s, rhs, surface_impedance=surface_impedance)
        )

    def diagnostics(self, *, passive_tol: float = 1.0e-10) -> dict[str, object]:
        """Return backend-independent matrix quality diagnostics.

        These checks are the production handoff between dense sampled kernels,
        ngsolve.bem projections, and Radia/HACApK backends: every backend must
        preserve the Hermitian/passive block structure before it is allowed to
        drive CLN fitting or motor-coupled solves.
        """

        if passive_tol < 0.0:
            raise ValueError("passive_tol must be non-negative")
        rmin = _min_hermitian_eigenvalue(self.resistance)
        lmin = _min_hermitian_eigenvalue(self.inductance)
        smin = _min_hermitian_eigenvalue(self.surface_mass)
        matrix_free = _is_matrix_free_operator(self.inductance)
        return {
            "n_modes": self.n_modes,
            "interaction_backend": self.interaction_backend,
            "basis_count": len(self.basis_names),
            "blocks": {name: [start, stop] for name, (start, stop) in self.blocks.items()},
            "resistance_hermitian_error": _relative_hermitian_error(self.resistance),
            "inductance_hermitian_error": _relative_hermitian_error(self.inductance),
            "surface_mass_hermitian_error": _relative_hermitian_error(self.surface_mass),
            "min_resistance_eigenvalue": rmin,
            "min_inductance_eigenvalue": lmin,
            "min_surface_mass_eigenvalue": smin,
            "inductance_matrix_free": matrix_free,
            "inductance_operator": self.inductance.stats() if matrix_free else None,
            "passive_blocks": (
                rmin >= -passive_tol
                and lmin >= -passive_tol
                and smin >= -passive_tol
            ),
        }


@dataclass(frozen=True)
class HCurlEddyCLNModel:
    """Passive CLN/EVRS descriptor produced by an HCurl Eddy Bubble VIM.

    The spatial reduction has already happened in the high-order HCurl parent
    space.  This object is the temporal/frequency coupling contract

    ``(R + s L + Zs(s) M_surface) c = P u``.

    For a coil current ``i`` whose incident vector potential is represented by
    ``P i``, Faraday excitation is ``u = -s i``.  A nonzero SIBC block remains
    a DtN termination and must be rationalized separately before ordinary
    state-space export.
    """

    resistance: np.ndarray
    inductance: object
    surface_mass: np.ndarray
    port_rhs: np.ndarray
    basis_names: tuple[str, ...] = ()
    blocks: dict[str, tuple[int, int]] | None = None

    def __post_init__(self) -> None:
        resistance = np.asarray(self.resistance)
        inductance = np.asarray(self.inductance)
        surface_mass = np.asarray(self.surface_mass)
        if resistance.ndim != 2 or resistance.shape[0] != resistance.shape[1]:
            raise ValueError("resistance must be square")
        if inductance.shape != resistance.shape:
            raise ValueError("inductance must match resistance")
        if surface_mass.shape != resistance.shape:
            raise ValueError("surface_mass must match resistance")
        if not (
            np.all(np.isfinite(resistance))
            and np.all(np.isfinite(inductance))
            and np.all(np.isfinite(surface_mass))
        ):
            raise ValueError("CLN matrices must contain only finite values")
        ports = _port_rhs_matrix(self.port_rhs, resistance.shape[0])
        names = _mode_names(self.basis_names or None, resistance.shape[0])
        blocks = {} if self.blocks is None else dict(self.blocks)
        for name, bounds in blocks.items():
            if len(bounds) != 2:
                raise ValueError(f"block {name!r} must contain (start, stop)")
            start, stop = (int(bounds[0]), int(bounds[1]))
            if start < 0 or stop < start or stop > resistance.shape[0]:
                raise ValueError(f"block {name!r} is out of range")
            blocks[name] = (start, stop)
        object.__setattr__(self, "resistance", np.array(resistance, copy=True))
        object.__setattr__(self, "inductance", np.array(inductance, copy=True))
        object.__setattr__(self, "surface_mass", np.array(surface_mass, copy=True))
        object.__setattr__(self, "port_rhs", np.array(ports, copy=True))
        object.__setattr__(self, "basis_names", names)
        object.__setattr__(self, "blocks", blocks)

    @property
    def state_order(self) -> int:
        return int(self.resistance.shape[0])

    @property
    def port_count(self) -> int:
        return int(self.port_rhs.shape[1])

    @property
    def has_sibc_termination(self) -> bool:
        return bool(np.linalg.norm(self.surface_mass) > 1.0e-14)

    def impedance(self, s, *, surface_impedance=0.0) -> np.ndarray:
        s = complex(s)
        if not (np.isfinite(s.real) and np.isfinite(s.imag)):
            raise ValueError("s must be finite")
        return (
            self.resistance
            + s * self.inductance
            + _surface_impedance_term(self.surface_mass, surface_impedance)
        )

    def solve(self, s, drive, *, surface_impedance=0.0) -> np.ndarray:
        """Solve for reduced current coefficients under generalized drive."""

        values = np.asarray(drive)
        vector_drive = values.ndim <= 1
        if values.ndim == 0:
            if self.port_count != 1:
                raise ValueError(
                    f"scalar drive is valid only for one port, not {self.port_count}"
                )
            values = values.reshape(1, 1)
        elif values.ndim == 1:
            values = values[:, np.newaxis]
        if values.ndim != 2 or values.shape[0] != self.port_count:
            raise ValueError(f"drive must have shape ({self.port_count}, n_cases)")
        solved = _solve_reduced_linear(
            self.impedance(s, surface_impedance=surface_impedance),
            self.port_rhs @ values,
        )
        return solved[:, 0] if vector_drive else solved

    def solve_vector_potential_drive(
        self,
        s,
        coil_current,
        *,
        surface_impedance=0.0,
    ) -> np.ndarray:
        """Solve ``(R+sL+ZsMs)c = -s P i`` for coil current ``i``."""

        return self.solve(
            s,
            -complex(s) * np.asarray(coil_current),
            surface_impedance=surface_impedance,
        )

    def port_admittance(self, s, *, surface_impedance=0.0) -> np.ndarray:
        z = self.impedance(s, surface_impedance=surface_impedance)
        response = _solve_reduced_linear(z, self.port_rhs)
        return self.port_rhs.conj().T @ response

    def derivative_input_state_space(self) -> dict[str, np.ndarray]:
        """Return ``c_dot=A c+B u``, ``y=C c`` for ``u=-i_dot``."""

        if self.has_sibc_termination:
            raise ValueError(
                "SIBC termination must be rationalized before state-space export"
            )
        a = -np.linalg.solve(self.inductance, self.resistance)
        b = np.linalg.solve(self.inductance, self.port_rhs)
        c = self.port_rhs.conj().T
        d = np.zeros(
            (self.port_count, self.port_count),
            dtype=np.result_type(a, b, c),
        )
        return {"A": a, "B": b, "C": c, "D": d}

    def diagnostics(self, *, passive_tol: float = 1.0e-10) -> dict[str, object]:
        rmin = _min_hermitian_eigenvalue(self.resistance)
        lmin = _min_hermitian_eigenvalue(self.inductance)
        smin = _min_hermitian_eigenvalue(self.surface_mass)
        return {
            "state_order": self.state_order,
            "port_count": self.port_count,
            "blocks": {
                name: [start, stop]
                for name, (start, stop) in (self.blocks or {}).items()
            },
            "input_convention": "u=-d(coil_current)/dt for vector-potential ports",
            "has_sibc_termination": self.has_sibc_termination,
            "finite_rl_state_space": not self.has_sibc_termination,
            "min_resistance_eigenvalue": rmin,
            "min_inductance_eigenvalue": lmin,
            "min_surface_mass_eigenvalue": smin,
            "passive": (
                rmin >= -passive_tol
                and lmin >= -passive_tol
                and smin >= -passive_tol
            ),
        }


def HCurlEddyCLNFromVIM(system: HybridVIMSystem, rhs) -> HCurlEddyCLNModel:
    """Create the CLN coupling descriptor from an HCurl Eddy Bubble VIM."""

    if not isinstance(system, HybridVIMSystem):
        raise TypeError("system must be a HybridVIMSystem")
    return HCurlEddyCLNModel(
        resistance=system.resistance,
        inductance=system.inductance,
        surface_mass=system.surface_mass,
        port_rhs=rhs,
        basis_names=system.basis_names,
        blocks=system.blocks,
    )


def ReducedPortAdmittance(
    system: HybridVIMSystem,
    s,
    rhs,
    *,
    surface_impedance=0.0,
) -> np.ndarray:
    """Evaluate ``B^* Z(s)^-1 B`` for a reduced eddy-current VIM system."""

    if not isinstance(system, HybridVIMSystem):
        raise TypeError("system must be a HybridVIMSystem")
    return system.port_admittance(s, rhs, surface_impedance=surface_impedance)


def ReducedPortImpedance(
    system: HybridVIMSystem,
    s,
    rhs,
    *,
    surface_impedance=0.0,
) -> np.ndarray:
    """Evaluate the inverse port admittance for a reduced VIM system."""

    if not isinstance(system, HybridVIMSystem):
        raise TypeError("system must be a HybridVIMSystem")
    return system.port_impedance(s, rhs, surface_impedance=surface_impedance)


@dataclass(frozen=True)
class SharedMeshMaterialModel:
    """Shared mesh/material registry for HCurl-VIM and HDiv-MMM branches.

    The object deliberately stores coefficients without interpreting their
    backend type.  Scalars, per-region dictionaries, and NGSolve coefficient
    functions can all be carried here; assembly code decides how to consume
    them.  The important production rule is that both branches receive this
    same registry instead of carrying independent permeability/conductivity
    copies.
    """

    mesh: object
    magnetic_regions: tuple[str, ...] | str | Iterable[str] = ()
    conductive_regions: tuple[str, ...] | str | Iterable[str] = ()
    mu: object | None = None
    nu: object | None = None
    sigma: object | None = None
    sibc: object | None = None
    metadata: object | None = None

    def __post_init__(self) -> None:
        if self.mesh is None:
            raise ValueError("mesh must not be None")
        object.__setattr__(
            self,
            "magnetic_regions",
            _label_tuple(self.magnetic_regions, "magnetic_regions"),
        )
        object.__setattr__(
            self,
            "conductive_regions",
            _label_tuple(self.conductive_regions, "conductive_regions"),
        )
        _check_positive_material_coeff(self.mu, "mu")
        _check_positive_material_coeff(self.nu, "nu")
        _check_positive_material_coeff(self.sigma, "sigma")

    @property
    def has_permeability(self) -> bool:
        return self.mu is not None

    @property
    def has_reluctivity(self) -> bool:
        return self.nu is not None

    @property
    def has_magnetic_law(self) -> bool:
        return self.has_permeability or self.has_reluctivity

    @property
    def has_conductivity(self) -> bool:
        return self.sigma is not None

    @property
    def has_sibc(self) -> bool:
        return self.sibc is not None

    def hdiv_mmm_coefficient(self):
        """Return the magnetic coefficient carried for the HDiv-MMM branch."""

        return self.nu if self.nu is not None else self.mu

    def hcurl_vim_coefficient(self):
        """Return the conductive/SIBC coefficient carried for the HCurl-VIM branch."""

        return self.sibc if self.sibc is not None else self.sigma

    def diagnostics(self) -> dict[str, int | bool]:
        """Return branch-visibility diagnostics for the shared registry."""

        return {
            "magnetic_region_count": len(self.magnetic_regions),
            "conductive_region_count": len(self.conductive_regions),
            "has_permeability": self.has_permeability,
            "has_reluctivity": self.has_reluctivity,
            "has_magnetic_law": self.has_magnetic_law,
            "has_conductivity": self.has_conductivity,
            "has_sibc": self.has_sibc,
        }


@dataclass(frozen=True)
class HCurlVIMHDivMMMSolution:
    """One-frequency result of the mixed HDiv-MMM / HCurl-VIM solve.

    Coefficients are retained together with physical sampled fields.  The
    volume EVRS coefficients are also lifted back to the high-order parent
    HCurl ``T`` coordinates when the response basis is available.
    """

    frequency_hz: float
    s: complex
    surface_impedance: object
    magnetization_coefficients: np.ndarray
    eddy_coefficients: np.ndarray
    parent_t_coefficients: np.ndarray | None
    parent_magnetization_coefficients: np.ndarray | None
    sampled_magnetization: np.ndarray
    eddy_block_names: tuple[str, ...]
    sampled_eddy_currents: tuple[np.ndarray, ...]
    eddy_sample_points: tuple[np.ndarray, ...]
    eddy_sample_weights: tuple[np.ndarray, ...]
    reduced_rhs: np.ndarray
    reduced_solution: np.ndarray
    port_response: np.ndarray
    average_joule_loss: np.ndarray
    residual_relative_norm: float
    solver_backend: str
    solver_diagnostics: dict[str, object] | None = None
    orthogonalized_rhs: np.ndarray | None = None
    orthogonalized_solution: np.ndarray | None = None
    mixed_galerkin_diagnostics: dict[str, object] | None = None
    residual_backward_error: float | None = None
    eddy_block_roles: tuple[str, ...] | None = None

    @property
    def n_excitations(self) -> int:
        return int(self.reduced_solution.shape[1])

    def _current_block_index(self, block: str) -> int:
        try:
            return self.eddy_block_names.index(block)
        except ValueError:
            if self.eddy_block_roles is not None:
                try:
                    return self.eddy_block_roles.index(block)
                except ValueError:
                    pass
            known = ", ".join(self.eddy_block_names)
            roles = (
                ""
                if self.eddy_block_roles is None
                else "; adjacency roles: " + ", ".join(self.eddy_block_roles)
            )
            raise KeyError(
                f"unknown eddy block {block!r}; known blocks: {known}{roles}"
            ) from None

    def current_samples(self, block: str) -> np.ndarray:
        """Return current samples shaped ``(n_excitation, n_point, 3)``."""

        return self.sampled_eddy_currents[self._current_block_index(block)]

    def eddy_flux_density(
        self,
        target_points,
        *,
        block: str | None = None,
        mu: float = MU0,
        kernel_epsilon: float = 0.0,
    ) -> np.ndarray:
        """Evaluate the retained eddy-current ``B`` field at target points.

        The return shape is ``(n_excitation, n_target, 3)``.  With ``block=None``
        the bulk EVRS, conductor-cycle bridge, and surface-Omega/SIBC fields are
        summed.  Passing a block name isolates that contribution.  This is the
        same quasi-static Biot--Savart reconstruction used to assemble the
        HDiv-MMM coupling.
        """

        if block is None:
            indices = range(len(self.eddy_block_names))
        else:
            indices = (self._current_block_index(block),)
        targets = _as_points(target_points, "target_points")
        total = np.zeros(
            (self.n_excitations, targets.shape[0], 3),
            dtype=self.reduced_solution.dtype,
        )
        for index in indices:
            currents = self.sampled_eddy_currents[index]
            role = (
                None
                if self.eddy_block_roles is None
                else self.eddy_block_roles[index]
            )
            sampled = SampledCurrentBasis(
                points=self.eddy_sample_points[index],
                weights=self.eddy_sample_weights[index],
                modes=currents,
                kind=(
                    "surface"
                    if role in ("sibc", "non_sibc_trace")
                    or self.eddy_block_names[index] == "surface"
                    else "volume"
                ),
                names=tuple(
                    f"excitation_{excitation}"
                    for excitation in range(self.n_excitations)
                ),
            )
            total += CurrentMagneticFluxDensitySamples(
                sampled,
                targets,
                mu=mu,
                kernel_epsilon=kernel_epsilon,
            )
        return total

    def diagnostics(self) -> dict[str, object]:
        """Return compact JSON-ready solve and reconstruction diagnostics."""

        if isinstance(self.surface_impedance, SurfaceImpedanceGram):
            surface_impedance_diagnostics = self.surface_impedance.diagnostics()
        else:
            value = complex(self.surface_impedance)
            surface_impedance_diagnostics = {
                "kind": "uniform",
                "real": float(value.real),
                "imag": float(value.imag),
            }

        return {
            "frequency_Hz": float(self.frequency_hz),
            "s": {"real": float(self.s.real), "imag": float(self.s.imag)},
            "surface_impedance": surface_impedance_diagnostics,
            "n_excitations": self.n_excitations,
            "magnetization_modes": int(self.magnetization_coefficients.shape[0]),
            "eddy_modes": int(self.eddy_coefficients.shape[0]),
            "parent_t_dofs": (
                None if self.parent_t_coefficients is None
                else int(self.parent_t_coefficients.shape[0])
            ),
            "parent_hdiv_dofs": (
                None if self.parent_magnetization_coefficients is None
                else int(self.parent_magnetization_coefficients.shape[0])
            ),
            "eddy_blocks": list(self.eddy_block_names),
            "eddy_block_roles": (
                None
                if self.eddy_block_roles is None
                else list(self.eddy_block_roles)
            ),
            "eddy_block_sample_counts": [
                int(values.shape[1]) for values in self.sampled_eddy_currents
            ],
            "sampled_magnetization_norm": float(np.linalg.norm(self.sampled_magnetization)),
            "sampled_eddy_current_norms": [
                float(np.linalg.norm(values)) for values in self.sampled_eddy_currents
            ],
            "average_joule_loss": [float(value) for value in self.average_joule_loss],
            "residual_relative_norm": float(self.residual_relative_norm),
            "residual_backward_error": (
                None
                if self.residual_backward_error is None
                else float(self.residual_backward_error)
            ),
            "solver_backend": self.solver_backend,
            "solver_diagnostics": (
                None
                if self.solver_diagnostics is None
                else dict(self.solver_diagnostics)
            ),
            "orthogonalized_modes": (
                None
                if self.orthogonalized_solution is None
                else int(self.orthogonalized_solution.shape[0])
            ),
            "mixed_galerkin": (
                None
                if self.mixed_galerkin_diagnostics is None
                else dict(self.mixed_galerkin_diagnostics)
            ),
        }


@dataclass(frozen=True)
class CoupledHDivEVRSSystem:
    """Rectangular HDiv-MMM / HCurl-VIM eddy-current coupling.

    This is the named container for the mixed Radia idea: HDiv-MMM owns the
    magnetic-material modes, while HCurl-VIM owns the eddy branch through a
    response-compressed HCurl parent space sampled as ``J = curl(T)`` currents.
    The coupling block is ``C_ij = int M_i dot B[J_j] dV``.

    ``material_model`` is optional metadata for the shared-mesh material
    registry used by both branches.  The material law is not folded into the
    coupling block itself: the full material operator is built first on the
    parent mesh/regions, then projected to the retained HCurl/HDiv coordinates.
    """

    magnetization_basis: SampledMagnetizationBasis
    eddy_basis: SampledCurrentBasis
    coupling: np.ndarray
    eddy_system: HybridVIMSystem | None = None
    material_model: object | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.magnetization_basis, SampledMagnetizationBasis):
            raise TypeError("magnetization_basis must be a SampledMagnetizationBasis")
        if not isinstance(self.eddy_basis, SampledCurrentBasis):
            raise TypeError("eddy_basis must be a SampledCurrentBasis")
        if self.eddy_system is not None and not isinstance(self.eddy_system, HybridVIMSystem):
            raise TypeError("eddy_system must be a HybridVIMSystem or None")
        coupling = np.asarray(self.coupling)
        expected = (self.magnetization_basis.n_modes, self.eddy_basis.n_modes)
        if coupling.shape != expected:
            raise ValueError(f"coupling must have shape {expected}")
        if not np.all(np.isfinite(coupling)):
            raise ValueError("coupling contains non-finite values")
        object.__setattr__(self, "coupling", coupling)

    @property
    def n_hdiv_modes(self) -> int:
        return self.magnetization_basis.n_modes

    @property
    def n_evrs_modes(self) -> int:
        return self.eddy_basis.n_modes

    @property
    def n_hdiv_mmm_modes(self) -> int:
        return self.n_hdiv_modes

    @property
    def n_hcurl_vim_modes(self) -> int:
        return self.n_evrs_modes

    def mixed_energy(self, magnetization_coeffs, eddy_coeffs):
        """Return ``m^* C j`` for the HDiv-MMM / HCurl-VIM coupling block."""

        m = np.asarray(magnetization_coeffs)
        j = np.asarray(eddy_coeffs)
        if m.shape != (self.n_hdiv_modes,):
            raise ValueError(f"magnetization_coeffs must have shape ({self.n_hdiv_modes},)")
        if j.shape != (self.n_evrs_modes,):
            raise ValueError(f"eddy_coeffs must have shape ({self.n_evrs_modes},)")
        return m.conj() @ self.coupling @ j

    def eddy_impedance(self, s, *, surface_impedance=0.0, eddy_operator=None) -> np.ndarray:
        """Return the HCurl-VIM eddy impedance block.

        Pass ``eddy_operator`` to use an already assembled matrix.  Otherwise
        the optional ``eddy_system`` is evaluated as ``R + s L + Zs M``.
        """

        if eddy_operator is not None:
            return _square_matrix(eddy_operator, self.n_evrs_modes, "eddy_operator")
        if self.eddy_system is None:
            raise ValueError("eddy_system is required unless eddy_operator is supplied")
        return self.eddy_system.impedance(s, surface_impedance=surface_impedance)

    def mixed_operator(
        self,
        magnetic_operator,
        s=None,
        *,
        eddy_operator=None,
        surface_impedance=0.0,
        coupling_scale=1.0,
        adjoint_coupling_scale=None,
    ) -> np.ndarray:
        """Assemble the coupled HDiv-MMM / HCurl-VIM block operator.

        The returned matrix has the block form

        ``[[A_m, alpha K], [beta K^*, Z_e]]``

        where ``A_m`` is the HDiv-MMM magnetic operator, ``Z_e`` is the
        HCurl-VIM eddy impedance, and ``K`` is the rectangular coupling block.
        By default ``beta = conj(alpha)`` so the off-diagonal coupling is
        reciprocal when the diagonal blocks are Hermitian.
        """

        magnetic = _square_matrix(magnetic_operator, self.n_hdiv_modes, "magnetic_operator")
        eddy = self.eddy_impedance(
            s,
            surface_impedance=surface_impedance,
            eddy_operator=eddy_operator,
        )
        if adjoint_coupling_scale is None:
            adjoint_coupling_scale = np.conjugate(coupling_scale)
        upper = coupling_scale * self.coupling
        lower = adjoint_coupling_scale * self.coupling.conj().T
        dtype = np.result_type(magnetic, eddy, upper, lower)
        out = np.zeros(
            (self.n_hdiv_modes + self.n_evrs_modes,
             self.n_hdiv_modes + self.n_evrs_modes),
            dtype=dtype,
        )
        mh = self.n_hdiv_modes
        out[:mh, :mh] = magnetic
        out[:mh, mh:] = upper
        out[mh:, :mh] = lower
        out[mh:, mh:] = eddy
        return out

    def solve(
        self,
        magnetic_operator,
        s=None,
        *,
        magnetic_rhs=None,
        eddy_rhs=None,
        eddy_operator=None,
        surface_impedance=0.0,
        coupling_scale=1.0,
        adjoint_coupling_scale=None,
        return_operator: bool = False,
    ) -> dict[str, np.ndarray]:
        """Solve the mixed HDiv-MMM / HCurl-VIM linear system."""

        op = self.mixed_operator(
            magnetic_operator,
            s,
            eddy_operator=eddy_operator,
            surface_impedance=surface_impedance,
            coupling_scale=coupling_scale,
            adjoint_coupling_scale=adjoint_coupling_scale,
        )
        mh = self.n_hdiv_modes
        ne = self.n_evrs_modes
        if magnetic_rhs is None and eddy_rhs is None:
            m_rhs = np.zeros((mh, 1), dtype=op.dtype)
            e_rhs = np.zeros((ne, 1), dtype=op.dtype)
        elif magnetic_rhs is None:
            e_rhs = _port_rhs_matrix(eddy_rhs, ne)
            m_rhs = np.zeros((mh, e_rhs.shape[1]), dtype=op.dtype)
        elif eddy_rhs is None:
            m_rhs = _port_rhs_matrix(magnetic_rhs, mh)
            e_rhs = np.zeros((ne, m_rhs.shape[1]), dtype=op.dtype)
        else:
            m_rhs = _port_rhs_matrix(magnetic_rhs, mh)
            e_rhs = _port_rhs_matrix(eddy_rhs, ne)
        if m_rhs.shape[1] != e_rhs.shape[1]:
            raise ValueError("magnetic_rhs and eddy_rhs must have the same number of columns")
        rhs = np.vstack([m_rhs, e_rhs]).astype(op.dtype, copy=False)
        sol = np.linalg.solve(op, rhs)
        result = {
            "magnetization": sol[:mh],
            "eddy": sol[mh:],
            "solution": sol,
        }
        if return_operator:
            result["operator"] = op
        return result

    def schur_magnetic_operator(
        self,
        magnetic_operator,
        s=None,
        *,
        eddy_operator=None,
        surface_impedance=0.0,
        coupling_scale=1.0,
        adjoint_coupling_scale=None,
    ) -> np.ndarray:
        """Return the Schur complement after eliminating HCurl-VIM eddy unknowns."""

        magnetic = _square_matrix(magnetic_operator, self.n_hdiv_modes, "magnetic_operator")
        eddy = self.eddy_impedance(
            s,
            surface_impedance=surface_impedance,
            eddy_operator=eddy_operator,
        )
        if adjoint_coupling_scale is None:
            adjoint_coupling_scale = np.conjugate(coupling_scale)
        upper = coupling_scale * self.coupling
        lower = adjoint_coupling_scale * self.coupling.conj().T
        return magnetic - upper @ np.linalg.solve(eddy, lower)

    def schur_eddy_operator(
        self,
        magnetic_operator,
        s=None,
        *,
        eddy_operator=None,
        surface_impedance=0.0,
        coupling_scale=1.0,
        adjoint_coupling_scale=None,
    ) -> np.ndarray:
        """Return the Schur complement after eliminating HDiv-MMM unknowns."""

        magnetic = _square_matrix(magnetic_operator, self.n_hdiv_modes, "magnetic_operator")
        eddy = self.eddy_impedance(
            s,
            surface_impedance=surface_impedance,
            eddy_operator=eddy_operator,
        )
        if adjoint_coupling_scale is None:
            adjoint_coupling_scale = np.conjugate(coupling_scale)
        upper = coupling_scale * self.coupling
        lower = adjoint_coupling_scale * self.coupling.conj().T
        return eddy - lower @ np.linalg.solve(magnetic, upper)

    def port_admittance(
        self,
        magnetic_operator,
        s,
        eddy_rhs,
        *,
        magnetic_rhs=None,
        eddy_operator=None,
        surface_impedance=0.0,
        coupling_scale=1.0,
        adjoint_coupling_scale=None,
    ) -> np.ndarray:
        """Return mixed-system port admittance for magnetic and eddy RHS columns."""

        mh = self.n_hdiv_modes
        ne = self.n_evrs_modes
        e_ports = _port_rhs_matrix(eddy_rhs, ne)
        if magnetic_rhs is None:
            m_ports = np.zeros((mh, e_ports.shape[1]), dtype=e_ports.dtype)
        else:
            m_ports = _port_rhs_matrix(magnetic_rhs, mh)
        if m_ports.shape[1] != e_ports.shape[1]:
            raise ValueError("magnetic_rhs and eddy_rhs must have the same number of columns")
        sol = self.solve(
            magnetic_operator,
            s,
            magnetic_rhs=m_ports,
            eddy_rhs=e_ports,
            eddy_operator=eddy_operator,
            surface_impedance=surface_impedance,
            coupling_scale=coupling_scale,
            adjoint_coupling_scale=adjoint_coupling_scale,
        )["solution"]
        ports = np.vstack([m_ports, e_ports])
        return ports.conj().T @ sol

    def diagnostics(self) -> dict[str, int | float | bool]:
        """Return compact mixed-system diagnostics."""

        return {
            "hdiv_modes": self.n_hdiv_modes,
            "hdiv_mmm_modes": self.n_hdiv_mmm_modes,
            "evrs_modes": self.n_evrs_modes,
            "hcurl_vim_modes": self.n_hcurl_vim_modes,
            "coupling_rows": int(self.coupling.shape[0]),
            "coupling_cols": int(self.coupling.shape[1]),
            "coupling_frobenius_norm": float(np.linalg.norm(self.coupling)),
            "has_eddy_system": self.eddy_system is not None,
            "has_shared_material_model": self.material_model is not None,
            "has_shared_mesh_material_model": isinstance(self.material_model, SharedMeshMaterialModel),
        }


@dataclass(frozen=True)
class CoupledHDivHybridVIMSystem:
    """HDiv-MMM coupled to a full hybrid HCurl-VIM eddy system.

    Unlike :class:`CoupledHDivEVRSSystem`, this container supports several
    current bases inside one eddy system, such as bulk EVRS, bridge-cycle
    currents, and surface-Omega/SIBC currents.  The coupling columns follow the
    mode order of the supplied :class:`HybridVIMSystem`.
    """

    magnetization_basis: SampledMagnetizationBasis
    eddy_system: HybridVIMSystem
    eddy_bases: tuple[SampledCurrentBasis, ...]
    coupling: np.ndarray
    material_model: object | None = None
    magnetic_operator: np.ndarray | None = None
    magnetic_rhs: np.ndarray | None = None
    eddy_rhs: np.ndarray | None = None
    response_basis: ResponseBasis | None = None
    eddy_bubbling: EddyBubbleDecomposition | None = None
    hdiv_reduction: HDivMMMReducedModel | None = None
    conductivity: float | None = None
    eddy_block_roles: dict[str, str] | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.magnetization_basis, SampledMagnetizationBasis):
            raise TypeError("magnetization_basis must be a SampledMagnetizationBasis")
        if not isinstance(self.eddy_system, HybridVIMSystem):
            raise TypeError("eddy_system must be a HybridVIMSystem")
        bases = tuple(self.eddy_bases)
        if not bases:
            raise ValueError("eddy_bases must not be empty")
        if not all(isinstance(basis, SampledCurrentBasis) for basis in bases):
            raise TypeError("all eddy_bases must be SampledCurrentBasis instances")
        if sum(basis.n_modes for basis in bases) != self.eddy_system.n_modes:
            raise ValueError("eddy_bases mode count must match eddy_system.n_modes")
        coupling = np.asarray(self.coupling)
        expected = (self.magnetization_basis.n_modes, self.eddy_system.n_modes)
        if coupling.shape != expected:
            raise ValueError(f"coupling must have shape {expected}")
        if not np.all(np.isfinite(coupling)):
            raise ValueError("coupling contains non-finite values")
        magnetic_operator = self.magnetic_operator
        if magnetic_operator is not None:
            if _is_matrix_free_operator(magnetic_operator):
                expected_magnetic = (
                    self.magnetization_basis.n_modes,
                    self.magnetization_basis.n_modes,
                )
                if tuple(magnetic_operator.shape) != expected_magnetic:
                    raise ValueError(
                        f"magnetic_operator must have shape {expected_magnetic}"
                    )
            else:
                magnetic_operator = _square_matrix(
                    magnetic_operator,
                    self.magnetization_basis.n_modes,
                    "magnetic_operator",
                )
        magnetic_rhs = self.magnetic_rhs
        if magnetic_rhs is not None:
            magnetic_rhs = _port_rhs_matrix(
                magnetic_rhs,
                self.magnetization_basis.n_modes,
            )
        eddy_rhs = self.eddy_rhs
        if eddy_rhs is not None:
            eddy_rhs = _port_rhs_matrix(eddy_rhs, self.eddy_system.n_modes)
        if self.response_basis is not None and not isinstance(self.response_basis, ResponseBasis):
            raise TypeError("response_basis must be a ResponseBasis or None")
        if self.eddy_bubbling is not None and not isinstance(
            self.eddy_bubbling,
            EddyBubbleDecomposition,
        ):
            raise TypeError("eddy_bubbling must be an EddyBubbleDecomposition or None")
        if self.hdiv_reduction is not None and not isinstance(
            self.hdiv_reduction,
            HDivMMMReducedModel,
        ):
            raise TypeError("hdiv_reduction must be an HDivMMMReducedModel or None")
        if (
            self.hdiv_reduction is not None
            and self.hdiv_reduction.n_modes != self.magnetization_basis.n_modes
        ):
            raise ValueError(
                "hdiv_reduction modes must match magnetization_basis modes"
            )
        if self.conductivity is not None and self.conductivity <= 0.0:
            raise ValueError("conductivity must be positive")
        roles = None
        if self.eddy_block_roles is not None:
            roles = {str(name): str(role) for name, role in self.eddy_block_roles.items()}
            if set(roles) != set(self.eddy_system.blocks):
                raise ValueError(
                    "eddy_block_roles must classify every HCurl-VIM block"
                )
            allowed_roles = {"bulk", "bridge", "sibc", "non_sibc_trace"}
            unknown = set(roles.values()) - allowed_roles
            if unknown:
                raise ValueError(
                    "unknown eddy block roles: " + ", ".join(sorted(unknown))
                )
        object.__setattr__(self, "eddy_bases", bases)
        object.__setattr__(self, "coupling", coupling)
        object.__setattr__(self, "magnetic_operator", magnetic_operator)
        object.__setattr__(self, "magnetic_rhs", magnetic_rhs)
        object.__setattr__(self, "eddy_rhs", eddy_rhs)
        object.__setattr__(self, "eddy_block_roles", roles)

    @property
    def n_hdiv_modes(self) -> int:
        return self.magnetization_basis.n_modes

    @property
    def n_hdiv_mmm_modes(self) -> int:
        return self.n_hdiv_modes

    @property
    def n_hcurl_vim_modes(self) -> int:
        return self.eddy_system.n_modes

    @property
    def n_evrs_modes(self) -> int:
        return self.n_hcurl_vim_modes

    def _resolved_magnetic_operator(self, magnetic_operator):
        if magnetic_operator is None:
            magnetic_operator = self.magnetic_operator
        if magnetic_operator is None:
            raise ValueError(
                "magnetic_operator is required; pass it to the builder or solve"
            )
        if _is_matrix_free_operator(magnetic_operator):
            expected = (self.n_hdiv_modes, self.n_hdiv_modes)
            if tuple(magnetic_operator.shape) != expected:
                raise ValueError(f"magnetic_operator must have shape {expected}")
            return magnetic_operator
        return _square_matrix(
            magnetic_operator, self.n_hdiv_modes, "magnetic_operator"
        )

    def _resolved_rhs(self, magnetic_rhs, eddy_rhs, dtype, *, require_excitation=False):
        if magnetic_rhs is None:
            magnetic_rhs = self.magnetic_rhs
        if eddy_rhs is None:
            eddy_rhs = self.eddy_rhs
        if require_excitation and magnetic_rhs is None and eddy_rhs is None:
            raise ValueError(
                "no excitation RHS is stored; pass magnetic_rhs or eddy_rhs"
            )

        mh = self.n_hdiv_modes
        ne = self.n_hcurl_vim_modes
        if magnetic_rhs is None and eddy_rhs is None:
            m_rhs = np.zeros((mh, 1), dtype=dtype)
            e_rhs = np.zeros((ne, 1), dtype=dtype)
        elif magnetic_rhs is None:
            e_rhs = _port_rhs_matrix(eddy_rhs, ne)
            m_rhs = np.zeros((mh, e_rhs.shape[1]), dtype=dtype)
        elif eddy_rhs is None:
            m_rhs = _port_rhs_matrix(magnetic_rhs, mh)
            e_rhs = np.zeros((ne, m_rhs.shape[1]), dtype=dtype)
        else:
            m_rhs = _port_rhs_matrix(magnetic_rhs, mh)
            e_rhs = _port_rhs_matrix(eddy_rhs, ne)
        if m_rhs.shape[1] != e_rhs.shape[1]:
            raise ValueError("magnetic_rhs and eddy_rhs must have the same number of columns")
        return (
            np.asarray(m_rhs, dtype=dtype),
            np.asarray(e_rhs, dtype=dtype),
        )

    def eddy_impedance(self, s, *, surface_impedance=0.0, eddy_operator=None) -> np.ndarray:
        """Return the full hybrid HCurl-VIM eddy impedance block."""

        if eddy_operator is not None:
            return _square_matrix(eddy_operator, self.n_hcurl_vim_modes, "eddy_operator")
        return self.eddy_system.impedance(s, surface_impedance=surface_impedance)

    def eddy_impedance_operator(
        self,
        s,
        *,
        surface_impedance=0.0,
        eddy_operator=None,
    ):
        """Return the HCurl block while retaining any HACApK BaseMatrix terms."""

        if eddy_operator is not None:
            return _square_matrix(
                eddy_operator,
                self.n_hcurl_vim_modes,
                "eddy_operator",
            )
        return self.eddy_system.impedance_operator(
            s,
            surface_impedance=surface_impedance,
        )

    def mixed_operator(
        self,
        magnetic_operator=None,
        s=None,
        *,
        eddy_operator=None,
        surface_impedance=0.0,
        coupling_scale=1.0,
        adjoint_coupling_scale=None,
    ):
        """Assemble ``[[A_m, alpha K], [beta K^*, Z_e]]`` matrix-free."""

        magnetic = self._resolved_magnetic_operator(magnetic_operator)
        eddy = self.eddy_impedance_operator(
            s,
            surface_impedance=surface_impedance,
            eddy_operator=eddy_operator,
        )
        if adjoint_coupling_scale is None:
            adjoint_coupling_scale = np.conjugate(coupling_scale)
        upper = coupling_scale * self.coupling
        lower = adjoint_coupling_scale * self.coupling.conj().T
        if _is_matrix_free_operator(magnetic) or _is_matrix_free_operator(eddy):
            if not _is_matrix_free_operator(eddy):
                eddy = ReducedBlockHMatrixOperator(eddy, ())
            return CoupledReducedOperator(magnetic, upper, lower, eddy)
        dtype = np.result_type(magnetic, eddy, upper, lower)
        out = np.zeros(
            (self.n_hdiv_modes + self.n_hcurl_vim_modes,
             self.n_hdiv_modes + self.n_hcurl_vim_modes),
            dtype=dtype,
        )
        mh = self.n_hdiv_modes
        out[:mh, :mh] = magnetic
        out[:mh, mh:] = upper
        out[mh:, :mh] = lower
        out[mh:, mh:] = eddy
        return out

    def mixed_galerkin_orthogonalization(
        self,
        keep_eddy_blocks,
        eliminate_eddy_blocks,
        magnetic_operator=None,
        s=None,
        *,
        eddy_operator=None,
        surface_impedance=0.0,
        coupling_scale=1.0,
        adjoint_coupling_scale=None,
    ) -> "MixedGalerkinHDivHybridVIMSystem":
        """Eliminate HCurl blocks in the complete HDiv/HCurl operator.

        HDiv-MMM modes are always retained.  Consequently, the eliminated
        HCurl correction contains both the retained eddy-current response and
        the magnetic coupling response.  The resulting Galerkin matrix is the
        exact Schur complement of the full coupled operator.
        """

        full_operator = self.mixed_operator(
            magnetic_operator,
            s,
            eddy_operator=eddy_operator,
            surface_impedance=surface_impedance,
            coupling_scale=coupling_scale,
            adjoint_coupling_scale=adjoint_coupling_scale,
        )
        if _is_matrix_free_operator(full_operator):
            full_operator = full_operator.to_dense()
        eddy_keep = self.eddy_system.block_indices(keep_eddy_blocks)
        eddy_eliminate = self.eddy_system.block_indices(eliminate_eddy_blocks)
        mh = self.n_hdiv_modes
        keep = np.concatenate(
            (np.arange(mh, dtype=int), mh + eddy_keep),
        )
        eliminate = mh + eddy_eliminate
        orthogonalization = _mixed_galerkin_orthogonalize_operator(
            full_operator,
            keep,
            eliminate,
        )

        def block_names(value) -> tuple[str, ...]:
            return (value,) if isinstance(value, str) else tuple(value)

        return MixedGalerkinHDivHybridVIMSystem(
            parent_system=self,
            full_operator=full_operator,
            orthogonalization=orthogonalization,
            keep_eddy_blocks=block_names(keep_eddy_blocks),
            eliminate_eddy_blocks=block_names(eliminate_eddy_blocks),
        )

    def adjacency_class_block_partition(self) -> tuple[tuple[str, ...], tuple[str, ...]]:
        """Return retained/eliminated blocks from neighboring-material roles."""

        if self.eddy_block_roles is None:
            raise ValueError(
                "eddy_block_roles are required for adjacency-class eddy bubbling"
            )
        ordered = tuple(
            name
            for name, _ in sorted(
                self.eddy_system.blocks.items(),
                key=lambda item: item[1][0],
            )
        )
        eliminate = tuple(
            name for name in ordered if self.eddy_block_roles[name] == "bulk"
        )
        keep = tuple(
            name for name in ordered if self.eddy_block_roles[name] != "bulk"
        )
        if not eliminate:
            raise ValueError("adjacency classification has no bulk eddy-bubble block")
        if not keep:
            raise ValueError("adjacency classification has no protected HCurl block")
        return keep, eliminate

    def eddy_bubble_mixed_galerkin_orthogonalization(
        self,
        magnetic_operator=None,
        s=None,
        *,
        eddy_operator=None,
        surface_impedance=0.0,
        coupling_scale=1.0,
        adjoint_coupling_scale=None,
    ) -> "MixedGalerkinHDivHybridVIMSystem":
        """Apply full coupled reduction using face-adjacency block roles."""

        keep, eliminate = self.adjacency_class_block_partition()
        return self.mixed_galerkin_orthogonalization(
            keep,
            eliminate,
            magnetic_operator,
            s,
            eddy_operator=eddy_operator,
            surface_impedance=surface_impedance,
            coupling_scale=coupling_scale,
            adjoint_coupling_scale=adjoint_coupling_scale,
        )

    def solve(
        self,
        magnetic_operator=None,
        s=None,
        *,
        magnetic_rhs=None,
        eddy_rhs=None,
        eddy_operator=None,
        surface_impedance=0.0,
        coupling_scale=1.0,
        adjoint_coupling_scale=None,
        solver: str = "gmres",
        tolerance: float = 1.0e-10,
        max_iterations: int | None = None,
        restart: int | None = None,
        return_operator: bool = False,
    ) -> dict[str, np.ndarray]:
        """Solve the coupled HDiv-MMM / hybrid HCurl-VIM system."""

        op = self.mixed_operator(
            magnetic_operator,
            s,
            eddy_operator=eddy_operator,
            surface_impedance=surface_impedance,
            coupling_scale=coupling_scale,
            adjoint_coupling_scale=adjoint_coupling_scale,
        )
        mh = self.n_hdiv_modes
        dtype = np.dtype(op.dtype) if _is_matrix_free_operator(op) else op.dtype
        m_rhs, e_rhs = self._resolved_rhs(magnetic_rhs, eddy_rhs, dtype)
        rhs = np.vstack([m_rhs, e_rhs]).astype(dtype, copy=False)
        if _is_matrix_free_operator(op):
            sol, solver_diagnostics = _solve_native_reduced_operator(
                op,
                rhs,
                solver=solver,
                tolerance=tolerance,
                max_iterations=max_iterations,
                restart=restart,
            )
        else:
            sol = _solve_reduced_linear(op, rhs)
            solver_diagnostics = {
                "backend": "native-dense-reduced-lu",
                "iterations": 0,
                "relative_residual_max": 0.0,
            }
        result = {
            "magnetization": sol[:mh],
            "eddy": sol[mh:],
            "solution": sol,
            "rhs": rhs,
            "solver_diagnostics": solver_diagnostics,
        }
        if return_operator:
            result["operator"] = op
        return result

    def solve_frequency(
        self,
        frequency_hz: float,
        *,
        magnetic_operator=None,
        magnetic_rhs=None,
        eddy_rhs=None,
        eddy_operator=None,
        surface_impedance=None,
        sigma: float | None = None,
        mu: float = MU0,
        coupling_scale=None,
        adjoint_coupling_scale=None,
        mixed_galerkin_keep_blocks=None,
        mixed_galerkin_eliminate_blocks=None,
    ) -> HCurlVIMHDivMMMSolution:
        """Solve one harmonic excitation and reconstruct all retained fields.

        The stored coupling block is the mutual-energy matrix
        ``K_ij = int M_i . B[J_j] dV``.  By default this physics driver applies
        the harmonic scales that follow from it and from the reciprocity
        ``int J_j . A_M[M_i] dV = K_ij``:

        * magnetic row  ``-K/mu``   (the currents' H-field, ``-<M_i, H_J>``),
        * eddy row      ``s K^H``   (the magnetization EMF, ``s <w_j, A_M[M]>``).

        The ``s`` factor is what makes a STATIC magnetization drive no eddy
        current; with unit scales a static M drives a spurious DC circulation
        limited only by ``Z(0) = R`` (measured on the sphere-alpha lane:
        coefficient norm ~2.4 vs the physical ~1e-5 at uHz).  Pass explicit
        ``coupling_scale`` / ``adjoint_coupling_scale`` for neutral-algebra
        uses (the raw containers keep their unit defaults).

        When both mixed-Galerkin block arguments are supplied, the elimination
        is performed on the complete HDiv-MMM/HCurl-VIM operator.  This updates
        the magnetic block, both coupling blocks, the RHS, and the HCurl field
        lift together.
        """

        frequency = float(frequency_hz)
        if not np.isfinite(frequency) or frequency <= 0.0:
            raise ValueError("frequency_hz must be positive")
        if mu <= 0.0:
            raise ValueError("mu must be positive")
        s = 1j * 2.0 * np.pi * frequency
        if coupling_scale is None:
            coupling_scale = -1.0 / mu
        if adjoint_coupling_scale is None:
            adjoint_coupling_scale = s
        if surface_impedance is None:
            conductivity = self.conductivity if sigma is None else sigma
            if conductivity is None and isinstance(
                self.material_model,
                SharedMeshMaterialModel,
            ) and np.isscalar(self.material_model.sigma):
                conductivity = float(self.material_model.sigma)
            if conductivity is None:
                raise ValueError(
                    "sigma or surface_impedance is required for the SIBC block"
                )
            surface_impedance = SkinImpedance(s, float(conductivity), mu)
        if not isinstance(surface_impedance, SurfaceImpedanceGram):
            surface_impedance = complex(surface_impedance)

        if magnetic_rhs is None:
            magnetic_rhs = self.magnetic_rhs
        if eddy_rhs is None:
            eddy_rhs = self.eddy_rhs
        self._resolved_rhs(
            magnetic_rhs,
            eddy_rhs,
            np.complex128,
            require_excitation=True,
        )
        use_mixed_galerkin = (
            mixed_galerkin_keep_blocks is not None
            or mixed_galerkin_eliminate_blocks is not None
        )
        if (
            mixed_galerkin_keep_blocks is None
            and mixed_galerkin_eliminate_blocks is not None
        ) or (
            mixed_galerkin_keep_blocks is not None
            and mixed_galerkin_eliminate_blocks is None
        ):
            raise ValueError(
                "mixed_galerkin_keep_blocks and "
                "mixed_galerkin_eliminate_blocks must be supplied together"
            )

        orthogonalized_rhs = None
        orthogonalized_solution = None
        mixed_galerkin_diagnostics = None
        if use_mixed_galerkin:
            reduction = self.mixed_galerkin_orthogonalization(
                mixed_galerkin_keep_blocks,
                mixed_galerkin_eliminate_blocks,
                magnetic_operator,
                s,
                eddy_operator=eddy_operator,
                surface_impedance=surface_impedance,
                coupling_scale=coupling_scale,
                adjoint_coupling_scale=adjoint_coupling_scale,
            )
            solved = reduction.solve(
                magnetic_rhs=magnetic_rhs,
                eddy_rhs=eddy_rhs,
                require_excitation=True,
                return_operator=True,
            )
            op = solved["operator"]
            rhs = solved["full_rhs"]
            orthogonalized_rhs = solved["reduced_rhs"]
            orthogonalized_solution = solved["reduced_solution"]
            mixed_galerkin_diagnostics = reduction.diagnostics()
        else:
            solved = self.solve(
                magnetic_operator,
                s,
                magnetic_rhs=magnetic_rhs,
                eddy_rhs=eddy_rhs,
                eddy_operator=eddy_operator,
                surface_impedance=surface_impedance,
                coupling_scale=coupling_scale,
                adjoint_coupling_scale=adjoint_coupling_scale,
                return_operator=True,
            )
            op = solved["operator"]
            rhs = solved["rhs"]
        coefficients = solved["solution"]
        magnetic_coefficients = solved["magnetization"]
        eddy_coefficients = solved["eddy"]

        sampled_magnetization = np.einsum(
            "mp,mnk->pnk",
            magnetic_coefficients,
            self.magnetization_basis.modes,
        )
        ordered_blocks = tuple(
            name
            for name, _ in sorted(
                self.eddy_system.blocks.items(),
                key=lambda item: item[1][0],
            )
        )
        if len(ordered_blocks) != len(self.eddy_bases):
            raise RuntimeError("eddy basis/block metadata is inconsistent")
        ordered_roles = (
            None
            if self.eddy_block_roles is None
            else tuple(self.eddy_block_roles[name] for name in ordered_blocks)
        )
        sampled_currents = []
        sample_points = []
        sample_weights = []
        for name, basis in zip(ordered_blocks, self.eddy_bases):
            block_coefficients = eddy_coefficients[self.eddy_system.block_slice(name), :]
            sampled_currents.append(
                np.einsum("mp,mnk->pnk", block_coefficients, basis.modes)
            )
            sample_points.append(basis.points)
            sample_weights.append(basis.weights)

        parent_t = None
        if self.response_basis is not None:
            first = self.eddy_system.block_slice(ordered_blocks[0])
            volume_coefficients = eddy_coefficients[first, :]
            if volume_coefficients.shape[0] != self.response_basis.rank:
                raise RuntimeError(
                    "EVRS response rank does not match the first hybrid VIM block"
                )
            parent_t = self.response_basis.vectors @ volume_coefficients
        parent_magnetization = None
        if self.hdiv_reduction is not None:
            parent_magnetization = self.hdiv_reduction.reconstruct_parent(
                magnetic_coefficients
            )

        applied = op @ coefficients
        residual = applied - rhs
        rhs_norm = float(np.linalg.norm(rhs))
        residual_relative = float(np.linalg.norm(residual) / max(rhs_norm, 1.0e-300))
        residual_backward = float(
            np.linalg.norm(residual)
            / max(
                float(np.linalg.norm(applied) + rhs_norm),
                1.0e-300,
            )
        )
        dissipative = (
            self.eddy_system.resistance
            + _surface_impedance_dissipative_term(
                self.eddy_system.surface_mass,
                surface_impedance,
            )
        )
        average_loss = 0.5 * np.real(
            np.einsum(
                "ip,ij,jp->p",
                eddy_coefficients.conj(),
                dissipative,
                eddy_coefficients,
            )
        )
        port_response = rhs.conj().T @ coefficients
        backend = solved.get("solver_diagnostics", {}).get(
            "backend",
            (
                "radia-cpp-dense"
                if _radia_cpp_kernel("_HybridVIMSolve") is not None
                else "numpy-linalg"
            ),
        )
        if use_mixed_galerkin:
            backend += "-mixed-galerkin"
        return HCurlVIMHDivMMMSolution(
            frequency_hz=frequency,
            s=s,
            surface_impedance=surface_impedance,
            magnetization_coefficients=magnetic_coefficients,
            eddy_coefficients=eddy_coefficients,
            parent_t_coefficients=parent_t,
            parent_magnetization_coefficients=parent_magnetization,
            sampled_magnetization=sampled_magnetization,
            eddy_block_names=ordered_blocks,
            sampled_eddy_currents=tuple(sampled_currents),
            eddy_sample_points=tuple(sample_points),
            eddy_sample_weights=tuple(sample_weights),
            reduced_rhs=rhs,
            reduced_solution=coefficients,
            port_response=port_response,
            average_joule_loss=np.asarray(average_loss),
            residual_relative_norm=residual_relative,
            solver_backend=backend,
            solver_diagnostics=solved.get("solver_diagnostics"),
            orthogonalized_rhs=orthogonalized_rhs,
            orthogonalized_solution=orthogonalized_solution,
            mixed_galerkin_diagnostics=mixed_galerkin_diagnostics,
            residual_backward_error=residual_backward,
            eddy_block_roles=ordered_roles,
        )

    def solve_frequency_eddy_bubbled(
        self,
        frequency_hz: float,
        **kwargs,
    ) -> HCurlVIMHDivMMMSolution:
        """Solve with adjacency-class-protected mixed Galerkin elimination."""

        forbidden = {
            "mixed_galerkin_keep_blocks",
            "mixed_galerkin_eliminate_blocks",
        } & set(kwargs)
        if forbidden:
            raise TypeError(
                "solve_frequency_eddy_bubbled selects mixed Galerkin blocks "
                "from eddy_block_roles"
            )
        keep, eliminate = self.adjacency_class_block_partition()
        return self.solve_frequency(
            frequency_hz,
            mixed_galerkin_keep_blocks=keep,
            mixed_galerkin_eliminate_blocks=eliminate,
            **kwargs,
        )

    def schur_magnetic_operator(
        self,
        magnetic_operator,
        s=None,
        *,
        eddy_operator=None,
        surface_impedance=0.0,
        coupling_scale=1.0,
        adjoint_coupling_scale=None,
    ) -> np.ndarray:
        """Return the magnetic Schur complement after eliminating eddies."""

        magnetic = _square_matrix(magnetic_operator, self.n_hdiv_modes, "magnetic_operator")
        eddy = self.eddy_impedance(
            s,
            surface_impedance=surface_impedance,
            eddy_operator=eddy_operator,
        )
        if adjoint_coupling_scale is None:
            adjoint_coupling_scale = np.conjugate(coupling_scale)
        upper = coupling_scale * self.coupling
        lower = adjoint_coupling_scale * self.coupling.conj().T
        return magnetic - upper @ np.linalg.solve(eddy, lower)

    def schur_eddy_operator(
        self,
        magnetic_operator,
        s=None,
        *,
        eddy_operator=None,
        surface_impedance=0.0,
        coupling_scale=1.0,
        adjoint_coupling_scale=None,
    ) -> np.ndarray:
        """Return the eddy Schur complement after eliminating HDiv-MMM."""

        magnetic = _square_matrix(magnetic_operator, self.n_hdiv_modes, "magnetic_operator")
        eddy = self.eddy_impedance(
            s,
            surface_impedance=surface_impedance,
            eddy_operator=eddy_operator,
        )
        if adjoint_coupling_scale is None:
            adjoint_coupling_scale = np.conjugate(coupling_scale)
        upper = coupling_scale * self.coupling
        lower = adjoint_coupling_scale * self.coupling.conj().T
        return eddy - lower @ np.linalg.solve(magnetic, upper)

    def port_admittance(
        self,
        magnetic_operator,
        s,
        eddy_rhs,
        *,
        magnetic_rhs=None,
        eddy_operator=None,
        surface_impedance=0.0,
        coupling_scale=1.0,
        adjoint_coupling_scale=None,
    ) -> np.ndarray:
        """Return mixed-system port admittance for hybrid eddy RHS columns."""

        mh = self.n_hdiv_modes
        ne = self.n_hcurl_vim_modes
        e_ports = _port_rhs_matrix(eddy_rhs, ne)
        if magnetic_rhs is None:
            m_ports = np.zeros((mh, e_ports.shape[1]), dtype=e_ports.dtype)
        else:
            m_ports = _port_rhs_matrix(magnetic_rhs, mh)
        if m_ports.shape[1] != e_ports.shape[1]:
            raise ValueError("magnetic_rhs and eddy_rhs must have the same number of columns")
        sol = self.solve(
            magnetic_operator,
            s,
            magnetic_rhs=m_ports,
            eddy_rhs=e_ports,
            eddy_operator=eddy_operator,
            surface_impedance=surface_impedance,
            coupling_scale=coupling_scale,
            adjoint_coupling_scale=adjoint_coupling_scale,
        )["solution"]
        ports = np.vstack([m_ports, e_ports])
        return ports.conj().T @ sol

    def diagnostics(self) -> dict[str, object]:
        """Return compact diagnostics for the hybrid mixed system."""

        return {
            "hdiv_modes": self.n_hdiv_modes,
            "hdiv_mmm_modes": self.n_hdiv_mmm_modes,
            "hcurl_vim_modes": self.n_hcurl_vim_modes,
            "evrs_modes": self.n_evrs_modes,
            "hybrid_blocks": self.eddy_system.diagnostics()["blocks"],
            "eddy_basis_count": len(self.eddy_bases),
            "eddy_basis_modes": [basis.n_modes for basis in self.eddy_bases],
            "eddy_block_roles": (
                None if self.eddy_block_roles is None else dict(self.eddy_block_roles)
            ),
            "coupling_rows": int(self.coupling.shape[0]),
            "coupling_cols": int(self.coupling.shape[1]),
            "coupling_frobenius_norm": float(np.linalg.norm(self.coupling)),
            "has_shared_material_model": self.material_model is not None,
            "has_shared_mesh_material_model": isinstance(self.material_model, SharedMeshMaterialModel),
            "has_magnetic_operator": self.magnetic_operator is not None,
            "has_magnetic_rhs": self.magnetic_rhs is not None,
            "has_eddy_rhs": self.eddy_rhs is not None,
            "has_response_basis": self.response_basis is not None,
            "response_basis": (
                None
                if self.response_basis is None
                else self.response_basis.diagnostics()
            ),
            "has_eddy_bubbling": self.eddy_bubbling is not None,
            "has_hdiv_reduction": self.hdiv_reduction is not None,
            "hdiv_reduction": (
                None
                if self.hdiv_reduction is None
                else self.hdiv_reduction.diagnostics()
            ),
            "conductivity": self.conductivity,
        }


@dataclass(frozen=True)
class MixedGalerkinHDivHybridVIMSystem:
    """Frequency-specific exact reduction of a coupled HDiv/HCurl system."""

    parent_system: CoupledHDivHybridVIMSystem
    full_operator: np.ndarray
    orthogonalization: MixedGalerkinOrthogonalization
    keep_eddy_blocks: tuple[str, ...]
    eliminate_eddy_blocks: tuple[str, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.parent_system, CoupledHDivHybridVIMSystem):
            raise TypeError("parent_system must be a CoupledHDivHybridVIMSystem")
        full = np.asarray(self.full_operator)
        expected = self.parent_system.n_hdiv_modes + self.parent_system.n_hcurl_vim_modes
        if full.shape != (expected, expected):
            raise ValueError(f"full_operator must have shape ({expected}, {expected})")
        if self.orthogonalization.trial_transform.shape[0] != expected:
            raise ValueError("orthogonalization does not match the full coupled operator")
        object.__setattr__(self, "full_operator", np.array(full, copy=True))
        object.__setattr__(self, "keep_eddy_blocks", tuple(self.keep_eddy_blocks))
        object.__setattr__(
            self,
            "eliminate_eddy_blocks",
            tuple(self.eliminate_eddy_blocks),
        )

    @property
    def reduced_operator(self) -> np.ndarray:
        return self.orthogonalization.reduced_operator

    @property
    def n_modes(self) -> int:
        return self.orthogonalization.rank

    @property
    def n_hdiv_modes(self) -> int:
        return self.parent_system.n_hdiv_modes

    @property
    def n_hcurl_retained_modes(self) -> int:
        return self.n_modes - self.n_hdiv_modes

    @property
    def n_hcurl_eliminated_modes(self) -> int:
        return int(self.orthogonalization.eliminate_indices.size)

    def project_rhs(
        self,
        *,
        magnetic_rhs=None,
        eddy_rhs=None,
        require_excitation: bool = False,
    ) -> dict[str, np.ndarray]:
        """Project a coupled RHS and build its eliminated-block affine lift."""

        m_rhs, e_rhs = self.parent_system._resolved_rhs(
            magnetic_rhs,
            eddy_rhs,
            self.full_operator.dtype,
            require_excitation=require_excitation,
        )
        full_rhs = np.vstack((m_rhs, e_rhs)).astype(
            self.full_operator.dtype,
            copy=False,
        )
        orth = self.orthogonalization
        reduced_rhs = orth.test_transform.T @ full_rhs
        particular = np.zeros_like(full_rhs)
        eliminated = orth.eliminate_indices
        z_ee = self.full_operator[np.ix_(eliminated, eliminated)]
        particular[eliminated, :] = np.linalg.solve(
            z_ee,
            full_rhs[eliminated, :],
        )
        return {
            "full_rhs": full_rhs,
            "reduced_rhs": reduced_rhs,
            "particular_solution": particular,
        }

    def solve(
        self,
        *,
        magnetic_rhs=None,
        eddy_rhs=None,
        require_excitation: bool = False,
        return_operator: bool = False,
    ) -> dict[str, np.ndarray | float]:
        """Solve the Schur system and reconstruct the exact full coefficients."""

        projected = self.project_rhs(
            magnetic_rhs=magnetic_rhs,
            eddy_rhs=eddy_rhs,
            require_excitation=require_excitation,
        )
        reduced_solution = _solve_reduced_linear(
            self.reduced_operator,
            projected["reduced_rhs"],
        )
        full_solution = (
            projected["particular_solution"]
            + self.orthogonalization.trial_transform @ reduced_solution
        )
        residual = self.full_operator @ full_solution - projected["full_rhs"]
        rhs_norm = max(float(np.linalg.norm(projected["full_rhs"])), 1.0e-300)
        projected_residual = (
            self.orthogonalization.test_transform.T @ residual
        )
        reduced_rhs_norm = max(
            float(np.linalg.norm(projected["reduced_rhs"])),
            1.0e-300,
        )
        mh = self.n_hdiv_modes
        result = {
            **projected,
            "reduced_solution": reduced_solution,
            "full_solution": full_solution,
            "solution": full_solution,
            "magnetization": full_solution[:mh, :],
            "eddy": full_solution[mh:, :],
            "full_residual_relative_norm": float(np.linalg.norm(residual) / rhs_norm),
            "projected_residual_relative_norm": float(
                np.linalg.norm(projected_residual) / reduced_rhs_norm
            ),
        }
        if return_operator:
            result["operator"] = self.full_operator
            result["reduced_operator"] = self.reduced_operator
        return result

    def diagnostics(self) -> dict[str, object]:
        info = self.orthogonalization.diagnostics()
        hdiv_reduction = self.parent_system.hdiv_reduction
        demag = None if hdiv_reduction is None else hdiv_reduction.diagnostics()
        return {
            **info,
            "kind": "mixed-galerkin-hdiv-hcurl",
            "full_coupled_schur": True,
            "hdiv_retained_modes": self.n_hdiv_modes,
            "hcurl_retained_modes": self.n_hcurl_retained_modes,
            "hcurl_eliminated_modes": self.n_hcurl_eliminated_modes,
            "keep_eddy_blocks": list(self.keep_eddy_blocks),
            "eliminate_eddy_blocks": list(self.eliminate_eddy_blocks),
            "hdiv_demag_backend": (
                None if demag is None else demag["demag_backend"]
            ),
            "hdiv_demag_hmatrix_backend": (
                None if demag is None else demag.get("demag_hmatrix_backend")
            ),
        }


def CoupleHDivMagnetizationToEVRS(
    magnetization_basis: SampledMagnetizationBasis,
    eddy_basis: SampledCurrentBasis,
    *,
    eddy_system: HybridVIMSystem | None = None,
    material_model: object | None = None,
    mu: float = MU0,
    kernel_epsilon: float = 0.0,
) -> CoupledHDivEVRSSystem:
    """Build the named HDiv-MMM / EVRS eddy-current coupling container."""

    coupling = MagnetizationCurrentCoupling(
        magnetization_basis,
        eddy_basis,
        mu=mu,
        kernel_epsilon=kernel_epsilon,
    )
    return CoupledHDivEVRSSystem(
        magnetization_basis=magnetization_basis,
        eddy_basis=eddy_basis,
        coupling=coupling,
        eddy_system=eddy_system,
        material_model=material_model,
    )


HCurlVIMHDivMMMSystem = CoupledHDivEVRSSystem


def CoupleHCurlVIMWithHDivMMM(
    magnetization_basis: SampledMagnetizationBasis,
    eddy_basis: SampledCurrentBasis,
    *,
    eddy_system: HybridVIMSystem | None = None,
    material_model: object | None = None,
    mu: float = MU0,
    kernel_epsilon: float = 0.0,
) -> HCurlVIMHDivMMMSystem:
    """Build the production-named HCurl-VIM / HDiv-MMM mixed system.

    Use this name when the intent is the coupled method: a shared mesh/material
    registry feeds an HDiv-MMM magnetic branch and a response-compressed
    HCurl-VIM eddy branch.  The returned object is the same concrete container
    as :func:`CoupleHDivMagnetizationToEVRS`.
    """

    return CoupleHDivMagnetizationToEVRS(
        magnetization_basis,
        eddy_basis,
        eddy_system=eddy_system,
        material_model=material_model,
        mu=mu,
        kernel_epsilon=kernel_epsilon,
    )


def CoupleHybridVIMWithHDivMMM(
    magnetization_basis: SampledMagnetizationBasis,
    eddy_system: HybridVIMSystem,
    eddy_bases,
    *,
    material_model: object | None = None,
    magnetic_operator=None,
    magnetic_rhs=None,
    eddy_rhs=None,
    response_basis: ResponseBasis | None = None,
    eddy_bubbling: EddyBubbleDecomposition | None = None,
    hdiv_reduction: HDivMMMReducedModel | None = None,
    conductivity: float | None = None,
    eddy_block_roles: dict[str, str] | None = None,
    mu: float = MU0,
    kernel_epsilon: float = 0.0,
) -> CoupledHDivHybridVIMSystem:
    """Couple a full hybrid HCurl-VIM eddy system to HDiv-MMM."""

    bases = tuple(eddy_bases)
    if not bases:
        raise ValueError("eddy_bases must not be empty")
    blocks = [
        MagnetizationCurrentCoupling(
            magnetization_basis,
            basis,
            mu=mu,
            kernel_epsilon=kernel_epsilon,
        )
        for basis in bases
    ]
    coupling = np.hstack(blocks)
    return CoupledHDivHybridVIMSystem(
        magnetization_basis=magnetization_basis,
        eddy_system=eddy_system,
        eddy_bases=bases,
        coupling=coupling,
        material_model=material_model,
        magnetic_operator=magnetic_operator,
        magnetic_rhs=magnetic_rhs,
        eddy_rhs=eddy_rhs,
        response_basis=response_basis,
        eddy_bubbling=eddy_bubbling,
        hdiv_reduction=hdiv_reduction,
        conductivity=conductivity,
        eddy_block_roles=eddy_block_roles,
    )


def CoupleEddyBubbleHCurlBasisWithHDivMMM(
    magnetization_basis: SampledMagnetizationBasis,
    eddy_basis: EddyBubbleHCurlBasis,
    *,
    eddy_system: HybridVIMSystem | None = None,
    sigma: float | None = None,
    material_model: object | None = None,
    mu: float = MU0,
    kernel_epsilon: float = 0.0,
) -> HCurlVIMHDivMMMSystem:
    """Couple a production eddy-bubbled HCurl basis to HDiv-MMM."""

    if not isinstance(eddy_basis, EddyBubbleHCurlBasis):
        raise TypeError("eddy_basis must be an EddyBubbleHCurlBasis")
    return eddy_basis.couple_hdiv_mmm(
        magnetization_basis,
        eddy_system=eddy_system,
        sigma=sigma,
        material_model=material_model,
        mu=mu,
        kernel_epsilon=kernel_epsilon,
    )


def _same_quadrature(left: SampledCurrentBasis, right: SampledCurrentBasis) -> bool:
    return (
        left.points.shape == right.points.shape
        and np.allclose(left.points, right.points)
        and np.allclose(left.weights, right.weights)
    )


def _mass_cross(left: SampledCurrentBasis, right: SampledCurrentBasis) -> np.ndarray:
    if not _same_quadrature(left, right):
        return np.zeros((left.n_modes, right.n_modes), dtype=np.result_type(left.modes, right.modes))
    return np.einsum(
        "aik,bik,i->ab", left.modes.conj(), right.modes, left.weights
    )


def _sampled_current_basis_diagnostics(basis: SampledCurrentBasis) -> dict[str, object]:
    mass = basis.mass_matrix()
    min_eigenvalue = 0.0
    if mass.size:
        hermitian_mass = 0.5 * (mass + mass.conj().T)
        min_eigenvalue = float(np.min(np.linalg.eigvalsh(hermitian_mass).real))
    return {
        "kind": basis.kind,
        "modes": int(basis.n_modes),
        "samples": int(basis.n_samples),
        "mass_trace": float(np.trace(mass).real) if mass.size else 0.0,
        "min_mass_eigenvalue": min_eigenvalue,
    }


def AssembleHybridVIM(
    *bases: SampledCurrentBasis,
    sigma: float | None = None,
    mu: float = MU0,
    kernel_epsilon: float | None = None,
    interaction=None,
) -> HybridVIMSystem:
    """Assemble reduced VIM matrices for volume-T and surface-Omega bases.

    ``sigma`` is applied to all volume-current bases.  Surface-current bases
    contribute to ``surface_mass``; multiply that block by ``SkinImpedance`` in
    :meth:`HybridVIMSystem.impedance`.

    ``interaction`` is an optional backend hook.  The pairwise contract calls
    ``interaction(left_basis, right_basis)`` for upper-triangular basis pairs.
    A matrix-free backend may additionally implement ``build_operator(bases)``
    and declare ``operator_scope`` as ``"full"`` or ``"cross"``.  The latter
    supplies all reciprocal off-diagonal blocks in one symmetric operator while
    diagonal pairs continue through the pairwise hook.  This is the production
    connection point for sampled HACApK and ``ngsolve.bem`` projected operators.
    """

    if not bases:
        raise ValueError("at least one basis is required")
    for basis in bases:
        if not isinstance(basis, SampledCurrentBasis):
            raise TypeError("all arguments must be SampledCurrentBasis objects")
    if sigma is not None and sigma <= 0.0:
        raise ValueError("sigma must be positive")
    if mu <= 0.0:
        raise ValueError("mu must be positive")
    if interaction is not None and not callable(interaction):
        raise TypeError("interaction must be callable")
    if interaction is None:
        if kernel_epsilon is None:
            kernel_epsilon = _default_kernel_epsilon(tuple(bases))
        if kernel_epsilon <= 0.0:
            raise ValueError("kernel_epsilon must be positive")

    total = sum(basis.n_modes for basis in bases)
    rmat = np.zeros((total, total), dtype=complex)
    lmat = np.zeros((total, total), dtype=complex)
    smat = np.zeros((total, total), dtype=complex)
    operator_blocks = []
    names: list[str] = []
    blocks: dict[str, tuple[int, int]] = {}

    offsets = _basis_offsets(tuple(bases))
    assembled_interaction_operator = None
    assembled_operator_scope = None
    if interaction is not None:
        operator_builder = getattr(interaction, "build_operator", None)
        if callable(operator_builder):
            assembled_interaction_operator = operator_builder(tuple(bases))
            assembled_operator_scope = str(
                getattr(interaction, "operator_scope", "full")
            ).lower()
            if assembled_operator_scope not in {"full", "cross"}:
                raise ValueError("interaction.operator_scope must be 'full' or 'cross'")
            if assembled_interaction_operator is not None:
                if not _is_matrix_free_operator(assembled_interaction_operator):
                    raise TypeError("interaction.build_operator must return a matrix-free operator")
                if tuple(assembled_interaction_operator.shape) != (total, total):
                    raise ValueError(
                        "interaction.build_operator returned shape "
                        f"{assembled_interaction_operator.shape}, expected {(total, total)}"
                    )
                if not getattr(assembled_interaction_operator, "is_hermitian", False):
                    raise ValueError("assembled interaction operator must declare is_hermitian=True")
    for i, basis in enumerate(bases):
        start, stop = offsets[i]
        prefix = basis.kind if basis.kind not in blocks else f"{basis.kind}{i}"
        blocks[prefix] = (start, stop)
        names.extend(basis.names)

    for ia, left in enumerate(bases):
        a0, a1 = offsets[ia]
        for ib in range(ia, len(bases)):
            right = bases[ib]
            b0, b1 = offsets[ib]
            assembled_pair = assembled_interaction_operator is not None and (
                assembled_operator_scope == "full"
                or (assembled_operator_scope == "cross" and ia != ib)
            )
            if not assembled_pair:
                if interaction is None:
                    block = _interaction_block(
                        left, right, mu=mu, kernel_epsilon=kernel_epsilon
                    )
                else:
                    block = _validate_interaction_block(
                        interaction(left, right), left, right
                    )
                if _is_matrix_free_operator(block):
                    if ia != ib:
                        raise NotImplementedError(
                            "rectangular pair operators must be supplied through "
                            "interaction.build_operator(bases)"
                        )
                    operator_blocks.append((a0, a1, block))
                else:
                    lmat[a0:a1, b0:b1] = block
                    if ia != ib:
                        lmat[b0:b1, a0:a1] = block.conj().T
            if left.kind == "volume" and right.kind == "volume":
                if sigma is None:
                    raise ValueError("sigma is required when volume bases are used")
                rmat[a0:a1, b0:b1] = _mass_cross(left, right) / sigma
                if ia != ib:
                    rmat[b0:b1, a0:a1] = rmat[a0:a1, b0:b1].conj().T
            if left.kind == "surface" and right.kind == "surface":
                smat[a0:a1, b0:b1] = _mass_cross(left, right)
                if ia != ib:
                    smat[b0:b1, a0:a1] = smat[a0:a1, b0:b1].conj().T

    rmat = 0.5 * (rmat + rmat.conj().T)
    lmat = 0.5 * (lmat + lmat.conj().T)
    smat = 0.5 * (smat + smat.conj().T)
    if assembled_interaction_operator is not None:
        operator_blocks.append((0, total, assembled_interaction_operator))
    inductance = (
        ReducedBlockHMatrixOperator(lmat, operator_blocks)
        if operator_blocks
        else lmat
    )
    return HybridVIMSystem(
        resistance=rmat,
        inductance=inductance,
        surface_mass=smat,
        basis_names=tuple(names),
        blocks=blocks,
        interaction_backend=_interaction_backend_name(interaction),
    )


def AssembleSurfaceImpedanceGram(
    system: HybridVIMSystem,
    surface_basis: SampledCurrentBasis,
    surface_impedance,
    *,
    block: str = "surface",
    label: str = "local-sibc",
) -> SurfaceImpedanceGram:
    """Assemble ``int_Gamma Zs(x) jt.jv dS`` on one surface block.

    ``surface_impedance`` may be a scalar, one complex value per surface
    quadrature sample, or a callable evaluated at ``surface_basis.points``.
    The returned full-system Gram can be passed anywhere a scalar
    ``surface_impedance`` is accepted, including Schur complements and mixed
    Galerkin orthogonalization.
    """

    if not isinstance(system, HybridVIMSystem):
        raise TypeError("system must be a HybridVIMSystem")
    if not isinstance(surface_basis, SampledCurrentBasis):
        raise TypeError("surface_basis must be a SampledCurrentBasis")
    if surface_basis.kind != "surface":
        raise ValueError("surface_basis must have kind='surface'")
    values = (
        surface_impedance(surface_basis.points)
        if callable(surface_impedance)
        else surface_impedance
    )
    if np.isscalar(values):
        values = np.full(surface_basis.n_samples, complex(values), dtype=complex)
    else:
        values = np.asarray(values, dtype=complex).reshape(-1)
    if values.shape != (surface_basis.n_samples,):
        raise ValueError(
            "surface_impedance must provide one value per surface sample; "
            f"expected {surface_basis.n_samples}, got {values.size}"
        )
    if not np.all(np.isfinite(values)):
        raise ValueError("surface_impedance contains non-finite sample values")

    weighted = np.einsum(
        "aik,bik,i,i->ab",
        surface_basis.modes.conj(),
        surface_basis.modes,
        surface_basis.weights,
        values,
    )
    target = system.block_slice(block)
    expected = target.stop - target.start
    if weighted.shape != (expected, expected):
        raise ValueError(
            f"surface basis has {surface_basis.n_modes} modes but block {block!r} "
            f"has size {expected}"
        )
    reference_mass = np.zeros_like(system.surface_mass, dtype=complex)
    reference_mass[target, target] = surface_basis.mass_matrix()
    mass_scale = max(
        float(np.linalg.norm(system.surface_mass)),
        np.finfo(float).tiny,
    )
    if np.linalg.norm(reference_mass - system.surface_mass) / mass_scale > 1.0e-10:
        raise ValueError(
            "surface_basis quadrature/modes do not match the requested hybrid "
            "system surface block"
        )
    matrix = np.zeros_like(system.resistance, dtype=complex)
    matrix[target, target] = weighted
    return SurfaceImpedanceGram(
        matrix=matrix,
        sample_values=values,
        surface_mass_reference=system.surface_mass,
        label=label,
    )


_LOCAL_ESIM_LUT_SCHEMA_VERSION = 1
_LOCAL_ESIM_LUT_INTERPOLATION = "bilinear-logf-logh-logabs-phase-v1"


def _canonical_json_sha256(value: object) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True)
class LocalESIMSurfaceLUT:
    """Persistent local-ESIM surface-impedance lookup table.

    The table stores one material/cell-model signature on a tensor grid in
    positive frequency and positive tangential-field magnitude.  Evaluation
    is bilinear in ``log(f)`` and ``log(|H_t|)``.  Extrapolation is forbidden:
    a production solve must either provide a table covering its operating
    range or explicitly run the cell problem to build a wider table.  The
    interpolated variables are log impedance magnitude and continuous phase,
    which exactly represents the square-root frequency scaling of the linear
    skin impedance.
    """

    frequency_nodes_hz: np.ndarray
    h_nodes_A_per_m: np.ndarray
    impedance_ohm: np.ndarray
    cell_iterations: np.ndarray
    model_signature: dict[str, object]
    schema_version: int = _LOCAL_ESIM_LUT_SCHEMA_VERSION
    _log_frequency_nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _log_h_nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _log_impedance_magnitude: np.ndarray = field(
        init=False,
        repr=False,
        compare=False,
    )
    _unwrapped_impedance_phase: np.ndarray = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        frequencies = np.asarray(self.frequency_nodes_hz, dtype=float).reshape(-1)
        fields = np.asarray(self.h_nodes_A_per_m, dtype=float).reshape(-1)
        impedance = np.asarray(self.impedance_ohm, dtype=complex)
        iterations = np.asarray(self.cell_iterations, dtype=int)
        if frequencies.size < 1 or fields.size < 2:
            raise ValueError(
                "local ESIM LUT requires at least one frequency and two field nodes"
            )
        if not np.all(np.isfinite(frequencies)) or np.any(frequencies <= 0.0):
            raise ValueError("local ESIM LUT frequencies must be finite and positive")
        if not np.all(np.isfinite(fields)) or np.any(fields <= 0.0):
            raise ValueError("local ESIM LUT field nodes must be finite and positive")
        if np.any(np.diff(frequencies) <= 0.0):
            raise ValueError("local ESIM LUT frequencies must be strictly increasing")
        if np.any(np.diff(fields) <= 0.0):
            raise ValueError("local ESIM LUT field nodes must be strictly increasing")
        expected_shape = (frequencies.size, fields.size)
        if impedance.shape != expected_shape:
            raise ValueError(
                "local ESIM LUT impedance must have shape "
                f"{expected_shape}, got {impedance.shape}"
            )
        if iterations.shape != expected_shape:
            raise ValueError(
                "local ESIM LUT cell_iterations must have shape "
                f"{expected_shape}, got {iterations.shape}"
            )
        if not np.all(np.isfinite(impedance)):
            raise ValueError("local ESIM LUT contains non-finite impedance")
        if np.any(np.abs(impedance) <= np.finfo(float).tiny):
            raise ValueError("local ESIM LUT impedance magnitude must be positive")
        if np.min(impedance.real) < -1.0e-12:
            raise ValueError("local ESIM LUT contains non-passive impedance")
        if np.any(iterations < 1):
            raise ValueError("local ESIM LUT cell iterations must be positive")
        signature = json.loads(json.dumps(self.model_signature, sort_keys=True))
        if not isinstance(signature, dict) or not signature:
            raise ValueError("local ESIM LUT model_signature must be a non-empty dict")
        schema_version = int(self.schema_version)
        if schema_version != _LOCAL_ESIM_LUT_SCHEMA_VERSION:
            raise ValueError(
                "unsupported local ESIM LUT schema version: "
                f"{schema_version}"
            )
        object.__setattr__(self, "frequency_nodes_hz", np.array(frequencies, copy=True))
        object.__setattr__(self, "h_nodes_A_per_m", np.array(fields, copy=True))
        object.__setattr__(self, "impedance_ohm", np.array(impedance, copy=True))
        object.__setattr__(self, "cell_iterations", np.array(iterations, copy=True))
        object.__setattr__(self, "model_signature", signature)
        object.__setattr__(self, "schema_version", schema_version)
        object.__setattr__(self, "_log_frequency_nodes", np.log(frequencies))
        object.__setattr__(self, "_log_h_nodes", np.log(fields))
        object.__setattr__(
            self,
            "_log_impedance_magnitude",
            np.log(np.abs(impedance)),
        )
        object.__setattr__(
            self,
            "_unwrapped_impedance_phase",
            np.unwrap(np.unwrap(np.angle(impedance), axis=1), axis=0),
        )

    @property
    def model_key(self) -> str:
        """Return the stable SHA-256 key for the material/cell configuration."""

        return _canonical_json_sha256(self.model_signature)

    def validate_model(self, model: "LocalESIMSurfaceModel") -> None:
        """Reject use with a different BH curve, conductivity, or cell setup."""

        expected = model.lut_signature()
        if expected != self.model_signature:
            raise ValueError(
                "local ESIM LUT model signature does not match the requested "
                "BH curve, conductivity, or cell settings"
            )

    def evaluate(self, tangential_field_magnitude, frequency_hz: float) -> np.ndarray:
        """Interpolate impedance without extrapolating beyond the LUT domain."""

        fields = np.asarray(tangential_field_magnitude, dtype=float).reshape(-1)
        if fields.size == 0:
            raise ValueError("tangential field samples must not be empty")
        if not np.all(np.isfinite(fields)) or np.any(fields <= 0.0):
            raise ValueError("LUT tangential field samples must be finite and positive")
        frequency = float(frequency_hz)
        if not np.isfinite(frequency) or frequency <= 0.0:
            raise ValueError("frequency_hz must be positive")
        h_min = float(self.h_nodes_A_per_m[0])
        h_max = float(self.h_nodes_A_per_m[-1])
        f_min = float(self.frequency_nodes_hz[0])
        f_max = float(self.frequency_nodes_hz[-1])
        range_rtol = 32.0 * np.finfo(float).eps
        if np.min(fields) < h_min * (1.0 - range_rtol) or np.max(fields) > h_max * (
            1.0 + range_rtol
        ):
            raise ValueError(
                "tangential field is outside the local ESIM LUT range "
                f"[{h_min:.6g}, {h_max:.6g}] A/m"
            )
        if frequency < f_min * (1.0 - range_rtol) or frequency > f_max * (
            1.0 + range_rtol
        ):
            raise ValueError(
                "frequency is outside the local ESIM LUT range "
                f"[{f_min:.6g}, {f_max:.6g}] Hz"
            )
        clipped_fields = np.clip(fields, h_min, h_max)
        clipped_frequency = float(np.clip(frequency, f_min, f_max))
        log_fields = np.log(clipped_fields)
        log_magnitude_by_frequency = np.empty(
            (self.frequency_nodes_hz.size, fields.size),
            dtype=float,
        )
        phase_by_frequency = np.empty_like(log_magnitude_by_frequency)
        for index in range(self.frequency_nodes_hz.size):
            log_magnitude_by_frequency[index] = np.interp(
                log_fields,
                self._log_h_nodes,
                self._log_impedance_magnitude[index],
            )
            phase_by_frequency[index] = np.interp(
                log_fields,
                self._log_h_nodes,
                self._unwrapped_impedance_phase[index],
            )
        if self.frequency_nodes_hz.size == 1:
            interpolated_log_magnitude = log_magnitude_by_frequency[0]
            interpolated_phase = phase_by_frequency[0]
        else:
            log_frequency = np.log(clipped_frequency)
            upper = int(np.searchsorted(self._log_frequency_nodes, log_frequency))
            upper = min(max(upper, 1), self.frequency_nodes_hz.size - 1)
            lower = upper - 1
            denominator = (
                self._log_frequency_nodes[upper]
                - self._log_frequency_nodes[lower]
            )
            fraction = (log_frequency - self._log_frequency_nodes[lower]) / denominator
            interpolated_log_magnitude = (
                (1.0 - fraction) * log_magnitude_by_frequency[lower]
                + fraction * log_magnitude_by_frequency[upper]
            )
            interpolated_phase = (
                (1.0 - fraction) * phase_by_frequency[lower]
                + fraction * phase_by_frequency[upper]
            )
        values = np.exp(interpolated_log_magnitude + 1j * interpolated_phase)
        if not np.all(np.isfinite(values)):
            raise RuntimeError("local ESIM LUT interpolation produced non-finite impedance")
        if np.min(values.real) < -1.0e-12:
            raise RuntimeError("local ESIM LUT interpolation produced non-passive impedance")
        return values

    def diagnostics(self) -> dict[str, object]:
        """Return JSON-ready provenance and coverage information."""

        return {
            "schema_version": int(self.schema_version),
            "interpolation": _LOCAL_ESIM_LUT_INTERPOLATION,
            "model_key_sha256": self.model_key,
            "frequency_count": int(self.frequency_nodes_hz.size),
            "field_count": int(self.h_nodes_A_per_m.size),
            "frequency_min_Hz": float(self.frequency_nodes_hz[0]),
            "frequency_max_Hz": float(self.frequency_nodes_hz[-1]),
            "h_min_A_per_m": float(self.h_nodes_A_per_m[0]),
            "h_max_A_per_m": float(self.h_nodes_A_per_m[-1]),
            "cell_iterations_max": int(np.max(self.cell_iterations)),
            "passive": bool(np.min(self.impedance_ohm.real) >= -1.0e-12),
        }

    def save(self, path) -> Path:
        """Atomically save this LUT as a pickle-free compressed NPZ file."""

        destination = Path(path)
        if destination.suffix.lower() != ".npz":
            raise ValueError("local ESIM LUT path must end in .npz")
        destination.parent.mkdir(parents=True, exist_ok=True)
        metadata = json.dumps(
            {
                "schema_version": int(self.schema_version),
                "interpolation": _LOCAL_ESIM_LUT_INTERPOLATION,
                "model_signature": self.model_signature,
                "model_key_sha256": self.model_key,
            },
            sort_keys=True,
            separators=(",", ":"),
        )
        temporary = destination.with_name(destination.name + ".tmp.npz")
        try:
            np.savez_compressed(
                temporary,
                frequency_nodes_hz=self.frequency_nodes_hz,
                h_nodes_A_per_m=self.h_nodes_A_per_m,
                impedance_real_ohm=self.impedance_ohm.real,
                impedance_imag_ohm=self.impedance_ohm.imag,
                cell_iterations=self.cell_iterations,
                metadata=np.asarray(metadata),
            )
            os.replace(temporary, destination)
        finally:
            if temporary.exists():
                temporary.unlink()
        return destination

    @classmethod
    def load(cls, path) -> "LocalESIMSurfaceLUT":
        """Load and validate a LUT written by :meth:`save`."""

        source = Path(path)
        with np.load(source, allow_pickle=False) as archive:
            required = {
                "frequency_nodes_hz",
                "h_nodes_A_per_m",
                "impedance_real_ohm",
                "impedance_imag_ohm",
                "cell_iterations",
                "metadata",
            }
            missing = required.difference(archive.files)
            if missing:
                raise ValueError(
                    "local ESIM LUT archive is missing: " + ", ".join(sorted(missing))
                )
            metadata = json.loads(str(archive["metadata"].item()))
            if not isinstance(metadata, dict):
                raise ValueError("local ESIM LUT metadata must be an object")
            if metadata.get("interpolation") != _LOCAL_ESIM_LUT_INTERPOLATION:
                raise ValueError("unsupported local ESIM LUT interpolation scheme")
            lut = cls(
                frequency_nodes_hz=archive["frequency_nodes_hz"],
                h_nodes_A_per_m=archive["h_nodes_A_per_m"],
                impedance_ohm=(
                    archive["impedance_real_ohm"]
                    + 1j * archive["impedance_imag_ohm"]
                ),
                cell_iterations=archive["cell_iterations"],
                model_signature=metadata.get("model_signature", {}),
                schema_version=metadata.get("schema_version", -1),
            )
        if metadata.get("model_key_sha256") != lut.model_key:
            raise ValueError("local ESIM LUT model signature hash is inconsistent")
        return lut


@dataclass(frozen=True)
class _LocalESIMImpedanceTable:
    frequency_hz: float
    h_nodes: np.ndarray
    z_nodes: np.ndarray
    cell_iterations_max: int
    skin_depth_initial_m: float

    def __post_init__(self) -> None:
        h_nodes = np.asarray(self.h_nodes, dtype=float).reshape(-1)
        z_nodes = np.asarray(self.z_nodes, dtype=complex).reshape(-1)
        if h_nodes.size < 2 or z_nodes.shape != h_nodes.shape:
            raise ValueError("local ESIM table requires matching node arrays")
        if not np.all(np.isfinite(h_nodes)) or not np.all(np.isfinite(z_nodes)):
            raise ValueError("local ESIM table contains non-finite values")
        if np.any(h_nodes <= 0.0) or np.any(np.diff(h_nodes) <= 0.0):
            raise ValueError("local ESIM table field nodes must be positive and increasing")
        object.__setattr__(self, "h_nodes", np.array(h_nodes, copy=True))
        object.__setattr__(self, "z_nodes", np.array(z_nodes, copy=True))


@dataclass(frozen=True)
class LocalESIMSurfaceModel:
    """Nonlinear one-dimensional ESIM cell model for one SIBC surface block."""

    bh_curve: np.ndarray
    sigma: float
    bins: int = 12
    n_nodes: int = 100
    num_skin_depths: float = 10.0
    cell_tolerance: float = 1.0e-5
    cell_max_iterations: int = 80
    cell_relaxation: float = 0.5
    h_floor: float = 1.0e-6
    table_range_margin: float = 2.0
    lut: LocalESIMSurfaceLUT | None = None
    _table_cache: dict[float, _LocalESIMImpedanceTable] = field(
        default_factory=dict,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        curve = np.asarray(self.bh_curve, dtype=float)
        if curve.ndim != 2 or curve.shape[1] != 2 or curve.shape[0] < 2:
            raise ValueError("bh_curve must have shape (n, 2) with n >= 2")
        if not np.all(np.isfinite(curve)):
            raise ValueError("bh_curve contains non-finite values")
        if np.any(curve[:, 0] < 0.0) or np.any(np.diff(curve[:, 0]) <= 0.0):
            raise ValueError("bh_curve H values must be non-negative and increasing")
        if np.any(curve[:, 1] < 0.0):
            raise ValueError("bh_curve B values must be non-negative")
        sigma = float(self.sigma)
        if not np.isfinite(sigma) or sigma <= 0.0:
            raise ValueError("sigma must be positive")
        bins = int(self.bins)
        n_nodes = int(self.n_nodes)
        max_iterations = int(self.cell_max_iterations)
        if bins < 2:
            raise ValueError("bins must be >= 2")
        if n_nodes < 3:
            raise ValueError("n_nodes must be >= 3")
        if max_iterations < 1:
            raise ValueError("cell_max_iterations must be >= 1")
        for name in (
            "num_skin_depths",
            "cell_tolerance",
            "h_floor",
            "table_range_margin",
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive")
            object.__setattr__(self, name, value)
        if self.table_range_margin <= 1.0:
            raise ValueError("table_range_margin must be greater than 1")
        relaxation = float(self.cell_relaxation)
        if not np.isfinite(relaxation) or not 0.0 < relaxation <= 1.0:
            raise ValueError("cell_relaxation must be in (0, 1]")
        object.__setattr__(self, "bh_curve", np.array(curve, copy=True))
        object.__setattr__(self, "sigma", sigma)
        object.__setattr__(self, "bins", bins)
        object.__setattr__(self, "n_nodes", n_nodes)
        object.__setattr__(self, "cell_max_iterations", max_iterations)
        object.__setattr__(self, "cell_relaxation", relaxation)
        if self.lut is not None:
            if not isinstance(self.lut, LocalESIMSurfaceLUT):
                raise TypeError("lut must be a LocalESIMSurfaceLUT")
            self.lut.validate_model(self)

    def lut_signature(self) -> dict[str, object]:
        """Return the stable material and cell-discretization signature."""

        curve = np.ascontiguousarray(self.bh_curve, dtype="<f8")
        return {
            "schema_version": _LOCAL_ESIM_LUT_SCHEMA_VERSION,
            "cell_model": "radia.esim_cell_problem.ESIMCellProblemSolver",
            "bh_curve_sha256": hashlib.sha256(curve.tobytes()).hexdigest(),
            "bh_curve_points": int(curve.shape[0]),
            "sigma_S_per_m": float(self.sigma),
            "n_nodes": int(self.n_nodes),
            "num_skin_depths": float(self.num_skin_depths),
            "cell_tolerance": float(self.cell_tolerance),
            "cell_max_iterations": int(self.cell_max_iterations),
            "cell_relaxation": float(self.cell_relaxation),
        }

    def with_lut(self, lut: LocalESIMSurfaceLUT) -> "LocalESIMSurfaceModel":
        """Return an equivalent model that evaluates a precomputed LUT."""

        if not isinstance(lut, LocalESIMSurfaceLUT):
            raise TypeError("lut must be a LocalESIMSurfaceLUT")
        lut.validate_model(self)
        return replace(self, lut=lut)

    def _cell_solver(self, frequency_hz: float):
        frequency = float(frequency_hz)
        if not np.isfinite(frequency) or frequency <= 0.0:
            raise ValueError("frequency_hz must be positive")
        from radia.esim_cell_problem import ESIMCellProblemSolver

        return ESIMCellProblemSolver(
            bh_curve=self.bh_curve,
            sigma=self.sigma,
            frequency=frequency,
            num_skin_depths=self.num_skin_depths,
            n_nodes=self.n_nodes,
        )

    def linear_impedance(self, frequency_hz: float) -> complex:
        """Return the low-field linear SIBC used to start the Karl iteration."""

        if self.lut is not None:
            self.lut.validate_model(self)
            return complex(
                self.lut.evaluate([float(self.lut.h_nodes_A_per_m[0])], frequency_hz)[0]
            )
        return complex(self._cell_solver(frequency_hz).get_linear_sibc())

    def _impedance_samples_cached(
        self,
        tangential_field_magnitude,
        frequency_hz: float,
        table: _LocalESIMImpedanceTable | None,
    ) -> tuple[
        np.ndarray,
        dict[str, object],
        _LocalESIMImpedanceTable | None,
    ]:
        """Interpolate a persistent table, extending it only when necessary."""

        fields = np.asarray(tangential_field_magnitude, dtype=float).reshape(-1)
        if fields.size == 0:
            raise ValueError("tangential field samples must not be empty")
        if not np.all(np.isfinite(fields)) or np.any(fields < 0.0):
            raise ValueError("tangential field magnitudes must be finite and non-negative")
        frequency = float(frequency_hz)
        if not np.isfinite(frequency) or frequency <= 0.0:
            raise ValueError("frequency_hz must be positive")
        effective = np.maximum(fields, self.h_floor)
        if self.lut is not None:
            self.lut.validate_model(self)
            values = self.lut.evaluate(effective, frequency)
            lut_info = self.lut.diagnostics()
            return values, {
                "cell_model": "precomputed-2d-local-esim-lut",
                "frequency_Hz": frequency,
                "bins": int(self.lut.h_nodes_A_per_m.size),
                "cell_iterations_max": int(lut_info["cell_iterations_max"]),
                "cell_converged": True,
                "cell_table_rebuilt": False,
                "cell_solve_count": 0,
                "lut_model_key_sha256": self.lut.model_key,
                "lut_frequency_count": int(self.lut.frequency_nodes_hz.size),
                "h_min_A_per_m": float(np.min(fields)),
                "h_max_A_per_m": float(np.max(fields)),
                "table_h_min_A_per_m": float(self.lut.h_nodes_A_per_m[0]),
                "table_h_max_A_per_m": float(self.lut.h_nodes_A_per_m[-1]),
                "z_real_min_ohm": float(np.min(values.real)),
                "z_real_max_ohm": float(np.max(values.real)),
                "z_imag_min_ohm": float(np.min(values.imag)),
                "z_imag_max_ohm": float(np.max(values.imag)),
            }, table
        required_nodes = max(2, min(self.bins, fields.size))
        if table is None:
            table = self._table_cache.get(frequency)
        table_valid = (
            table is not None
            and np.isclose(table.frequency_hz, frequency, rtol=0.0, atol=0.0)
            and table.h_nodes.size >= required_nodes
            and float(np.min(effective)) >= float(table.h_nodes[0])
            and float(np.max(effective)) <= float(table.h_nodes[-1])
        )
        rebuilt = not table_valid
        cell_solve_count = 0
        if rebuilt:
            h_min = max(
                float(np.min(effective)) / self.table_range_margin,
                self.h_floor,
            )
            h_max = max(
                float(np.max(effective)) * self.table_range_margin,
                h_min * (1.0 + 1.0e-6),
            )
            if table is not None and np.isclose(
                table.frequency_hz,
                frequency,
                rtol=0.0,
                atol=0.0,
            ):
                h_min = min(h_min, float(table.h_nodes[0]))
                h_max = max(h_max, float(table.h_nodes[-1]))
            h_nodes = np.geomspace(
                h_min,
                h_max,
                required_nodes,
            )
            solver = self._cell_solver(frequency)
            z_nodes = []
            iterations = []
            for h_value in h_nodes:
                result = solver.solve(
                    float(h_value),
                    tol=self.cell_tolerance,
                    max_iter=self.cell_max_iterations,
                    relaxation=self.cell_relaxation,
                )
                if not bool(result["converged"]):
                    raise RuntimeError(
                        "ESIM cell problem did not converge at "
                        f"H_t={float(h_value):.6g} A/m"
                    )
                z_nodes.append(complex(result["Z"]))
                iterations.append(int(result["iterations"]))
            table = _LocalESIMImpedanceTable(
                frequency_hz=frequency,
                h_nodes=h_nodes,
                z_nodes=np.asarray(z_nodes, dtype=complex),
                cell_iterations_max=max(iterations),
                skin_depth_initial_m=float(solver.delta_initial),
            )
            self._table_cache[frequency] = table
            cell_solve_count = int(h_nodes.size)
        assert table is not None
        h_nodes = table.h_nodes
        z_nodes = table.z_nodes
        values = np.interp(effective, h_nodes, z_nodes.real) + 1j * np.interp(
            effective,
            h_nodes,
            z_nodes.imag,
        )
        if not np.all(np.isfinite(values)):
            raise RuntimeError("ESIM cell interpolation produced non-finite impedance")
        if np.min(values.real) < -1.0e-12:
            raise RuntimeError("ESIM cell interpolation produced a non-passive impedance")
        return values, {
            "cell_model": "semi-infinite-geometrically-graded-1d-esim",
            "frequency_Hz": frequency,
            "bins": int(h_nodes.size),
            "cell_iterations_max": int(table.cell_iterations_max),
            "cell_converged": True,
            "cell_table_rebuilt": rebuilt,
            "cell_solve_count": cell_solve_count,
            "skin_depth_initial_m": float(table.skin_depth_initial_m),
            "h_nodes_A_per_m": [float(value) for value in h_nodes],
            "z_nodes_ohm": [
                {"real": float(value.real), "imag": float(value.imag)}
                for value in z_nodes
            ],
            "h_min_A_per_m": float(np.min(fields)),
            "h_max_A_per_m": float(np.max(fields)),
            "table_h_min_A_per_m": float(h_nodes[0]),
            "table_h_max_A_per_m": float(h_nodes[-1]),
            "z_real_min_ohm": float(np.min(values.real)),
            "z_real_max_ohm": float(np.max(values.real)),
            "z_imag_min_ohm": float(np.min(values.imag)),
            "z_imag_max_ohm": float(np.max(values.imag)),
        }, table

    def impedance_samples(
        self,
        tangential_field_magnitude,
        frequency_hz: float,
    ) -> tuple[np.ndarray, dict[str, object]]:
        """Evaluate an attached LUT or build and cache a runtime cell table."""

        values, diagnostics, _ = self._impedance_samples_cached(
            tangential_field_magnitude,
            frequency_hz,
            None,
        )
        return values, diagnostics


def BuildLocalESIMSurfaceLUT(
    model: LocalESIMSurfaceModel,
    frequencies_hz,
    h_nodes_A_per_m,
    *,
    output_path=None,
) -> LocalESIMSurfaceLUT:
    """Precompute a reusable ``Z_s(f, |H_t|)`` database from cell problems.

    All nodes are solved directly; interpolation is used only later by
    :class:`LocalESIMSurfaceLUT`.  The returned table can be saved explicitly
    or atomically written by passing ``output_path``.
    """

    if not isinstance(model, LocalESIMSurfaceModel):
        raise TypeError("model must be a LocalESIMSurfaceModel")
    frequencies = np.asarray(frequencies_hz, dtype=float).reshape(-1)
    fields = np.asarray(h_nodes_A_per_m, dtype=float).reshape(-1)
    if frequencies.size < 1 or not np.all(np.isfinite(frequencies)):
        raise ValueError("frequencies_hz must contain finite positive nodes")
    if fields.size < 2 or not np.all(np.isfinite(fields)):
        raise ValueError("h_nodes_A_per_m must contain at least two finite nodes")
    if np.any(frequencies <= 0.0) or np.any(np.diff(frequencies) <= 0.0):
        raise ValueError("frequencies_hz must be positive and strictly increasing")
    if np.any(fields <= 0.0) or np.any(np.diff(fields) <= 0.0):
        raise ValueError("h_nodes_A_per_m must be positive and strictly increasing")
    impedance = np.empty((frequencies.size, fields.size), dtype=complex)
    iterations = np.empty((frequencies.size, fields.size), dtype=int)
    for frequency_index, frequency in enumerate(frequencies):
        solver = model._cell_solver(float(frequency))
        for field_index, h_value in enumerate(fields):
            result = solver.solve(
                float(h_value),
                tol=model.cell_tolerance,
                max_iter=model.cell_max_iterations,
                relaxation=model.cell_relaxation,
            )
            if not bool(result["converged"]):
                raise RuntimeError(
                    "ESIM cell problem did not converge while building LUT at "
                    f"f={float(frequency):.6g} Hz, H_t={float(h_value):.6g} A/m"
                )
            value = complex(result["Z"])
            if not np.isfinite(value) or value.real < -1.0e-12:
                raise RuntimeError(
                    "ESIM cell problem produced invalid impedance while building LUT at "
                    f"f={float(frequency):.6g} Hz, H_t={float(h_value):.6g} A/m"
                )
            impedance[frequency_index, field_index] = value
            iterations[frequency_index, field_index] = int(result["iterations"])
    lut = LocalESIMSurfaceLUT(
        frequency_nodes_hz=frequencies,
        h_nodes_A_per_m=fields,
        impedance_ohm=impedance,
        cell_iterations=iterations,
        model_signature=model.lut_signature(),
    )
    if output_path is not None:
        lut.save(output_path)
    return lut


def ValidateLocalESIMSurfaceLUT(
    model: LocalESIMSurfaceModel,
    lut: LocalESIMSurfaceLUT,
    *,
    frequencies_hz=None,
    h_values_A_per_m=None,
) -> dict[str, object]:
    """Compare LUT interpolation with direct cell solves at validation points.

    By default, validation uses every LUT frequency node plus each geometric
    frequency midpoint, crossed with every geometric field midpoint.  The
    result is JSON-ready and does not modify the LUT or the model cache.
    """

    if not isinstance(model, LocalESIMSurfaceModel):
        raise TypeError("model must be a LocalESIMSurfaceModel")
    if not isinstance(lut, LocalESIMSurfaceLUT):
        raise TypeError("lut must be a LocalESIMSurfaceLUT")
    lut.validate_model(model)
    if frequencies_hz is None:
        frequencies = np.asarray(lut.frequency_nodes_hz, dtype=float)
        if frequencies.size > 1:
            frequency_midpoints = np.sqrt(frequencies[:-1] * frequencies[1:])
            frequencies = np.sort(np.concatenate((frequencies, frequency_midpoints)))
    else:
        frequencies = np.asarray(frequencies_hz, dtype=float).reshape(-1)
    if h_values_A_per_m is None:
        fields = np.sqrt(lut.h_nodes_A_per_m[:-1] * lut.h_nodes_A_per_m[1:])
    else:
        fields = np.asarray(h_values_A_per_m, dtype=float).reshape(-1)
    if frequencies.size < 1 or fields.size < 1:
        raise ValueError("LUT validation requires frequency and field samples")
    if not np.all(np.isfinite(frequencies)) or np.any(frequencies <= 0.0):
        raise ValueError("validation frequencies must be finite and positive")
    if not np.all(np.isfinite(fields)) or np.any(fields <= 0.0):
        raise ValueError("validation fields must be finite and positive")

    relative_errors = []
    direct_impedances = []
    cell_iterations = []
    for frequency in frequencies:
        predicted = lut.evaluate(fields, float(frequency))
        solver = model._cell_solver(float(frequency))
        direct = np.empty(fields.size, dtype=complex)
        for index, h_value in enumerate(fields):
            result = solver.solve(
                float(h_value),
                tol=model.cell_tolerance,
                max_iter=model.cell_max_iterations,
                relaxation=model.cell_relaxation,
            )
            if not bool(result["converged"]):
                raise RuntimeError(
                    "ESIM cell problem did not converge while validating LUT at "
                    f"f={float(frequency):.6g} Hz, H_t={float(h_value):.6g} A/m"
                )
            direct[index] = complex(result["Z"])
            cell_iterations.append(int(result["iterations"]))
        denominator = np.maximum(np.abs(direct), np.finfo(float).tiny)
        relative_errors.extend(np.abs(predicted - direct) / denominator)
        direct_impedances.extend(direct)
    errors = np.asarray(relative_errors, dtype=float)
    direct_values = np.asarray(direct_impedances, dtype=complex)
    return {
        "model_key_sha256": lut.model_key,
        "interpolation": _LOCAL_ESIM_LUT_INTERPOLATION,
        "sample_count": int(errors.size),
        "frequency_sample_count": int(frequencies.size),
        "field_sample_count": int(fields.size),
        "max_relative_error": float(np.max(errors)),
        "rms_relative_error": float(np.sqrt(np.mean(errors**2))),
        "p95_relative_error": float(np.percentile(errors, 95.0)),
        "cell_iterations_max": int(max(cell_iterations)),
        "direct_cells_passive": bool(np.min(direct_values.real) >= -1.0e-12),
    }


@dataclass(frozen=True)
class LocalESIMSurfaceSolution:
    """Converged nonlinear local-ESIM solution of a reduced HCurl-VIM system."""

    coefficients: np.ndarray
    surface_impedance: SurfaceImpedanceGram
    tangential_field_magnitude: np.ndarray
    frequency_hz: float
    converged: bool
    iterations: int
    relative_impedance_change: float
    history: tuple[dict[str, object], ...]
    linear_solve_diagnostics: dict[str, object]

    def __post_init__(self) -> None:
        coefficients = np.asarray(self.coefficients, dtype=complex).reshape(-1)
        fields = np.asarray(self.tangential_field_magnitude, dtype=float).reshape(-1)
        if not np.all(np.isfinite(coefficients)) or not np.all(np.isfinite(fields)):
            raise ValueError("local ESIM solution contains non-finite values")
        object.__setattr__(self, "coefficients", np.array(coefficients, copy=True))
        object.__setattr__(
            self,
            "tangential_field_magnitude",
            np.array(fields, copy=True),
        )
        object.__setattr__(self, "history", tuple(dict(row) for row in self.history))
        object.__setattr__(
            self,
            "linear_solve_diagnostics",
            dict(self.linear_solve_diagnostics),
        )

    def diagnostics(self) -> dict[str, object]:
        """Return JSON-ready nonlinear convergence and passivity diagnostics."""

        return {
            "frequency_Hz": float(self.frequency_hz),
            "converged": bool(self.converged),
            "iterations": int(self.iterations),
            "relative_impedance_change": float(self.relative_impedance_change),
            "h_min_A_per_m": float(np.min(self.tangential_field_magnitude)),
            "h_max_A_per_m": float(np.max(self.tangential_field_magnitude)),
            "surface_impedance": self.surface_impedance.diagnostics(),
            "linear_solve": dict(self.linear_solve_diagnostics),
            "history": [dict(row) for row in self.history],
        }


def SolveLocalESIMSurfaceVIM(
    system: HybridVIMSystem,
    surface_basis: SampledCurrentBasis,
    rhs,
    model: LocalESIMSurfaceModel,
    frequency_hz: float,
    *,
    block: str = "surface",
    outer_tolerance: float = 1.0e-3,
    outer_max_iterations: int = 20,
    outer_relaxation: float = 0.5,
    field_amplitude_scale: float = 1.0,
    initial_surface_impedance=None,
    solver: str = "cocr",
    solve_tolerance: float = 1.0e-10,
    solve_max_iterations: int | None = None,
    mixed_galerkin_keep_blocks=None,
    mixed_galerkin_eliminate_blocks=None,
    raise_on_nonconvergence: bool = True,
) -> LocalESIMSurfaceSolution:
    """Solve one physical excitation with a local nonlinear ESIM SIBC block.

    Nonlinear superposition is not valid, so ``rhs`` must contain exactly one
    excitation.  The surface-Omega current is the tangential magnetic-field
    trace (up to ``field_amplitude_scale`` for the caller's phasor convention).
    """

    if not isinstance(system, HybridVIMSystem):
        raise TypeError("system must be a HybridVIMSystem")
    if not isinstance(surface_basis, SampledCurrentBasis):
        raise TypeError("surface_basis must be a SampledCurrentBasis")
    if not isinstance(model, LocalESIMSurfaceModel):
        raise TypeError("model must be a LocalESIMSurfaceModel")
    frequency = float(frequency_hz)
    if not np.isfinite(frequency) or frequency <= 0.0:
        raise ValueError("frequency_hz must be positive")
    values = np.asarray(rhs)
    if values.ndim == 2 and values.shape[1] == 1:
        values = values[:, 0]
    if values.ndim != 1 or values.shape != (system.n_modes,):
        raise ValueError(
            "rhs must contain one nonlinear excitation with shape "
            f"({system.n_modes},)"
        )
    if not np.all(np.isfinite(values)):
        raise ValueError("rhs contains non-finite values")
    outer_tolerance = float(outer_tolerance)
    outer_relaxation = float(outer_relaxation)
    field_amplitude_scale = float(field_amplitude_scale)
    outer_max_iterations = int(outer_max_iterations)
    if not np.isfinite(outer_tolerance) or outer_tolerance <= 0.0:
        raise ValueError("outer_tolerance must be positive")
    if outer_max_iterations < 1:
        raise ValueError("outer_max_iterations must be >= 1")
    if not np.isfinite(outer_relaxation) or not 0.0 < outer_relaxation <= 1.0:
        raise ValueError("outer_relaxation must be in (0, 1]")
    if not np.isfinite(field_amplitude_scale) or field_amplitude_scale <= 0.0:
        raise ValueError("field_amplitude_scale must be positive")
    use_mixed = (
        mixed_galerkin_keep_blocks is not None
        or mixed_galerkin_eliminate_blocks is not None
    )
    if (mixed_galerkin_keep_blocks is None) != (
        mixed_galerkin_eliminate_blocks is None
    ):
        raise ValueError(
            "mixed_galerkin_keep_blocks and mixed_galerkin_eliminate_blocks "
            "must be supplied together"
        )
    surface_slice = system.block_slice(block)
    if surface_slice.stop - surface_slice.start != surface_basis.n_modes:
        raise ValueError("surface basis mode count does not match the selected block")
    s = 1j * 2.0 * np.pi * frequency
    if initial_surface_impedance is None:
        initial_surface_impedance = model.linear_impedance(frequency)
    if isinstance(initial_surface_impedance, SurfaceImpedanceGram):
        _surface_impedance_term(system.surface_mass, initial_surface_impedance)
        current_gram = initial_surface_impedance
    elif np.isscalar(initial_surface_impedance):
        current_gram = AssembleSurfaceImpedanceGram(
            system,
            surface_basis,
            complex(initial_surface_impedance),
            block=block,
            label="local-esim-initial-linear",
        )
    else:
        raise TypeError(
            "initial_surface_impedance must be scalar or SurfaceImpedanceGram"
        )

    history = []
    coefficients = None
    fields = None
    solved_gram = current_gram
    solve_info = {}
    impedance_table = None
    relative_change = float("inf")
    converged = False
    for iteration in range(1, outer_max_iterations + 1):
        solved_gram = current_gram
        if use_mixed:
            matrix = system.impedance(s, surface_impedance=solved_gram)
            orthogonalization = system.mixed_galerkin_orthogonalization(
                mixed_galerkin_keep_blocks,
                mixed_galerkin_eliminate_blocks,
                s,
                surface_impedance=solved_gram,
            )
            coefficients = orthogonalization.solve(matrix, values)
            solve_info = {
                "backend": "dense-mixed-galerkin-schur",
                **orthogonalization.diagnostics(),
            }
        else:
            coefficients, solve_info = system.solve(
                s,
                values,
                surface_impedance=solved_gram,
                solver=solver,
                tolerance=solve_tolerance,
                max_iterations=solve_max_iterations,
                return_diagnostics=True,
            )
        surface_coefficients = coefficients[surface_slice]
        surface_current = np.einsum(
            "m,mik->ik",
            surface_coefficients,
            surface_basis.modes,
        )
        fields = field_amplitude_scale * np.sqrt(
            np.sum(np.abs(surface_current) ** 2, axis=1)
        )
        target_values, cell_info, impedance_table = model._impedance_samples_cached(
            fields,
            frequency,
            impedance_table,
        )
        target_gram = AssembleSurfaceImpedanceGram(
            system,
            surface_basis,
            target_values,
            block=block,
            label="local-esim",
        )
        denominator = max(
            float(np.linalg.norm(target_values)),
            np.finfo(float).tiny,
        )
        relative_change = float(
            np.linalg.norm(target_values - solved_gram.sample_values) / denominator
        )
        history.append(
            {
                "iteration": iteration,
                "relative_impedance_change": relative_change,
                **cell_info,
            }
        )
        if relative_change <= outer_tolerance:
            converged = True
            break
        relaxed_values = (
            (1.0 - outer_relaxation) * solved_gram.sample_values
            + outer_relaxation * target_values
        )
        current_gram = AssembleSurfaceImpedanceGram(
            system,
            surface_basis,
            relaxed_values,
            block=block,
            label="local-esim-relaxed",
        )

    if not converged and raise_on_nonconvergence:
        raise RuntimeError(
            "local ESIM outer iteration did not converge: "
            f"relative change={relative_change:.3e} after {outer_max_iterations} iterations"
        )
    assert coefficients is not None and fields is not None
    return LocalESIMSurfaceSolution(
        coefficients=coefficients,
        surface_impedance=solved_gram,
        tangential_field_magnitude=fields,
        frequency_hz=frequency,
        converged=converged,
        iterations=iteration,
        relative_impedance_change=relative_change,
        history=tuple(history),
        linear_solve_diagnostics=solve_info,
    )


def _rhs_matrix_for_basis(basis: SampledCurrentBasis, vector_potentials) -> np.ndarray:
    ports = tuple(vector_potentials)
    if not ports:
        raise ValueError("vector_potentials must not be empty")
    values = [port(basis.points) if callable(port) else port for port in ports]
    return np.column_stack([ExternalVectorPotentialRHS(basis, value) for value in values])


def _topology_sibc_boundary_labels(topology: EddyMeshTopology) -> tuple[str, ...]:
    labels = {
        label
        for face in topology.sibc_faces
        for label in face.boundary_labels
    }
    return tuple(sorted(labels))


@dataclass(frozen=True)
class TopologyAwareHybridVIM:
    """Production-facing tri-block HCurl-VIM assembly result.

    The wrapped system has block names:

    ``volume``
        EVRS bulk current basis sampled from ``curl(T)``.
    ``volume1``
        conductor-graph bridge-cycle current basis.
    ``surface``
        surface-Omega/SIBC current basis on air-touching conductor faces.
    """

    system: HybridVIMSystem
    volume_basis: SampledCurrentBasis
    bridge_cycle_basis: SampledCurrentBasis
    surface_basis: SampledCurrentBasis
    topology: EddyMeshTopology
    conductor_graph: EddyConductorGraph
    dof_policy: EddyDofPolicy
    reduction_plan: EddyReductionPlan
    rhs: np.ndarray | None = None
    parent_order: int | None = None
    parent_order_ledger: EddyParentOrderLedger | None = None
    response_basis: ResponseBasis | None = None
    conductivity: float | None = None
    sibc_applicability: EddySIBCApplicability | None = None

    def diagnostics(self) -> dict[str, object]:
        """Return complete production diagnostics for the tri-block VIM."""

        return {
            "system": self.system.diagnostics(),
            "volume_basis": _sampled_current_basis_diagnostics(self.volume_basis),
            "bridge_cycle_basis": _sampled_current_basis_diagnostics(self.bridge_cycle_basis),
            "surface_basis": _sampled_current_basis_diagnostics(self.surface_basis),
            "topology": self.topology.diagnostics(),
            "conductor_graph": self.conductor_graph.diagnostics(),
            "dof_policy": self.dof_policy.diagnostics(),
            "reduction_plan": self.reduction_plan.diagnostics(),
            "eddy_bubbling": self.eddy_bubble_decomposition().diagnostics(),
            "has_rhs": self.rhs is not None,
            "has_response_basis": self.response_basis is not None,
            "conductivity": self.conductivity,
            "surface_model": (
                "explicit"
                if self.sibc_applicability is None
                else self.sibc_applicability.selected_model
            ),
            "sibc_applicability": (
                None
                if self.sibc_applicability is None
                else self.sibc_applicability.diagnostics()
            ),
        }

    def eddy_bubble_decomposition(self) -> EddyBubbleDecomposition:
        """Return the topology-aware eddy-bubbling class split."""

        return EddyBubbleDecomposition(
            policy=self.dof_policy,
            plan=self.reduction_plan,
            topology=self.topology,
            conductor_graph=self.conductor_graph,
            parent_order=self.parent_order,
            parent_order_ledger=self.parent_order_ledger,
        )

    def port_admittance(self, s, *, surface_impedance=0.0) -> np.ndarray:
        """Evaluate port admittance using the stored RHS."""

        if self.rhs is None:
            raise ValueError("rhs was not assembled; pass port_vector_potentials")
        return self.system.port_admittance(
            s,
            self.rhs,
            surface_impedance=surface_impedance,
        )

    def surface_impedance_gram(
        self,
        surface_impedance,
        *,
        label: str = "local-sibc",
    ) -> SurfaceImpedanceGram:
        """Assemble local linear/ESIM values on the air-facing SIBC block."""

        return AssembleSurfaceImpedanceGram(
            self.system,
            self.surface_basis,
            surface_impedance,
            block="surface",
            label=label,
        )

    def solve_local_esim(
        self,
        model: LocalESIMSurfaceModel,
        frequency_hz: float,
        rhs=None,
        **kwargs,
    ) -> LocalESIMSurfaceSolution:
        """Run the production local-ESIM Karl loop on the air-facing block."""

        if rhs is None:
            if self.rhs is None:
                raise ValueError("rhs was not assembled; supply one excitation")
            rhs = self.rhs
        return SolveLocalESIMSurfaceVIM(
            self.system,
            self.surface_basis,
            rhs,
            model,
            frequency_hz,
            block="surface",
            **kwargs,
        )

    def cln_model(self) -> HCurlEddyCLNModel:
        """Return the CLN/EVRS descriptor for the assembled physical ports."""

        if self.rhs is None:
            raise ValueError("rhs was not assembled; pass port_vector_potentials")
        return HCurlEddyCLNFromVIM(self.system, self.rhs)

    def couple_hdiv_mmm(
        self,
        magnetization_basis: SampledMagnetizationBasis,
        *,
        material_model: object | None = None,
        magnetic_operator=None,
        magnetic_rhs=None,
        hdiv_reduction: HDivMMMReducedModel | None = None,
        mu: float = MU0,
        kernel_epsilon: float = 0.0,
    ) -> CoupledHDivHybridVIMSystem:
        """Couple this tri-block HCurl-VIM system to an HDiv-MMM branch."""

        return CoupleHybridVIMWithHDivMMM(
            magnetization_basis,
            self.system,
            (self.volume_basis, self.bridge_cycle_basis, self.surface_basis),
            material_model=material_model,
            magnetic_operator=magnetic_operator,
            magnetic_rhs=magnetic_rhs,
            eddy_rhs=self.rhs,
            response_basis=self.response_basis,
            eddy_bubbling=self.eddy_bubble_decomposition(),
            hdiv_reduction=hdiv_reduction,
            conductivity=self.conductivity,
            eddy_block_roles={
                "volume": "bulk",
                "volume1": "bridge",
                "surface": "sibc",
            },
            mu=mu,
            kernel_epsilon=kernel_epsilon,
        )


def NgsolveTopologyAwareHybridVIM(
    mesh,
    fes,
    response_vectors,
    surface_grad_modes,
    *,
    sigma: float,
    conductive_materials,
    air_materials=("air", "vacuum"),
    volume_materials=None,
    surface_boundaries=None,
    intorder: int = 2,
    current_gram_rtol: float = 1.0e-10,
    geometry_intorder: int | None = None,
    kernel_epsilon: float | None = None,
    topology: EddyMeshTopology | None = None,
    free_dofs=None,
    volume_names=None,
    bridge_names=None,
    surface_names=None,
    non_sibc_trace_modes: int | None = None,
    parent_order: int | None = None,
    parent_order_ledger: EddyParentOrderLedger | None = None,
    sibc_applicability: EddySIBCApplicability | None = None,
    port_vector_potentials=None,
    interaction=None,
) -> TopologyAwareHybridVIM:
    """Build the production tri-block topology-aware HCurl VIM system.

    This is the high-level entry point for the current research path:

    ``EVRS curl(T) bulk`` + ``bridge-cycle VIM`` + ``surface-Omega/SIBC``.

    It deliberately assembles a VIM system.  The optional ``interaction`` hook
    only replaces the Laplace/VIM interaction backend; it does not change the
    formulation into a BEM method.
    """

    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    if sibc_applicability is not None:
        if not isinstance(sibc_applicability, EddySIBCApplicability):
            raise TypeError("sibc_applicability must be an EddySIBCApplicability")
        if not np.isclose(sibc_applicability.sigma, sigma):
            raise ValueError("sibc_applicability.sigma must match sigma")
    surface_grad_modes = tuple(surface_grad_modes)
    if sibc_applicability is not None and not sibc_applicability.sibc_applicable:
        surface_grad_modes = ()
        surface_names = None
    if topology is None:
        topology = ClassifyNgsolveEddyTopology(
            mesh,
            conductive_materials,
            air_materials=air_materials,
        )
    if volume_materials is None:
        volume_materials = conductive_materials
    if surface_boundaries is None:
        labels = _topology_sibc_boundary_labels(topology)
        surface_boundaries = labels if labels else None

    response_basis = response_vectors if isinstance(response_vectors, ResponseBasis) else None
    vectors = response_basis.vectors if response_basis is not None else response_vectors
    volume = NgsolveHCurlCurlBasis(
        mesh,
        fes,
        vectors,
        intorder=intorder,
        materials=volume_materials,
        names=volume_names,
    )
    if response_basis is not None:
        response_basis, volume = CompressHCurlResponseInCurrentGram(
            response_basis,
            volume,
            rtol=current_gram_rtol,
        )
        vectors = response_basis.vectors
    graph = topology.conductor_graph()
    if bridge_names is None:
        bridge_names = [f"bridge_cycle_{i}" for i in range(graph.cycle_rank)]
    bridge = NgsolveBridgeCycleCurrentBasis(
        mesh,
        topology,
        geometry_intorder=geometry_intorder,
        names=bridge_names,
    )
    surface = NgsolveSurfaceOmegaBasis(
        mesh,
        surface_grad_modes,
        intorder=intorder,
        boundaries=surface_boundaries,
        names=surface_names,
    )
    if interaction is None:
        sampled_epsilon = (
            _default_kernel_epsilon((volume, bridge, surface))
            if kernel_epsilon is None
            else float(kernel_epsilon)
        )
        if int(mesh.dim) == 3:
            from ._hcurl_tet_interaction import (
                NgsolveHCurlCellVolumeInteraction,
            )

            volume_interaction = NgsolveHCurlCellVolumeInteraction(
                mesh,
                fes,
                vectors,
                volume,
                materials=volume_materials,
                matrix_free=True,
            )
            interaction = HACApKSampledLaplaceInteraction(
                kernel_epsilon=sampled_epsilon,
                cross_only=True,
                diagonal_interaction=volume_interaction,
                diagonal_bases=(volume,),
            )
        elif int(mesh.dim) == 2:
            from ._hcurl_planar_interaction import (
                NgsolveHCurlPlanarVolumeInteraction,
            )

            volume_interaction = NgsolveHCurlPlanarVolumeInteraction(
                mesh,
                fes,
                vectors,
                volume,
                materials=volume_materials,
                matrix_free=True,
            )
            interaction = HACApKSampledPlanarLogInteraction(
                kernel_epsilon=sampled_epsilon,
                reference_length=1.0,
                cross_only=True,
                diagonal_interaction=volume_interaction,
                diagonal_bases=(volume,),
            )
        else:
            raise ValueError("HCurl hybrid VIM supports only 2-D or 3-D meshes")
    system = AssembleHybridVIM(
        volume,
        bridge,
        surface,
        sigma=sigma,
        kernel_epsilon=kernel_epsilon,
        interaction=interaction,
    )
    dof_policy = NgsolveEddyDofPolicy(
        mesh,
        fes,
        topology,
        free_dofs=free_dofs,
    )
    reduction_plan = dof_policy.reduction_plan(
        evrs_rank=volume.n_modes,
        surface_modes=surface.n_modes,
        non_sibc_trace_modes=non_sibc_trace_modes,
        loop_bridge_modes=graph.cycle_rank,
        bridge_strategy="cycle-basis",
    )
    rhs = None
    if port_vector_potentials is not None:
        rhs = system.block_rhs(
            volume=_rhs_matrix_for_basis(volume, port_vector_potentials),
            volume1=_rhs_matrix_for_basis(bridge, port_vector_potentials),
            surface=_rhs_matrix_for_basis(surface, port_vector_potentials),
        )
    return TopologyAwareHybridVIM(
        system=system,
        volume_basis=volume,
        bridge_cycle_basis=bridge,
        surface_basis=surface,
        topology=topology,
        conductor_graph=graph,
        dof_policy=dof_policy,
        reduction_plan=reduction_plan,
        rhs=rhs,
        parent_order=(
            _ngsolve_fes_order(fes) if parent_order is None else parent_order
        ),
        parent_order_ledger=parent_order_ledger,
        response_basis=response_basis,
        conductivity=float(sigma),
        sibc_applicability=sibc_applicability,
    )


def NgsolveEddyBubbleHybridVIM(
    mesh,
    fes,
    stiffness,
    mass,
    ports,
    surface_grad_modes,
    *,
    steps: int,
    sigma: float,
    conductive_materials,
    air_materials=("air", "vacuum"),
    volume_materials=None,
    surface_boundaries=None,
    intorder: int = 2,
    geometry_intorder: int | None = None,
    kernel_epsilon: float | None = None,
    topology: EddyMeshTopology | None = None,
    free_dofs=None,
    condense: bool = False,
    response_backend: str = "auto",
    inverse: str = "sparsecholesky",
    rtol: float = 1.0e-10,
    current_gram_rtol: float = 1.0e-10,
    volume_names=None,
    bridge_names=None,
    surface_names=None,
    non_sibc_trace_modes: int | None = None,
    parent_order: int | None = None,
    parent_order_ledger: EddyParentOrderLedger | None = None,
    sibc_applicability: EddySIBCApplicability | None = None,
    port_vector_potentials=None,
    interaction=None,
) -> TopologyAwareHybridVIM:
    """Build the full production eddy-bubbled hybrid HCurl-VIM system.

    This is the one-call station/stator-side builder.  It generates the EVRS
    response basis from the high-order HCurl parent matrices, then assembles the
    topology-aware tri-block VIM system:

    ``bulk EVRS curl(T)`` + ``bridge-cycle current`` + ``surface-Omega/SIBC``.
    """

    if steps < 1:
        raise ValueError("steps must be >= 1")
    if topology is None:
        topology = ClassifyNgsolveEddyTopology(
            mesh,
            conductive_materials,
            air_materials=air_materials,
        )
    if free_dofs is None:
        free_dofs = fes.FreeDofs(bool(condense))
    response = _ngsolve_response_basis_for_eddy_bubbling(
        stiffness,
        mass,
        ports,
        steps=steps,
        free_dofs=free_dofs,
        condense=condense,
        response_backend=response_backend,
        inverse=inverse,
        rtol=rtol,
    )
    if volume_names is None:
        p = _ngsolve_fes_order(fes) if parent_order is None else parent_order
        prefix = f"evrs_p{p}_n{steps}" if p is not None else f"evrs_n{steps}"
        volume_names = [f"{prefix}_{i}" for i in range(response.rank)]
    return NgsolveTopologyAwareHybridVIM(
        mesh,
        fes,
        response,
        surface_grad_modes,
        sigma=sigma,
        conductive_materials=conductive_materials,
        air_materials=air_materials,
        volume_materials=volume_materials,
        surface_boundaries=surface_boundaries,
        intorder=intorder,
        current_gram_rtol=current_gram_rtol,
        geometry_intorder=geometry_intorder,
        kernel_epsilon=kernel_epsilon,
        topology=topology,
        free_dofs=free_dofs,
        volume_names=volume_names,
        bridge_names=bridge_names,
        surface_names=surface_names,
        non_sibc_trace_modes=non_sibc_trace_modes,
        parent_order=parent_order,
        parent_order_ledger=parent_order_ledger,
        sibc_applicability=sibc_applicability,
        port_vector_potentials=port_vector_potentials,
        interaction=interaction,
    )


def NgsolveHCurlVIMHDivMMM(
    mesh,
    hcurl_fes,
    stiffness,
    mass,
    ports,
    surface_grad_modes,
    magnetization_basis: SampledMagnetizationBasis | HDivMMMReducedModel,
    *,
    steps: int,
    sigma: float,
    conductive_materials,
    air_materials=("air", "vacuum"),
    volume_materials=None,
    surface_boundaries=None,
    intorder: int = 2,
    geometry_intorder: int | None = None,
    kernel_epsilon: float | None = None,
    topology: EddyMeshTopology | None = None,
    free_dofs=None,
    condense: bool = False,
    response_backend: str = "auto",
    inverse: str = "sparsecholesky",
    rtol: float = 1.0e-10,
    current_gram_rtol: float = 1.0e-10,
    volume_names=None,
    bridge_names=None,
    surface_names=None,
    non_sibc_trace_modes: int | None = None,
    parent_order: int | None = None,
    parent_order_ledger: EddyParentOrderLedger | None = None,
    sibc_applicability: EddySIBCApplicability | None = None,
    port_vector_potentials=None,
    interaction=None,
    material_model: object | None = None,
    magnetic_operator=None,
    magnetic_rhs=None,
    mu: float = MU0,
    coupling_kernel_epsilon: float = 0.0,
) -> CoupledHDivHybridVIMSystem:
    """Build the one-call HCurl-VIM / HDiv-MMM production mixed system."""

    if getattr(hcurl_fes, "mesh", mesh) is not mesh:
        raise ValueError("hcurl_fes must be built on the shared production mesh")
    hdiv_reduction = None
    if isinstance(magnetization_basis, HDivMMMReducedModel):
        hdiv_reduction = magnetization_basis
        if hdiv_reduction.mesh is not mesh:
            raise ValueError(
                "the HDiv-MMM reduction and HCurl-VIM branch must share mesh"
            )
        magnetization_basis = hdiv_reduction.magnetization_basis
        if magnetic_operator is None:
            magnetic_operator = hdiv_reduction.magnetic_operator
        if magnetic_rhs is None:
            magnetic_rhs = hdiv_reduction.magnetic_rhs
    elif not isinstance(magnetization_basis, SampledMagnetizationBasis):
        raise TypeError(
            "magnetization_basis must be a SampledMagnetizationBasis or "
            "HDivMMMReducedModel"
        )
    if (
        isinstance(material_model, SharedMeshMaterialModel)
        and material_model.mesh is not mesh
    ):
        raise ValueError("material_model must reference the shared production mesh")

    hybrid = NgsolveEddyBubbleHybridVIM(
        mesh,
        hcurl_fes,
        stiffness,
        mass,
        ports,
        surface_grad_modes,
        steps=steps,
        sigma=sigma,
        conductive_materials=conductive_materials,
        air_materials=air_materials,
        volume_materials=volume_materials,
        surface_boundaries=surface_boundaries,
        intorder=intorder,
        geometry_intorder=geometry_intorder,
        kernel_epsilon=kernel_epsilon,
        topology=topology,
        free_dofs=free_dofs,
        condense=condense,
        response_backend=response_backend,
        inverse=inverse,
        rtol=rtol,
        current_gram_rtol=current_gram_rtol,
        volume_names=volume_names,
        bridge_names=bridge_names,
        surface_names=surface_names,
        non_sibc_trace_modes=non_sibc_trace_modes,
        parent_order=parent_order,
        parent_order_ledger=parent_order_ledger,
        sibc_applicability=sibc_applicability,
        port_vector_potentials=port_vector_potentials,
        interaction=interaction,
    )
    return hybrid.couple_hdiv_mmm(
        magnetization_basis,
        material_model=material_model,
        magnetic_operator=magnetic_operator,
        magnetic_rhs=magnetic_rhs,
        hdiv_reduction=hdiv_reduction,
        mu=mu,
        kernel_epsilon=coupling_kernel_epsilon,
    )


def NgsolveBDMEddyBubbleVIM(
    mesh,
    hcurl_fes,
    stiffness,
    mass,
    ports,
    surface_grad_modes,
    *,
    hdiv_order: int = 1,
    mu_r: float,
    external_fields=None,
    training_fields=None,
    external_names=None,
    training_names=None,
    hdiv_port_weights=None,
    hdiv_normalize_ports: bool = True,
    hdiv_pod_rtol: float = 1.0e-10,
    hdiv_max_modes=None,
    hdiv_solve_tol: float = 1.0e-10,
    hdiv_solve_maxit: int = 5000,
    hdiv_inverse: str = "sparsecholesky",
    hdiv_intorder: int = 2,
    magnetic_materials=None,
    hdiv_mass=None,
    demag_operator=None,
    demag_intorder=None,
    demag_eps: float = 1.0e-7,
    demag_leafsize: int = 16,
    demag_eta: float = 2.0,
    demag_far_quad: int = 3,
    demag_ho_far_factor: float = 2.0,
    steps: int,
    sigma: float,
    conductive_materials,
    air_materials=("air", "vacuum"),
    volume_materials=None,
    surface_boundaries=None,
    intorder: int = 2,
    geometry_intorder: int | None = None,
    kernel_epsilon: float | None = None,
    topology: EddyMeshTopology | None = None,
    free_dofs=None,
    condense: bool = False,
    response_backend: str = "auto",
    inverse: str = "sparsecholesky",
    rtol: float = 1.0e-10,
    current_gram_rtol: float = 1.0e-10,
    volume_names=None,
    bridge_names=None,
    surface_names=None,
    non_sibc_trace_modes: int | None = None,
    parent_order: int | None = None,
    parent_order_ledger: EddyParentOrderLedger | None = None,
    sibc_applicability: EddySIBCApplicability | None = None,
    port_vector_potentials=None,
    interaction=None,
    material_model: SharedMeshMaterialModel | None = None,
    mu: float = MU0,
    coupling_kernel_epsilon: float = 0.0,
) -> CoupledHDivHybridVIMSystem:
    """Build the production BDM-MMM plus eddy-bubble HCurl-VIM system.

    One shared NGSolve mesh feeds both branches.  The magnetic parent is always
    bare ``HDiv(order=hdiv_order)`` (BDM on simplex cells); the eddy parent is
    the supplied high-order HCurl space and is reduced class by class into bulk
    EVRS, conductor-cycle bridge, and air-facing SIBC modes.  Physical applied-H
    responses are protected before POD compression.
    """

    if magnetic_materials is None:
        magnetic_materials = conductive_materials
    if material_model is None:
        material_model = SharedMeshMaterialModel(
            mesh=mesh,
            magnetic_regions=magnetic_materials,
            conductive_regions=conductive_materials,
            mu=float(mu) * float(mu_r),
            sigma=sigma,
            sibc=(
                "half-space"
                if sibc_applicability is None or sibc_applicability.sibc_applicable
                else None
            ),
            metadata={
                "hdiv_family": "BDM",
                "hdiv_order": int(hdiv_order),
                "eddy_reduction": "topology-aware-eddy-bubble",
                "eddy_surface_model": (
                    "explicit"
                    if sibc_applicability is None
                    else sibc_applicability.selected_model
                ),
            },
        )
    elif not isinstance(material_model, SharedMeshMaterialModel):
        raise TypeError("material_model must be a SharedMeshMaterialModel or None")
    elif material_model.mesh is not mesh:
        raise ValueError("material_model must reference the shared production mesh")

    hdiv_reduction = NgsolveBDMHDivMMMResponseReduction(
        mesh,
        order=hdiv_order,
        mu_r=mu_r,
        external_fields=external_fields,
        training_fields=training_fields,
        external_names=external_names,
        training_names=training_names,
        port_weights=hdiv_port_weights,
        normalize_ports=hdiv_normalize_ports,
        pod_rtol=hdiv_pod_rtol,
        max_modes=hdiv_max_modes,
        solve_tol=hdiv_solve_tol,
        solve_maxit=hdiv_solve_maxit,
        inverse=hdiv_inverse,
        intorder=hdiv_intorder,
        materials=magnetic_materials,
        mass=hdiv_mass,
        demag_operator=demag_operator,
        demag_intorder=demag_intorder,
        demag_eps=demag_eps,
        demag_leafsize=demag_leafsize,
        demag_eta=demag_eta,
        demag_far_quad=demag_far_quad,
        demag_ho_far_factor=demag_ho_far_factor,
    )
    return NgsolveHCurlVIMHDivMMM(
        mesh,
        hcurl_fes,
        stiffness,
        mass,
        ports,
        surface_grad_modes,
        hdiv_reduction,
        steps=steps,
        sigma=sigma,
        conductive_materials=conductive_materials,
        air_materials=air_materials,
        volume_materials=volume_materials,
        surface_boundaries=surface_boundaries,
        intorder=intorder,
        geometry_intorder=geometry_intorder,
        kernel_epsilon=kernel_epsilon,
        topology=topology,
        free_dofs=free_dofs,
        condense=condense,
        response_backend=response_backend,
        inverse=inverse,
        rtol=rtol,
        current_gram_rtol=current_gram_rtol,
        volume_names=volume_names,
        bridge_names=bridge_names,
        surface_names=surface_names,
        non_sibc_trace_modes=non_sibc_trace_modes,
        parent_order=parent_order,
        parent_order_ledger=parent_order_ledger,
        sibc_applicability=sibc_applicability,
        port_vector_potentials=port_vector_potentials,
        interaction=interaction,
        material_model=material_model,
        mu=mu,
        coupling_kernel_epsilon=coupling_kernel_epsilon,
    )


def SkinImpedance(s, sigma: float, mu: float = MU0):
    """Half-space SIBC impedance ``Zs = sqrt(mu*s/sigma)``.

    The principal square-root branch gives positive real part for
    ``s = 1j*omega`` with ``omega > 0``.
    """

    s = complex(s)
    sigma = float(sigma)
    mu = float(mu)
    if not (np.isfinite(s.real) and np.isfinite(s.imag)):
        raise ValueError("s must be finite")
    if not np.isfinite(sigma) or sigma <= 0.0:
        raise ValueError("sigma must be positive")
    if not np.isfinite(mu) or mu <= 0.0:
        raise ValueError("mu must be positive")
    func = _radia_cpp_kernel("_SkinImpedance")
    if func is not None:
        return func(s, sigma, mu)
    return np.sqrt(mu * s / sigma)


def SIBCAdmittanceTail(s, surface_measure: float, sigma: float, mu: float = MU0):
    """Return the leading SIBC admittance tail ``S sqrt(sigma/(mu s))``.

    ``surface_measure`` is perimeter for 2-D per-unit-length reductions and
    surface area for 3-D bodies, matching the IGTE digest convention.
    """

    s = complex(s)
    surface_measure = float(surface_measure)
    sigma = float(sigma)
    mu = float(mu)
    if not (np.isfinite(s.real) and np.isfinite(s.imag)):
        raise ValueError("s must be finite")
    if s == 0.0:
        raise ValueError("SIBC admittance tail is undefined at s=0")
    if not np.isfinite(surface_measure) or surface_measure <= 0.0:
        raise ValueError("surface_measure must be positive")
    if not np.isfinite(sigma) or sigma <= 0.0:
        raise ValueError("sigma must be positive")
    if not np.isfinite(mu) or mu <= 0.0:
        raise ValueError("mu must be positive")
    func = _radia_cpp_kernel("_SIBCAdmittanceTail")
    if func is not None:
        return func(s, surface_measure, sigma, mu)
    return surface_measure * np.sqrt(sigma / (mu * s))


def SIBCSchurTerminationImpedance(s, k_sibc: float, d: float = 0.0):
    """Return the scalar Schur/Warburg surface-port impedance.

    The supplement writes the Schur block as
    ``z(s) = (s + d) / (K_SIBC sqrt(s))``.  Its inverse tends to the SIBC
    admittance tail ``K_SIBC / sqrt(s)`` at high frequency while vanishing at
    DC when ``d > 0``.
    """

    s = complex(s)
    k_sibc = float(k_sibc)
    d = float(d)
    if not (np.isfinite(s.real) and np.isfinite(s.imag)):
        raise ValueError("s must be finite")
    if not np.isfinite(k_sibc) or k_sibc <= 0.0:
        raise ValueError("k_sibc must be positive")
    if not np.isfinite(d) or d < 0.0:
        raise ValueError("d must be non-negative")
    if s == 0.0:
        if d == 0.0:
            return 0.0j
        raise ValueError("SIBC termination impedance has a pole at s=0 when d>0")
    func = _radia_cpp_kernel("_SIBCSchurTerminationImpedance")
    if func is not None:
        return func(s, k_sibc, d)
    return (s + d) / (k_sibc * np.sqrt(s))


def SIBCSchurTerminationAdmittance(s, k_sibc: float, d: float = 0.0):
    """Return the inverse of :func:`SIBCSchurTerminationImpedance`."""

    s = complex(s)
    k_sibc = float(k_sibc)
    d = float(d)
    if not (np.isfinite(s.real) and np.isfinite(s.imag)):
        raise ValueError("s must be finite")
    if not np.isfinite(k_sibc) or k_sibc <= 0.0:
        raise ValueError("k_sibc must be positive")
    if not np.isfinite(d) or d < 0.0:
        raise ValueError("d must be non-negative")
    if s == 0.0:
        if d > 0.0:
            return 0.0j
        raise ValueError("SIBC termination admittance has a pole at s=0 when d=0")
    func = _radia_cpp_kernel("_SIBCSchurTerminationAdmittance")
    if func is not None:
        return func(s, k_sibc, d)
    return 1.0 / SIBCSchurTerminationImpedance(s, k_sibc, d=d)


def ExternalVectorPotentialRHS(basis: SampledCurrentBasis, vector_potential) -> np.ndarray:
    """Project an external vector potential onto sampled current modes.

    The returned entry is ``int mode_i dot A_ext dV`` or ``dS`` depending on
    the basis kind.  ``vector_potential`` can be a constant 3-vector or samples
    with shape ``(n, 3)``.
    """

    aext = np.asarray(vector_potential)
    if aext.shape == (3,):
        aext = np.broadcast_to(aext, basis.points.shape)
    if aext.shape != basis.points.shape:
        raise ValueError(
            f"vector_potential must have shape (3,) or {basis.points.shape}"
        )
    if not np.all(np.isfinite(aext)):
        raise ValueError("vector_potential contains non-finite values")
    return np.einsum("aik,ik,i->a", basis.modes.conj(), aext, basis.weights)


__all__ = [
    "MU0",
    "SkinDepth",
    "EddySIBCApplicability",
    "EddyTracePolynomialDim",
    "EddyParentOrderLedger",
    "SampledCurrentBasis",
    "SampledMagnetizationBasis",
    "VolumeCurrentBasis",
    "MagnetizationBasis",
    "SurfaceOmegaBasis",
    "EddyFaceTopology",
    "EddyConductorGraphEdge",
    "EddyConductorCycle",
    "EddyConductorGraph",
    "EddyMeshTopology",
    "EddyDofPolicy",
    "EddyReductionPlan",
    "EddyBubbleDecomposition",
    "EddyBubbleHCurlBasis",
    "HCurlCellFamilyInventory",
    "EddyBubbleReduction",
    "NgsolveHCurlCellFamilies",
    "ClassifyNgsolveEddyTopology",
    "NgsolveEddyDofPolicy",
    "NgsolveEddyBubbleReduction",
    "NgsolveEddyBubbleHCurlBasis",
    "NgsolveBridgeCycleCurrentBasis",
    "SampleNgsolveVectorCFs",
    "NgsolveVolumeCurrentBasis",
    "NgsolveMagnetizationBasis",
    "NgsolveHDivMagnetizationBasis",
    "HDivMultipolePortSet",
    "PlanarHarmonicPortSet",
    "NgsolveHDivRegularSolidHarmonicPorts",
    "NgsolvePlanarHarmonicPorts",
    "NgsolveHDivExternalFieldRHS",
    "HDivMMMReducedModel",
    "NgsolveHDivMMMReduction",
    "NgsolveHDivMMMResponseReduction",
    "NgsolveBDMHDivMMMResponseReduction",
    "PlanarHDivMMMReducedSolution",
    "PlanarHDivMMMReducedModel",
    "NgsolvePlanarHDivMMMResponseReduction",
    "NgsolveHCurlCurlBasis",
    "NgsolveSurfaceOmegaBasis",
    "NgsolveMatrixToDense",
    "NgsolveVectorToArray",
    "NgsolveCouplingDofMasks",
    "ResponseBasis",
    "EVRSBasis",
    "CompressHCurlResponseInCurrentGram",
    "BlockKrylovBasis",
    "NgsolveBlockKrylovBasis",
    "NgsolveOperatorBlockKrylovBasis",
    "NgsolveStaticCondensedBlockKrylovBasis",
    "SampledLaplaceInteraction",
    "HACApKSampledLaplaceInteraction",
    "HACApKSampledPlanarLogInteraction",
    "SampledHACApKOperator",
    "ReducedInteractionMatrix",
    "NGSolveProjectedOperator",
    "NGSolveProjectedInteraction",
    "CurrentMagneticFluxDensitySamples",
    "MagnetizationCurrentCoupling",
    "EVRSTMethodAlgebra",
    "ReducedPortAdmittance",
    "ReducedPortImpedance",
    "SharedMeshMaterialModel",
    "HCurlVIMHDivMMMSolution",
    "CoupledHDivEVRSSystem",
    "CoupledHDivHybridVIMSystem",
    "MixedGalerkinHDivHybridVIMSystem",
    "HCurlVIMHDivMMMSystem",
    "CoupleHDivMagnetizationToEVRS",
    "CoupleHCurlVIMWithHDivMMM",
    "CoupleHybridVIMWithHDivMMM",
    "CoupleEddyBubbleHCurlBasisWithHDivMMM",
    "MixedGalerkinOrthogonalization",
    "SurfaceImpedanceGram",
    "LocalESIMSurfaceLUT",
    "LocalESIMSurfaceModel",
    "LocalESIMSurfaceSolution",
    "BuildLocalESIMSurfaceLUT",
    "ValidateLocalESIMSurfaceLUT",
    "HybridVIMSystem",
    "HCurlEddyCLNModel",
    "HCurlEddyCLNFromVIM",
    "AssembleHybridVIM",
    "AssembleSurfaceImpedanceGram",
    "SolveLocalESIMSurfaceVIM",
    "TopologyAwareHybridVIM",
    "NgsolveTopologyAwareHybridVIM",
    "NgsolveEddyBubbleHybridVIM",
    "NgsolveHCurlVIMHDivMMM",
    "NgsolveBDMEddyBubbleVIM",
    "SkinImpedance",
    "SIBCAdmittanceTail",
    "SIBCSchurTerminationImpedance",
    "SIBCSchurTerminationAdmittance",
    "ExternalVectorPotentialRHS",
]
