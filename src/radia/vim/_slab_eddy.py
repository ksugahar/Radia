"""Readable finite-slab eddy-current validation for the HDiv/HCurl lane.

For a conducting slab in a tangential harmonic magnetic field, the coupled
HDiv magnetization and HCurl eddy-current equations reduce exactly to the
one-dimensional magnetic-diffusion equation.  This module assembles that
reduced equation with P1 Galerkin elements.  The closed form is kept as an
independent validation oracle, not as the numerical implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np


MU0 = 4.0e-7 * math.pi
_HALF_SPACE_SIBC_MIN_THICKNESS_TO_SKIN_DEPTH = 6.0


@dataclass(frozen=True)
class ConductiveSlab:
    """SI-unit material and geometry for a parallel-field conducting slab."""

    thickness_m: float
    relative_permeability: float
    conductivity_s_per_m: float
    area_m2: float = 1.0
    surface_field_a_per_m: complex = 1.0 + 0.0j

    def __post_init__(self) -> None:
        for name in ("thickness_m", "relative_permeability", "conductivity_s_per_m", "area_m2"):
            if not math.isfinite(float(getattr(self, name))) or float(getattr(self, name)) <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
        if not (math.isfinite(complex(self.surface_field_a_per_m).real)
                and math.isfinite(complex(self.surface_field_a_per_m).imag)):
            raise ValueError("surface_field_a_per_m must be finite")
        if abs(complex(self.surface_field_a_per_m)) == 0.0:
            raise ValueError("surface_field_a_per_m must be nonzero")


def _positive_frequency(frequency_hz: float) -> float:
    value = float(frequency_hz)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("frequency_hz must be finite and positive")
    return value


def _skin_depth(slab: ConductiveSlab, frequency_hz: float) -> float:
    omega = 2.0 * math.pi * _positive_frequency(frequency_hz)
    permeability = MU0 * slab.relative_permeability
    return math.sqrt(2.0 / (omega * permeability * slab.conductivity_s_per_m))


def slab_surface_model(
    slab: ConductiveSlab,
    frequency_hz: float,
    *,
    requested: str = "auto",
    min_thickness_to_skin_depth: float = _HALF_SPACE_SIBC_MIN_THICKNESS_TO_SKIN_DEPTH,
) -> dict[str, object]:
    """Select a volumetric model or validate a half-space SIBC request.

    A half-space SIBC is rejected unless the two slab surfaces are separated
    by at least ``min_thickness_to_skin_depth`` skin depths.  Thin conductors
    retain their finite-thickness interaction and must use the volume model.
    """

    if requested not in {"auto", "volumetric", "half_space_sibc"}:
        raise ValueError("requested must be auto, volumetric, or half_space_sibc")
    threshold = float(min_thickness_to_skin_depth)
    if not math.isfinite(threshold) or threshold <= 0.0:
        raise ValueError("min_thickness_to_skin_depth must be finite and positive")
    skin_depth_m = _skin_depth(slab, frequency_hz)
    ratio = slab.thickness_m / skin_depth_m
    sibc_valid = ratio >= threshold
    if requested == "half_space_sibc" and not sibc_valid:
        raise ValueError(
            "half_space_sibc is invalid because thickness/skin_depth "
            f"is {ratio:.6g}, below the required {threshold:.6g}"
        )
    selected = "half_space_sibc" if requested == "auto" and sibc_valid else requested
    if selected == "auto":
        selected = "volumetric"
    return {
        "requested": requested,
        "selected": selected,
        "skin_depth_m": skin_depth_m,
        "thickness_to_skin_depth_ratio": ratio,
        "half_space_sibc_min_ratio": threshold,
        "finite_thickness_separation_verified": sibc_valid,
    }


def parallel_slab_exact_response(slab: ConductiveSlab, frequency_hz: float) -> complex:
    """Return the exact complex effective relative permeability of the slab."""

    delta = _skin_depth(slab, frequency_hz)
    gamma = (1.0 + 1.0j) / delta
    half_argument = gamma * slab.thickness_m / 2.0
    return slab.relative_permeability * np.tanh(half_argument) / half_argument


def solve_parallel_slab_reduced(
    slab: ConductiveSlab,
    frequency_hz: float,
    *,
    elements: int = 128,
    surface_model: str = "volumetric",
) -> dict[str, object]:
    """Solve the exact 1-D HDiv-MMM/HCurl-eddy reduction with P1 FEM."""

    if isinstance(elements, bool) or int(elements) != elements or int(elements) < 2:
        raise ValueError("elements must be an integer >= 2")
    elements = int(elements)
    model = slab_surface_model(slab, frequency_hz, requested=surface_model)
    if model["selected"] != "volumetric":
        raise ValueError("the P1 reduced solve is a volumetric finite-thickness model")

    omega = 2.0 * math.pi * _positive_frequency(frequency_hz)
    permeability = MU0 * slab.relative_permeability
    gamma_squared = 1.0j * omega * permeability * slab.conductivity_s_per_m
    h = slab.thickness_m / elements

    diagonal = np.full(elements + 1, 2.0 / h + gamma_squared * (2.0 * h / 3.0), dtype=complex)
    diagonal[[0, -1]] = 1.0 / h + gamma_squared * (h / 3.0)
    off_diagonal = np.full(elements, -1.0 / h + gamma_squared * (h / 6.0), dtype=complex)
    boundary_value = complex(slab.surface_field_a_per_m)
    rhs = np.zeros(elements - 1, dtype=complex)
    # Accumulate the two boundary contributions separately.  For the minimum
    # two-element mesh they both address the same interior degree of freedom.
    rhs[0] -= off_diagonal[0] * boundary_value
    rhs[-1] -= off_diagonal[-1] * boundary_value
    try:
        from scipy.sparse import diags
        from scipy.sparse.linalg import spsolve

        free_matrix = diags(
            (off_diagonal[1:-1], diagonal[1:-1], off_diagonal[1:-1]),
            (-1, 0, 1),
            shape=(elements - 1, elements - 1),
            format="csr",
        )
        interior = spsolve(free_matrix, rhs)
    except ImportError:  # NumPy is a base dependency; SciPy is an acceleration.
        free_matrix = np.diag(diagonal[1:-1])
        if elements > 2:
            free_matrix += np.diag(off_diagonal[1:-1], 1)
            free_matrix += np.diag(off_diagonal[1:-1], -1)
        interior = np.linalg.solve(free_matrix, rhs)
    field = np.empty(elements + 1, dtype=complex)
    field[[0, -1]] = boundary_value
    field[1:-1] = interior

    residual_vector = free_matrix @ interior - rhs
    residual = float(np.linalg.norm(residual_vector) / max(np.linalg.norm(rhs), np.finfo(float).tiny))
    integral_field = h * (0.5 * field[0] + field[1:-1].sum() + 0.5 * field[-1])
    response = slab.relative_permeability * integral_field / (slab.thickness_m * boundary_value)
    current_density = np.diff(field) / h
    joule_loss_w = float(
        0.5 * slab.area_m2 * h * np.sum(np.abs(current_density) ** 2) / slab.conductivity_s_per_m
    )

    return {
        "lane_id": "hdiv_mmm_hcurl_eddy_bubble",
        "dimensional_reduction": "exact_parallel_slab_hdiv_hcurl_reduction",
        "element_family": "P1",
        "elements": elements,
        "frequency_hz": float(frequency_hz),
        "skin_depth_m": model["skin_depth_m"],
        "elements_per_skin_depth": float(model["skin_depth_m"]) / h,
        "surface_model": model,
        "effective_relative_permeability": complex(response),
        "joule_loss_w": joule_loss_w,
        "normalized_algebraic_residual": residual,
    }


def refine_parallel_slab_reduced(
    slab: ConductiveSlab,
    frequencies_hz: Iterable[float],
    *,
    levels: Iterable[int] = (64, 96, 128),
) -> dict[str, object]:
    """Run a deterministic independent mesh-refinement ledger."""

    mesh_levels = tuple(int(level) for level in levels)
    if len(mesh_levels) < 2 or sorted(set(mesh_levels)) != list(mesh_levels) or mesh_levels[0] < 2:
        raise ValueError("levels must contain at least two strictly increasing element counts")
    frequencies = tuple(_positive_frequency(value) for value in frequencies_hz)
    if not frequencies:
        raise ValueError("frequencies_hz must not be empty")

    runs = []
    for level in mesh_levels:
        rows = [solve_parallel_slab_reduced(slab, frequency, elements=level) for frequency in frequencies]
        runs.append({"elements": level, "rows": rows})
    previous = runs[-2]["rows"]
    final = runs[-1]["rows"]
    changes = [
        abs(complex(current["effective_relative_permeability"])
            - complex(prior["effective_relative_permeability"]))
        / max(abs(complex(current["effective_relative_permeability"])), np.finfo(float).tiny)
        for prior, current in zip(previous, final)
    ]
    return {
        "lane_id": "hdiv_mmm_hcurl_eddy_bubble",
        "dimensional_reduction": "exact_parallel_slab_hdiv_hcurl_reduction",
        "levels": list(mesh_levels),
        "runs": runs,
        "final_rows": final,
        "final_max_relative_change": float(max(changes)),
        "final_relative_changes": [float(value) for value in changes],
    }
