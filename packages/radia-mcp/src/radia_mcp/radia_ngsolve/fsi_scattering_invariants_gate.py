"""Conservation and reciprocity gate for lossless acoustic FSI scattering."""
from __future__ import annotations

import math


def fsi_scattering_invariants_gate(
    *,
    reciprocity_relative_error: float,
    optical_theorem_relative_error: float,
    bem_dtn_relative_error: float,
    max_solver_residual: float,
    lossless_material: bool,
    time_convention: str,
    max_invariant_error: float = 0.05,
    max_bem_dtn_error: float = 0.05,
    max_solver_residual_allowed: float = 1.0e-8,
) -> dict:
    """Require reciprocal, energy-conserving P1 FEM/BEM scattering evidence."""

    scalars = (
        reciprocity_relative_error,
        optical_theorem_relative_error,
        bem_dtn_relative_error,
        max_solver_residual,
        max_invariant_error,
        max_bem_dtn_error,
        max_solver_residual_allowed,
    )
    if not all(math.isfinite(float(value)) for value in scalars):
        raise ValueError("all errors, residuals, and tolerances must be finite")
    if any(float(value) < 0.0 for value in scalars):
        raise ValueError("all errors, residuals, and tolerances must be nonnegative")

    convention = str(time_convention).strip().lower().replace(" ", "").replace("*", "")
    recognized_convention = convention in {"exp(-iomegat)", "exp(-iwt)", "e^-iwt"}
    checks = {
        "lossless_material_declared": lossless_material is True,
        "outgoing_time_convention_explicit": recognized_convention,
        "farfield_reciprocity": reciprocity_relative_error <= max_invariant_error,
        "optical_theorem_energy_closure": optical_theorem_relative_error <= max_invariant_error,
        "p1_bem_high_order_dtn_agreement": bem_dtn_relative_error <= max_bem_dtn_error,
        "coupled_linear_system_converged": max_solver_residual <= max_solver_residual_allowed,
    }
    return {
        "policy": "fsi_scattering_invariants_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "reciprocity_relative_error": float(reciprocity_relative_error),
            "optical_theorem_relative_error": float(optical_theorem_relative_error),
            "bem_dtn_relative_error": float(bem_dtn_relative_error),
            "max_solver_residual": float(max_solver_residual),
        },
        "tolerances": {
            "max_invariant_error": float(max_invariant_error),
            "max_bem_dtn_error": float(max_bem_dtn_error),
            "max_solver_residual_allowed": float(max_solver_residual_allowed),
        },
        "conventions": {
            "time_dependence": str(time_convention),
            "reciprocity_transform": "f(obs=b,inc=a)=f(obs=-a,inc=-b)",
            "optical_theorem": "sigma_scat=(4*pi/k)*Im(f_forward)",
        },
        "notes": [
            "The optical-theorem check is valid only for a lossless coupled scatterer.",
            "Dense P1 BEM and high-order DtN are independent exterior closures on the same P1 interface.",
        ],
    }
