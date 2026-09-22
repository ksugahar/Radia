"""H1 mixed total/reduced scalar-potential route for electromagnets.

This is deliberately a small adapter around NGSolve-owned finite-element
spaces and the Kelvin solver.  It fixes the physical partition and source
trace contract shared by every static-electromagnet acceptance calculation.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


MIXED_TOTAL_REDUCED_OMEGA = "mixed_total_reduced_omega"
MIXED_TOTAL_REDUCED_OMEGA_LABEL = "H1 mixed total/reduced Omega"


def _unique_names(names: Iterable[str], *, field: str) -> tuple[str, ...]:
    values = tuple(str(name) for name in names)
    if not values:
        raise ValueError(f"{field} must name at least one mesh material")
    if len(values) != len(set(values)):
        raise ValueError(f"{field} must not contain duplicate material names")
    if any(not name for name in values):
        raise ValueError(f"{field} must not contain an empty material name")
    return values


@dataclass(frozen=True)
class StaticElectromagnetMixedDomain:
    """Explicit source/reduced and total-potential material partition.

    The physical air holding the CoilBuilder source is the reduced region;
    iron and Kelvin exterior are total-potential regions.  The Kelvin gauge is
    a ``BBND`` point constraint, never an ordinary surface label.
    """

    reduced_materials: tuple[str, ...]
    total_materials: tuple[str, ...]
    nonlinear_materials: tuple[str, ...]
    reduced_total_interface: str = "iron_air_interface"
    kelvin_interface: str = "kelvin_int"
    ground_boundary: str = "GND"
    kelvin_materials: tuple[str, ...] = ("kelvin",)

    def __post_init__(self) -> None:
        reduced = _unique_names(self.reduced_materials, field="reduced_materials")
        total = _unique_names(self.total_materials, field="total_materials")
        nonlinear = tuple(str(name) for name in self.nonlinear_materials)
        kelvin = _unique_names(self.kelvin_materials, field="kelvin_materials")
        if len(nonlinear) != len(set(nonlinear)):
            raise ValueError("nonlinear_materials must not contain duplicates")
        if not set(nonlinear) <= set(total):
            raise ValueError("nonlinear_materials must be contained in total_materials")
        if set(reduced) & set(total):
            raise ValueError("reduced_materials and total_materials must be disjoint")
        if not set(kelvin) <= set(total):
            raise ValueError("kelvin_materials must be contained in total_materials")
        for name in (
            "reduced_total_interface",
            "kelvin_interface",
            "ground_boundary",
        ):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"{name} must be a non-empty boundary name")

    def validate_mesh_labels(
        self,
        materials: Iterable[str],
        boundaries: Iterable[str],
        bbboundaries: Iterable[str] = (),
    ) -> None:
        """Reject a mesh that cannot represent this physical split."""
        material_set = {str(name) for name in materials}
        expected = set(self.reduced_materials) | set(self.total_materials)
        if material_set != expected:
            raise ValueError(
                "mixed total/reduced Omega requires an exhaustive declared "
                f"material partition; undeclared_mesh_materials="
                f"{sorted(material_set - expected)}, "
                f"declared_but_absent_materials={sorted(expected - material_set)}"
            )
        boundary_set = {str(name) for name in boundaries}
        missing_boundaries = sorted(
            {self.reduced_total_interface, self.kelvin_interface} - boundary_set
        )
        missing_bbboundaries = sorted(
            {self.ground_boundary} - {str(name) for name in bbboundaries}
        )
        if missing_boundaries or missing_bbboundaries:
            raise ValueError(
                "mixed total/reduced Omega requires explicit source/total and "
                "Kelvin trace BNDs plus its point GND constraint; "
                f"missing_boundaries={missing_boundaries}, "
                f"missing_bbboundaries={missing_bbboundaries}"
            )

    def as_dict(self) -> dict[str, object]:
        return {
            "formulation": MIXED_TOTAL_REDUCED_OMEGA_LABEL,
            "reduced_materials": list(self.reduced_materials),
            "total_materials": list(self.total_materials),
            "nonlinear_materials": list(self.nonlinear_materials),
            "reduced_total_interface": self.reduced_total_interface,
            "kelvin_interface": self.kelvin_interface,
            "ground_boundary": self.ground_boundary,
            "kelvin_materials": list(self.kelvin_materials),
        }


def solve_static_electromagnet_mixed_total_reduced_omega(
    mesh,
    source_h,
    domain: StaticElectromagnetMixedDomain,
    kelvin_radius: float,
    kelvin_offset,
    *,
    order: int,
    linear_mu_r_by_material: dict[str, float] | None = None,
    bh_table=None,
    source_trace_tolerance: float | None = None,
    source_potential_contract: str = "surface_trace",
    source_projection_order: int | None = None,
    nonlinear_tolerance: float = 2.0e-5,
    nonlinear_max_iterations: int = 80,
    nonlinear_relaxation: float = 0.3,
    kelvin_source_h=None,
    nonlinear_anderson_depth: int = 0,
    nonlinear_anderson_transform: str = "log",
    nonlinear_mu_r_initial=1000.0,
    nonlinear_observation_points=None,
    nonlinear_material_update_order: int | None = None,
    nonlinear_material_log_state_initial=None,
    nonlinear_material_sampling: str = "element_centroid",
    nonlinear_bh_interpolation: str = "pchip",
    nonlinear_method: str = "picard",
    nonlinear_residual_tolerance: float = 1e-8,
    nonlinear_progress_callback=None,
    inverse: str = "pardiso",
    bonus_intorder: int = 4,
) -> dict[str, object]:
    """Solve one static electromagnet through the required H1 formulation.

    ``nonlinear_method="newton"`` selects quadrature-based PCHIP Newton for
    orders one and two, with residual backtracking. It does not use a projected
    material state or Anderson mixing. ``nonlinear_residual_tolerance`` bounds
    its free-DOF equation residual; ``nonlinear_tolerance`` bounds field change.

    The nonlinear loop is the Picard iteration of
    :func:`radia.kelvin_solver.solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin`:
    ``nonlinear_anderson_depth`` enables its constrained Anderson mixing,
    ``nonlinear_mu_r_initial`` is a scalar or the per-element warm start of an
    earlier order-one ``nonlinear_stats["mu_r_elements"]``, and
    ``nonlinear_observation_points`` records the per-iteration field change at
    the points where the result is consumed.  A non-converged loop raises
    :class:`radia.kelvin_solver.MixedOmegaPicardNotConverged` with that state.
    For a P1 diagnostic, ``nonlinear_material_sampling="integration_point"``
    evaluates the B(H) secant at volume quadrature points. It requires plain
    Picard (relaxation=1, Anderson depth=0), and reports a separate returned-field
    constitutive defect rather than treating iterate convergence as accuracy.
    ``nonlinear_bh_interpolation="linear_spline"`` uses NGSolve's compact
    piecewise-linear table lookup; its interpolation error must be checked
    against the original B-H curve for a validation comparison.
    Response order two requires an explicit
    ``nonlinear_material_update_order=1``; its positive log-permeability field
    is a separate spatial material state. Resume it with the complete
    ``nonlinear_stats["material_restart_state"]`` mapping so mesh, B-H table,
    material selector, orders, and active DOFs are checked before solving.

    ``source_potential_contract="total_hodge"`` is the general CoilBuilder
    route.  It retains the non-exact harmonic/cut component of a linked source
    inside total-potential iron instead of forcing it into a scalar trace.
    ``"surface_trace"`` is the strict scalar-only contract for simply connected
    interfaces.  Fixed permanent magnetization may use ``"global_physical"``.
    The total-Hodge projection defaults to the response order: a higher-order
    lift is not generally representable in a lower-order total space, even in
    vacuum. Explicit over-order projections remain available for convergence
    studies and carry a discretization warning. Other trace contracts retain
    their at-least-order-two projection default.
    """
    if int(order) < 1:
        raise ValueError("order must be positive")
    if source_projection_order is None:
        source_projection_order = (
            int(order) if source_potential_contract == "total_hodge"
            else max(2, int(order)))
    if int(source_projection_order) < 1:
        raise ValueError("source_projection_order must be positive")
    if float(kelvin_radius) <= 0.0:
        raise ValueError("kelvin_radius must be positive")
    if linear_mu_r_by_material is None and bh_table is None:
        raise ValueError(
            "supply linear_mu_r_by_material or bh_table (or both for hybrid materials) for "
            "mixed total/reduced Omega"
        )
    if source_potential_contract not in {
            "surface_trace", "total_hodge", "global_physical"}:
        raise ValueError(
            "source_potential_contract must be 'surface_trace', "
            "'total_hodge', or 'global_physical'"
        )
    domain.validate_mesh_labels(
        mesh.GetMaterials(), mesh.GetBoundaries(), mesh.GetBBBoundaries()
    )

    from radia.kelvin_solver import (
        project_source_total_hodge,
        project_source_physical_potential,
        project_source_interface_potential,
        solve_magnetostatic_mixed_total_reduced_omega_kelvin,
        solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin,
    )

    if source_potential_contract == "surface_trace":
        source_trace = project_source_interface_potential(
            mesh,
            source_h,
            domain.reduced_total_interface,
            order=int(source_projection_order),
            relative_tolerance=source_trace_tolerance,
        )
        kelvin_trace = project_source_interface_potential(
            mesh,
            source_h,
            domain.kelvin_interface,
            order=int(source_projection_order),
            relative_tolerance=source_trace_tolerance,
        )
        source_potential = source_trace["potential"]
        kelvin_source_potential = kelvin_trace["potential"]
        source_diagnostics = {
            "contract": source_potential_contract,
            "projection_order": int(source_projection_order),
            "iron_air_relative_tangential_residual": float(
                source_trace["relative_tangential_residual"]
            ),
            "kelvin_relative_tangential_residual": float(
                kelvin_trace["relative_tangential_residual"]
            ),
            "relative_tolerance": source_trace_tolerance,
        }
        total_source_h = None
        total_source_materials = ()
    elif source_potential_contract == "total_hodge":
        total_source_materials = tuple(
            name for name in domain.total_materials
            if name not in domain.kelvin_materials
        )
        source_hodge = project_source_total_hodge(
            mesh,
            source_h,
            total_source_materials,
            order=int(source_projection_order),
            bonus_intorder=bonus_intorder,
        )
        # The exact pulled-back exterior source needs no interface trace at
        # all, so the projection is skipped rather than computed and dropped.
        kelvin_trace = (
            None if kelvin_source_h is not None
            else project_source_interface_potential(
                mesh,
                source_h,
                domain.kelvin_interface,
                order=int(source_projection_order),
                relative_tolerance=source_trace_tolerance,
            )
        )
        source_potential = source_hodge["potential"]
        kelvin_source_potential = (
            None if kelvin_trace is None else kelvin_trace["potential"])
        total_source_h = source_hodge["harmonic_field"]
        source_diagnostics = {
            "contract": source_potential_contract,
            "projection_order": int(source_projection_order),
            "total_source_materials": list(total_source_materials),
            "iron_relative_harmonic_norm": float(
                source_hodge["relative_harmonic_norm"]
            ),
            "projection_bonus_intorder": source_hodge["bonus_intorder"],
            "kelvin_exterior_source": (
                "exact pulled-back field" if kelvin_trace is None
                else "projected interface trace"),
            "relative_tolerance": source_trace_tolerance,
        }
        if kelvin_trace is not None:
            source_diagnostics["kelvin_relative_tangential_residual"] = float(
                kelvin_trace["relative_tangential_residual"])
    else:
        physical_materials = tuple(
            name for name in domain.reduced_materials + domain.total_materials
            if name not in domain.kelvin_materials
        )
        source_volume = project_source_physical_potential(
            mesh,
            source_h,
            physical_materials,
            order=int(source_projection_order),
            relative_tolerance=source_trace_tolerance,
        )
        source_potential = source_volume["potential"]
        kelvin_source_potential = source_potential
        total_source_h = None
        total_source_materials = ()
        source_diagnostics = {
            "contract": source_potential_contract,
            "projection_order": int(source_projection_order),
            "physical_materials": list(physical_materials),
            "relative_volume_residual": float(
                source_volume["relative_volume_residual"]
            ),
            "relative_tolerance": source_trace_tolerance,
        }
    common = {
        "reduced_materials": domain.reduced_materials,
        "total_materials": domain.total_materials,
        "interface_boundary": domain.reduced_total_interface,
        "order": int(order),
        "dirichlet_bbbnd": domain.ground_boundary,
        "bonus_intorder": int(bonus_intorder),
        "inverse": inverse,
        "kelvin_mats": domain.kelvin_materials,
        "kelvin_match_exact": True,
        "kelvin_interface_boundary": domain.kelvin_interface,
        "kelvin_source_potential": (None if kelvin_source_h is not None
                                    else kelvin_source_potential),
        "kelvin_source_h": kelvin_source_h,
        "total_source_h": total_source_h,
        "total_source_materials": total_source_materials,
    }
    if bh_table is None:
        result = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh,
            source_h,
            source_potential,
            float(kelvin_radius),
            kelvin_offset,
            mu_r_by_material=dict(linear_mu_r_by_material),
            **common,
        )
    else:
        if not domain.nonlinear_materials:
            raise ValueError("bh_table requires declared nonlinear_materials")
        if nonlinear_method not in ("picard", "newton"):
            raise ValueError("nonlinear_method must be 'picard' or 'newton'")
        nonlinear_solver = solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin
        iteration_options = dict(
            relaxation=float(nonlinear_relaxation),
            anderson_depth=int(nonlinear_anderson_depth),
            anderson_transform=str(nonlinear_anderson_transform),
            material_update_order=nonlinear_material_update_order,
            material_log_state_initial=nonlinear_material_log_state_initial,
            material_sampling=nonlinear_material_sampling,
            bh_interpolation=nonlinear_bh_interpolation)
        if nonlinear_method == "newton":
            from .mixed_omega_newton import solve_magnetostatic_mixed_total_reduced_omega_newton_kelvin
            if (nonlinear_anderson_depth != 0
                    or nonlinear_material_update_order is not None
                    or nonlinear_material_log_state_initial is not None
                    or nonlinear_bh_interpolation != "pchip"):
                raise ValueError("Newton requires PCHIP without Anderson or projected material state")
            nonlinear_solver = solve_magnetostatic_mixed_total_reduced_omega_newton_kelvin
            iteration_options = dict(residual_tolerance=float(nonlinear_residual_tolerance))
        result = nonlinear_solver(
            mesh,
            source_h,
            source_potential,
            float(kelvin_radius),
            kelvin_offset,
            bh_table=bh_table,
            mu_r_by_material=linear_mu_r_by_material,
            nonlinear_materials=domain.nonlinear_materials,
            tolerance=float(nonlinear_tolerance),
            max_iterations=int(nonlinear_max_iterations),
            mu_r_initial=nonlinear_mu_r_initial,
            observation_points=nonlinear_observation_points,
            progress_callback=nonlinear_progress_callback,
            **iteration_options,
            **common,
        )
    trace_gate_applied = source_trace_tolerance is not None and (
        source_potential_contract != "total_hodge" or kelvin_trace is not None)
    source_diagnostics["gate_enabled"] = trace_gate_applied
    source_diagnostics["acceptance"] = "passed" if trace_gate_applied else "not_evaluated"
    source_diagnostics["gate_scope"] = (
        "kelvin_trace_only" if source_potential_contract == "total_hodge"
        else source_potential_contract)
    if source_potential_contract == "total_hodge":
        source_diagnostics["response_order"] = int(order)
        source_diagnostics["lift_order_compatible"] = int(source_projection_order) <= int(order)
        source_diagnostics["split_invariance_verified"] = False
        if int(source_projection_order) > int(order):
            source_diagnostics["discretization_warning"] = (
                "source lift exceeds the response space; a smaller source projection "
                "residual does not imply a more accurate physical field. Run a vacuum "
                "split-invariance test or increase response order.")
    result["static_electromagnet_contract"] = domain.as_dict()
    result["static_electromagnet_contract"]["source_trace"] = source_diagnostics
    return result
