"""H1 mixed total/reduced Omega route for electromagnets.

The formulation is the total/reduced scalar-potential split: a total potential
in the iron, a reduced potential carrying the source field in the air, and an
interface condition between them.  It is due to Simkin and Trowbridge, who
introduced the total scalar potential for this purpose in 1979 [1] and gave
the nonlinear three-dimensional treatment in 1980 [2].  The cancellation the
split exists to avoid -- a large source field minus a nearly equal gradient,
inside high-permeability iron -- is the reason a single reduced potential is
not used everywhere.

    [1] J. Simkin and C. W. Trowbridge, "On the use of the total scalar
        potential in the numerical solution of field problems in
        electromagnetics", Int. J. Numer. Methods Eng. 14(3), 423-440, 1979.
        doi:10.1002/nme.1620140308  (bibliography key ``simkin1979use``)
    [2] J. Simkin and C. W. Trowbridge, "Three-dimensional nonlinear
        electromagnetic field computations, using scalar potentials", IEE
        Proc. B 127(6), 368-374, 1980.  doi:10.1049/ip-b.1980.0052
        (bibliography key ``simkin1980three``)

This module itself is deliberately a small adapter around NGSolve-owned
finite-element spaces and the Kelvin solver.  It fixes the physical partition
and source trace contract shared by every static-electromagnet acceptance
calculation.
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
    reduced_source_load: str = "volume",
    total_source_load: str = "volume",
    source_representation: str = "exact",
    source_interpolation_order: int = 2,
) -> dict[str, object]:
    """Solve one static electromagnet through the required H1 formulation.

    ``nonlinear_method="newton"`` selects quadrature-based PCHIP Newton for
    orders one and two, with residual backtracking. It requires an explicit
    ``nonlinear_material_sampling="integration_point"``; centroid sampling is
    a different discrete material law and is rejected. It does not use a projected
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

    ``reduced_source_load`` / ``total_source_load`` select how the coil
    source enters the load: ``"volume"`` (default) integrates ``H_s`` over the
    region, ``"surface_flux"`` uses the equal boundary normal-flux form
    (``div H_s = 0``) and evaluates the source on faces only.
    ``total_source_load="surface_flux"`` requires the ``total_hodge``
    contract and ``source_trace_tolerance``, which then also gates the iron
    boundary tangential residual: a linked source that the surface form cannot
    represent fails instead of losing its harmonic part.

    With ``bh_table`` the reduced surface load is exact (air stays at mu0).
    The iron term ``mu(H) (H_s + grad Phi_s) . grad v`` cannot move to the
    boundary because ``mu`` varies, so ``total_source_load="surface_flux"``
    drops that harmonic remainder instead: once the gate shows the source is
    exact on the iron, the remainder is zero in the continuous problem and
    only the Hodge projection error of ``Phi_s``.  The result records it as
    ``iron_harmonic_remainder="dropped"``.

    ``source_representation="nodal"`` replaces the coil field in every load
    (reduced load, source traces, iron Hodge projection and iron load) by its
    nodal H1 interpolant of ``source_interpolation_order`` on the region each
    load integrates over: the region's volume for a volume load, its enclosing
    boundary for a surface-flux load.  The reported reduced field keeps the
    exact source.  The interpolant adds an ``O(h^{p+1})`` source error that
    grows where the mesh does not resolve the coil's proximity, so it must be
    validated per mesh.  Linear only.
    """
    from radia.kelvin_solver import SOURCE_LOADS

    for name, value in (("reduced_source_load", reduced_source_load),
                        ("total_source_load", total_source_load)):
        if value not in SOURCE_LOADS:
            raise ValueError(f"{name} must be one of {SOURCE_LOADS}; got {value!r}")
    if total_source_load == "surface_flux" and source_potential_contract != "total_hodge":
        raise ValueError("total_source_load='surface_flux' requires source_potential_contract='total_hodge'")
    if source_representation not in ("exact", "nodal"):
        raise ValueError("source_representation must be 'exact' or 'nodal'")
    if source_representation == "nodal" and bh_table is not None:
        raise ValueError("source_representation='nodal' is implemented for the linear solve only")
    if source_representation == "nodal" and source_potential_contract == "global_physical":
        raise ValueError("source_representation='nodal' does not cover the global_physical contract")
    if int(order) < 1:
        raise ValueError("order must be positive")
    if (bh_table is not None and nonlinear_method == "newton"
            and nonlinear_material_sampling != "integration_point"):
        raise ValueError(
            "Newton requires nonlinear_material_sampling='integration_point'; "
            "element_centroid is a different material discretization")
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

    reduced_load_h = None
    total_load_h = source_h
    interpolation = {}
    if source_representation == "nodal":
        from radia.kelvin_solver import interpolate_source_field, material_set_boundary_signs

        def load_region(materials, load):
            if load == "volume":
                return mesh.Materials("|".join(materials))
            return mesh.Boundaries("|".join(sorted(
                material_set_boundary_signs(mesh, materials))))

        reduced_interpolant = interpolate_source_field(
            mesh, source_h, load_region(tuple(domain.reduced_materials), reduced_source_load),
            order=int(source_interpolation_order))
        reduced_load_h = reduced_interpolant["field"]
        interpolation["reduced"] = {key: reduced_interpolant[key]
                                    for key in ("nodes", "order", "timings_seconds")}
        if source_potential_contract == "total_hodge":
            iron_materials = tuple(name for name in domain.total_materials
                                   if name not in domain.kelvin_materials)
            total_interpolant = interpolate_source_field(
                mesh, source_h, load_region(iron_materials, total_source_load),
                order=int(source_interpolation_order))
            total_load_h = total_interpolant["field"]
            interpolation["total"] = {key: total_interpolant[key]
                                      for key in ("nodes", "order", "timings_seconds")}
    trace_source_h = source_h if reduced_load_h is None else reduced_load_h

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
            trace_source_h,
            domain.reduced_total_interface,
            order=int(source_projection_order),
            relative_tolerance=source_trace_tolerance,
        )
        kelvin_trace = project_source_interface_potential(
            mesh,
            trace_source_h,
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
        total_source_potential = None
        total_source_materials = ()
    elif source_potential_contract == "total_hodge":
        total_source_materials = tuple(
            name for name in domain.total_materials
            if name not in domain.kelvin_materials
        )
        source_hodge = project_source_total_hodge(
            mesh,
            total_load_h,
            total_source_materials,
            order=int(source_projection_order),
            bonus_intorder=bonus_intorder,
            source_load=total_source_load,
            tangential_tolerance=(source_trace_tolerance
                                  if total_source_load == "surface_flux" else None),
        )
        # The exact pulled-back exterior source needs no interface trace at
        # all, so the projection is skipped rather than computed and dropped.
        kelvin_trace = (
            None if kelvin_source_h is not None
            else project_source_interface_potential(
                mesh,
                trace_source_h,
                domain.kelvin_interface,
                order=int(source_projection_order),
                relative_tolerance=source_trace_tolerance,
            )
        )
        source_potential = source_hodge["potential"]
        kelvin_source_potential = (
            None if kelvin_trace is None else kelvin_trace["potential"])
        if total_source_load == "surface_flux" and bh_table is not None:
            # mu(H) varies in the iron: keep the gated exact part only.
            total_source_h = None
            total_source_potential = None
            total_source_materials = ()
        elif total_source_load == "surface_flux":
            total_source_h = total_load_h
            total_source_potential = source_hodge["potential"]
        else:
            total_source_h = source_hodge["harmonic_field"]
            total_source_potential = None
        source_diagnostics = {
            "contract": source_potential_contract,
            "projection_order": int(source_projection_order),
            "total_source_materials": list(total_source_materials),
            "iron_relative_harmonic_norm": (
                None if source_hodge["relative_harmonic_norm"] is None
                else float(source_hodge["relative_harmonic_norm"])),
            "iron_relative_tangential_residual": source_hodge.get(
                "relative_tangential_residual"),
            "reduced_source_load": reduced_source_load,
            "total_source_load": total_source_load,
            "iron_harmonic_remainder": (
                "dropped" if total_source_load == "surface_flux" and bh_table is not None
                else "boundary_identity" if total_source_load == "surface_flux"
                else "volume"),
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
        total_source_potential = None
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
            reduced_source_load=reduced_source_load,
            total_source_load=total_source_load,
            total_source_potential=total_source_potential,
            load_source_h=reduced_load_h,
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
        if reduced_source_load != "volume":
            # Picard takes it explicitly; Newton forwards it to its linear solve.
            iteration_options["reduced_source_load"] = reduced_source_load
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
    source_diagnostics["source_representation"] = source_representation
    if interpolation:
        source_diagnostics["interpolation"] = interpolation
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
