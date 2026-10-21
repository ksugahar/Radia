"""CPU thermal FE operators for IH, independent of MATLAB and Simulink.

The caller owns ngsolve.TaskManager(). NGSolve owns spaces, assembly and
mapped evaluation. Integration adapters consume these operators without
copying the numerical implementation.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol
import math
import numpy as np


def nonlinear_thermal_quadrature(*, basis, gradient, volume_weights, density,
                                boundary_basis, boundary_weights, convection,
                                constant_coefficients, kelvin, cp, kappa,
                                emissivity):
    """Validate NGSolve-owned sample maps for the explicit nonlinear IH route.

    The caller must evaluate its own space at the declared mapped quadrature
    points. This function does not infer element numbering, mass lumping,
    gradients, rotation or a nodal representation from an exported matrix.
    Returned ``native_config`` is a fragment, not a complete EM configuration.
    """
    from .ih_thermal_coupling import ThermalLaw, ThermalQuadrature
    return ThermalQuadrature(basis, gradient, volume_weights, density,
        boundary_basis, boundary_weights, convection, constant_coefficients,
        ThermalLaw(kelvin, cp, kappa, emissivity))


class ThermalAssemblyOptions(Protocol):
    axisymmetric_thermal_vol: str | None
    density_kg_per_m3: float
    heat_capacity_J_per_kgK: float
    thermal_conductivity_W_per_mK: float
    thermal_order: int
    workpiece_label: str


@dataclass(frozen=True, slots=True)
class ThermalOperators:
    n_temperature: int
    heat_dofs: np.ndarray
    unit_heat_density_W_per_m3: np.ndarray
    heat_cell_weights_m3: np.ndarray
    temperature_cell_weights_J_per_K: np.ndarray
    heat_to_temperature_projection: list[float]
    mass_row_ptr: list[int]
    mass_col: list[int]
    mass_value: list[float]
    stiffness_row_ptr: list[int]
    stiffness_col: list[int]
    stiffness_value: list[float]
    convection_row_ptr: list[int]
    convection_col: list[int]
    convection_value: list[float]
    heat_power_W: float
    full_heat_density_W_per_m3: np.ndarray
    constant_coefficients: np.ndarray | None = None
    temperature_evaluation: dict[str, Any] | None = None


def _coo_values(matrix: Any) -> dict[tuple[int, int], float]:
    rows, columns, values = matrix.COO()
    result: dict[tuple[int, int], float] = {}
    for row, column, value in zip(rows, columns, values):
        key = (int(row), int(column))
        result[key] = result.get(key, 0.0) + float(value)
    return result


def _aligned_csr(
    n: int, *matrices: dict[tuple[int, int], float]
) -> tuple[list[int], list[int], list[list[float]]]:
    keys = set().union(*(matrix.keys() for matrix in matrices))
    keys.update((index, index) for index in range(n))
    by_row: list[list[int]] = [[] for _ in range(n)]
    for row, column in keys:
        if not (0 <= row < n and 0 <= column < n):
            raise ValueError("assembled matrix contains an out-of-range entry")
        by_row[row].append(column)
    row_ptr = [0]
    columns: list[int] = []
    for row_columns in by_row:
        columns.extend(sorted(set(row_columns)))
        row_ptr.append(len(columns))
    aligned: list[list[float]] = []
    for matrix in matrices:
        aligned.append(
            [
                float(matrix.get((row, column), 0.0))
                for row, row_columns in enumerate(by_row)
                for column in sorted(set(row_columns))
            ]
        )
    return row_ptr, columns, aligned


def _row_sums(n: int, matrix: dict[tuple[int, int], float]) -> np.ndarray:
    result = np.zeros(n, dtype=float)
    for (row, _), value in matrix.items():
        result[row] += value
    return result


def _matvec(n: int, matrix: dict[tuple[int, int], float], values: np.ndarray) -> np.ndarray:
    result = np.zeros(n, dtype=float)
    for (row, column), value in matrix.items():
        result[row] += value * values[column]
    return result


def _assemble_thermal_operators(
    workpiece: Path,
    heat_flux_W_per_m2: np.ndarray,
    options: ThermalAssemblyOptions,
) -> ThermalOperators:
    from ngsolve import CF, H1, BilinearForm, Mesh, ds, dx, grad, x

    mesh = Mesh(str(workpiece))
    expected_dim = 2 if options.axisymmetric_thermal_vol else 3
    if mesh.dim != expected_dim or mesh.ne <= 0:
        raise ValueError(f"the IH thermal workpiece must be a {expected_dim}D mesh")
    weight = 2 * math.pi * x if expected_dim == 2 else CF(1)
    convection_region = options.workpiece_label if expected_dim == 2 else ".*"
    if expected_dim == 2:
        from radia.panels.calc_heat_axisym import _reject_axis_boundary_roles
        if any(float(vertex.point[0]) < 0 for vertex in mesh.vertices):
            raise ValueError("axisymmetric thermal radii must be nonnegative")
        _reject_axis_boundary_roles(mesh, [("native heat/convection", [options.workpiece_label])])
    materials = sorted(set(str(name) for name in mesh.GetMaterials()))
    if materials != ["workpiece"]:
        raise ValueError(
            "the IH thermal mesh must contain exactly one material named "
            f"'workpiece'; got {materials}"
        )
    if options.workpiece_label not in set(str(name) for name in mesh.GetBoundaries()):
        raise ValueError(f"workpiece boundary {options.workpiece_label!r} is absent from the mesh")

    if options.thermal_order == 2:
        return _assemble_p2_thermal(mesh, heat_flux_W_per_m2, options, weight, convection_region)
    fes = H1(mesh, order=options.thermal_order)
    if heat_flux_W_per_m2.size != fes.ndof:
        raise ValueError("the unit-current heat field and thermal H1 space have different sizes")
    u, v = fes.TnT()
    rho_cp = options.density_kg_per_m3 * options.heat_capacity_J_per_kgK
    mass = BilinearForm(fes, symmetric=True)
    mass += weight * CF(rho_cp) * u * v * dx
    mass.Assemble()
    stiffness = BilinearForm(fes, symmetric=True)
    stiffness += weight * CF(options.thermal_conductivity_W_per_mK) * grad(u) * grad(v) * dx
    stiffness.Assemble()
    heat_mass = BilinearForm(fes, symmetric=True, check_unused=False)
    heat_mass += weight * u * v * ds(options.workpiece_label)
    heat_mass.Assemble()
    convection = BilinearForm(fes, symmetric=True, check_unused=False)
    convection += weight * u * v * ds(convection_region)
    convection.Assemble()

    mass_values = _coo_values(mass.mat)
    stiffness_values = _coo_values(stiffness.mat)
    heat_mass_values = _coo_values(heat_mass.mat)
    convection_values = _coo_values(convection.mat)
    row_ptr, columns, aligned = _aligned_csr(
        fes.ndof, mass_values, stiffness_values, convection_values
    )
    mass_csr, stiffness_csr, convection_csr = aligned

    capacity_weights = _row_sums(fes.ndof, mass_values)
    volume_weights = capacity_weights / rho_cp
    surface_weights = _row_sums(fes.ndof, heat_mass_values)
    scale = max(1.0, float(np.max(np.abs(surface_weights))))
    heat_dofs = np.flatnonzero(surface_weights > 1.0e-14 * scale)
    if heat_dofs.size == 0:
        raise RuntimeError("the workpiece SIBC boundary has no active H1 heat DOFs")
    if np.any(capacity_weights <= 0.0) or np.any(volume_weights <= 0.0):
        raise RuntimeError("the assembled H1 mass matrix has non-positive lumped weights")

    heat_load = _matvec(fes.ndof, heat_mass_values, heat_flux_W_per_m2)
    outside = np.ones(fes.ndof, dtype=bool)
    outside[heat_dofs] = False
    if np.any(outside) and np.linalg.norm(heat_load[outside], ord=np.inf) > 1.0e-11 * max(
        1.0, np.linalg.norm(heat_load, ord=np.inf)
    ):
        raise RuntimeError("the SIBC boundary load leaked into non-boundary heat DOFs")
    unit_density = heat_load[heat_dofs] / volume_weights[heat_dofs]
    if not np.all(np.isfinite(unit_density)) or np.min(unit_density) < -1.0e-10:
        raise RuntimeError("the assembled unit-current heat density is invalid")
    unit_density = np.maximum(unit_density, 0.0)
    heat_power = float(np.dot(volume_weights[heat_dofs], unit_density))
    if not math.isfinite(heat_power) or heat_power <= 0.0:
        raise RuntimeError("the unit-current electromagnetic solve produced no heat")

    n_heat = int(heat_dofs.size)
    projection = [0.0] * (fes.ndof * n_heat)
    for heat_index, dof in enumerate(heat_dofs):
        projection[int(dof) * n_heat + heat_index] = float(volume_weights[dof])
    full_density = np.zeros(fes.ndof, dtype=float)
    full_density[heat_dofs] = unit_density
    return ThermalOperators(
        n_temperature=fes.ndof,
        heat_dofs=heat_dofs,
        unit_heat_density_W_per_m3=unit_density,
        heat_cell_weights_m3=volume_weights[heat_dofs],
        temperature_cell_weights_J_per_K=capacity_weights,
        heat_to_temperature_projection=projection,
        mass_row_ptr=row_ptr,
        mass_col=columns,
        mass_value=mass_csr,
        stiffness_row_ptr=row_ptr,
        stiffness_col=columns,
        stiffness_value=stiffness_csr,
        convection_row_ptr=row_ptr,
        convection_col=columns,
        convection_value=convection_csr,
        heat_power_W=heat_power,
        full_heat_density_W_per_m3=full_density,
    )


def _assemble_p2_thermal(mesh, flux, options, weight, convection_region):
    """Keep the P1 surface source and P2 temperature spaces independent.

    Native temperatures are NGSolve coefficients, not nodal kelvin values.
    NGSolve owns the mixed load, constant function and mapped evaluation.
    """
    import ngsolve as ng

    heat_fes, fes = ng.H1(mesh, order=1), ng.H1(mesh, order=2)
    if flux.size != heat_fes.ndof or not np.all(np.isfinite(flux)) or np.min(flux) < 0:
        raise ValueError("P2 thermal assembly requires a finite nonnegative P1 surface field")
    u, v = fes.TnT()
    rho_cp = options.density_kg_per_m3 * options.heat_capacity_J_per_kgK
    forms = []
    for integrand in (weight*rho_cp*u*v*ng.dx,
                      weight*options.thermal_conductivity_W_per_mK*ng.grad(u)*ng.grad(v)*ng.dx,
                      weight*u*v*ng.ds(convection_region)):
        form = ng.BilinearForm(fes, symmetric=True, check_unused=False)
        form += integrand
        form.Assemble()
        forms.append(_coo_values(form.mat))
    constant = ng.GridFunction(fes)
    constant.Set(ng.CF(1))
    c = constant.vec.FV().NumPy().copy()
    capacity = _matvec(fes.ndof, forms[0], c)
    mixed = ng.BilinearForm(trialspace=heat_fes, testspace=fes)
    mixed += weight*heat_fes.TrialFunction()*fes.TestFunction()*ng.ds(options.workpiece_label)
    mixed.Assemble()
    mixed_values = _coo_values(mixed.mat)
    area = ng.LinearForm(heat_fes)
    area += weight*heat_fes.TestFunction()*ng.ds(options.workpiece_label)
    area.Assemble()
    surface_weights = area.vec.FV().NumPy().copy()
    heat_dofs = np.flatnonzero(surface_weights > 1e-14*max(1., np.max(surface_weights)))
    volume = float(np.dot(c, capacity)/rho_cp)
    length = volume / float(np.sum(surface_weights))
    if not math.isfinite(length) or length <= 0 or not heat_dofs.size:
        raise ValueError("P2 thermal geometry has no positive heated area or volume")
    projection = np.zeros((fes.ndof, len(heat_dofs)))
    indices = {int(d): i for i, d in enumerate(heat_dofs)}
    for (row, col), value in mixed_values.items():
        if col in indices:
            projection[row, indices[col]] += length*value
    density = flux[heat_dofs]/length
    power = float(np.dot(surface_weights[heat_dofs], flux[heat_dofs]))
    if not math.isfinite(power) or power <= 0:
        raise ValueError("P2 thermal source must have positive power")
    # Sample the actual FE function at mapped quadrature points, never min/max
    # of modal coefficients. Retain a sparse evaluation operator for MATLAB.
    rules = {el.type: ng.IntegrationRule(el.type, 4) for el in mesh.Elements(ng.VOL)}
    points = mesh.MapToAllElements(rules, ng.VOL)
    rows, cols, values = [], [], []
    probe = ng.GridFunction(fes)
    for col in range(fes.ndof):
        probe.vec[:] = 0
        probe.vec[col] = 1
        evaluated = np.asarray(probe(points)).reshape(-1)
        active = np.flatnonzero(np.abs(evaluated) > 1e-14)
        rows.extend(active.tolist()); cols.extend([col]*len(active)); values.extend(evaluated[active].tolist())
    evaluation = dict(n_samples=len(points), rows=rows, cols=cols, values=values,
                      sampling="mapped-volume-quadrature-order-4")
    row_ptr, columns, aligned = _aligned_csr(fes.ndof, *forms)
    return ThermalOperators(fes.ndof, heat_dofs, density, surface_weights[heat_dofs]*length,
        capacity, projection.ravel().tolist(), row_ptr, columns, aligned[0],
        row_ptr, columns, aligned[1], row_ptr, columns, aligned[2], power,
        flux/length, c, evaluation)


def _export_fields(
    workpiece: Path,
    heat_flux_W_per_m2: np.ndarray,
    full_heat_density_W_per_m3: np.ndarray,
    run_dir: Path,
) -> list[Path]:
    from ngsolve import H1, GridFunction, Mesh

    from radia.gmsh_post_export import GmshPostExport

    mesh = Mesh(str(workpiece))
    fes = H1(mesh, order=1)
    qsurf = GridFunction(fes)
    qsurf.vec.FV().NumPy()[:] = heat_flux_W_per_m2
    density = GridFunction(fes)
    density.vec.FV().NumPy()[:] = full_heat_density_W_per_m3

    volume_path = run_dir / "ih_heat_density.msh"
    volume = GmshPostExport(mesh)
    volume.add_scalar_field("heat_density_W_per_m3_at_1A", density)
    volume.write(str(volume_path))
    surface_path = run_dir / "ih_surface_heat_flux.msh"
    # The scalar trace is carried by the 2D H1 grid function. The current
    # exporter has no 1D boundary topology route; retain the meridian domain
    # instead of silently writing a zero-element boundary artifact.
    surface = GmshPostExport(mesh, boundary=(mesh.dim == 3))
    surface.add_scalar_field("heat_flux_W_per_m2_at_1A", qsurf)
    surface.write(str(surface_path))
    return [volume_path.resolve(), surface_path.resolve()]
