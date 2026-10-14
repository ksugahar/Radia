"""Replayable transient heat execution for dimension-2 Netgen ``.vol`` meshes."""

from __future__ import annotations

import csv
import io
import json
import math
import time
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from typing import Any

import numpy as np

from .vol2d_circuit import _runtime_vol_path
from .vol2d_postprocess import _canonical, _export_entry, _package_version, _sha
from .vol2d_scalar import (
    _SAFE_NAME,
    _dense_matrix,
    _dof_indices,
    _finite,
    _gmsh_export,
    _prepare,
)


TRANSIENT_HEAT_SCHEMA = "radia.vol2d-transient-heat.v1"


def _time_axis(value: Any) -> np.ndarray:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) < 2:
        raise ValueError("time_s must contain at least two samples")
    times = np.asarray([_finite(item, "time_s") for item in value], dtype=float)
    if times[0] < 0.0 or np.any(np.diff(times) <= 0.0):
        raise ValueError("time_s must be nonnegative and strictly increasing")
    if times.size > 2000:
        raise ValueError("time_s is limited to 2000 samples")
    return times


def _history(
    raw: Any,
    names: Sequence[str],
    times: np.ndarray,
    defaults: Mapping[str, float],
    label: str,
) -> dict[str, np.ndarray]:
    if raw is None:
        raw = {}
    if not isinstance(raw, Mapping):
        raise ValueError(f"{label} must be an object keyed by name")
    unknown = sorted({str(key) for key in raw} - set(names))
    if unknown:
        raise ValueError(f"{label} contains unknown names: {unknown}")
    result: dict[str, np.ndarray] = {}
    for name in names:
        value = raw.get(name)
        if value is None:
            result[name] = np.full(times.size, float(defaults[name]), dtype=float)
            continue
        if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
            raise ValueError(f"{label}[{name}] must contain len(time_s) samples")
        row = np.asarray([_finite(item, f"{label}[{name}]") for item in value], dtype=float)
        if row.shape != times.shape:
            raise ValueError(f"{label}[{name}] must contain len(time_s) samples")
        result[name] = row
    return result


def _initial_temperature(raw: Any, material_names: Sequence[str]) -> dict[str, float]:
    if isinstance(raw, Mapping):
        supplied = {str(key) for key in raw}
        expected = set(material_names)
        if supplied != expected:
            raise ValueError(
                "initial_temperature_k must cover materials exactly: "
                f"missing={sorted(expected-supplied)}, unknown={sorted(supplied-expected)}"
            )
        return {
            name: _finite(raw[name], f"initial_temperature_k[{name}]")
            for name in material_names
        }
    value = _finite(raw, "initial_temperature_k")
    return {name: value for name in material_names}


def _csv_summary(times: np.ndarray, minima: list[float], maxima: list[float], energies: list[float], residuals: list[float]) -> str:
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(["time_s", "minimum_temperature_k", "maximum_temperature_k", "thermal_l2_energy_j_k", "residual_inf"])
    for row in zip(times, minima, maxima, energies, residuals):
        writer.writerow(row)
    return output.getvalue()


def solve_vol2d_transient_heat(request: Mapping[str, Any]) -> dict[str, Any]:
    """Solve ``rho*c*dT/dt - div(k grad T) = q`` with a theta method.

    Sources and prescribed temperatures may be sampled on an irregular time
    axis.  Material conductivity and volumetric heat capacity remain setup-time
    contracts.  The implementation stores the full accepted state history and
    verifies the semi-discrete residual at every step.
    """

    if not isinstance(request, Mapping):
        raise ValueError("request must be an object")
    started = time.perf_counter()
    prepared_request = dict(request)
    prepared_request["physics"] = "transient_heat"
    prepared, mesh_view, materials, boundaries = _prepare(prepared_request)
    times = _time_axis(request.get("time_s"))
    theta = _finite(request.get("theta", 1.0), "theta")
    if theta < 0.5 or theta > 1.0:
        raise ValueError("theta must be in [0.5, 1]")
    basename = str(request.get("export_basename", "vol2d_transient_heat")).strip()
    if not _SAFE_NAME.fullmatch(basename):
        raise ValueError("export_basename must be a portable filename stem")

    material_names = list(mesh_view.contract()["material_names"])
    initial_by_material = _initial_temperature(
        request.get("initial_temperature_k", 293.15), material_names
    )
    source_defaults = {
        name: float(materials["materials"][name]["volumetric_source_si"])
        for name in material_names
    }
    source_history = _history(
        request.get("volumetric_source_history_w_per_m3"),
        material_names,
        times,
        source_defaults,
        "volumetric_source_history_w_per_m3",
    )
    dirichlet_names = list(boundaries["dirichlet_values"])
    dirichlet_defaults = {
        name: float(boundaries["dirichlet_values"][name][0])
        for name in dirichlet_names
    }
    dirichlet_history = _history(
        request.get("dirichlet_history_k"),
        dirichlet_names,
        times,
        dirichlet_defaults,
        "dirichlet_history_k",
    )

    from ngsolve import BND, BilinearForm, GridFunction, H1, LinearForm, Mesh, grad, ds, dx, x  # type: ignore

    runtime = _runtime_vol_path(str(request["vol_text"]), mesh_view.content_sha256)
    mesh = Mesh(str(runtime))
    if prepared["family_contract"]["curved_geometry"]:
        mesh.Curve(prepared["family_contract"]["order"])
    fes = H1(
        mesh,
        order=prepared["family_contract"]["order"],
        dirichlet="|".join(dirichlet_names),
    )
    trial, test = fes.TnT()
    conductivity_x = mesh.MaterialCF({
        name: row["coefficient_si"][0] for name, row in materials["materials"].items()
    })
    conductivity_y = mesh.MaterialCF({
        name: row["coefficient_si"][1] for name, row in materials["materials"].items()
    })
    capacity = mesh.MaterialCF({
        name: row["volumetric_heat_capacity_j_per_m3_k"]
        for name, row in materials["materials"].items()
    })
    if prepared["formulation"] == "planar":
        measure = prepared["model_depth_m"]
        measure_contract = {"kind": "planar_depth", "model_depth_m": measure}
    else:
        measure = 2.0 * math.pi * x
        measure_contract = {"kind": "axisymmetric_full_revolution", "weight": "2*pi*r"}

    stiffness_form = BilinearForm(fes, symmetric=True)
    stiffness_form += (
        conductivity_x * grad(trial)[0] * grad(test)[0]
        + conductivity_y * grad(trial)[1] * grad(test)[1]
    ) * measure * dx
    robin_load = LinearForm(fes)
    for name, row in boundaries["robin_boundaries"].items():
        region = mesh.Boundaries(name)
        h = row["transfer_w_per_m2_k"]
        stiffness_form += h * trial * test * measure * ds(definedon=region)
        robin_load += h * row["ambient_k"] * test * measure * ds(definedon=region)
    capacity_form = BilinearForm(fes, symmetric=True)
    capacity_form += capacity * trial * test * measure * dx
    stiffness_form.Assemble()
    capacity_form.Assemble()
    robin_load.Assemble()
    stiffness = _dense_matrix(stiffness_form.mat, complex_matrix=False)
    mass = _dense_matrix(capacity_form.mat, complex_matrix=False)
    free = _dof_indices(fes.FreeDofs())
    fixed = np.asarray(sorted(set(range(fes.ndof)) - set(free.tolist())), dtype=int)
    if len(free) == 0:
        raise ValueError("transient heat solve has no free degrees of freedom")

    def prescribed(index: int) -> np.ndarray:
        values = {name: float(row[index]) for name, row in dirichlet_history.items()}
        lift = GridFunction(fes)
        if values:
            lift.Set(mesh.BoundaryCF(values, default=0.0), BND)
        return np.asarray(lift.vec.FV().NumPy(), dtype=float).copy()

    def load(index: int) -> np.ndarray:
        values = {name: float(row[index]) for name, row in source_history.items()}
        form = LinearForm(fes)
        form += mesh.MaterialCF(values) * test * measure * dx
        form.Assemble()
        result = np.asarray(form.vec.FV().NumPy(), dtype=float).copy()
        result += np.asarray(robin_load.vec.FV().NumPy(), dtype=float)
        return result

    initial = GridFunction(fes)
    initial.Set(mesh.MaterialCF(initial_by_material))
    state = np.asarray(initial.vec.FV().NumPy(), dtype=float).copy()
    boundary = prescribed(0)
    state[fixed] = boundary[fixed]
    states = [state.copy()]
    residuals = [0.0]
    factorization_count = 0
    loads = [load(index) for index in range(times.size)]
    assembly_s = time.perf_counter() - started

    solve_started = time.perf_counter()
    for index, dt in enumerate(np.diff(times), start=1):
        previous = state
        lhs = mass / dt + theta * stiffness
        rhs = (
            (mass / dt - (1.0 - theta) * stiffness) @ previous
            + theta * loads[index]
            + (1.0 - theta) * loads[index - 1]
        )
        boundary = prescribed(index)
        reduced_rhs = rhs[free]
        if fixed.size:
            reduced_rhs = reduced_rhs - lhs[np.ix_(free, fixed)] @ boundary[fixed]
        try:
            candidate = np.linalg.solve(lhs[np.ix_(free, free)], reduced_rhs)
        except np.linalg.LinAlgError as exc:
            raise ValueError(f"transient heat step {index} matrix is singular") from exc
        factorization_count += 1
        state = boundary
        state[free] = candidate
        derivative = (state - previous) / dt
        weighted = theta * state + (1.0 - theta) * previous
        weighted_load = theta * loads[index] + (1.0 - theta) * loads[index - 1]
        residual = mass @ derivative + stiffness @ weighted - weighted_load
        states.append(state.copy())
        residuals.append(float(np.linalg.norm(residual[free], ord=np.inf)))
    solve_s = time.perf_counter() - solve_started

    post_started = time.perf_counter()
    energies = [float(0.5 * value @ mass @ value) for value in states]
    minima = [float(np.min(value)) for value in states]
    maxima = [float(np.max(value)) for value in states]
    zero_external = (
        all(np.max(np.abs(row)) == 0.0 for row in source_history.values())
        and all(np.max(np.abs(row)) == 0.0 for row in dirichlet_history.values())
        and all(abs(row["ambient_k"]) == 0.0 for row in boundaries["robin_boundaries"].values())
    )
    passive_decay = None
    if zero_external:
        passive_decay = all(
            right <= left + 1.0e-10 * max(1.0, abs(left))
            for left, right in zip(energies, energies[1:])
        )
    final = GridFunction(fes)
    final.vec.FV().NumPy()[:] = states[-1]
    field = -grad(final)
    request_contract = {
        "schema": "radia.vol2d-transient-heat-request.v1",
        "mesh_contract_sha256": mesh_view.contract()["contract_sha256"],
        "material_contract_sha256": materials["contract_sha256"],
        "boundary_contract_sha256": boundaries["contract_sha256"],
        "element_family": prepared["element_family"],
        "formulation": prepared["formulation"],
        "measure": measure_contract,
        "time_s": times.tolist(),
        "theta": theta,
        "initial_temperature_k": initial_by_material,
        "source_history": {name: row.tolist() for name, row in source_history.items()},
        "dirichlet_history": {name: row.tolist() for name, row in dirichlet_history.items()},
    }
    request_sha = _sha(request_contract)
    state_rows = [row.tolist() for row in states]
    result_contract = {
        "schema": TRANSIENT_HEAT_SCHEMA,
        "status": "solved",
        "request_contract": request_contract,
        "request_contract_sha256": request_sha,
        "field_state_history": state_rows,
        "field_state_history_sha256": _sha(state_rows),
        "thermal_l2_energy_history_j_k": energies,
        "minimum_temperature_history_k": minima,
        "maximum_temperature_history_k": maxima,
        "step_residual_inf": residuals,
        "maximum_step_residual_inf": max(residuals),
        "zero_external_power": zero_external,
        "passive_energy_decay": passive_decay,
        "factorization_count": factorization_count,
        "generated_vol_git_required": False,
    }
    gmsh = _gmsh_export(
        mesh,
        final,
        field,
        basename=basename,
        request_sha256=request_sha,
        complex_field=False,
    )
    csv_content = _csv_summary(times, minima, maxima, energies, residuals)
    json_content = _canonical(result_contract)
    exports = {
        "json": _export_entry(json_content, f"{basename}.json", "application/json"),
        "csv": _export_entry(csv_content, f"{basename}.csv", "text/csv"),
        "gmsh_msh": _export_entry(gmsh["gmsh_msh"], f"{basename}.msh", "model/mesh"),
        "gmsh_geo": _export_entry(gmsh["gmsh_geo"], f"{basename}.geo", "text/plain"),
        "gmsh_geo_opt": _export_entry(gmsh["gmsh_geo_opt"], f"{basename}.geo.opt", "text/plain"),
        "gmsh_msh_opt": _export_entry(gmsh["gmsh_msh_opt"], f"{basename}.msh.opt", "text/plain"),
    }
    result_contract["export_content_sha256"] = {
        name: row["sha256"] for name, row in exports.items()
    }
    post_s = time.perf_counter() - post_started
    total = assembly_s + solve_s + post_s
    return {
        "schema": TRANSIENT_HEAT_SCHEMA,
        "status": "solved",
        "operation": "transient_heat",
        "execution_version": {
            "radia_mcp": _package_version(),
            "ngsolve": getattr(__import__("ngsolve"), "__version__", "unknown"),
        },
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "timing_s": {
            "prepare_and_assemble": assembly_s,
            "factorize_and_solve": solve_s,
            "postprocess_and_export": post_s,
            "total": total,
        },
        "result_contract": result_contract,
        "exports": exports,
    }
