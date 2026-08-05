"""Assemble native Simulink IH operators directly from geometry files.

The workpiece is a checked Netgen ``.vol`` volume mesh.  A STEP coil selects
PEEC; a checked ``.vol`` coil selects BEM-A.  The electromagnetic solve is
performed once at unit current with the existing Radia BEM-SIBC workflow.  For
the single scalar current port and linear material law this is an exact response
basis: the field is linear in current and the distributed loss scales with
``current**2``.  It is not a LUT or a lumped thermal surrogate.

The thermal mass, conductivity, convection, and heat-load projection operators
are assembled with NGSolve H1 forms on the workpiece mesh.  The resulting JSON
is consumed directly by Radia's Level-2 MATLAB S-Functions and standalone MEX
object handles; Python is not called per Simulink time step.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import io
import json
import math
import platform
import shutil
import sys
import time
import traceback
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import numpy as np

CONFIG_SCHEMA = "radia.ih.simulink.native_sfunction.v1"
ASSEMBLY_SCHEMA = "radia.ih.operator_assembly.v1"
RESULT_SCHEMA = "radia.simulink.application_run.v1"
_WORKPIECE_CONTRACT = "ih_workpiece_v1.json"
_COIL_CONTRACT = "ih_coil_bema_v1.json"


@dataclass(frozen=True, slots=True)
class IHOperatorAssemblyOptions:
    """Physical and discretization settings not inferable from geometry."""

    frequency_hz: float = 7000.0
    coil_conductivity_S_per_m: float = 5.8e7
    workpiece_conductivity_S_per_m: float = 5.0e6
    workpiece_relative_permeability: float = 100.0
    density_kg_per_m3: float = 7800.0
    heat_capacity_J_per_kgK: float = 467.0
    thermal_conductivity_W_per_mK: float = 46.6
    convection_W_per_m2K: float = 10.0
    initial_temperature_K: float = 293.15
    sample_time_s: float = 0.5
    workpiece_label: str = "sibc"
    coil_body_label: str = "body"
    coil_source_label: str = "source"
    coil_sink_label: str = "sink"
    peec_perimeter_filaments: int = 16
    peec_proximity: bool = True
    coupling_mode: str = "weak"
    workpiece_bem_backend: str = "intree-dense"
    thermal_order: int = 1

    def checked(self) -> "IHOperatorAssemblyOptions":
        positive = {
            "frequency_hz": self.frequency_hz,
            "coil_conductivity_S_per_m": self.coil_conductivity_S_per_m,
            "workpiece_conductivity_S_per_m": self.workpiece_conductivity_S_per_m,
            "workpiece_relative_permeability": self.workpiece_relative_permeability,
            "density_kg_per_m3": self.density_kg_per_m3,
            "heat_capacity_J_per_kgK": self.heat_capacity_J_per_kgK,
            "thermal_conductivity_W_per_mK": self.thermal_conductivity_W_per_mK,
            "initial_temperature_K": self.initial_temperature_K,
            "sample_time_s": self.sample_time_s,
        }
        for name, value in positive.items():
            if not math.isfinite(float(value)) or float(value) <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
        if not math.isfinite(float(self.convection_W_per_m2K)) or self.convection_W_per_m2K < 0.0:
            raise ValueError("convection_W_per_m2K must be finite and non-negative")
        for name in (
            "workpiece_label",
            "coil_body_label",
            "coil_source_label",
            "coil_sink_label",
        ):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"{name} must not be empty")
        canonical_labels = {
            "workpiece_label": "sibc",
            "coil_body_label": "body",
            "coil_source_label": "source",
            "coil_sink_label": "sink",
        }
        for name, expected in canonical_labels.items():
            if getattr(self, name) != expected:
                raise ValueError(
                    f"{name} is fixed by the versioned IH label contract as "
                    f"{expected!r}; relabel the .vol before assembly"
                )
        if self.peec_perimeter_filaments <= 0:
            raise ValueError("peec_perimeter_filaments must be positive")
        if self.coupling_mode not in {"weak", "strong"}:
            raise ValueError("coupling_mode must be weak or strong")
        if self.workpiece_bem_backend not in {"intree-dense", "hacapk"}:
            raise ValueError("workpiece_bem_backend must be intree-dense or hacapk")
        if self.thermal_order != 1:
            raise ValueError(
                "thermal_order must be 1 because the current qsurf.sol contract "
                "is an H1 P1 field"
            )
        return self


@dataclass(frozen=True, slots=True)
class UnitCurrentResult:
    heat_flux_W_per_m2: np.ndarray
    solver_payload: dict[str, Any]
    qsurf_solution: Path
    field_mesh: Path


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


class _Tee(io.TextIOBase):
    def __init__(self, *streams: io.TextIOBase) -> None:
        self._streams = streams

    def write(self, text: str) -> int:
        for stream in self._streams:
            stream.write(text)
            stream.flush()
        return len(text)

    def flush(self) -> None:
        for stream in self._streams:
            stream.flush()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _radia_version() -> str:
    try:
        return importlib.metadata.version("radia")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _geometry_stem(path: Path) -> str:
    name = path.name
    lower = name.lower()
    for suffix in (".vol.gz", ".step", ".stp", ".vol"):
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def default_output_path(workpiece: str | Path, coil: str | Path) -> Path:
    workpiece_path = Path(workpiece).resolve()
    coil_path = Path(coil).resolve()
    return workpiece_path.with_name(
        f"{_geometry_stem(workpiece_path)}_{_geometry_stem(coil_path)}_ih_native.json"
    )


def default_run_directory(output: str | Path) -> Path:
    output_path = Path(output).resolve()
    return output_path.with_name(output_path.stem + "_artifacts")


def _geometry_kind(path: Path) -> str:
    lower = path.name.lower()
    if lower.endswith((".step", ".stp")):
        return "peec"
    if lower.endswith((".vol", ".vol.gz")):
        return "bem-a"
    raise ValueError(
        f"coil geometry must end in .step/.stp (PEEC) or .vol/.vol.gz " f"(BEM-A): {path}"
    )


def _checked_geometry_paths(workpiece: str | Path, coil: str | Path) -> tuple[Path, Path, str]:
    workpiece_path = Path(workpiece).resolve()
    coil_path = Path(coil).resolve()
    if not workpiece_path.is_file():
        raise FileNotFoundError(f"workpiece mesh not found: {workpiece_path}")
    if not workpiece_path.name.lower().endswith((".vol", ".vol.gz")):
        raise ValueError(f"workpiece geometry must be a .vol/.vol.gz mesh: {workpiece_path}")
    if not coil_path.is_file():
        raise FileNotFoundError(f"coil geometry not found: {coil_path}")
    if workpiece_path == coil_path:
        raise ValueError("workpiece and coil geometry must be different files")
    return workpiece_path, coil_path, _geometry_kind(coil_path)


def _contract_path(name: str) -> Path:
    path = Path(__file__).resolve().parent / "contracts" / name
    if not path.is_file():
        raise FileNotFoundError(f"packaged IH label contract is missing: {path}")
    return path


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _materialize_solver_vol(source: Path, run_dir: Path, role: str) -> tuple[Path, list[Path]]:
    """Return an NGSolve-safe .vol, normalizing compressed CRLF input."""

    if not source.name.lower().endswith(".vol.gz"):
        return source, []

    materialized = (run_dir / f"{role}.solver.vol").resolve()
    with gzip.open(source, "rb") as compressed, materialized.open("wb") as target:
        for line in compressed:
            target.write(line.rstrip(b"\r\n") + b"\n")
    if materialized.stat().st_size == 0:
        raise RuntimeError(f"compressed mesh is empty: {source}")

    artifacts = [materialized]
    uncompressed_name = Path(str(source)[: -len(".gz")])
    companion = Path(str(uncompressed_name) + ".json")
    if companion.is_file():
        copied_companion = Path(str(materialized) + ".json")
        shutil.copyfile(companion, copied_companion)
        artifacts.append(copied_companion.resolve())
    return materialized, artifacts


def _check_vol(
    vol_path: Path,
    contract_path: Path,
    report_path: Path,
    *,
    workpiece: bool,
    workpiece_label: str,
) -> dict[str, Any]:
    try:
        from cubit_mesh_export.check import check_consistency
    except ImportError as exc:
        raise RuntimeError(
            "cubit-mesh-export is required for the mandatory check-vol gate; "
            "install radia[cubit]"
        ) from exc

    keyword: dict[str, Any] = {
        "contract": contract_path,
        "strict_labels": True,
    }
    if workpiece:
        keyword.update(
            conductive_materials=("workpiece",),
            sibc_boundaries=(workpiece_label,),
            require_all_sibc_labeled=True,
        )
    report = check_consistency(vol_path, **keyword)
    _write_json(report_path, report)
    if not report.get("passed", False):
        warnings = report.get("warnings", [])
        detail = "; ".join(str(item) for item in warnings) or "unknown finding"
        raise RuntimeError(f"check-vol failed for {vol_path}: {detail}")
    return report


def _solve_unit_current(
    workpiece: Path,
    coil: Path,
    backend: str,
    options: IHOperatorAssemblyOptions,
    run_dir: Path,
) -> UnitCurrentResult:
    from ngsolve import H1, GridFunction, Mesh

    from radia.panels.calc_inductance import build_argparser, run_inductance

    if backend == "peec" and options.coupling_mode == "strong" and options.peec_proximity:
        raise ValueError(
            "strong PEEC coupling cannot be combined with PEEC proximity; "
            "pass --no-peec-proximity or use weak coupling"
        )

    field_path = run_dir / "ih_fields.msh"
    argv = [
        "--coil-solver",
        backend,
        "--vol",
        str(workpiece),
        "--wp-label",
        options.workpiece_label,
        "--sigma",
        format(options.workpiece_conductivity_S_per_m, ".17g"),
        "--mu-r",
        format(options.workpiece_relative_permeability, ".17g"),
        "--frequency",
        format(options.frequency_hz, ".17g"),
        "--current",
        "1",
        "--coil-sigma",
        format(options.coil_conductivity_S_per_m, ".17g"),
        "--coupling-mode",
        options.coupling_mode,
        "--wp-bem-backend",
        options.workpiece_bem_backend,
        "--h1-order",
        "1",
        "--wp-loop-dof",
        "auto",
        "--impedance-model",
        "sibc",
        "--msh-output",
        str(field_path),
    ]
    if backend == "peec":
        argv.extend(
            [
                "--coil-step",
                str(coil),
                "--peec-n-peri",
                str(options.peec_perimeter_filaments),
                "--peec-proximity" if options.peec_proximity else "--no-peec-proximity",
            ]
        )
    else:
        argv.extend(
            [
                "--coil-vol",
                str(coil),
                "--coil-source-name",
                options.coil_source_label,
                "--coil-sink-name",
                options.coil_sink_label,
            ]
        )

    arguments = build_argparser().parse_args(argv)
    payload = run_inductance(arguments)
    if payload.get("status") == "error" or payload.get("error"):
        raise RuntimeError(str(payload.get("error", "IH electromagnetic solve failed")))
    qsurf_solution = Path(str(payload.get("qsurf_sol", ""))).resolve()
    if not qsurf_solution.is_file():
        raise RuntimeError(
            "the IH electromagnetic solve did not produce the required qsurf.sol field"
        )
    field_mesh = Path(str(payload.get("msh_file", ""))).resolve()
    if not field_mesh.is_file():
        raise RuntimeError("the IH electromagnetic solve did not produce the required GMSH field")
    mesh = Mesh(str(workpiece))
    fes = H1(mesh, order=1)
    qsurf = GridFunction(fes)
    qsurf.Load(str(qsurf_solution))
    values = np.asarray(qsurf.vec.FV().NumPy(), dtype=float).copy()
    if values.size != fes.ndof or not np.all(np.isfinite(values)):
        raise RuntimeError("qsurf.sol does not match the workpiece H1 P1 space")
    return UnitCurrentResult(values, payload, qsurf_solution, field_mesh)


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
    options: IHOperatorAssemblyOptions,
) -> ThermalOperators:
    from ngsolve import CF, H1, BilinearForm, Mesh, TaskManager, ds, dx, grad

    mesh = Mesh(str(workpiece))
    if mesh.dim != 3 or mesh.ne <= 0:
        raise ValueError("the IH thermal workpiece must be a 3D volume .vol mesh")
    materials = sorted(set(str(name) for name in mesh.GetMaterials()))
    if materials != ["workpiece"]:
        raise ValueError(
            "the IH thermal mesh must contain exactly one material named "
            f"'workpiece'; got {materials}"
        )
    if options.workpiece_label not in set(str(name) for name in mesh.GetBoundaries()):
        raise ValueError(f"workpiece boundary {options.workpiece_label!r} is absent from the mesh")

    fes = H1(mesh, order=options.thermal_order)
    if heat_flux_W_per_m2.size != fes.ndof:
        raise ValueError("the unit-current heat field and thermal H1 space have different sizes")
    u, v = fes.TnT()
    rho_cp = options.density_kg_per_m3 * options.heat_capacity_J_per_kgK
    with TaskManager():
        mass = BilinearForm(fes, symmetric=True)
        mass += CF(rho_cp) * u * v * dx
        mass.Assemble()
        stiffness = BilinearForm(fes, symmetric=True)
        stiffness += CF(options.thermal_conductivity_W_per_mK) * grad(u) * grad(v) * dx
        stiffness.Assemble()
        heat_mass = BilinearForm(fes, symmetric=True, check_unused=False)
        heat_mass += u * v * ds(options.workpiece_label)
        heat_mass.Assemble()
        convection = BilinearForm(fes, symmetric=True, check_unused=False)
        convection += u * v * ds
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
    if np.linalg.norm(heat_load[outside], ord=np.inf) > 1.0e-11 * max(
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
    surface = GmshPostExport(mesh, boundary=True)
    surface.add_scalar_field("heat_flux_W_per_m2_at_1A", qsurf)
    surface.write(str(surface_path))
    return [volume_path.resolve(), surface_path.resolve()]


def _native_config(
    workpiece: Path,
    coil: Path,
    solver_workpiece: Path,
    solver_coil: Path,
    backend: str,
    options: IHOperatorAssemblyOptions,
    electromagnetic: UnitCurrentResult,
    thermal: ThermalOperators,
    reports: list[Path],
    contracts: list[Path],
    gmsh_files: list[Path],
) -> dict[str, Any]:
    genus = int(electromagnetic.solver_payload.get("wp_genus", 0))
    if genus >= 1 and not electromagnetic.solver_payload.get("wp_loop_dof", False):
        reason = electromagnetic.solver_payload.get(
            "wp_loop_dof_skip_reason", "no cohomology mode was reported"
        )
        raise RuntimeError(
            f"workpiece genus={genus} requires a solved cohomology loop DOF; "
            f"the electromagnetic solve skipped it: {reason}. Use the "
            "intree-dense workpiece backend and a supported genus-1 mesh."
        )
    solver_power = float(electromagnetic.solver_payload.get("P_wp_W", math.nan))
    if not math.isfinite(solver_power) or solver_power <= 0.0:
        raise RuntimeError("the electromagnetic result does not contain positive P_wp_W")
    relative_power_error = abs(thermal.heat_power_W - solver_power) / solver_power
    if relative_power_error > 1.0e-6:
        raise RuntimeError(
            "the assembled thermal source does not preserve electromagnetic power: "
            f"relative error {relative_power_error:.3e}"
        )

    n_temperature = thermal.n_temperature
    n_heat = int(thermal.heat_dofs.size)
    config: dict[str, Any] = {
        "schema": CONFIG_SCHEMA,
        "assembly_schema": ASSEMBLY_SCHEMA,
        "backend": "matlab-level2+radia-mex-handles",
        "python_fallback": False,
        "operator_assembly": "preassembled",
        "operator_basis": "exact-single-current-linear-response",
        "distributed_field": True,
        "surrogate": False,
        "release_channel": "preview",
        "n_eddy_unknown": 1,
        "n_heat": n_heat,
        "n_temperature": n_temperature,
        "bh_mode": "linear",
        "eddy_matrix_real": [1.0],
        "eddy_matrix_imag": [0.0],
        "eddy_rhs_real": [1.0],
        "eddy_rhs_imag": [0.0],
        "heat_projection": thermal.unit_heat_density_W_per_m3.tolist(),
        "heat_cell_weights": thermal.heat_cell_weights_m3.tolist(),
        "heat_to_temperature_projection": thermal.heat_to_temperature_projection,
        "mass_row_ptr": thermal.mass_row_ptr,
        "mass_col": thermal.mass_col,
        "mass_value": thermal.mass_value,
        "stiffness_row_ptr": thermal.stiffness_row_ptr,
        "stiffness_col": thermal.stiffness_col,
        "stiffness_value": thermal.stiffness_value,
        "convection_row_ptr": thermal.convection_row_ptr,
        "convection_col": thermal.convection_col,
        "convection_value": thermal.convection_value,
        "temperature_cell_weights": thermal.temperature_cell_weights_J_per_K.tolist(),
        "initial_temperature_K": [options.initial_temperature_K] * n_temperature,
        "sample_time_s": options.sample_time_s,
        "thermal_tolerance": 1.0e-10,
        "thermal_max_iterations": 1000,
        "convection_W_per_m2K": options.convection_W_per_m2K,
        "rotation_mode": "none",
        "angle_origin_rad": 0.0,
        "eddy_solver": backend,
        "eddy_method": (
            "PEEC + BEM-SIBC unit-current response"
            if backend == "peec"
            else "BEM-A + BEM-SIBC unit-current response"
        ),
        "linear_solver": options.workpiece_bem_backend,
        "thermal_solver": "fem",
        "thermal_mesh_type": "3D volume",
        "current_change_recomputes_eddy": False,
        "temperature_change_recomputes_eddy": False,
        "temperature_coordinate_system": "workpiece",
        "rotation_transport": "none",
        "dt_order": "eddy;transport(theta_prev,theta_now);thermal",
        "vol_check_required": True,
        "vol_check_reports": [str(path) for path in reports],
        "workpiece_vol_label_contract": str(contracts[0]),
        "geometry": {
            "workpiece_vol": str(workpiece),
            "coil_file": str(coil),
            "coil_backend": backend,
            "workpiece_sha256": _sha256(workpiece),
            "coil_sha256": _sha256(coil),
            "solver_workpiece_vol": str(solver_workpiece),
            "solver_coil_file": str(solver_coil),
            "gzip_materialized": bool(solver_workpiece != workpiece or solver_coil != coil),
        },
        "physical_parameters": asdict(options),
        "unit_current": {
            "reference_current_A": 1.0,
            "electromagnetic_power_W": solver_power,
            "thermal_source_power_W": thermal.heat_power_W,
            "relative_power_error": relative_power_error,
            "qsurf_solution": str(electromagnetic.qsurf_solution),
        },
        "artifacts": {
            "gmsh_format": "msh-v4.1",
            "gmsh": [str(path) for path in gmsh_files],
        },
        "generated_at_utc": _utc_now(),
        "runtime_radia_version": _radia_version(),
        "runtime_python": platform.python_version(),
        "runtime_platform": platform.platform(),
    }
    if backend == "bem-a":
        config["coil_vol_label_contract"] = str(contracts[1])
    return config


def _result_payload(
    *,
    status: str,
    started_at: str,
    elapsed_s: float,
    output: Path,
    run_dir: Path,
    workpiece: Path,
    coil: Path,
    error: str | None,
    artifacts: list[Path],
    primary_power_W: float | None,
    config_current: bool,
) -> dict[str, Any]:
    return {
        "radia_result": {
            "schema": RESULT_SCHEMA,
            "application": "ih-operator-assembly",
            "backend": "python-headless-cli",
            "status": status,
            "returncode": 0 if status == "passed" else 1,
            "error": error,
            "executed_at_utc": started_at,
            "completed_at_utc": _utc_now(),
            "elapsed_s": round(elapsed_s, 6),
            "runtime_radia_version": _radia_version(),
            "runtime_python": platform.python_version(),
            "runtime_platform": platform.platform(),
            "config": str(output) if config_current and output.is_file() else None,
            "config_sha256": (_sha256(output) if config_current and output.is_file() else None),
            "run_dir": str(run_dir),
            "log": str(run_dir / "run.log"),
            "geometry": {
                "workpiece": str(workpiece),
                "coil": str(coil),
            },
            "primary": {"key": "P_wp_W", "value": primary_power_W},
            "artifacts": {
                "gmsh_policy": "required",
                "gmsh_format": "msh-v4.1",
                "gmsh": [str(path) for path in artifacts if path.suffix == ".msh"],
                "all": [str(path) for path in artifacts],
            },
        }
    }


def assemble_ih_operators(
    workpiece: str | Path,
    coil: str | Path,
    *,
    output: str | Path | None = None,
    run_dir: str | Path | None = None,
    options: IHOperatorAssemblyOptions | None = None,
) -> dict[str, Any]:
    """Generate a checked native IH configuration from two geometry files."""

    options = (options or IHOperatorAssemblyOptions()).checked()
    workpiece_path, coil_path, backend = _checked_geometry_paths(workpiece, coil)
    output_path = (
        Path(output).resolve()
        if output is not None
        else default_output_path(workpiece_path, coil_path)
    )
    if output_path.suffix.lower() != ".json":
        raise ValueError("the built-in IH operator assembler writes a .json config")
    run_path = (
        Path(run_dir).resolve() if run_dir is not None else default_run_directory(output_path)
    )
    run_path.mkdir(parents=True, exist_ok=True)
    log_path = run_path / "run.log"
    result_path = run_path / "result.json"
    started_at = _utc_now()
    started = time.monotonic()
    artifacts: list[Path] = []
    primary_power: float | None = None

    with log_path.open("w", encoding="utf-8") as log:
        tee_out = _Tee(sys.stdout, log)
        tee_err = _Tee(sys.stderr, log)
        try:
            with redirect_stdout(tee_out), redirect_stderr(tee_err):
                print(f"workpiece: {workpiece_path}")
                print(f"coil ({backend}): {coil_path}")
                print(f"output: {output_path}")
                workpiece_contract = _contract_path(_WORKPIECE_CONTRACT)
                workpiece_report = run_path / "workpiece.vol-check.json"
                _check_vol(
                    workpiece_path,
                    workpiece_contract,
                    workpiece_report,
                    workpiece=True,
                    workpiece_label=options.workpiece_label,
                )
                reports = [workpiece_report.resolve()]
                contracts = [workpiece_contract.resolve()]
                if backend == "bem-a":
                    coil_contract = _contract_path(_COIL_CONTRACT)
                    coil_report = run_path / "coil.vol-check.json"
                    _check_vol(
                        coil_path,
                        coil_contract,
                        coil_report,
                        workpiece=False,
                        workpiece_label=options.workpiece_label,
                    )
                    reports.append(coil_report.resolve())
                    contracts.append(coil_contract.resolve())
                artifacts.extend(reports)

                solver_workpiece, workpiece_solver_artifacts = _materialize_solver_vol(
                    workpiece_path, run_path, "workpiece"
                )
                artifacts.extend(workpiece_solver_artifacts)
                if backend == "bem-a":
                    solver_coil, coil_solver_artifacts = _materialize_solver_vol(
                        coil_path, run_path, "coil"
                    )
                    artifacts.extend(coil_solver_artifacts)
                else:
                    solver_coil = coil_path

                electromagnetic = _solve_unit_current(
                    solver_workpiece, solver_coil, backend, options, run_path
                )
                electromagnetic_result = run_path / "electromagnetic_result.json"
                _write_json(electromagnetic_result, electromagnetic.solver_payload)
                artifacts.append(electromagnetic_result.resolve())
                thermal = _assemble_thermal_operators(
                    solver_workpiece,
                    electromagnetic.heat_flux_W_per_m2,
                    options,
                )
                thermal_gmsh_files = _export_fields(
                    solver_workpiece,
                    electromagnetic.heat_flux_W_per_m2,
                    thermal.full_heat_density_W_per_m3,
                    run_path,
                )
                gmsh_files = [electromagnetic.field_mesh, *thermal_gmsh_files]
                artifacts.extend(gmsh_files)
                artifacts.append(electromagnetic.qsurf_solution)
                config = _native_config(
                    workpiece_path,
                    coil_path,
                    solver_workpiece,
                    solver_coil,
                    backend,
                    options,
                    electromagnetic,
                    thermal,
                    reports,
                    contracts,
                    gmsh_files,
                )
                _write_json(output_path, config)
                artifacts.append(output_path)
                primary_power = float(config["unit_current"]["electromagnetic_power_W"])
                print(
                    f"assembled n_temperature={thermal.n_temperature} "
                    f"n_heat={thermal.heat_dofs.size} "
                    f"P_wp(1 A)={primary_power:.9g} W"
                )
        except Exception as exc:
            log.write("\n" + traceback.format_exc())
            payload = _result_payload(
                status="failed",
                started_at=started_at,
                elapsed_s=time.monotonic() - started,
                output=output_path,
                run_dir=run_path,
                workpiece=workpiece_path,
                coil=coil_path,
                error=str(exc),
                artifacts=artifacts,
                primary_power_W=primary_power,
                config_current=False,
            )
            _write_json(result_path, payload)
            raise

    artifacts.append(log_path)
    payload = _result_payload(
        status="passed",
        started_at=started_at,
        elapsed_s=time.monotonic() - started,
        output=output_path,
        run_dir=run_path,
        workpiece=workpiece_path,
        coil=coil_path,
        error=None,
        artifacts=artifacts,
        primary_power_W=primary_power,
        config_current=True,
    )
    _write_json(result_path, payload)
    return config


def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workpiece", type=Path, help="checked workpiece .vol/.vol.gz")
    parser.add_argument("coil", type=Path, help="coil .step/.stp (PEEC) or .vol/.vol.gz (BEM-A)")
    parser.add_argument("-o", "--output", type=Path, default=None, help="native IH JSON config")
    parser.add_argument(
        "--run-dir", type=Path, default=None, help="persistent report/artifact directory"
    )
    parser.add_argument("--frequency", type=float, default=7000.0, dest="frequency_hz")
    parser.add_argument("--coil-sigma", type=float, default=5.8e7)
    parser.add_argument("--workpiece-sigma", type=float, default=5.0e6)
    parser.add_argument("--workpiece-mu-r", type=float, default=100.0)
    parser.add_argument("--density", type=float, default=7800.0)
    parser.add_argument("--heat-capacity", type=float, default=467.0)
    parser.add_argument("--thermal-conductivity", type=float, default=46.6)
    parser.add_argument("--convection", type=float, default=10.0)
    parser.add_argument("--initial-temperature-K", type=float, default=293.15)
    parser.add_argument("--sample-time", type=float, default=0.5)
    parser.add_argument(
        "--workpiece-label",
        choices=("sibc",),
        default="sibc",
        help="fixed by ih_workpiece_v1: sibc",
    )
    parser.add_argument(
        "--coil-body-label",
        choices=("body",),
        default="body",
        help="fixed by ih_coil_bema_v1: body",
    )
    parser.add_argument(
        "--coil-source-label",
        choices=("source",),
        default="source",
        help="fixed by ih_coil_bema_v1: source",
    )
    parser.add_argument(
        "--coil-sink-label",
        choices=("sink",),
        default="sink",
        help="fixed by ih_coil_bema_v1: sink",
    )
    parser.add_argument("--peec-n-peri", type=int, default=16)
    parser.add_argument(
        "--peec-proximity",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--coupling-mode", choices=("weak", "strong"), default="weak")
    parser.add_argument(
        "--workpiece-bem-backend",
        choices=("intree-dense", "hacapk"),
        default="intree-dense",
    )
    return parser


def _options_from_args(args: argparse.Namespace) -> IHOperatorAssemblyOptions:
    return IHOperatorAssemblyOptions(
        frequency_hz=args.frequency_hz,
        coil_conductivity_S_per_m=args.coil_sigma,
        workpiece_conductivity_S_per_m=args.workpiece_sigma,
        workpiece_relative_permeability=args.workpiece_mu_r,
        density_kg_per_m3=args.density,
        heat_capacity_J_per_kgK=args.heat_capacity,
        thermal_conductivity_W_per_mK=args.thermal_conductivity,
        convection_W_per_m2K=args.convection,
        initial_temperature_K=args.initial_temperature_K,
        sample_time_s=args.sample_time,
        workpiece_label=args.workpiece_label,
        coil_body_label=args.coil_body_label,
        coil_source_label=args.coil_source_label,
        coil_sink_label=args.coil_sink_label,
        peec_perimeter_filaments=args.peec_n_peri,
        peec_proximity=args.peec_proximity,
        coupling_mode=args.coupling_mode,
        workpiece_bem_backend=args.workpiece_bem_backend,
    )


def main(argv: Sequence[str] | None = None) -> int:
    args = build_argparser().parse_args(argv)
    output = args.output or default_output_path(args.workpiece, args.coil)
    try:
        assemble_ih_operators(
            args.workpiece,
            args.coil,
            output=output,
            run_dir=args.run_dir,
            options=_options_from_args(args),
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(Path(output).resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
