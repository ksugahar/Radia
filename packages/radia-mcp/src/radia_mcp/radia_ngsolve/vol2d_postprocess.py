"""Solved-field postprocessing for dimension-2 Netgen ``.vol`` models.

One linear field operator and one Cholesky factor are reused for every current
row.  Planar H1 and axisymmetric H1Henrotte fields then share a compact result
contract for point probes, oriented paths, material integrals, and Gmsh 4.1
handoff.  The MCP surface returns content and hashes; it never writes to a
caller-selected path.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import re
import time
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import numpy as np

from .vol2d_dynamics import _open_space, _prepare_request, assemble_vol2d_dynamics


POSTPROCESS_SCHEMA = "radia.vol2d-postprocess-analysis.v1"
REPLAY_SCHEMA = "radia.vol2d-postprocess-replay.v1"
_MAX_SWEEP_ROWS = 32
_MAX_POINT_PROBES = 128
_MAX_PATH_SAMPLES = 2048
_SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")


def _package_version() -> str:
    try:
        return version("radia-mcp")
    except PackageNotFoundError:
        return "source-tree"


def _canonical(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, separators=(",", ": ")) + "\n"


def _sha(value: Any) -> str:
    if isinstance(value, str):
        payload = value.encode("utf-8")
    else:
        payload = json.dumps(
            value, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _finite(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be finite")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def _positive(value: Any, label: str) -> float:
    result = _finite(value, label)
    if result <= 0.0:
        raise ValueError(f"{label} must be positive")
    return result


def _sequence(value: Any, label: str, *, nonempty: bool = True) -> Sequence[Any]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{label} must be a sequence")
    if nonempty and not value:
        raise ValueError(f"{label} must be non-empty")
    return value


def _real_scalar(value: Any, label: str) -> float:
    result = complex(value)
    if not math.isfinite(result.real) or not math.isfinite(result.imag):
        raise ValueError(f"{label} is not finite")
    if abs(result.imag) > 1.0e-10 * max(1.0, abs(result.real)):
        raise ValueError(f"{label} unexpectedly contains an imaginary component")
    return float(result.real)


def _normalize_current_rows(raw: Any, branch_count: int) -> list[list[float]]:
    rows = _sequence(raw, "current_rows_a")
    if len(rows) > _MAX_SWEEP_ROWS:
        raise ValueError(f"current_rows_a is limited to {_MAX_SWEEP_ROWS} rows")
    normalized: list[list[float]] = []
    for row_index, row in enumerate(rows):
        values = _sequence(row, f"current_rows_a[{row_index}]")
        if len(values) != branch_count:
            raise ValueError(
                f"current_rows_a[{row_index}] must contain {branch_count} values"
            )
        normalized.append(
            [_finite(value, f"current_rows_a[{row_index}][{index}]") for index, value in enumerate(values)]
        )
    identities = [_sha(row) for row in normalized]
    if len(set(identities)) != len(identities):
        raise ValueError("current_rows_a must not contain duplicate operating points")
    return normalized


def _named_points(raw: Any, coordinate_names: tuple[str, str]) -> list[dict[str, Any]]:
    rows = _sequence(raw, "point_probes")
    if len(rows) > _MAX_POINT_PROBES:
        raise ValueError(f"point_probes is limited to {_MAX_POINT_PROBES} rows")
    result: list[dict[str, Any]] = []
    names: set[str] = set()
    for index, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise ValueError(f"point_probes[{index}] must be an object")
        name = str(row.get("name", "")).strip()
        if not name or name in names:
            raise ValueError("point probe names must be non-empty and unique")
        coordinates = _sequence(row.get("coordinates_m"), f"point_probes[{index}].coordinates_m")
        if len(coordinates) != 2:
            raise ValueError("point probe coordinates_m must contain exactly two values")
        point = [
            _finite(coordinates[0], f"point_probes[{index}].coordinates_m[0]"),
            _finite(coordinates[1], f"point_probes[{index}].coordinates_m[1]"),
        ]
        names.add(name)
        result.append(
            {"name": name, "coordinate_names": list(coordinate_names), "coordinates_m": point}
        )
    return result


def _named_paths(raw: Any, coordinate_names: tuple[str, str]) -> list[dict[str, Any]]:
    rows = _sequence(raw, "path_probes")
    result: list[dict[str, Any]] = []
    names: set[str] = set()
    sample_count = 0
    for index, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise ValueError(f"path_probes[{index}] must be an object")
        name = str(row.get("name", "")).strip()
        if not name or name in names:
            raise ValueError("path probe names must be non-empty and unique")
        raw_points = _sequence(row.get("points_m"), f"path_probes[{index}].points_m")
        if len(raw_points) < 2:
            raise ValueError("each path probe needs at least two points")
        points: list[list[float]] = []
        for point_index, raw_point in enumerate(raw_points):
            coordinates = _sequence(
                raw_point, f"path_probes[{index}].points_m[{point_index}]"
            )
            if len(coordinates) != 2:
                raise ValueError("path points must contain exactly two coordinates")
            point = [
                _finite(coordinates[0], "path coordinate"),
                _finite(coordinates[1], "path coordinate"),
            ]
            if points and point == points[-1]:
                raise ValueError("a path must not contain adjacent duplicate points")
            points.append(point)
        samples_per_segment = int(row.get("samples_per_segment", 4))
        if samples_per_segment < 1 or samples_per_segment > 64:
            raise ValueError("samples_per_segment must be in [1, 64]")
        sample_count += samples_per_segment * (len(points) - 1)
        names.add(name)
        result.append(
            {
                "name": name,
                "coordinate_names": list(coordinate_names),
                "points_m": points,
                "samples_per_segment": samples_per_segment,
                "closed": points[0] == points[-1],
            }
        )
    if sample_count > _MAX_PATH_SAMPLES:
        raise ValueError(f"path sampling is limited to {_MAX_PATH_SAMPLES} midpoint samples")
    return result


def _linear_permeability(materials: Mapping[str, Any]) -> dict[str, float]:
    result: dict[str, float] = {}
    for name, material in materials["materials"].items():
        if material["kind"] != "linear":
            raise ValueError(
                "postprocessing sweeps require linear materials; use the nonlinear solver for one operating point"
            )
        result[name] = float(material["permeability_h_per_m"])
    return result


def _coordinate_contract(formulation: str) -> dict[str, Any]:
    if formulation == "planar":
        return {
            "coordinate_frame": "global_cartesian_xy",
            "coordinate_names": ["x", "y"],
            "potential_name": "A_z",
            "field_components": ["B_x", "B_y"],
            "path_orientation": "ordered_xy_with_left_normal",
        }
    return {
        "coordinate_frame": "meridional_rz",
        "coordinate_names": ["r", "z"],
        "potential_name": "A_phi",
        "field_components": ["B_r", "B_z"],
        "path_orientation": "ordered_rz_with_left_normal",
    }


def _field_coefficient(gfu: Any, formulation: str, *, axis_safe: bool = False) -> Any:
    from ngsolve import CoefficientFunction, IfPos, grad, x  # type: ignore

    if formulation == "planar":
        return CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    bz = grad(gfu)[0] + gfu / x
    if axis_safe:
        bz = IfPos(x - 1.0e-12, bz, 2.0 * grad(gfu)[0])
    return CoefficientFunction((-grad(gfu)[1], bz))


def _point_material(mesh: Any, mip: Any, material_names: Sequence[str]) -> str:
    for name in material_names:
        value = _real_scalar(mesh.MaterialCF({name: 1.0}, default=0.0)(mip), "material indicator")
        if value > 0.5:
            return name
    raise ValueError("probe point is outside all named mesh materials")


def _sample_point(
    mesh: Any,
    gfu: Any,
    field: Any,
    coordinates: Sequence[float],
    *,
    formulation: str,
    material_names: Sequence[str],
    permeability: Mapping[str, float],
) -> dict[str, Any]:
    first, second = float(coordinates[0]), float(coordinates[1])
    if formulation == "axisymmetric_henrotte" and first <= 1.0e-12:
        raise ValueError("axisymmetric point/path probes require r > 1e-12 m")
    try:
        mip = mesh(first, second)
        potential = _real_scalar(gfu(mip), "potential")
        raw_field = field(mip)
        components = [
            _real_scalar(raw_field[0], "field component 0"),
            _real_scalar(raw_field[1], "field component 1"),
        ]
        material = _point_material(mesh, mip, material_names)
    except Exception as exc:
        if isinstance(exc, ValueError):
            raise
        raise ValueError(f"probe point {list(coordinates)} is outside the solved mesh") from exc
    mu = float(permeability[material])
    return {
        "coordinates_m": [first, second],
        "material": material,
        "potential_wb_per_m": potential,
        "b_t": components,
        "b_magnitude_t": float(math.hypot(*components)),
        "h_a_per_m": [component / mu for component in components],
    }


def _sample_path(
    mesh: Any,
    gfu: Any,
    field: Any,
    path: Mapping[str, Any],
    *,
    formulation: str,
    material_names: Sequence[str],
    permeability: Mapping[str, float],
) -> dict[str, Any]:
    samples: list[dict[str, Any]] = []
    line_b = 0.0
    line_h = 0.0
    length = 0.0
    nsub = int(path["samples_per_segment"])
    for segment_index, (left, right) in enumerate(
        zip(path["points_m"], path["points_m"][1:])
    ):
        delta = np.asarray(right, dtype=float) - np.asarray(left, dtype=float)
        segment_length = float(np.linalg.norm(delta))
        if segment_length <= 0.0:
            raise ValueError("path segment length must be positive")
        tangent = delta / segment_length
        normal = np.asarray((-tangent[1], tangent[0]), dtype=float)
        ds = segment_length / nsub
        for local_index in range(nsub):
            fraction = (local_index + 0.5) / nsub
            coordinates = np.asarray(left, dtype=float) + fraction * delta
            sampled = _sample_point(
                mesh,
                gfu,
                field,
                coordinates,
                formulation=formulation,
                material_names=material_names,
                permeability=permeability,
            )
            b = np.asarray(sampled["b_t"], dtype=float)
            h = np.asarray(sampled["h_a_per_m"], dtype=float)
            b_tangent = float(b @ tangent)
            b_normal = float(b @ normal)
            h_tangent = float(h @ tangent)
            line_b += b_tangent * ds
            line_h += h_tangent * ds
            length += ds
            samples.append(
                {
                    **sampled,
                    "segment_index": segment_index,
                    "sample_index": local_index,
                    "tangent": tangent.tolist(),
                    "left_normal": normal.tolist(),
                    "ds_m": ds,
                    "b_tangent_t": b_tangent,
                    "b_left_normal_t": b_normal,
                    "h_tangent_a_per_m": h_tangent,
                }
            )
    return {
        "name": path["name"],
        "closed": bool(path["closed"]),
        "length_m": length,
        "b_tangent_line_integral_t_m": line_b,
        "h_tangent_line_integral_a": line_h,
        "samples": samples,
    }


def _region_integrals(
    mesh: Any,
    field: Any,
    names: Sequence[str],
    permeability: Mapping[str, float],
    *,
    formulation: str,
    model_depth_m: float | None,
) -> list[dict[str, Any]]:
    from ngsolve import Integrate, InnerProduct, x  # type: ignore

    rows: list[dict[str, Any]] = []
    for name in names:
        indicator = mesh.MaterialCF({name: 1.0}, default=0.0)
        density = 0.5 * InnerProduct(field, field) / float(permeability[name])
        if formulation == "planar":
            weight = float(model_depth_m)
            measure = _real_scalar(Integrate(indicator * weight, mesh), "planar volume")
            energy = _real_scalar(Integrate(indicator * density * weight, mesh), "magnetic energy")
            measure_name = "volume_m3"
        else:
            weight = 2.0 * math.pi * x
            measure = _real_scalar(Integrate(indicator * weight, mesh), "axisymmetric volume")
            energy = _real_scalar(Integrate(indicator * density * weight, mesh), "magnetic energy")
            measure_name = "full_revolution_volume_m3"
        rows.append(
            {
                "material": name,
                measure_name: measure,
                "magnetic_energy_j": energy,
            }
        )
    return rows


def _export_csv(branch_order: Sequence[str], sweep_rows: Sequence[Mapping[str, Any]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(
        [
            "row_index",
            *[f"current_a:{name}" for name in branch_order],
            "magnetic_energy_j",
            "maximum_probe_b_t",
            "field_residual_inf",
        ]
    )
    for row in sweep_rows:
        writer.writerow(
            [
                row["row_index"],
                *[format(value, ".17g") for value in row["branch_current_a"]],
                format(row["total_magnetic_energy_j"], ".17g"),
                format(row["maximum_point_probe_b_t"], ".17g"),
                format(row["field_residual_inf"], ".17g"),
            ]
        )
    return output.getvalue()


def _export_gmsh(
    mesh: Any,
    gfu: Any,
    field: Any,
    *,
    basename: str,
    request_sha256: str,
    potential_name: str,
) -> dict[str, str]:
    from ngsolve import CoefficientFunction, InnerProduct, sqrt  # type: ignore
    from radia.gmsh_post_export import GmshPostExport

    root = Path(r"C:\temp") / "radia_mcp_vol2d_postprocess" / request_sha256[:20]
    root.mkdir(parents=True, exist_ok=True)
    msh = root / f"{basename}.msh"
    post = GmshPostExport(mesh)
    post.add_scalar_field(potential_name, gfu)
    post.add_scalar_field("B_magnitude", sqrt(InnerProduct(field, field)))
    post.add_vector_field(
        "B_vector",
        CoefficientFunction((field[0], field[1], 0.0)),
    )
    post.write(str(msh))
    geo = msh.with_suffix(".geo")
    geo_opt = Path(str(geo) + ".opt")
    msh_opt = Path(str(msh) + ".opt")
    paths = {"gmsh_msh": msh, "gmsh_geo": geo, "gmsh_geo_opt": geo_opt, "gmsh_msh_opt": msh_opt}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise ValueError(f"Gmsh exporter did not create required companions: {missing}")
    return {name: path.read_text(encoding="utf-8") for name, path in paths.items()}


def _export_entry(content: str, filename: str, media_type: str) -> dict[str, Any]:
    return {
        "filename": filename,
        "media_type": media_type,
        "sha256": _sha(content),
        "bytes": len(content.encode("utf-8")),
        "content": content,
    }


def solve_vol2d_postprocess(request: Mapping[str, Any]) -> dict[str, Any]:
    """Solve a linear current sweep and emit one portable postprocess artifact."""

    if not isinstance(request, Mapping):
        raise ValueError("request must be an object")
    assembly_started = time.perf_counter()
    operators = assemble_vol2d_dynamics(request)
    prepared, mesh_view, material_contract = _prepare_request(request)
    mesh, fes, _ = _open_space(prepared, mesh_view)
    permeability = _linear_permeability(material_contract)
    assembly_s = time.perf_counter() - assembly_started

    formulation = prepared["formulation"]
    coordinates = _coordinate_contract(formulation)
    coordinate_names = tuple(coordinates["coordinate_names"])
    point_probes = _named_points(request.get("point_probes"), coordinate_names)
    path_probes = _named_paths(request.get("path_probes"), coordinate_names)
    region_names = [str(value) for value in _sequence(request.get("region_materials"), "region_materials")]
    if len(set(region_names)) != len(region_names):
        raise ValueError("region_materials must not contain duplicates")
    unknown_regions = sorted(set(region_names) - set(permeability))
    if unknown_regions:
        raise ValueError(f"region_materials are absent from .vol: {unknown_regions}")
    if formulation == "planar":
        model_depth_m: float | None = _positive(request.get("model_depth_m"), "model_depth_m")
    else:
        if request.get("model_depth_m") is not None:
            raise ValueError("axisymmetric full-revolution postprocessing must not specify model_depth_m")
        model_depth_m = None

    branch_order = list(operators["assembly"]["branch_order"])
    current_rows = _normalize_current_rows(request.get("current_rows_a"), len(branch_order))
    basename = str(request.get("export_basename", "vol2d_postprocess")).strip()
    if not _SAFE_NAME.fullmatch(basename):
        raise ValueError("export_basename must be a portable filename stem")
    gmsh_row = int(request.get("gmsh_sweep_row", len(current_rows) - 1))
    if gmsh_row < 0 or gmsh_row >= len(current_rows):
        raise ValueError("gmsh_sweep_row is outside current_rows_a")

    request_contract = {
        "schema": "radia.vol2d-postprocess-request.v1",
        "mesh_contract_sha256": operators["assembly"]["mesh_contract"]["contract_sha256"],
        "material_contract_sha256": material_contract["contract_sha256"],
        "operator_sha256": operators["operator_sha256"],
        "formulation": formulation,
        "element_family": prepared["element_family"],
        "coordinate_contract": coordinates,
        "branch_order": branch_order,
        "branch_materials": list(operators["assembly"]["branch_materials"]),
        "branch_turns": list(operators["assembly"]["branch_turns"]),
        "current_rows_a": current_rows,
        "point_probes": point_probes,
        "path_probes": path_probes,
        "region_materials": region_names,
        "model_depth_m": model_depth_m,
        "gmsh_sweep_row": gmsh_row,
        "export_basename": basename,
    }
    request_sha = _sha(request_contract)
    for key in ("expected_mesh_contract_sha256", "expected_material_contract_sha256", "expected_operator_sha256"):
        if request.get(key) is not None:
            actual_key = key.removeprefix("expected_")
            if str(request[key]) != str(request_contract[actual_key]):
                raise ValueError(f"{key} does not match the solved request")

    factor_started = time.perf_counter()
    stiffness = np.asarray(operators["assembly"]["field_matrix"], dtype=float)
    source = np.asarray(operators["assembly"]["source_matrix"], dtype=float)
    try:
        factor = np.linalg.cholesky(0.5 * (stiffness + stiffness.T))
    except np.linalg.LinAlgError as exc:
        raise ValueError("postprocess field operator is not positive definite") from exc
    factor_sha = _sha(
        {"operator_sha256": operators["operator_sha256"], "cholesky": factor.tolist()}
    )
    factor_s = time.perf_counter() - factor_started

    from ngsolve import GridFunction  # type: ignore

    gfu = GridFunction(fes)
    free = np.asarray(operators["assembly"]["free_dof_indices_0based"], dtype=int)
    field = _field_coefficient(gfu, formulation)
    states: list[list[float]] = []
    sweep_rows: list[dict[str, Any]] = []
    solve_s = 0.0
    post_s = 0.0
    for row_index, currents in enumerate(current_rows):
        solve_started = time.perf_counter()
        rhs = source @ np.asarray(currents, dtype=float)
        state = np.linalg.solve(factor.T, np.linalg.solve(factor, rhs))
        residual = stiffness @ state - rhs
        vector = gfu.vec.FV().NumPy()
        vector[:] = 0.0
        vector[free] = state
        solve_s += time.perf_counter() - solve_started

        post_started = time.perf_counter()
        probes = [
            {
                "name": probe["name"],
                **_sample_point(
                    mesh,
                    gfu,
                    field,
                    probe["coordinates_m"],
                    formulation=formulation,
                    material_names=list(permeability),
                    permeability=permeability,
                ),
            }
            for probe in point_probes
        ]
        paths = [
            _sample_path(
                mesh,
                gfu,
                field,
                path,
                formulation=formulation,
                material_names=list(permeability),
                permeability=permeability,
            )
            for path in path_probes
        ]
        regions = _region_integrals(
            mesh,
            field,
            region_names,
            permeability,
            formulation=formulation,
            model_depth_m=model_depth_m,
        )
        total_energy = float(sum(region["magnetic_energy_j"] for region in regions))
        maximum_probe = float(max(probe["b_magnitude_t"] for probe in probes))
        state_list = state.tolist()
        states.append(state_list)
        sweep_rows.append(
            {
                "row_index": row_index,
                "branch_current_a": currents,
                "field_state_sha256": _sha(state_list),
                "field_residual_inf": float(np.linalg.norm(residual, ord=np.inf)),
                "point_probes": probes,
                "path_probes": paths,
                "region_integrals": regions,
                "total_magnetic_energy_j": total_energy,
                "maximum_point_probe_b_t": maximum_probe,
            }
        )
        post_s += time.perf_counter() - post_started

    export_started = time.perf_counter()
    selected = np.asarray(states[gmsh_row], dtype=float)
    vector = gfu.vec.FV().NumPy()
    vector[:] = 0.0
    vector[free] = selected
    gmsh_field = _field_coefficient(gfu, formulation, axis_safe=True)
    gmsh_contents = _export_gmsh(
        mesh,
        gfu,
        gmsh_field,
        basename=basename,
        request_sha256=request_sha,
        potential_name=coordinates["potential_name"],
    )
    csv_content = _export_csv(branch_order, sweep_rows)
    export_s = time.perf_counter() - export_started

    state_table_sha = _sha(states)
    result_table_sha = _sha(sweep_rows)
    export_hashes = {name: _sha(content) for name, content in gmsh_contents.items()}
    export_hashes["csv"] = _sha(csv_content)
    result_contract = {
        "schema": POSTPROCESS_SCHEMA,
        "status": "solved",
        "operation": "solve",
        "execution_version": {
            "radia_mcp": _package_version(),
            "ngsolve": getattr(__import__("ngsolve"), "__version__", "unknown"),
        },
        "request_contract": request_contract,
        "request_contract_sha256": request_sha,
        "mesh_contract_sha256": request_contract["mesh_contract_sha256"],
        "material_contract_sha256": request_contract["material_contract_sha256"],
        "operator_sha256": request_contract["operator_sha256"],
        "factorization_sha256": factor_sha,
        "field_state_rows": states,
        "field_state_table_sha256": state_table_sha,
        "sweep_rows": sweep_rows,
        "result_table_sha256": result_table_sha,
        "mesh_rebuild_count": 0,
        "operator_rebuild_count": 0,
        "factorization_count": 1,
        "solve_count": len(current_rows),
        "gmsh_sweep_row": gmsh_row,
        "export_content_sha256": export_hashes,
        "generated_vol_git_required": False,
    }
    json_content = _canonical(result_contract)
    exports = {
        "json": _export_entry(json_content, f"{basename}.json", "application/json"),
        "csv": _export_entry(csv_content, f"{basename}.csv", "text/csv"),
        "gmsh_msh": _export_entry(gmsh_contents["gmsh_msh"], f"{basename}.msh", "model/mesh"),
        "gmsh_geo": _export_entry(gmsh_contents["gmsh_geo"], f"{basename}.geo", "text/plain"),
        "gmsh_geo_opt": _export_entry(gmsh_contents["gmsh_geo_opt"], f"{basename}.geo.opt", "text/plain"),
        "gmsh_msh_opt": _export_entry(gmsh_contents["gmsh_msh_opt"], f"{basename}.msh.opt", "text/plain"),
    }
    result_contract["canonical_json_sha256"] = exports["json"]["sha256"]
    return {
        "execution_version": result_contract["execution_version"],
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "schema": POSTPROCESS_SCHEMA,
        "status": "solved",
        "operation": "solve",
        "timing_s": {
            "assembly": assembly_s,
            "factorization": factor_s,
            "solve_rows": solve_s,
            "postprocess_and_export": post_s + export_s,
            "total": assembly_s + factor_s + solve_s + post_s + export_s,
        },
        "result_contract": result_contract,
        "exports": exports,
    }


def postprocess_replay_gate(artifact: Mapping[str, Any]) -> dict[str, Any]:
    """Recompute deterministic identities without trusting stored pass flags."""

    if not isinstance(artifact, Mapping):
        raise ValueError("replay_artifact must be an object")
    contract = artifact.get("result_contract")
    exports = artifact.get("exports")
    if not isinstance(contract, Mapping) or not isinstance(exports, Mapping):
        raise ValueError("replay artifact needs result_contract and exports objects")
    checks: dict[str, bool] = {
        "schema": contract.get("schema") == POSTPROCESS_SCHEMA,
        "status": contract.get("status") == "solved",
        "request_contract_sha256": _sha(contract.get("request_contract"))
        == contract.get("request_contract_sha256"),
        "field_state_table_sha256": _sha(contract.get("field_state_rows"))
        == contract.get("field_state_table_sha256"),
        "result_table_sha256": _sha(contract.get("sweep_rows"))
        == contract.get("result_table_sha256"),
        "single_factorization": contract.get("factorization_count") == 1,
        "no_mesh_rebuild": contract.get("mesh_rebuild_count") == 0,
        "no_operator_rebuild": contract.get("operator_rebuild_count") == 0,
        "solve_count_matches_rows": contract.get("solve_count")
        == len(contract.get("sweep_rows", [])),
    }
    content_hashes = contract.get("export_content_sha256")
    if not isinstance(content_hashes, Mapping):
        raise ValueError("result_contract.export_content_sha256 must be an object")
    for name in ("csv", "gmsh_msh", "gmsh_geo", "gmsh_geo_opt", "gmsh_msh_opt"):
        entry = exports.get(name)
        checks[f"{name}_content_sha256"] = (
            isinstance(entry, Mapping)
            and _sha(str(entry.get("content", ""))) == entry.get("sha256")
            and entry.get("sha256") == content_hashes.get(name)
        )
    json_entry = exports.get("json")
    stored_contract = dict(contract)
    canonical_sha = stored_contract.pop("canonical_json_sha256", None)
    expected_json = _canonical(stored_contract)
    checks["canonical_json_content"] = (
        isinstance(json_entry, Mapping)
        and json_entry.get("content") == expected_json
        and _sha(expected_json) == json_entry.get("sha256")
        and canonical_sha == json_entry.get("sha256")
    )
    checks["gmsh_v41"] = (
        isinstance(exports.get("gmsh_msh"), Mapping)
        and "$MeshFormat\n4.1 " in str(exports["gmsh_msh"].get("content", ""))
    )
    basename = str(contract.get("request_contract", {}).get("export_basename", ""))
    checks["geo_merges_msh"] = (
        isinstance(exports.get("gmsh_geo"), Mapping)
        and f'Merge "{basename}.msh";' in str(exports["gmsh_geo"].get("content", ""))
    )
    checks["geo_opt_exact_sidecar"] = (
        isinstance(exports.get("gmsh_geo_opt"), Mapping)
        and str(exports["gmsh_geo_opt"].get("filename", "")).endswith(".geo.opt")
    )
    passed = all(checks.values())
    return {
        "schema": REPLAY_SCHEMA,
        "status": "accepted" if passed else "rejected",
        "checks": checks,
        "pass": passed,
        "request_contract_sha256": contract.get("request_contract_sha256"),
        "field_state_table_sha256": contract.get("field_state_table_sha256"),
        "result_table_sha256": contract.get("result_table_sha256"),
    }


def analyze_vol2d_postprocess(request: Mapping[str, Any]) -> dict[str, Any]:
    """Dispatch the closed-world solve or replay operation."""

    if not isinstance(request, Mapping):
        raise ValueError("request must be an object")
    operation = str(request.get("operation", "solve"))
    if operation == "solve":
        return solve_vol2d_postprocess(request)
    if operation == "replay_gate":
        return postprocess_replay_gate(request.get("replay_artifact"))
    raise ValueError("operation must be solve or replay_gate")
