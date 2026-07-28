"""MCP-safe periodic-sector and AGE phase-rotation execution contracts."""

from __future__ import annotations

import hashlib
import json
import math
import time
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from typing import Any

import numpy as np

from .airgap_machine import (
    airgap_coupling,
    airgap_factorize,
    airgap_solve,
    airgap_torque,
)
from .solve import machine_symmetry_sector
from .vol2d_circuit import MU0, _runtime_vol_path, parse_netgen_2d_vol


SCHEMA = "radia.age-periodic-motion-analysis.v1"
SECTOR_SCHEMA = "radia.age-periodic-sector.v1"


def _package_version() -> str:
    try:
        return version("radia-mcp")
    except PackageNotFoundError:
        return "source-tree"


def _sha(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _finite(value: Any, label: str) -> float:
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


def _positive_int(value: Any, label: str, *, maximum: int | None = None) -> int:
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a positive integer") from exc
    if result <= 0 or float(result) != _finite(value, label):
        raise ValueError(f"{label} must be a positive integer")
    if maximum is not None and result > maximum:
        raise ValueError(f"{label} must not exceed {maximum}")
    return result


def _matrix_sha(matrix: Any, metadata: Mapping[str, Any]) -> str:
    rows, columns, values = matrix.COO()
    digest = hashlib.sha256(
        json.dumps(dict(metadata), sort_keys=True, separators=(",", ":")).encode(
            "utf-8"
        )
    )
    for array, dtype in (
        (np.asarray(rows), "<i8"),
        (np.asarray(columns), "<i8"),
        (np.asarray(values), "<c16"),
    ):
        normalized = np.ascontiguousarray(array, dtype=dtype)
        digest.update(str(normalized.shape).encode("ascii"))
        digest.update(normalized.tobytes())
    return digest.hexdigest()


def _angle_grid(raw: Any, period: float) -> dict[str, Any]:
    if (
        not isinstance(raw, Sequence)
        or isinstance(raw, (str, bytes))
        or len(raw) < 4
        or len(raw) > 361
    ):
        raise ValueError("rotor_angles_rad must contain 4..361 samples")
    angles = [_finite(value, "rotor_angles_rad") for value in raw]
    if any(right <= left for left, right in zip(angles, angles[1:])):
        raise ValueError("rotor_angles_rad must be strictly increasing")
    steps = np.diff(np.asarray(angles))
    spacing = float(steps[0])
    if not np.allclose(steps, spacing, rtol=1.0e-10, atol=1.0e-14):
        raise ValueError("rotor_angles_rad must use one uniform spacing")
    if not math.isclose(angles[0], 0.0, rel_tol=0.0, abs_tol=1.0e-14):
        raise ValueError("rotor angle grid must start at zero")
    if angles[-1] >= period - 1.0e-14:
        raise ValueError("rotor angle grid must exclude the repeated period endpoint")
    if not math.isclose(
        spacing * len(angles), period, rel_tol=1.0e-9, abs_tol=1.0e-12
    ):
        raise ValueError("rotor angle grid must uniformly cover exactly one torque period")
    result = {
        "angle_basis": "mechanical_rad",
        "endpoint_policy": "exclude_repeated_period_endpoint",
        "angles_rad": angles,
        "spacing_rad": spacing,
        "period_rad": period,
    }
    result["angle_grid_sha256"] = _sha(result)
    return result


def normalize_periodic_sector(raw: Any) -> dict[str, Any]:
    """Require one internally consistent slot/pole sector and phase convention."""

    if not isinstance(raw, Mapping):
        raise ValueError("periodic_sector must be an object")
    slots = _positive_int(raw.get("slots"), "slots")
    poles = _positive_int(raw.get("poles"), "poles")
    expected = machine_symmetry_sector(slots, poles)
    sectors = _positive_int(
        raw.get("sector_count", raw.get("sectors")), "sector_count"
    )
    angle = _positive(raw.get("sector_angle_deg"), "sector_angle_deg")
    boundary = str(raw.get("boundary", "")).strip().lower()
    phase = _finite(raw.get("boundary_phase"), "boundary_phase")
    expected_phase = -1.0 if expected["boundary"] == "anti-periodic" else 1.0
    if sectors != expected["sectors"]:
        raise ValueError("sector_count does not match gcd(slots, poles)")
    if not math.isclose(angle, expected["sector_angle_deg"], abs_tol=1.0e-12):
        raise ValueError("sector_angle_deg does not match the slot/pole symmetry")
    if boundary != expected["boundary"]:
        raise ValueError("periodic boundary sign does not match poles per sector")
    if phase != expected_phase:
        raise ValueError("boundary_phase must be +1 periodic or -1 anti-periodic")
    result = {
        "schema": SECTOR_SCHEMA,
        **expected,
        "boundary_phase": expected_phase,
        "angle_basis": "mechanical_deg",
        "whole_machine_multiplier": sectors,
    }
    result["periodicity_contract_sha256"] = _sha(result)
    return result


def _normalize_materials(mesh_contract: Mapping[str, Any], raw: Any) -> dict[str, Any]:
    names = set(mesh_contract["material_names"])
    if not isinstance(raw, Mapping) or set(map(str, raw)) != names:
        raise ValueError("materials must cover the .vol materials exactly")
    rows: dict[str, dict[str, float]] = {}
    for name in sorted(names):
        value = raw[name]
        if not isinstance(value, Mapping):
            raise ValueError(f"materials[{name}] must be an object")
        mu_r = _positive(value.get("relative_permeability", 1.0), "relative_permeability")
        sigma = _finite(value.get("conductivity_s_per_m", 0.0), "conductivity_s_per_m")
        if sigma < 0.0:
            raise ValueError("conductivity_s_per_m must be nonnegative")
        rows[name] = {
            "relative_permeability": mu_r,
            "conductivity_s_per_m": sigma,
        }
    result = {
        "schema": "radia.age-material-contract.v1",
        "mesh_contract_sha256": mesh_contract["contract_sha256"],
        "materials": rows,
    }
    result["material_contract_sha256"] = _sha(result)
    return result


def _boundary_radius(mesh_view: Any, name: str) -> list[float]:
    numbers = {
        number for number, boundary_name in mesh_view.boundary_names.items() if boundary_name == name
    }
    return [
        math.hypot(mesh_view.points[node - 1][0], mesh_view.points[node - 1][1])
        for edge in mesh_view.boundary_edges
        if edge.boundary_number in numbers
        for node in edge.nodes
    ]


def _component_contract(mesh_view: Any, names: Mapping[str, str]) -> dict[str, Any]:
    node_to_cells: dict[int, set[int]] = {}
    for cell_index, cell in enumerate(mesh_view.cells):
        for node in cell.nodes:
            node_to_cells.setdefault(node, set()).add(cell_index)
    unseen = set(range(len(mesh_view.cells)))
    components: list[set[int]] = []
    while unseen:
        first = min(unseen)
        unseen.remove(first)
        pending = [first]
        component: set[int] = set()
        while pending:
            cell_index = pending.pop()
            component.add(cell_index)
            neighbours = {
                other
                for node in mesh_view.cells[cell_index].nodes
                for other in node_to_cells[node]
            }
            new = neighbours & unseen
            unseen -= new
            pending.extend(sorted(new, reverse=True))
        components.append(component)
    if len(components) != 2:
        raise ValueError("AGE .vol must contain exactly two disconnected FE regions")

    cell_component = {
        cell_index: component_index
        for component_index, component in enumerate(components)
        for cell_index in component
    }
    node_components = {
        node: {cell_component[cell_index] for cell_index in cell_indices}
        for node, cell_indices in node_to_cells.items()
    }

    def boundary_component(boundary_name: str) -> int:
        numbers = {
            number
            for number, value in mesh_view.boundary_names.items()
            if value == boundary_name
        }
        memberships = {
            component
            for edge in mesh_view.boundary_edges
            if edge.boundary_number in numbers
            for node in edge.nodes
            for component in node_components.get(node, set())
        }
        if len(memberships) != 1:
            raise ValueError(f"{boundary_name} must belong to exactly one FE region")
        return memberships.pop()

    rotor_component = boundary_component(names["rotor_ring"])
    stator_component = boundary_component(names["stator_ring"])
    if rotor_component == stator_component:
        raise ValueError("rotor_ring and stator_ring must lie on disconnected FE regions")

    material_components: dict[str, set[int]] = {}
    for cell_index, cell in enumerate(mesh_view.cells):
        material_components.setdefault(mesh_view.material_name(cell.material_number), set()).add(
            cell_component[cell_index]
        )
    for role, component in (
        ("rotor_material", rotor_component),
        ("stator_material", stator_component),
    ):
        material = names[role]
        if material_components.get(material) != {component}:
            raise ValueError(f"{role} must be confined to its named FE region")
    return {
        "component_count": 2,
        "rotor_component": rotor_component,
        "stator_component": stator_component,
        "component_cell_counts": [len(component) for component in components],
        "topology": "two_disconnected_fe_regions_coupled_by_unmeshed_age",
    }


def _normalize_gap(mesh_view: Any, raw: Any) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise ValueError("airgap must be an object")
    ri = _positive(raw.get("inner_radius_m"), "inner_radius_m")
    ro = _positive(raw.get("outer_radius_m"), "outer_radius_m")
    if ro <= ri:
        raise ValueError("outer_radius_m must exceed inner_radius_m")
    names = {
        key: str(raw.get(key, "")).strip()
        for key in (
            "rotor_ring",
            "stator_ring",
            "rotor_inner",
            "outer",
            "rotor_material",
            "stator_material",
        )
    }
    boundary_names = {names[key] for key in ("rotor_ring", "stator_ring", "rotor_inner", "outer")}
    if any(not name for name in names.values()) or len(boundary_names) != 4:
        raise ValueError("airgap names must identify four boundaries and two materials")
    if names["rotor_material"] == names["stator_material"]:
        raise ValueError("rotor_material and stator_material must differ")
    mesh_contract = mesh_view.contract()
    missing = sorted(boundary_names - set(mesh_contract["boundary_names"]))
    if missing:
        raise ValueError(f"airgap boundaries are absent from .vol: {missing}")
    missing_materials = sorted(
        {names["rotor_material"], names["stator_material"]}
        - set(mesh_contract["material_names"])
    )
    if missing_materials:
        raise ValueError(f"airgap materials are absent from .vol: {missing_materials}")
    for key, radius in (("rotor_ring", ri), ("stator_ring", ro)):
        values = _boundary_radius(mesh_view, names[key])
        if not values:
            raise ValueError(f"{key} has no boundary edges")
        relative = max(abs(value - radius) for value in values) / radius
        if relative > 1.0e-8:
            raise ValueError(f"{key} nodes do not lie on the declared radius")
    harmonics_raw = raw.get("harmonics")
    if not isinstance(harmonics_raw, Sequence) or isinstance(harmonics_raw, (str, bytes)):
        raise ValueError("harmonics must be a non-empty sequence")
    harmonics = sorted({_positive_int(value, "harmonic", maximum=64) for value in harmonics_raw})
    if not harmonics or len(harmonics) != len(harmonics_raw) or len(harmonics) > 16:
        raise ValueError("harmonics must contain 1..16 unique positive orders")
    result = {
        "schema": "radia.age-gap-contract.v1",
        "center_m": [0.0, 0.0],
        "inner_radius_m": ri,
        "outer_radius_m": ro,
        **names,
        "component_contract": _component_contract(mesh_view, names),
        "harmonics": harmonics,
    }
    result["gap_contract_sha256"] = _sha(result)
    return result


def _normalize_excitation(harmonics: Sequence[int], raw: Any) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise ValueError("excitation must be an object keyed by harmonic")
    if set(raw) != {str(value) for value in harmonics}:
        raise ValueError("excitation keys must match retained harmonics exactly")
    rows = {}
    for harmonic in harmonics:
        row = raw[str(harmonic)]
        if not isinstance(row, Mapping):
            raise ValueError("each excitation harmonic must be an object")
        rows[str(harmonic)] = {
            "rotor_amplitude": _finite(row.get("rotor_amplitude"), "rotor_amplitude"),
            "stator_amplitude": _finite(row.get("stator_amplitude"), "stator_amplitude"),
        }
    if not any(
        abs(row["rotor_amplitude"] * row["stator_amplitude"]) > 0.0
        for row in rows.values()
    ):
        raise ValueError("at least one harmonic needs nonzero rotor and stator amplitudes")
    result = {
        "schema": "radia.age-excitation-contract.v1",
        "rotation": "rotor_harmonic_phase_only",
        "harmonics": rows,
    }
    result["excitation_sha256"] = _sha(result)
    return result


def _harmonic_cf(harmonics: Sequence[int], excitation: Mapping[str, Any], angle: float, *, rotor: bool) -> Any:
    from ngsolve import atan2, cos, sin, x, y  # type: ignore

    theta = atan2(y, x)
    result = 0.0
    for harmonic in harmonics:
        row = excitation["harmonics"][str(harmonic)]
        amplitude = row["rotor_amplitude" if rotor else "stator_amplitude"]
        phase = harmonic * angle if rotor else 0.0
        result = result + amplitude * (
            math.cos(phase) * cos(harmonic * theta)
            + math.sin(phase) * sin(harmonic * theta)
        )
    return result


def solve_age_periodic_motion(request: Mapping[str, Any]) -> dict[str, Any]:
    """Run one no-remesh AGE angle sweep from a dimension-2 ``.vol``."""

    started = time.perf_counter()
    mesh_view = parse_netgen_2d_vol(
        request.get("vol_text"), source_name=str(request.get("source_name", "generated.vol"))
    )
    mesh_contract = mesh_view.contract()
    gap = _normalize_gap(mesh_view, request.get("airgap"))
    materials = _normalize_materials(mesh_contract, request.get("materials"))
    periodicity = normalize_periodic_sector(request.get("periodic_sector"))
    excitation = _normalize_excitation(gap["harmonics"], request.get("excitation"))
    harmonic_gcd = math.gcd(*gap["harmonics"])
    period = 2.0 * math.pi / harmonic_gcd
    angle_grid = _angle_grid(request.get("rotor_angles_rad"), period)
    axial_length = _positive(request.get("axial_length_m"), "axial_length_m")
    frequency = _finite(request.get("frequency_hz", 0.0), "frequency_hz")
    if frequency < 0.0:
        raise ValueError("frequency_hz must be nonnegative")
    order = _positive_int(request.get("element_order", 3), "element_order", maximum=6)
    prepare_s = time.perf_counter() - started

    from ngsolve import BilinearForm, H1, Mesh, TaskManager, dx, grad  # type: ignore

    mesh = Mesh(str(_runtime_vol_path(request["vol_text"], mesh_view.content_sha256)))
    complex_space = frequency > 0.0
    dirichlet = f"{gap['rotor_inner']}|{gap['outer']}"
    fes = H1(mesh, order=order, complex=complex_space, dirichlet=dirichlet)
    trial, test = fes.TnT()
    inv_mu_r = mesh.MaterialCF(
        {name: 1.0 / row["relative_permeability"] for name, row in materials["materials"].items()}
    )
    sigma = mesh.MaterialCF(
        {name: row["conductivity_s_per_m"] for name, row in materials["materials"].items()}
    )
    form = BilinearForm(fes)
    form += inv_mu_r * grad(trial) * grad(test) * dx
    if frequency > 0.0:
        form += 1j * (2.0 * math.pi * frequency) * MU0 * sigma * trial * test * dx
    assemble_started = time.perf_counter()
    with TaskManager():
        form.Assemble()
    operator_sha = _matrix_sha(
        form.mat,
        {
            "mesh": mesh_contract["contract_sha256"],
            "material": materials["material_contract_sha256"],
            "frequency_hz": frequency,
            "order": order,
        },
    )
    coupling = airgap_coupling(
        fes,
        gap["inner_radius_m"],
        gap["outer_radius_m"],
        gap["rotor_ring"],
        gap["stator_ring"],
        gap["harmonics"],
    )
    factor = airgap_factorize(form.mat, coupling, fes.FreeDofs())
    factor_sha = _sha(
        {
            "operator_sha256": operator_sha,
            "gap_contract_sha256": gap["gap_contract_sha256"],
            "small_inverse": [
                [[float(value.real), float(value.imag)] for value in row]
                for row in np.asarray(factor["smallinv"], dtype=complex)
            ],
        }
    )
    assemble_s = time.perf_counter() - assemble_started

    solve_started = time.perf_counter()
    rows = []
    with TaskManager():
        for angle in angle_grid["angles_rad"]:
            boundary_cf = mesh.BoundaryCF(
                {
                    gap["rotor_inner"]: _harmonic_cf(
                        gap["harmonics"], excitation, angle, rotor=True
                    ),
                    gap["outer"]: _harmonic_cf(
                        gap["harmonics"], excitation, angle, rotor=False
                    ),
                }
            )
            field = airgap_solve(factor, fes, dirichlet_cf=boundary_cf)
            torque = float(airgap_torque(coupling, field, axial_length=axial_length))
            vector = np.ascontiguousarray(field.vec.FV().NumPy(), dtype="<c16")
            rows.append(
                {
                    "rotor_angle_rad": angle,
                    "torque_nm": torque,
                    "field_state_sha256": hashlib.sha256(vector.tobytes()).hexdigest(),
                    "mesh_contract_sha256": mesh_contract["contract_sha256"],
                    "material_contract_sha256": materials["material_contract_sha256"],
                    "operator_sha256": operator_sha,
                    "age_factorization_sha256": factor_sha,
                    "excitation_sha256": excitation["excitation_sha256"],
                }
            )
        closure_cf = mesh.BoundaryCF(
            {
                gap["rotor_inner"]: _harmonic_cf(
                    gap["harmonics"], excitation, period, rotor=True
                ),
                gap["outer"]: _harmonic_cf(
                    gap["harmonics"], excitation, period, rotor=False
                ),
            }
        )
        closure_field = airgap_solve(factor, fes, dirichlet_cf=closure_cf)
        closure_torque = float(
            airgap_torque(coupling, closure_field, axial_length=axial_length)
        )
    solve_s = time.perf_counter() - solve_started

    post_started = time.perf_counter()
    torques = np.asarray([row["torque_nm"] for row in rows])
    scale = max(float(np.max(np.abs(torques))), 1.0e-30)
    variation = float(np.ptp(torques))
    closure_relative = float(abs(closure_torque - torques[0]) / scale)
    sign_reversal = bool(torques.min() < 0.0 < torques.max())
    if scale <= 1.0e-12:
        raise ValueError("AGE torque sweep is indeterminate because torque is near zero")
    if variation <= 1.0e-12 * scale:
        raise ValueError("AGE torque sweep is indeterminate because torque is constant")
    if frequency == 0.0 and not sign_reversal:
        raise ValueError("static AGE torque sweep must demonstrate phase-sign reversal")
    if closure_relative > 1.0e-8:
        raise ValueError("AGE torque does not close over the declared mechanical period")
    output_sha = _sha(rows)
    post_s = time.perf_counter() - post_started
    return {
        "execution_version": {
            "radia_mcp": _package_version(),
            "ngsolve": getattr(__import__("ngsolve"), "__version__", "unknown"),
        },
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "schema": SCHEMA,
        "status": "solved",
        "operation": "solve",
        "mesh_contract": mesh_contract,
        "material_contract": materials,
        "gap_contract": gap,
        "periodicity_contract": periodicity,
        "excitation_contract": excitation,
        "angle_grid": angle_grid,
        "frequency_hz": frequency,
        "axial_length_m": axial_length,
        "operator_sha256": operator_sha,
        "age_factorization_sha256": factor_sha,
        "angle_grid_sha256": angle_grid["angle_grid_sha256"],
        "excitation_sha256": excitation["excitation_sha256"],
        "periodicity_contract_sha256": periodicity["periodicity_contract_sha256"],
        "torque_output_sha256": output_sha,
        "mesh_reused_all_angles": True,
        "operator_reused_all_angles": True,
        "factorization_reused_all_angles": True,
        "rotation_method": "rotor_harmonic_phase_only_no_remesh",
        "torque_rows": rows,
        "torque_summary": {
            "minimum_nm": float(torques.min()),
            "maximum_nm": float(torques.max()),
            "peak_to_peak_nm": variation,
            "closure_relative_error": closure_relative,
            "phase_sign_reversal_observed": sign_reversal,
            "phase_sign_reversal_required": frequency == 0.0,
            "full_machine_scope": True,
        },
        "result_output_schema_id": "radia.age-torque-angle-table.v1",
        "result_output_columns": ["rotor_angle_rad", "torque_nm"],
        "result_output_units": {"rotor_angle_rad": "rad", "torque_nm": "N_m"},
        "timing_breakdown_s": {
            "prepare": prepare_s,
            "assemble_and_factor": assemble_s,
            "angle_sweep_and_postprocess": solve_s + post_s,
            "total": time.perf_counter() - started,
        },
    }


def age_sector_torque_gate(request: Mapping[str, Any]) -> dict[str, Any]:
    """Validate one periodic-sector torque table and whole-machine scaling."""

    periodicity = normalize_periodic_sector(request.get("periodic_sector"))
    rows = request.get("rows")
    if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes)) or len(rows) < 4:
        raise ValueError("sector torque rows must contain at least four samples")
    parsed = []
    identities: dict[str, set[str]] = {
        key: set()
        for key in (
            "mesh_contract_sha256",
            "material_contract_sha256",
            "operator_sha256",
            "age_factorization_sha256",
            "excitation_sha256",
        )
    }
    for index, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise ValueError(f"rows[{index}] must be an object")
        angle = _finite(row.get("rotor_angle_rad"), "rotor_angle_rad")
        sector_torque = _finite(row.get("sector_torque_nm"), "sector_torque_nm")
        full_torque = _finite(row.get("full_machine_torque_nm"), "full_machine_torque_nm")
        if not math.isclose(
            full_torque,
            periodicity["whole_machine_multiplier"] * sector_torque,
            rel_tol=1.0e-10,
            abs_tol=1.0e-12,
        ):
            raise ValueError("full-machine torque must equal sector multiplier times sector torque")
        for key in identities:
            value = str(row.get(key, ""))
            if not value:
                raise ValueError(f"rows[{index}] requires {key}")
            identities[key].add(value)
        parsed.append({"rotor_angle_rad": angle, "sector_torque_nm": sector_torque, "full_machine_torque_nm": full_torque})
    if any(len(values) != 1 for values in identities.values()):
        raise ValueError("sector torque rows must reuse mesh/material/operator/factor/excitation identities")
    angles = [row["rotor_angle_rad"] for row in parsed]
    if any(right <= left for left, right in zip(angles, angles[1:])):
        raise ValueError("sector torque angles must be strictly increasing")
    result = {
        "schema": "radia.age-sector-torque-gate.v1",
        "status": "ok",
        "periodicity_contract": periodicity,
        "rows": parsed,
        "identity_reused_without_remesh": True,
        "result_output_schema_id": "radia.age-sector-torque-table.v1",
        "result_output_units": {"angle": "rad", "torque": "N_m"},
    }
    result["result_sha256"] = _sha(result)
    return result


def analyze_age_periodic_motion(request: Mapping[str, Any]) -> dict[str, Any]:
    operation = str(request.get("operation", "solve"))
    if operation == "solve":
        return solve_age_periodic_motion(request)
    if operation == "periodic_sector_gate":
        return age_sector_torque_gate(request)
    raise ValueError("operation must be solve or periodic_sector_gate")
