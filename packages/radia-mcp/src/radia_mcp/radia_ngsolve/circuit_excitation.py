"""Typed circuit and AGE application contracts for 2D magnetic solvers.

The finite-element formulation should not know whether a user entered a
current through a GUI or through Python.  This module turns that application
intent into explicit source and constraint records shared by planar 2D and
axisymmetric formulations.

Series and parallel circuits are deliberately different:

* a series circuit prescribes the same branch current in every assigned region;
* a parallel circuit keeps branch currents unknown, shares one terminal-voltage
  unknown, and adds a total-current constraint.

The module also compiles rotary and linear Air-Gap Element (AGE) motion into
Fourier phase factors.  Neither AGE path rebuilds the rotor/stator FE meshes.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

from .airgap_element import annular_rotation_phase, planar_translation_phase


SCHEMA = "radia.circuit-age-application.v1"
_PROBLEM_KINDS = {"planar_2d", "axisymmetric"}
_ELEMENT_FAMILIES = {"P1", "P2", "Q1", "Q2", "P2_curved", "Q2_curved"}
_CONNECTIONS = {"series", "parallel"}


def _finite(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be a finite number")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a finite number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} must be a finite number")
    return result


def _phasor(value: Any, label: str) -> complex:
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        result = complex(_finite(value, label), 0.0)
    elif (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and len(value) == 2
    ):
        result = complex(
            _finite(value[0], f"{label}[0]"),
            _finite(value[1], f"{label}[1]"),
        )
    else:
        raise ValueError(f"{label} must be a real value or [real, imag]")
    return result


def _pair(value: complex) -> list[float]:
    return [float(value.real), float(value.imag)]


def _name(value: Any, label: str) -> str:
    result = str(value or "").strip()
    if not result:
        raise ValueError(f"{label} must be a non-empty string")
    return result


def _compile_motion(motion: Any, problem_kind: str) -> dict[str, Any] | None:
    if motion in (None, {}):
        return None
    if not isinstance(motion, Mapping):
        raise ValueError("motion must be an object")
    kind = str(motion.get("kind", "")).strip().lower()
    if kind not in {"annular_age", "planar_age"}:
        raise ValueError("motion.kind must be annular_age or planar_age")
    if problem_kind != "planar_2d":
        raise ValueError("AGE motion belongs to planar_2d cross-sections")

    if kind == "annular_age":
        position = _finite(motion.get("position_rad", 0.0), "motion.position_rad")
        modes = motion.get("harmonics", [])
        if not isinstance(modes, Sequence) or isinstance(modes, (str, bytes)):
            raise ValueError("motion.harmonics must be a sequence")
        normalized: list[int] = []
        for value in modes:
            if isinstance(value, bool) or int(value) != value or int(value) <= 0:
                raise ValueError("annular AGE harmonics must be positive integers")
            normalized.append(int(value))
        if not normalized or len(set(normalized)) != len(normalized):
            raise ValueError("annular AGE harmonics must be non-empty and unique")
        factors = {
            str(mode): _pair(annular_rotation_phase(mode, position))
            for mode in normalized
        }
        return {
            "kind": kind,
            "position_rad": position,
            "trace_phase_convention": "A_rotor,n(theta-delta) = exp(-i*n*delta) A_rotor,n(theta)",
            "phase_factors": factors,
            "coupling_operator": "annular Fourier DtN (AGE)",
            "mechanical_observable": "mesh-independent harmonic torque",
            "mesh_rebuild": False,
        }

    position = _finite(motion.get("position_m", 0.0), "motion.position_m")
    modes = motion.get("wavenumbers_per_m", [])
    if not isinstance(modes, Sequence) or isinstance(modes, (str, bytes)):
        raise ValueError("motion.wavenumbers_per_m must be a sequence")
    normalized_float = [_finite(value, "AGE wavenumber") for value in modes]
    if (
        not normalized_float
        or any(value <= 0.0 for value in normalized_float)
        or len(set(normalized_float)) != len(normalized_float)
    ):
        raise ValueError("planar AGE wavenumbers must be positive and unique")
    factors = {
        format(wavenumber, ".17g"): _pair(
            planar_translation_phase(wavenumber, position)
        )
        for wavenumber in normalized_float
    }
    return {
        "kind": kind,
        "position_m": position,
        "trace_phase_convention": "A_moving,k(x-d) = exp(-i*k*d) A_moving,k(x)",
        "phase_factors": factors,
        "coupling_operator": "planar Fourier DtN (AGE)",
        "mechanical_observable": "mesh-independent harmonic thrust",
        "mesh_rebuild": False,
    }


def compile_circuit_age_application(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Compile application-level current and AGE intent into solver contracts.

    Complex quantities use ``[real, imag]`` so the returned object is directly
    JSON serializable.  The compiler does not silently divide a prescribed
    parallel current equally between regions; the field-circuit system owns
    that current split.
    """

    if not isinstance(payload, Mapping):
        raise ValueError("application payload must be an object")
    problem_kind = str(payload.get("problem_kind", "")).strip()
    if problem_kind not in _PROBLEM_KINDS:
        raise ValueError("problem_kind must be planar_2d or axisymmetric")
    element_family = str(payload.get("element_family", "")).strip()
    if element_family not in _ELEMENT_FAMILIES:
        raise ValueError(
            "element_family must be P1, P2, Q1, Q2, P2_curved, or Q2_curved"
        )

    circuit_rows = payload.get("circuits", [])
    region_rows = payload.get("regions", [])
    if not isinstance(circuit_rows, Sequence) or isinstance(circuit_rows, (str, bytes)):
        raise ValueError("circuits must be a sequence")
    if not isinstance(region_rows, Sequence) or isinstance(region_rows, (str, bytes)):
        raise ValueError("regions must be a sequence")
    if not circuit_rows or not region_rows:
        raise ValueError("at least one circuit and one assigned region are required")

    circuits: dict[str, dict[str, Any]] = {}
    for index, row in enumerate(circuit_rows):
        if not isinstance(row, Mapping):
            raise ValueError(f"circuits[{index}] must be an object")
        name = _name(row.get("name"), f"circuits[{index}].name")
        if name in circuits:
            raise ValueError(f"duplicate circuit name: {name}")
        connection = str(row.get("connection", "")).strip().lower()
        if connection not in _CONNECTIONS:
            raise ValueError(f"circuit {name} connection must be series or parallel")
        current = _phasor(row.get("current_a", 0.0), f"circuit {name} current_a")
        frequency_hz = _finite(
            row.get("frequency_hz", 0.0), f"circuit {name} frequency_hz"
        )
        if frequency_hz < 0.0:
            raise ValueError(f"circuit {name} frequency_hz must be nonnegative")
        circuits[name] = {
            "name": name,
            "connection": connection,
            "current_a": current,
            "frequency_hz": frequency_hz,
            "regions": [],
        }

    compiled_regions: list[dict[str, Any]] = []
    region_names: set[str] = set()
    for index, row in enumerate(region_rows):
        if not isinstance(row, Mapping):
            raise ValueError(f"regions[{index}] must be an object")
        name = _name(row.get("name"), f"regions[{index}].name")
        if name in region_names:
            raise ValueError(f"duplicate region name: {name}")
        region_names.add(name)
        circuit_name = _name(row.get("circuit"), f"region {name} circuit")
        if circuit_name not in circuits:
            raise ValueError(f"region {name} references undefined circuit {circuit_name}")
        turns_raw = row.get("turns", 1)
        if isinstance(turns_raw, bool) or int(turns_raw) != turns_raw or int(turns_raw) == 0:
            raise ValueError(f"region {name} turns must be a non-zero integer")
        turns = int(turns_raw)
        area = _finite(row.get("area_m2"), f"region {name} area_m2")
        if area <= 0.0:
            raise ValueError(f"region {name} area_m2 must be positive")
        conductivity = _finite(
            row.get("conductivity_s_per_m", 0.0),
            f"region {name} conductivity_s_per_m",
        )
        if conductivity < 0.0:
            raise ValueError(f"region {name} conductivity_s_per_m must be nonnegative")
        circuit = circuits[circuit_name]
        circuit["regions"].append(name)
        compiled = {
            "name": name,
            "circuit": circuit_name,
            "turns": turns,
            "area_m2": area,
            "conductivity_s_per_m": conductivity,
            "field_source_coefficient_per_m2": turns / area,
        }
        if circuit["connection"] == "series":
            compiled.update(
                {
                    "branch_current_status": "prescribed",
                    "branch_current_a": _pair(circuit["current_a"]),
                    "impressed_current_density_a_per_m2": _pair(
                        turns * circuit["current_a"] / area
                    ),
                }
            )
        else:
            compiled.update(
                {
                    "branch_current_status": "field_circuit_unknown",
                    "branch_current_unknown": f"I_branch:{name}",
                    "impressed_current_density_a_per_m2": None,
                }
            )
        compiled_regions.append(compiled)

    compiled_circuits: list[dict[str, Any]] = []
    for circuit in circuits.values():
        if not circuit["regions"]:
            raise ValueError(f"circuit {circuit['name']} has no assigned regions")
        base = {
            "name": circuit["name"],
            "connection": circuit["connection"],
            "prescribed_total_current_a": _pair(circuit["current_a"]),
            "frequency_hz": circuit["frequency_hz"],
            "assigned_regions": list(circuit["regions"]),
            "observables": [
                "current_a",
                "terminal_voltage_v",
                "flux_linkage_wb_turn",
                "joule_loss_w",
            ],
        }
        if circuit["connection"] == "series":
            base.update(
                {
                    "current_constraint": "I_branch(region) = I_circuit for every assigned region",
                    "voltage_reduction": "V_circuit = sum(V_region) with signed turns",
                    "additional_circuit_unknowns": [],
                }
            )
        else:
            branch_unknowns = [f"I_branch:{name}" for name in circuit["regions"]]
            base.update(
                {
                    "current_constraint": "sum(I_branch(region)) = I_circuit",
                    "voltage_reduction": "one common terminal-voltage/voltage-gradient unknown",
                    "additional_circuit_unknowns": [
                        f"V_common:{circuit['name']}",
                        *branch_unknowns,
                    ],
                    "equal_current_split_assumed": False,
                }
            )
        compiled_circuits.append(base)

    motion = _compile_motion(payload.get("motion"), problem_kind)
    return {
        "schema": SCHEMA,
        "status": "compiled",
        "problem_kind": problem_kind,
        "element_family": element_family,
        "element_order_contract": "circuit topology is independent of P1/P2/Q1/Q2 interpolation order",
        "geometry_contract": (
            "curved geometry is retained for P2_curved/Q2_curved"
            if element_family.endswith("_curved")
            else "straight-sided geometry"
        ),
        "circuits": compiled_circuits,
        "regions": compiled_regions,
        "motion": motion,
        "assembly_contract": {
            "magnetic_source": "J_source = source_coefficient * branch_current",
            "transient_voltage": "v = R i + d(lambda)/dt",
            "parallel_closure": "augment the field system; never pre-divide total current",
            "kelvin_open_boundary_compatible": True,
            "eddy_transient_compatible": True,
        },
    }
