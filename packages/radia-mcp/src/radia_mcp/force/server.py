"""MCP Server: radia_mcp.force — shared electromagnetic-force layer.

Usage:
    mcp-server-force
    mcp-server-force --selftest

The server imports without Radia, NumPy, or NGSolve.  Numerical sample
integration lazily loads ``radia.force`` when one of the calculation tools is
called; the knowledge and discovery tools remain standalone.
"""

from __future__ import annotations

import json
import sys
from typing import Any

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool, register_topics_tool
from ..common.tool_group import CoarseToolRegistry
from ..common.lazy_call import lazy_callable
_method_selection_gate = lazy_callable(".gates", "electromagnetic_force_method_selection_gate", __package__)
_action_reaction_gate = lazy_callable(".gates", "force_action_reaction_gate", __package__)
_method_agreement_gate = lazy_callable(".gates", "force_torque_method_agreement_gate", __package__)
_weight_equilibrium_gate = lazy_callable(".gates", "force_weight_equilibrium_gate", __package__)
from .knowledge import (
    TOPICS,
    get_force_extras,
    get_force_knowledge,
    get_force_methods,
    get_force_recipe,
    get_force_validation,
)

mcp = FastMCP("mcp-server-force")
_validation = CoarseToolRegistry(mcp, namespace="force")


def _load_force_api():
    try:
        from radia import force as force_api
    except (ImportError, OSError) as exc:
        raise RuntimeError(
            "Numerical Force tools require the Radia package. Install "
            "'radia-mcp[radia]' or 'radia'; knowledge tools remain available."
        ) from exc
    return force_api


def _json_result(payload: dict[str, Any]) -> str:
    return json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True)


def _phasor_rows(real_rows, imag_rows, name: str):
    if len(real_rows) != len(imag_rows):
        raise ValueError(f"{name} real and imaginary tables must have the same length")
    rows = []
    for index, (real_row, imag_row) in enumerate(zip(real_rows, imag_rows)):
        if len(real_row) != 3 or len(imag_row) != 3:
            raise ValueError(f"{name}[{index}] real and imaginary rows must have length 3")
        rows.append([
            complex(float(real_value), float(imag_value))
            for real_value, imag_value in zip(real_row, imag_row)
        ])
    return rows


def _phasor_values(real_values, imag_values, name: str):
    if len(real_values) != len(imag_values):
        raise ValueError(f"{name} real and imaginary tables must have the same length")
    return [
        complex(float(real_value), float(imag_value))
        for real_value, imag_value in zip(real_values, imag_values)
    ]


def _force_result(
    force_api,
    force_n,
    torque_nm,
    *,
    method: str,
    frame: str,
    pivot_m,
    field_convention: str = "static",
    amplitude: str | None = None,
    sample_count: int | None = None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    payload = force_api.force_torque_result(
        force_n,
        torque_nm,
        method=method,
        frame=frame,
        pivot_m=pivot_m,
        field_convention=field_convention,
        amplitude=amplitude,
    )
    payload["status"] = "ok"
    if sample_count is not None:
        payload["sample_count"] = int(sample_count)
    if extra:
        payload.update(extra)
    return payload


def _error_payload(method: str, status: str, exc: Exception) -> dict[str, Any]:
    return {
        "schema": "radia.force-result/v1",
        "status": status,
        "method": method,
        "error": str(exc),
    }


@mcp.tool()
def force(topic: str = "overview") -> str:
    """Common EM-force overview and routing; call `force_topics` for topics."""

    return get_force_knowledge(topic)


@mcp.tool()
def force_methods(section: str = "all") -> str:
    """Unified theory: Maxwell stress, eggshell, Arkkio, nodal, and Lorentz."""

    return get_force_methods(section)


@mcp.tool()
def force_recipe(section: str = "method_choice") -> str:
    """Practical method choice, high-order setup, pitfalls, and examples."""

    return get_force_recipe(section)


@mcp.tool()
def force_extras(section: str = "all") -> str:
    """PM, energy/coenergy, shape-derivative, Lorentz, and Meissner guidance."""

    return get_force_extras(section)


@mcp.tool()
def force_validation_guide(section: str = "all") -> str:
    """Force method map, cross-checks, eggshell guidance, and evidence contract."""

    return get_force_validation(section)


@mcp.tool()
def force_result(
    force_n: list[float] | None,
    torque_nm: list[float] | None,
    method: str,
    frame: str = "global_cartesian",
    pivot_m: list[float] | None = None,
    field_convention: str = "static",
    amplitude: str | None = None,
    dimensionality: str = "3d",
    per_unit_depth: bool = False,
) -> str:
    """Normalize solver-owned resultants to ``radia.force-result/v1``."""

    try:
        force_api = _load_force_api()
        payload = force_api.force_torque_result(
            force_n,
            torque_nm,
            method=method,
            frame=frame,
            pivot_m=pivot_m,
            field_convention=field_convention,
            amplitude=amplitude,
            dimensionality=dimensionality,
            per_unit_depth=per_unit_depth,
        )
        payload["status"] = "ok"
    except RuntimeError as exc:
        payload = _error_payload(method, "unavailable", exc)
    except (TypeError, ValueError) as exc:
        payload = _error_payload(method, "invalid_input", exc)
    return _json_result(payload)


@mcp.tool()
def force_lorentz(
    current_density_a_per_m2: list[list[float]],
    magnetic_flux_density_t: list[list[float]],
    volume_weights_m3: list[float],
    sample_points_m: list[list[float]] | None = None,
    pivot_m: list[float] | None = None,
    frame: str = "global_cartesian",
) -> str:
    """Integrate static Lorentz force and optional torque from SI samples.

    Each J and B row is `[x, y, z]`; each nonnegative volume weight includes
    the physical Jacobian and any symmetry-sector factor.  Supply sample points
    to also return torque in N m about `pivot_m`.
    """

    try:
        force_api = _load_force_api()
        if sample_points_m is None:
            force_n = force_api.integrate_lorentz_force(
                current_density_a_per_m2,
                magnetic_flux_density_t,
                volume_weights_m3,
            )
            torque_nm = None
            result_pivot = None
        else:
            force_n, torque_nm = force_api.integrate_lorentz_force_and_torque(
                current_density_a_per_m2,
                magnetic_flux_density_t,
                volume_weights_m3,
                sample_points_m,
                pivot_m=pivot_m,
            )
            result_pivot = pivot_m or [0.0, 0.0, 0.0]
        payload = _force_result(
            force_api,
            force_n,
            torque_nm,
            method="lorentz_body_force",
            frame=frame,
            pivot_m=result_pivot,
            sample_count=len(current_density_a_per_m2),
        )
    except RuntimeError as exc:
        payload = _error_payload("lorentz_body_force", "unavailable", exc)
    except (TypeError, ValueError) as exc:
        payload = _error_payload("lorentz_body_force", "invalid_input", exc)
    return _json_result(payload)


@mcp.tool()
def force_maxwell_surface(
    magnetic_flux_density_t: list[list[float]],
    outward_normals: list[list[float]],
    area_weights_m2: list[float],
    permeability_h_per_m: float = 1.2566370614359173e-6,
    sample_points_m: list[list[float]] | None = None,
    pivot_m: list[float] | None = None,
    frame: str = "global_cartesian",
) -> str:
    """Integrate static Maxwell surface force and optional torque.

    Each B and normal row is `[x, y, z]`.  Normals point outward from the body;
    area weights include physical surface Jacobians.  Supply sample points to
    also return torque in N m about `pivot_m`.
    """

    try:
        force_api = _load_force_api()
        if sample_points_m is None:
            force_n = force_api.integrate_maxwell_surface_force(
                magnetic_flux_density_t,
                outward_normals,
                area_weights_m2,
                permeability_H_per_m=permeability_h_per_m,
            )
            torque_nm = None
            result_pivot = None
        else:
            force_n, torque_nm = force_api.integrate_maxwell_surface_force_and_torque(
                magnetic_flux_density_t,
                outward_normals,
                area_weights_m2,
                sample_points_m,
                pivot_m=pivot_m,
                permeability_H_per_m=permeability_h_per_m,
            )
            result_pivot = pivot_m or [0.0, 0.0, 0.0]
        payload = _force_result(
            force_api,
            force_n,
            torque_nm,
            method="maxwell_surface_stress_air",
            frame=frame,
            pivot_m=result_pivot,
            sample_count=len(magnetic_flux_density_t),
            extra={"permeability_H_per_m": float(permeability_h_per_m)},
        )
    except RuntimeError as exc:
        payload = _error_payload("maxwell_surface_stress_air", "unavailable", exc)
    except (TypeError, ValueError) as exc:
        payload = _error_payload("maxwell_surface_stress_air", "invalid_input", exc)
    return _json_result(payload)


@mcp.tool()
def force_time_average_lorentz(
    current_density_real_a_per_m2: list[list[float]],
    current_density_imag_a_per_m2: list[list[float]],
    magnetic_flux_density_real_t: list[list[float]],
    magnetic_flux_density_imag_t: list[list[float]],
    volume_weights_m3: list[float],
    amplitude: str = "peak",
    sample_points_m: list[list[float]] | None = None,
    pivot_m: list[float] | None = None,
    frame: str = "global_cartesian",
) -> str:
    """Integrate cycle-averaged phasor Lorentz force and optional torque."""

    method = "time_average_lorentz_body_force"
    try:
        force_api = _load_force_api()
        current = _phasor_rows(
            current_density_real_a_per_m2,
            current_density_imag_a_per_m2,
            "current_density",
        )
        field = _phasor_rows(
            magnetic_flux_density_real_t,
            magnetic_flux_density_imag_t,
            "magnetic_flux_density",
        )
        if sample_points_m is None:
            force_n = force_api.integrate_time_average_lorentz_force(
                current,
                field,
                volume_weights_m3,
                amplitude=amplitude,
            )
            torque_nm = None
            result_pivot = None
        else:
            force_n, torque_nm = (
                force_api.integrate_time_average_lorentz_force_and_torque(
                    current,
                    field,
                    volume_weights_m3,
                    sample_points_m,
                    pivot_m=pivot_m,
                    amplitude=amplitude,
                )
            )
            result_pivot = pivot_m or [0.0, 0.0, 0.0]
        payload = _force_result(
            force_api,
            force_n,
            torque_nm,
            method=method,
            frame=frame,
            pivot_m=result_pivot,
            field_convention="time_average_phasor",
            amplitude=amplitude,
            sample_count=len(current),
        )
    except RuntimeError as exc:
        payload = _error_payload(method, "unavailable", exc)
    except (TypeError, ValueError) as exc:
        payload = _error_payload(method, "invalid_input", exc)
    return _json_result(payload)


@mcp.tool()
def force_time_average_maxwell_surface(
    magnetic_flux_density_real_t: list[list[float]],
    magnetic_flux_density_imag_t: list[list[float]],
    outward_normals: list[list[float]],
    area_weights_m2: list[float],
    amplitude: str = "peak",
    permeability_h_per_m: float = 1.2566370614359173e-6,
    sample_points_m: list[list[float]] | None = None,
    pivot_m: list[float] | None = None,
    frame: str = "global_cartesian",
) -> str:
    """Integrate cycle-averaged phasor Maxwell force and optional torque."""

    method = "time_average_maxwell_surface_stress_air"
    try:
        force_api = _load_force_api()
        field = _phasor_rows(
            magnetic_flux_density_real_t,
            magnetic_flux_density_imag_t,
            "magnetic_flux_density",
        )
        if sample_points_m is None:
            force_n = force_api.integrate_time_average_maxwell_surface_force(
                field,
                outward_normals,
                area_weights_m2,
                permeability_H_per_m=permeability_h_per_m,
                amplitude=amplitude,
            )
            torque_nm = None
            result_pivot = None
        else:
            force_n, torque_nm = (
                force_api.integrate_time_average_maxwell_surface_force_and_torque(
                    field,
                    outward_normals,
                    area_weights_m2,
                    sample_points_m,
                    pivot_m=pivot_m,
                    permeability_H_per_m=permeability_h_per_m,
                    amplitude=amplitude,
                )
            )
            result_pivot = pivot_m or [0.0, 0.0, 0.0]
        payload = _force_result(
            force_api,
            force_n,
            torque_nm,
            method=method,
            frame=frame,
            pivot_m=result_pivot,
            field_convention="time_average_phasor",
            amplitude=amplitude,
            sample_count=len(field),
            extra={"permeability_H_per_m": float(permeability_h_per_m)},
        )
    except RuntimeError as exc:
        payload = _error_payload(method, "unavailable", exc)
    except (TypeError, ValueError) as exc:
        payload = _error_payload(method, "invalid_input", exc)
    return _json_result(payload)


@mcp.tool()
def force_virtual_work(
    positions_m: list[float],
    energy_j: list[float],
    energy_kind: str = "coenergy",
) -> str:
    """Differentiate energy/coenergy versus displacement into force samples."""

    try:
        force_api = _load_force_api()
        force_n = force_api.virtual_work_force_from_displacement_samples(
            positions_m,
            energy_j,
            energy_kind=energy_kind,
        )
        payload = {
            "schema": "radia.force-sweep/v1",
            "status": "ok",
            "method": "virtual_work",
            "energy_kind": energy_kind,
            "positions_m": [float(value) for value in positions_m],
            "energy_J": [float(value) for value in energy_j],
            "force_N": force_n.tolist(),
            "sample_count": len(positions_m),
        }
    except RuntimeError as exc:
        payload = _error_payload("virtual_work", "unavailable", exc)
    except (TypeError, ValueError) as exc:
        payload = _error_payload("virtual_work", "invalid_input", exc)
    return _json_result(payload)


@mcp.tool()
def force_coenergy_torque(
    angles_rad: list[float],
    coenergy_j: list[float],
    periodic: bool = False,
    period_rad: float = 6.283185307179586,
) -> str:
    """Differentiate fixed-current coenergy versus angle into torque samples."""

    try:
        force_api = _load_force_api()
        torque_nm = force_api.coenergy_torque_from_angle_samples(
            angles_rad,
            coenergy_j,
            periodic=periodic,
            period_rad=period_rad,
        )
        payload = {
            "schema": "radia.torque-sweep/v1",
            "status": "ok",
            "method": "coenergy_virtual_work",
            "angles_rad": [float(value) for value in angles_rad],
            "coenergy_J": [float(value) for value in coenergy_j],
            "torque_Nm": torque_nm.tolist(),
            "periodic": bool(periodic),
            "period_rad": float(period_rad),
            "sample_count": len(angles_rad),
        }
    except RuntimeError as exc:
        payload = _error_payload("coenergy_virtual_work", "unavailable", exc)
    except (TypeError, ValueError) as exc:
        payload = _error_payload("coenergy_virtual_work", "invalid_input", exc)
    return _json_result(payload)


@mcp.tool()
def force_air_gap_torque(
    magnetic_flux_density_radial_t: float,
    magnetic_flux_density_tangential_t: float,
    radius_m: float,
    axial_length_m: float = 1.0,
    angle_rad: float = 6.283185307179586,
    permeability_h_per_m: float = 1.2566370614359173e-6,
    frame: str = "global_cartesian",
) -> str:
    """Compute uniform cylindrical air-gap Maxwell shear torque."""

    method = "air_gap_maxwell_shear"
    try:
        force_api = _load_force_api()
        torque_z = force_api.air_gap_shear_torque(
            magnetic_flux_density_radial_t,
            magnetic_flux_density_tangential_t,
            radius_m,
            axial_length_m=axial_length_m,
            angle_rad=angle_rad,
            permeability_H_per_m=permeability_h_per_m,
        )
        payload = _force_result(
            force_api,
            None,
            [0.0, 0.0, torque_z],
            method=method,
            frame=frame,
            pivot_m=[0.0, 0.0, 0.0],
            extra={
                "radius_m": float(radius_m),
                "axial_length_m": float(axial_length_m),
                "angle_rad": float(angle_rad),
                "permeability_H_per_m": float(permeability_h_per_m),
            },
        )
    except RuntimeError as exc:
        payload = _error_payload(method, "unavailable", exc)
    except (TypeError, ValueError) as exc:
        payload = _error_payload(method, "invalid_input", exc)
    return _json_result(payload)


@mcp.tool()
def force_air_gap_torque_samples(
    angles_rad: list[float],
    magnetic_flux_density_radial_t: list[float],
    magnetic_flux_density_tangential_t: list[float],
    radius_m: float,
    axial_length_m: float = 1.0,
    periodic: bool = True,
    period_rad: float = 6.283185307179586,
    permeability_h_per_m: float = 1.2566370614359173e-6,
    frame: str = "global_cartesian",
) -> str:
    """Integrate sampled cylindrical air-gap Maxwell shear torque."""

    method = "air_gap_maxwell_shear_samples"
    try:
        force_api = _load_force_api()
        summary = force_api.air_gap_shear_torque_from_angle_samples(
            angles_rad,
            magnetic_flux_density_radial_t,
            magnetic_flux_density_tangential_t,
            radius_m,
            axial_length_m=axial_length_m,
            periodic=periodic,
            period_rad=period_rad,
            permeability_H_per_m=permeability_h_per_m,
        )
        payload = _force_result(
            force_api,
            None,
            [0.0, 0.0, summary["torque_Nm"]],
            method=method,
            frame=frame,
            pivot_m=[0.0, 0.0, 0.0],
            sample_count=summary["n_samples"],
            extra={"air_gap_integration": summary},
        )
    except RuntimeError as exc:
        payload = _error_payload(method, "unavailable", exc)
    except (TypeError, ValueError) as exc:
        payload = _error_payload(method, "invalid_input", exc)
    return _json_result(payload)


@mcp.tool()
def force_time_average_air_gap_torque_samples(
    angles_rad: list[float],
    magnetic_flux_density_radial_real_t: list[float],
    magnetic_flux_density_radial_imag_t: list[float],
    magnetic_flux_density_tangential_real_t: list[float],
    magnetic_flux_density_tangential_imag_t: list[float],
    radius_m: float,
    axial_length_m: float = 1.0,
    periodic: bool = True,
    period_rad: float = 6.283185307179586,
    amplitude: str = "peak",
    permeability_h_per_m: float = 1.2566370614359173e-6,
    frame: str = "global_cartesian",
) -> str:
    """Integrate sampled peak/RMS phasor air-gap Maxwell shear torque."""

    method = "time_average_air_gap_maxwell_shear_samples"
    try:
        force_api = _load_force_api()
        radial = _phasor_values(
            magnetic_flux_density_radial_real_t,
            magnetic_flux_density_radial_imag_t,
            "radial flux density",
        )
        tangential = _phasor_values(
            magnetic_flux_density_tangential_real_t,
            magnetic_flux_density_tangential_imag_t,
            "tangential flux density",
        )
        summary = force_api.time_average_air_gap_shear_torque_from_angle_samples(
            angles_rad,
            radial,
            tangential,
            radius_m,
            axial_length_m=axial_length_m,
            periodic=periodic,
            period_rad=period_rad,
            permeability_H_per_m=permeability_h_per_m,
            amplitude=amplitude,
        )
        payload = _force_result(
            force_api,
            None,
            [0.0, 0.0, summary["torque_Nm"]],
            method=method,
            frame=frame,
            pivot_m=[0.0, 0.0, 0.0],
            field_convention="time_average_phasor",
            amplitude=amplitude,
            sample_count=summary["n_samples"],
            extra={"air_gap_integration": summary},
        )
    except RuntimeError as exc:
        payload = _error_payload(method, "unavailable", exc)
    except (TypeError, ValueError) as exc:
        payload = _error_payload(method, "invalid_input", exc)
    return _json_result(payload)


@_validation.tool()
def force_method_selection_gate(
    target_kind: str,
    requested_method: str,
    relative_permeability: float = 1.0,
    weighted_stress_available: bool = False,
    virtual_work_samples_available: bool = False,
    contour_clearance_mesh_layers: int = 0,
) -> str:
    """Select and gate a robust primary force-extraction method."""

    try:
        payload = _method_selection_gate(
            target_kind,
            requested_method,
            relative_permeability=relative_permeability,
            weighted_stress_available=weighted_stress_available,
            virtual_work_samples_available=virtual_work_samples_available,
            contour_clearance_mesh_layers=contour_clearance_mesh_layers,
        )
    except (TypeError, ValueError) as exc:
        payload = _error_payload("method_selection", "invalid_input", exc)
    return _json_result(payload)


@mcp.tool()
def force_method_agreement_gate(
    primary: dict,
    independent: dict,
    maximum_force_relative_difference: float = 0.05,
    maximum_torque_relative_difference: float = 0.05,
) -> str:
    """Gate two independent force/torque result records for agreement."""

    try:
        payload = _method_agreement_gate(
            primary,
            independent,
            maximum_force_relative_difference=maximum_force_relative_difference,
            maximum_torque_relative_difference=maximum_torque_relative_difference,
        )
    except (TypeError, ValueError) as exc:
        payload = _error_payload("method_agreement", "invalid_input", exc)
    return _json_result(payload)


@_validation.tool()
def force_action_reaction_gate(
    force_a_n: list[float],
    force_b_n: list[float],
    torque_a_nm: list[float] | None = None,
    torque_b_nm: list[float] | None = None,
    maximum_force_relative_residual: float = 0.01,
    maximum_torque_relative_residual: float = 0.01,
) -> str:
    """Gate Newton action-reaction closure for force and common-pivot torque."""

    try:
        payload = _action_reaction_gate(
            force_a_n,
            force_b_n,
            torque_a_Nm=torque_a_nm,
            torque_b_Nm=torque_b_nm,
            maximum_force_relative_residual=maximum_force_relative_residual,
            maximum_torque_relative_residual=maximum_torque_relative_residual,
        )
    except (TypeError, ValueError) as exc:
        payload = _error_payload("action_reaction", "invalid_input", exc)
    return _json_result(payload)


@_validation.tool()
def force_weight_equilibrium_gate(
    force_n: list[float],
    mass_kg: float,
    lift_axis: int = 2,
    gravity_m_per_s2: float = 9.80665,
    maximum_relative_residual: float = 0.02,
) -> str:
    """Gate levitation or bearing lift against weight."""

    try:
        payload = _weight_equilibrium_gate(
            force_n,
            mass_kg,
            lift_axis=lift_axis,
            gravity_m_per_s2=gravity_m_per_s2,
            maximum_relative_residual=maximum_relative_residual,
        )
    except (TypeError, ValueError) as exc:
        payload = _error_payload("weight_equilibrium", "invalid_input", exc)
    return _json_result(payload)


_validation.install()

register_status_tool(
    mcp,
    server_name="mcp-server-force",
    description=(
        "Shared electromagnetic-force and torque method selection, "
        "solver-independent sample integration, and validation guidance"
    ),
    subpackage="radia_mcp.force",
    related_servers=["motor", "maglev", "differential-forms", "radia-ngsolve"],
    optional_deps=["radia"],
)


register_topics_tool(
    mcp,
    server_name="mcp-server-force",
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        assert "common electromagnetic-force layer" in force("overview")
        assert "EM-force method" in force_recipe("method_choice")
        assert len(force_validation_guide("method_map")) > 100
        print("Force MCP server self-test:")
        print(f"  topics: {len(TOPICS)}")
        print(f"  tools:  {len(mcp._tool_manager._tools)}")
        print("  PASSED")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
