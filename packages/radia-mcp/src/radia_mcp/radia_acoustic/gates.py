"""Cheap pre-solve gates for FSI and convolution quadrature."""

from __future__ import annotations

import cmath
import math
import sys


def fsi_preflight_gate(*, wavenumber: float, c_longitudinal: float = 2.0,
                       c_transverse: float = 1.0, solid_density: float = 1.5,
                       fluid_density: float = 1.0, order: int = 1,
                       boundary: str = "gamma", radius_deviation: float = 0.0,
                       radius_tolerance: float = 0.03) -> dict:
    checks = {
        "positive_wavenumber": wavenumber > 0,
        "positive_wave_speeds": c_longitudinal > 0 and c_transverse > 0,
        "positive_densities": solid_density > 0 and fluid_density > 0,
        "positive_lame_lambda": c_longitudinal**2 > 2 * c_transverse**2,
        "valid_order": isinstance(order, int) and order >= 1,
        "boundary_named": bool(boundary.strip()),
        "spherical_dtn_geometry": 0 <= radius_deviation <= radius_tolerance,
    }
    return {
        "schema": "radia-mcp.radia-acoustic-fsi-preflight/v1",
        "ok": all(checks.values()), "checks": checks,
        "formulation": "NGSolve VectorH1 elasticity + spherical Helmholtz DtN",
        "next": "call radia.acoustics.fsi.fsi_dtn_solve only when ok is true",
    }


def cq_grid_gate(*, num_time: int, time_step: float, sound_speed: float = 1.0,
                 method: str = "BDF2") -> dict:
    method = method.upper()
    basic = isinstance(num_time, int) and num_time >= 4 and time_step > 0 and sound_speed > 0 and method in {"BDF1", "BDF2"}
    if not basic:
        return {"schema": "radia-mcp.radia-acoustic-cq-grid/v1", "ok": False,
                "errors": ["require num_time>=4, time_step>0, sound_speed>0, method BDF1/BDF2"]}
    rho = (sys.float_info.epsilon ** 0.5) ** (1.0 / num_time)
    zeta = [rho * cmath.exp(-2j * math.pi * n / num_time) for n in range(num_time)]
    if method == "BDF1":
        delta = [1.0 - value for value in zeta]
    else:
        delta = [1.5 - 2.0 * value + 0.5 * value * value for value in zeta]
    s = [value / time_step for value in delta]
    kappa = [1j * value / sound_speed for value in s]

    def _finite(value: complex) -> bool:
        return math.isfinite(value.real) and math.isfinite(value.imag)

    checks = {
        "laplace_nodes_right_half_plane": all(value.real > 0 for value in s),
        "finite_nodes": all(_finite(value) for value in (*s, *kappa)),
        "conjugate_pair_symmetry": max(
            abs(s[index] - s[-index].conjugate())
            for index in range(1, num_time)
        ) < 1e-10,
    }
    return {
        "schema": "radia-mcp.radia-acoustic-cq-grid/v1", "ok": all(checks.values()),
        "checks": checks, "method": method, "num_time": num_time,
        "time_step": time_step, "sound_speed": sound_speed, "cq_radius": float(rho),
        "min_real_s": min(value.real for value in s),
        "max_abs_kappa": max(abs(value) for value in kappa),
        "convention": "s=delta(zeta)/dt; kappa=i*s/c",
    }
