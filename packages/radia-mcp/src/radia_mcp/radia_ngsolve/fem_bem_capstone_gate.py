"""Solver-neutral first-order FEM/BEM capstone-suite gate."""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any


def fem_bem_capstone_suite_gate(payload: Mapping[str, Any]) -> dict[str, object]:
    cases = payload.get("cases")
    capabilities = payload.get("capabilities")
    if not isinstance(cases, list) or len(cases) != 10 or not all(isinstance(row, Mapping) for row in cases):
        raise ValueError("cases must contain exactly ten mappings")
    if not isinstance(capabilities, Mapping):
        raise ValueError("capabilities must be a mapping")
    indexed = {str(row.get("id")): row for row in cases}
    expected = [f"GYP-{index:03d}" for index in range(91, 101)]
    if sorted(indexed) != expected:
        raise ValueError(f"case ids must be exactly {expected}")

    def details(case_id: str) -> Mapping[str, Any]:
        value = indexed[case_id].get("details")
        return value if isinstance(value, Mapping) else {}

    def finite(case_id: str, name: str) -> float:
        try:
            value = float(details(case_id)[name])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{case_id}.{name} must be finite") from exc
        if not math.isfinite(value):
            raise ValueError(f"{case_id}.{name} must be finite")
        return value

    low = details("GYP-096").get("lowFrequencyValue")
    if not isinstance(low, Mapping):
        raise ValueError("GYP-096.lowFrequencyValue must be a complex-value mapping")
    low_real = float(low.get("real"))
    low_imag = float(low.get("imag"))
    h1_error = finite("GYP-100", "h1Error")
    hcurl_error = finite("GYP-100", "hcurlError")
    checks = {
        "all_source_cases_passed": all(row.get("passed") is True and not row.get("failures") for row in cases),
        "single_tet_mesh_identity": details("GYP-091").get("meshVertices") == 4
        and details("GYP-091").get("meshElements") == 1,
        "h1_matrix_matches_reference": finite("GYP-092", "error") <= 1.0e-12,
        "hcurl_matrix_matches_reference": finite("GYP-093", "error") <= 1.0e-12,
        "laplace_single_layer_available": details("GYP-094").get("hasLaplaceSL") is True,
        "helmholtz_single_layer_available": details("GYP-095").get("hasHelmholtzSL") is True,
        "low_frequency_kernel_reaches_laplace_limit": abs(low_real - 1.0 / (4.0 * math.pi)) <= 1.0e-12
        and abs(low_imag) <= 1.0e-8,
        "analytic_sphere_observables_positive": finite("GYP-097", "meanPotential") > 0.0
        and finite("GYP-098", "meanAmplitude") > 0.0,
        "p1_trace_has_four_boundary_rows": details("GYP-099").get("traceRows") == 4,
        "combined_h1_hcurl_capstone_matches": h1_error <= 1.0e-12 and hcurl_error <= 1.0e-12,
        "reference_capability_shape_matches": capabilities.get("ok") is True
        and capabilities.get("mesh_vertices") == 4
        and capabilities.get("mesh_elements") == 1
        and capabilities.get("h1_dofs") == 4
        and capabilities.get("hcurl_dofs") == 6
        and capabilities.get("has_bem") is True
        and capabilities.get("has_laplace_sl") is True
        and capabilities.get("has_helmholtz_sl") is True,
    }
    return {
        "policy": "fem_bem_capstone_suite_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "case_count": len(cases),
            "h1_combined_error": h1_error,
            "hcurl_combined_error": hcurl_error,
            "low_frequency_real_error": abs(low_real - 1.0 / (4.0 * math.pi)),
            "low_frequency_imaginary_magnitude": abs(low_imag),
            "trace_row_count": details("GYP-099").get("traceRows"),
        },
        "lesson": (
            "A readable FEM/BEM capstone should close local H1 and HCurl matrices, low-frequency kernel "
            "limits, scalar BEM availability, and the volume-to-boundary P1 trace against one consistent mesh contract."
        ),
    }
