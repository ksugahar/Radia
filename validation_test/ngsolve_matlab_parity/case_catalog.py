"""Deterministic 100-case catalog for NGSolve Python/MATLAB parity."""

from __future__ import annotations

from typing import Any


DOF_LIMIT = 1_000_000

MESH_SPECS: tuple[dict[str, Any], ...] = (
    {"id": "square_coarse", "dimension": 2, "generator": "unit_square", "maxh": 0.45},
    {"id": "square_fine", "dimension": 2, "generator": "unit_square", "maxh": 0.25},
    {"id": "cube_coarse", "dimension": 3, "generator": "unit_cube", "maxh": 0.75},
    {"id": "cube_fine", "dimension": 3, "generator": "unit_cube", "maxh": 0.50},
)


def _operator_specs() -> tuple[dict[str, Any], ...]:
    specs: list[dict[str, Any]] = []
    for order in range(1, 5):
        specs.append(_operator("h1", order, "mass"))
    for order in range(1, 5):
        specs.append(_operator("h1", order, "stiffness", dirichlet=".*"))
    for order in range(1, 4):
        specs.append(_operator("hcurl", order, "mass"))
    for order in range(1, 4):
        specs.append(_operator("hcurl", order, "curlcurl", solve=False))
    for order in range(1, 4):
        specs.append(_operator("hdiv", order, "mass"))
    for order in range(1, 4):
        specs.append(_operator("hdiv", order, "divdiv", solve=False))

    specs.extend(
        (
            _operator("h1", 3, "mass", weight=0.75, variant="weighted075"),
            _operator(
                "h1", 3, "stiffness", weight=2.0,
                dirichlet=".*", variant="weighted200",
            ),
            _operator("hcurl", 2, "mass", weight=1.25, variant="weighted125"),
            _operator("hdiv", 2, "mass", weight=0.5, variant="weighted050"),
            _operator(
                "hdiv", 2, "divdiv", weight=1.5,
                solve=False, variant="weighted150",
            ),
        )
    )
    assert len(specs) == 25
    return tuple(specs)


def _operator(
    space: str,
    order: int,
    form: str,
    *,
    weight: float = 1.0,
    dirichlet: str = "",
    solve: bool = True,
    variant: str = "base",
) -> dict[str, Any]:
    return {
        "space": space,
        "order": order,
        "form": form,
        "weight": weight,
        "dirichlet": dirichlet,
        "solve": solve,
        "variant": variant,
    }


OPERATOR_SPECS = _operator_specs()


def build_case_catalog() -> list[dict[str, Any]]:
    """Return four meshes times 25 operators in a stable order."""

    cases: list[dict[str, Any]] = []
    for mesh in MESH_SPECS:
        for operator in OPERATOR_SPECS:
            number = len(cases) + 1
            variant = "" if operator["variant"] == "base" else f"_{operator['variant']}"
            case_id = (
                f"{number:03d}_{mesh['id']}_{operator['space']}_"
                f"p{operator['order']}_{operator['form']}{variant}"
            )
            cases.append(
                {
                    "case_id": case_id,
                    "oracle_key": f"case_{number:03d}",
                    "mesh_id": mesh["id"],
                    "dimension": mesh["dimension"],
                    **operator,
                }
            )
    assert len(cases) == 100
    return cases

