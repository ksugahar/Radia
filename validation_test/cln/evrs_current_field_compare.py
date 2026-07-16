"""Compare reduced EVRS currents against a full p=6 HCurl parent solve.

This validation lane checks the field-level questions behind the p-order
discussion:

1. L2 / energy norm of ``J = curl(T)`` against a full p=6 parent solve.
2. Joule-loss error ``int |J|^2 / sigma dV``.
3. Local current error near the notched-box re-entrant corner.
4. p=4 reduced versus p=6 reduced versus full p=6.

The solve is the same dimensionless HCurl parent model used to generate the
EVRS basis, ``(K + shift M) T = b``.  It is a desktop smoke for the reduction
mechanism, not a production motor benchmark.

Run from the repository root:

    python validation_test/cln/evrs_current_field_compare.py
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
SRC = REPO_ROOT / "src"
DEFAULT_OUTPUT = HERE / "evrs_current_field_compare_smoke.json"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import radia.vim as vim
from evrs_sibc_mixed_schur import (  # noqa: E402
    _make_skin_mesh,
    _parse_floats,
    _parse_ints,
    _port_samples,
    _response_basis,
)


def _assemble_parent_on_mesh(mesh, order: int, condense: bool):
    import ngsolve as ng

    fes = ng.HCurl(mesh, order=order, nograds=True)
    u, v = fes.TnT()

    stiffness = ng.BilinearForm(fes, condense=condense)
    stiffness += ng.curl(u) * ng.curl(v) * ng.dx + 0.05 * u * v * ng.dx
    mass = ng.BilinearForm(fes)
    mass += u * v * ng.dx

    ports = []
    for cf in (
        ng.CoefficientFunction((-ng.y, ng.x, 0.0)),
        ng.CoefficientFunction((0.0, -ng.z, ng.y)),
    ):
        port = ng.LinearForm(fes)
        port += cf * v * ng.dx
        ports.append(port)

    with ng.TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        for port in ports:
            port.Assemble()

    return fes, stiffness, mass, ports


def _free_mask(fes, condense: bool) -> np.ndarray:
    return np.asarray(fes.FreeDofs(condense), dtype=bool)


def _dense_parent_arrays(stiffness, mass, ports):
    k = vim.NgsolveMatrixToDense(stiffness, dtype=np.complex128)
    m = vim.NgsolveMatrixToDense(mass, dtype=np.complex128)
    b = np.column_stack(
        [vim.NgsolveVectorToArray(port, dtype=np.complex128) for port in ports]
    )
    return k, m, b


def _solve_full(k: np.ndarray, m: np.ndarray, b: np.ndarray, free: np.ndarray, shift: float):
    active = np.flatnonzero(free)
    a = k + float(shift) * m
    af = a[np.ix_(active, active)]
    bf = b[active, :]
    solved = np.linalg.solve(af, bf)
    out = np.zeros_like(b, dtype=np.complex128)
    out[active, :] = solved
    return out


def _solve_reduced(
    k: np.ndarray,
    m: np.ndarray,
    b: np.ndarray,
    free: np.ndarray,
    response: vim.ResponseBasis,
    shift: float,
):
    active = np.flatnonzero(free)
    q = np.asarray(response.vectors, dtype=np.complex128)[active, :]
    a = k + float(shift) * m
    af = a[np.ix_(active, active)]
    bf = b[active, :]
    ar = q.conj().T @ af @ q
    br = q.conj().T @ bf
    yr = np.linalg.solve(ar, br)
    solved = q @ yr
    out = np.zeros_like(b, dtype=np.complex128)
    out[active, :] = solved
    return out


def _sample_current(mesh, fes, coeffs: np.ndarray, *, intorder: int, prefix: str):
    arr = np.asarray(coeffs)
    if np.iscomplexobj(arr):
        imag_scale = float(np.max(np.abs(arr.imag))) if arr.size else 0.0
        real_scale = max(float(np.max(np.abs(arr.real))) if arr.size else 0.0, 1.0)
        if imag_scale > 1.0e-12 * real_scale:
            raise ValueError("NGSolve current sampling currently expects real coefficients")
        arr = arr.real
    return vim.NgsolveHCurlCurlBasis(
        mesh,
        fes,
        arr,
        intorder=intorder,
        materials="cond",
        names=[f"{prefix}_port{i}" for i in range(coeffs.shape[1])],
    )


def _weighted_norms(modes: np.ndarray, weights: np.ndarray) -> np.ndarray:
    values = np.einsum("psj,psj,s->p", modes.conj(), modes, weights)
    return np.maximum(values.real, 0.0)


def _corner_mask(points: np.ndarray, geometry: str, radius: float) -> np.ndarray:
    if geometry == "notched-box":
        line_distance = np.sqrt((points[:, 0] - 0.45) ** 2 + (points[:, 1] - 0.45) ** 2)
        return line_distance <= radius
    if geometry == "l-prism":
        line_distance = np.sqrt((points[:, 0] - 0.45) ** 2 + (points[:, 1] - 0.45) ** 2)
        return line_distance <= radius
    return np.linalg.norm(points, axis=1) <= radius


def _field_metrics(candidate, reference, *, sigma: float, geometry: str, corner_radius: float):
    if not np.allclose(candidate.points, reference.points, rtol=1.0e-10, atol=1.0e-12):
        raise RuntimeError("candidate and reference current samples are not colocated")
    if not np.allclose(candidate.weights, reference.weights, rtol=1.0e-10, atol=1.0e-12):
        raise RuntimeError("candidate and reference current weights are not identical")

    weights = reference.weights
    diff = candidate.modes - reference.modes
    ref_norm2 = _weighted_norms(reference.modes, weights)
    diff_norm2 = _weighted_norms(diff, weights)
    cand_norm2 = _weighted_norms(candidate.modes, weights)
    rel_l2 = np.sqrt(diff_norm2 / np.maximum(ref_norm2, 1.0e-300))
    loss_ref = ref_norm2 / sigma
    loss_candidate = cand_norm2 / sigma
    loss_rel = np.abs(loss_candidate - loss_ref) / np.maximum(np.abs(loss_ref), 1.0e-300)

    mask = _corner_mask(reference.points, geometry, corner_radius)
    if not np.any(mask):
        local_rel = np.full(reference.n_modes, np.nan)
        local_peak_rel = np.full(reference.n_modes, np.nan)
        local_sample_count = 0
    else:
        local_weights = weights[mask]
        local_ref = reference.modes[:, mask, :]
        local_diff = diff[:, mask, :]
        local_ref_norm2 = _weighted_norms(local_ref, local_weights)
        local_diff_norm2 = _weighted_norms(local_diff, local_weights)
        local_rel = np.sqrt(local_diff_norm2 / np.maximum(local_ref_norm2, 1.0e-300))

        ref_mag = np.linalg.norm(local_ref, axis=2)
        cand_mag = np.linalg.norm(candidate.modes[:, mask, :], axis=2)
        local_peak_rel = np.max(np.abs(cand_mag - ref_mag), axis=1) / np.maximum(
            np.max(ref_mag, axis=1),
            1.0e-300,
        )
        local_sample_count = int(np.count_nonzero(mask))

    return {
        "ports": [
            {
                "port": int(i),
                "relative_current_l2": float(rel_l2[i]),
                "relative_energy_norm": float(rel_l2[i]),
                "joule_loss_reference": float(loss_ref[i]),
                "joule_loss_candidate": float(loss_candidate[i]),
                "relative_joule_loss_error": float(loss_rel[i]),
                "corner_relative_current_l2": float(local_rel[i]),
                "corner_peak_magnitude_relative_error": float(local_peak_rel[i]),
            }
            for i in range(reference.n_modes)
        ],
        "max_relative_current_l2": float(np.nanmax(rel_l2)),
        "max_relative_energy_norm": float(np.nanmax(rel_l2)),
        "max_relative_joule_loss_error": float(np.nanmax(loss_rel)),
        "max_corner_relative_current_l2": float(np.nanmax(local_rel)),
        "max_corner_peak_magnitude_relative_error": float(np.nanmax(local_peak_rel)),
        "corner_sample_count": local_sample_count,
    }


def _case_result(
    mesh,
    *,
    order: int,
    steps: int,
    condense: bool,
    shift: float,
    sigma: float,
    intorder: int,
    rtol: float,
    reference_current,
    geometry: str,
    corner_radius: float,
):
    t0 = time.perf_counter()
    fes, stiffness, mass, ports = _assemble_parent_on_mesh(mesh, order, condense)
    k, m, b = _dense_parent_arrays(stiffness, mass, ports)
    free = _free_mask(fes, condense)
    response = _response_basis(
        fes,
        stiffness,
        mass,
        ports,
        steps=steps,
        condense=condense,
        rtol=rtol,
    )
    coeffs = _solve_reduced(k, m, b, free, response, shift)
    current = _sample_current(
        mesh,
        fes,
        coeffs,
        intorder=intorder,
        prefix=f"p{order}_n{steps}",
    )
    metrics = _field_metrics(
        current,
        reference_current,
        sigma=sigma,
        geometry=geometry,
        corner_radius=corner_radius,
    )
    elapsed = time.perf_counter() - t0
    info = response.diagnostics()
    return {
        "kind": "reduced",
        "order": order,
        "krylov_steps": steps,
        "shift": float(shift),
        "condensed": condense,
        "ndof": int(info["ndof"]),
        "active_dofs": int(info["active_dofs"]),
        "rank": int(info["rank"]),
        "compression_ratio": float(info["compression_ratio"]),
        "response_seconds_lab_smoke": float(elapsed),
        "current_samples": int(current.n_samples),
        **metrics,
    }


def run_compare(args: argparse.Namespace) -> dict[str, object]:
    import ngsolve as ng

    mesh = _make_skin_mesh(args.maxh, args.geometry)
    reference_condense = args.reference_order >= args.condense_from
    ref_fes, ref_stiffness, ref_mass, ref_ports = _assemble_parent_on_mesh(
        mesh,
        args.reference_order,
        reference_condense,
    )
    ref_k, ref_m, ref_b = _dense_parent_arrays(ref_stiffness, ref_mass, ref_ports)
    ref_free = _free_mask(ref_fes, reference_condense)

    cases = []
    for shift in args.shifts:
        ref_coeffs = _solve_full(ref_k, ref_m, ref_b, ref_free, shift)
        ref_current = _sample_current(
            mesh,
            ref_fes,
            ref_coeffs,
            intorder=args.intorder,
            prefix=f"full_p{args.reference_order}",
        )
        cases.append(
            {
                "kind": "full-reference",
                "order": args.reference_order,
                "krylov_steps": None,
                "shift": float(shift),
                "condensed": reference_condense,
                "ndof": int(ref_fes.ndof),
                "active_dofs": int(np.count_nonzero(ref_free)),
                "rank": int(np.count_nonzero(ref_free)),
                "compression_ratio": 1.0,
                "current_samples": int(ref_current.n_samples),
                "max_relative_current_l2": 0.0,
                "max_relative_energy_norm": 0.0,
                "max_relative_joule_loss_error": 0.0,
                "max_corner_relative_current_l2": 0.0,
                "max_corner_peak_magnitude_relative_error": 0.0,
                "corner_sample_count": int(
                    np.count_nonzero(
                        _corner_mask(ref_current.points, args.geometry, args.corner_radius)
                    )
                ),
            }
        )

        for order in args.orders:
            condense = order >= args.condense_from
            for steps in args.steps:
                cases.append(
                    _case_result(
                        mesh,
                        order=order,
                        steps=steps,
                        condense=condense,
                        shift=shift,
                        sigma=args.sigma,
                        intorder=args.intorder,
                        rtol=args.rtol,
                        reference_current=ref_current,
                        geometry=args.geometry,
                        corner_radius=args.corner_radius,
                    )
                )

    return {
        "schema": "radia.validation.evrs_current_field_compare.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "validation_host": platform.node(),
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "ngsolve_version": getattr(ng, "__version__", "unknown"),
        },
        "note": (
            "LAB/desktop smoke only.  The parent solve uses the dimensionless "
            "(K + shift M) HCurl model; timings are not benchmark data."
        ),
        "configuration": {
            "geometry": args.geometry,
            "orders": args.orders,
            "krylov_steps": args.steps,
            "reference_order": args.reference_order,
            "shifts": args.shifts,
            "maxh": args.maxh,
            "sigma": args.sigma,
            "intorder": args.intorder,
            "corner_radius": args.corner_radius,
            "condense_from": args.condense_from,
            "rtol": args.rtol,
        },
        "cases": cases,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", choices=("box", "bar", "notched-box", "l-prism"), default="notched-box")
    parser.add_argument("--orders", type=_parse_ints, default=[4, 6])
    parser.add_argument("--steps", type=_parse_ints, default=[22])
    parser.add_argument("--reference-order", type=int, default=6)
    parser.add_argument("--shifts", type=_parse_floats, default=[1.0])
    parser.add_argument("--maxh", type=float, default=2.0)
    parser.add_argument("--sigma", type=float, default=5.8e7)
    parser.add_argument("--intorder", type=int, default=2)
    parser.add_argument("--corner-radius", type=float, default=0.35)
    parser.add_argument("--condense-from", type=int, default=7)
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)

    if args.reference_order < 1:
        parser.error("--reference-order must be positive")
    if args.corner_radius <= 0.0:
        parser.error("--corner-radius must be positive")
    if any(order < 1 for order in args.orders):
        parser.error("--orders must contain positive integers")

    result = run_compare(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    print("EVRS current-field comparison smoke")
    print(f"  output: {args.output}")
    print("  kind      p   n  shift  rank/active  J_L2     loss     corner_L2  corner_peak")
    for case in result["cases"]:
        n = "-" if case["krylov_steps"] is None else str(case["krylov_steps"])
        print(
            f"  {case['kind']:<9} "
            f"{case['order']:>2} {n:>3} "
            f"{case['shift']:>6.3g} "
            f"{case['rank']:>4}/{case['active_dofs']:<5} "
            f"{case['max_relative_current_l2']:.3e} "
            f"{case['max_relative_joule_loss_error']:.3e} "
            f"{case['max_corner_relative_current_l2']:.3e} "
            f"{case['max_corner_peak_magnitude_relative_error']:.3e}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
