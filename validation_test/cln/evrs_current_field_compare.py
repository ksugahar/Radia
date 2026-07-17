"""Compare reduced EVRS/VIM currents against a full high-order HCurl-FEM solve.

This validation lane checks the field-level questions behind the p-order
discussion:

1. L2 / energy norm of ``J = curl(T)`` against a full p=6 parent solve.
2. Joule-loss error ``int |J|^2 / sigma dV``.
3. Local current error near a re-entrant corner or curved surface layer.
4. p=3--6 reduced versus full p=6 on one shared physical mesh.
5. Curved surface area, volume, tangential trace, and cycle-bridge geometry.

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


def _response_basis(
    fes,
    stiffness,
    mass,
    ports,
    *,
    steps: int,
    condense: bool,
    inverse: str,
    rtol: float,
):
    if condense:
        return vim.NgsolveStaticCondensedBlockKrylovBasis(
            stiffness,
            mass,
            ports,
            steps=steps,
            free_dofs=fes.FreeDofs(True),
            inverse=inverse,
            rtol=rtol,
        )
    return vim.NgsolveOperatorBlockKrylovBasis(
        stiffness,
        mass,
        ports,
        steps=steps,
        free_dofs=fes.FreeDofs(False),
        inverse=inverse,
        rtol=rtol,
    )


def _free_mask(fes, condense: bool) -> np.ndarray:
    return np.asarray(fes.FreeDofs(condense), dtype=bool)


def _shifted_sparse_matrix(stiffness, mass, shift: float):
    k = stiffness.mat
    m = mass.mat
    operator = k.CreateMatrix()
    operator.AsVector().data = k.AsVector() + float(shift) * m.AsVector()
    return operator


def _solve_full(stiffness, mass, ports, free_dofs, shift: float, inverse: str):
    t0 = time.perf_counter()
    operator = _shifted_sparse_matrix(stiffness, mass, shift)
    solver = operator.Inverse(freedofs=free_dofs, inverse=inverse)
    columns = []
    for port in ports:
        solution = operator.CreateColVector()
        solution.data = solver * port.vec
        columns.append(np.array(solution.FV().NumPy(), copy=True))
    return np.column_stack(columns), time.perf_counter() - t0


def _apply_matrix(matrix, vectors: np.ndarray) -> np.ndarray:
    operator = matrix.mat
    columns = []
    for column in vectors.T:
        source = operator.CreateColVector()
        source.FV().NumPy()[:] = column
        target = source.CreateVector()
        target.data = operator * source
        columns.append(np.array(target.FV().NumPy(), copy=True))
    return np.column_stack(columns)


def _project_reduced(stiffness, mass, ports, response: vim.ResponseBasis):
    q = np.asarray(response.vectors)
    kq = _apply_matrix(stiffness, q)
    mq = _apply_matrix(mass, q)
    b = np.column_stack([vim.NgsolveVectorToArray(port) for port in ports])
    kr = q.conj().T @ kq
    mr = q.conj().T @ mq
    br = q.conj().T @ b
    return 0.5 * (kr + kr.conj().T), 0.5 * (mr + mr.conj().T), br


def _solve_reduced(
    response: vim.ResponseBasis,
    kr: np.ndarray,
    mr: np.ndarray,
    br: np.ndarray,
    shift: float,
    repeats: int,
):
    t0 = time.perf_counter()
    operator = kr + float(shift) * mr
    for _ in range(repeats):
        reduced = np.linalg.solve(operator, br)
    coefficients = np.asarray(response.vectors) @ reduced
    return coefficients, (time.perf_counter() - t0) / repeats


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
    if geometry == "sphere":
        return np.linalg.norm(points, axis=1) >= 1.0 - radius
    return np.linalg.norm(points, axis=1) <= radius


def _geometry_diagnostics(mesh, geometry: str, intorder: int) -> dict[str, object]:
    import ngsolve as ng

    area = float(
        ng.Integrate(
            1.0,
            mesh,
            definedon=mesh.Boundaries("skin"),
            order=intorder,
        )
    )
    volume = float(ng.Integrate(1.0, mesh, order=intorder))
    surface = vim.NgsolveSurfaceOmegaBasis(
        mesh,
        (ng.CoefficientFunction((1.0, 0.0, 0.0)),),
        intorder=intorder,
        boundaries="skin",
    )
    _, _, normals = vim.SampleNgsolveVectorCFs(
        mesh,
        (ng.specialcf.normal(3),),
        vb="BND",
        intorder=intorder,
        boundaries="skin",
    )
    tangential_defect = float(
        np.max(np.abs(np.einsum("ij,ij->i", surface.modes[0], normals[0])))
    )
    topology = vim.ClassifyNgsolveEddyTopology(mesh, conductive_materials="cond")
    bridge = vim.NgsolveBridgeCycleCurrentBasis(
        mesh,
        topology,
        geometry_intorder=intorder,
    )
    exact_area = 4.0 * np.pi if geometry == "sphere" else None
    exact_volume = 4.0 * np.pi / 3.0 if geometry == "sphere" else None
    return {
        "mesh_curve_order": int(mesh.GetCurveOrder()),
        "elements": int(mesh.ne),
        "surface_area_fem": area,
        "surface_area_vim_samples": float(surface.weights.sum()),
        "surface_area_vim_to_fem_relative_error": float(
            abs(surface.weights.sum() - area) / area
        ),
        "volume_fem": volume,
        "surface_tangential_defect": tangential_defect,
        "conductor_graph_edges": topology.conductor_graph().edge_count,
        "conductor_graph_cycle_rank": topology.conductor_graph().cycle_rank,
        "bridge_dual_volume_sum": float(bridge.weights.sum()),
        "exact_surface_area": exact_area,
        "surface_area_relative_error_to_exact": (
            None if exact_area is None else float(abs(area - exact_area) / exact_area)
        ),
        "exact_volume": exact_volume,
        "volume_relative_error_to_exact": (
            None if exact_volume is None else float(abs(volume - exact_volume) / exact_volume)
        ),
    }


def _refine_reentrant_corner(
    mesh,
    geometry: str,
    levels: int,
    radius: float,
) -> list[dict[str, int | float]]:
    import ngsolve as ng

    if levels == 0:
        return []
    if geometry not in {"notched-box", "l-prism"}:
        raise ValueError("corner refinement requires notched-box or l-prism geometry")
    history = []
    for level in range(levels):
        candidates = []
        for element in mesh.Elements(ng.VOL):
            vertices = np.vstack(
                [
                    np.asarray(mesh.vertices[vertex.nr].point, dtype=float)
                    for vertex in element.vertices
                ]
            )
            centroid = np.mean(vertices, axis=0)
            distance = float(np.linalg.norm(centroid[:2] - (0.45, 0.45)))
            candidates.append((element.nr, distance))
        level_radius = radius / (2**level)
        marked = [nr for nr, distance in candidates if distance <= level_radius]
        if not marked:
            marked = [min(candidates, key=lambda item: item[1])[0]]
        marked_set = set(marked)
        elements_before = int(mesh.ne)
        for element_nr, _ in candidates:
            mesh.SetRefinementFlag(
                ng.ElementId(ng.VOL, element_nr),
                element_nr in marked_set,
            )
        mesh.Refine()
        history.append(
            {
                "level": level + 1,
                "radius": float(level_radius),
                "marked_parent_elements": len(marked),
                "elements_before": elements_before,
                "elements_after": int(mesh.ne),
            }
        )
    return history


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
    fes,
    response,
    kr,
    mr,
    br,
    *,
    order: int,
    steps: int,
    condense: bool,
    shift: float,
    sigma: float,
    intorder: int,
    reference_current,
    full_solve_seconds: float,
    parent_assembly_seconds: float,
    offline_basis_seconds: float,
    online_repeats: int,
    geometry: str,
    corner_radius: float,
):
    coeffs, online_seconds = _solve_reduced(
        response,
        kr,
        mr,
        br,
        shift,
        online_repeats,
    )
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
    info = response.diagnostics()
    saved_seconds = full_solve_seconds - online_seconds
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
        "parent_assembly_seconds": float(parent_assembly_seconds),
        "offline_basis_seconds": float(offline_basis_seconds),
        "online_reduced_solve_seconds": float(online_seconds),
        "online_timing_repeats": int(online_repeats),
        "full_reference_solve_seconds": float(full_solve_seconds),
        "online_speedup": float(full_solve_seconds / max(online_seconds, 1.0e-300)),
        "break_even_online_cases": (
            None
            if saved_seconds <= 0.0
            else float(offline_basis_seconds / saved_seconds)
        ),
        "current_samples": int(current.n_samples),
        **metrics,
    }


def run_compare(args: argparse.Namespace) -> dict[str, object]:
    import ngsolve as ng

    started = time.perf_counter()
    mesh = _make_skin_mesh(
        args.maxh,
        args.geometry,
        args.curve_order,
        args.corner_edge_maxh,
    )
    refinement_history = _refine_reentrant_corner(
        mesh,
        args.geometry,
        args.corner_refinements,
        args.refinement_radius,
    )
    geometry = _geometry_diagnostics(mesh, args.geometry, args.geometry_intorder)
    reference_condense = args.reference_order >= args.condense_from
    if reference_condense:
        raise ValueError(
            "full FEM reference must be uncondensed; set --condense-from above "
            "--reference-order"
        )
    t0 = time.perf_counter()
    ref_fes, ref_stiffness, ref_mass, ref_ports = _assemble_parent_on_mesh(
        mesh,
        args.reference_order,
        reference_condense,
    )
    reference_assembly_seconds = time.perf_counter() - t0
    ref_free = _free_mask(ref_fes, reference_condense)
    ref_free_dofs = ref_fes.FreeDofs(reference_condense)

    cases = []
    reference_currents = {}
    reference_solve_seconds = {}
    for shift in args.shifts:
        ref_coeffs, solve_seconds = _solve_full(
            ref_stiffness,
            ref_mass,
            ref_ports,
            ref_free_dofs,
            shift,
            args.inverse,
        )
        ref_current = _sample_current(
            mesh,
            ref_fes,
            ref_coeffs,
            intorder=args.intorder,
            prefix=f"full_p{args.reference_order}",
        )
        reference_currents[float(shift)] = ref_current
        reference_solve_seconds[float(shift)] = solve_seconds
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
                "parent_assembly_seconds": float(reference_assembly_seconds),
                "full_reference_solve_seconds": float(solve_seconds),
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
        t0 = time.perf_counter()
        fes, stiffness, mass, ports = _assemble_parent_on_mesh(mesh, order, condense)
        parent_assembly_seconds = time.perf_counter() - t0
        for steps in args.steps:
            t0 = time.perf_counter()
            response = _response_basis(
                fes,
                stiffness,
                mass,
                ports,
                steps=steps,
                condense=condense,
                inverse=args.inverse,
                rtol=args.rtol,
            )
            kr, mr, br = _project_reduced(stiffness, mass, ports, response)
            offline_basis_seconds = time.perf_counter() - t0
            for shift in args.shifts:
                cases.append(
                    _case_result(
                        mesh,
                        fes,
                        response,
                        kr,
                        mr,
                        br,
                        order=order,
                        steps=steps,
                        condense=condense,
                        shift=shift,
                        sigma=args.sigma,
                        intorder=args.intorder,
                        reference_current=reference_currents[float(shift)],
                        full_solve_seconds=reference_solve_seconds[float(shift)],
                        parent_assembly_seconds=parent_assembly_seconds,
                        offline_basis_seconds=offline_basis_seconds,
                        online_repeats=args.online_repeats,
                        geometry=args.geometry,
                        corner_radius=args.corner_radius,
                    )
                )

    return {
        "schema": "radia.validation.evrs_current_field_compare.v2",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "validation_host": platform.node(),
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "ngsolve_version": getattr(ng, "__version__", "unknown"),
        },
        "note": (
            "The parent solve uses the dimensionless (K + shift M) HCurl model. "
            "Timing fields are wall-clock observations on validation_host; publish "
            "them only when this file was generated on a designated compute host."
        ),
        "configuration": {
            "geometry": args.geometry,
            "orders": args.orders,
            "krylov_steps": args.steps,
            "reference_order": args.reference_order,
            "shifts": args.shifts,
            "maxh": args.maxh,
            "curve_order": args.curve_order,
            "corner_refinements": args.corner_refinements,
            "refinement_radius": args.refinement_radius,
            "corner_edge_maxh": args.corner_edge_maxh,
            "sigma": args.sigma,
            "intorder": args.intorder,
            "geometry_intorder": args.geometry_intorder,
            "corner_radius": args.corner_radius,
            "condense_from": args.condense_from,
            "inverse": args.inverse,
            "online_repeats": args.online_repeats,
            "rtol": args.rtol,
        },
        "geometry_diagnostics": geometry,
        "refinement_history": refinement_history,
        "total_wall_seconds": float(time.perf_counter() - started),
        "cases": cases,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--geometry",
        choices=("box", "bar", "notched-box", "l-prism", "sphere"),
        default="notched-box",
    )
    parser.add_argument("--orders", type=_parse_ints, default=[3, 4, 5, 6])
    parser.add_argument("--steps", type=_parse_ints, default=[4, 8, 12, 16, 22])
    parser.add_argument("--reference-order", type=int, default=6)
    parser.add_argument("--shifts", type=_parse_floats, default=[1.0])
    parser.add_argument("--maxh", type=float, default=2.0)
    parser.add_argument("--curve-order", type=int, default=1)
    parser.add_argument("--corner-refinements", type=int, default=0)
    parser.add_argument("--refinement-radius", type=float, default=0.18)
    parser.add_argument("--corner-edge-maxh", type=float, default=None)
    parser.add_argument("--sigma", type=float, default=5.8e7)
    parser.add_argument("--intorder", type=int, default=2)
    parser.add_argument("--geometry-intorder", type=int, default=8)
    parser.add_argument("--corner-radius", type=float, default=0.35)
    parser.add_argument("--condense-from", type=int, default=7)
    parser.add_argument("--inverse", default="sparsecholesky")
    parser.add_argument("--online-repeats", type=int, default=200)
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)

    if args.reference_order < 1:
        parser.error("--reference-order must be positive")
    if args.corner_radius <= 0.0:
        parser.error("--corner-radius must be positive")
    if args.curve_order < 1:
        parser.error("--curve-order must be positive")
    if args.corner_refinements < 0:
        parser.error("--corner-refinements must be non-negative")
    if args.refinement_radius <= 0.0:
        parser.error("--refinement-radius must be positive")
    if args.corner_edge_maxh is not None and args.corner_edge_maxh <= 0.0:
        parser.error("--corner-edge-maxh must be positive")
    if args.geometry_intorder < 0:
        parser.error("--geometry-intorder must be non-negative")
    if args.online_repeats < 1:
        parser.error("--online-repeats must be positive")
    if any(order < 1 for order in args.orders):
        parser.error("--orders must contain positive integers")

    result = run_compare(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    print("EVRS/VIM current-field comparison against HCurl-FEM")
    print(f"  output: {args.output}")
    print(
        "  kind      p   n  shift  rank/active  J_L2     loss     "
        "local_L2   local_peak  speedup"
    )
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
            f"{case['max_corner_peak_magnitude_relative_error']:.3e} "
            f"{case.get('online_speedup', 1.0):.2e}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
