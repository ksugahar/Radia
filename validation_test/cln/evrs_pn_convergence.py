"""EVRS p-by-Krylov-depth convergence smoke for high-order HCurl eddy bases.

This validation lane is intentionally small.  It is not a benchmark and the
numbers should not be used as timing claims.  The purpose is to keep the
research observable concrete:

    HCurl(p) parent -> EVRS_n -> curl(T) current basis -> reduced VIM
    -> port admittance Y_r(s) = B_r^* Z_r(s)^-1 B_r.

Run from the repository root:

    python validation_test/cln/evrs_pn_convergence.py
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
DEFAULT_OUTPUT = HERE / "evrs_pn_convergence_smoke.json"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import radia.vim as vim


def _parse_ints(text: str) -> list[int]:
    values = [int(part.strip()) for part in text.split(",") if part.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected at least one integer")
    if any(value < 1 for value in values):
        raise argparse.ArgumentTypeError("values must be positive")
    return values


def _complex_matrix_parts(matrix: np.ndarray) -> dict[str, list[list[float]]]:
    return {
        "real": np.asarray(matrix.real, dtype=float).tolist(),
        "imag": np.asarray(matrix.imag, dtype=float).tolist(),
    }


def _relative_frobenius_error(matrix: np.ndarray, reference: np.ndarray) -> float:
    denom = float(np.linalg.norm(reference))
    if denom == 0.0:
        return float(np.linalg.norm(matrix - reference))
    return float(np.linalg.norm(matrix - reference) / denom)


def _sampled_basis_diagnostics(basis: vim.SampledCurrentBasis) -> dict[str, object]:
    mass = basis.mass_matrix()
    min_eigenvalue = 0.0
    if mass.size:
        min_eigenvalue = float(np.min(np.linalg.eigvalsh(0.5 * (mass + mass.conj().T)).real))
    return {
        "kind": basis.kind,
        "modes": int(basis.n_modes),
        "samples": int(basis.n_samples),
        "mass_trace": float(np.trace(mass).real) if mass.size else 0.0,
        "min_mass_eigenvalue": min_eigenvalue,
    }


def _make_mesh(maxh: float):
    import ngsolve as ng
    import netgen.occ as occ

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("cond")
    return ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=maxh))


def _assemble_parent(order: int, maxh: float, condense: bool):
    import ngsolve as ng

    mesh = _make_mesh(maxh)
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

    return mesh, fes, stiffness, mass, ports


def _port_samples(points: np.ndarray) -> list[np.ndarray]:
    return [
        np.column_stack(
            (
                -points[:, 1],
                points[:, 0],
                np.zeros(points.shape[0]),
            )
        ),
        np.column_stack(
            (
                np.zeros(points.shape[0]),
                -points[:, 2],
                points[:, 1],
            )
        ),
    ]


def _response_basis(fes, stiffness, mass, ports, *, steps: int, condense: bool, rtol: float):
    if condense:
        return vim.NgsolveStaticCondensedBlockKrylovBasis(
            stiffness,
            mass,
            ports,
            steps=steps,
            free_dofs=fes.FreeDofs(True),
            rtol=rtol,
        )
    return vim.NgsolveBlockKrylovBasis(
        stiffness,
        mass,
        ports,
        steps=steps,
        free_dofs=fes.FreeDofs(False),
        rtol=rtol,
    )


def _case_result(
    mesh,
    fes,
    stiffness,
    mass,
    ports,
    *,
    order: int,
    steps: int,
    condense: bool,
    sigma: float,
    omega: float,
    intorder: int,
    kernel_epsilon: float,
    rtol: float,
) -> dict[str, object]:
    t0 = time.perf_counter()
    response = _response_basis(
        fes,
        stiffness,
        mass,
        ports,
        steps=steps,
        condense=condense,
        rtol=rtol,
    )
    response_seconds = time.perf_counter() - t0

    names = [f"EVRS_p{order}_n{steps}_{i}" for i in range(response.rank)]
    volume = vim.NgsolveHCurlCurlBasis(
        mesh,
        fes,
        response.vectors,
        intorder=intorder,
        materials="cond",
        names=names,
    )
    system = vim.AssembleHybridVIM(
        volume,
        sigma=sigma,
        kernel_epsilon=kernel_epsilon,
    )
    rhs = np.column_stack(
        [vim.ExternalVectorPotentialRHS(volume, samples) for samples in _port_samples(volume.points)]
    )
    admittance = vim.ReducedPortAdmittance(system, 1j * omega, rhs)
    info = response.diagnostics()
    topology = vim.ClassifyNgsolveEddyTopology(mesh, conductive_materials="cond")
    dof_policy = vim.NgsolveEddyDofPolicy(mesh, fes, topology)
    conductor_graph = topology.conductor_graph()
    bridge_cycle = vim.NgsolveBridgeCycleCurrentBasis(mesh, topology)
    bridge_cycle_system = vim.AssembleHybridVIM(
        bridge_cycle,
        sigma=sigma,
        kernel_epsilon=kernel_epsilon,
    )
    reduction_plan = dof_policy.reduction_plan(
        evrs_rank=response.rank,
        surface_modes=0,
        loop_bridge_modes=conductor_graph.cycle_rank,
        bridge_strategy="cycle-basis",
    )
    return {
        "order": order,
        "krylov_steps": steps,
        "condensed": condense,
        "ndof": int(info["ndof"]),
        "active_dofs": int(info["active_dofs"]),
        "rank": int(info["rank"]),
        "inactive_dofs": int(info["inactive_dofs"]),
        "eddy_visible_dofs": int(info["eddy_visible_dofs"]),
        "eddy_invisible_dofs": int(info["eddy_invisible_dofs"]),
        "compression_ratio": float(info["compression_ratio"]),
        "response_seconds_lab_smoke": float(response_seconds),
        "n_samples": volume.n_samples,
        "topology_diagnostics": topology.diagnostics(),
        "conductor_graph_diagnostics": conductor_graph.diagnostics(),
        "bridge_cycle_basis_diagnostics": _sampled_basis_diagnostics(bridge_cycle),
        "bridge_cycle_vim_diagnostics": bridge_cycle_system.diagnostics(),
        "dof_policy_diagnostics": dof_policy.diagnostics(),
        "reduction_plan_diagnostics": reduction_plan.diagnostics(),
        "system_diagnostics": system.diagnostics(),
        "admittance": _complex_matrix_parts(admittance),
        "_admittance_array": admittance,
    }


def run_sweep(args: argparse.Namespace) -> dict[str, object]:
    import ngsolve as ng

    rows: list[dict[str, object]] = []
    parent_cache = {}
    for order in args.orders:
        condense = order >= args.condense_from
        parent_cache[order] = _assemble_parent(order, args.maxh, condense)
        mesh, fes, stiffness, mass, ports = parent_cache[order]
        for steps in args.steps:
            rows.append(
                _case_result(
                    mesh,
                    fes,
                    stiffness,
                    mass,
                    ports,
                    order=order,
                    steps=steps,
                    condense=condense,
                    sigma=args.sigma,
                    omega=args.omega,
                    intorder=args.intorder,
                    kernel_epsilon=args.kernel_epsilon,
                    rtol=args.rtol,
                )
            )

    reference_order = max(args.orders)
    reference_steps = max(args.steps)
    reference = next(
        row["_admittance_array"]
        for row in rows
        if row["order"] == reference_order and row["krylov_steps"] == reference_steps
    )
    for row in rows:
        row["relative_error_to_reference"] = _relative_frobenius_error(
            row["_admittance_array"],
            reference,
        )
        del row["_admittance_array"]

    return {
        "schema": "radia.validation.evrs_pn_convergence.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "validation_host": platform.node(),
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "ngsolve_version": getattr(ng, "__version__", "unknown"),
        },
        "note": (
            "LAB/desktop smoke only: timings are recorded to diagnose the run, "
            "not as benchmark data."
        ),
        "configuration": {
            "orders": args.orders,
            "krylov_steps": args.steps,
            "reference": {"order": reference_order, "krylov_steps": reference_steps},
            "maxh": args.maxh,
            "sigma": args.sigma,
            "omega_rad_per_s": args.omega,
            "intorder": args.intorder,
            "kernel_epsilon": args.kernel_epsilon,
            "condense_from": args.condense_from,
            "rtol": args.rtol,
        },
        "rows": rows,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orders", type=_parse_ints, default=[2, 4, 6, 8])
    parser.add_argument("--steps", type=_parse_ints, default=[1, 2, 3])
    parser.add_argument("--maxh", type=float, default=3.0)
    parser.add_argument("--sigma", type=float, default=5.8e7)
    parser.add_argument("--omega", type=float, default=2.0 * np.pi * 10_000.0)
    parser.add_argument("--intorder", type=int, default=2)
    parser.add_argument("--kernel-epsilon", type=float, default=0.12)
    parser.add_argument("--condense-from", type=int, default=7)
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)

    result = run_sweep(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    print("EVRS p-by-n convergence smoke")
    print(f"  output: {args.output}")
    print("  p  n  rank/active       compression     relerr")
    for row in result["rows"]:
        print(
            f"  {row['order']:>1}  {row['krylov_steps']:>1}  "
            f"{row['rank']:>4}/{row['active_dofs']:<5}  "
            f"{row['compression_ratio']:.3e}  "
            f"{row['relative_error_to_reference']:.3e}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
