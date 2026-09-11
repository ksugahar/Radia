"""Independent 3-D HCurl-VIM reference for the ECB thin-plate force model.

The production ``radia.maglev.ecb.lorentz`` kernel restricts the reaction
magnetic flux density to ``v z_hat`` and reconstructs an in-plane current from
its curl.  This lane instead builds a three-dimensional divergence-free HCurl
current basis, couples it with the open-boundary Newton-potential VIM operator,
and solves the frequency-domain current response to the same point-dipole
vector potential.  The two routes share only the incident field and the common
sampled Lorentz integrator.

The default ``smoke`` profile is intentionally small enough for development.
The ``full`` profile adds an h-refinement row and is the release evidence lane.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import socket
import sys
import time
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
REPO = next(path for path in HERE.parents if (path / "src" / "radia").is_dir())
if str(REPO / "src") not in sys.path:
    sys.path.insert(0, str(REPO / "src"))

import radia
from radia import vim
from radia.maglev.ecb import compute_lorentz_force_torque_via_hcurl_vim

DEFAULT_OUTPUT = HERE / "ecb_foster_lorentz_3d_reference_summary.json"
MU0 = 4.0e-7 * math.pi
PLATE_M = (0.10, 0.06, 0.002)
SIGMA = 3.5e7
M_PM = 10.0
Z_PM = 0.020
FREQUENCIES_HZ = (50.0, 500.0, 5000.0)
THREE_D_CONVERGENCE_RTOL = 0.005
SCALAR_MESH_CONVERGENCE_RTOL = 0.001
SCALAR_REJECTION_MIN_RELATIVE_DIFFERENCE = 0.5


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _package_version(name):
    try:
        return version(name)
    except PackageNotFoundError:
        return "unknown"


def _lift_relative_differences(left, right):
    return [
        abs(left_case["force_N"][2] / right_case["force_N"][2] - 1.0)
        for left_case, right_case in zip(left["cases"], right["cases"])
    ]


def dipole_fields(points, *, x_pm=0.0):
    """Return analytic peak-phasor ``(A, B)`` samples for a z dipole."""

    points = np.asarray(points, dtype=float)
    relative = points - np.asarray((x_pm, 0.0, Z_PM))
    radius2 = np.einsum("ij,ij->i", relative, relative)
    radius = np.sqrt(radius2)
    if np.any(radius == 0.0):
        raise ValueError("dipole source must remain outside the conductor")
    factor = MU0 * M_PM / (4.0 * math.pi)
    vector_potential = factor * np.column_stack(
        (-relative[:, 1], relative[:, 0], np.zeros(len(points)))
    ) / radius[:, None] ** 3
    moment = np.asarray((0.0, 0.0, M_PM))
    magnetic_flux_density = MU0 / (4.0 * math.pi) * (
        3.0
        * relative
        * (relative @ moment)[:, None]
        / radius[:, None] ** 5
        - moment[None, :] / radius[:, None] ** 3
    )
    return vector_potential, magnetic_flux_density


def _dipole_vector_potential_cf(ng, *, x_pm=0.0):
    relative = ng.CF((ng.x - x_pm, ng.y, ng.z - Z_PM))
    radius2 = relative * relative + 1.0e-30
    factor = MU0 * M_PM / (4.0 * math.pi)
    return factor * ng.CF((-relative[1], relative[0], 0.0)) / radius2 ** 1.5


def _dipole_flux_density_cf(ng, *, x_pm=0.0):
    relative = ng.CF((ng.x - x_pm, ng.y, ng.z - Z_PM))
    radius2 = relative * relative + 1.0e-30
    dot = M_PM * relative[2]
    return MU0 / (4.0 * math.pi) * (
        3.0 * dot * relative / radius2 ** 2.5
        - ng.CF((0.0, 0.0, M_PM)) / radius2 ** 1.5
    )


def _ngsolve_csr(form, size):
    rows, columns, values = form.mat.COO()
    return sp.csr_matrix(
        (np.asarray(values), (np.asarray(rows), np.asarray(columns))),
        shape=(size, size),
    )


def run_thin_plate_direct(
    *, mesh_counts=(24, 14, 2), x_pm=0.0, dirichlet_label="back|left|front|right"
):
    """Solve the scalar thin-plate PDE directly, without a Foster truncation."""

    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh

    started = time.perf_counter()
    lx, ly, lz = PLATE_M
    nx, ny, nz = mesh_counts
    base = MakeStructured3DMesh(
        hexes=True,
        nx=nx,
        ny=ny,
        nz=nz,
        mapping=lambda x, y, z: (
            -0.5 * lx + lx * x,
            -0.5 * ly + ly * y,
            -lz + lz * z,
        ),
    )
    mesh = ng.Mesh(base.ngmesh)
    fes = ng.H1(mesh, order=2, dirichlet=dirichlet_label)
    trial, test = fes.TnT()
    stiffness_form = ng.BilinearForm(fes, symmetric=True)
    stiffness_form += ng.grad(trial) * ng.grad(test) * ng.dx
    mass_form = ng.BilinearForm(fes, symmetric=True)
    mass_form += trial * test * ng.dx
    source_b = _dipole_flux_density_cf(ng, x_pm=x_pm)
    drive_form = ng.LinearForm(fes)
    drive_form += source_b[2] * test * ng.dx
    with ng.TaskManager():
        stiffness_form.Assemble()
        mass_form.Assemble()
        drive_form.Assemble()
    stiffness = _ngsolve_csr(stiffness_form, fes.ndof)
    mass = _ngsolve_csr(mass_form, fes.ndof)
    free = np.asarray(
        [bool(fes.FreeDofs()[index]) for index in range(fes.ndof)], dtype=bool
    )
    drive = np.asarray(drive_form.vec.FV().NumPy())
    cases = []
    for frequency_hz in FREQUENCIES_HZ:
        s = 2.0j * math.pi * frequency_hz
        s_mu_sigma = s * MU0 * SIGMA
        system = (stiffness + s_mu_sigma * mass)[free][:, free].tocsc()
        values = np.zeros(fes.ndof, dtype=complex)
        values[free] = spla.spsolve(system, -s_mu_sigma * drive[free])
        real_field = ng.GridFunction(fes)
        real_field.vec.FV().NumPy()[:] = values.real
        gradient = ng.grad(real_field)
        current = ng.CF((gradient[1] / MU0, -gradient[0] / MU0, 0.0))
        density = ng.CF(
            (
                current[1] * source_b[2],
                -current[0] * source_b[2],
                current[0] * source_b[1] - current[1] * source_b[0],
            )
        )
        force = [
            float(0.5 * ng.Integrate(density[index], mesh, order=6))
            for index in range(3)
        ]
        cases.append({"frequency_hz": frequency_hz, "force_N": force})
    return {
        "mesh_counts": list(mesh_counts),
        "dirichlet_label": dirichlet_label,
        "elements": int(mesh.ne),
        "dofs": int(fes.ndof),
        "cases": cases,
        "elapsed_seconds": time.perf_counter() - started,
    }


def _build_hcurl_vim(*, maxh_m, order, steps, outer_quad, x_pm=0.0):
    import ngsolve as ng
    from netgen import occ

    lx, ly, lz = PLATE_M
    plate = occ.Box(
        occ.Pnt(-0.5 * lx, -0.5 * ly, -lz),
        occ.Pnt(0.5 * lx, 0.5 * ly, 0.0),
    )
    plate.mat("plate")
    plate.maxh = maxh_m
    for face in plate.faces:
        face.name = "plate_air"
    mesh = ng.Mesh(occ.OCCGeometry(plate).GenerateMesh(maxh=maxh_m))
    fes = ng.HCurl(mesh, order=order, nograds=True)
    trial, test = fes.TnT()
    stiffness = ng.BilinearForm(fes)
    stiffness += (ng.curl(trial) * ng.curl(test) + 100.0 * trial * test) * ng.dx
    metric = ng.BilinearForm(fes)
    metric += trial * test * ng.dx
    source = _dipole_vector_potential_cf(ng, x_pm=x_pm)
    ports = (
        vim.NgsolveHCurlVectorPotentialPort(fes, source, materials="plate"),
    )
    with ng.TaskManager():
        stiffness.Assemble()
        metric.Assemble()
        for port in ports:
            port.Assemble()
        basis = vim.NgsolveEddyBubbleHCurlBasis(
            mesh,
            fes,
            stiffness,
            metric,
            ports,
            steps=steps,
            conductive_materials="plate",
            volume_materials="plate",
            intorder=max(order + 2, 5),
            parent_order=order,
            current_gram_rtol=1.0e-10,
        )
        interaction = basis.tet_volume_interaction(
            mesh,
            fes,
            degree=max(order - 1, 0),
            projection_quad=max(order + 2, 5),
            outer_quad=outer_quad,
            projection_tolerance=1.0e-9,
            materials="plate",
        )
    system = basis.assemble_vim(sigma=SIGMA, interaction=interaction)
    external_a, external_b = dipole_fields(basis.current_basis.points, x_pm=x_pm)
    rhs = vim.ExternalVectorPotentialRHS(basis.current_basis, external_a)
    model = vim.HCurlEddyCLNFromVIM(system, rhs)
    return ng, mesh, fes, basis, interaction, model, external_b


def run_3d_case(*, maxh_m=0.012, order=3, steps=3, outer_quad=5, x_pm=0.0):
    """Solve the independent 3-D current response and integrate Lorentz force."""

    started = time.perf_counter()
    _ng, mesh, fes, basis, interaction, model, external_b = _build_hcurl_vim(
        maxh_m=maxh_m,
        order=order,
        steps=steps,
        outer_quad=outer_quad,
        x_pm=x_pm,
    )
    cases = []
    for frequency_hz in FREQUENCIES_HZ:
        s = 2.0j * math.pi * frequency_hz
        force_raw, torque_raw = compute_lorentz_force_torque_via_hcurl_vim(
                model,
                basis.current_basis,
                external_b,
                s,
        )
        force = np.asarray(force_raw, dtype=float)
        torque = np.asarray(torque_raw, dtype=float)
        cases.append(
            {
                "frequency_hz": frequency_hz,
                "force_N": np.asarray(force, dtype=float).tolist(),
                "torque_Nm": torque.tolist(),
                "horizontal_to_vertical_ratio": float(
                    np.linalg.norm(force[:2])
                    / max(abs(force[2]), np.finfo(float).tiny)
                ),
                "torque_to_force_length_ratio": float(
                    np.linalg.norm(torque)
                    / max(
                        abs(force[2]) * max(PLATE_M[:2]),
                        np.finfo(float).tiny,
                    )
                ),
            }
        )
    return {
        "maxh_m": float(maxh_m),
        "order": int(order),
        "steps": int(steps),
        "outer_quad": int(outer_quad),
        "mesh_elements": int(mesh.ne),
        "parent_ndof": int(fes.ndof),
        "reduced_rank": int(basis.rank),
        "sample_count": int(basis.current_basis.n_samples),
        "cases": cases,
        "passive": bool(model.diagnostics()["passive"]),
        "interaction": interaction.diagnostics(),
        "elapsed_seconds": time.perf_counter() - started,
    }


def run(profile="smoke"):
    # The first two rows isolate reduced-rank convergence on the same mesh.
    # The full lane adds a finer mesh at the converged rank.
    print("3-D HCurl-VIM: coarse mesh, three reduction steps", flush=True)
    rows = [run_3d_case(maxh_m=0.012, order=3, steps=3, outer_quad=5)]
    print("3-D HCurl-VIM: coarse mesh, six reduction steps", flush=True)
    rows.append(run_3d_case(maxh_m=0.012, order=3, steps=6, outer_quad=5))
    if profile == "full":
        print("3-D HCurl-VIM: refined mesh, six reduction steps", flush=True)
        rows.append(run_3d_case(maxh_m=0.008, order=3, steps=6, outer_quad=6))

    # The legacy row reproduces the released all-boundary Dirichlet model.  The
    # lateral-only rows give the thin-plate ansatz its more favourable boundary
    # treatment; the comparison still rejects it quantitatively.
    print("Scalar direct solve: legacy all-boundary condition", flush=True)
    legacy_thin = run_thin_plate_direct(dirichlet_label=".*")
    print("Scalar direct solve: lateral boundary, coarse mesh", flush=True)
    thin_rows = [run_thin_plate_direct()]
    if profile == "full":
        print("Scalar direct solve: lateral boundary, refined mesh", flush=True)
        thin_rows.append(run_thin_plate_direct(mesh_counts=(36, 22, 3)))

    reference_3d = rows[-1]
    reference_thin = thin_rows[-1]
    comparisons = []
    for three_d, thin in zip(reference_3d["cases"], reference_thin["cases"]):
        comparisons.append(
            {
                "frequency_hz": three_d["frequency_hz"],
                "three_d_force_N": three_d["force_N"],
                "thin_plate_force_N": thin["force_N"],
                "lift_relative_difference": abs(
                    thin["force_N"][2] / three_d["force_N"][2] - 1.0
                ),
            }
        )
    rank_errors = _lift_relative_differences(rows[0], rows[1])
    mesh_errors = (
        _lift_relative_differences(rows[1], rows[2])
        if profile == "full"
        else None
    )
    scalar_mesh_errors = (
        _lift_relative_differences(thin_rows[0], thin_rows[1])
        if profile == "full"
        else None
    )
    scalar_model_errors = [row["lift_relative_difference"] for row in comparisons]
    legacy_model_errors = _lift_relative_differences(legacy_thin, reference_3d)
    image_bound = 0.5 * 3.0 * MU0 * M_PM**2 / (32.0 * math.pi * Z_PM**4)
    checks = {
        "all_models_passive": all(row["passive"] for row in rows),
        "centred_horizontal_force_small": all(
            case["horizontal_to_vertical_ratio"] < 0.01
            for row in rows
            for case in row["cases"]
        ),
        "lift_is_repulsive": all(
            case["force_N"][2] < 0.0 for row in rows for case in row["cases"]
        ),
        "centred_torque_small": all(
            case["torque_to_force_length_ratio"] < 0.01
            for row in rows
            for case in row["cases"]
        ),
        "thin_plate_lift_is_repulsive": all(
            case["force_N"][2] < 0.0
            for row in thin_rows
            for case in row["cases"]
        ),
        "three_d_rank_converged": max(rank_errors) < THREE_D_CONVERGENCE_RTOL,
        "three_d_mesh_converged": (
            mesh_errors is None or max(mesh_errors) < THREE_D_CONVERGENCE_RTOL
        ),
        "scalar_mesh_converged": (
            scalar_mesh_errors is None
            or max(scalar_mesh_errors) < SCALAR_MESH_CONVERGENCE_RTOL
        ),
        "three_d_lift_below_image_bound": all(
            -case["force_N"][2] < image_bound
            for row in rows
            for case in row["cases"]
        ),
        "scalar_model_rejected_for_quantitative_3d_force": min(
            scalar_model_errors
        ) > SCALAR_REJECTION_MIN_RELATIVE_DIFFERENCE,
        "legacy_all_boundary_model_rejected": min(
            legacy_model_errors
        ) > SCALAR_REJECTION_MIN_RELATIVE_DIFFERENCE,
    }
    payload = {
        "schema": "radia.maglev.ecb-foster-lorentz-3d-reference.v1",
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "runtime": {
            "host": socket.gethostname(),
            "python_version": platform.python_version(),
            "radia_version": getattr(radia, "__version__", "unknown"),
            "ngsolve_version": _package_version("ngsolve"),
            "numpy_version": _package_version("numpy"),
            "scipy_version": _package_version("scipy"),
        },
        "source_sha256": {
            "validation_script": _sha256(__file__),
            "lorentz_module": _sha256(REPO / "src" / "radia" / "maglev" / "ecb" / "lorentz.py"),
        },
        "problem": {
            "plate_m": list(PLATE_M),
            "sigma_S_per_m": SIGMA,
            "dipole_moment_Am2": M_PM,
            "dipole_height_m": Z_PM,
            "frequencies_hz": list(FREQUENCIES_HZ),
            "peak_phasor_image_lift_bound_N": image_bound,
        },
        "profile": profile,
        "rows": rows,
        "legacy_all_boundary_scalar_row": legacy_thin,
        "thin_plate_direct_rows": thin_rows,
        "comparisons": comparisons,
        "convergence": {
            "three_d_rank_lift_relative_differences": rank_errors,
            "three_d_mesh_lift_relative_differences": mesh_errors,
            "scalar_mesh_lift_relative_differences": scalar_mesh_errors,
            "three_d_rtol": THREE_D_CONVERGENCE_RTOL,
            "scalar_mesh_rtol": SCALAR_MESH_CONVERGENCE_RTOL,
        },
        "conclusion": {
            "scalar_model_status": "rejected-for-quantitative-3d-force",
            "recommended_backend": "hcurl_vim",
            "foster_scope": "fast-reduced-solve-of-the-scalar-local-reaction-model-only",
            "minimum_scalar_vs_3d_lift_relative_difference": min(
                scalar_model_errors
            ),
            "experiment_status": "not-run-no-measurement-dataset-supplied",
        },
        "checks": checks,
        "pass": all(checks.values()),
    }
    return payload


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--profile", choices=("smoke", "full"), default="smoke")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    payload = run(args.profile)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "pass": payload["pass"],
                      "rows": payload["rows"]}, indent=2))
    return 0 if payload["pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
