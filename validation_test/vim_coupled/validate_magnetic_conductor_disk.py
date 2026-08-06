"""Cross-validate a thin magnetic conductor without a tracked mesh artifact.

The quick profile checks four independent responsibilities:

* axisymmetric Q2 ``radia.axifem`` reference;
* full 3-D HCurl A-form reference;
* mapped-HEX BDM1 HDiv-MMM static response;
* HDiv-MMM plus HCurl eddy-bubble mixed mechanics and frequency routing.

With ``--include-bdm2-gate``, the runner also verifies that the known-unsafe
mapped/non-affine HEX BDM2 material lane fails loudly.  The legacy
``--include-bdm2-negative`` spelling remains an alias for replay compatibility.

The full profile adds the expensive HCurl p/h and BDM1 h ladders used by the
curated adjudication artifact.  Cubit always runs headless, and generated
``.vol`` files remain under ``C:\\temp``.

Run from the repository root:

    python validation_test/vim_coupled/validate_magnetic_conductor_disk.py
    python validation_test/vim_coupled/validate_magnetic_conductor_disk.py --profile full
    python validation_test/vim_coupled/validate_magnetic_conductor_disk.py --direct-h-ladder
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import netgen.occ as occ
import ngsolve as ng
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "src"
AXIFEM_VERIFY = REPO / "validation_test" / "axifem" / "research" / "verification"
GENERATOR = HERE / "build_magnetic_conductor_disk_hex.py"
DEFAULT_OUTPUT = Path(r"C:\temp\magnetic_conductor_disk_live.json")

for path in (SRC, AXIFEM_VERIFY):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import radia  # noqa: E402
from radia import vim  # noqa: E402
from radia.axifem import (  # noqa: E402
    AxiHenrotteSigmaMassBFI,
    AxiHenrotteStiffnessBFI,
    H1Henrotte,
)
import test_hiruma_disk_q1 as axifem_disk  # noqa: E402


MU0 = 4.0e-7 * math.pi
RADIUS_M = 0.010
THICKNESS_M = 0.0005
MU_R = 100.0
SIGMA_S_PER_M = 1.0e7
H0_A_PER_M = 1.0


def _package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _git_head() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO,
        check=False,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip() if completed.returncode == 0 else "unknown"


def _git_dirty() -> bool:
    completed = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=REPO,
        check=False,
        capture_output=True,
        text=True,
    )
    return completed.returncode != 0 or bool(completed.stdout.strip())


def _source_fingerprints() -> dict[str, str]:
    paths = (
        SRC / "radia" / "__init__.py",
        SRC / "radia" / "axifem.pyd",
        SRC / "radia" / "_radia_pybind.pyd",
        SRC / "radia" / "vim" / "_eddy_hybrid.py",
        SRC / "radia" / "vim" / "_hcurl_tet_interaction.py",
        Path(__file__).resolve(),
        GENERATOR,
    )
    return {
        path.relative_to(REPO).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in paths
    }


def _symmetric_csr(matrix, size: int) -> sp.csr_matrix:
    rows, columns, values = matrix.COO()
    result = sp.csr_matrix(
        (np.asarray(values, dtype=float), (rows, columns)),
        shape=(size, size),
    )
    return 0.5 * (result + result.T)


def solve_axisymmetric_q2(
    frequency_hz: float,
    *,
    fine: bool,
) -> dict[str, object]:
    started = time.perf_counter()
    mesh_parameters = (
        (160, 32, 32, 32) if fine else (80, 16, 24, 24)
    )
    nr_disk, nz_disk, nr_air, nz_air = mesh_parameters
    old_thickness = axifem_disk.T_DISK
    old_sigma = axifem_disk.SIGMA_CU
    try:
        axifem_disk.T_DISK = THICKNESS_M
        axifem_disk.SIGMA_CU = SIGMA_S_PER_M
        mesh = axifem_disk.make_structured_disk_quad_mesh(
            nr_disk, nz_disk, nr_air, nz_air, 0.5, 0.5
        )
    finally:
        axifem_disk.T_DISK = old_thickness
        axifem_disk.SIGMA_CU = old_sigma

    fes = H1Henrotte(mesh, order=2, dirichlet="axis|right|top|bot")
    permeability = mesh.MaterialCF(
        {"conductor": MU0 * MU_R},
        default=MU0,
    )
    conductivity = mesh.MaterialCF(
        {"conductor": SIGMA_S_PER_M},
        default=0.0,
    )
    stiffness = ng.BilinearForm(fes, symmetric=True)
    stiffness += AxiHenrotteStiffnessBFI(permeability)
    mass = ng.BilinearForm(fes, symmetric=True)
    mass += AxiHenrotteSigmaMassBFI(conductivity)
    with ng.TaskManager():
        stiffness.Assemble()
        mass.Assemble()

    system = _symmetric_csr(stiffness.mat, fes.ndof)
    system = system + 2j * math.pi * frequency_hz * _symmetric_csr(
        mass.mat, fes.ndof
    )
    free = np.asarray(
        [index for index in range(fes.ndof) if fes.FreeDofs()[index]],
        dtype=int,
    )
    boundary = np.asarray(
        [index for index in range(fes.ndof) if not fes.FreeDofs()[index]],
        dtype=int,
    )
    values = np.zeros(fes.ndof, dtype=complex)
    prescribed = ng.GridFunction(fes)
    prescribed.Set(
        0.5 * MU0 * H0_A_PER_M * ng.x,
        definedon=mesh.Boundaries("axis|right|top|bot"),
    )
    values[boundary] = prescribed.vec.FV().NumPy()[boundary]
    values[free] = spla.spsolve(
        system[free[:, None], free[None, :]],
        -system[free[:, None], boundary[None, :]] @ values[boundary],
    )

    conductor = mesh.Materials("conductor")
    volume = float(ng.Integrate(2.0 * math.pi * ng.x, mesh, definedon=conductor))
    averages = []
    for component in (values.real, values.imag):
        grid_function = ng.GridFunction(fes)
        grid_function.vec.FV().NumPy()[:] = component
        bz = ng.grad(grid_function)[0] + grid_function / ng.x
        averages.append(
            float(
                ng.Integrate(
                    2.0 * math.pi * ng.x * bz,
                    mesh,
                    definedon=conductor,
                    order=4,
                )
                / volume
            )
        )
    normalized = complex(*averages) / (MU0 * H0_A_PER_M)
    return {
        "frequency_hz": float(frequency_hz),
        "normalized_Bz": [float(normalized.real), float(normalized.imag)],
        "elements": int(mesh.ne),
        "dofs": int(fes.ndof),
        "mesh_parameters": {
            "nr_disk": nr_disk,
            "nz_disk": nz_disk,
            "nr_air": nr_air,
            "nz_air": nz_air,
            "outer_radius_m": 0.5,
            "outer_half_height_m": 0.5,
        },
        "wall_s": time.perf_counter() - started,
    }


def _full_3d_mesh(disk_maxh_m: float) -> ng.Mesh:
    disk = occ.Cylinder(
        occ.Pnt(0.0, 0.0, -0.5 * THICKNESS_M),
        occ.Z,
        RADIUS_M,
        THICKNESS_M,
    )
    disk.mat("conductor")
    disk.maxh = disk_maxh_m
    for face in disk.faces:
        face.name = "interface"
    box = occ.Box(
        occ.Pnt(-0.04, -0.04, -0.02),
        occ.Pnt(0.04, 0.04, 0.02),
    )
    box.mat("air")
    box.maxh = 0.010
    for face in box.faces:
        face.name = "outer"
    air = box - disk
    air.mat("air")
    return ng.Mesh(occ.OCCGeometry(occ.Glue([air, disk])).GenerateMesh(maxh=0.010))


def solve_full_3d_hcurl(
    frequency_hz: float,
    *,
    order: int,
    disk_maxh_m: float,
) -> dict[str, object]:
    started = time.perf_counter()
    mesh = _full_3d_mesh(disk_maxh_m)
    mesh_s = time.perf_counter() - started
    fes = ng.HCurl(
        mesh,
        order=order,
        dirichlet="outer",
        complex=True,
        nograds=True,
    )
    trial, test = fes.TnT()
    inv_mu = mesh.MaterialCF(
        {"conductor": 1.0 / (MU0 * MU_R)},
        default=1.0 / MU0,
    )
    sigma = mesh.MaterialCF({"conductor": SIGMA_S_PER_M}, default=0.0)
    form = ng.BilinearForm(fes, symmetric=False)
    form += inv_mu * ng.curl(trial) * ng.curl(test) * ng.dx
    form += (
        1j * 2.0 * math.pi * frequency_hz * sigma * trial * test * ng.dx
    )
    with ng.TaskManager():
        form.Assemble()
    assembly_s = time.perf_counter() - started - mesh_s

    external_a = 0.5 * MU0 * H0_A_PER_M * ng.CoefficientFunction(
        (-ng.y, ng.x, 0.0)
    )
    solution = ng.GridFunction(fes)
    solution.Set(external_a, definedon=mesh.Boundaries("outer"))
    initial_residual = -form.mat * solution.vec
    solve_started = time.perf_counter()
    with ng.TaskManager():
        inverse = form.mat.Inverse(fes.FreeDofs(), inverse="umfpack")
        solution.vec.data += inverse * initial_residual
    solve_s = time.perf_counter() - solve_started

    conductor = mesh.Materials("conductor")
    volume = float(ng.Integrate(1.0, mesh, definedon=conductor))
    average_bz = complex(
        ng.Integrate(
            ng.curl(solution)[2],
            mesh,
            definedon=conductor,
            order=4,
        )
        / volume
    )
    algebraic_residual = solution.vec.CreateVector()
    algebraic_residual.data = form.mat * solution.vec
    for index, is_free in enumerate(fes.FreeDofs()):
        if not is_free:
            algebraic_residual[index] = 0.0
    relative_residual = float(
        ng.Norm(algebraic_residual) / max(ng.Norm(initial_residual), 1.0e-300)
    )
    normalized = average_bz / (MU0 * H0_A_PER_M)
    return {
        "frequency_hz": float(frequency_hz),
        "order": int(order),
        "disk_maxh_m": float(disk_maxh_m),
        "tetrahedra": int(mesh.ne),
        "dofs": int(fes.ndof),
        "normalized_Bz": [float(normalized.real), float(normalized.imag)],
        "relative_residual": relative_residual,
        "timing_s": {
            "mesh": mesh_s,
            "assembly": assembly_s,
            "solve": solve_s,
            "total": time.perf_counter() - started,
        },
    }


def _generate_hex(size_mm: float, tag: str) -> tuple[Path, dict[str, object]]:
    output = Path(r"C:\temp") / f"magnetic_conductor_disk_{tag}.vol"
    command = [
        sys.executable,
        str(GENERATOR),
        "--size-mm",
        f"{size_mm:g}",
        "--output",
        str(output),
    ]
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=REPO,
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=120.0,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"headless Cubit mesh generation failed ({completed.returncode})"
        )
    metadata = None
    for line in reversed(completed.stdout.splitlines()):
        try:
            candidate = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(candidate, dict) and candidate.get("output"):
            metadata = candidate
            break
    if metadata is None:
        raise RuntimeError("headless Cubit output did not contain mesh metadata")
    metadata["generation_wall_s"] = time.perf_counter() - started
    metadata["headless"] = True
    return output, metadata


def solve_hdiv_static(mesh_path: Path, *, order: int) -> dict[str, object]:
    started = time.perf_counter()
    mesh = ng.Mesh(str(mesh_path))
    with ng.TaskManager():
        result = vim.Solve(
            mesh,
            mu_r=MU_R,
            H_ext=ng.CoefficientFunction((0.0, 0.0, H0_A_PER_M)),
            order=order,
        )
    normalized = MU_R / (MU_R - 1.0) * float(result["M_avg"][2]) / H0_A_PER_M
    return {
        "order": int(order),
        "cell_family": "hex",
        "elements": int(result["n_el"]),
        "dofs": int(result["ndof"]),
        "normalized_Bz": float(normalized),
        "demag": float(result["demag"]),
        "iterations": int(result["iters"]),
        "timing_s": {
            "charge_basis": float(result["charge_basis_wall_s"]),
            "charge_gram": float(result["charge_gram_wall_s"]),
            "solve": float(result["solve_wall_s"]),
            "total": time.perf_counter() - started,
        },
    }


def probe_mapped_hex_bdm2_material_gate(mesh_path: Path) -> dict[str, object]:
    """Record the expected production rejection without promoting bad physics."""
    started = time.perf_counter()
    try:
        unexpected_result = solve_hdiv_static(mesh_path, order=2)
    except NotImplementedError as exc:
        message = str(exc)
        expected = (
            "mapped/non-affine HEX BDM2 material solve" in message
            and "separate volume/surface charge quadrature" in message
            and "mapped HEX BDM1" in message
        )
        return {
            "status": "rejected_as_expected" if expected else "unexpected_rejection",
            "expected_gate": bool(expected),
            "error_type": type(exc).__name__,
            "message": message,
            "wall_s": time.perf_counter() - started,
        }
    return {
        "status": "unexpectedly_accepted",
        "expected_gate": False,
        "result": unexpected_result,
        "wall_s": time.perf_counter() - started,
    }


def solve_coupled_smoke(
    mesh_path: Path,
    *,
    maxh_m: float,
    frequency_hz: float = 100.0,
    interaction_mode: str = "sampled",
) -> dict[str, object]:
    started = time.perf_counter()
    interaction_mode = str(interaction_mode).strip().lower()
    if interaction_mode not in {"sampled", "direct-q2-hex"}:
        raise ValueError("interaction_mode must be 'sampled' or 'direct-q2-hex'")
    mesh = ng.Mesh(str(mesh_path))
    fes = ng.HCurl(mesh, order=1, nograds=True)
    trial, test = fes.TnT()
    stiffness = ng.BilinearForm(fes)
    stiffness += (ng.curl(trial) * ng.curl(test) + 100.0 * trial * test) * ng.dx
    mass = ng.BilinearForm(fes)
    mass += trial * test * ng.dx
    base = 0.5 * MU0 * H0_A_PER_M * ng.CoefficientFunction(
        (-ng.y, ng.x, 0.0)
    )
    radial = (ng.x * ng.x + ng.y * ng.y) / (RADIUS_M * RADIUS_M)
    axial = (2.0 * ng.z / THICKNESS_M) ** 2
    ports = tuple(
        vim.NgsolveHCurlVectorPotentialPort(fes, potential, materials="conductor")
        for potential in (base, radial * base, radial * radial * base, axial * base)
    )
    with ng.TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        for port in ports:
            port.Assemble()

    def external_a(points):
        values = np.asarray(points)
        return 0.5 * MU0 * H0_A_PER_M * np.column_stack(
            (-values[:, 1], values[:, 0], np.zeros(values.shape[0]))
        )

    gate = vim.EddySIBCApplicability(
        frequency_hz=frequency_hz,
        sigma=SIGMA_S_PER_M,
        characteristic_thickness_m=THICKNESS_M,
        mu=MU0 * MU_R,
    )
    interaction_options = {}
    if interaction_mode == "sampled":
        interaction_options["interaction"] = vim.HACApKSampledLaplaceInteraction(
            mu=MU0,
            kernel_epsilon=0.25 * maxh_m,
            cross_only=False,
        )
    with ng.TaskManager():
        mixed = vim.NgsolveBDMEddyBubbleVIM(
            mesh,
            fes,
            stiffness,
            mass,
            ports,
            (
                ng.CoefficientFunction((1.0, 0.0, 0.0)),
                ng.CoefficientFunction((0.0, 1.0, 0.0)),
                ng.CoefficientFunction((0.0, 0.0, 1.0)),
            ),
            hdiv_order=1,
            mu_r=MU_R,
            external_fields=(ng.CoefficientFunction((0.0, 0.0, H0_A_PER_M)),),
            external_names=("uniform_Hz",),
            training_fields=vim.NgsolveHDivRegularSolidHarmonicPorts(
                mesh, max_degree=2
            ),
            hdiv_pod_rtol=1.0e-10,
            hdiv_max_modes=12,
            hdiv_solve_tol=1.0e-10,
            hdiv_intorder=3,
            magnetic_materials="conductor",
            demag_eps=0.25 * maxh_m,
            demag_eta=2.0,
            steps=3,
            sigma=SIGMA_S_PER_M,
            conductive_materials="conductor",
            surface_boundaries="|".join(mesh.GetBoundaries()),
            intorder=4,
            kernel_epsilon=0.25 * maxh_m,
            response_backend="operator",
            inverse="sparsecholesky",
            rtol=1.0e-11,
            current_gram_rtol=1.0e-11,
            parent_order=1,
            parent_order_ledger=vim.EddyParentOrderLedger(
                bulk_degree=1,
                bridge_trace_degree=0,
                surface_current_degree=0,
            ),
            sibc_applicability=gate,
            port_vector_potentials=(external_a,),
            coupling_kernel_epsilon=0.25 * maxh_m,
            hcurl_interaction_max_subtets=2048,
            **interaction_options,
        )
    solution = mixed.solve_frequency(
        frequency_hz,
        solver="dense",
        surface_impedance=vim.SkinImpedance(
            2j * math.pi * frequency_hz,
            SIGMA_S_PER_M,
            MU0 * MU_R,
        ),
    )
    stale_error = None
    try:
        mixed.solve_frequency(10_000.0, solver="dense", surface_impedance=0.0)
    except ValueError as exc:
        stale_error = str(exc)
    weights = np.asarray(mixed.magnetization_basis.weights)
    mz = np.asarray(solution.sampled_magnetization)[0, :, 2]
    average_mz = np.sum(weights * mz) / np.sum(weights)
    normalized = MU_R / (MU_R - 1.0) * average_mz / H0_A_PER_M
    eddy_diagnostics = mixed.eddy_system.diagnostics()
    diagonal_backends = sorted(
        {
            str(value)
            for value in _nested_values(
                eddy_diagnostics["inductance_operator"],
                "scalar_gram_backend",
            )
        }
    )
    return {
        "interaction": (
            "explicit sampled HACApK smoke; not the accuracy oracle"
            if interaction_mode == "sampled"
            else "production auto HCurl diagonal with direct Q2 HEX reference density"
        ),
        "interaction_mode": interaction_mode,
        "interaction_backend": eddy_diagnostics["interaction_backend"],
        "inductance_matrix_free": eddy_diagnostics["inductance_matrix_free"],
        "inductance_operator": eddy_diagnostics["inductance_operator"],
        "hcurl_diagonal_backends": diagonal_backends,
        "frequency_hz": float(frequency_hz),
        "route": mixed.sibc_applicability.selected_model,
        "normalized_Bz": [float(normalized.real), float(normalized.imag)],
        "hexes": int(mesh.ne),
        "hcurl_parent_dofs": int(fes.ndof),
        "hdiv_parent_dofs": int(mixed.hdiv_reduction.parent_ndof),
        "hdiv_modes": int(mixed.n_hdiv_mmm_modes),
        "eddy_modes": int(mixed.n_hcurl_vim_modes),
        "relative_residual": float(solution.residual_relative_norm),
        "joule_loss": float(solution.average_joule_loss[0]),
        "stale_10khz_route_rejected": stale_error is not None,
        "stale_route_error": stale_error,
        "wall_s": time.perf_counter() - started,
    }


def _complex(row: dict[str, object]) -> complex:
    values = row["normalized_Bz"]
    return complex(float(values[0]), float(values[1]))


def _nested_values(value, key: str) -> list[object]:
    found = []
    if isinstance(value, dict):
        if key in value:
            found.append(value[key])
        for child in value.values():
            found.extend(_nested_values(child, key))
    elif isinstance(value, (list, tuple)):
        for child in value:
            found.extend(_nested_values(child, key))
    return found


def run(
    profile: str,
    include_bdm2_negative: bool,
    coupled_h_ladder: bool,
    direct_h_ladder: bool = False,
) -> dict[str, object]:
    started = time.perf_counter()
    source_start = {
        "head": _git_head(),
        "dirty": _git_dirty(),
        "fingerprints": _source_fingerprints(),
    }
    frequencies = (100.0,) if profile == "quick" else (100.0, 10_000.0)
    fine_axisymmetric_reference = (
        profile == "full" or coupled_h_ladder or direct_h_ladder
    )
    axisymmetric = [
        solve_axisymmetric_q2(
            frequency,
            fine=fine_axisymmetric_reference,
        )
        for frequency in frequencies
    ]
    reference_by_frequency = {
        row["frequency_hz"]: _complex(row) for row in axisymmetric
    }
    hcurl_cases = (
        ((100.0, 2, 0.002),)
        if profile == "quick"
        else (
            (100.0, 2, 0.002),
            (10_000.0, 1, 0.002),
            (10_000.0, 2, 0.002),
            (10_000.0, 3, 0.002),
            (10_000.0, 3, 0.001),
        )
    )
    hcurl = []
    for frequency_hz, order, maxh_m in hcurl_cases:
        row = solve_full_3d_hcurl(
            frequency_hz,
            order=order,
            disk_maxh_m=maxh_m,
        )
        row["reference_relative_error"] = float(
            abs(_complex(row) - reference_by_frequency[frequency_hz])
            / abs(reference_by_frequency[frequency_hz])
        )
        hcurl.append(row)

    sizes_mm = (4.0,) if profile == "quick" else (2.0, 1.0, 0.5)
    bdm1 = []
    meshes = {}
    mesh_metadata_by_size = {}
    for size_mm in sizes_mm:
        tag = f"{profile}_h{size_mm:g}".replace(".", "p")
        mesh_path, mesh_metadata = _generate_hex(size_mm, tag)
        meshes[size_mm] = mesh_path
        mesh_metadata_by_size[size_mm] = mesh_metadata
        row = solve_hdiv_static(mesh_path, order=1)
        row["size_mm"] = size_mm
        row["mesh_generation"] = mesh_metadata
        row["reference_relative_error"] = float(
            abs(
                complex(float(row["normalized_Bz"]), 0.0)
                - reference_by_frequency[100.0]
            )
            / abs(reference_by_frequency[100.0])
        )
        bdm1.append(row)

    bdm2 = []
    if include_bdm2_negative:
        for size_mm in (2.0, 1.0):
            mesh_path = meshes.get(size_mm)
            if mesh_path is None:
                mesh_path, _ = _generate_hex(size_mm, f"bdm2_h{size_mm:g}")
            row = probe_mapped_hex_bdm2_material_gate(mesh_path)
            row["size_mm"] = size_mm
            bdm2.append(row)

    coupled_sizes_mm = (
        (4.0, 2.0, 1.0)
        if profile == "full" or coupled_h_ladder or direct_h_ladder
        else (4.0,)
    )
    coupled_ladder = []
    for size_mm in coupled_sizes_mm:
        coupled_mesh = meshes.get(size_mm)
        coupled_metadata = mesh_metadata_by_size.get(size_mm)
        if coupled_mesh is None or coupled_metadata is None:
            tag = f"{profile}_coupled_h{size_mm:g}".replace(".", "p")
            coupled_mesh, coupled_metadata = _generate_hex(size_mm, tag)
        row = solve_coupled_smoke(
            coupled_mesh,
            maxh_m=size_mm * 1.0e-3,
        )
        row["size_mm"] = size_mm
        row["validation_role"] = (
            "h_refinement_ladder_member"
            if len(coupled_sizes_mm) > 1
            else "single_point_mechanics_smoke"
        )
        row["mesh_generation"] = coupled_metadata
        row["reference_relative_error"] = float(
            abs(_complex(row) - reference_by_frequency[100.0])
            / abs(reference_by_frequency[100.0])
        )
        coupled_ladder.append(row)
    coupled = coupled_ladder[0]
    direct_ladder = []
    direct_count = len(coupled_ladder) if direct_h_ladder else 1
    for sampled_row in coupled_ladder[:direct_count]:
        direct_row = solve_coupled_smoke(
            Path(sampled_row["mesh_generation"]["output"]),
            maxh_m=float(sampled_row["size_mm"]) * 1.0e-3,
            interaction_mode="direct-q2-hex",
        )
        direct_row["size_mm"] = sampled_row["size_mm"]
        direct_row["mesh_generation"] = sampled_row["mesh_generation"]
        direct_row["reference_relative_error"] = float(
            abs(_complex(direct_row) - reference_by_frequency[100.0])
            / abs(reference_by_frequency[100.0])
        )
        direct_row["sampled_observable_relative_difference"] = float(
            abs(_complex(direct_row) - _complex(sampled_row))
            / abs(_complex(sampled_row))
        )
        direct_row["sampled_joule_relative_difference"] = float(
            abs(direct_row["joule_loss"] - sampled_row["joule_loss"])
            / max(abs(sampled_row["joule_loss"]), np.finfo(float).tiny)
        )
        direct_ladder.append(direct_row)
    direct_coupled = direct_ladder[0]

    p_rows = [
        row
        for row in hcurl
        if row["frequency_hz"] == 10_000.0 and row["disk_maxh_m"] == 0.002
    ]
    checks = {
        "full_3d_hcurl_residuals_below_1e-8": max(
            row["relative_residual"] for row in hcurl
        )
        < 1.0e-8,
        "full_3d_hcurl_100hz_error_below_2pct": hcurl[0][
            "reference_relative_error"
        ]
        < 0.02,
        "coupled_residuals_below_1e-10": max(
            row["relative_residual"] for row in coupled_ladder
        )
        < 1.0e-10,
        "coupled_joule_nonnegative": min(
            row["joule_loss"] for row in coupled_ladder
        )
        >= 0.0,
        "stale_frequency_route_rejected": all(
            row["stale_10khz_route_rejected"] for row in coupled_ladder
        ),
        "production_direct_hex_residual_below_1e-10": max(
            row["relative_residual"] for row in direct_ladder
        )
        < 1.0e-10,
        "production_direct_hex_joule_nonnegative": min(
            row["joule_loss"] for row in direct_ladder
        )
        >= 0.0,
        "production_direct_hex_backend_observed": all(
            "direct-q2-hex-reference-density" in row["hcurl_diagonal_backends"]
            for row in direct_ladder
        ),
        "production_direct_hex_matches_sampled_observable_below_1e-4": max(
            row["sampled_observable_relative_difference"] for row in direct_ladder
        )
        < 1.0e-4,
        "production_direct_hex_matches_sampled_joule_below_1e-4": max(
            row["sampled_joule_relative_difference"] for row in direct_ladder
        )
        < 1.0e-4,
    }
    if profile == "full":
        checks.update(
            {
                "full_3d_hcurl_10khz_p_error_decreases": [
                    row["reference_relative_error"] for row in p_rows
                ]
                == sorted(
                    (row["reference_relative_error"] for row in p_rows),
                    reverse=True,
                ),
                "full_3d_hcurl_10khz_fine_p3_error_below_2pct": hcurl[-1][
                    "reference_relative_error"
                ]
                < 0.02,
                "bdm1_hex_h_error_decreases": [
                    row["reference_relative_error"] for row in bdm1
                ]
                == sorted(
                    (row["reference_relative_error"] for row in bdm1),
                    reverse=True,
                ),
                "bdm1_hex_fine_error_below_2pct": bdm1[-1][
                    "reference_relative_error"
                ]
                < 0.02,
            }
        )
    if include_bdm2_negative:
        checks["mapped_hex_bdm2_material_gate_verified"] = all(
            row["status"] == "rejected_as_expected" and row["expected_gate"]
            for row in bdm2
        )
    if len(coupled_ladder) > 1:
        coupled_errors = [
            row["reference_relative_error"] for row in coupled_ladder
        ]
        checks["coupled_sampled_h_error_decreases"] = coupled_errors == sorted(
            coupled_errors,
            reverse=True,
        )
        checks["coupled_sampled_fine_error_below_2pct"] = (
            coupled_errors[-1] < 0.02
        )
    if len(direct_ladder) > 1:
        direct_errors = [row["reference_relative_error"] for row in direct_ladder]
        checks["production_direct_hex_h_error_decreases"] = (
            direct_errors == sorted(direct_errors, reverse=True)
        )
        checks["production_direct_hex_fine_error_below_2pct"] = (
            direct_errors[-1] < 0.02
        )
    source_end = {
        "head": _git_head(),
        "dirty": _git_dirty(),
        "fingerprints": _source_fingerprints(),
    }
    source_head_stable = source_start["head"] == source_end["head"]
    checks["source_fingerprints_stable_during_run"] = (
        source_start["fingerprints"] == source_end["fingerprints"]
    )

    return {
        "schema": "cae-ai-lab.solver-run.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "tool_versions": {
            "radia": getattr(radia, "__version__", _package_version("radia-ngsolve")),
            "radia_source_head": source_start["head"],
            "radia_source_head_end": source_end["head"],
            "radia_source_dirty": source_start["dirty"],
            "python": platform.python_version(),
            "ngsolve": _package_version("ngsolve"),
            "numpy": _package_version("numpy"),
            "scipy": _package_version("scipy"),
        },
        "profile": profile,
        "source_state": {
            "start": source_start,
            "end": source_end,
            "head_stable_during_run": source_head_stable,
            "execution_fingerprints_stable_during_run": checks[
                "source_fingerprints_stable_during_run"
            ],
            "identity_policy": (
                "exact hashes of the validation driver, mesh generator, Python "
                "solver modules, and loaded native extensions are mandatory; "
                "the repository HEAD is diagnostic because unrelated commits may "
                "land in the shared worktree during a long validation"
            ),
        },
        "identity": {
            "geometry": "solid circular disk",
            "radius_m": RADIUS_M,
            "thickness_m": THICKNESS_M,
            "mu_r": MU_R,
            "sigma_S_per_m": SIGMA_S_PER_M,
            "excitation_Hz_A_per_m": H0_A_PER_M,
            "observable": "volume-average complex Bz divided by mu0*H0",
        },
        "axisymmetric_q2_reference": axisymmetric,
        "full_3d_hcurl_A_form": hcurl,
        "hdiv_mmm_static_hex_bdm1": bdm1,
        "mapped_hex_bdm2_material_gate": bdm2,
        "hdiv_mmm_static_hex_bdm2_negative": bdm2,
        "coupled_hdiv_mmm_hcurl_eddy_bubble_smoke": coupled,
        "coupled_hdiv_mmm_hcurl_eddy_bubble_h_ladder": coupled_ladder,
        "coupled_hdiv_mmm_hcurl_eddy_bubble_direct_hex": direct_coupled,
        "coupled_hdiv_mmm_hcurl_eddy_bubble_direct_hex_h_ladder": direct_ladder,
        "checks": checks,
        "pass": all(checks.values()),
        "timing_s": {"total": time.perf_counter() - started},
        "claim_boundary": {
            "established": (
                "cross-formulation references, BDM1 h convergence in full profile, "
                "sampled coupled h convergence when requested, production direct-Q2 "
                "HEX execution/parity"
                + (
                    " and h convergence"
                    if len(direct_ladder) > 1
                    else ""
                )
                + ", mixed mechanics/routing"
                + (
                    ", and mapped-HEX BDM2 fail-loud material gate"
                    if include_bdm2_negative
                    else ""
                )
            ),
            "not_established": (
                "a cancellation-preserving mapped-HEX BDM2 material operator or universal solver superiority"
                if len(direct_ladder) > 1
                else "a direct-Q2 HEX h-convergence ladder, a cancellation-preserving mapped-HEX BDM2 material operator, or universal solver superiority"
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=("quick", "full"), default="quick")
    parser.add_argument(
        "--include-bdm2-gate",
        "--include-bdm2-negative",
        dest="include_bdm2_negative",
        action="store_true",
        help="verify that mapped/non-affine HEX BDM2 material solves fail loudly",
    )
    parser.add_argument("--coupled-h-ladder", action="store_true")
    parser.add_argument("--direct-h-ladder", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run(
        args.profile,
        args.include_bdm2_negative,
        args.coupled_h_ladder,
        args.direct_h_ladder,
    )
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    validation_passed = bool(result["pass"])
    hcurl_rows = result["full_3d_hcurl_A_form"]
    bdm1_rows = result["hdiv_mmm_static_hex_bdm1"]
    coupled_rows = result["coupled_hdiv_mmm_hcurl_eddy_bubble_h_ladder"]
    direct_rows = result[
        "coupled_hdiv_mmm_hcurl_eddy_bubble_direct_hex_h_ladder"
    ]
    accepted_relative_errors = [
        float(hcurl_rows[0]["reference_relative_error"]),
    ]
    if args.profile == "full":
        accepted_relative_errors.extend(
            (
                float(hcurl_rows[-1]["reference_relative_error"]),
                float(bdm1_rows[-1]["reference_relative_error"]),
            )
        )
    if len(coupled_rows) > 1:
        accepted_relative_errors.append(
            float(coupled_rows[-1]["reference_relative_error"])
        )
    if len(direct_rows) > 1:
        accepted_relative_errors.append(
            float(direct_rows[-1]["reference_relative_error"])
        )
    command = (
        "python validation_test/vim_coupled/validate_magnetic_conductor_disk.py "
        f"--profile {args.profile} --output {output}"
    )
    if args.include_bdm2_negative:
        command += " --include-bdm2-gate"
    if args.coupled_h_ladder:
        command += " --coupled-h-ladder"
    if args.direct_h_ladder:
        command += " --direct-h-ladder"
    result.update(
        {
            "case": "thin magnetic-conductor HDiv-MMM plus HCurl validation",
            "solver": "radia-ngsolve",
            "source_artifact": str(GENERATOR),
            "run": {
                "command": command,
                "workdir": str(REPO),
                "exit_code": 0 if validation_passed else 2,
                "duration_s": float(result["timing_s"]["total"]),
            },
            "result_files": [str(output)],
            "tolerances": {
                "max_rel": 0.02,
                "max_abs": 1.0e-8,
                "full_3d_reference_max_rel": 0.02,
                "algebraic_residual_max_abs": 1.0e-8,
            },
            "errors": {
                "max_rel": max(accepted_relative_errors),
                "max_abs": float(
                    max(
                        *(row["relative_residual"] for row in hcurl_rows),
                        *(row["relative_residual"] for row in coupled_rows),
                        *(row["relative_residual"] for row in direct_rows),
                    )
                ),
                "all_hcurl_levels_max_rel": float(
                    max(row["reference_relative_error"] for row in hcurl_rows)
                ),
                "all_coupled_levels_max_rel": float(
                    max(
                        row["reference_relative_error"]
                        for row in coupled_rows
                    )
                ),
                "all_direct_hex_levels_max_rel": float(
                    max(row["reference_relative_error"] for row in direct_rows)
                ),
                "direct_sampled_observable_max_rel": float(
                    max(
                        row["sampled_observable_relative_difference"]
                        for row in direct_rows
                    )
                ),
                "direct_sampled_joule_max_rel": float(
                    max(
                        row["sampled_joule_relative_difference"]
                        for row in direct_rows
                    )
                ),
            },
            "timing_breakdown_s": {
                "axisymmetric_reference": float(
                    sum(
                        row["wall_s"]
                        for row in result["axisymmetric_q2_reference"]
                    )
                ),
                "full_3d_hcurl": float(
                    sum(row["timing_s"]["total"] for row in hcurl_rows)
                ),
                "hdiv_mesh_and_static": float(
                    sum(
                        row["mesh_generation"]["generation_wall_s"]
                        + row["timing_s"]["total"]
                        for row in bdm1_rows
                    )
                ),
                "coupled_mesh_and_solve": float(
                    sum(
                        row["mesh_generation"]["generation_wall_s"]
                        + row["wall_s"]
                        for row in coupled_rows
                    )
                ),
                "coupled_direct_hex": float(
                    sum(row["wall_s"] for row in direct_rows)
                ),
            },
            "verification": {
                "method": (
                    "axisymmetric Q2 versus full 3-D HCurl reference, "
                    "HDiv mapped-HEX response, direct-Q2 HEX versus sampled "
                    "coupled parity, mixed residual/passivity, and frequency-route "
                    "rejection"
                ),
                "command": command,
            },
        }
    )
    result["checks"].update(
        {
            "ran_to_completion": True,
            "result_files_exist": True,
            "validation_passed": validation_passed,
        }
    )
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(output),
                "pass": result["pass"],
                "checks": result["checks"],
                "timing_s": result["timing_s"],
            },
            indent=2,
        )
    )
    return 0 if result["pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
