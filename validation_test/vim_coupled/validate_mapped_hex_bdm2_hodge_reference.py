"""Separate the mapped-HEX BDM2 space from its open-boundary charge kernel.

The same two-cell mapped HEX body is evaluated in two independent ways:

* a finite-domain H1 Omega formulation, ``N_omega = C.T K^-1 C``;
* the production-style volume/surface charge BEM diagnostic.

``N_omega`` is a discrete Hodge projection, so its generalized spectrum must
lie in ``[0, 1]`` relative to the HDiv mass.  A bounded H1 spectrum together
with an out-of-range charge spectrum localizes the defect to the mapped charge
operator rather than to the BDM2 approximation space.  The finite outer box is
not promoted as an open-boundary accuracy oracle.

No mesh artifact is tracked; both meshes are regenerated in memory.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import ngsolve as ng
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
from netgen.meshing import (
    Element2D,
    Element3D,
    FaceDescriptor,
    Mesh,
    MeshPoint,
    Pnt,
)
from ngsolve.meshes import MakeStructured3DMesh


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "src"
DEFAULT_OUTPUT = HERE / "results_mapped_hex_bdm2_hodge_reference.json"

if str(SRC) not in __import__("sys").path:
    __import__("sys").path.insert(0, str(SRC))

import radia  # noqa: E402
from radia import vim as vim_api  # noqa: E402
from radia.vim import _vim as vim_core  # noqa: E402


NX, NY, NZ = 6, 3, 3
BODY_I = {2, 3}
BODY_J = {1}
BODY_K = {1}


def _package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _git_head() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
    ).strip()


def _git_dirty() -> bool:
    return bool(
        subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=REPO, text=True
        ).strip()
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _source_fingerprints() -> dict[str, str]:
    paths = (
        Path(__file__).resolve(),
        SRC / "radia" / "vim" / "_vim.py",
        SRC / "radia" / "vim" / "_eddy_hybrid.py",
        SRC / "radia" / "_radia_pybind.pyd",
    )
    return {
        path.relative_to(REPO).as_posix(): _sha256(path)
        for path in paths
    }


def physical_map(x: float, y: float, z: float) -> tuple[float, float, float]:
    """Continuous trilinear-cell map with a strong non-affine body region."""
    base = 0.12 * (np.asarray((x, y, z), dtype=float) - 0.5)
    warp = np.asarray(
        (
            0.010 * np.sin(2.0 * np.pi * y) * np.sin(2.0 * np.pi * z),
            0.006 * np.sin(2.0 * np.pi * x) * np.sin(2.0 * np.pi * z),
            0.006 * np.sin(2.0 * np.pi * x) * np.sin(2.0 * np.pi * y),
        )
    )
    return tuple(float(value) for value in base + warp)


def make_air_body_mesh() -> ng.Mesh:
    """Build a conforming 54-HEX outer box with a central two-HEX body."""
    netmesh = Mesh(dim=3)
    netmesh.SetMaterial(1, "air")
    netmesh.SetMaterial(2, "body")
    points = []
    for i in range(NX + 1):
        for j in range(NY + 1):
            for k in range(NZ + 1):
                points.append(
                    netmesh.Add(
                        MeshPoint(Pnt(*physical_map(i / NX, j / NY, k / NZ)))
                    )
                )

    def offset(i: int, j: int, k: int) -> int:
        return i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k

    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                base = offset(i, j, k)
                baseup = base + (NY + 1) * (NZ + 1)
                local = (
                    base,
                    base + 1,
                    base + (NZ + 1) + 1,
                    base + (NZ + 1),
                    baseup,
                    baseup + 1,
                    baseup + (NZ + 1) + 1,
                    baseup + (NZ + 1),
                )
                material = (
                    2 if i in BODY_I and j in BODY_J and k in BODY_K else 1
                )
                netmesh.Add(Element3D(material, [points[index] for index in local]))

    netmesh.Add(FaceDescriptor(surfnr=0, domin=1, bc=1))

    def add_quad(indices: tuple[int, int, int, int]) -> None:
        netmesh.Add(Element2D(1, [points[index] for index in indices]))

    for j in range(NY):
        for k in range(NZ):
            add_quad(
                (
                    offset(0, j, k),
                    offset(0, j + 1, k),
                    offset(0, j + 1, k + 1),
                    offset(0, j, k + 1),
                )
            )
            add_quad(
                (
                    offset(NX, j, k),
                    offset(NX, j, k + 1),
                    offset(NX, j + 1, k + 1),
                    offset(NX, j + 1, k),
                )
            )
    for i in range(NX):
        for k in range(NZ):
            add_quad(
                (
                    offset(i, 0, k),
                    offset(i, 0, k + 1),
                    offset(i + 1, 0, k + 1),
                    offset(i + 1, 0, k),
                )
            )
            add_quad(
                (
                    offset(i, NY, k),
                    offset(i + 1, NY, k),
                    offset(i + 1, NY, k + 1),
                    offset(i, NY, k + 1),
                )
            )
    for i in range(NX):
        for j in range(NY):
            add_quad(
                (
                    offset(i, j, 0),
                    offset(i + 1, j, 0),
                    offset(i + 1, j + 1, 0),
                    offset(i, j + 1, 0),
                )
            )
            add_quad(
                (
                    offset(i, j, NZ),
                    offset(i, j + 1, NZ),
                    offset(i + 1, j + 1, NZ),
                    offset(i + 1, j, NZ),
                )
            )
    netmesh.SetBCName(0, "outer")
    netmesh.Compress()
    return ng.Mesh(netmesh)


def body_local_map(x: float, y: float, z: float) -> tuple[float, float, float]:
    return physical_map(
        (2.0 + 2.0 * x) / NX,
        (1.0 + y) / NY,
        (1.0 + z) / NZ,
    )


def make_body_mesh() -> ng.Mesh:
    return MakeStructured3DMesh(
        hexes=True,
        nx=2,
        ny=1,
        nz=1,
        mapping=body_local_map,
    )


def _generalized_spectrum(operator: np.ndarray, mass: np.ndarray) -> dict[str, object]:
    eigenvalues = sla.eigh(operator, mass, eigvals_only=True)
    return {
        "minimum": float(eigenvalues[0]),
        "maximum": float(eigenvalues[-1]),
        "outside_count": int(
            np.count_nonzero(
                (eigenvalues < -1.0e-8) | (eigenvalues > 1.0 + 1.0e-5)
            )
        ),
    }


def solve_h1_hodge_reference(mesh: ng.Mesh, *, h1_order: int) -> dict[str, object]:
    body = mesh.Materials("body")
    hdiv = ng.HDiv(mesh, order=2, definedon=body)
    h1 = ng.H1(mesh, order=h1_order, dirichlet="outer")
    hodge = vim_core.H1HodgeDemagOperator(
        hdiv,
        h1,
        definedon=body,
        boundary_contract="finite-dirichlet-validation",
    )

    active_hdiv = np.flatnonzero(np.asarray(hdiv.FreeDofs(), dtype=bool))
    free_h1 = np.flatnonzero(np.asarray(h1.FreeDofs(), dtype=bool))
    operator = np.asarray(hodge.mat.ToDense())[np.ix_(active_hdiv, active_hdiv)]
    operator = 0.5 * (operator + operator.T)
    mass_rows, mass_columns, mass_values = hodge.mass.COO()
    mass = sp.csr_matrix(
        (mass_values, (mass_rows, mass_columns)),
        shape=(hdiv.ndof, hdiv.ndof),
    )
    mass_dense = mass[active_hdiv, :][:, active_hdiv].toarray()
    return {
        "formulation": "finite-domain H1 Omega discrete Hodge projection",
        "backend": hodge.Diagnostics(),
        "h1_order": int(h1_order),
        "h1_dofs": int(h1.ndof),
        "h1_free_dofs": int(len(free_h1)),
        "hdiv_dofs": int(hdiv.ndof),
        "hdiv_active_dofs": int(len(active_hdiv)),
        "spectrum": _generalized_spectrum(operator, mass_dense),
    }


def solve_charge_bem_diagnostic(mesh: ng.Mesh) -> dict[str, object]:
    hdiv = ng.HDiv(mesh, order=2)
    charge_map, gram, mass, mass_ngsolve = vim_core._build_charge_gram_hex(
        hdiv,
        glout_n=5,
        glin_n=7,
        near_grade=0.5,
        far_inner=1.0,
        eps=1.0e-14,
        leafsize=4096,
    )
    _, gram, mass = vim_core._configure_cpp_operator(
        charge_map, gram, mass, mass_ngsolve
    )
    operator = np.asarray(gram.demag_matrix().ToDense())
    operator = 0.5 * (operator + operator.T)
    return {
        "formulation": "mapped volume/surface charge BEM diagnostic",
        "hdiv_dofs": int(hdiv.ndof),
        "spectrum": _generalized_spectrum(
            operator, sp.csr_matrix(mass).toarray()
        ),
        "affinity": vim_core._hex_mapping_affinity_report(mesh),
    }


def solve_h1_hcurl_mixed_reference(mesh: ng.Mesh) -> dict[str, object]:
    """Exercise mapped BDM2 response reduction inside the mixed solver.

    This is a shared-mesh mechanics and contraction gate.  Its sampled
    Laplace interaction and finite Dirichlet H1 exterior are deliberately not
    promoted as an open-boundary accuracy reference.
    """

    started = time.perf_counter()
    body = mesh.Materials("body")
    h1 = ng.H1(mesh, order=3, dirichlet="outer")
    hcurl = ng.HCurl(mesh, order=2, nograds=True)
    trial, test = hcurl.TnT()
    stiffness = ng.BilinearForm(hcurl)
    stiffness += (
        ng.curl(trial) * ng.curl(test) + 0.05 * trial * test
    ) * ng.dx
    mass = ng.BilinearForm(hcurl)
    mass += trial * test * ng.dx
    port = ng.LinearForm(hcurl)
    port += ng.CF((-ng.y, ng.x, 0.0)) * test * ng.dx("body")

    stiffness.Assemble()
    mass.Assemble()
    port.Assemble()
    mixed = vim_api.NgsolveBDMEddyBubbleVIM(
        mesh,
        hcurl,
        stiffness,
        mass,
        port,
        (),
        hdiv_order=2,
        mu_r=100.0,
        external_fields=(ng.CF((0.0, 0.0, 1.0)),),
        external_names=("uniform_Hz",),
        training_fields=(ng.CF((ng.x, 0.0, 0.0)),),
        training_names=("linear_Hx",),
        hdiv_max_modes=2,
        hdiv_intorder=4,
        hdiv_solve_tol=1.0e-10,
        magnetic_materials="body",
        hdiv_definedon=body,
        demag_operator_factory=lambda hdiv: vim_api.H1HodgeDemagOperator(
            hdiv,
            h1,
            definedon=body,
            boundary_contract="finite-dirichlet-validation",
        ),
        steps=2,
        sigma=1.0e6,
        conductive_materials="body",
        air_materials=("air",),
        response_backend="dense",
        intorder=3,
        kernel_epsilon=0.002,
        coupling_kernel_epsilon=0.002,
        interaction=vim_api.HACApKSampledLaplaceInteraction(
            mu=4.0e-7 * np.pi,
            kernel_epsilon=0.002,
            cross_only=False,
        ),
    )
    reduction = mixed.hdiv_reduction
    hdiv = reduction.fes
    solution = mixed.solve_frequency(100.0)
    generation = reduction.basis_generation
    reduced_spectrum = sla.eigh(
        0.5 * (reduction.demag + reduction.demag.T),
        0.5 * (reduction.mass + reduction.mass.T),
        eigvals_only=True,
    )
    return {
        "formulation": "finite-domain H1-Hodge BDM2 plus HCurl eddy-bubble",
        "validation_role": "shared-mesh mixed-mechanics gate, not open-boundary accuracy",
        "frequency_hz": 100.0,
        "sigma_s_per_m": 1.0e6,
        "mu_r": 100.0,
        "mesh_elements": int(mesh.ne),
        "hdiv_parent_dofs": int(hdiv.ndof),
        "hdiv_active_dofs": int(sum(hdiv.FreeDofs())),
        "hdiv_modes": int(reduction.n_modes),
        "hcurl_parent_dofs": int(hcurl.ndof),
        "hcurl_modes": int(mixed.n_hcurl_vim_modes),
        "demag_backend": reduction.diagnostics()["demag_backend"],
        "snapshot_backend": generation["snapshot_backend"],
        "snapshot_iterations": generation["snapshot_iterations"],
        "max_snapshot_relative_residual": float(
            generation["max_snapshot_relative_residual"]
        ),
        "reduced_demag_generalized_spectrum": {
            "minimum": float(reduced_spectrum[0]),
            "maximum": float(reduced_spectrum[-1]),
        },
        "hcurl_interaction_backend": mixed.eddy_system.interaction_backend,
        "mixed_relative_residual": float(solution.residual_relative_norm),
        "average_joule_loss": float(solution.average_joule_loss[0]),
        "magnetization_coefficient_norm": float(
            np.linalg.norm(solution.magnetization_coefficients)
        ),
        "eddy_coefficient_norm": float(np.linalg.norm(solution.eddy_coefficients)),
        "wall_s": time.perf_counter() - started,
    }


def run() -> dict[str, object]:
    started = time.perf_counter()
    source_start = _source_fingerprints()
    head_start = _git_head()
    outer = make_air_body_mesh()
    body = make_body_mesh()
    with ng.TaskManager():
        hodge = [
            solve_h1_hodge_reference(outer, h1_order=order)
            for order in (2, 3)
        ]
        mixed = solve_h1_hcurl_mixed_reference(outer)
        charge = solve_charge_bem_diagnostic(body)
    source_end = _source_fingerprints()
    hodge_spectra = [row["spectrum"] for row in hodge]
    checks = {
        "source_fingerprints_stable_during_run": source_start == source_end,
        "body_is_nonaffine": charge["affinity"]["nonaffine_cell_count"] == 2,
        "h1_hodge_active_dofs_match_body_bdm2": all(
            row["hdiv_active_dofs"] == charge["hdiv_dofs"] for row in hodge
        ),
        "h1_hodge_spectra_are_contractions": all(
            row["minimum"] >= -1.0e-8
            and row["maximum"] <= 1.0 + 1.0e-5
            and row["outside_count"] == 0
            for row in hodge_spectra
        ),
        "h1_hodge_public_operator_path_verified": all(
            row["backend"]["operator"] == "C.T @ K^-1 @ C"
            and row["backend"]["boundary_contract"]
            == "finite-dirichlet-validation"
            for row in hodge
        ),
        "charge_bem_spectrum_violation_reproduced": (
            charge["spectrum"]["maximum"] > 1.1
            and charge["spectrum"]["outside_count"] > 0
        ),
        "mapped_bdm2_response_uses_generic_ngsolve_cg": (
            mixed["demag_backend"] == "H1HodgeDemagOperator"
            and mixed["snapshot_backend"] == "ngsolve-mass-preconditioned-cg"
            and mixed["max_snapshot_relative_residual"] < 1.0e-8
        ),
        "mapped_bdm2_hcurl_mixed_solve_is_finite_and_converged": (
            mixed["hdiv_active_dofs"] == 207
            and mixed["reduced_demag_generalized_spectrum"]["minimum"]
            >= -1.0e-8
            and mixed["reduced_demag_generalized_spectrum"]["maximum"]
            <= 1.0 + 1.0e-5
            and mixed["mixed_relative_residual"] < 1.0e-10
            and mixed["average_joule_loss"] > 0.0
            and mixed["magnetization_coefficient_norm"] > 0.0
            and mixed["eddy_coefficient_norm"] > 0.0
        ),
    }
    return {
        "schema": "radia.validation.mapped-hex-bdm2-hodge-reference.v2",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "tool_versions": {
            "radia": getattr(radia, "__version__", _package_version("radia-ngsolve")),
            "radia_source_head": head_start,
            "radia_source_head_end": _git_head(),
            "radia_source_dirty": _git_dirty(),
            "python": platform.python_version(),
            "ngsolve": _package_version("ngsolve"),
            "numpy": _package_version("numpy"),
            "scipy": _package_version("scipy"),
        },
        "source_fingerprints": source_start,
        "identity": {
            "outer_mesh": "54 mapped HEX, air plus central body",
            "body_mesh": "same central two mapped HEX regenerated independently",
            "hdiv_family": "BDM2",
            "body_active_dofs": 207,
            "meshes_tracked": False,
        },
        "h1_hodge_reference": hodge,
        "h1_hcurl_mixed_reference": mixed,
        "charge_bem_diagnostic": charge,
        "checks": checks,
        "pass": all(checks.values()),
        "timing_s": {"total": time.perf_counter() - started},
        "claim_boundary": {
            "established": (
                "the mapped BDM2 HDiv space and NGSolve mass/coupling admit a "
                "bounded discrete Hodge projection; the out-of-range spectrum "
                "is localized to the current mapped volume/surface charge operator; "
                "the same BDM2 response reduction feeds a converged shared-mesh "
                "HCurl eddy-bubble mixed solve"
            ),
            "not_established": (
                "an open-boundary H1 accuracy reference, a repaired composite "
                "charge operator, a production open-boundary replacement for "
                "vim.Solve, or universal solver superiority"
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "pass": result["pass"], "checks": result["checks"]}, indent=2))
    return 0 if result["pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
