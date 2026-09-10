"""Reference solution for the ECB plate force of radia.maglev.ecb.lorentz.

``compute_lorentz_force_via_foster`` solves, in a truncated Dirichlet
eigenbasis (the Foster expansion),

    (-Delta + s mu sigma) v = -s mu sigma B_z,   v = 0 on "outer",

for the z component ``v`` of the reaction field, and builds the eddy current as
``J = (1/mu) curl(v z)``.  This lane solves the SAME model directly (no
eigen-truncation) on a mirror-symmetric structured hex mesh and reconstructs
the current the same way.  Force is <F> = 0.5 int Re(J) x B dV with the full
dipole field (B_x, B_y, B_z).

The reference must meet what a centred z-dipole over a plate has to show: zero
horizontal force (mirror symmetry), a repulsive lift (conductor Fz < 0), mesh
convergence, a lift rising with frequency, and a lift below the infinite
perfect-conductor image bound.  The shipped kernel is then held to the
reference: the same symmetry and sign, the same lift within 2.5 % at 50 and
500 Hz with 200 modes, and -- at 5 kHz, where 200 modes stop far short of the
skin depth -- monotone convergence to the reference as the basis grows.

Until 2026-09-11 the kernel built the current as -omega sigma Im(v).  On this
lane it gave zero lift for a centred magnet and a horizontal force of 1906 N at
5 kHz against the 11.7 N bound; ``history`` in the summary keeps that record.

Solver-heavy validation: run on an idle mdx or hibino host.  The recorded host
is written into the summary.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import socket
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "src" / "radia").is_dir())
if str(REPO / "src") not in sys.path:
    sys.path.insert(0, str(REPO / "src"))

import scipy.sparse as sp  # noqa: E402
import scipy.sparse.linalg as spla  # noqa: E402
from ngsolve import (H1, BilinearForm, GridFunction, Integrate,  # noqa: E402
                     LinearForm, Mesh, TaskManager, dx, grad)
from ngsolve import x as xC  # noqa: E402
from ngsolve import y as yC  # noqa: E402
from ngsolve import z as zC  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

import radia  # noqa: E402
from radia.maglev.ecb.lorentz import compute_lorentz_force_via_foster  # noqa: E402
from radia.maglev.mixed_galerkin.alpha import _dirichlet_eigenmodes  # noqa: E402

DEFAULT_OUTPUT = HERE / "ecb_foster_lorentz_reference_summary.json"
MU0 = 4.0e-7 * math.pi
PLATE_M = (0.20, 0.08, 0.010)          # top face at z = 0
SIGMA = 3.5e7
M_PM = 10.0                            # A m^2, z-oriented
Z_PM = 0.02
MESHES = ((30, 12, 3), (40, 16, 4))
FREQUENCIES_HZ = (50.0, 500.0, 5000.0)
OFFSETS_M = (0.0, 0.03)
KERNEL_EIGENMODES = 200
KERNEL_AGREEMENT_FREQUENCIES_HZ = (50.0, 500.0)
KERNEL_LIFT_RTOL = 0.025
MODE_STUDY_MESH = (30, 12, 3)
MODE_STUDY_FREQUENCY_HZ = 5000.0
MODE_STUDY_COUNTS = (200, 400, 800, 1600)
MODE_STUDY_FINAL_RTOL = 0.01

HISTORY = {
    "before": "2026-09-11",
    "current_reconstruction": "J_y = -omega sigma Im(v)",
    "drive_projection": "M_free @ Bz_free (boundary values of B_z dropped)",
    "centred_horizontal_force_N_50_500_5000Hz": [-31.57, -1194.09, -1906.27],
    "centred_vertical_force_N": 0.0,
    "image_bound_N": 11.72,
}


def plate_mesh(nx, ny, nz):
    lx, ly, lz = PLATE_M
    mesh = MakeStructured3DMesh(
        hexes=True, nx=nx, ny=ny, nz=nz,
        mapping=lambda x, y, z: (-lx / 2 + lx * x, -ly / 2 + ly * y, -lz + lz * z),
    )
    ngmesh = mesh.ngmesh
    for index in range(len(set(mesh.GetBoundaries()))):
        ngmesh.SetBCName(index, "outer")
    mesh = Mesh(ngmesh)
    if set(mesh.GetBoundaries()) != {"outer"}:
        raise RuntimeError(f"boundary relabel failed: {sorted(set(mesh.GetBoundaries()))}")
    return mesh


def dipole_field(x_pm):
    r2 = (xC - x_pm) ** 2 + yC ** 2 + (zC - Z_PM) ** 2 + 1e-30
    k = MU0 / (4 * math.pi)
    dz = zC - Z_PM
    bx = k * 3 * M_PM * (xC - x_pm) * dz / r2 ** 2.5
    by = k * 3 * M_PM * yC * dz / r2 ** 2.5
    bz = k * (3 * M_PM * dz * dz / r2 ** 2.5 - M_PM / r2 ** 1.5)
    return bx, by, bz


def _csr(form, ndof):
    rows, cols, vals = form.mat.COO()
    return sp.csr_matrix((np.asarray(vals), (np.asarray(rows), np.asarray(cols))),
                         shape=(ndof, ndof))


def solve_reaction_field(mesh, fes, bz, s):
    """Direct Galerkin solve of the kernel's scalar problem -- no eigen truncation."""
    smus = s * MU0 * SIGMA
    u, w = fes.TnT()
    a = BilinearForm(fes); a += grad(u) * grad(w) * dx; a.Assemble()
    m = BilinearForm(fes); m += u * w * dx; m.Assemble()
    f = LinearForm(fes); f += bz * w * dx; f.Assemble()
    free = np.array([fes.FreeDofs()[i] for i in range(fes.ndof)], dtype=bool)
    system = (_csr(a, fes.ndof) + smus * _csr(m, fes.ndof))[free][:, free].tocsc()
    v = np.zeros(fes.ndof, dtype=complex)
    v[free] = spla.spsolve(system, -smus * np.asarray(f.vec.FV().NumPy())[free])
    return v


def reference_force(mesh, fes, v, x_pm):
    """<F> on the conductor from J = (1/mu) curl(v z), full dipole B."""
    bx, by, bz = dipole_field(x_pm)
    v_re = GridFunction(fes)
    v_re.vec.FV().NumPy()[:] = v.real
    gradient = grad(v_re)
    jx = gradient[1] / MU0
    jy = -gradient[0] / MU0
    return [
        float(0.5 * Integrate(jy * bz, mesh)),
        float(0.5 * Integrate(-jx * bz, mesh)),
        float(0.5 * Integrate(jx * by - jy * bx, mesh)),
    ]


def mode_study():
    """Kernel lift error against the direct solve as the Foster basis grows."""
    mesh = plate_mesh(*MODE_STUDY_MESH)
    fes = H1(mesh, order=2, dirichlet="outer")
    s = 2j * math.pi * MODE_STUDY_FREQUENCY_HZ
    _bx, _by, bz = dipole_field(0.0)
    reference = reference_force(mesh, fes, solve_reaction_field(mesh, fes, bz, s), 0.0)[2]
    rows = []
    for count in MODE_STUDY_COUNTS:
        lam, vecs, _mass, free, _fes, _volume = _dirichlet_eigenmodes(mesh, count, "outer")
        lift = compute_lorentz_force_via_foster(
            mesh, lam, vecs, free, SIGMA, MU0, s, M_PM, Z_PM, 0.0)[2]
        rows.append({
            "modes": int(len(lam)),
            "lambda_max_per_m2": float(lam[-1]),
            "lambda_max_over_s_mu_sigma": float(lam[-1]) / abs(s * MU0 * SIGMA),
            "kernel_lift_N": float(lift),
            "relative_error": abs(float(lift) / reference - 1.0),
        })
        print(f"mode study {count:5d} modes: lambda_max/|s mu sigma| = "
              f"{rows[-1]['lambda_max_over_s_mu_sigma']:.3f}, lift error "
              f"{rows[-1]['relative_error']:.3e}", flush=True)
    return {"mesh": list(MODE_STUDY_MESH), "frequency_hz": MODE_STUDY_FREQUENCY_HZ,
            "reference_lift_N": reference, "rows": rows}


def run(output):
    lift_bound = 0.5 * 3 * MU0 * M_PM ** 2 / (32 * math.pi * Z_PM ** 4)
    cases = []
    with TaskManager():
        for nx, ny, nz in MESHES:
            mesh = plate_mesh(nx, ny, nz)
            fes = H1(mesh, order=2, dirichlet="outer")
            lam, vecs, _mass, free, _fes, _volume = _dirichlet_eigenmodes(
                mesh, KERNEL_EIGENMODES, "outer")
            for frequency in FREQUENCIES_HZ:
                s = 2j * math.pi * frequency
                for x_pm in OFFSETS_M:
                    _bx, _by, bz = dipole_field(x_pm)
                    v = solve_reaction_field(mesh, fes, bz, s)
                    reference = reference_force(mesh, fes, v, x_pm)
                    kernel = [float(value) for value in compute_lorentz_force_via_foster(
                        mesh, lam, vecs, free, SIGMA, MU0, s, M_PM, Z_PM, x_pm)]
                    cases.append({
                        "mesh": [nx, ny, nz],
                        "ndof": fes.ndof,
                        "frequency_hz": frequency,
                        "x_pm_m": x_pm,
                        "reference_force_N": reference,
                        "kernel_force_N": kernel,
                        "kernel_lift_relative_error": abs(kernel[2] / reference[2] - 1.0),
                    })
                    print(f"mesh={nx}x{ny}x{nz} f={frequency:7.1f} x_pm={x_pm:+.3f} | "
                          f"reference F=({reference[0]:+.3e},{reference[1]:+.3e},"
                          f"{reference[2]:+.3e}) | kernel F=({kernel[0]:+.3e},"
                          f"{kernel[1]:+.3e},{kernel[2]:+.3e})", flush=True)
        study = mode_study()

    def pick(mesh, frequency, x_pm):
        return next(c for c in cases if c["mesh"] == list(mesh)
                    and c["frequency_hz"] == frequency and c["x_pm_m"] == x_pm)

    coarse, fine = MESHES
    centred = [pick(fine, f, 0.0) for f in FREQUENCIES_HZ]
    lifts = [-c["reference_force_N"][2] for c in centred]
    mesh_change = [
        abs(pick(fine, f, 0.0)["reference_force_N"][2]
            / pick(coarse, f, 0.0)["reference_force_N"][2] - 1.0)
        for f in FREQUENCIES_HZ
    ]
    agreement = [c["kernel_lift_relative_error"] for c in cases
                 if c["frequency_hz"] in KERNEL_AGREEMENT_FREQUENCIES_HZ]
    study_errors = [row["relative_error"] for row in study["rows"]]

    def horizontal_vanishes(force):
        return max(abs(force[0]), abs(force[1])) < 1e-9 * abs(force[2])

    checks = {
        "reference_horizontal_force_vanishes_when_centred": all(
            horizontal_vanishes(c["reference_force_N"]) for c in centred),
        "reference_lift_repels_the_magnet": all(
            c["reference_force_N"][2] < 0.0 for c in cases),
        "reference_mesh_converged_within_1pct": max(mesh_change) < 0.01,
        "reference_lift_rises_with_frequency": all(
            a < b for a, b in zip(lifts, lifts[1:])),
        "reference_lift_below_image_bound": max(lifts) < lift_bound,
    }
    kernel_checks = {
        "kernel_horizontal_force_vanishes_when_centred": all(
            horizontal_vanishes(c["kernel_force_N"]) for c in centred),
        "kernel_lift_repels_the_magnet": all(
            c["kernel_force_N"][2] < 0.0 for c in cases),
        "kernel_lift_matches_reference_at_50_and_500Hz": max(agreement) < KERNEL_LIFT_RTOL,
        "kernel_converges_with_modes_at_5kHz": (
            all(a > b for a, b in zip(study_errors, study_errors[1:]))
            and study_errors[-1] < MODE_STUDY_FINAL_RTOL),
    }
    import ngsolve
    payload = {
        "schema": "radia.maglev.ecb-foster-lorentz-reference.v2",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "runtime": {
            "radia_version": radia.__version__,
            "ngsolve_version": ngsolve.__version__,
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "host": socket.gethostname(),
        },
        "problem": {
            "plate_m": list(PLATE_M), "sigma_S_per_m": SIGMA,
            "dipole_moment_Am2": M_PM, "dipole_height_m": Z_PM,
            "meshes": [list(m) for m in MESHES],
            "frequencies_hz": list(FREQUENCIES_HZ), "offsets_m": list(OFFSETS_M),
            "kernel_eigenmodes": KERNEL_EIGENMODES,
            "kernel_agreement_frequencies_hz": list(KERNEL_AGREEMENT_FREQUENCIES_HZ),
            "kernel_lift_rtol": KERNEL_LIFT_RTOL,
            "mode_study_final_rtol": MODE_STUDY_FINAL_RTOL,
            "image_lift_bound_N": lift_bound,
        },
        "cases": cases,
        "reference_mesh_relative_change": mesh_change,
        "mode_study": study,
        "checks": checks,
        "kernel_checks": kernel_checks,
        "reference_passed": all(checks.values()),
        "kernel_passed": all(kernel_checks.values()),
        "history": HISTORY,
    }
    output.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    passed = payload["reference_passed"] and payload["kernel_passed"]
    print(f"[{'OK' if passed else 'FAIL'}] kernel 50/500 Hz lift error "
          f"{max(agreement):.3e}; 5 kHz mode study {study_errors}; wrote {output}")
    return payload


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    payload = run(args.output)
    return 0 if payload["reference_passed"] and payload["kernel_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
