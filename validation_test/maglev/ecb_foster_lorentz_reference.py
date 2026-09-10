"""Reference solution for the ECB plate force of radia.maglev.ecb.lorentz.

``compute_lorentz_force_via_foster`` had no caller in the repository and had
never been run end to end.  This lane solves the SAME scalar model it uses,

    (-Delta + s mu sigma) v = -s mu sigma B_z,   v = 0 on "outer",

directly (no eigen-expansion) on a mirror-symmetric structured hex mesh, and
reconstructs the eddy current two ways:

  curl    J = (1/mu) curl(v z) = (1/mu) (d_y v, -d_x v, 0)   -- the reference
  kernel  J = (0, -omega sigma Im v, 0)                       -- what lorentz.py does

``v`` is the z component of the reaction field (tesla), so only the curl route
gives a current density in A/m^2 with the parity of a real eddy current.  Force
is <F> = 0.5 int Re(J) x B dV with the full dipole field (B_x, B_y, B_z).

A centred z-dipole over the plate must show: zero horizontal force (mirror
symmetry), a repulsive lift (conductor force Fz < 0), mesh convergence, a lift
that rises with frequency, and a lift below the infinite perfect-conductor
image bound.  The shipped kernel is run alongside and its violations recorded.

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
                    kernel_fx, kernel_fz = compute_lorentz_force_via_foster(
                        mesh, lam, vecs, free, SIGMA, MU0, s, M_PM, Z_PM, x_pm)
                    cases.append({
                        "mesh": [nx, ny, nz],
                        "ndof": fes.ndof,
                        "frequency_hz": frequency,
                        "x_pm_m": x_pm,
                        "reference_force_N": reference,
                        "kernel_force_xz_N": [float(kernel_fx), float(kernel_fz)],
                    })
                    print(f"mesh={nx}x{ny}x{nz} f={frequency:7.1f} x_pm={x_pm:+.3f} | "
                          f"reference F=({reference[0]:+.3e},{reference[1]:+.3e},"
                          f"{reference[2]:+.3e}) | kernel (Fx,Fz)=({kernel_fx:+.3e},"
                          f"{kernel_fz:+.3e})", flush=True)

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
    checks = {
        "reference_horizontal_force_vanishes_when_centred": all(
            max(abs(c["reference_force_N"][0]), abs(c["reference_force_N"][1]))
            < 1e-9 * abs(c["reference_force_N"][2]) for c in centred),
        "reference_lift_repels_the_magnet": all(
            c["reference_force_N"][2] < 0.0 for c in cases),
        "reference_mesh_converged_within_1pct": max(mesh_change) < 0.01,
        "reference_lift_rises_with_frequency": all(
            a < b for a, b in zip(lifts, lifts[1:])),
        "reference_lift_below_image_bound": max(lifts) < lift_bound,
    }
    kernel_record = {
        "centred_horizontal_force_N": [c["kernel_force_xz_N"][0] for c in centred],
        "centred_vertical_force_N": [c["kernel_force_xz_N"][1] for c in centred],
        "violates_mirror_symmetry": any(
            abs(c["kernel_force_xz_N"][0]) > 1.0 for c in centred),
        "produces_no_lift_when_centred": all(
            abs(c["kernel_force_xz_N"][1]) < 1e-9 for c in centred),
        "exceeds_image_bound": any(
            abs(c["kernel_force_xz_N"][0]) > lift_bound for c in centred),
        "cause": ("J_y = -omega sigma Im(v) treats the reaction-field z component "
                  "v [T] as a vector potential; J must be (1/mu) curl(v z)."),
    }
    import ngsolve
    payload = {
        "schema": "radia.maglev.ecb-foster-lorentz-reference.v1",
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
            "image_lift_bound_N": lift_bound,
        },
        "cases": cases,
        "reference_mesh_relative_change": mesh_change,
        "checks": checks,
        "reference_passed": all(checks.values()),
        "shipped_kernel": kernel_record,
    }
    output.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    print(f"[{'OK' if payload['reference_passed'] else 'FAIL'}] wrote {output}")
    return payload


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    payload = run(args.output)
    return 0 if payload["reference_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
