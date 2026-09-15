"""Diagnose the missing electric reciprocity term on a genus-zero cylinder.

The source is a closed 720-segment, one-ampere peak ring at r=30 mm.
The workpiece is r=25 mm, h=25 mm, sigma=5.8e7 S/m, mu_r=1 at 1 kHz.
Uses explicit in-memory Netgen surface meshes, not production VOL exports.

For exp(+i omega t), total H=-grad(phi), incident H=-grad(psi):
  delta_L = int(phi n.B_inc) / I**2
            + Zs/(i omega I**2) int(grad_s(psi).grad_s(phi)).
The second term is from -E cross H_inc in electromagnetic reciprocity.
It is bilinear, NOT a conjugated norm of the solved field. Therefore the
reaction power is independent of the surface-loss evaluation.

This diagnostic does NOT validate nonlinear ESIM, multiply-connected
surfaces, finite-source FEM equivalence, or the production application path.
The old term is retained in the output solely to reproduce the defect.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import ngsolve as ng
from netgen.occ import Cylinder, Pnt, Vec, OCCGeometry, Glue
import radia.bem_sibc_solver as sibc
import radia.biot_savart as bs
import radia._radia_pybind as native


def solve(maxh):
    body = Cylinder(Pnt(0, 0, -.0125), Vec(0, 0, 1), r=.025, h=.025)
    mesh = ng.Mesh(OCCGeometry(Glue(list(body.faces))).GenerateMesh(maxh=maxh))
    pts = np.array([list(v.point) for v in mesh.vertices])
    tri = np.array([[v.nr for v in el.vertices] for el in mesh.Elements(ng.BND)])
    cross = np.cross(pts[tri[:, 1]] - pts[tri[:, 0]], pts[tri[:, 2]] - pts[tri[:, 0]])
    areas = np.linalg.norm(cross, axis=1) / 2
    normals = cross / (2 * areas[:, None])
    cents = pts[tri].mean(axis=1)
    signed_volume = np.sum(areas * np.sum(cents * normals, axis=1)) / 3
    if signed_volume <= 0:
        raise ValueError("Expected outward-oriented surface")
    angle = np.linspace(0, 2 * np.pi, 721)
    loop = np.stack([.03 * np.cos(angle), .03 * np.sin(angle), np.zeros_like(angle)], axis=1)
    segments = np.stack([loop[:-1], loop[1:]], axis=1)
    h_inc = bs.h_segments_batch(segments, pts)
    h_centres = bs.h_segments_batch(segments, cents)
    psi, residual = sibc.compute_phi_inc_surface_poisson(pts, tri, h_inc, max_grad_residual=.2)
    omega = 2 * np.pi * 1000
    zs = (1 + 1j) * np.sqrt(omega * 4e-7 * np.pi / (2 * 5.8e7))
    solver = sibc.ScalarBIESIBCSolver(
        mesh, order=1, assemble_dense=True, use_intree_bem=True,
        intree_geom_order=1, intree_singular_n_q=6, intree_regular_quad_degree=7)
    result = solver.solve(psi.astype(complex), Z_s=zs, omega=omega)
    phi = result["phi_vec"]  # result['phi'] contains only the real part!
    magnetic = np.sum(phi[tri].mean(axis=1) * np.sum(normals * h_centres, axis=1) * areas) * 4e-7 * np.pi
    electric = zs / (1j * omega) * (psi @ solver.K @ phi)
    surface = result["P_density"] * result["area"]
    reaction_old = -.5 * omega * magnetic.imag
    reaction = -.5 * omega * (magnetic + electric).imag
    return dict(maxh=maxh, nv=mesh.nv, triangles=len(tri), P_surface=surface,
                P_reaction_magnetic_only=reaction_old, P_reaction_complete=reaction,
                P_electric_term=-.5 * omega * electric.imag,
                old_loss_to_reaction_ratio=surface / reaction_old,
                complete_power_relative_error=abs(surface-reaction)/surface,
                phi_projection_relative_residual=residual)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ng.SetNumThreads(4)
    rows = []
    with ng.TaskManager():
        for maxh in (.006, .003, .002):
            row = solve(maxh)
            rows.append(row)
            print(json.dumps(row, allow_nan=False), flush=True)
    files = [Path(__file__), Path(sibc.__file__), Path(bs.__file__), Path(native.__file__)]
    passed = bool(rows[-1]["complete_power_relative_error"] < .01)
    report = dict(scope=__doc__, host=platform.node(), python=sys.version,
                  ngsolve=ng.__version__, production_path_fixed=False,
                  source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                  rows=rows, finest_power_balance_tolerance=.01, passed=passed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")
    if not passed:
        raise RuntimeError("Complete reciprocity diagnostic failed; inspect JSON")


if __name__ == "__main__":
    main()
