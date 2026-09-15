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
from radia.workpiece_surface import _complete_sibc_reaction, delta_L_telegen_phiB
from radia.bem_loop_extension import solve_loop_extended, A_from_filaments


def solve(maxh, bore=0.0, coil_a=0.0):
    body = Cylinder(Pnt(0, 0, -.0125), Vec(0, 0, 1), r=.025, h=.025)
    if bore:
        body = body - Cylinder(Pnt(0, 0, -.013), Vec(0, 0, 1), r=bore, h=.026)
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
    # Degree-five disk cubature for the same uniform-J circular FEM coil.
    offsets = [(0., 0., .25)] + [
        (coil_a*np.sqrt(2/3)*np.cos(t), coil_a*np.sqrt(2/3)*np.sin(t), .125)
        for t in np.arange(6)*np.pi/3]
    paths, currents = [], []
    for dr, dz, weight in offsets:
        loop = np.stack([(.03+dr)*np.cos(angle), (.03+dr)*np.sin(angle),
                         np.full_like(angle, dz)], axis=1)
        paths.append(np.stack([loop[:-1], loop[1:]], axis=1))
        currents.append(weight)
    h_inc = sum(w*bs.h_segments_batch(p, pts) for p,w in zip(paths,currents))
    psi, residual = sibc.compute_phi_inc_surface_poisson(pts, tri, h_inc, max_grad_residual=.2)
    omega = 2 * np.pi * 1000
    zs = (1 + 1j) * np.sqrt(omega * 4e-7 * np.pi / (2 * 5.8e7))
    solver = sibc.ScalarBIESIBCSolver(
        mesh, order=1, assemble_dense=True, use_intree_bem=True,
        intree_geom_order=1, intree_singular_n_q=6, intree_regular_quad_degree=7)
    result = solver.solve(psi.astype(complex), Z_s=zs, omega=omega)
    phi = result["phi_vec"]  # result['phi'] contains only the real part!
    magnetic = delta_L_telegen_phiB(paths, currents, cents, areas, normals,
                                   phi[tri].mean(axis=1), I_port=1.)
    complete = _complete_sibc_reaction(magnetic, psi, phi, solver.K, zs, omega, 1.)
    surface = result["P_density"] * result["area"]
    if bore:
        extended = solve_loop_extended(solver, psi, zs, omega,
            lambda p: A_from_filaments(p, paths, currents))
        surface = extended["P_total"]
        complete = extended["reaction_integral"]
        h_total = extended["H_t_tri"]
        heat, projected = sibc._project_sibc_surface_heat(
            solver.fes, extended["phi_u"], zs, element_heat=extended["q_tri"])
    else:
        heat, projected = sibc._project_sibc_surface_heat(solver.fes, phi, zs)
        real, imag = ng.GridFunction(solver.fes), ng.GridFunction(solver.fes)
        real.vec.FV().NumPy()[:] = phi.real
        imag.vec.FV().NumPy()[:] = phi.imag
        h_total = np.stack([-(np.array(ng.Integrate(ng.grad(real)[j], mesh, ng.BND, element_wise=True))
            + 1j*np.array(ng.Integrate(ng.grad(imag)[j], mesh, ng.BND, element_wise=True)))/areas
            for j in range(3)], axis=1)
    z_samples = np.linspace(-.010, .010, 11)
    hz = []
    for zi in z_samples:
        selected = (abs(cents[:,2]-zi) <= .001) & (np.linalg.norm(cents[:,:2],axis=1)>.024) & (abs(normals[:,2])<.1)
        if not np.any(selected):
            raise ValueError("Surface too coarse for local-field bins")
        hz.append(np.sum(areas[selected]*h_total[selected,2])/sum(areas[selected]))
    hz = np.array(hz)
    reaction_old = -.5 * omega * magnetic.imag
    reaction = -.5 * omega * complete.imag
    return dict(maxh=maxh, bore=bore, coil_a=coil_a, nv=mesh.nv, triangles=len(tri), P_surface=surface,
                P_reaction_magnetic_only=reaction_old, P_reaction_complete=reaction,
                P_added_reciprocity_terms=reaction-reaction_old,
                heat_projection_relative_error=abs(projected-surface)/surface,
                z=z_samples.tolist(), H_z_real=hz.real.tolist(), H_z_imag=hz.imag.tolist(),
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
                  ngsolve=ng.__version__, workpiece_stage_exercised=False,
                  source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                  rows=rows, finest_power_balance_tolerance=.01, passed=passed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")
    if not passed:
        raise RuntimeError("Complete reciprocity diagnostic failed; inspect JSON")


if __name__ == "__main__":
    main()
