"""Self-authored thin-ring circuit check for the production scalar loop closure.

Run with a built Radia/NGSolve environment, e.g. --levels 32 64 128.
These levels have 512, 2048, 8192 flat triangles. The finest level requires
dense operators and substantial memory/time. This is an analytic approximation
at finite wire thickness, not a finite-element comparison or an error bound.
The production closure is the only solver path used by this script.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import netgen.meshing as nm
import ngsolve as ng
import numpy as np

from radia.bem_loop_extension import solve_loop_extended
from radia.bem_sibc_solver import ScalarBIESIBCSolver


def ring_mesh(nphi):
    if nphi < 16 or nphi % 4:
        raise ValueError("Each level must be a multiple of four, at least 16")
    ntheta = nphi // 4
    t, p = np.meshgrid(np.arange(ntheta)*2*np.pi/ntheta,
                       np.arange(nphi)*2*np.pi/nphi, indexing="ij")
    radius, wire = .03, .003
    points = np.stack([(radius+wire*np.cos(t))*np.cos(p),
                       (radius+wire*np.cos(t))*np.sin(p),
                       wire*np.sin(t)], axis=-1).reshape(-1, 3)
    triangles = []
    for k in range(ntheta):
        for j in range(nphi):
            a, b = k*nphi+j, k*nphi+(j+1) % nphi
            c, d = ((k+1) % ntheta)*nphi+j, ((k+1) % ntheta)*nphi+(j+1) % nphi
            triangles.extend([(a, b, d), (a, d, c)])
    # This analytic parameterization has outward orientation.
    mesh = nm.Mesh(dim=3)
    face = mesh.Add(nm.FaceDescriptor(surfnr=1, domin=1, bc=1))
    vertices = [mesh.Add(nm.MeshPoint(nm.Pnt(*p))) for p in points]
    for tri in triangles:
        mesh.Add(nm.Element2D(face, [vertices[i] for i in tri]))
    return ng.Mesh(mesh), points, np.asarray(triangles, dtype=np.int64)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--levels", nargs="+", type=int, default=[32])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ng.SetNumThreads(1)
    mu0, radius, wire, frequency, sigma = 4e-7*math.pi, .03, .003, 5e4, 5.8e7
    omega = 2*math.pi*frequency
    zs = (1+1j)*math.sqrt(omega*mu0/(2*sigma))
    inductance = mu0*radius*(math.log(8*radius/wire)-2)
    alpha_analytic = -1j*omega*mu0*math.pi*radius**2/(zs*radius/wire+1j*omega*inductance)

    def a_inc(points):
        return .5*mu0*np.stack([-points[:, 1], points[:, 0],
                                np.zeros(len(points))], axis=1)

    rows = []
    for level in args.levels:
        with ng.TaskManager():
            mesh, points, triangles = ring_mesh(level)
            solver = ScalarBIESIBCSolver(mesh, order=1, assemble_dense=True,
                use_intree_bem=True, intree_geom_order=1,
                intree_singular_n_q=6, intree_regular_quad_degree=7)
            cases = []
            for offset in [0., .25*wire]:
                angles = np.arange(96)*2*np.pi/96
                carrier = np.c_[(radius+offset)*np.cos(angles),
                                (radius+offset)*np.sin(angles), np.zeros(96)]
                out = solve_loop_extended(solver, -points[:, 2].astype(complex),
                    zs, omega, a_inc, section_anchor=(radius, 0), carrier_ring=carrier)
                cases.append(dict(carrier_offset_m=offset,
                    alpha=[out["alpha"].real, out["alpha"].imag],
                    analytic_alpha_relative_error=abs(out["alpha"]/alpha_analytic-1),
                    heat_power_w=out["P_total"], reaction_power_w=out["P_reaction"],
                    power_balance_relative=abs(out["P_reaction"]/out["P_total"]-1),
                    linear_residual=out["linear_residual_rel"],
                    faraday_residual=out["faraday_residual_rel"],
                    seam_jump_deviation=out["theta_jump_max_deviation"],
                    work_diagnostics=out["loop_work_diagnostics"]))
            identity = hashlib.sha256(points.tobytes()+triangles.tobytes()).hexdigest()
            rows.append(dict(triangles=len(triangles), vertices=len(points),
                             geometry_sha256=identity, cases=cases))
        result = dict(schema="radia.scalar-loop-analytic-ring.v1",
            model=dict(radius_m=radius, wire_radius_m=wire, frequency_hz=frequency,
                       conductivity_s_per_m=sigma, axial_field_a_per_m=1),
            analytic_alpha=[alpha_analytic.real, alpha_analytic.imag], rows=rows,
            qualification="Finite-thickness thin-ring approximation; no rigorous bound. "
                          "Observed assembly time includes shared-host load; no speed claim.")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
