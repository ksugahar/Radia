"""Face-aligned prescribed impedance on a self-authored genus-1 ring.

Use --levels 32 64 128 for h-halving. Only the adopted production solver
is used. Jump-edge field differences are diagnostics, not interface traces.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ngsolve as ng
import numpy as np

from radia.bem_loop_extension import solve_loop_extended
from radia.bem_sibc_solver import ScalarBIESIBCSolver
from radia.surface_impedance import PanelSurfaceImpedance
from run_loop_work_ring import ring_mesh


def jump_edge_diagnostic(triangles, values, field):
    """Compare adjacent face averages across each impedance jump edge."""
    edges = {}
    differences = []
    for face, tri in enumerate(triangles):
        for a, b in zip(tri, np.roll(tri, -1)):
            key = tuple(sorted((int(a), int(b))))
            if key in edges:
                other = edges[key]
                if values[face] != values[other]:
                    differences.append(float(np.linalg.norm(field[face]-field[other])))
            else:
                edges[key] = face
    return dict(jump_edges=len(differences),
                maximum_adjacent_face_difference_a_per_m=max(differences, default=0.),
                qualification="Adjacent face averages; no carrier line trace or continuity assertion.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--levels", nargs="+", type=int, default=[32])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ng.SetNumThreads(1)
    mu0, omega = 4e-7*np.pi, 2*np.pi*50000
    z = (1+1j)*np.sqrt(omega*mu0/(2*5.8e7))

    def a_inc(p):
        return .5*mu0*np.c_[-p[:, 1], p[:, 0], np.zeros(len(p))]

    rows = []
    for level in args.levels:
        with ng.TaskManager():
            mesh, points, triangles = ring_mesh(level)
            bem = ScalarBIESIBCSolver(mesh, order=1, assemble_dense=True,
                use_intree_bem=True, intree_geom_order=1,
                intree_singular_n_q=6, intree_regular_quad_degree=7)
            cell = np.arange(len(triangles)) // (2*level)
            values = np.where((cell >= level/16) & (cell < 3*level/16), .25*z, 3*z)
            area = .5*np.linalg.norm(np.cross(points[triangles[:, 1]]-points[triangles[:, 0]],
                                             points[triangles[:, 2]]-points[triangles[:, 0]]), axis=1)
            cases = []
            for offset in [0., .00075]:
                angles = np.arange(96)*2*np.pi/96
                carrier = np.c_[(.03+offset)*np.cos(angles), (.03+offset)*np.sin(angles), np.zeros(96)]
                controls = dict(section_anchor=(.03, 0), carrier_ring=carrier)
                nonuniform = solve_loop_extended(bem, -points[:, 2].astype(complex),
                    PanelSurfaceImpedance(values), omega, a_inc, **controls)
                uniform = solve_loop_extended(bem, -points[:, 2].astype(complex),
                    np.average(values, weights=area), omega, a_inc, **controls)
                np.testing.assert_array_equal(nonuniform["Z_s_per_panel"], values)
                cases.append(dict(carrier_offset_m=offset,
                    alpha=[nonuniform["alpha"].real, nonuniform["alpha"].imag],
                    heat_power_w=nonuniform["P_total"], uniform_mean_heat_power_w=uniform["P_total"],
                    delta_power_w=nonuniform["P_total"]-uniform["P_total"],
                    reaction_power_w=nonuniform["P_reaction"],
                    power_balance_relative=abs(nonuniform["P_reaction"]/nonuniform["P_total"]-1),
                    linear_residual=nonuniform["linear_residual_rel"],
                    faraday_residual=nonuniform["faraday_residual_rel"],
                    seam_jump_deviation=nonuniform["theta_jump_max_deviation"],
                    jump_diagnostic=jump_edge_diagnostic(triangles, values, nonuniform["H_t_tri"])))
            rows.append(dict(triangles=len(triangles), vertices=len(points), cases=cases))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(dict(schema="radia.panel-loop-ring.v1", rows=rows,
            law="Inner minor-angle half .25 Zref; outer half 3 Zref, aligned to faces.",
            qualification="Self-authored prescribed linear weak-coupling case; no nonlinear or finite-element claim."), indent=2)+"\n")


if __name__ == "__main__":
    main()
