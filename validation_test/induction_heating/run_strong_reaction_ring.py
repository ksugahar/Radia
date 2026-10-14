"""Self-authored two-filament/torus scalar strong-coupling checks.

Prescribed circuit matrices isolate body reaction from CAD coil fitting.
No thermal solve or external solver reference is used.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import ngsolve as ng
from run_loop_work_ring import ring_mesh
from radia.bem_complete_reaction import electric_incident_vertex_load, solve_complete_body
from radia.bem_loop_extension import A_from_filaments
from radia.bem_sibc_solver import SurfacePoissonPhiInc
from radia.biot_savart import h_segments_batch
from radia.peec_coupled_bem_solver import CoupledPEECBEMSolver


def run(levels=(16, 32)):
    angle = np.arange(96)*2*np.pi/96
    paths = []
    for radius, height in [(.06, .005), (.065, -.005)]:
        p = np.c_[radius*np.cos(angle), radius*np.sin(angle), np.full(96, height)]
        paths.append(list(zip(p, np.roll(p, -1, axis=0))))
    omega = 2*np.pi*50000
    z = (1+1j)*np.sqrt(omega*4e-7*np.pi/(2*5.8e7))
    rows = []
    with ng.TaskManager():
        for resolution in levels:
            mesh, points, tri = ring_mesh(resolution)
            solver = CoupledPEECBEMSolver(paths, np.diag([.002, .003]),
                np.array([[3e-7, 1e-7], [1e-7, 3.5e-7]]), mesh)
            result = solver.solve(z, omega, max_iter=40, tol=1e-9,
                                  relax=.8, loop_dof=True)
            assert result['body_power_balance_relative_error'] < .1
            assert result['coupling_residual'] <= 1e-9
            assert result['body_residual'] <= 1e-6
            assert np.isclose(result['port_power_W'], result['body_reaction_power_W']
                              +result['coil_loss_W'], rtol=1e-12)
            assert np.isclose(.5*result['Delta_R'], result['body_reaction_power_W']
                              +result['coil_loss_change_W'], rtol=1e-12)
            poisson = SurfacePoissonPhiInc(points, tri)
            centers = points[tri].mean(axis=1)
            hs = np.stack([h_segments_batch(path, points) for path in paths])
            aa = np.stack([A_from_filaments(centers, [path], [1.]) for path in paths])
            reaction = np.empty((2, 2), complex)
            for k in range(2):
                phi, _ = poisson(hs[k], max_grad_residual=.1)
                body = solve_complete_body(solver.wp_solver, poisson, phi, z, omega,
                    lambda x, k=k: A_from_filaments(x, [paths[k]], [1.]), loop_dof=True)
                load = electric_incident_vertex_load(poisson, z, body['field'])
                reaction[:, k] = (1j*omega*np.einsum('ktd,td,t->k', aa,
                    body['current'], poisson._areas)+np.einsum('kvd,vd->k', hs, load))
            symmetry = float(np.linalg.norm(reaction-reaction.T)/np.linalg.norm(reaction))
            row = {k:v for k,v in result.items() if isinstance(v, (float, int, bool, str))}
            row.update(triangles=len(tri), reaction_reciprocity_relative_error=symmetry,
                       currents=[[v.real, v.imag] for v in result['I_f']],
                       alpha=[result['wp_loop_alpha'].real, result['wp_loop_alpha'].imag])
            if resolution == levels[0]:
                scaled = solver.solve(z, omega, I_port=5., max_iter=40, tol=1e-9,
                                      relax=.8, loop_dof=True)
                for key in ('P_total', 'body_reaction_power_W', 'port_power_W',
                            'coil_loss_W', 'coil_loss_change_W'):
                    assert np.isclose(scaled[key], 25*result[key], rtol=1e-9, atol=1e-20)
                np.testing.assert_allclose(scaled['I_f'], 5*result['I_f'], rtol=1e-9)
                assert np.isclose(scaled['wp_loop_alpha'], 5*result['wp_loop_alpha'], rtol=1e-9)
                row['terminal_current_scaling_passed'] = True
            rows.append(row)
            print(json.dumps(row), flush=True)
    return dict(schema='radia.scalar-strong-reaction.v1',
                qualification='Self-authored prescribed circuit and flat torus; finite-mesh power and reciprocity defects, no general accuracy claim.', rows=rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    ng.SetNumThreads(1)
    args.output.write_text(json.dumps(run(), indent=2)+'\n', encoding='utf-8')
