"""Self-authored fixed-panel strong torus check with a prescribed circuit.

Coil fitting and thermal coupling are outside this electromagnetic example.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import ngsolve as ng
from run_loop_work_ring import ring_mesh
from radia.peec_coupled_bem_solver import CoupledPEECBEMSolver
from radia.surface_impedance import PanelSurfaceImpedance


def run(levels=(16, 32)):
    angle = np.arange(96)*2*np.pi/96
    paths = []
    for radius, height in [(.06, .005), (.065, -.005)]:
        points = np.c_[radius*np.cos(angle), radius*np.sin(angle), np.full(96, height)]
        paths.append(list(zip(points, np.roll(points, -1, axis=0))))
    omega = 2*np.pi*50000
    z = (1+1j)*np.sqrt(omega*4e-7*np.pi/(2*5.8e7))
    rows = []
    with ng.TaskManager():
        for resolution in levels:
            mesh, points, tri = ring_mesh(resolution)
            solver = CoupledPEECBEMSolver(paths, np.diag([.002, .003]),
                np.array([[3e-7, 1e-7], [1e-7, 3.5e-7]]), mesh)
            values = np.where(points[tri].mean(axis=1)[:,2]>=0, 3*z, .25*z)
            results = {}
            for name, impedance, current in [
                    ('scalar', z, 1.),
                    ('uniform', PanelSurfaceImpedance(np.full(len(tri), z)), 1.),
                    ('panel', PanelSurfaceImpedance(values), 1.),
                    ('rewrite', PanelSurfaceImpedance(1.1*values), 1.),
                    ('scaled', PanelSurfaceImpedance(values), 5.)]:
                result = solver.solve(impedance, omega, I_port=current,
                    max_iter=40, tol=1e-9, relax=.8, loop_dof=True)
                assert result['body_power_balance_relative_error'] < .1
                assert result['coupling_residual'] <= 1e-9 and result['body_residual'] <= 1e-6
                assert np.isclose(result['port_power_W'], result['body_reaction_power_W']
                                  +result['coil_loss_W'], rtol=1e-12)
                assert np.isclose(.5*result['Delta_R']*current**2,
                    result['body_reaction_power_W']+result['coil_loss_change_W'], rtol=1e-12)
                if isinstance(impedance, PanelSurfaceImpedance):
                    np.testing.assert_array_equal(result['Z_s_per_panel'], impedance.values)
                    np.testing.assert_allclose(result['wp_q_tri'], .5*impedance.values.real
                        *np.sum(abs(result['wp_H_t_tri'])**2, axis=1), rtol=1e-12)
                results[name] = result
                row = {k:v for k,v in result.items() if isinstance(v,(float,int,bool,str))}
                row.update(case=name, triangles=len(tri), current_a=current,
                    alpha=[result['wp_loop_alpha'].real,result['wp_loop_alpha'].imag])
                rows.append(row)
            for key in ('P_total','L_total','R_total','body_reaction_power_W'):
                assert np.isclose(results['scalar'][key],results['uniform'][key],rtol=1e-10)
            assert abs(results['panel']['P_total']/results['scalar']['P_total']-1)>1e-4
            assert abs(results['rewrite']['P_total']/results['panel']['P_total']-1)>1e-4
            assert np.isclose(results['scaled']['P_total'],25*results['panel']['P_total'],rtol=1e-10)
    return dict(schema='radia.fixed-panel-strong.v1',
        qualification='Self-authored prescribed two-filament circuit and flat torus. Finite-mesh EM checks; no coil-fitting, thermal or general physical-accuracy certification.', rows=rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    ng.SetNumThreads(1)
    args.output.write_text(json.dumps(run(),indent=2)+'\n',encoding='utf-8')
