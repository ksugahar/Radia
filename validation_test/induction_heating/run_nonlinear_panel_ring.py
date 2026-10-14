"""Finite-cell nonlinear panel EM on a self-authored torus/circuit.

This checks discrete material/EM consistency, not general material accuracy,
coil CAD fitting, mesh convergence or thermal coupling.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import ngsolve as ng
from run_loop_work_ring import ring_mesh
from radia.em_material import EMMaterial
from radia.esim_cell_problem import require_esim_converged
from radia.esim_panel_evaluator import PanelESIMEvaluator, solve_panel_esim
from radia.peec_coupled_bem_solver import CoupledPEECBEMSolver


def run():
    frequency=50000.; omega=2*np.pi*frequency
    h=np.r_[0.,np.geomspace(.001,10000,160)]
    bh=np.c_[h,4e-7*np.pi*(h+100*np.tanh(h/5))].tolist()
    angle=np.arange(96)*2*np.pi/96; paths=[]
    for radius,height in [(.06,.005),(.065,-.005)]:
        p=np.c_[radius*np.cos(angle),radius*np.sin(angle),np.full(96,height)]
        paths.append(list(zip(p,np.roll(p,-1,axis=0))))
    rows=[]; solutions=[]
    with ng.TaskManager():
        mesh,points,tri=ring_mesh(16)
        solver=CoupledPEECBEMSolver(paths,np.diag([.002,.003]),
            np.array([[3e-7,1e-7],[1e-7,3.5e-7]]),mesh)
        for current,mode,seed_factor in [(1.,'direct',1.),(2.,'direct',1.),
                                        (1.,'table',1.),(1.,'direct',2.)]:
            cell=EMMaterial('self',5.8e7,1,bh_curve=bh).create_esim_solver(
                frequency,.003,geometry='cylinder')
            seed=complex(require_esim_converged(cell.solve(5.))['Z'])
            ev=PanelESIMEvaluator(cell,mode=mode,interpolation_tol=1e-4)
            def solve_em(z):
                return solver.solve(z,omega,I_port=current,max_iter=40,
                    tol=1e-8,relax=.8,loop_dof=True)
            z,state,record=solve_panel_esim(ev,solve_em,
                np.full(len(tri),seed_factor*seed),tolerance=1e-3,max_iter=30,
                relaxation=.7,inner_tolerance=1e-8)
            np.testing.assert_array_equal(z.values,state['Z_s_per_panel'])
            assert record['esim_panel_evaluation']['certification_cell_calls']==len(tri)
            np.testing.assert_allclose(state['wp_q_tri'],.5*z.values.real*
                np.sum(abs(state['wp_H_t_tri'])**2,axis=1),rtol=1e-12)
            solutions.append((z.values,state))
            rows.append(dict(current_A=current,evaluator=mode,seed_factor=seed_factor,
                triangles=len(tri),power_W=state['P_total'],
                Delta_L_H=state['Delta_L'],Delta_R_ohm=state['Delta_R'],
                body_residual=state['body_residual'],coil_residual=state['coupling_residual'],
                faraday_residual=state['wp_loop_faraday_residual'],
                body_balance=state['body_power_balance_relative_error'],
                constitutive_error=record['esim_fixed_point_relative_error'],
                outer_iterations=record['esim_iterations'],
                inner_iterations_total=record['esim_inner_iterations_total'],
                final_inner_iterations=record['esim_inner_history'][-1]['iterations'],
                initial_Z=[float((seed_factor*seed).real),float((seed_factor*seed).imag)],
                final_Z_real_range=[float(z.values.real.min()),float(z.values.real.max())],
                final_H_range=[float(np.linalg.norm(state['wp_H_t_tri'],axis=1).min()),
                               float(np.linalg.norm(state['wp_H_t_tri'],axis=1).max())],
                evaluation=dict(**record['esim_panel_evaluation'],seed_cell_calls=1),
                port_quantity_kind=record['port_quantity_kind']))
            print(f'Accepted current={current} mode={mode} seed={seed_factor}',flush=True)
    np.testing.assert_allclose(solutions[0][0],solutions[2][0],rtol=.003,atol=1e-15)
    np.testing.assert_allclose(solutions[0][0],solutions[3][0],rtol=.001,atol=1e-15)
    assert abs(solutions[1][1]['P_total']/(4*solutions[0][1]['P_total'])-1)>1e-3
    return dict(schema='radia.nonlinear-panel-complete-state.v1',
        qualification='Self-authored saturating B(H), prescribed two-filament circuit, flat 128-face torus. Discrete EM/constitutive consistency only; no material accuracy, CAD fitting, mesh convergence or thermal certification.',
        seed_relative_difference=float(np.max(abs(solutions[0][0]-solutions[3][0])/abs(solutions[0][0]))),
        table_relative_difference=float(np.max(abs(solutions[0][0]-solutions[2][0])/abs(solutions[0][0]))),
        two_amp_power_over_four_times_one_amp=float(solutions[1][1]['P_total']/(4*solutions[0][1]['P_total'])),
        material=dict(law='B=mu0*(H+100*tanh(H/5))',sigma_S_per_m=5.8e7,
            frequency_Hz=frequency,cell_geometry='cylinder',cell_radius_m=.003,
            cell_nodes=200,table_H_min_A_per_m=.001,table_H_max_A_per_m=10000,
            table_positive_rows=160),rows=rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();ng.SetNumThreads(1)
    args.output.write_text(json.dumps(run(),indent=2)+'\n',encoding='utf-8')
