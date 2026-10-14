"""Table/direct finite-cell ESIM checks on three self-authored material laws.

Prescribed two-filament circuit; weak sphere and strong one-handle body.
This is discrete wiring/certification evidence, not material accuracy or
mesh convergence. Input VOLs must be checked before invoking this study.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import ngsolve as ng
from radia.panels import calc_inductance as ci
from radia.em_material import EMMaterial
from radia.esim_panel_evaluator import PanelESIMEvaluator


def run(directory):
    import radia.coil_from_cad as cad, radia.peec_bundle as bundle
    angle=np.arange(96)*2*np.pi/96; paths=[]
    for radius,height in [(.06,.005),(.065,-.005)]:
        p=np.c_[radius*np.cos(angle),radius*np.sin(angle),np.full(96,height)]
        paths.append(list(zip(p,np.roll(p,-1,axis=0))))
    resistance=np.diag([.002,.003]); inductance=np.array([[3e-7,1e-7],[1e-7,3.5e-7]])
    cad.filaments_from_step=lambda *a,**k:dict(filament_paths=paths,solver=None,seg_of_filament=None)
    bundle.build_loop_bundle_impedance=lambda *a,**k:(resistance,inductance)
    ci._peec_cross_section_model=lambda *a,**k:None
    ci._peec_internal_impedance_per_filament=lambda *a,**k:None
    def prescribed_coil(args):
        omega=2*np.pi*args.frequency; matrix=resistance+1j*omega*inductance
        unit=np.linalg.solve(matrix,np.ones(2)); z=1/unit.sum()
        return dict(source_type='filament',paths=paths,I_fil=args.current*unit/unit.sum(),
            L_coil=z.imag/omega,R_coil=z.real,t_coil_topology_s=0.,t_coil_solve_s=0.,n_filaments=2)
    ci._solve_coil_peec=prescribed_coil
    ci._export_coil_msh_viz=lambda *a,**k:{}
    h=np.r_[0.,np.geomspace(.001,10000,240)]
    laws=dict(near_linear=.2*h/(1+h/1000),sharp_knee=60*np.tanh(h/3),
              nonmonotone_mu=30*h**3/(8**3+h**3))
    rows=[]
    for name,magnetization in laws.items():
        bh=np.c_[h,4e-7*np.pi*(h+magnetization)]
        bhpath=directory/(name+'.txt');np.savetxt(bhpath,bh)
        for route,filename in [('weak-genus-0','analytic-sphere.vol'),('strong-genus-1','analytic-structured-torus.vol')]:
            pair=[]
            for mode in ['direct','table']:
                options=['--no-peec-proximity','--coil-solver','peec','--coil-step','prescribed-circuit.step',
                    '--vol',str(directory/filename),'--wp-label','sibc','--frequency','50000',
                    '--current','1','--sigma','5.8e7','--mu-r','1','--half-thickness','.003',
                    '--wp-bem-backend','intree-dense','--impedance-model','esim','--esim-per-panel',
                    '--bh-file',str(bhpath),'--esim-max-iter','60','--esim-relax','.5','--esim-tol','.001']
                # Omitting the option exercises the actual CLI default for table.
                if mode=='direct': options += ['--esim-panel-evaluator','direct']
                if route.startswith('strong'): options += ['--coupling-mode','strong','--coupling-max-iter','40',
                    '--coupling-tol','1e-8','--coupling-relax','.8']
                args=ci.build_argparser().parse_args(options)
                result=ci.run_inductance(args);assert result['status']=='ok',result
                z=np.array(result['esim_per_panel_Z_s_real'])+1j*np.array(result['esim_per_panel_Z_s_imag'])
                field=np.array(result['esim_per_panel_H_t'])
                cell=EMMaterial('self',5.8e7,1,bh_curve=bh.tolist()).create_esim_solver(50000,.003,geometry='cylinder')
                target=PanelESIMEvaluator(cell,mode='direct').certify(field)
                mismatch=float(np.max(abs(target-z)/abs(z)));assert mismatch<=.001
                evaluation=result['esim_panel_evaluation']
                assert evaluation['mode']==mode and evaluation['certification_cell_calls']==len(z)
                history=result['esim_history']
                pair.append((z,field,result['P_wp_W']))
                rows.append(dict(material=name,route=route,evaluator=mode,faces=len(z),
                    power_W=result['P_wp_W'],constitutive_error=mismatch,
                    outer_iterations=result['esim_iterations'],outer_history=history,
                    direct_final_cells=evaluation['certification_cell_calls'],
                    table_cells=evaluation['table_cell_calls'],
                    final_Z_real_range=[float(z.real.min()),float(z.real.max())],
                    final_H_range=[float(field.min()),float(field.max())]))
                print(f'Accepted {name} {route} {mode}',flush=True)
            comparisons={}
            for label,a,b in zip(['Z','H','power'],pair[0],pair[1]):
                error=float(np.max(abs(np.asarray(a)-b)/np.maximum(abs(np.asarray(a)),1e-30)))
                assert error<=.003,(name,route,label,error)
                comparisons[label]=error
            rows[-1]['table_direct_relative_difference']=comparisons
    return dict(schema='radia.esim-table-material-laws.v1',tolerance=.003,
        certification_tolerance=.001,
        qualification='Self-authored B(H); weak sphere and strong genus-1 prescribed two-filament circuit. Discrete material/EM certification only; no mesh-convergence, material accuracy or thermal claim.',
        laws=dict(near_linear='M=0.2H/(1+H/1000)',sharp_knee='M=60*tanh(H/3)',
                  nonmonotone_mu='M=30H^3/(8^3+H^3)'),rows=rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();ng.SetNumThreads(1)
    args.output.write_text(json.dumps(run(args.directory),indent=2)+'\n',encoding='utf-8')
