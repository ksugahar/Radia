"""Independent CPU cell/direct-equation check; not independent BEM assembly."""
from pathlib import Path
import argparse,json,sys,types
import numpy as np
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[2]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('job',type=Path);args=parser.parse_args()
    package=types.ModuleType('radia');package.__path__=[str(ROOT/'src/radia')];sys.modules['radia']=package
    from radia.esim_cell_problem import ESIMFiniteSlabSolver
    cfg=loadmat(args.job/'input.mat',simplify_cells=True)['cfg']
    gpu=loadmat(args.job/'accepted.mat',simplify_cells=True)
    dense=lambda x:x.toarray() if hasattr(x,'toarray') else np.asarray(x)
    a=dense(cfg['eddy']);p=dense(cfg['loss_projection']);mapping=dense(cfg['operator_impedance_map'])
    temp=np.array(cfg['initial_temperature_K'],dtype=float);history=[temp.copy()];heats=[np.zeros(4)];zs=[np.zeros(4,dtype=complex)]
    knots=cfg['bh_knots'];axis=cfg['temperature_axis_K'];mu0=4*np.pi*1e-7
    maximum_residual=0.;maximum_field=0.
    for step in range(5):
        local=dense(cfg['temperature_evaluation'])@temp
        solvers=[]
        for e,t in enumerate(local):
            f=(t-axis[0])/(axis[1]-axis[0]);assert 0<=f<=1
            b=(1-f)*cfg['bh_values_T'][e,0,:]+f*cfg['bh_values_T'][e,1,:]
            rho=np.interp(t,axis,cfg['resistivity_ohm_m'][e])
            initial=b[1]/knots[1];sat=b[-1]/mu0-knots[-1]
            geom=cfg['esim_cell_parameters'][e]
            solver=ESIMFiniteSlabSolver(half_thickness=geom[3],sigma=1/rho,frequency=geom[4]/(2*np.pi),n_nodes=int(geom[0]),mu_r=initial/mu0,geometry='slab')
            # Fixed validation adapter selects the declared piecewise-linear
            # amplitude law; original CPU finite-difference solve is unchanged.
            def law(h,b=b,initial=initial,sat=sat):
                if h<1e-10:return initial
                return (np.interp(h,knots,b) if h<=knots[-1] else mu0*(h+sat))/h
            solver._get_mu=law
            solvers.append(solver)
        def cells(fields):
            values=[]
            for e,h in enumerate(fields):
                g=cfg['esim_cell_parameters'][e]
                r=solvers[e].solve(float(h),tol=g[5],max_iter=int(g[1]),relaxation=g[6])
                assert r['converged'];values.append(r['Z'])
            return np.array(values)
        fields=cfg['initial_field_A_per_m'].copy();rhs=cfg['rhs_per_amp']*1e4
        for iteration in range(300):
            z=cells(fields);matrix=a+(mapping@z).reshape(a.shape)
            x=np.linalg.solve(matrix,rhs)
            new=np.maximum(np.linalg.norm((p@x).reshape(4,3),axis=1),cfg['field_floor_A_per_m'])
            z=cells(new);matrix=a+(mapping@z).reshape(a.shape)
            residual=np.linalg.norm(matrix@x-rhs)/np.linalg.norm(rhs)
            change=np.linalg.norm(new-fields)/max(np.linalg.norm(new),1e-30)
            if iteration and residual<=cfg['constitutive_tolerance'] and change<=cfg['constitutive_tolerance']:break
            fields=.5*(fields+new)
        else:raise RuntimeError('Independent CPU constitutive iteration failed')
        maximum_residual=max(maximum_residual,residual);maximum_field=max(maximum_field,float(new.max()))
        heat=np.real(z)*np.sum(np.abs((p@x).reshape(4,3))**2,axis=1)*cfg['loss_weights'][::3]
        k=dense(cfg['stiffness']);c=dense(cfg['convection'])
        rhs_t=cfg['sample_time_s']*(dense(cfg['heat_projection'])@heat-k@temp+c@(300-temp))
        temp=temp+np.linalg.solve(dense(cfg['thermal_system']),rhs_t).real
        history.append(temp.copy());heats.append(heat);zs.append(z)
    relative=lambda x,y:float(np.linalg.norm(x-y)/max(np.linalg.norm(y),1e-30))
    errors=dict(temperature=relative(np.array(history),gpu['t']),heat=relative(np.array(heats),gpu['q']),impedance=relative(np.array(zs),gpu['z']),electromagnetic=relative(x,gpu['em']))
    result=dict(accepted=all(np.isfinite(v) and v<=1e-7 for v in errors.values()),scope='Independent CPU ESIM finite differences and dense complete equations; shared native geometry assembly',relative_errors=errors,true_relative_residual=maximum_residual,max_field_A_per_m=maximum_field)
    (args.job/'cpu-result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result));assert result['accepted']

if __name__=='__main__':main()
