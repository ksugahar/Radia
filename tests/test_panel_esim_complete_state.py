"""Native complete-state material iteration contracts; no external fixtures."""
import numpy as np
import pytest
import ngsolve as ng
from radia.esim_panel_evaluator import PanelESIMEvaluator, solve_panel_esim
from tests.test_bem_strong_panel_impedance import setup_cases, make_solver, compare


class Law:
    def __init__(self, z, nonlinear=False):
        self.z, self.nonlinear = z, nonlinear
    def solve(self, h):
        return dict(Z=self.z*(1+.15*np.log1p(h)) if self.nonlinear else self.z,
                    converged=True)


def evaluate(law):
    return PanelESIMEvaluator(law, configuration=lambda: {'analytic_law': 1})


@pytest.mark.parametrize('model', ['surface', 'filament'])
@pytest.mark.parametrize('genus', [0, 1])
def test_constant_law_and_nonlinear_seed_independence(setup_cases, model, genus):
    omega=2*np.pi*50000; z=(1+1j)*np.sqrt(omega*4e-7*np.pi/(2*5.8e7))
    with ng.TaskManager():
        solver=make_solver(setup_cases, model, genus)
        def solve_em(values):
            return solver.solve(values,omega,max_iter=40,tol=1e-8,relax=.8,loop_dof=bool(genus))
        fixed=solve_em(z)
        impedance,state,record=solve_panel_esim(evaluate(Law(z)),solve_em,
            np.full(len(solver.wp_tris),z),tolerance=1e-5,max_iter=2,
            inner_tolerance=1e-8)
        compare(fixed,state)
        assert record['esim_panel_evaluation']['certification_cell_calls']==len(solver.wp_tris)
        certified=[]
        for seed in (.5,2.):
            impedance,state,record=solve_panel_esim(evaluate(Law(z,True)),solve_em,
                np.full(len(solver.wp_tris),seed*z),tolerance=1e-5,max_iter=30,
                relaxation=.8,inner_tolerance=1e-8)
            fields=np.linalg.norm(state['wp_H_t_tri'],axis=1)
            target=z*(1+.15*np.log1p(np.maximum(fields,1e-3)))
            assert np.max(abs(target-impedance.values)/abs(impedance.values))<=1e-5
            np.testing.assert_array_equal(state['Z_s_per_panel'],impedance.values)
            assert record['esim_inner_history'][-1]['final_certification']
            assert np.ptp(impedance.values.real)>1e-4*abs(z)
            certified.append(impedance.values)
        np.testing.assert_allclose(*certified,rtol=1e-5,atol=1e-15)


def test_final_direct_certificate_rejects_stale_field():
    law=Law(1+1j,True); evaluator=evaluate(law); calls=[]
    def changing_em(z):
        calls.append(z.values.copy())
        return dict(wp_H_t_tri=np.array([[0.,0.,0.]]) if len(calls)==1 else np.array([[100.,0.,0.]]))
    seed=np.array([law.solve(1e-3)['Z']])
    with pytest.raises(RuntimeError,match='direct ESIM'):
        solve_panel_esim(evaluator,changing_em,seed,tolerance=1e-3,max_iter=2)
    assert len(calls)==2 and evaluator.counts['certification_cell_calls']==1


def test_material_budget_and_inner_tolerance_fail_loud():
    with pytest.raises(ValueError,match='Inner coupling tolerance'):
        solve_panel_esim(evaluate(Law(1+1j)),None,[1+1j],tolerance=1e-3,
                         max_iter=2,inner_tolerance=1e-3)
    with pytest.raises(RuntimeError,match='did not converge'):
        solve_panel_esim(evaluate(Law(2+2j)),lambda z:dict(wp_H_t_tri=np.ones((1,3))),
                         [1+1j],tolerance=1e-6,max_iter=1)


@pytest.mark.parametrize('genus', [0, 1])
def test_surface_coil_actual_current_law_and_secant_power(setup_cases, genus):
    omega=2*np.pi*50000; z=(1+1j)*np.sqrt(omega*4e-7*np.pi/(2*5.8e7))
    with ng.TaskManager():
        solver=make_solver(setup_cases,'surface',genus)
        def solve_em(values):
            return solver.solve(values,omega,max_iter=40,tol=1e-8,relax=.8,loop_dof=bool(genus))
        def material(scale, callback):
            return solve_panel_esim(evaluate(Law(z,True)),callback,
                np.full(len(solver.wp_tris),z),tolerance=1e-5,max_iter=30,
                relaxation=.8,inner_tolerance=1e-8,field_scale=scale)
        _,one,_=material(1.,solve_em)
        z_two,two,_=material(2.,solve_em)
        def independently_scaled(values):
            result=solve_em(values)
            for key in ('P_total','body_reaction_power_W','port_power_W','coil_loss_W',
                        'coil_loss_air_W','coil_loss_change_W','wp_q_tri'):
                result[key]=4*result[key]
            result['wp_H_t_tri']=2*result['wp_H_t_tri']
            return result
        z_actual,actual,_=material(1.,independently_scaled)
        np.testing.assert_allclose(z_two.values,z_actual.values,rtol=1e-10,atol=1e-18)
        np.testing.assert_allclose(4*two['P_total'],actual['P_total'],rtol=1e-10)
        assert abs(two['P_total']/one['P_total']-1)>1e-3


@pytest.mark.parametrize('nonlinear', [False,True])
def test_weak_genus_one_law_uses_complete_loop_field(setup_cases,tmp_path,monkeypatch,nonlinear):
    from radia.panels import calc_inductance as ci
    import em_material
    import radia.esim_panel_evaluator as evaluation
    mesh=setup_cases[0][1]
    mesh.ngmesh.SetBCName(0,'sibc')
    omega=2*np.pi*50000;z=(1+1j)*np.sqrt(omega*4e-7*np.pi/(2*5.8e7))
    bh=tmp_path/'self-law.txt';bh.write_text('0 0\n100 0.0001256637061435917\n',encoding='utf-8')
    path=str(tmp_path/'in-memory-body.vol')
    # Paired field sidecars identify this exact generated surface; the test
    # passes its in-memory mesh directly instead of exercising a VOL loader.
    mesh.ngmesh.Save(path)
    original_mesh=ng.Mesh
    monkeypatch.setattr(ng,'Mesh',lambda value,*a,**k:mesh if value==path else original_mesh(value,*a,**k))
    monkeypatch.setattr(em_material.EMMaterial,'create_esim_solver',lambda *a,**k:Law(z,nonlinear))
    original_evaluator=evaluation.PanelESIMEvaluator
    monkeypatch.setattr(evaluation,'PanelESIMEvaluator',lambda cell,**k:
        original_evaluator(cell,configuration=lambda:dict(analytic_law=1),**k))
    args=ci.build_argparser().parse_args(['--coil-solver','peec','--coil-step','prescribed.step',
        '--vol',path,'--wp-label','sibc','--frequency','50000','--current','1',
        '--sigma','5.8e7','--mu-r','1','--wp-bem-backend','intree-dense',
        '--impedance-model','esim','--esim-per-panel','--bh-file',str(bh),
        '--esim-tol','1e-5','--esim-max-iter','30','--esim-relax','.8'])
    source=dict(source_type='filament',paths=setup_cases[2],I_fil=np.array([.6,.4]))
    with ng.TaskManager():
        result=ci._solve_workpiece_weak_coupled(args,source)
        assert result['wp_loop_dof'] and result['wp_loop_faraday_residual_rel']<=1e-6
        field=np.array(result['esim_per_panel_H_t'])
        values=np.array(result['esim_per_panel_Z_s_real'])+1j*np.array(result['esim_per_panel_Z_s_imag'])
        target=np.array([Law(z,nonlinear).solve(max(h,1e-3))['Z'] for h in field])
        assert np.max(abs(target-values)/abs(values))<=1e-5
        assert result['esim_panel_evaluation']['certification_cell_calls']==len(field)
        assert result['wp_power_balance_relative_error']<.1
        if not nonlinear:
            args.impedance_model='sibc';args.esim_per_panel=False
            fixed=ci._solve_workpiece_weak_coupled(args,source)
            for key in ('P_wp','H_t_rms','delta_L_nH','delta_R_mOhm'):
                np.testing.assert_allclose(result[key],fixed[key],rtol=1e-10,atol=1e-18)


def test_reflected_port_gate_rejects_a_stale_air_loss_change():
    def inconsistent_state(z):
        return dict(wp_H_t_tri=np.ones((1,3)),Z_s_per_panel=z.values.copy(),
            converged=True,coupling_residual=0.,body_residual=0.,iterations=2,
            body_power_balance_relative_error=0.,body_reaction_power_W=.1,
            coil_loss_W=1.,coil_loss_air_W=1.,coil_loss_change_W=.5,
            port_power_W=1.1)
    with pytest.raises(RuntimeError,match='reflected port'):
        solve_panel_esim(evaluate(Law(1+1j)),inconsistent_state,[1+1j],
            tolerance=1e-3,max_iter=2,inner_tolerance=1e-8)
