import numpy as np
import ngsolve as ng
from ngsolve.meshes import MakeStructured3DMesh
from types import SimpleNamespace

from radia.hcurl_topology_optimization import (
    assemble_ngsolve_hcurl_resistance_shape_tangents,
    HCurlJouleLoadCase,
    HCurlConductivityInterpolation,
    HCurlMultiFrequencyJouleLinearization,
    linearize_and_solve_hcurl_sheet_joule_lp,
    linearize_and_solve_hcurl_activation_sheet_joule_lp,
    linearize_hcurl_joule_loss_from_ngsolve,
    linearize_hcurl_multifrequency_joule_loss_from_ngsolve,
    linearize_hcurl_multifrequency_activation_joule_loss_from_ngsolve,
    optimize_hcurl_eddy_bubble_hex_sheet,
    solve_hcurl_joule_lp,
)
from radia import vim


def test_ngsolve_hcurl_resistance_piola_scaling_derivative():
    mesh=MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1)
    fes=ng.HCurl(mesh,order=1)
    space=ng.VectorH1(mesh,order=1)
    scaling=ng.GridFunction(space);scaling.Set(ng.CF((ng.x,ng.y,ng.z)))
    vectors=np.random.default_rng(72).normal(size=(fes.ndof,2))
    with ng.TaskManager():
        result=assemble_ngsolve_hcurl_resistance_shape_tangents(
            fes,vectors,[scaling],conductivity=5.8e7)
    np.testing.assert_allclose(result.jacobian[0],-result.matrix,
        rtol=2e-12,atol=2e-20)


def test_hcurl_joule_modal_lp_obeys_bounds_move_and_volume():
    result=solve_hcurl_joule_lp([.5,.5],[-2.,1.],move_limit=.2,
        lower_bounds=0.,upper_bounds=1.,volume_gradient=[1.,1.],volume_limit=0.)
    assert np.all(np.abs(result.delta)<=.2+1e-12)
    assert np.sum(result.delta)<=1e-12
    assert result.design[0]>=result.design[1]


def test_gettrafo_hcurl_joule_complex_adjoint_closes_scaling_tangent():
    mesh=MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1)
    fes=ng.HCurl(mesh,order=2,nograds=True)
    field=ng.GridFunction(fes);field.Set(ng.CF((-ng.y,ng.x,0)))
    vectors=field.vec.FV().NumPy().copy()
    basis=vim.NgsolveHCurlCurlBasis(mesh,fes,vectors,intorder=4)
    interaction=vim.NgsolveHCurlCellVolumeInteraction(mesh,fes,vectors,basis,
        degree=0,projection_quad=4,outer_quad=4,projection_tolerance=1e-10)
    space=ng.VectorH1(mesh,order=1)
    scaling=ng.GridFunction(space);scaling.Set(ng.CF((ng.x,ng.y,ng.z)))
    with ng.TaskManager():
        result=linearize_hcurl_joule_loss_from_ngsolve(mesh=mesh,fes=fes,
            vectors=vectors[:,None],interaction=interaction,
            deformation_modes=[scaling],frequency_hz=100.,rhs=[1.],
            conductivity=5.8e7)
    R=result.resistance;L=interaction.matrix.to_dense();omega=2*np.pi*100.;step=2e-6
    def objective(scale):
        Rs=R/scale;Ls=scale*L;x=np.linalg.solve(Rs+1j*omega*Ls,np.ones(1))
        return .5*np.real(np.vdot(x,Rs@x))
    validation=(objective(1+step)-objective(1-step))/(2*step)
    np.testing.assert_allclose(result.gradient,[validation],rtol=2e-9,atol=2e-11)


def test_multifrequency_joule_objective_and_gradient_are_weighted_adjoint_sum():
    mesh=MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1)
    fes=ng.HCurl(mesh,order=2,nograds=True)
    field=ng.GridFunction(fes);field.Set(ng.CF((-ng.y,ng.x,0)))
    vectors=field.vec.FV().NumPy().copy()
    basis=vim.NgsolveHCurlCurlBasis(mesh,fes,vectors,intorder=4)
    interaction=vim.NgsolveHCurlCellVolumeInteraction(mesh,fes,vectors,basis,
        degree=0,projection_quad=4,outer_quad=4,projection_tolerance=1e-10)
    space=ng.VectorH1(mesh,order=1)
    scaling=ng.GridFunction(space);scaling.Set(ng.CF((ng.x,ng.y,ng.z)))
    cases=(HCurlJouleLoadCase(50.,np.array([1.+.2j]),.25),
           HCurlJouleLoadCase(500.,np.array([.7-.1j]),1.5))
    with ng.TaskManager():
        result=linearize_hcurl_multifrequency_joule_loss_from_ngsolve(
            mesh=mesh,fes=fes,vectors=vectors[:,None],interaction=interaction,
            deformation_modes=[scaling],load_cases=cases,conductivity=5.8e7)
        singles=[linearize_hcurl_joule_loss_from_ngsolve(
            mesh=mesh,fes=fes,vectors=vectors[:,None],interaction=interaction,
            deformation_modes=[scaling],frequency_hz=case.frequency_hz,
            rhs=case.rhs,conductivity=5.8e7) for case in cases]
    expected_objective=sum(case.weight*x.objective for case,x in zip(cases,singles))
    expected_gradient=sum(case.weight*x.gradient for case,x in zip(cases,singles))
    np.testing.assert_allclose(result.objective,expected_objective,rtol=2e-13)
    np.testing.assert_allclose(result.gradient,expected_gradient,rtol=2e-12,atol=1e-18)


def test_sheet_lp_adapter_maps_analytic_modes_without_inventing_topology_gradient(monkeypatch):
    import radia.hcurl_topology_optimization as ht
    linearization=HCurlMultiFrequencyJouleLinearization((),3.,np.array([-2.]),np.array([1.]))
    monkeypatch.setattr(ht,"linearize_hcurl_multifrequency_joule_loss_from_ngsolve",
        lambda **kwargs:linearization)
    state=SimpleNamespace(normal_displacement=np.zeros(2),thickness=np.ones(2),
        activation=np.ones(2))
    mapping=np.array([[1.,-1.,0.,0.,0.,0.]])
    step=linearize_and_solve_hcurl_sheet_joule_lp(state=state,mesh=None,fes=None,
        vectors=None,interaction=None,deformation_modes=[object()],load_cases=[object()],
        design_mode_jacobian=mapping,area=np.ones(2),volume_max=2.,
        displacement_move=.1,thickness_move=0.,activation_move=0.,
        thickness_bounds=(1.,1.))
    np.testing.assert_allclose(step.update.normal_displacement,[.1,-.1])
    np.testing.assert_allclose(step.update.thickness,state.thickness)
    np.testing.assert_allclose(step.update.activation,state.activation)


def test_hcurl_hex_driver_rebuilds_step_inputs_after_mesh_update(monkeypatch):
    import radia.hcurl_topology_optimization as ht
    import radia.sheet_metal_optimization as sm
    states=[SimpleNamespace(mesh=SimpleNamespace(name="initial")),
            SimpleNamespace(mesh=SimpleNamespace(name="remeshed"))]
    seen=[]
    monkeypatch.setattr(ht,"linearize_and_solve_hcurl_sheet_joule_lp",
        lambda **kwargs:seen.append((kwargs["state"],kwargs["mesh"],kwargs["token"])) or "step")
    def fake_driver(initial_state,*,linearize_step,**kwargs):
        assert initial_state is states[0] and kwargs["sentinel"]==7
        return [linearize_step(state) for state in states]
    monkeypatch.setattr(sm,"optimize_hex_sheet_topology",fake_driver)
    result=optimize_hcurl_eddy_bubble_hex_sheet(states[0],
        build_step_inputs=lambda state:{"token":state.mesh.name},sentinel=7)
    assert result==["step","step"]
    assert [(mesh.name,token) for _,mesh,token in seen]==[
        ("initial","initial"),("remeshed","remeshed")]


def test_eddy_bubble_activation_adjoint_matches_objective_regression_fd():
    mesh=MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1)
    fes=ng.HCurl(mesh,order=2,nograds=True)
    field=ng.GridFunction(fes);field.Set(ng.CF((-ng.y,ng.x,0)))
    vectors=field.vec.FV().NumPy().copy()
    basis=vim.NgsolveHCurlCurlBasis(mesh,fes,vectors,intorder=4)
    interaction=vim.NgsolveHCurlCellVolumeInteraction(mesh,fes,vectors,basis,
        degree=0,projection_quad=4,outer_quad=4,projection_tolerance=1e-10)
    rho=np.full(mesh.ne,.67)
    law=HCurlConductivityInterpolation(5.8e7,2.0e5,2.3)
    cases=(HCurlJouleLoadCase(240.,np.array([1.+.15j]),1.),)
    with ng.TaskManager():
        result=linearize_hcurl_multifrequency_activation_joule_loss_from_ngsolve(
            mesh=mesh,fes=fes,vectors=vectors[:,None],interaction=interaction,
            activation=rho,load_cases=cases,conductivity=law,
            inductance_power=1.4)
    sigma,_=law.evaluate(rho)
    # Structured tetrahedra share the same activation here, so this regression
    # can perturb all cells together without reassembling FE plumbing.
    Rbase=result.cases[0].resistance*sigma[0]
    Lbase=interaction.matrix.to_dense(); omega=2*np.pi*cases[0].frequency_hz
    def objective(value):
        sig,_=law.evaluate(np.array([value]))
        R=Rbase/sig[0];L=value**(2*1.4)*Lbase
        x=np.linalg.solve(R+1j*omega*L,cases[0].rhs)
        return .5*np.real(np.vdot(x,R@x))
    step=2e-6
    fd=(objective(rho[0]+step)-objective(rho[0]-step))/(2*step)
    np.testing.assert_allclose(np.sum(result.gradient),fd,rtol=3e-8,atol=2e-12)


def test_eddy_bubble_cellwise_activation_gradient_keeps_element_locality():
    mesh=MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1)
    fes=ng.HCurl(mesh,order=2,nograds=True)
    field=ng.GridFunction(fes);field.Set(ng.CF((-ng.y,ng.x,0)))
    vectors=field.vec.FV().NumPy().copy()
    basis=vim.NgsolveHCurlCurlBasis(mesh,fes,vectors,intorder=4)
    interaction=vim.NgsolveHCurlCellVolumeInteraction(mesh,fes,vectors,basis,
        degree=0,projection_quad=4,outer_quad=4,projection_tolerance=1e-10)
    law=HCurlConductivityInterpolation(5.8e7,3.e5,2.)
    case=(HCurlJouleLoadCase(120.,np.array([1.-.1j])),)
    rho=np.linspace(.45,.8,mesh.ne)
    def evaluate(values):
        return linearize_hcurl_multifrequency_activation_joule_loss_from_ngsolve(
            mesh=mesh,fes=fes,vectors=vectors[:,None],interaction=interaction,
            activation=values,load_cases=case,conductivity=law,
            inductance_power=1.2)
    with ng.TaskManager():
        result=evaluate(rho)
        step=2e-6; fd=[]
        for cell in (0,mesh.ne-1):
            plus=rho.copy();minus=rho.copy();plus[cell]+=step;minus[cell]-=step
            fd.append((evaluate(plus).objective-evaluate(minus).objective)/(2*step))
    np.testing.assert_allclose(result.gradient[[0,-1]],fd,rtol=5e-8,atol=2e-12)


def test_activation_sheet_lp_updates_only_material_topology(monkeypatch):
    import radia.hcurl_topology_optimization as ht
    linearization=HCurlMultiFrequencyJouleLinearization(
        (),1.,np.array([-2.,1.]),np.array([1.]))
    monkeypatch.setattr(ht,
        "linearize_hcurl_multifrequency_activation_joule_loss_from_ngsolve",
        lambda **kwargs:linearization)
    state=SimpleNamespace(normal_displacement=np.zeros(2),
        thickness=np.ones(2),activation=np.array([.5,.5]))
    step=linearize_and_solve_hcurl_activation_sheet_joule_lp(
        state=state,mesh=None,fes=None,vectors=None,interaction=None,
        load_cases=None,conductivity=None,area=np.ones(2),volume_max=1.,
        activation_move=.2)
    np.testing.assert_allclose(step.update.normal_displacement,
        state.normal_displacement,atol=3e-16)
    np.testing.assert_allclose(step.update.thickness,state.thickness,atol=3e-16)
    assert step.update.activation[0]>.5 and step.update.activation[1]<.5
