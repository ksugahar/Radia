from pathlib import Path

import numpy as np
import ngsolve as ng

from radia.magnetic_shield_optimization import (MagneticShieldDesign,
    linearize_and_step_production_shield_from_ngsolve,
    linearize_and_step_production_shield_streaming,
    linearize_and_step_shield, linearize_and_step_shield_from_ngsolve,
    linearize_shield_response, rms_response_gradient, solve_shield)
from radia.topology_optimization import VIMOperatorLinearization


def test_streaming_contractions_serialize_shared_charge_gram_modes():
    """Direction workers must not reenter one native ChargeGram instance."""
    source=(Path(__file__).resolve().parents[1]/"src"/"core"/
            "rad_hacapk_hdiv.cpp").read_text(encoding="utf-8")
    start=source.index("RadHACApKChargeGram::DirectionalDerivativeContractions(")
    end=source.index(
        "RadHACApKChargeGram::DirectionalDerivativeContractionsMany(", start)
    implementation=source[start:end]
    assert "ParallelFor" not in implementation
    assert "for(int kk=0;kk<nDirections;++kk)" in implementation


def test_parallel_field_shield_improves_with_thickness():
    with ng.TaskManager():
        thin=solve_shield(MagneticShieldDesign(thickness=.001,nx=3,ny=3))
        thick=solve_shield(MagneticShieldDesign(thickness=.003,nx=3,ny=3))
    assert thin.shielding_factor>1
    assert thick.rms_parallel_h<thin.rms_parallel_h
    assert thick.shielding_factor>thin.shielding_factor


def test_rms_gradient_uses_supplied_analytic_vim_jacobian():
    response=np.array([3.,4.]); jac=np.array([[1.,0.,2.],[0.,2.,1.]])
    observed=rms_response_gradient(response,jac)
    expected=(response@jac)/(2*np.sqrt(np.mean(response**2)))
    np.testing.assert_allclose(observed,expected)


def test_shield_operator_tangent_closes_into_probe_response():
    A=np.array([[2.,.2],[.2,1.5]]); b=np.array([1.,.5]); dA=np.array([[[.1,0],[0,-.05]]]); db=np.array([[.02,-.01]])
    operator=VIMOperatorLinearization(A,b,dA,db); C=np.array([[1.,0.],[0.,1.]])
    model=linearize_shield_response(operator,C,incident_response=[3.,4.])
    eps=1e-7; shifted=C@np.linalg.solve(A+eps*dA[0],b+eps*db[0])+np.array([3.,4.])
    np.testing.assert_allclose((shifted-model.response)/eps,model.response_jacobian[:,0],rtol=2e-6,atol=2e-8)


def test_analytic_geometry_to_shield_lp_is_closed_without_fd_in_optimizer():
    cells=[np.array([[0.,0,0],[1.,0,0],[0,1.,0],[0,0,1.]]),
           np.array([[2.,0,0],[3.,0,0],[2,1.,0],[2,0,1.]])]
    modes=[]
    for c in cells:
        modes.append(np.stack([.02*c,np.tile([0,.01,0],(4,1)),np.zeros_like(c)]))
    points=np.array([c.mean(axis=0) for c in cells]); weights=np.array([1/6,1/6])
    point_modes=np.stack([[m[k].mean(axis=0) for m in modes] for k in range(3)])
    result=linearize_and_step_shield(points=points,weights=weights,
        displacement_modes=point_modes,relative_weight_derivatives=np.array([[.06,.06],[0,0],[0,0]]),
        self_cell_types=["tet","tet"],self_nodes=cells,self_node_displacements=modes,
        mass=np.eye(2),charge_map=np.eye(2),applied_coefficients=np.array([2.,1.]),inv_chi=.01,
        dmass=np.zeros((3,2,2)),dcharge_map=np.zeros((3,2,2)),response_matrix=np.eye(2),
        incident_response=np.array([3.,3.]),normal_displacement=[0.],thickness=[.002],activation=[1.],
        cell_areas=[.01],volume_max=2e-5,displacement_move=1e-3,thickness_move=2e-4,
        thickness_bounds=(.001,.003))
    assert result.charge_gram.matrix[0,0]>0
    assert result.operator.matrix_jacobian.shape==(3,2,2)
    assert result.response.response_jacobian.shape==(2,3)
    assert result.rms_gradient.shape==(3,)
    assert result.update.delta.shape==(3,)


def test_real_gettrafo_hdiv_mass_to_reduced_shield_lp_closure():
    from ngsolve.meshes import MakeStructured3DMesh
    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=1); vf=ng.VectorH1(mesh,order=1)
    modes=[]
    for field in ((.01*ng.x,0,0),(0,.01*ng.y,0),(0,0,.01*ng.z)):
        gf=ng.GridFunction(vf); gf.Set(ng.CF(field)); modes.append(gf)
    B=np.zeros((1,fes.ndof)); B[0,0]=1
    with ng.TaskManager():
        result=linearize_and_step_shield_from_ngsolve(fes=fes,deformation_modes=modes,
            charge_map=B,weights=[1.],relative_weight_derivatives=np.full((3,1),.01),
            applied_coefficients=np.ones(fes.ndof),inv_chi=.01,
            response_matrix=np.ones((1,fes.ndof)),incident_response=[2.],
            normal_displacement=[0.],thickness=[.002],activation=[1.],cell_areas=[.01],
            volume_max=2e-5,displacement_move=1e-3,thickness_move=2e-4,
            thickness_bounds=(.001,.003))
    assert result.operator.matrix.shape==(fes.ndof,fes.ndof)
    assert result.update.delta.shape==(3,)


def test_real_production_hex_gettrafo_to_shield_lp_closure():
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis_hex,build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1,
        mapping=lambda x,y,z:(x+.03*y*z,y+.02*x*z,.1*z+.01*x*y))
    fes=ng.HDiv(mesh,order=1); vf=ng.VectorH1(mesh,order=1)
    modes=[]
    for field in ((.01*ng.x,0,0),(0,.01*ng.y,0),(0,0,.01*ng.z)):
        gf=ng.GridFunction(vf);gf.Set(ng.CF(field));modes.append(gf)
    with ng.TaskManager():
        cb=_charge_basis_hex(fes,cob_quad=3)
        B,gram,_=build_charge_gram(fes,eps=1e-10,leafsize=256,eta=2.)
        result=linearize_and_step_production_shield_from_ngsolve(fes=fes,
            deformation_modes=modes,charge_basis=cb,charge_gram=gram,
            charge_map=B,applied_coefficients=np.ones(fes.ndof),inv_chi=.01,
            response_matrix=np.ones((1,fes.ndof)),incident_response=[2.],
            normal_displacement=[0.],thickness=[.002],activation=[1.],
            cell_areas=[.01],volume_max=2e-5,displacement_move=1e-3,
            thickness_move=2e-4,thickness_bounds=(.001,.003),family="hex")
        streaming=linearize_and_step_production_shield_streaming(fes=fes,
            deformation_modes=modes,charge_basis=cb,charge_gram=gram,
            charge_map=B,applied_coefficients=np.ones(fes.ndof),inv_chi=.01,
            response_matrix=np.ones((1,fes.ndof)),incident_response=[2.],
            normal_displacement=[0.],thickness=[.002],activation=[1.],
            cell_areas=[.01],volume_max=2e-5,displacement_move=1e-3,
            thickness_move=2e-4,thickness_bounds=(.001,.003),family="hex",
            eps=1e-12,leaf=256,eta=2.)
    assert result.charge_gram.jacobian.shape[0]==3
    assert result.operator.matrix_jacobian.shape==(3,fes.ndof,fes.ndof)
    assert result.response.response_jacobian.shape==(1,3)
    assert result.update.delta.shape==(3,)
    np.testing.assert_allclose(streaming.adjoint.gradient,result.rms_gradient,
        rtol=2e-7,atol=2e-9)
    assert streaming.adjoint.peak_directional_operators==0
    assert streaming.update.delta.shape==(3,)
