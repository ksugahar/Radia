"""Required native contracts for exact and proposal-only HDiv responses."""

import numpy as np
import pytest

from radia.topology_optimization import (
    linearize_hdiv_mmm_element_generation,
    hdiv_mmm_block_insertion_response,
    hdiv_mmm_removal_group_responses,
    grow_hdiv_mmm_by_superposition,
)


def test_native_bdm1_directional_schur_matches_full_schur_projection():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram

    mesh=MakeStructured3DMesh(hexes=True,nx=2,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=1,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(
            fes,eps=1e-10,leafsize=256,eta=2.0,
            internal_interfaces=True)
    rng=np.random.default_rng(20260823)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(3,fes.ndof))
    active=np.array([True,False])
    common=dict(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        solve_tolerance=1e-11,candidate_batch_size=128)
    full=linearize_hdiv_mmm_element_generation(**common)
    directional=linearize_hdiv_mmm_element_generation(
        **common,candidate_selector=lambda elements,*_:elements,
        candidate_direction_reduction=True)

    assert directional.candidate_directional_reduction
    assert directional.candidate_coupling_rank<=1
    np.testing.assert_array_equal(directional.candidate_dof_offsets,[0,1])
    assert directional.candidate_ritz_directions is not None
    direction=directional.candidate_ritz_directions[0]
    np.testing.assert_allclose(np.linalg.norm(direction),1.0,rtol=1e-13)
    expected_schur=float(direction@full.reduced_schur_complement@direction)
    expected_rhs=float(direction@full.reduced_schur_rhs)
    expected_response=full.reduced_response_matrix@direction
    np.testing.assert_allclose(
        directional.reduced_schur_complement,[[expected_schur]],
        rtol=2e-9,atol=2e-11)
    np.testing.assert_allclose(
        directional.reduced_schur_rhs,[expected_rhs],
        rtol=2e-9,atol=2e-11)
    np.testing.assert_allclose(
        directional.reduced_response_matrix[:,0],expected_response,
        rtol=2e-9,atol=2e-11)
    np.testing.assert_allclose(
        directional.candidate_response_delta[:,0],
        expected_response*(expected_rhs/expected_schur),
        rtol=3e-9,atol=3e-11)
    directional_bundle=hdiv_mmm_block_insertion_response(
        directional,directional.candidate_elements)
    assert directional_bundle.candidate_directional_reduction
    with pytest.raises(ValueError,match="full-block candidate linearization"):
        hdiv_mmm_removal_group_responses(
            directional,[directional.candidate_elements],
            full_elements=directional.candidate_elements)
    partial=linearize_hdiv_mmm_element_generation(
        **common,candidate_selector=lambda elements,*_:elements,
        candidate_direction_reduction=True,
        screen_adjoint_rows=np.array([0,2],dtype=np.int64),
        screen_response_band=np.ones(3))
    assert len(partial.adjoint_iterations)==2
    partial_reference=linearize_hdiv_mmm_element_generation(
        **common,active_state=partial.state)
    partial_direction=partial.candidate_ritz_directions[0]
    partial_expected=(
        partial_reference.reduced_response_matrix@partial_direction)
    np.testing.assert_allclose(
        partial.reduced_response_matrix[[0,2],0],
        partial_expected[[0,2]],
        rtol=3e-9,atol=3e-11)
    assert np.all(np.isfinite(partial.reduced_response_matrix))


def test_hdiv_mmm_nonlinear_transform_may_reduce_response_dimension():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram

    mesh=MakeStructured3DMesh(hexes=False,nx=2,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(
            fes,eps=1e-10,leafsize=256,eta=2.0,
            internal_interfaces=True)
    rng=np.random.default_rng(20260831)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(2,fes.ndof))
    active=np.zeros(mesh.ne,dtype=bool);active[0]=True
    initial=linearize_hdiv_mmm_element_generation(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        solve_tolerance=1e-11)
    inserted=initial.response+initial.candidate_response_delta[:,0]
    transform=lambda values:np.array([values[0]+0.25*values[1]**2])
    target=transform(inserted)
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))

    result=grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        element_volumes=volumes,response_target=target,
        response_band=[1e-8],volume_max=float(np.sum(volumes))+1e-14,
        maximum_batch_elements=1,max_iterations=1,
        solve_tolerance=1e-11,response_transform=transform)

    assert result.converged
    np.testing.assert_allclose(result.objective_response,target,rtol=0,atol=4e-12)
