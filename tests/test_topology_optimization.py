import numpy as np
import pytest

from radia.topology_optimization import (
    affine_cell_self_energy_shape_derivative,
    assemble_ngsolve_hdiv_shape_tangents,
    assemble_ngsolve_hdiv_linear_form_shape_tangents,
    assemble_ngsolve_hdiv_mass_shape_contractions,
    linearize_laplace_charge_gram,
    production_hex_volume_self_block_derivatives,
    production_hex_face_self_block_derivatives,
    production_wedge_volume_self_block_derivatives,
    production_wedge_face_self_block_derivatives,
    production_wedge_charge_gram_derivatives,
    production_tet_volume_self_block_derivatives,
    production_tet_face_self_block_derivatives,
    linearize_laplace_pair_gram, linearize_vim_operator,
    VIMLinearization,
    linearize_vim_system,
    optimize_vim_lp,
    solve_lp_update,
    write_cubit_density_journal,
    linearize_production_vim_from_ngsolve,
    linearize_production_vim_matrix_free_from_ngsolve,
    sample_production_gettrafo_displacements,
    production_vim_functional_shape_jacobian_streaming,
    ElementInsertionResponse,
    ShapeLinearization,
    finite_element_insertion_response,
    hdiv_mmm_all_single_removal_responses,
    hdiv_mmm_block_insertion_response,
    linearize_hdiv_mmm_element_generation,
    grow_hdiv_mmm_by_superposition,
    ngsolve_boundary_growth_candidates,
    ngsolve_boundary_removal_candidates,
    ngsolve_growth_topology,
    ngsolve_discontinuous_element_dof_blocks,
    select_collaborative_element_batch,
    select_tsvd_element_candidates,
    select_tsvd_exact_block_batch,
    solve_hdiv_mmm_active_elements,
    solve_element_generation_lp,
    solve_shape_lp,
)


def test_empty_insertion_side_tsvd_front_is_a_clean_no_update():
    update=select_tsvd_exact_block_batch(
        current_response=[2.0],response_target=[0.0],response_band=[1.0],
        candidate_elements=[3],candidate_volumes=[1.0],
        proposal_elements=[],representative_elements=[],
        evaluate_bundle_response=lambda _:np.array([0.0]),
        volume_budget=1.0)
    assert update.selected_elements.size==0
    assert update.evaluated_bundles==0
    assert update.predicted_max_band_ratio==pytest.approx(2.0)
    assert "no insertion-side" in update.status

    deferred_removal=select_tsvd_exact_block_batch(
        current_response=[2.0],response_target=[0.0],response_band=[1.0],
        candidate_elements=[],candidate_volumes=[],
        proposal_elements=[3],representative_elements=[3],
        evaluate_bundle_response=lambda _:np.array([0.0]),
        volume_budget=1.0)
    assert deferred_removal.selected_elements.size==0
    assert deferred_removal.evaluated_bundles==0


def test_screened_adjoint_correction_is_lifted_through_candidate_row_space():
    import radia.topology_optimization as topopt
    skeleton=np.array([[1.0,2.0,-1.0],[0.5,-1.0,3.0]])
    interpolation=np.array([[1.0,0.0],[0.0,1.0],[2.0,-.5],[-1.0,3.0]])
    direct=interpolation@skeleton
    sampled_correction=np.array([[.2,-.1,.4],[-.3,.5,.1]])
    expected=direct+interpolation@sampled_correction
    partial=direct.copy();partial[[0,1]]+=sampled_correction
    actual=topopt._interpolate_screened_response_correction(
        direct,partial,[0,1],np.ones(4))
    np.testing.assert_allclose(actual,expected,rtol=0,atol=2e-14)


def _element_centroids(mesh):
    import ngsolve as ng
    return np.asarray([
        np.mean([np.asarray(mesh[vertex].point,dtype=float)
                 for vertex in element.vertices],axis=0)
        for element in mesh.Elements(ng.VOL)])


def test_growth_topology_rejects_cavity_and_disconnected_iron():
    from ngsolve.meshes import MakeStructured3DMesh
    mesh=MakeStructured3DMesh(hexes=True,nx=3,ny=3,nz=3)
    centers=_element_centroids(mesh)
    center=int(np.argmin(np.linalg.norm(centers-.5,axis=1)))
    cavity=np.ones(mesh.ne,dtype=bool);cavity[center]=False
    report=ngsolve_growth_topology(mesh,cavity)
    assert report.iron_connected and not report.inactive_reaches_exterior
    np.testing.assert_array_equal(report.enclosed_inactive_elements,[center])
    assert not report.valid

    corners=np.zeros(mesh.ne,dtype=bool)
    corners[int(np.argmin(np.sum(centers,axis=1)))]=True
    corners[int(np.argmax(np.sum(centers,axis=1)))]=True
    report=ngsolve_growth_topology(mesh,corners)
    assert len(report.iron_components)==2 and not report.iron_connected
    assert report.inactive_reaches_exterior and not report.valid

    slab=centers[:,0]<.34
    report=ngsolve_growth_topology(mesh,slab)
    assert report.valid


def test_growth_candidates_respect_fixed_air_and_column_predecessor():
    from ngsolve.meshes import MakeStructured3DMesh
    mesh=MakeStructured3DMesh(hexes=True,nx=3,ny=1,nz=1)
    centers=_element_centroids(mesh);order=np.argsort(centers[:,0])
    active=np.zeros(mesh.ne,dtype=bool);active[order[0]]=True
    predecessor=np.full(mesh.ne,-1,dtype=np.int64)
    predecessor[order[1]]=order[0];predecessor[order[2]]=order[1]
    candidates=ngsolve_boundary_growth_candidates(
        mesh,active,predecessor_elements=predecessor)
    np.testing.assert_array_equal(candidates,[order[1]])
    descendants=ngsolve_boundary_growth_candidates(
        mesh,active,predecessor_elements=predecessor,
        include_predecessor_descendants=True)
    np.testing.assert_array_equal(descendants,order[1:])
    fixed=np.zeros(mesh.ne,dtype=bool);fixed[order[1]]=True
    assert ngsolve_boundary_growth_candidates(
        mesh,active,fixed_inactive_elements=fixed,
        predecessor_elements=predecessor).size==0


def test_binary_generation_lp_enforces_predecessor_selection():
    update=solve_element_generation_lp(
        current_response=[0.0],response_target=[0.0],response_band=[1.0],
        candidate_response_delta=np.zeros((1,2)),candidate_volumes=[1.0,1.0],
        volume_budget=2.0,maximum_new_elements=2,
        candidate_objective_change=[0.0,-1.0],predecessor_pairs=[(1,0)])
    np.testing.assert_array_equal(update.selected,[True,True])


def test_binary_generation_lp_enforces_mutually_exclusive_terminal_states():
    update=solve_element_generation_lp(
        current_response=[2.0],response_target=[0.0],response_band=[1.0],
        candidate_response_delta=[[-1.1,-1.2]],candidate_volumes=[1.0,2.0],
        volume_budget=3.0,maximum_new_elements=2,
        candidate_objective_change=[-1.0,-1.0],
        candidate_exclusion_groups=[7,7])
    assert np.count_nonzero(update.selected)==1
    assert update.predicted_max_band_ratio>=0.8-1e-10


def test_binary_generation_lp_separates_net_volume_from_flip_trust_region():
    unrestricted=solve_element_generation_lp(
        current_response=[2.0,3.0],response_target=[0.0,0.0],
        response_band=[1.0,1.0],
        candidate_response_delta=[[-2.0,0.0],[0.0,-3.0]],
        candidate_volumes=[1.0,1.0],volume_budget=0.0,
        candidate_volume_change=[1.0,-1.0])
    np.testing.assert_array_equal(unrestricted.selected,[True,True])
    assert unrestricted.predicted_max_band_ratio==pytest.approx(0.0)

    trusted=solve_element_generation_lp(
        current_response=[2.0,3.0],response_target=[0.0,0.0],
        response_band=[1.0,1.0],
        candidate_response_delta=[[-2.0,0.0],[0.0,-3.0]],
        candidate_volumes=[1.0,1.0],volume_budget=0.0,
        candidate_volume_change=[1.0,-1.0],
        maximum_changed_volume=1.0)
    np.testing.assert_array_equal(trusted.selected,[False,True])
    assert trusted.predicted_max_band_ratio==pytest.approx(2.0)


def test_ngsolve_hdiv_linear_form_shape_tangent_includes_spatial_kernel_motion():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    mesh=MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=1);space=ng.VectorH1(mesh,order=1)
    velocity=ng.GridFunction(space)
    velocity.Set(ng.CF((.1*ng.x,.02*ng.y,-.03*ng.z)))
    coefficient=ng.CF((ng.x*ng.x,ng.y*ng.z,ng.z+.2*ng.x))
    with ng.TaskManager():
        _,analytic=assemble_ngsolve_hdiv_linear_form_shape_tangents(
            fes,coefficient,[velocity])
        values=[];epsilon=2e-6
        for sign in (1.0,-1.0):
            deformation=ng.GridFunction(space)
            deformation.vec.data=sign*epsilon*velocity.vec
            mesh.SetDeformation(deformation)
            form=ng.LinearForm(fes)
            form+=ng.InnerProduct(coefficient,fes.TestFunction())*ng.dx
            form.Assemble();values.append(form.vec.FV().NumPy().copy())
            mesh.UnsetDeformation()
    fd=(values[0]-values[1])/(2*epsilon)
    np.testing.assert_allclose(analytic[0],fd,rtol=2e-8,atol=2e-10)


def test_finite_filament_hdiv_load_shape_tangent_matches_rebuild():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.biot_savart import h_segments_cf
    mesh=MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=1);space=ng.VectorH1(mesh,order=1)
    velocity=ng.GridFunction(space)
    velocity.Set(ng.CF((.04*ng.x-.01*ng.y,.02*ng.y,-.03*ng.z)))
    coefficient=h_segments_cf([
        ((-0.3,-0.2,1.4),(1.2,-0.2,1.4)),
        ((1.2,-0.2,1.4),(1.2,1.3,1.4)),
    ],current=1200.0)
    with ng.TaskManager():
        _,analytic=assemble_ngsolve_hdiv_linear_form_shape_tangents(
            fes,coefficient,[velocity],bonus_intorder=10)
        values=[];epsilon=2e-6
        for sign in (1.0,-1.0):
            deformation=ng.GridFunction(space)
            deformation.vec.data=sign*epsilon*velocity.vec
            mesh.SetDeformation(deformation)
            form=ng.LinearForm(fes)
            form+=ng.InnerProduct(coefficient,fes.TestFunction())*ng.dx(
                bonus_intorder=10)
            form.Assemble();values.append(form.vec.FV().NumPy().copy())
            mesh.UnsetDeformation()
    fd=(values[0]-values[1])/(2*epsilon)
    np.testing.assert_allclose(analytic[0],fd,rtol=3e-8,atol=2e-9)


def test_ngsolve_hdiv_mass_shape_contractions_match_sparse_directional_matrices():
    import ngsolve as ng
    import scipy.sparse as sp
    from ngsolve.meshes import MakeStructured3DMesh
    mesh=MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=1);space=ng.VectorH1(mesh,order=1)
    modes=[]
    for coefficient in ((.08,-.03,.02),(-.01,.04,.06)):
        mode=ng.GridFunction(space)
        mode.Set(ng.CF((coefficient[0]*ng.x,coefficient[1]*ng.y,
                        coefficient[2]*ng.z)))
        modes.append(mode)
    rng=np.random.default_rng(20260731)
    right=rng.normal(size=fes.ndof);left=rng.normal(size=(3,fes.ndof))
    with ng.TaskManager():
        _,dmass,_=assemble_ngsolve_hdiv_shape_tangents(
            fes,modes,sp.eye(fes.ndof),sparse=True)
        contractions=assemble_ngsolve_hdiv_mass_shape_contractions(
            fes,modes,left,right)
    expected=np.asarray([[row@(matrix@right) for matrix in dmass]
                         for row in left])
    np.testing.assert_allclose(contractions,expected,rtol=3e-13,atol=3e-13)


def test_hex_charge_basis_lattice_follows_live_ngsolve_deformation():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis_hex
    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=1)
    space=ng.VectorH1(mesh,order=1)
    deformation=ng.GridFunction(space)
    deformation.Set(ng.CF((.07*ng.x-.02*ng.y,.03*ng.y,.04*ng.z+.01*ng.x)))
    with ng.TaskManager():
        undeformed=_charge_basis_hex(fes,cob_quad=3,materialize_mass=False)
        mesh.SetDeformation(deformation)
        deformed=_charge_basis_hex(fes,cob_quad=3,materialize_mass=False)
        mesh.UnsetDeformation()
    x0=np.asarray(undeformed["cell_nodes"]).reshape(-1,3)
    x1=np.asarray(deformed["cell_nodes"]).reshape(-1,3)
    expected=x0+np.column_stack((.07*x0[:,0]-.02*x0[:,1],
                                 .03*x0[:,1],
                                 .04*x0[:,2]+.01*x0[:,0]))
    np.testing.assert_allclose(x1,expected,rtol=2e-13,atol=2e-13)


def test_tet_charge_basis_vertices_follow_live_ngsolve_deformation():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis
    mesh=MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=1)
    space=ng.VectorH1(mesh,order=1)
    deformation=ng.GridFunction(space)
    deformation.Set(ng.CF((.07*ng.x-.02*ng.y,.03*ng.y,
                           .04*ng.z+.01*ng.x)))
    with ng.TaskManager():
        undeformed=_charge_basis(fes,quad=4,materialize_mass=False)
        mesh.SetDeformation(deformation)
        deformed=_charge_basis(fes,quad=4,materialize_mass=False)
        mesh.UnsetDeformation()
    for key in ("vV","bV"):
        x0=np.asarray(undeformed[key]).reshape(-1,3)
        x1=np.asarray(deformed[key]).reshape(-1,3)
        expected=x0+np.column_stack((.07*x0[:,0]-.02*x0[:,1],
                                     .03*x0[:,1],
                                     .04*x0[:,2]+.01*x0[:,0]))
        np.testing.assert_allclose(x1,expected,rtol=2e-13,atol=2e-13)


def test_multi_functional_streaming_shape_jacobian_matches_rebuilt_geometry():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram,_charge_basis_hex
    gradient=np.array([[.04,-.02,.01],[.01,.03,-.015],[-.02,.005,.025]])
    shift=np.array([.01,-.006,.004]);inv_chi=.2
    applied=ng.CF((.2,-.1,1.0))
    observed=ng.CF((ng.x*ng.x+.1*ng.z,.3*ng.y,ng.z+.2*ng.x))
    def build(epsilon):
        def mapping(x,y,z):
            point=np.array([x,y,z]);return tuple(point+epsilon*(gradient@point+shift))
        mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1,mapping=mapping)
        fes=ng.HDiv(mesh,order=1)
        with ng.TaskManager():
            basis=_charge_basis_hex(fes,cob_quad=3)
            B,gram,_=build_charge_gram(fes,eps=1e-12,leafsize=256,eta=2.0)
            rhs,_=assemble_ngsolve_hdiv_linear_form_shape_tangents(
                fes,applied,(),bonus_intorder=4)
            row,_=assemble_ngsolve_hdiv_linear_form_shape_tangents(
                fes,observed,(),bonus_intorder=4)
        solved=gram.solve_configured_linear_material_auto_prec_many(inv_chi,
            np.ascontiguousarray(rhs[None,:]),tol=1e-12,maxit=5000,
            mass_riesz=True)["m"][0]
        return float(row@solved),(mesh,fes,basis,B,gram,rhs,row)
    value,data=build(0.0);mesh,fes,basis,B,gram,rhs,row=data
    space=ng.VectorH1(mesh,order=1);mode=ng.GridFunction(space)
    mode.Set(ng.CF(tuple(gradient@np.array([ng.x,ng.y,ng.z],dtype=object)+shift)))
    with ng.TaskManager():
        _,drhs=assemble_ngsolve_hdiv_linear_form_shape_tangents(
            fes,applied,[mode],bonus_intorder=4)
        _,drow=assemble_ngsolve_hdiv_linear_form_shape_tangents(
            fes,observed,[mode],bonus_intorder=4)
        analytic=production_vim_functional_shape_jacobian_streaming(
            fes=fes,deformation_modes=[mode],charge_basis=basis,
            charge_gram=gram,charge_map=B,inv_chi=inv_chi,rhs=rhs,
            response_matrix=row[None,:],rhs_jacobian=drhs,
            dresponse_matrix=drow[:,None,:],family="hex",solve_tolerance=1e-12)
    epsilon=2e-6;plus,_=build(epsilon);minus,_=build(-epsilon)
    fd=(plus-minus)/(2*epsilon)
    np.testing.assert_allclose(analytic.response,[value],rtol=2e-12,atol=2e-12)
    np.testing.assert_allclose(analytic.response_jacobian[0,0],fd,rtol=2e-4,atol=2e-7)


def test_finite_element_insertion_response_is_exact_for_full_strength_block():
    rng=np.random.default_rng(20260730)
    R=rng.normal(size=(7,7)); Aaa=R.T@R+4*np.eye(7)
    Aae=.08*rng.normal(size=(7,3));
    Ree=rng.normal(size=(3,3)); Aee=Ree.T@Ree+3*np.eye(3)
    ba=rng.normal(size=7);be=rng.normal(size=3)
    Ca=rng.normal(size=(4,7));Ce=rng.normal(size=(4,3))
    ma=np.linalg.solve(Aaa,ba)
    adjoint=np.linalg.solve(Aaa,Ca.T)
    result=finite_element_insertion_response(
        solve_active=lambda rhs:np.linalg.solve(Aaa,rhs),active_state=ma,
        active_to_candidate=Aae,candidate_matrix=Aee,candidate_rhs=be,
        active_response_matrix=Ca,candidate_response_matrix=Ce,
        active_adjoint=adjoint)
    enlarged=np.block([[Aaa,Aae],[Aae.T,Aee]])
    full=np.linalg.solve(enlarged,np.r_[ba,be])
    expected=(np.c_[Ca,Ce]@full)-(Ca@ma)
    assert isinstance(result,ElementInsertionResponse)
    np.testing.assert_allclose(result.active_state_delta,full[:7]-ma,rtol=2e-13,atol=2e-14)
    np.testing.assert_allclose(result.candidate_state,full[7:],rtol=2e-13,atol=2e-14)
    np.testing.assert_allclose(result.response_delta,expected,rtol=2e-13,atol=2e-14)


def test_whole_element_generation_lp_has_no_gray_material():
    update=solve_element_generation_lp([0.0],[2.0],[.05],
        np.array([[1.0,1.0,.2]]),[1.0,1.0,.1],volume_budget=2.0,
        maximum_new_elements=2,whole_elements=True)
    assert update.selected.tolist()==[True,True,False]
    assert set(update.weights.tolist())<={0.0,1.0}
    np.testing.assert_allclose(update.predicted_response,[2.0],atol=1e-12)
    assert update.added_volume==2.0 and update.predicted_max_band_ratio<1e-10


def test_native_hdiv_mmm_boundary_insertions_match_exact_active_resolves():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram
    mesh=MakeStructured3DMesh(hexes=False,nx=2,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(fes,eps=1e-10,leafsize=256,
            eta=2.0,internal_interfaces=True)
    assert hasattr(gram,"reduce_configured_candidate_schur")
    rng=np.random.default_rng(9182)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(2,fes.ndof))
    active=np.zeros(mesh.ne,dtype=bool);active[0]=True
    result=linearize_hdiv_mmm_element_generation(charge_gram=gram,fes=fes,
        inv_chi=.2,rhs=rhs,response_matrix=response_matrix,
        active_elements=active,solve_tolerance=1e-11,candidate_batch_size=4)
    probes=np.ascontiguousarray(rng.normal(size=(3,fes.ndof)))
    for respect_constraints in (False,True):
        batched=gram.apply_configured_linear_material_operator_many(
            .2,probes,respect_constraints=respect_constraints)
        scalar=np.stack([
            gram.apply_configured_linear_material_operator(
                .2,row,respect_constraints=respect_constraints)
            for row in probes])
        np.testing.assert_allclose(batched,scalar,rtol=2e-13,atol=2e-13)
        assert np.asarray(batched).flags.c_contiguous
    blocks=ngsolve_discontinuous_element_dof_blocks(fes)
    active_dofs=blocks[0]
    candidate_blocks=tuple(blocks[int(element)]
                           for element in result.candidate_elements)
    candidate_dofs=np.concatenate(candidate_blocks).astype(np.int32)
    block_offsets=np.asarray(np.r_[0,np.cumsum(
        [len(block) for block in candidate_blocks])],dtype=np.int32)
    packed=np.asarray(gram.configured_linear_material_element_blocks(
        .2,candidate_dofs,block_offsets))
    packed_offset=0
    for block in candidate_blocks:
        basis=np.zeros((len(block),fes.ndof),dtype=float)
        basis[np.arange(len(block)),block]=1.0
        applied=np.asarray(
            gram.apply_configured_linear_material_operator_many(
                .2,np.ascontiguousarray(basis),
                respect_constraints=False))
        expected=applied[:,block].T
        next_offset=packed_offset+len(block)**2
        local=packed[packed_offset:next_offset].reshape(
            len(block),len(block))
        packed_offset=next_offset
        np.testing.assert_allclose(local,expected,rtol=3e-13,atol=3e-13)
    assert result.candidate_response_delta.shape==(2,len(result.candidate_elements))
    assert result.available_candidate_count==len(result.candidate_elements)
    assert set(result.native_reduction_timings)=={
        "operator_s","solve_s","contraction_s"}
    assert all(value>0.0 for value in result.native_reduction_timings.values())
    for column,element in enumerate(result.candidate_elements):
        keep=np.r_[active_dofs,blocks[int(element)]]
        constrained=np.ones(fes.ndof,dtype=bool);constrained[keep]=False
        gram.set_configured_constraints(
            np.flatnonzero(constrained).astype(np.int32))
        exact_rhs=rhs.copy();exact_rhs[constrained]=0.0
        exact=gram.solve_configured_linear_material_auto_prec_many(.2,
            np.ascontiguousarray(exact_rhs[None,:]),tol=1e-11,maxit=5000,
            mass_riesz=True)["m"][0]
        delta=response_matrix@exact-result.response
        np.testing.assert_allclose(result.candidate_response_delta[:,column],
            delta,rtol=2e-12,atol=2e-12)

    selected=result.candidate_elements[:1]
    screened=linearize_hdiv_mmm_element_generation(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        candidate_elements=result.candidate_elements,solve_tolerance=1e-11,
        candidate_batch_size=4,
        candidate_selector=lambda elements,delta,state,response:selected)
    assert screened.available_candidate_count==len(result.candidate_elements)
    np.testing.assert_array_equal(screened.candidate_elements,selected)
    np.testing.assert_allclose(screened.state,result.state,rtol=2e-13,atol=2e-13)
    np.testing.assert_allclose(screened.candidate_response_delta[:,0],
                               result.candidate_response_delta[:,0],
                               rtol=2e-12,atol=2e-12)


def test_large_hdiv_mmm_candidate_screen_uses_one_shared_hmatrix_batch():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram

    mesh=MakeStructured3DMesh(hexes=True,nx=5,ny=5,nz=1)
    # Production Lego growth uses broken mapped-HEX BDM1 (36 DOFs/cell).
    fes=ng.HDiv(mesh,order=1,discontinuous=True)
    with ng.TaskManager():
        _,native,mass=build_charge_gram(
            fes,eps=1e-10,leafsize=256,eta=2.0,
            internal_interfaces=True)

    class CountingChargeGram:
        def __init__(self,inner):
            self.inner=inner
            self.operator_batches=[]

        def __getattr__(self,name):
            return getattr(self.inner,name)

        def apply_configured_linear_material_operator_many(
                self,inv_chi,rows,*,respect_constraints=True):
            self.operator_batches.append(
                (int(np.asarray(rows).shape[0]),bool(respect_constraints)))
            return self.inner.apply_configured_linear_material_operator_many(
                inv_chi,rows,respect_constraints=respect_constraints)

    gram=CountingChargeGram(native)
    active=np.zeros(mesh.ne,dtype=bool)
    active[[6,7,8,11,12,13,16,17,18]]=True
    candidates=ngsolve_boundary_growth_candidates(mesh,active)
    assert len(candidates)==12
    blocks=ngsolve_discontinuous_element_dof_blocks(fes)
    probe=blocks[int(candidates[0])]
    local=np.asarray(native.configured_linear_material_element_blocks(
        .2,np.asarray(probe,dtype=np.int32),
        np.asarray([0,len(probe)],dtype=np.int32))).reshape(
            len(probe),len(probe))
    basis=np.zeros((len(probe),fes.ndof),dtype=float)
    basis[np.arange(len(probe)),probe]=1.0
    applied=np.asarray(native.apply_configured_linear_material_operator_many(
        .2,np.ascontiguousarray(basis),respect_constraints=False))
    np.testing.assert_allclose(local,applied[:,probe].T,
                               rtol=3e-13,atol=3e-13)
    rng=np.random.default_rng(20260807)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(2,fes.ndof))
    active_state,_,_=solve_hdiv_mmm_active_elements(
        charge_gram=native,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        solve_tolerance=1e-11)
    screened=linearize_hdiv_mmm_element_generation(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        candidate_elements=candidates,solve_tolerance=1e-11,
        candidate_batch_size=64,
        candidate_selector=lambda elements,delta,state,response:[],
        active_state=active_state,screen_with_adjoint=False)
    assert screened.candidate_elements.size==0
    assert screened.state_iterations==0
    assert screened.adjoint_iterations==tuple()
    # The reused state is the only unconstrained global H-matrix apply.  Exact
    # candidate A_ee blocks come from their local charge supports; production
    # proposal adjoints are deferred to a rejected-proposal fallback.
    unconstrained=[count for count,respect in gram.operator_batches
                   if not respect]
    assert unconstrained==[1]
    gram.operator_batches.clear()
    working=linearize_hdiv_mmm_element_generation(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        candidate_elements=candidates,solve_tolerance=1e-11,
        candidate_batch_size=64,
        candidate_selector=lambda elements,delta,state,response:[],
        active_state=active_state,screen_with_adjoint=True,
        screen_adjoint_rows=[1])
    assert working.state_iterations==0
    assert len(working.adjoint_iterations)==1
    unconstrained=[count for count,respect in gram.operator_batches
                   if not respect]
    assert unconstrained==[2]


def test_native_hdiv_mmm_block_schur_bundle_matches_exact_active_resolve():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=3,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(fes,eps=1e-10,leafsize=256,
            eta=2.0,internal_interfaces=True)
    rng=np.random.default_rng(20260731)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(3,fes.ndof))
    active=np.zeros(mesh.ne,dtype=bool);active[1]=True
    linearization=linearize_hdiv_mmm_element_generation(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        solve_tolerance=1e-11,candidate_batch_size=4)
    selected=linearization.candidate_elements
    assert selected.size==2
    block=hdiv_mmm_block_insertion_response(linearization,selected)
    exact_state,exact_response,_=solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,
        active_elements=np.ones(mesh.ne,dtype=bool),solve_tolerance=1e-11)
    np.testing.assert_allclose(
        linearization.response+block.response_delta,exact_response,
        rtol=3e-12,atol=3e-12)
    assert block.schur_complement.shape[0]==sum(
        len(linearization.candidate_dof_blocks[index])
        for index in range(len(selected)))
    assert np.all(np.isfinite(exact_state))


def test_native_hdiv_mmm_single_removal_oracle_reuses_one_schur_inverse():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=3,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(fes,eps=1e-10,leafsize=256,
            eta=2.0,internal_interfaces=True)
    rng=np.random.default_rng(20260821)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(4,fes.ndof))
    retained=np.zeros(mesh.ne,dtype=bool);retained[1]=True
    linearization=linearize_hdiv_mmm_element_generation(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=retained,
        solve_tolerance=1e-11,candidate_batch_size=8)
    candidates=linearization.candidate_elements
    oracle=hdiv_mmm_all_single_removal_responses(
        linearization,candidates)
    full=hdiv_mmm_block_insertion_response(linearization,candidates)
    np.testing.assert_allclose(
        oracle.full_response_delta,full.response_delta,
        rtol=3e-12,atol=3e-12)
    for column,removed in enumerate(candidates):
        kept=candidates[candidates!=removed]
        expected=hdiv_mmm_block_insertion_response(
            linearization,kept).response_delta
        np.testing.assert_allclose(
            oracle.removed_response_delta[:,column],expected,
            rtol=3e-12,atol=3e-12)
        np.testing.assert_allclose(
            oracle.full_response_delta-
            oracle.positive_material_response[:,column],expected,
            rtol=3e-12,atol=3e-12)
    assert oracle.positive_material_response.flags.c_contiguous


def test_native_bdm1_candidate_schur_tsvd_solves_only_charge_coupling_rank():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=2,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=1,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(fes,eps=1e-10,leafsize=256,
            eta=2.0,internal_interfaces=True)
    rng=np.random.default_rng(20260803)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(2,fes.ndof))
    active=np.array([True,False])
    linearization=linearize_hdiv_mmm_element_generation(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        solve_tolerance=1e-12,candidate_batch_size=128)
    candidate_width=len(linearization.candidate_dof_blocks[0])
    assert candidate_width==36
    assert 0<linearization.candidate_coupling_rank<candidate_width
    assert len(linearization.schur_iterations)==\
        linearization.candidate_coupling_rank
    assert linearization.candidate_coupling_relative_truncation_error<1e-12
    candidate_dofs=np.concatenate(
        linearization.candidate_dof_blocks).astype(np.int32)
    raw=gram.reduce_configured_candidate_schur(
        .2,candidate_dofs,rhs,np.zeros(fes.ndof),response_matrix,
        np.zeros_like(response_matrix),tol=1e-10,maxit=5000,
        solve_batch_size=128,mass_riesz=True)
    assert len(raw["iters"])==candidate_width
    assert len(raw["coupling_mode_iters"])==raw["coupling_rank"]
    exact_state,exact_response,_=solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,
        active_elements=np.ones(mesh.ne,dtype=bool),solve_tolerance=1e-12)
    np.testing.assert_allclose(
        linearization.response+linearization.candidate_response_delta[:,0],
        exact_response,rtol=3e-10,atol=3e-10)
    assert np.all(np.isfinite(exact_state))


def test_native_candidate_schur_reports_zero_coupling_rank_and_stable_iters():
    import ngsolve as ng
    import pytest
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram

    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,_=build_charge_gram(
            fes,eps=1e-10,leafsize=256,eta=2.0,
            internal_interfaces=True)
    gram.set_configured_constraints(
        np.arange(fes.ndof,dtype=np.int32),preserve_existing=False)
    candidate=ngsolve_discontinuous_element_dof_blocks(fes)[0]
    kwargs=dict(
        inv_chi=.2,candidate_dofs=candidate,rhs=np.ones(fes.ndof),
        state=np.zeros(fes.ndof),response_matrix=np.zeros((1,fes.ndof)),
        adjoints=np.zeros((1,fes.ndof)),maxit=20,
        solve_batch_size=64,mass_riesz=True)
    raw=gram.reduce_configured_candidate_schur(tol=1e-10,**kwargs)
    assert raw["coupling_rank"]==0
    assert raw["coupling_mode_iters"]==[]
    assert raw["iters"]==[0]*len(candidate)
    with pytest.raises(RuntimeError,match="invalid material or solver parameters"):
        gram.reduce_configured_candidate_schur(tol=1.0,**kwargs)


def test_hdiv_generation_rejects_unconverged_native_multi_rhs():
    import pytest
    from ngsolve.meshes import MakeStructured3DMesh

    mesh=MakeStructured3DMesh(hexes=True,nx=2,ny=1,nz=1)
    import ngsolve as ng
    fes=ng.HDiv(mesh,order=0,discontinuous=True)

    class UnconvergedGram:
        def set_configured_constraints(self,*args,**kwargs):
            return None

        def solve_configured_linear_material_auto_prec_many(
                self,inv_chi,rows,**kwargs):
            return {"m":np.zeros_like(rows),"iters":[kwargs["maxit"]]*len(rows)}

        def apply_configured_linear_material_operator_many(
                self,inv_chi,rows,**kwargs):
            return np.zeros_like(rows)

    with pytest.raises(RuntimeError,match="did not meet its checked relative residual"):
        linearize_hdiv_mmm_element_generation(
            charge_gram=UnconvergedGram(),fes=fes,inv_chi=.2,
            rhs=np.ones(fes.ndof),response_matrix=np.zeros((1,fes.ndof)),
            active_elements=np.array([True,False]),solve_max_iterations=3)


def test_collaborative_search_finds_pair_when_neither_singleton_improves():
    responses={
        (10,):np.array([0.0,0.0]),
        (11,):np.array([0.0,0.0]),
        (12,):np.array([0.1,0.0]),
        (10,11):np.array([1.0,-1.0]),
        (10,12):np.array([0.1,0.0]),
        (11,12):np.array([0.1,0.0]),
    }
    update=select_collaborative_element_batch(
        current_response=np.zeros(2),response_target=np.array([1.0,-1.0]),
        response_band=np.full(2,.05),candidate_elements=np.array([10,11,12]),
        candidate_volumes=np.ones(3),
        evaluate_bundle_response=lambda bundle:responses[tuple(bundle)],
        volume_budget=2.0,maximum_new_elements=2,
        candidate_limit=3,beam_width=4)
    np.testing.assert_array_equal(update.selected_elements,[10,11])
    np.testing.assert_allclose(update.predicted_response,[1.0,-1.0])
    assert update.predicted_max_band_ratio==0.0
    assert update.evaluated_bundles==6
    assert "block-Schur" in update.status


def test_collaborative_search_uses_smallest_bundle_capturing_most_improvement():
    responses={
        (10,):np.array([5.2]),
        (11,):np.array([9.0]),
        (10,11):np.array([5.0]),
    }
    update=select_collaborative_element_batch(
        current_response=np.array([10.0]),response_target=np.array([0.0]),
        response_band=np.array([1.0]),candidate_elements=np.array([10,11]),
        candidate_volumes=np.ones(2),
        evaluate_bundle_response=lambda bundle:responses[tuple(bundle)],
        volume_budget=2.0,maximum_new_elements=2,
        candidate_limit=2,beam_width=2,improvement_capture=.9)
    # The pair improves the normalized error from 10 to 5, but element 10
    # alone captures 96% of that reduction with half the added volume.
    np.testing.assert_array_equal(update.selected_elements,[10])
    assert update.predicted_max_band_ratio==5.2


def test_all_candidate_aca_qr_tsvd_determines_binary_cardinality():
    current=np.array([3.0,3.0])
    target=np.zeros(2);band=np.ones(2)
    elements=np.arange(10,dtype=np.int64)+20
    # Ten candidates enter one rank-two response matrix.  Only the first two
    # are required; cardinality follows the TSVD target solve, not a fixed cap.
    delta=np.array([
        [-2.9,0.0,-1.1,-.8,-.6,-.4,-.3,-.2,-.1,-.05],
        [0.0,-2.9,-1.1,-.7,-.5,-.4,-.25,-.15,-.08,-.04]])
    update=select_tsvd_element_candidates(
        current_response=current,response_target=target,response_band=band,
        candidate_elements=elements,candidate_response_delta=delta,
        candidate_volumes=np.ones(10),volume_budget=10.0,
        relative_tolerance=1e-10,improvement_capture=.9)
    assert update.aca_rank==2 and update.numerical_rank==2
    np.testing.assert_array_equal(update.selected_elements,[20,21])
    assert update.predicted_max_band_ratio<0.11
    assert update.relative_truncation_error<1e-12
    assert update.linearized_reachable
    assert update.linearized_reachability_max_band_ratio<1e-12


def test_all_candidate_tsvd_reports_unreachable_response_component():
    update=select_tsvd_element_candidates(
        current_response=[0.0,0.0],response_target=[1.0,1.0],
        response_band=[.1,.1],candidate_elements=[20],
        candidate_response_delta=[[1.0],[0.0]],candidate_volumes=[1.0],
        volume_budget=1.0,relative_tolerance=1e-10,
        improvement_capture=1.0)
    assert update.numerical_rank==1
    assert not update.linearized_reachable
    assert update.linearized_reachability_max_band_ratio==pytest.approx(10.0)
    np.testing.assert_allclose(
        update.linearized_reachability_residual,[0.0,-1.0],atol=1e-12)


def test_tsvd_secondary_cost_breaks_physics_tie_without_filtering_columns():
    update=select_tsvd_element_candidates(
        current_response=[1.0],response_target=[0.0],response_band=[1.0],
        candidate_elements=[10,11],candidate_response_delta=[[-1.0,-1.0]],
        candidate_volumes=[1.0,1.0],volume_budget=2.0,
        maximum_changed_volume=1.0,candidate_secondary_cost=[1.0,-1.0],
        relative_tolerance=1e-10,improvement_capture=1.0)
    np.testing.assert_array_equal(update.selected_elements,[11])
    np.testing.assert_allclose(update.predicted_response,[0.0],atol=1e-12)


def test_tsvd_magnetization_sign_selects_addition_or_removal():
    elements=np.array([10,11],dtype=np.int64)
    volumes=np.ones(2)
    add=select_tsvd_element_candidates(
        current_response=[0.0],response_target=[1.0],response_band=[.1],
        candidate_elements=elements,candidate_response_delta=[[1.0,.2]],
        candidate_volumes=volumes,volume_budget=2.0,
        candidate_material_active=[False,True],relative_tolerance=1e-10,
        improvement_capture=1.0)
    np.testing.assert_array_equal(add.selected_elements,[10])
    np.testing.assert_array_equal(add.selected_directions,[1])
    remove=select_tsvd_element_candidates(
        current_response=[1.0],response_target=[0.0],response_band=[.1],
        candidate_elements=elements,candidate_response_delta=[[.2,1.0]],
        candidate_volumes=volumes,volume_budget=0.0,
        candidate_material_active=[False,True],relative_tolerance=1e-10,
        improvement_capture=1.0)
    np.testing.assert_array_equal(remove.selected_elements,[11])
    np.testing.assert_array_equal(remove.selected_directions,[-1])
    assert remove.added_volume==-1.0


def test_abe_murata_tsvd_reports_equivalent_material_volume_predictor():
    response=np.array([[2.0,0.1,0.3],[0.2,1.0,0.4]])
    target=np.array([1.9,-0.8])
    volumes=np.array([3.0,5.0,7.0])
    update=select_tsvd_element_candidates(
        current_response=[0.0,0.0],response_target=target,
        response_band=[1.0,1.0],candidate_elements=[40,41,42],
        candidate_response_delta=response,candidate_volumes=volumes,
        volume_budget=3.0,candidate_material_active=[False,True,True],
        relative_tolerance=1e-10,
        improvement_capture=1.0)
    diagnostics=update.abe_murata_diagnostics
    assert diagnostics is not None and diagnostics.retained_rank==2
    np.testing.assert_array_equal(diagnostics.candidate_elements,[40,41,42])
    np.testing.assert_array_equal(
        diagnostics.candidate_material_active,[False,True,True])
    expected_fractions=np.linalg.pinv(response,rcond=1e-10)@target
    np.testing.assert_allclose(
        diagnostics.signed_material_fractions,expected_fractions,atol=1e-12)
    np.testing.assert_allclose(
        diagnostics.equivalent_volume_changes,
        expected_fractions*volumes,atol=1e-12)
    np.testing.assert_allclose(
        diagnostics.projected_normalized_correction,[1.9,-0.8],atol=1e-12)
    reconstructed=(
        np.asarray(update.abe_murata_diagnostics.signed_material_fractions))
    np.testing.assert_allclose(
        response@reconstructed,target,atol=1e-12)
    assert diagnostics.signed_material_fractions[0]>0.0
    assert np.all(diagnostics.signed_material_fractions[1:]<0.0)
    np.testing.assert_array_equal(
        update.representative_directions,[1,-1])


def test_clustered_tsvd_removal_front_preserves_spatial_modes():
    import radia.topology_optimization as topopt
    elements=np.arange(10,16,dtype=np.int64)
    # Three spatial clusters contain locally near-parallel columns.  A global
    # response-only rank decision may keep only the first cluster; the native
    # tree front must retain at least one column from all three.
    delta=np.array([
        [1.0,.9, .05,.04, .01,.02],
        [.02,.01, 1.0,.8, .04,.03],
        [.01,.02, .03,.04, 1.0,.9]])
    labels=np.array([0,0,1,1,2,2])
    front=topopt._clustered_tsvd_candidate_front(
        candidate_elements=elements,candidate_response_delta=delta,
        response_band=np.ones(3),cluster_labels=labels,
        signed_coefficients=-np.array([3.,2.,5.,4.,7.,6.]),
        selected_elements=[11],representative_elements=[10],
        relative_tolerance=1e-6,front_limit=4)
    assert 11 in front and 10 in front
    assert len(front)<=4
    np.testing.assert_array_equal(
        np.unique(labels[np.searchsorted(elements,front)]),[0,1,2])


def test_adjoint_corrected_removal_response_uses_local_schur_action():
    import radia.topology_optimization as topopt
    A=np.array([[4.,1.,0.],[1.,3.,.5],[0.,.5,2.]])
    state=np.array([.2,-.4,.7])
    C=np.array([[1.,2.,-.5],[-.3,.1,1.2]])
    adjoints=np.linalg.solve(A,C.T)
    block=np.array([1,2],dtype=np.int32)
    class LocalBlockGram:
        @staticmethod
        def configured_linear_material_element_blocks(inv_chi,dofs,offsets):
            del inv_chi,offsets
            return A[np.ix_(dofs,dofs)].reshape(-1)
    observed=topopt._adjoint_corrected_removal_material_response(
        charge_gram=LocalBlockGram(),inv_chi=.2,dof_blocks=(block,),
        state=state,response_matrix=C,
        screen_context={"adjoints":adjoints,"adjoint_rows":np.array([0,1])},
        response_band=np.ones(2))
    expected=adjoints[block].T@(A[np.ix_(block,block)]@state[block])
    np.testing.assert_allclose(observed[:,0],expected,rtol=2e-15,atol=2e-15)


def test_native_candidate_blocks_follow_charge_hmatrix_cluster_tree():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=4,ny=2,nz=1)
    fes=ng.HDiv(mesh,order=1,discontinuous=True)
    with ng.TaskManager():
        _,gram,_=build_charge_gram(
            fes,eps=1e-9,leafsize=4,eta=1.0,
            internal_interfaces=True)
    blocks=ngsolve_discontinuous_element_dof_blocks(fes)
    packed=np.concatenate(blocks).astype(np.int32)
    offsets=np.asarray(np.r_[0,np.cumsum([len(block) for block in blocks])],
                       dtype=np.int32)
    membership=gram.configured_linear_material_candidate_clusters(
        packed,offsets,4)
    labels=np.asarray(membership["labels"],dtype=np.int64)
    assert labels.shape==(mesh.ne,)
    assert 1<int(membership["n_cluster"])<=4
    assert np.unique(labels).size>1
    assert np.all((labels>=0)&(labels<int(membership["n_cluster"])))


def test_tsvd_keeps_cumulative_terminal_depths_mutually_exclusive():
    update=select_tsvd_element_candidates(
        current_response=[2.0],response_target=[0.0],response_band=[1.0],
        candidate_elements=[10,11],candidate_response_delta=[[1.1,1.2]],
        candidate_volumes=[1.0,2.0],volume_budget=0.0,
        candidate_material_active=[True,True],
        candidate_volume_changes=[1.0,2.0],
        candidate_exclusion_groups=[3,3],relative_tolerance=1e-10,
        improvement_capture=1.0)
    assert len(update.selected_elements)==1
    np.testing.assert_array_equal(update.selected_directions,[-1])


def test_boundary_removal_is_outermost_and_preserves_seed():
    from ngsolve.meshes import MakeStructured3DMesh
    mesh=MakeStructured3DMesh(hexes=True,nx=3,ny=1,nz=1)
    active=np.ones(mesh.ne,dtype=bool)
    fixed=np.zeros(mesh.ne,dtype=bool);fixed[0]=True
    predecessor=np.array([-1,0,1],dtype=np.int64)
    np.testing.assert_array_equal(ngsolve_boundary_removal_candidates(
        mesh,active,fixed_active_elements=fixed,
        predecessor_elements=predecessor),[2])
    active[2]=False
    np.testing.assert_array_equal(ngsolve_boundary_removal_candidates(
        mesh,active,fixed_active_elements=fixed,
        predecessor_elements=predecessor),[1])


def test_hdiv_mmm_generation_driver_accepts_purely_collaborative_pair():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=3,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(fes,eps=1e-10,leafsize=256,
            eta=2.0,internal_interfaces=True)
    rng=np.random.default_rng(4471)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    zero_response=np.zeros((1,fes.ndof))
    masks=[]
    for selected in ((1,),(0,1),(1,2),(0,1,2)):
        mask=np.zeros(mesh.ne,dtype=bool);mask[list(selected)]=True;masks.append(mask)
    states=[solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=zero_response,active_elements=mask,
        solve_tolerance=1e-11)[0] for mask in masks]
    state_deltas=np.vstack((states[1]-states[0],states[2]-states[0],
                            states[3]-states[0]))
    response_row=np.linalg.lstsq(
        state_deltas,np.array([0.0,0.0,1.0]),rcond=None)[0][None,:]
    realized=state_deltas@response_row[0]
    np.testing.assert_allclose(realized,[0.0,0.0,1.0],atol=2e-11)
    current=float((response_row@states[0]).item())
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))
    result=grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_row,active_elements=masks[0],
        element_volumes=volumes,response_target=[current+1.0],
        response_band=[1e-8],volume_max=float(np.sum(volumes))+1e-14,
        maximum_batch_elements=2,max_iterations=1,
        solve_tolerance=1e-11)
    assert result.converged and len(result.history)==1
    assert result.stop_reason=="target_met"
    np.testing.assert_array_equal(result.history[0].added_elements,[0,2])
    assert result.history[0].selection_model==\
        "all-candidate-aca-qr-tsvd-exact-conditional"
    assert result.history[0].collaborative_bundles_evaluated==3
    assert result.history[0].superposed_max_band_ratio>=9e7
    np.testing.assert_allclose(result.response,[current+1.0],atol=3e-11)


def test_hdiv_mmm_graph_front_finds_connected_tsvd_seed_collaboration(
        monkeypatch):
    import ngsolve as ng
    import radia.topology_optimization as topopt
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram

    mesh=MakeStructured3DMesh(hexes=True,nx=3,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(
            fes,eps=1e-10,leafsize=256,eta=2.0,
            internal_interfaces=True)
    rng=np.random.default_rng(20260811)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    zero_response=np.zeros((1,fes.ndof))
    masks=[]
    for selected in ((1,),(0,1),(1,2),(0,1,2)):
        mask=np.zeros(mesh.ne,dtype=bool);mask[list(selected)]=True
        masks.append(mask)
    states=[topopt.solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=zero_response,active_elements=mask,
        solve_tolerance=1e-11)[0] for mask in masks]
    state_deltas=np.vstack((states[1]-states[0],states[2]-states[0],
                            states[3]-states[0]))
    response_row=np.linalg.lstsq(
        state_deltas,np.array([0.0,0.0,1.0]),rcond=None)[0][None,:]
    current=float((response_row@states[0]).item())

    def force_empty_global_but_keep_qr_seeds(**kwargs):
        elements=np.asarray(kwargs["candidate_elements"],dtype=np.int64)
        material=np.asarray(kwargs["candidate_material_active"],dtype=bool)
        seeds=np.sort(elements[~material])
        coefficients=np.ones(len(elements))
        return topopt.TSVDElementCandidateSelection(
            selected_elements=np.empty(0,dtype=np.int64),
            selected_directions=np.empty(0,dtype=np.int8),
            representative_elements=seeds,
            representative_directions=np.ones(len(seeds),dtype=np.int8),
            predicted_response=np.asarray(kwargs["current_response"],dtype=float),
            predicted_max_band_ratio=float(np.max(np.abs((
                np.asarray(kwargs["current_response"])-
                np.asarray(kwargs["response_target"]))/
                np.asarray(kwargs["response_band"])))),
            added_volume=0.0,numerical_rank=1,aca_rank=1,
            singular_values=np.ones(1),signed_coefficients=coefficients,
            relative_truncation_error=0.0,
            status="forced empty global proposal with QR seeds")

    monkeypatch.setattr(topopt,"select_tsvd_element_candidates",
                        force_empty_global_but_keep_qr_seeds)
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))
    result=topopt.grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_row,active_elements=masks[0],
        element_volumes=volumes,response_target=[current+1.0],
        response_band=[1e-8],volume_max=float(np.sum(volumes))+1e-14,
        maximum_batch_elements=2,max_iterations=1,solve_tolerance=1e-11,
        graph_front_maximum_components=2,
        graph_front_proposal_limit=8,
        graph_front_response_novelty_weight=0.5)
    assert result.converged and len(result.history)==1
    np.testing.assert_array_equal(result.active_elements,masks[-1])
    assert result.history[0].selection_model==(
        "aca-qr-tsvd-connected-graph-front-full-resolve")
    assert result.history[0].graph_front_proposals_evaluated>=3
    assert result.graph_front_diagnostics
    diagnostics=result.graph_front_diagnostics[0]
    assert diagnostics.novelty_weight==0.5
    assert diagnostics.pool_proposal_count>=3
    assert diagnostics.selected_response_rank==1


def test_hdiv_mmm_exact_beam_crosses_one_worsening_lego_state(monkeypatch):
    import ngsolve as ng
    import radia.topology_optimization as topopt
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram

    mesh=MakeStructured3DMesh(hexes=True,nx=3,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(
            fes,eps=1e-10,leafsize=256,eta=2.0,
            internal_interfaces=True)
    rng=np.random.default_rng(20260812)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    zero_response=np.zeros((1,fes.ndof))
    masks=[]
    for selected in ((0,),(0,1),(0,1,2)):
        mask=np.zeros(mesh.ne,dtype=bool);mask[list(selected)]=True
        masks.append(mask)
    states=[topopt.solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=zero_response,active_elements=mask,
        solve_tolerance=1e-11)[0] for mask in masks]
    state_deltas=np.vstack((states[1]-states[0],states[2]-states[0]))
    response_row=np.linalg.lstsq(
        state_deltas,np.array([-0.2,1.0]),rcond=None)[0][None,:]
    realized=state_deltas@response_row[0]
    np.testing.assert_allclose(realized,[-0.2,1.0],atol=3e-11)
    current=float((response_row@states[0]).item())

    def force_next_addition(**kwargs):
        elements=np.asarray(kwargs["candidate_elements"],dtype=np.int64)
        material=np.asarray(kwargs["candidate_material_active"],dtype=bool)
        additions=np.sort(elements[~material])
        assert additions.size
        wanted=int(additions[0])
        coefficients=np.zeros(len(elements))
        coefficients[np.flatnonzero(elements==wanted)[0]]=1.0
        return topopt.TSVDElementCandidateSelection(
            selected_elements=np.array([wanted],dtype=np.int64),
            selected_directions=np.array([1],dtype=np.int8),
            representative_elements=np.array([wanted],dtype=np.int64),
            representative_directions=np.array([1],dtype=np.int8),
            predicted_response=np.asarray(kwargs["response_target"],dtype=float),
            predicted_max_band_ratio=0.0,
            added_volume=float(np.asarray(kwargs["candidate_volumes"])[
                np.flatnonzero(elements==wanted)[0]]),
            numerical_rank=1,aca_rank=1,singular_values=np.ones(1),
            signed_coefficients=coefficients,relative_truncation_error=0.0,
            status="forced next outward addition")

    monkeypatch.setattr(topopt,"select_tsvd_element_candidates",
                        force_next_addition)
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))
    result=topopt.grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_row,active_elements=masks[0],
        element_volumes=volumes,response_target=[current+1.0],
        response_band=[1e-8],volume_max=float(np.sum(volumes))+1e-14,
        fixed_active_elements=masks[0],
        predecessor_elements=np.array([-1,0,1],dtype=np.int64),
        max_iterations=1,solve_tolerance=1e-11,
        exact_beam_width=2,exact_beam_depth=2,
        exact_beam_barrier_fraction=0.25)
    assert result.converged and len(result.history)==1
    np.testing.assert_array_equal(result.active_elements,masks[-1])
    np.testing.assert_array_equal(result.history[0].added_elements,[1,2])
    assert result.history[0].removed_elements.size==0
    assert result.history[0].nonmonotone_search_depth==1
    np.testing.assert_allclose(result.response,[current+1.0],atol=3e-11)
    assert [trial.depth for trial in result.exact_search_trace]==[1,2]
    assert result.exact_search_trace[0].parent_max_band_ratio==(
        result.exact_search_trace[0].incumbent_max_band_ratio)
    assert result.exact_search_trace[0].max_band_ratio>(
        result.exact_search_trace[0].parent_max_band_ratio)
    np.testing.assert_array_equal(
        result.exact_search_trace[0].added_elements,[1])
    np.testing.assert_array_equal(
        result.exact_search_trace[1].added_elements,[1,2])
    assert all(trial.solve_iterations>0
               for trial in result.exact_search_trace)

    blocked=topopt.grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_row,active_elements=masks[0],
        element_volumes=volumes,response_target=[current+1.0],
        response_band=[1e-8],volume_max=float(np.sum(volumes))+1e-14,
        fixed_active_elements=masks[0],
        predecessor_elements=np.array([-1,0,1],dtype=np.int64),
        max_iterations=1,solve_tolerance=1e-11,
        exact_beam_width=2,exact_beam_depth=2,
        exact_beam_barrier_fraction=0.1)
    assert not blocked.converged and len(blocked.history)==0
    np.testing.assert_array_equal(blocked.active_elements,masks[0])
    assert blocked.stop_reason=="exact_nonmonotone_beam_exhausted"
    assert [trial.depth for trial in blocked.exact_search_trace]==[1]

    guard_calls=[]
    guarded=topopt.grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_row,active_elements=masks[0],
        element_volumes=volumes,response_target=[current+1.0],
        response_band=[1e-8],volume_max=float(np.sum(volumes))+1e-14,
        fixed_active_elements=masks[0],
        predecessor_elements=np.array([-1,0,1],dtype=np.int64),
        max_iterations=1,solve_tolerance=1e-11,
        exact_beam_width=2,exact_beam_depth=2,
        exact_beam_barrier_fraction=0.25,
        exact_response_validator=lambda raw,objective: (
            guard_calls.append((raw.copy(),objective.copy())) or
            bool(raw[0]>=current-1e-10)))
    assert not guarded.converged and len(guarded.history)==0
    assert guarded.stop_reason=="exact_nonmonotone_beam_exhausted"
    assert guarded.exact_search_trace==()
    assert len(guard_calls)>=2
    np.testing.assert_allclose(guard_calls[0][0],[current],atol=3e-11)
    np.testing.assert_allclose(guard_calls[0][1],[current],atol=3e-11)

    depth_one=topopt.grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_row,active_elements=masks[0],
        element_volumes=volumes,response_target=[current+1.0],
        response_band=[1e-8],volume_max=float(np.sum(volumes))+1e-14,
        fixed_active_elements=masks[0],
        predecessor_elements=np.array([-1,0,1],dtype=np.int64),
        max_iterations=1,solve_tolerance=1e-11,
        exact_beam_width=2,exact_beam_depth=1,
        exact_beam_barrier_fraction=0.25)
    assert not depth_one.converged and len(depth_one.history)==0
    assert depth_one.stop_reason=="exact_nonmonotone_beam_exhausted"
    assert [trial.depth for trial in depth_one.exact_search_trace]==[1]


def test_hdiv_mmm_generation_checks_global_addition_proposal_before_dense_schur(
        monkeypatch):
    import ngsolve as ng
    import radia.topology_optimization as topopt
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram

    mesh=MakeStructured3DMesh(hexes=True,nx=3,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(
            fes,eps=1e-10,leafsize=256,eta=2.0,
            internal_interfaces=True)
    active=np.array([False,True,False])
    target_active=np.array([True,True,False])
    rng=np.random.default_rng(20260806)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(2,fes.ndof))
    _,target,_=solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=target_active,
        solve_tolerance=1e-11)

    def force_global_addition(**kwargs):
        elements=np.asarray(kwargs["candidate_elements"],dtype=np.int64)
        material=np.asarray(kwargs["candidate_material_active"],dtype=bool)
        wanted=int(elements[(elements==0)&~material][0])
        coefficients=np.zeros(len(elements))
        coefficients[np.flatnonzero(elements==wanted)[0]]=1.0
        return topopt.TSVDElementCandidateSelection(
            selected_elements=np.array([wanted],dtype=np.int64),
            selected_directions=np.array([1],dtype=np.int8),
            representative_elements=np.array([wanted],dtype=np.int64),
            representative_directions=np.array([1],dtype=np.int8),
            predicted_response=np.asarray(kwargs["response_target"],dtype=float),
            predicted_max_band_ratio=0.0,added_volume=1.0,
            numerical_rank=1,aca_rank=1,singular_values=np.ones(1),
            signed_coefficients=coefficients,relative_truncation_error=0.0,
            status="forced global addition")

    monkeypatch.setattr(
        topopt,"select_tsvd_element_candidates",force_global_addition)
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))
    result=topopt.grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        element_volumes=volumes,response_target=target,
        response_band=np.full(2,1e-8),volume_max=float(np.sum(volumes)),
        fixed_active_elements=np.array([False,True,False]),
        maximum_batch_elements=1,max_iterations=1,solve_tolerance=1e-11)
    assert result.converged and len(result.history)==1
    np.testing.assert_array_equal(result.active_elements,target_active)
    np.testing.assert_array_equal(result.history[0].added_elements,[0])
    assert result.history[0].selection_model==(
        "all-candidate-aca-qr-tsvd-direct-full-resolve")
    assert result.history[0].collaborative_bundles_evaluated==1
    assert result.history[0].candidate_coupling_rank==0
    assert result.history[0].native_reduction_timings["solve_s"]==0.0


def test_hdiv_mmm_generation_uses_qr_representative_adjoint_rows(monkeypatch):
    import ngsolve as ng
    import radia.topology_optimization as topopt
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram

    mesh=MakeStructured3DMesh(hexes=True,nx=4,ny=4,nz=4)
    centers=_element_centroids(mesh)
    active=np.all((centers>.25)&(centers<.75),axis=1)
    candidates=ngsolve_boundary_growth_candidates(mesh,active)
    assert candidates.size>8
    wanted=int(candidates[0])
    target_active=active.copy();target_active[wanted]=True
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(
            fes,eps=1e-10,leafsize=256,eta=2.0,
            internal_interfaces=True)
    rng=np.random.default_rng(20260808)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(4,fes.ndof))
    _,target,_=solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=target_active,
        solve_tolerance=1e-11)
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))

    def force_wanted(**kwargs):
        elements=np.asarray(kwargs["candidate_elements"],dtype=np.int64)
        coefficients=np.zeros(len(elements))
        coefficients[np.flatnonzero(elements==wanted)[0]]=1.0
        return topopt.TSVDElementCandidateSelection(
            selected_elements=np.array([wanted],dtype=np.int64),
            selected_directions=np.array([1],dtype=np.int8),
            representative_elements=np.array([wanted],dtype=np.int64),
            representative_directions=np.array([1],dtype=np.int8),
            predicted_response=np.asarray(kwargs["response_target"],dtype=float),
            predicted_max_band_ratio=0.0,added_volume=float(volumes[wanted]),
            numerical_rank=1,aca_rank=1,singular_values=np.ones(1),
            signed_coefficients=coefficients,relative_truncation_error=0.0,
            status="forced QR-row proposal")

    monkeypatch.setattr(topopt,"select_tsvd_element_candidates",force_wanted)
    result=topopt.grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        element_volumes=volumes,response_target=target,
        response_band=np.full(4,1e-8),
        volume_max=float(volumes@active+volumes[wanted]),
        fixed_active_elements=active,max_iterations=1,
        solve_tolerance=1e-11,proposal_adjoint_count=4)
    assert result.converged and len(result.history)==1
    np.testing.assert_array_equal(result.active_elements,target_active)
    assert result.history[0].response_adjoint_count==4
    assert result.history[0].screened_candidate_count==0


def test_hdiv_mmm_generation_shrinks_all_candidate_volume_trust_region(
        monkeypatch):
    import ngsolve as ng
    import radia.topology_optimization as topopt
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram

    mesh=MakeStructured3DMesh(hexes=True,nx=3,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(
            fes,eps=1e-10,leafsize=256,eta=2.0,
            internal_interfaces=True)
    active=np.array([False,True,False])
    target_active=np.array([True,True,False])
    rng=np.random.default_rng(20260807)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(2,fes.ndof))
    _,target,_=solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=target_active,
        solve_tolerance=1e-11)
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))
    first_budget=float(volumes[0]+volumes[2])

    def force_budget_dependent_batch(**kwargs):
        elements=np.asarray(kwargs["candidate_elements"],dtype=np.int64)
        material=np.asarray(kwargs["candidate_material_active"],dtype=bool)
        additions=elements[~material]
        wanted=int(additions[additions==0][0])
        trust=kwargs.get("maximum_changed_volume")
        trust=first_budget if trust is None else float(trust)
        selected=(additions if trust>
                  .75*first_budget else np.array([wanted],dtype=np.int64))
        coefficients=np.zeros(len(elements))
        coefficients[np.isin(elements,selected)]=1.0
        return topopt.TSVDElementCandidateSelection(
            selected_elements=np.asarray(selected,dtype=np.int64),
            selected_directions=np.ones(len(selected),dtype=np.int8),
            representative_elements=np.asarray(selected,dtype=np.int64),
            representative_directions=np.ones(len(selected),dtype=np.int8),
            predicted_response=np.asarray(kwargs["response_target"],dtype=float),
            predicted_max_band_ratio=0.0,
            added_volume=float(np.sum(volumes[selected])),
            numerical_rank=1,aca_rank=1,singular_values=np.ones(1),
            signed_coefficients=coefficients,relative_truncation_error=0.0,
            status="forced volume-trust proposal")

    monkeypatch.setattr(
        topopt,"select_tsvd_element_candidates",force_budget_dependent_batch)
    result=topopt.grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        element_volumes=volumes,response_target=target,
        response_band=np.full(2,1e-8),volume_max=float(np.sum(volumes)),
        fixed_active_elements=np.array([False,True,False]),
        max_iterations=1,solve_tolerance=1e-11,
        minimum_model_agreement=.999,proposal_trust_region_trials=3)
    assert result.converged and len(result.history)==1
    np.testing.assert_array_equal(result.active_elements,target_active)
    np.testing.assert_array_equal(result.history[0].added_elements,[0])
    assert result.history[0].selection_model==(
        "all-candidate-aca-qr-tsvd-adaptive-trust-full-resolve")
    assert result.history[0].collaborative_bundles_evaluated==2


def test_hdiv_mmm_generation_driver_commits_exact_whole_element_batch():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram
    mesh=MakeStructured3DMesh(hexes=False,nx=2,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(fes,eps=1e-10,leafsize=256,
            eta=2.0,internal_interfaces=True)
    rng=np.random.default_rng(5521)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(2,fes.ndof))
    active=np.zeros(mesh.ne,dtype=bool);active[0]=True
    initial=linearize_hdiv_mmm_element_generation(charge_gram=gram,fes=fes,
        inv_chi=.2,rhs=rhs,response_matrix=response_matrix,
        active_elements=active,solve_tolerance=1e-11)
    wanted=int(initial.candidate_elements[0])
    target=initial.response+initial.candidate_response_delta[:,0]
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))
    result=grow_hdiv_mmm_by_superposition(charge_gram=gram,fes=fes,
        inv_chi=.2,rhs=rhs,response_matrix=response_matrix,
        active_elements=active,element_volumes=volumes,
        response_target=target,response_band=np.full(2,1e-8),
        volume_max=volumes[0]+volumes[wanted]+1e-14,
        maximum_batch_elements=1,max_iterations=2,solve_tolerance=1e-11,
        graph_interface_weight=0.02)
    assert result.converged and len(result.history)==1
    assert result.stop_reason=="target_met"
    assert "aca-qr-tsvd" in result.history[0].selection_model
    assert result.history[0].added_elements.tolist()==[wanted]
    assert result.active_elements[wanted] and np.count_nonzero(result.active_elements)==2
    np.testing.assert_allclose(result.response,target,rtol=0,atol=3e-12)


def test_hdiv_mmm_generation_removes_negative_magnetization_candidate():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=3,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(fes,eps=1e-10,leafsize=256,
            eta=2.0,internal_interfaces=True)
    rng=np.random.default_rng(20260730)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(2,fes.ndof))
    current=np.array([True,True,False])
    target_active=np.array([True,False,False])
    _,target,_=solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=target_active,
        solve_tolerance=1e-11)
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))
    result=grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=current,
        element_volumes=volumes,response_target=target,
        response_band=np.full(2,1e-8),volume_max=float(np.sum(volumes)),
        fixed_active_elements=np.array([True,False,False]),
        predecessor_elements=np.array([-1,0,1]),max_iterations=1,
        solve_tolerance=1e-11)
    assert result.converged and len(result.history)==1
    np.testing.assert_array_equal(result.history[0].added_elements,[])
    np.testing.assert_array_equal(result.history[0].removed_elements,[1])
    assert result.history[0].selection_model==\
        "all-candidate-aca-qr-tsvd-direct-full-resolve"
    assert result.history[0].candidate_coupling_rank==0
    assert result.history[0].native_reduction_timings["solve_s"]==0.0
    diagnostics=result.history[0].abe_murata_diagnostics
    assert diagnostics is not None
    assert diagnostics.retained_rank==result.history[0].tsvd_rank
    assert diagnostics.candidate_elements.size==(
        result.history[0].addition_candidate_count+
        result.history[0].removal_candidate_count)
    assert np.any(diagnostics.equivalent_volume_changes<0.0)
    np.testing.assert_array_equal(result.active_elements,target_active)
    np.testing.assert_allclose(result.response,target,rtol=0,atol=4e-11)


def test_hdiv_mmm_projected_adjoint_removal_is_not_double_transformed(
        monkeypatch):
    import ngsolve as ng
    import radia.topology_optimization as topopt
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=3,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(
            fes,eps=1e-10,leafsize=256,eta=2.0,
            internal_interfaces=True)
    rng=np.random.default_rng(20260809)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    raw_rows=rng.normal(size=(3,fes.ndof))
    projection=np.array([[1.,-.2,.3],[.1,.7,-.4]])
    current=np.array([True,True,False])
    target_active=np.array([True,False,False])
    _,target_raw,_=solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=raw_rows,active_elements=target_active,
        solve_tolerance=1e-11)
    target=projection@target_raw
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))
    def force_removal(**kwargs):
        elements=np.asarray(kwargs["candidate_elements"],dtype=np.int64)
        material=np.asarray(kwargs["candidate_material_active"],dtype=bool)
        wanted=int(elements[(elements==1)&material][0])
        coefficients=np.zeros(len(elements))
        coefficients[np.flatnonzero(elements==wanted)[0]]=-1.0
        return topopt.TSVDElementCandidateSelection(
            selected_elements=np.array([wanted],dtype=np.int64),
            selected_directions=np.array([-1],dtype=np.int8),
            representative_elements=np.array([wanted],dtype=np.int64),
            representative_directions=np.array([-1],dtype=np.int8),
            predicted_response=np.asarray(kwargs["response_target"],dtype=float),
            predicted_max_band_ratio=0.0,added_volume=-float(volumes[wanted]),
            numerical_rank=1,aca_rank=1,singular_values=np.ones(1),
            signed_coefficients=coefficients,relative_truncation_error=0.0,
            status="forced projected deletion")
    monkeypatch.setattr(topopt,"select_tsvd_element_candidates",force_removal)
    result=topopt.grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=raw_rows,active_elements=current,
        element_volumes=volumes,response_target=target,
        response_band=np.full(2,1e-8),volume_max=float(np.sum(volumes)),
        fixed_active_elements=np.array([True,False,False]),
        predecessor_elements=np.array([-1,0,1]),max_iterations=1,
        solve_tolerance=1e-11,
        response_transform=lambda values:projection@values,
        response_transform_jacobian=lambda values:projection)
    assert result.converged and len(result.history)==1
    np.testing.assert_array_equal(result.active_elements,target_active)
    np.testing.assert_array_equal(result.history[0].removed_elements,[1])
    np.testing.assert_allclose(result.objective_response,target,
                               rtol=0,atol=5e-11)


def test_hdiv_mmm_generation_exactly_checks_alternate_removal_after_bad_tsvd(
        monkeypatch):
    import ngsolve as ng
    import radia.topology_optimization as topopt
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram

    mesh=MakeStructured3DMesh(hexes=True,nx=2,ny=2,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(fes,eps=1e-10,leafsize=256,
            eta=2.0,internal_interfaces=True)
    active=np.array([True,True,True,False])
    fixed_active=np.array([True,False,False,False])
    good_removal=1
    bad_removal=2
    target_active=active.copy();target_active[good_removal]=False
    rng=np.random.default_rng(20260805)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(3,fes.ndof))
    _,target,_=solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=target_active,
        solve_tolerance=1e-11)

    def deliberately_bad_signed_proposal(**kwargs):
        elements=np.asarray(kwargs["candidate_elements"],dtype=np.int64)
        material=np.asarray(kwargs["candidate_material_active"],dtype=bool)
        assert set(elements[material])=={good_removal,bad_removal}
        coefficients=np.zeros(len(elements))
        coefficients[np.flatnonzero(elements==bad_removal)[0]]=-1.0
        return topopt.TSVDElementCandidateSelection(
            selected_elements=np.array([bad_removal],dtype=np.int64),
            selected_directions=np.array([-1],dtype=np.int8),
            representative_elements=np.array([good_removal],dtype=np.int64),
            representative_directions=np.array([-1],dtype=np.int8),
            predicted_response=np.asarray(kwargs["response_target"],dtype=float),
            predicted_max_band_ratio=0.0,added_volume=-1.0,
            numerical_rank=2,aca_rank=2,singular_values=np.ones(2),
            signed_coefficients=coefficients,relative_truncation_error=0.0,
            status="deliberately bad deletion proposal")

    monkeypatch.setattr(
        topopt,"select_tsvd_element_candidates",deliberately_bad_signed_proposal)
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))
    result=topopt.grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        element_volumes=volumes,response_target=target,
        response_band=np.full(3,1e-8),volume_max=float(np.sum(volumes)),
        fixed_active_elements=fixed_active,max_iterations=1,
        solve_tolerance=1e-11,graph_front_proposal_limit=0)
    assert result.converged and len(result.history)==1
    np.testing.assert_array_equal(result.active_elements,target_active)
    np.testing.assert_array_equal(result.history[0].added_elements,[])
    np.testing.assert_array_equal(result.history[0].removed_elements,[good_removal])
    assert result.history[0].selection_model==(
        "signed-magnetization-aca-qr-tsvd-alternate-removal-exact")
    assert result.history[0].collaborative_bundles_evaluated==2


def test_hdiv_mmm_generation_supports_removal_only_front():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=2,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(fes,eps=1e-10,leafsize=256,
            eta=2.0,internal_interfaces=True)
    rng=np.random.default_rng(20260731)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(2,fes.ndof))
    current=np.ones(2,dtype=bool);target_active=np.array([True,False])
    _,target,_=solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=target_active,
        solve_tolerance=1e-11)
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))
    result=grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=current,
        element_volumes=volumes,response_target=target,
        response_band=np.full(2,1e-8),volume_max=float(np.sum(volumes)),
        fixed_active_elements=np.array([True,False]),
        predecessor_elements=np.array([-1,0]),max_iterations=1,
        solve_tolerance=1e-11)
    assert result.converged and len(result.history)==1
    np.testing.assert_array_equal(result.history[0].removed_elements,[1])
    assert result.history[0].addition_candidate_count==0
    assert result.history[0].removal_candidate_count==1
    assert result.history[0].selection_model==\
        "signed-magnetization-aca-qr-tsvd-conditional-exact"
    assert result.history[0].collaborative_bundles_evaluated>=1
    np.testing.assert_array_equal(result.active_elements,target_active)


def test_hdiv_mmm_generation_removes_through_thickness_group_as_one_move():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=2,ny=1,nz=2)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(fes,eps=1e-10,leafsize=256,
            eta=2.0,internal_interfaces=True)
    centers=_element_centroids(mesh)
    target_group=np.flatnonzero(centers[:,0]>.5)
    assert target_group.size==2
    current=np.ones(mesh.ne,dtype=bool)
    target_active=current.copy();target_active[target_group]=False
    fixed_active=~np.isin(np.arange(mesh.ne),target_group)
    coupling=np.full(mesh.ne,-1,dtype=np.int64);coupling[target_group]=0
    rng=np.random.default_rng(20260804)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(2,fes.ndof))
    _,target,_=solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=target_active,
        solve_tolerance=1e-11)
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))
    result=grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=current,
        element_volumes=volumes,response_target=target,
        response_band=np.full(2,1e-8),volume_max=float(np.sum(volumes)),
        fixed_active_elements=fixed_active,
        removal_coupling_groups=coupling,max_iterations=1,
        solve_tolerance=1e-11)
    assert result.converged and len(result.history)==1
    np.testing.assert_array_equal(
        np.sort(result.history[0].removed_elements),target_group)
    np.testing.assert_array_equal(result.active_elements,target_active)


def test_hdiv_mmm_generation_removes_complete_coupled_predecessor_chain():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram

    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=3)
    centers=_element_centroids(mesh)
    order=np.argsort(centers[:,2])
    outer,middle,root=map(int,order)
    predecessor=np.full(mesh.ne,-1,dtype=np.int64)
    predecessor[outer]=middle;predecessor[middle]=root
    active=np.ones(mesh.ne,dtype=bool)
    target_active=active.copy();target_active[[outer,middle]]=False
    fixed_active=np.zeros(mesh.ne,dtype=bool);fixed_active[root]=True
    coupling=np.full(mesh.ne,-1,dtype=np.int64)
    coupling[[outer,middle]]=0
    raw=ngsolve_boundary_removal_candidates(
        mesh,active,fixed_active_elements=fixed_active,
        predecessor_elements=predecessor)
    np.testing.assert_array_equal(raw,[outer])
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(
            fes,eps=1e-10,leafsize=256,eta=2.0,
            internal_interfaces=True)
    rng=np.random.default_rng(20260809)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(2,fes.ndof))
    _,target,_=solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=target_active,
        solve_tolerance=1e-11)
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))
    result=grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        element_volumes=volumes,response_target=target,
        response_band=np.full(2,1e-8),volume_max=float(np.sum(volumes)),
        fixed_active_elements=fixed_active,
        predecessor_elements=predecessor,
        removal_coupling_groups=coupling,max_iterations=1,
        solve_tolerance=1e-11)
    assert result.converged and len(result.history)==1
    np.testing.assert_array_equal(
        np.sort(result.history[0].removed_elements),
        np.sort([outer,middle]))
    np.testing.assert_array_equal(result.active_elements,target_active)


def test_hdiv_mmm_generation_rejects_group_that_would_remove_all_iron():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram

    mesh=MakeStructured3DMesh(hexes=True,nx=2,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(
            fes,eps=1e-10,leafsize=256,eta=2.0,
            internal_interfaces=True)
    active=np.ones(mesh.ne,dtype=bool)
    rhs=np.asarray(mass@np.ones(fes.ndof))
    response_matrix=np.ones((1,fes.ndof))
    _,response,_=solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        solve_tolerance=1e-11)
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))
    result=grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        element_volumes=volumes,response_target=response+1.0,
        response_band=[1e-8],volume_max=float(np.sum(volumes)),
        removal_coupling_groups=np.zeros(mesh.ne,dtype=np.int64),
        max_iterations=1,solve_tolerance=1e-11)
    assert not result.converged
    assert result.stop_reason=="no_growth_candidates"
    assert len(result.history)==0


def test_hdiv_mmm_generation_recalibrates_linear_coil_source_after_batch():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram
    mesh=MakeStructured3DMesh(hexes=False,nx=2,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(fes,eps=1e-10,leafsize=256,
            eta=2.0,internal_interfaces=True)
    rng=np.random.default_rng(8015)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    response_matrix=rng.normal(size=(2,fes.ndof))
    active=np.zeros(mesh.ne,dtype=bool);active[0]=True
    initial=linearize_hdiv_mmm_element_generation(charge_gram=gram,fes=fes,
        inv_chi=.2,rhs=rhs,response_matrix=response_matrix,
        active_elements=active,solve_tolerance=1e-11)
    inserted=initial.response+initial.candidate_response_delta[:,0]
    raw_target=1.7*inserted
    transform=lambda values:np.array([values[0],values[1]**2])
    target=transform(raw_target)
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))
    result=grow_hdiv_mmm_by_superposition(charge_gram=gram,fes=fes,
        inv_chi=.2,rhs=rhs,response_matrix=response_matrix,
        active_elements=active,element_volumes=volumes,
        response_target=target,response_band=np.full(2,1e-8),
        volume_max=float(np.sum(volumes))+1e-14,
        maximum_batch_elements=1,max_iterations=2,solve_tolerance=1e-11,
        source_calibration_rows=[0],
        source_calibration_target=[raw_target[0]],
        response_transform=transform)
    assert result.converged and len(result.history)==1
    assert result.stop_reason=="target_met"
    np.testing.assert_allclose(result.response,raw_target,rtol=0,atol=4e-12)
    np.testing.assert_allclose(result.objective_response,target,rtol=0,atol=4e-12)
    np.testing.assert_allclose(result.source_scale,1.7,rtol=0,atol=2e-11)
    np.testing.assert_allclose(result.history[0].source_scale,1.7,
                               rtol=0,atol=2e-11)


def test_positive_minimax_source_scale_has_piecewise_analytic_gradient():
    from radia.topology_optimization import (
        _positive_minimax_source_scale_and_gradient,
    )

    response=np.array([1.0,1.0,0.4])
    target=np.array([1.0,3.0,0.6])
    band=np.array([1.0,0.1,2.0])
    scale,gradient=_positive_minimax_source_scale_and_gradient(
        response,target,band)
    np.testing.assert_allclose(scale,31.0/11.0,rtol=0,atol=2e-14)
    ratios=np.abs((scale*response-target)/band)
    np.testing.assert_allclose(ratios[:2],[20.0/11.0]*2,
                               rtol=0,atol=2e-14)
    direction=np.array([0.3,-0.2,0.1])
    step=2e-7
    upper,_=_positive_minimax_source_scale_and_gradient(
        response+step*direction,target,band)
    lower,_=_positive_minimax_source_scale_and_gradient(
        response-step*direction,target,band)
    regression=(upper-lower)/(2.0*step)
    np.testing.assert_allclose(
        gradient@direction,regression,rtol=2e-7,atol=2e-9)


def test_hdiv_mmm_generation_contracts_raw_rows_before_metric_adjoints():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import build_charge_gram
    mesh=MakeStructured3DMesh(hexes=False,nx=2,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=0,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(fes,eps=1e-10,leafsize=256,
            eta=2.0,internal_interfaces=True)
    rng=np.random.default_rng(20260802)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    raw_rows=rng.normal(size=(3,fes.ndof))
    active=np.zeros(mesh.ne,dtype=bool);active[0]=True
    initial=linearize_hdiv_mmm_element_generation(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=raw_rows,active_elements=active,
        solve_tolerance=1e-11)
    target_raw=initial.response+initial.candidate_response_delta[:,0]
    projection=np.array([[1.0,-.25,.5]])
    target=projection@target_raw
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))
    result=grow_hdiv_mmm_by_superposition(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=raw_rows,active_elements=active,
        element_volumes=volumes,response_target=target,
        response_band=[1e-8],volume_max=float(np.sum(volumes))+1e-14,
        maximum_batch_elements=1,max_iterations=1,solve_tolerance=1e-11,
        response_transform=lambda values:projection@values,
        response_transform_jacobian=lambda values:projection)
    assert result.converged and len(result.history)==1
    assert result.objective_response.shape==(1,)
    np.testing.assert_allclose(result.response,target_raw,rtol=0,atol=3e-11)
    np.testing.assert_allclose(result.objective_response,target,rtol=0,atol=3e-11)


def test_ngsolve_boundary_growth_candidates_only_share_complete_facets():
    from ngsolve.meshes import MakeStructured3DMesh
    mesh=MakeStructured3DMesh(hexes=True,nx=3,ny=1,nz=1)
    active=np.array([True,False,False])
    candidates=ngsolve_boundary_growth_candidates(mesh,active)
    assert candidates.tolist()==[1]
    active[:2]=True
    assert ngsolve_boundary_growth_candidates(mesh,active).tolist()==[2]


def test_shape_lp_restores_field_band_with_real_boundary_parameters():
    linearization=ShapeLinearization(objective=0.0,
        objective_gradient=np.array([-1.0,0.2]),response=np.array([2.0]),
        response_jacobian=np.array([[1.0,0.0]]),response_target=np.array([0.0]),
        response_band=np.array([1.0]))
    update=solve_shape_lp(np.zeros(2),linearization,move_limit=[.1,.2],
        parameter_bounds=([-1.0,-1.0],[1.0,1.0]))
    # Feasibility restoration wins over the objective, which alone would move q0 positive.
    np.testing.assert_allclose(update.parameters,[-.1,-.2],atol=1e-10)
    assert update.restoration and abs(update.predicted_max_band_ratio-1.9)<3e-9


def test_shape_lp_keeps_feasible_response_and_obeys_curvature():
    linearization=ShapeLinearization(objective=1.0,
        objective_gradient=np.array([-1.0,0.0,1.0]),response=np.array([0.0]),
        response_jacobian=np.array([[1.0,0.0,0.0]]),response_target=np.array([0.0]),
        response_band=np.array([.05]))
    L=np.array([[1.0,-2.0,1.0]])
    update=solve_shape_lp(np.zeros(3),linearization,move_limit=.1,
        laplacian=L,curvature_limit=.02)
    assert update.predicted_max_band_ratio<=1.0+1e-10
    assert abs((L@update.parameters).item())<=.02+1e-10
    assert np.all(np.abs(update.delta)<=.1+1e-12)


def test_shape_lp_rejects_nonfinite_auxiliary_constraints():
    import pytest

    linearization=ShapeLinearization(objective=0.0,
        objective_gradient=np.array([1.0]),response=np.empty(0),
        response_jacobian=np.zeros((0,1)),response_target=np.empty(0),
        response_band=np.empty(0))
    with pytest.raises(ValueError,match="parameter bounds"):
        solve_shape_lp([0.0],linearization,move_limit=0.1,
            parameter_bounds=([np.nan],[1.0]))
    with pytest.raises(ValueError,match="A_ub requires b_ub"):
        solve_shape_lp([0.0],linearization,move_limit=0.1,A_ub=[[1.0]])
    with pytest.raises(ValueError,match="curvature_limit"):
        solve_shape_lp([0.0],linearization,move_limit=0.1,
            laplacian=[[1.0]],curvature_limit=np.inf)


def test_ngsolve_gettrafo_production_hex_closes_full_vim_scaling_tangent():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis_hex,build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1,
        mapping=lambda x,y,z:(1.1*x+.03*y*z,.9*y+.02*x*z,1.2*z+.04*x*y))
    fes=ng.HDiv(mesh,order=1)
    with ng.TaskManager():
        cb=_charge_basis_hex(fes,cob_quad=3)
        B,gram,_=build_charge_gram(fes,eps=1e-10,leafsize=256,eta=2.0)
        mode=ng.GridFunction(ng.VectorH1(mesh,order=1))
        mode.Set(ng.CF((ng.x,ng.y,ng.z)))
        result=linearize_production_vim_from_ngsolve(fes=fes,
            deformation_modes=[mode],charge_basis=cb,charge_gram=gram,
            charge_map=B,applied_coefficients=np.ones(fes.ndof),inv_chi=.2,
            family="hex")
        _,mf_charge,mf=linearize_production_vim_matrix_free_from_ngsolve(
            fes=fes,deformation_modes=[mode],charge_basis=cb,charge_gram=gram,
            charge_map=B,inv_chi=.2,family="hex",eps=1e-12,leaf=256,eta=2.0)
    np.testing.assert_allclose(result.charge_gram.jacobian[0],
        -result.charge_gram.matrix,rtol=4e-10,atol=4e-13)
    np.testing.assert_allclose(result.operator.matrix_jacobian[0],
        -result.operator.matrix,rtol=2e-9,atol=2e-11)
    np.testing.assert_allclose(result.operator.rhs_jacobian[0],
        -result.operator.rhs,rtol=2e-9,atol=2e-11)
    probe=np.linspace(-.3,.6,fes.ndof)
    np.testing.assert_allclose(mf.matvec(probe),result.operator.matrix@probe,
        rtol=2e-11,atol=2e-12)
    np.testing.assert_allclose(mf.directional_matvec(0,probe),
        result.operator.matrix_jacobian[0]@probe,rtol=2e-10,atol=2e-11)
    assert not hasattr(mf_charge,"matrix") or mf_charge.matrix is gram


def test_gettrafo_sampler_routes_tet_and_mixed_wedge_face_lattices():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis,_charge_basis_wedge
    for family,mesh,builder in (
        ("tet",MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1),
         lambda f:_charge_basis(f,4)),
        ("wedge",MakeStructured3DMesh(prism=True,nx=1,ny=1,nz=1),
         _charge_basis_wedge)):
        fes=ng.HDiv(mesh,order=1);vf=ng.VectorH1(mesh,order=1)
        mode=ng.GridFunction(vf);mode.Set(ng.CF((.01*ng.x,.02*ng.y,.03*ng.z)))
        with ng.TaskManager(): cb=builder(fes)
        sampled=sample_production_gettrafo_displacements(
            fes,[mode],cb,family=family)
        assert sampled.family==family and all(x.shape[0]==1 for x in sampled.cell)
        if family=="tet":
            assert all(x.shape[1:]==(4,3) for x in sampled.cell)
            assert all(x.shape[1:]==(3,3) for x in sampled.face)
        else:
            assert all(x.shape[1:]==(18,3) for x in sampled.cell)
            assert {x.shape[1] for x in sampled.face}=={6,9}


def test_gettrafo_sampler_uses_reference_vertices_and_restores_live_deformation():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis
    mesh=MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=1);space=ng.VectorH1(mesh,order=1)
    current=ng.GridFunction(space)
    current.Set(ng.CF((.07*ng.x-.02*ng.y,.03*ng.y,.04*ng.z)))
    mode=ng.GridFunction(space)
    mode.Set(ng.CF((.11*ng.x+.01*ng.y,-.05*ng.y,.08*ng.z)))
    with ng.TaskManager():
        mesh.SetDeformation(current)
        try:
            cb=_charge_basis(fes,quad=4,materialize_mass=False)
            sampled=sample_production_gettrafo_displacements(
                fes,[mode],cb,family="tet")
            assert mesh.deformation is not None
            np.testing.assert_allclose(
                mesh.deformation.vec.FV().NumPy(),
                current.vec.FV().NumPy(),rtol=0.0,atol=0.0)
        finally:
            mesh.UnsetDeformation()
    reference=np.asarray(cb["reference_vV"],dtype=float)
    physical=np.asarray(cb["vV"],dtype=float)
    assert np.max(np.abs(physical-reference))>1e-3
    expected=np.stack((.11*reference[...,0]+.01*reference[...,1],
                       -.05*reference[...,1],.08*reference[...,2]),axis=-1)
    np.testing.assert_allclose(np.asarray(sampled.cell)[:,0],expected,
                               rtol=2e-13,atol=2e-13)


def test_production_wedge_self_block_python_boundaries_preserve_host_mode_order():
    class Gram:
        def wedge_volume_self_block_directional_derivative(self,host,velocity):
            return np.full((2,2),host+np.asarray(velocity)[0,0])
        def wedge_face_self_block_directional_derivative(self,host,velocity):
            return np.full((3,3),host+np.asarray(velocity)[0,0])
    volume=production_wedge_volume_self_block_derivatives(
        Gram(),[np.stack([np.zeros((18,3)),np.ones((18,3))])])
    faces=production_wedge_face_self_block_derivatives(
        Gram(),[np.stack([np.zeros((6,3)),np.ones((6,3))]),
                np.stack([2*np.ones((9,3)),3*np.ones((9,3))])])
    assert volume[0].shape==(2,2,2)
    assert np.all(volume[0][1]==1)
    assert faces[0].shape==(2,3,3) and faces[1].shape==(2,3,3)
    assert np.all(faces[1][0]==3)


def test_production_wedge_full_gram_python_boundary_preserves_mode_order():
    class Gram:
        def wedge_charge_gram_directional_derivative(self,cells,faces):
            return np.eye(2)*(np.asarray(cells)[0,0,0]+np.asarray(faces)[0,0,0])
    cells=np.zeros((2,1,18,3));faces=np.zeros((2,3,9,3))
    cells[1,0,0,0]=2;faces[1,0,0,0]=3
    out=production_wedge_charge_gram_derivatives(Gram(),cells,faces)
    assert out.shape==(2,2,2) and np.array_equal(out[1],5*np.eye(2))


def test_native_production_wedge_self_block_derivative_invariants():
    import pytest
    ng=pytest.importorskip("ngsolve")
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _build_charge_gram_wedge, _charge_basis_wedge
    try:
        mesh=MakeStructured3DMesh(nx=1,ny=1,nz=1,prism=True)
    except TypeError:
        pytest.skip("this NGSolve build cannot generate a prism mesh")
    with ng.TaskManager():
        fes=ng.HDiv(mesh,order=1)
        cb=_charge_basis_wedge(fes)
        _,gram,_,_=_build_charge_gram_wedge(fes,eps=1e-14,leafsize=256)
    def value_block(kind,host):
        ids=np.flatnonzero((np.asarray(cb["kind"])==kind)&(np.asarray(cb["host"])==host))
        return np.asarray([[gram.entry(int(i),int(j)) for j in ids] for i in ids])
    cell=np.asarray(cb["cell_nodes"]).reshape(-1,18,3)[0]
    dc=gram.wedge_volume_self_block_directional_derivative(0,cell)
    tc=gram.wedge_volume_self_block_directional_derivative(0,np.ones((18,3)))
    assert np.linalg.norm(dc+value_block(0,0))/np.linalg.norm(value_block(0,0))<2e-10
    assert np.linalg.norm(tc)<2e-12 and np.array_equal(dc,dc.T)
    face_nodes=np.asarray(cb["face_nodes"]).reshape(-1,9,3)
    for host,ft in enumerate(np.asarray(cb["face_type"])):
        nn=6 if ft==0 else 9; nodes=face_nodes[host,:nn]
        df=gram.wedge_face_self_block_directional_derivative(host,nodes)
        tf=gram.wedge_face_self_block_directional_derivative(host,np.ones((nn,3)))
        block=value_block(1,host)
        assert np.linalg.norm(df+block)/np.linalg.norm(block)<2e-10
        assert np.linalg.norm(tf)<2e-12 and np.array_equal(df,df.T)


def test_native_production_wedge_self_block_derivative_matches_general_affine_fd():
    import pytest
    ng=pytest.importorskip("ngsolve")
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _build_charge_gram_wedge, _charge_basis_wedge
    A=np.array([[.17,-.09,.04],[.06,-.13,.08],[-.03,.05,.11]])
    base=np.array([[1.0,.12,-.04],[.03,.91,.08],[-.06,.05,1.13]])
    shift=np.array([.02,-.01,.03])
    def make(eps):
        def mapping(x,y,z):
            q=np.array([x,y,z],dtype=float)
            r=base@q+eps*(A@q+shift)
            return tuple(r)
        try:
            mesh=MakeStructured3DMesh(nx=1,ny=1,nz=1,prism=True,mapping=mapping)
        except TypeError:
            pytest.skip("this NGSolve build cannot generate a prism mesh")
        with ng.TaskManager():
            fes=ng.HDiv(mesh,order=1)
            cb=_charge_basis_wedge(fes)
            _,gram,_,_=_build_charge_gram_wedge(fes,eps=1e-14,leafsize=256)
        return cb,gram
    step=2e-6
    cb0,g0=make(0.0); cbp,gp=make(step); cbm,gm=make(-step)
    kind=np.asarray(cb0["kind"]); hosts=np.asarray(cb0["host"])
    def block(gram,k,h):
        ids=np.flatnonzero((kind==k)&(hosts==h))
        return np.asarray([[gram.entry(int(i),int(j)) for j in ids] for i in ids])
    cells=np.asarray(cb0["cell_nodes"]).reshape(-1,18,3)
    faces=np.asarray(cb0["face_nodes"]).reshape(-1,9,3)
    types=np.asarray(cb0["face_type"])
    cases=[("volume",0,0,cells[0])]
    for ft,name in ((0,"tri"),(1,"quad")):
        h=int(np.flatnonzero(types==ft)[0]); nn=6 if ft==0 else 9
        cases.append((name,1,h,faces[h,:nn]))
    for name,k,h,nodes in cases:
        reference=nodes@np.linalg.inv(base).T
        velocity=reference@A.T+shift
        analytic=(g0.wedge_volume_self_block_directional_derivative(h,velocity) if k==0
                  else g0.wedge_face_self_block_directional_derivative(h,velocity))
        fd=(block(gp,k,h)-block(gm,k,h))/(2*step)
        relative=np.linalg.norm(analytic-fd)/np.linalg.norm(fd)
        assert relative<2e-7,(name,relative)


def test_native_production_wedge_full_gram_derivative_translation_scale_and_local_fd():
    import pytest
    ng=pytest.importorskip("ngsolve")
    from ngsolve.meshes import MakeStructured3DMesh
    from radia import _radia_pybind as rp
    from radia.vim._vim import (_charge_basis_wedge,_g01,_SYM5_TET,_SYM5_TRI,
                                _f64_buffer,_i32_buffer)
    try:mesh=MakeStructured3DMesh(nx=2,ny=1,nz=1,prism=True)
    except TypeError:pytest.skip("this NGSolve build cannot generate a prism mesh")
    with ng.TaskManager():cb=_charge_basis_wedge(ng.HDiv(mesh,order=1))
    cells=np.asarray(cb["cell_nodes"],dtype=float).reshape(-1,18,3)
    faces=np.asarray(cb["face_nodes"],dtype=float).reshape(-1,9,3).copy()
    types=np.asarray(cb["face_type"],dtype=np.int32)
    glo,gwo=_g01(6);gli,gwi=_g01(5)
    tri=ng.IntegrationRule(ng.ET.TRIG,5)
    ftp=np.asarray([(ip.point[0],ip.point[1]) for ip in tri]);ftw=np.asarray([ip.weight for ip in tri])
    def make(c,f,near=.6,far=1.5):
        return rp._ChargeGramHMatrix(
            wedge_cell_nodes=_f64_buffer(c),face_nodes=_f64_buffer(f),face_type=_i32_buffer(types),
            n_el=int(cb["n_el"]),n_bf=int(cb["n_bf"]),charge_host=_i32_buffer(cb["host"]),
            charge_kind=_i32_buffer(cb["kind"]),charge_expo=_i32_buffer(cb["expo"]),
            sym_tet_pts=_f64_buffer(_SYM5_TET[0]),sym_tet_w=_f64_buffer(_SYM5_TET[1]),
            sym_tri_pts=_f64_buffer(_SYM5_TRI[0]),sym_tri_w=_f64_buffer(_SYM5_TRI[1]),
            field_tri_pts=_f64_buffer(ftp),field_tri_w=_f64_buffer(ftw),
            gl_out=_f64_buffer(glo),gw_out=_f64_buffer(gwo),gl_in=_f64_buffer(gli),gw_in=_f64_buffer(gwi),
            far_tet_pts=_f64_buffer(_SYM5_TET[0]),far_tet_w=_f64_buffer(_SYM5_TET[1]),
            far_tri_pts=_f64_buffer(_SYM5_TRI[0]),far_tri_w=_f64_buffer(_SYM5_TRI[1]),
            near_grade=near,far_inner_factor=far,image_masks=np.empty(0,np.int32),
            image_signs=np.empty(0),eps=1e-14,leaf=256,eta=2.,build=False)
    n=len(cb["kind"])
    def dense(g):return np.asarray([[g.entry(i,j) for j in range(n)] for i in range(n)])
    g0=make(cells,faces);G=dense(g0)
    ones_c=np.ones_like(cells);ones_f=np.ones_like(faces)
    dt=np.asarray(g0.wedge_charge_gram_directional_derivative(ones_c,ones_f))
    ds=np.asarray(g0.wedge_charge_gram_directional_derivative(cells,faces))
    dop=g0.directional_derivative_operator("wedge",cells,faces,
        eps=1e-12,leaf=256,eta=2.0)
    probe=np.linspace(-.4,.7,n)
    np.testing.assert_allclose(dop.matvec_sym(probe),ds@probe,rtol=2e-12,atol=2e-12)
    right=probe[::-1].copy()
    contractions=g0.directional_derivative_contractions("wedge",
        np.ascontiguousarray(np.stack([ones_c,cells])),
        np.ascontiguousarray(np.stack([ones_f,faces])),probe,right)
    contractions_many=g0.directional_derivative_contractions_many("wedge",
        np.ascontiguousarray(np.stack([ones_c,cells])),
        np.ascontiguousarray(np.stack([ones_f,faces])),
        np.ascontiguousarray(np.stack([probe,2*probe-right])),right)
    np.testing.assert_allclose(contractions,[probe@dt@right,probe@ds@right],
        rtol=3e-13,atol=3e-13)
    np.testing.assert_allclose(contractions_many,
        [[probe@dt@right,probe@ds@right],
         [(2*probe-right)@dt@right,(2*probe-right)@ds@right]],
        rtol=3e-13,atol=3e-13)
    assert np.linalg.norm(dt)<3e-11
    assert np.linalg.norm(ds+G)/np.linalg.norm(G)<3e-10
    assert np.array_equal(dt,dt.T) and np.array_equal(ds,ds.T)
    # A continuous piecewise-affine hat field is localized to the interior
    # x-plane.  It agrees exactly on duplicate cell/face nodes and preserves
    # affine hosts, just like an H1/GetTrafo P1 deformation mode.
    def local_velocity(x):
        h=np.maximum(0.,1.-2.*np.abs(x[...,0]-.5))
        return h[...,None]*np.array([.031,-.027,.023])
    vc,vf=local_velocity(cells),local_velocity(faces)
    gf=make(cells,faces,near=0.,far=0.)
    analytic=np.asarray(gf.wedge_charge_gram_directional_derivative(vc,vf))
    step=2e-6
    Gp=dense(make(cells+step*vc,faces+step*vf,near=0.,far=0.));Gm=dense(make(cells-step*vc,faces-step*vf,near=0.,far=0.))
    fd=(Gp-Gm)/(2*step)
    # Touching-but-separately-integrated host pairs have a discrete nearest-site
    # quadrature dispatch and are not a differentiable FD oracle.  Compare the
    # localized non-self derivative on the farthest cell pair; all self blocks
    # are independently locked above, while translation/scale cover full dG.
    centers=cells.mean(axis=1);dist=np.linalg.norm(centers[:,None]-centers[None,:],axis=2)
    ha,hb=np.unravel_index(np.argmax(dist),dist.shape)
    kind=np.asarray(cb["kind"]);host=np.asarray(cb["host"])
    ia=np.flatnonzero((kind==0)&(host==ha));ib=np.flatnonzero((kind==0)&(host==hb))
    aa=analytic[np.ix_(ia,ib)];ff=fd[np.ix_(ia,ib)]
    relative=np.linalg.norm(aa-ff)/np.linalg.norm(ff)
    assert relative<3e-6,(relative,ha,hb,np.linalg.norm(aa),np.linalg.norm(ff))


def test_affine_self_term_derivative_closes_tet_hex_wedge_diagonal():
    cells={
        "tet":np.array([[0.,0,0],[1.1,0,0],[.1,.9,0],[.2,.1,1.2]]),
        "hex":np.array([[0.,0,0],[1.,0,0],[1.,1,0],[0,1,0],
                        [0,0,1],[1.,0,1],[1.,1,1],[0,1,1.]]),
        "wedge":np.array([[0.,0,0],[1.,0,0],[0,1.,0],[0,0,1],[1.,0,1],[0,1.,1]]),
    }
    rng=np.random.default_rng(741)
    for kind,nodes in cells.items():
        direction=rng.normal(size=nodes.shape)
        value,jac=affine_cell_self_energy_shape_derivative(
            kind,nodes,np.stack([np.ones_like(nodes),nodes,direction]))
        assert value>0 and np.all(np.isfinite(jac))
        assert jac[0]==0.0  # rigid translation is removed analytically
        np.testing.assert_allclose(jac[1],5*value,rtol=3e-4,atol=2e-8)
        eps=2e-6
        plus=affine_cell_self_energy_shape_derivative(
            kind,nodes+eps*direction,np.zeros((0,*nodes.shape)))[0]
        minus=affine_cell_self_energy_shape_derivative(
            kind,nodes-eps*direction,np.zeros((0,*nodes.shape)))[0]
        np.testing.assert_allclose((plus-minus)/(2*eps),jac[2],rtol=2e-3,atol=2e-6)


def test_full_charge_gram_combines_nonself_and_analytic_self_tangents():
    cells=[np.array([[0.,0,0],[1.,0,0],[0,1.,0],[0,0,1.]]),
           np.array([[2.,.1,0],[3.1,.1,0],[2.,1.2,0],[2.,.1,.8]])]
    velocity=[np.array([[0,0,0],[.1,0,0],[0,.02,0],[0,0,-.03]]),
              np.array([[.03,0,0],[-.02,0,0],[0,.04,0],[0,0,.01]])]
    def geometry(cells_now):
        points=np.array([c.mean(axis=0) for c in cells_now])
        weights=[]
        for c in cells_now:
            J=np.column_stack((c[1]-c[0],c[2]-c[0],c[3]-c[0]))
            weights.append(abs(np.linalg.det(J))/6)
        return points,np.array(weights)
    points,weights=geometry(cells)
    point_velocity=np.array([[v.mean(axis=0) for v in velocity]])
    eps=2e-6
    plus_cells=[c+eps*v for c,v in zip(cells,velocity)]
    minus_cells=[c-eps*v for c,v in zip(cells,velocity)]
    plus_points,plus_weights=geometry(plus_cells); minus_points,minus_weights=geometry(minus_cells)
    rel=((plus_weights-minus_weights)/(2*eps))/weights
    model=linearize_laplace_charge_gram(points,weights,point_velocity,
        relative_weight_derivatives=rel[None,:],self_cell_types=["tet","tet"],
        self_nodes=cells,self_node_displacements=[np.array([v]) for v in velocity])
    plus=linearize_laplace_charge_gram(plus_points,plus_weights,np.zeros((0,2,3)),
        self_cell_types=["tet","tet"],self_nodes=plus_cells,
        self_node_displacements=[np.zeros((0,4,3)),np.zeros((0,4,3))]).matrix
    minus=linearize_laplace_charge_gram(minus_points,minus_weights,np.zeros((0,2,3)),
        self_cell_types=["tet","tet"],self_nodes=minus_cells,
        self_node_displacements=[np.zeros((0,4,3)),np.zeros((0,4,3))]).matrix
    assert np.all(np.diag(model.matrix)>0)
    np.testing.assert_allclose((plus-minus)/(2*eps),model.jacobian[0],rtol=2e-3,atol=2e-6)


def test_ngsolve_hdiv_mass_shape_tangent_uses_piola_weak_form():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _csr
    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=1); vf=ng.VectorH1(mesh,order=1)
    velocity=ng.GridFunction(vf); velocity.Set(ng.CF((.07*ng.x,-.03*ng.y,.02*ng.z)))
    with ng.TaskManager():
        mass,dmass,dB=assemble_ngsolve_hdiv_shape_tangents(fes,[velocity],np.eye(fes.ndof))
        sparse_mass,sparse_dmass,sparse_dB=assemble_ngsolve_hdiv_shape_tangents(
            fes,[velocity],np.eye(fes.ndof),sparse=True)
        eps=2e-6; trial=ng.GridFunction(vf); trial.vec.data=eps*velocity.vec
        mesh.SetDeformation(trial)
        u,v=fes.TnT(); shifted=ng.BilinearForm(fes); shifted+=u*v*ng.dx; shifted.Assemble()
        plus=_csr(shifted).toarray(); mesh.UnsetDeformation()
        trial.vec.data=-eps*velocity.vec; mesh.SetDeformation(trial)
        shifted=ng.BilinearForm(fes); shifted+=u*v*ng.dx; shifted.Assemble()
        minus=_csr(shifted).toarray(); mesh.UnsetDeformation()
    np.testing.assert_allclose((plus-minus)/(2*eps),dmass[0],rtol=2e-6,atol=2e-8)
    assert mass.shape==(fes.ndof,fes.ndof) and np.count_nonzero(dB)==0
    import scipy.sparse as sp
    assert sp.isspmatrix_csr(sparse_mass)
    assert len(sparse_dmass)==1 and sp.isspmatrix_csr(sparse_dmass[0])
    assert len(sparse_dB)==1 and sp.isspmatrix_csr(sparse_dB[0])
    np.testing.assert_allclose(sparse_mass.toarray(),mass,rtol=0,atol=0)
    np.testing.assert_allclose(sparse_dmass[0].toarray(),dmass[0],rtol=0,atol=0)


def test_production_hex_self_block_python_boundary_preserves_host_mode_order():
    class Gram:
        def hex_volume_self_block_directional_derivative(self,host,velocity):
            return np.full((2,2),host+np.sum(velocity))
    modes=[np.zeros((2,27,3)),np.ones((2,27,3))]
    blocks=production_hex_volume_self_block_derivatives(Gram(),modes)
    assert len(blocks)==2 and blocks[0].shape==(2,2,2)
    np.testing.assert_allclose(blocks[0],0)
    np.testing.assert_allclose(blocks[1],82)


def test_production_hex_face_self_block_python_boundary_preserves_host_mode_order():
    class Gram:
        def hex_face_self_block_directional_derivative(self,host,velocity):
            return np.full((3,3),host+np.sum(velocity))
    modes=[np.zeros((2,9,3)),np.ones((2,9,3))]
    blocks=production_hex_face_self_block_derivatives(Gram(),modes)
    assert len(blocks)==2 and blocks[0].shape==(2,3,3)
    np.testing.assert_allclose(blocks[0],0)
    np.testing.assert_allclose(blocks[1],28)


def test_production_tet_self_block_python_boundaries_preserve_host_mode_order():
    class Gram:
        def tet_volume_self_block_directional_derivative(self,host,velocity): return np.full((2,2),host+np.sum(velocity))
        def tet_face_self_block_directional_derivative(self,host,velocity): return np.full((3,3),host+np.sum(velocity))
    volume=production_tet_volume_self_block_derivatives(Gram(),[np.ones((2,4,3))])
    face=production_tet_face_self_block_derivatives(Gram(),[np.ones((2,3,3))])
    assert volume[0].shape==(2,2,2) and face[0].shape==(2,3,3)
    np.testing.assert_allclose(volume[0],12);np.testing.assert_allclose(face[0],9)


def test_native_production_hex_volume_self_block_derivative_invariants():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis_hex,build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1,
        mapping=lambda x,y,z:(x,y,z+.12*x*y))
    fes=ng.HDiv(mesh,order=1)
    with ng.TaskManager():
        cb=_charge_basis_hex(fes,cob_quad=3)
        # ACA recompression at +/-epsilon is not differentiable; this
        # regression compares against the production dense block path.
        _,gram,_=build_charge_gram(fes,eps=1e-10,leafsize=256,eta=2.0)
    nodes=np.asarray(cb["cell_nodes"]).reshape(-1,27,3)[0]
    translation=gram.hex_volume_self_block_directional_derivative(0,np.ones((27,3)))
    scaling=gram.hex_volume_self_block_directional_derivative(0,nodes)
    n=scaling.shape[0]
    value=np.array([[gram.entry(i,j) for j in range(n)] for i in range(n)])
    np.testing.assert_allclose(translation,0,atol=2e-17)
    np.testing.assert_allclose(scaling,scaling.T,rtol=0,atol=0)
    # Reference Piola charges have fixed measure; uniform physical scaling
    # therefore differentiates the Laplace kernel with homogeneity -1.
    np.testing.assert_allclose(scaling,-value,rtol=3e-13,atol=3e-15)


def test_native_production_tet_self_block_derivatives_match_scaling_and_fd():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis,build_charge_gram

    gradient=np.array([[.11,-.04,.03],[.02,.08,-.05],[-.01,.06,.09]])
    offset=np.array([.013,-.021,.017]); epsilon=2e-6
    def build(step):
        mesh=MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1,
            mapping=lambda x,y,z:tuple(np.array([x,y,z])+step*epsilon*(gradient@np.array([x,y,z])+offset)))
        fes=ng.HDiv(mesh,order=1)
        with ng.TaskManager():
            cb=_charge_basis(fes,4)
            _,gram,_=build_charge_gram(fes,eps=1e-10,leafsize=256,eta=2.0)
        return cb,gram
    cb,gram=build(0); _,plus=build(1); _,minus=build(-1)
    kinds=np.asarray(cb["kind"]);hosts=np.asarray(cb["host"])
    for kind,nodes,key in ((0,np.asarray(cb["vV"]),"tet_volume_self_block_directional_derivative"),
                           (1,np.asarray(cb["bV"]),"tet_face_self_block_directional_derivative")):
        method=getattr(gram,key)
        for host,vertices in enumerate(nodes):
            ids=np.flatnonzero((kinds==kind)&(hosts==host))
            value=np.array([[gram.entry(int(i),int(j)) for j in ids] for i in ids])
            value_plus=np.array([[plus.entry(int(i),int(j)) for j in ids] for i in ids])
            value_minus=np.array([[minus.entry(int(i),int(j)) for j in ids] for i in ids])
            scaling=method(host,vertices)
            translation=method(host,np.ones_like(vertices))
            velocity=vertices@gradient.T+offset
            derivative=method(host,velocity)
            fd=(value_plus-value_minus)/(2*epsilon)
            np.testing.assert_allclose(translation,0,atol=8e-17)
            np.testing.assert_allclose(scaling,scaling.T,rtol=0,atol=0)
            # Flat TET ChargeGram stores physical charge monomials and both
            # physical measures: volume-volume and surface-surface blocks are
            # homogeneous of degree five and three respectively.  The Piola
            # B-map derivatives close the final B.T@G@B product separately.
            np.testing.assert_allclose(scaling,(5 if kind==0 else 3)*value,rtol=8e-12,atol=4e-15)
            np.testing.assert_allclose(derivative,fd,rtol=2e-7,atol=3e-11)

    # Reverse the physical orientation so the same regression covers both
    # signs of det(E); the analytic moments are physical (unsigned), while
    # the affine coordinate map retains the production orientation convention.
    reflection=np.diag([-1.,1.,1.])
    def build_reflected(step):
        mesh=MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1,
            mapping=lambda x,y,z:tuple((np.eye(3)+step*epsilon*gradient)@(reflection@np.array([x,y,z]))+step*epsilon*offset))
        fes=ng.HDiv(mesh,order=1)
        with ng.TaskManager():
            local_cb=_charge_basis(fes,4)
            _,local_gram,_=build_charge_gram(fes,eps=1e-10,leafsize=16,eta=2.0)
        return local_cb,local_gram
    rcb,rgram=build_reflected(0);_,rplus=build_reflected(1);_,rminus=build_reflected(-1)
    rk=np.asarray(rcb["kind"]);rh=np.asarray(rcb["host"]);vertices=np.asarray(rcb["vV"])[0]
    ids=np.flatnonzero((rk==0)&(rh==0));value=np.array([[rgram.entry(int(i),int(j)) for j in ids] for i in ids])
    derivative=rgram.tet_volume_self_block_directional_derivative(0,vertices@gradient.T+offset)
    fd=np.array([[(rplus.entry(int(i),int(j))-rminus.entry(int(i),int(j)))/(2*epsilon) for j in ids] for i in ids])
    np.testing.assert_allclose(rgram.tet_volume_self_block_directional_derivative(0,vertices),5*value,rtol=8e-12,atol=4e-15)
    np.testing.assert_allclose(derivative,fd,rtol=2e-7,atol=3e-11)


def test_native_production_tet_complete_gram_and_piola_product_derivative():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.topology_optimization import production_tet_charge_gram_derivatives
    from radia.vim._vim import _charge_basis,build_charge_gram

    gradient=np.array([[.073,-.031,.019],[.014,.052,-.027],[-.022,.041,.064]])
    offset=np.array([.009,-.015,.011]); epsilon=1e-6
    def build(step, scaling=False, translation=False):
        def mapping(x,y,z):
            point=np.array([x,y,z])
            if scaling: velocity=point
            elif translation: velocity=np.ones(3)
            else: velocity=gradient@point+offset
            return tuple(point+step*epsilon*velocity)
        mesh=MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1,mapping=mapping)
        fes=ng.HDiv(mesh,order=1)
        with ng.TaskManager():
            cb=_charge_basis(fes,4)
            # Keep this derivative regression on one dense leaf.  Rebuilding
            # ACA factors at +/-epsilon is not a differentiable reference.
            _,gram,_=build_charge_gram(fes,eps=1e-10,leafsize=256,eta=2.0)
        n=len(cb["host"])
        G=np.array([[gram.entry(i,j) for j in range(n)] for i in range(n)])
        return cb,gram,G

    cb,gram,G=build(0); cbp,_,Gp=build(1); cbm,_,Gm=build(-1)
    cells=np.asarray(cb["vV"]); faces=np.asarray(cb["bV"])
    cell_v=cells@gradient.T+offset; face_v=faces@gradient.T+offset
    dG,dB=production_tet_charge_gram_derivatives(
        gram,cell_v[None,...],face_v[None,...],cb["B"])
    dG=dG[0]; dB=dB[0].toarray(); B=cb["B"].toarray()
    dop=gram.directional_derivative_operator("tet",cell_v,face_v,
        eps=1e-12,leaf=256,eta=2.0)
    probe=np.linspace(-.6,.8,dG.shape[0])
    np.testing.assert_allclose(dop.matvec_sym(probe),dG@probe,rtol=2e-12,atol=2e-12)
    right=probe[::-1].copy()
    contraction=gram.directional_derivative_contractions("tet",
        np.ascontiguousarray(cell_v[None,...]),np.ascontiguousarray(face_v[None,...]),
        probe,right)
    contraction_many=gram.directional_derivative_contractions_many("tet",
        np.ascontiguousarray(cell_v[None,...]),np.ascontiguousarray(face_v[None,...]),
        np.ascontiguousarray(np.stack([probe,2*probe-right])),right)
    np.testing.assert_allclose(contraction,[probe@dG@right],rtol=3e-13,atol=3e-13)
    np.testing.assert_allclose(contraction_many,
        [[probe@dG@right],[(2*probe-right)@dG@right]],rtol=3e-13,atol=3e-13)
    fdG=(Gp-Gm)/(2*epsilon)
    np.testing.assert_allclose(dG,dG.T,rtol=0,atol=0)
    np.testing.assert_allclose(dG,fdG,rtol=4e-7,atol=8e-11)
    N=B.T@G@B
    dN=dB.T@G@B+B.T@dG@B+B.T@G@dB
    Np=cbp["B"].toarray().T@Gp@cbp["B"].toarray()
    Nm=cbm["B"].toarray().T@Gm@cbm["B"].toarray()
    np.testing.assert_allclose(dN,(Np-Nm)/(2*epsilon),rtol=7e-7,atol=2e-10)

    # Translation leaves every raw block and Piola row unchanged.
    _,gt,_=build(0,translation=True)
    dGt,rt=production_tet_charge_gram_derivatives(
        gt,np.ones_like(cells)[None,...],np.ones_like(faces)[None,...])
    np.testing.assert_allclose(dGt,0,atol=2e-16)
    np.testing.assert_allclose(rt,0,atol=2e-16)

    # Under uniform scaling raw VV/VF/FF blocks have degrees 5/4/3,
    # while dB/B is -3/-2.  Every block of B.T@G@B is therefore degree -1.
    dGs,dBs=production_tet_charge_gram_derivatives(
        gram,cells[None,...],faces[None,...],cb["B"])
    dBs=dBs[0].toarray()
    np.testing.assert_allclose(dBs.T@G@B+B.T@dGs[0]@B+B.T@G@dBs,-N,
                               rtol=2e-11,atol=2e-13)


def test_native_tet_cluster_contractions_many_matches_across_thread_counts():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis,build_charge_gram

    mesh=MakeStructured3DMesh(hexes=False,nx=2,ny=2,nz=1,
        mapping=lambda x,y,z:(x+.03*y*z,y+.02*x*z,z+.01*x*y))
    fes=ng.HDiv(mesh,order=1)
    ng.SetNumThreads(4)
    with ng.TaskManager():
        cb=_charge_basis(fes,4)
        _,gram,_=build_charge_gram(
            fes,eps=1e-7,leafsize=8,eta=2.0)
    assert gram.stats()["n_lowrank"]>0
    cells=np.asarray(cb["vV"]);faces=np.asarray(cb["bV"])
    velocity=lambda x:np.stack((
        .02*x[...,0]-.01*x[...,1],
        .015*x[...,1]+.005*x[...,2],
        -.012*x[...,0]+.008*x[...,2]),axis=-1)
    cell_velocity=velocity(cells);face_velocity=velocity(faces)
    cell_modes=np.ascontiguousarray(np.stack((
        cell_velocity,1.7*cell_velocity)))
    face_modes=np.ascontiguousarray(np.stack((
        face_velocity,1.7*face_velocity)))
    n=int(gram.stats()["n_dof"])
    right=np.ascontiguousarray(np.cos(np.arange(n)))
    left=np.ascontiguousarray(np.stack((
        np.linspace(-.4,.7,n),np.sin(np.arange(n)))))
    ng.SetNumThreads(1)
    with ng.TaskManager():
        serial=np.asarray(gram.directional_derivative_contractions_many(
            "tet",cell_modes,face_modes,left,right))
    ng.SetNumThreads(4)
    with ng.TaskManager():
        parallel=np.asarray(gram.directional_derivative_contractions_many(
            "tet",cell_modes,face_modes,left,right))
    np.testing.assert_allclose(parallel,serial,rtol=2e-12,atol=2e-12)


def test_native_tet_configured_field_rows_directional_derivative_matches_fd():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis,build_charge_gram

    gradient=np.array([[.073,-.031,.019],[.014,.052,-.027],[-.022,.041,.064]])
    offset=np.array([.009,-.015,.011])
    observations=np.array([[.31,.27,1e-3],[.5,.4,1.7],[-.2,.8,.5]])
    weights=np.array([
        [[.2,-.3,.7],[.1,.4,-.2],[-.5,.2,.3]],
        [[.7,.1,-.4],[.3,-.2,.6],[.1,.5,-.3]]])
    image_masks=(2,4,6);image_signs=(1.,-1.,-1.)

    def build(step):
        def mapping(x,y,z):
            point=np.array([x,y,z])
            return tuple(point+step*(gradient@point+offset))
        mesh=MakeStructured3DMesh(
            hexes=False,nx=1,ny=1,nz=1,mapping=mapping)
        fes=ng.HDiv(mesh,order=1)
        with ng.TaskManager():
            charge_basis=_charge_basis(fes,4)
            _,gram,_=build_charge_gram(
                fes,eps=1e-10,leafsize=256,eta=2.0,
                image_masks=image_masks,image_signs=image_signs)
        return charge_basis,gram

    charge_basis,gram=build(0.0)
    cells=np.asarray(charge_basis["vV"])
    faces=np.asarray(charge_basis["bV"])
    analytic=np.asarray(
        gram.configured_field_functional_rows_directional_derivative(
            observations,weights,(cells@gradient.T+offset)[None,...],
            (faces@gradient.T+offset)[None,...]))[0]
    step=2e-6
    _,plus=build(step);_,minus=build(-step)
    finite_difference=(
        np.asarray(plus.configured_field_functional_rows(observations,weights))
        -np.asarray(minus.configured_field_functional_rows(observations,weights))
        )/(2*step)

    assert analytic.flags.c_contiguous and np.all(np.isfinite(analytic))
    np.testing.assert_allclose(
        analytic,finite_difference,rtol=3e-6,atol=2e-8)


@pytest.mark.parametrize("outside_y",[-0.25,-5e-14])
def test_native_tet_configured_field_derivative_is_finite_on_coplanar_extension(
        outside_y):
    """Differentiate the smooth exterior limit where a panel plane hits a probe."""
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis,build_charge_gram

    translation=np.array([0.0,0.0,0.04])
    observations=np.array([[0.25,outside_y,0.0]],dtype=float)
    weights=np.array([[[0.0,0.0,1.0]],[[0.3,-0.2,0.4]]],dtype=float)

    def build(step):
        def mapping(x,y,z):
            point=np.array([x,y,z])
            return tuple(point+step*translation)
        mesh=MakeStructured3DMesh(
            hexes=False,nx=1,ny=1,nz=1,mapping=mapping)
        fes=ng.HDiv(mesh,order=1)
        with ng.TaskManager():
            charge_basis=_charge_basis(fes,4)
            _,gram,_=build_charge_gram(
                fes,eps=1e-10,leafsize=256,eta=2.0)
        return charge_basis,gram

    charge_basis,gram=build(0.0)
    cells=np.asarray(charge_basis["vV"])
    faces=np.asarray(charge_basis["bV"])
    cell_velocity=np.broadcast_to(translation,cells.shape).copy()
    face_velocity=np.broadcast_to(translation,faces.shape).copy()
    analytic=np.asarray(
        gram.configured_field_functional_rows_directional_derivative(
            observations,weights,cell_velocity[None,...],
            face_velocity[None,...]))[0]
    step=2e-6
    _,plus=build(step);_,minus=build(-step)
    finite_difference=(
        np.asarray(plus.configured_field_functional_rows(observations,weights))
        -np.asarray(minus.configured_field_functional_rows(observations,weights))
        )/(2*step)

    assert analytic.flags.c_contiguous and np.all(np.isfinite(analytic))
    assert np.all(np.isfinite(finite_difference))
    if abs(outside_y)>1e-10:
        np.testing.assert_allclose(
            analytic,finite_difference,rtol=8e-6,atol=3e-8)


def test_tet_streaming_shape_jacobian_uses_exact_configured_field_derivative():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis,build_charge_gram

    gradient=np.array([[.04,-.02,.01],[.01,.03,-.015],[-.02,.005,.025]])
    shift=np.array([.01,-.006,.004]);inv_chi=.2
    applied=ng.CF((.2,-.1,1.0))
    observations=np.array([[.3,.2,2e-3],[.7,.4,1.4],[-.1,.6,.3]])
    weights=np.array([
        [[0.,0.,.5],[0.,0.,.3],[0.,0.,.2]],
        [[.2,-.1,.4],[-.3,.5,.2],[.1,.2,-.4]]])

    def build(epsilon):
        def mapping(x,y,z):
            point=np.array([x,y,z])
            return tuple(point+epsilon*(gradient@point+shift))
        mesh=MakeStructured3DMesh(
            hexes=False,nx=1,ny=1,nz=1,mapping=mapping)
        fes=ng.HDiv(mesh,order=1)
        with ng.TaskManager():
            basis=_charge_basis(fes,4)
            B,gram,_=build_charge_gram(
                fes,eps=1e-12,leafsize=256,eta=2.0)
            rhs,_=assemble_ngsolve_hdiv_linear_form_shape_tangents(
                fes,applied,(),bonus_intorder=4)
        C=np.asarray(gram.configured_field_functional_rows(
            observations,weights))
        state=np.asarray(gram.solve_configured_linear_material_auto_prec_many(
            inv_chi,np.ascontiguousarray(rhs[None,:]),tol=1e-12,
            maxit=5000,mass_riesz=True)["m"])[0]
        return C@state,(mesh,fes,basis,B,gram,rhs,C,state)

    value,data=build(0.0);mesh,fes,basis,B,gram,rhs,C,state=data
    space=ng.VectorH1(mesh,order=1);mode=ng.GridFunction(space)
    mode.Set(ng.CF(tuple(
        gradient@np.array([ng.x,ng.y,ng.z],dtype=object)+shift)))
    with ng.TaskManager():
        _,drhs=assemble_ngsolve_hdiv_linear_form_shape_tangents(
            fes,applied,[mode],bonus_intorder=4)
        analytic=production_vim_functional_shape_jacobian_streaming(
            fes=fes,deformation_modes=[mode],charge_basis=basis,
            charge_gram=gram,charge_map=B,inv_chi=inv_chi,rhs=rhs,
            response_matrix=C,rhs_jacobian=drhs,
            response_observations=observations,response_weights=weights,
            family="tet",solve_tolerance=1e-12,
            mass_riesz=False,
            cluster_coarse_size=8,cluster_deflation_size=2,
            recycle_size=2)
        reused=production_vim_functional_shape_jacobian_streaming(
            fes=fes,deformation_modes=[mode],charge_basis=basis,
            charge_gram=gram,charge_map=B,inv_chi=inv_chi,rhs=rhs,
            response_matrix=C,rhs_jacobian=drhs,
            response_observations=observations,response_weights=weights,
            family="tet",solve_tolerance=1e-12,
            mass_riesz=False,
            cluster_coarse_size=8,cluster_deflation_size=2,
            recycle_size=2,state=state,state_iterations=123)
    epsilon=2e-6;plus,_=build(epsilon);minus,_=build(-epsilon)
    finite_difference=(plus-minus)/(2*epsilon)

    np.testing.assert_allclose(analytic.response,value,rtol=2e-12,atol=2e-12)
    np.testing.assert_allclose(
        analytic.response_jacobian[:,0],finite_difference,
        rtol=3e-4,atol=3e-7)
    np.testing.assert_allclose(reused.state,state,rtol=0,atol=0)
    np.testing.assert_allclose(reused.response,analytic.response,
                               rtol=2e-12,atol=2e-12)
    np.testing.assert_allclose(reused.response_jacobian,
                               analytic.response_jacobian,
                               rtol=2e-11,atol=2e-11)
    assert reused.state_iterations==123


def test_native_tet_ima_directional_derivative_matches_geometry_regression():
    """The matrix-free dG contraction differentiates the full IMA fold."""
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.topology_optimization import production_tet_charge_gram_derivatives
    from radia.vim._vim import _charge_basis,build_charge_gram

    gradient=np.array([[.07,-.03,.02],[0.,.05,0.],[0.,0.,.04]])
    offset=np.array([.01,0.,0.])
    epsilon=1e-6
    image_masks=(2,4,6)
    image_signs=(1.,-1.,-1.)

    def build(step):
        def mapping(x,y,z):
            point=np.array([x,y,z])
            return tuple(point+step*epsilon*(gradient@point+offset))
        mesh=MakeStructured3DMesh(
            hexes=False,nx=1,ny=1,nz=1,mapping=mapping)
        fes=ng.HDiv(mesh,order=1)
        with ng.TaskManager():
            charge_basis=_charge_basis(fes,4)
            _,gram,_=build_charge_gram(
                fes,eps=1e-10,leafsize=256,eta=2.0,
                image_masks=image_masks,image_signs=image_signs)
        count=len(charge_basis["host"])
        dense=np.array([[gram.entry(i,j) for j in range(count)]
                        for i in range(count)])
        return charge_basis,gram,dense

    charge_basis,gram,_=build(0)
    _,_,plus=build(1)
    _,_,minus=build(-1)
    cells=np.asarray(charge_basis["vV"])
    faces=np.asarray(charge_basis["bV"])
    derivative,_=production_tet_charge_gram_derivatives(
        gram,(cells@gradient.T+offset)[None,...],
        (faces@gradient.T+offset)[None,...],charge_basis["B"])
    finite_difference=(plus-minus)/(2*epsilon)

    np.testing.assert_allclose(derivative[0],derivative[0].T,rtol=0,atol=0)
    np.testing.assert_allclose(
        derivative[0],finite_difference,rtol=8e-7,atol=2e-10)


def test_native_production_hex_face_self_block_derivative_invariants():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis_hex,build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1,
        mapping=lambda x,y,z:(x,y,z+.12*x*y))
    fes=ng.HDiv(mesh,order=1)
    with ng.TaskManager():
        cb=_charge_basis_hex(fes,cob_quad=3)
        _,gram,_=build_charge_gram(fes,eps=1e-10,leafsize=16,eta=2.0)
    face_nodes=np.asarray(cb["face_nodes"]).reshape(-1,9,3)
    host=4  # warped z-face: production value and derivative share radial-Duffy
    translation=gram.hex_face_self_block_directional_derivative(host,np.ones((9,3)))
    scaling=gram.hex_face_self_block_directional_derivative(host,face_nodes[host])
    offset=8+4*host
    value=np.array([[gram.entry(offset+i,offset+j) for j in range(4)] for i in range(4)])
    np.testing.assert_allclose(translation,0,atol=5e-17)
    np.testing.assert_allclose(scaling,scaling.T,rtol=0,atol=0)
    np.testing.assert_allclose(scaling,-value,rtol=3e-13,atol=3e-15)


def test_native_affine_hex_face_self_block_derivative_matches_production_value_and_fd():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis_hex,build_charge_gram

    velocity_gradient=np.array([[.13,-.07,.04],[.02,.09,-.05],[-.03,.06,.11]])
    velocity_offset=np.array([.01,-.02,.03])
    epsilon=2e-6

    def build(scale):
        mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1,
            mapping=lambda x,y,z:tuple(np.array([x,y,z])+scale*epsilon*(
                velocity_gradient@np.array([x,y,z])+velocity_offset)))
        fes=ng.HDiv(mesh,order=1)
        with ng.TaskManager():
            cb=_charge_basis_hex(fes,cob_quad=3)
            # Keep +/-epsilon on the same dense numerical path; ACA factor
            # recompression is not a differentiable finite-difference oracle.
            _,gram,_=build_charge_gram(fes,eps=1e-10,leafsize=256,eta=2.0)
        return cb,gram

    cb,gram=build(0); _,plus=build(1); _,minus=build(-1)
    face_nodes=np.asarray(cb["face_nodes"]).reshape(-1,9,3)
    for host,nodes in enumerate(face_nodes):
        offset=8+4*host
        value=np.array([[gram.entry(offset+i,offset+j) for j in range(4)] for i in range(4)])
        scaling=gram.hex_face_self_block_directional_derivative(host,nodes)
        translation=gram.hex_face_self_block_directional_derivative(host,np.ones((9,3)))
        velocity=nodes@velocity_gradient.T+velocity_offset
        derivative=gram.hex_face_self_block_directional_derivative(host,velocity)
        value_plus=np.array([[plus.entry(offset+i,offset+j) for j in range(4)] for i in range(4)])
        value_minus=np.array([[minus.entry(offset+i,offset+j) for j in range(4)] for i in range(4)])
        validation_difference=(value_plus-value_minus)/(2*epsilon)
        np.testing.assert_allclose(scaling,-value,rtol=2e-13,atol=3e-15)
        np.testing.assert_allclose(translation,0,atol=5e-17)
        np.testing.assert_allclose(derivative,derivative.T,rtol=0,atol=0)
        np.testing.assert_allclose(derivative,validation_difference,rtol=2e-8,atol=1e-10)


def test_native_complete_hex_charge_gram_directional_derivative():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis_hex,build_charge_gram
    gradient=np.array([[.13,-.07,.04],[.02,.09,-.05],[-.03,.06,.11]])
    offset=np.array([.01,-.02,.03]); epsilon=2e-6
    # Fully generic Q2 warp keeps every closest-corner/near decision away
    # from symmetry ties, so the validation difference follows one fixed
    # production quadrature branch on both sides.
    def base(x,y,z): return np.array([
        x+.071*y*z+.013*x*y,
        y+.037*x*z+.009*y*z,
        z+.053*x*y+.011*x*z,
    ])
    def build(scale):
        mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1,
            mapping=lambda x,y,z:tuple(base(x,y,z)+scale*epsilon*(gradient@base(x,y,z)+offset)))
        fes=ng.HDiv(mesh,order=1)
        with ng.TaskManager():
            cb=_charge_basis_hex(fes,cob_quad=3)
            # Keep this regression on a single dense leaf.  Differentiating
            # ACA recompression factors is outside the ChargeGram kernel API.
            _,gram,_=build_charge_gram(fes,eps=1e-10,leafsize=256,eta=2.0)
        return cb,gram
    cb,gram=build(0); _,plus=build(1); _,minus=build(-1)
    cells=np.asarray(cb["cell_nodes"]).reshape(-1,27,3)
    faces=np.asarray(cb["face_nodes"]).reshape(-1,9,3)
    zeros_c=np.zeros_like(cells); zeros_f=np.zeros_like(faces)
    translation=gram.hex_charge_gram_directional_derivative(
        zeros_c+np.array([.2,-.1,.3]),zeros_f+np.array([.2,-.1,.3]))
    scaling=gram.hex_charge_gram_directional_derivative(cells,faces)
    derivative=gram.hex_charge_gram_directional_derivative(cells@gradient.T+offset,faces@gradient.T+offset)
    left=np.linspace(-.3,.8,derivative.shape[0])
    right=np.linspace(.7,-.2,derivative.shape[0])
    contractions=gram.directional_derivative_contractions(
        "hex",
        np.ascontiguousarray(np.stack([cells, cells@gradient.T+offset])),
        np.ascontiguousarray(np.stack([faces, faces@gradient.T+offset])),
        np.ascontiguousarray(left),np.ascontiguousarray(right))
    left_many=np.ascontiguousarray(np.stack([left,2.0*left-right]))
    contractions_many=gram.directional_derivative_contractions_many(
        "hex",
        np.ascontiguousarray(np.stack([cells,cells@gradient.T+offset])),
        np.ascontiguousarray(np.stack([faces,faces@gradient.T+offset])),
        left_many,np.ascontiguousarray(right))
    derivative_operator=gram.directional_derivative_operator(
        "hex",cells@gradient.T+offset,faces@gradient.T+offset,
        eps=1e-12,leaf=256,eta=2.0)
    n=derivative.shape[0]
    value=np.array([[gram.entry(i,j) for j in range(n)] for i in range(n)])
    value_plus=np.array([[plus.entry(i,j) for j in range(n)] for i in range(n)])
    value_minus=np.array([[minus.entry(i,j) for j in range(n)] for i in range(n)])
    np.testing.assert_allclose(translation,0,atol=2e-15)
    np.testing.assert_allclose(scaling,-value,rtol=3e-10,atol=3e-13)
    np.testing.assert_allclose(derivative,derivative.T,rtol=0,atol=0)
    np.testing.assert_allclose(contractions,
        [left@scaling@right,left@derivative@right],rtol=3e-12,atol=3e-13)
    np.testing.assert_allclose(contractions_many,
        [[left@scaling@right,left@derivative@right],
         [(2*left-right)@scaling@right,
          (2*left-right)@derivative@right]],rtol=3e-12,atol=3e-13)
    operator_dense=np.array([[derivative_operator.entry(i,j) for j in range(n)] for i in range(n)])
    np.testing.assert_allclose(operator_dense,derivative,rtol=2e-13,atol=2e-15)
    probe=np.linspace(-.7,.9,n)
    np.testing.assert_allclose(derivative_operator.matvec_sym(probe),derivative@probe,
        rtol=2e-12,atol=2e-13)
    right=probe[::-1].copy()
    contractions=gram.directional_derivative_contractions("hex",
        np.ascontiguousarray(np.stack([cells@gradient.T+offset,cells])),
        np.ascontiguousarray(np.stack([faces@gradient.T+offset,faces])),probe,right)
    np.testing.assert_allclose(contractions,
        [probe@derivative@right,probe@scaling@right],rtol=3e-13,atol=3e-13)
    np.testing.assert_allclose(derivative,(value_plus-value_minus)/(2*epsilon),rtol=3e-7,atol=2e-10)


def test_native_hex_cluster_leaf_contraction_matches_analytic_dense_derivative():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis_hex,build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=2,ny=2,nz=1,
        mapping=lambda x,y,z:(x+.031*y*z,y+.019*x*z,z+.027*x*y))
    fes=ng.HDiv(mesh,order=1)
    with ng.TaskManager():
        cb=_charge_basis_hex(fes,cob_quad=3)
        _,gram,_=build_charge_gram(fes,eps=1e-8,leafsize=8,eta=1.0)
    assert gram.stats()["n_lowrank"]>0
    cells=np.asarray(cb["cell_nodes"]).reshape(-1,27,3)
    faces=np.asarray(cb["face_nodes"]).reshape(-1,9,3)
    def velocity(x):
        weight=np.maximum(0.,1.-2.*np.abs(x[...,0]-.5))
        return weight[...,None]*np.array([.021,-.017,.013])
    vc,vf=velocity(cells),velocity(faces)
    dense=np.asarray(gram.hex_charge_gram_directional_derivative(vc,vf))
    left=np.linspace(-.4,.7,dense.shape[0]);right=np.cos(np.arange(dense.shape[0]))
    observed=gram.directional_derivative_contractions("hex",vc[None],vf[None],left,right)[0]
    left_many=np.ascontiguousarray(np.stack((left,2*left-right)))
    observed_many=np.asarray(gram.directional_derivative_contractions_many(
        "hex",np.ascontiguousarray(np.stack((vc,2*vc))),
        np.ascontiguousarray(np.stack((vf,2*vf))),left_many,right))
    np.testing.assert_allclose(observed,left@dense@right,rtol=2e-6,atol=2e-9)
    np.testing.assert_allclose(observed_many,
        [[left@dense@right,2*left@dense@right],
         [(2*left-right)@dense@right,2*(2*left-right)@dense@right]],
        rtol=2e-6,atol=2e-9)


def test_vim_linearization_matches_analytic_two_cell_system():
    A=np.array([[3.0,-1.0],[-1.0,2.0]])
    b=np.array([1.0,0.5]); C=np.array([[1.0,2.0]])
    dA=np.array([[[0.4,0.0],[0.0,0.0]],[[0.0,0.0],[0.0,0.3]]])
    result=linearize_vim_system(A,b,C,dA)
    epsilon=1e-7
    for cell in range(2):
        shifted=np.linalg.solve(A+epsilon*dA[cell],b)
        observed=(C@shifted-result.response)/epsilon
        np.testing.assert_allclose(observed,result.response_jacobian[:,cell],rtol=2e-6,atol=2e-8)


def test_laplace_pair_shape_derivative_matches_validation_difference():
    points=np.array([[0.,0.,0.],[1.,0.,0.],[0.,2.,0.]])
    weights=np.array([1.,1.5,.75]); velocity=np.array([[[.1,0,0],[0,.2,0],[-.1,0,0]]])
    rel=np.array([[.03,-.02,.01]])
    gram,derivative=linearize_laplace_pair_gram(points,weights,velocity,rel)
    eps=1e-7
    moved,_=linearize_laplace_pair_gram(points+eps*velocity[0],weights*(1+eps*rel[0]),velocity*0)
    np.testing.assert_allclose((moved-gram)/eps,derivative[0],rtol=2e-6,atol=2e-9)


def test_vim_operator_product_rule_matches_validation_difference():
    M=np.array([[2.,.2],[.2,1.5]]); B=np.array([[1.,-1.],[.5,.25]])
    G=np.array([[.8,.1],[.1,.6]]); h=np.array([2.,-1.]); inv_chi=.1
    dM=np.array([[[.1,.02],[.02,-.03]]]); dB=np.array([[[.03,0],[-.01,.02]]]); dG=np.array([[[.04,.01],[.01,-.02]]])
    lin=linearize_vim_operator(M,B,G,h,inv_chi=inv_chi,dmass=dM,dcharge_map=dB,dcharge_gram=dG)
    eps=1e-7; shifted=linearize_vim_operator(M+eps*dM[0],B+eps*dB[0],G+eps*dG[0],h,
        inv_chi=inv_chi,dmass=np.zeros_like(dM),dcharge_map=np.zeros_like(dB),dcharge_gram=np.zeros_like(dG))
    np.testing.assert_allclose((shifted.matrix-lin.matrix)/eps,lin.matrix_jacobian[0],rtol=2e-6,atol=2e-9)
    np.testing.assert_allclose((shifted.rhs-lin.rhs)/eps,lin.rhs_jacobian[0],rtol=2e-6,atol=2e-9)


def test_matrix_free_vim_directional_action_matches_dense_product_rule():
    from radia.topology_optimization import linearize_vim_operator_matrix_free
    class SymOperator:
        def __init__(self,matrix): self.matrix=np.asarray(matrix)
        def matvec_sym(self,x): return self.matrix@x
    M=np.array([[2.,.2],[.2,1.3]])
    B=np.array([[1.,.3],[-.2,.8],[.4,-.1]])
    G=np.array([[1.2,.1,.05],[.1,.9,-.03],[.05,-.03,.7]])
    dM=np.array([[[.1,.02],[.02,-.03]]])
    dB=np.array([[[.03,0],[-.01,.02],[.04,-.02]]])
    dG=np.array([[[.04,.01,0],[.01,-.02,.03],[0,.03,.01]]])
    op=linearize_vim_operator_matrix_free(M,B,SymOperator(G),inv_chi=2.5,
        dmass=dM,dcharge_map=dB,dcharge_gram=(SymOperator(dG[0]),))
    x=np.array([.7,-.4])
    A=2.5*M+B.T@G@B
    dA=2.5*dM[0]+dB[0].T@G@B+B.T@dG[0]@B+B.T@G@dB[0]
    np.testing.assert_allclose(op.matvec(x),A@x,rtol=2e-15,atol=2e-15)
    np.testing.assert_allclose(op.directional_matvec(0,x),dA@x,rtol=2e-15,atol=2e-15)
    assert op.as_scipy_linear_operator(0).shape==(2,2)


def test_lp_update_obeys_volume_and_move_limit():
    result=solve_lp_update([0.5,0.5,0.5],[-3.0,-1.0,2.0],[1.0,1.0,1.0],1.5,move_limit=0.2)
    assert np.all(np.abs(result.delta)<=0.2+1e-12)
    assert np.sum(result.density)<=1.5+1e-12
    assert result.density[0]>=result.density[1]>=result.density[2]


def test_cubit_density_journal_has_deterministic_blocks(tmp_path):
    path=tmp_path/"density.jou"
    info=write_cubit_density_journal(path,[11,12,13,14],[0.9,0.1,0.7,0.2])
    text=path.read_text(encoding="utf-8")
    assert info["solid_count"]==2 and info["void_count"]==2
    assert "group 'radia_topopt_solid' add hex 11 13" in text
    assert "block 1002 hex in group 'radia_topopt_void'" in text


def test_sequential_vim_lp_reaches_volume_constrained_material_layout():
    def linearize(density):
        A=np.eye(2); b=np.array([1.0,1.0]); C=np.eye(2)
        result=linearize_vim_system(A,b,C,np.zeros((3,2,2)))
        return VIMLinearization(
            result.state,
            np.array([density[0]-density[2],density[1]]),
            result.state_jacobian,
            np.array([[-1.0,0.0,1.0],[0.0,1.0,0.0]]),
        )
    result=optimize_vim_lp([0.5]*3,[1.0]*3,0.5,linearize,objective_weights=[1.0,0.0],move_limit=0.25,max_iterations=5)
    assert np.sum(result.density)<=1.5+1e-12
    assert result.density[0]>=result.density[2]
