import numpy as np
import pytest
from types import SimpleNamespace

from radia.sheet_metal_optimization import (apply_ngsolve_mesh_route,
    backtrack_ngsolve_deformation, backtrack_ngsolve_target_deformation,
    combine_deformation_modes, elastic_normal_deformation_modes,
    optimize_topology_preserving_shape, relative_gettrafo_displacements,
    route_mesh_update, sample_trafo_quality,
    sample_affine_gettrafo_cells, solve_sheet_metal_lp, local_trust_region)
from radia.sheet_metal_optimization import (
    ShapeModelEvaluation, TopologyPreservingShapeState)
from radia.topology_optimization import ShapeLinearization


def test_real_elastic_normal_modes_are_gettrafo_ready_without_gray_material():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh

    mesh = MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=2)
    with ng.TaskManager():
        modes = elastic_normal_deformation_modes(
            mesh, [ng.CF(.01)], movable_boundaries="top",
            fixed_boundaries="bottom")
    assert len(modes) == 1
    np.testing.assert_allclose(
        modes[0](mesh(.5, .5, 1.0)), [0, 0, .01], atol=2e-14)
    np.testing.assert_allclose(
        modes[0](mesh(.5, .5, 0.0)), [0, 0, 0], atol=2e-14)
    deformation = combine_deformation_modes(modes, [2.0])
    relative = relative_gettrafo_displacements(mesh, deformation)
    assert relative.shape == (mesh.ne,)
    assert np.all(relative > 0)
    accepted = backtrack_ngsolve_deformation(mesh, deformation, relative)
    assert accepted.accepted and mesh.deformation is not None
    mesh.UnsetDeformation()


def test_shape_deformation_inputs_reject_nonfinite_values():
    with pytest.raises(ValueError, match="coefficients must be finite"):
        combine_deformation_modes([object()], [np.nan])
    with pytest.raises(ValueError, match="moduli"):
        elastic_normal_deformation_modes(
            object(), [1.0], movable_boundaries="top", fixed_boundaries="",
            shear_modulus=np.nan)


def test_topology_preserving_shape_driver_resolves_every_trial(monkeypatch):
    import radia.sheet_metal_optimization as sm

    class Mesh:
        deformation = None

        def SetDeformation(self, value):
            self.deformation = value

        def UnsetDeformation(self):
            self.deformation = None

    mesh = Mesh()
    calls = {"solve": 0}
    monkeypatch.setattr(
        sm, "sample_trafo_quality",
        lambda mesh, **kwargs: (np.ones(2), np.ones(2)))
    monkeypatch.setattr(
        sm, "relative_gettrafo_displacements",
        lambda mesh, deformation: np.full(2, .02))
    initial = TopologyPreservingShapeState(
        mesh, {"q": 0.0}, np.array([0.]), np.array([0.]),
        ShapeModelEvaluation(1.0, np.empty(0)))

    def linearize(state):
        q = float(state.parameters[0])
        return ShapeLinearization(
            (q - 1.)**2, np.array([2 * (q - 1.)]), np.empty(0),
            np.zeros((0, 1)), np.empty(0), np.empty(0))

    def rebuild(active_mesh, parameters, route):
        calls["solve"] += 1
        return {"q": float(parameters[0]), "route": route}

    result = optimize_topology_preserving_shape(
        initial, linearize_step=linearize,
        deformation_factory=lambda mesh, reference, candidate: float(
            candidate[0] - reference[0]),
        rebuild_model=rebuild,
        evaluate_model=lambda model: ShapeModelEvaluation(
            (model["q"] - 1.)**2, np.empty(0)),
        move_limit=.25, parameter_bounds=([-1.], [1.]), max_iterations=4)
    np.testing.assert_allclose(result.state.parameters, [1.0], atol=1e-12)
    assert len(result.history) == 4 and calls["solve"] == 4
    assert all(item.route == "ngsolve_deform" for item in result.history)
    assert all(item.nonlinear_resolves == 1 for item in result.history)
    objectives = [item.objective_after for item in result.history]
    assert objectives == sorted(objectives, reverse=True)
    assert objectives[-1] == 0.0


def test_topology_preserving_shape_clears_deformation_after_solver_error(
        monkeypatch):
    import radia.sheet_metal_optimization as sm

    class Mesh:
        deformation = None

        def SetDeformation(self, value):
            self.deformation = value

        def UnsetDeformation(self):
            self.deformation = None

    mesh = Mesh()
    monkeypatch.setattr(
        sm, "sample_trafo_quality",
        lambda mesh, **kwargs: (np.ones(1), np.ones(1)))
    monkeypatch.setattr(
        sm, "relative_gettrafo_displacements",
        lambda mesh, deformation: np.full(1, 0.01))
    initial = TopologyPreservingShapeState(
        mesh, object(), np.array([0.0]), np.array([0.0]),
        ShapeModelEvaluation(1.0, np.empty(0)))
    linearization = ShapeLinearization(
        1.0, np.array([-1.0]), np.empty(0), np.zeros((0, 1)),
        np.empty(0), np.empty(0))

    with pytest.raises(RuntimeError, match="trial solve failed"):
        optimize_topology_preserving_shape(
            initial, linearize_step=lambda state: linearization,
            deformation_factory=lambda mesh, reference, candidate: object(),
            rebuild_model=lambda mesh, parameters, route: (_ for _ in ()).throw(
                RuntimeError("trial solve failed")),
            evaluate_model=lambda model: None, move_limit=0.1,
            max_iterations=1)
    assert mesh.deformation is None


def test_hex_topology_driver_prefers_deformation_then_uses_cubit(monkeypatch,tmp_path):
    import radia.sheet_metal_optimization as sm
    mesh0=SimpleNamespace(name="initial")
    state=sm.HexSheetTopologyState(mesh0,{"iteration":-1},np.zeros(1),
        np.ones(1)*.002,np.ones(1),10.)
    calls={"iteration":0,"cubit":0}
    def linearize(current):
        k=calls["iteration"]; calls["iteration"]+=1
        rho=np.array([1. if k==0 else .2])
        update=sm.SheetMetalUpdate(np.array([.01*(k+1)]),np.array([.002]),rho,
            np.zeros(3),"ok")
        return SimpleNamespace(update=update)
    def accept(mesh,factory,current,target,relative,**kwargs):
        topology=bool(kwargs["topology_changed"])
        decision=sm.MeshUpdateDecision("cubit_rebuild" if topology else "ngsolve_deform",
            np.empty(0,dtype=int),("test",),.8,2.,float(np.max(relative)))
        return sm.DeformationAcceptance(True,1.,decision,1),object()
    monkeypatch.setattr(sm,"backtrack_ngsolve_target_deformation",accept)
    class Backend:
        def rebuild(self,request):
            calls["cubit"]+=1
            assert request.journal_path.parent==tmp_path
            return SimpleNamespace(name="cubit")
    result=sm.optimize_hex_sheet_topology(state,linearize_step=linearize,
        deformation_factory=lambda mesh,u:object(),
        rebuild_model=lambda mesh,u,t,r,route:{"mesh":mesh,"route":route,"k":calls["iteration"]},
        evaluate_objective=lambda model:10.-model["k"],element_sizes=[.1],
        cubit_backend=Backend(),cubit_work_directory=tmp_path,max_iterations=2,
        design_tolerance=0.,objective_tolerance=0.)
    assert [x.route for x in result.history]==["ngsolve_deform","cubit_rebuild"]
    assert calls["cubit"]==1 and result.state.mesh.name=="cubit"
    np.testing.assert_allclose(result.state.activation,[.2])


def test_hex_topology_driver_batches_sparse_cubit_changes(monkeypatch,tmp_path):
    import radia.sheet_metal_optimization as sm
    n=10
    state=sm.HexSheetTopologyState(SimpleNamespace(name="initial"),{},np.zeros(n),
        np.full(n,.002),np.ones(n),10.)
    calls={"iteration":0,"cubit":0}
    def linearize(current):
        calls["iteration"]+=1
        rho=np.ones(n);rho[0]=.2
        return SimpleNamespace(update=sm.SheetMetalUpdate(
            np.zeros(n),np.full(n,.002),rho,np.zeros(3),"ok"))
    def accept(mesh,factory,current,target,relative,**kwargs):
        route="cubit_rebuild" if kwargs["topology_changed"] else "ngsolve_deform"
        return sm.DeformationAcceptance(True,1.,sm.MeshUpdateDecision(
            route,np.empty(0,dtype=int),("test",),.8,2.,0.),1),object()
    monkeypatch.setattr(sm,"backtrack_ngsolve_target_deformation",accept)
    class Backend:
        def rebuild(self,request):
            calls["cubit"]+=1
            return SimpleNamespace(name="cubit")
    result=sm.optimize_hex_sheet_topology(state,linearize_step=linearize,
        deformation_factory=lambda mesh,u:object(),cubit_backend=Backend(),
        cubit_work_directory=tmp_path,element_sizes=np.ones(n),max_iterations=3,
        cubit_batch_interval=3,cubit_batch_fraction=.5,
        rebuild_model=lambda mesh,u,t,r,route:{"route":route},
        evaluate_objective=lambda model:10.-calls["iteration"],
        design_tolerance=0.,objective_tolerance=0.)
    assert [item.route for item in result.history]==[
        "ngsolve_deform","ngsolve_deform","cubit_rebuild"]
    assert [item.pending_topology_changes for item in result.history]==[1,1,1]
    assert calls["cubit"]==1


def test_cubit_pending_change_cancels_inside_hysteresis_band(monkeypatch,tmp_path):
    import radia.sheet_metal_optimization as sm
    calls={"iteration":0,"cubit":0}
    state=sm.HexSheetTopologyState(SimpleNamespace(),{},np.zeros(2),
        np.full(2,.002),np.ones(2),1.)
    def linearize(current):
        rho=np.ones(2); rho[0]=.2 if calls["iteration"]==0 else .5
        calls["iteration"]+=1
        return SimpleNamespace(update=sm.SheetMetalUpdate(
            np.zeros(2),np.full(2,.002),rho,np.zeros(6),"ok"))
    def accept(mesh,factory,current,target,relative,**kwargs):
        decision=sm.MeshUpdateDecision("ngsolve_deform",np.empty(0,int),
            ("test",),1.,1.,0.)
        return sm.DeformationAcceptance(True,1.,decision,1),object()
    monkeypatch.setattr(sm,"backtrack_ngsolve_target_deformation",accept)
    class Backend:
        def rebuild(self,request):
            calls["cubit"]+=1
            return SimpleNamespace()
    result=sm.optimize_hex_sheet_topology(state,linearize_step=linearize,
        deformation_factory=lambda mesh,u:object(),rebuild_model=lambda *args:{},
        evaluate_objective=lambda model:1.,element_sizes=np.ones(2),
        cubit_backend=Backend(),cubit_work_directory=tmp_path,max_iterations=2,
        cubit_batch_fraction=1.,cubit_batch_interval=5,
        design_tolerance=0.,objective_tolerance=0.)
    assert [item.pending_topology_changes for item in result.history]==[1,0]
    assert calls["cubit"]==0


def test_sheet_metal_lp_obeys_moves_thickness_curvature_and_volume():
    n=3; u=np.zeros(n); t=np.ones(n); rho=np.ones(n); area=np.ones(n)
    L=np.array([[1.,-2.,1.]])
    trust=local_trust_region([1,2,1],fraction=.1)
    update=solve_sheet_metal_lp(u,t,rho,np.r_[-np.ones(n),np.ones(n),np.ones(n)],area,
        volume_max=2.7,displacement_move=trust,thickness_move=.2,activation_move=.1,
        thickness_bounds=(.8,1.2),laplacian=L,curvature_limit=.05)
    assert np.all(np.abs(update.normal_displacement) <= trust+1e-12)
    assert np.max(np.abs(L@update.normal_displacement)) <= .05+1e-12
    assert np.all((update.thickness>=.8)&(update.thickness<=1.2))
    linearized_volume=np.sum(area*(t*rho + rho*(update.thickness-t) + t*(update.activation-rho)))
    assert linearized_volume <= 2.7+1e-10


def test_mesh_update_routes_deform_refine_and_rebuild():
    deform=route_mesh_update([.9,.8],[2,3],[.02,.04])
    refine=route_mesh_update([.9,.3],[2,10],[.02,.3])
    rebuild=route_mesh_update([.9,.8],[2,3],[.02,.6])
    topology=route_mesh_update([.9],[2],[.02],topology_changed=True)
    assert deform.route=="ngsolve_deform"
    assert refine.route=="ngsolve_refine" and refine.refine_elements.tolist()==[1]
    assert rebuild.route=="cubit_rebuild"
    assert topology.route=="cubit_rebuild"


def test_ngsolve_route_executor_keeps_cubit_rebuild_external():
    class Mesh:
        deformation=None
        def SetDeformation(self,value): self.deformation=value
    mesh=Mesh(); decision=route_mesh_update([.9],[2],[.02])
    assert apply_ngsolve_mesh_route(mesh,"gf",decision)["action"]=="deformation_set"
    assert mesh.deformation=="gf"
    rebuild=route_mesh_update([.9],[2],[.02],topology_changed=True)
    assert apply_ngsolve_mesh_route(mesh,"gf",rebuild)["action"]=="cubit_rebuild_required"


def test_real_ngsolve_trafo_quality_and_backtracking():
    import ngsolve as ng
    from netgen.geom2d import unit_square
    mesh=ng.Mesh(unit_square.GenerateMesh(maxh=.35))
    fes=ng.VectorH1(mesh,order=1); deformation=ng.GridFunction(fes)
    deformation.Set(ng.CF((0.08*ng.x,0.0)))
    base_det,base_cond=sample_trafo_quality(mesh)
    assert base_det.size==mesh.ne and np.all(np.isfinite(base_cond))
    accepted=backtrack_ngsolve_deformation(mesh,deformation,np.full(mesh.ne,.08))
    assert accepted.accepted and accepted.scale==1.0
    assert accepted.decision.route=="ngsolve_deform"
    mesh.UnsetDeformation()


def test_real_gettrafo_displacement_reaches_affine_self_kernel_contract():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.topology_optimization import affine_cell_self_energy_shape_derivative
    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1,
        mapping=lambda x,y,z:(2*x,3*y,4*z))
    fes=ng.VectorH1(mesh,order=1); deformation=ng.GridFunction(fes)
    deformation.Set(ng.CF((.01*ng.x,.02*ng.y,.03*ng.z)))
    sampled=sample_affine_gettrafo_cells(mesh,[deformation])
    assert sampled.cell_types==("hex",)
    assert sampled.nodes[0].shape==(8,3)
    np.testing.assert_allclose(sampled.centroids,[[1,1.5,2]],atol=1e-14)
    value,jac=affine_cell_self_energy_shape_derivative(
        sampled.cell_types[0],sampled.nodes[0],sampled.node_displacements[0])
    assert value>0 and jac.shape==(1,) and np.isfinite(jac[0])


def test_real_hex_absolute_target_deformation_path():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1)
    space=ng.VectorH1(mesh,order=1)
    def factory(active_mesh,normal):
        field=ng.GridFunction(space)
        field.Set(ng.CF((float(normal[0])*ng.x,0,0)))
        return field
    accepted,deformation=backtrack_ngsolve_target_deformation(
        mesh,factory,np.array([0.]),np.array([.04]),np.array([.04]))
    assert accepted.accepted and accepted.decision.route=="ngsolve_deform"
    assert deformation is not None and mesh.deformation is not None
    mesh.UnsetDeformation()
