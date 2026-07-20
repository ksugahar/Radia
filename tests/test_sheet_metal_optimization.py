import numpy as np
from types import SimpleNamespace

from radia.sheet_metal_optimization import (apply_ngsolve_mesh_route,
    backtrack_ngsolve_deformation, backtrack_ngsolve_target_deformation,
    route_mesh_update, sample_trafo_quality,
    sample_affine_gettrafo_cells, solve_sheet_metal_lp, local_trust_region)


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
