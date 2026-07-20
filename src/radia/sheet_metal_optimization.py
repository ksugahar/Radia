"""Manufacturing-aware sheet-metal updates for Radia-VIM topology workflows.

NGSolve owns the deformation and element transformations.  This module builds
the bounded LP update and decides when deformation is still safe, when local
refinement is required, and when a topology-changing Cubit rebuild is required.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import subprocess

import numpy as np


@dataclass(frozen=True)
class SheetMetalUpdate:
    normal_displacement: np.ndarray
    thickness: np.ndarray
    activation: np.ndarray
    delta: np.ndarray
    status: str


@dataclass(frozen=True)
class MeshUpdateDecision:
    route: str
    refine_elements: np.ndarray
    reasons: tuple[str, ...]
    minimum_jacobian: float
    maximum_condition: float
    maximum_relative_displacement: float


@dataclass(frozen=True)
class DeformationAcceptance:
    accepted: bool
    scale: float
    decision: MeshUpdateDecision
    trials: int


@dataclass(frozen=True)
class AffineGetTrafoCells:
    cell_types: tuple[str, ...]
    nodes: tuple[np.ndarray, ...]
    node_displacements: tuple[np.ndarray, ...]
    centroids: np.ndarray
    centroid_displacements: np.ndarray


@dataclass(frozen=True)
class HexSheetTopologyState:
    mesh: object
    model: object
    normal_displacement: np.ndarray
    thickness: np.ndarray
    activation: np.ndarray
    objective: float


@dataclass(frozen=True)
class HexSheetTopologyIteration:
    iteration: int
    objective_before: float
    objective_after: float
    accepted_scale: float
    route: str
    minimum_jacobian: float
    maximum_condition: float
    topology_changed: bool


@dataclass(frozen=True)
class HexSheetTopologyResult:
    state: HexSheetTopologyState
    history: tuple[HexSheetTopologyIteration, ...]
    converged: bool


@dataclass(frozen=True)
class CubitHexRemeshRequest:
    iteration: int
    normal_displacement: np.ndarray
    thickness: np.ndarray
    activation: np.ndarray
    journal_path: Path
    mesh_path: Path


@dataclass(frozen=True)
class CubitHexRemeshBackend:
    """Batch Cubit adapter; CAD/journal policy remains application-owned."""
    executable: str
    work_directory: Path
    write_journal: object
    load_mesh: object
    arguments: tuple[str, ...] = ("-batch",)

    def rebuild(self, request: CubitHexRemeshRequest):
        work=Path(self.work_directory); work.mkdir(parents=True,exist_ok=True)
        self.write_journal(request)
        journal=Path(request.journal_path); output=Path(request.mesh_path)
        if not journal.is_file(): raise RuntimeError(f"Cubit journal was not created: {journal}")
        completed=subprocess.run([self.executable,*self.arguments,str(journal)],
            cwd=work,text=True,capture_output=True,check=False)
        if completed.returncode!=0:
            tail="\n".join((completed.stdout+"\n"+completed.stderr).splitlines()[-30:])
            raise RuntimeError(f"Cubit remesh failed ({completed.returncode}):\n{tail}")
        if not output.is_file(): raise RuntimeError(f"Cubit did not create mesh: {output}")
        return self.load_mesh(output)


def sample_affine_gettrafo_cells(mesh, displacement_modes) -> AffineGetTrafoCells:
    """Sample affine TET/HEX/WEDGE geometry and deformation via ``GetTrafo``.

    NGSolve owns the reference-to-physical map and deformation evaluation.
    This adapter only packages mapped corner values for the C++ singular
    self-kernel; it does not reproduce element shapes, orientations, or Piola
    transforms in Python.
    """
    import ngsolve as ng
    modes=tuple(displacement_modes)
    references={
        "ET.TET":("tet",[(0,0,0),(1,0,0),(0,1,0),(0,0,1)]),
        "ET.HEX":("hex",[(0,0,0),(0,0,1),(0,1,1),(0,1,0),
                          (1,0,0),(1,0,1),(1,1,1),(1,1,0)]),
        "ET.PRISM":("wedge",[(0,0,0),(1,0,0),(0,1,0),
                              (0,0,1),(1,0,1),(0,1,1)]),
    }
    kinds=[]; all_nodes=[]; all_modes=[]
    for element in mesh.Elements(ng.VOL):
        key=str(element.type)
        if key not in references:
            raise NotImplementedError(f"affine GetTrafo sampling supports TET/HEX/WEDGE, got {key}")
        kind,ref=references[key]; trafo=mesh.GetTrafo(element)
        # Batched ``trafo(rule)`` is a mesh-point array on some NGSolve
        # versions, not physical coordinates.  Mapping each IntegrationPoint
        # is the stable public GetTrafo path.
        mapped_nodes=[]
        for p in ref:
            rule=ng.IntegrationRule([p],[1.0])  # keep owner alive while the IP is mapped
            ip=next(iter(rule))
            mapped_nodes.append(np.asarray(trafo(ip).point,dtype=float).copy())
        nodes=np.array(mapped_nodes)
        sampled=[]
        for field in modes:
            sampled.append(np.array([field(mesh(*point)) for point in nodes],dtype=float))
        kinds.append(kind); all_nodes.append(nodes)
        all_modes.append(np.asarray(sampled,dtype=float).reshape(len(modes),len(ref),3))
    centroids=np.array([x.mean(axis=0) for x in all_nodes])
    centroid_modes=np.stack([[v[k].mean(axis=0) for v in all_modes]
                             for k in range(len(modes))]) if modes else np.zeros((0,len(all_nodes),3))
    return AffineGetTrafoCells(tuple(kinds),tuple(all_nodes),tuple(all_modes),centroids,centroid_modes)


def solve_sheet_metal_lp(
    normal_displacement,
    thickness,
    activation,
    objective_gradient,
    cell_areas,
    *,
    volume_max,
    displacement_move,
    thickness_move,
    activation_move=0.2,
    thickness_bounds,
    laplacian=None,
    curvature_limit=None,
    A_ub=None,
    b_ub=None,
) -> SheetMetalUpdate:
    """Solve one LP over normal displacement, thickness, and material activation.

    Volume is linearized as ``sum(area * thickness * activation)`` at the
    current design.  ``laplacian`` constrains the updated normal displacement
    through ``abs(L*u) <= curvature_limit``.
    """
    from scipy.optimize import linprog

    u=np.asarray(normal_displacement,dtype=float).reshape(-1)
    t=np.asarray(thickness,dtype=float).reshape(-1)
    rho=np.asarray(activation,dtype=float).reshape(-1)
    area=np.asarray(cell_areas,dtype=float).reshape(-1)
    n=u.size
    if not (t.size==rho.size==area.size==n) or n==0: raise ValueError("sheet cell vectors must be non-empty and equal length")
    if np.any(area<=0) or np.any(t<=0) or np.any((rho<0)|(rho>1)): raise ValueError("invalid area, thickness, or activation")
    gradient=np.asarray(objective_gradient,dtype=float).reshape(-1)
    if gradient.size!=3*n: raise ValueError("objective_gradient must have 3*n entries ordered [normal, thickness, activation]")
    tmin,tmax=map(float,thickness_bounds)
    if not 0<tmin<=tmax: raise ValueError("invalid thickness_bounds")
    current=np.r_[u,t,rho]
    u_move=np.broadcast_to(np.asarray(displacement_move,dtype=float),u.shape)
    if np.any(u_move<=0): raise ValueError("displacement_move must be positive")
    lower=np.r_[u-u_move,np.maximum(tmin,t-float(thickness_move)),np.maximum(0,rho-float(activation_move))]
    upper=np.r_[u+u_move,np.minimum(tmax,t+float(thickness_move)),np.minimum(1,rho+float(activation_move))]
    # First-order product t*rho about the current design.
    volume_row=np.r_[np.zeros(n),area*rho,area*t]
    volume_rhs=float(volume_max)+float(np.sum(area*t*rho))
    rows=[volume_row]; limits=[volume_rhs]
    if laplacian is not None:
        L=np.atleast_2d(np.asarray(laplacian,dtype=float))
        if L.shape[1]!=n or curvature_limit is None: raise ValueError("laplacian requires n columns and curvature_limit")
        pad=np.zeros((L.shape[0],2*n)); rows.extend([np.c_[L,pad],np.c_[-L,pad]])
        limit=np.broadcast_to(np.asarray(curvature_limit,dtype=float),L.shape[0])
        limits.extend([limit,limit])
    if A_ub is not None:
        extra=np.atleast_2d(np.asarray(A_ub,dtype=float)); rhs=np.asarray(b_ub,dtype=float).reshape(-1)
        if extra.shape!=(rhs.size,3*n): raise ValueError("A_ub/b_ub shape mismatch")
        rows.append(extra); limits.append(rhs)
    result=linprog(gradient,A_ub=np.vstack(rows),b_ub=np.concatenate([np.atleast_1d(x) for x in limits]),
                   bounds=list(zip(lower,upper)),method="highs")
    if not result.success: raise RuntimeError(f"sheet-metal LP failed: {result.message}")
    x=result.x
    return SheetMetalUpdate(x[:n],x[n:2*n],x[2*n:],x-current,str(result.message))


def local_trust_region(element_sizes, *, fraction=0.1, minimum=None, maximum=None):
    """Return per-design-cell displacement bounds proportional to local size."""
    sizes=np.asarray(element_sizes,dtype=float).reshape(-1)
    if sizes.size==0 or np.any(sizes<=0) or not 0<fraction<1: raise ValueError("invalid element sizes or trust fraction")
    bounds=fraction*sizes
    if minimum is not None: bounds=np.maximum(bounds,float(minimum))
    if maximum is not None: bounds=np.minimum(bounds,float(maximum))
    return bounds


def route_mesh_update(jacobian_determinants, jacobian_conditions, relative_displacements, *,
                      refine_threshold=0.25, rebuild_threshold=0.5,
                      minimum_jacobian=0.2, maximum_condition=20.0,
                      topology_changed=False) -> MeshUpdateDecision:
    """Select NGSolve deformation, local refinement, or Cubit reconstruction."""
    det=np.asarray(jacobian_determinants,dtype=float).reshape(-1)
    cond=np.asarray(jacobian_conditions,dtype=float).reshape(-1)
    disp=np.asarray(relative_displacements,dtype=float).reshape(-1)
    if not (det.size==cond.size==disp.size) or det.size==0: raise ValueError("quality vectors must be non-empty and equal length")
    if not np.all(np.isfinite(np.r_[det,cond,disp])): raise ValueError("quality vectors must be finite")
    reasons=[]
    bad=(det<=minimum_jacobian)|(cond>=maximum_condition)|(disp>=refine_threshold)
    if topology_changed: reasons.append("material topology changed")
    if np.any(det<=0): reasons.append("inverted element")
    if np.any(disp>=rebuild_threshold): reasons.append("deformation exceeded rebuild threshold")
    if np.any(cond>=2*maximum_condition): reasons.append("severe Jacobian distortion")
    if reasons: route="cubit_rebuild"
    elif np.any(bad): route="ngsolve_refine"; reasons.append("local deformation quality threshold exceeded")
    else: route="ngsolve_deform"; reasons.append("deformation remains inside quality limits")
    return MeshUpdateDecision(route,np.flatnonzero(bad),tuple(reasons),float(det.min()),float(cond.max()),float(disp.max()))


def apply_ngsolve_mesh_route(mesh, deformation, decision: MeshUpdateDecision):
    """Apply an NGSolve-owned deformation/refinement route.

    A Cubit rebuild is deliberately returned as an external action because it
    changes CAD topology and must pass the Cubit Sculpt validation gates.
    """
    if decision.route=="ngsolve_deform":
        mesh.SetDeformation(deformation)
        return {"action":"deformation_set","element_count":0}
    if decision.route=="ngsolve_refine":
        from ngsolve import ElementId, VOL
        if getattr(mesh,"deformation",None) is not None: mesh.UnsetDeformation()
        for index in decision.refine_elements: mesh.SetRefinementFlag(ElementId(VOL,int(index)),True)
        mesh.Refine()
        return {"action":"locally_refined","element_count":int(decision.refine_elements.size)}
    if decision.route=="cubit_rebuild":
        return {"action":"cubit_rebuild_required","element_count":int(decision.refine_elements.size)}
    raise ValueError(f"unknown mesh route: {decision.route}")


def sample_trafo_quality(mesh, *, integration_order=2, reference_determinants=None):
    """Sample every volume-element Jacobian through NGSolve ``GetTrafo``."""
    import ngsolve as ng
    determinants=[]; conditions=[]
    for element in mesh.Elements(ng.VOL):
        trafo=mesh.GetTrafo(element)
        rule=ng.IntegrationRule(element.type,int(integration_order))
        local_det=[]; local_cond=[]
        for ip in rule:
            mapped=trafo(ip)
            J=np.asarray(mapped.jacobi,dtype=float)
            local_det.append(float(np.linalg.det(J)))
            local_cond.append(float(np.linalg.cond(J)))
        determinants.append(min(local_det)); conditions.append(max(local_cond))
    determinants=np.asarray(determinants); conditions=np.asarray(conditions)
    if reference_determinants is not None:
        reference=np.asarray(reference_determinants,dtype=float).reshape(-1)
        if reference.shape!=determinants.shape or np.any(reference==0): raise ValueError("reference determinant shape/value mismatch")
        determinants=determinants/reference
    return determinants,conditions


def backtrack_ngsolve_deformation(mesh, deformation, relative_displacements, *,
                                  minimum_scale=1/64, contraction=0.5,
                                  minimum_jacobian_ratio=0.2, maximum_condition=20.0,
                                  refine_threshold=0.25, rebuild_threshold=0.5,
                                  integration_order=2) -> DeformationAcceptance:
    """Trial a VectorH1 displacement and backtrack before a Trafo becomes invalid.

    The accepted deformation remains installed only for ``ngsolve_deform``.
    Refinement and Cubit routes leave the reference mesh undeformed so the
    caller can prolong/reproject the mesh-independent design field first.
    """
    import ngsolve as ng
    if not 0<contraction<1 or not 0<minimum_scale<=1: raise ValueError("invalid backtracking parameters")
    if getattr(mesh,"deformation",None) is not None: mesh.UnsetDeformation()
    reference,_=sample_trafo_quality(mesh,integration_order=integration_order)
    trial=ng.GridFunction(deformation.space)
    scale=1.0; trials=0; last=None
    while scale>=minimum_scale:
        trials+=1; trial.vec.data=scale*deformation.vec
        mesh.SetDeformation(trial)
        ratios,conditions=sample_trafo_quality(mesh,integration_order=integration_order,reference_determinants=reference)
        last=route_mesh_update(ratios,conditions,scale*np.asarray(relative_displacements),
            refine_threshold=refine_threshold,rebuild_threshold=rebuild_threshold,
            minimum_jacobian=minimum_jacobian_ratio,maximum_condition=maximum_condition)
        unsafe=np.any(ratios<=0) or np.any(ratios<=minimum_jacobian_ratio) or np.any(conditions>=2*maximum_condition)
        if not unsafe:
            if last.route!="ngsolve_deform": mesh.UnsetDeformation()
            return DeformationAcceptance(True,scale,last,trials)
        mesh.UnsetDeformation(); scale*=contraction
    assert last is not None
    return DeformationAcceptance(False,0.0,last,trials)


def backtrack_ngsolve_target_deformation(mesh, deformation_factory,
        current_normal, target_normal, relative_displacements, *,
        minimum_scale=1/64, contraction=0.5, minimum_jacobian_ratio=0.2,
        maximum_condition=20.0, refine_threshold=0.25,
        rebuild_threshold=0.5, integration_order=2,
        topology_changed=False):
    """Backtrack an absolute sheet shape while preserving the current shape."""
    if not 0<contraction<1 or not 0<minimum_scale<=1: raise ValueError("invalid backtracking parameters")
    current=np.asarray(current_normal,dtype=float); target=np.asarray(target_normal,dtype=float)
    if current.shape!=target.shape: raise ValueError("current/target normal shape mismatch")
    if getattr(mesh,"deformation",None) is not None: mesh.UnsetDeformation()
    reference,_=sample_trafo_quality(mesh,integration_order=integration_order)
    scale=1.0; trials=0; last=None; accepted_deformation=None
    while scale>=minimum_scale:
        trials+=1; candidate=current+scale*(target-current)
        deformation=deformation_factory(mesh,candidate)
        mesh.SetDeformation(deformation)
        ratios,conditions=sample_trafo_quality(mesh,integration_order=integration_order,
                                               reference_determinants=reference)
        last=route_mesh_update(ratios,conditions,scale*np.asarray(relative_displacements),
            refine_threshold=refine_threshold,rebuild_threshold=rebuild_threshold,
            minimum_jacobian=minimum_jacobian_ratio,maximum_condition=maximum_condition,
            topology_changed=topology_changed)
        unsafe=np.any(ratios<=minimum_jacobian_ratio) or np.any(conditions>=2*maximum_condition)
        if not unsafe:
            accepted_deformation=deformation
            if last.route!="ngsolve_deform": mesh.UnsetDeformation()
            return DeformationAcceptance(True,scale,last,trials),accepted_deformation
        mesh.UnsetDeformation(); scale*=contraction
    assert last is not None
    return DeformationAcceptance(False,0.0,last,trials),None


def optimize_hex_sheet_topology(initial_state: HexSheetTopologyState, *,
        linearize_step, deformation_factory, rebuild_model, evaluate_objective,
        element_sizes, cubit_backend: CubitHexRemeshBackend,
        cubit_work_directory, max_iterations=20, objective_tolerance=1e-4,
        design_tolerance=1e-3, activation_threshold=0.5,
        minimum_scale=1/64, contraction=0.5, minimum_jacobian_ratio=0.2,
        maximum_condition=20.0, refine_threshold=0.25,
        rebuild_threshold=0.5, integration_order=2,
        allow_ngsolve_refine=False):
    """Run HEX-sheet topology iterations with deform-first/Cubit fallback.

    ``linearize_step(state)`` returns a sheet LP result with ``.update``.
    ``rebuild_model(mesh,u,t,rho,route)`` rebuilds FESpace/ChargeGram data and
    returns the next model.  This keeps solver plumbing outside the driver and
    makes the remesher backend replaceable without changing optimization.
    """
    state=initial_state; history=[]; converged=False
    sizes=np.asarray(element_sizes,dtype=float).reshape(-1)
    work=Path(cubit_work_directory); work.mkdir(parents=True,exist_ok=True)
    for iteration in range(int(max_iterations)):
        step=linearize_step(state); update=step.update
        old_u=np.asarray(state.normal_displacement); old_t=np.asarray(state.thickness); old_r=np.asarray(state.activation)
        new_u=np.asarray(update.normal_displacement); new_t=np.asarray(update.thickness); new_r=np.asarray(update.activation)
        if not (old_u.shape==old_t.shape==old_r.shape==new_u.shape==new_t.shape==new_r.shape==sizes.shape):
            raise ValueError("HEX sheet design/element-size shape mismatch")
        relative=np.abs(new_u-old_u)/sizes
        topology_changed=bool(np.any((old_r>=activation_threshold)!=(new_r>=activation_threshold)))
        acceptance,_=backtrack_ngsolve_target_deformation(state.mesh,deformation_factory,
            old_u,new_u,relative,minimum_scale=minimum_scale,contraction=contraction,
            minimum_jacobian_ratio=minimum_jacobian_ratio,
            maximum_condition=maximum_condition,refine_threshold=refine_threshold,
            rebuild_threshold=rebuild_threshold,integration_order=integration_order,
            topology_changed=topology_changed)
        scale=acceptance.scale if acceptance.accepted else 0.0
        decision=acceptance.decision
        use_cubit=(not acceptance.accepted or decision.route=="cubit_rebuild"
                   or (decision.route=="ngsolve_refine" and not allow_ngsolve_refine))
        if use_cubit:
            scale=1.0 if not acceptance.accepted else scale
            u=old_u+scale*(new_u-old_u); t=old_t+scale*(new_t-old_t); r=old_r+scale*(new_r-old_r)
            request=CubitHexRemeshRequest(iteration,u,t,r,
                work/f"hex_topopt_{iteration:04d}.jou",work/f"hex_topopt_{iteration:04d}.vol")
            mesh=cubit_backend.rebuild(request); route="cubit_rebuild"
        else:
            u=old_u+scale*(new_u-old_u); t=old_t+scale*(new_t-old_t); r=old_r+scale*(new_r-old_r)
            route=decision.route; mesh=state.mesh
            if route=="ngsolve_refine":
                apply_ngsolve_mesh_route(mesh,deformation_factory(mesh,u),decision)
        model=rebuild_model(mesh,u,t,r,route)
        objective=float(evaluate_objective(model))
        if not np.isfinite(objective): raise RuntimeError("non-finite topology objective")
        change=float(max(np.max(np.abs(u-old_u)),np.max(np.abs(t-old_t)),np.max(np.abs(r-old_r))))
        history.append(HexSheetTopologyIteration(iteration,float(state.objective),objective,
            float(scale),route,decision.minimum_jacobian,decision.maximum_condition,topology_changed))
        previous=float(state.objective)
        state=HexSheetTopologyState(mesh,model,u,t,r,objective)
        relative_objective=abs(objective-previous)/max(1.0,abs(previous))
        if change<=design_tolerance and relative_objective<=objective_tolerance:
            converged=True; break
    return HexSheetTopologyResult(state,tuple(history),converged)


__all__=["SheetMetalUpdate","MeshUpdateDecision","DeformationAcceptance",
         "AffineGetTrafoCells","HexSheetTopologyState","HexSheetTopologyIteration",
         "HexSheetTopologyResult","CubitHexRemeshRequest","CubitHexRemeshBackend",
         "sample_affine_gettrafo_cells","solve_sheet_metal_lp","local_trust_region",
         "route_mesh_update","apply_ngsolve_mesh_route","sample_trafo_quality",
         "backtrack_ngsolve_deformation","backtrack_ngsolve_target_deformation",
         "optimize_hex_sheet_topology"]
