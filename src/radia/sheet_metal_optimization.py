"""Manufacturing-aware sheet-metal updates for Radia-VIM topology workflows.

NGSolve owns the deformation and element transformations.  This module builds
the bounded LP update and decides when deformation is still safe, when local
refinement is required, and when a topology-changing Cubit rebuild is required.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
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
    pending_topology_changes: int


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
class ShapeModelEvaluation:
    """Objective and physical responses from one fully solved shape."""
    objective: float
    response: np.ndarray


@dataclass(frozen=True)
class TopologyPreservingShapeState:
    """Real iron geometry; ``reference_parameters`` own the current mesh."""
    mesh: object
    model: object
    reference_parameters: np.ndarray
    parameters: np.ndarray
    evaluation: ShapeModelEvaluation


@dataclass(frozen=True)
class TopologyPreservingShapeIteration:
    iteration: int
    objective_before: float
    objective_after: float
    max_band_ratio_before: float
    max_band_ratio_after: float
    accepted_scale: float
    maximum_parameter_change: float
    route: str
    minimum_jacobian: float
    maximum_condition: float
    nonlinear_resolves: int
    remesh_attempted: bool = False
    remesh_accepted: bool = False
    remesh_reason: str = ""


@dataclass(frozen=True)
class TopologyPreservingShapeResult:
    state: TopologyPreservingShapeState
    history: tuple[TopologyPreservingShapeIteration, ...]
    converged: bool


@dataclass(frozen=True)
class CubitShapeRemeshRequest:
    iteration: int
    shape_parameters: np.ndarray
    journal_path: Path
    mesh_path: Path
    source_mesh: object = None
    source_deformation: object = None


@dataclass(frozen=True)
class CubitShapeRemeshResult:
    """Mesh plus the machine-readable Cubit/Sculpt validation report."""
    mesh: object
    report: object = None


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


@dataclass(frozen=True)
class CubitSculptShapeRemeshBackend:
    """Anisotropic Sculpt checkpoint backend for long topology-fixed magnets.

    ``sculpt`` is normally ``radia_mcp.cubit.server.cubit_stl_to_vol``.  It is
    injected so the numerical Radia package keeps the Cubit execution layer
    optional and testable.  The accepted GetTrafo exterior is compressed by
    ``coordinate_scale`` before isotropic Sculpt and expanded afterwards.
    """
    sculpt: object
    boundary_classifier: object
    coordinate_scale: tuple[float, float, float] = (1.0, 0.16, 1.0)
    sculpt_size: float = 0.039
    closure_tolerance: float = 0.03
    material_name: str = "iron"
    required_boundaries: tuple[str, ...] = ("edge", "fixed", "sym_z")
    timeout_s: int = 900
    gq_iters: int = 0
    gq_threshold: float = 0.2
    mesh_check: object = None
    tetrahedralize_for_analysis: bool = False

    def rebuild(self, request: CubitShapeRemeshRequest):
        from .topopt_cad import (exact_surface_stl_from_mesh,
                                 nearest_boundary_label_classifier,
                                 relabel_straight_mesh,
                                 rescale_netgen_vol_points)
        import ngsolve as ng

        if request.source_mesh is None:
            raise ValueError(
                "Sculpt shape rebuild requires request.source_mesh")
        scale = np.asarray(self.coordinate_scale, dtype=float).reshape(-1)
        if (scale.shape != (3,) or not np.all(np.isfinite(scale))
                or np.any(scale <= 0.0)):
            raise ValueError(
                "Sculpt coordinate_scale must have three positive entries")
        size = float(self.sculpt_size)
        closure = float(self.closure_tolerance)
        if (not np.isfinite(size) or size <= 0.0
                or not np.isfinite(closure) or closure < 0.0):
            raise ValueError("invalid Sculpt size or closure tolerance")
        output = Path(request.mesh_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        stem = output.with_suffix("")
        surface_path = stem.with_name(stem.name + "_scaled.stl")
        scaled_vol = stem.with_name(stem.name + "_scaled.vol")
        scaled_msh = stem.with_name(stem.name + "_scaled.msh")
        physical_unlabeled = stem.with_name(stem.name + "_physical_raw.vol")
        surface = exact_surface_stl_from_mesh(
            request.source_mesh, surface_path,
            deformation=request.source_deformation,
            coordinate_scale=scale)
        raw = self.sculpt(
            stl_path=str(surface_path), scheme="hex", size=size,
            closure_tolerance=closure, out_vol=str(scaled_vol),
            out_msh=str(scaled_msh), timeout_s=int(self.timeout_s),
            gq_iters=int(self.gq_iters),
            gq_threshold=float(self.gq_threshold))
        report = json.loads(raw) if isinstance(raw, str) else dict(raw)
        if report.get("status") != "ok" or not scaled_vol.is_file():
            raise RuntimeError(
                "Cubit Sculpt did not pass its scaled-mesh gates: "
                + json.dumps({"status": report.get("status"),
                              "gates": report.get("gates")},
                             ensure_ascii=False))
        rescale_netgen_vol_points(
            scaled_vol, physical_unlabeled, 1.0 / scale)
        raw_mesh = ng.Mesh(str(physical_unlabeled))
        transferred_classifier = nearest_boundary_label_classifier(
            request.source_mesh, deformation=request.source_deformation,
            fallback=self.boundary_classifier)
        labeled = relabel_straight_mesh(
            raw_mesh, transferred_classifier,
            material_name=self.material_name)
        sculpt_elements = tuple(labeled.Elements(ng.VOL))
        sculpt_families = {
            {4: "tet", 6: "wedge", 8: "hex"}.get(len(element.vertices),
                                                     "unsupported")
            for element in sculpt_elements
        }
        if bool(self.tetrahedralize_for_analysis):
            if sculpt_families != {"hex"}:
                raise RuntimeError(
                    "tetrahedralize_for_analysis requires a pure Sculpt HEX "
                    f"mesh, got {sorted(sculpt_families)}")
            # Netgen owns the conforming six-TET split, including consistent
            # diagonals on shared and boundary quadrilaterals.  Sculpt still
            # owns the accepted exterior; this is only the one-time handoff
            # to the topology-fixed exact TET HDiv-MMM/Trafo path.
            labeled.ngmesh.Split2Tets()
            labeled = ng.Mesh(labeled.ngmesh)
        labeled.ngmesh.Save(str(output))
        mesh = ng.Mesh(str(output))
        analysis_elements = tuple(mesh.Elements(ng.VOL))
        analysis_families = sorted({
            {4: "tet", 6: "wedge", 8: "hex"}.get(len(element.vertices),
                                                     "unsupported")
            for element in analysis_elements
        })
        boundaries = sorted(set(map(str, mesh.GetBoundaries())))
        materials = sorted(set(map(str, mesh.GetMaterials())))
        missing_boundaries = sorted(
            set(map(str, self.required_boundaries)) - set(boundaries))
        labels_ok = (not missing_boundaries
                     and self.material_name in materials)
        physical_volume = float(ng.Integrate(1.0, mesh))
        physical_closure = abs(
            physical_volume - float(surface["physical_volume"])) / max(
                float(surface["physical_volume"]), 1.0e-300)
        physical_closure_ok = bool(physical_closure <= closure)
        external_check = None
        check_ok = True
        if self.mesh_check is not None:
            external_check = self.mesh_check(str(output))
            if isinstance(external_check, str):
                external_check = json.loads(external_check)
            check_ok = bool(external_check.get("passed", False))
        original_gates = dict(report.get("gates") or {})
        gates = {
            **original_gates,
            "closure_ok": bool(
                original_gates.get("closure_ok") is True
                and physical_closure_ok),
            "no_inverted_elements": bool(
                original_gates.get("no_inverted_elements") is True),
            "boundary_faces_ok": bool(
                original_gates.get("boundary_faces_ok") is True
                and labels_ok and check_ok),
            "labels_ok": labels_ok,
            "external_check_ok": check_ok,
        }
        combined = {
            "schema": "radia.cubit-sculpt-shape-remesh/v1",
            "status": "ok" if all(gates.values()) else "gate_failed",
            "gates": gates,
            "surface": surface,
            "coordinate_scale": scale.tolist(),
            "physical_target_sizes": (size / scale).tolist(),
            "sculpt": report,
            "sculpt_element_families": sorted(sculpt_families),
            "sculpt_element_count": len(sculpt_elements),
            "tetrahedralized_for_analysis": bool(
                self.tetrahedralize_for_analysis),
            "analysis_element_families": analysis_families,
            "analysis_element_count": len(analysis_elements),
            "physical_vol": str(output),
            "physical_volume": physical_volume,
            "physical_closure": physical_closure,
            "boundaries": boundaries,
            "materials": materials,
            "missing_boundaries": missing_boundaries,
            "external_check": external_check,
        }
        manifest = stem.with_name(stem.name + "_sculpt.json")
        manifest.write_text(
            json.dumps(combined, indent=2, ensure_ascii=False),
            encoding="utf-8")
        combined["manifest"] = str(manifest)
        if combined["status"] != "ok":
            raise RuntimeError(
                "physical Sculpt mesh failed postprocessing gates: "
                + json.dumps(gates, ensure_ascii=False))
        return CubitShapeRemeshResult(mesh, combined)


def elastic_normal_deformation_modes(mesh, scalar_boundary_modes, *,
        movable_boundaries, fixed_boundaries, order=1, shear_modulus=1.0,
        bulk_modulus=0.1, inverse="sparsecholesky"):
    """Extend boundary-normal shape modes into the volume by elasticity.

    NGSolve owns boundary normals, H1 orientation, weak assembly, and the
    volume-field solve.  Caller owns ``ngsolve.TaskManager``.
    """
    import ngsolve as ng

    modes = tuple(scalar_boundary_modes)
    if not modes:
        raise ValueError("at least one scalar boundary mode is required")
    order_value = float(order)
    shear = float(shear_modulus)
    bulk = float(bulk_modulus)
    if (not np.isfinite(order_value) or order_value < 1.0
            or order_value != np.floor(order_value)):
        raise ValueError("elastic extension order must be a positive integer")
    if (not np.isfinite(shear) or not np.isfinite(bulk)
            or shear <= 0.0 or bulk < 0.0):
        raise ValueError(
            "elastic extension moduli must satisfy shear>0 and bulk>=0")
    movable = str(movable_boundaries)
    fixed = str(fixed_boundaries)
    if not movable:
        raise ValueError("movable_boundaries must be named")
    dirichlet = movable if not fixed else movable + "|" + fixed
    fes = ng.VectorH1(mesh, order=int(order_value), dirichlet=dirichlet)
    u, v = fes.TnT()
    form = ng.BilinearForm(fes)
    form += (2.0 * shear * ng.InnerProduct(
        ng.Sym(ng.Grad(u)), ng.Sym(ng.Grad(v)))
        + bulk * ng.div(u) * ng.div(v)) * ng.dx
    form.Assemble()
    inverse_operator = form.mat.Inverse(fes.FreeDofs(), inverse=str(inverse))
    normal = ng.specialcf.normal(mesh.dim)
    boundary = mesh.Boundaries(movable)
    result = []
    for scalar in modes:
        field = ng.GridFunction(fes)
        field.Set(scalar * normal, definedon=boundary)
        residual = -form.mat * field.vec
        field.vec.data += inverse_operator * residual
        result.append(field)
    return tuple(result)


def combine_deformation_modes(displacement_modes, coefficients):
    """Form one GetTrafo field from a superposed set of shape modes."""
    import ngsolve as ng

    modes = tuple(displacement_modes)
    coeff = np.asarray(coefficients, dtype=float).reshape(-1)
    if len(modes) != coeff.size or not modes:
        raise ValueError(
            "displacement modes and coefficients must be non-empty and match")
    if not np.all(np.isfinite(coeff)):
        raise ValueError("deformation coefficients must be finite")
    space = modes[0].space
    if any(mode.space is not space for mode in modes):
        raise ValueError("all displacement modes must share one NGSolve space")
    field = ng.GridFunction(space)
    field.vec[:] = 0.0
    for value, mode in zip(coeff, modes):
        if value != 0.0:
            field.vec.data += float(value) * mode.vec
    return field


def relative_gettrafo_displacements(mesh, deformation):
    """Maximum nodal displacement divided by cell diameter, per element."""
    sampled = sample_affine_gettrafo_cells(mesh, (deformation,))
    relative = []
    for nodes, values in zip(sampled.nodes, sampled.node_displacements):
        diameter = float(np.max(np.linalg.norm(
            nodes[:, None, :] - nodes[None, :, :], axis=2)))
        if not np.isfinite(diameter) or diameter <= 0.0:
            raise ValueError(
                "degenerate cell encountered while scaling deformation")
        relative.append(
            float(np.max(np.linalg.norm(values[0], axis=1))) / diameter)
    return np.asarray(relative, dtype=float)


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


def reference_aware_condition_limit(reference_conditions, *, requested=20.0,
                                    margin=1.25):
    """Return an absolute condition limit that admits the base mesh."""
    conditions = np.asarray(reference_conditions, dtype=float).reshape(-1)
    if (conditions.size == 0 or not np.all(np.isfinite(conditions))
            or np.any(conditions < 1.0)):
        raise ValueError(
            "reference conditions must be finite and at least one")
    if not np.isfinite(requested) or float(requested) < 1.0:
        raise ValueError(
            "requested condition limit must be finite and at least one")
    if not np.isfinite(margin) or float(margin) <= 1.0:
        raise ValueError(
            "condition margin must be finite and greater than one")
    return float(max(
        float(requested), float(margin) * float(np.max(conditions))))


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


def _shape_band_ratio(response, target, band):
    values = np.asarray(response, dtype=float).reshape(-1)
    if values.size == 0:
        return 0.0
    target = np.asarray(target, dtype=float).reshape(-1)
    band = np.asarray(band, dtype=float).reshape(-1)
    if (target.shape != values.shape or band.shape != values.shape
            or np.any(band <= 0)
            or not np.all(np.isfinite(np.r_[values, target, band]))):
        raise ValueError(
            "shape evaluation response/target/band must match and be finite "
            "with positive bands")
    return float(np.max(np.abs((values - target) / band)))


def _accept_perturbative_shape_trial(
        current, trial, linearization, delta, scale, *, armijo,
        objective_tolerance, band_tolerance):
    """Accept only a fully re-solved physical shape trial."""
    current_ratio = _shape_band_ratio(
        current.response, linearization.response_target,
        linearization.response_band)
    trial_ratio = _shape_band_ratio(
        trial.response, linearization.response_target,
        linearization.response_band)
    predicted_response = (
        np.asarray(linearization.response, dtype=float)
        + float(scale)
        * np.asarray(linearization.response_jacobian, dtype=float) @ delta)
    predicted_ratio = _shape_band_ratio(
        predicted_response, linearization.response_target,
        linearization.response_band)
    if current_ratio > 1.0 + band_tolerance:
        expected = max(0.0, current_ratio - predicted_ratio)
        required = float(armijo) * max(expected, band_tolerance)
        return (
            trial_ratio <= current_ratio - required,
            current_ratio,
            trial_ratio,
        )
    if trial_ratio > 1.0 + band_tolerance:
        return False, current_ratio, trial_ratio
    predicted_change = (
        float(scale)
        * float(np.asarray(linearization.objective_gradient, dtype=float) @ delta))
    allowed = float(current.objective) + float(objective_tolerance)
    if predicted_change < 0.0:
        allowed += float(armijo) * predicted_change
    return float(trial.objective) <= allowed, current_ratio, trial_ratio


def _normalize_cubit_shape_result(value):
    """Accept the new report-bearing result without breaking old backends."""
    if isinstance(value, CubitShapeRemeshResult):
        return value
    return CubitShapeRemeshResult(value, None)


def _default_cubit_shape_gate(result):
    """Require the three solver-facing Sculpt gates when a report is present."""
    report = result.report
    if report is None:
        return True, "legacy backend supplied no mesh report"
    if not isinstance(report, dict):
        return False, "Cubit remesh report is not a mapping"
    if report.get("status") not in (None, "ok"):
        return False, f"Cubit remesh status is {report.get('status')!r}"
    gates = report.get("gates")
    if gates is None:
        return False, "Cubit remesh report has no gates"
    required = ("closure_ok", "no_inverted_elements", "boundary_faces_ok")
    failed = [name for name in required if gates.get(name) is not True]
    if failed:
        return False, "failed Cubit remesh gates: " + ", ".join(failed)
    return True, "Cubit remesh gates passed"


def _cubit_shape_equivalence(trial, remeshed, linearization, *,
                             response_tolerance, objective_tolerance):
    """Check that remeshing preserved the already accepted physical shape."""
    if response_tolerance is not None:
        tolerance = float(response_tolerance)
        if not np.isfinite(tolerance) or tolerance < 0.0:
            raise ValueError(
                "cubit_response_equivalence_tolerance must be nonnegative")
        response_delta = np.asarray(remeshed.response, dtype=float) - np.asarray(
            trial.response, dtype=float)
        if response_delta.size:
            ratio = float(np.max(np.abs(
                response_delta / np.asarray(
                    linearization.response_band, dtype=float))))
            if not np.isfinite(ratio) or ratio > tolerance:
                return False, (
                    "Sculpt response changed by "
                    f"{ratio:.6g} response-band units")
    if objective_tolerance is not None:
        tolerance = float(objective_tolerance)
        if not np.isfinite(tolerance) or tolerance < 0.0:
            raise ValueError(
                "cubit_objective_equivalence_tolerance must be nonnegative")
        relative = abs(float(remeshed.objective) - float(trial.objective)) / max(
            abs(float(trial.objective)), 1.0e-300)
        if not np.isfinite(relative) or relative > tolerance:
            return False, (
                "Sculpt objective changed by relative "
                f"{relative:.6g}")
    return True, "Sculpt physical equivalence gates passed"


def optimize_topology_preserving_shape(
        initial_state: TopologyPreservingShapeState, *, linearize_step,
        deformation_factory, rebuild_model, evaluate_model, move_limit,
        parameter_bounds=None, laplacian=None, curvature_limit=None,
        A_ub=None, b_ub=None, cubit_backend=None,
        cubit_work_directory=None, max_iterations=20,
        cubit_batch_interval=None, cubit_batch_parameter_change=None,
        cubit_remesh_gate=None,
        cubit_response_equivalence_tolerance=None,
        cubit_objective_equivalence_tolerance=None,
        parameter_tolerance=1e-4, objective_tolerance=1e-10,
        armijo=1e-4, band_tolerance=1e-8, minimum_scale=1/64,
        contraction=0.5, minimum_jacobian_ratio=0.2,
        maximum_condition=20.0, refine_threshold=0.25,
        rebuild_threshold=0.5, integration_order=2,
        iteration_callback=None):
    """Run clay-like, topology-preserving HDiv-MMM shape optimization.

    Every trial is fully re-solved.  Safe steps stay as an NGSolve
    deformation; a quality-limit crossing requests one application-owned
    Cubit rebuild without changing the accepted iron topology.  Optional
    interval/displacement batching checkpoints several accepted GetTrafo
    steps in one Sculpt rebuild.  A scheduled checkpoint that fails its mesh
    or physical-equivalence gates leaves the accepted GetTrafo shape intact;
    a quality-mandated rebuild still backtracks.
    """
    from .topology_optimization import solve_shape_lp

    max_iterations_value = float(max_iterations)
    if (not np.isfinite(max_iterations_value) or max_iterations_value < 1.0
            or max_iterations_value != np.floor(max_iterations_value)):
        raise ValueError("max_iterations must be a positive integer")
    for value, name in (
            (parameter_tolerance, "parameter_tolerance"),
            (objective_tolerance, "objective_tolerance"),
            (band_tolerance, "band_tolerance")):
        if not np.isfinite(value) or float(value) < 0.0:
            raise ValueError(f"{name} must be finite and nonnegative")
    if not (0 < contraction < 1 and 0 < minimum_scale <= 1):
        raise ValueError(
            "invalid topology-preserving shape backtracking controls")
    if not (0 < armijo <= 1):
        raise ValueError("armijo must be in (0,1]")
    if cubit_batch_interval is not None:
        interval_value = float(cubit_batch_interval)
        if (not np.isfinite(interval_value) or interval_value < 1.0
                or interval_value != np.floor(interval_value)):
            raise ValueError("cubit_batch_interval must be a positive integer")
        cubit_batch_interval = int(interval_value)
    if cubit_batch_parameter_change is not None:
        cubit_batch_parameter_change = float(cubit_batch_parameter_change)
        if (not np.isfinite(cubit_batch_parameter_change)
                or cubit_batch_parameter_change <= 0.0):
            raise ValueError(
                "cubit_batch_parameter_change must be finite and positive")
    if ((cubit_batch_interval is not None
         or cubit_batch_parameter_change is not None)
            and cubit_backend is None):
        raise ValueError("Cubit batching controls require cubit_backend")
    state = initial_state
    history = []
    converged = False
    q = np.asarray(state.parameters, dtype=float).reshape(-1)
    qref = np.asarray(state.reference_parameters, dtype=float).reshape(-1)
    if (q.size == 0 or qref.shape != q.shape
            or not np.all(np.isfinite(q)) or not np.all(np.isfinite(qref))):
        raise ValueError(
            "shape state requires matching finite non-empty parameter vectors")
    work = (None if cubit_work_directory is None
            else Path(cubit_work_directory))
    if work is not None:
        work.mkdir(parents=True, exist_ok=True)
    accepted_since_cubit = 0

    for iteration in range(int(max_iterations_value)):
        linearization = linearize_step(state)
        update = solve_shape_lp(
            q, linearization, move_limit=move_limit,
            parameter_bounds=parameter_bounds, laplacian=laplacian,
            curvature_limit=curvature_limit, A_ub=A_ub, b_ub=b_ub)
        if float(np.max(np.abs(update.delta))) <= float(parameter_tolerance):
            converged = True
            break
        mesh = state.mesh
        if getattr(mesh, "deformation", None) is not None:
            mesh.UnsetDeformation()
        reference_determinants, _ = sample_trafo_quality(
            mesh, integration_order=integration_order)
        scale = 1.0
        accepted = None
        nonlinear_resolves = 0
        while scale >= minimum_scale:
            candidate = q + scale * update.delta
            deformation = deformation_factory(mesh, qref, candidate)
            relative = relative_gettrafo_displacements(mesh, deformation)
            mesh.SetDeformation(deformation)
            try:
                ratios, conditions = sample_trafo_quality(
                    mesh, integration_order=integration_order,
                    reference_determinants=reference_determinants)
                decision = route_mesh_update(
                    ratios, conditions, relative,
                    refine_threshold=refine_threshold,
                    rebuild_threshold=rebuild_threshold,
                    minimum_jacobian=minimum_jacobian_ratio,
                    maximum_condition=maximum_condition,
                    topology_changed=False)
                unsafe = (np.any(ratios <= minimum_jacobian_ratio)
                          or np.any(conditions >= 2 * maximum_condition))
                needs_cubit = decision.route != "ngsolve_deform"
                if unsafe or (needs_cubit and cubit_backend is None):
                    mesh.UnsetDeformation()
                    scale *= contraction
                    continue
                trial_model = rebuild_model(
                    mesh, candidate, "ngsolve_deform")
                trial_evaluation = evaluate_model(trial_model)
                nonlinear_resolves += 1
                ok, before_ratio, after_ratio = (
                    _accept_perturbative_shape_trial(
                        state.evaluation, trial_evaluation, linearization,
                        update.delta, scale, armijo=armijo,
                        objective_tolerance=objective_tolerance,
                        band_tolerance=band_tolerance))
            except Exception:
                if getattr(mesh, "deformation", None) is not None:
                    mesh.UnsetDeformation()
                raise
            if not ok:
                mesh.UnsetDeformation()
                scale *= contraction
                continue
            route = "ngsolve_deform"
            next_mesh = mesh
            next_model = trial_model
            next_evaluation = trial_evaluation
            next_reference = qref.copy()
            scheduled_cubit = bool(
                cubit_backend is not None and not needs_cubit and (
                    (cubit_batch_interval is not None
                     and accepted_since_cubit + 1 >= cubit_batch_interval)
                    or (cubit_batch_parameter_change is not None
                        and float(np.max(np.abs(candidate - qref)))
                        >= cubit_batch_parameter_change)))
            attempt_cubit = bool(needs_cubit or scheduled_cubit)
            remesh_attempted = False
            remesh_accepted = False
            remesh_reason = ""
            if attempt_cubit:
                mesh.UnsetDeformation()
                if work is None:
                    raise ValueError(
                        "cubit_work_directory is required with cubit_backend")
                request = CubitShapeRemeshRequest(
                    iteration, candidate.copy(),
                    work / f"shape_{iteration:04d}.jou",
                    work / f"shape_{iteration:04d}.vol",
                    source_mesh=mesh, source_deformation=deformation)
                remesh_attempted = True
                try:
                    result = _normalize_cubit_shape_result(
                        cubit_backend.rebuild(request))
                    gate = (cubit_remesh_gate or
                            _default_cubit_shape_gate)(result)
                    if isinstance(gate, tuple):
                        mesh_gate_ok, remesh_reason = bool(gate[0]), str(gate[1])
                    else:
                        mesh_gate_ok = bool(gate)
                        remesh_reason = (
                            "application Cubit remesh gate passed"
                            if mesh_gate_ok else
                            "application Cubit remesh gate failed")
                except Exception as exc:
                    if not scheduled_cubit:
                        raise
                    mesh_gate_ok = False
                    remesh_reason = f"scheduled Sculpt failed: {exc}"
                if mesh_gate_ok:
                    remeshed_mesh = result.mesh
                    remeshed_model = rebuild_model(
                        remeshed_mesh, candidate, "cubit_rebuild")
                    remeshed_evaluation = evaluate_model(remeshed_model)
                    nonlinear_resolves += 1
                    remesh_ok, remesh_before_ratio, remesh_after_ratio = (
                        _accept_perturbative_shape_trial(
                            state.evaluation, remeshed_evaluation,
                            linearization, update.delta, scale, armijo=armijo,
                            objective_tolerance=objective_tolerance,
                            band_tolerance=band_tolerance))
                    equivalent, equivalence_reason = _cubit_shape_equivalence(
                        trial_evaluation, remeshed_evaluation, linearization,
                        response_tolerance=
                            cubit_response_equivalence_tolerance,
                        objective_tolerance=
                            cubit_objective_equivalence_tolerance)
                    remesh_ok = bool(remesh_ok and equivalent)
                    remesh_reason = remesh_reason + "; " + equivalence_reason
                    if remesh_ok:
                        next_mesh = remeshed_mesh
                        next_model = remeshed_model
                        next_evaluation = remeshed_evaluation
                        next_reference = candidate.copy()
                        before_ratio = remesh_before_ratio
                        after_ratio = remesh_after_ratio
                        route = "cubit_rebuild"
                        remesh_accepted = True
                if not remesh_accepted:
                    if needs_cubit:
                        scale *= contraction
                        continue
                    # A scheduled checkpoint is an optimization accelerator,
                    # not a reason to discard an already accepted physical
                    # GetTrafo step.  Restore that exact deformation/model.
                    mesh.SetDeformation(deformation)
            accepted = (
                candidate, next_mesh, next_model, next_evaluation,
                next_reference, decision, before_ratio, after_ratio,
                scale, route, remesh_attempted, remesh_accepted,
                remesh_reason)
            break
        if accepted is None:
            current_deformation = deformation_factory(mesh, qref, q)
            if np.max(np.abs(q - qref)) > 0.0:
                mesh.SetDeformation(current_deformation)
            break
        (candidate, next_mesh, next_model, next_evaluation, next_reference,
         decision, before_ratio, after_ratio, accepted_scale, route,
         remesh_attempted, remesh_accepted, remesh_reason) = accepted
        change = float(np.max(np.abs(candidate - q)))
        objective_before = float(state.evaluation.objective)
        state = TopologyPreservingShapeState(
            next_mesh, next_model, next_reference, candidate.copy(),
            next_evaluation)
        history.append(TopologyPreservingShapeIteration(
            iteration, objective_before, float(next_evaluation.objective),
            before_ratio, after_ratio, float(accepted_scale), change, route,
            decision.minimum_jacobian, decision.maximum_condition,
            nonlinear_resolves, remesh_attempted, remesh_accepted,
            remesh_reason))
        if iteration_callback is not None:
            iteration_callback(history[-1], state)
        q = candidate.copy()
        qref = next_reference.copy()
        accepted_since_cubit = (
            0 if remesh_accepted else accepted_since_cubit + 1)
        if change <= float(parameter_tolerance):
            converged = True
            break
    return TopologyPreservingShapeResult(state, tuple(history), converged)


def optimize_hex_sheet_topology(initial_state: HexSheetTopologyState, *,
        linearize_step, deformation_factory, rebuild_model, evaluate_objective,
        element_sizes, cubit_backend: CubitHexRemeshBackend,
        cubit_work_directory, max_iterations=20, objective_tolerance=1e-4,
        design_tolerance=1e-3, activation_threshold=0.5,
        activation_remove_threshold=0.35, activation_restore_threshold=0.65,
        cubit_batch_interval=5, cubit_batch_fraction=0.05,
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
    if not 0<=activation_remove_threshold<activation_restore_threshold<=1: raise ValueError("activation hysteresis thresholds must satisfy 0 <= remove < restore <= 1")
    if int(cubit_batch_interval)<1 or not 0<cubit_batch_fraction<=1: raise ValueError("invalid Cubit batching controls")
    state=initial_state; history=[]; converged=False
    sizes=np.asarray(element_sizes,dtype=float).reshape(-1)
    committed_topology=np.asarray(state.activation)>=float(activation_threshold)
    pending_topology=np.zeros_like(committed_topology,dtype=bool);last_cubit_iteration=0
    work=Path(cubit_work_directory); work.mkdir(parents=True,exist_ok=True)
    for iteration in range(int(max_iterations)):
        step=linearize_step(state); update=step.update
        old_u=np.asarray(state.normal_displacement); old_t=np.asarray(state.thickness); old_r=np.asarray(state.activation)
        new_u=np.asarray(update.normal_displacement); new_t=np.asarray(update.thickness); new_r=np.asarray(update.activation)
        if not (old_u.shape==old_t.shape==old_r.shape==new_u.shape==new_t.shape==new_r.shape==sizes.shape):
            raise ValueError("HEX sheet design/element-size shape mismatch")
        relative=np.abs(new_u-old_u)/sizes
        desired_topology=committed_topology.copy()
        desired_topology[committed_topology & (new_r<=activation_remove_threshold)]=False
        desired_topology[(~committed_topology) & (new_r>=activation_restore_threshold)]=True
        # Cancel a queued change if the continuous design returns through the
        # hysteresis band before the batch is committed.
        pending_topology=(desired_topology!=committed_topology)
        pending_count=int(np.count_nonzero(pending_topology))
        topology_changed=bool(pending_count>0 and (pending_count/max(1,pending_topology.size)>=cubit_batch_fraction or iteration-last_cubit_iteration+1>=int(cubit_batch_interval)))
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
            committed_topology=desired_topology.copy();pending_topology[:]=False;last_cubit_iteration=iteration+1
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
            float(scale),route,decision.minimum_jacobian,decision.maximum_condition,topology_changed,pending_count))
        previous=float(state.objective)
        state=HexSheetTopologyState(mesh,model,u,t,r,objective)
        relative_objective=abs(objective-previous)/max(1.0,abs(previous))
        if change<=design_tolerance and relative_objective<=objective_tolerance:
            converged=True; break
    return HexSheetTopologyResult(state,tuple(history),converged)


__all__=["SheetMetalUpdate","MeshUpdateDecision","DeformationAcceptance",
         "AffineGetTrafoCells","HexSheetTopologyState","HexSheetTopologyIteration",
         "HexSheetTopologyResult","CubitHexRemeshRequest","CubitHexRemeshBackend",
         "CubitSculptShapeRemeshBackend",
         "ShapeModelEvaluation","TopologyPreservingShapeState",
         "TopologyPreservingShapeIteration","TopologyPreservingShapeResult",
         "CubitShapeRemeshRequest","CubitShapeRemeshResult",
         "elastic_normal_deformation_modes",
         "combine_deformation_modes","relative_gettrafo_displacements",
         "reference_aware_condition_limit",
         "optimize_topology_preserving_shape",
         "sample_affine_gettrafo_cells","solve_sheet_metal_lp","local_trust_region",
         "route_mesh_update","apply_ngsolve_mesh_route","sample_trafo_quality",
         "backtrack_ngsolve_deformation","backtrack_ngsolve_target_deformation",
         "optimize_hex_sheet_topology"]
