"""Binary graph-front operators for the HDiv-MMM Lego search.

This internal module contains no FE basis reconstruction and no design finite
difference.  It operates only on exact/analytic candidate response columns and
on the element-face graph supplied by NGSolve.  The physical HDiv-MMM solve,
Schur reduction, particle adjoint, and final topology gate remain authoritative.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class GraphFrontProposal:
    candidate_indices: np.ndarray
    predicted_response: np.ndarray
    predicted_max_band_ratio: float
    regularization_change: float
    connected_components: int


@dataclass(frozen=True)
class GraphFrontTrustUpdate:
    budget_before: int
    budget_after: int
    interface_weight_before: float
    interface_weight_after: float
    agreement_ratio: float
    action: str


@dataclass(frozen=True)
class GraphFrontDiversityDiagnostics:
    proposal_count: int
    numerical_rank: int
    duplicate_pair_fraction: float
    maximum_absolute_correlation: float
    minimum_subspace_novelty: float


@dataclass(frozen=True)
class GraphFrontBeamResult:
    proposals: tuple[GraphFrontProposal, ...]
    pool_diagnostics: GraphFrontDiversityDiagnostics
    selected_diagnostics: GraphFrontDiversityDiagnostics


def best_admissible_singleton(representatives, predicted_ratios, *,
                              current_ratio, is_valid,
                              improvement_tolerance=1.0e-8):
    """Return the globally best improving physical singleton.

    Graph seeds are a beam-search heuristic and must not restrict this safety
    lane.  ``is_valid`` is the authoritative volume/topology/connectivity gate.
    The input order is deliberately ignored when comparing predicted ratios.
    """
    values=tuple(map(int,representatives))
    ratios=np.asarray(predicted_ratios,dtype=float).reshape(-1)
    current=float(current_ratio);tolerance=float(improvement_tolerance)
    if (ratios.shape!=(len(values),) or not np.isfinite(current) or
            tolerance<0.0 or not np.all(np.isfinite(ratios))):
        raise ValueError("singleton screening inputs are invalid")
    improving=[]
    for representative,predicted in zip(values,ratios):
        if (predicted<current-tolerance and
                bool(is_valid(int(representative)))):
            improving.append((float(predicted),int(representative)))
    return None if not improving else min(improving)[1]


def minimax_driving_potential(current_response, response_target,
                              response_band, candidate_response_delta, *,
                              tie_tolerance=1.0e-10):
    """Return a positive-for-improvement subgradient contraction.

    The candidate columns are the *physical signed move* ``r(x+move)-r(x)``.
    Only this scalar contraction is filtered.  The response columns themselves
    must remain untouched in every LP/MILP and exact acceptance calculation.
    """
    current=np.asarray(current_response,dtype=float).reshape(-1)
    target=np.asarray(response_target,dtype=float).reshape(-1)
    band=np.asarray(response_band,dtype=float).reshape(-1)
    delta=np.asarray(candidate_response_delta,dtype=float)
    if (current.shape!=target.shape or band.shape!=current.shape or
            current.size==0 or np.any(band<=0.0) or
            delta.ndim!=2 or delta.shape[0]!=current.size or
            not np.all(np.isfinite(np.r_[current,target,band,delta.ravel()]))):
        raise ValueError("minimax driving-potential arrays are incompatible")
    normalized=(current-target)/band
    maximum=float(np.max(np.abs(normalized)))
    if maximum==0.0:
        return np.zeros(delta.shape[1],dtype=float),np.zeros_like(normalized)
    tolerance=max(float(tie_tolerance),float(tie_tolerance)*maximum)
    active=np.flatnonzero(np.abs(normalized)>=maximum-tolerance)
    subgradient=np.zeros_like(normalized)
    subgradient[active]=np.sign(normalized[active])/len(active)
    driving=-(subgradient[:,None]*(delta/band[:,None])).sum(axis=0)
    return np.asarray(driving,dtype=float),subgradient


def candidate_face_adjacency(candidate_members: Sequence[Iterable[int]],
                             element_adjacency: Sequence[Iterable[int]], *,
                             exclusion_groups=None):
    """Lift an NGSolve element-face graph to finite add/remove pack nodes."""
    members=tuple(set(map(int,value)) for value in candidate_members)
    nc=len(members)
    exclusions=(np.full(nc,-1,dtype=np.int64) if exclusion_groups is None else
                np.asarray(exclusion_groups,dtype=np.int64).reshape(-1))
    if exclusions.shape!=(nc,):
        raise ValueError("candidate exclusion groups must match candidate count")
    owner={}
    for candidate,cells in enumerate(members):
        if not cells:
            raise ValueError("each graph-front candidate needs an element pack")
        for cell in cells:
            if cell<0 or cell>=len(element_adjacency):
                raise ValueError("candidate pack contains an invalid element")
            owner.setdefault(cell,[]).append(candidate)
    adjacency=[set() for _ in range(nc)]
    for left,cells in enumerate(members):
        neighbours=set()
        for cell in cells:
            neighbours.update(map(int,element_adjacency[cell]))
        for cell in neighbours:
            for right in owner.get(cell,()):
                if left==right:
                    continue
                # Cumulative depths in one terminal family are alternatives,
                # not neighbouring moves in a connected cluster.
                if exclusions[left]>=0 and exclusions[left]==exclusions[right]:
                    continue
                adjacency[left].add(int(right));adjacency[right].add(int(left))
    return tuple(np.asarray(sorted(row),dtype=np.int64) for row in adjacency)


def binary_graph_interface_energy(active_elements, element_adjacency, *,
                                  exterior_degree=None, edge_weights=None):
    """Weighted binary cut energy, including an optional exterior-air node.

    ``edge_weights[(min(e,f),max(e,f))]`` may contain NGSolve facet measures.
    Unit weights give the dimensionless graph-cut version.
    ``exterior_degree[e]`` counts/weights exposed
    superset facets connected to exterior air.
    """
    active=np.asarray(active_elements,dtype=bool).reshape(-1)
    if len(element_adjacency)!=active.size:
        raise ValueError("active set and element graph have incompatible sizes")
    exterior=(np.zeros(active.size,dtype=float) if exterior_degree is None else
              np.asarray(exterior_degree,dtype=float).reshape(-1))
    if exterior.shape!=active.shape or np.any(exterior<0.0):
        raise ValueError("exterior graph weights must match the active set")
    energy=float(exterior@active.astype(float))
    for left,row in enumerate(element_adjacency):
        for right_value in row:
            right=int(right_value)
            if right<=left:
                continue
            weight=(1.0 if edge_weights is None else float(
                edge_weights.get((left,right),edge_weights.get((right,left),1.0))))
            if weight<0.0 or not np.isfinite(weight):
                raise ValueError("graph interface weights must be finite/nonnegative")
            energy+=weight*float(active[left]!=active[right])
    return float(energy)


def ngsolve_facet_measure_graph(mesh):
    """Return the volume graph with NGSolve-integrated facet measures.

    A discontinuous element-boundary integral on order-zero ``FacetFESpace``
    gives one mass diagonal per facet.  Interior facets occur in two volume
    element boundaries, so their diagonal is divided by two.  This delegates
    curved mappings, quadrature, topology, and mixed TET/HEX/WEDGE geometry to
    NGSolve instead of reconstructing polygon areas in Python.
    """
    import ngsolve as ng

    elements=tuple(mesh.Elements(ng.VOL));owners={}
    for index,element in enumerate(elements):
        for facet in element.facets:
            owners.setdefault(int(facet.nr),[]).append(index)
    if not owners:
        raise ValueError("NGSolve facet graph needs volume elements")
    fes=ng.FacetFESpace(mesh,order=0)
    if int(fes.ndof)!=int(mesh.nfacet):
        raise RuntimeError("order-zero FacetFESpace must own one DOF per facet")
    facet_dof={}
    for element in elements:
        facets=tuple(int(value.nr) for value in element.facets)
        dofs=tuple(int(value) for value in fes.GetDofNrs(element))
        if len(facets)!=len(dofs):
            raise RuntimeError("NGSolve facet/DOF incidence is inconsistent")
        for facet,dof in zip(facets,dofs):
            previous=facet_dof.setdefault(facet,dof)
            if previous!=dof:
                raise RuntimeError("NGSolve facet DOF changes across neighbours")
    u,v=fes.TnT();mass=ng.BilinearForm(fes)
    mass+=u*v*ng.dx(element_boundary=True)
    mass.Assemble()
    adjacency=[set() for _ in elements]
    exterior=np.zeros(len(elements),dtype=float);weights={}
    for facet,cells in owners.items():
        if len(cells) not in (1,2):
            raise RuntimeError("non-manifold volume facet has invalid owner count")
        diagonal=float(mass.mat[facet_dof[facet],facet_dof[facet]])
        area=diagonal/len(cells)
        if not np.isfinite(area) or area<=0.0:
            raise RuntimeError("NGSolve returned a nonpositive facet measure")
        if len(cells)==1:
            exterior[cells[0]]+=area
        else:
            left,right=map(int,cells)
            adjacency[left].add(right);adjacency[right].add(left)
            weights[(min(left,right),max(left,right))]=area
    return (tuple(np.asarray(sorted(row),dtype=np.int64)
                  for row in adjacency),exterior,weights)


def terminal_l1_curvature_energy(active_elements, radial_coordinates,
                                 longitudinal_coordinates,
                                 designable_elements, *, total_length):
    """Dimensionless L1 second-difference energy of both terminal fronts."""
    active=np.asarray(active_elements,dtype=bool).reshape(-1)
    radial=np.asarray(radial_coordinates,dtype=float).reshape(-1)
    longitudinal=np.asarray(longitudinal_coordinates,dtype=float).reshape(-1)
    designable=np.asarray(designable_elements,dtype=bool).reshape(-1)
    if not (active.shape==radial.shape==longitudinal.shape==designable.shape):
        raise ValueError("terminal-front geometry must match active elements")
    length=float(total_length)
    if not np.isfinite(length) or length<=0.0:
        raise ValueError("terminal-front total length must be positive")
    rounded=np.round(radial,10)
    radii=np.unique(rounded[designable])
    if radii.size<3:
        return 0.0
    pitch=float(np.median(np.diff(np.sort(radii))))
    if pitch<=0.0:
        raise ValueError("terminal-front radial stations must be distinct")
    energy=0.0
    for entrance in (True,False):
        depths=[]
        for radius in np.sort(radii):
            mask=(designable&(rounded==radius)&
                  ((longitudinal<0.5*length) if entrance else
                   (longitudinal>=0.5*length))&active)
            if not np.any(mask):
                depths.append(0.0)
            elif entrance:
                depths.append(0.5*length-float(np.min(longitudinal[mask])))
            else:
                depths.append(float(np.max(longitudinal[mask]))-0.5*length)
        second=np.diff(np.asarray(depths,dtype=float),n=2)
        energy+=float(np.sum(np.abs(second))/pitch)
    return energy


def _selected_components(selected, adjacency):
    selected=set(map(int,selected));count=0
    while selected:
        count+=1;stack=[selected.pop()]
        while stack:
            for neighbour in adjacency[stack.pop()]:
                neighbour=int(neighbour)
                if neighbour in selected:
                    selected.remove(neighbour);stack.append(neighbour)
    return count


def _response_subspace_novelty(units, index, selected):
    """Return the residual norm after projection on selected directions."""
    vector=units[:,int(index)]
    if np.linalg.norm(vector)<=64.0*np.finfo(float).eps:
        return 0.0
    if not selected:
        return 1.0
    basis=np.column_stack([units[:,int(value)] for value in selected])
    left,singular,_=np.linalg.svd(basis,full_matrices=False)
    if singular.size==0 or singular[0]==0.0:
        return 1.0
    rank=int(np.count_nonzero(
        singular>1.0e-12*singular[0]))
    if rank==0:
        return 1.0
    projection=left[:,:rank].T@vector
    return float(np.sqrt(max(0.0,1.0-float(projection@projection))))


def graph_front_response_diversity(proposals, current_response,
                                   response_band, *,
                                   relative_tolerance=1.0e-3,
                                   duplicate_correlation=0.995):
    """Measure rank and directional duplication of graph-front responses."""
    values=tuple(proposals)
    current=np.asarray(current_response,dtype=float).reshape(-1)
    band=np.asarray(response_band,dtype=float).reshape(-1)
    tolerance=float(relative_tolerance)
    duplicate=float(duplicate_correlation)
    if (band.shape!=current.shape or np.any(band<=0.0) or
            not 0.0<tolerance<1.0 or not 0.0<=duplicate<=1.0):
        raise ValueError("graph-front diversity settings are invalid")
    if not values:
        return GraphFrontDiversityDiagnostics(0,0,0.0,0.0,0.0)
    columns=np.column_stack([
        (np.asarray(item.predicted_response,dtype=float).reshape(-1)-current)
        /band for item in values])
    if columns.shape[0]!=current.size or not np.all(np.isfinite(columns)):
        raise ValueError("graph-front proposal response is invalid")
    singular=np.linalg.svd(columns,compute_uv=False)
    rank=(0 if singular.size==0 or singular[0]==0.0 else int(
        np.count_nonzero(singular>tolerance*singular[0])))
    norms=np.linalg.norm(columns,axis=0)
    units=np.zeros_like(columns)
    nonzero=norms>64.0*np.finfo(float).eps
    units[:,nonzero]=columns[:,nonzero]/norms[nonzero]
    correlation=np.abs(units.T@units)
    pair=np.triu_indices(len(values),k=1)
    pair_values=correlation[pair]
    duplicate_fraction=(0.0 if pair_values.size==0 else float(
        np.mean(pair_values>=duplicate)))
    maximum=(0.0 if pair_values.size==0 else float(np.max(pair_values)))
    sequential_novelty=[
        _response_subspace_novelty(units,index,tuple(range(index)))
        for index in range(1,len(values))]
    minimum_novelty=(1.0 if not sequential_novelty else float(
        np.min(sequential_novelty)))
    return GraphFrontDiversityDiagnostics(
        len(values),rank,duplicate_fraction,maximum,minimum_novelty)


def _select_response_diverse_proposals(proposals, current_response,
                                       response_target, response_band, *, proposal_limit,
                                       novelty_weight):
    ranked=tuple(proposals);limit=min(len(ranked),int(proposal_limit))
    weight=float(novelty_weight)
    if limit<=0:
        return ()
    if weight<=0.0:
        return ranked[:limit]
    current=np.asarray(current_response,dtype=float).reshape(-1)
    target=np.asarray(response_target,dtype=float).reshape(-1)
    band=np.asarray(response_band,dtype=float).reshape(-1)
    current_ratio=float(np.max(np.abs((current-target)/band)))
    improvements=np.asarray([
        current_ratio-float(item.predicted_max_band_ratio)
        for item in ranked],dtype=float)
    improvement_scale=max(float(np.max(np.maximum(improvements,0.0))),
                          np.finfo(float).eps)
    quality=np.clip(improvements/improvement_scale,0.0,1.0)
    columns=np.column_stack([
        (np.asarray(item.predicted_response,dtype=float).reshape(-1)-current)
        /band for item in ranked])
    norms=np.linalg.norm(columns,axis=0)
    units=np.zeros_like(columns)
    nonzero=norms>64.0*np.finfo(float).eps
    units[:,nonzero]=columns[:,nonzero]/norms[nonzero]
    selected=[0];available=set(range(1,len(ranked)))
    while available and len(selected)<limit:
        def merit(index):
            novelty=_response_subspace_novelty(units,index,selected)
            score=(1.0-weight)*quality[index]+weight*novelty
            return (score,-float(ranked[index].predicted_max_band_ratio),
                    -len(ranked[index].candidate_indices),-index)
        chosen=max(available,key=merit)
        selected.append(chosen);available.remove(chosen)
    return tuple(ranked[index] for index in selected)


def connected_graph_front_beam(*, current_response, response_target,
        response_band, candidate_response_delta, adjacency, seed_indices,
        exclusion_groups=None, maximum_size, maximum_components=2,
        beam_width=64, proposal_limit=8,
        response_novelty_weight=0.0, return_result=False,
        regularization_change: Callable[[tuple[int, ...]],float] | None=None,
        is_valid: Callable[[tuple[int, ...]],bool] | None=None):
    """Rank connected, QR-seeded binary front clusters without exact solves.

    Every component must begin at an ACA/QR/TSVD seed and may then grow only
    through a face-neighbour node.  This replaces arbitrary-cardinality macro
    bundles while retaining more than one physical end through
    ``maximum_components``.  Returned responses use the raw analytic columns.
    """
    current=np.asarray(current_response,dtype=float).reshape(-1)
    target=np.asarray(response_target,dtype=float).reshape(-1)
    band=np.asarray(response_band,dtype=float).reshape(-1)
    delta=np.asarray(candidate_response_delta,dtype=float)
    nc=delta.shape[1] if delta.ndim==2 else -1
    if (current.shape!=target.shape or band.shape!=current.shape or
            delta.shape!=(current.size,nc) or nc<1 or len(adjacency)!=nc or
            np.any(band<=0.0)):
        raise ValueError("connected graph-front response arrays are incompatible")
    seeds=np.unique(np.asarray(seed_indices,dtype=np.int64).reshape(-1))
    if seeds.size==0 or np.any(seeds<0) or np.any(seeds>=nc):
        raise ValueError("connected graph-front needs valid QR/TSVD seeds")
    exclusions=(np.full(nc,-1,dtype=np.int64) if exclusion_groups is None else
                np.asarray(exclusion_groups,dtype=np.int64).reshape(-1))
    if exclusions.shape!=(nc,):
        raise ValueError("connected graph-front exclusions must match candidates")
    maximum_size=min(nc,int(maximum_size));maximum_components=int(maximum_components)
    beam_width=int(beam_width);proposal_limit=int(proposal_limit)
    novelty_weight=float(response_novelty_weight)
    if min(maximum_size,maximum_components,beam_width,proposal_limit)<1:
        raise ValueError("connected graph-front limits must be positive")
    if not 0.0<=novelty_weight<=1.0:
        raise ValueError("response_novelty_weight must lie in [0, 1]")
    ratio=lambda value:float(np.max(np.abs((value-target)/band)))

    def compatible(bundle):
        labels=[int(exclusions[value]) for value in bundle
                if int(exclusions[value])>=0]
        return len(labels)==len(set(labels))

    def make(bundle):
        bundle=tuple(sorted(map(int,bundle)))
        response=current+np.sum(delta[:,bundle],axis=1)
        regularization=(0.0 if regularization_change is None else
                        float(regularization_change(bundle)))
        if not np.isfinite(regularization):
            raise ValueError("graph-front regularization must be finite")
        components=_selected_components(bundle,adjacency)
        return (ratio(response)+regularization,bundle,response,
                regularization,components)

    beam=[make((int(seed),)) for seed in seeds]
    beam=sorted(beam,key=lambda value:(value[0],value[1]))[:beam_width]
    proposals=[];seen=set()
    for depth in range(1,maximum_size+1):
        for score,bundle,response,regularization,components in beam:
            if bundle in seen:
                continue
            seen.add(bundle)
            if is_valid is None or bool(is_valid(bundle)):
                proposals.append(GraphFrontProposal(
                    np.asarray(bundle,dtype=np.int64),response,
                    ratio(response),regularization,components))
        if depth==maximum_size:
            break
        expanded={}
        for _,bundle,_,_,components in beam:
            members=set(bundle);front=set()
            for value in bundle:
                front.update(map(int,adjacency[value]))
            # A new disconnected component may start only at a response-space
            # seed.  All later nodes in that component are face-grown.
            if components<maximum_components:
                front.update(map(int,seeds))
            for value in front-members:
                candidate=tuple(sorted(bundle+(int(value),)))
                if not compatible(candidate):
                    continue
                item=make(candidate)
                if item[4]>maximum_components:
                    continue
                previous=expanded.get(candidate)
                if previous is None or item[0]<previous[0]:
                    expanded[candidate]=item
        beam=sorted(expanded.values(),key=lambda value:(value[0],value[1]))[
            :beam_width]
        if not beam:
            break
    proposals.sort(key=lambda value:(
        value.predicted_max_band_ratio+value.regularization_change,
        len(value.candidate_indices),tuple(value.candidate_indices)))
    pool=tuple(proposals)
    selected=_select_response_diverse_proposals(
        pool,current,target,band,proposal_limit=proposal_limit,
        novelty_weight=novelty_weight)
    if not return_result:
        return selected
    return GraphFrontBeamResult(
        selected,
        graph_front_response_diversity(pool,current,band),
        graph_front_response_diversity(selected,current,band))


def update_graph_front_trust(*, budget, minimum_budget, maximum_budget,
                             predicted_ratio, actual_ratio, current_ratio,
                             selected_size, interface_weight,
                             low_agreement=0.25, high_agreement=0.75):
    """Update discrete front radius from exact/predicted model agreement."""
    before=int(budget);minimum=int(minimum_budget);maximum=int(maximum_budget)
    selected_size=int(selected_size);weight=float(interface_weight)
    predicted_reduction=float(current_ratio)-float(predicted_ratio)
    actual_reduction=float(current_ratio)-float(actual_ratio)
    agreement=(actual_reduction/predicted_reduction
               if predicted_reduction>1.0e-12 else 0.0)
    if predicted_reduction<=1.0e-12 or agreement<float(low_agreement):
        after=max(minimum,int(np.floor(0.5*before)))
        next_weight=1.5*weight;action="shrink"
    elif agreement>float(high_agreement) and selected_size>=before:
        after=min(maximum,max(before+1,int(np.ceil(1.5*before))))
        next_weight=weight/1.25;action="expand"
    else:
        after=before;next_weight=weight;action="hold"
    return GraphFrontTrustUpdate(before,after,weight,float(next_weight),
                                 float(agreement),action)


__all__=[
    "GraphFrontProposal","GraphFrontTrustUpdate",
    "GraphFrontDiversityDiagnostics","GraphFrontBeamResult",
    "best_admissible_singleton",
    "minimax_driving_potential","candidate_face_adjacency",
    "binary_graph_interface_energy",
    "ngsolve_facet_measure_graph",
    "terminal_l1_curvature_energy","graph_front_response_diversity",
    "connected_graph_front_beam",
    "update_graph_front_trust",
]
