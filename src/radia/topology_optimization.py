"""Cubit material topology optimization driven by a linearized VIM system."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class VIMLinearization:
    state: np.ndarray
    response: np.ndarray
    state_jacobian: np.ndarray
    response_jacobian: np.ndarray


@dataclass(frozen=True)
class VIMOperatorLinearization:
    matrix: np.ndarray
    rhs: np.ndarray
    matrix_jacobian: np.ndarray
    rhs_jacobian: np.ndarray


@dataclass(frozen=True)
class ChargeGramLinearization:
    matrix: np.ndarray
    jacobian: np.ndarray


@dataclass(frozen=True)
class ProductionGetTrafoDisplacements:
    """Production ChargeGram node velocities sampled from NGSolve fields."""
    family: str
    cell: tuple[np.ndarray, ...]
    face: tuple[np.ndarray, ...]


@dataclass(frozen=True)
class ProductionVIMLinearization:
    """NGSolve/production ChargeGram geometry tangent before response solve."""
    geometry: ProductionGetTrafoDisplacements
    charge_gram: ChargeGramLinearization
    operator: VIMOperatorLinearization


@dataclass(frozen=True)
class ChargeGramDirectionalOperators:
    """Base ChargeGram and analytic directional H-matrix operators."""
    matrix: object
    jacobian: tuple[object, ...]


@dataclass(frozen=True)
class VIMMatrixFreeLinearization:
    """VIM operator and its full-product-rule directional actions."""
    mass: object
    charge_map: object
    charge_gram: object
    dmass: tuple[object, ...]
    dcharge_map: tuple[object, ...]
    dcharge_gram: tuple[object, ...]
    inv_chi: float

    @property
    def shape(self):
        n=int(self.charge_map.shape[1])
        return (n,n)

    def matvec(self, x):
        x=np.asarray(x,dtype=float).reshape(-1)
        bx=np.asarray(self.charge_map@x).reshape(-1)
        return (self.inv_chi*np.asarray(self.mass@x).reshape(-1)
                +np.asarray(self.charge_map.T@self.charge_gram.matvec_sym(bx)).reshape(-1))

    def directional_matvec(self, mode, x):
        """Apply dA[mode] without materializing dG or dA."""
        x=np.asarray(x,dtype=float).reshape(-1); k=int(mode)
        B=self.charge_map; dB=self.dcharge_map[k]
        bx=np.asarray(B@x).reshape(-1); dbx=np.asarray(dB@x).reshape(-1)
        left=np.asarray(dB.T@self.charge_gram.matvec_sym(bx)).reshape(-1)
        middle=np.asarray(B.T@self.dcharge_gram[k].matvec_sym(bx)).reshape(-1)
        right=np.asarray(B.T@self.charge_gram.matvec_sym(dbx)).reshape(-1)
        return self.inv_chi*np.asarray(self.dmass[k]@x).reshape(-1)+left+middle+right

    def as_scipy_linear_operator(self, mode=None):
        from scipy.sparse.linalg import LinearOperator
        action=self.matvec if mode is None else lambda x:self.directional_matvec(mode,x)
        return LinearOperator(self.shape,matvec=action,rmatvec=action,dtype=float)


@dataclass(frozen=True)
class VIMAdjointGradient:
    state: np.ndarray
    response: np.ndarray
    objective: float
    gradient: np.ndarray
    state_solver_info: int
    adjoint_solver_info: int
    peak_directional_operators: int = 1


def sample_production_gettrafo_displacements(fes, displacement_modes, charge_basis,
                                              *, family=None):
    """Evaluate deformation modes on the exact nodes stored by ChargeGram.

    NGSolve remains responsible for evaluating each vector field.  This
    adapter only reshapes the already constructed production geometry arrays;
    it does not evaluate reference finite-element shapes or orientations.
    """
    modes=tuple(displacement_modes); cb=charge_basis; mesh=fes.mesh
    if family is None:
        if "face_type" in cb: family="wedge"
        elif "vV" in cb: family="tet"
        elif "cell_nodes" in cb: family="hex"
        else: raise ValueError("cannot infer ChargeGram element family")
    family=str(family).lower()
    if family=="tet":
        cell_nodes=tuple(np.asarray(x,dtype=float) for x in cb["vV"])
        face_nodes=tuple(np.asarray(x,dtype=float) for x in cb["bV"])
    elif family=="hex":
        ncell=len(cb.get("cell_charges", ())) or int(cb.get("n_el", 0))
        raw=np.asarray(cb["cell_nodes"],dtype=float).reshape(-1,27,3)
        cell_nodes=tuple(raw)
        face_nodes=tuple(np.asarray(cb["face_nodes"],dtype=float).reshape(-1,9,3))
        if ncell and len(cell_nodes)!=ncell: raise ValueError("HEX cell-node count mismatch")
    elif family=="wedge":
        cell_nodes=tuple(np.asarray(cb["cell_nodes"],dtype=float).reshape(-1,18,3))
        padded=np.asarray(cb["face_nodes"],dtype=float).reshape(-1,9,3)
        types=np.asarray(cb["face_type"],dtype=int)
        if len(padded)!=len(types): raise ValueError("WEDGE face-node count mismatch")
        face_nodes=tuple(padded[i,:6] if types[i]==0 else padded[i] for i in range(len(types)))
    else:
        raise NotImplementedError(f"production GetTrafo sampling does not support {family}")

    def sample(hosts):
        result=[]
        for nodes in hosts:
            values=[]
            for field in modes:
                values.append(np.array([field(mesh(*point)) for point in nodes],dtype=float))
            result.append(np.asarray(values,dtype=float).reshape(len(modes),len(nodes),3))
        return tuple(result)
    return ProductionGetTrafoDisplacements(family,sample(cell_nodes),sample(face_nodes))


def _materialize_charge_gram(charge_gram):
    n=int(charge_gram.ndof())
    return np.array([[charge_gram.entry(i,j) for j in range(n)] for i in range(n)])


def linearize_production_charge_gram(charge_gram, geometry, *, charge_map=None):
    """Dispatch complete production ``dG`` for TET/HEX/WEDGE.

    ``geometry`` must be sampled by
    :func:`sample_production_gettrafo_displacements`.  The production C++
    kernels own all singular and non-singular integration.  No finite
    difference is used here.
    """
    family=geometry.family
    if not geometry.cell and not geometry.face: raise ValueError("empty production geometry")
    q=(geometry.cell or geometry.face)[0].shape[0]
    cells=np.stack(geometry.cell,axis=1) if geometry.cell else None
    if family=="wedge":
        # The native mixed-face contract is padded to nine Q2 slots.  Tri
        # faces use the first six slots; unused velocities remain exactly zero.
        faces=np.zeros((q,len(geometry.face),9,3))
        for host,modes in enumerate(geometry.face): faces[:,host,:modes.shape[1]]=modes
    else:
        faces=np.stack(geometry.face,axis=1) if geometry.face else None
    if family=="tet":
        dG,dB=production_tet_charge_gram_derivatives(
            charge_gram,cells,faces,charge_map=charge_map)
        if charge_map is None: dB=None
        else: dB=np.asarray([x.toarray() if hasattr(x,"toarray") else x for x in dB])
    elif family=="hex":
        dG=np.asarray([charge_gram.hex_charge_gram_directional_derivative(
            np.ascontiguousarray(c),np.ascontiguousarray(f)) for c,f in zip(cells,faces)])
        dB=None
    elif family=="wedge":
        dG=np.asarray([charge_gram.wedge_charge_gram_directional_derivative(
            np.ascontiguousarray(c),np.ascontiguousarray(f)) for c,f in zip(cells,faces)])
        dB=None
    else: raise NotImplementedError(f"production ChargeGram derivative does not support {family}")
    return ChargeGramLinearization(_materialize_charge_gram(charge_gram),dG),dB


def linearize_production_charge_gram_matrix_free(charge_gram, geometry, *,
                                                  eps=1e-10, leaf=64, eta=2.0):
    """Build analytic directional ChargeGram H-matrices without dense ``dG``.

    Each returned derivative owns an ACA/H-matrix representation and exposes
    only entry probes and symmetric matvec.  Finite differences are not used.
    """
    family=str(geometry.family).lower()
    if family not in {"tet","hex","wedge"}:
        raise NotImplementedError(f"production ChargeGram derivative does not support {family}")
    if not geometry.cell and not geometry.face: raise ValueError("empty production geometry")
    q=(geometry.cell or geometry.face)[0].shape[0]
    cells=np.stack(geometry.cell,axis=1) if geometry.cell else np.empty((q,0,0,3))
    if family=="wedge":
        faces=np.zeros((q,len(geometry.face),9,3))
        for host,modes in enumerate(geometry.face): faces[:,host,:modes.shape[1]]=modes
    else:
        faces=np.stack(geometry.face,axis=1) if geometry.face else np.empty((q,0,0,3))
    operators=[]
    for cell_velocity,face_velocity in zip(cells,faces):
        operators.append(charge_gram.directional_derivative_operator(
            family,np.ascontiguousarray(cell_velocity),np.ascontiguousarray(face_velocity),
            eps=float(eps),leaf=int(leaf),eta=float(eta)))
    return ChargeGramDirectionalOperators(charge_gram,tuple(operators))


def linearize_vim_operator_matrix_free(mass, charge_map, charge_gram,
                                       *, inv_chi, dmass, dcharge_map,
                                       dcharge_gram):
    """Create matrix-free ``A``/``dA`` actions using the complete product rule."""
    import scipy.sparse as sp
    B=sp.csr_matrix(charge_map); M=sp.csr_matrix(mass)
    dM=tuple(sp.csr_matrix(x) for x in dmass)
    dB=tuple(sp.csr_matrix(x) for x in dcharge_map)
    dG=tuple(dcharge_gram)
    if not (len(dM)==len(dB)==len(dG)): raise ValueError("VIM derivative mode count mismatch")
    if M.shape!=(B.shape[1],B.shape[1]): raise ValueError("mass/charge-map shape mismatch")
    if any(x.shape!=M.shape for x in dM) or any(x.shape!=B.shape for x in dB):
        raise ValueError("VIM derivative shape mismatch")
    return VIMMatrixFreeLinearization(M,B,charge_gram,dM,dB,dG,float(inv_chi))


def assemble_ngsolve_hdiv_shape_tangents(fes, displacement_modes, charge_map):
    """Assemble analytic HDiv mass and Piola-exact charge-map tangents.

    NGSolve evaluates the HDiv trial/test functions and ``Grad(V)``.  Under
    the contravariant Piola transport,
    ``dM = integral u.(DV+DV.T-div(V)I).v``.  Radia's reference-charge map
    has exactly cancelling physical Jacobians, hence its material derivative
    is zero; topology/orientation remains owned by the already assembled map.
    The caller owns ``ngsolve.TaskManager``.
    """
    import ngsolve as ng
    import scipy.sparse as sp
    def csr(matrix):
        rows,cols,values=matrix.COO()
        return sp.csr_matrix((np.asarray(values),(np.asarray(rows),np.asarray(cols))),
                             shape=(matrix.height,matrix.width))
    modes=tuple(displacement_modes)
    u,v=fes.TnT(); base=ng.BilinearForm(fes); base+=u*v*ng.dx; base.Assemble()
    mass=csr(base.mat); tangents=[]
    for velocity in modes:
        DV=ng.Grad(velocity)
        form=ng.BilinearForm(fes)
        form+=ng.InnerProduct(u,(DV+DV.trans-ng.div(velocity)*ng.Id(fes.mesh.dim))*v)*ng.dx
        form.Assemble(); tangents.append(csr(form.mat).toarray())
    B=sp.csr_matrix(charge_map)
    return mass.toarray(),np.asarray(tangents),np.zeros((len(modes),*B.shape))


def linearize_production_vim_from_ngsolve(*, fes, deformation_modes,
                                           charge_basis, charge_gram, charge_map,
                                           applied_coefficients, inv_chi,
                                           dapplied_coefficients=None,
                                           family=None):
    """Close GetTrafo fields through production ``dM,dB,dG`` and VIM ``dA``.

    The caller owns ``ngsolve.TaskManager``.  NGSolve assembles the HDiv mass
    tangent; native Radia kernels differentiate the complete production
    ChargeGram.  Flat TET charge-row Piola rates are applied to the already
    oriented NGSolve charge map, while HEX/WEDGE reference charge maps have
    zero material derivative.
    """
    modes=tuple(deformation_modes)
    geometry=sample_production_gettrafo_displacements(
        fes,modes,charge_basis,family=family)
    mass,dmass,dcharge=assemble_ngsolve_hdiv_shape_tangents(
        fes,modes,charge_map)
    charge,dcharge_native=linearize_production_charge_gram(
        charge_gram,geometry,charge_map=charge_map)
    if dcharge_native is not None: dcharge=dcharge_native
    B=charge_map.toarray() if hasattr(charge_map,"toarray") else np.asarray(charge_map)
    operator=linearize_vim_operator(mass,B,charge.matrix,applied_coefficients,
        inv_chi=inv_chi,dmass=dmass,dcharge_map=dcharge,
        dcharge_gram=charge.jacobian,
        dapplied_coefficients=dapplied_coefficients)
    return ProductionVIMLinearization(geometry,charge,operator)


def linearize_production_vim_matrix_free_from_ngsolve(*, fes, deformation_modes,
            charge_basis, charge_gram, charge_map, inv_chi, family=None,
            eps=1e-10, leaf=64, eta=2.0):
    """Close GetTrafo through matrix-free ``A`` and analytic ``dA`` actions.

    The base and every directional ChargeGram remain H-matrices.  NGSolve
    owns FE assembly/evaluation and the caller owns its ``TaskManager``.
    """
    import scipy.sparse as sp
    modes=tuple(deformation_modes)
    geometry=sample_production_gettrafo_displacements(
        fes,modes,charge_basis,family=family)
    mass,dmass,dcharge=assemble_ngsolve_hdiv_shape_tangents(
        fes,modes,charge_map)
    charge=linearize_production_charge_gram_matrix_free(
        charge_gram,geometry,eps=eps,leaf=leaf,eta=eta)
    B=sp.csr_matrix(charge_map)
    if geometry.family=="tet":
        cells=np.stack(geometry.cell,axis=1); faces=np.stack(geometry.face,axis=1)
        dcharge=tuple(sp.diags(np.asarray(
            charge_gram.tet_charge_map_row_directional_rates(
                np.ascontiguousarray(c),np.ascontiguousarray(f))))@B
            for c,f in zip(cells,faces))
    else:
        dcharge=tuple(sp.csr_matrix(x) for x in dcharge)
    operator=linearize_vim_operator_matrix_free(
        mass,B,charge_gram,inv_chi=inv_chi,dmass=dmass,
        dcharge_map=dcharge,dcharge_gram=charge.jacobian)
    return geometry,charge,operator


def production_vim_rms_adjoint_gradient_streaming(*, fes, deformation_modes,
            charge_basis, charge_gram, charge_map, applied_coefficients,
            inv_chi, response_matrix, family=None, incident_response=None,
            dapplied_coefficients=None, dresponse_matrix=None,
            eps=1e-10, leaf=64, eta=2.0, solve_tolerance=1e-10,
            solve_max_iterations=None):
    """Compute an RMS response gradient with one live derivative H-matrix.

    For each GetTrafo mode, the native analytic ``dG`` H-matrix is built,
    contracted as ``lambda.T @ (db - dA @ state)``, and immediately released.
    Neither dense ``dG``/``dA`` nor a tuple of directional H-matrices exists.
    The caller owns ``ngsolve.TaskManager``.
    """
    import gc
    import scipy.sparse as sp
    from scipy.sparse.linalg import cg
    modes=tuple(deformation_modes); q=len(modes)
    if q==0: raise ValueError("at least one deformation mode is required")
    B=sp.csr_matrix(charge_map); h=np.asarray(applied_coefficients,dtype=float).reshape(-1)
    C=np.atleast_2d(np.asarray(response_matrix,dtype=float))
    if B.shape[1]!=h.size or C.shape[1]!=h.size: raise ValueError("VIM response shape mismatch")
    dh=np.zeros((q,h.size)) if dapplied_coefficients is None else np.asarray(dapplied_coefficients,dtype=float)
    dC=np.zeros((q,*C.shape)) if dresponse_matrix is None else np.asarray(dresponse_matrix,dtype=float)
    if dh.shape!=(q,h.size) or dC.shape!=(q,*C.shape): raise ValueError("response derivative shape mismatch")

    # Geometry mass is mode-independent; assembling it once also keeps all FE
    # orientation/Piola plumbing in NGSolve.
    mass,_,_=assemble_ngsolve_hdiv_shape_tangents(fes,(),B)
    M=sp.csr_matrix(mass)
    base=linearize_vim_operator_matrix_free(M,B,charge_gram,inv_chi=inv_chi,
        dmass=(),dcharge_map=(),dcharge_gram=())
    A=base.as_scipy_linear_operator(); rhs=np.asarray(M@h).reshape(-1)
    state,state_info=cg(A,rhs,rtol=solve_tolerance,atol=0.0,maxiter=solve_max_iterations)
    if state_info!=0: raise RuntimeError(f"matrix-free VIM state solve failed (info={state_info})")
    response=C@state
    if incident_response is not None:
        incident=np.asarray(incident_response,dtype=float).reshape(-1)
        if incident.shape!=response.shape: raise ValueError("incident_response shape mismatch")
        response=response+incident
    objective=float(np.sqrt(np.mean(response**2)))
    if objective==0: raise ValueError("RMS response derivative is undefined at zero")
    response_weight=response/(response.size*objective)
    adjoint,adjoint_info=cg(A,C.T@response_weight,rtol=solve_tolerance,
                            atol=0.0,maxiter=solve_max_iterations)
    if adjoint_info!=0: raise RuntimeError(f"matrix-free VIM adjoint solve failed (info={adjoint_info})")

    gradient=np.empty(q)
    for k,mode in enumerate(modes):
        geometry=sample_production_gettrafo_displacements(
            fes,(mode,),charge_basis,family=family)
        _,dmass,dcharge=assemble_ngsolve_hdiv_shape_tangents(fes,(mode,),B)
        dM=sp.csr_matrix(dmass[0])
        if geometry.family=="tet":
            cv=np.stack(geometry.cell,axis=1)[0]; fv=np.stack(geometry.face,axis=1)[0]
            rates=np.asarray(charge_gram.tet_charge_map_row_directional_rates(
                np.ascontiguousarray(cv),np.ascontiguousarray(fv)))
            dB=sp.diags(rates)@B
        else:
            dB=sp.csr_matrix(dcharge[0])
        directional=linearize_production_charge_gram_matrix_free(
            charge_gram,geometry,eps=eps,leaf=leaf,eta=eta).jacobian[0]
        bx=np.asarray(B@state).reshape(-1); dbx=np.asarray(dB@state).reshape(-1)
        dA_state=(float(inv_chi)*np.asarray(dM@state).reshape(-1)
            +np.asarray(dB.T@charge_gram.matvec_sym(bx)).reshape(-1)
            +np.asarray(B.T@directional.matvec_sym(bx)).reshape(-1)
            +np.asarray(B.T@charge_gram.matvec_sym(dbx)).reshape(-1))
        db=np.asarray(dM@h+M@dh[k]).reshape(-1)
        gradient[k]=float(adjoint@(db-dA_state)+response_weight@(dC[k]@state))
        del directional,geometry,dM,dB
        gc.collect()
    return VIMAdjointGradient(state,response,objective,gradient,
                              int(state_info),int(adjoint_info),1)


def linearize_laplace_pair_gram(points, weights, displacement_modes,
                                relative_weight_derivatives=None):
    """Analytic geometry derivative of the off-diagonal Laplace pair Gram.

    Self/singular panel terms are intentionally excluded; production callers
    must add their analytic self-panel value and derivative separately.
    """
    x=np.asarray(points,dtype=float); w=np.asarray(weights,dtype=float).reshape(-1)
    velocity=np.asarray(displacement_modes,dtype=float)
    if x.ndim!=2 or x.shape[0]!=w.size or velocity.ndim!=3 or velocity.shape[1:]!=x.shape:
        raise ValueError("expected points (n,d), weights (n), displacement_modes (q,n,d)")
    if np.any(w<=0): raise ValueError("weights must be positive")
    q=velocity.shape[0]
    rel=np.zeros((q,w.size)) if relative_weight_derivatives is None else np.asarray(relative_weight_derivatives,dtype=float)
    if rel.shape!=(q,w.size): raise ValueError("relative_weight_derivatives must have shape (q,n)")
    delta=x[:,None,:]-x[None,:,:]; distance=np.linalg.norm(delta,axis=2)
    if np.any(distance[np.triu_indices(w.size,1)]==0): raise ValueError("distinct sample points are required")
    inverse=np.zeros_like(distance); mask=~np.eye(w.size,dtype=bool); inverse[mask]=1/distance[mask]
    gram=(w[:,None]*w[None,:])*inverse/(4*np.pi)
    derivative=np.empty((q,w.size,w.size))
    for k in range(q):
        dv=velocity[k,:,None,:]-velocity[k,None,:,:]
        radial=np.einsum("ijd,ijd->ij",delta,dv)
        kernel_term=np.zeros_like(distance); kernel_term[mask]=-radial[mask]/distance[mask]**3
        derivative[k]=(w[:,None]*w[None,:])/(4*np.pi)*(kernel_term+inverse*(rel[k,:,None]+rel[k,None,:]))
        np.fill_diagonal(derivative[k],0)
    return gram,derivative


def affine_cell_self_energy_shape_derivative(cell_type, nodes, displacement_modes):
    """C++ analytic singular self-term derivative for an affine cell.

    This is the constant physical volume-charge block.  ``displacement_modes``
    are nodal GetTrafo velocities; finite differences are not used here.
    """
    from . import _radia_pybind as _core
    result=_core._AffineCellSelfEnergyShapeDerivative(
        str(cell_type).lower(),np.ascontiguousarray(nodes,dtype=float),
        np.ascontiguousarray(displacement_modes,dtype=float))
    return float(result["value"]),np.asarray(result["derivative"],dtype=float)


def production_hex_volume_self_block_derivatives(charge_gram, host_node_displacements):
    """Differentiate production HEX polynomial volume-charge self blocks.

    ``host_node_displacements[h]`` has shape ``(q,27,3)`` in the exact Q2
    GetTrafo lattice ordering used to construct ``charge_gram``.  Returned
    blocks are row-major NumPy arrays grouped as ``result[host][mode]``.
    """
    result=[]; expected_modes=None
    for host,modes in enumerate(host_node_displacements):
        modes=np.asarray(modes,dtype=float)
        if modes.ndim!=3 or modes.shape[1:]!=(27,3):
            raise ValueError("each HEX host displacement must have shape (q,27,3)")
        if expected_modes is None: expected_modes=modes.shape[0]
        elif modes.shape[0]!=expected_modes: raise ValueError("all HEX hosts must have the same mode count")
        result.append(np.stack([
            np.asarray(charge_gram.hex_volume_self_block_directional_derivative(
                host,np.ascontiguousarray(mode)),dtype=float)
            for mode in modes]))
    return tuple(result)


def production_hex_face_self_block_derivatives(charge_gram, host_node_displacements):
    """Differentiate production HEX quad-face polynomial self blocks.

    Each host displacement array has shape ``(q,9,3)`` in the Q2 face
    GetTrafo lattice ordering used by the production ChargeGram.
    """
    result=[]; expected_modes=None
    for host,modes in enumerate(host_node_displacements):
        modes=np.asarray(modes,dtype=float)
        if modes.ndim!=3 or modes.shape[1:]!=(9,3):
            raise ValueError("each HEX face displacement must have shape (q,9,3)")
        if expected_modes is None: expected_modes=modes.shape[0]
        elif modes.shape[0]!=expected_modes: raise ValueError("all HEX faces must have the same mode count")
        result.append(np.stack([
            np.asarray(charge_gram.hex_face_self_block_directional_derivative(
                host,np.ascontiguousarray(mode)),dtype=float)
            for mode in modes]))
    return tuple(result)


def production_tet_volume_self_block_derivatives(charge_gram, host_vertex_displacements):
    """Differentiate flat production TET volume self blocks analytically."""
    result=[]; expected=None
    for host,modes in enumerate(host_vertex_displacements):
        modes=np.asarray(modes,dtype=float)
        if modes.ndim!=3 or modes.shape[1:]!=(4,3): raise ValueError("each TET host displacement must have shape (q,4,3)")
        if expected is None: expected=modes.shape[0]
        elif modes.shape[0]!=expected: raise ValueError("all TET hosts must have the same mode count")
        result.append(np.stack([np.asarray(charge_gram.tet_volume_self_block_directional_derivative(host,np.ascontiguousarray(mode)),dtype=float) for mode in modes]))
    return tuple(result)


def production_tet_face_self_block_derivatives(charge_gram, host_vertex_displacements):
    """Differentiate flat production TET triangular-face self blocks analytically."""
    result=[]; expected=None
    for host,modes in enumerate(host_vertex_displacements):
        modes=np.asarray(modes,dtype=float)
        if modes.ndim!=3 or modes.shape[1:]!=(3,3): raise ValueError("each TET face displacement must have shape (q,3,3)")
        if expected is None: expected=modes.shape[0]
        elif modes.shape[0]!=expected: raise ValueError("all TET faces must have the same mode count")
        result.append(np.stack([np.asarray(charge_gram.tet_face_self_block_directional_derivative(host,np.ascontiguousarray(mode)),dtype=float) for mode in modes]))
    return tuple(result)


def production_tet_charge_gram_derivatives(charge_gram, cell_vertex_displacements,
                                           face_vertex_displacements, charge_map=None):
    """Differentiate the complete flat production TET ChargeGram analytically.

    Inputs have shapes ``(q,ncell,4,3)`` and ``(q,nface,3,3)``.  The returned
    dense ``dG`` is the exact derivative of production's symmetric row-major
    Gram.  If the current NGSolve-owned sparse/dense ``charge_map`` is supplied,
    ``dB`` is also returned.  Flat TET monomials use fixed host reference
    coordinates, so Piola transport changes each row only by inverse volume or
    surface-area measure; element orientation remains exactly as stored in B.
    """
    cells=np.asarray(cell_vertex_displacements,dtype=float)
    faces=np.asarray(face_vertex_displacements,dtype=float)
    if cells.ndim!=4 or cells.shape[2:]!=(4,3):
        raise ValueError("cell_vertex_displacements must have shape (q,ncell,4,3)")
    if faces.ndim!=4 or faces.shape[2:]!=(3,3):
        raise ValueError("face_vertex_displacements must have shape (q,nface,3,3)")
    if cells.shape[0]!=faces.shape[0]: raise ValueError("cell and face mode counts must match")
    dG=[]; rates=[]
    for c,f in zip(cells,faces):
        cc=np.ascontiguousarray(c); ff=np.ascontiguousarray(f)
        dG.append(np.asarray(charge_gram.tet_charge_gram_directional_derivative(cc,ff),dtype=float))
        rates.append(np.asarray(charge_gram.tet_charge_map_row_directional_rates(cc,ff),dtype=float))
    dG=np.asarray(dG); rates=np.asarray(rates)
    if charge_map is None: return dG,rates
    try:
        import scipy.sparse as sp
        B=sp.csr_matrix(charge_map)
        dB=tuple(sp.diags(r)@B for r in rates)
    except ImportError:
        B=np.asarray(charge_map,dtype=float)
        dB=rates[:,:,None]*B[None,:,:]
    return dG,dB


def production_wedge_volume_self_block_derivatives(charge_gram, host_node_displacements):
    """Differentiate production WEDGE volume self blocks (18-node lattice)."""
    result=[]; expected_modes=None
    for host,modes in enumerate(host_node_displacements):
        modes=np.asarray(modes,dtype=float)
        if modes.ndim!=3 or modes.shape[1:]!=(18,3):
            raise ValueError("each WEDGE host displacement must have shape (q,18,3)")
        if expected_modes is None: expected_modes=modes.shape[0]
        elif modes.shape[0]!=expected_modes: raise ValueError("all WEDGE hosts must have the same mode count")
        result.append(np.stack([np.asarray(charge_gram.wedge_volume_self_block_directional_derivative(host,np.ascontiguousarray(mode)),dtype=float) for mode in modes]))
    return tuple(result)


def production_wedge_face_self_block_derivatives(charge_gram, host_node_displacements):
    """Differentiate mixed WEDGE face blocks (six-node tri or nine-node quad)."""
    result=[]; expected_modes=None
    for host,modes in enumerate(host_node_displacements):
        modes=np.asarray(modes,dtype=float)
        if modes.ndim!=3 or modes.shape[1] not in (6,9) or modes.shape[2]!=3:
            raise ValueError("each WEDGE face displacement must have shape (q,6,3) or (q,9,3)")
        if expected_modes is None: expected_modes=modes.shape[0]
        elif modes.shape[0]!=expected_modes: raise ValueError("all WEDGE faces must have the same mode count")
        result.append(np.stack([np.asarray(charge_gram.wedge_face_self_block_directional_derivative(host,np.ascontiguousarray(mode)),dtype=float) for mode in modes]))
    return tuple(result)


def production_wedge_charge_gram_derivatives(charge_gram, cell_node_displacements,
                                              face_node_displacements):
    """Differentiate the complete production WEDGE ChargeGram for each mode.

    Triangular face modes occupy the first six of the common nine-node face
    slots; their final three rows must be zero.  The returned array is
    ``(q,ncharge,ncharge)`` and C-contiguous/row-major.
    """
    cells=np.asarray(cell_node_displacements,dtype=float)
    faces=np.asarray(face_node_displacements,dtype=float)
    if cells.ndim!=4 or cells.shape[2:]!=(18,3):
        raise ValueError("cell_node_displacements must have shape (q,ncell,18,3)")
    if faces.ndim!=4 or faces.shape[2:]!=(9,3):
        raise ValueError("face_node_displacements must have shape (q,nface,9,3)")
    if cells.shape[0]!=faces.shape[0]:raise ValueError("cell and face mode counts must match")
    return np.asarray([np.asarray(charge_gram.wedge_charge_gram_directional_derivative(
        np.ascontiguousarray(c),np.ascontiguousarray(f)),dtype=float)
        for c,f in zip(cells,faces)])


def linearize_laplace_charge_gram(points, weights, displacement_modes, *,
                                  relative_weight_derivatives=None,
                                  self_cell_types=None, self_nodes=None,
                                  self_node_displacements=None):
    """Assemble the full analytic Laplace Gram tangent, including its diagonal.

    Off-diagonal pairs use the analytic kernel derivative.  Every diagonal is
    supplied by the C++ singular affine-cell kernel.  The node displacement
    arrays are the values of caller-owned NGSolve/GetTrafo deformation modes;
    this function deliberately does not reconstruct FE bases or Piola maps.
    """
    gram,jac=linearize_laplace_pair_gram(
        points,weights,displacement_modes,relative_weight_derivatives)
    n=gram.shape[0]; q=jac.shape[0]
    if self_cell_types is None or self_nodes is None or self_node_displacements is None:
        raise ValueError("analytic self_cell_types, self_nodes, and self_node_displacements are required")
    if len(self_cell_types)!=n or len(self_nodes)!=n or len(self_node_displacements)!=n:
        raise ValueError("one analytic self-cell description is required per Gram diagonal")
    for i,(kind,nodes,modes) in enumerate(zip(self_cell_types,self_nodes,self_node_displacements)):
        modes=np.asarray(modes,dtype=float)
        if modes.ndim!=3 or modes.shape[0]!=q:
            raise ValueError("each self-node displacement array must have shape (q,nnode,3)")
        value,derivative=affine_cell_self_energy_shape_derivative(kind,nodes,modes)
        gram[i,i]=value; jac[:,i,i]=derivative
    return ChargeGramLinearization(gram,jac)


def linearize_vim_operator(mass, charge_map, charge_gram, applied_coefficients, *,
                           inv_chi, dmass, dcharge_map, dcharge_gram,
                           dapplied_coefficients=None) -> VIMOperatorLinearization:
    """Differentiate ``A=inv_chi*M+B.T*G*B`` and ``b=M*h`` analytically."""
    M=np.asarray(mass); B=np.asarray(charge_map); G=np.asarray(charge_gram); h=np.asarray(applied_coefficients).reshape(-1)
    dM=np.asarray(dmass); dB=np.asarray(dcharge_map); dG=np.asarray(dcharge_gram)
    if M.ndim!=2 or M.shape[0]!=M.shape[1] or B.ndim!=2 or B.shape[1]!=M.shape[0] or G.shape!=(B.shape[0],B.shape[0]) or h.size!=M.shape[0]:
        raise ValueError("incompatible VIM operator shapes")
    q=dM.shape[0]
    if dM.shape!=(q,*M.shape) or dB.shape!=(q,*B.shape) or dG.shape!=(q,*G.shape): raise ValueError("VIM derivative shape mismatch")
    dh=np.zeros((q,h.size),dtype=np.result_type(M,h)) if dapplied_coefficients is None else np.asarray(dapplied_coefficients)
    if dh.shape!=(q,h.size): raise ValueError("dapplied_coefficients shape mismatch")
    A=float(inv_chi)*M+B.T@G@B; rhs=M@h
    dA=np.empty((q,*A.shape),dtype=np.result_type(M,B,G)); db=np.empty((q,h.size),dtype=np.result_type(M,h))
    for k in range(q):
        dA[k]=float(inv_chi)*dM[k]+dB[k].T@G@B+B.T@dG[k]@B+B.T@G@dB[k]
        db[k]=dM[k]@h+M@dh[k]
    return VIMOperatorLinearization(A,rhs,dA,db)


def linearize_vim_system(A, b, C, dA, db=None, dC=None) -> VIMLinearization:
    """Analytically linearize ``A(rho)m=b(rho), y=C(rho)m`` by design cell."""
    A=np.asarray(A); b=np.asarray(b); C=np.asarray(C); dA=np.asarray(dA)
    if A.ndim!=2 or A.shape[0]!=A.shape[1]: raise ValueError("A must be square")
    n=A.shape[0]
    if b.shape not in {(n,),(n,1)}: raise ValueError("b must match A")
    b=b.reshape(n); C=np.atleast_2d(C)
    if C.shape[1]!=n: raise ValueError("C must have one column per state")
    if dA.ndim!=3 or dA.shape[1:]!=(n,n): raise ValueError("dA must have shape (cells,n,n)")
    cells=dA.shape[0]
    db=np.zeros((cells,n),dtype=np.result_type(A,b)) if db is None else np.asarray(db)
    if db.shape!=(cells,n): raise ValueError("db must have shape (cells,n)")
    dC=np.zeros((cells,*C.shape),dtype=np.result_type(A,C)) if dC is None else np.asarray(dC)
    if dC.shape!=(cells,*C.shape): raise ValueError("dC must have shape (cells,outputs,n)")
    state=np.linalg.solve(A,b)
    rhs=db-np.einsum("cij,j->ci",dA,state)
    state_jacobian=np.linalg.solve(A,rhs.T).T
    response=C@state
    response_jacobian=np.einsum("oi,ci->oc",C,state_jacobian)+np.einsum("coi,i->oc",dC,state)
    return VIMLinearization(state,response,state_jacobian,response_jacobian)


@dataclass(frozen=True)
class LPUpdate:
    density: np.ndarray
    delta: np.ndarray
    predicted_objective: float
    status: str
    iterations: int


@dataclass(frozen=True)
class TopologyOptimizationResult:
    density: np.ndarray
    history: tuple[dict, ...]
    converged: bool


def solve_lp_update(density, objective_gradient, cell_volumes, volume_max, *,
                    move_limit=0.2, A_ub=None, b_ub=None) -> LPUpdate:
    """Solve one bounded linear material-distribution update with HiGHS."""
    from scipy.optimize import linprog
    density=np.asarray(density,dtype=float).reshape(-1)
    gradient=np.asarray(objective_gradient,dtype=float).reshape(-1)
    volumes=np.asarray(cell_volumes,dtype=float).reshape(-1)
    if not (density.size==gradient.size==volumes.size): raise ValueError("cell vectors must have equal length")
    if np.any((density<0)|(density>1)) or np.any(volumes<=0): raise ValueError("density must be in [0,1] and volumes positive")
    if not (0<move_limit<=1): raise ValueError("move_limit must be in (0,1]")
    lower=np.maximum(0,density-move_limit); upper=np.minimum(1,density+move_limit)
    rows=[volumes]; limits=[float(volume_max)]
    if A_ub is not None:
        extra=np.atleast_2d(np.asarray(A_ub,dtype=float)); rhs=np.asarray(b_ub,dtype=float).reshape(-1)
        if extra.shape!=(rhs.size,density.size): raise ValueError("A_ub/b_ub shape mismatch")
        rows.extend(extra); limits.extend(rhs)
    result=linprog(gradient,A_ub=np.vstack(rows),b_ub=np.asarray(limits),bounds=list(zip(lower,upper)),method="highs")
    if not result.success: raise RuntimeError(f"topology LP failed: {result.message}")
    return LPUpdate(result.x,result.x-density,float(gradient@result.x),str(result.message),int(result.nit))


def optimize_vim_lp(initial_density, cell_volumes, volume_fraction, linearize, *,
                    objective_weights, move_limit=0.2, max_iterations=30,
                    density_tolerance=1e-3):
    """Run sequential VIM linearization and bounded LP material updates."""
    density=np.asarray(initial_density,dtype=float).reshape(-1)
    volumes=np.asarray(cell_volumes,dtype=float).reshape(-1)
    weights=np.asarray(objective_weights,dtype=float).reshape(-1)
    if density.size!=volumes.size: raise ValueError("initial_density and cell_volumes must match")
    if not (0<volume_fraction<=1): raise ValueError("volume_fraction must be in (0,1]")
    history=[]; converged=False; volume_max=float(volume_fraction*np.sum(volumes))
    for iteration in range(int(max_iterations)):
        model=linearize(density.copy())
        jacobian=np.asarray(model.response_jacobian)
        response=np.asarray(model.response).reshape(-1)
        if jacobian.shape!=(weights.size,density.size): raise ValueError("linearized response shape mismatch")
        gradient=np.real(weights@jacobian)
        update=solve_lp_update(density,gradient,volumes,volume_max,move_limit=move_limit)
        change=float(np.max(np.abs(update.delta)))
        history.append({"iteration":iteration,"objective":float(np.real(weights@response)),
                        "volume":float(volumes@update.density),"max_density_change":change})
        density=update.density
        if change<=density_tolerance: converged=True; break
    return TopologyOptimizationResult(density,tuple(history),converged)


def write_cubit_density_journal(path, element_ids, density, *, threshold=0.5,
                                solid_block=1001, void_block=1002):
    """Write deterministic Cubit block assignment commands for a density field."""
    ids=np.asarray(element_ids,dtype=np.int64).reshape(-1); rho=np.asarray(density,dtype=float).reshape(-1)
    if ids.size!=rho.size or ids.size==0: raise ValueError("element_ids and density must be non-empty and equal length")
    if len(np.unique(ids))!=ids.size or np.any(ids<=0): raise ValueError("element IDs must be unique positive integers")
    if np.any((rho<0)|(rho>1)): raise ValueError("density must be in [0,1]")
    solid=ids[rho>=threshold]; void=ids[rho<threshold]
    lines=["# Radia VIM linearized topology material assignment", "set echo off"]
    for name,values,block in (("radia_topopt_solid",solid,solid_block),("radia_topopt_void",void,void_block)):
        lines.append(f"group '{name}' add hex {' '.join(map(str,values))}" if values.size else f"group '{name}'")
        lines.append(f"block {int(block)} hex in group '{name}'")
    destination=Path(path); destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text("\n".join(lines)+"\n",encoding="utf-8")
    return {"path":str(destination),"solid_count":int(solid.size),"void_count":int(void.size),"threshold":float(threshold)}


__all__=["VIMLinearization","VIMOperatorLinearization","ChargeGramLinearization",
         "ChargeGramDirectionalOperators","VIMMatrixFreeLinearization","VIMAdjointGradient",
         "ProductionGetTrafoDisplacements","ProductionVIMLinearization","LPUpdate",
         "TopologyOptimizationResult","sample_production_gettrafo_displacements",
         "assemble_ngsolve_hdiv_shape_tangents","linearize_production_charge_gram",
         "linearize_production_charge_gram_matrix_free","linearize_vim_operator_matrix_free",
         "linearize_production_vim_matrix_free_from_ngsolve",
         "production_vim_rms_adjoint_gradient_streaming",
         "linearize_production_vim_from_ngsolve","linearize_laplace_pair_gram",
         "affine_cell_self_energy_shape_derivative",
         "production_hex_volume_self_block_derivatives","production_hex_face_self_block_derivatives",
         "production_tet_volume_self_block_derivatives","production_tet_face_self_block_derivatives",
         "production_tet_charge_gram_derivatives",
         "production_wedge_volume_self_block_derivatives","production_wedge_face_self_block_derivatives",
         "production_wedge_charge_gram_derivatives","linearize_laplace_charge_gram",
         "linearize_vim_operator","linearize_vim_system","solve_lp_update",
         "optimize_vim_lp","write_cubit_density_journal"]
