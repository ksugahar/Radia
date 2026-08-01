"""Cubit material topology optimization driven by a linearized VIM system."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
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


@dataclass(frozen=True)
class VIMFunctionalShapeJacobian:
    """State, physical responses, and GetTrafo Jacobian for many functionals."""
    state: np.ndarray
    response: np.ndarray
    response_jacobian: np.ndarray
    state_iterations: int
    adjoint_iterations: tuple[int, ...]


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


def assemble_ngsolve_hdiv_shape_tangents(fes, displacement_modes, charge_map,
                                          *, sparse=False):
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
        form.Assemble()
        tangent=csr(form.mat)
        tangents.append(tangent if sparse else tangent.toarray())
    B=sp.csr_matrix(charge_map)
    if sparse:
        return mass,tuple(tangents),tuple(
            sp.csr_matrix(B.shape,dtype=float) for _ in modes)
    return mass.toarray(),np.asarray(tangents),np.zeros((len(modes),*B.shape))


def assemble_ngsolve_hdiv_linear_form_shape_tangents(fes, coefficient,
        displacement_modes, *, bonus_intorder=None, base_values=None):
    """Assemble ``int coefficient.u`` and its analytic material derivatives.

    Contravariant Piola transport cancels the volume Jacobian, leaving
    ``coefficient.(DV*u)``.  A spatially varying coefficient contributes its
    material derivative ``(Grad(coefficient)*V).u`` as well.  This closes both
    fixed coil fields and dipole-reciprocity observation loads without finite
    differences.  Caller owns ``ngsolve.TaskManager``.
    """
    import ngsolve as ng
    modes=tuple(displacement_modes);v=fes.TestFunction()
    measure=(ng.dx if bonus_intorder is None else
             ng.dx(bonus_intorder=int(bonus_intorder)))
    if base_values is None:
        base=ng.LinearForm(fes);base+=ng.InnerProduct(coefficient,v)*measure;base.Assemble()
        values=np.asarray(base.vec.FV().NumPy(),dtype=float).copy()
    else:
        values=np.asarray(base_values,dtype=float).reshape(-1).copy()
        if values.shape!=(fes.ndof,):
            raise ValueError("base_values must match fes.ndof")
    tangents=[]
    try:
        coordinates=(ng.x,ng.y,ng.z)[:fes.mesh.dim]
        partials=[coefficient.Diff(coordinate) for coordinate in coordinates]
        gradient=ng.CoefficientFunction(tuple(
            partials[column][row] for row in range(fes.mesh.dim)
            for column in range(fes.mesh.dim)),dims=(fes.mesh.dim,fes.mesh.dim))
    except Exception as error:
        raise ValueError("coefficient must expose an NGSolve spatial gradient") from error
    # GridFunction modes produced by the production elastic extension share a
    # VectorH1 space.  The tangent is linear in that velocity, so assemble one
    # rectangular operator and apply it to every mode.  The former loop built
    # the same element quadrature separately for every mode (and was especially
    # costly for a many-segment analytic coil CoefficientFunction).
    mode_space=(getattr(modes[0],"space",None) if modes else None)
    # With automatic quadrature, NGSolve can choose a different rule for the
    # equivalent mixed BilinearForm than for the base LinearForm.  Batch only
    # when the caller fixes the bonus order, which makes the two discrete
    # functionals identical.  Otherwise retain the exact direct assembly.
    batched=(bonus_intorder is not None and mode_space is not None and
             all(getattr(mode,"space",None) is mode_space for mode in modes) and
             all(hasattr(mode,"vec") for mode in modes))
    if batched:
        velocity=mode_space.TrialFunction()
        form=ng.BilinearForm(trialspace=mode_space,testspace=fes,
                             check_unused=False)
        form+=ng.InnerProduct(
            gradient*velocity+ng.Grad(velocity).trans*coefficient,v)*measure
        form.Assemble()
        work=ng.GridFunction(fes).vec
        for mode in modes:
            work.data=form.mat*mode.vec
            tangents.append(np.asarray(work.FV().NumPy(),dtype=float).copy())
    else:
        for velocity in modes:
            DV=ng.Grad(velocity)
            form=ng.LinearForm(fes)
            form+=ng.InnerProduct(gradient*velocity+DV.trans*coefficient,v)*measure
            form.Assemble()
            tangents.append(np.asarray(form.vec.FV().NumPy(),dtype=float).copy())
    return values,np.asarray(tangents,dtype=float).reshape(len(modes),fes.ndof)


def assemble_ngsolve_hdiv_mass_shape_contractions(fes, displacement_modes,
        left_vectors, right_vector, *, bonus_intorder=None):
    """Assemble all ``left.T @ dM[mode] @ right`` contractions once.

    The HDiv mass shape derivative is trilinear in the left state, right state,
    and VectorH1 deformation.  With both HDiv states fixed it is a LinearForm on
    the shared deformation space, so one assembly per left state evaluates all
    modes and avoids building one full sparse HDiv matrix per mode.  NGSolve
    remains the owner of Piola, orientation, quadrature, and local/global
    transformations.  The caller owns ``ngsolve.TaskManager``.
    """
    import ngsolve as ng
    modes=tuple(displacement_modes)
    left=np.atleast_2d(np.asarray(left_vectors,dtype=float))
    right=np.asarray(right_vector,dtype=float).reshape(-1)
    if not modes: return np.empty((left.shape[0],0),dtype=float)
    if left.shape[1]!=fes.ndof or right.shape!=(fes.ndof,):
        raise ValueError("mass shape contraction vectors must match fes.ndof")
    mode_space=getattr(modes[0],"space",None)
    if (mode_space is None or any(getattr(mode,"space",None) is not mode_space
                                  for mode in modes) or
            any(not hasattr(mode,"vec") for mode in modes)):
        raise ValueError("mass shape contractions require GridFunction modes in one shared space")
    # NGSolve 6.2.2604's cross-space LinearForm path does not apply the same
    # local HDiv transform on tensor-product cells as its HDiv BilinearForm.
    # Keep the optimized contraction path on the affine TET production model;
    # retain the already verified sparse dM assembly for HEX/WEDGE parity.
    if any(len(element.vertices)!=4 for element in fes.mesh.Elements(ng.VOL)):
        import scipy.sparse as sp
        _,dmass,_=assemble_ngsolve_hdiv_shape_tangents(
            fes,modes,sp.csr_matrix((0,fes.ndof)),sparse=True)
        return np.asarray([[row@(matrix@right) for matrix in dmass]
                           for row in left],dtype=float)
    right_field=ng.GridFunction(fes)
    right_field.vec.FV().NumPy()[:]=right
    velocity=mode_space.TestFunction();DV=ng.Grad(velocity)
    measure=(ng.dx if bonus_intorder is None else
             ng.dx(bonus_intorder=int(bonus_intorder)))
    result=np.empty((left.shape[0],len(modes)),dtype=float)
    left_field=ng.GridFunction(fes)
    for row,coefficients in enumerate(left):
        left_field.vec.FV().NumPy()[:]=coefficients
        form=ng.LinearForm(mode_space)
        form+=ng.InnerProduct(left_field,
            (DV+DV.trans-ng.div(velocity)*ng.Id(fes.mesh.dim))*right_field)*measure
        form.Assemble()
        for column,mode in enumerate(modes):
            result[row,column]=float(form.vec.InnerProduct(mode.vec))
    return result


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
        fes,modes,charge_map,sparse=True)
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
    mass,_,_=assemble_ngsolve_hdiv_shape_tangents(fes,(),B,sparse=True)
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

    geometry=sample_production_gettrafo_displacements(fes,modes,charge_basis,family=family)
    cells=np.stack(geometry.cell,axis=1)
    if geometry.family=="wedge":
        faces=np.zeros((q,len(geometry.face),9,3))
        for host,values in enumerate(geometry.face): faces[:,host,:values.shape[1]]=values
    else: faces=np.stack(geometry.face,axis=1)
    bx=np.asarray(B@state).reshape(-1); badjoint=np.asarray(B@adjoint).reshape(-1)
    gram_terms=np.asarray(charge_gram.directional_derivative_contractions(
        geometry.family,np.ascontiguousarray(cells),np.ascontiguousarray(faces),
        np.ascontiguousarray(badjoint),np.ascontiguousarray(bx)))
    Gbx=np.asarray(charge_gram.matvec_sym(bx)).reshape(-1)
    gradient=np.empty(q)
    for k,mode in enumerate(modes):
        _,dmass,dcharge=assemble_ngsolve_hdiv_shape_tangents(
            fes,(mode,),B,sparse=True)
        dM=sp.csr_matrix(dmass[0])
        if geometry.family=="tet":
            rates=np.asarray(charge_gram.tet_charge_map_row_directional_rates(
                np.ascontiguousarray(cells[k]),np.ascontiguousarray(faces[k])))
            dB=sp.diags(rates)@B
        else:
            dB=sp.csr_matrix(dcharge[0])
        dbx=np.asarray(dB@state).reshape(-1);dbadjoint=np.asarray(dB@adjoint).reshape(-1)
        db=np.asarray(dM@h+M@dh[k]).reshape(-1)
        gradient[k]=float(adjoint@db-float(inv_chi)*adjoint@np.asarray(dM@state).reshape(-1)
            -dbadjoint@Gbx-badjoint@np.asarray(charge_gram.matvec_sym(dbx)).reshape(-1)
            -gram_terms[k]+response_weight@(dC[k]@state))
    return VIMAdjointGradient(state,response,objective,gradient,
                              int(state_info),int(adjoint_info),0)


def production_vim_functional_shape_jacobian_streaming(*, fes,
        deformation_modes, charge_basis, charge_gram, charge_map,
        inv_chi, rhs, response_matrix, rhs_jacobian=None,
        dresponse_matrix=None, response_observations=None,
        response_weights=None, family=None, incident_response=None,
        dincident_response=None, solve_tolerance=1e-9,
        solve_max_iterations=5000, mass_riesz=True):
    """Differentiate many linear accelerator-field functionals at once.

    The state and all response adjoints share one row-major native solve.
    Directional ``dG`` terms stay on the H-matrix cluster tree; for each
    response, native support-pruned contractions evaluate every deformation
    mode without materializing a derivative matrix.  For flat TET geometry,
    ``response_observations`` and row-major vector ``response_weights`` route
    ``dC`` through the exact native configured-field derivative.  An explicit
    ``dresponse_matrix`` remains available for other element families.
    """
    import scipy.sparse as sp
    modes=tuple(deformation_modes);q=len(modes)
    if q==0: raise ValueError("at least one deformation mode is required")
    B=sp.csr_matrix(charge_map);n=B.shape[1]
    b=np.asarray(rhs,dtype=float).reshape(-1)
    C=np.atleast_2d(np.asarray(response_matrix,dtype=float));nout=C.shape[0]
    if b.shape!=(n,) or C.shape[1]!=n:
        raise ValueError("functional shape RHS/response matrix mismatch")
    db=(np.zeros((q,n)) if rhs_jacobian is None else
        np.asarray(rhs_jacobian,dtype=float))
    native_response=(response_observations is not None or
                     response_weights is not None)
    if native_response and (response_observations is None or
                            response_weights is None):
        raise ValueError("response_observations and response_weights must be supplied together")
    if native_response and dresponse_matrix is not None:
        raise ValueError("choose native response sampling or dresponse_matrix, not both")
    dC=(None if native_response else
        (np.zeros((q,nout,n)) if dresponse_matrix is None else
         np.asarray(dresponse_matrix,dtype=float)))
    if db.shape!=(q,n) or (dC is not None and dC.shape!=(q,nout,n)):
        raise ValueError("functional shape RHS/response derivatives mismatch")
    if native_response:
        response_observations=np.asarray(response_observations,dtype=float)
        response_weights=np.asarray(response_weights,dtype=float)
        if (response_observations.ndim!=2 or
                response_observations.shape[1:]!=(3,) or
                response_weights.shape!=(nout,len(response_observations),3)):
            raise ValueError("native response sampling shape mismatch")
    incident=(np.zeros(nout) if incident_response is None else
              np.asarray(incident_response,dtype=float).reshape(-1))
    dincident=(np.zeros((q,nout)) if dincident_response is None else
               np.asarray(dincident_response,dtype=float))
    if incident.shape!=(nout,) or dincident.shape!=(q,nout):
        raise ValueError("functional incident response derivatives mismatch")

    charge_gram.restore_geometry_mass_matrix()
    right_hand_sides=np.ascontiguousarray(np.vstack((b,C)),dtype=np.float64)
    solved=charge_gram.solve_configured_linear_material_auto_prec_many(
        float(inv_chi),right_hand_sides,tol=float(solve_tolerance),
        maxit=int(solve_max_iterations),cluster_coarse_size=0,
        cluster_deflation_size=0,recycle_size=0,
        mass_riesz=bool(mass_riesz))
    solutions=np.asarray(solved["m"],dtype=float)
    iterations=tuple(int(value) for value in solved["iters"])
    if solutions.shape!=(nout+1,n) or len(iterations)!=nout+1:
        raise RuntimeError("native functional shape solve returned invalid shapes")
    state=solutions[0];adjoints=solutions[1:]
    response=C@state+incident

    geometry=sample_production_gettrafo_displacements(
        fes,modes,charge_basis,family=family)
    cells=np.stack(geometry.cell,axis=1)
    if geometry.family=="wedge":
        faces=np.zeros((q,len(geometry.face),9,3))
        for host,values in enumerate(geometry.face):
            faces[:,host,:values.shape[1]]=values
    else: faces=np.stack(geometry.face,axis=1)
    bx=np.asarray(B@state).reshape(-1)
    Gbx=np.asarray(charge_gram.matvec_sym(bx)).reshape(-1)
    badjoints=[np.asarray(B@value).reshape(-1) for value in adjoints]
    cell_velocities=np.ascontiguousarray(cells)
    face_velocities=np.ascontiguousarray(faces)
    if native_response:
        if geometry.family!="tet":
            raise ValueError("native configured-field response derivatives require TET geometry")
        native=getattr(
            charge_gram,
            "configured_field_functional_rows_directional_derivative",None)
        if native is None:
            raise RuntimeError(
                "the native configured-field row directional derivative is required")
        dC=np.asarray(native(
            np.ascontiguousarray(response_observations,dtype=np.float64),
            np.ascontiguousarray(response_weights,dtype=np.float64),
            cell_velocities,face_velocities),dtype=float)
        if dC.shape!=(q,nout,n) or not np.all(np.isfinite(dC)):
            raise RuntimeError("native configured-field response derivative returned invalid data")
    left_matrix=np.ascontiguousarray(np.stack(badjoints),dtype=float)
    if hasattr(charge_gram,"directional_derivative_contractions_many"):
        gram_terms=np.asarray(
            charge_gram.directional_derivative_contractions_many(
                geometry.family,cell_velocities,face_velocities,left_matrix,
                np.ascontiguousarray(bx)),dtype=float)
    else:  # Compatibility with an older downloaded native binary.
        gram_terms=np.stack([
            np.asarray(charge_gram.directional_derivative_contractions(
                geometry.family,cell_velocities,face_velocities,
                np.ascontiguousarray(left),np.ascontiguousarray(bx)),dtype=float)
            for left in badjoints],axis=0)
    mass_terms=assemble_ngsolve_hdiv_mass_shape_contractions(
        fes,modes,adjoints,state)
    jacobian=np.empty((nout,q),dtype=float)
    for k in range(q):
        if geometry.family=="tet":
            rates=np.asarray(charge_gram.tet_charge_map_row_directional_rates(
                np.ascontiguousarray(cells[k]),np.ascontiguousarray(faces[k])))
            dbx=rates*bx
            dbadjoints=left_matrix*rates[None,:]
        else:
            dbx=np.zeros_like(bx)
            dbadjoints=np.zeros_like(left_matrix)
        Gdbx=np.asarray(charge_gram.matvec_sym(dbx)).reshape(-1)
        jacobian[:,k]=(adjoints@db[k]-float(inv_chi)*mass_terms[:,k]
            -dbadjoints@Gbx-left_matrix@Gdbx-gram_terms[:,k]
            +dC[k]@state+dincident[k])
    return VIMFunctionalShapeJacobian(state,response,jacobian,
        int(iterations[0]),tuple(iterations[1:]))


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


@dataclass(frozen=True)
class ShapeLinearization:
    """One topology-preserving free-boundary linearization.

    ``response`` is compared with ``response_target`` using the absolute
    half-width ``response_band``.  The Jacobian columns correspond to the
    caller's GetTrafo displacement modes.  An empty response vector denotes
    unconstrained objective-only shape optimization.
    """
    objective: float
    objective_gradient: np.ndarray
    response: np.ndarray
    response_jacobian: np.ndarray
    response_target: np.ndarray
    response_band: np.ndarray


@dataclass(frozen=True)
class ShapeLPUpdate:
    """Lexicographic trust-region LP update for a real material boundary."""
    parameters: np.ndarray
    delta: np.ndarray
    predicted_objective_change: float
    predicted_max_band_ratio: float
    restoration: bool
    status: str
    iterations: int


@dataclass(frozen=True)
class ElementInsertionResponse:
    """Exact response of one finite element block added to a linear system."""
    candidate_state: np.ndarray
    active_state_delta: np.ndarray
    response_delta: np.ndarray
    schur_complement: np.ndarray


@dataclass(frozen=True)
class ElementGenerationLPUpdate:
    """Whole-element growth batch selected from superposed insertion responses."""
    selected: np.ndarray
    weights: np.ndarray
    predicted_response: np.ndarray
    predicted_max_band_ratio: float
    predicted_objective_change: float
    added_volume: float
    status: str


@dataclass(frozen=True)
class HDivMMMElementGenerationLinearization:
    """Matrix-free finite insertion responses for boundary-growth candidates."""
    state: np.ndarray
    response: np.ndarray
    candidate_elements: np.ndarray
    candidate_dof_blocks: tuple[np.ndarray, ...]
    candidate_response_delta: np.ndarray
    candidate_states: tuple[np.ndarray, ...]
    candidate_dof_offsets: np.ndarray
    reduced_schur_complement: np.ndarray
    reduced_schur_rhs: np.ndarray
    reduced_response_matrix: np.ndarray
    available_candidate_count: int
    native_reduction_timings: dict
    state_iterations: int
    adjoint_iterations: tuple[int, ...]
    schur_iterations: tuple[int, ...]


@dataclass(frozen=True)
class HDivMMMBlockInsertionResponse:
    """Exact fixed-active-set response of a finite candidate-element bundle."""
    selected_elements: np.ndarray
    candidate_state: np.ndarray
    response_delta: np.ndarray
    schur_complement: np.ndarray


@dataclass(frozen=True)
class CollaborativeElementBatchUpdate:
    """Best exact block-Schur bundle found by a discrete candidate search."""
    selected_elements: np.ndarray
    predicted_response: np.ndarray
    predicted_max_band_ratio: float
    added_volume: float
    evaluated_bundles: int
    status: str


@dataclass(frozen=True)
class TSVDElementCandidateSelection:
    """Global whole-element proposal obtained from every candidate response."""
    selected_elements: np.ndarray
    selected_directions: np.ndarray
    representative_elements: np.ndarray
    representative_directions: np.ndarray
    predicted_response: np.ndarray
    predicted_max_band_ratio: float
    added_volume: float
    numerical_rank: int
    aca_rank: int
    singular_values: np.ndarray
    signed_coefficients: np.ndarray
    relative_truncation_error: float
    status: str


@dataclass(frozen=True)
class HDivMMMGenerationIteration:
    iteration: int
    candidate_count: int
    added_elements: np.ndarray
    predicted_max_band_ratio: float
    actual_max_band_ratio: float
    active_element_count: int
    active_volume: float
    batch_trials: int
    solve_iterations: int
    source_scale: float = 1.0
    superposed_max_band_ratio: float | None = None
    selection_model: str = "superposed-lp"
    collaborative_bundles_evaluated: int = 0
    screened_candidate_count: int = 0
    native_reduction_timings: dict | None = None
    batch_limit_before: int = 1
    batch_limit_after: int = 1
    model_agreement_ratio: float = 1.0
    tsvd_rank: int = 0
    tsvd_aca_rank: int = 0
    tsvd_relative_truncation_error: float = 0.0
    tsvd_selected_count: int = 0
    removed_elements: np.ndarray | None = None
    addition_candidate_count: int = 0
    removal_candidate_count: int = 0


@dataclass(frozen=True)
class HDivMMMGenerationResult:
    active_elements: np.ndarray
    state: np.ndarray
    response: np.ndarray
    history: tuple[HDivMMMGenerationIteration, ...]
    converged: bool
    source_scale: float = 1.0
    objective_response: np.ndarray | None = None
    stop_reason: str = "unknown"


@dataclass(frozen=True)
class GrowthTopologyReport:
    """Connectivity gate for a whole-element iron-growth state."""
    iron_components: tuple[np.ndarray, ...]
    inactive_components: tuple[np.ndarray, ...]
    enclosed_inactive_elements: np.ndarray
    iron_connected: bool
    inactive_reaches_exterior: bool
    valid: bool


def finite_element_insertion_response(*, solve_active, active_state,
        active_to_candidate, candidate_matrix, candidate_rhs,
        active_response_matrix, candidate_response_matrix,
        active_adjoint=None) -> ElementInsertionResponse:
    """Return the exact one-element Schur/Woodbury response.

    The enlarged symmetric system is ``[[Aaa,Aae],[Aae.T,Aee]]``.  No density
    perturbation or finite difference is involved: the candidate block is
    inserted at full material strength.  ``solve_active`` may be a sparse,
    H-matrix, or batched native solve.  Supplying the already computed active
    adjoints reuses them across every boundary candidate.
    """
    ma=np.asarray(active_state,dtype=float).reshape(-1)
    Aae=np.asarray(active_to_candidate,dtype=float)
    Aee=np.asarray(candidate_matrix,dtype=float)
    be=np.asarray(candidate_rhs,dtype=float).reshape(-1)
    Ca=np.atleast_2d(np.asarray(active_response_matrix,dtype=float))
    Ce=np.atleast_2d(np.asarray(candidate_response_matrix,dtype=float))
    if Aae.ndim!=2 or Aae.shape[0]!=ma.size:
        raise ValueError("active_to_candidate must have one row per active state DOF")
    ne=Aae.shape[1]
    if Aee.shape!=(ne,ne) or be.shape!=(ne,):
        raise ValueError("candidate matrix/RHS shape mismatch")
    if Ca.shape[1]!=ma.size or Ce.shape!=(Ca.shape[0],ne):
        raise ValueError("candidate response matrices are incompatible")
    W=np.asarray(solve_active(Aae),dtype=float)
    if W.shape!=Aae.shape:
        raise ValueError("solve_active must preserve a two-dimensional RHS shape")
    schur=Aee-Aae.T@W
    residual=be-Aae.T@ma
    try:
        me=np.linalg.solve(schur,residual)
    except np.linalg.LinAlgError as error:
        raise RuntimeError("candidate element Schur complement is singular") from error
    dma=-W@me
    if active_adjoint is None:
        adj=np.asarray(solve_active(Ca.T),dtype=float)
    else:
        adj=np.asarray(active_adjoint,dtype=float)
    if adj.shape!=(ma.size,Ca.shape[0]):
        raise ValueError("active_adjoint must have shape (n_active,n_response)")
    # Ca*(-Aaa^-1*Aae*me) + Ce*me, evaluated from shared adjoints.
    dy=(Ce-adj.T@Aae)@me
    return ElementInsertionResponse(me,dma,np.asarray(dy).reshape(-1),schur)


def _ngsolve_volume_facet_graph(mesh):
    """Return NGSolve-owned volume adjacency and exterior-touching cells."""
    import ngsolve as ng
    elements=tuple(mesh.Elements(ng.VOL));owners={}
    for index,element in enumerate(elements):
        for facet in element.facets:
            owners.setdefault(int(facet.nr),[]).append(index)
    adjacency=[set() for _ in elements];exterior=set()
    for cells in owners.values():
        if len(cells)==1:
            exterior.add(cells[0])
        elif len(cells)==2:
            left,right=cells;adjacency[left].add(right);adjacency[right].add(left)
        else:
            raise RuntimeError("non-manifold volume facet has more than two owners")
    return tuple(np.asarray(sorted(row),dtype=np.int64) for row in adjacency),\
        np.asarray(sorted(exterior),dtype=np.int64)


def ngsolve_growth_topology(mesh, active_elements) -> GrowthTopologyReport:
    """Check connected iron and absence of an enclosed inactive component.

    An inactive component is admissible when it touches a boundary facet of
    the design superset: those components all connect through the unmeshed
    exterior air.  A component that touches no exterior facet is a generated
    internal cavity and fails loudly.  Prescribed apertures excluded from the
    superset are boundary facets and therefore remain admissible.
    """
    adjacency,exterior=_ngsolve_volume_facet_graph(mesh)
    active=np.asarray(active_elements,dtype=bool).reshape(-1)
    if active.shape!=(len(adjacency),) or not np.any(active):
        raise ValueError("active_elements must select a non-empty volume subset")

    def components(mask):
        remaining=set(int(v) for v in np.flatnonzero(mask));result=[]
        while remaining:
            seed=remaining.pop();stack=[seed];component=[seed]
            while stack:
                current=stack.pop()
                for neighbour in adjacency[current]:
                    neighbour=int(neighbour)
                    if neighbour in remaining:
                        remaining.remove(neighbour);stack.append(neighbour)
                        component.append(neighbour)
            result.append(np.asarray(sorted(component),dtype=np.int64))
        return tuple(result)

    iron=components(active);inactive=components(~active)
    exterior_set=set(int(v) for v in exterior)
    enclosed=[component for component in inactive
              if not any(int(cell) in exterior_set for cell in component)]
    enclosed_cells=(np.concatenate(enclosed) if enclosed else
                    np.empty(0,dtype=np.int64))
    iron_connected=len(iron)==1
    air_ok=enclosed_cells.size==0
    return GrowthTopologyReport(iron,inactive,np.sort(enclosed_cells),
                                iron_connected,air_ok,
                                bool(iron_connected and air_ok))


def ngsolve_boundary_growth_candidates(mesh, active_elements, *,
        fixed_inactive_elements=None, predecessor_elements=None,
        include_predecessor_descendants=False):
    """Return inactive volume elements sharing a complete facet with iron.

    Restricting generation to this set preserves iron connectivity when the
    initial seed is connected, but facet adjacency alone does *not* prevent a
    growing shell from enclosing a void.  Use :func:`ngsolve_growth_topology`
    as the acceptance gate.  ``predecessor_elements[e]`` optionally identifies
    the inward cell of an outward growth column; candidate ``e`` is eligible
    only after that predecessor is active.  NGSolve owns facet topology for
    TET/HEX/WEDGE; no vertex-order table is reimplemented here.
    With ``include_predecessor_descendants=True``, all later cells in an
    anchored growth column are also returned.  A caller selecting more than
    one layer in a batch must impose ``outer <= inner`` binary constraints.
    """
    import ngsolve as ng
    active=np.asarray(active_elements,dtype=bool).reshape(-1)
    elements=tuple(mesh.Elements(ng.VOL))
    if active.shape!=(len(elements),):
        raise ValueError("active_elements must have one flag per volume element")
    if not np.any(active):
        raise ValueError("boundary growth requires a non-empty iron seed")
    fixed=(np.zeros_like(active) if fixed_inactive_elements is None else
           np.asarray(fixed_inactive_elements,dtype=bool).reshape(-1))
    if fixed.shape!=active.shape or np.any(active&fixed):
        raise ValueError("fixed inactive elements must match the mesh and cannot be active")
    predecessors=(None if predecessor_elements is None else
                  np.asarray(predecessor_elements,dtype=np.int64).reshape(-1))
    if predecessors is not None:
        if predecessors.shape!=active.shape or np.any(predecessors>=len(active)) or np.any(predecessors<-1):
            raise ValueError("predecessor_elements must contain mesh indices or -1")
    facet_active=set()
    facet_inactive={}
    for index,element in enumerate(elements):
        facets=tuple(int(f.nr) for f in element.facets)
        if active[index]:
            facet_active.update(facets)
        else:
            for facet in facets: facet_inactive.setdefault(facet,[]).append(index)
    candidates=set()
    for facet in facet_active:
        candidates.update(facet_inactive.get(facet,()))
    candidates={cell for cell in candidates if not fixed[cell]}
    if predecessors is not None:
        candidates={cell for cell in candidates
                    if predecessors[cell]<0 or active[predecessors[cell]]}
        if include_predecessor_descendants:
            anchored=set(candidates)
            changed=True
            while changed:
                changed=False
                for cell in np.flatnonzero(~active&~fixed):
                    predecessor=int(predecessors[cell])
                    if (cell not in anchored and predecessor>=0
                            and (active[predecessor] or predecessor in anchored)):
                        anchored.add(int(cell));changed=True
            candidates=anchored
    return np.asarray(sorted(candidates),dtype=np.int64)


def ngsolve_boundary_removal_candidates(mesh, active_elements, *,
        fixed_active_elements=None, predecessor_elements=None):
    """Return removable outer iron cells while preserving the no-hole topology.

    A cell must lie on the current iron boundary, must not belong to the fixed
    initial yoke/pole, and cannot be the predecessor of another active growth
    cell.  Each removal is independently checked by the same iron/exterior-air
    connectivity gate used for accepted batches.
    """
    active=np.asarray(active_elements,dtype=bool).reshape(-1)
    adjacency,exterior=_ngsolve_volume_facet_graph(mesh)
    if active.shape!=(len(adjacency),) or not np.any(active):
        raise ValueError("active_elements must select a non-empty mesh subset")
    fixed=(np.zeros_like(active) if fixed_active_elements is None else
           np.asarray(fixed_active_elements,dtype=bool).reshape(-1))
    if fixed.shape!=active.shape or np.any(fixed&~active):
        raise ValueError("fixed_active_elements must be an active mesh subset")
    predecessors=(None if predecessor_elements is None else
                  np.asarray(predecessor_elements,dtype=np.int64).reshape(-1))
    if predecessors is not None and predecessors.shape!=active.shape:
        raise ValueError("predecessor_elements must match active_elements")
    has_active_child=np.zeros_like(active)
    if predecessors is not None:
        for child in np.flatnonzero(active):
            predecessor=int(predecessors[child])
            if predecessor>=0: has_active_child[predecessor]=True
    exterior_set=set(map(int,exterior));result=[]
    for element in np.flatnonzero(active&~fixed):
        element=int(element)
        if (element not in exterior_set and
                all(active[int(neighbour)] for neighbour in adjacency[element])):
            continue
        if has_active_child[element]:
            continue
        trial=active.copy();trial[element]=False
        if np.any(trial) and ngsolve_growth_topology(mesh,trial).valid:
            result.append(element)
    return np.asarray(result,dtype=np.int64)


def ngsolve_discontinuous_element_dof_blocks(fes):
    """Return one disjoint HDiv DOF block per volume element.

    Element generation requires a broken space so constraining an inactive
    element represents exactly zero magnetization without clamping its active
    neighbour.  This fail-loud check prevents accidental use of conforming HDiv.
    """
    import ngsolve as ng
    blocks=[]; owner={}
    for index,element in enumerate(fes.mesh.Elements(ng.VOL)):
        block=np.asarray([int(dof) for dof in fes.GetDofNrs(element)
                          if int(dof)>=0],dtype=np.int32)
        if block.size==0 or np.unique(block).size!=block.size:
            raise ValueError("each HDiv element must own a non-empty unique DOF block")
        for dof in block:
            if int(dof) in owner:
                raise ValueError(
                    "HDiv-MMM element generation requires a discontinuous HDiv space; "
                    f"DOF {int(dof)} is shared by elements {owner[int(dof)]} and {index}")
            owner[int(dof)]=index
        blocks.append(block)
    if len(owner)!=int(fes.ndof):
        missing=sorted(set(range(int(fes.ndof)))-set(owner))
        raise ValueError(
            "HDiv-MMM element generation requires every DOF to belong to one volume element "
            f"(unowned count={len(missing)})")
    return tuple(blocks)


def linearize_hdiv_mmm_element_generation(*, charge_gram, fes, inv_chi,
        rhs, response_matrix, active_elements, candidate_elements=None,
        incident_response=None, solve_tolerance=1e-9,
        solve_max_iterations=5000, candidate_batch_size=64,
        mass_riesz=True,
        candidate_selector=None) -> HDivMMMElementGenerationLinearization:
    """Close boundary element generation on the configured H-matrix operator.

    State and response adjoints are solved once on the active iron.  For every
    selected candidate RT/BDM block, a fused native call constructs full-strength
    columns of ``A``, performs the constrained active solves, and returns only
    the reduced Schur/RHS/response arrays.  When ``candidate_selector`` is
    supplied, all available candidates are first ranked using the local
    full-strength insertion that neglects active relaxation.  The callback is
    called as ``selector(elements, approximate_delta, state, response)`` and
    returns the element indices that receive the exact fused reduction.  The
    full dense ``A`` and ``N`` matrices are never materialized.
    """
    blocks=ngsolve_discontinuous_element_dof_blocks(fes)
    active_el=np.asarray(active_elements,dtype=bool).reshape(-1)
    if active_el.shape!=(len(blocks),) or not np.any(active_el):
        raise ValueError("active_elements must select a non-empty subset of the mesh")
    candidates=(ngsolve_boundary_growth_candidates(fes.mesh,active_el)
                if candidate_elements is None else
                np.asarray(candidate_elements,dtype=np.int64).reshape(-1))
    if candidates.size==0:
        raise ValueError("no boundary-growth candidate elements are available")
    if np.unique(candidates).size!=candidates.size or np.any(candidates<0) or np.any(candidates>=len(blocks)):
        raise ValueError("candidate element indices must be unique and in range")
    if np.any(active_el[candidates]):
        raise ValueError("candidate elements must currently be inactive")
    n=int(fes.ndof); rhs_full=np.asarray(rhs,dtype=float).reshape(-1)
    C=np.atleast_2d(np.asarray(response_matrix,dtype=float))
    if rhs_full.shape!=(n,) or C.shape[1]!=n:
        raise ValueError("HDiv-MMM generation RHS/response matrix shape mismatch")
    if not np.isfinite(inv_chi) or float(inv_chi)<0:
        raise ValueError("inv_chi must be finite and nonnegative")
    batch_size=int(candidate_batch_size)
    if batch_size<1: raise ValueError("candidate_batch_size must be positive")
    active_dofs=np.concatenate([blocks[k] for k in np.flatnonzero(active_el)]).astype(np.int32)
    active_mask=np.zeros(n,dtype=bool); active_mask[active_dofs]=True
    inactive_dofs=np.flatnonzero(~active_mask).astype(np.int32)
    charge_gram.set_configured_constraints(inactive_dofs,preserve_existing=False)

    def solve_many(rows):
        data=np.ascontiguousarray(rows,dtype=np.float64)
        result=charge_gram.solve_configured_linear_material_auto_prec_many(
            float(inv_chi),data,tol=float(solve_tolerance),
            maxit=int(solve_max_iterations),cluster_coarse_size=0,
            cluster_deflation_size=0,recycle_size=0,
            mass_riesz=bool(mass_riesz))
        solutions=np.asarray(result["m"],dtype=float)
        iterations=tuple(int(x) for x in result["iters"])
        if solutions.shape!=data.shape or len(iterations)!=data.shape[0]:
            raise RuntimeError("native HDiv-MMM multi-RHS solve returned invalid shapes")
        return solutions,iterations

    base_rhs=np.vstack((rhs_full,C)).copy();base_rhs[:,inactive_dofs]=0.0
    solved,base_iterations=solve_many(base_rhs)
    state=solved[0];adjoints=solved[1:].T
    response=C@state
    if incident_response is not None:
        incident=np.asarray(incident_response,dtype=float).reshape(-1)
        if incident.shape!=response.shape:
            raise ValueError("incident_response shape mismatch")
        response=response+incident

    available_candidate_count=int(candidates.size)

    def apply_operator_many(dofs):
        dofs=np.asarray(dofs,dtype=np.int32).reshape(-1)
        basis=np.zeros((len(dofs),n),dtype=np.float64)
        basis[np.arange(len(dofs)),dofs]=1.0
        if hasattr(charge_gram,
                   "apply_configured_linear_material_operator_many"):
            applied=np.asarray(
                charge_gram.apply_configured_linear_material_operator_many(
                    float(inv_chi),np.ascontiguousarray(basis),
                    respect_constraints=False),dtype=float)
        else:  # Compatibility with an older downloaded native binary.
            applied=np.stack([
                np.asarray(charge_gram.apply_configured_linear_material_operator(
                    float(inv_chi),row,respect_constraints=False),dtype=float)
                for row in basis])
        if applied.shape!=basis.shape:
            raise RuntimeError("native HDiv-MMM operator batch has invalid shape")
        return applied

    if candidate_selector is not None:
        # The screening model retains the complete physical A_ee block and
        # source/response coupling of each full material element, but omits the
        # expensive A_aa^-1 relaxation.  Candidate groups are discarded as soon
        # as their small local response has been formed, bounding Python memory.
        approximate=[];begin=0
        while begin<len(candidates):
            end=begin;count=0
            while end<len(candidates):
                width=len(blocks[int(candidates[end])])
                if end>begin and count+width>batch_size:
                    break
                count+=width;end+=1
                if count>=batch_size:
                    break
            group_blocks=tuple(blocks[int(k)] for k in candidates[begin:end])
            group_dofs=np.concatenate(group_blocks).astype(np.int32)
            applied=apply_operator_many(group_dofs)
            local_offsets=np.r_[0,np.cumsum([len(block) for block in group_blocks])]
            for local_column,block in enumerate(group_blocks):
                local=np.arange(local_offsets[local_column],
                                local_offsets[local_column+1])
                columns=applied[local].T
                active_coupling=columns[active_dofs,:]
                local_matrix=0.5*(columns[block,:]+columns[block,:].T)
                local_rhs=(rhs_full[block]
                           -active_coupling.T@state[active_dofs])
                local_response=(C[:,block]
                    -adjoints[active_dofs].T@active_coupling)
                try:
                    local_state=np.linalg.solve(local_matrix,local_rhs)
                    approximate.append(local_response@local_state)
                except np.linalg.LinAlgError:
                    approximate.append(np.full(C.shape[0],np.nan))
            begin=end
        approximate_delta=np.stack(approximate,axis=1)
        selected=np.asarray(candidate_selector(
            candidates,approximate_delta,state,response),dtype=np.int64).reshape(-1)
        if selected.size==0 or np.unique(selected).size!=selected.size:
            raise ValueError("candidate_selector must return a non-empty unique subset")
        available=set(int(k) for k in candidates)
        if any(int(k) not in available for k in selected):
            raise ValueError("candidate_selector returned an unavailable element")
        candidates=selected

    candidate_blocks=tuple(blocks[int(k)] for k in candidates)
    candidate_dofs=np.concatenate(candidate_blocks).astype(np.int32)
    native_timings={"operator_s":0.0,"solve_s":0.0,"contraction_s":0.0}
    if hasattr(charge_gram,"reduce_configured_candidate_schur"):
        reduced=charge_gram.reduce_configured_candidate_schur(
            float(inv_chi),candidate_dofs,np.ascontiguousarray(rhs_full),
            np.ascontiguousarray(state),np.ascontiguousarray(C),
            np.ascontiguousarray(adjoints.T),tol=float(solve_tolerance),
            maxit=int(solve_max_iterations),solve_batch_size=batch_size,
            mass_riesz=bool(mass_riesz))
        reduced_schur=np.asarray(reduced["schur"],dtype=float)
        reduced_rhs=np.asarray(reduced["rhs"],dtype=float).reshape(-1)
        reduced_response=np.asarray(reduced["response"],dtype=float)
        schur_iterations=tuple(int(x) for x in reduced["iters"])
        native_timings={key:float(reduced[key]) for key in native_timings}
        nc=len(candidate_dofs)
        if (reduced_schur.shape!=(nc,nc) or reduced_rhs.shape!=(nc,)
                or reduced_response.shape!=(C.shape[0],nc)
                or len(schur_iterations)!=nc):
            raise RuntimeError("native candidate Schur reduction returned invalid shapes")
    else:  # Compatibility with an older downloaded native binary.
        operator_columns=[]
        for begin in range(0,len(candidate_dofs),batch_size):
            operator_columns.extend(apply_operator_many(
                candidate_dofs[begin:begin+batch_size]))
        W_columns=[];schur_iterations=[]
        for begin in range(0,len(operator_columns),batch_size):
            columns=operator_columns[begin:begin+batch_size]
            coupling=np.zeros((len(columns),n),dtype=np.float64)
            coupling[:,active_dofs]=np.stack(columns,axis=1)[active_dofs].T
            W,iters=solve_many(coupling)
            W_columns.extend(W);schur_iterations.extend(iters)
        operator_matrix=np.stack(operator_columns,axis=1)
        active_solutions=np.stack(W_columns,axis=1)[active_dofs,:]
        active_coupling=operator_matrix[active_dofs,:]
        reduced_schur=(operator_matrix[candidate_dofs,:]
                       -active_coupling.T@active_solutions)
        reduced_schur=0.5*(reduced_schur+reduced_schur.T)
        reduced_rhs=(rhs_full[candidate_dofs]
                     -active_coupling.T@state[active_dofs])
        reduced_response=(C[:,candidate_dofs]
                          -adjoints[active_dofs].T@active_coupling)

    offsets=np.r_[0,np.cumsum([len(block) for block in candidate_blocks])]
    response_deltas=[];candidate_states=[]
    for column in range(len(candidate_blocks)):
        local=np.arange(offsets[column],offsets[column+1])
        schur=reduced_schur[np.ix_(local,local)]
        try: me=np.linalg.solve(schur,reduced_rhs[local])
        except np.linalg.LinAlgError as error:
            raise RuntimeError("candidate HDiv element Schur complement is singular") from error
        candidate_states.append(me)
        response_deltas.append(reduced_response[:,local]@me)
    delta=np.stack(response_deltas,axis=1)
    return HDivMMMElementGenerationLinearization(state,response,candidates,
        candidate_blocks,delta,tuple(candidate_states),offsets.astype(np.int32),
        np.ascontiguousarray(reduced_schur),np.ascontiguousarray(reduced_rhs),
        np.ascontiguousarray(reduced_response),available_candidate_count,
        native_timings,int(base_iterations[0]),
        tuple(base_iterations[1:]),tuple(schur_iterations))


def hdiv_mmm_block_insertion_response(linearization,
        selected_elements) -> HDivMMMBlockInsertionResponse:
    """Evaluate a full-strength multi-element insertion from one reduced Schur matrix.

    ``linearize_hdiv_mmm_element_generation`` has already eliminated every
    active DOF.  This routine therefore solves only the candidate DOFs in the
    requested bundle, while retaining all candidate-candidate interactions.
    It is algebraically identical to enlarging and solving the active system,
    but does not repeat an H-matrix solve during combinatorial search.
    """
    selected=np.asarray(selected_elements,dtype=np.int64).reshape(-1)
    if selected.size==0 or np.unique(selected).size!=selected.size:
        raise ValueError("selected_elements must contain a non-empty unique bundle")
    candidates=np.asarray(linearization.candidate_elements,dtype=np.int64)
    lookup={int(element):column for column,element in enumerate(candidates)}
    if any(int(element) not in lookup for element in selected):
        raise ValueError("selected_elements must belong to the current candidate set")
    selected=np.sort(selected)
    offsets=np.asarray(linearization.candidate_dof_offsets,dtype=np.int64)
    local=np.concatenate([
        np.arange(offsets[lookup[int(element)]],
                  offsets[lookup[int(element)]+1],dtype=np.int64)
        for element in selected])
    schur=np.asarray(linearization.reduced_schur_complement,dtype=float)[
        np.ix_(local,local)]
    rhs=np.asarray(linearization.reduced_schur_rhs,dtype=float)[local]
    try:
        state=np.linalg.solve(schur,rhs)
    except np.linalg.LinAlgError as error:
        raise RuntimeError("candidate bundle Schur complement is singular") from error
    response=(np.asarray(linearization.reduced_response_matrix,dtype=float)[:,local]
              @state)
    return HDivMMMBlockInsertionResponse(
        selected,state,np.asarray(response).reshape(-1),schur)


def select_collaborative_element_batch(*, current_response, response_target,
        response_band, candidate_elements, candidate_volumes,
        evaluate_bundle_response, volume_budget, maximum_new_elements,
        active_elements=None, predecessor_elements=None,
        candidate_limit=32, beam_width=48,
        seed_bundles=(), ratio_tolerance=1e-12,
        improvement_capture=0.9
        ) -> CollaborativeElementBatchUpdate:
    """Find synergistic whole-element bundles with exact pair/block responses.

    Every admissible singleton and pair in the screened pool is evaluated,
    even when neither singleton improves the current solution.  Larger bundles
    are constructed by a bounded beam search.  ``evaluate_bundle_response``
    owns the physical block-Schur solve and any source recalibration or
    nonlinear response transform.  Predecessor closure enforces outward-only
    growth throughout the search.
    """
    current=np.asarray(current_response,dtype=float).reshape(-1)
    target=np.asarray(response_target,dtype=float).reshape(-1)
    band=np.asarray(response_band,dtype=float).reshape(-1)
    candidates=np.asarray(candidate_elements,dtype=np.int64).reshape(-1)
    volumes=np.asarray(candidate_volumes,dtype=float).reshape(-1)
    if current.shape!=target.shape or target.shape!=band.shape or target.size==0:
        raise ValueError("current, target, and band responses must have one shape")
    if np.any(band<=0.0) or candidates.shape!=volumes.shape or np.any(volumes<=0.0):
        raise ValueError("candidate volumes and response bands must be positive")
    if (candidates.size==0 or np.unique(candidates).size!=candidates.size
            or int(maximum_new_elements)<1 or int(candidate_limit)<2
            or int(beam_width)<1):
        raise ValueError("collaborative search dimensions are invalid")
    volume_budget=float(volume_budget)
    if not np.isfinite(volume_budget) or volume_budget<0.0:
        raise ValueError("volume_budget must be finite and nonnegative")
    improvement_capture=float(improvement_capture)
    if not 0.0<improvement_capture<=1.0:
        raise ValueError("improvement_capture must lie in (0,1]")
    lookup={int(element):column for column,element in enumerate(candidates)}
    active=(None if active_elements is None else
            np.asarray(active_elements,dtype=bool).reshape(-1))
    predecessors=(None if predecessor_elements is None else
                  np.asarray(predecessor_elements,dtype=np.int64).reshape(-1))
    if predecessors is not None:
        if active is None or predecessors.shape!=active.shape:
            raise ValueError("predecessor search requires matching active_elements")

    def ratio(values):
        return float(np.max(np.abs((np.asarray(values)-target)/band)))

    def valid(bundle):
        if not bundle or len(bundle)>int(maximum_new_elements): return False
        if any(element not in lookup for element in bundle): return False
        if sum(volumes[lookup[element]] for element in bundle)>volume_budget+1e-14:
            return False
        if predecessors is None: return True
        selected=set(bundle)
        for element in bundle:
            predecessor=int(predecessors[element])
            if predecessor>=0 and not active[predecessor] and predecessor not in selected:
                return False
        return True

    cache={}
    evaluated=0
    def evaluate(bundle):
        nonlocal evaluated
        key=tuple(sorted(int(element) for element in bundle))
        if key in cache: return cache[key]
        if not valid(key):
            cache[key]=None;return None
        evaluated_response=evaluate_bundle_response(
            np.asarray(key,dtype=np.int64))
        if evaluated_response is None:
            cache[key]=None;return None
        values=np.asarray(evaluated_response,dtype=float).reshape(-1)
        if values.shape!=target.shape or not np.all(np.isfinite(values)):
            raise ValueError("evaluate_bundle_response returned an invalid response")
        evaluated+=1
        item=(ratio(values),float(sum(volumes[lookup[e]] for e in key)),key,values)
        cache[key]=item
        return item

    singles=[]
    for element in candidates:
        item=evaluate((int(element),))
        if item is not None: singles.append(item)
    singles.sort(key=lambda item:(item[0],item[1],item[2]))

    # Exact singleton scores only screen very large fronts.  Add every missing
    # predecessor so no outward column is cut by screening.  Candidate fronts
    # up to candidate_limit are exhaustive at the pair level.
    if len(candidates)<=int(candidate_limit):
        screened=set(int(element) for element in candidates)
    else:
        screened=set(item[2][0] for item in singles[:int(candidate_limit)])
        if predecessors is not None:
            pending=list(screened)
            while pending:
                predecessor=int(predecessors[pending.pop()])
                if predecessor in lookup and predecessor not in screened:
                    screened.add(predecessor);pending.append(predecessor)
    screened=tuple(sorted(screened))

    items=list(singles)
    for seed in seed_bundles:
        item=evaluate(tuple(int(element) for element in seed))
        if item is not None: items.append(item)
    pair_items=[]
    if int(maximum_new_elements)>=2:
        for bundle in combinations(screened,2):
            item=evaluate(bundle)
            if item is not None:
                items.append(item);pair_items.append(item)
    pair_items.sort(key=lambda item:(item[0],item[1],item[2]))
    beam=pair_items[:int(beam_width)]
    for size in range(3,int(maximum_new_elements)+1):
        expanded={}
        for item in beam:
            base=set(item[2])
            for element in screened:
                if element in base: continue
                candidate=tuple(sorted(base|{element}))
                if len(candidate)!=size or candidate in expanded: continue
                result=evaluate(candidate)
                if result is not None: expanded[candidate]=result
        beam=sorted(expanded.values(),key=lambda item:(item[0],item[1],item[2]))[
            :int(beam_width)]
        items.extend(beam)
        if not beam: break

    current_ratio=ratio(current)
    improving=[item for item in items if item[0]<current_ratio-float(ratio_tolerance)]
    if not improving:
        return CollaborativeElementBatchUpdate(
            np.empty(0,dtype=np.int64),current,current_ratio,0.0,evaluated,
            "no improving exact block-Schur bundle")
    # A pure minimum-ratio criterion almost always consumes the entire allowed
    # batch for a tiny final gain.  Prefer the smallest-volume bundle that
    # captures most of the best predicted reduction; once any bundle reaches
    # the target band, choose the smallest target-reaching bundle instead.
    target_reaching=[item for item in improving if item[0]<=1.0+ratio_tolerance]
    if target_reaching:
        best=min(target_reaching,key=lambda item:(item[1],item[0],item[2]))
    else:
        best_reduction=max(current_ratio-item[0] for item in improving)
        efficient=[item for item in improving
                   if current_ratio-item[0]>=improvement_capture*best_reduction]
        best=min(efficient,key=lambda item:(item[1],item[0],item[2]))
    return CollaborativeElementBatchUpdate(
        np.asarray(best[2],dtype=np.int64),best[3],best[0],best[1],evaluated,
        "exact block-Schur collaborative bundle")


def select_tsvd_element_candidates(*, current_response, response_target,
        response_band, candidate_elements, candidate_response_delta,
        candidate_volumes, volume_budget, active_elements=None,
        predecessor_elements=None, relative_tolerance=1e-3,
        improvement_capture=0.9, ratio_tolerance=1e-12,
        candidate_volume_changes=None, candidate_material_active=None
        ) -> TSVDElementCandidateSelection:
    """Select a global binary proposal from the TSVD of *all* candidates.

    Rows are normalized by their engineering bands before the decomposition,
    so a field row with large SI units cannot hide an optics row with small SI
    units.  TSVD removes response directions below ``relative_tolerance``.  A
    first 0--1 minimax solve finds the best truncated response under the volume
    and predecessor constraints; a second solve chooses the minimum-volume
    set that captures ``improvement_capture`` of that best reduction.  No
    element-count constraint is imposed here: proposal cardinality is a result
    of the numerical response rank, target, volume budget, and binary LP.
    When ``candidate_material_active`` is supplied, the retained TSVD
    pseudoinverse estimates signed material coefficients: positive coefficients
    activate currently inactive candidates and negative coefficients deactivate
    currently active candidates.  Infeasible opposite signs are discarded.

    The returned proposal is still only a global screening decision.  Its
    candidate block is subsequently evaluated with the exact Schur complement
    and accepted only after a complete active-system re-solve.
    """
    current=np.asarray(current_response,dtype=float).reshape(-1)
    target=np.asarray(response_target,dtype=float).reshape(-1)
    band=np.asarray(response_band,dtype=float).reshape(-1)
    elements=np.asarray(candidate_elements,dtype=np.int64).reshape(-1)
    delta=np.asarray(candidate_response_delta,dtype=float)
    volumes=np.asarray(candidate_volumes,dtype=float).reshape(-1)
    volume_changes=(volumes.copy() if candidate_volume_changes is None else
                    np.asarray(candidate_volume_changes,dtype=float).reshape(-1))
    material_active=(None if candidate_material_active is None else
                     np.asarray(candidate_material_active,dtype=bool).reshape(-1))
    if (current.shape!=target.shape or target.shape!=band.shape or
            target.size==0 or np.any(band<=0.0)):
        raise ValueError("TSVD response, target, and band vectors must match")
    if (delta.shape!=(target.size,elements.size) or
            volumes.shape!=elements.shape or volume_changes.shape!=elements.shape or
            elements.size==0 or (material_active is not None and
                                  material_active.shape!=elements.shape) or
            np.unique(elements).size!=elements.size or np.any(volumes<=0.0)):
        raise ValueError("TSVD candidate arrays have incompatible shapes")
    relative_tolerance=float(relative_tolerance)
    improvement_capture=float(improvement_capture)
    if (not 0.0<relative_tolerance<1.0 or
            not 0.0<improvement_capture<=1.0 or
            not np.all(np.isfinite(delta))):
        raise ValueError("TSVD tolerances and candidate responses must be finite")

    normalized=np.ascontiguousarray(delta/band[:,None])
    # The repository's canonical ACA+ -> thin QR -> small TSVD kernel owns the
    # factorization.  It treats every candidate as a source column while
    # sampling only the entries required by ACA; no Python SVD implementation
    # or candidate-count heuristic is maintained here.
    from .stream_function import aca_tsvd
    factor=aca_tsvd(
        normalized.shape[0],normalized.shape[1],
        lambda row,column:normalized[row,column],
        modes=min(normalized.shape),kmax=min(normalized.shape),
        aca_eps=relative_tolerance,method="aca_qr_tsvd")
    U=np.asarray(factor.U,dtype=float)
    singular=np.asarray(factor.S,dtype=float)
    V=np.asarray(factor.V,dtype=float)
    residual_scale=max(1.0,float(np.linalg.norm((current-target)/band)))
    if (singular.size==0 or not np.isfinite(singular[0]) or
            singular[0]<=relative_tolerance*residual_scale):
        return TSVDElementCandidateSelection(
            np.empty(0,dtype=np.int64),np.empty(0,dtype=np.int8),
            np.empty(0,dtype=np.int64),np.empty(0,dtype=np.int8),current.copy(),
            float(np.max(np.abs((current-target)/band))),0.0,0,
            int(factor.k_aca),np.asarray(singular,dtype=float),
            np.zeros(elements.size),0.0,
            "all normalized candidate responses are zero")
    rank=max(1,int(np.count_nonzero(
        singular>=relative_tolerance*singular[0])))
    retained=(U[:,:rank]*singular[:rank])@V[:,:rank].T
    correction=(target-current)/band
    coefficients=V[:,:rank]@((U[:,:rank].T@correction)/singular[:rank])
    truncated_delta=band[:,None]*retained
    discarded=float(np.linalg.norm(normalized-retained))
    total=float(np.linalg.norm(normalized))
    relative_error=discarded/total if total>0.0 else 0.0
    # Column-pivoted QR of V^T is an interpolative skeleton of the complete
    # candidate set.  Keep those representatives in the exact Schur front even
    # when the binary TSVD proposal does not activate them; they preserve
    # distinct response modes and allow exact pair collaboration to overturn
    # the additive screening model.
    if material_active is None:
        feasible=np.ones(elements.size,dtype=bool)
        directions=np.ones(elements.size,dtype=np.int8)
    else:
        coefficient_floor=relative_tolerance*max(
            1.0e-300,float(np.max(np.abs(coefficients))))
        feasible=((~material_active)&(coefficients>coefficient_floor))|\
                 (material_active&(coefficients<-coefficient_floor))
        directions=np.where(material_active,-1,1).astype(np.int8)
    feasible_columns=np.flatnonzero(feasible)
    if feasible_columns.size==0:
        return TSVDElementCandidateSelection(
            np.empty(0,dtype=np.int64),np.empty(0,dtype=np.int8),
            np.empty(0,dtype=np.int64),np.empty(0,dtype=np.int8),current.copy(),
            float(np.max(np.abs((current-target)/band))),0.0,rank,
            int(factor.k_aca),np.asarray(singular,dtype=float),coefficients,
            relative_error,"TSVD magnetization signs admit no feasible boundary move")
    from scipy.linalg import qr
    feasible_v=V[feasible_columns,:rank]
    _,_,pivots=qr(feasible_v.T,mode="economic",pivoting=True)
    representative_columns=set(int(column) for column in pivots[:rank])
    for mode in range(rank):
        representative_columns.add(int(np.argmax(feasible_v[:,mode])))
        representative_columns.add(int(np.argmin(feasible_v[:,mode])))
    representative_columns=feasible_columns[np.asarray(
        sorted(representative_columns),dtype=np.int64)]
    representatives=elements[representative_columns]
    representative_directions=directions[representative_columns]

    move_elements=elements[feasible_columns]
    move_directions=directions[feasible_columns]
    move_delta=truncated_delta[:,feasible_columns]*move_directions[None,:]
    move_volumes=volumes[feasible_columns]
    move_volume_changes=volume_changes[feasible_columns]*move_directions

    predecessor_pairs=None
    if predecessor_elements is not None:
        predecessors=np.asarray(predecessor_elements,dtype=np.int64).reshape(-1)
        active=(None if active_elements is None else
                np.asarray(active_elements,dtype=bool).reshape(-1))
        if active is None or predecessors.shape!=active.shape:
            raise ValueError("TSVD predecessor selection requires active_elements")
        lookup={int(element):column for column,element in enumerate(move_elements)}
        predecessor_pairs=[(column,lookup[int(predecessors[element])])
            for column,element in enumerate(move_elements)
            if move_directions[column]>0 and
               int(predecessors[element]) in lookup]

    best=solve_element_generation_lp(
        current,target,band,move_delta,move_volumes,
        volume_budget=float(volume_budget),maximum_new_elements=None,
        whole_elements=True,predecessor_pairs=predecessor_pairs,
        candidate_volume_change=move_volume_changes)
    current_ratio=float(np.max(np.abs((current-target)/band)))
    if (not np.any(best.selected) or
            best.predicted_max_band_ratio>=current_ratio-float(ratio_tolerance)):
        return TSVDElementCandidateSelection(
            np.empty(0,dtype=np.int64),np.empty(0,dtype=np.int8),
            representatives,representative_directions,current.copy(),current_ratio,0.0,
            rank,int(factor.k_aca),np.asarray(singular,dtype=float),coefficients,
            relative_error,
            "TSVD global binary model found no improving insertion")
    capture_ratio=(current_ratio-improvement_capture*
                   (current_ratio-best.predicted_max_band_ratio))
    compact=solve_element_generation_lp(
        current,target,band,move_delta,move_volumes,
        volume_budget=float(volume_budget),maximum_new_elements=None,
        whole_elements=True,predecessor_pairs=predecessor_pairs,
        predicted_ratio_cap=capture_ratio,
        candidate_volume_change=move_volume_changes)
    selected_columns=feasible_columns[np.asarray(compact.selected,dtype=bool)]
    selected=elements[selected_columns]
    return TSVDElementCandidateSelection(
        np.asarray(selected,dtype=np.int64),directions[selected_columns],
        np.asarray(representatives,dtype=np.int64),representative_directions,
        compact.predicted_response,
        compact.predicted_max_band_ratio,compact.added_volume,rank,
        int(factor.k_aca),np.asarray(singular,dtype=float),coefficients,relative_error,
        "all-candidate band-normalized TSVD plus binary LP")


def select_tsvd_exact_block_batch(*, current_response, response_target,
        response_band, candidate_elements, candidate_volumes,
        proposal_elements, representative_elements, evaluate_bundle_response,
        volume_budget, active_elements=None, predecessor_elements=None,
        maximum_new_elements=None, improvement_capture=0.9,
        ratio_tolerance=1e-12, bundle_is_valid=None
        ) -> CollaborativeElementBatchUpdate:
    """Refine one TSVD proposal by conditional exact Schur selection.

    The former bounded beam repeatedly expanded candidate combinations up to a
    caller-supplied cardinality.  Here cardinality comes from the global TSVD
    proposal.  Exact singleton and signed/QR representative-pair seeds are
    evaluated, then the best seed is enlarged one element at a time by its
    conditional exact-Schur improvement.  A backward path is used only when a
    higher-order TSVD proposal improves although no singleton/pair does.  This
    preserves pair-only collaboration without an exponential subset search.
    """
    current=np.asarray(current_response,dtype=float).reshape(-1)
    target=np.asarray(response_target,dtype=float).reshape(-1)
    band=np.asarray(response_band,dtype=float).reshape(-1)
    candidates=np.asarray(candidate_elements,dtype=np.int64).reshape(-1)
    volumes=np.asarray(candidate_volumes,dtype=float).reshape(-1)
    proposal=np.asarray(proposal_elements,dtype=np.int64).reshape(-1)
    representatives=np.asarray(representative_elements,dtype=np.int64).reshape(-1)
    if (current.shape!=target.shape or target.shape!=band.shape or
            np.any(band<=0.0) or candidates.shape!=volumes.shape or
            np.any(volumes<=0.0)):
        raise ValueError("exact TSVD block selection arrays are invalid")
    lookup={int(element):column for column,element in enumerate(candidates)}
    front=tuple(sorted(set(map(int,proposal))|set(map(int,representatives))))
    if not front or any(element not in lookup for element in front):
        raise ValueError("TSVD exact front must be a non-empty candidate subset")
    maximum=(None if maximum_new_elements is None else int(maximum_new_elements))
    if maximum is not None and maximum<1:
        raise ValueError("maximum_new_elements must be positive or None")
    active=(None if active_elements is None else
            np.asarray(active_elements,dtype=bool).reshape(-1))
    predecessors=(None if predecessor_elements is None else
                  np.asarray(predecessor_elements,dtype=np.int64).reshape(-1))
    if predecessors is not None and (active is None or predecessors.shape!=active.shape):
        raise ValueError("exact TSVD predecessor selection requires active_elements")

    def ratio(values):
        return float(np.max(np.abs((np.asarray(values)-target)/band)))
    def structurally_valid(bundle):
        selected=set(bundle)
        return bool(selected) and all(element in lookup for element in selected)
    def eligible(item):
        if item is None: return False
        selected=set(item[2])
        if item[1]>volume_budget+1e-14:
            return False
        if maximum is not None and len(selected)>maximum:
            return False
        if predecessors is not None:
            for element in selected:
                predecessor=int(predecessors[element])
                if predecessor>=0 and not active[predecessor] and predecessor not in selected:
                    return False
        if (bundle_is_valid is not None and not bundle_is_valid(
                np.asarray(sorted(selected),dtype=np.int64))):
            return False
        return True
    cache={};evaluated=0
    def evaluate(bundle):
        nonlocal evaluated
        key=tuple(sorted(set(map(int,bundle))))
        if key in cache: return cache[key]
        if not structurally_valid(key):
            cache[key]=None;return None
        values=evaluate_bundle_response(np.asarray(key,dtype=np.int64))
        if values is None:
            cache[key]=None;return None
        values=np.asarray(values,dtype=float).reshape(-1)
        if values.shape!=target.shape or not np.all(np.isfinite(values)):
            raise ValueError("evaluate_bundle_response returned an invalid response")
        evaluated+=1
        item=(ratio(values),float(sum(volumes[lookup[e]] for e in key)),key,values)
        cache[key]=item
        return item

    items=[];full_items=[]
    for seed in (front,tuple(sorted(set(map(int,proposal))))):
        item=evaluate(seed)
        if eligible(item):
            items.append(item);full_items.append(item)
    representative_front=tuple(sorted(set(map(int,representatives))))
    for element in front:
        item=evaluate((element,))
        if eligible(item): items.append(item)
    for pair in combinations(representative_front,2):
        item=evaluate(pair)
        if eligible(item): items.append(item)

    current_ratio=ratio(current)
    small_improving=[item for item in items if len(item[2])<=2 and
        item[0]<current_ratio-float(ratio_tolerance)]
    if small_improving:
        path=min(small_improving,key=lambda item:(item[0],item[1],item[2]))
        while len(path[2])<len(front):
            children=[];base=set(path[2])
            for element in front:
                if element in base: continue
                child=evaluate(tuple(sorted(base|{element})))
                if eligible(child):
                    children.append(child);items.append(child)
            if not children: break
            best_child=min(children,key=lambda item:(item[0],item[1],item[2]))
            if best_child[0]>=path[0]-float(ratio_tolerance): break
            path=best_child
            if path[0]<=1.0+ratio_tolerance: break
    else:
        # Rare higher-order synergy fallback: if the complete TSVD proposal is
        # improving while no singleton/pair is, prune one element per step.
        improving_full=[item for item in full_items
                        if item[0]<current_ratio-float(ratio_tolerance)]
        if improving_full:
            path=min(improving_full,key=lambda item:(item[0],item[1],item[2]))[2]
            while len(path)>1:
                children=[]
                for removed in path:
                    child=evaluate(tuple(element for element in path
                                         if element!=removed))
                    if child is not None: children.append(child)
                if not children: break
                items.extend(item for item in children if eligible(item))
                path=min(children,key=lambda item:(item[0],item[1],item[2]))[2]

    improving=[item for item in items
               if item[0]<current_ratio-float(ratio_tolerance)]
    if not improving:
        return CollaborativeElementBatchUpdate(
            np.empty(0,dtype=np.int64),current,current_ratio,0.0,evaluated,
            "no improving exact TSVD/Schur bundle")
    target_reaching=[item for item in improving if item[0]<=1.0+ratio_tolerance]
    if target_reaching:
        best=min(target_reaching,key=lambda item:(item[1],item[0],item[2]))
    else:
        best_reduction=max(current_ratio-item[0] for item in improving)
        efficient=[item for item in improving if
            current_ratio-item[0]>=float(improvement_capture)*best_reduction]
        best=min(efficient,key=lambda item:(item[1],item[0],item[2]))
    return CollaborativeElementBatchUpdate(
        np.asarray(best[2],dtype=np.int64),best[3],best[0],best[1],evaluated,
        "all-candidate TSVD plus conditional exact block-Schur selection")


def solve_hdiv_mmm_active_elements(*, charge_gram, fes, inv_chi, rhs,
        response_matrix, active_elements, incident_response=None,
        solve_tolerance=1e-9, solve_max_iterations=5000, mass_riesz=True):
    """Solve one exact whole-element active iron set on the fixed superset mesh."""
    blocks=ngsolve_discontinuous_element_dof_blocks(fes)
    active=np.asarray(active_elements,dtype=bool).reshape(-1)
    if active.shape!=(len(blocks),) or not np.any(active):
        raise ValueError("active_elements must select a non-empty subset")
    n=int(fes.ndof); rhs_full=np.asarray(rhs,dtype=float).reshape(-1)
    C=np.atleast_2d(np.asarray(response_matrix,dtype=float))
    if rhs_full.shape!=(n,) or C.shape[1]!=n:
        raise ValueError("active HDiv-MMM RHS/response matrix shape mismatch")
    active_dofs=np.concatenate([blocks[k] for k in np.flatnonzero(active)]).astype(np.int32)
    mask=np.ones(n,dtype=bool);mask[active_dofs]=False
    charge_gram.set_configured_constraints(np.flatnonzero(mask).astype(np.int32),
                                           preserve_existing=False)
    active_rhs=rhs_full.copy();active_rhs[mask]=0.0
    result=charge_gram.solve_configured_linear_material_auto_prec_many(
        float(inv_chi),np.ascontiguousarray(active_rhs[None,:]),
        tol=float(solve_tolerance),maxit=int(solve_max_iterations),
        cluster_coarse_size=0,cluster_deflation_size=0,recycle_size=0,
        mass_riesz=bool(mass_riesz))
    state=np.asarray(result["m"],dtype=float)[0]
    response=C@state
    if incident_response is not None:
        incident=np.asarray(incident_response,dtype=float).reshape(-1)
        if incident.shape!=response.shape: raise ValueError("incident_response shape mismatch")
        response=response+incident
    return state,response,int(result["iters"][0])


def grow_hdiv_mmm_by_superposition(*, charge_gram, fes, inv_chi, rhs,
        response_matrix, active_elements, element_volumes,
        response_target, response_band, volume_max,
        incident_response=None, maximum_batch_elements=None,
        max_iterations=30, ratio_tolerance=1e-8,
        solve_tolerance=1e-9, solve_max_iterations=5000,
        candidate_batch_size=64, mass_riesz=True,
        fixed_inactive_elements=None,
        fixed_active_elements=None,
        predecessor_elements=None,
        source_calibration_rows=None,
        source_calibration_target=None,
        response_transform=None,
        include_predecessor_descendants=False,
        exact_candidate_limit=64,
        batch_improvement_capture=0.9,
        tsvd_relative_tolerance=1e-3,
        iteration_callback=None) -> HDivMMMGenerationResult:
    """Grow connected whole iron elements by Schur superposition and 0-1 LP.

    Each iteration cheaply evaluates every face-adjacent inactive addition and
    every topology-safe active boundary removal,
    factors the band-normalized global response by the canonical native
    ACA+--thin-QR--TSVD kernel, and solves a binary minimum-volume proposal on
    the retained response subspace.  Only those proposed DOFs enter the fused
    native block-Schur reduction and conditional singleton/pair/forward
    selection.  Signed singular-mode representative pairs are checked even
    when neither singleton improves, so cooperative insertions remain visible.
    A batch is committed only after an exact active-set HDiv-MMM re-solve.  The
    Positive recovered material coefficients add iron; negative coefficients
    remove only non-fixed outer iron.  Production imposes no cardinality cap:
    the TSVD rank, target, volume, and
    predecessor constraints determine how many elements are proposed.
    ``maximum_batch_elements`` is only an optional regression safety cap.
    Among block-Schur bundles the selector chooses the smallest volume that
    captures most of the best predicted reduction.  Thus TSVD-compressed
    superposition is a global proposal oracle, block Schur resolves local
    collaboration, and neither substitutes for the accepted physical solve.
    When ``source_calibration_rows`` is supplied, the RHS and incident response
    define a reference coil current.  The linear source amplitude is then
    recalibrated after every candidate insertion so the mean response on those
    rows equals the corresponding mean target.  This eliminates coil current
    analytically while the 0-1 LP remains responsible only for iron geometry.
    ``response_transform`` may map the raw linear field response to nonlinear
    design metrics.  It is evaluated on every exact one-element Schur response
    and exact accepted batch, so the LP sees finite metric changes rather than
    a finite-difference derivative.
    """
    active=np.asarray(active_elements,dtype=bool).reshape(-1).copy()
    volumes=np.asarray(element_volumes,dtype=float).reshape(-1)
    target=np.asarray(response_target,dtype=float).reshape(-1)
    band=np.asarray(response_band,dtype=float).reshape(-1)
    if active.shape!=volumes.shape or np.any(volumes<=0.0):
        raise ValueError("element_volumes must be positive and match active_elements")
    if target.shape!=band.shape or target.size==0 or np.any(band<=0.0):
        raise ValueError("generation target/band vectors are invalid")
    maximum_cap=(None if maximum_batch_elements is None else
                 int(maximum_batch_elements))
    if maximum_cap is not None and maximum_cap<1:
        raise ValueError("maximum_batch_elements must be positive or None")
    if maximum_cap is not None and int(exact_candidate_limit)<maximum_cap:
        raise ValueError("exact_candidate_limit must cover maximum_batch_elements")
    if (not 0.0<float(batch_improvement_capture)<=1.0 or
            not 0.0<float(tsvd_relative_tolerance)<1.0):
        raise ValueError("TSVD selection parameters are invalid")
    calibration=(None if source_calibration_rows is None else
                 np.asarray(source_calibration_rows,dtype=np.int64).reshape(-1))
    if calibration is not None:
        if calibration.size==0 or np.any(calibration<0):
            raise ValueError("source_calibration_rows must index a non-empty response subset")
        if source_calibration_target is None:
            if np.any(calibration>=target.size):
                raise ValueError("source calibration target is required for raw rows outside the design response")
            calibration_values=target[calibration]
        else:
            calibration_values=np.asarray(
                source_calibration_target,dtype=float).reshape(-1)
            if calibration_values.shape!=calibration.shape:
                raise ValueError("source_calibration_target must match source_calibration_rows")
        calibration_target=float(np.mean(calibration_values))
        if not np.isfinite(calibration_target) or calibration_target==0.0:
            raise ValueError("source calibration target mean must be finite and nonzero")

    def transform_response(raw_response):
        raw=np.asarray(raw_response,dtype=float).reshape(-1)
        transformed=(raw if response_transform is None else
                     np.asarray(response_transform(raw),dtype=float).reshape(-1))
        if transformed.shape!=target.shape or not np.all(np.isfinite(transformed)):
            raise ValueError("response_transform must return one finite design-response vector")
        return transformed

    def calibrate_source(base_state,base_response):
        state_value=np.asarray(base_state,dtype=float)
        response_value=np.asarray(base_response,dtype=float)
        if calibration is None:
            return state_value,response_value,1.0
        if np.any(calibration>=response_value.size):
            raise ValueError("source_calibration_rows index outside the raw response")
        denominator=float(np.mean(response_value[calibration]))
        if not np.isfinite(denominator) or denominator==0.0:
            raise RuntimeError("source calibration response mean is zero or invalid")
        scale=calibration_target/denominator
        if not np.isfinite(scale) or scale<=0.0:
            raise RuntimeError("source calibration requires a positive finite source scale")
        return state_value*scale,response_value*scale,float(scale)
    fixed=(np.zeros_like(active) if fixed_inactive_elements is None else
           np.asarray(fixed_inactive_elements,dtype=bool).reshape(-1))
    if fixed.shape!=active.shape or np.any(active&fixed):
        raise ValueError("fixed inactive elements must match the mesh and cannot be active")
    fixed_active=(np.zeros_like(active) if fixed_active_elements is None else
                  np.asarray(fixed_active_elements,dtype=bool).reshape(-1))
    if fixed_active.shape!=active.shape or np.any(fixed_active&~active):
        raise ValueError("fixed active elements must match and belong to the initial iron")
    predecessors=(None if predecessor_elements is None else
                  np.asarray(predecessor_elements,dtype=np.int64).reshape(-1))
    if predecessors is not None and (predecessors.shape!=active.shape
            or np.any(predecessors>=len(active)) or np.any(predecessors<-1)):
        raise ValueError("predecessor_elements must contain element indices or -1")
    if include_predecessor_descendants and predecessors is None:
        raise ValueError("predecessor descendants require predecessor_elements")
    topology=ngsolve_growth_topology(fes.mesh,active)
    if not topology.valid:
        raise ValueError("initial iron seed must be connected and contain no enclosed inactive cavity")
    state,response,solve_iterations=solve_hdiv_mmm_active_elements(
        charge_gram=charge_gram,fes=fes,inv_chi=inv_chi,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        incident_response=incident_response,solve_tolerance=solve_tolerance,
        solve_max_iterations=solve_max_iterations,mass_riesz=mass_riesz)
    state,response,source_scale=calibrate_source(state,response)
    element_blocks=ngsolve_discontinuous_element_dof_blocks(fes)
    ratio=lambda values:float(np.max(np.abs((np.asarray(values)-target)/band)))
    objective_response=transform_response(response)
    current_ratio=ratio(objective_response);history=[]
    converged=current_ratio<=1.0+ratio_tolerance
    stop_reason="target_met" if converged else "max_iterations"
    for iteration in range(int(max_iterations)):
        if converged:
            stop_reason="target_met";break
        remaining=float(volume_max)-float(volumes@active)
        candidates=ngsolve_boundary_growth_candidates(
            fes.mesh,active,fixed_inactive_elements=fixed,
            predecessor_elements=predecessors,
            include_predecessor_descendants=include_predecessor_descendants)
        removal_candidates=ngsolve_boundary_removal_candidates(
            fes.mesh,active,fixed_active_elements=fixed_active,
            predecessor_elements=predecessors)
        if candidates.size==0 and removal_candidates.size==0:
            stop_reason="no_growth_candidates";break
        if (removal_candidates.size==0 and
                remaining<min(volumes[candidates],default=np.inf)-1e-14):
            stop_reason="volume_budget_exhausted";break
        if candidates.size==0:
            material=[]
            for element in removal_candidates:
                block=element_blocks[int(element)]
                local=np.asarray(response_matrix,dtype=float)[:,block]@state[block]
                try:
                    _,removed,_=calibrate_source(state,response-local)
                    material.append(objective_response-transform_response(removed))
                except (RuntimeError,ValueError):
                    material.append(np.zeros_like(objective_response))
            removal_tsvd=select_tsvd_element_candidates(
                current_response=objective_response,response_target=target,
                response_band=band,candidate_elements=removal_candidates,
                candidate_response_delta=np.column_stack(material),
                candidate_volumes=volumes[removal_candidates],
                volume_budget=max(0.0,remaining),
                candidate_material_active=np.ones(len(removal_candidates),dtype=bool),
                relative_tolerance=float(tsvd_relative_tolerance),
                improvement_capture=float(batch_improvement_capture),
                ratio_tolerance=ratio_tolerance)
            selected_remove=np.asarray(removal_tsvd.selected_elements,dtype=np.int64)[
                np.asarray(removal_tsvd.selected_directions)<0]
            if selected_remove.size==0:
                stop_reason="no_improving_removal_candidate";break
            # A removal-only TSVD proposal is an additive magnetization model,
            # not an exact downdate of the active VIM system.  A collaborative
            # bundle can therefore be rejected even though one of its cells is
            # genuinely useful.  Evaluate the full proposal, nested smaller
            # proposals, and a bounded representative singleton front with
            # complete active-system re-solves.  This is the removal analogue
            # of conditional exact-Schur selection on the insertion path.
            material_matrix=np.column_stack(material)
            removal_lookup={int(element):column for column,element in
                            enumerate(removal_candidates)}
            coefficients=np.asarray(removal_tsvd.signed_coefficients,dtype=float)
            ordered=np.asarray(sorted(
                (int(element) for element in selected_remove),
                key=lambda element:abs(coefficients[removal_lookup[element]]),
                reverse=True),dtype=np.int64)
            attempts=[];seen=set()
            def add_attempt(elements):
                bundle=tuple(sorted(int(element) for element in
                                    np.asarray(elements,dtype=np.int64).reshape(-1)))
                if bundle and bundle not in seen:
                    seen.add(bundle);attempts.append(np.asarray(bundle,dtype=np.int64))
            add_attempt(ordered)
            size=len(ordered)
            while size>1:
                size=max(1,(size+1)//2);add_attempt(ordered[:size])
            single_ratio=np.asarray([
                ratio(objective_response-material_matrix[:,column])
                for column in range(len(removal_candidates))])
            ranked=sorted((int(element) for element in removal_candidates),
                          key=lambda element:single_ratio[
                              removal_lookup[element]])
            tsvd_representatives=[int(x) for x in
                                  removal_tsvd.representative_elements]
            front_limit=min(int(exact_candidate_limit),
                            max(len(tsvd_representatives),
                                2*int(removal_tsvd.numerical_rank)+2))
            # ``ordered`` can contain the entire TSVD bundle.  It is already
            # covered by the full/nested bundle trials above and must not
            # silently turn the singleton fallback back into one full solve
            # per candidate.  Rank the union of TSVD/QR representatives,
            # coefficient-leading cells, and predicted best singletons, then
            # enforce the advertised exact-front bound.
            representative_pool=set(tsvd_representatives)
            representative_pool.update(int(x) for x in ordered[:front_limit])
            representative_pool.update(ranked[:front_limit])
            representative=set(sorted(
                representative_pool,
                key=lambda value:single_ratio[removal_lookup[value]])[
                    :front_limit])
            for element in sorted(representative,
                                  key=lambda value:single_ratio[
                                      removal_lookup[value]]):
                add_attempt([element])

            best=None;batch_trials=0
            for bundle in attempts:
                trial_active=active.copy();trial_active[bundle]=False
                if not ngsolve_growth_topology(fes.mesh,trial_active).valid:
                    continue
                batch_trials+=1
                trial_state,trial_response,trial_iterations=\
                    solve_hdiv_mmm_active_elements(
                        charge_gram=charge_gram,fes=fes,inv_chi=inv_chi,rhs=rhs,
                        response_matrix=response_matrix,active_elements=trial_active,
                        incident_response=incident_response,
                        solve_tolerance=solve_tolerance,
                        solve_max_iterations=solve_max_iterations,
                        mass_riesz=mass_riesz)
                try:
                    trial_state,trial_response,trial_scale=calibrate_source(
                        trial_state,trial_response)
                    trial_objective=transform_response(trial_response)
                except (RuntimeError,ValueError):
                    continue
                trial_ratio=ratio(trial_objective)
                if (trial_ratio<current_ratio-ratio_tolerance and
                        (best is None or trial_ratio<best[0])):
                    best=(trial_ratio,bundle,trial_active,trial_state,
                          trial_response,trial_iterations,trial_scale,
                          trial_objective)
            if best is None:
                stop_reason="conditional_exact_rejected_removal_front";break
            (actual,selected_remove,trial_active,trial_state,trial_response,
             trial_iterations,trial_scale,trial_objective)=best
            selected_columns=np.asarray([
                removal_lookup[int(element)] for element in selected_remove],
                dtype=np.int64)
            predicted=ratio(objective_response-
                            np.sum(material_matrix[:,selected_columns],axis=1))
            agreement=((current_ratio-actual)/(current_ratio-predicted)
                       if predicted<current_ratio-ratio_tolerance else 0.0)
            active=trial_active;state=trial_state;response=trial_response
            solve_iterations=trial_iterations;source_scale=trial_scale
            objective_response=trial_objective;current_ratio=actual
            row=HDivMMMGenerationIteration(
                iteration,len(removal_candidates),np.empty(0,dtype=np.int64),
                predicted,current_ratio,int(np.count_nonzero(active)),
                float(volumes@active),batch_trials,solve_iterations,source_scale,
                predicted,"signed-magnetization-aca-qr-tsvd-conditional-exact",
                batch_trials,len(representative),{},len(ordered),
                len(selected_remove),float(agreement),
                int(removal_tsvd.numerical_rank),int(removal_tsvd.aca_rank),
                float(removal_tsvd.relative_truncation_error),
                int(len(removal_tsvd.selected_elements)),selected_remove,0,
                int(len(removal_candidates)))
            history.append(row)
            if iteration_callback is not None: iteration_callback(row)
            converged=current_ratio<=1.0+ratio_tolerance
            if converged: stop_reason="target_met"
            continue

        tsvd_proposal=None
        def select_exact_candidates(elements,approximate_delta,
                                    approximate_state,approximate_response):
            nonlocal tsvd_proposal
            elements=np.asarray(elements,dtype=np.int64).reshape(-1)
            _,base_calibrated,_=calibrate_source(
                approximate_state,approximate_response)
            base_objective=transform_response(base_calibrated)
            effective=[];valid=[]
            for column in range(elements.size):
                delta=np.asarray(approximate_delta[:,column],dtype=float)
                if not np.all(np.isfinite(delta)):
                    continue
                try:
                    _,inserted,_=calibrate_source(
                        approximate_state,approximate_response+delta)
                    objective=transform_response(inserted)
                except (RuntimeError,ValueError):
                    continue
                valid.append(column);effective.append(objective-base_objective)
            if not valid:
                raise RuntimeError("all candidate insertion responses are invalid")
            valid=np.asarray(valid,dtype=np.int64)
            valid_elements=elements[valid]
            material_effective=list(effective)
            material_elements=list(map(int,valid_elements))
            material_is_active=[False]*len(material_elements)
            for element in removal_candidates:
                block=element_blocks[int(element)]
                local_material=(np.asarray(response_matrix,dtype=float)[:,block]
                                @np.asarray(approximate_state)[block])
                try:
                    _,removed,_=calibrate_source(
                        approximate_state,approximate_response-local_material)
                    removed_objective=transform_response(removed)
                except (RuntimeError,ValueError):
                    continue
                # TSVD columns represent positive material.  Removing the
                # active cell is the negative of this column; the recovered
                # coefficient sign chooses which operation is feasible.
                material_effective.append(base_objective-removed_objective)
                material_elements.append(int(element))
                material_is_active.append(True)
            material_elements=np.asarray(material_elements,dtype=np.int64)
            material_matrix=np.column_stack(material_effective)
            material_is_active=np.asarray(material_is_active,dtype=bool)
            tsvd_proposal=select_tsvd_element_candidates(
                current_response=base_objective,response_target=target,
                response_band=band,candidate_elements=material_elements,
                candidate_response_delta=material_matrix,
                candidate_volumes=volumes[material_elements],
                volume_budget=max(0.0,remaining),active_elements=active,
                predecessor_elements=predecessors,
                candidate_material_active=material_is_active,
                relative_tolerance=float(tsvd_relative_tolerance),
                improvement_capture=float(batch_improvement_capture),
                ratio_tolerance=ratio_tolerance)
            selected=np.asarray(tsvd_proposal.selected_elements,dtype=np.int64)[
                np.asarray(tsvd_proposal.selected_directions)>0]
            representatives=np.asarray(
                tsvd_proposal.representative_elements,dtype=np.int64)[
                np.asarray(tsvd_proposal.representative_directions)>0]
            front=np.union1d(selected,representatives).astype(np.int64)
            if valid_elements.size<=int(exact_candidate_limit):
                return valid_elements
            if front.size:
                return front
            # A no-improvement TSVD proposal still receives one exact Schur
            # probe.  When the zero-rank set is small, keep it whole so a pair
            # that improves only jointly is not discarded.  The limit controls
            # exact look-ahead work, not accepted batch cardinality.
            effective_matrix=np.column_stack(effective)
            score=np.max(np.abs((base_objective[:,None]+effective_matrix-
                                 target[:,None])/band[:,None]),axis=0)
            return valid_elements[np.asarray([int(np.argmin(score))])]

        lin=linearize_hdiv_mmm_element_generation(charge_gram=charge_gram,
            fes=fes,inv_chi=inv_chi,rhs=rhs,response_matrix=response_matrix,
            active_elements=active,candidate_elements=candidates,
            incident_response=incident_response,
            solve_tolerance=solve_tolerance,
            solve_max_iterations=solve_max_iterations,
            candidate_batch_size=candidate_batch_size,mass_riesz=mass_riesz,
            candidate_selector=select_exact_candidates)
        _,lp_raw_response,linear_scale=calibrate_source(lin.state,lin.response)
        lp_response=transform_response(lp_raw_response)
        state=lin.state*linear_scale
        response=lp_raw_response
        source_scale=linear_scale
        objective_response=lp_response
        current_ratio=ratio(objective_response)
        if tsvd_proposal is None:
            raise RuntimeError("all-candidate TSVD selector did not run")
        move_elements=np.asarray(tsvd_proposal.selected_elements,dtype=np.int64)
        move_directions=np.asarray(tsvd_proposal.selected_directions,dtype=np.int8)
        proposed_additions=move_elements[move_directions>0]
        proposed_removals=move_elements[move_directions<0]
        representative_elements=np.asarray(
            tsvd_proposal.representative_elements,dtype=np.int64)
        representative_directions=np.asarray(
            tsvd_proposal.representative_directions,dtype=np.int8)
        representative_additions=representative_elements[
            representative_directions>0]
        batch_limit_before=int(len(np.union1d(
            move_elements,representative_elements)))
        removed_elements=np.empty(0,dtype=np.int64)
        mixed_accepted=False

        # A negative recovered material coefficient means deletion of an
        # already active boundary cell.  Mixed add/remove proposals are tested
        # by one complete re-solve; no finite difference or gray material is
        # introduced.  If the signed low-rank model is inaccurate, fall back
        # to the exact-Schur addition front for this iteration.
        if proposed_removals.size:
            trial_active=active.copy()
            trial_active[proposed_additions]=True
            trial_active[proposed_removals]=False
            if ngsolve_growth_topology(fes.mesh,trial_active).valid:
                trial_state,trial_response,trial_iterations=\
                    solve_hdiv_mmm_active_elements(
                        charge_gram=charge_gram,fes=fes,inv_chi=inv_chi,
                        rhs=rhs,response_matrix=response_matrix,
                        active_elements=trial_active,
                        incident_response=incident_response,
                        solve_tolerance=solve_tolerance,
                        solve_max_iterations=solve_max_iterations,
                        mass_riesz=mass_riesz)
                try:
                    trial_state,trial_response,trial_scale=calibrate_source(
                        trial_state,trial_response)
                    trial_objective=transform_response(trial_response)
                    actual=ratio(trial_objective)
                except (RuntimeError,ValueError):
                    actual=np.inf
                if actual<current_ratio-ratio_tolerance:
                    selected=proposed_additions
                    removed_elements=proposed_removals
                    predicted_ratio=tsvd_proposal.predicted_max_band_ratio
                    exact_evaluated=1
                    selection_model=\
                        "signed-magnetization-aca-qr-tsvd-full-resolve"
                    mixed_accepted=True

        def evaluate_bundle(elements):
            block=hdiv_mmm_block_insertion_response(lin,elements)
            inserted=lin.response+block.response_delta
            try:
                _,calibrated,_=calibrate_source(lin.state,inserted)
                return transform_response(calibrated)
            except (RuntimeError,ValueError):
                return None

        if not mixed_accepted:
            exact=select_tsvd_exact_block_batch(
                current_response=lp_response,response_target=target,
                response_band=band,candidate_elements=lin.candidate_elements,
                candidate_volumes=volumes[lin.candidate_elements],
                proposal_elements=(proposed_additions if proposed_additions.size
                                   else representative_additions),
                representative_elements=(
                    lin.candidate_elements
                    if lin.available_candidate_count<=int(exact_candidate_limit)
                    else representative_additions),
                evaluate_bundle_response=evaluate_bundle,
                volume_budget=max(0.0,remaining),active_elements=active,
                predecessor_elements=predecessors,
                maximum_new_elements=maximum_cap,
                improvement_capture=float(batch_improvement_capture),
                ratio_tolerance=ratio_tolerance,
                bundle_is_valid=lambda elements:ngsolve_growth_topology(
                    fes.mesh,np.isin(np.arange(len(active)),elements)|active).valid)
            selected=exact.selected_elements
            predicted_ratio=exact.predicted_max_band_ratio
            exact_evaluated=exact.evaluated_bundles
            selection_model="all-candidate-aca-qr-tsvd-exact-conditional"
            if selected.size==0 or predicted_ratio>=current_ratio-ratio_tolerance:
                stop_reason="no_improving_exact_bundle";break
            # Block Schur is exact for the current active set, but acceptance
            # remains tied to a fresh full solve and the topology gate.
            trial_active=active.copy();trial_active[selected]=True
            trial_topology=ngsolve_growth_topology(fes.mesh,trial_active)
            if not trial_topology.valid:
                stop_reason="topology_gate_rejected_exact_bundle";break
            trial_state,trial_response,trial_iterations=\
                solve_hdiv_mmm_active_elements(
                    charge_gram=charge_gram,fes=fes,inv_chi=inv_chi,rhs=rhs,
                    response_matrix=response_matrix,active_elements=trial_active,
                    incident_response=incident_response,
                    solve_tolerance=solve_tolerance,
                    solve_max_iterations=solve_max_iterations,
                    mass_riesz=mass_riesz)
            try:
                trial_state,trial_response,trial_scale=calibrate_source(
                    trial_state,trial_response)
                trial_objective=transform_response(trial_response)
            except (RuntimeError,ValueError):
                stop_reason="source_calibration_rejected_exact_bundle";break
            actual=ratio(trial_objective)
            if actual>=current_ratio-ratio_tolerance:
                stop_reason="full_solve_rejected_exact_bundle";break
        predicted_improvement=current_ratio-predicted_ratio
        actual_improvement=current_ratio-actual
        agreement=(actual_improvement/predicted_improvement
                   if predicted_improvement>ratio_tolerance else 0.0)
        active=trial_active;state=trial_state;response=trial_response
        solve_iterations=trial_iterations;source_scale=trial_scale
        objective_response=trial_objective;current_ratio=actual
        batch_limit_after=len(selected)+len(removed_elements)
        history.append(HDivMMMGenerationIteration(iteration,
            lin.available_candidate_count+len(removal_candidates),
            np.asarray(selected,dtype=np.int64),predicted_ratio,
            current_ratio,int(np.count_nonzero(active)),float(volumes@active),
            1,solve_iterations,source_scale,
            tsvd_proposal.predicted_max_band_ratio,selection_model,
            exact_evaluated,
            len(lin.candidate_elements)+len(removal_candidates),
            dict(lin.native_reduction_timings),
            int(batch_limit_before),int(batch_limit_after),float(agreement),
            int(tsvd_proposal.numerical_rank),
            int(tsvd_proposal.aca_rank),
            float(tsvd_proposal.relative_truncation_error),
            int(len(tsvd_proposal.selected_elements)),
            np.asarray(removed_elements,dtype=np.int64),
            int(lin.available_candidate_count),int(len(removal_candidates))))
        if iteration_callback is not None:
            iteration_callback(history[-1])
        converged=current_ratio<=1.0+ratio_tolerance
        if converged: stop_reason="target_met"
    return HDivMMMGenerationResult(
        active,state,response,tuple(history),converged,source_scale,
        objective_response,stop_reason)


def solve_element_generation_lp(current_response, response_target,
        response_band, candidate_response_delta, candidate_volumes, *,
        volume_budget, maximum_new_elements=None,
        candidate_objective_change=None, whole_elements=True,
        relative_mip_gap=0.0, predecessor_pairs=None,
        predicted_ratio_cap=None,
        candidate_volume_change=None) -> ElementGenerationLPUpdate:
    """Select a small full-strength element-growth batch by a 0-1 LP.

    Single-element Schur responses are superposed only for the selection model;
    the committed batch must subsequently be solved as one exact HDiv-MMM
    problem.  With ``whole_elements=True`` (the default), HiGHS receives binary
    element variables, so no gray material can enter the physical model.
    """
    from scipy.optimize import Bounds, LinearConstraint, milp

    y=np.asarray(current_response,dtype=float).reshape(-1)
    target=np.asarray(response_target,dtype=float).reshape(-1)
    band=np.asarray(response_band,dtype=float).reshape(-1)
    D=np.asarray(candidate_response_delta,dtype=float)
    volumes=np.asarray(candidate_volumes,dtype=float).reshape(-1)
    volume_change=(volumes.copy() if candidate_volume_change is None else
                   np.asarray(candidate_volume_change,dtype=float).reshape(-1))
    if y.size==0 or target.shape!=y.shape or band.shape!=y.shape or np.any(band<=0):
        raise ValueError("element-generation response/target/band vectors are invalid")
    if D.ndim!=2 or D.shape[0]!=y.size or D.shape[1]!=volumes.size or volumes.size==0:
        raise ValueError("candidate_response_delta must have shape (n_response,n_candidate)")
    if (np.any(volumes<=0) or volume_change.shape!=volumes.shape or
            not np.all(np.isfinite(volume_change)) or
            not np.isfinite(volume_budget) or volume_budget<0):
        raise ValueError("candidate volumes must be positive and volume_budget nonnegative")
    if not np.all(np.isfinite(np.r_[y,target,band,D.ravel(),volumes])):
        raise ValueError("element-generation LP inputs must be finite")
    nc=volumes.size
    maximum=nc if maximum_new_elements is None else int(maximum_new_elements)
    if maximum<0: raise ValueError("maximum_new_elements must be nonnegative")
    if predicted_ratio_cap is not None:
        predicted_ratio_cap=float(predicted_ratio_cap)
        if not np.isfinite(predicted_ratio_cap) or predicted_ratio_cap<0.0:
            raise ValueError("predicted_ratio_cap must be finite and nonnegative")
    change=(volumes.copy() if candidate_objective_change is None else
            np.asarray(candidate_objective_change,dtype=float).reshape(-1))
    if change.shape!=(nc,) or not np.all(np.isfinite(change)):
        raise ValueError("candidate_objective_change must match candidate count")
    normalized=(y-target)/band; normalized_delta=D/band[:,None]
    rows=[np.c_[normalized_delta,-np.ones(y.size)],
          np.c_[-normalized_delta,-np.ones(y.size)],
          np.r_[volume_change,0.0][None,:],
          np.r_[np.ones(nc),0.0][None,:]]
    upper_parts=[-normalized,normalized,np.array([float(volume_budget)]),
                 np.array([float(maximum)])]
    if predecessor_pairs is not None:
        pairs=np.asarray(predecessor_pairs,dtype=np.int64)
        if pairs.size:
            pairs=pairs.reshape(-1,2)
            if np.any(pairs<0) or np.any(pairs>=nc) or np.any(pairs[:,0]==pairs[:,1]):
                raise ValueError("predecessor_pairs must contain distinct candidate indices")
            predecessor_rows=np.zeros((len(pairs),nc+1))
            predecessor_rows[np.arange(len(pairs)),pairs[:,0]]=1.0
            predecessor_rows[np.arange(len(pairs)),pairs[:,1]]=-1.0
            rows.append(predecessor_rows);upper_parts.append(np.zeros(len(pairs)))
    A=np.vstack(rows);upper=np.concatenate(upper_parts)
    constraint=LinearConstraint(A,np.full(A.shape[0],-np.inf),upper)
    bounds=Bounds(np.zeros(nc+1),np.r_[np.ones(nc),np.inf])
    integrality=np.r_[np.full(nc,1 if whole_elements else 0,dtype=int),0]
    options={"mip_rel_gap":float(relative_mip_gap)} if whole_elements else None
    phase1=milp(np.r_[np.zeros(nc),1.0],integrality=integrality,
        bounds=bounds,constraints=constraint,options=options)
    if not phase1.success:
        raise RuntimeError(f"element-generation restoration LP failed: {phase1.message}")
    best=max(0.0,float(phase1.x[-1]))
    # A normalized tiny volume tie-break keeps the batch compact when the
    # supplied objective has equal optima.
    scale=max(1.0,float(np.max(np.abs(change))))
    objective=np.r_[change+1e-10*scale*volumes/float(np.max(volumes)),0.0]
    phase2=None
    for relative_slack in (1e-9,1e-7,1e-5):
        if predicted_ratio_cap is None:
            cap=(max(1.0,best) if best<=1.0 else
                 best+relative_slack*max(1.0,best))
        else:
            cap=max(best,predicted_ratio_cap)+relative_slack*max(
                1.0,best,predicted_ratio_cap)
        phase2_bounds=Bounds(np.zeros(nc+1),np.r_[np.ones(nc),cap])
        phase2=milp(objective,integrality=integrality,bounds=phase2_bounds,
            constraints=constraint,options=options)
        if phase2.success:
            break
    if phase2 is not None and phase2.success:
        weights=np.asarray(phase2.x[:nc],dtype=float)
        status=str(phase2.message)
    else:
        # The minimax MILP is already a valid whole-element solution.  Some
        # HiGHS builds reject the second lexicographic phase when its t-bound
        # is within feasibility tolerance of the phase-one optimum; retaining
        # phase one is safer than turning a volume tie-break into a hard error.
        weights=np.asarray(phase1.x[:nc],dtype=float)
        status=(str(phase1.message)+"; objective tie-break unavailable: "
                +str(phase2.message if phase2 is not None else "not attempted"))
    selected=(weights>0.5) if whole_elements else (weights>0.0)
    predicted=y+D@weights
    ratio=float(np.max(np.abs((predicted-target)/band)))
    return ElementGenerationLPUpdate(selected,weights,predicted,ratio,
        float(change@weights),float(volume_change@weights),status)


def _shape_linearization_arrays(parameters, linearization):
    q=np.asarray(parameters,dtype=float).reshape(-1)
    g=np.asarray(linearization.objective_gradient,dtype=float).reshape(-1)
    y=np.asarray(linearization.response,dtype=float).reshape(-1)
    J=np.asarray(linearization.response_jacobian,dtype=float)
    target=np.asarray(linearization.response_target,dtype=float).reshape(-1)
    band=np.asarray(linearization.response_band,dtype=float).reshape(-1)
    if q.size==0 or g.shape!=q.shape:
        raise ValueError("shape parameters and objective gradient must be non-empty matching vectors")
    if y.size==0:
        if J.size and J.shape!=(0,q.size):
            raise ValueError("empty shape response requires a (0,n) Jacobian")
        J=np.zeros((0,q.size),dtype=float)
    elif J.shape!=(y.size,q.size):
        raise ValueError("shape response Jacobian must have shape (n_response,n_parameter)")
    if target.shape!=y.shape or band.shape!=y.shape:
        raise ValueError("shape response, target, and band vectors must match")
    if np.any(band<=0.0) or not np.all(np.isfinite(np.r_[q,g,y,J.ravel(),target,band])):
        raise ValueError("shape LP inputs must be finite and response bands positive")
    return q,g,y,J,target,band


def solve_shape_lp(parameters, linearization: ShapeLinearization, *, move_limit,
                   parameter_bounds=None, laplacian=None,
                   curvature_limit=None, A_ub=None, b_ub=None) -> ShapeLPUpdate:
    """Solve a topology-preserving GetTrafo shape step.

    The first LP minimizes the worst normalized response-band violation.  A
    second LP minimizes the physical objective without losing that optimum (or
    leaving the feasible response bands).  This is the shape analogue of a
    topology-optimization SLP, but every candidate remains an ordinary iron
    domain: there is no density, ersatz permeability, or gray material.
    """
    from scipy.optimize import linprog

    q,g,y,J,target,band=_shape_linearization_arrays(parameters,linearization)
    n=q.size
    move=np.broadcast_to(np.asarray(move_limit,dtype=float),q.shape).copy()
    if np.any(move<=0.0) or not np.all(np.isfinite(move)):
        raise ValueError("shape move_limit must contain positive finite values")
    if parameter_bounds is None:
        lower=np.full(n,-np.inf); upper=np.full(n,np.inf)
    else:
        if len(parameter_bounds)!=2:
            raise ValueError("parameter_bounds must be (lower,upper)")
        lower=np.broadcast_to(np.asarray(parameter_bounds[0],dtype=float),q.shape).copy()
        upper=np.broadcast_to(np.asarray(parameter_bounds[1],dtype=float),q.shape).copy()
        if np.any(lower>upper) or np.any(q<lower) or np.any(q>upper):
            raise ValueError("invalid shape parameter bounds or current parameters outside them")
    lower=np.maximum(lower,q-move); upper=np.minimum(upper,q+move)
    if np.any(lower>upper):
        raise ValueError("shape trust region does not intersect the parameter bounds")

    static_rows=[]; static_rhs=[]
    if laplacian is not None:
        L=np.atleast_2d(np.asarray(laplacian,dtype=float))
        if L.shape[1]!=n or curvature_limit is None:
            raise ValueError("shape laplacian requires n columns and curvature_limit")
        limit=np.broadcast_to(np.asarray(curvature_limit,dtype=float),(L.shape[0],))
        if np.any(limit<0.0): raise ValueError("curvature_limit must be nonnegative")
        static_rows.extend([L,-L]); static_rhs.extend([limit,limit])
    if A_ub is not None:
        extra=np.atleast_2d(np.asarray(A_ub,dtype=float))
        rhs=np.asarray(b_ub,dtype=float).reshape(-1)
        if extra.shape!=(rhs.size,n):
            raise ValueError("shape A_ub/b_ub mismatch")
        static_rows.append(extra); static_rhs.append(rhs)
    elif b_ub is not None:
        raise ValueError("shape b_ub requires A_ub")

    def _static_with_t(rows):
        return [np.c_[row,np.zeros((row.shape[0],1))] for row in rows]

    current_ratio=(float(np.max(np.abs((y-target)/band))) if y.size else 0.0)
    total_iterations=0
    if y.size:
        normalized=(y-target)/band
        normalized_jac=J/band[:,None]
        # normalized + Jn*(x-q) <= t and its negative counterpart.
        response_rows=[np.c_[normalized_jac,-np.ones(y.size)],
                       np.c_[-normalized_jac,-np.ones(y.size)]]
        response_rhs=[-normalized+normalized_jac@q,
                      normalized-normalized_jac@q]
        rows=response_rows+_static_with_t(static_rows)
        rhs=response_rhs+static_rhs
        phase1=linprog(np.r_[np.zeros(n),1.0],A_ub=np.vstack(rows),
            b_ub=np.concatenate(rhs),bounds=list(zip(lower,upper))+[(0.0,None)],
            method="highs")
        if not phase1.success:
            raise RuntimeError(f"shape restoration LP failed: {phase1.message}")
        total_iterations+=int(phase1.nit)
        best_ratio=max(0.0,float(phase1.x[-1]))
        ratio_cap=max(1.0,best_ratio) if best_ratio<=1.0 else best_ratio+1e-9*max(1.0,best_ratio)
        phase2=linprog(np.r_[g,0.0],A_ub=np.vstack(rows),b_ub=np.concatenate(rhs),
            bounds=list(zip(lower,upper))+[(0.0,ratio_cap)],method="highs")
        if not phase2.success:
            raise RuntimeError(f"shape objective LP failed: {phase2.message}")
        total_iterations+=int(phase2.nit); x=phase2.x[:n]
        predicted_ratio=float(np.max(np.abs(normalized+normalized_jac@(x-q))))
        status=str(phase2.message)
    else:
        rows=static_rows or None; rhs=static_rhs or None
        phase=linprog(g,A_ub=(None if rows is None else np.vstack(rows)),
            b_ub=(None if rhs is None else np.concatenate(rhs)),
            bounds=list(zip(lower,upper)),method="highs")
        if not phase.success:
            raise RuntimeError(f"shape objective LP failed: {phase.message}")
        total_iterations=int(phase.nit); x=phase.x; predicted_ratio=0.0
        status=str(phase.message)
    delta=np.asarray(x-q,dtype=float)
    return ShapeLPUpdate(np.asarray(x,dtype=float),delta,float(g@delta),
        predicted_ratio,bool(current_ratio>1.0),status,total_iterations)


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
         "VIMFunctionalShapeJacobian",
         "ProductionGetTrafoDisplacements","ProductionVIMLinearization","LPUpdate",
         "ShapeLinearization","ShapeLPUpdate","solve_shape_lp",
         "ElementInsertionResponse","ElementGenerationLPUpdate",
         "HDivMMMElementGenerationLinearization",
         "TSVDElementCandidateSelection",
         "HDivMMMGenerationIteration","HDivMMMGenerationResult",
         "GrowthTopologyReport","ngsolve_growth_topology",
         "finite_element_insertion_response","ngsolve_boundary_growth_candidates",
         "ngsolve_boundary_removal_candidates",
         "ngsolve_discontinuous_element_dof_blocks",
         "linearize_hdiv_mmm_element_generation","solve_hdiv_mmm_active_elements",
         "select_tsvd_element_candidates",
         "select_tsvd_exact_block_batch",
         "grow_hdiv_mmm_by_superposition","solve_element_generation_lp",
         "TopologyOptimizationResult","sample_production_gettrafo_displacements",
         "assemble_ngsolve_hdiv_shape_tangents","linearize_production_charge_gram",
         "assemble_ngsolve_hdiv_linear_form_shape_tangents",
         "assemble_ngsolve_hdiv_mass_shape_contractions",
         "linearize_production_charge_gram_matrix_free","linearize_vim_operator_matrix_free",
         "linearize_production_vim_matrix_free_from_ngsolve",
         "production_vim_rms_adjoint_gradient_streaming",
         "production_vim_functional_shape_jacobian_streaming",
         "linearize_production_vim_from_ngsolve","linearize_laplace_pair_gram",
         "affine_cell_self_energy_shape_derivative",
         "production_hex_volume_self_block_derivatives","production_hex_face_self_block_derivatives",
         "production_tet_volume_self_block_derivatives","production_tet_face_self_block_derivatives",
         "production_tet_charge_gram_derivatives",
         "production_wedge_volume_self_block_derivatives","production_wedge_face_self_block_derivatives",
         "production_wedge_charge_gram_derivatives","linearize_laplace_charge_gram",
         "linearize_vim_operator","linearize_vim_system","solve_lp_update",
         "optimize_vim_lp","write_cubit_density_journal"]
