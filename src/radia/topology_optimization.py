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
    timings_s: dict[str, float] | None = None


@dataclass(frozen=True)
class VIMStateShapeJacobian:
    """Matrix-free analytic GetTrafo derivative of the HDiv-VIM state.

    ``state_jacobian[k]`` is ``dm/dq_k``.  The derivative H-matrix is used
    only as an operator while forming ``db/dq-dA/dq*m``; no dense ``dG`` or
    dense VIM operator derivative is materialized.
    """

    state: np.ndarray
    state_jacobian: np.ndarray
    state_iterations: int
    tangent_iterations: tuple[int, ...]
    timings_s: dict | None = None


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
    reference_sampling=False
    if family=="tet":
        reference_cells=cb.get("reference_vV")
        reference_faces=cb.get("reference_bV")
        reference_sampling=(reference_cells is not None and
                            reference_faces is not None)
        cell_nodes=tuple(np.asarray(x,dtype=float) for x in
                         (reference_cells if reference_sampling else cb["vV"]))
        face_nodes=tuple(np.asarray(x,dtype=float) for x in
                         (reference_faces if reference_sampling else cb["bV"]))
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
            mapped=[mesh(*point) for point in nodes]
            values=[]
            for field in modes:
                values.append(np.array([field(point) for point in mapped],dtype=float))
            result.append(np.asarray(values,dtype=float).reshape(len(modes),len(nodes),3))
        return tuple(result)
    current_deformation=(getattr(mesh,"deformation",None)
                         if reference_sampling else None)
    if current_deformation is not None:
        mesh.UnsetDeformation()
    try:
        cells=sample(cell_nodes)
        faces=sample(face_nodes)
    finally:
        if current_deformation is not None:
            mesh.SetDeformation(current_deformation)
    return ProductionGetTrafoDisplacements(family,cells,faces)


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


def production_vim_state_shape_jacobian_streaming(*, fes,
        deformation_modes, charge_basis, charge_gram, charge_map,
        inv_chi, rhs, rhs_jacobian, family=None,
        solve_tolerance=1e-9, solve_max_iterations=5000,
        mass_riesz=True, derivative_eps=1e-8, derivative_leaf=32,
        derivative_eta=2.0, cluster_coarse_size=0,
        cluster_deflation_size=0, recycle_size=0,
        state=None, state_iterations=None):
    """Differentiate the complete HDiv-VIM state without dense ``dA``.

    This is the forward counterpart of
    :func:`production_vim_functional_shape_jacobian_streaming`.  It is used
    when a downstream nonlinear operation, such as closed-orbit tracking,
    needs the field tangent itself rather than a fixed list of adjoint
    contractions.  ``dM``, ``dB``, ``dG`` and ``drhs`` all participate in
    ``A dm = drhs-dA m``.  Each ``dG m`` action stays on a directional
    H-matrix; optimizer finite differences are not used.

    The caller owns ``ngsolve.TaskManager``.
    """
    import time
    import scipy.sparse as sp

    total_started = time.perf_counter()
    modes = tuple(deformation_modes)
    q = len(modes)
    if q == 0:
        raise ValueError("at least one deformation mode is required")
    B = sp.csr_matrix(charge_map)
    n = B.shape[1]
    b = np.asarray(rhs, dtype=float).reshape(-1)
    db = np.asarray(rhs_jacobian, dtype=float)
    if b.shape != (n,) or db.shape != (q, n):
        raise ValueError("state shape RHS/Jacobian mismatch")
    reused_state = (None if state is None else
                    np.asarray(state, dtype=float).reshape(-1))
    if reused_state is not None and (
            reused_state.shape != (n,)
            or not np.all(np.isfinite(reused_state))):
        raise ValueError("reused state must be finite and match fes.ndof")
    derivative_eps = float(derivative_eps)
    derivative_leaf = int(derivative_leaf)
    derivative_eta = float(derivative_eta)
    if (not np.isfinite(derivative_eps) or derivative_eps <= 0.0
            or derivative_leaf < 1 or not np.isfinite(derivative_eta)
            or derivative_eta <= 0.0):
        raise ValueError("invalid directional H-matrix controls")

    geometry_started = time.perf_counter()
    geometry = sample_production_gettrafo_displacements(
        fes, modes, charge_basis, family=family)
    cells = np.stack(geometry.cell, axis=1)
    if geometry.family == "wedge":
        faces = np.zeros((q, len(geometry.face), 9, 3))
        for host, values in enumerate(geometry.face):
            faces[:, host, :values.shape[1]] = values
    else:
        faces = np.stack(geometry.face, axis=1)
    geometry_s = time.perf_counter()-geometry_started

    charge_gram.restore_geometry_mass_matrix()
    state_solve_started = time.perf_counter()
    if reused_state is None:
        solved = charge_gram.solve_configured_linear_material_auto_prec_many(
            float(inv_chi), np.ascontiguousarray(b[None, :]),
            tol=float(solve_tolerance), maxit=int(solve_max_iterations),
            cluster_coarse_size=int(cluster_coarse_size),
            cluster_deflation_size=int(cluster_deflation_size),
            recycle_size=0, mass_riesz=bool(mass_riesz))
        state_value = np.asarray(solved["m"], dtype=float)[0]
        state_count = int(solved["iters"][0])
    else:
        state_value = reused_state
        state_count = (0 if state_iterations is None else
                       int(state_iterations))
    state_solve_s = time.perf_counter()-state_solve_started

    tangent_started = time.perf_counter()
    _, dmass, dcharge = assemble_ngsolve_hdiv_shape_tangents(
        fes, modes, B, sparse=True)
    bx = np.asarray(B@state_value).reshape(-1)
    Gbx = np.asarray(charge_gram.matvec_sym(bx)).reshape(-1)
    tangent_rhs = np.empty((q, n), dtype=float)
    derivative_stats = []
    for k in range(q):
        dM = sp.csr_matrix(dmass[k])
        if geometry.family == "tet":
            rates = np.asarray(
                charge_gram.tet_charge_map_row_directional_rates(
                    np.ascontiguousarray(cells[k]),
                    np.ascontiguousarray(faces[k])), dtype=float)
            dB = sp.diags(rates)@B
        else:
            dB = sp.csr_matrix(dcharge[k])
        dbx = np.asarray(dB@state_value).reshape(-1)
        derivative = charge_gram.directional_derivative_operator(
            geometry.family, np.ascontiguousarray(cells[k]),
            np.ascontiguousarray(faces[k]), derivative_eps,
            derivative_leaf, derivative_eta)
        dGbx = np.asarray(derivative.matvec_sym(bx)).reshape(-1)
        derivative_stats.append(dict(derivative.stats))
        dA_state = (
            float(inv_chi)*np.asarray(dM@state_value).reshape(-1)
            + np.asarray(dB.T@Gbx).reshape(-1)
            + np.asarray(B.T@(
                dGbx + np.asarray(charge_gram.matvec_sym(dbx)).reshape(-1)
            )).reshape(-1)
        )
        tangent_rhs[k] = db[k]-dA_state
    tangent_rhs_s = time.perf_counter()-tangent_started

    tangent_solve_started = time.perf_counter()
    solved = charge_gram.solve_configured_linear_material_auto_prec_many(
        float(inv_chi), np.ascontiguousarray(tangent_rhs),
        tol=float(solve_tolerance), maxit=int(solve_max_iterations),
        cluster_coarse_size=int(cluster_coarse_size),
        cluster_deflation_size=int(cluster_deflation_size),
        recycle_size=min(int(recycle_size), q),
        mass_riesz=bool(mass_riesz))
    state_jacobian = np.asarray(solved["m"], dtype=float)
    tangent_iterations = tuple(int(value) for value in solved["iters"])
    if state_jacobian.shape != (q, n):
        raise RuntimeError("native VIM tangent solve returned invalid shape")
    tangent_solve_s = time.perf_counter()-tangent_solve_started
    timings = {
        "geometry": geometry_s,
        "state_solve": state_solve_s,
        "tangent_rhs": tangent_rhs_s,
        "tangent_solve": tangent_solve_s,
        "directional_hmatrix_stats": derivative_stats,
        "total": time.perf_counter()-total_started,
    }
    return VIMStateShapeJacobian(
        state_value, state_jacobian, state_count,
        tangent_iterations, timings)


def production_vim_functional_shape_jacobian_streaming(*, fes,
        deformation_modes, charge_basis, charge_gram, charge_map,
        inv_chi, rhs, response_matrix, rhs_jacobian=None,
        dresponse_matrix=None, response_observations=None,
        response_weights=None, family=None, incident_response=None,
        dincident_response=None, solve_tolerance=1e-9,
        solve_max_iterations=5000, mass_riesz=True,
        cluster_coarse_size=0, cluster_deflation_size=0, recycle_size=0,
        state=None, state_iterations=None):
    """Differentiate many linear accelerator-field functionals at once.

    The state and all response adjoints share one row-major native solve.
    Directional ``dG`` terms stay on the H-matrix cluster tree; for each
    response, native support-pruned contractions evaluate every deformation
    mode without materializing a derivative matrix.  For flat TET geometry,
    ``response_observations`` and row-major vector ``response_weights`` route
    ``dC`` through the exact native configured-field derivative.  An explicit
    ``dresponse_matrix`` remains available for other element families.
    Optional cluster coarse, deflation, and recycle dimensions are forwarded
    to the shared native multi-RHS solve; their zero defaults preserve the
    measured cluster-free baseline.
    A caller that already owns a converged state may pass ``state`` and its
    optional ``state_iterations``.  Only response adjoints are then solved at
    ``solve_tolerance``; the high-accuracy state is reused exactly.
    """
    import time
    import scipy.sparse as sp
    total_started=time.perf_counter()
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

    # Geometry sampling is independent of the state/adjoint solutions.  Do it
    # before the expensive multi-RHS solve so an invalid resumed GetTrafo
    # mapping fails immediately instead of after a long Krylov run.
    geometry_started=time.perf_counter()
    geometry=sample_production_gettrafo_displacements(
        fes,modes,charge_basis,family=family)
    geometry_s=time.perf_counter()-geometry_started
    cells=np.stack(geometry.cell,axis=1)
    if geometry.family=="wedge":
        faces=np.zeros((q,len(geometry.face),9,3))
        for host,values in enumerate(geometry.face):
            faces[:,host,:values.shape[1]]=values
    else: faces=np.stack(geometry.face,axis=1)

    reused_state=None if state is None else np.asarray(
        state,dtype=float).reshape(-1)
    if reused_state is not None and (
            reused_state.shape!=(n,) or not np.all(np.isfinite(reused_state))):
        raise ValueError("reused functional shape state must be finite and match fes.ndof")
    charge_gram.restore_geometry_mass_matrix()
    right_hand_sides=np.ascontiguousarray(
        C if reused_state is not None else np.vstack((b,C)),dtype=np.float64)
    cluster_coarse_size=int(cluster_coarse_size)
    cluster_deflation_size=int(cluster_deflation_size)
    recycle_size=int(recycle_size)
    if min(cluster_coarse_size,cluster_deflation_size,recycle_size)<0:
        raise ValueError("cluster solver dimensions must be non-negative")
    solve_started=time.perf_counter()
    solved=charge_gram.solve_configured_linear_material_auto_prec_many(
        float(inv_chi),right_hand_sides,tol=float(solve_tolerance),
        maxit=int(solve_max_iterations),
        cluster_coarse_size=cluster_coarse_size,
        cluster_deflation_size=cluster_deflation_size,
        recycle_size=min(recycle_size,len(right_hand_sides)),
        mass_riesz=bool(mass_riesz))
    solve_s=time.perf_counter()-solve_started
    solutions=np.asarray(solved["m"],dtype=float)
    iterations=tuple(int(value) for value in solved["iters"])
    expected_solutions=nout if reused_state is not None else nout+1
    if solutions.shape!=(expected_solutions,n) or len(iterations)!=expected_solutions:
        raise RuntimeError("native functional shape solve returned invalid shapes")
    if reused_state is None:
        state=solutions[0];adjoints=solutions[1:]
        state_iteration_count=int(iterations[0])
        adjoint_iteration_counts=tuple(iterations[1:])
    else:
        state=reused_state;adjoints=solutions
        state_iteration_count=(0 if state_iterations is None else
                               int(state_iterations))
        adjoint_iteration_counts=tuple(iterations)
    response=C@state+incident

    bx=np.asarray(B@state).reshape(-1)
    Gbx=np.asarray(charge_gram.matvec_sym(bx)).reshape(-1)
    badjoints=[np.asarray(B@value).reshape(-1) for value in adjoints]
    cell_velocities=np.ascontiguousarray(cells)
    face_velocities=np.ascontiguousarray(faces)
    response_derivative_started=time.perf_counter()
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
    response_derivative_s=time.perf_counter()-response_derivative_started
    left_matrix=np.ascontiguousarray(np.stack(badjoints),dtype=float)
    gram_started=time.perf_counter()
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
    gram_s=time.perf_counter()-gram_started
    mass_started=time.perf_counter()
    mass_terms=assemble_ngsolve_hdiv_mass_shape_contractions(
        fes,modes,adjoints,state)
    mass_s=time.perf_counter()-mass_started
    post_started=time.perf_counter()
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
    post_s=time.perf_counter()-post_started
    timings={
        "geometry":geometry_s,
        "solve":solve_s,
        "response_derivative":response_derivative_s,
        "gram_contractions":gram_s,
        "mass_contractions":mass_s,
        "postprocess":post_s,
        "total":time.perf_counter()-total_started,
    }
    return VIMFunctionalShapeJacobian(state,response,jacobian,
        state_iteration_count,adjoint_iteration_counts,timings)


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
    candidate_coupling_rank: int
    candidate_coupling_relative_truncation_error: float


@dataclass(frozen=True)
class HDivMMMBlockInsertionResponse:
    """Exact fixed-active-set response of a finite candidate-element bundle."""
    selected_elements: np.ndarray
    candidate_state: np.ndarray
    response_delta: np.ndarray
    schur_complement: np.ndarray


@dataclass(frozen=True)
class HDivMMMSingleRemovalResponses:
    """All exact singleton deletions from one reduced candidate Schur solve.

    ``positive_material_response[:, j]`` is the response of candidate element
    ``j`` that disappears when that element is removed from the full selected
    set.  ``removed_response_delta[:, j]`` is the corresponding response,
    relative to the retained active base, after the deletion.  Both are exact
    for the already-reduced fixed active set; no new H-matrix solve is used.
    """

    candidate_elements: np.ndarray
    full_candidate_state: np.ndarray
    full_response_delta: np.ndarray
    positive_material_response: np.ndarray
    removed_response_delta: np.ndarray


@dataclass(frozen=True)
class HDivMMMRemovalGroupResponses:
    """Exact deletion responses for arbitrary groups from one Schur inverse."""

    candidate_elements: np.ndarray
    removal_groups: tuple[np.ndarray, ...]
    full_candidate_state: np.ndarray
    full_response_delta: np.ndarray
    positive_material_response: np.ndarray
    removed_response_delta: np.ndarray


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
class AbeMurataEquivalentMaterialDiagnostics:
    """Continuous DUCAS/TSVD predictor behind a binary material proposal.

    Murata and Abe compute a magnetizing-current distribution for the error
    field, convert it to an equivalent ferromagnetic volume, update the mesh,
    and repeat the field solve.  Here one source column is one complete
    candidate element (or one coupled removal group).  The continuous TSVD
    coefficients are therefore fractions of those material columns and
    ``equivalent_volume_changes`` is their signed physical-volume analogue.

    The response rows are divided by their engineering bands before the SVD,
    because an accelerator objective mixes field and transfer-map units.
    Consequently ``normalized_mode_field_strengths`` is Abe's per-mode field
    strength applied to this dimensionless, band-normalized error field.
    The values are a proposal diagnostic; committed geometry remains binary
    and is accepted only after an exact active-system re-solve.
    """
    candidate_elements: np.ndarray
    candidate_material_active: np.ndarray
    retained_rank: int
    singular_values: np.ndarray
    normalized_mode_field_strengths: np.ndarray
    mode_material_amplitudes: np.ndarray
    signed_material_fractions: np.ndarray
    equivalent_volume_changes: np.ndarray
    projected_normalized_correction: np.ndarray


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
    linearized_reachability_residual: np.ndarray | None = None
    linearized_reachability_max_band_ratio: float = float("inf")
    linearized_reachability_relative_residual: float = float("inf")
    linearized_reachable: bool = False
    abe_murata_diagnostics: AbeMurataEquivalentMaterialDiagnostics | None = None


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
    candidate_coupling_rank: int = 0
    candidate_coupling_relative_truncation_error: float = 0.0
    candidate_schur_max_iterations: int = 0
    response_adjoint_count: int = 0
    linearized_reachability_residual: np.ndarray | None = None
    linearized_reachability_max_band_ratio: float = float("inf")
    linearized_reachability_relative_residual: float = float("inf")
    linearized_reachable: bool = False
    material_trust_volume_before: float | None = None
    material_trust_volume_after: float | None = None
    material_changed_volume: float = 0.0
    graph_front_proposals_evaluated: int = 0
    nonmonotone_search_depth: int = 0
    abe_murata_diagnostics: AbeMurataEquivalentMaterialDiagnostics | None = None


@dataclass(frozen=True)
class HDivMMMGraphFrontDiagnostics:
    iteration: int
    search_depth: int
    candidate_count: int
    novelty_weight: float
    pool_proposal_count: int
    pool_response_rank: int
    pool_duplicate_pair_fraction: float
    pool_maximum_absolute_correlation: float
    selected_proposal_count: int
    selected_response_rank: int
    selected_duplicate_pair_fraction: float
    selected_minimum_subspace_novelty: float


@dataclass(frozen=True)
class HDivMMMExactSearchTrial:
    """One topology-valid beam state evaluated by a complete active solve."""

    depth: int
    parent_max_band_ratio: float
    incumbent_max_band_ratio: float
    max_band_ratio: float
    added_elements: np.ndarray
    removed_elements: np.ndarray
    solve_iterations: int
    path: tuple[tuple[tuple[int, int], ...], ...]


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
    exact_search_trace: tuple[HDivMMMExactSearchTrial, ...] = ()
    graph_front_diagnostics: tuple[HDivMMMGraphFrontDiagnostics, ...] = ()


@dataclass(frozen=True)
class _HDivMMMExactBeamState:
    active: np.ndarray
    state: np.ndarray
    response: np.ndarray
    objective_response: np.ndarray
    ratio: float
    source_scale: float
    solve_iterations: int
    depth: int
    path: tuple[tuple[tuple[int, int], ...], ...]


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
        mass_riesz=True, cluster_coarse_size=0,
        cluster_deflation_size=0, recycle_size=0,
        candidate_selector=None, active_state=None,
        screen_with_adjoint=True,
        screen_adjoint_rows=None,
        candidate_screen_context=None) -> HDivMMMElementGenerationLinearization:
    """Close boundary element generation on the configured H-matrix operator.

    State and response adjoints are solved once on the active iron.  For every
    selected candidate RT/BDM block, a fused native call constructs full-strength
    columns of ``A``, performs the constrained active solves, and returns only
    the reduced Schur/RHS/response arrays.  When ``candidate_selector`` is
    supplied, the native kernel contracts each exact element-local block from
    its sparse mass/charge supports after one shared state/adjoint H-matrix
    application.  Older binaries use a one-direction Ritz screen on production
    fronts and retain complete blocks on tiny fronts.  The callback is called as
    ``selector(elements, approximate_delta, state, response)`` and returns the
    element indices that receive the exact fused reduction.  An empty selection
    is permitted when the global signed proposal will instead be checked by one
    complete active-set solve.  The full dense ``A`` and ``N`` matrices are
    never materialized.
    """
    blocks=ngsolve_discontinuous_element_dof_blocks(fes)
    if (candidate_screen_context is not None and
            not isinstance(candidate_screen_context,dict)):
        raise ValueError("candidate_screen_context must be a dictionary")
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
    coarse_size=int(cluster_coarse_size)
    deflation_size=int(cluster_deflation_size)
    recycle=int(recycle_size)
    if coarse_size<0 or deflation_size<0 or recycle<0:
        raise ValueError("cluster-tree solver sizes must be nonnegative")
    active_dofs=np.concatenate([blocks[k] for k in np.flatnonzero(active_el)]).astype(np.int32)
    active_mask=np.zeros(n,dtype=bool); active_mask[active_dofs]=True
    inactive_dofs=np.flatnonzero(~active_mask).astype(np.int32)
    charge_gram.set_configured_constraints(inactive_dofs,preserve_existing=False)

    def solve_many(rows):
        data=np.ascontiguousarray(rows,dtype=np.float64)
        if data.ndim!=2 or data.shape[1]!=n or not np.all(np.isfinite(data)):
            raise ValueError(
                "HDiv-MMM multi-RHS rows must be a finite (nrhs, ndof) array")
        # The state and all response adjoints share the same active matrix.
        # Mass-Riesz is the production default for broken BDM.  The optional
        # H-matrix cluster coarse space is explicit: on the end-pack system it
        # increased six-RHS iteration counts by roughly an order of magnitude.
        use_cluster=(data.shape[0]>1 and coarse_size>0 and
                     n>=max(1024,4*coarse_size))
        result=charge_gram.solve_configured_linear_material_auto_prec_many(
            float(inv_chi),data,tol=float(solve_tolerance),
            maxit=int(solve_max_iterations),
            cluster_coarse_size=(coarse_size if use_cluster else 0),
            cluster_deflation_size=(deflation_size if use_cluster else 0),
            recycle_size=(min(recycle,data.shape[0]) if use_cluster else 0),
            # The native solver deliberately treats mass-Riesz and the
            # cluster-tree two-level preconditioner as alternative policies.
            mass_riesz=(False if use_cluster else bool(mass_riesz)))
        solutions=np.asarray(result["m"],dtype=float)
        iterations=tuple(int(x) for x in result["iters"])
        if (solutions.shape!=data.shape or len(iterations)!=data.shape[0]
                or not np.all(np.isfinite(solutions))):
            raise RuntimeError("native HDiv-MMM multi-RHS solve returned invalid shapes")
        if hasattr(charge_gram,
                   "apply_configured_linear_material_operator_many"):
            applied=np.asarray(
                charge_gram.apply_configured_linear_material_operator_many(
                    float(inv_chi),np.ascontiguousarray(solutions),
                    respect_constraints=True),dtype=float)
        else:
            applied=np.stack([
                np.asarray(charge_gram.apply_configured_linear_material_operator(
                    float(inv_chi),row,respect_constraints=True),dtype=float)
                for row in solutions])
        if applied.shape!=data.shape or not np.all(np.isfinite(applied)):
            raise RuntimeError(
                "native HDiv-MMM convergence check returned invalid shapes")
        denominator=np.maximum(
            np.linalg.norm(data,axis=1),np.finfo(float).tiny)
        relative=np.linalg.norm(data-applied,axis=1)/denominator
        residual_limit=max(5.0*float(solve_tolerance),1.0e-12)
        failed=np.flatnonzero(relative>residual_limit)
        if failed.size:
            evidence=", ".join(
                f"row {int(row)}={relative[row]:.3e}"
                for row in failed[:8])
            raise RuntimeError(
                "native HDiv-MMM multi-RHS solve did not meet its checked "
                f"relative residual ({evidence}; limit={residual_limit:.3e})")
        return solutions,iterations

    provided_state=(None if active_state is None else
                    np.asarray(active_state,dtype=float).reshape(-1))
    if (provided_state is not None and
            (provided_state.shape!=(n,) or
             not np.all(np.isfinite(provided_state)) or
             np.any(provided_state[inactive_dofs]!=0.0))):
        raise ValueError(
            "active_state must be finite, configured-size, and zero outside active DOFs")
    screen_with_adjoint=bool(screen_with_adjoint)
    if screen_with_adjoint:
        adjoint_rows=(np.arange(C.shape[0],dtype=np.int64)
            if screen_adjoint_rows is None else
            np.asarray(screen_adjoint_rows,dtype=np.int64).reshape(-1))
        if (np.unique(adjoint_rows).size!=adjoint_rows.size or
                np.any(adjoint_rows<0) or np.any(adjoint_rows>=C.shape[0])):
            raise ValueError("screen_adjoint_rows must be unique response rows")
        if adjoint_rows.size==0:
            screen_with_adjoint=False
    if not screen_with_adjoint:
        adjoint_rows=np.empty(0,dtype=np.int64)
    if provided_state is None:
        base_rhs=(np.vstack((rhs_full,C[adjoint_rows]))
                  if screen_with_adjoint else
                  rhs_full[None,:].copy())
        base_rhs[:,inactive_dofs]=0.0
        solved,base_iterations=solve_many(base_rhs)
        state=solved[0]
        adjoints=(solved[1:].T if screen_with_adjoint else None)
        state_iterations=int(base_iterations[0])
        adjoint_iterations=(tuple(base_iterations[1:])
                            if screen_with_adjoint else tuple())
    else:
        state=provided_state.copy();state_iterations=0
        if screen_with_adjoint:
            adjoint_rhs=C[adjoint_rows].copy()
            adjoint_rhs[:,inactive_dofs]=0.0
            solved_adjoint,adjoint_iterations=solve_many(adjoint_rhs)
            adjoints=solved_adjoint.T
        else:
            adjoints=None;adjoint_iterations=tuple()
    response=C@state
    if incident_response is not None:
        incident=np.asarray(incident_response,dtype=float).reshape(-1)
        if incident.shape!=response.shape:
            raise ValueError("incident_response shape mismatch")
        response=response+incident
    if candidate_screen_context is not None:
        candidate_screen_context.clear()
        candidate_screen_context.update(
            state=np.asarray(state,dtype=float).copy(),
            adjoints=(None if adjoints is None else
                      np.asarray(adjoints,dtype=float).copy()),
            adjoint_rows=np.asarray(adjoint_rows,dtype=np.int64).copy())

    available_candidate_count=int(candidates.size)

    def apply_operator_rows(rows,*,respect_constraints=False):
        rows=np.ascontiguousarray(rows,dtype=np.float64)
        if rows.ndim!=2 or rows.shape[1]!=n or not np.all(np.isfinite(rows)):
            raise ValueError(
                "HDiv-MMM operator rows must be a finite (nrhs, ndof) array")
        if hasattr(charge_gram,
                   "apply_configured_linear_material_operator_many"):
            applied=np.asarray(
                charge_gram.apply_configured_linear_material_operator_many(
                    float(inv_chi),rows,
                    respect_constraints=bool(respect_constraints)),dtype=float)
        else:  # Compatibility with an older downloaded native binary.
            applied=np.stack([
                np.asarray(charge_gram.apply_configured_linear_material_operator(
                    float(inv_chi),row,
                    respect_constraints=bool(respect_constraints)),dtype=float)
                for row in rows])
        if applied.shape!=rows.shape:
            raise RuntimeError("native HDiv-MMM operator batch has invalid shape")
        return applied

    def apply_operator_many(dofs):
        dofs=np.asarray(dofs,dtype=np.int32).reshape(-1)
        basis=np.zeros((len(dofs),n),dtype=np.float64)
        basis[np.arange(len(dofs)),dofs]=1.0
        return apply_operator_rows(basis,respect_constraints=False)

    def contracted_local_response(block,active_coupling=None,
                                  unconstrained_base=None):
        local=C[:,block].copy()
        if adjoints is None:
            return local
        if unconstrained_base is not None:
            local[adjoint_rows]-=unconstrained_base[1:,block]
        else:
            local[adjoint_rows]-=(
                adjoints[active_dofs].T@active_coupling)
        return local

    if candidate_selector is not None:
        approximate=[]
        if hasattr(charge_gram,
                   "configured_linear_material_element_blocks"):
            # The active/candidate couplings for every element follow from one
            # state/adjoint batch.  C++ then assembles each exact A_ee from its
            # local sparse mass and charge supports, without any candidate-wise
            # global H-matrix application.
            unconstrained_rows=(state[None,:] if adjoints is None else
                np.vstack((state,adjoints.T)))
            unconstrained_base=apply_operator_rows(
                np.ascontiguousarray(unconstrained_rows),
                respect_constraints=False)
            candidate_blocks=tuple(blocks[int(k)] for k in candidates)
            candidate_dofs=np.concatenate(candidate_blocks).astype(np.int32)
            offsets=np.asarray(np.r_[0,np.cumsum(
                [len(block) for block in candidate_blocks])],dtype=np.int32)
            packed=np.asarray(
                charge_gram.configured_linear_material_element_blocks(
                    float(inv_chi),candidate_dofs,offsets),dtype=float)
            packed_offset=0
            for block in candidate_blocks:
                width=len(block);next_offset=packed_offset+width*width
                local_matrix=packed[packed_offset:next_offset].reshape(
                    width,width)
                packed_offset=next_offset
                local_rhs=rhs_full[block]-unconstrained_base[0,block]
                local_response=contracted_local_response(
                    block,unconstrained_base=unconstrained_base)
                try:
                    local_state=np.linalg.solve(local_matrix,local_rhs)
                    approximate.append(local_response@local_state)
                except np.linalg.LinAlgError:
                    approximate.append(np.full(C.shape[0],np.nan))
            if packed_offset!=packed.size:
                raise RuntimeError(
                    "native HDiv-MMM local block packing is inconsistent")
        elif len(candidates)<=8:
            # On a tiny front the complete element-local block is both cheap
            # and useful for exact collaborative regression problems.
            begin=0
            while begin<len(candidates):
                end=begin;count=0
                while end<len(candidates):
                    width=len(blocks[int(candidates[end])])
                    if end>begin and count+width>batch_size:
                        break
                    count+=width;end+=1
                    if count>=batch_size:
                        break
                group_blocks=tuple(
                    blocks[int(k)] for k in candidates[begin:end])
                group_dofs=np.concatenate(group_blocks).astype(np.int32)
                applied=apply_operator_many(group_dofs)
                local_offsets=np.r_[0,np.cumsum(
                    [len(block) for block in group_blocks])]
                for local_column,block in enumerate(group_blocks):
                    local=np.arange(local_offsets[local_column],
                                    local_offsets[local_column+1])
                    columns=applied[local].T
                    active_coupling=columns[active_dofs,:]
                    local_matrix=0.5*(
                        columns[block,:]+columns[block,:].T)
                    local_rhs=(rhs_full[block]
                               -active_coupling.T@state[active_dofs])
                    local_response=contracted_local_response(
                        block,active_coupling=active_coupling)
                    try:
                        local_state=np.linalg.solve(local_matrix,local_rhs)
                        approximate.append(local_response@local_state)
                    except np.linalg.LinAlgError:
                        approximate.append(np.full(C.shape[0],np.nan))
                begin=end
        else:
            # Compatibility with downloaded binaries predating the local-block
            # kernel: production fronts are screened by one Ritz direction per
            # element, not by extracting every A[:,e_j] column.  Because the
            # active state and adjoints vanish on inactive DOFs, a single
            # unconstrained application gives the insertion residual
            # r_e=b_e-A_ea*m_a and the adjoint contraction for every candidate.
            # The element-local steepest direction is scaled by its exact
            # H-matrix Rayleigh quotient.  Thus an order-p BDM block costs one
            # shared-H-matrix RHS rather than one RHS per local DOF.  Every
            # proposed signed batch is still checked by a full active-set solve.
            unconstrained_rows=(state[None,:] if adjoints is None else
                np.vstack((state,adjoints.T)))
            unconstrained_base=apply_operator_rows(
                np.ascontiguousarray(unconstrained_rows),
                respect_constraints=False)
            for begin in range(0,len(candidates),batch_size):
                group=candidates[begin:begin+batch_size]
                group_blocks=tuple(blocks[int(k)] for k in group)
                directions=np.zeros((len(group_blocks),n),dtype=np.float64)
                residuals=[];responses=[]
                for local_column,block in enumerate(group_blocks):
                    residual=(rhs_full[block]
                              -unconstrained_base[0,block])
                    residuals.append(residual)
                    responses.append(contracted_local_response(
                        block,unconstrained_base=unconstrained_base))
                    directions[local_column,block]=residual
                applied=apply_operator_rows(
                    directions,respect_constraints=False)
                for local_column,block in enumerate(group_blocks):
                    residual=residuals[local_column]
                    numerator=float(residual@residual)
                    denominator=float(
                        residual@applied[local_column,block])
                    scale_floor=np.finfo(float).eps*max(1.0,numerator)
                    if (not np.isfinite(denominator) or
                            denominator<=scale_floor or numerator==0.0):
                        approximate.append(np.full(C.shape[0],np.nan))
                        continue
                    local_state=(numerator/denominator)*residual
                    approximate.append(
                        responses[local_column]@local_state)
        approximate_delta=np.stack(approximate,axis=1)
        selected=np.asarray(candidate_selector(
            candidates,approximate_delta,state,response),dtype=np.int64).reshape(-1)
        if np.unique(selected).size!=selected.size:
            raise ValueError("candidate_selector must return a unique subset")
        available=set(int(k) for k in candidates)
        if any(int(k) not in available for k in selected):
            raise ValueError("candidate_selector returned an unavailable element")
        candidates=selected

    if candidates.size==0:
        empty=np.empty(0,dtype=np.int64)
        return HDivMMMElementGenerationLinearization(
            state,response,empty,tuple(),
            np.empty((C.shape[0],0),dtype=float),tuple(),
            np.zeros(1,dtype=np.int32),np.empty((0,0),dtype=float),
            np.empty(0,dtype=float),
            np.empty((C.shape[0],0),dtype=float),available_candidate_count,
            {"operator_s":0.0,"solve_s":0.0,"contraction_s":0.0},
            int(state_iterations),tuple(adjoint_iterations),tuple(),0,0.0)

    if adjoints is None or adjoint_rows.size!=C.shape[0]:
        raise RuntimeError(
            "partial-adjoint global screening must return an empty exact Schur front")

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
        schur_iterations=tuple(int(x) for x in reduced.get(
            "coupling_mode_iters",reduced["iters"]))
        nc=len(candidate_dofs)
        coupling_rank=int(reduced.get("coupling_rank",nc))
        coupling_error=float(reduced.get(
            "coupling_relative_truncation_error",0.0))
        native_timings={key:float(reduced[key]) for key in native_timings}
        if (reduced_schur.shape!=(nc,nc) or reduced_rhs.shape!=(nc,)
                or reduced_response.shape!=(C.shape[0],nc)
                or coupling_rank<0 or coupling_rank>nc
                or len(schur_iterations)!=coupling_rank
                or not np.isfinite(coupling_error)):
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
        coupling_rank=len(candidate_dofs)
        coupling_error=0.0

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
        native_timings,int(state_iterations),
        tuple(adjoint_iterations),tuple(schur_iterations),
        int(coupling_rank),float(coupling_error))


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


def hdiv_mmm_removal_group_responses(linearization,removal_groups,*,
        full_elements=None) -> HDivMMMRemovalGroupResponses:
    """Evaluate arbitrary element-group deletions with one Schur factorization.

    Let ``S`` be the candidate Schur complement, ``x=S^-1 b``, and
    ``Q=S^-1``.  Constraining one element block ``e`` to zero changes the
    candidate state by

    ``delta x = -Q[:,e] Q[e,e]^-1 x[e]``.

    Therefore all requested removal-group responses follow after one solve
    with right-hand sides ``[b, I]`` plus one local block solve per group.
    This replaces repeatedly factoring the nearly identical
    ``S_without_group`` matrices and is the intended exact binary oracle after
    ACA--QR--TSVD.
    """
    if not isinstance(linearization,HDivMMMElementGenerationLinearization):
        raise TypeError(
            "linearization must be HDivMMMElementGenerationLinearization")
    available=np.asarray(
        linearization.candidate_elements,dtype=np.int64).reshape(-1)
    selected=(available.copy() if full_elements is None else
              np.asarray(full_elements,dtype=np.int64).reshape(-1))
    if (selected.size==0 or np.unique(selected).size!=selected.size):
        raise ValueError("full_elements must contain a non-empty unique set")
    lookup={int(element):column for column,element in enumerate(available)}
    if any(int(element) not in lookup for element in selected):
        raise ValueError("full_elements contains an unavailable candidate")
    try:
        groups=tuple(np.asarray(group,dtype=np.int64).reshape(-1)
                     for group in removal_groups)
    except TypeError as error:
        raise ValueError("removal_groups must be an iterable of groups") from error
    selected_set=set(map(int,selected))
    if (not groups or any(group.size==0 or
            np.unique(group).size!=group.size or
            any(int(element) not in selected_set for element in group)
            for group in groups)):
        raise ValueError(
            "removal_groups must contain non-empty unique subsets of "
            "full_elements")
    offsets=np.asarray(
        linearization.candidate_dof_offsets,dtype=np.int64).reshape(-1)
    if offsets.shape!=(available.size+1,):
        raise ValueError("candidate DOF offsets are inconsistent")
    source_blocks=[np.arange(
        int(offsets[lookup[int(element)]]),
        int(offsets[lookup[int(element)]+1]),dtype=np.int64)
        for element in selected]
    if any(block.size==0 for block in source_blocks):
        raise ValueError("every selected candidate must own a non-empty block")
    source=np.concatenate(source_blocks)
    schur_all=np.asarray(
        linearization.reduced_schur_complement,dtype=float)
    rhs_all=np.asarray(
        linearization.reduced_schur_rhs,dtype=float).reshape(-1)
    response_all=np.asarray(
        linearization.reduced_response_matrix,dtype=float)
    if (schur_all.shape!=(offsets[-1],offsets[-1]) or
            rhs_all.shape!=(offsets[-1],) or
            response_all.ndim!=2 or response_all.shape[1]!=offsets[-1]):
        raise ValueError("reduced candidate Schur arrays are inconsistent")
    schur=np.ascontiguousarray(schur_all[np.ix_(source,source)])
    rhs=np.ascontiguousarray(rhs_all[source])
    response=np.ascontiguousarray(response_all[:,source])
    width=source.size
    try:
        solved=np.linalg.solve(
            schur,np.column_stack((rhs,np.eye(width,dtype=float))))
    except np.linalg.LinAlgError as error:
        raise RuntimeError(
            "candidate Schur complement is singular") from error
    state=np.asarray(solved[:,0],dtype=float)
    inverse=np.asarray(solved[:,1:],dtype=float)
    full_delta=response@state
    local_offsets=np.r_[0,np.cumsum([block.size for block in source_blocks])]
    selected_lookup={int(element):column
                     for column,element in enumerate(selected)}
    positive=[]
    for group in groups:
        local=np.concatenate([
            np.arange(int(local_offsets[selected_lookup[int(element)]]),
                      int(local_offsets[selected_lookup[int(element)]+1]))
            for element in group])
        inverse_local=inverse[np.ix_(local,local)]
        try:
            multiplier=np.linalg.solve(inverse_local,state[local])
        except np.linalg.LinAlgError as error:
            raise RuntimeError(
                "candidate inverse principal block is singular") from error
        constrained_correction=inverse[:,local]@multiplier
        positive.append(response@constrained_correction)
    positive=np.ascontiguousarray(np.column_stack(positive))
    removed=np.ascontiguousarray(full_delta[:,None]-positive)
    return HDivMMMRemovalGroupResponses(
        candidate_elements=np.ascontiguousarray(selected),
        removal_groups=tuple(np.ascontiguousarray(group.copy())
                             for group in groups),
        full_candidate_state=np.ascontiguousarray(state),
        full_response_delta=np.ascontiguousarray(full_delta),
        positive_material_response=positive,
        removed_response_delta=removed)


def hdiv_mmm_all_single_removal_responses(linearization,
        full_elements=None) -> HDivMMMSingleRemovalResponses:
    """Evaluate every one-element deletion with one dense Schur factorization."""
    available=np.asarray(
        linearization.candidate_elements,dtype=np.int64).reshape(-1)
    selected=(available.copy() if full_elements is None else
              np.asarray(full_elements,dtype=np.int64).reshape(-1))
    grouped=hdiv_mmm_removal_group_responses(
        linearization,tuple(np.asarray([element],dtype=np.int64)
                            for element in selected),
        full_elements=selected)
    return HDivMMMSingleRemovalResponses(
        candidate_elements=grouped.candidate_elements,
        full_candidate_state=grouped.full_candidate_state,
        full_response_delta=grouped.full_response_delta,
        positive_material_response=grouped.positive_material_response,
        removed_response_delta=grouped.removed_response_delta)


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
        candidate_volume_changes=None, candidate_material_active=None,
        candidate_exclusion_groups=None, maximum_changed_volume=None,
        maximum_changed_elements=None, candidate_secondary_cost=None
        ) -> TSVDElementCandidateSelection:
    """Select a global binary proposal from the TSVD of *all* candidates.

    Rows are normalized by their engineering bands before the decomposition,
    so a field row with large SI units cannot hide an optics row with small SI
    units.  TSVD removes response directions below ``relative_tolerance``.  A
    first 0--1 minimax solve finds the best truncated response under the volume
    and predecessor constraints; a second solve chooses the minimum-volume
    set that captures ``improvement_capture`` of that best reduction.
    ``maximum_changed_volume`` is the TOBS/SAIP trust radius on total material
    flipped, distinct from the signed net ``volume_budget``.  It prevents a
    large addition and a large removal from cancelling in the volume row while
    producing a physically remote, poorly predicted proposal.
    When ``candidate_material_active`` is supplied, the retained TSVD
    pseudoinverse estimates signed material coefficients: positive coefficients
    activate currently inactive candidates and negative coefficients deactivate
    currently active candidates.  Infeasible opposite signs are discarded.

    ``candidate_exclusion_groups`` labels mutually exclusive alternatives.
    For example, retreating one end-contour column by one, two, or three packs
    defines three alternative terminal states rather than three independent
    material moves.  Negative labels leave candidates unconstrained.

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
    exclusion_groups=(None if candidate_exclusion_groups is None else
                      np.asarray(candidate_exclusion_groups,dtype=np.int64).reshape(-1))
    secondary_cost=(None if candidate_secondary_cost is None else
                    np.asarray(candidate_secondary_cost,dtype=float).reshape(-1))
    if (current.shape!=target.shape or target.shape!=band.shape or
            target.size==0 or np.any(band<=0.0)):
        raise ValueError("TSVD response, target, and band vectors must match")
    if (delta.shape!=(target.size,elements.size) or
            volumes.shape!=elements.shape or volume_changes.shape!=elements.shape or
            elements.size==0 or (material_active is not None and
                                  material_active.shape!=elements.shape) or
            (exclusion_groups is not None and
             exclusion_groups.shape!=elements.shape) or
            (secondary_cost is not None and
             secondary_cost.shape!=elements.shape) or
            np.unique(elements).size!=elements.size or np.any(volumes<=0.0)):
        raise ValueError("TSVD candidate arrays have incompatible shapes")
    relative_tolerance=float(relative_tolerance)
    improvement_capture=float(improvement_capture)
    if (not 0.0<relative_tolerance<1.0 or
            not 0.0<improvement_capture<=1.0 or
            not np.all(np.isfinite(delta)) or
            (secondary_cost is not None and
             not np.all(np.isfinite(secondary_cost)))):
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
    correction=(target-current)/band
    mode_projection=U.T@correction
    normalized_mode_field_strengths=(mode_projection/
        np.sqrt(float(correction.size)))
    diagnostic_material_active=(
        np.zeros(elements.size,dtype=bool) if material_active is None else
        material_active.copy())

    def abe_murata_diagnostics(retained_rank):
        retained_rank=int(retained_rank)
        mode_material_amplitudes=np.zeros_like(mode_projection)
        if retained_rank:
            mode_material_amplitudes[:retained_rank]=(
                mode_projection[:retained_rank]/singular[:retained_rank])
            signed_material_fractions=(
                V[:,:retained_rank]@mode_material_amplitudes[:retained_rank])
            projected=(
                U[:,:retained_rank]@mode_projection[:retained_rank])
        else:
            signed_material_fractions=np.zeros(elements.size,dtype=float)
            projected=np.zeros_like(correction)
        return AbeMurataEquivalentMaterialDiagnostics(
            candidate_elements=elements.copy(),
            candidate_material_active=diagnostic_material_active.copy(),
            retained_rank=retained_rank,
            singular_values=np.asarray(singular,dtype=float).copy(),
            normalized_mode_field_strengths=np.asarray(
                normalized_mode_field_strengths,dtype=float).copy(),
            mode_material_amplitudes=np.asarray(
                mode_material_amplitudes,dtype=float).copy(),
            signed_material_fractions=np.asarray(
                signed_material_fractions,dtype=float).copy(),
            equivalent_volume_changes=np.asarray(
                signed_material_fractions*volumes,dtype=float),
            projected_normalized_correction=np.asarray(
                projected,dtype=float).copy())

    residual_scale=max(1.0,float(np.linalg.norm((current-target)/band)))
    if (singular.size==0 or not np.isfinite(singular[0]) or
            singular[0]<=relative_tolerance*residual_scale):
        diagnostics=abe_murata_diagnostics(0)
        reachability_residual=current-target
        reachability_ratio=float(np.max(np.abs(reachability_residual/band)))
        reachability_relative=float(np.linalg.norm(reachability_residual/band)/
            max(np.finfo(float).tiny,
                float(np.linalg.norm((target-current)/band))))
        return TSVDElementCandidateSelection(
            np.empty(0,dtype=np.int64),np.empty(0,dtype=np.int8),
            np.empty(0,dtype=np.int64),np.empty(0,dtype=np.int8),current.copy(),
            float(np.max(np.abs((current-target)/band))),0.0,0,
            int(factor.k_aca),np.asarray(singular,dtype=float),
            np.zeros(elements.size),0.0,
            "all normalized candidate responses are zero",
            reachability_residual,reachability_ratio,reachability_relative,
            reachability_ratio<=1.0+float(ratio_tolerance),
            abe_murata_diagnostics=diagnostics)
    rank=max(1,int(np.count_nonzero(
        singular>=relative_tolerance*singular[0])))
    retained=(U[:,:rank]*singular[:rank])@V[:,:rank].T
    projected_correction=U[:,:rank]@(U[:,:rank].T@correction)
    normalized_reachability_residual=projected_correction-correction
    reachability_residual=band*normalized_reachability_residual
    reachability_ratio=float(np.max(np.abs(normalized_reachability_residual)))
    reachability_relative=float(np.linalg.norm(
        normalized_reachability_residual)/max(
            np.finfo(float).tiny,float(np.linalg.norm(correction))))
    linearized_reachable=reachability_ratio<=1.0+float(ratio_tolerance)
    diagnostics=abe_murata_diagnostics(rank)
    coefficients=diagnostics.signed_material_fractions
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
            relative_error,"TSVD magnetization signs admit no feasible boundary move",
            reachability_residual,reachability_ratio,reachability_relative,
            linearized_reachable,abe_murata_diagnostics=diagnostics)
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
    move_exclusion_groups=(None if exclusion_groups is None else
                           exclusion_groups[feasible_columns])
    move_secondary_cost=(None if secondary_cost is None else
                         secondary_cost[feasible_columns])

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
        volume_budget=float(volume_budget),
        maximum_new_elements=maximum_changed_elements,
        whole_elements=True,predecessor_pairs=predecessor_pairs,
        candidate_objective_change=move_secondary_cost,
        candidate_volume_change=move_volume_changes,
        candidate_exclusion_groups=move_exclusion_groups,
        maximum_changed_volume=maximum_changed_volume)
    current_ratio=float(np.max(np.abs((current-target)/band)))
    if (not np.any(best.selected) or
            best.predicted_max_band_ratio>=current_ratio-float(ratio_tolerance)):
        return TSVDElementCandidateSelection(
            np.empty(0,dtype=np.int64),np.empty(0,dtype=np.int8),
            representatives,representative_directions,current.copy(),current_ratio,0.0,
            rank,int(factor.k_aca),np.asarray(singular,dtype=float),coefficients,
            relative_error,
            "TSVD global binary model found no improving insertion",
            reachability_residual,reachability_ratio,reachability_relative,
            linearized_reachable,abe_murata_diagnostics=diagnostics)
    capture_ratio=(current_ratio-improvement_capture*
                   (current_ratio-best.predicted_max_band_ratio))
    compact=solve_element_generation_lp(
        current,target,band,move_delta,move_volumes,
        volume_budget=float(volume_budget),
        maximum_new_elements=maximum_changed_elements,
        whole_elements=True,predecessor_pairs=predecessor_pairs,
        candidate_objective_change=move_secondary_cost,
        predicted_ratio_cap=capture_ratio,
        candidate_volume_change=move_volume_changes,
        candidate_exclusion_groups=move_exclusion_groups,
        maximum_changed_volume=maximum_changed_volume)
    selected_columns=feasible_columns[np.asarray(compact.selected,dtype=bool)]
    selected=elements[selected_columns]
    return TSVDElementCandidateSelection(
        np.asarray(selected,dtype=np.int64),directions[selected_columns],
        np.asarray(representatives,dtype=np.int64),representative_directions,
        compact.predicted_response,
        compact.predicted_max_band_ratio,compact.added_volume,rank,
        int(factor.k_aca),np.asarray(singular,dtype=float),coefficients,relative_error,
        "all-candidate band-normalized TSVD plus binary LP",
        reachability_residual,reachability_ratio,reachability_relative,
        linearized_reachable,abe_murata_diagnostics=diagnostics)


def _interpolate_screened_response_correction(
        direct_delta,partial_delta,selected_rows,response_band):
    """Lift a few exact adjoint-row corrections to the candidate row space."""
    direct=np.asarray(direct_delta,dtype=float)
    partial=np.asarray(partial_delta,dtype=float)
    rows=np.asarray(selected_rows,dtype=np.int64).reshape(-1)
    band=np.asarray(response_band,dtype=float).reshape(-1)
    if (direct.ndim!=2 or partial.shape!=direct.shape or
            band.shape!=(direct.shape[0],) or np.any(band<=0.0) or
            rows.size==0 or np.unique(rows).size!=rows.size or
            np.any(rows<0) or np.any(rows>=direct.shape[0])):
        raise ValueError("screened response-correction arrays are incompatible")
    normalized=direct/band[:,None]
    skeleton=normalized[rows]
    interpolation=normalized@np.linalg.pinv(skeleton,rcond=1e-10)
    sampled_correction=(partial[rows]-direct[rows])/band[rows,None]
    return np.ascontiguousarray(
        direct+band[:,None]*(interpolation@sampled_correction),dtype=float)


def _configured_candidate_cluster_labels(charge_gram,dof_blocks,
                                         requested_clusters):
    """Map complete element blocks to the native charge H-matrix tree."""
    blocks=tuple(np.asarray(block,dtype=np.int32).reshape(-1)
                 for block in dof_blocks)
    if not blocks or any(block.size==0 for block in blocks):
        raise ValueError("candidate cluster blocks must be non-empty")
    requested=int(requested_clusters)
    if requested<1:
        raise ValueError("requested candidate clusters must be positive")
    packed=np.concatenate(blocks).astype(np.int32,copy=False)
    offsets=np.asarray(np.r_[0,np.cumsum([len(block) for block in blocks])],
                       dtype=np.int32)
    if hasattr(charge_gram,
               "configured_linear_material_candidate_clusters"):
        result=charge_gram.configured_linear_material_candidate_clusters(
            packed,offsets,requested)
        labels=np.asarray(result["labels"],dtype=np.int64).reshape(-1)
        count=int(result["n_cluster"])
        if (labels.shape!=(len(blocks),) or count<1 or
                np.any(labels<0) or np.any(labels>=count)):
            raise RuntimeError(
                "native candidate cluster membership is inconsistent")
        return labels,count
    # Compatibility for a downloaded binary that predates the native tree
    # surface.  This preserves correctness (the complete signed proposal still
    # receives a physical solve) but deliberately provides no spatial speedup.
    return np.zeros(len(blocks),dtype=np.int64),1


def _clustered_tsvd_candidate_front(*,candidate_elements,
        candidate_response_delta,response_band,cluster_labels,
        signed_coefficients,selected_elements=(),representative_elements=(),
        relative_tolerance=1e-3,front_limit=64):
    """Retain signed TSVD modes in every occupied native H-tree cluster.

    The global ACA--QR--TSVD proposal owns batch cardinality.  This helper
    builds only its conditional exact-physics fallback: selected/global QR
    columns are preserved, then each native spatial cluster contributes its
    first local QR skeleton column before higher local modes compete for the
    remaining front.  Hence a small exact front cannot collapse onto one pole
    station merely because neighbouring response columns are nearly parallel.
    """
    elements=np.asarray(candidate_elements,dtype=np.int64).reshape(-1)
    delta=np.asarray(candidate_response_delta,dtype=float)
    band=np.asarray(response_band,dtype=float).reshape(-1)
    labels=np.asarray(cluster_labels,dtype=np.int64).reshape(-1)
    coefficients=np.asarray(signed_coefficients,dtype=float).reshape(-1)
    selected=np.asarray(selected_elements,dtype=np.int64).reshape(-1)
    representatives=np.asarray(
        representative_elements,dtype=np.int64).reshape(-1)
    limit=int(front_limit);tolerance=float(relative_tolerance)
    if (delta.shape!=(band.size,elements.size) or labels.shape!=elements.shape
            or coefficients.shape!=elements.shape or np.any(band<=0.0)
            or np.any(labels<0) or np.unique(elements).size!=elements.size
            or limit<1 or not 0.0<tolerance<1.0
            or not np.all(np.isfinite(delta))):
        raise ValueError("clustered TSVD candidate arrays are incompatible")
    lookup={int(element):column for column,element in enumerate(elements)}
    mandatory=[]
    for value in np.r_[selected,representatives]:
        element=int(value)
        if element in lookup and element not in mandatory:
            mandatory.append(element)

    from .stream_function import aca_tsvd
    from scipy.linalg import qr
    first=[];extras=[]
    normalized=np.ascontiguousarray(delta/band[:,None])
    for cluster in np.unique(labels):
        columns=np.flatnonzero(labels==cluster)
        local=normalized[:,columns]
        if columns.size==1:
            ranked=[(int(columns[0]),1.0)]
        else:
            factor=aca_tsvd(
                local.shape[0],local.shape[1],
                lambda row,column:local[row,column],
                modes=min(local.shape),kmax=min(local.shape),
                aca_eps=tolerance,method="aca_qr_tsvd")
            singular=np.asarray(factor.S,dtype=float)
            V=np.asarray(factor.V,dtype=float)
            if singular.size==0 or not np.isfinite(singular[0]) or singular[0]<=0.0:
                order=np.argsort(-np.abs(coefficients[columns]),kind="stable")
                ranked=[(int(columns[index]),0.0) for index in order[:1]]
            else:
                rank=max(1,int(np.count_nonzero(
                    singular>=tolerance*singular[0])))
                local_v=V[:,:rank]
                _,_,pivots=qr(local_v.T,mode="economic",pivoting=True)
                leverage=np.linalg.norm(local_v,axis=1)
                ranked=[(int(columns[int(index)]),float(leverage[int(index)]))
                        for index in pivots[:rank]]
        if ranked:
            column,leverage=ranked[0]
            first.append((int(cluster),column,leverage,
                          abs(float(coefficients[column]))))
            for column,leverage in ranked[1:]:
                extras.append((int(cluster),column,leverage,
                               abs(float(coefficients[column]))))

    # Preserve spatial coverage before accepting second/third modes from any
    # one cluster.  If mandatory global modes already exceed the advertised
    # limit they remain intact: the limit bounds optional look-ahead only.
    front=list(mandatory)
    first.sort(key=lambda item:(-item[3],-item[2],item[0],
                                int(elements[item[1]])))
    extras.sort(key=lambda item:(-item[2],-item[3],item[0],
                                 int(elements[item[1]])))
    for _,column,_,_ in first+extras:
        element=int(elements[column])
        if element in front:
            continue
        if len(front)>=limit:
            break
        front.append(element)
    return np.asarray(front,dtype=np.int64)


def _adjoint_corrected_removal_material_response(*,charge_gram,inv_chi,
        dof_blocks,state,response_matrix,screen_context,response_band):
    """Approximate exact active-set deletion with an adjoint Schur contraction.

    For a removed block e, the exact response change is
    ``-C A^-1[:,e] S_e x_e`` with ``S_e=(A^-1[e,e])^-1``.  The native local
    block ``A_ee`` is the zero-extra-solve approximation to ``S_e``.  Exact
    state/adjoint H-matrix solves therefore give
    ``lambda_e.T A_ee x_e`` on sampled response rows; the usual candidate-row
    interpolation lifts that correction to the complete design response.
    """
    blocks=tuple(np.asarray(block,dtype=np.int32).reshape(-1)
                 for block in dof_blocks)
    values=np.asarray(state,dtype=float).reshape(-1)
    C=np.asarray(response_matrix,dtype=float)
    band=np.asarray(response_band,dtype=float).reshape(-1)
    if (not blocks or C.ndim!=2 or C.shape[1]!=values.size or
            band.shape!=(C.shape[0],)):
        raise ValueError("removal response contraction arrays are incompatible")
    direct=np.column_stack([C[:,block]@values[block] for block in blocks])
    context={} if screen_context is None else screen_context
    adjoints=context.get("adjoints")
    rows=np.asarray(context.get("adjoint_rows",()),dtype=np.int64).reshape(-1)
    if adjoints is None or rows.size==0:
        return np.ascontiguousarray(direct,dtype=float)
    adjoints=np.asarray(adjoints,dtype=float)
    if (adjoints.shape!=(values.size,rows.size) or
            np.any(rows<0) or np.any(rows>=C.shape[0])):
        raise RuntimeError("removal screen adjoints are inconsistent")
    packed=np.concatenate(blocks).astype(np.int32,copy=False)
    offsets=np.asarray(np.r_[0,np.cumsum([len(block) for block in blocks])],
                       dtype=np.int32)
    packed_matrices=np.asarray(
        charge_gram.configured_linear_material_element_blocks(
            float(inv_chi),packed,offsets),dtype=float).reshape(-1)
    partial=direct.copy();offset=0
    for column,block in enumerate(blocks):
        width=len(block);next_offset=offset+width*width
        local=packed_matrices[offset:next_offset].reshape(width,width)
        if local.size!=width*width:
            raise RuntimeError("native removal block packing is inconsistent")
        partial[rows,column]=adjoints[block].T@(local@values[block])
        offset=next_offset
    if offset!=packed_matrices.size:
        raise RuntimeError("native removal block packing has trailing values")
    if rows.size==C.shape[0]:
        return np.ascontiguousarray(partial,dtype=float)
    return _interpolate_screened_response_correction(
        direct,partial,rows,band)


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
    if candidates.size==0 or not front:
        current_ratio=float(np.max(np.abs((current-target)/band)))
        return CollaborativeElementBatchUpdate(
            np.empty(0,dtype=np.int64),current,current_ratio,0.0,0,
            "TSVD proposed no insertion-side exact front")
    if any(element not in lookup for element in front):
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
        solve_tolerance=1e-9, solve_max_iterations=5000, mass_riesz=True,
        cluster_coarse_size=0, cluster_deflation_size=0, recycle_size=0):
    """Solve one exact whole-element active iron set on the fixed superset mesh.

    The native Krylov routine uses periodically refreshed true residuals, but
    its historical result dictionary exposes only the iteration count.  Apply
    the same configured H-matrix operator once more here and reject a state
    that reached ``solve_max_iterations`` without satisfying the requested
    relative residual.  Candidate Schur models must never be accepted against
    an unconverged physical re-solve.
    """
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
    cluster_coarse_size=int(cluster_coarse_size)
    cluster_deflation_size=int(cluster_deflation_size)
    recycle_size=int(recycle_size)
    if min(cluster_coarse_size,cluster_deflation_size,recycle_size)<0:
        raise ValueError("cluster solve sizes must be nonnegative")
    result=charge_gram.solve_configured_linear_material_auto_prec_many(
        float(inv_chi),np.ascontiguousarray(active_rhs[None,:]),
        tol=float(solve_tolerance),maxit=int(solve_max_iterations),
        cluster_coarse_size=cluster_coarse_size,
        cluster_deflation_size=cluster_deflation_size,
        recycle_size=recycle_size,
        mass_riesz=bool(mass_riesz))
    state=np.asarray(result["m"],dtype=float)[0]
    applied=np.asarray(
        charge_gram.apply_configured_linear_material_operator(
            float(inv_chi),np.ascontiguousarray(state),
            respect_constraints=True),dtype=float).reshape(-1)
    rhs_norm=float(np.linalg.norm(active_rhs))
    residual_norm=float(np.linalg.norm(applied-active_rhs))
    relative_residual=(residual_norm/rhs_norm if rhs_norm>0.0 else
                       residual_norm)
    residual_limit=max(5.0*float(solve_tolerance),1.0e-12)
    if (not np.isfinite(relative_residual) or
            relative_residual>residual_limit):
        iterations=int(np.asarray(result["iters"]).reshape(-1)[0])
        raise RuntimeError(
            "HDiv-MMM active-set solve failed the true-residual gate: "
            f"relative residual {relative_residual:.6e} exceeds "
            f"{residual_limit:.6e} after {iterations} iterations")
    response=C@state
    if incident_response is not None:
        incident=np.asarray(incident_response,dtype=float).reshape(-1)
        if incident.shape!=response.shape: raise ValueError("incident_response shape mismatch")
        response=response+incident
    return state,response,int(result["iters"][0])


def _positive_minimax_source_scale_and_gradient(response,target,band):
    """Solve ``min_{alpha>0} max(abs(alpha*y-t)/b)`` analytically.

    The scalar convex objective is the upper envelope of two affine lines per
    observation.  An interior minimum is therefore an intersection of an
    active negative- and positive-slope line.  Enumerating those intersections
    is exact for the small source-calibration row set and also identifies the
    active pair whose implicit derivative gives ``d alpha / d response``.
    """
    values=np.asarray(response,dtype=float).reshape(-1)
    targets=np.asarray(target,dtype=float).reshape(-1)
    bands=np.asarray(band,dtype=float).reshape(-1)
    if (values.shape!=targets.shape or values.shape!=bands.shape or
            values.size==0 or np.any(~np.isfinite(values)) or
            np.any(~np.isfinite(targets)) or
            np.any(~np.isfinite(bands)) or np.any(bands<=0.0)):
        raise ValueError(
            "minimax source calibration arrays must be finite, non-empty, "
            "shape matched, and have positive bands")
    rows=np.r_[np.arange(values.size,dtype=np.int64),
               np.arange(values.size,dtype=np.int64)]
    signs=np.r_[np.ones(values.size),-np.ones(values.size)]
    slopes=signs*values[rows]/bands[rows]
    intercepts=-signs*targets[rows]/bands[rows]
    best_alpha=None;best_value=np.inf
    slope_scale=max(1.0,float(np.max(np.abs(slopes))))
    slope_tolerance=64.0*np.finfo(float).eps*slope_scale
    for left in range(len(slopes)):
        for right in range(left+1,len(slopes)):
            denominator=float(slopes[left]-slopes[right])
            if abs(denominator)<=slope_tolerance:
                continue
            alpha=float(
                (intercepts[right]-intercepts[left])/denominator)
            if not np.isfinite(alpha) or alpha<=0.0:
                continue
            value=float(np.max(slopes*alpha+intercepts))
            tie_tolerance=64.0*np.finfo(float).eps*max(
                1.0,abs(value),abs(best_value) if np.isfinite(best_value)
                else 1.0)
            if (value<best_value-tie_tolerance or
                    (abs(value-best_value)<=tie_tolerance and
                     (best_alpha is None or alpha<best_alpha))):
                best_alpha=alpha;best_value=value
    if best_alpha is None:
        raise RuntimeError(
            "source calibration has no finite positive minimax scale")
    boundary_value=float(np.max(intercepts))
    objective_tolerance=128.0*np.finfo(float).eps*max(
        1.0,abs(boundary_value),abs(best_value))
    if boundary_value<best_value-objective_tolerance:
        raise RuntimeError(
            "positive source calibration is minimized only at zero scale")
    line_values=slopes*best_alpha+intercepts
    active_tolerance=512.0*np.finfo(float).eps*max(
        1.0,abs(best_value))
    active=np.flatnonzero(
        np.abs(line_values-best_value)<=active_tolerance)
    if active.size<2:
        raise RuntimeError(
            "minimax source calibration did not expose an active pair")
    left=int(active[np.argmin(slopes[active])])
    right=int(active[np.argmax(slopes[active])])
    denominator=float(slopes[left]-slopes[right])
    if (slopes[left]>slope_tolerance or
            slopes[right]<-slope_tolerance or
            abs(denominator)<=slope_tolerance):
        raise RuntimeError(
            "minimax source calibration active slopes do not bracket zero")
    numerator=float(
        signs[left]*targets[rows[left]]/bands[rows[left]]-
        signs[right]*targets[rows[right]]/bands[rows[right]])
    alpha=numerator/denominator
    if (not np.isfinite(alpha) or alpha<=0.0 or
            not np.isclose(alpha,best_alpha,rtol=2e-11,
                           atol=2e-13*max(1.0,abs(best_alpha)))):
        raise RuntimeError(
            "minimax source calibration active-pair reconstruction failed")
    gradient=np.zeros(values.size,dtype=float)
    np.add.at(gradient,rows[left],
              -alpha*signs[left]/(bands[rows[left]]*denominator))
    np.add.at(gradient,rows[right],
              alpha*signs[right]/(bands[rows[right]]*denominator))
    return float(alpha),gradient


def grow_hdiv_mmm_by_superposition(*, charge_gram, fes, inv_chi, rhs,
        response_matrix, active_elements, element_volumes,
        response_target, response_band, volume_max,
        incident_response=None, maximum_batch_elements=None,
        max_iterations=30, ratio_tolerance=1e-8,
        solve_tolerance=1e-9, solve_max_iterations=5000,
        candidate_batch_size=64, mass_riesz=True,
        cluster_coarse_size=0, cluster_deflation_size=0,
        recycle_size=0,
        fixed_inactive_elements=None,
        fixed_active_elements=None,
        removal_coupling_groups=None,
        predecessor_elements=None,
        active_set_validator=None,
        source_calibration_rows=None,
        source_calibration_target=None,
        source_calibration_band=None,
        source_calibration_norm="mean",
        response_transform=None,
        response_transform_jacobian=None,
        exact_response_validator=None,
        include_predecessor_descendants=False,
        exact_candidate_limit=64,
        batch_improvement_capture=0.9,
        tsvd_relative_tolerance=1e-3,
        proposal_adjoint_count=4,
        proposal_solve_tolerance=1e-5,
        minimum_model_agreement=0.1,
        proposal_trust_region_trials=4,
        initial_material_move_fraction=None,
        maximum_material_move_fraction=0.25,
        graph_interface_weight=0.0,
        graph_front_maximum_components=2,
        graph_front_beam_width=64,
        graph_front_proposal_limit=8,
        graph_front_response_novelty_weight=0.0,
        exact_beam_width=0,
        exact_beam_depth=0,
        exact_beam_barrier_fraction=0.25,
        removal_cluster_count=0,
        iteration_callback=None) -> HDivMMMGenerationResult:
    """Grow connected whole iron elements by Schur superposition and 0-1 LP.

    Each iteration cheaply evaluates every face-adjacent inactive addition and
    every topology-safe active boundary removal,
    factors the band-normalized global response by the canonical native
    ACA+--thin-QR--TSVD kernel, and solves a binary minimum-volume proposal on
    the retained response subspace.  The complete signed proposal is checked
    first by one exact active-set HDiv-MMM re-solve.  Only a rejected/empty
    proposal enters the bounded adjoint block-Schur fallback, where signed
    singular-mode representatives preserve cooperative insertions.  A batch is
    committed only after an exact active-set HDiv-MMM re-solve.  The
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
    recalibrated after every candidate insertion.  ``source_calibration_norm``
    may retain the legacy ``"mean"`` match, or select ``"linf"`` to solve the
    exact positive one-variable Chebyshev problem in the supplied
    ``source_calibration_band``.  The active pair of affine residuals supplies
    its piecewise-analytic scale gradient.  This eliminates coil current while
    the 0-1 LP remains responsible only for iron geometry.
    ``response_transform`` may map the raw linear field response to nonlinear
    design metrics.  It is evaluated on every exact one-element Schur response
    and exact accepted batch, so the LP sees finite metric changes rather than
    a finite-difference derivative.  ``response_transform_jacobian`` may return
    its analytic Jacobian at a calibrated raw response.  The driver then
    contracts the raw observation rows to the design-metric dimension *before*
    the active-system adjoint solves.  Source-recalibration differentiation is
    included exactly, so an 80-row field profile can require only four optics
    adjoints while every committed material move still receives a full raw
    field re-solve.  The production ``proposal_adjoint_count`` keeps only the
    candidate-response row space is first exposed by the unbiased direct
    material-superposition screen.  Pivoted QR then chooses at most
    ``proposal_adjoint_count`` representative response rows (four by default),
    including the current worst minimax row.  Their exact active-relaxation
    adjoint corrections are interpolated across the complete candidate-response
    row space.  This avoids both the biased one-worst-row screen and an
    unbounded solve of every response adjoint.  A
    rejected or poorly predicted proposal shrinks a discrete volume trust
    region and resolves the binary LP over the same complete candidate set.
    Thus production does not fall back to an unbounded all-row adjoint solve.
    Deletion columns use an adjoint-contracted local Schur action instead of
    raw element observation alone.  If a signed proposal is rejected, actual
    charge H-matrix cluster membership and per-cluster ACA--QR--TSVD skeletons
    define the bounded exact-physics removal front.  ``removal_cluster_count``
    may prescribe the number of native tree nodes; zero chooses it from the
    retained global response rank.
    Tiny fronts retain the full adjoint model automatically.  Proposal adjoints
    use the separate, deliberately looser ``proposal_solve_tolerance`` because
    the accepted active-set solve still uses ``solve_tolerance`` and is the sole
    physical acceptance gate.  ``minimum_model_agreement`` is the usual ratio
    of actual to predicted minimax reduction; a smaller ratio triggers at most
    ``proposal_trust_region_trials`` geometrically smaller whole-element
    proposals, never gray material or a fixed element-count guess.

    ``removal_coupling_groups`` assigns a nonnegative group id to elements
    that must be deleted together (``-1`` leaves an element independent).
    This is useful for a sheet/pole Lego whose through-thickness cells must
    share one terminal edge.  TSVD columns, volume accounting, topology gates,
    and committed full solves all use the complete group; no member can leave
    a partial-thickness notch.

    ``active_set_validator`` may impose an application-specific binary geometry
    rule after the generic iron/air connectivity gate.  It receives the full
    boolean active mask and must return truth for an admissible shape.  This is
    intended for exact rules such as a one-cell Lipschitz bound on neighbouring
    pole-end columns; it is never a density penalty or finite-difference model.

    ``exact_response_validator`` may impose an application-specific physics
    guard after source calibration and a complete active-set solve.  It receives
    copies of the calibrated raw response and transformed design response and
    must return one boolean.  Guard-rejected states are neither accepted nor
    retained as nonmonotone beam parents.  Proposal linearization remains
    unchanged; the complete solve is the single source of truth for the guard.

    ``initial_material_move_fraction`` enables a TOBS/SAIP-style trust region
    on the *total physical volume flipped* by one signed add/remove proposal.
    This is deliberately independent of the net iron-volume constraint: a
    large addition and a large deletion may have zero net volume but remain a
    remote, poorly predicted binary state.  Accepted physical re-solves expand,
    hold, or shrink the radius from their actual/predicted reduction ratio.

    ``graph_interface_weight`` adds a physical facet-area regularizer to the
    lexicographic binary master after the best response-band reduction has been
    fixed.  Analytic signed Schur response columns are never filtered or
    averaged across iterations.

    The same face graph also supplies a connected-front challenge when
    ``graph_front_proposal_limit`` is positive.  ACA/QR/TSVD representatives
    seed each component; every subsequent move must touch that component.
    Its maximum cardinality is derived from the current physical-volume trust
    radius (or, without one, from the retained response rank), never from a
    fixed ``try N elements`` batch.
    ``graph_front_response_novelty_weight`` greedily retains response-space
    independent proposals after the best predicted proposal.  Zero preserves
    pure minimax ranking; positive values trade a bounded amount of predicted
    quality for band-normalized subspace novelty before any exact solve.

    A positive ``exact_beam_width`` and ``exact_beam_depth`` enable a bounded
    nonmonotone Lego search.  Rejected graph-front states that have already
    passed a complete HDiv-MMM solve may be retained inside
    ``exact_beam_barrier_fraction`` of the best incumbent.  They are
    relinearized as ordinary binary states, but are never returned unless a
    descendant improves the incumbent.  Thus the beam can cross a shallow
    discrete barrier without accepting a worse final design or using a design
    finite difference.
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
    minimum_model_agreement=float(minimum_model_agreement)
    trust_region_trials=int(proposal_trust_region_trials)
    if (not np.isfinite(minimum_model_agreement) or
            not 0.0<=minimum_model_agreement<=1.0):
        raise ValueError("minimum_model_agreement must lie in [0, 1]")
    if trust_region_trials<1:
        raise ValueError("proposal_trust_region_trials must be positive")
    removal_cluster_count=int(removal_cluster_count)
    if removal_cluster_count<0:
        raise ValueError("removal_cluster_count must be nonnegative")
    if maximum_cap is not None and int(exact_candidate_limit)<maximum_cap:
        raise ValueError("exact_candidate_limit must cover maximum_batch_elements")
    if (not 0.0<float(batch_improvement_capture)<=1.0 or
            not 0.0<float(tsvd_relative_tolerance)<1.0):
        raise ValueError("TSVD selection parameters are invalid")
    move_fraction=(None if initial_material_move_fraction is None else
                   float(initial_material_move_fraction))
    maximum_move_fraction=float(maximum_material_move_fraction)
    if ((move_fraction is not None and
         (not np.isfinite(move_fraction) or not 0.0<move_fraction<=1.0)) or
            not np.isfinite(maximum_move_fraction) or
            not 0.0<maximum_move_fraction<=1.0 or
            (move_fraction is not None and
             maximum_move_fraction<move_fraction)):
        raise ValueError(
            "material move fractions must satisfy 0 < initial <= maximum <= 1")
    graph_interface_weight=float(graph_interface_weight)
    graph_front_maximum_components=int(graph_front_maximum_components)
    graph_front_beam_width=int(graph_front_beam_width)
    graph_front_proposal_limit=int(graph_front_proposal_limit)
    graph_front_response_novelty_weight=float(
        graph_front_response_novelty_weight)
    exact_beam_width=int(exact_beam_width)
    exact_beam_depth=int(exact_beam_depth)
    exact_beam_barrier_fraction=float(exact_beam_barrier_fraction)
    if (not np.isfinite(graph_interface_weight) or
            graph_interface_weight<0.0 or
            graph_front_maximum_components<1 or
            graph_front_beam_width<1 or graph_front_proposal_limit<0 or
            not 0.0<=graph_front_response_novelty_weight<=1.0 or
            exact_beam_width<0 or exact_beam_depth<0 or
            not np.isfinite(exact_beam_barrier_fraction) or
            exact_beam_barrier_fraction<0.0):
        raise ValueError("graph-front parameters are invalid")
    if exact_beam_width==0 or exact_beam_depth==0:
        exact_beam_width=0;exact_beam_depth=0
    calibration=(None if source_calibration_rows is None else
                 np.asarray(source_calibration_rows,dtype=np.int64).reshape(-1))
    calibration_norm=str(source_calibration_norm).lower()
    if calibration_norm not in ("mean","linf"):
        raise ValueError(
            "source_calibration_norm must be 'mean' or 'linf'")
    if response_transform_jacobian is not None and response_transform is None:
        raise ValueError(
            "response_transform_jacobian requires response_transform")
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
        if np.any(~np.isfinite(calibration_values)):
            raise ValueError("source calibration target must be finite")
        calibration_target=float(np.mean(calibration_values))
        if (calibration_norm=="mean" and
                (not np.isfinite(calibration_target) or
                 calibration_target==0.0)):
            raise ValueError(
                "mean source calibration target must have a finite, "
                "nonzero mean")
        if source_calibration_band is None:
            calibration_band=np.ones(calibration.size,dtype=float)
        else:
            calibration_band=np.asarray(
                source_calibration_band,dtype=float).reshape(-1)
            if (calibration_band.shape!=calibration.shape or
                    np.any(~np.isfinite(calibration_band)) or
                    np.any(calibration_band<=0.0)):
                raise ValueError(
                    "source_calibration_band must be positive and match "
                    "source_calibration_rows")
        if calibration_norm=="linf" and source_calibration_band is None:
            raise ValueError(
                "linf source calibration requires source_calibration_band")

    def source_scale_and_gradient(base_response):
        values=np.asarray(base_response,dtype=float).reshape(-1)
        if calibration is None:
            return 1.0,np.zeros(values.size,dtype=float)
        if np.any(calibration>=values.size):
            raise ValueError(
                "source_calibration_rows index outside the raw response")
        selected=values[calibration]
        gradient=np.zeros(values.size,dtype=float)
        if calibration_norm=="mean":
            denominator=float(np.mean(selected))
            if not np.isfinite(denominator) or denominator==0.0:
                raise RuntimeError(
                    "source calibration response mean is zero or invalid")
            scale=calibration_target/denominator
            local_gradient=np.full(
                calibration.size,
                -scale/(float(calibration.size)*denominator),dtype=float)
        else:
            scale,local_gradient=_positive_minimax_source_scale_and_gradient(
                selected,calibration_values,calibration_band)
        if not np.isfinite(scale) or scale<=0.0:
            raise RuntimeError(
                "source calibration requires a positive finite source scale")
        np.add.at(gradient,calibration,local_gradient)
        return float(scale),gradient

    def transform_response(raw_response):
        raw=np.asarray(raw_response,dtype=float).reshape(-1)
        transformed=(raw if response_transform is None else
                     np.asarray(response_transform(raw),dtype=float).reshape(-1))
        if transformed.shape!=target.shape or not np.all(np.isfinite(transformed)):
            raise ValueError("response_transform must return one finite design-response vector")
        return transformed

    def transform_jacobian(calibrated_response):
        if response_transform_jacobian is None:
            return None
        raw=np.asarray(calibrated_response,dtype=float).reshape(-1)
        jacobian=np.asarray(response_transform_jacobian(raw),dtype=float)
        if (jacobian.shape!=(target.size,raw.size) or
                not np.all(np.isfinite(jacobian))):
            raise ValueError(
                "response_transform_jacobian must return a finite matrix "
                "with shape (n_design_response,n_raw_response)")
        return jacobian

    def exact_response_is_valid(raw_response,design_response):
        if exact_response_validator is None:
            return True
        verdict=exact_response_validator(
            np.asarray(raw_response,dtype=float).copy(),
            np.asarray(design_response,dtype=float).copy())
        result=np.asarray(verdict)
        if result.ndim!=0:
            raise ValueError(
                "exact_response_validator must return one boolean")
        return bool(result)

    def calibrate_source(base_state,base_response):
        state_value=np.asarray(base_state,dtype=float)
        response_value=np.asarray(base_response,dtype=float)
        if calibration is None:
            return state_value,response_value,1.0
        scale,_=source_scale_and_gradient(response_value)
        return state_value*scale,response_value*scale,float(scale)
    fixed=(np.zeros_like(active) if fixed_inactive_elements is None else
           np.asarray(fixed_inactive_elements,dtype=bool).reshape(-1))
    if fixed.shape!=active.shape or np.any(active&fixed):
        raise ValueError("fixed inactive elements must match the mesh and cannot be active")
    fixed_active=(np.zeros_like(active) if fixed_active_elements is None else
                  np.asarray(fixed_active_elements,dtype=bool).reshape(-1))
    if fixed_active.shape!=active.shape or np.any(fixed_active&~active):
        raise ValueError("fixed active elements must match and belong to the initial iron")
    removal_groups=(np.full(active.shape,-1,dtype=np.int64)
                    if removal_coupling_groups is None else
                    np.asarray(removal_coupling_groups,dtype=np.int64).reshape(-1))
    if (removal_groups.shape!=active.shape or np.any(removal_groups<-1)):
        raise ValueError(
            "removal_coupling_groups must contain group ids or -1")
    movable=~(fixed|fixed_active)
    if not np.any(movable):
        raise ValueError("HDiv-MMM generation needs movable material elements")
    movable_volume=float(np.sum(volumes[movable]))
    minimum_material_move=float(np.min(volumes[movable]))
    maximum_material_move=max(
        minimum_material_move,maximum_move_fraction*movable_volume)
    material_trust_volume=(None if move_fraction is None else min(
        maximum_material_move,max(
            minimum_material_move,move_fraction*movable_volume)))
    graph_enabled=bool(
        graph_front_proposal_limit>0 or graph_interface_weight>0.0)
    if graph_enabled:
        from ._topopt_graph import ngsolve_facet_measure_graph
        (element_graph,element_exterior,
         element_interface_weights)=ngsolve_facet_measure_graph(fes.mesh)
    else:
        element_graph=None;element_exterior=None
        element_interface_weights=None
    graph_front_data=None
    graph_front_budget=None

    def update_material_trust(before,changed,agreement):
        if before is None:
            return None
        before=float(before);changed=float(changed);agreement=float(agreement)
        if agreement<0.25:
            return max(minimum_material_move,0.5*before)
        if agreement>0.75 and changed>=0.8*before:
            return min(maximum_material_move,max(
                before+minimum_material_move,1.5*before))
        return before

    def valid_active_set(mask):
        values=np.asarray(mask,dtype=bool).reshape(-1)
        if values.shape!=active.shape:
            return False
        if not np.any(values):
            return False
        if predecessors is not None:
            dependent=np.flatnonzero(values&(predecessors>=0))
            if (dependent.size and
                    np.any(~values[predecessors[dependent]])):
                return False
        if not ngsolve_growth_topology(fes.mesh,values).valid:
            return False
        if active_set_validator is None:
            return True
        verdict=active_set_validator(values.copy())
        result=np.asarray(verdict)
        if result.ndim!=0:
            raise ValueError("active_set_validator must return one boolean")
        return bool(result)

    def grouped_removal_candidates(raw_candidates,current_active):
        raw=np.asarray(raw_candidates,dtype=np.int64).reshape(-1)
        raw_set=set(int(value) for value in raw)
        members={}
        for element in raw:
            group=int(removal_groups[int(element)])
            if group<0:
                block=np.asarray([int(element)],dtype=np.int64)
            else:
                block=np.flatnonzero(removal_groups==group).astype(np.int64)
            representative=int(np.min(block))
            if representative in members:
                continue
            if (not np.all(current_active[block]) or np.any(fixed_active[block])
                    or not any(int(value) in raw_set for value in block)):
                continue
            trial=current_active.copy();trial[block]=False
            if predecessors is not None:
                # A coupled pack may deliberately contain an entire
                # predecessor chain (for example every movable depth cell of
                # one pole-end station).  Its inner members are not legal
                # singleton removals, so requiring every member in ``raw``
                # incorrectly deletes the whole removal design space.  Accept
                # the chain when at least one exposed member is raw and no
                # surviving active cell depends on a removed predecessor.
                surviving=np.flatnonzero(trial)
                if np.any(np.isin(predecessors[surviving],block)):
                    continue
            if not valid_active_set(trial):
                continue
            members[representative]=block
        representatives=np.asarray(sorted(members),dtype=np.int64)
        return representatives,members
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
    if active_set_validator is not None and not valid_active_set(active):
        raise ValueError("initial iron seed violates active_set_validator")
    state,response,solve_iterations=solve_hdiv_mmm_active_elements(
        charge_gram=charge_gram,fes=fes,inv_chi=inv_chi,rhs=rhs,
        response_matrix=response_matrix,active_elements=active,
        incident_response=incident_response,solve_tolerance=solve_tolerance,
        solve_max_iterations=solve_max_iterations,mass_riesz=mass_riesz,
        cluster_coarse_size=cluster_coarse_size,
        cluster_deflation_size=cluster_deflation_size,
        recycle_size=recycle_size)
    state,response,source_scale=calibrate_source(state,response)
    element_blocks=ngsolve_discontinuous_element_dof_blocks(fes)
    ratio=lambda values:float(np.max(np.abs((np.asarray(values)-target)/band)))
    objective_response=transform_response(response)
    if not exact_response_is_valid(response,objective_response):
        raise ValueError("initial exact response violates exact_response_validator")
    current_ratio=ratio(objective_response);history=[]
    exact_search_trace=[]
    graph_front_diagnostics=[]
    converged=current_ratio<=1.0+ratio_tolerance
    stop_reason="target_met" if converged else "max_iterations"
    search_depth=0
    search_path=()
    pending_exact_states=[]
    visited_exact_states={np.packbits(active).tobytes()}

    def exact_snapshot(*,active_value,state_value,response_value,
                       objective_value,ratio_value,scale_value,
                       iterations_value,depth_value,path_value):
        return _HDivMMMExactBeamState(
            np.asarray(active_value,dtype=bool).copy(),
            np.asarray(state_value,dtype=float).copy(),
            np.asarray(response_value,dtype=float).copy(),
            np.asarray(objective_value,dtype=float).copy(),float(ratio_value),
            float(scale_value),int(iterations_value),int(depth_value),
            tuple(path_value))

    incumbent=exact_snapshot(
        active_value=active,state_value=state,response_value=response,
        objective_value=objective_response,ratio_value=current_ratio,
        scale_value=source_scale,iterations_value=solve_iterations,
        depth_value=0,path_value=())

    def record_exact_trial(snapshot):
        exact_search_trace.append(HDivMMMExactSearchTrial(
            depth=int(snapshot.depth),
            parent_max_band_ratio=float(current_ratio),
            incumbent_max_band_ratio=float(incumbent.ratio),
            max_band_ratio=float(snapshot.ratio),
            added_elements=np.flatnonzero(
                snapshot.active & ~incumbent.active),
            removed_elements=np.flatnonzero(
                incumbent.active & ~snapshot.active),
            solve_iterations=int(snapshot.solve_iterations),
            path=snapshot.path))
        return snapshot

    def next_nonmonotone_state(trials):
        """Queue exact, topology-valid barrier states and pop one beam node."""
        nonlocal pending_exact_states
        if not exact_beam_width:
            return None
        if search_depth<exact_beam_depth:
            barrier=incumbent.ratio*(1.0+exact_beam_barrier_fraction)
            for trial in trials:
                key=np.packbits(trial.active).tobytes()
                if (trial.depth<=exact_beam_depth and
                        trial.ratio<=barrier+ratio_tolerance and
                        key not in visited_exact_states):
                    visited_exact_states.add(key)
                    pending_exact_states.append(trial)

        # Preserve the best exact ratio at every depth, then use normalized
        # response-space novelty.  This prevents a beam from spending all its
        # width on near-identical Lego states from one terminal station.
        pruned=[]
        for depth in sorted({value.depth for value in pending_exact_states}):
            pool=sorted((value for value in pending_exact_states
                         if value.depth==depth),
                        key=lambda value:(value.ratio,value.path))
            if not pool:
                continue
            chosen=[pool.pop(0)]
            while pool and len(chosen)<exact_beam_width:
                def diversity(value):
                    novelty=min(float(np.linalg.norm(
                        (value.objective_response-other.objective_response)/band))
                        for other in chosen)
                    return (novelty/max(value.ratio,1.0),-value.ratio)
                index=max(range(len(pool)),key=lambda value:diversity(pool[value]))
                chosen.append(pool.pop(index))
            pruned.extend(chosen)
        pending_exact_states=pruned
        if not pending_exact_states:
            return None
        pending_exact_states.sort(
            key=lambda value:(value.depth,value.ratio,value.path))
        return pending_exact_states.pop(0)

    maximum_improvements=max(0,int(max_iterations))
    expansion_limit=maximum_improvements*(
        1+exact_beam_width*exact_beam_depth)
    expansion=0
    while expansion<expansion_limit and len(history)<maximum_improvements:
        iteration=len(history);expansion+=1
        if converged:
            stop_reason="target_met";break
        if exact_beam_width and search_depth>=exact_beam_depth:
            # Depth is a strict solve budget: a state at the maximum depth may
            # be retained as a diagnostic barrier node, but it must never be
            # expanded into an unrequested depth+1 active solve.  Prefer any
            # shallower queued branch before declaring the beam exhausted.
            branch=next_nonmonotone_state(())
            if branch is None:
                stop_reason="exact_nonmonotone_beam_exhausted"
                break
            active=branch.active;state=branch.state
            response=branch.response
            objective_response=branch.objective_response
            current_ratio=branch.ratio
            source_scale=branch.source_scale
            solve_iterations=branch.solve_iterations
            search_depth=branch.depth;search_path=branch.path
            stop_reason="exact_nonmonotone_beam_in_progress"
            continue
        material_trust_before=material_trust_volume
        graph_front_data=None
        exact_trial_states=[]
        remaining=float(volume_max)-float(volumes@active)
        candidates=ngsolve_boundary_growth_candidates(
            fes.mesh,active,fixed_inactive_elements=fixed,
            predecessor_elements=predecessors,
            include_predecessor_descendants=include_predecessor_descendants)
        raw_removal_candidates=ngsolve_boundary_removal_candidates(
            fes.mesh,active,fixed_active_elements=fixed_active,
            predecessor_elements=predecessors)
        removal_candidates,removal_members=grouped_removal_candidates(
            raw_removal_candidates,active)
        if candidates.size==0 and removal_candidates.size==0:
            stop_reason="no_growth_candidates";break
        if (removal_candidates.size==0 and
                remaining<min(volumes[candidates],default=np.inf)-1e-14):
            stop_reason="volume_budget_exhausted";break
        if candidates.size==0:
            material=[]
            for element in removal_candidates:
                block=np.concatenate([
                    element_blocks[int(member)]
                    for member in removal_members[int(element)]])
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
                candidate_volumes=np.asarray([
                    np.sum(volumes[removal_members[int(element)]])
                    for element in removal_candidates]),
                volume_budget=max(0.0,remaining),
                candidate_material_active=np.ones(len(removal_candidates),dtype=bool),
                relative_tolerance=float(tsvd_relative_tolerance),
                improvement_capture=float(batch_improvement_capture),
                ratio_tolerance=ratio_tolerance,
                maximum_changed_volume=material_trust_volume,
                maximum_changed_elements=maximum_cap)
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
                expanded=np.unique(np.concatenate([
                    removal_members[int(element)] for element in bundle]))
                trial_active=active.copy();trial_active[expanded]=False
                if not valid_active_set(trial_active):
                    continue
                batch_trials+=1
                trial_state,trial_response,trial_iterations=\
                    solve_hdiv_mmm_active_elements(
                        charge_gram=charge_gram,fes=fes,inv_chi=inv_chi,rhs=rhs,
                        response_matrix=response_matrix,active_elements=trial_active,
                        incident_response=incident_response,
                        solve_tolerance=solve_tolerance,
                        solve_max_iterations=solve_max_iterations,
                        mass_riesz=mass_riesz,
                        cluster_coarse_size=cluster_coarse_size,
                        cluster_deflation_size=cluster_deflation_size,
                        recycle_size=recycle_size)
                try:
                    trial_state,trial_response,trial_scale=calibrate_source(
                        trial_state,trial_response)
                    trial_objective=transform_response(trial_response)
                    if not exact_response_is_valid(
                            trial_response,trial_objective):
                        continue
                except (RuntimeError,ValueError):
                    continue
                trial_ratio=ratio(trial_objective)
                exact_trial_states.append(record_exact_trial(exact_snapshot(
                    active_value=trial_active,state_value=trial_state,
                    response_value=trial_response,
                    objective_value=trial_objective,ratio_value=trial_ratio,
                    scale_value=trial_scale,
                    iterations_value=trial_iterations,
                    depth_value=search_depth+1,
                    path_value=search_path+(tuple(
                        (int(element),-1) for element in bundle),))))
                if (trial_ratio<current_ratio-ratio_tolerance and
                        (best is None or trial_ratio<best[0])):
                    best=(trial_ratio,bundle,expanded,trial_active,trial_state,
                          trial_response,trial_iterations,trial_scale,
                          trial_objective)
                    # The globally compressed TSVD proposal is the first
                    # attempt.  Exact solves are an acceptance gate, not an
                    # exhaustive local optimizer: once an ordered proposal
                    # improves, commit it and let the next outer iteration
                    # relinearize every candidate.  Alternatives are paid for
                    # only when all earlier global/nested proposals fail.
                    break
            if best is None:
                branch=next_nonmonotone_state(exact_trial_states)
                if branch is not None:
                    active=branch.active;state=branch.state
                    response=branch.response
                    objective_response=branch.objective_response
                    current_ratio=branch.ratio
                    source_scale=branch.source_scale
                    solve_iterations=branch.solve_iterations
                    search_depth=branch.depth;search_path=branch.path
                    stop_reason="exact_nonmonotone_beam_in_progress"
                    continue
                stop_reason=("exact_nonmonotone_beam_exhausted" if
                    exact_beam_width else
                    "conditional_exact_rejected_removal_front")
                break
            (actual,selected_remove,removed_physical,trial_active,trial_state,trial_response,
             trial_iterations,trial_scale,trial_objective)=best
            if (search_depth>0 and
                    actual>=incumbent.ratio-ratio_tolerance):
                branch=next_nonmonotone_state(exact_trial_states)
                if branch is not None:
                    active=branch.active;state=branch.state
                    response=branch.response
                    objective_response=branch.objective_response
                    current_ratio=branch.ratio
                    source_scale=branch.source_scale
                    solve_iterations=branch.solve_iterations
                    search_depth=branch.depth;search_path=branch.path
                    stop_reason="exact_nonmonotone_beam_in_progress"
                    continue
                stop_reason="exact_nonmonotone_beam_exhausted"
                break
            accepted_search_depth=search_depth
            selected_columns=np.asarray([
                removal_lookup[int(element)] for element in selected_remove],
                dtype=np.int64)
            predicted=ratio(objective_response-
                            np.sum(material_matrix[:,selected_columns],axis=1))
            agreement=((current_ratio-actual)/(current_ratio-predicted)
                       if predicted<current_ratio-ratio_tolerance else 0.0)
            accepted_added=np.flatnonzero(
                trial_active & ~incumbent.active)
            accepted_removed=np.flatnonzero(
                incumbent.active & ~trial_active)
            changed_volume=float(
                np.sum(volumes[accepted_added])+np.sum(
                    volumes[accepted_removed]))
            material_trust_volume=update_material_trust(
                material_trust_before,changed_volume,agreement)
            active=trial_active;state=trial_state;response=trial_response
            solve_iterations=trial_iterations;source_scale=trial_scale
            objective_response=trial_objective;current_ratio=actual
            row=HDivMMMGenerationIteration(
                iteration,len(removal_candidates),accepted_added,
                predicted,current_ratio,int(np.count_nonzero(active)),
                float(volumes@active),batch_trials,solve_iterations,source_scale,
                predicted,"signed-magnetization-aca-qr-tsvd-conditional-exact",
                batch_trials,len(representative),{},len(ordered),
                len(selected_remove),float(agreement),
                int(removal_tsvd.numerical_rank),int(removal_tsvd.aca_rank),
                float(removal_tsvd.relative_truncation_error),
                int(len(removal_tsvd.selected_elements)),accepted_removed,0,
                int(len(removal_candidates)),
                linearized_reachability_residual=(None if
                    removal_tsvd.linearized_reachability_residual is None else
                    np.asarray(removal_tsvd.linearized_reachability_residual,
                               dtype=float)),
                linearized_reachability_max_band_ratio=float(
                    removal_tsvd.linearized_reachability_max_band_ratio),
                linearized_reachability_relative_residual=float(
                    removal_tsvd.linearized_reachability_relative_residual),
                linearized_reachable=bool(removal_tsvd.linearized_reachable),
                material_trust_volume_before=material_trust_before,
                material_trust_volume_after=material_trust_volume,
                material_changed_volume=changed_volume,
                nonmonotone_search_depth=int(accepted_search_depth),
                abe_murata_diagnostics=(
                    removal_tsvd.abe_murata_diagnostics))
            history.append(row)
            if iteration_callback is not None: iteration_callback(row)
            incumbent=exact_snapshot(
                active_value=active,state_value=state,response_value=response,
                objective_value=objective_response,ratio_value=current_ratio,
                scale_value=source_scale,iterations_value=solve_iterations,
                depth_value=0,path_value=())
            search_depth=0;search_path=();pending_exact_states=[]
            visited_exact_states={np.packbits(active).tobytes()}
            converged=current_ratio<=1.0+ratio_tolerance
            if converged: stop_reason="target_met"
            continue

        # Contract field observations to the actual design metrics before the
        # state/adjoint batch.  The complete derivative of the selected mean or
        # Chebyshev source elimination is composed with the analytic optics
        # Jacobian.  The projection is only a proposal model.  Accepted active
        # sets are still solved and scored with the full raw response above.
        objective_projection=None
        linear_response_matrix=np.asarray(response_matrix,dtype=float)
        linear_incident_response=incident_response
        if response_transform_jacobian is not None:
            if not np.isfinite(source_scale) or source_scale<=0.0:
                raise RuntimeError(
                    "response-row contraction requires a positive source scale")
            raw_base=np.asarray(response,dtype=float)/float(source_scale)
            calibrated_jacobian=transform_jacobian(response)
            calibration_jacobian=float(source_scale)*np.eye(raw_base.size)
            if calibration is not None:
                projection_scale,scale_gradient=source_scale_and_gradient(
                    raw_base)
                if not np.isclose(
                        projection_scale,float(source_scale),rtol=2e-12,
                        atol=2e-14*max(1.0,abs(float(source_scale)))):
                    raise RuntimeError(
                        "source-calibrated response projection scale drifted")
                calibration_jacobian+=np.outer(
                    raw_base,scale_gradient)
            objective_projection=np.ascontiguousarray(
                calibrated_jacobian@calibration_jacobian,dtype=float)
            linear_response_matrix=np.ascontiguousarray(
                objective_projection@linear_response_matrix,dtype=float)
            if incident_response is not None:
                linear_incident_response=np.ascontiguousarray(
                    objective_projection@np.asarray(
                        incident_response,dtype=float).reshape(-1),dtype=float)

        tsvd_proposal=None
        tsvd_material_data=None
        removal_material_response={}
        removal_cluster_front=np.empty(0,dtype=np.int64)
        fallback_addition_elements=np.empty(0,dtype=np.int64)
        screen_context={}
        working_adjoint_count=int(proposal_adjoint_count)
        if working_adjoint_count<0:
            raise ValueError("proposal_adjoint_count must be nonnegative")
        proposal_tolerance=max(
            float(solve_tolerance),float(proposal_solve_tolerance))
        if not np.isfinite(proposal_tolerance) or proposal_tolerance<=0.0:
            raise ValueError("proposal_solve_tolerance must be positive and finite")
        if len(candidates)<=8:
            proposal_tolerance=float(solve_tolerance)
            proposal_adjoint_rows=np.arange(target.size,dtype=np.int64)
        else:
            proposal_adjoint_rows=np.empty(0,dtype=np.int64)
        use_proposal_adjoint=proposal_adjoint_rows.size>0
        def select_exact_candidates(elements,approximate_delta,
                                    approximate_state,approximate_response):
            nonlocal tsvd_proposal,tsvd_material_data
            nonlocal fallback_addition_elements,removal_cluster_front
            nonlocal graph_front_data
            elements=np.asarray(elements,dtype=np.int64).reshape(-1)
            if objective_projection is None:
                _,base_calibrated,_=calibrate_source(
                    approximate_state,approximate_response)
                base_objective=transform_response(base_calibrated)
            else:
                base_objective=np.asarray(objective_response,dtype=float)
            effective=[];valid=[]
            for column in range(elements.size):
                delta=np.asarray(approximate_delta[:,column],dtype=float)
                if not np.all(np.isfinite(delta)):
                    continue
                if objective_projection is None:
                    try:
                        _,inserted,_=calibrate_source(
                            approximate_state,approximate_response+delta)
                        objective=transform_response(inserted)
                    except (RuntimeError,ValueError):
                        continue
                    delta=objective-base_objective
                valid.append(column);effective.append(delta)
            if not valid:
                raise RuntimeError("all candidate insertion responses are invalid")
            valid=np.asarray(valid,dtype=np.int64)
            valid_elements=elements[valid]
            material_effective=list(effective)
            material_elements=list(map(int,valid_elements))
            material_is_active=[False]*len(material_elements)
            material_volumes=list(map(float,volumes[valid_elements]))
            material_members=[np.asarray([int(element)],dtype=np.int64)
                              for element in valid_elements]
            removal_blocks=tuple(np.concatenate([
                element_blocks[int(member)]
                for member in removal_members[int(element)]])
                for element in removal_candidates)
            removal_linear_material=(None if not removal_blocks else
                _adjoint_corrected_removal_material_response(
                    charge_gram=charge_gram,inv_chi=inv_chi,
                    dof_blocks=removal_blocks,state=approximate_state,
                    response_matrix=linear_response_matrix,
                    screen_context=screen_context,response_band=band))
            for removal_column,element in enumerate(removal_candidates):
                local_material=removal_linear_material[:,removal_column]
                if objective_projection is None:
                    try:
                        _,removed,_=calibrate_source(
                            approximate_state,
                            approximate_response-local_material)
                        material_response=base_objective-transform_response(removed)
                    except (RuntimeError,ValueError):
                        continue
                else:
                    # ``removal_linear_material`` was assembled with
                    # ``linear_response_matrix`` and is therefore already in
                    # the projected transfer-map response space.
                    material_response=local_material
                removal_material_response[int(element)]=np.asarray(
                    material_response,dtype=float).copy()
                # TSVD columns represent positive material.  Removing the
                # active cell is the negative of this column; the recovered
                # coefficient sign chooses which operation is feasible.
                material_effective.append(material_response)
                material_elements.append(int(element))
                material_is_active.append(True)
                material_volumes.append(float(np.sum(
                    volumes[removal_members[int(element)]])))
                material_members.append(np.asarray(
                    removal_members[int(element)],dtype=np.int64))
            material_elements=np.asarray(material_elements,dtype=np.int64)
            material_matrix=np.column_stack(material_effective)
            material_is_active=np.asarray(material_is_active,dtype=bool)
            material_volumes=np.asarray(material_volumes,dtype=float)
            secondary_cost=None
            if graph_enabled:
                from ._topopt_graph import (
                    binary_graph_interface_energy,
                    candidate_face_adjacency,
                )
                directions=np.where(material_is_active,-1,1).astype(np.int8)
                signed_matrix=material_matrix*directions[None,:]
                candidate_graph=candidate_face_adjacency(
                    material_members,element_graph)
                base_interface=binary_graph_interface_energy(
                    active,element_graph,exterior_degree=element_exterior,
                    edge_weights=element_interface_weights)
                interface_delta=[]
                for direction,members in zip(directions,material_members):
                    trial=active.copy();trial[members]=(direction>0)
                    interface_delta.append(binary_graph_interface_energy(
                        trial,element_graph,
                        exterior_degree=element_exterior,
                        edge_weights=element_interface_weights)-base_interface)
                interface_delta=np.asarray(interface_delta,dtype=float)
                interface_scale=max(
                    1.0,float(np.max(np.abs(interface_delta))))
                volume_scale=max(
                    1.0e-300,float(np.max(material_volumes)))
                secondary_cost=(
                    0.05*material_volumes/volume_scale
                    +graph_interface_weight*interface_delta/interface_scale)
            tsvd_material_data=(
                np.asarray(base_objective,dtype=float).copy(),
                material_elements.copy(),material_matrix.copy(),
                material_volumes.copy(),material_is_active.copy())
            tsvd_proposal=select_tsvd_element_candidates(
                current_response=base_objective,response_target=target,
                response_band=band,candidate_elements=material_elements,
                candidate_response_delta=material_matrix,
                candidate_volumes=material_volumes,
                volume_budget=max(0.0,remaining),active_elements=active,
                predecessor_elements=predecessors,
                candidate_material_active=material_is_active,
                relative_tolerance=float(tsvd_relative_tolerance),
                improvement_capture=float(batch_improvement_capture),
                ratio_tolerance=ratio_tolerance,
                maximum_changed_volume=material_trust_volume,
                maximum_changed_elements=maximum_cap,
                candidate_secondary_cost=secondary_cost)
            if graph_enabled and graph_front_proposal_limit:
                graph_front_data=dict(
                    base=np.asarray(base_objective,dtype=float).copy(),
                    elements=material_elements.copy(),
                    directions=np.where(
                        material_is_active,-1,1).astype(np.int8),
                    signed_matrix=(material_matrix*np.where(
                        material_is_active,-1,1)[None,:]),
                    volumes=material_volumes.copy(),
                    members=tuple(value.copy() for value in material_members),
                    adjacency=candidate_graph,
                    secondary=(None if secondary_cost is None else
                               np.asarray(secondary_cost,dtype=float).copy()))
            removal_material_columns=np.flatnonzero(material_is_active)
            if removal_material_columns.size:
                requested=(removal_cluster_count if removal_cluster_count>0
                    else max(8,2*int(tsvd_proposal.numerical_rank)+2))
                cluster_labels,_=_configured_candidate_cluster_labels(
                    charge_gram,removal_blocks,
                    min(len(removal_blocks),requested))
                proposal_elements=np.asarray(
                    tsvd_proposal.selected_elements,dtype=np.int64)
                proposal_directions=np.asarray(
                    tsvd_proposal.selected_directions,dtype=np.int8)
                representative_elements=np.asarray(
                    tsvd_proposal.representative_elements,dtype=np.int64)
                representative_directions=np.asarray(
                    tsvd_proposal.representative_directions,dtype=np.int8)
                removal_cluster_front=_clustered_tsvd_candidate_front(
                    candidate_elements=material_elements[
                        removal_material_columns],
                    candidate_response_delta=material_matrix[:,
                        removal_material_columns],
                    response_band=band,cluster_labels=cluster_labels,
                    signed_coefficients=np.asarray(
                        tsvd_proposal.signed_coefficients,dtype=float)[
                            removal_material_columns],
                    selected_elements=proposal_elements[
                        proposal_directions<0],
                    representative_elements=representative_elements[
                        representative_directions<0],
                    relative_tolerance=float(tsvd_relative_tolerance),
                    front_limit=int(exact_candidate_limit))
            selected=np.asarray(tsvd_proposal.selected_elements,dtype=np.int64)[
                np.asarray(tsvd_proposal.selected_directions)>0]
            representatives=np.asarray(
                tsvd_proposal.representative_elements,dtype=np.int64)[
                np.asarray(tsvd_proposal.representative_directions)>0]
            positive_probe=np.union1d(selected,representatives).astype(np.int64)
            effective_matrix=np.column_stack(effective)
            addition_score=np.max(np.abs((
                base_objective[:,None]+effective_matrix-target[:,None]) /
                band[:,None]),axis=0)
            score_by_element={int(element):float(addition_score[column])
                              for column,element in enumerate(valid_elements)}
            ranked_additions=valid_elements[np.argsort(
                addition_score,kind="stable")]
            # This front is reached only after the global signed proposal has
            # been checked by one complete active-set solve.  Keep its dense
            # block-Schur challenge deliberately small: reducing hundreds of
            # BDM blocks forms a multi-gigabyte candidate-active matrix and
            # defeats the all-candidate low-rank screen.
            fallback_limit=min(16,int(exact_candidate_limit),max(
                2,len(positive_probe),
                2*int(tsvd_proposal.numerical_rank)+2))
            fallback_pool=np.union1d(
                positive_probe,ranked_additions[:fallback_limit])
            fallback_addition_elements=np.asarray(sorted(
                (int(element) for element in fallback_pool),
                key=lambda element:(score_by_element.get(element,np.inf),
                                    element))[:fallback_limit],dtype=np.int64)
            front=np.union1d(selected,representatives).astype(np.int64)
            if proposal_adjoint_rows.size!=target.size:
                # The production global screen solves only the minimax working
                # set of active-relaxation adjoints.  Its nonempty signed
                # proposal is verified by a full active-set solve below.  If it
                # is empty or rejected, the bounded fallback front is recomputed
                # with the complete adjoint Schur model.
                return np.empty(0,dtype=np.int64)
            if np.asarray(tsvd_proposal.selected_elements).size:
                # The global ACA/QR/TSVD+MILP proposal already determines the
                # complete signed batch.  Verify that batch by one physical
                # active-set solve before paying for any dense candidate Schur
                # matrix.  The latter is only a rejection fallback on a small
                # representative front; it is never the all-candidate screen.
                return np.empty(0,dtype=np.int64)
            # Exhaustive exact Schur search remains useful for tiny regression
            # fronts, but scaling it to every BDM1 element defeats the global
            # ACA/QR/TSVD screen.  Above eight elements retain only the signed
            # TSVD/QR representatives, bounded by the explicit exact-front
            # limit; accepted batch cardinality is still decided by the LP.
            exhaustive_limit=min(8,int(exact_candidate_limit))
            if valid_elements.size<=exhaustive_limit:
                return valid_elements
            if front.size:
                selected_set=set(int(value) for value in selected)
                # The global TSVD/MILP proposal determines batch cardinality;
                # an exact-front safety limit must never truncate that proposal
                # (which would turn a cooperative move back into arbitrary
                # "try N elements").  The limit bounds only extra QR
                # representatives used to challenge the proposal.
                proposed=np.asarray(sorted(selected_set),dtype=np.int64)
                extras=np.asarray(sorted(
                    int(value) for value in front
                    if int(value) not in selected_set),dtype=np.int64)
                extra_limit=max(0,int(exact_candidate_limit)-len(proposed))
                return np.r_[proposed,extras[:extra_limit]].astype(
                    np.int64,copy=False)
            # A no-improvement TSVD proposal still receives one exact Schur
            # probe.  When the zero-rank set is small, keep it whole so a pair
            # that improves only jointly is not discarded.  The limit controls
            # exact look-ahead work, not accepted batch cardinality.
            score=np.max(np.abs((base_objective[:,None]+effective_matrix-
                                 target[:,None])/band[:,None]),axis=0)
            return valid_elements[np.asarray([int(np.argmin(score))])]

        common_linearization=dict(charge_gram=charge_gram,fes=fes,
            inv_chi=inv_chi,rhs=rhs,response_matrix=linear_response_matrix,
            active_elements=active,candidate_elements=candidates,
            incident_response=linear_incident_response,
            solve_tolerance=proposal_tolerance,
            solve_max_iterations=solve_max_iterations,
            candidate_batch_size=candidate_batch_size,mass_riesz=mass_riesz,
            cluster_coarse_size=cluster_coarse_size,
            cluster_deflation_size=cluster_deflation_size,
            recycle_size=recycle_size,
            active_state=np.asarray(state,dtype=float)/float(source_scale),
            candidate_screen_context=screen_context)
        if len(candidates)>8 and working_adjoint_count>0:
            direct_capture={}
            def capture_direct_candidates(elements,approximate_delta,
                                          approximate_state,
                                          approximate_response):
                direct_capture["elements"]=np.asarray(
                    elements,dtype=np.int64).copy()
                direct_capture["delta"]=np.asarray(
                    approximate_delta,dtype=float).copy()
                return np.empty(0,dtype=np.int64)
            linearize_hdiv_mmm_element_generation(
                **common_linearization,
                candidate_selector=capture_direct_candidates,
                screen_with_adjoint=False)
            direct_delta=direct_capture.get("delta")
            direct_elements=direct_capture.get("elements")
            if (direct_delta is None or direct_delta.shape[0]!=target.size or
                    direct_elements is None):
                raise RuntimeError(
                    "candidate-response QR screen returned incompatible rows")
            from scipy.linalg import qr
            normalized_direct=direct_delta/band[:,None]
            _,_,row_pivots=qr(
                normalized_direct.T,mode="economic",pivoting=True)
            worst=int(np.argmax(np.abs(
                (objective_response-target)/band)))
            row_order=[worst]+[
                int(row) for row in row_pivots if int(row)!=worst]
            proposal_adjoint_rows=np.asarray(
                row_order[:min(working_adjoint_count,target.size)],
                dtype=np.int64)
            use_proposal_adjoint=proposal_adjoint_rows.size>0
            def select_interpolated_candidates(elements,approximate_delta,
                                               approximate_state,
                                               approximate_response):
                elements=np.asarray(elements,dtype=np.int64)
                if not np.array_equal(elements,direct_elements):
                    raise RuntimeError(
                        "candidate-response QR screen changed candidate order")
                interpolated=_interpolate_screened_response_correction(
                    direct_delta,approximate_delta,proposal_adjoint_rows,band)
                return select_exact_candidates(
                    elements,interpolated,approximate_state,
                    approximate_response)
            lin=linearize_hdiv_mmm_element_generation(
                **common_linearization,
                candidate_selector=select_interpolated_candidates,
                screen_with_adjoint=use_proposal_adjoint,
                screen_adjoint_rows=proposal_adjoint_rows)
        else:
            use_proposal_adjoint=proposal_adjoint_rows.size>0
            lin=linearize_hdiv_mmm_element_generation(
                **common_linearization,
                candidate_selector=select_exact_candidates,
                screen_with_adjoint=use_proposal_adjoint,
                screen_adjoint_rows=proposal_adjoint_rows)
        if objective_projection is None:
            _,lp_raw_response,linear_scale=calibrate_source(
                lin.state,lin.response)
            lp_response=transform_response(lp_raw_response)
            state=lin.state*linear_scale
            response=lp_raw_response
            source_scale=linear_scale
            objective_response=lp_response
        else:
            # The projected base response lacks the nonlinear affine offset;
            # retain the last exact raw solve and use the new unscaled state
            # only for the reduced candidate system.
            state=lin.state*source_scale
            lp_response=np.asarray(objective_response,dtype=float)
        current_ratio=ratio(objective_response)
        if tsvd_proposal is None:
            raise RuntimeError("all-candidate TSVD selector did not run")
        move_elements=np.asarray(tsvd_proposal.selected_elements,dtype=np.int64)
        move_directions=np.asarray(tsvd_proposal.selected_directions,dtype=np.int8)
        proposed_additions=move_elements[move_directions>0]
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
        graph_proposal_count=0

        # The global signed proposal -- additions, removals, or both -- is
        # tested by one complete re-solve.  Poor model agreement shrinks a
        # volume trust region and resolves the binary LP over every original
        # candidate.  This is not an arbitrary ``try k elements`` rule: the LP
        # still determines membership and cardinality, while the physical
        # re-solve controls only the admissible size of the material move.
        # No finite difference or gray material is introduced.
        if move_elements.size or graph_front_data is not None:
            exact_evaluated=0
            trust_proposals=[tsvd_proposal]
            if trust_region_trials>1 and tsvd_material_data is not None:
                (material_base,material_elements,material_matrix,
                 material_volumes,material_is_active)=tsvd_material_data
                volume_lookup={int(element):float(material_volumes[column])
                    for column,element in enumerate(material_elements)}
                proposed_changed_volume=float(sum(
                    volume_lookup[int(element)]
                    for element in move_elements))
                minimum_volume=float(np.min(material_volumes))
                initial_trust=(proposed_changed_volume if
                    material_trust_volume is None else
                    float(material_trust_volume))
                seen={tuple(sorted(zip(
                    map(int,move_elements),map(int,move_directions))))}
                for shrink in range(1,trust_region_trials):
                    trust_budget=max(
                        minimum_volume,initial_trust*(0.5**shrink))
                    if trust_budget>=proposed_changed_volume-ratio_tolerance:
                        continue
                    alternate=select_tsvd_element_candidates(
                        current_response=material_base,
                        response_target=target,response_band=band,
                        candidate_elements=material_elements,
                        candidate_response_delta=material_matrix,
                        candidate_volumes=material_volumes,
                        volume_budget=max(0.0,remaining),
                        active_elements=active,
                        predecessor_elements=predecessors,
                        candidate_material_active=material_is_active,
                        relative_tolerance=float(tsvd_relative_tolerance),
                        improvement_capture=float(batch_improvement_capture),
                        ratio_tolerance=ratio_tolerance,
                        maximum_changed_volume=trust_budget,
                        maximum_changed_elements=maximum_cap)
                    alternate_elements=np.asarray(
                        alternate.selected_elements,dtype=np.int64)
                    alternate_directions=np.asarray(
                        alternate.selected_directions,dtype=np.int8)
                    key=tuple(sorted(zip(map(int,alternate_elements),
                                         map(int,alternate_directions))))
                    if alternate_elements.size and key not in seen:
                        seen.add(key);trust_proposals.append(alternate)
                    if trust_budget<=minimum_volume+ratio_tolerance:
                        break

            if graph_front_data is not None:
                from ._topopt_graph import (
                    connected_graph_front_beam,minimax_driving_potential)
                graph_elements=graph_front_data["elements"]
                graph_directions=graph_front_data["directions"]
                graph_volumes=graph_front_data["volumes"]
                graph_members=graph_front_data["members"]
                graph_lookup={(int(element),int(direction)):column
                    for column,(element,direction) in enumerate(zip(
                        graph_elements,graph_directions))}
                seed_indices=[]
                for element,direction in zip(
                        np.r_[tsvd_proposal.selected_elements,
                              tsvd_proposal.representative_elements],
                        np.r_[tsvd_proposal.selected_directions,
                              tsvd_proposal.representative_directions]):
                    column=graph_lookup.get((int(element),int(direction)))
                    if column is not None and column not in seed_indices:
                        seed_indices.append(column)
                if not seed_indices:
                    raw_drive,_=minimax_driving_potential(
                        graph_front_data["base"],target,band,
                        graph_front_data["signed_matrix"])
                    seed_indices=[int(np.argmax(raw_drive))]

                physical_max=len(graph_elements)
                if material_trust_before is not None:
                    cumulative=np.cumsum(np.sort(graph_volumes))
                    physical_max=max(1,int(np.count_nonzero(
                        cumulative<=material_trust_before+ratio_tolerance)))
                if maximum_cap is not None:
                    physical_max=min(physical_max,maximum_cap)
                rank_budget=max(1,int(tsvd_proposal.numerical_rank)+1,
                    int(len(tsvd_proposal.selected_elements)))
                if graph_front_budget is None:
                    graph_front_budget=min(physical_max,rank_budget)
                else:
                    graph_front_budget=min(physical_max,
                                           max(1,int(graph_front_budget)))

                def graph_bundle_mask(bundle):
                    trial=active.copy()
                    for column in bundle:
                        trial[graph_members[int(column)]]=(
                            graph_directions[int(column)]>0)
                    return trial

                def graph_bundle_valid(bundle):
                    columns=np.asarray(bundle,dtype=np.int64)
                    if (material_trust_before is not None and
                            float(np.sum(graph_volumes[columns]))>
                            material_trust_before+ratio_tolerance):
                        return False
                    trial=graph_bundle_mask(bundle)
                    return (float(volumes@trial)<=float(volume_max)+1e-14 and
                            valid_active_set(trial))

                secondary=graph_front_data["secondary"]
                secondary_scale=(1.0 if secondary is None else
                    max(1.0,float(np.max(np.abs(secondary)))))
                def graph_regularization(bundle):
                    if secondary is None:
                        return 0.0
                    # Strictly a near-tie regularizer; raw response columns
                    # and the exact physical acceptance ratio remain unchanged.
                    return (1.0e-9/secondary_scale)*float(
                        np.sum(secondary[np.asarray(bundle,dtype=np.int64)]))

                graph_result=connected_graph_front_beam(
                    current_response=graph_front_data["base"],
                    response_target=target,response_band=band,
                    candidate_response_delta=graph_front_data["signed_matrix"],
                    adjacency=graph_front_data["adjacency"],
                    seed_indices=seed_indices,
                    maximum_size=graph_front_budget,
                    maximum_components=graph_front_maximum_components,
                    beam_width=graph_front_beam_width,
                    proposal_limit=graph_front_proposal_limit,
                    response_novelty_weight=(
                        graph_front_response_novelty_weight),
                    return_result=True,
                    regularization_change=graph_regularization,
                    is_valid=graph_bundle_valid)
                graph_proposals=graph_result.proposals
                pool_diag=graph_result.pool_diagnostics
                selected_diag=graph_result.selected_diagnostics
                graph_front_diagnostics.append(
                    HDivMMMGraphFrontDiagnostics(
                        int(iteration),int(search_depth),
                        int(len(graph_elements)),
                        float(graph_front_response_novelty_weight),
                        int(pool_diag.proposal_count),
                        int(pool_diag.numerical_rank),
                        float(pool_diag.duplicate_pair_fraction),
                        float(pool_diag.maximum_absolute_correlation),
                        int(selected_diag.proposal_count),
                        int(selected_diag.numerical_rank),
                        float(selected_diag.duplicate_pair_fraction),
                        float(selected_diag.minimum_subspace_novelty)))
                seen_proposals={tuple(sorted(zip(
                    map(int,proposal.selected_elements),
                    map(int,proposal.selected_directions))))
                    for proposal in trust_proposals}
                for graph_proposal in graph_proposals:
                    columns=np.asarray(
                        graph_proposal.candidate_indices,dtype=np.int64)
                    key=tuple(sorted(zip(
                        map(int,graph_elements[columns]),
                        map(int,graph_directions[columns]))))
                    if not key or key in seen_proposals:
                        continue
                    seen_proposals.add(key)
                    trust_proposals.append(TSVDElementCandidateSelection(
                        graph_elements[columns].copy(),
                        graph_directions[columns].copy(),
                        np.asarray(tsvd_proposal.representative_elements,
                                   dtype=np.int64).copy(),
                        np.asarray(tsvd_proposal.representative_directions,
                                   dtype=np.int8).copy(),
                        np.asarray(graph_proposal.predicted_response,
                                   dtype=float).copy(),
                        float(graph_proposal.predicted_max_band_ratio),
                        float(np.sum(graph_volumes[columns]*
                                     graph_directions[columns])),
                        int(tsvd_proposal.numerical_rank),
                        int(tsvd_proposal.aca_rank),
                        np.asarray(tsvd_proposal.singular_values,
                                   dtype=float).copy(),
                        np.asarray(tsvd_proposal.signed_coefficients,
                                   dtype=float).copy(),
                        float(tsvd_proposal.relative_truncation_error),
                        "connected ACA/QR/TSVD-seeded graph-front",
                        tsvd_proposal.linearized_reachability_residual,
                        tsvd_proposal.linearized_reachability_max_band_ratio,
                        tsvd_proposal.linearized_reachability_relative_residual,
                        tsvd_proposal.linearized_reachable))
                    graph_proposal_count+=1

            best_trial=None
            meaningful_ratio_tolerance=(
                ratio_tolerance*max(1.0,float(current_ratio)))
            for proposal_index,proposal in enumerate(trust_proposals):
                trial_elements=np.asarray(
                    proposal.selected_elements,dtype=np.int64)
                trial_directions=np.asarray(
                    proposal.selected_directions,dtype=np.int8)
                if trial_elements.size==0:
                    continue
                trial_additions=trial_elements[trial_directions>0]
                trial_removals=trial_elements[trial_directions<0]
                trial_active=active.copy();trial_active[trial_additions]=True
                expanded_removals=(np.empty(0,dtype=np.int64)
                    if trial_removals.size==0 else np.unique(np.concatenate([
                        removal_members[int(element)]
                        for element in trial_removals])))
                trial_active[expanded_removals]=False
                if not valid_active_set(trial_active):
                    continue
                exact_evaluated+=1
                trial_state,trial_response,trial_iterations=\
                    solve_hdiv_mmm_active_elements(
                        charge_gram=charge_gram,fes=fes,inv_chi=inv_chi,
                        rhs=rhs,response_matrix=response_matrix,
                        active_elements=trial_active,
                        incident_response=incident_response,
                        solve_tolerance=solve_tolerance,
                        solve_max_iterations=solve_max_iterations,
                        mass_riesz=mass_riesz,
                        cluster_coarse_size=cluster_coarse_size,
                        cluster_deflation_size=cluster_deflation_size,
                        recycle_size=recycle_size)
                try:
                    trial_state,trial_response,trial_scale=calibrate_source(
                        trial_state,trial_response)
                    trial_objective=transform_response(trial_response)
                    if not exact_response_is_valid(
                            trial_response,trial_objective):
                        continue
                    actual=ratio(trial_objective)
                except (RuntimeError,ValueError):
                    continue
                predicted_ratio=float(proposal.predicted_max_band_ratio)
                exact_trial_states.append(record_exact_trial(exact_snapshot(
                    active_value=trial_active,state_value=trial_state,
                    response_value=trial_response,
                    objective_value=trial_objective,ratio_value=actual,
                    scale_value=trial_scale,
                    iterations_value=trial_iterations,
                    depth_value=search_depth+1,
                    path_value=search_path+(tuple(sorted(zip(
                        map(int,trial_elements),map(int,trial_directions)))),))))
                predicted_reduction=current_ratio-predicted_ratio
                model_agreement=((current_ratio-actual)/predicted_reduction
                    if predicted_reduction>ratio_tolerance else 0.0)
                if (actual<current_ratio-meaningful_ratio_tolerance and
                        (best_trial is None or actual<best_trial[0])):
                    best_trial=(
                        actual,trial_additions,expanded_removals,trial_active,
                        trial_state,trial_response,trial_iterations,trial_scale,
                        trial_objective,predicted_ratio,proposal_index,
                        model_agreement)
                # A locally well-predicted singleton must not hide a
                # collaborative bundle that is needed to enter the response
                # band.  The proposal list is deliberately bounded, so keep
                # resolving it until one proposal satisfies every band; only
                # then is an early exit lossless for the feasibility goal.
                if actual<=1.0+ratio_tolerance:
                    break
            if best_trial is not None:
                (actual,selected,removed_elements,trial_active,trial_state,
                 trial_response,trial_iterations,trial_scale,trial_objective,
                 predicted_ratio,proposal_index,_)=best_trial
                accepted_proposal=trust_proposals[int(proposal_index)]
                if accepted_proposal.status.startswith("connected"):
                    selection_model=(
                        "aca-qr-tsvd-connected-graph-front-full-resolve")
                else:
                    selection_model=(
                        "all-candidate-aca-qr-tsvd-direct-full-resolve"
                        if proposal_index==0 else
                        "all-candidate-aca-qr-tsvd-adaptive-trust-full-resolve")
                mixed_accepted=True

            if not mixed_accepted:
                # The additive magnetization deletion column is only a cheap
                # proposal; it is not an exact active-system downdate.  Do not
                # return to every exposed deletion after a rejected global
                # move.  The native charge-tree clusters plus local
                # ACA--QR--TSVD retain spatially distinct response modes in a
                # bounded removal front; the minimax LP and trust-region
                # checks operate only on that compressed front.
                available_removals=np.asarray([
                    int(element) for element in removal_cluster_front
                    if int(element) in removal_material_response],
                    dtype=np.int64)
                best_alternate=None
                if available_removals.size:
                    removal_matrix=np.column_stack([
                        removal_material_response[int(element)]
                        for element in available_removals])
                    removal_volumes=np.asarray([
                        np.sum(volumes[removal_members[int(element)]])
                        for element in available_removals],dtype=float)
                    unrestricted=solve_element_generation_lp(
                        objective_response,target,band,-removal_matrix,
                        removal_volumes,volume_budget=max(0.0,remaining),
                        maximum_new_elements=maximum_cap,
                        whole_elements=True,
                        candidate_volume_change=-removal_volumes,
                        maximum_changed_volume=material_trust_before)
                    initial_count=int(np.count_nonzero(unrestricted.selected))
                    if initial_count==0:
                        initial_count=1
                    counts=[];count=initial_count
                    for _ in range(trust_region_trials):
                        count=max(1,int(count))
                        if count not in counts:
                            counts.append(count)
                        if count==1:
                            break
                        count=max(1,int(np.ceil(0.5*count)))
                    seen_removal_batches=set()
                    for count in counts:
                        proposal=(unrestricted if
                            count==initial_count and
                            np.count_nonzero(unrestricted.selected) else
                            solve_element_generation_lp(
                                objective_response,target,band,
                                -removal_matrix,removal_volumes,
                                volume_budget=max(0.0,remaining),
                                maximum_new_elements=count,
                                whole_elements=True,
                                candidate_volume_change=-removal_volumes,
                                maximum_changed_volume=material_trust_before))
                        selected_groups=available_removals[
                            np.asarray(proposal.selected,dtype=bool)]
                        if selected_groups.size==0:
                            continue
                        expanded=np.unique(np.concatenate([
                            removal_members[int(element)]
                            for element in selected_groups]))
                        batch_key=tuple(map(int,expanded))
                        if batch_key in seen_removal_batches:
                            continue
                        seen_removal_batches.add(batch_key)
                        alternate_active=active.copy()
                        alternate_active[expanded]=False
                        if not valid_active_set(alternate_active):
                            continue
                        exact_evaluated+=1
                        (alternate_state,alternate_response,
                         alternate_iterations)=solve_hdiv_mmm_active_elements(
                            charge_gram=charge_gram,fes=fes,
                            inv_chi=inv_chi,rhs=rhs,
                            response_matrix=response_matrix,
                            active_elements=alternate_active,
                            incident_response=incident_response,
                            solve_tolerance=solve_tolerance,
                            solve_max_iterations=solve_max_iterations,
                            mass_riesz=mass_riesz,
                            cluster_coarse_size=cluster_coarse_size,
                            cluster_deflation_size=cluster_deflation_size,
                            recycle_size=recycle_size)
                        try:
                            (alternate_state,alternate_response,
                             alternate_scale)=calibrate_source(
                                alternate_state,alternate_response)
                            alternate_objective=transform_response(
                                alternate_response)
                            if not exact_response_is_valid(
                                    alternate_response,alternate_objective):
                                continue
                            alternate_ratio=ratio(alternate_objective)
                        except (RuntimeError,ValueError):
                            continue
                        predicted=float(proposal.predicted_max_band_ratio)
                        predicted_reduction=current_ratio-predicted
                        model_agreement=(
                            (current_ratio-alternate_ratio)/predicted_reduction
                            if predicted_reduction>ratio_tolerance else 0.0)
                        if (alternate_ratio<current_ratio-ratio_tolerance and
                                (best_alternate is None or
                                 alternate_ratio<best_alternate[0])):
                            best_alternate=(alternate_ratio,expanded,
                                alternate_active,alternate_state,
                                alternate_response,alternate_iterations,
                                alternate_scale,alternate_objective,predicted,
                                selected_groups,model_agreement)
                        if (alternate_ratio<current_ratio-ratio_tolerance and
                                model_agreement>=minimum_model_agreement):
                            break
                if best_alternate is not None:
                    (actual,removed_elements,trial_active,trial_state,
                     trial_response,trial_iterations,trial_scale,
                     trial_objective,predicted_ratio,selected_groups,
                     _)=best_alternate
                    selected=np.empty(0,dtype=np.int64)
                    selection_model=(
                        "signed-magnetization-aca-qr-tsvd-"
                        "alternate-removal-exact"
                        if selected_groups.size==1 else
                        "clustered-removal-aca-qr-tsvd-adaptive-trust-full-resolve")
                    mixed_accepted=True

        if not mixed_accepted and len(lin.candidate_elements)==0:
            # A production partial-adjoint screen has already challenged its
            # global proposal with adaptively shrunken physical moves.  Do not
            # silently expand it to every response adjoint: that defeated the
            # matrix-free search on the 63k-DOF end-pack model.  The exact
            # block-Schur fallback remains available for tiny/full-adjoint
            # fronts where its cost is explicitly bounded.
            if proposal_adjoint_rows.size!=target.size:
                branch=next_nonmonotone_state(exact_trial_states)
                if branch is not None:
                    active=branch.active;state=branch.state
                    response=branch.response
                    objective_response=branch.objective_response
                    current_ratio=branch.ratio
                    source_scale=branch.source_scale
                    solve_iterations=branch.solve_iterations
                    search_depth=branch.depth;search_path=branch.path
                    stop_reason="exact_nonmonotone_beam_in_progress"
                    continue
                stop_reason=("exact_nonmonotone_beam_exhausted" if
                    exact_beam_width else
                    "adaptive_trust_region_proposals_rejected")
                break
            if fallback_addition_elements.size==0:
                stop_reason=(
                    "no_insertion_front_after_removal_rejection")
                break
            def fallback_addition_front(elements,approximate_delta,
                                        approximate_state,
                                        approximate_response):
                elements=np.asarray(elements,dtype=np.int64).reshape(-1)
                preferred=np.intersect1d(
                    elements,fallback_addition_elements).astype(np.int64)
                if preferred.size:
                    return preferred[:min(16,int(exact_candidate_limit))]
                return np.empty(0,dtype=np.int64)
            lin=linearize_hdiv_mmm_element_generation(
                charge_gram=charge_gram,fes=fes,inv_chi=inv_chi,rhs=rhs,
                response_matrix=linear_response_matrix,
                active_elements=active,candidate_elements=candidates,
                incident_response=linear_incident_response,
                solve_tolerance=solve_tolerance,
                solve_max_iterations=solve_max_iterations,
                candidate_batch_size=candidate_batch_size,
                mass_riesz=mass_riesz,
                cluster_coarse_size=cluster_coarse_size,
                cluster_deflation_size=cluster_deflation_size,
                recycle_size=recycle_size,
                candidate_selector=fallback_addition_front,
                active_state=np.asarray(state,dtype=float)/float(source_scale),
                screen_with_adjoint=True)
            if len(lin.candidate_elements)==0:
                stop_reason="full_solve_rejected_removal_only_tsvd";break

        def evaluate_bundle(elements):
            block=hdiv_mmm_block_insertion_response(lin,elements)
            if objective_projection is not None:
                return objective_response+block.response_delta
            inserted=lin.response+block.response_delta
            try:
                _,calibrated,_=calibrate_source(lin.state,inserted)
                return transform_response(calibrated)
            except (RuntimeError,ValueError):
                return None

        if not mixed_accepted:
            exact_candidates=np.asarray(
                lin.candidate_elements,dtype=np.int64).reshape(-1)
            exact_proposal=np.intersect1d(
                proposed_additions,exact_candidates).astype(np.int64)
            exact_representatives=np.intersect1d(
                (exact_candidates
                 if lin.available_candidate_count<=int(exact_candidate_limit)
                 else representative_additions),
                exact_candidates).astype(np.int64)
            exact=select_tsvd_exact_block_batch(
                current_response=lp_response,response_target=target,
                response_band=band,candidate_elements=exact_candidates,
                candidate_volumes=volumes[exact_candidates],
                proposal_elements=(exact_proposal if exact_proposal.size
                                   else exact_representatives),
                representative_elements=exact_representatives,
                evaluate_bundle_response=evaluate_bundle,
                volume_budget=(max(0.0,remaining) if
                    material_trust_before is None else min(
                        max(0.0,remaining),material_trust_before)),
                active_elements=active,
                predecessor_elements=predecessors,
                maximum_new_elements=maximum_cap,
                improvement_capture=float(batch_improvement_capture),
                ratio_tolerance=ratio_tolerance,
                bundle_is_valid=lambda elements:valid_active_set(
                    np.isin(np.arange(len(active)),elements)|active))
            selected=exact.selected_elements
            predicted_ratio=exact.predicted_max_band_ratio
            exact_evaluated=exact.evaluated_bundles
            selection_model="all-candidate-aca-qr-tsvd-exact-conditional"
            if selected.size==0 or predicted_ratio>=current_ratio-ratio_tolerance:
                branch=next_nonmonotone_state(exact_trial_states)
                if branch is not None:
                    active=branch.active;state=branch.state
                    response=branch.response
                    objective_response=branch.objective_response
                    current_ratio=branch.ratio
                    source_scale=branch.source_scale
                    solve_iterations=branch.solve_iterations
                    search_depth=branch.depth;search_path=branch.path
                    stop_reason="exact_nonmonotone_beam_in_progress"
                    continue
                stop_reason=("exact_nonmonotone_beam_exhausted" if
                    exact_beam_width else "no_improving_exact_bundle")
                break
            # Block Schur is exact for the current active set, but acceptance
            # remains tied to a fresh full solve and the topology gate.
            trial_active=active.copy();trial_active[selected]=True
            trial_topology=ngsolve_growth_topology(fes.mesh,trial_active)
            if not trial_topology.valid or not valid_active_set(trial_active):
                stop_reason="topology_gate_rejected_exact_bundle";break
            trial_state,trial_response,trial_iterations=\
                solve_hdiv_mmm_active_elements(
                    charge_gram=charge_gram,fes=fes,inv_chi=inv_chi,rhs=rhs,
                    response_matrix=response_matrix,active_elements=trial_active,
                    incident_response=incident_response,
                    solve_tolerance=solve_tolerance,
                    solve_max_iterations=solve_max_iterations,
                    mass_riesz=mass_riesz,
                    cluster_coarse_size=cluster_coarse_size,
                    cluster_deflation_size=cluster_deflation_size,
                    recycle_size=recycle_size)
            try:
                trial_state,trial_response,trial_scale=calibrate_source(
                    trial_state,trial_response)
                trial_objective=transform_response(trial_response)
            except (RuntimeError,ValueError):
                stop_reason="source_calibration_rejected_exact_bundle";break
            if not exact_response_is_valid(trial_response,trial_objective):
                stop_reason="response_guard_rejected_exact_bundle"
                branch=next_nonmonotone_state(exact_trial_states)
                if branch is not None:
                    active=branch.active;state=branch.state
                    response=branch.response
                    objective_response=branch.objective_response
                    current_ratio=branch.ratio
                    source_scale=branch.source_scale
                    solve_iterations=branch.solve_iterations
                    search_depth=branch.depth;search_path=branch.path
                    stop_reason="exact_nonmonotone_beam_in_progress"
                    continue
                break
            actual=ratio(trial_objective)
            exact_trial_states.append(record_exact_trial(exact_snapshot(
                active_value=trial_active,state_value=trial_state,
                response_value=trial_response,
                objective_value=trial_objective,ratio_value=actual,
                scale_value=trial_scale,iterations_value=trial_iterations,
                depth_value=search_depth+1,
                path_value=search_path+(tuple(
                    (int(element),1) for element in selected),))))
            if actual>=current_ratio-ratio_tolerance:
                branch=next_nonmonotone_state(exact_trial_states)
                if branch is not None:
                    active=branch.active;state=branch.state
                    response=branch.response
                    objective_response=branch.objective_response
                    current_ratio=branch.ratio
                    source_scale=branch.source_scale
                    solve_iterations=branch.solve_iterations
                    search_depth=branch.depth;search_path=branch.path
                    stop_reason="exact_nonmonotone_beam_in_progress"
                    continue
                stop_reason=("exact_nonmonotone_beam_exhausted" if
                    exact_beam_width else "full_solve_rejected_exact_bundle")
                break
        if (search_depth>0 and
                actual>=incumbent.ratio-ratio_tolerance):
            branch=next_nonmonotone_state(exact_trial_states)
            if branch is not None:
                active=branch.active;state=branch.state
                response=branch.response
                objective_response=branch.objective_response
                current_ratio=branch.ratio
                source_scale=branch.source_scale
                solve_iterations=branch.solve_iterations
                search_depth=branch.depth;search_path=branch.path
                stop_reason="exact_nonmonotone_beam_in_progress"
                continue
            stop_reason="exact_nonmonotone_beam_exhausted"
            break

        accepted_search_depth=search_depth
        predicted_improvement=current_ratio-predicted_ratio
        actual_improvement=current_ratio-actual
        agreement=(actual_improvement/predicted_improvement
                   if predicted_improvement>ratio_tolerance else 0.0)
        accepted_added=np.flatnonzero(trial_active & ~incumbent.active)
        accepted_removed=np.flatnonzero(incumbent.active & ~trial_active)
        changed_volume=float(
            np.sum(volumes[accepted_added])+np.sum(volumes[accepted_removed]))
        if graph_front_data is not None and graph_front_budget is not None:
            from ._topopt_graph import update_graph_front_trust
            graph_trust=update_graph_front_trust(
                budget=graph_front_budget,minimum_budget=1,
                maximum_budget=max(1,int(physical_max)),
                predicted_ratio=predicted_ratio,actual_ratio=actual,
                current_ratio=current_ratio,
                selected_size=max(
                    1,int(len(accepted_added)+len(accepted_removed))),
                interface_weight=graph_interface_weight)
            graph_front_budget=graph_trust.budget_after
            graph_interface_weight=graph_trust.interface_weight_after
        material_trust_volume=update_material_trust(
            material_trust_before,changed_volume,agreement)
        active=trial_active;state=trial_state;response=trial_response
        solve_iterations=trial_iterations;source_scale=trial_scale
        objective_response=trial_objective;current_ratio=actual
        batch_limit_after=len(accepted_added)+len(accepted_removed)
        history.append(HDivMMMGenerationIteration(iteration,
            lin.available_candidate_count+len(removal_candidates),
            accepted_added,predicted_ratio,
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
            accepted_removed,
            int(lin.available_candidate_count),int(len(removal_candidates)),
            int(lin.candidate_coupling_rank),
            float(lin.candidate_coupling_relative_truncation_error),
            int(max(lin.schur_iterations,default=0)),
            int(len(lin.adjoint_iterations)),
            (None if tsvd_proposal.linearized_reachability_residual is None else
             np.asarray(tsvd_proposal.linearized_reachability_residual,dtype=float)),
            float(tsvd_proposal.linearized_reachability_max_band_ratio),
            float(tsvd_proposal.linearized_reachability_relative_residual),
             bool(tsvd_proposal.linearized_reachable),
             material_trust_volume_before=material_trust_before,
             material_trust_volume_after=material_trust_volume,
             material_changed_volume=changed_volume,
             graph_front_proposals_evaluated=int(graph_proposal_count),
             nonmonotone_search_depth=int(accepted_search_depth),
             abe_murata_diagnostics=(
                 tsvd_proposal.abe_murata_diagnostics)))
        if iteration_callback is not None:
            iteration_callback(history[-1])
        incumbent=exact_snapshot(
            active_value=active,state_value=state,response_value=response,
            objective_value=objective_response,ratio_value=current_ratio,
            scale_value=source_scale,iterations_value=solve_iterations,
            depth_value=0,path_value=())
        search_depth=0;search_path=();pending_exact_states=[]
        visited_exact_states={np.packbits(active).tobytes()}
        converged=current_ratio<=1.0+ratio_tolerance
        if converged: stop_reason="target_met"
    if exact_beam_width:
        # Exploratory barrier states are never a deliverable design.  Return
        # only the best fully solved incumbent, even when the expansion budget
        # ends while the beam is away from it.
        active=incumbent.active;state=incumbent.state
        response=incumbent.response
        objective_response=incumbent.objective_response
        current_ratio=incumbent.ratio;source_scale=incumbent.source_scale
        converged=current_ratio<=1.0+ratio_tolerance
        if stop_reason=="exact_nonmonotone_beam_in_progress":
            stop_reason="exact_nonmonotone_beam_budget_exhausted"
    return HDivMMMGenerationResult(
        active,state,response,tuple(history),converged,source_scale,
        objective_response,stop_reason,tuple(exact_search_trace),
        tuple(graph_front_diagnostics))


def solve_element_generation_lp(current_response, response_target,
        response_band, candidate_response_delta, candidate_volumes, *,
        volume_budget, maximum_new_elements=None,
        candidate_objective_change=None, whole_elements=True,
        relative_mip_gap=0.0, predecessor_pairs=None,
        predicted_ratio_cap=None,
        candidate_volume_change=None,
        candidate_exclusion_groups=None,
        maximum_changed_volume=None) -> ElementGenerationLPUpdate:
    """Select a small full-strength element-growth batch by a 0-1 LP.

    Single-element Schur responses are superposed only for the selection model;
    the committed batch must subsequently be solved as one exact HDiv-MMM
    problem.  With ``whole_elements=True`` (the default), HiGHS receives binary
    element variables, so no gray material can enter the physical model.
    ``candidate_volume_change`` owns the signed net-volume constraint, whereas
    ``maximum_changed_volume`` bounds the positive total flipped volume.  The
    latter is the discrete trust region required when additions and removals
    coexist.
    """
    from scipy.optimize import Bounds, LinearConstraint, milp

    y=np.asarray(current_response,dtype=float).reshape(-1)
    target=np.asarray(response_target,dtype=float).reshape(-1)
    band=np.asarray(response_band,dtype=float).reshape(-1)
    D=np.asarray(candidate_response_delta,dtype=float)
    volumes=np.asarray(candidate_volumes,dtype=float).reshape(-1)
    volume_change=(volumes.copy() if candidate_volume_change is None else
                   np.asarray(candidate_volume_change,dtype=float).reshape(-1))
    exclusion_groups=(None if candidate_exclusion_groups is None else
                      np.asarray(candidate_exclusion_groups,dtype=np.int64).reshape(-1))
    if y.size==0 or target.shape!=y.shape or band.shape!=y.shape or np.any(band<=0):
        raise ValueError("element-generation response/target/band vectors are invalid")
    if D.ndim!=2 or D.shape[0]!=y.size or D.shape[1]!=volumes.size or volumes.size==0:
        raise ValueError("candidate_response_delta must have shape (n_response,n_candidate)")
    if (np.any(volumes<=0) or volume_change.shape!=volumes.shape or
            (exclusion_groups is not None and
             exclusion_groups.shape!=volumes.shape) or
            not np.all(np.isfinite(volume_change)) or
            not np.isfinite(volume_budget) or volume_budget<0):
        raise ValueError("candidate volumes must be positive and volume_budget nonnegative")
    if not np.all(np.isfinite(np.r_[y,target,band,D.ravel(),volumes])):
        raise ValueError("element-generation LP inputs must be finite")
    nc=volumes.size
    maximum=nc if maximum_new_elements is None else int(maximum_new_elements)
    if maximum<0: raise ValueError("maximum_new_elements must be nonnegative")
    changed_volume=(None if maximum_changed_volume is None else
                    float(maximum_changed_volume))
    if (changed_volume is not None and
            (not np.isfinite(changed_volume) or changed_volume<0.0)):
        raise ValueError(
            "maximum_changed_volume must be finite and nonnegative")
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
    if changed_volume is not None:
        rows.append(np.r_[volumes,0.0][None,:])
        upper_parts.append(np.array([changed_volume]))
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
    if exclusion_groups is not None:
        labels=np.unique(exclusion_groups[exclusion_groups>=0])
        if labels.size:
            exclusive_rows=np.zeros((len(labels),nc+1))
            for row,label in enumerate(labels):
                exclusive_rows[row,np.flatnonzero(
                    exclusion_groups==label)]=1.0
            rows.append(exclusive_rows);upper_parts.append(np.ones(len(labels)))
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
        if (np.any(np.isnan(lower)) or np.any(np.isnan(upper))
                or np.any(lower>upper) or np.any(q<lower) or np.any(q>upper)):
            raise ValueError("invalid shape parameter bounds or current parameters outside them")
    lower=np.maximum(lower,q-move); upper=np.minimum(upper,q+move)
    if np.any(lower>upper):
        raise ValueError("shape trust region does not intersect the parameter bounds")

    static_rows=[]; static_rhs=[]
    if laplacian is not None:
        L=np.atleast_2d(np.asarray(laplacian,dtype=float))
        if L.shape[1]!=n or curvature_limit is None or not np.all(np.isfinite(L)):
            raise ValueError("shape laplacian requires n columns and curvature_limit")
        limit=np.broadcast_to(np.asarray(curvature_limit,dtype=float),(L.shape[0],))
        if np.any(limit<0.0) or not np.all(np.isfinite(limit)):
            raise ValueError("curvature_limit must be finite and nonnegative")
        static_rows.extend([L,-L]); static_rhs.extend([limit,limit])
    if A_ub is not None:
        if b_ub is None:
            raise ValueError("shape A_ub requires b_ub")
        extra=np.atleast_2d(np.asarray(A_ub,dtype=float))
        rhs=np.asarray(b_ub,dtype=float).reshape(-1)
        if (extra.shape!=(rhs.size,n) or not np.all(np.isfinite(extra))
                or not np.all(np.isfinite(rhs))):
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


__all__=["VIMLinearization","VIMStateShapeJacobian",
         "VIMOperatorLinearization","ChargeGramLinearization",
         "ChargeGramDirectionalOperators","VIMMatrixFreeLinearization","VIMAdjointGradient",
         "VIMFunctionalShapeJacobian",
         "ProductionGetTrafoDisplacements","ProductionVIMLinearization","LPUpdate",
         "ShapeLinearization","ShapeLPUpdate","solve_shape_lp",
         "ElementInsertionResponse","ElementGenerationLPUpdate",
         "HDivMMMElementGenerationLinearization",
         "HDivMMMSingleRemovalResponses",
         "HDivMMMRemovalGroupResponses",
         "AbeMurataEquivalentMaterialDiagnostics",
         "TSVDElementCandidateSelection",
         "HDivMMMGenerationIteration","HDivMMMGraphFrontDiagnostics",
         "HDivMMMExactSearchTrial",
         "HDivMMMGenerationResult",
         "GrowthTopologyReport","ngsolve_growth_topology",
         "finite_element_insertion_response","ngsolve_boundary_growth_candidates",
         "ngsolve_boundary_removal_candidates",
         "ngsolve_discontinuous_element_dof_blocks",
         "linearize_hdiv_mmm_element_generation","solve_hdiv_mmm_active_elements",
         "hdiv_mmm_block_insertion_response",
         "hdiv_mmm_all_single_removal_responses",
         "hdiv_mmm_removal_group_responses",
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
         "production_vim_state_shape_jacobian_streaming",
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
