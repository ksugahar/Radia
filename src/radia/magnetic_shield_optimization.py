"""Thin-sheet magnetic-shield application adapter for Radia-VIM optimization."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .sheet_metal_optimization import solve_sheet_metal_lp
from .topology_optimization import (VIMOperatorLinearization,
    assemble_ngsolve_hdiv_shape_tangents, linearize_laplace_charge_gram,
    linearize_production_vim_from_ngsolve, linearize_vim_operator,
    linearize_vim_system, production_vim_rms_adjoint_gradient_streaming)


@dataclass(frozen=True)
class MagneticShieldDesign:
    width: float=0.1
    height: float=0.1
    thickness: float=0.002
    mu_r: float=1000.0
    nx: int=4
    ny: int=4
    nz: int=1
    applied_h: tuple[float,float,float]=(1.0e4,0.0,0.0)
    probes: tuple[tuple[float,float,float],...]=((0.0,0.0,0.01),(0.0,0.0,0.02),(0.02,0.0,0.01))

    def __post_init__(self):
        if min(self.width,self.height,self.thickness,self.mu_r)<=0: raise ValueError("shield dimensions and mu_r must be positive")
        if min(self.nx,self.ny,self.nz)<1 or not self.probes: raise ValueError("mesh subdivisions and probes must be non-empty")


@dataclass(frozen=True)
class MagneticShieldResult:
    mesh: object
    vim_result: dict
    total_h: np.ndarray
    rms_parallel_h: float
    shielding_factor: float


@dataclass(frozen=True)
class MagneticShieldLPStep:
    charge_gram: object
    operator: VIMOperatorLinearization
    response: object
    rms_gradient: np.ndarray
    update: object


@dataclass(frozen=True)
class StreamingMagneticShieldLPStep:
    adjoint: object
    update: object


def build_plate_mesh(design: MagneticShieldDesign):
    """Build the initial structured HEX neutral-sheet volume."""
    from ngsolve.meshes import MakeStructured3DMesh
    return MakeStructured3DMesh(hexes=True,nx=design.nx,ny=design.ny,nz=design.nz,
        mapping=lambda X,Y,Z:(design.width*(X-.5),design.height*(Y-.5),design.thickness*(Z-.5)))


def solve_shield(design: MagneticShieldDesign, *, order=1, tolerance=1e-9):
    """Solve one plate and measure the total field at protected probes.

    The caller owns ``ngsolve.TaskManager`` in accordance with the repository
    FE policy.
    """
    import ngsolve as ng
    import radia.vim as vim
    mesh=build_plate_mesh(design); applied=np.asarray(design.applied_h,dtype=float)
    result=vim.Solve(mesh,mu_r=design.mu_r,H_ext=ng.CF(tuple(applied)),order=order,tol=tolerance)
    scattered=np.asarray(vim.FieldFromSolution(result,np.asarray(design.probes),algorithm="direct"),dtype=float)
    total=scattered+applied
    direction=applied/np.linalg.norm(applied)
    parallel=total@direction; rms=float(np.sqrt(np.mean(parallel**2)))
    return MagneticShieldResult(mesh,result,total,rms,float(np.linalg.norm(applied)/rms))


def rms_response_gradient(response, response_jacobian):
    """Map an analytic VIM probe-field Jacobian to the RMS objective gradient."""
    values=np.asarray(response,dtype=float).reshape(-1)
    jac=np.asarray(response_jacobian,dtype=float)
    if jac.ndim!=2 or jac.shape[0]!=values.size or values.size==0: raise ValueError("response/Jacobian shape mismatch")
    rms=float(np.sqrt(np.mean(values**2)))
    if rms==0: raise ValueError("RMS response derivative is undefined at zero")
    return (values@jac)/(values.size*rms)


def linearize_shield_response(operator: VIMOperatorLinearization, response_matrix,
                              *, dresponse_matrix=None, incident_response=None):
    """Close the analytic VIM operator tangent into protected-zone responses."""
    C=np.asarray(response_matrix)
    model=linearize_vim_system(operator.matrix,operator.rhs,C,operator.matrix_jacobian,
        db=operator.rhs_jacobian,dC=dresponse_matrix)
    if incident_response is None: return model
    incident=np.asarray(incident_response).reshape(-1)
    if incident.shape!=model.response.shape: raise ValueError("incident_response shape mismatch")
    return type(model)(model.state,model.response+incident,model.state_jacobian,model.response_jacobian)


def sheet_lp_step(normal_displacement, thickness, activation, cell_areas, *,
                  parallel_field, analytic_response_jacobian, volume_max,
                  displacement_move, thickness_move, thickness_bounds,
                  laplacian=None, curvature_limit=None):
    """Perform one no-finite-difference magnetic-shield sheet LP update."""
    gradient=rms_response_gradient(parallel_field,analytic_response_jacobian)
    return solve_sheet_metal_lp(normal_displacement,thickness,activation,gradient,cell_areas,
        volume_max=volume_max,displacement_move=displacement_move,thickness_move=thickness_move,
        thickness_bounds=thickness_bounds,laplacian=laplacian,curvature_limit=curvature_limit)


def linearize_and_step_shield(*, points, weights, displacement_modes,
                              relative_weight_derivatives, self_cell_types,
                              self_nodes, self_node_displacements,
                              mass, charge_map, applied_coefficients, inv_chi,
                              dmass, dcharge_map, response_matrix,
                              normal_displacement, thickness, activation,
                              cell_areas, volume_max, displacement_move,
                              thickness_move, thickness_bounds,
                              dresponse_matrix=None, incident_response=None,
                              laplacian=None, curvature_limit=None):
    """Close GetTrafo geometry tangents through the VIM response and sheet LP.

    ``dmass`` and ``dcharge_map`` must come from NGSolve assembly.  The only
    singular operation performed here is the C++ analytic self-panel kernel;
    no finite difference is used by this optimization path.
    """
    charge=linearize_laplace_charge_gram(points,weights,displacement_modes,
        relative_weight_derivatives=relative_weight_derivatives,
        self_cell_types=self_cell_types,self_nodes=self_nodes,
        self_node_displacements=self_node_displacements)
    operator=linearize_vim_operator(mass,charge_map,charge.matrix,applied_coefficients,
        inv_chi=inv_chi,dmass=dmass,dcharge_map=dcharge_map,dcharge_gram=charge.jacobian)
    response=linearize_shield_response(operator,response_matrix,
        dresponse_matrix=dresponse_matrix,incident_response=incident_response)
    gradient=rms_response_gradient(response.response,response.response_jacobian)
    update=solve_sheet_metal_lp(normal_displacement,thickness,activation,gradient,cell_areas,
        volume_max=volume_max,displacement_move=displacement_move,
        thickness_move=thickness_move,thickness_bounds=thickness_bounds,
        laplacian=laplacian,curvature_limit=curvature_limit)
    return MagneticShieldLPStep(charge,operator,response,gradient,update)


def linearize_and_step_shield_from_ngsolve(*, fes, deformation_modes,
                                           charge_map, weights,
                                           relative_weight_derivatives, **kwargs):
    """NGSolve-native entry from GetTrafo modes through the sheet LP.

    Geometry corners and mode values come from ``GetTrafo``; HDiv ``M,dM``
    come from NGSolve weak-form assembly; the Piola-exact reference charge map
    has ``dB=0``.  The caller owns ``ngsolve.TaskManager``.
    """
    from .sheet_metal_optimization import sample_affine_gettrafo_cells
    sampled=sample_affine_gettrafo_cells(fes.mesh,deformation_modes)
    if len(getattr(charge_map,"shape",()))!=2 or charge_map.shape[0]!=len(sampled.nodes):
        raise NotImplementedError(
            "the current affine Python closure supports one constant volume-charge row per cell; "
            "full BDM polynomial cell/face self blocks are not yet wired")
    mass,dmass,dcharge=assemble_ngsolve_hdiv_shape_tangents(
        fes,deformation_modes,charge_map)
    return linearize_and_step_shield(points=sampled.centroids,weights=weights,
        displacement_modes=sampled.centroid_displacements,
        relative_weight_derivatives=relative_weight_derivatives,
        self_cell_types=sampled.cell_types,self_nodes=sampled.nodes,
        self_node_displacements=sampled.node_displacements,
        mass=mass,charge_map=charge_map,dmass=dmass,dcharge_map=dcharge,**kwargs)


def linearize_and_step_production_shield_from_ngsolve(*, fes, deformation_modes,
        charge_basis, charge_gram, charge_map, applied_coefficients, inv_chi,
        response_matrix, normal_displacement, thickness, activation, cell_areas,
        volume_max, displacement_move, thickness_move, thickness_bounds,
        dapplied_coefficients=None, dresponse_matrix=None,
        incident_response=None, laplacian=None, curvature_limit=None,
        family=None):
    """Production TET/HEX/WEDGE GetTrafo-to-sheet-LP closure.

    Unlike the legacy constant-cell adapter above, this path differentiates
    every polynomial volume/face ChargeGram block, including non-self pairs.
    All geometry derivatives are analytic; the caller owns TaskManager.
    """
    production=linearize_production_vim_from_ngsolve(fes=fes,
        deformation_modes=deformation_modes,charge_basis=charge_basis,
        charge_gram=charge_gram,charge_map=charge_map,
        applied_coefficients=applied_coefficients,inv_chi=inv_chi,
        dapplied_coefficients=dapplied_coefficients,family=family)
    response=linearize_shield_response(production.operator,response_matrix,
        dresponse_matrix=dresponse_matrix,incident_response=incident_response)
    gradient=rms_response_gradient(response.response,response.response_jacobian)
    update=solve_sheet_metal_lp(normal_displacement,thickness,activation,gradient,
        cell_areas,volume_max=volume_max,displacement_move=displacement_move,
        thickness_move=thickness_move,thickness_bounds=thickness_bounds,
        laplacian=laplacian,curvature_limit=curvature_limit)
    return MagneticShieldLPStep(production.charge_gram,production.operator,
        response,gradient,update)


def linearize_and_step_production_shield_streaming(*, fes, deformation_modes,
        charge_basis, charge_gram, charge_map, applied_coefficients, inv_chi,
        response_matrix, normal_displacement, thickness, activation, cell_areas,
        volume_max, displacement_move, thickness_move, thickness_bounds,
        dapplied_coefficients=None, dresponse_matrix=None,
        incident_response=None, laplacian=None, curvature_limit=None,
        family=None, eps=1e-10, leaf=64, eta=2.0,
        solve_tolerance=1e-10, solve_max_iterations=None):
    """Large-design-variable shield LP step using streaming adjoints."""
    adjoint=production_vim_rms_adjoint_gradient_streaming(
        fes=fes,deformation_modes=deformation_modes,charge_basis=charge_basis,
        charge_gram=charge_gram,charge_map=charge_map,
        applied_coefficients=applied_coefficients,inv_chi=inv_chi,
        response_matrix=response_matrix,family=family,
        incident_response=incident_response,
        dapplied_coefficients=dapplied_coefficients,
        dresponse_matrix=dresponse_matrix,eps=eps,leaf=leaf,eta=eta,
        solve_tolerance=solve_tolerance,
        solve_max_iterations=solve_max_iterations)
    update=solve_sheet_metal_lp(normal_displacement,thickness,activation,
        adjoint.gradient,cell_areas,volume_max=volume_max,
        displacement_move=displacement_move,thickness_move=thickness_move,
        thickness_bounds=thickness_bounds,laplacian=laplacian,
        curvature_limit=curvature_limit)
    return StreamingMagneticShieldLPStep(adjoint,update)


__all__=["MagneticShieldDesign","MagneticShieldResult","MagneticShieldLPStep",
         "StreamingMagneticShieldLPStep","build_plate_mesh","solve_shield",
         "rms_response_gradient","linearize_shield_response","sheet_lp_step",
         "linearize_and_step_shield","linearize_and_step_shield_from_ngsolve",
         "linearize_and_step_production_shield_from_ngsolve",
         "linearize_and_step_production_shield_streaming"]
