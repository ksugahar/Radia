"""Per-element density adjoint for HDiv-MMM topology optimization.

Application-layer module for isochronous-magnet (and general accelerator
pole/yoke) density topology optimization on the HDiv-MMM forward engine
(design record: ``docs/hdiv_vim/TOPOLOGY_OPTIMIZATION_ISOCHRONOUS.md``).
It composes :class:`radia.vim.DemagOperator` (build-once geometry operator
``N = B^T G B``) with the per-element design variable ``s_e = 1/chi_e``
(an ``L2(order=0)`` mass weight) into the adjoint gradient route:

    (M_s + N) m      = f_state        state solve (native H-matrix Jacobi-PCG)
    (M_s + N) lambda = f_adjoint      adjoint (same operator, same factor)
    J(s)             = f_adjoint^T m  (the design-dependent objective part)
    dJ/ds_e          = -lambda^T M^(e) m
                     = -Integrate(lambda . m, element_wise=True)[e]

The last identity is exact because ``A(s) = sum_e s_e M^(e) + N`` is affine in
``s`` and the element mass contraction of two HDiv fields IS the element-wise
integral of their inner product -- all sensitivities come from ONE
``Integrate`` call, with no per-mode machinery and no dense Jacobian.

The field functional ``F^T c`` (weighted point evaluations of the
demagnetizing field, e.g. orbit-quadrature ``dB_z/dx`` for the focusing
integral ``gL``) is realized by kernel reciprocity: the adjoint load is the
mass projection of the analytic H field of a point-DIPOLE array at the
evaluation points.  The dense magnetization-to-field matrix ``F`` is never
formed; derivative functionals enter as +/- point pairs
(:func:`gradient_pair_points`).

Verified 2026-07-28 on the promoted unit-ball gate (270 tets, log-uniform
``s`` in [1e-2, 1]): adjoint gradient matches
central finite differences to 8.1e-10 (directional, all elements) and 3.9e-7
worst-case per element (the smallest-gradient element, FD noise floor); the
reciprocity load matches the independent C++ analytic charge evaluator
(``vim.FieldFromSolution``) to 1.1e-10 at ``bonus_intorder=10``.  Locked by
``tests/test_isochronous_topopt.py``.

The caller wraps all NGSolve work (construction and solves) in
``with TaskManager():`` per the repository TaskManager policy.  The LP/OC
design-update drivers live in :mod:`radia.topology_optimization`; this module
supplies the objective/gradient callback layer above the same operators.
"""
import time
from dataclasses import dataclass
from math import pi

import numpy as np
import ngsolve as ng
from ngsolve.krylovspace import CGSolver

from .vim import DemagOperator, FieldFromSolution

MU0 = 4.0e-7 * pi
#: Validated ersatz-void susceptibility floor: with the s-weighted mass-Riesz
#: preconditioner, CG iteration counts DECREASE as chi_void -> 0 (measured
#: 74 -> 46 from 1e-1 to 1e-6), and the void response stays exactly the
#: physical chi*H over five decades.  See the design document, Sec. 3.
CHI_MIN = 1.0e-6

__all__ = [
    "MU0", "CHI_MIN", "AdjointGradientResult", "DensityAdjointVIM",
    "FunctionalLinearization", "DensityDesignResult", "HelmholtzFilter",
    "HeavisideProjection", "density_discreteness",
    "iron_only_verification_ready",
    "IronOnlyVerification", "density_to_s",
    "density_gradient_from_s_gradient", "gradient_pair_points",
    "dipole_array_field_cf", "field_functional_load", "uniform_field_load",
    "demag_field_from_solution", "demag_field_evaluator",
    "orbit_arc_points", "optimize_density",
    "SectorOrbit", "SectorLinearOptics", "track_sector_orbit",
    "sector_linear_optics", "isochronous_increment_targets",
    "isochronous_total_field_bands", "isochronous_profile_metrics",
    "CombinedFunctionLinearOptics", "combined_function_linear_optics",
    "CombinedFunctionTransferMap", "combined_function_transfer_map",
    "combined_function_transfer_map_from_field_response",
    "TransferMapReachability", "transfer_map_reachability",
    "TransferMapTargetDesign", "design_combined_function_transfer_map_target",
    "CombinedFunctionExitMetrics", "combined_function_exit_metrics",
    "combined_function_exit_metrics_from_field_response",
    "StraightenedBendValidation", "straightened_bend_validation",
    "AchromaticGradientDesign", "design_achromatic_gradient_profile",
    "iron_only_mesh", "verify_design_iron_only",
]


# --------------------------------------------------------------------------
# density <-> s mapping
# --------------------------------------------------------------------------
def density_to_s(density, chi_iron, chi_min=CHI_MIN, penalty=1.0):
    """Map a per-element density ``rho in [0, 1]`` to ``s = 1/chi(rho)``.

    SIMP-style interpolation in susceptibility, ``chi(rho) = chi_min +
    (chi_iron - chi_min) rho^penalty``.  ``penalty=1`` (default) is the
    linear map; ``penalty=3`` makes intermediate densities material-
    inefficient and drives designs toward 0/1 (use for runs whose result
    will be thresholded/manufactured).  The filtering layer of the design
    loop operates on ``rho`` before this map.  ``chi_min`` defaults to the
    validated ersatz floor :data:`CHI_MIN`.
    """
    rho = np.asarray(density, dtype=float)
    chi_iron = float(chi_iron)
    chi_min = float(chi_min)
    penalty = float(penalty)
    if not (np.isfinite(chi_iron) and np.isfinite(chi_min)
            and chi_iron > chi_min > 0.0):
        raise ValueError(
            "density_to_s: need chi_iron > chi_min > 0 (got chi_iron=%r, chi_min=%r)"
            % (chi_iron, chi_min))
    if not np.isfinite(penalty) or penalty < 1.0:
        raise ValueError("density_to_s: penalty must be >= 1")
    if not np.all(np.isfinite(rho)):
        raise ValueError("density_to_s: density must be finite")
    if np.any(rho < -1e-9) or np.any(rho > 1.0 + 1e-9):
        raise ValueError(
            "density_to_s: density must lie in [0, 1] (got min=%r, max=%r)"
            % (float(rho.min()), float(rho.max())))
    chi = chi_min + (chi_iron - chi_min) * np.clip(rho, 0.0, 1.0) ** penalty
    return 1.0 / chi


def density_gradient_from_s_gradient(density, s_gradient, chi_iron,
                                     chi_min=CHI_MIN, penalty=1.0):
    """Chain rule ``dJ/drho_e = dJ/ds_e * ds/drho_e`` for :func:`density_to_s`.

    ``ds/drho = -(chi_iron - chi_min) penalty rho^(penalty-1) / chi(rho)^2``.
    """
    rho = np.asarray(density, dtype=float)
    grad_s = np.asarray(s_gradient, dtype=float)
    if rho.shape != grad_s.shape:
        raise ValueError("density_gradient_from_s_gradient: shape mismatch %r vs %r"
                         % (rho.shape, grad_s.shape))
    chi_iron = float(chi_iron)
    chi_min = float(chi_min)
    penalty = float(penalty)
    if not (np.isfinite(chi_iron) and np.isfinite(chi_min)
            and chi_iron > chi_min > 0.0):
        raise ValueError("density_gradient_from_s_gradient: need "
                         "chi_iron > chi_min > 0")
    if not np.isfinite(penalty) or penalty < 1.0:
        raise ValueError("density_gradient_from_s_gradient: penalty must be >= 1")
    if (not np.all(np.isfinite(rho)) or not np.all(np.isfinite(grad_s))
            or np.any(rho < 0.0) or np.any(rho > 1.0)):
        raise ValueError("density_gradient_from_s_gradient: density and "
                         "gradient must be finite, with density in [0, 1]")
    chi = chi_min + (chi_iron - chi_min) * rho ** penalty
    dchi = ((chi_iron - chi_min) * penalty * rho ** (penalty - 1.0)
            if penalty != 1.0 else (chi_iron - chi_min) * np.ones_like(rho))
    return grad_s * (-dchi / (chi * chi))


# --------------------------------------------------------------------------
# field-functional loads (F^T c by kernel reciprocity)
# --------------------------------------------------------------------------
def gradient_pair_points(points, weights, delta, axis=0, direction=None):
    """Realize ``sum_i w_i d(field)/d(direction) (x_i)`` as +/- evaluation pairs.

    Returns ``(points_out, weights_out)`` with ``2 N`` rows: each input point
    splits into ``x_i +/- (delta/2) d_i`` carrying weights ``+/- w_i/delta``
    (central finite difference of the point functional; ``delta`` is a
    physical stencil length on the orbit, not a numerical epsilon).
    ``direction`` is a single (3,) vector or per-point (N, 3) vectors
    (normalized here); when omitted the coordinate ``axis`` is used --
    per-point radial directions realize the ``d B_z/d r`` orbit functionals.
    """
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    wts = np.asarray(weights, dtype=float).ravel()
    if len(pts) != len(wts):
        raise ValueError("gradient_pair_points: %d points but %d weights"
                         % (len(pts), len(wts)))
    if not delta > 0.0:
        raise ValueError("gradient_pair_points: delta must be positive")
    if direction is None:
        if axis not in (0, 1, 2):
            raise ValueError("gradient_pair_points: axis must be 0, 1, or 2")
        direction = np.zeros(3)
        direction[axis] = 1.0
    d = np.asarray(direction, dtype=float)
    d = np.broadcast_to(d.reshape(-1, 3), pts.shape).copy()
    norms = np.linalg.norm(d, axis=1)
    if np.any(norms <= 0.0):
        raise ValueError("gradient_pair_points: zero direction vector")
    d /= norms[:, None]
    shift = 0.5 * float(delta) * d
    return (np.concatenate([pts + shift, pts - shift]),
            np.concatenate([wts / delta, -wts / delta]))


def orbit_arc_points(radius, z, n, span=(0.0, 2.0 * pi)):
    """Equispaced points on a circular arc (radius, height ``z``), with the
    per-point radial unit vectors for :func:`gradient_pair_points`.

    Returns ``(points (n,3), radial_directions (n,3))``.  A full circle
    (span of 2 pi) omits the duplicate endpoint.
    """
    lo, hi = float(span[0]), float(span[1])
    full = abs((hi - lo) - 2.0 * pi) < 1e-12
    theta = np.linspace(lo, hi, int(n), endpoint=not full)
    pts = np.stack([radius * np.cos(theta), radius * np.sin(theta),
                    float(z) * np.ones_like(theta)], axis=1)
    radial = np.stack([np.cos(theta), np.sin(theta),
                       np.zeros_like(theta)], axis=1)
    return pts, radial


@dataclass
class SectorOrbit:
    """RK4 reference trajectory through one azimuthal magnet sector."""
    position: np.ndarray
    tangent: np.ndarray
    path_length: float
    exit_radius: float
    exit_angle: float


@dataclass
class SectorLinearOptics:
    """Tracked radial/vertical first-order maps about a reference orbit."""
    reference: SectorOrbit
    radial_matrix: np.ndarray
    vertical_matrix: np.ndarray
    radial_determinant: float
    vertical_determinant: float


@dataclass(frozen=True)
class CombinedFunctionLinearOptics:
    """Piecewise-constant combined-function map and analytic sensitivities.

    ``curvature`` is ``h=B/B_rho`` and ``normalized_gradient`` is
    ``k1=(dB_y/dx)/B_rho``.  The adopted accelerator convention is

    ``eta'' + (h**2+k1) eta = h`` and ``y'' - k1 y = 0``.

    The endpoint Jacobian differentiates with respect to caller-supplied
    design parameters through exact matrix-exponential Frechet derivatives;
    no momentum or design finite difference is used.
    """
    s: np.ndarray
    dispersion: np.ndarray
    radial_matrix: np.ndarray
    vertical_matrix: np.ndarray
    endpoint_jacobian: np.ndarray
    radial_matrix_jacobian: np.ndarray
    vertical_matrix_jacobian: np.ndarray
    radial_trace: float
    vertical_trace: float
    radial_stable: bool
    vertical_stable: bool


@dataclass(frozen=True)
class CombinedFunctionTransferMap:
    """First-order bend map and analytic field/design sensitivities.

    The state ordering is ``(x, x', y, y', ell, delta)``.  The returned
    6-by-6 matrix contains the horizontal and vertical betatron blocks, the
    horizontal dispersion column, and the geometric path-length row.
    ``matrix_jacobian[p]`` is the derivative of the complete matrix with
    respect to design parameter ``p``.  The magnetic model has no RF or
    velocity-slip input, so ``R56`` is the geometric contribution only.

    ``response`` packs the entries requested through ``response_entries``;
    this is the row-major numerical contract consumed by the HDiv-MMM
    add/remove LP.  The default entries are the two complete 2-by-2 blocks and
    ``R16, R26, R51, R52, R56``.  No particle or design finite difference is
    used.
    """
    matrix: np.ndarray
    matrix_jacobian: np.ndarray
    response: np.ndarray
    response_jacobian: np.ndarray
    response_entries: tuple[tuple[int, int], ...]
    optics: CombinedFunctionLinearOptics


@dataclass(frozen=True)
class TransferMapReachability:
    """TSVD reachability certificate for one linearized topology step."""
    numerical_rank: int
    singular_values: np.ndarray
    parameter_step: np.ndarray
    predicted_response: np.ndarray
    residual: np.ndarray
    max_normalized_residual: float
    reachable: bool


@dataclass(frozen=True)
class TransferMapTargetDesign:
    """Ideal-optics target subsequently realized by HDiv-MMM topology."""
    matrix: np.ndarray
    segment_lengths: np.ndarray
    curvature: np.ndarray
    normalized_gradient: np.ndarray
    transfer_map: CombinedFunctionTransferMap
    maximum_scaled_residual: float
    iterations: int
    status: str


@dataclass(frozen=True)
class CombinedFunctionExitMetrics:
    """Four exit responses used by the five-particle acceptance contract.

    ``reference_orbit_error`` is the central-particle displacement/slope from
    the prescribed design orbit.  ``response`` is ordered as
    ``(x0, psi0, eta, eta_prime)`` and ``response_jacobian`` has one column
    per caller-supplied field/design parameter.  The Jacobian is assembled
    from matrix-exponential Frechet derivatives; finite differences are not
    part of the optimization path.
    """
    s: np.ndarray
    reference_orbit_error: np.ndarray
    x0_m: float
    psi0_rad: float
    eta_m: float
    eta_prime_rad: float
    response: np.ndarray
    response_jacobian: np.ndarray
    downstream_drift_m: float
    optics: CombinedFunctionLinearOptics


@dataclass(frozen=True)
class StraightenedBendValidation:
    """Field-driven reference orbit and optics in a straightened bend frame.

    The electromagnetic mesh uses longitudinal coordinate ``s`` while the
    returned Cartesian orbit integrates the realized curvature into the
    physical bend.  This is the standard straightened Frenet representation:
    it avoids pretending that the straight HEX mesh itself is a curved global
    geometry while still checking the field, bend angle, transfer maps, and
    endpoint dispersion produced by that mesh.
    """
    s: np.ndarray
    bz: np.ndarray
    gradient: np.ndarray
    bend_angle: float
    position: np.ndarray
    tangent: np.ndarray
    optics: CombinedFunctionLinearOptics


@dataclass(frozen=True)
class AchromaticGradientDesign:
    """Stable longitudinal gradient target with zero endpoint dispersion."""
    segment_lengths: np.ndarray
    curvature: np.ndarray
    normalized_gradient: np.ndarray
    optics: CombinedFunctionLinearOptics
    objective: float
    iterations: int
    status: str


def _orbit_field_value(field, point):
    value = np.asarray(field(np.asarray(point, dtype=float)), dtype=float)
    if value.shape == (1, 3):
        value = value[0]
    if value.shape != (3,) or not np.all(np.isfinite(value)):
        raise ValueError("track_sector_orbit: field must return one finite 3-vector")
    return value


def track_sector_orbit(field, radius, span, inverse_rigidity, *, n_steps=720,
                       radial_offset=0.0, vertical_offset=0.0,
                       radial_slope=0.0, vertical_slope=0.0):
    """Track ``dr/ds=u, du/ds=(q/p) u x B`` through a nominal circular arc.

    ``field(point)`` returns magnetic flux density in tesla and
    ``inverse_rigidity=q/p`` in inverse tesla-metre, including its sign.
    The integration length is ``radius * abs(span[1]-span[0])``.  Offsets and
    slopes are validation perturbations for :func:`sector_linear_optics`; they
    are not topology sensitivities and never enter the optimization gradient.
    """
    radius = float(radius)
    inverse_rigidity = float(inverse_rigidity)
    n_steps = int(n_steps)
    lo, hi = map(float, span)
    if not radius > 0.0 or not np.isfinite(inverse_rigidity) or n_steps < 4:
        raise ValueError("track_sector_orbit: need radius>0, finite rigidity, n_steps>=4")
    direction = 1.0 if hi >= lo else -1.0
    er = np.array([np.cos(lo), np.sin(lo), 0.0])
    et = direction * np.array([-np.sin(lo), np.cos(lo), 0.0])
    ez = np.array([0.0, 0.0, 1.0])
    position = (radius + float(radial_offset)) * er + float(vertical_offset) * ez
    tangent = et + float(radial_slope) * er + float(vertical_slope) * ez
    tangent /= np.linalg.norm(tangent)
    length = radius * abs(hi - lo)
    ds = length / n_steps
    positions = np.empty((n_steps + 1, 3), dtype=float)
    tangents = np.empty_like(positions)
    positions[0], tangents[0] = position, tangent

    def rhs(r, u):
        return u, inverse_rigidity * np.cross(u, _orbit_field_value(field, r))

    for step in range(n_steps):
        k1r, k1u = rhs(position, tangent)
        k2r, k2u = rhs(position + 0.5*ds*k1r,
                       tangent + 0.5*ds*k1u)
        k3r, k3u = rhs(position + 0.5*ds*k2r,
                       tangent + 0.5*ds*k2u)
        k4r, k4u = rhs(position + ds*k3r, tangent + ds*k3u)
        position = position + ds*(k1r + 2*k2r + 2*k3r + k4r)/6.0
        tangent = tangent + ds*(k1u + 2*k2u + 2*k3u + k4u)/6.0
        tangent /= np.linalg.norm(tangent)
        positions[step+1], tangents[step+1] = position, tangent
    return SectorOrbit(
        position=positions, tangent=tangents, path_length=length,
        exit_radius=float(np.hypot(position[0], position[1])),
        exit_angle=float(np.arctan2(position[1], position[0])))


def sector_linear_optics(field, radius, span, inverse_rigidity, *,
                         n_steps=720, position_step=1e-5,
                         slope_step=1e-5):
    """Track symmetric orbit pencils and return radial/vertical 2x2 maps."""
    if position_step <= 0.0 or slope_step <= 0.0:
        raise ValueError("sector_linear_optics: perturbation steps must be positive")
    kwargs = dict(n_steps=n_steps)
    reference = track_sector_orbit(field, radius, span, inverse_rigidity, **kwargs)
    end = reference.position[-1]
    er = np.array([end[0], end[1], 0.0])
    er /= np.linalg.norm(er)
    et = np.array([-er[1], er[0], 0.0])

    def difference(kind, step):
        plus = track_sector_orbit(
            field, radius, span, inverse_rigidity, **kwargs, **{kind: step})
        minus = track_sector_orbit(
            field, radius, span, inverse_rigidity, **kwargs, **{kind: -step})
        def exit_plane_state(orbit):
            r = orbit.position[-1].copy()
            u = orbit.tangent[-1].copy()
            correction = -((r-end)@et)/(u@et)
            r += correction*u
            u += correction*inverse_rigidity*np.cross(
                u, _orbit_field_value(field, r))
            u /= np.linalg.norm(u)
            return r, u
        rp, up = exit_plane_state(plus)
        rm, um = exit_plane_state(minus)
        dr = (rp-rm)/(2*step)
        du = (up-um)/(2*step)
        return dr, du

    drx, dux = difference("radial_offset", position_step)
    drxp, duxp = difference("radial_slope", slope_step)
    drz, duz = difference("vertical_offset", position_step)
    drzp, duzp = difference("vertical_slope", slope_step)
    radial = np.array([[drx@er, drxp@er], [dux@er, duxp@er]])
    vertical = np.array([[drz[2], drzp[2]], [duz[2], duzp[2]]])
    return SectorLinearOptics(
        reference=reference, radial_matrix=radial, vertical_matrix=vertical,
        radial_determinant=float(np.linalg.det(radial)),
        vertical_determinant=float(np.linalg.det(vertical)))


def combined_function_linear_optics(curvature, normalized_gradient,
                                    segment_lengths, *, eta0=0.0,
                                    eta_prime0=0.0,
                                    curvature_jacobian=None,
                                    gradient_jacobian=None,
                                    stability_tolerance=1e-10):
    """Propagate dispersion and betatron maps through a combined-function bend.

    Inputs are one constant value per longitudinal segment.  Optional
    Jacobians have shape ``(n_segment, n_parameter)``.  The propagators and all
    requested derivatives are evaluated analytically with
    :func:`scipy.linalg.expm_frechet`.  This provides the optics link used to
    chain HDiv-MMM field/gradient response Jacobians into an endpoint-
    dispersion constraint without finite differences in the optimizer.
    """
    from scipy.linalg import expm, expm_frechet

    h=np.asarray(curvature,dtype=float).reshape(-1)
    k=np.asarray(normalized_gradient,dtype=float).reshape(-1)
    ds=np.asarray(segment_lengths,dtype=float).reshape(-1)
    if h.size==0 or k.shape!=h.shape or ds.shape!=h.shape:
        raise ValueError("combined-function optics requires matching non-empty segment arrays")
    if (np.any(ds<=0.0) or not np.all(np.isfinite(np.r_[h,k,ds,eta0,eta_prime0]))):
        raise ValueError("combined-function optics inputs must be finite and lengths positive")

    supplied=(curvature_jacobian is not None or gradient_jacobian is not None)
    if supplied:
        raw=(curvature_jacobian if curvature_jacobian is not None
             else gradient_jacobian)
        raw=np.asarray(raw,dtype=float)
        if raw.ndim!=2 or raw.shape[0]!=h.size:
            raise ValueError("combined-function Jacobians need shape (n_segment,n_parameter)")
        n_parameter=raw.shape[1]
        dh=(np.zeros((h.size,n_parameter)) if curvature_jacobian is None else
            np.asarray(curvature_jacobian,dtype=float))
        dk=(np.zeros((h.size,n_parameter)) if gradient_jacobian is None else
            np.asarray(gradient_jacobian,dtype=float))
        if dh.shape!=(h.size,n_parameter) or dk.shape!=(h.size,n_parameter):
            raise ValueError("curvature and gradient Jacobians must have matching shapes")
        if not np.all(np.isfinite(np.r_[dh.ravel(),dk.ravel()])):
            raise ValueError("combined-function Jacobians must be finite")
    else:
        n_parameter=0
        dh=np.zeros((h.size,0));dk=np.zeros((h.size,0))

    z=np.array([float(eta0),float(eta_prime0),1.0])
    dz=np.zeros((3,n_parameter))
    radial=np.eye(2); vertical=np.eye(2)
    dradial=np.zeros((n_parameter,2,2));dvertical=np.zeros_like(dradial)
    history=[z[:2].copy()];stations=[0.0]
    for segment,(hi,ki,length) in enumerate(zip(h,k,ds)):
        kx=hi*hi+ki
        generator=np.array([[0.0,1.0,0.0],[-kx,0.0,hi],[0.0,0.0,0.0]])
        scaled=generator*length
        propagator=expm(scaled)
        old_z=z; old_dz=dz; old_radial=radial;old_dradial=dradial
        old_vertical=vertical;old_dvertical=dvertical
        dprop=[];dvertical_prop=[]
        vertical_generator=np.array([[0.0,1.0],[ki,0.0]])
        vertical_propagator=expm(vertical_generator*length)
        for parameter in range(n_parameter):
            dhi=dh[segment,parameter];dki=dk[segment,parameter]
            if dhi == 0.0 and dki == 0.0:
                dprop.append(np.zeros((3,3)))
                dvertical_prop.append(np.zeros((2,2)))
                continue
            derivative=np.zeros((3,3))
            derivative[1,0]=-(2.0*hi*dhi+dki)
            derivative[1,2]=dhi
            dprop.append(expm_frechet(
                scaled,derivative*length,compute_expm=False))
            vertical_derivative=np.array([[0.0,0.0],[dki,0.0]])
            dvertical_prop.append(expm_frechet(
                vertical_generator*length,vertical_derivative*length,
                compute_expm=False))
        z=propagator@old_z
        radial=propagator[:2,:2]@old_radial
        vertical=vertical_propagator@old_vertical
        for parameter in range(n_parameter):
            dz[:,parameter]=(propagator@old_dz[:,parameter]
                             +dprop[parameter]@old_z)
            dradial[parameter]=(propagator[:2,:2]@old_dradial[parameter]
                                +dprop[parameter][:2,:2]@old_radial)
            dvertical[parameter]=(vertical_propagator@old_dvertical[parameter]
                                  +dvertical_prop[parameter]@old_vertical)
        history.append(z[:2].copy());stations.append(stations[-1]+length)
    radial_trace=float(np.trace(radial));vertical_trace=float(np.trace(vertical))
    tolerance=float(stability_tolerance)
    if not np.isfinite(tolerance) or tolerance<0.0:
        raise ValueError("stability_tolerance must be finite and nonnegative")
    return CombinedFunctionLinearOptics(
        np.asarray(stations),np.asarray(history),radial,vertical,dz[:2],
        dradial,dvertical,radial_trace,vertical_trace,
        bool(abs(radial_trace)<2.0-tolerance),
        bool(abs(vertical_trace)<2.0-tolerance))


_DEFAULT_TRANSFER_RESPONSE_ENTRIES = (
    (0, 0), (0, 1), (0, 5),
    (1, 0), (1, 1), (1, 5),
    (2, 2), (2, 3),
    (3, 2), (3, 3),
    (4, 0), (4, 1), (4, 5),
)


def combined_function_transfer_map(
        curvature, normalized_gradient, segment_lengths, *,
        curvature_jacobian=None, gradient_jacobian=None,
        response_entries=None, stability_tolerance=1e-10):
    """Return the 6-by-6 first-order map of a combined-function bend.

    The horizontal-longitudinal generator for one constant segment is

    ``x' = px``, ``px' = -(h**2+k1)x + h*delta``,
    ``ell' = h*x``, ``delta' = 0``.

    Exact matrix exponentials are multiplied in traversal order.  Their
    derivatives use ``expm_frechet`` and therefore form an analytic chain from
    HDiv-MMM field/gradient response rows to every selected transfer-matrix
    entry.  ``response_entries`` uses zero-based ``(row, column)`` pairs in the
    standard state order ``(x,x',y,y',ell,delta)``.
    """
    from scipy.linalg import expm, expm_frechet

    h = np.asarray(curvature, dtype=float).reshape(-1)
    k = np.asarray(normalized_gradient, dtype=float).reshape(-1)
    ds = np.asarray(segment_lengths, dtype=float).reshape(-1)
    if h.size == 0 or k.shape != h.shape or ds.shape != h.shape:
        raise ValueError(
            "combined-function transfer map requires matching non-empty "
            "segment arrays")
    if np.any(ds <= 0.0) or not np.all(np.isfinite(np.r_[h, k, ds])):
        raise ValueError(
            "combined-function transfer-map inputs must be finite and "
            "lengths positive")

    supplied = (curvature_jacobian is not None
                or gradient_jacobian is not None)
    if supplied:
        raw = (curvature_jacobian if curvature_jacobian is not None
               else gradient_jacobian)
        raw = np.asarray(raw, dtype=float)
        if raw.ndim != 2 or raw.shape[0] != h.size:
            raise ValueError(
                "combined-function transfer-map Jacobians need shape "
                "(n_segment,n_parameter)")
        n_parameter = raw.shape[1]
        dh = (np.zeros((h.size, n_parameter))
              if curvature_jacobian is None
              else np.asarray(curvature_jacobian, dtype=float))
        dk = (np.zeros((h.size, n_parameter))
              if gradient_jacobian is None
              else np.asarray(gradient_jacobian, dtype=float))
        if (dh.shape != (h.size, n_parameter)
                or dk.shape != (h.size, n_parameter)
                or not np.all(np.isfinite(np.r_[dh.ravel(), dk.ravel()]))):
            raise ValueError(
                "curvature and gradient transfer-map Jacobians must have "
                "matching finite shapes")
    else:
        n_parameter = 0
        dh = np.zeros((h.size, 0))
        dk = np.zeros((h.size, 0))

    entries = (_DEFAULT_TRANSFER_RESPONSE_ENTRIES
               if response_entries is None else
               tuple(tuple(int(value) for value in pair)
                     for pair in response_entries))
    if (not entries or any(len(pair) != 2 for pair in entries)
            or any(row < 0 or row >= 6 or column < 0 or column >= 6
                   for row, column in entries)
            or len(set(entries)) != len(entries)):
        raise ValueError(
            "response_entries must contain unique 6-by-6 matrix indices")

    # Local order: (x, x', ell, delta).  This augments the usual dispersion
    # propagator by the geometric path-length equation ell' = h*x.
    horizontal = np.eye(4)
    dhorizontal = np.zeros((n_parameter, 4, 4))
    for segment, (hi, ki, length) in enumerate(zip(h, k, ds)):
        generator = np.array([
            [0.0, 1.0, 0.0, 0.0],
            [-(hi * hi + ki), 0.0, 0.0, hi],
            [hi, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ])
        scaled = generator * length
        propagator = expm(scaled)
        old_horizontal = horizontal
        old_dhorizontal = dhorizontal.copy()
        horizontal = propagator @ old_horizontal
        for parameter in range(n_parameter):
            dhi = dh[segment, parameter]
            dki = dk[segment, parameter]
            if dhi == 0.0 and dki == 0.0:
                dhorizontal[parameter] = (
                    propagator @ old_dhorizontal[parameter])
                continue
            derivative = np.zeros((4, 4))
            derivative[1, 0] = -(2.0 * hi * dhi + dki)
            derivative[1, 3] = dhi
            derivative[2, 0] = dhi
            dpropagator = expm_frechet(
                scaled, derivative * length, compute_expm=False)
            dhorizontal[parameter] = (
                propagator @ old_dhorizontal[parameter]
                + dpropagator @ old_horizontal)

    optics = combined_function_linear_optics(
        h, k, ds,
        curvature_jacobian=(dh if supplied else None),
        gradient_jacobian=(dk if supplied else None),
        stability_tolerance=stability_tolerance)
    matrix = np.eye(6)
    matrix_jacobian = np.zeros((n_parameter, 6, 6))
    horizontal_indices = np.array([0, 1, 4, 5], dtype=np.int64)
    matrix[np.ix_(horizontal_indices, horizontal_indices)] = horizontal
    matrix[2:4, 2:4] = optics.vertical_matrix
    for parameter in range(n_parameter):
        matrix_jacobian[parameter][
            np.ix_(horizontal_indices, horizontal_indices)
        ] = dhorizontal[parameter]
        matrix_jacobian[parameter, 2:4, 2:4] = (
            optics.vertical_matrix_jacobian[parameter])
    response = np.asarray(
        [matrix[row, column] for row, column in entries], dtype=float)
    response_jacobian = np.asarray(
        [matrix_jacobian[:, row, column] for row, column in entries],
        dtype=float).reshape(len(entries), n_parameter)
    return CombinedFunctionTransferMap(
        matrix=matrix,
        matrix_jacobian=matrix_jacobian,
        response=response,
        response_jacobian=response_jacobian,
        response_entries=entries,
        optics=optics)


def combined_function_transfer_map_from_field_response(
        field_response, segment_lengths, magnetic_rigidity, *,
        field_response_jacobian=None, curvature_sign=1.0,
        gradient_sign=1.0, response_entries=None,
        stability_tolerance=1e-10):
    """Fuse row-major HDiv-MMM ``[B..., dB/dx...]`` rows into a bend map."""
    values = np.asarray(field_response, dtype=float).reshape(-1)
    lengths = np.asarray(segment_lengths, dtype=float).reshape(-1)
    rigidity = float(magnetic_rigidity)
    curvature_sign = float(curvature_sign)
    gradient_sign = float(gradient_sign)
    if (lengths.size == 0 or values.shape != (2 * lengths.size,)
            or not np.isfinite(rigidity) or rigidity == 0.0
            or not np.all(np.isfinite(
                np.r_[values, lengths, curvature_sign, gradient_sign]))):
        raise ValueError(
            "field response must contain finite B/G rows for every segment "
            "and magnetic_rigidity must be finite and nonzero")
    curvature_jacobian = None
    gradient_jacobian = None
    if field_response_jacobian is not None:
        jacobian = np.asarray(field_response_jacobian, dtype=float)
        if (jacobian.ndim != 2 or jacobian.shape[0] != values.size
                or not np.all(np.isfinite(jacobian))):
            raise ValueError(
                "field_response_jacobian needs finite shape "
                "(2*n_segment,n_parameter)")
        curvature_jacobian = (
            curvature_sign * jacobian[:lengths.size] / rigidity)
        gradient_jacobian = (
            gradient_sign * jacobian[lengths.size:] / rigidity)
    return combined_function_transfer_map(
        curvature_sign * values[:lengths.size] / rigidity,
        gradient_sign * values[lengths.size:] / rigidity,
        lengths,
        curvature_jacobian=curvature_jacobian,
        gradient_jacobian=gradient_jacobian,
        response_entries=response_entries,
        stability_tolerance=stability_tolerance)


def design_combined_function_transfer_map_target(*, length, bend_angle,
        radial_matrix, vertical_matrix, n_segments=16,
        normalized_gradient_limit=None, initial_normalized_gradient=None,
        endpoint_tolerance=1e-9, symplectic_tolerance=1e-8,
        max_iterations=1000):
    """Calculate a realizable achromatic map before electromagnetic design.

    The caller supplies stable, symplectic horizontal and vertical 2-by-2
    betatron blocks.  A piecewise-constant combined-function profile is then
    found whose full map has those blocks and ``R16=R26=0`` while preserving
    the specified total bend angle.  Only this ideal one-dimensional optics
    problem is solved here; the returned gradient profile is *not* prescribed
    to the magnet.  HDiv-MMM topology optimization instead matches the
    returned 6-by-6 matrix through
    :func:`combined_function_transfer_map_from_field_response`.

    The nonlinear least-squares Jacobian is assembled from exact matrix-
    exponential Frechet derivatives.  Finite differences are not used.
    """
    from scipy.optimize import least_squares

    length = float(length)
    angle = float(bend_angle)
    count = int(n_segments)
    tolerance = float(endpoint_tolerance)
    symplectic = float(symplectic_tolerance)
    radial = np.asarray(radial_matrix, dtype=float)
    vertical = np.asarray(vertical_matrix, dtype=float)
    if (not np.isfinite(length) or length <= 0.0
            or not np.isfinite(angle) or angle == 0.0 or count < 4
            or radial.shape != (2, 2) or vertical.shape != (2, 2)
            or not np.all(np.isfinite(np.r_[radial.ravel(), vertical.ravel()]))
            or not np.isfinite(tolerance) or tolerance <= 0.0
            or not np.isfinite(symplectic) or symplectic <= 0.0
            or int(max_iterations) <= 0):
        raise ValueError("invalid transfer-map target specification")
    determinants = np.array([np.linalg.det(radial), np.linalg.det(vertical)])
    if np.any(np.abs(determinants - 1.0) > symplectic):
        raise ValueError(
            "target radial and vertical blocks must be symplectic "
            "(determinant one)")

    lengths = np.full(count, length / count)
    curvature = np.full(count, angle / length)
    identity = np.eye(count)
    entries = (
        (0, 0), (0, 1), (1, 0), (1, 1), (0, 5), (1, 5),
        (2, 2), (2, 3), (3, 2), (3, 3),
    )
    target = np.r_[radial.ravel(), 0.0, 0.0, vertical.ravel()]
    # Balance dimensionally different matrix entries at one-percent relative
    # accuracy while keeping zero diagonal and dispersion targets observable.
    # A pure length nondimensionalization makes R16/R26 too cheap and can trap
    # the ideal-optics solve on a non-achromatic local branch.
    def block_scale(block):
        return np.array([
            max(1e-2, 1e-2 * abs(block[0, 0])),
            max(1e-3 * length, 1e-2 * abs(block[0, 1])),
            max(1e-3 / length, 1e-2 * abs(block[1, 0])),
            max(1e-2, 1e-2 * abs(block[1, 1])),
        ])
    scale = np.r_[block_scale(radial),
                  1e-4 * length, 1e-4,
                  block_scale(vertical)]

    def evaluate(values):
        return combined_function_transfer_map(
            curvature, values, lengths,
            gradient_jacobian=identity, response_entries=entries)

    def residual(values):
        return (evaluate(values).response - target) / scale

    def residual_jacobian(values):
        return evaluate(values).response_jacobian / scale[:, None]

    if initial_normalized_gradient is None:
        initial = np.zeros(count)
    else:
        initial = np.asarray(
            initial_normalized_gradient, dtype=float).reshape(-1)
        if initial.shape != (count,) or not np.all(np.isfinite(initial)):
            raise ValueError(
                "initial_normalized_gradient must match n_segments")
    if normalized_gradient_limit is None:
        lower = np.full(count, -np.inf)
        upper = np.full(count, np.inf)
    else:
        limit = float(normalized_gradient_limit)
        if not np.isfinite(limit) or limit <= 0.0:
            raise ValueError(
                "normalized_gradient_limit must be positive and finite")
        lower = np.full(count, -limit)
        upper = np.full(count, limit)
        initial = np.clip(initial, lower, upper)
    result = least_squares(
        residual, initial, jac=residual_jacobian,
        bounds=(lower, upper), max_nfev=int(max_iterations),
        xtol=1e-12, ftol=1e-12, gtol=1e-12)
    ideal = combined_function_transfer_map(
        curvature, result.x, lengths)
    requested = np.eye(6)
    requested[:2, :2] = radial
    requested[:2, 5] = 0.0
    requested[2:4, 2:4] = vertical
    # A closed-dispersion symplectic map fixes R51 and R52 to zero.  R56 is
    # not specified by the endpoint betatron blocks; retain the geometric
    # value calculated by the realizable ideal profile.
    requested[4, 0] = 0.0
    requested[4, 1] = 0.0
    requested[4, 5] = ideal.matrix[4, 5]
    maximum = float(np.max(np.abs(residual(result.x))))
    if (not result.success or maximum > tolerance
            or not ideal.optics.radial_stable
            or not ideal.optics.vertical_stable):
        raise RuntimeError(
            "transfer-map target design failed: %s; scaled residual=%.3e; "
            "traces=(%.6g,%.6g)" % (
                result.message, maximum, ideal.optics.radial_trace,
                ideal.optics.vertical_trace))
    return TransferMapTargetDesign(
        matrix=requested,
        segment_lengths=lengths,
        curvature=curvature,
        normalized_gradient=np.asarray(result.x, dtype=float),
        transfer_map=ideal,
        maximum_scaled_residual=maximum,
        iterations=int(result.nfev),
        status=str(result.message))


def transfer_map_reachability(current_response, response_jacobian,
                              target_response, response_band, *,
                              relative_tolerance=1e-10,
                              acceptance_ratio=1.0):
    """Project a target transfer map onto one linearized design space.

    This is a necessary reachability gate for topology optimization.  It solves
    the band-scaled, unconstrained least-squares problem by TSVD and reports the
    residual orthogonal to the candidate-response column space.  A failed gate
    means the current design variables cannot reproduce the requested map even
    before 0/1, volume, predecessor, and connectivity constraints are imposed.
    """
    response = np.asarray(current_response, dtype=float).reshape(-1)
    target = np.asarray(target_response, dtype=float).reshape(-1)
    band = np.asarray(response_band, dtype=float).reshape(-1)
    jacobian = np.asarray(response_jacobian, dtype=float)
    tolerance = float(relative_tolerance)
    acceptance = float(acceptance_ratio)
    if (response.size == 0 or target.shape != response.shape
            or band.shape != response.shape or np.any(band <= 0.0)
            or jacobian.ndim != 2 or jacobian.shape[0] != response.size
            or not np.all(np.isfinite(
                np.r_[response, target, band, jacobian.ravel()]))
            or not np.isfinite(tolerance) or tolerance < 0.0
            or not np.isfinite(acceptance) or acceptance < 0.0):
        raise ValueError("invalid transfer-map reachability inputs")
    scaled = jacobian / band[:, None]
    rhs = (target - response) / band
    u, singular_values, vh = np.linalg.svd(scaled, full_matrices=False)
    threshold = (tolerance * singular_values[0]
                 if singular_values.size else 0.0)
    rank = int(np.count_nonzero(singular_values > threshold))
    if rank:
        parameter_step = (vh[:rank].T
                          @ ((u[:, :rank].T @ rhs)
                             / singular_values[:rank]))
    else:
        parameter_step = np.zeros(jacobian.shape[1], dtype=float)
    predicted = response + jacobian @ parameter_step
    residual = predicted - target
    ratio = float(np.max(np.abs(residual / band)))
    return TransferMapReachability(
        numerical_rank=rank,
        singular_values=singular_values,
        parameter_step=np.asarray(parameter_step, dtype=float),
        predicted_response=np.asarray(predicted, dtype=float),
        residual=np.asarray(residual, dtype=float),
        max_normalized_residual=ratio,
        reachable=bool(ratio <= acceptance))


def combined_function_exit_metrics(
        curvature, normalized_gradient, segment_lengths, *,
        reference_curvature, downstream_drift=0.0,
        curvature_jacobian=None, gradient_jacobian=None):
    """Return analytic ``x0, psi0, eta, eta_prime`` exit responses.

    The prescribed design orbit has curvature ``reference_curvature``.  The
    central-orbit error ``q`` produced by the realized field obeys

    ``q'' + (h**2 + k1) q = h_reference - h``.

    Positive momentum error drives dispersion outward through ``+h``.  A
    positive field-curvature error bends the central particle inward, so its
    forcing has the opposite sign.

    Dispersion obeys the usual
    ``eta'' + (h**2 + k1) eta = h`` equation and is delegated to
    :func:`combined_function_linear_optics`.  Both states are transported
    through an optional straight downstream drift to the common score plane.
    Optional curvature/gradient Jacobians have shape
    ``(n_segment, n_parameter)`` and are propagated analytically.

    For a field convention ``Bz = -B0 + x*G`` with positive design curvature,
    callers normally use ``h = -Bz/B_rho`` and ``k1 = -G/B_rho``.
    """
    from scipy.linalg import expm, expm_frechet

    h = np.asarray(curvature, dtype=float).reshape(-1)
    k = np.asarray(normalized_gradient, dtype=float).reshape(-1)
    ds = np.asarray(segment_lengths, dtype=float).reshape(-1)
    href = np.asarray(reference_curvature, dtype=float)
    if href.ndim == 0:
        href = np.full(h.shape, float(href))
    else:
        href = href.reshape(-1)
    drift = float(downstream_drift)
    if h.size == 0 or k.shape != h.shape or ds.shape != h.shape:
        raise ValueError(
            "combined-function exit metrics require matching non-empty "
            "segment arrays")
    if href.shape != h.shape:
        raise ValueError(
            "reference_curvature must be scalar or match the segment arrays")
    if (np.any(ds <= 0.0) or drift < 0.0
            or not np.all(np.isfinite(np.r_[h, k, ds, href, drift]))):
        raise ValueError(
            "combined-function exit metrics inputs must be finite, with "
            "positive lengths and nonnegative downstream drift")

    optics = combined_function_linear_optics(
        h, k, ds,
        curvature_jacobian=curvature_jacobian,
        gradient_jacobian=gradient_jacobian)

    supplied = (curvature_jacobian is not None
                or gradient_jacobian is not None)
    if supplied:
        raw = (curvature_jacobian if curvature_jacobian is not None
               else gradient_jacobian)
        raw = np.asarray(raw, dtype=float)
        n_parameter = raw.shape[1]
        dh = (np.zeros((h.size, n_parameter))
              if curvature_jacobian is None
              else np.asarray(curvature_jacobian, dtype=float))
        dk = (np.zeros((h.size, n_parameter))
              if gradient_jacobian is None
              else np.asarray(gradient_jacobian, dtype=float))
    else:
        n_parameter = 0
        dh = np.zeros((h.size, 0))
        dk = np.zeros((h.size, 0))

    state = np.array([0.0, 0.0, 1.0])
    derivative_state = np.zeros((3, n_parameter))
    history = [state[:2].copy()]
    stations = [0.0]
    for segment, (hi, ki, hrefi, length) in enumerate(zip(h, k, href, ds)):
        generator = np.array([
            [0.0, 1.0, 0.0],
            [-(hi * hi + ki), 0.0, hrefi - hi],
            [0.0, 0.0, 0.0],
        ])
        scaled = generator * length
        propagator = expm(scaled)
        old_state = state
        old_derivative_state = derivative_state
        state = propagator @ old_state
        for parameter in range(n_parameter):
            dhi = dh[segment, parameter]
            dki = dk[segment, parameter]
            derivative_generator = np.zeros((3, 3))
            derivative_generator[1, 0] = -(2.0 * hi * dhi + dki)
            derivative_generator[1, 2] = -dhi
            derivative_propagator = expm_frechet(
                scaled, derivative_generator * length, compute_expm=False)
            derivative_state[:, parameter] = (
                propagator @ old_derivative_state[:, parameter]
                + derivative_propagator @ old_state)
        history.append(state[:2].copy())
        stations.append(stations[-1] + length)

    orbit_score = np.array([state[0] + drift * state[1], state[1]])
    orbit_jacobian = derivative_state[:2].copy()
    orbit_jacobian[0] += drift * orbit_jacobian[1]
    dispersion_score = np.array([
        optics.dispersion[-1, 0]
        + drift * optics.dispersion[-1, 1],
        optics.dispersion[-1, 1],
    ])
    dispersion_jacobian = optics.endpoint_jacobian.copy()
    dispersion_jacobian[0] += drift * dispersion_jacobian[1]
    response = np.r_[orbit_score, dispersion_score]
    response_jacobian = np.vstack((orbit_jacobian, dispersion_jacobian))
    return CombinedFunctionExitMetrics(
        s=np.asarray(stations),
        reference_orbit_error=np.asarray(history),
        x0_m=float(response[0]),
        psi0_rad=float(response[1]),
        eta_m=float(response[2]),
        eta_prime_rad=float(response[3]),
        response=response,
        response_jacobian=response_jacobian,
        downstream_drift_m=drift,
        optics=optics,
    )


def combined_function_exit_metrics_from_field_response(
        field_response, segment_lengths, magnetic_rigidity, *,
        reference_curvature, downstream_drift=0.0,
        field_response_jacobian=None, curvature_sign=1.0,
        gradient_sign=1.0):
    """Fuse HDiv-MMM field rows directly into the four optics responses.

    ``field_response`` is the row-major vector
    ``[B_0 ... B_(n-1), G_0 ... G_(n-1)]`` produced by the HDiv-MMM response
    matrix, where ``G=dB/dx``.  Its optional Jacobian therefore has shape
    ``(2*n_segment, n_parameter)``.  Sign arguments make the electromagnetic
    coordinate convention explicit; no silent field-axis assumption is made.
    """
    values = np.asarray(field_response, dtype=float).reshape(-1)
    lengths = np.asarray(segment_lengths, dtype=float).reshape(-1)
    rigidity = float(magnetic_rigidity)
    curvature_sign = float(curvature_sign)
    gradient_sign = float(gradient_sign)
    if (lengths.size == 0 or values.shape != (2 * lengths.size,)
            or not np.isfinite(rigidity) or rigidity == 0.0
            or not np.all(np.isfinite(
                np.r_[values, lengths, curvature_sign, gradient_sign]))):
        raise ValueError(
            "field response must contain finite B/G rows for every segment "
            "and magnetic_rigidity must be finite and nonzero")
    jacobian = None
    curvature_jacobian = None
    gradient_jacobian = None
    if field_response_jacobian is not None:
        jacobian = np.asarray(field_response_jacobian, dtype=float)
        if (jacobian.ndim != 2
                or jacobian.shape[0] != values.size
                or not np.all(np.isfinite(jacobian))):
            raise ValueError(
                "field_response_jacobian needs finite shape "
                "(2*n_segment, n_parameter)")
        curvature_jacobian = (
            curvature_sign * jacobian[:lengths.size] / rigidity)
        gradient_jacobian = (
            gradient_sign * jacobian[lengths.size:] / rigidity)
    return combined_function_exit_metrics(
        curvature_sign * values[:lengths.size] / rigidity,
        gradient_sign * values[lengths.size:] / rigidity,
        lengths,
        reference_curvature=reference_curvature,
        downstream_drift=downstream_drift,
        curvature_jacobian=curvature_jacobian,
        gradient_jacobian=gradient_jacobian,
    )


def straightened_bend_validation(s, bz, gradient, magnetic_rigidity, *,
                                 eta0=0.0, eta_prime0=0.0,
                                 initial_angle=None):
    """Reconstruct a physical bend orbit and linear optics from field samples.

    ``s`` is the straightened longitudinal coordinate of the electromagnetic
    model.  Trapezoidal cell averages of ``bz`` and ``gradient=dBz/dx`` drive
    the same combined-function equations as :func:`combined_function_linear_optics`.
    The reference orbit is integrated exactly within every constant-curvature
    cell.  When ``initial_angle`` is omitted, the orbit is centred about the
    longitudinal chord by choosing half the realized bend angle on each side.

    This is a post-solve field/optics check.  It does not calculate or
    approximate a topology or shape derivative.
    """
    s=np.asarray(s,dtype=float).reshape(-1)
    bz=np.asarray(bz,dtype=float).reshape(-1)
    gradient=np.asarray(gradient,dtype=float).reshape(-1)
    magnetic_rigidity=float(magnetic_rigidity)
    if (s.size<2 or bz.shape!=s.shape or gradient.shape!=s.shape
            or not np.all(np.isfinite(np.r_[s,bz,gradient,magnetic_rigidity]))
            or magnetic_rigidity==0.0 or np.any(np.diff(s)<=0.0)):
        raise ValueError("straightened bend inputs must be finite matching "
                         "arrays, with increasing s and nonzero rigidity")
    ds=np.diff(s)
    curvature=0.5*(bz[:-1]+bz[1:])/magnetic_rigidity
    normalized_gradient=0.5*(gradient[:-1]+gradient[1:])/magnetic_rigidity
    optics=combined_function_linear_optics(
        curvature,normalized_gradient,ds,eta0=float(eta0),
        eta_prime0=float(eta_prime0))
    increments=curvature*ds
    bend_angle=float(np.sum(increments))
    angle=np.empty(s.size,dtype=float)
    angle[0]=(-0.5*bend_angle if initial_angle is None
              else float(initial_angle))
    angle[1:]=angle[0]+np.cumsum(increments)
    position=np.zeros((s.size,3),dtype=float)
    for index,(length,h) in enumerate(zip(ds,curvature)):
        before=angle[index];after=angle[index+1]
        if abs(h)>1.0e-14:
            position[index+1,0]=(position[index,0]
                +(np.cos(before)-np.cos(after))/h)
            position[index+1,1]=(position[index,1]
                +(np.sin(after)-np.sin(before))/h)
        else:
            position[index+1,0]=position[index,0]+length*np.sin(before)
            position[index+1,1]=position[index,1]+length*np.cos(before)
    tangent=np.column_stack((np.sin(angle),np.cos(angle),
                             np.zeros_like(angle)))
    return StraightenedBendValidation(
        s=s,bz=bz,gradient=gradient,bend_angle=bend_angle,
        position=position,tangent=tangent,optics=optics)


def design_achromatic_gradient_profile(*, length, bend_angle,
        n_segments=4, normalized_gradient_limit=None,
        stability_margin=0.2, smoothness=2e-3,
        endpoint_tolerance=1e-10, max_iterations=1000,
        initial_normalized_gradient=None, gradient_sign_pattern=None):
    """Design a stable constant-curvature bend with ``eta_in=eta_out=0``.

    The incident dispersion and slope are both zero.  Only ``eta_out`` is
    constrained, matching the achromatic boundary requested by the design;
    ``eta'_out`` is reported but left free.  The longitudinal normalized
    gradient is optimized with analytic endpoint and transfer-map derivatives.
    The bend-angle integral is exact because the curvature is fixed to
    ``bend_angle/length`` in every segment.
    """
    from scipy.optimize import minimize

    length=float(length);angle=float(bend_angle);n=int(n_segments)
    margin=float(stability_margin);smooth=float(smoothness)
    if (not np.isfinite(length) or length<=0.0 or not np.isfinite(angle)
            or angle==0.0 or n<4 or not 0.0<margin<2.0 or smooth<0.0):
        raise ValueError("invalid achromatic-gradient design specification")
    ds=np.full(n,length/n);h=np.full(n,angle/length)
    # Dimensionless deterministic seed found from the four-cell alternating-
    # gradient family.  Interpolation keeps the same physical profile for a
    # finer segmentation; the 1/L^2 scaling follows Hill's equation.
    seed4=np.array([29.55146698,7.53520068,45.44214309,-36.32797540])
    centers=(np.arange(n)+0.5)/n
    if initial_normalized_gradient is None:
        seed=np.interp(centers,(np.arange(4)+0.5)/4,seed4)/length**2
    else:
        seed=np.asarray(initial_normalized_gradient,dtype=float).reshape(-1)
        if seed.shape!=(n,) or not np.all(np.isfinite(seed)):
            raise ValueError("initial_normalized_gradient must match n_segments")
    scale=max(1.0,1.0/length**2)
    identity=np.eye(n)

    def evaluate(values):
        return combined_function_linear_optics(
            h,values,ds,gradient_jacobian=identity)

    difference=np.zeros((max(0,n-1),n))
    for row in range(n-1):
        difference[row,row]=-1.0;difference[row,row+1]=1.0
    def objective(values):
        return float(0.5*np.dot(values,values)/scale**2
                     +0.5*smooth*np.dot(difference@values,difference@values)/scale**2)
    def objective_jacobian(values):
        return ((values+smooth*difference.T@(difference@values))/scale**2)
    limit=2.0-margin
    def equality(values):
        return np.array([evaluate(values).dispersion[-1,0]])
    def equality_jacobian(values):
        return evaluate(values).endpoint_jacobian[0:1]
    def stability(values):
        optics=evaluate(values)
        return np.array([limit-optics.radial_trace,limit+optics.radial_trace,
                         limit-optics.vertical_trace,limit+optics.vertical_trace])
    def stability_jacobian(values):
        optics=evaluate(values)
        tx=np.trace(optics.radial_matrix_jacobian,axis1=1,axis2=2)
        ty=np.trace(optics.vertical_matrix_jacobian,axis1=1,axis2=2)
        return np.vstack((-tx,tx,-ty,ty))
    bound=(None if normalized_gradient_limit is None else
           float(normalized_gradient_limit))
    if bound is not None and (not np.isfinite(bound) or bound<=0.0):
        raise ValueError("normalized_gradient_limit must be positive and finite")
    lower=np.full(n,-np.inf if bound is None else -bound)
    upper=np.full(n,np.inf if bound is None else bound)
    if gradient_sign_pattern is not None:
        signs=np.asarray(gradient_sign_pattern,dtype=int).reshape(-1)
        if signs.shape!=(n,) or np.any(~np.isin(signs,[-1,0,1])):
            raise ValueError("gradient_sign_pattern entries must be -1, 0, or +1")
        lower[signs>0]=0.0;upper[signs<0]=0.0
    seed=np.clip(seed,lower,upper)
    bounds=[((None if not np.isfinite(lo) else float(lo)),
             (None if not np.isfinite(hi) else float(hi)))
            for lo,hi in zip(lower,upper)]
    result=minimize(objective,seed,jac=objective_jacobian,method="SLSQP",
        bounds=bounds,constraints=(
            {"type":"eq","fun":equality,"jac":equality_jacobian},
            {"type":"ineq","fun":stability,"jac":stability_jacobian}),
        options={"ftol":1e-12,"maxiter":int(max_iterations),"disp":False})
    optics=evaluate(result.x)
    residual=abs(float(optics.dispersion[-1,0]))
    if (not result.success or residual>float(endpoint_tolerance)
            or not optics.radial_stable or not optics.vertical_stable):
        raise RuntimeError(
            "achromatic gradient design failed: %s; |eta_out|=%.3e, traces=(%.6g,%.6g)"
            %(result.message,residual,optics.radial_trace,optics.vertical_trace))
    return AchromaticGradientDesign(ds,h,np.asarray(result.x),optics,
        float(result.fun),int(result.nit),str(result.message))


def isochronous_profile_metrics(radii, bz, gamma):
    """Relative field-shape and revolution-period errors for ``B~gamma``."""
    radii = np.asarray(radii, dtype=float).ravel()
    bz = np.asarray(bz, dtype=float).ravel()
    gamma = np.asarray(gamma, dtype=float).ravel()
    if radii.size < 2 or bz.shape != radii.shape or gamma.shape != radii.shape:
        raise ValueError("isochronous_profile_metrics: matching vectors of length >=2 required")
    if (not np.all(np.isfinite(radii)) or not np.all(np.isfinite(bz))
            or not np.all(np.isfinite(gamma)) or np.any(bz == 0.0)
            or np.any(gamma <= 0.0)):
        raise ValueError("isochronous_profile_metrics: invalid profile")
    reference = radii.size//2
    normalized = (bz/bz[reference])/(gamma/gamma[reference])
    field_error = normalized-1.0
    period_error = 1.0/normalized-1.0
    return dict(
        normalized_field=normalized, field_error=field_error,
        period_error=period_error,
        max_abs_field_error=float(np.max(np.abs(field_error))),
        rms_field_error=float(np.sqrt(np.mean(field_error**2))),
        max_abs_period_error=float(np.max(np.abs(period_error))),
        rms_period_error=float(np.sqrt(np.mean(period_error**2))))


def isochronous_increment_targets(gamma, reference_increment, external_bz,
                                  reference_index=None):
    """Targets for a design-dependent field increment under a fixed field.

    The isochronous law applies to total ``Bz``.  If the optimizer controls
    only ``delta_Bz`` while a fixed coil contributes ``external_bz``, the
    correct target is

    ``(external_bz + reference_increment) * gamma/gamma_ref - external_bz``.

    Applying ``gamma/gamma_ref`` to the increment alone is incorrect because
    the additive fixed field does not cancel in a ratio.
    """
    gamma = np.asarray(gamma, dtype=float).ravel()
    if gamma.size < 2 or not np.all(np.isfinite(gamma)) or np.any(gamma <= 0.0):
        raise ValueError("isochronous_increment_targets: gamma must be a "
                         "positive finite vector of length >=2")
    reference = gamma.size//2 if reference_index is None else int(reference_index)
    if reference < 0 or reference >= gamma.size:
        raise ValueError("isochronous_increment_targets: reference_index out of range")
    reference_increment = float(reference_increment)
    external_bz = float(external_bz)
    if not np.isfinite(reference_increment) or not np.isfinite(external_bz):
        raise ValueError("isochronous_increment_targets: field values must be finite")
    return ((external_bz + reference_increment) * gamma / gamma[reference]
            - external_bz)


def isochronous_total_field_bands(increment_targets, external_bz,
                                  relative_band):
    """Absolute increment-functional bands defined relative to total ``Bz``."""
    increments = np.asarray(increment_targets, dtype=float).ravel()
    external_bz = float(external_bz)
    relative_band = float(relative_band)
    if (increments.size < 1 or not np.all(np.isfinite(increments))
            or not np.isfinite(external_bz)
            or not np.isfinite(relative_band) or relative_band <= 0.0):
        raise ValueError("isochronous_total_field_bands: need finite targets, "
                         "field, and positive relative_band")
    bands = relative_band * np.abs(external_bz + increments)
    if np.any(bands == 0.0):
        raise ValueError("isochronous_total_field_bands: total target is zero")
    return bands


def dipole_array_field_cf(points, moments):
    """Analytic H field CoefficientFunction of a point-dipole array.

    ``H(y) = (1/4pi) sum_i [3 (m_i . r) r/|r|^5 - m_i/|r|^3]``, ``r = y - x_i``.
    Singular at the dipole points -- integrate it only over regions that
    exclude them (:func:`field_functional_load` enforces this).
    """
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    mom = np.asarray(moments, dtype=float).reshape(-1, 3)
    if len(pts) != len(mom):
        raise ValueError("dipole_array_field_cf: %d points but %d moments"
                         % (len(pts), len(mom)))
    if not len(pts):
        raise ValueError("dipole_array_field_cf: empty dipole array")
    total = None
    for (px, py, pz), (mx, my, mz) in zip(pts, mom):
        rx, ry, rz = ng.x - px, ng.y - py, ng.z - pz
        r2 = rx * rx + ry * ry + rz * rz
        r1 = ng.sqrt(r2)
        r3 = r2 * r1
        r5 = r2 * r3
        mdotr = mx * rx + my * ry + mz * rz
        term = ng.CoefficientFunction((3.0 * mdotr * rx / r5 - mx / r3,
                                       3.0 * mdotr * ry / r5 - my / r3,
                                       3.0 * mdotr * rz / r5 - mz / r3))
        total = term if total is None else total + term
    return total * (1.0 / (4.0 * pi))


def _require_points_outside(mesh, pts, who):
    """Fail loud when an evaluation point lies inside the meshed iron."""
    mapped = mesh(pts[:, 0], pts[:, 1], pts[:, 2])
    if not (isinstance(mapped, np.ndarray) and mapped.dtype.names
            and "nr" in mapped.dtype.names):
        raise RuntimeError(
            "%s: NGSolve's vectorized MeshPoint array API is required for the "
            "inside-mesh check" % who)
    inside = np.flatnonzero(np.asarray(mapped["nr"]) >= 0)
    if len(inside):
        raise ValueError(
            "%s: %d evaluation/dipole point(s) lie INSIDE the meshed body "
            "(first: %s) -- the dipole reciprocity load is singular there; "
            "orbit points must be outside the iron mesh"
            % (who, len(inside), pts[inside[0]].tolist()))


def field_functional_load(fes, points, weights, axis=2, scale=MU0,
                          bonus_intorder=10):
    """Assembled adjoint load ``f`` with ``f^T m = scale * sum_i w_i H_d[m]_axis(x_i)``.

    Kernel reciprocity: the load is the mass projection of the analytic field
    of dipoles ``w_i e_axis`` at ``x_i``.  With ``scale=MU0`` (default),
    ``f^T m`` is the magnetization-dependent part of
    ``sum_i w_i B_axis(x_i)`` in Tesla; the design-independent coil part
    ``scale * sum_i w_i H_ext_axis(x_i)`` is the caller's constant.
    ``bonus_intorder=10`` reproduces the independent C++ charge evaluator to
    ~1e-10 relative on the verification geometry (points >= half an element
    size outside the body); the returned value is exact FD-consistent for the
    adjoint regardless of quadrature, since J is DEFINED as ``f^T m``.
    """
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    wts = np.asarray(weights, dtype=float).ravel()
    if len(pts) != len(wts):
        raise ValueError("field_functional_load: %d points but %d weights"
                         % (len(pts), len(wts)))
    if axis not in (0, 1, 2):
        raise ValueError("field_functional_load: axis must be 0, 1, or 2")
    _require_points_outside(fes.mesh, pts, "field_functional_load")
    moments = np.zeros((len(wts), 3))
    moments[:, axis] = wts
    cf = dipole_array_field_cf(pts, moments)
    v = fes.TestFunction()
    f = ng.LinearForm(fes)
    f += float(scale) * ng.InnerProduct(cf, v) * ng.dx(
        bonus_intorder=int(bonus_intorder))
    f.Assemble()
    return f


def uniform_field_load(fes, H_ext):
    """Assembled state load for a uniform external field ``H_ext`` (A/m)."""
    H = tuple(float(c) for c in H_ext)
    if len(H) != 3:
        raise ValueError("uniform_field_load: H_ext must have 3 components")
    v = fes.TestFunction()
    f = ng.LinearForm(fes)
    f += ng.CoefficientFunction(H) * v * ng.dx
    f.Assemble()
    return f


def demag_field_from_solution(demag, gfM, points, algorithm="direct"):
    """Demagnetizing field H_d (A/m) of a solved magnetization at ``points``.

    Independent verification/reporting path: evaluates the exact analytic
    field of the BDM charge representation through the C++ evaluator of
    ``vim.FieldFromSolution``.  This is the ONE sanctioned coupling point to
    the vim result-dict internals (the same keys ``vim.Solve`` stores).
    """
    if not isinstance(demag, DemagOperator):
        raise TypeError("demag_field_from_solution: demag must be a vim.DemagOperator")
    if gfM.space is not demag.space:
        raise ValueError("demag_field_from_solution: gfM does not live on the "
                         "operator's HDiv space")
    res = {"gfM": gfM, "order": int(demag.space.globalorder),
           "_charge_gram": demag._G,
           "_m_coefficients": np.ascontiguousarray(
               gfM.vec.FV().NumPy(), dtype=np.float64)}
    return FieldFromSolution(res, points, algorithm=algorithm)


def demag_field_evaluator(demag, gfM, algorithm="direct"):
    """Return a reusable callable for one solved demagnetizing field.

    The callable retains the immutable native source evaluator after its first
    query.  Use this for RK4 orbit tracking, which evaluates the same final
    magnetization thousands of times.
    """
    if not isinstance(demag, DemagOperator):
        raise TypeError("demag_field_evaluator: demag must be a vim.DemagOperator")
    if gfM.space is not demag.space:
        raise ValueError("demag_field_evaluator: gfM does not live on the "
                         "operator's HDiv space")
    result = {"gfM": gfM, "order": int(demag.space.globalorder),
              "_charge_gram": demag._G,
              "_m_coefficients": np.ascontiguousarray(
                  gfM.vec.FV().NumPy(), dtype=np.float64)}

    def evaluate(points):
        values = np.asarray(
            FieldFromSolution(result, points, algorithm=algorithm), dtype=float)
        return values[0] if np.asarray(points).shape == (3,) else values

    return evaluate


# --------------------------------------------------------------------------
# density filter
# --------------------------------------------------------------------------
class HelmholtzFilter:
    """Element-density Helmholtz filter ``rho_f = D_V^-1 B^T K^-1 B rho``.

    ``K = (u v + radius^2 grad u . grad v)`` on H1 order 1 (SPD, factored
    once); ``B_je = int_e phi_j``; ``D_V`` = element volumes.  ``radius`` is
    the minimum-feature length (a manufacturability input).  The K^-1 core is
    self-adjoint, so the gradient chain is the transpose map
    ``dJ/drho = B^T K^-1 B (D_V^-1 dJ/drho_f)`` -- :meth:`chain` -- exact to
    solver precision (locked by the filter FD test).

    The caller wraps construction and every apply in ``with TaskManager():``.
    """

    def __init__(self, mesh, radius):
        if not radius > 0.0:
            raise ValueError("HelmholtzFilter: radius must be positive")
        self.mesh = mesh
        self.radius = float(radius)
        fes = ng.H1(mesh, order=1)
        u, v = fes.TnT()
        a = ng.BilinearForm(fes, symmetric=True)
        a += (u * v + self.radius ** 2 * ng.grad(u) * ng.grad(v)) * ng.dx
        a.Assemble()
        self._fes = fes
        self._inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        self._volumes = np.asarray(
            ng.Integrate(ng.CoefficientFunction(1.0), mesh,
                         element_wise=True), dtype=float)
        self._l2 = ng.L2(mesh, order=0)
        if self._l2.ndof != int(mesh.ne):
            raise RuntimeError("HelmholtzFilter: L2(order=0) dof count %d != "
                               "element count %d" % (self._l2.ndof, mesh.ne))
        self._rho_gf = ng.GridFunction(self._l2)
        self._u_gf = ng.GridFunction(fes)

    def _bt_kinv_b(self, values):
        values = np.ascontiguousarray(values, dtype=float).ravel()
        if values.size != self._volumes.size:
            raise ValueError("HelmholtzFilter: vector has %d entries, mesh "
                             "has %d elements" % (values.size,
                                                  self._volumes.size))
        self._rho_gf.vec.FV().NumPy()[:] = values
        f = ng.LinearForm(self._fes)
        f += self._rho_gf * self._fes.TestFunction() * ng.dx
        f.Assemble()
        self._u_gf.vec.data = self._inv * f.vec
        return np.asarray(ng.Integrate(self._u_gf, self.mesh,
                                       element_wise=True), dtype=float)

    def apply(self, density):
        """Filtered element densities (same shape, smoothed over ``radius``)."""
        return self._bt_kinv_b(density) / self._volumes

    def chain(self, gradient_filtered):
        """Transpose map: gradient w.r.t. raw density from the filtered one."""
        return self._bt_kinv_b(
            np.asarray(gradient_filtered, dtype=float) / self._volumes)


class HeavisideProjection:
    """Smooth density projection with an exact analytic chain rule.

    ``beta`` controls sharpness and ``eta`` is the transition density.  Use a
    continuation sequence (for example 1, 2, 4, 8) across design phases; a
    large beta from a flat start suppresses useful topology gradients.
    """

    def __init__(self, beta, eta=0.5):
        self.beta = float(beta)
        self.eta = float(eta)
        if not np.isfinite(self.beta) or self.beta <= 0.0:
            raise ValueError("HeavisideProjection: beta must be positive")
        if not np.isfinite(self.eta) or not 0.0 < self.eta < 1.0:
            raise ValueError("HeavisideProjection: eta must lie in (0, 1)")
        self._denominator = (np.tanh(self.beta*self.eta)
                             + np.tanh(self.beta*(1.0-self.eta)))

    def apply(self, density):
        rho = np.asarray(density, dtype=float)
        if (not np.all(np.isfinite(rho)) or np.any(rho < 0.0)
                or np.any(rho > 1.0)):
            raise ValueError("HeavisideProjection: density must lie in [0, 1]")
        return ((np.tanh(self.beta*self.eta)
                 + np.tanh(self.beta*(rho-self.eta))) / self._denominator)

    def derivative(self, density):
        rho = np.asarray(density, dtype=float)
        if (not np.all(np.isfinite(rho)) or np.any(rho < 0.0)
                or np.any(rho > 1.0)):
            raise ValueError("HeavisideProjection: density must lie in [0, 1]")
        return (self.beta / self._denominator
                / np.cosh(self.beta*(rho-self.eta))**2)

    def chain(self, density, gradient_projected):
        gradient = np.asarray(gradient_projected, dtype=float)
        rho = np.asarray(density, dtype=float)
        if gradient.shape != rho.shape:
            raise ValueError("HeavisideProjection.chain: shape mismatch")
        return gradient*self.derivative(rho)


def restore_projected_volume(density, element_volumes, volume_fraction, *,
                             density_filter=None, density_projection=None,
                             volume_tolerance=1e-12,
                             density_tolerance=1e-12,
                             max_iterations=80):
    """Make a continuation start feasible after changing its projection.

    Increasing a Heaviside ``beta`` can increase the projected material
    volume even though the raw design is unchanged.  This routine finds the
    smallest uniform downward shift

    ``rho_feasible = clip(rho - shift, 0, 1)``

    whose *filtered and projected* volume satisfies ``volume_fraction``.
    The scalar correction preserves the ordering of all non-clipped design
    variables, so a continuation phase does not invent a new topology merely
    to repair feasibility.  It is intended for phase transitions and restart
    recovery, not as a replacement for the analytic volume row used by
    :func:`optimize_density` within a phase.

    Returns ``(rho_feasible, diagnostics)``.  ``density_filter`` and
    ``density_projection`` follow the same contracts as
    :func:`optimize_density`.
    """
    rho = np.asarray(density, dtype=float).copy()
    volumes = np.asarray(element_volumes, dtype=float).reshape(-1)
    if rho.ndim != 1 or rho.shape != volumes.shape or rho.size == 0:
        raise ValueError(
            "restore_projected_volume: density and element_volumes must "
            "be non-empty vectors with identical shapes")
    if (not np.all(np.isfinite(rho)) or np.any(rho < 0.0)
            or np.any(rho > 1.0)):
        raise ValueError(
            "restore_projected_volume: density must be finite in [0, 1]")
    if (not np.all(np.isfinite(volumes)) or np.any(volumes <= 0.0)):
        raise ValueError(
            "restore_projected_volume: element_volumes must be finite and "
            "positive")
    volume_fraction = float(volume_fraction)
    if not np.isfinite(volume_fraction) or not 0.0 < volume_fraction <= 1.0:
        raise ValueError(
            "restore_projected_volume: volume_fraction must be in (0, 1]")
    volume_tolerance = float(volume_tolerance)
    density_tolerance = float(density_tolerance)
    max_iterations_i = int(max_iterations)
    if (not np.isfinite(volume_tolerance) or volume_tolerance <= 0.0
            or not np.isfinite(density_tolerance)
            or density_tolerance <= 0.0):
        raise ValueError(
            "restore_projected_volume: tolerances must be positive and finite")
    if max_iterations_i < 1 or max_iterations_i != max_iterations:
        raise ValueError(
            "restore_projected_volume: max_iterations must be a positive integer")

    volume_limit = float(volume_fraction * volumes.sum())
    absolute_tolerance = volume_tolerance * max(1.0, abs(volume_limit))

    def transformed(candidate):
        if density_filter is None:
            material = candidate
        else:
            material = np.clip(density_filter.apply(candidate), 0.0, 1.0)
        if density_projection is not None:
            material = density_projection.apply(material)
        return np.asarray(material, dtype=float)

    volume_before = float(volumes @ transformed(rho))
    if volume_before <= volume_limit + absolute_tolerance:
        return rho, dict(
            changed=False, shift=0.0, iterations=0,
            volume_before=volume_before, volume_after=volume_before,
            volume_limit=volume_limit,
            relative_excess_before=(volume_before-volume_limit)
            / max(abs(volume_limit), 1e-300))

    # shift=1 maps every admissible raw density to zero and must therefore
    # bracket a feasible point for the standard Helmholtz/projection maps.
    lo, hi = 0.0, 1.0
    volume_hi = float(volumes @ transformed(np.zeros_like(rho)))
    if volume_hi > volume_limit + absolute_tolerance:
        raise RuntimeError(
            "restore_projected_volume: zero density is not volume-feasible; "
            "the supplied filter/projection does not preserve void")
    iterations = 0
    while iterations < max_iterations_i and hi-lo > density_tolerance:
        mid = 0.5*(lo+hi)
        candidate = np.clip(rho-mid, 0.0, 1.0)
        volume_mid = float(volumes @ transformed(candidate))
        iterations += 1
        if volume_mid > volume_limit:
            lo = mid
        else:
            hi, volume_hi = mid, volume_mid
        if (volume_limit-volume_hi >= 0.0
                and volume_limit-volume_hi <= absolute_tolerance):
            break
    feasible = np.clip(rho-hi, 0.0, 1.0)
    volume_after = float(volumes @ transformed(feasible))
    if volume_after > volume_limit + absolute_tolerance:
        raise RuntimeError(
            "restore_projected_volume: scalar feasibility search did not "
            "reach the requested volume budget")
    return feasible, dict(
        changed=True, shift=float(hi), iterations=iterations,
        volume_before=volume_before, volume_after=volume_after,
        volume_limit=volume_limit,
        relative_excess_before=(volume_before-volume_limit)
        / max(abs(volume_limit), 1e-300))


def density_discreteness(density, *, lower=0.1, upper=0.9):
    """Return stable grayness diagnostics for a continuous density design."""
    rho = np.asarray(density, dtype=float).reshape(-1)
    if (rho.size == 0 or not np.all(np.isfinite(rho))
            or np.any(rho < 0.0) or np.any(rho > 1.0)):
        raise ValueError("density_discreteness: density must be non-empty in [0,1]")
    lower, upper = float(lower), float(upper)
    if not 0.0 <= lower < upper <= 1.0:
        raise ValueError("density_discreteness: need 0 <= lower < upper <= 1")
    intermediate = (rho > lower) & (rho < upper)
    return {
        "intermediate_fraction": float(np.mean(intermediate)),
        "binary_distance_mean": float(np.mean(np.minimum(rho, 1.0-rho))),
        "binary_distance_max": float(np.max(np.minimum(rho, 1.0-rho))),
    }


def iron_only_verification_ready(density, *, maximum_intermediate_fraction=0.05,
                                 lower=0.1, upper=0.9):
    """Whether thresholded exact-void verification is scientifically useful."""
    limit = float(maximum_intermediate_fraction)
    if not np.isfinite(limit) or not 0.0 <= limit <= 1.0:
        raise ValueError("maximum_intermediate_fraction must lie in [0,1]")
    metrics = density_discreteness(density, lower=lower, upper=upper)
    return metrics["intermediate_fraction"] <= limit, metrics


# --------------------------------------------------------------------------
# adjoint solver
# --------------------------------------------------------------------------
@dataclass
class AdjointGradientResult:
    """State + adjoint solution pair with objective and per-element gradient."""
    objective: float          # f_adjoint^T m  (design-dependent part of J)
    gradient: np.ndarray      # dJ/ds_e, shape (n_el,)
    gfM: "ng.GridFunction"    # state magnetization m
    gfLambda: "ng.GridFunction"  # adjoint field lambda
    state_iterations: int
    adjoint_iterations: int


@dataclass
class FunctionalLinearization:
    """State + one adjoint per functional, all on one configured operator.

    ``values[k] = f_k^T m`` and ``jacobians[k, e] = d(f_k^T m)/ds_e`` for the
    functional loads passed to :meth:`DensityAdjointVIM.linearize` --
    functional 0 is conventionally the objective, the rest are constraints.
    """
    values: np.ndarray        # (n_loads,)
    jacobians: np.ndarray     # (n_loads, n_el)
    gfM: "ng.GridFunction"
    gfLambdas: tuple          # one adjoint GridFunction per load
    state_iterations: int
    adjoint_iterations: tuple


class DensityAdjointVIM:
    """Build-once HDiv-MMM operator with per-element density adjoint gradients.

    ``fes`` is the HDiv space on the design mesh (order 0, 1, or 2 where the
    selected element topology supports it).  Material
    topology with element-wise iron/void interfaces must use
    ``ng.HDiv(..., discontinuous=True)`` and
    ``internal_interfaces=True`` so the ChargeGram contains the physical
    jump of ``M.n`` on every element facet.  A conforming HDiv space is the
    body-fitted single-material path and cannot create a new internal boundary
    merely by changing an element's constitutive coefficient.  The geometry
    operator ``N`` is built once (``demag=None``) or shared from an existing
    :class:`radia.vim.DemagOperator`; every design iterate then costs one
    weighted-mass assembly and the state/adjoint CG solves.  The caller
    wraps construction and every solve in ``with TaskManager():``.
    """

    def __init__(self, fes, demag=None, **gram_kwargs):
        if demag is None:
            demag = DemagOperator(fes, **gram_kwargs)
        else:
            if gram_kwargs:
                raise ValueError(
                    "DensityAdjointVIM: gram kwargs are only accepted when the "
                    "operator is built here (demag=None)")
            if demag.space is not fes:
                raise ValueError(
                    "DensityAdjointVIM: the supplied DemagOperator lives on a "
                    "different HDiv space")
        self.fes = fes
        self.mesh = fes.mesh
        self.demag = demag
        self.n_el = int(self.mesh.ne)
        self._l2 = ng.L2(self.mesh, order=0)
        if self._l2.ndof != self.n_el:
            raise RuntimeError(
                "DensityAdjointVIM: L2(order=0) dof count %d != element count %d"
                % (self._l2.ndof, self.n_el))
        self._s_gf = ng.GridFunction(self._l2)
        self._volumes = None

    @property
    def element_volumes(self):
        """Per-element volumes (the LP volume-budget row), shape (n_el,)."""
        if self._volumes is None:
            self._volumes = np.asarray(
                ng.Integrate(ng.CoefficientFunction(1.0), self.mesh,
                             element_wise=True), dtype=float)
        return self._volumes

    # -------------------------------------------------------------- internals
    def _system(self, s, *, need_ngsolve_solver=True):
        s = np.ascontiguousarray(s, dtype=float).ravel()
        if s.size != self.n_el:
            raise ValueError("DensityAdjointVIM: s has %d entries, mesh has %d "
                             "elements" % (s.size, self.n_el))
        if not np.all(s > 0.0):
            raise ValueError(
                "DensityAdjointVIM: s = 1/chi must be strictly positive "
                "(min=%r); map densities through density_to_s" % float(s.min()))
        self._s_gf.vec.FV().NumPy()[:] = s
        u, v = self.fes.TnT()
        mass = ng.BilinearForm(self.fes, symmetric=True)
        mass += self._s_gf * u * v * ng.dx
        mass.Assemble()
        if not need_ngsolve_solver:
            return mass, None, None
        A = mass.mat + self.demag.mat
        pre = mass.mat.Inverse(self.fes.FreeDofs(), inverse="sparsecholesky")
        return mass, A, pre

    def _cg(self, A, pre, rhs_vec, gf, tol, maxiter, warm_vec=None):
        # krylovspace's relative tol is anchored to the FIRST residual of the
        # run, which for a warm start is already tiny -- that would make warm
        # restarts HARDER, not cheaper.  Anchor the stopping criterion to the
        # x0=0 first residual ||rhs||_pre instead, so warm and cold solves
        # share one absolute criterion (cold behavior is unchanged).
        work = rhs_vec.CreateVector()
        work.data = pre * rhs_vec
        reference = abs(ng.InnerProduct(work, rhs_vec)) ** 0.5
        if reference == 0.0:
            gf.vec[:] = 0.0
            return 0
        inv = CGSolver(A, pre=pre, tol=None, atol=tol * reference,
                       maxiter=int(maxiter), printrates=False)
        if warm_vec is not None:
            gf.vec.data = warm_vec
            inv.Solve(rhs=rhs_vec, sol=gf.vec, initialize=False)
        else:
            inv.Solve(rhs=rhs_vec, sol=gf.vec, initialize=True)
        iters = int(inv.iterations)
        if iters >= int(maxiter):
            raise RuntimeError(
                "DensityAdjointVIM: CG did not converge within %d iterations "
                "(tol=%g); raise maxiter or inspect the design state"
                % (int(maxiter), tol))
        return iters

    def _native_solve_many(self, mass, rhs_vecs, tol, maxiter,
                           warm_vecs=None, mass_riesz=True,
                           cluster_tree=False, cluster_coarse_size=64,
                           cluster_deflation_size=8, recycle_size=8):
        """Solve shared-matrix VIM systems in the configured C++ kernel.

        The weighted HDiv mass is registered once.  H-matrix application,
        Krylov updates, and true-residual stopping tests remain inside C++.
        ``mass_riesz=False`` uses the inexpensive exact system diagonal;
        ``True`` reuses the persistent PARDISO mass factor and applies it to
        all right-hand sides in one phase-33 call.  With
        ``cluster_tree=True``, all right-hand sides cross one
        row-major native boundary; the preserved H-matrix cluster tree supplies
        aggregate ``D^-1 B^T`` modes to a balanced two-level preconditioner,
        and converged columns form a small Ritz recycle space for later columns.
        The operator and true-residual PCG convergence contract are unchanged.
        """
        rhs_vecs = list(rhs_vecs)
        use_many = bool(cluster_tree or (mass_riesz and len(rhs_vecs) > 1))
        gram = getattr(self.demag, "_G", None)
        solve_name = (
            "solve_configured_linear_material_auto_prec_many"
            if use_many else
            ("solve_configured_linear_material_mass_riesz"
             if mass_riesz else
             "solve_configured_linear_material_auto_prec"))
        required = ("configure_mass_matrix_ngsolve", solve_name)
        if gram is None or any(not hasattr(gram, name) for name in required):
            raise RuntimeError(
                "DensityAdjointVIM: native H-matrix CG is unavailable "
                "on this DemagOperator; pass solver='ngsolve-cg' only for "
                "the explicit reference path")
        gram.configure_mass_matrix_ngsolve(mass.mat)
        if warm_vecs is None:
            warm_vecs = [None] * len(rhs_vecs)
        else:
            warm_vecs = list(warm_vecs)
            if len(warm_vecs) != len(rhs_vecs):
                raise ValueError(
                    "DensityAdjointVIM: native warm-start count mismatch")
        if use_many:
            rhs_matrix = np.ascontiguousarray(np.stack([
                np.asarray(rhs.FV().NumPy(), dtype=float)
                for rhs in rhs_vecs]), dtype=float)
            if all(value is None for value in warm_vecs):
                x0 = None
            else:
                x0 = np.ascontiguousarray(np.stack([
                    np.zeros(rhs_matrix.shape[1], dtype=float)
                    if value is None else
                    np.asarray(value.FV().NumPy(), dtype=float)
                    for value in warm_vecs]), dtype=float)
            result = gram.solve_configured_linear_material_auto_prec_many(
                1.0, rhs_matrix, float(tol), int(maxiter),
                int(cluster_coarse_size if cluster_tree else 0),
                int(cluster_deflation_size if cluster_tree else 0),
                int(recycle_size if cluster_tree else 0),
                mass_riesz=bool(mass_riesz), x0=x0)
            solution = np.asarray(result["m"], dtype=float)
            iterations = [int(value) for value in result["iters"]]
            if solution.shape != rhs_matrix.shape or len(iterations) != len(rhs_vecs):
                raise RuntimeError(
                    "DensityAdjointVIM: invalid native multi-RHS result shape")
            if any(value >= int(maxiter) for value in iterations):
                raise RuntimeError(
                    "DensityAdjointVIM: native batched %s CG did not "
                    "converge within %d iterations (tol=%g)"
                    % ("mass-Riesz" if mass_riesz else "cluster-tree",
                       int(maxiter), tol))
            fields = []
            for row in solution:
                gf = ng.GridFunction(self.fes)
                gf.vec.FV().NumPy()[:] = row
                fields.append(gf)
            return fields, iterations

        fields, iterations = [], []
        for rhs, warm_vec in zip(rhs_vecs, warm_vecs):
            rhs_array = np.ascontiguousarray(rhs.FV().NumPy(), dtype=float)
            x0 = None if warm_vec is None else np.ascontiguousarray(
                warm_vec.FV().NumPy(), dtype=float)
            if mass_riesz:
                result = gram.solve_configured_linear_material_mass_riesz(
                    1.0, rhs_array, float(tol), int(maxiter), True, x0=x0)
            else:
                result = gram.solve_configured_linear_material_auto_prec(
                    1.0, rhs_array, float(tol), int(maxiter), x0=x0)
            iters = int(result["iters"])
            if iters >= int(maxiter):
                raise RuntimeError(
                    "DensityAdjointVIM: native %s CG did not converge "
                    "within %d iterations (tol=%g)"
                    % ("mass-Riesz" if mass_riesz else "Jacobi",
                       int(maxiter), tol))
            gf = ng.GridFunction(self.fes)
            gf.vec.FV().NumPy()[:] = np.asarray(result["m"], dtype=float)
            fields.append(gf)
            iterations.append(iters)
        return fields, iterations

    # ------------------------------------------------------------ public API
    def solve(self, s, load, tol=1e-12, maxiter=20000, warm=None,
              solver="native-jacobi"):
        """One SPD solve ``(M_s + N) x = load.vec``.

        ``warm`` is an optional GridFunction whose vector seeds CG (warm start
        across design iterates).  Returns ``(GridFunction, iterations)``.
        """
        mass, A, pre = self._system(
            s, need_ngsolve_solver=(solver == "ngsolve-cg"))
        if solver in {"native", "native-jacobi", "native-batch", "native-cluster"}:
            fields, iterations = self._native_solve_many(
                mass, [load.vec], tol, maxiter,
                warm_vecs=None if warm is None else [warm.vec],
                mass_riesz=(solver == "native"),
                cluster_tree=(solver in {"native-batch", "native-cluster"}),
                cluster_coarse_size=(64 if solver == "native-cluster" else 0),
                cluster_deflation_size=(8 if solver == "native-cluster" else 0))
            return fields[0], iterations[0]
        if solver != "ngsolve-cg":
            raise ValueError(
                "DensityAdjointVIM.solve: solver must be 'native-batch', 'native-cluster', "
                "'native-jacobi', 'native', or 'ngsolve-cg'")
        gf = ng.GridFunction(self.fes)
        iters = self._cg(A, pre, load.vec, gf, tol, maxiter,
                         warm_vec=None if warm is None else warm.vec)
        return gf, iters

    def linearize(self, s, state_load, functional_loads,
                  tol=1e-12, maxiter=20000, warm=None,
                  solver="native-jacobi"):
        """State + one adjoint per functional, sharing one configured matrix.

        The self-adjoint operator serves every solve; ``warm`` is an optional
        :class:`FunctionalLinearization` from the previous design iterate
        whose fields seed all CG solves.  Sensitivities come from one
        element-wise ``Integrate`` per functional:
        ``jacobians[k] = -Integrate(lambda_k . m, element_wise=True)``.
        """
        loads = list(functional_loads)
        if not loads:
            raise ValueError("DensityAdjointVIM.linearize: no functional loads")
        if warm is not None and len(warm.gfLambdas) != len(loads):
            raise ValueError(
                "DensityAdjointVIM.linearize: warm start carries %d adjoints "
                "for %d loads" % (len(warm.gfLambdas), len(loads)))
        mass, A, pre = self._system(
            s, need_ngsolve_solver=(solver == "ngsolve-cg"))
        if solver in {"native", "native-jacobi", "native-batch", "native-cluster"}:
            rhs = [state_load.vec] + [load.vec for load in loads]
            warm_vecs = None if warm is None else [warm.gfM.vec] + [
                value.vec for value in warm.gfLambdas]
            fields, native_iterations = self._native_solve_many(
                mass, rhs, tol, maxiter, warm_vecs=warm_vecs,
                mass_riesz=(solver == "native"),
                cluster_tree=(solver in {"native-batch", "native-cluster"}),
                cluster_coarse_size=(64 if solver == "native-cluster" else 0),
                cluster_deflation_size=(8 if solver == "native-cluster" else 0))
            gfm, gfls = fields[0], fields[1:]
            it_m, its = native_iterations[0], native_iterations[1:]
        elif solver == "ngsolve-cg":
            gfm = ng.GridFunction(self.fes)
            it_m = self._cg(A, pre, state_load.vec, gfm, tol, maxiter,
                            warm_vec=None if warm is None else warm.gfM.vec)
            gfls, its = [], []
            for k, load in enumerate(loads):
                gfl = ng.GridFunction(self.fes)
                wv = None if warm is None else warm.gfLambdas[k].vec
                its.append(self._cg(A, pre, load.vec, gfl, tol, maxiter,
                                    warm_vec=wv))
                gfls.append(gfl)
        else:
            raise ValueError(
                "DensityAdjointVIM.linearize: solver must be "
                "'native-batch', 'native-cluster', 'native-jacobi', 'native', or 'ngsolve-cg'")
        values = np.array([float(ng.InnerProduct(load.vec, gfm.vec))
                           for load in loads])
        jacobians = np.stack([
            -np.asarray(ng.Integrate(gfl * gfm, self.mesh,
                                     element_wise=True), dtype=float)
            for gfl in gfls])
        return FunctionalLinearization(
            values=values, jacobians=jacobians, gfM=gfm,
            gfLambdas=tuple(gfls), state_iterations=it_m,
            adjoint_iterations=tuple(its))

    def objective_and_gradient(self, s, state_load, adjoint_load,
                               tol=1e-12, maxiter=20000, warm=None,
                               solver="native-jacobi"):
        """State + adjoint solve and the full per-element gradient.

        One-functional convenience over :meth:`linearize` (same configured
        operator).  ``warm`` is an optional
        :class:`AdjointGradientResult` from the previous design iterate; its
        fields seed both CG solves.  The objective is the design-dependent
        part ``adjoint_load^T m``; the gradient is
        ``dJ/ds_e = -Integrate(lambda . m, element_wise=True)``.
        """
        warm_lin = None
        if warm is not None:
            warm_lin = FunctionalLinearization(
                values=np.array([warm.objective]),
                jacobians=warm.gradient.reshape(1, -1), gfM=warm.gfM,
                gfLambdas=(warm.gfLambda,),
                state_iterations=warm.state_iterations,
                adjoint_iterations=(warm.adjoint_iterations,))
        lin = self.linearize(s, state_load, [adjoint_load], tol=tol,
                             maxiter=maxiter, warm=warm_lin, solver=solver)
        return AdjointGradientResult(
            objective=float(lin.values[0]), gradient=lin.jacobians[0],
            gfM=lin.gfM, gfLambda=lin.gfLambdas[0],
            state_iterations=lin.state_iterations,
            adjoint_iterations=lin.adjoint_iterations[0])


# --------------------------------------------------------------------------
# constrained design loop (trust-region SLP)
# --------------------------------------------------------------------------
# Ascent ACCEPTANCE zone above the engineering band.  A hard cap at exactly
# 1.0*band starves the acceptance headroom to zero as the iterate approaches
# an active band edge: the LP targets <= band, the true (quadratic) response
# overshoots by more than the vanishing headroom at ANY move, and the trial
# backtracking exhausts move_min -- the 194-tet sector lane rejection-died at
# 5 iterates / +6.7 % (2026-08-10) against its measured +16.1 % golden.  A
# fixed overshoot zone keeps the headroom nonvanishing: ascent steps may LAND
# in (band, (1+zeta)*band] and remain in ascent mode (J stays monotone; the
# 0.9*|viol| rows already pull the next prediction back inside the band).
# DEEP Chebyshev restoration starts only ABOVE the zone -- the 40k beta=8
# stranding at 1.235 bands (see the deep trigger) is above 1.1 bands, so its
# fix is preserved, while the golden lane's accepted peak (1.06 bands) is
# inside.  0.1 deliberately threads between those two measured points.
_BAND_ACCEPT_OVERSHOOT = 0.1


@dataclass
class DensityDesignResult:
    """Converged/terminal state of :func:`optimize_density`."""
    density: np.ndarray       # raw (unfiltered) element densities in [0, 1]
    history: tuple            # one dict per ACCEPTED iterate
    converged: bool           # move limit collapsed below move_min
    final_move: float
    solves: int               # linearizations evaluated (incl. rejected trials)


def _accept_deep_restoration(violation, violation_new, band, volume_ok):
    """Chebyshev restoration acceptance for an infeasible SLP iterate.

    A non-worst row may rise while remaining below the new maximum.  Forcing
    every row to stay within 5 % of its own previous value can deadlock a
    genuine minimax step when the active worst row changes; only the maximum
    normalized violation (with an L1 tie-break) defines restoration progress.
    """
    if not volume_ok:
        return False
    normalized = np.abs(np.asarray(violation, dtype=float)) / band
    normalized_new = np.abs(np.asarray(violation_new, dtype=float)) / band
    worst = float(np.max(normalized))
    worst_new = float(np.max(normalized_new))
    total = float(np.sum(normalized))
    total_new = float(np.sum(normalized_new))
    progress = (worst_new < 0.995 * worst
                or (worst_new < worst and total_new < 0.995 * total))
    return bool(progress or np.all(normalized_new <= 1.25))


def _solve_minimax_lp_update(density, gradients, violation, band,
                             cell_volumes, raw_volume_max, *, move_limit,
                             A_ub=None, b_ub=None):
    """Solve the exact linearized Chebyshev restoration subproblem.

    The final LP variable is the common normalized violation cap ``t``.
    Unlike a high-order smooth approximation, this epigraph formulation
    continues to balance rows when two orbit radii exchange the active-worst
    role near a minimax point.
    """
    from scipy.optimize import linprog
    from .topology_optimization import LPUpdate

    rho = np.asarray(density, dtype=float).reshape(-1)
    gradients = np.asarray(gradients, dtype=float)
    violation = np.asarray(violation, dtype=float).reshape(-1)
    band = np.asarray(band, dtype=float).reshape(-1)
    volumes = np.asarray(cell_volumes, dtype=float).reshape(-1)
    if gradients.shape != (violation.size, rho.size):
        raise ValueError("minimax gradients must have shape (rows, density)")
    if band.shape != violation.shape or np.any(band <= 0.0):
        raise ValueError("minimax bands must be positive and match violations")
    if volumes.shape != rho.shape or np.any(volumes <= 0.0):
        raise ValueError("minimax cell volumes must match density and be positive")

    rows, limits = [], []

    def add_row(coeff_density, coeff_t, rhs):
        row = np.r_[np.asarray(coeff_density, dtype=float), float(coeff_t)]
        scale = 1.0 / max(float(np.max(np.abs(row))), 1e-300)
        rows.append(row*scale)
        limits.append(float(rhs)*scale)

    add_row(volumes, 0.0, raw_volume_max)
    if A_ub is not None:
        extra = np.atleast_2d(np.asarray(A_ub, dtype=float))
        rhs = np.asarray(b_ub, dtype=float).reshape(-1)
        if extra.shape != (rhs.size, rho.size):
            raise ValueError("minimax A_ub/b_ub shape mismatch")
        for row, limit in zip(extra, rhs):
            add_row(row, 0.0, limit)
    for gradient, value, width in zip(gradients, violation, band):
        scaled_gradient = gradient/width
        base = float(scaled_gradient @ rho)
        normalized_value = float(value/width)
        # value + gradient @ (x-rho) <= width*t
        add_row(scaled_gradient, -1.0, base-normalized_value)
        # -(value + gradient @ (x-rho)) <= width*t
        add_row(-scaled_gradient, -1.0, normalized_value-base)

    lower = np.maximum(0.0, rho-move_limit)
    upper = np.minimum(1.0, rho+move_limit)
    current_cap = float(np.max(np.abs(violation)/band))
    bounds = list(zip(lower, upper)) + [(0.0, current_cap)]
    objective = np.r_[np.zeros(rho.size), 1.0]
    result = linprog(objective, A_ub=np.asarray(rows),
                     b_ub=np.asarray(limits), bounds=bounds, method="highs")
    if not result.success:
        raise RuntimeError("topology minimax LP failed: %s" % result.message)

    # The minimax epigraph has only a handful of response rows and thousands
    # of density variables.  A one-stage LP therefore has a huge nullspace;
    # HiGHS may return an arbitrary bound-saturated topology with the same
    # predicted t, invalidating the local model.  A second sparse LP fixes t
    # at its optimum and minimizes the volume-weighted L1 material movement.
    # This is the lexicographic LP analogue of a proximal trust-region step.
    from scipy import sparse
    n = rho.size
    base = sparse.csr_matrix(np.asarray(rows))
    zero_base_u = sparse.csr_matrix((base.shape[0], n))
    stage2_rows = [sparse.hstack([base, zero_base_u], format="csr")]
    stage2_limits = [np.asarray(limits, dtype=float)]
    t_tolerance = 1e-8*max(1.0, abs(float(result.x[-1])))
    t_cap = sparse.csr_matrix((
        [1.0], ([0], [n])), shape=(1, 2*n+1))
    stage2_rows.append(t_cap)
    stage2_limits.append(np.asarray([float(result.x[-1])+t_tolerance]))
    identity = sparse.identity(n, format="csr")
    zero_t = sparse.csr_matrix((n, 1))
    # x-rho <= u and -(x-rho) <= u
    stage2_rows.append(sparse.hstack(
        [identity, zero_t, -identity], format="csr"))
    stage2_limits.append(rho)
    stage2_rows.append(sparse.hstack(
        [-identity, zero_t, -identity], format="csr"))
    stage2_limits.append(-rho)
    stage2_objective = np.r_[np.zeros(n+1), volumes/volumes.sum()]
    stage2_bounds = bounds + [(0.0, move_limit)]*n
    result2 = linprog(
        stage2_objective,
        A_ub=sparse.vstack(stage2_rows, format="csr"),
        b_ub=np.concatenate(stage2_limits), bounds=stage2_bounds,
        method="highs")
    if not result2.success:
        raise RuntimeError(
            "topology minimax proximal LP failed: %s" % result2.message)
    new_density = np.asarray(result2.x[:n], dtype=float)
    return LPUpdate(
        new_density, new_density-rho, float(result2.x[n]),
        str(result2.message), int(result.nit)+int(result2.nit))


def optimize_density(problem, state_load, objective_load, constraint_loads=(),
                     targets=(), *, chi_iron, volume_fraction,
                     density_filter=None, density_projection=None,
                     initial_density=None, penalty=1.0,
                     move_limit=0.1, max_iterations=30,
                     band_relative=5e-3, band_floor=None,
                     move_min=1e-3, objective_slack=1e-6,
                     tol=1e-10, cg_maxiter=20000,
                     linear_solver="native-jacobi", callback=None,
                     checkpoint_callback=None, initial_warm=None,
                     evaluation_callback=None, trial_callback=None):
    """MAXIMIZE ``J = objective_load^T m`` under linear-functional equality
    bands, an iron volume budget, and box/move limits (trust-region SLP).

    Constraints are engineering bands: ``|constraint_k - target_k| <=
    band_k`` with ``band_k = band_relative * |target_k|`` (or explicit
    ``band_floor``, absolute, per constraint).  Each iterate costs one
    weighted-mass assembly and ``1 + n_constraints``
    warm-started CG solves (:meth:`DensityAdjointVIM.linearize`), one
    element-wise Integrate per functional, and a milliseconds HiGHS LP
    (:func:`radia.topology_optimization.solve_lp_update` with the constraint
    rows in its ``A_ub`` slot, normalized to O(1) -- Tesla-scale rows sit
    below HiGHS's absolute feasibility tolerance and read as noise).

    Two-phase trust-region SLP.  While every violation is inside the
    acceptance zone (``(1 + _BAND_ACCEPT_OVERSHOOT) * band``), the LP
    maximizes J with rows targeting the band itself and acceptance keeps
    the ascent MONOTONE; accepted states may overshoot into the zone, from
    which the ``0.9 * |viol|`` rows pull the next prediction back inside
    the band.  While any violation sits ABOVE the zone (e.g. an infeasible
    profile start), the loop switches to RESTORATION: the LP objective
    becomes the exact linearized minimax band violation (J free), the rows
    hold a common epigraph cap, and acceptance requires the worst
    normalized violation to decrease.
    Rejected trials halve the step along the SAME analytic LP direction.
    When a constraint is active at the optimum the design rides its band
    edge inside the zone rather than rejection-dying against a hard cap or
    cycling between objective ascent and restoration.

    ``callback(entry)`` receives each accepted history dict.  The caller
    wraps the whole call in ``with TaskManager():``.  Informal per-iterate
    timings are recorded in the history; benchmark-grade timings belong on
    the quiet compute hosts, per the repository benchmark policy.
    """
    from .topology_optimization import LPUpdate, solve_lp_update

    constraint_loads = list(constraint_loads)
    targets = np.asarray(targets, dtype=float).ravel()
    if len(constraint_loads) != targets.size:
        raise ValueError("optimize_density: %d constraint loads but %d targets"
                         % (len(constraint_loads), targets.size))
    if not 0.0 < volume_fraction <= 1.0:
        raise ValueError("optimize_density: volume_fraction must be in (0, 1]")
    move_limit = float(move_limit)
    move_min = float(move_min)
    if (not np.isfinite(move_limit) or not np.isfinite(move_min)
            or not 0.0 < move_min <= move_limit <= 1.0):
        raise ValueError("optimize_density: need 0 < move_min <= "
                         "move_limit <= 1")
    max_iterations_i = int(max_iterations)
    if max_iterations_i < 1 or max_iterations_i != max_iterations:
        raise ValueError("optimize_density: max_iterations must be a positive integer")
    objective_slack = float(objective_slack)
    if not np.isfinite(objective_slack) or objective_slack < 0.0:
        raise ValueError("optimize_density: objective_slack must be non-negative")
    if not np.all(np.isfinite(targets)):
        raise ValueError("optimize_density: targets must be finite")
    if linear_solver not in {
            "native-batch", "native-cluster", "native-jacobi", "native", "ngsolve-cg"}:
        raise ValueError(
            "optimize_density: linear_solver must be 'native-batch', 'native-cluster', "
            "'native-jacobi', 'native', or 'ngsolve-cg'")
    volumes = problem.element_volumes
    volume_max = float(volume_fraction * volumes.sum())
    if initial_density is None:
        rho = np.full(problem.n_el, float(volume_fraction))
    else:
        rho = np.asarray(initial_density, dtype=float).copy()
        if rho.shape != (problem.n_el,):
            raise ValueError("optimize_density: initial_density shape %r != "
                             "(%d,)" % (rho.shape, problem.n_el))
    if (not np.all(np.isfinite(rho)) or np.any(rho < 0.0)
            or np.any(rho > 1.0)):
        raise ValueError("optimize_density: initial_density must be finite "
                         "and lie in [0, 1]")

    def transform_density(rho_vec):
        if density_filter is not None:
            rho_f_raw = density_filter.apply(rho_vec)
            rho_f = np.clip(rho_f_raw, 0.0, 1.0)
            unclipped = (rho_f_raw > 0.0) & (rho_f_raw < 1.0)
        else:
            rho_f = rho_vec
            unclipped = None
        if density_projection is not None:
            rho_material = density_projection.apply(rho_f)
        else:
            rho_material = rho_f
        return rho_f, rho_material, unclipped

    _rho_f_initial, rho_material_initial, _unclipped_initial = (
        transform_density(rho))
    initial_volume = float(volumes @ rho_material_initial)
    if initial_volume > volume_max + 1e-12 * max(1.0, abs(volume_max)):
        raise ValueError(
            "optimize_density: initial_density exceeds the projected "
            "iron volume budget")
    if band_floor is None:
        if not np.isfinite(band_relative) or band_relative <= 0.0:
            raise ValueError("optimize_density: band_relative must be positive")
        band = float(band_relative) * np.maximum(np.abs(targets), 1e-300)
    else:
        band = np.broadcast_to(np.asarray(band_floor, dtype=float),
                               targets.shape).astype(float)
    if targets.size and not np.all(band > 0.0):
        raise ValueError("optimize_density: constraint bands must be positive")
    loads = [objective_load] + constraint_loads

    evaluation_index = 0
    transformed_volume = (density_filter is not None
                          or density_projection is not None)
    raw_volume_max = float(volumes.sum())

    def evaluate(rho_vec, warm):
        nonlocal evaluation_index
        t_evaluate = time.perf_counter()
        # The P1 Helmholtz filter under/overshoots at bang-bang transitions
        # (measured -1.2e-2 on a coarse ball); clip the FILTERED density to
        # [0, 1] with the exact piecewise chain rule (zero derivative on
        # clipped elements) before the material map.
        rho_f, rho_material, unclipped = transform_density(rho_vec)
        lin = problem.linearize(density_to_s(
            rho_material, chi_iron, penalty=penalty),
                                state_load, loads, tol=tol,
                                maxiter=cg_maxiter, warm=warm,
                                solver=linear_solver)

        def to_raw(g_material):
            g_rf = np.asarray(g_material, dtype=float)
            if density_projection is not None:
                g_rf = density_projection.chain(rho_f, g_rf)
            if density_filter is None:
                return g_rf
            return density_filter.chain(np.where(unclipped, g_rf, 0.0))

        def to_rho(g_s):
            return to_raw(density_gradient_from_s_gradient(
                rho_material, g_s, chi_iron, penalty=penalty))

        gJ = to_rho(lin.jacobians[0])
        gks = [to_rho(lin.jacobians[1 + k])
               for k in range(len(constraint_loads))]
        material_volume = float(volumes @ rho_material)
        volume_gradient = to_raw(volumes)
        if evaluation_callback is not None:
            evaluation_callback(dict(
                evaluation=evaluation_index,
                elapsed_s=time.perf_counter() - t_evaluate,
                warm_start=warm is not None,
                state_iterations=int(lin.state_iterations),
                adjoint_iterations=[int(v) for v in lin.adjoint_iterations],
                values=np.asarray(lin.values, dtype=float).tolist(),
                iron_volume=material_volume))
        evaluation_index += 1
        return lin, gJ, gks, material_volume, volume_gradient

    lin, gJ, gks, material_volume, volume_gradient = evaluate(
        rho, initial_warm)
    n_solves = 1
    history = []
    move = move_limit
    converged = False
    for iteration in range(max_iterations_i):
        t_iter = time.perf_counter()
        J = float(lin.values[0])
        viol = lin.values[1:] - targets
        # Lexicographic SLP.  DEEP restoration owns everything ABOVE the
        # acceptance zone; the objective competes for material inside it.
        # Ascent in the former (band, 1.25*band] transition zone stranded the
        # 40k beta=8 study at 1.235 bands even though the minimax direction
        # was still productive -- 1.235 sits above the 1.1 zone, so that case
        # still deep-restores.  States inside (band, 1.1*band] stay in ascent
        # (monotone J) and the 0.9*|viol| rows below pull them back into the
        # band; see _BAND_ACCEPT_OVERSHOOT for the band-edge rejection-death
        # this zone prevents.
        deep = bool(targets.size) and bool(np.any(
            np.abs(viol) > (1.0 + _BAND_ACCEPT_OVERSHOOT) * band))

        def lp_rows(bands_eff):
            # rows normalized to O(1) for HiGHS's absolute tolerances
            A_ub, b_ub = [], []
            for G, v, b in zip(gks, viol, bands_eff):
                scale = 1.0 / max(np.abs(G).max(), 1e-300)
                base = float(G @ rho)
                A_ub.append(G * scale)
                b_ub.append((b - v + base) * scale)
                A_ub.append(-G * scale)
                b_ub.append((b + v - base) * scale)
            if transformed_volume:
                scale = 1.0 / max(np.abs(volume_gradient).max(), 1e-300)
                A_ub.append(volume_gradient * scale)
                b_ub.append((volume_max - material_volume
                             + float(volume_gradient @ rho)) * scale)
            if not A_ub:
                return None, None
            return np.array(A_ub), np.array(b_ub)

        if deep:
            lp_objective = None
        else:
            lp_objective = -gJ / max(np.abs(gJ).max(), 1e-300)

        accepted = False
        trials = 0
        band_mode = "deep" if deep else "band"
        deep_full_update = None
        deep_full_move = None
        while move >= move_min:
            trials += 1
            if deep:
                if deep_full_update is None:
                    if transformed_volume:
                        scale = 1.0 / max(
                            np.abs(volume_gradient).max(), 1e-300)
                        A_ub = np.asarray([volume_gradient*scale])
                        b_ub = np.asarray([(
                            volume_max-material_volume
                            + float(volume_gradient @ rho))*scale])
                    else:
                        A_ub = b_ub = None
                    deep_full_update = _solve_minimax_lp_update(
                        rho, gks, viol, band, volumes,
                        raw_volume_max if transformed_volume else volume_max,
                        move_limit=move, A_ub=A_ub, b_ub=b_ub)
                    deep_full_move = float(move)
                # Backtrack on one analytic minimax direction.  Re-solving
                # this highly underdetermined LP at every smaller move changes
                # the active set and can replace descent by ascent.  The fixed
                # direction was independently checked against directional FD
                # on the 40k beta=8 checkpoint (predicted/actual 1.867/1.853
                # bands at scale 1/4).
                line_scale = float(move/deep_full_move)
                delta = line_scale*deep_full_update.delta
                predicted_cap = float(np.max(
                    np.abs(viol + np.asarray(gks) @ delta)/band))
                update = LPUpdate(
                    density=rho+delta, delta=delta,
                    predicted_objective=predicted_cap,
                    status=deep_full_update.status,
                    iterations=deep_full_update.iterations)
            else:
                line_scale = 1.0
                bands_eff = np.where(np.abs(viol) > band,
                                     np.maximum(band, 0.9 * np.abs(viol)),
                                     band)
                band_mode = ("band" if not targets.size
                             or np.all(np.abs(viol) <= band) else "restore")
                A_ub, b_ub = lp_rows(bands_eff)
                try:
                    update = solve_lp_update(rho, lp_objective, volumes,
                                             raw_volume_max if transformed_volume
                                             else volume_max, move_limit=move,
                                             A_ub=A_ub, b_ub=b_ub)
                except RuntimeError:
                    bands_eff = np.maximum(band, np.abs(viol))
                    band_mode = "hold"
                    A_ub, b_ub = lp_rows(bands_eff)
                    update = solve_lp_update(rho, lp_objective, volumes,
                                             raw_volume_max if transformed_volume
                                             else volume_max, move_limit=move,
                                             A_ub=A_ub, b_ub=b_ub)
            try:
                (lin_new, gJ_new, gks_new, material_volume_new,
                 volume_gradient_new) = evaluate(update.density, warm=lin)
                n_solves += 1
            except RuntimeError as exc:
                if "CG did not converge" not in str(exc):
                    raise
                # A large material step can make the trial system much harder
                # than the accepted state.  This is a trust-region rejection,
                # not a reason to lose the checkpointed continuation.  Count
                # the attempted linearization and retry from the same analytic
                # sensitivity with a smaller move; no finite difference or
                # approximate accepted state is introduced.
                n_solves += 1
                move *= 0.5
                continue
            J_new = float(lin_new.values[0])
            viol_new = np.abs(lin_new.values[1:] - targets)
            volume_ok = material_volume_new <= (
                volume_max + 1e-10 * max(1.0, abs(volume_max)))
            if trial_callback is not None:
                trial_callback(dict(
                    move=float(move), trials=int(trials),
                    predicted_max_violation_over_band=(
                        float(update.predicted_objective) if deep else None),
                    line_search_scale=float(line_scale),
                    current_max_violation_over_band=(
                        float(np.max(np.abs(viol)/band)) if targets.size
                        else 0.0),
                    actual_max_violation_over_band=(
                        float(np.max(viol_new/band)) if targets.size else 0.0),
                    max_density_change=float(np.max(np.abs(update.delta))),
                    weighted_l1_density_change=float(
                        volumes @ np.abs(update.delta)),
                    volume_ok=bool(volume_ok)))
            if deep:
                # Accept Chebyshev feasibility progress.  A smaller L1 total
                # is not sufficient: it can sacrifice the worst orbit radius
                # while improving easier rows and never reach the all-row cap.
                # Individual easier rows are allowed to rise below the new
                # Chebyshev maximum.  A per-row 5 % guard deadlocked the 40k
                # sector continuation even while the true maximum improved
                # 8.20 -> 6.91 bands; the active worst row was simply moving
                # from the outer to the inner orbit radius.
                ok = _accept_deep_restoration(
                    viol, viol_new, band, volume_ok)
            else:
                ok_J = J_new >= J - objective_slack * abs(J)
                # Ascent may land inside the fixed acceptance zone above the
                # band (see _BAND_ACCEPT_OVERSHOOT): a hard 1.0*band cap
                # rejection-died at the active edge, while the former
                # 1.25-band cap fed accepted states to restoration and
                # cycled.  The zone is NOT a restoration trigger -- the
                # 0.9*|viol| rows pull the next prediction back into the
                # band from inside it.
                ok_g = np.all(viol_new <= (
                    1.0 + _BAND_ACCEPT_OVERSHOOT + 1e-6) * band)
                ok = volume_ok and ok_J and bool(ok_g)
            if ok:
                accepted = True
                break
            move *= 0.5
        if not accepted:
            converged = True
            break
        change = float(np.max(np.abs(update.delta)))
        rho, lin, gJ, gks = update.density, lin_new, gJ_new, gks_new
        material_volume, volume_gradient = (
            material_volume_new, volume_gradient_new)
        move = min(float(move_limit), 1.5 * move)
        entry = dict(iteration=iteration, objective=float(lin.values[0]),
                     constraints=lin.values[1:].tolist(),
                     violation=viol_new.tolist(), band=band.tolist(),
                     volume=material_volume,
                     design_volume=float(volumes @ rho),
                     max_density_change=change,
                     move=move, trials=trials, band_mode=band_mode,
                     t_iter_s=time.perf_counter() - t_iter,
                     state_iterations=lin.state_iterations,
                     adjoint_iterations=list(lin.adjoint_iterations))
        if targets.size:
            entry["max_violation_over_band"] = float(np.max(viol_new / band))
        if density_filter is not None:
            rho_report = np.clip(density_filter.apply(rho), 0.0, 1.0)
        else:
            rho_report = rho
        if density_projection is not None:
            rho_report = density_projection.apply(rho_report)
        entry.update(density_discreteness(rho_report))
        history.append(entry)
        if callback is not None:
            callback(entry)
        if checkpoint_callback is not None:
            checkpoint_callback(entry, rho.copy())
    return DensityDesignResult(density=rho, history=tuple(history),
                               converged=converged, final_move=move,
                               solves=n_solves)


# --------------------------------------------------------------------------
# Stage-3 verification protocol: exact-void iron-only re-evaluation
# --------------------------------------------------------------------------
# Boundary triangles of a netgen tet, ordered so the right-hand-rule normal
# points OUT of the element -- netgen stores surface elements outward for
# FaceDescriptor(domin=1, domout=0) (probed on an OCC sphere mesh).  The
# opposite handedness flips every surface charge and turns the demag solve
# into runaway magnetization (measured <Mz> 18-34 instead of ~2.2).
_TET_BOUNDARY_FACES = ((0, 3, 1), (1, 3, 2), (2, 3, 0), (0, 1, 2))


def iron_only_mesh(mesh, keep, *, boundary_classifier=None,
                   tetrahedralize_hex=False):
    """New straight TET/HEX/WEDGE mesh from the kept VOL elements.

    Exact void removal: vertices of the kept set are copied, kept volume
    elements are re-added as one ``iron`` material, and facets owned by exactly
    one kept element become the new exterior surface.  Surface orientation is
    determined geometrically from the owning cell, avoiding topology-specific
    Python vertex tables.  ``boundary_classifier(center, outward_normal)`` may
    return a boundary name such as ``"pole"`` or ``"fixed"`` for the later
    GetTrafo elasticity solve; its default is ``"outer"``.
    """
    import netgen.meshing as nm

    keep = np.asarray(keep, dtype=bool).ravel()
    if mesh.dim != 3:
        raise NotImplementedError("iron_only_mesh: 3D meshes only")
    if mesh.GetCurveOrder() >= 2:
        raise NotImplementedError(
            "iron_only_mesh: curved parent meshes are not supported (the "
            "extraction copies vertices only); pass the straight design mesh")
    src = mesh.ngmesh
    points = list(src.Points())
    elements = list(src.Elements3D())
    if keep.size != len(elements):
        raise ValueError("iron_only_mesh: mask has %d entries, mesh has %d "
                         "volume elements" % (keep.size, len(elements)))
    if not keep.any() or keep.all():
        raise ValueError("iron_only_mesh: the kept set must be a proper "
                         "non-empty subset (got %d of %d elements)"
                         % (int(keep.sum()), keep.size))
    # the mask lives in NGSolve VOL numbering -- verify it matches the
    # netgen element list order before using it
    for ng_el, nm_el in zip(mesh.Elements(ng.VOL), elements):
        if (tuple(sorted(v.nr for v in ng_el.vertices))
                != tuple(sorted(v.nr - 1 for v in nm_el.vertices))):
            raise RuntimeError(
                "iron_only_mesh: NGSolve VOL element numbering does not "
                "match the netgen Elements3D order on this mesh")
    new = nm.Mesh(3)
    new.SetMaterial(1, "iron")
    pmap = {}

    def pid(nr):
        if nr not in pmap:
            p = points[nr - 1].p
            pmap[nr] = new.Add(nm.MeshPoint(nm.Pnt(p[0], p[1], p[2])))
        return pmap[nr]

    ng_elements=tuple(mesh.Elements(ng.VOL));facets={}
    supported={4:"TET",6:"WEDGE",8:"HEX"};family=None
    for index in np.flatnonzero(keep):
        vs=[v.nr for v in elements[index].vertices]
        current=supported.get(len(vs))
        if current is None:
            raise NotImplementedError(
                "iron_only_mesh: only straight TET/HEX/WEDGE are supported")
        family=current if family is None else family
        if family!=current:
            raise NotImplementedError("iron_only_mesh: mixed element families are not supported")
        new.Add(nm.Element3D(1,[pid(v) for v in vs]))
        cell_coordinates=np.asarray([points[v-1].p for v in vs],dtype=float)
        cell_center=np.mean(cell_coordinates,axis=0)
        for facet in ng_elements[index].facets:
            face_vertices=[int(vertex.nr)+1 for vertex in mesh.faces[facet.nr].vertices]
            face_coordinates=np.asarray([points[v-1].p for v in face_vertices],dtype=float)
            face_center=np.mean(face_coordinates,axis=0)
            normal=np.cross(face_coordinates[1]-face_coordinates[0],
                            face_coordinates[2]-face_coordinates[0])
            if float(normal@(face_center-cell_center))<0.0:
                face_vertices=[face_vertices[0]]+list(reversed(face_vertices[1:]))
                face_coordinates=np.asarray([points[v-1].p for v in face_vertices],dtype=float)
                normal=np.cross(face_coordinates[1]-face_coordinates[0],
                                face_coordinates[2]-face_coordinates[0])
            norm=float(np.linalg.norm(normal))
            if norm==0.0:
                raise RuntimeError("iron_only_mesh: degenerate boundary facet")
            facets.setdefault(tuple(sorted(face_vertices)),[]).append(
                (tuple(face_vertices),face_center,normal/norm))
    boundary=[]
    for occurrences in facets.values():
        if len(occurrences)==1:
            boundary.append(occurrences[0])
        elif len(occurrences)!=2:
            raise RuntimeError("iron_only_mesh: facet shared by %d kept elements"
                               %len(occurrences))
    classified=[]
    for vertices,center,normal in boundary:
        name=("outer" if boundary_classifier is None else
              str(boundary_classifier(center.copy(),normal.copy())))
        if not name:
            raise ValueError("iron_only_mesh: boundary classifier returned an empty name")
        classified.append((name,vertices))
    names=sorted({name for name,_ in classified})
    descriptors={}
    for bc,name in enumerate(names,start=1):
        descriptors[name]=new.Add(nm.FaceDescriptor(
            surfnr=bc,domin=1,domout=0,bc=bc))
        new.SetBCName(bc-1,name)
    for name,vertices in classified:
        new.Add(nm.Element2D(descriptors[name],[pid(v) for v in vertices]))
    if tetrahedralize_hex:
        if family != "HEX":
            raise ValueError(
                "iron_only_mesh: tetrahedralize_hex requires a pure HEX input")
        # Netgen owns the conforming HEX -> six-TET split, including the
        # matching diagonal on every shared/boundary quadrilateral.  Calling
        # it only after exact-void extraction makes this a one-time handoff
        # from the fixed-grid Schur/TSVD stage to the topology-fixed Trafo
        # stage; it is not a design-variable interpolation or finite
        # difference approximation.
        new.Split2Tets()
    return ng.Mesh(new)


@dataclass
class IronOnlyVerification:
    """Exact-void re-evaluation of a thresholded design (Stage-3 protocol).

    ``values_*[k]`` are the functional values (state solve + load inner
    product) for the supplied builders; ``bands[k] = (embedded - iron_only)
    / |iron_only|`` quantifies the ersatz-void error of the embedded model
    honestly, per functional, at matched 0/1 representation.
    """
    keep: np.ndarray               # thresholded iron mask (parent numbering)
    iron_mesh: "ng.Mesh"
    values_embedded: np.ndarray    # 0/1 ersatz void on the parent mesh
    values_iron_only: np.ndarray   # exact void on the extracted mesh
    bands: np.ndarray
    embedded_iterations: int
    iron_only_iterations: int
    embedded_solution: object
    iron_problem: object
    iron_solution: object


def verify_design_iron_only(problem, density, state_load_builder,
                            functional_builders, *, chi_iron, threshold=0.5,
                            density_filter=None, chi_min=CHI_MIN,
                            tol=1e-10, cg_maxiter=20000, gram_kwargs=None,
                            linear_solver="native-jacobi"):
    """Stage-3 final-verification protocol for a converged density design.

    Thresholds the (filtered) density at ``threshold``, evaluates every
    functional on (a) the PARENT mesh with the 0/1 ersatz void
    (``chi_min``) and (b) a NEW iron-only mesh with the void REMOVED
    (:func:`iron_only_mesh` -- the exact-void gold standard), and reports
    the per-functional ersatz band.  Loads are mesh-bound, so the caller
    passes BUILDERS: ``state_load_builder(fes)`` and
    ``functional_builders = [f(fes) -> assembled load, ...]`` are invoked
    on both spaces (orbit points and weights are geometry-fixed).

    The comparison is at MATCHED 0/1 representation: it quantifies the
    ersatz-void/interface-charge error alone.  The separate gap between
    the CONTINUOUS design objective and the thresholded one (large for
    gray designs; drive it down with ``penalty``/projection before
    manufacturing conclusions) is the caller's report from the design
    history.  Never report final design numbers from the embedded model --
    this function exists to replace them.
    """
    density = np.asarray(density, dtype=float).ravel()
    if density.shape != (problem.n_el,):
        raise ValueError("verify_design_iron_only: density shape %r != (%d,)"
                         % (density.shape, problem.n_el))
    rho_f = (np.clip(density_filter.apply(density), 0.0, 1.0)
             if density_filter is not None else density)
    keep = rho_f >= float(threshold)
    loads = [b(problem.fes) for b in functional_builders]
    s_bin = np.where(keep, 1.0 / chi_iron, 1.0 / chi_min)
    gf_emb, it_emb = problem.solve(s_bin, state_load_builder(problem.fes),
                                   tol=tol, maxiter=cg_maxiter,
                                   solver=linear_solver)
    values_emb = np.array([float(ng.InnerProduct(load.vec, gf_emb.vec))
                           for load in loads])
    mesh_iron = iron_only_mesh(problem.mesh, keep)
    internal_interfaces = bool(
        getattr(problem.demag, "internal_interfaces", False))
    fes_iron = ng.HDiv(
        mesh_iron, order=int(problem.fes.globalorder),
        discontinuous=internal_interfaces)
    iron_gram_kwargs = dict(gram_kwargs or {})
    if internal_interfaces:
        supplied = iron_gram_kwargs.setdefault("internal_interfaces", True)
        if not supplied:
            raise ValueError(
                "verify_design_iron_only: the parent problem uses explicit "
                "internal interface charges but gram_kwargs disables them")
    problem_iron = DensityAdjointVIM(fes_iron, **iron_gram_kwargs)
    loads_iron = [b(fes_iron) for b in functional_builders]
    gf_iron, it_iron = problem_iron.solve(
        np.full(problem_iron.n_el, 1.0 / chi_iron),
        state_load_builder(fes_iron), tol=tol, maxiter=cg_maxiter,
        solver=linear_solver)
    values_iron = np.array([float(ng.InnerProduct(load.vec, gf_iron.vec))
                            for load in loads_iron])
    bands = (values_emb - values_iron) / np.maximum(np.abs(values_iron),
                                                    1e-300)
    return IronOnlyVerification(
        keep=keep, iron_mesh=mesh_iron, values_embedded=values_emb,
        values_iron_only=values_iron, bands=bands,
        embedded_iterations=it_emb, iron_only_iterations=it_iron,
        embedded_solution=gf_emb, iron_problem=problem_iron,
        iron_solution=gf_iron)
