"""lorentz.py -- plate eddy-current force from a Foster expansion of the reaction field.

A z-oriented magnetic dipole of peak moment ``m`` (time-harmonic excitation,
``s = j omega``) sits at ``(x_pm, 0, z_pm)`` above a conducting plate.  Inside
the plate the z component of the reaction field, ``v = B_r,z`` in tesla,
satisfies

    (-Delta + s mu sigma) v = -s mu sigma B_z^source,   v = 0 on "outer",

which is diagonal in the Dirichlet eigenbasis of -Delta (the Foster expansion):

    v = sum_n c_n phi_n,   c_n = -s mu sigma <B_z, phi_n>_M / (lambda_n + s mu sigma).

The eddy current is the in-plane curl of that field,

    J = (1/mu) curl(v z) = (1/mu) (d_y v, -d_x v, 0),

and the cycle-averaged force on the conductor is

    <F_conductor> = (1/2) integral Re(J) x B dV,   <F_PM> = -<F_conductor>,

with the full real source field ``B = (B_x, B_y, B_z)``: the dipole phase is
the phase reference.  Keeping only the z component of the reaction field in the
curl is the usual thin-plate stream-function approximation.
``validation_test/maglev/ecb_foster_lorentz_reference.py`` locks this module
against a direct solve of the same scalar model.  The independent three-
dimensional HCurl-VIM lane in
``validation_test/maglev/ecb_foster_lorentz_3d_reference.py`` shows that this
local-reaction approximation is not quantitatively reliable for plate force.
Use the HCurl-VIM adapter below for quantitative three-dimensional work.

This is an AC-excitation model.  It has no translation velocity, so it computes
no motional drag.

Before 2026-09-11 the current was built as ``J_y = -omega sigma Im(v)``, which
reads the reaction field [T] as a vector potential: a current density in
A/m^3 with the wrong parity, zero lift for a centred magnet, and a horizontal
force of ~1900 N at 5 kHz against an 11.7 N physical bound.  The drive was also
projected as ``M_free @ Bz_free``, dropping the boundary values of B_z that the
drive -- unlike the eigenmodes -- does not vanish on.
"""
from __future__ import annotations

import math
import warnings

import numpy as np

MU_0 = 4 * math.pi * 1e-7

_SCALAR_MODEL_WARNING = (
    "The scalar local-reaction Foster model is not a quantitatively validated "
    "3-D plate-force model. It reproduces its own scalar PDE, but an independent "
    "HCurl-VIM comparison shows large force differences. Use "
    "compute_lorentz_force_via_hcurl_vim for quantitative 3-D work."
)


def pm_field_z_dipole(x, y, z, m, x_pm=0.0, y_pm=0.0, z_pm=0.0):
    """B_z(x,y,z) from a magnetic dipole m at (x_pm, y_pm, z_pm) aligned with z."""
    dx = x - x_pm
    dy = y - y_pm
    dz = z - z_pm
    r2 = dx*dx + dy*dy + dz*dz
    r = np.sqrt(r2 + 1e-30)
    r5 = r**5
    return (MU_0 / (4*math.pi)) * (3*m*dz*dz/r5 - m/(r**3))


def pm_field_xz_dipole(x, y, z, m, x_pm=0.0, y_pm=0.0, z_pm=0.0):
    """B_x(x,y,z) component from z-aligned dipole (for force in x direction)."""
    dx = x - x_pm
    dy = y - y_pm
    dz = z - z_pm
    r2 = dx*dx + dy*dy + dz*dz
    r5 = (r2 + 1e-30)**(5/2.0)
    return (MU_0 / (4*math.pi)) * 3*m*dx*dz/r5


def _dipole_source_field(m_pm, z_pm_center, x_pm):
    """(B_x, B_y, B_z) CoefficientFunctions of a z dipole at (x_pm, 0, z_pm)."""
    from ngsolve import x as xC
    from ngsolve import y as yC
    from ngsolve import z as zC

    dz = zC - z_pm_center
    r2 = (xC - x_pm)**2 + yC**2 + dz**2 + 1e-30
    k = MU_0 / (4*math.pi)
    bx = k * 3*m_pm*(xC - x_pm)*dz / r2**2.5
    by = k * 3*m_pm*yC*dz / r2**2.5
    bz = k * (3*m_pm*dz*dz / r2**2.5 - m_pm / r2**1.5)
    return bx, by, bz


def compute_lorentz_force_via_foster(
    mesh,
    lam,
    vecs_free,
    free_mask,
    sigma,
    mu,
    s,
    m_pm,
    z_pm_center,
    x_pm,
    *,
    dirichlet_label="outer",
    warn_model_limit=True,
):
    """Cycle-averaged force [N] on the conductor, as ``(F_x, F_y, F_z)``.

    Parameters
    ----------
    mesh : ngsolve.Mesh
        Conductor mesh.
    lam : ndarray (n_eigen,)
        Eigenvalues of -Lap on the conductor with Dirichlet BC on "outer".
    vecs_free : ndarray (n_free, n_eigen)
        M-normalized eigenvectors on the free DOFs of ``H1(order=2)``, e.g. from
        ``radia.maglev.mixed_galerkin.alpha._dirichlet_eigenmodes``.
    free_mask : ndarray (n_dof,) bool
        Mask for free (non-Dirichlet) DOFs.
    sigma : float
        Conductivity [S/m].
    mu : float
        Conductor permeability [H/m].
    s : complex
        Laplace variable of the excitation, ``j omega``.
    m_pm : float
        Peak dipole moment [A m^2], oriented along +z.
    z_pm_center : float
        z-coordinate of the dipole.
    x_pm : float
        x-coordinate of the dipole.
    dirichlet_label : str
        Boundary region on which the scalar reaction field vanishes.  The
        default preserves the historical all-boundary model.  A thin plate can
        instead name only its lateral boundary, but that boundary correction
        does not resolve the independently measured 3-D model discrepancy.

    Returns
    -------
    (F_x, F_y, F_z) : tuple of float
        Force on the conductor.  For a dipole above the plate F_z < 0 (the
        magnet is repelled); a centred dipole has F_x = F_y = 0 by symmetry.

    Notes
    -----
    Accuracy within the scalar PDE is set by the Foster truncation, and it
    depends on frequency.  On
    a 30x12x3 plate mesh the lift error against a direct solve of the same
    model was 1.7 % at 500 Hz and 15 % at 5 kHz with 200 modes, falling to
    0.08 % and 0.6 % with 1600: at high frequency the basis must reach
    eigenvalues comparable to ``|s mu sigma|`` to resolve the skin depth.  Grow
    ``lam`` / ``vecs_free`` until the force stops changing.  Prefer
    :func:`compute_lorentz_force_via_foster_verified`, which checks the
    truncation against a direct scalar solve and falls back rather than
    silently accepting a high-frequency error.

    This convergence does not establish three-dimensional physical accuracy.
    The independent HCurl-VIM lane
    ``validation_test/maglev/ecb_foster_lorentz_3d_reference.py`` rejects this
    local scalar model for quantitative force prediction on its plate case.
    """
    import scipy.sparse as sp
    from ngsolve import H1, BilinearForm, GridFunction, Integrate, dx, grad

    if warn_model_limit:
        warnings.warn(_SCALAR_MODEL_WARNING, RuntimeWarning, stacklevel=2)

    fes = H1(mesh, order=2, dirichlet=dirichlet_label)
    bx_cf, by_cf, bz_cf = _dipole_source_field(m_pm, z_pm_center, x_pm)

    bz_gfu = GridFunction(fes)
    bz_gfu.Set(bz_cf)
    bz_full = np.array(bz_gfu.vec.FV().NumPy())

    u, w = fes.TnT()
    m_form = BilinearForm(fes, symmetric=True)
    m_form += u * w * dx
    m_form.Assemble()
    rows_m, cols_m, vals_m = m_form.mat.COO()
    mass = sp.csr_matrix(
        (np.asarray(vals_m), (np.asarray(rows_m), np.asarray(cols_m))),
        shape=(fes.ndof, fes.ndof),
    )

    # The eigenmodes vanish on "outer" but the drive does not: project with the
    # free rows against ALL mass-matrix columns.
    projection = vecs_free.T @ (mass[free_mask, :] @ bz_full)
    s_mu_sigma = s * mu * sigma
    coefficients = -s_mu_sigma * projection / (lam + s_mu_sigma)

    v_full = np.zeros(fes.ndof, dtype=complex)
    v_full[free_mask] = vecs_free @ coefficients
    v_re = GridFunction(fes)
    v_re.vec.FV().NumPy()[:] = v_full.real

    # Re(J) = (1/mu) curl(Re(v) z); B is real, so <F> = 0.5 int Re(J) x B.
    gradient = grad(v_re)
    jx = gradient[1] / mu
    jy = -gradient[0] / mu
    force_x = 0.5 * Integrate(jy * bz_cf, mesh)
    force_y = 0.5 * Integrate(-jx * bz_cf, mesh)
    force_z = 0.5 * Integrate(jx * by_cf - jy * bx_cf, mesh)
    return float(force_x), float(force_y), float(force_z)


def compute_lorentz_force_via_direct_scalar(
    mesh,
    sigma,
    mu,
    s,
    m_pm,
    z_pm_center,
    x_pm,
    *,
    dirichlet_label="outer",
):
    """Directly solve the same scalar local-reaction PDE and return force.

    This is a truncation oracle for the Foster expansion, not an independent
    physical model.  It is intentionally exposed so high-frequency callers can
    distinguish an insufficient eigenbasis from the much larger 3-D-model
    discrepancy documented by the HCurl-VIM validation lane.
    """

    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    from ngsolve import H1, BilinearForm, GridFunction, Integrate, LinearForm, dx, grad

    fes = H1(mesh, order=2, dirichlet=dirichlet_label)
    bx_cf, by_cf, bz_cf = _dipole_source_field(m_pm, z_pm_center, x_pm)
    trial, test = fes.TnT()
    stiffness_form = BilinearForm(fes, symmetric=True)
    stiffness_form += grad(trial) * grad(test) * dx
    mass_form = BilinearForm(fes, symmetric=True)
    mass_form += trial * test * dx
    drive_form = LinearForm(fes)
    drive_form += bz_cf * test * dx
    stiffness_form.Assemble()
    mass_form.Assemble()
    drive_form.Assemble()

    def csr(form):
        rows, columns, values = form.mat.COO()
        return sp.csr_matrix(
            (np.asarray(values), (np.asarray(rows), np.asarray(columns))),
            shape=(fes.ndof, fes.ndof),
        )

    free = np.asarray(
        [bool(fes.FreeDofs()[index]) for index in range(fes.ndof)], dtype=bool
    )
    s_mu_sigma = complex(s) * float(mu) * float(sigma)
    system = (csr(stiffness_form) + s_mu_sigma * csr(mass_form))[free][:, free]
    values = np.zeros(fes.ndof, dtype=complex)
    drive = np.asarray(drive_form.vec.FV().NumPy())
    values[free] = spla.spsolve(system.tocsc(), -s_mu_sigma * drive[free])
    reaction_real = GridFunction(fes)
    reaction_real.vec.FV().NumPy()[:] = values.real
    gradient = grad(reaction_real)
    jx = gradient[1] / mu
    jy = -gradient[0] / mu
    return (
        float(0.5 * Integrate(jy * bz_cf, mesh)),
        float(0.5 * Integrate(-jx * bz_cf, mesh)),
        float(0.5 * Integrate(jx * by_cf - jy * bx_cf, mesh)),
    )


def compute_lorentz_force_via_foster_verified(
    mesh,
    lam,
    vecs_free,
    free_mask,
    sigma,
    mu,
    s,
    m_pm,
    z_pm_center,
    x_pm,
    *,
    dirichlet_label="outer",
    rtol=0.01,
    atol_N=1.0e-9,
    fallback_to_direct=True,
):
    """Check Foster truncation against the direct scalar solve.

    The returned dictionary states whether the supplied eigenbasis met the
    requested tolerance and which backend supplied ``force_N``.  When the
    tolerance is missed, the default is a direct scalar fallback; setting
    ``fallback_to_direct=False`` raises instead.  Neither outcome upgrades the
    scalar PDE to a validated 3-D plate model.
    """

    tolerance = float(rtol)
    absolute = float(atol_N)
    if not math.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("rtol must be finite and >= 0")
    if not math.isfinite(absolute) or absolute < 0.0:
        raise ValueError("atol_N must be finite and >= 0")
    warnings.warn(_SCALAR_MODEL_WARNING, RuntimeWarning, stacklevel=2)
    foster = np.asarray(
        compute_lorentz_force_via_foster(
            mesh,
            lam,
            vecs_free,
            free_mask,
            sigma,
            mu,
            s,
            m_pm,
            z_pm_center,
            x_pm,
            dirichlet_label=dirichlet_label,
            warn_model_limit=False,
        ),
        dtype=float,
    )
    direct = np.asarray(
        compute_lorentz_force_via_direct_scalar(
            mesh,
            sigma,
            mu,
            s,
            m_pm,
            z_pm_center,
            x_pm,
            dirichlet_label=dirichlet_label,
        ),
        dtype=float,
    )
    difference = float(np.linalg.norm(foster - direct))
    scale = max(float(np.linalg.norm(direct)), absolute)
    relative = difference / scale if scale > 0.0 else 0.0
    converged = difference <= absolute + tolerance * float(np.linalg.norm(direct))
    if not converged and not fallback_to_direct:
        raise RuntimeError(
            "Foster force did not converge to the direct scalar solve: "
            f"relative difference {relative:.6g} exceeds rtol={tolerance:.6g}"
        )
    selected = foster if converged else direct
    return {
        "schema": "radia.maglev-foster-force-verification/v1",
        "force_N": selected.tolist(),
        "backend": "foster" if converged else "direct_scalar_fallback",
        "foster_force_N": foster.tolist(),
        "direct_scalar_force_N": direct.tolist(),
        "foster_mode_count": int(np.asarray(lam).size),
        "absolute_difference_N": difference,
        "relative_difference": relative,
        "rtol": tolerance,
        "atol_N": absolute,
        "foster_converged": bool(converged),
        "model_scope": "scalar_local_reaction_only",
        "validated_against_3d": False,
    }


def compute_lorentz_force_via_hcurl_vim(
    model,
    current_basis,
    magnetic_flux_density_phasor_T,
    s,
    drive=1.0,
):
    """Return the quantitative 3-D HCurl-VIM conductor force in newtons.

    ``model`` is a :class:`radia.vim.HCurlEddyCLNModel`; ``current_basis`` is
    its sampled three-dimensional divergence-free current basis, and ``B`` is
    the incident-field phasor sampled at the same quadrature points.  The
    open-boundary VIM interaction and model construction remain solver-owned;
    this adapter centralizes the current reconstruction and common peak-phasor
    Lorentz convention.
    """

    force, _torque = compute_lorentz_force_torque_via_hcurl_vim(
        model,
        current_basis,
        magnetic_flux_density_phasor_T,
        s,
        drive,
    )
    return force


def compute_lorentz_force_torque_via_hcurl_vim(
    model,
    current_basis,
    magnetic_flux_density_phasor_T,
    s,
    drive=1.0,
    *,
    pivot_m=None,
):
    """Return quantitative 3-D HCurl-VIM force and torque.

    Torque is evaluated from the same sampled Lorentz density about
    ``pivot_m`` (the global origin by default), so force and torque cannot drift
    to different quadrature or phasor conventions.
    """

    from radia.force import integrate_time_average_lorentz_force_and_torque

    coefficients = np.asarray(model.solve_vector_potential_drive(s, drive))
    if coefficients.ndim != 1:
        raise ValueError("model coefficients must be one-dimensional")
    modes = np.asarray(current_basis.modes)
    if modes.ndim != 3 or modes.shape[0] != coefficients.size or modes.shape[2] != 3:
        raise ValueError("current_basis.modes must have shape (n_modes, n_samples, 3)")
    current = np.einsum("a,aik->ik", coefficients, modes)
    force, torque = integrate_time_average_lorentz_force_and_torque(
        current,
        magnetic_flux_density_phasor_T,
        current_basis.weights,
        current_basis.points,
        pivot_m=pivot_m,
        amplitude="peak",
    )
    return (
        tuple(float(value) for value in force),
        tuple(float(value) for value in torque),
    )


def compute_lorentz_force_result_via_hcurl_vim(
    model,
    current_basis,
    magnetic_flux_density_phasor_T,
    s,
    drive=1.0,
    *,
    frame="global_cartesian",
    pivot_m=None,
):
    """Return the 3-D HCurl-VIM conductor/source force pair.

    The source force is the reaction on the exciting permanent magnet or coil.
    It is reported explicitly so MagLev consumers use the same action/reaction
    and peak-phasor conventions as Motor and the rest of ``radia.force``.
    """

    from radia.force import force_torque_result

    conductor_force_raw, conductor_torque_raw = (
        compute_lorentz_force_torque_via_hcurl_vim(
            model,
            current_basis,
            magnetic_flux_density_phasor_T,
            s,
            drive,
            pivot_m=pivot_m,
        )
    )
    conductor_force = list(conductor_force_raw)
    conductor_torque = list(conductor_torque_raw)
    source_force = [-component for component in conductor_force]
    source_torque = [-component for component in conductor_torque]
    common = {
        "method": "hcurl_vim_time_average_lorentz_body_force",
        "frame": frame,
        "pivot_m": pivot_m,
        "field_convention": "time_average_phasor",
        "amplitude": "peak",
    }
    return {
        "schema": "radia.maglev-force-pair/v1",
        "conductor": force_torque_result(
            conductor_force, conductor_torque, **common
        ),
        "source_pm": force_torque_result(source_force, source_torque, **common),
        "action_reaction_residual_N": [
            conductor_force[index] + source_force[index] for index in range(3)
        ],
        "action_reaction_torque_residual_Nm": [
            conductor_torque[index] + source_torque[index] for index in range(3)
        ],
        "validated_against_3d": True,
    }


def compute_lorentz_force_result_via_foster(
    mesh,
    lam,
    vecs_free,
    free_mask,
    sigma,
    mu,
    s,
    m_pm,
    z_pm_center,
    x_pm,
    *,
    frame="global_cartesian",
    pivot_m=None,
):
    """Return scalar-model conductor/PM forces using the shared contract.

    This compatibility adapter records the peak-phasor convention and both
    sides of Newton's third-law pair.  It emits the scalar-model warning; use
    :func:`compute_lorentz_force_result_via_hcurl_vim` for quantitative 3-D
    prediction.
    """

    from radia.force import force_torque_result

    conductor_force = list(compute_lorentz_force_via_foster(
        mesh,
        lam,
        vecs_free,
        free_mask,
        sigma,
        mu,
        s,
        m_pm,
        z_pm_center,
        x_pm,
    ))
    source_force = [-component for component in conductor_force]
    common = {
        "method": "time_average_lorentz_body_force",
        "frame": frame,
        "pivot_m": pivot_m,
        "field_convention": "time_average_phasor",
        "amplitude": "peak",
    }
    return {
        "schema": "radia.maglev-force-pair/v1",
        "conductor": force_torque_result(conductor_force, None, **common),
        "source_pm": force_torque_result(source_force, None, **common),
        "action_reaction_residual_N": [
            conductor_force[index] + source_force[index] for index in range(3)
        ],
        "validated_against_3d": False,
    }
