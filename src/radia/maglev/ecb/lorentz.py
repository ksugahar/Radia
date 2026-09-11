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
against a direct solve of the same model and against what a centred dipole must
show: zero horizontal force, a repulsive lift, and a lift below the infinite
perfect-conductor image bound.

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

import numpy as np

MU_0 = 4 * math.pi * 1e-7


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
    mesh, lam, vecs_free, free_mask, sigma, mu, s, m_pm, z_pm_center, x_pm
):
    """Cycle-averaged force [N] on the conductor, as ``(F_x, F_y, F_z)``.

    Parameters
    ----------
    mesh : ngsolve.Mesh
        Conductor mesh whose whole boundary is labelled "outer".
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

    Returns
    -------
    (F_x, F_y, F_z) : tuple of float
        Force on the conductor.  For a dipole above the plate F_z < 0 (the
        magnet is repelled); a centred dipole has F_x = F_y = 0 by symmetry.

    Notes
    -----
    Accuracy is set by the Foster truncation, and it depends on frequency.  On
    a 30x12x3 plate mesh the lift error against a direct solve of the same
    model was 1.7 % at 500 Hz and 15 % at 5 kHz with 200 modes, falling to
    0.08 % and 0.6 % with 1600: at high frequency the basis must reach
    eigenvalues comparable to ``|s mu sigma|`` to resolve the skin depth.  Grow
    ``lam`` / ``vecs_free`` until the force stops changing.
    """
    import scipy.sparse as sp
    from ngsolve import H1, BilinearForm, GridFunction, Integrate, dx, grad

    fes = H1(mesh, order=2, dirichlet="outer")
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
    """Return conductor/PM Lorentz forces using the shared result contract.

    The Foster/NGSolve solve remains application-owned.  This adapter records
    the peak-phasor convention and both sides of Newton's third-law pair.
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
    }
