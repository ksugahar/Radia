"""lorentz.py -- Linear ECB drag/lift via Lorentz integral F = integral J x B dV.

For a PM (cube Nd or dipole) at position (x_pm, 0, z_pm) moving along
+x at velocity v, this computes the time-averaged drag and lift forces
on the plate using the Lorentz integral over the conductor volume:

    <F_conductor> = (1/2) Re[ integral conj(J(s)) x B_ext(r) dV ]
    <F_PM> = -<F_conductor>   (Newton's 3rd)

J(r, s) is reconstructed from the Foster eigenmode expansion of the
scalar A_z model.  For the full vector eddy-current formulation see
calc_fem_kelvin.py in src/radia/panels/.

Why Lorentz instead of Maxwell stress: in our framework Lorentz has
a SHORTER error chain than Maxwell stress (no Biot-Savart of induced
field, no surface-evaluation artifacts).  See package README for
details.
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


def compute_lorentz_force_via_foster(
    mesh, lam, vecs_free, free_mask, sigma, mu, s, m_pm, z_pm_center, x_pm
):
    """Time-averaged drag (F_x) and lift (F_z) on the conductor.

    Parameters
    ----------
    mesh : ngsolve.Mesh
    lam : ndarray (n_eigen,)
        Eigenvalues of -Lap on the conductor with Dirichlet BC.
    vecs_free : ndarray (n_free, n_eigen)
        M-normalized eigenvectors on free DOFs.
    free_mask : ndarray (n_dof,) bool
        Mask for free (non-Dirichlet) DOFs.
    sigma : float
        Conductivity.
    mu : float
        Permeability.
    s : complex
        Laplace variable (typically j*omega).
    m_pm : float
        PM magnetic moment magnitude (A m^2).
    z_pm_center : float
        z-coordinate of PM dipole center.
    x_pm : float
        x-coordinate of PM dipole center.

    Returns
    -------
    F_x_conductor : float
        Time-averaged x-force on the conductor (drag opposes motion).
    F_z_conductor : float
        Time-averaged z-force on the conductor (lift on PM is opposite).
    """
    import scipy.sparse as sp
    from ngsolve import (
        H1,
        BilinearForm,
        GridFunction,
        Integrate,
        dx,
    )
    from ngsolve import x as xC
    from ngsolve import y as yC
    from ngsolve import z as zC

    fes = H1(mesh, order=2, dirichlet="outer")
    Bz_cf = (MU_0 / (4*math.pi)) * (
        3*m_pm*(zC - z_pm_center)*(zC - z_pm_center) /
        ((xC-x_pm)**2 + yC**2 + (zC-z_pm_center)**2 + 1e-30)**(5/2.0)
        - m_pm /
        ((xC-x_pm)**2 + yC**2 + (zC-z_pm_center)**2 + 1e-30)**(3/2.0)
    )
    Bx_cf = (MU_0 / (4*math.pi)) * (
        3*m_pm*(xC-x_pm)*(zC - z_pm_center) /
        ((xC-x_pm)**2 + yC**2 + (zC-z_pm_center)**2 + 1e-30)**(5/2.0)
    )

    Bz_gfu = GridFunction(fes)
    Bz_gfu.Set(Bz_cf)
    Bz_vec_full = np.array(Bz_gfu.vec.FV().NumPy())
    Bz_vec_free = Bz_vec_full[free_mask]

    u, v = fes.TnT()
    m_form = BilinearForm(fes, symmetric=True)
    m_form += u * v * dx
    m_form.Assemble()
    rows_m, cols_m, vals_m = m_form.mat.COO()
    ndof = fes.ndof
    M = sp.csr_matrix(
        (np.asarray(vals_m), (np.asarray(rows_m), np.asarray(cols_m))),
        shape=(ndof, ndof),
    )
    M_free = M[free_mask][:, free_mask]

    M_Bz = M_free.dot(Bz_vec_free)
    proj_Bz = vecs_free.T @ M_Bz

    s_mu_sigma = s * mu * sigma
    c_n = -s_mu_sigma * proj_Bz / (lam + s_mu_sigma)

    v_vec_free = vecs_free @ c_n
    v_vec_full = np.zeros(fes.ndof, dtype=complex)
    v_vec_full[free_mask] = v_vec_free

    v_gfu_re = GridFunction(fes)
    v_gfu_im = GridFunction(fes)
    v_gfu_re.vec.FV().NumPy()[:] = v_vec_full.real
    v_gfu_im.vec.FV().NumPy()[:] = v_vec_full.imag

    omega = s.imag
    Jy_re = -omega * sigma * v_gfu_im   # Re[J_y] = omega sigma Im[v]

    Fx_conductor = -0.5 * Integrate(Jy_re * Bz_cf, mesh)
    Fz_conductor = +0.5 * Integrate(Jy_re * Bx_cf, mesh)

    return float(Fx_conductor), float(Fz_conductor)


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

    The underlying Foster/NGSolve solve remains application-owned. This
    adapter records the peak-phasor convention and both sides of Newton's
    third-law pair without changing the legacy tuple-returning function.
    """

    from radia.force import force_torque_result

    force_x, force_z = compute_lorentz_force_via_foster(
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
    )
    conductor_force = [force_x, 0.0, force_z]
    source_force = [-force_x, 0.0, -force_z]
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
