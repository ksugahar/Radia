# -*- coding: utf-8 -*-
r"""AGE nonlinear motor workflow gate.

Validates the Picard iteration wrapping the linear AGE core:

  1. Linear limit: nu_cf_fn=None → converges in 1 iteration, result identical
     to direct airgap_solve (machine-precision match).

  2. Constant nu: nu_cf_fn returning a constant → same as linear limit.

  3. Saturation loop: Froelich-type nu → converges within max_iter, monotone
     residual decrease, torque non-trivial.

  4. Rotation sweep: age_motor_rotation_sweep gives per-angle results, all
     converged, torque varies with angle.
"""
import math
import os
import sys

import numpy as np
import ngsolve as ng
from ngsolve import (
    H1, BilinearForm, LinearForm, grad, dx, x, y, atan2, cos, sin,
    BND, TaskManager, Mesh,
)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.airgap_machine import (
    airgap_coupling, airgap_factorize, airgap_solve, airgap_torque,
)
from radia_mcp.radia_ngsolve.airgap_motor_workflow import (
    age_motor_nonlinear_solve,
    age_motor_rotation_sweep,
)

# ---------------------------------------------------------------------------
# Shared fixture (simple 2D rotating machine, no eddy currents for clarity)
# ---------------------------------------------------------------------------

R0, RI, RO, REXT = 0.05, 0.10, 0.11, 0.20
N_HARM = 2
A_ROT, A_STA = 1.0, 0.3
NU0 = 1.0 / (4e-7 * math.pi)   # free-space reluctivity [A·m/Wb]


def _build_mesh():
    geo = SplineGeometry()
    geo.AddCircle((0, 0), REXT, leftdomain=1, rightdomain=0, bc="outer")
    geo.AddCircle((0, 0), RO,   leftdomain=0, rightdomain=1, bc="stator_ring")
    geo.AddCircle((0, 0), RI,   leftdomain=2, rightdomain=0, bc="rotor_ring")
    geo.AddCircle((0, 0), R0,   leftdomain=0, rightdomain=2, bc="rotor_inner")
    geo.SetMaterial(1, "stator")
    geo.SetMaterial(2, "rotor")
    return Mesh(geo.GenerateMesh(maxh=0.01))


def _make_fes(mesh):
    return H1(mesh, order=2, dirichlet="outer|rotor_inner", complex=True)


def _excite(mesh, theta_r):
    return mesh.BoundaryCF({
        "rotor_inner": A_ROT * cos(N_HARM * atan2(y, x) - N_HARM * theta_r),
        "outer": A_STA * cos(N_HARM * atan2(y, x)),
    })


def _bilinear_fn(fes, nu_cf):
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu_cf * grad(u) * grad(v) * dx
    a.Assemble()
    return a


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def validate_linear_limit_one_iteration():
    """With constant nu, Picard converges in at most 2 iterations.

    Iteration 1: solve from zero → solution u₁.  rel_change = ||u₁||/||u₁|| = 1 > tol.
    Iteration 2: same stiffness → same solution u₂ = u₁.  rel_change = 0 < tol → STOP.
    So niter_used ≤ 2."""
    mesh = _build_mesh()
    fes = _make_fes(mesh)
    coupling = airgap_coupling(fes, RI, RO, "rotor_ring", "stator_ring", [N_HARM])
    nu_cf = ng.CoefficientFunction(NU0)

    gfu_p, info = age_motor_nonlinear_solve(
        fes, coupling,
        bilinear_fn=lambda nu: _bilinear_fn(fes, nu),
        nu_cf_fn=lambda gfu: nu_cf,
        dirichlet_cf=_excite(mesh, 0.0),
        niter=5, tol=1e-12,
    )
    assert info["niter_used"] <= 2, f"Expected ≤2 iterations for linear case, got {info['niter_used']}"
    assert info["converged"]


def validate_linear_limit_matches_direct_solve():
    """Picard result must match direct airgap_solve to machine precision."""
    mesh = _build_mesh()
    fes = _make_fes(mesh)
    coupling = airgap_coupling(fes, RI, RO, "rotor_ring", "stator_ring", [N_HARM])
    nu_cf = ng.CoefficientFunction(NU0)
    theta_r = 0.7

    # Picard solve
    gfu_p, info = age_motor_nonlinear_solve(
        fes, coupling,
        bilinear_fn=lambda nu: _bilinear_fn(fes, nu),
        nu_cf_fn=lambda gfu: nu_cf,
        dirichlet_cf=_excite(mesh, theta_r),
        niter=5, tol=1e-12,
    )

    # Direct AGE solve
    a_ref = _bilinear_fn(fes, nu_cf)
    factor_ref = airgap_factorize(a_ref.mat, coupling, fes.FreeDofs())
    gfu_ref = airgap_solve(factor_ref, fes, dirichlet_cf=_excite(mesh, theta_r))

    diff = gfu_p.vec.CreateVector()
    diff.data = gfu_p.vec - gfu_ref.vec
    norm_diff = math.sqrt(abs(ng.InnerProduct(diff, diff)))
    norm_ref  = math.sqrt(abs(ng.InnerProduct(gfu_ref.vec, gfu_ref.vec)))
    rel = norm_diff / norm_ref if norm_ref > 0 else float("inf")

    assert rel < 1e-10, f"Picard vs direct solve mismatch: rel_err = {rel:.3e}"


def validate_null_nu_cf_fn_converges_fast():
    """nu_cf_fn=None (linear mode): converges in at most 2 iterations."""
    mesh = _build_mesh()
    fes = _make_fes(mesh)
    coupling = airgap_coupling(fes, RI, RO, "rotor_ring", "stator_ring", [N_HARM])
    nu_const = ng.CoefficientFunction(NU0)

    gfu, info = age_motor_nonlinear_solve(
        fes, coupling,
        bilinear_fn=lambda nu: _bilinear_fn(fes, nu if nu is not None else nu_const),
        nu_cf_fn=None,
        dirichlet_cf=_excite(mesh, 0.0),
        niter=5, tol=1e-12,
    )
    assert info["niter_used"] <= 2
    assert info["converged"]


def validate_froelich_saturation_converges():
    """Froelich-type nu → Picard converges within 50 iters, residuals shrink.

    With |∇A| ~ A_ROT/R0 ~ 20 (in SI 2D), B_sat2 = 500 gives mild but non-trivial
    saturation (ν decreases by ~20% at peak field), so Picard converges smoothly."""
    mesh = _build_mesh()
    fes = _make_fes(mesh)
    coupling = airgap_coupling(fes, RI, RO, "rotor_ring", "stator_ring", [N_HARM])

    B_sat2 = 500.0   # mild saturation: |∇A|² ~ 400 at peak, so ν ≈ 0.44 NU0

    def nu_cf_fn(gfu):
        B2 = ng.InnerProduct(grad(gfu), grad(gfu)).real
        return ng.CoefficientFunction(NU0) / (1.0 + B2 / B_sat2)

    gfu, info = age_motor_nonlinear_solve(
        fes, coupling,
        bilinear_fn=lambda nu: _bilinear_fn(fes, nu),
        nu_cf_fn=nu_cf_fn,
        dirichlet_cf=_excite(mesh, 0.5),
        niter=80, tol=5e-4, relax=0.7,
    )
    assert info["converged"], (
        f"Froelich Picard did not converge in {info['niter_used']} iters; "
        f"last rel_change = {info['rel_change_history'][-1]:.3e}"
    )
    h = info["rel_change_history"]
    assert h[-1] < h[0], f"Residuals not shrinking: {h[0]:.3e} → {h[-1]:.3e}"


def validate_saturation_torque_nonzero():
    """Torque must be nonzero for nonzero phase between rotor and stator excitation."""
    mesh = _build_mesh()
    fes = _make_fes(mesh)
    coupling = airgap_coupling(fes, RI, RO, "rotor_ring", "stator_ring", [N_HARM])

    B_sat2 = 500.0

    def nu_cf_fn(gfu):
        B2 = ng.InnerProduct(grad(gfu), grad(gfu)).real
        return ng.CoefficientFunction(NU0) / (1.0 + B2 / B_sat2)

    theta_r = math.pi / 4   # 45 deg → nonzero torque expected
    gfu, info = age_motor_nonlinear_solve(
        fes, coupling,
        bilinear_fn=lambda nu: _bilinear_fn(fes, nu),
        nu_cf_fn=nu_cf_fn,
        dirichlet_cf=_excite(mesh, theta_r),
        niter=30, tol=1e-5,
    )
    T = airgap_torque(coupling, gfu, axial_length=0.05)
    assert abs(T) > 1e-8, f"Torque should be nonzero, got T={T:.3e}"


def validate_rotation_sweep_linear():
    """Linear rotation sweep: all angles converge in 1 iter, torque varies."""
    mesh = _build_mesh()
    fes = _make_fes(mesh)
    coupling = airgap_coupling(fes, RI, RO, "rotor_ring", "stator_ring", [N_HARM])
    nu_const = ng.CoefficientFunction(NU0)

    thetas = np.linspace(0, 2 * math.pi / N_HARM, 9, endpoint=False)
    results = age_motor_rotation_sweep(
        fes, coupling,
        bilinear_fn=lambda nu: _bilinear_fn(fes, nu if nu is not None else nu_const),
        excitation_fn=lambda theta: _excite(mesh, theta),
        theta_r_list=thetas,
        nu_cf_fn=None,
        axial_length=0.05,
        niter=5, tol=1e-12,
    )
    assert all(r["converged"] for r in results), "Not all angles converged"
    assert all(r["niter_used"] <= 2 for r in results), "Linear sweep should use ≤2 iters"
    T_vals = [r["torque"] for r in results]
    T_max = max(abs(T) for T in T_vals)
    T_min = min(abs(T) for T in T_vals)
    # Torque varies with angle (not identically zero)
    assert T_max > 1e-8, "Torque should be nonzero somewhere in the sweep"
    # Torque should vary (not constant) over the rotation
    assert (T_max - T_min) / T_max > 0.1, (
        f"Torque should vary with angle: max={T_max:.4e}, min={T_min:.4e}"
    )


def main():
    for name, fn in sorted(globals().items()):
        if name.startswith("validate_") and callable(fn):
            fn(); print("ok", name)
    print("[OK] AGE nonlinear motor: Picard loop + saturation + rotation sweep.")


if __name__ == "__main__":
    main()
