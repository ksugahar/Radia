"""Exact symbolic reference verification; imports no Radia/FE solver.

This independently checks the declared analytic function against the ODE,
boundary data and impedance identity. It does not solve a discrete problem or
certify the measured convergence table. Scientific/reference review is separate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform

import sympy as sp


def verify():
    x = sp.Symbol("x", real=True)
    a, sigma, mu, omega, H0 = sp.symbols("a sigma mu omega H0", positive=True)
    rho = 1 / sigma
    k = sp.sqrt(sp.I * omega * mu * sigma)
    field = H0 * sp.cosh(k * x) / sp.cosh(k * a)
    identities = {}

    def zero(name, expression):
        residual = sp.simplify(expression)
        if residual != 0:
            raise AssertionError(f"{name}: actual residual={residual}; expected exact zero")
        identities[name] = {"residual": str(residual), "expected": "0", "passed": True}

    zero("dispersion_relation", k**2 - sp.I * omega * mu * sigma)
    zero("diffusion_ODE", sp.diff(field, x, 2) - sp.I * omega * mu * sigma * field)
    zero("centre_zero_flux", sp.diff(field, x).subs(x, 0))
    zero("surface_imposed_H", field.subs(x, a) - H0)
    impedance = rho * sp.diff(field, x).subs(x, a) / H0
    zero("surface_impedance", impedance - rho * k * sp.tanh(k * a))

    # Verify the product rule/ODE identity before integration, without using
    # the production assembly or its computed heating/impedance as an oracle.
    F = sp.Function("F")(x)
    density = rho * sp.diff(F, x) * sp.conjugate(sp.diff(F, x)) + sp.I * omega * mu * F * sp.conjugate(F)
    residual = sp.diff(rho * sp.diff(F, x) * sp.conjugate(F), x) - density
    zero("complex_power_flux_identity", residual.subs(sp.diff(F, x, 2), sp.I * omega * mu * sigma * F))

    conjugate = sp.conjugate(field)
    wrong_residual = sp.diff(conjugate, x, 2) - sp.I * omega * mu * sigma * conjugate
    zero("conjugate_has_opposite_ODE_sign", wrong_residual + 2 * sp.I * omega * mu * sigma * conjugate)
    wrong_surface = sp.simplify(wrong_residual.subs(x, a))
    expected_nonzero = -2 * sp.I * omega * mu * sigma * H0
    zero("wrong_phase_surface_residual_formula", wrong_surface - expected_nonzero)
    if wrong_surface == 0 or expected_nonzero.is_zero is not False:
        raise AssertionError(f"Conjugate negative control must fail: actual={wrong_surface}; expected nonzero")
    return {
        "schema": "radia.symbolic-reference-verification.v1",
        "bibliography_key": "haus1989electromagnetic",
        "method": "Exact symbolic differentiation/substitution and algebra; no numerical solver import",
        "assumptions": ["x real; a,sigma,mu,omega,H0 positive and real", "exp(+j omega t)",
                        "k principal sqrt(j omega mu sigma), with positive real part; cosh(k*a) nonzero"],
        "boundary_specialization": "Symmetric slab centre Neumann and imposed surface H; derived from the cited equation, not a claim about the book's worked-example boundary conditions",
        "identities": identities,
        "negative_control": {"reference": "conjugate analytic profile", "correct_ODE_surface_residual": str(wrong_surface),
                             "expected_nonzero": str(expected_nonzero), "killed": True},
        "software": {"sympy": sp.__version__, "python": platform.python_version()},
        "passed": True, "radia_solver_imported": False, "convergence_campaign_executed": False,
        "review_required": "Symbolic identities do not authenticate reference applicability, user approval or measured convergence",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    report = verify()
    folder = Path(__file__).resolve().parent
    report["artifacts"] = {name: hashlib.sha256((folder / name).read_bytes()).hexdigest()
                           for name in ("verify_reference.py", "reference-verification.md")}
    args.report.write_bytes((json.dumps(report, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    print(json.dumps({"passed": True, "identities": len(report["identities"]), "negative_control_killed": True}))
