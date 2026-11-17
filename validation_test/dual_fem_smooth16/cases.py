"""Predeclared manufactured cases for the smooth periodic A/Omega rate study.

This table is fixed before any solve.  Each family has a 2D member on the unit
2-torus and a genuinely z-dependent 3D member on the unit 3-torus.  Fields are
pure trigonometric, so phi* and B* = curl A* have zero mean.

Truth construction (both dimensions, no material jump, no boundary):
    H* = -grad phi*,  B* = curl A*  (2D: rot psi* = (d_y psi*, -d_x psi*)),
    S  = B* + mu grad phi*,  so that  B* = mu H* + S  and  curl H* = 0, div B* = 0.

Descriptive rate rules, declared before any solve.  They are not acceptance
criteria, and equal convergence order is a hypothesis under test:
  * h ladder (fixed p): "equal h-rate" if the observed field-norm orders of A and
    Omega on the finest pair differ by at most 0.15.
  * p ladder (fixed mesh): fit ln(field error) against p over the non-floored levels
    with p >= 2; "equal p-decay" if the fitted decay constants differ by at most
    10% of the larger one.
  * Case verdict: "equal" if both hold, "mixed" if one holds, "different" if
    neither holds, "undetermined" if a ladder has fewer than two usable levels.
  * Floor: a level is excluded from fits when its relative field error is below
    1e-11, or when its readout-order kappa check fails (see smooth16.py).
Fixed-weight estimates: the midpoint (A+Omega)/2 is primary; 2:1 (2A+Omega)/3 is
the predeclared comparison.  No weight is fitted.
"""

from math import pi

U = 2 * pi

# Ladders (fixed before any solve).
P_LADDER = {2: {"n": 4, "orders": (1, 2, 3, 4, 5, 6)},
            3: {"n": 3, "orders": (1, 2, 3, 4)}}
H_LADDER = {2: {"order": 2, "n": (4, 8, 16, 32)},
            3: {"order": 2, "n": (3, 6, 12)}}
ASSEMBLY_BONUS = 8

# Reference-first blend test (user procedure, 2026-11-08), fixed before any solve:
#   1. Q_ref = midpoint of the converged bracket [Q_A, Q_Omega] at REFERENCE, with
#      its half-width as uncertainty.  Q_A and Q_Omega are the energy functionals
#      C -/+ 2J, which bound Q* for any conforming discrete field, so this is the
#      corrected interval.  The last p step gives the drift of each formulation;
#      the h-ladder finest level is a secondary check.
#   2. Freeze Q_ref.  Then, at COARSE, compare Q_A, Q_Omega, (Q_A+Q_Omega)/2 and
#      (2Q_A+Q_Omega)/3 against it.  A blend is "closer" only if the difference of
#      absolute errors exceeds twice the reference half-width.  Otherwise the
#      result is "inconclusive".
#   3. The analytic Q* is an independent check on Q_ref only.
COARSE = {2: {"n": 4, "p": 2}, 3: {"n": 3, "p": 2}}
REFERENCE = {2: {"n": 4, "p": 6}, 3: {"n": 3, "p": 4}}
READOUT_EXTRA = (12, 18)  # readout order = 2p + extra, two levels for kappa stability


def _families(ng):
    sin, cos, exp = ng.sin, ng.cos, ng.exp
    x, y, z = ng.x, ng.y, ng.z
    return [
        dict(key="F1", label="constant mu, separable modes",
             d2=dict(mu=ng.CF(1.0), phi=sin(U * x), psi=cos(U * y)),
             d3=dict(mu=ng.CF(1.0), phi=sin(U * (x + z)),
                     A=(sin(U * z), ng.CF(0.0), cos(U * y)))),
        dict(key="F2", label="constant mu, oblique modes",
             d2=dict(mu=ng.CF(1.0), phi=sin(U * (x + y)), psi=cos(U * (x - 2 * y))),
             d3=dict(mu=ng.CF(1.0), phi=sin(U * (x + y + z)),
                     A=(cos(U * (y - z)), sin(U * (x + 2 * z)), ng.CF(0.0)))),
        dict(key="F3", label="low contrast (about 2) axis-aligned mu",
             d2=dict(mu=exp(0.35 * sin(U * x)), phi=sin(U * y), psi=sin(U * x)),
             d3=dict(mu=exp(0.35 * sin(U * z)), phi=sin(U * y) * cos(U * z),
                     A=(ng.CF(0.0), sin(U * x), cos(U * (x + z))))),
        dict(key="F4", label="high contrast (about 20) oblique mu",
             d2=dict(mu=exp(1.5 * sin(U * (x + y))), phi=sin(U * x) * cos(U * y),
                     psi=sin(U * (2 * x + y))),
             d3=dict(mu=exp(1.5 * sin(U * (x + y + z))), phi=sin(U * x) * cos(U * z),
                     A=(sin(U * (y + 2 * z)), ng.CF(0.0), cos(U * (2 * x + y))))),
        dict(key="F5", label="product-form mu (contrast 3), mixed modes",
             d2=dict(mu=1 + 0.5 * cos(U * x) * cos(U * y), phi=cos(U * (2 * x + y)),
                     psi=sin(U * x) * sin(U * y)),
             d3=dict(mu=1 + 0.5 * cos(U * x) * cos(U * y) * cos(U * z),
                     phi=cos(U * (2 * x + y - z)),
                     A=(sin(U * y) * sin(U * z), cos(U * x) * sin(U * z), ng.CF(0.0)))),
        dict(key="F6", label="two-mode mu (contrast about 16), two-mode phi",
             d2=dict(mu=exp(0.7 * sin(U * x) + 0.7 * cos(U * y)),
                     phi=sin(U * (x - y)) + 0.3 * cos(2 * U * y), psi=cos(U * (x + y))),
             d3=dict(mu=exp(0.7 * sin(U * x) + 0.7 * cos(U * (y + z))),
                     phi=sin(U * (x - y)) + 0.3 * cos(2 * U * z),
                     A=(ng.CF(0.0), cos(U * (x + z)), sin(U * (x - y))))),
        dict(key="F7", label="short-wave phi, contrast 3 mu",
             d2=dict(mu=2 + sin(U * (x + 2 * y)), phi=sin(U * (2 * x + 3 * y)),
                     psi=cos(U * (x - y))),
             d3=dict(mu=2 + sin(U * (x + 2 * y + z)), phi=sin(U * (2 * x + 3 * y + z)),
                     A=(cos(U * (y - z)), ng.CF(0.0), sin(U * (x + y))))),
        dict(key="F8", label="short-wave high-contrast (about 20) mu",
             d2=dict(mu=exp(1.5 * cos(2 * U * x)), phi=sin(U * y) + 0.5 * sin(U * x),
                     psi=cos(U * (x + 2 * y))),
             d3=dict(mu=exp(1.5 * cos(2 * U * z)), phi=sin(U * y) + 0.5 * sin(U * (x + z)),
                     A=(cos(U * (y + 2 * z)), sin(U * (z + x)), cos(U * (x + 2 * y))))),
    ]


def cases(ng):
    """Return the 16 cases as (case_id, dim, family_label, fields)."""
    out = []
    for family in _families(ng):
        out.append((family["key"] + "-2D", 2, family["label"], family["d2"]))
        out.append((family["key"] + "-3D", 3, family["label"], family["d3"]))
    return out
