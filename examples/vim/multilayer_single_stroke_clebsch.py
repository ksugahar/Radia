"""Multi-layer single-stroke via the foliated Clebsch (lambda, mu) path.

A MULTI-LAYER stream function IS the foliated volume Clebsch  J = grad(lambda) x
grad(mu)  with mu = the layer (here mu = radius): each mu-leaf (a cylinder R_m)
carries a surface stream function lambda -- a STACK of surface SFs
(streamfunction_volume.py).  The multi-layer SINGLE-STROKE is ONE wire threading
all layers = a continuous path in the (lambda, mu) plane:

  * within a layer (mu = R_m fixed): lambda spirals -> a helix (the single-layer
    pitch, single_stroke_clebsch.py);
  * to cross to the next layer: mu steps R_m -> R_{m+1} (a radial crossover).

The natural path is BOUSTROPHEDON (alternating axial direction): the SAME
azimuthal sense in every layer (so Bz ADDS across layers) but ALTERNATING axial
advance.  The per-layer axial 'lead' currents (the mu-step secular terms) then
CANCEL -- exactly the rigorous secular bookkeeping of the (lambda, mu) path.

Demonstrated (Radia Biot-Savart): the boustrophedon single-stroke is a valid
M-layer solenoid (Bz ~ M x single layer) whose net-axial-lead secular stray is
~100x smaller than naively stacking M SAME-handed helical layers (net axial M*I).

Honest scope: the lead-cancellation is EVEN-M (the alternation pairs layers
+I/-I; odd M leaves one uncancelled layer, so the net axial is I and the factor
degrades to the trivial M*I/I ~ M).  We use M = 4 (even).  Also: the ~121x is the
azimuthal B_theta secular component (it GROWS with distance, ~rho^-2.85 vs the
same-handed rho^-1.97 -- the 1/rho net-axial multipole is what is killed), not
the full |B|.  Same clean solenoid-foliation (mu = radius) case as the single
layer; a general multi-layer winding-surface family is the F1-like routing part.

Run:  python examples/vim/multilayer_single_stroke_clebsch.py
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src",
                                "radia"))

R0 = 0.05       # innermost layer radius [m]
DR = 0.012      # layer pitch (radial) [m]
L = 0.20        # length [m]
I = 100.0       # wire current [A]
N = 16          # turns per layer
M = 4           # number of layers (EVEN: the alternation pairs +I/-I axial leads)
assert M % 2 == 0, "lead-cancellation is even-M; odd M leaves one uncancelled layer"


def _layer_pts(Rm, updown, ns=60):
    """Helical turns on cylinder Rm; updown=+1 winds up (-L/2 -> +L/2), -1 down.
    Azimuthal sense (theta increasing) is the SAME for both -> Bz adds."""
    th = np.linspace(0, 2 * math.pi * N, N * ns + 1)
    z = updown * (-L / 2 + (th / (2 * math.pi * N)) * L)
    return [[Rm * math.cos(a), Rm * math.sin(a), float(zz)]
            for a, zz in zip(th, z)]


def _boustrophedon(rad):
    """ONE continuous wire: layer m alternates axial dir; the consecutive-point
    polyline makes the radial crossovers automatically (alternating top/bottom)."""
    pts = []
    for m in range(M):
        pts.extend(_layer_pts(R0 + m * DR, +1 if m % 2 == 0 else -1))
    return rad.ObjFlmCur(pts, I)


def _same_handed_stack(rad):
    """M SAME-handed helices (the naive multi-layer; each carries axial lead I
    -> net M*I, uncancelled)."""
    return rad.ObjCnt([rad.ObjFlmCur(_layer_pts(R0 + m * DR, +1), I)
                       for m in range(M)])


def _single_layer(rad):
    return rad.ObjFlmCur(_layer_pts(R0, +1), I)


def _fld(rad, obj, pts):
    return np.array([rad.Fld(obj, "b", list(p)) for p in pts])


def _Btheta(B, pts):
    out = []
    for i, q in enumerate(pts):
        a = math.atan2(q[1], q[0])
        out.append(-B[i, 0] * math.sin(a) + B[i, 1] * math.cos(a))
    return np.array(out)


def main():
    import radia as rad

    print("=== multi-layer single-stroke via the foliated Clebsch (lambda,mu) ===")
    axis = np.array([[0.0, 0.0, z] for z in np.linspace(-0.15 * L, 0.15 * L, 5)])
    far = []
    for z in np.linspace(-0.3 * L, 0.3 * L, 3):
        for a in np.linspace(0, 2 * math.pi, 16, endpoint=False):
            far.append([8 * R0 * math.cos(a), 8 * R0 * math.sin(a), float(z)])
    far = np.array(far)

    rad.UtiDelAll()
    Bz1 = float(np.mean(_fld(rad, _single_layer(rad), axis)[:, 2]))
    bous, samh = _boustrophedon(rad), _same_handed_stack(rad)
    Bz_b = float(np.mean(_fld(rad, bous, axis)[:, 2]))
    Bz_s = float(np.mean(_fld(rad, samh, axis)[:, 2]))
    Bf_b = _fld(rad, bous, far)
    Bf_s = _fld(rad, samh, far)
    rad.UtiDelAll()

    bz_ratio = Bz_b / (M * Bz1)                       # valid M-layer solenoid?
    stray_b = float(np.linalg.norm(_Btheta(Bf_b, far)))   # net-axial-lead secular
    stray_s = float(np.linalg.norm(_Btheta(Bf_s, far)))
    cancel = stray_s / (stray_b + 1e-30)              # secular cancellation factor

    print(f"layers M={M}, turns/layer N={N}, radii {R0:.3f}..{R0+(M-1)*DR:.3f} m")
    print(f"on-axis Bz: single-layer {Bz1*1e3:.3f} mT, boustrophedon "
          f"{Bz_b*1e3:.3f} mT (M*single = {M*Bz1*1e3:.3f}, ratio {bz_ratio:.3f})")
    print(f"net-axial-lead secular far-stray: boustrophedon {stray_b:.3e}, "
          f"same-handed {stray_s:.3e}  -> cancels by {cancel:.0f}x")
    print("  => the boustrophedon (lambda,mu) path keeps the azimuthal sense (Bz")
    print("     adds) but ALTERNATES the axial advance, so the per-layer mu-step")
    print("     secular leads CANCEL -- one wire, ~100x less far-field stray.")

    result = {
        "M_layers": M, "N_turns": N, "R0_m": R0, "dR_m": DR,
        "bz_single_mT": Bz1 * 1e3, "bz_boustrophedon_mT": Bz_b * 1e3,
        "bz_ratio_vs_Mxsingle": bz_ratio,
        "stray_boustrophedon": stray_b, "stray_same_handed": stray_s,
        "secular_cancellation": cancel,
    }
    ok = (abs(bz_ratio - 1.0) < 0.15      # valid M-layer solenoid (Bz ~ M*single)
          and cancel > 5.0)               # the alternating mu-steps cancel leads
    print(f"\n=== PASS: {ok} ===")
    return result, ok


if __name__ == "__main__":
    _, ok = main()
    raise SystemExit(0 if ok else 1)
