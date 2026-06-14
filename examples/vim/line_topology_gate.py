"""Stage-B line-topology gate: do the current lines CLOSE (windable) or not?

The SECOND topological extraction gate (the first is B1 helicity_diagnostic.py;
F3 helical_current_no_clebsch.py demonstrated the phenomenon for the trivial
single-turn azimuthal line).  Even a helicity-0 current gives clean equal-current
WIRES only if its field lines CLOSE: a current line is windable as one loop iff it
returns to itself.  On a winding surface the line's ROTATION NUMBER iota = dphi/dtheta
decides this:

  * iota = p/q RATIONAL  -> the line closes after q toroidal turns -> a windable
    q-fold loop  (q = 1 is the azimuthal ring; q > 1 is a (p,q) torus winding);
  * iota IRRATIONAL      -> the line never closes; it covers an invariant curve
    densely -> NO finite equal-current loop family can be cut from it.

This packages a reusable diagnostic -- closure_defect(Jf, seed, ...) traces ONE
line and returns, per toroidal turn, how far it is from its seed; the MINIMUM over
the first Q_max turns (relative to the minor radius) is the gate score.  Small ->
the line closes (extractable); bounded-away-from-zero -> it does not.

This is the ROBUST closure metric (cf. F3's "Clebsch ring closes, ABC line drifts
143 ring-lengths"), NOT a chaotic-Poincare 2-D-fill test -- the ABC flow has mixed
KAM islands + chaotic sea, so a single-seed 2-D-fill score is unreliable; closure
defect of a prescribed-rotation-number winding is clean and seed-robust.  The gate
takes ANY callable J(p), so it applies to a mesh-evaluated current too; the torus
winding with an exact constant iota is the analytically-clean demonstration vehicle.

Run:  python examples/vim/line_topology_gate.py
"""
from __future__ import annotations

import math

import numpy as np

R_MAJOR = 1.0     # torus major radius
A_MINOR = 0.3     # torus minor radius (the transverse closure scale)
TWO_PI = 2.0 * math.pi
GOLDEN_INV = 2.0 / (1.0 + math.sqrt(5.0))     # 1/phi, the "most irrational" iota


def torus_winding_field(iota, R0=R_MAJOR, a=A_MINOR):
    """Tangent field of a torus winding with EXACT constant rotation number iota =
    dphi/dtheta.  J = (R0 + a cos phi) e_theta + iota a e_phi (toroidal + poloidal),
    tangent to the (R0,a) torus everywhere."""
    def J(p):
        px, py, pz = p
        th = math.atan2(py, px)
        eR = np.array([math.cos(th), math.sin(th), 0.0])
        rho = p - R0 * eR
        phi = math.atan2(pz, float(np.dot(rho, eR)))
        e_th = np.array([-math.sin(th), math.cos(th), 0.0])
        e_ph = np.array([-math.sin(phi) * math.cos(th),
                         -math.sin(phi) * math.sin(th), math.cos(phi)])
        return (R0 + a * math.cos(phi)) * e_th + iota * a * e_ph
    return J


def closure_defect(Jf, seed, max_turns=6, steps_per_turn=400):
    """RK4-trace a field line; report |p - seed| at each completed toroidal turn.
    Returns (defects: list[(turn, dist)], min_rel: min dist / minor radius)."""
    p = np.array(seed, float)
    p0 = p.copy()
    dt = TWO_PI / steps_per_turn
    th_prev = math.atan2(p[1], p[0])
    total_th = 0.0
    mark = 1
    defects = []
    for _ in range(max_turns * steps_per_turn * 2):
        k1 = Jf(p)
        k2 = Jf(p + 0.5 * dt * k1)
        k3 = Jf(p + 0.5 * dt * k2)
        k4 = Jf(p + dt * k3)
        p = p + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        th = math.atan2(p[1], p[0])
        dth = th - th_prev
        if dth < -math.pi:
            dth += TWO_PI
        if dth > math.pi:
            dth -= TWO_PI
        total_th += dth
        th_prev = th
        if total_th >= mark * TWO_PI:
            defects.append((mark, float(np.linalg.norm(p - p0))))
            mark += 1
            if mark > max_turns:
                break
    min_rel = min(d for _, d in defects) / A_MINOR if defects else float("inf")
    return defects, min_rel


def gate(Jf, seed, max_turns=6, thresh=0.15):
    """Line-topology verdict.  closes (windable) iff some turn returns within
    thresh*minor-radius of the seed."""
    defects, min_rel = closure_defect(Jf, seed, max_turns=max_turns)
    q = min((d for d in defects if d[1] / A_MINOR < thresh),
            key=lambda d: d[1], default=(None, None))[0]
    return {"min_rel_defect": min_rel, "closes_at_turn": q,
            "closed": bool(min_rel < thresh)}


def main():
    print("=== Stage-B line-topology gate (rotation-number closure defect) ===")
    seed = (R_MAJOR + A_MINOR, 0.0, 0.0)        # on the outer equator of the torus

    cases = [("iota = 1/5  (rational)", 0.2),
             ("iota = 2/5  (rational)", 0.4),
             ("iota = 1/phi (irrational)", GOLDEN_INV)]
    res = {}
    for tag, iota in cases:
        g = gate(torus_winding_field(iota), seed)
        res[tag] = g
        verdict = (f"CLOSES at turn {g['closes_at_turn']} -> windable"
                   if g["closed"] else "never closes -> NOT windable")
        print(f"  {tag:26s}: min closure-defect {g['min_rel_defect']:.3f} (minor-radii)"
              f"  -> {verdict}")

    rational_ok = all(res[t]["closed"] for t, _ in cases[:2])
    irrational_open = not res[cases[2][0]]["closed"]
    margin = res[cases[2][0]]["min_rel_defect"] / max(
        res[cases[0][0]]["min_rel_defect"], 1e-9)
    print(f"  separation margin (irrational/rational min-defect) = {margin:.1f}x")
    print("  => rational-rotation lines close after q turns (windable q-fold loops);")
    print("     irrational lines never close -> no finite equal-current loop family.")

    ok = (rational_ok and irrational_open and margin > 3.0)
    print(f"\n=== PASS: {ok} ===")
    return res, float(margin), ok


if __name__ == "__main__":
    *_, ok = main()
    raise SystemExit(0 if ok else 1)
