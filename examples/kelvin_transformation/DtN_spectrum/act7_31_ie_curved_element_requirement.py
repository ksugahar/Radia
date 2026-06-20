# -*- coding: utf-8 -*-
"""
act7_31_ie_curved_element_requirement.py  (Act 7 -- why the IE truncation surface MUST be curved)
================================================================================================
The infinite element attaches a radial decay layer to a truncation surface Gamma and extends it
along the surface NORMAL.  So the IE is *more* sensitive to a faceted (flat-element) Gamma than an
ordinary volume FE: a facet gets the surface NORMAL wrong, and the decay rays then point the wrong
way.  This file quantifies the geometry error of a faceted vs a curved (high-order isoparametric)
truncation surface, on the cleanest separable case (a circle, the 2-D sphere), and shows:

  - faceted (order-1) surface:   NORMAL error ~ O(1/N)   (the IE-critical metric -- decay-ray tilt),
                                 perimeter error ~ O(1/N^2).
  - curved (order-p) surface:    both errors fall ~ O(1/N^p)-class -- removed at modest N.

So the IE truncation surface MUST use curved (isoparametric) elements -- this is the SAME principle
as the lab's `mesh.Curve(order)` policy ("polygon approximation of circles loses ~2% area -> ~9%
force"), but BITES HARDER for the IE because the decay direction is the surface normal (O(1/N), not
O(1/N^2)).  In the 3-D build the IE attaches to NGSolve curved (`mesh.Curve(p)`) / Cubit high-order
surface elements; the vector ends (act7_30) carry the Piola map of the curved element.

Pure numpy.  (Geometry-fidelity measurement; the modal IE artifacts act7_25-30 are geometry-EXACT --
analytic harmonics = a perfectly curved sphere -- so curved-element support is a 3-D-discretization
requirement, recorded + quantified here and in the Gate-2/3 spec.)
"""
import os
import json
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np

N_FAIL = 0


def check(name, cond, detail=""):
    global N_FAIL
    if not cond:
        N_FAIL += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}{('  -- ' + detail) if detail else ''}")


def perimeter_error(N, p, nq=24):
    """|(order-p interpolated circle perimeter) - 2 pi| for N elements."""
    si = np.linspace(0.0, 1.0, p + 1)
    gx, gw = np.polynomial.legendre.leggauss(nq)
    sg = 0.5 * (gx + 1.0)
    wg = 0.5 * gw
    dth = 2.0 * np.pi / N
    per = 0.0
    for e in range(N):
        th = e * dth + si * dth
        cx = np.polyfit(si, np.cos(th), p)
        cy = np.polyfit(si, np.sin(th), p)
        xp = np.polyval(np.polyder(cx), sg)
        yp = np.polyval(np.polyder(cy), sg)
        per += np.sum(wg * np.hypot(xp, yp))   # element arc length: x(s) on s in [0,1] spans the element
    return abs(per - 2.0 * np.pi)


def normal_max_error(N, p, nq=24):
    """max angle between the (order-p) curve normal and the TRUE radial normal (the IE decay direction)."""
    si = np.linspace(0.0, 1.0, p + 1)
    gx, _ = np.polynomial.legendre.leggauss(nq)
    sg = 0.5 * (gx + 1.0)
    dth = 2.0 * np.pi / N
    maxdev = 0.0
    for e in range(N):
        th = e * dth + si * dth
        cx = np.polyfit(si, np.cos(th), p)
        cy = np.polyfit(si, np.sin(th), p)
        xp = np.polyval(np.polyder(cx), sg)
        yp = np.polyval(np.polyder(cy), sg)
        nrm = np.hypot(xp, yp)
        cnx, cny = yp / nrm, -xp / nrm                 # curve unit normal (rotate tangent)
        thg = e * dth + sg * dth                        # true angle at the same points
        dot = np.abs(cnx * np.cos(thg) + cny * np.sin(thg))
        maxdev = max(maxdev, float(np.arccos(np.clip(dot, 0.0, 1.0)).max()))
    return maxdev


def order_fit(Ns, errs):
    return -float(np.polyfit(np.log(Ns), np.log(errs), 1)[0])


print("=" * 92)
print(" act7_31 : the IE truncation surface MUST be curved -- faceted vs curved geometry error")
print("=" * 92)

Ns = np.array([4, 8, 16, 32, 64])
RESULTS = {"N": Ns.tolist(), "perimeter": {}, "normal": {}}

print("\n  NORMAL error (the IE-critical metric -- decay rays follow the surface normal):")
print("   p     " + "  ".join(f"N={n}" for n in Ns) + "      order")
for p in (1, 2, 3):
    e = np.array([normal_max_error(n, p) for n in Ns])
    RESULTS["normal"][f"p={p}"] = e.tolist()
    print(f"   {p}   " + "  ".join(f"{x:.1e}" for x in e) + f"     ~{order_fit(Ns, e):.1f}")

print("\n  PERIMETER error:")
print("   p     " + "  ".join(f"N={n}" for n in Ns) + "      order")
for p in (1, 2, 3):
    e = np.array([perimeter_error(n, p) for n in Ns])
    RESULTS["perimeter"][f"p={p}"] = e.tolist()
    print(f"   {p}   " + "  ".join(f"{x:.1e}" for x in e) + f"     ~{order_fit(Ns, e):.1f}")

# --- the IE-critical findings ---
nrm1 = np.array([normal_max_error(n, 1) for n in Ns])
nrm3 = np.array([normal_max_error(n, 3) for n in Ns])
per1 = np.array([perimeter_error(n, 1) for n in Ns])

check("faceted (p=1) NORMAL error is O(1/N) (decay-ray tilt -- the IE-critical, slow rate)",
      0.7 < order_fit(Ns, nrm1) < 1.3, f"order {order_fit(Ns, nrm1):.2f}")
check("faceted (p=1) perimeter error is O(1/N^2) (faster than the normal -- so NORMAL dominates IE)",
      1.6 < order_fit(Ns, per1) < 2.4, f"order {order_fit(Ns, per1):.2f}")
check("curved (p=3) NORMAL error converges much faster than faceted (order p=3 > order p=1 + 1.5)",
      order_fit(Ns, nrm3) > order_fit(Ns, nrm1) + 1.5,
      f"p3 order {order_fit(Ns, nrm3):.2f} vs p1 {order_fit(Ns, nrm1):.2f}")
check("curved (p=3, N=8) NORMAL error << faceted (p=1, N=64) -- curving wins at far coarser N",
      nrm3[1] < nrm1[-1], f"p3@N8 {nrm3[1]:.1e} vs p1@N64 {nrm1[-1]:.1e}")

print("\n" + "-" * 92)
print(" CURVED-ELEMENT REQUIREMENT (for the IE):")
print("   - the IE decay layer follows the surface NORMAL -> a faceted Gamma tilts the rays at O(1/N)")
print("     (worse than the O(1/N^2) area error) -> the IE is MORE curved-element-sensitive than a")
print("     volume FE. Curving (order p) removes it ~O(1/N^p).")
print("   - => the 3-D IE MUST attach to CURVED (isoparametric) surface elements: NGSolve mesh.Curve(p)")
print("     / Cubit high-order export; the vector ends (act7_30) carry the Piola map of the curved cell.")
print("   - same principle as the lab mesh.Curve(order) policy, but it bites harder for the IE.")
print("   - (the modal artifacts act7_25-30 are geometry-EXACT; this is a 3-D-build requirement.)")
print("-" * 92)

RESULTS["normal_order_p1"] = order_fit(Ns, nrm1)
RESULTS["perimeter_order_p1"] = order_fit(Ns, per1)
RESULTS["normal_order_p3"] = order_fit(Ns, nrm3)
RESULTS["n_fail"] = N_FAIL
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "act7_31_ie_curved_element_requirement.json")
with open(out, "w") as f:
    json.dump(RESULTS, f, indent=2)
print(f"\n  wrote {os.path.basename(out)}")

print("\n" + "=" * 92)
print(" ALL CHECKS PASSED" if N_FAIL == 0 else f" {N_FAIL} CHECK(S) FAILED")
print("=" * 92)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
