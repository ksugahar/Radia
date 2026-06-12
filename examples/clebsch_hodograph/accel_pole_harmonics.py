"""(A) Pole shape -> multipole harmonics: the equipotential design lever.

RESEARCH example (track A foundation).  The iron pole is an equipotential
Phi = const.  For a quadrupole the IDEAL pole is the hyperbola x^2-y^2 = c
(= the n=2 equipotential), giving a PURE quad.  A pole that DEVIATES from
that equipotential (a shim / a shaped tip) forces higher harmonics to keep
Phi = const along the iron face -- those are the spurious FIELD harmonics.

This quantifies the central design lever "deviation from the equipotential
= harmonics" with NO FEM: the current-free aperture field is expanded in
the quad-allowed interior multipoles

    Phi = sum_{n in {2,6,10,...}} a_n r^n cos(n theta)        (n = 2 mod 4)

and the best Phi = const along the (given) pole face is found by the
smallest singular vector of the collocation matrix.  The ideal hyperbola
-> a_6 = a_10 = 0 (machine); a shimmed pole -> a_6 grows ~linearly with the
shim amplitude.

The 2-D cross-section is the principle; the magnet ENDS are the frontier
(the equipotential SURFACE curves in 3-D -> the end chamfer; see
docs/clebsch_hodograph/DESIGN_METHODOLOGY.md sec 3.2).
"""
import math

import numpy as np

ALLOWED = (2, 6, 10)            # quad-allowed normal harmonics (n = 2 mod 4)


def _pole_face(c=0.5, y_half=0.8, n_pts=120, shim=0.0):
    """The +x quadrupole pole face x^2 - y^2 = c (one of the 4 poles; quad
    symmetry handles the rest), optionally shimmed by a Gaussian bump at the
    pole tip (y=0) of amplitude `shim` (deviation from the equipotential)."""
    y = np.linspace(-y_half, y_half, n_pts)
    x = np.sqrt(c + y * y)
    x = x + shim * np.exp(-(y / (0.35 * y_half)) ** 2)
    return x, y


def harmonics_from_pole(c=0.5, y_half=0.8, shim=0.0, allowed=ALLOWED):
    """Best Phi=const multipole fit on the (possibly shimmed) pole face.
    Returns {n: a_n / a_2} for n in `allowed` (the spurious content) + the
    relative equipotential residual."""
    x, y = _pole_face(c, y_half, shim=shim)
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    # columns: r^n cos(n th) for each allowed n, plus -1 for the constant Phi_p
    cols = [r ** n * np.cos(n * th) for n in allowed] + [-np.ones_like(r)]
    A = np.column_stack(cols)
    _, S, Vt = np.linalg.svd(A, full_matrices=False)
    v = Vt[-1]                              # smallest singular vector = best Phi=const
    a = {n: v[i] for i, n in enumerate(allowed)}
    a2 = a[allowed[0]]
    rel = {n: a[n] / a2 for n in allowed}   # normalized to the main (n=2)
    return rel, float(S[-1] / S[0])


def solve():
    """Verify: ideal hyperbola -> pure quad; shim -> a_6 grows ~linearly."""
    ideal, res_ideal = harmonics_from_pole(shim=0.0)
    spurious_ideal = max(abs(ideal[n]) for n in ALLOWED if n != 2)

    # shim sweep: a_6 should grow ~linearly with the shim amplitude
    shims = [0.0, 0.01, 0.02, 0.04]
    a6 = []
    for s in shims:
        rel, _ = harmonics_from_pole(shim=s)
        a6.append(abs(rel[6]))
    # linearity: a6/shim approximately constant for the nonzero shims
    slopes = [a6[i] / shims[i] for i in range(1, len(shims))]
    slope_spread = (max(slopes) - min(slopes)) / np.mean(slopes)

    return {
        "ideal_spurious_rel": float(spurious_ideal),   # ~0: pure quad
        "ideal_residual": float(res_ideal),            # ~0: exact equipotential
        "a6_at_shim": {f"{s:.2f}": float(v) for s, v in zip(shims, a6)},
        "a6_slope_spread": float(slope_spread),        # ~0: linear in shim
    }


def main():
    r = solve()
    print("(A) Pole shape -> multipole harmonics (equipotential design lever)\n")
    print(f"  ideal hyperbola pole: spurious |a_n/a_2| = "
          f"{r['ideal_spurious_rel']:.2e}  (-> pure quad)")
    print(f"                        equipotential residual = "
          f"{r['ideal_residual']:.2e}\n")
    print("  shim (deviation from equipotential) -> normal sextupole a_6/a_2:")
    for s, v in r["a6_at_shim"].items():
        print(f"    shim={s}: a_6/a_2 = {v:.4e}")
    print(f"  -> a_6 grows ~linearly with the shim "
          f"(slope spread {r['a6_slope_spread']:.1%})\n")
    print("  => 'iron face off the equipotential' is the harmonic source; the")
    print("     3-D ENDS extend this to the curved equipotential surface.")


if __name__ == "__main__":
    main()
