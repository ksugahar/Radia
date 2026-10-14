"""
10_edge_corner_dofs.py -- must face, edge and corner get their own surface DOF?

The tensor envelope psi_1 = F(x) F(y) F(z), F = 1 - g, g the planar layer,
carries the face layers, the edge products g g and the corner product g g g
with FIXED relative weights (they are the expansion of one product).  The
talk's cube row stops at 0.33 % in the transition band for that reason, or
so the summary slide says.  This script measures it against the exact
heat-content admittance of the box (_references/box_heat_content.py) by
adding CONFORMING, LOCALISED surface functions with their own amplitudes:

    d(x)  = g(x; t) - g(x; a t)      a layer bump: 0 on the faces and 0 in the bulk
    psi_F = sum_axes  d F F          face-profile correction    (codim 1)
    psi_E = sum_pairs d d F          edge-localised amplitude   (codim 2)
    psi_C = d d d                    corner-localised amplitude (codim 3)
    psi_E2                           psi_E at rate 2a (edge shape)

plus the odd Dirichlet sine modes as the bulk.  Everything is a tensor product
of 1-D factors, so every 3-D integral is a product of 1-D Gauss sums on a rule
graded to the skin depth (checked against the closed forms of cube3d/03 to
1e-13).  At low frequency the bumps and F are all ~ x(L - x); the family is
then linearly dependent while the Galerkin value is not, so the Jacobi-scaled
system is solved by truncated-SVD least squares (rcond 1e-12; 1e-10 and 1e-13
give the same numbers to five digits).

Earlier attempts in this directory split the envelope into NON-conforming
pieces (08: f(x) alone is not zero on the y and z faces; asymptotically
rank-deficient) or used non-separable wedge bumps (09: diverged).  The bumps
here are conforming and localised, which is why the family behaves.

Result (2026-09-04, exact reference, 1 Hz .. 1 GHz, 46 points):
    cube, envelope only:  0.71 % (2 unknowns) .. 0.21 % (126 unknowns);
                          the floor sits in the transition band 6-60 kHz and
                          the asymptotic-band error (0.045 % at 1 MHz) does not
                          move with the bulk rank
    + corner DOF:         0.16 %       + face-profile DOF: 0.20 % (nothing)
    + EDGE DOF:           0.014 %      (15x; 1 MHz 0.002 %, 100 MHz 0.0002 %)
    + edge + corner:      0.0125 %, and 0.030 % with only 11 unknowns
    square:               0.34 % -> 0.006 % with the corner (= 2-D edge) DOF
So the edge amplitude is what the single coefficient gets wrong, at every
frequency; what remains after it is a transition-band residual that falls
with the bulk rank (0.108 % at 4 unknowns -> 0.0125 % at 128), i.e. the
coupling with the interior modes.

Is an INDEPENDENT edge DOF needed, or would a better fixed envelope do?  The
free 3-DOF solve gives amplitude ratios xi_E/xi_1 -> 0.800 and xi_C/xi_1 ->
0.457 (real) above ~1 MHz, but complex and swinging through the transition
band (e.g. 0.50 - 0.19 j and 0.39 - 0.59 j at 4 kHz).  Freezing the ratios at
their 1 GHz values into ONE surface function reproduces the free result above
1 MHz exactly, and in the transition band only when the bulk is rich (n <= 9:
0.020 % vs 0.0125 %); with a small bulk it fails there (n <= 3: 0.21 % vs
0.030 %; n <= 1: 0.74 % vs 0.11 %).  In the asymptotic band a fixed envelope
with the right ratios suffices; in the transition band the edge and corner
amplitudes must be free, because they couple to the interior modes.
"""
import cmath
import itertools
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _references.box_heat_content import Y_exact, self_test  # noqa: E402
from _references.cube3d_foster import Y_DC_cube3d  # noqa: E402

SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
L = 5e-3
MS = MU * SIGMA
PI = math.pi
ALPHA = 2.0


# --------------------------------------------------------------------------
# 1-D machinery
# --------------------------------------------------------------------------
def graded_rule(L=L, a_min_frac=1e-7, ratio=1.5, max_panel_frac=1 / 24, nodes=10):
    """Symmetric Gauss-Legendre rule on [0, L], panels graded geometrically
    towards both ends down to a_min_frac L (0.5 nm: the 1 GHz layer is 2 um)."""
    xg, wg = np.polynomial.legendre.leggauss(nodes)
    edges = [0.0, a_min_frac * L]
    while edges[-1] < L / 2:
        step = min(edges[-1] * (ratio - 1.0), max_panel_frac * L)
        edges.append(min(edges[-1] + step, L / 2))
    half = np.array(edges)
    full = np.concatenate([half, L - half[-2::-1]])
    x, w = [], []
    for a, b in zip(full[:-1], full[1:]):
        x.append(0.5 * (b - a) * xg + 0.5 * (b + a))
        w.append(0.5 * (b - a) * wg)
    return np.concatenate(x), np.concatenate(w)


def layer(x, t, L=L):
    """g = cosh(t(x - L/2)) / cosh(tL/2) and g', overflow-free for Re t > 0."""
    e1 = np.exp(-t * x)
    e2 = np.exp(-t * (L - x))
    den = 1.0 + np.exp(-t * L)
    return (e1 + e2) / den, t * (e2 - e1) / den


class Factors:
    """The 1-D factor functions on the rule and their Gram, stiffness and load
    tables: M[p, q] = int u_p u_q, S[p, q] = int u_p' u_q', b[p] = int u_p
    (bilinear, no conjugation: the system is complex symmetric)."""

    def __init__(self, t, x, w, alphas=(ALPHA, 2 * ALPHA), n_sines=(1, 3, 5, 7, 9)):
        self.names, vals, ders = [], [], []
        g, dg = layer(x, t)
        self.names.append("F")
        vals.append(1.0 - g)
        ders.append(-dg)
        for a in alphas:
            ga, dga = layer(x, a * t)
            self.names.append(f"d{a:g}")
            vals.append(g - ga)
            ders.append(dg - dga)
        for n in n_sines:
            k = n * PI / L
            self.names.append(f"s{n}")
            vals.append(np.sin(k * x) + 0j)
            ders.append(k * np.cos(k * x) + 0j)
        V, dV = np.array(vals), np.array(ders)
        self.M = (V * w[None, :]) @ V.T
        self.S = (dV * w[None, :]) @ dV.T
        self.b = V @ w
        self.idx = {n: i for i, n in enumerate(self.names)}


def make_basis(D, kind, n_bulk_max, alpha=ALPHA):
    """[(name, [(coef, factor-name tuple), ...])] for the D-box.  kind is a
    string of flags: '1' envelope, 'F' face, 'E' edge, 'C' corner, 'G' edge at
    rate 2 alpha."""
    a, a2 = f"d{alpha:g}", f"d{2 * alpha:g}"
    axes = list(range(D))
    basis = []
    odd = [n for n in (1, 3, 5, 7, 9) if n <= n_bulk_max]
    modes = sorted(itertools.product(odd, repeat=D), key=lambda m: sum(n * n for n in m))
    for m in modes:
        basis.append(("bulk" + "".join(str(n) for n in m), [(1.0, tuple(f"s{n}" for n in m))]))
    if "1" in kind:
        basis.append(("psi_1", [(1.0, tuple("F" for _ in axes))]))
    if "F" in kind:
        basis.append(("psi_F", [(1.0, tuple(a if j == i else "F" for j in axes)) for i in axes]))
    if "E" in kind and D >= 2:
        basis.append(("psi_E", [(1.0, tuple(a if j in pair else "F" for j in axes))
                                for pair in itertools.combinations(axes, 2)]))
    if "C" in kind:
        basis.append(("psi_C", [(1.0, tuple(a for _ in axes))]))
    if "G" in kind and D >= 2:
        basis.append(("psi_E2", [(1.0, tuple(a2 if j in pair else "F" for j in axes))
                                 for pair in itertools.combinations(axes, 2)]))
    return basis


def assemble(basis, fac, D, sMS):
    n = len(basis)
    K = np.zeros((n, n), dtype=complex)
    Mm = np.zeros((n, n), dtype=complex)
    load = np.zeros(n, dtype=complex)
    ix = fac.idx
    for i, (_, ti) in enumerate(basis):
        for ci, fi in ti:
            load[i] += ci * np.prod([fac.b[ix[f]] for f in fi])
        for j in range(i, n):
            tj = basis[j][1]
            mij = kij = 0j
            for ci, fi in ti:
                for cj, fj in tj:
                    m1 = [fac.M[ix[fi[ax]], ix[fj[ax]]] for ax in range(D)]
                    s1 = [fac.S[ix[fi[ax]], ix[fj[ax]]] for ax in range(D)]
                    mij += ci * cj * np.prod(m1)
                    for ax in range(D):
                        others = np.prod([m1[k] for k in range(D) if k != ax]) if D > 1 else 1.0
                        kij += ci * cj * s1[ax] * others
            Mm[i, j] = Mm[j, i] = mij
            K[i, j] = K[j, i] = kij
    return K + sMS * Mm, -sMS * load, load


def Y_galerkin(s, D, kind, n_bulk_max, rule, alpha=ALPHA, surface=None):
    """Mixed Galerkin admittance of the D-box with the given basis.  `surface`
    replaces the surface part of the basis by an explicit list (frozen-ratio
    test).  Returns (Y, n_unknowns, {name: coefficient})."""
    x, w = rule
    sMS = s * MS
    fac = Factors(cmath.sqrt(sMS), x, w, alphas=(alpha, 2 * alpha))
    basis = make_basis(D, kind, n_bulk_max, alpha)
    if surface is not None:
        basis = [b for b in basis if b[0].startswith("bulk")] + list(surface)
    A, rhs, load = assemble(basis, fac, D, sMS)
    dsc = 1.0 / np.sqrt(np.abs(np.diag(A)))
    As = dsc[:, None] * A * dsc[None, :]
    xi = dsc * np.linalg.lstsq(As, dsc * rhs, rcond=1e-12)[0]
    Y = Y_DC_cube3d(L, SIGMA) * L ** (D - 3) * (1.0 + (xi @ load) / L**D)
    return Y, len(basis), {b[0]: complex(c) for b, c in zip(basis, xi)}


def frozen_ratio_surface(D, rule, f_hz=1e9, alpha=ALPHA):
    """ONE surface function psi_1 + r_E psi_E + r_C psi_C with the amplitude
    ratios frozen at their asymptotic (1 GHz) values from the free solve: a
    better fixed envelope, no independent edge or corner DOF.  Whether it
    matches the free 3-DOF surface tells whether independence is needed."""
    _, _, xi = Y_galerkin(2j * PI * f_hz, D, "1EC", 9, rule, alpha)
    r_E = xi["psi_E"] / xi["psi_1"]
    r_C = xi["psi_C"] / xi["psi_1"]
    a = f"d{alpha:g}"
    axes = list(range(D))
    terms = [(1.0, tuple("F" for _ in axes))]
    terms += [(r_E, tuple(a if j in pair else "F" for j in axes)) for pair in itertools.combinations(axes, 2)]
    terms += [(r_C, tuple(a for _ in axes))]
    return [("psi_frozen", terms)], r_E, r_C


CONFIGS = [
    ("envelope", "1", 1), ("envelope", "1", 3), ("envelope", "1", 5), ("envelope", "1", 9),
    ("envelope + face", "1F", 9),
    ("envelope + corner", "1C", 9),
    ("envelope + edge", "1E", 9),
    ("envelope + edge + corner", "1EC", 1), ("envelope + edge + corner", "1EC", 3),
    ("envelope + edge + corner", "1EC", 5), ("envelope + edge + corner", "1EC", 9),
    ("envelope + face + edge + edge2 + corner", "1FEGC", 9),
]
FREQS = np.logspace(0, 9, 46)
REPORT_AT = (1e3, 1e4, 1e5, 1e6, 1e8)


def _row(label, kind, nb, n, ys, ex):
    err = np.abs(np.array(ys) - ex) / np.abs(ex)
    i = int(np.argmax(err))
    return {
        "surface": label, "flags": kind, "bulk_n_max": nb, "n_unknowns": int(n),
        "max_error_pct": float(100 * err[i]), "at_f_Hz": float(FREQS[i]),
        "error_pct_at": {f"{f:.0e}": float(100 * err[int(np.argmin(np.abs(FREQS - f)))]) for f in REPORT_AT},
    }


def _pair(z):
    return [float(z.real), float(z.imag)]


def sweep(D, rule):
    ex = np.array([Y_exact(2j * PI * f, D, L, SIGMA) for f in FREQS])
    rows, ratios = [], []
    for label, kind, nb in CONFIGS:
        sols = [Y_galerkin(2j * PI * f, D, kind, nb, rule) for f in FREQS]
        rows.append(_row(label, kind, nb, sols[0][1], [y for y, _, _ in sols], ex))
        if kind == "1EC" and nb == 9:
            for f, (_, _, xi) in zip(FREQS, sols):
                ratios.append({"f_Hz": float(f), "xi_1": _pair(xi["psi_1"]),
                               "xi_E_over_xi_1": _pair(xi["psi_E"] / xi["psi_1"]),
                               "xi_C_over_xi_1": _pair(xi["psi_C"] / xi["psi_1"])})
    surface, r_E, r_C = frozen_ratio_surface(D, rule)
    for nb in (1, 3, 9):
        sols = [Y_galerkin(2j * PI * f, D, "1", nb, rule, surface=surface) for f in FREQS]
        rows.append(_row("frozen-ratio single DOF (ratios from 1 GHz)", "frozen", nb,
                         sols[0][1], [y for y, _, _ in sols], ex))
    return ex, rows, {"r_E_at_1GHz": _pair(r_E), "r_C_at_1GHz": _pair(r_C),
                      "free_3dof_bulk_n_max_9": ratios}


def summary() -> dict:
    """The numbers this script stands for.  Read by emit_results.py."""
    slab = self_test(L, SIGMA, MU)
    rule = graded_rule()
    y10k = Y_exact(2j * PI * 1e4, 3, L, SIGMA)
    out = {
        "case": "box_edge_corner_dofs",
        "body": f"cube and square, L = {L} m, sigma = {SIGMA:.3g} S/m, copper",
        "metric": ("|Y_galerkin - Y_exact| / |Y_exact| over 1 Hz .. 1 GHz (46 points); "
                   "Y_exact = heat-content integral (_references/box_heat_content.py)"),
        "reference": {
            "slab_self_test_max_abs_diff": float(slab),
            "cube_abs_Y_at_10kHz": float(abs(y10k)),
            "aitken_foster_N799_at_10kHz": 3.431919,
            "mixed_rank20_closed_Kss_at_10kHz": 3.443338,
        },
        "bump_rate_alpha": ALPHA,
        "quadrature_nodes_per_axis": int(len(rule[0])),
    }
    for D, name in ((3, "cube"), (2, "square")):
        _, rows, ratios = sweep(D, rule)
        out[name] = rows
        out[name + "_amplitude_ratios"] = ratios
    return out


def main():
    r = summary()
    print("Reference: slab self-test %.1e; cube |Y|(10 kHz) exact %.6f, Aitken Foster %.6f, mixed rank-20 %.6f"
          % (r["reference"]["slab_self_test_max_abs_diff"], r["reference"]["cube_abs_Y_at_10kHz"],
             r["reference"]["aitken_foster_N799_at_10kHz"], r["reference"]["mixed_rank20_closed_Kss_at_10kHz"]))
    for name in ("cube", "square"):
        print(f"\n{name}: max relative error over 1 Hz .. 1 GHz (%)")
        print(f"{'surface DOFs':42s} {'bulk':>5s} {'n':>4s} {'max':>8s} {'at f':>7s} " + " ".join(f"{k:>7s}" for k in r[name][0]["error_pct_at"]))
        for row in r[name]:
            print(f"{row['surface'][:42]:42s} n<={row['bulk_n_max']:<2d} {row['n_unknowns']:4d} {row['max_error_pct']:8.4f} {row['at_f_Hz']:7.0e} "
                  + " ".join(f"{v:7.4f}" for v in row["error_pct_at"].values()))
        rat = r[name + "_amplitude_ratios"]
        print(f"  free 3-DOF amplitude ratios (bulk n<=9); at 1 GHz r_E = {rat['r_E_at_1GHz'][0]:.4f}{rat['r_E_at_1GHz'][1]:+.4f}j, "
              f"r_C = {rat['r_C_at_1GHz'][0]:.4f}{rat['r_C_at_1GHz'][1]:+.4f}j")
        for q in rat["free_3dof_bulk_n_max_9"][15::3]:
            print(f"    f = {q['f_Hz']:7.0e} Hz: xi_E/xi_1 = {q['xi_E_over_xi_1'][0]:7.4f}{q['xi_E_over_xi_1'][1]:+8.4f}j   "
                  f"xi_C/xi_1 = {q['xi_C_over_xi_1'][0]:7.4f}{q['xi_C_over_xi_1'][1]:+8.4f}j")


if __name__ == "__main__":
    main()
