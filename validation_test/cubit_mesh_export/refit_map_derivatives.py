"""Higher derivatives of the curved element map: geometric refit vs Netgen.

The convergence theory of curved finite elements (Ciarlet-Raviart, Lenoir)
assumes that the element map F satisfies |D^k F| <= C h^k.  The geometric
refit changes F's edge and face coefficients, including tangential
re-parametrisation, so it could raise these constants while bringing the
boundary closer to the CAD.  This script measures them on exported .vol
files of the same Cubit mesh with Netgen's coefficients and with the refit.

F is a polynomial in the reference coordinates (total degree p on TET,
degree p per variable on HEX).  It is sampled on a reference grid through
NGSolve's vectorised map (the assembly path), recovered exactly by least
squares (the residual is reported), and differentiated analytically.  Per
element and order k = 2..p the quantity is max over the samples of the
Frobenius norm of D^k F divided by h_K^k, h_K the corner diameter.

    python refit_map_derivatives.py --netgen DIR --refit DIR

DIR holds {tet,hex}_o{p}.vol as written by geometric_refit_benchmark.py.
The result JSON is written next to this script.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import itertools
import json
import math
import platform
from pathlib import Path

import numpy as np
from ngsolve import HEX, TET, VOL, CoefficientFunction, IntegrationRule, Mesh, x, y, z

KINDS = {"tet": TET, "hex": HEX}
ORDERS = (2, 3, 4, 5)


def exponents(kind: str, p: int) -> list[tuple[int, int, int]]:
    rng = range(p + 1)
    if kind == "tet":
        return [e for e in itertools.product(rng, rng, rng) if sum(e) <= p]
    return list(itertools.product(rng, rng, rng))


def reference_points(kind: str, p: int) -> np.ndarray:
    m = 2 * p + 3
    s = (np.arange(m) + 0.5) / m
    pts = [(a, b, c) for a in s for b in s for c in s
           if kind == "hex" or a + b + c < 1.0]
    return np.array(pts)


def design(pts: np.ndarray, exps, deriv=(0, 0, 0)) -> np.ndarray:
    """Columns: d^deriv of xi^a eta^b zeta^c at the points."""
    cols = []
    for e in exps:
        col = np.ones(len(pts))
        for axis in range(3):
            n, k = e[axis], deriv[axis]
            if k > n:
                col = np.zeros(len(pts))
                break
            # Monomials in t = 2 xi - 1 keep the fit well conditioned;
            # d/dxi = 2 d/dt.
            coef = math.perm(n, k) * 2.0 ** k
            col = col * coef * (2.0 * pts[:, axis] - 1.0) ** (n - k)
        cols.append(col)
    return np.stack(cols, axis=1)


def element_measures(mesh: Mesh, kind: str, p: int) -> dict:
    exps = exponents(kind, p)
    pts = reference_points(kind, p)
    ir = IntegrationRule(points=[tuple(q) for q in pts], weights=[1.0] * len(pts))
    xyz = np.array(CoefficientFunction((x, y, z))(mesh.MapToAllElements({KINDS[kind]: ir}, VOL)))
    ne = mesh.ne
    xyz = xyz.reshape(ne, len(pts), 3)
    V = design(pts, exps)
    Vpinv = np.linalg.pinv(V)
    derivs = {k: [d for d in itertools.product(range(k + 1), repeat=3) if sum(d) == k]
              for k in range(2, p + 1)}
    dmats = {k: [(d, design(pts, exps, d)) for d in derivs[k]] for k in derivs}
    corners = np.array([[mesh[v].point for v in el.vertices] for el in mesh.Elements(VOL)])
    h = np.max(np.linalg.norm(corners[:, :, None, :] - corners[:, None, :, :], axis=3), axis=(1, 2))
    fit_res = 0.0
    ratio = {k: np.zeros(ne) for k in derivs}
    for i in range(ne):
        C = Vpinv @ xyz[i]                       # monomial coefficients, (nmono, 3)
        fit_res = max(fit_res, float(np.abs(V @ C - xyz[i]).max()) / h[i])
        for k, mats in dmats.items():
            # Frobenius norm of the k-th derivative tensor: sum over ordered
            # index tuples = sum over multi-indices times multinomial count.
            total = np.zeros(len(pts))
            for d, D in mats:
                mult = math.factorial(k) // (math.factorial(d[0]) * math.factorial(d[1]) * math.factorial(d[2]))
                g = D @ C
                total += mult * np.sum(g * g, axis=1)
            ratio[k][i] = math.sqrt(total.max()) / h[i] ** k
    return {"h": h, "ratio": ratio, "fit_residual_rel": fit_res}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--netgen", type=Path, required=True)
    ap.add_argument("--refit", type=Path, required=True)
    args = ap.parse_args()
    rows = []
    for kind in KINDS:
        for p in ORDERS:
            base = element_measures(Mesh(str(args.netgen / f"{kind}_o{p}.vol")), kind, p)
            new = element_measures(Mesh(str(args.refit / f"{kind}_o{p}.vol")), kind, p)
            for k in base["ratio"]:
                b, n = base["ratio"][k], new["ratio"][k]
                curved = b > 1e-3 * b.max()       # elements Netgen visibly curves
                q = n[curved] / b[curved]
                row = dict(kind=kind, order=p, k=k, elements=int(len(b)),
                           curved_elements=int(curved.sum()),
                           netgen_max=float(b.max()), refit_max=float(n.max()),
                           ratio_median=float(np.median(q)), ratio_p95=float(np.percentile(q, 95)),
                           ratio_max=float(q.max()),
                           fit_residual_rel=max(base["fit_residual_rel"], new["fit_residual_rel"]))
                rows.append(row)
                print(f"{kind} p={p} k={k}: max |D^kF|/h^k netgen {row['netgen_max']:.3e} "
                      f"refit {row['refit_max']:.3e}; refit/netgen per element median "
                      f"{row['ratio_median']:.2f} p95 {row['ratio_p95']:.2f} max {row['ratio_max']:.2f} "
                      f"(fit residual {row['fit_residual_rel']:.1e})", flush=True)
    result = {
        "protocol": "same Cubit unit-sphere meshes (TET size 0.65, HEX sphere scheme 0.4) exported "
                    "with Netgen coefficients and with the geometric refit; D^k F of each element "
                    "recovered from the vectorised NGSolve map and normalised by h_K^k",
        "environment": {"host": platform.node(), "python": platform.python_version(),
                        "ngsolve": importlib.metadata.version("ngsolve"),
                        "netgen_dir": str(args.netgen), "refit_dir": str(args.refit)},
        "rows": rows,
    }
    target = Path(__file__).with_name("refit_map_derivatives_results.json")
    target.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {target}")


if __name__ == "__main__":
    main()
