"""Validation-class d-q cross-saturation flux-map example.

This is an example/validation run rather than a pytest test.  A synthetic but
physics-shaped magnetic coenergy map is used as a stand-in for a nonlinear
machine solve:

* flux linkages are the gradient of coenergy
* the incremental inductance matrix is the coenergy Hessian
* reciprocity requires L_dq = L_qd
* saturation makes the small-signal diagonal inductances roll off with current

Run:

    python examples/electric_machine/validation_cross_saturation_flux_map.py
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.solve import dq_flux_torque, incremental_inductance_matrix  # noqa: E402


OUT_JSON = HERE / "validation_cross_saturation_flux_map_summary.json"
DI = 1.0e-3

PARAMS = {
    "lambda_m": 0.100,
    "Ld_sat": 3.0e-3,
    "Lq_sat": 9.0e-3,
    "Ldrop": 4.0e-3,
    "Iref": 25.0,
    "Id_scale": 18.0,
    "Iq_scale": 22.0,
    "coupling": 0.35,
    "pole_pairs": 4,
}

OPERATING_POINTS = [
    (-1.0, 1.0),
    (-5.0, 5.0),
    (-10.0, 10.0),
    (-15.0, 15.0),
    (-20.0, 20.0),
    (-25.0, 25.0),
    (-30.0, 30.0),
    (-20.0, 5.0),
    (-5.0, 20.0),
]

DIAGONAL_PATH = [(-5.0, 5.0), (-10.0, 10.0), (-15.0, 15.0), (-20.0, 20.0), (-25.0, 25.0), (-30.0, 30.0)]


def _rho(id_: float, iq: float) -> float:
    id_s = PARAMS["Id_scale"]
    iq_s = PARAMS["Iq_scale"]
    c = PARAMS["coupling"]
    return (id_ / id_s) ** 2 + (iq / iq_s) ** 2 + 2.0 * c * id_ * iq / (id_s * iq_s)


def _rho_grad(id_: float, iq: float) -> tuple[float, float]:
    id_s = PARAMS["Id_scale"]
    iq_s = PARAMS["Iq_scale"]
    c = PARAMS["coupling"]
    return (
        2.0 * id_ / id_s ** 2 + 2.0 * c * iq / (id_s * iq_s),
        2.0 * iq / iq_s ** 2 + 2.0 * c * id_ / (id_s * iq_s),
    )


def _rho_hessian() -> list[list[float]]:
    id_s = PARAMS["Id_scale"]
    iq_s = PARAMS["Iq_scale"]
    c = PARAMS["coupling"]
    return [
        [2.0 / id_s ** 2, 2.0 * c / (id_s * iq_s)],
        [2.0 * c / (id_s * iq_s), 2.0 / iq_s ** 2],
    ]


def coenergy(currents: tuple[float, float] | list[float]) -> float:
    id_, iq = currents
    a = 0.5 * PARAMS["Ldrop"] * PARAMS["Iref"] ** 2
    return (
        PARAMS["lambda_m"] * id_
        + 0.5 * PARAMS["Ld_sat"] * id_ ** 2
        + 0.5 * PARAMS["Lq_sat"] * iq ** 2
        + a * math.log1p(_rho(id_, iq))
    )


def flux_linkages(currents: tuple[float, float] | list[float]) -> list[float]:
    id_, iq = currents
    a = 0.5 * PARAMS["Ldrop"] * PARAMS["Iref"] ** 2
    rho = _rho(id_, iq)
    gd, gq = _rho_grad(id_, iq)
    scale = a / (1.0 + rho)
    return [
        PARAMS["lambda_m"] + PARAMS["Ld_sat"] * id_ + scale * gd,
        PARAMS["Lq_sat"] * iq + scale * gq,
    ]


def analytic_incremental_matrix(currents: tuple[float, float] | list[float]) -> list[list[float]]:
    id_, iq = currents
    a = 0.5 * PARAMS["Ldrop"] * PARAMS["Iref"] ** 2
    rho = _rho(id_, iq)
    grad = _rho_grad(id_, iq)
    hess = _rho_hessian()
    mat = [[PARAMS["Ld_sat"], 0.0], [0.0, PARAMS["Lq_sat"]]]
    for r in range(2):
        for c in range(2):
            mat[r][c] += a * (
                hess[r][c] / (1.0 + rho)
                - grad[r] * grad[c] / (1.0 + rho) ** 2
            )
    return mat


def _max_matrix_abs_diff(a: list[list[float]], b: list[list[float]]) -> float:
    return max(abs(a[i][j] - b[i][j]) for i in range(2) for j in range(2))


def _max_matrix_rel_diff(a: list[list[float]], b: list[list[float]]) -> float:
    return max(abs(a[i][j] - b[i][j]) / max(abs(b[i][j]), 1.0e-15) for i in range(2) for j in range(2))


def _record_point(currents: tuple[float, float]) -> dict:
    fd = incremental_inductance_matrix(flux_linkages, list(currents), DI)
    exact = analytic_incremental_matrix(currents)
    lam_d, lam_q = flux_linkages(currents)
    id_, iq = currents
    det = fd[0][0] * fd[1][1] - fd[0][1] * fd[1][0]
    return {
        "id": id_,
        "iq": iq,
        "coenergy": coenergy(currents),
        "lambda_d": lam_d,
        "lambda_q": lam_q,
        "torque": dq_flux_torque(PARAMS["pole_pairs"], lam_d, lam_q, id_, iq),
        "rho": _rho(id_, iq),
        "incremental_matrix_fd": fd,
        "incremental_matrix_exact": exact,
        "symmetry_abs_error": abs(fd[0][1] - fd[1][0]),
        "exact_max_abs_error": _max_matrix_abs_diff(fd, exact),
        "exact_max_rel_error": _max_matrix_rel_diff(fd, exact),
        "determinant": det,
        "cross_ratio": fd[0][1] / math.sqrt(fd[0][0] * fd[1][1]),
    }


def _validate(records: list[dict]) -> dict:
    max_sym = max(rec["symmetry_abs_error"] for rec in records)
    max_abs = max(rec["exact_max_abs_error"] for rec in records)
    max_rel = max(rec["exact_max_rel_error"] for rec in records)
    min_det = min(rec["determinant"] for rec in records)
    by_point = {(rec["id"], rec["iq"]): rec for rec in records}
    path = [by_point[p] for p in DIAGONAL_PATH]
    ldd = [rec["incremental_matrix_fd"][0][0] for rec in path]
    lqq = [rec["incremental_matrix_fd"][1][1] for rec in path]
    ldq = [rec["incremental_matrix_fd"][0][1] for rec in path]
    torque = [rec["torque"] for rec in path]
    checks = {
        "max_symmetry_abs_error": max_sym,
        "max_exact_matrix_abs_error": max_abs,
        "max_exact_matrix_rel_error": max_rel,
        "min_incremental_matrix_determinant": min_det,
        "Ldd_rolloff_ratio_last_over_first": ldd[-1] / ldd[0],
        "Lqq_rolloff_ratio_last_over_first": lqq[-1] / lqq[0],
        "max_cross_inductance": max(abs(x) for x in ldq),
        "last_cross_below_peak": abs(ldq[-1]) < max(abs(x) for x in ldq),
        "torque_monotone_on_diagonal_path": all(b > a for a, b in zip(torque, torque[1:])),
    }
    assert max_sym < 1.0e-9
    assert max_abs < 1.0e-9
    assert max_rel < 1.0e-7
    assert min_det > 0.0
    assert checks["Ldd_rolloff_ratio_last_over_first"] < 0.40
    assert checks["Lqq_rolloff_ratio_last_over_first"] < 0.75
    assert checks["max_cross_inductance"] > 1.0e-3
    assert checks["last_cross_below_peak"]
    assert checks["torque_monotone_on_diagonal_path"]
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    records = [_record_point(p) for p in OPERATING_POINTS]
    checks = _validate(records)
    summary = {
        "kind": "cross_saturation_flux_map_validation",
        "validation_class": True,
        "finite_difference_step": DI,
        "parameters": PARAMS,
        "operating_points": records,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("[cross-saturation d-q flux map]")
    for rec in records:
        mat = rec["incremental_matrix_fd"]
        print(
            f"  id={rec['id']:6.1f} A iq={rec['iq']:6.1f} A  "
            f"Ldd={mat[0][0]*1e3:8.4f} mH  Lqq={mat[1][1]*1e3:8.4f} mH  "
            f"Ldq={mat[0][1]*1e3:8.4f} mH  "
            f"T={rec['torque']:9.5f} Nm  "
            f"sym={rec['symmetry_abs_error']:.3e}"
        )
    print(
        "[checks] "
        f"Ldd rolloff={checks['Ldd_rolloff_ratio_last_over_first']:.6f}, "
        f"Lqq rolloff={checks['Lqq_rolloff_ratio_last_over_first']:.6f}, "
        f"max Hessian rel err={checks['max_exact_matrix_rel_error']:.3e}"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
