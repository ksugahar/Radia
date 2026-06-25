"""Harmonic-subspace HDiv-VIM solve on a Cubit curved-hex sphere.

This is the self-consistent companion to ``hdiv_cubit_hex_sphere_pconv.py``.
It uses the same Cubit-generated curved hex sphere meshes, but now computes the
demagnetizing response coefficient from the HDiv-VIM surface single-layer energy:

    D_h = <sigma, V sigma> / <M, M>
    M_h = chi / (1 + chi D_h) * H_ext,h.

For a permeable sphere in a solid-harmonic applied field of degree ``l``, the
analytic demag coefficient is ``D_l = l/(2*l+1)``:

    l=1 dipole/uniform field:       D_1 = 1/3
    l=2 quadrupole field:           D_2 = 2/5

Why a harmonic subspace instead of the full HDiv space?  High-order HDiv contains
charge-free solenoidal modes.  Those modes are essential for de Rham robustness,
but a pure energy solve over the whole HDiv space can excite non-physical modes
unless the curl-free material relation is enforced as well.  For a sphere in a
single harmonic field, the physical solution is known to stay in the same
solid-harmonic subspace, so this scalar Rayleigh solve is the correct
self-consistent benchmark.

The important geometry point remains the same: Cubit's high-order ``.vol`` carries
curved hex nodes that multipole-moment MMM's flat per-cell surface formulation does not use.

Run:
    python hdiv_cubit_hex_sphere_harmonic_vim.py
"""

from __future__ import annotations

import argparse
import json
import math
import os
from dataclasses import dataclass
from typing import Callable

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import ngsolve as ng
from ngsolve import CoefficientFunction, Integrate, Mesh, TaskManager, VOL
from ngsolve.bem import SingleLayerPotentialOperator

from radia.vim._field import _project_pointdata_to_hdiv, reconstruct_field_polynomial


HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DEFAULT_MESH_DIR = os.path.join(REPO, "examples", "cubit_mesh_export", "hex_sphere_highorder")
OUT_JSON = os.path.join(HERE, "hdiv_cubit_hex_sphere_harmonic_vim.json")

RADIUS = 0.05
MU_R = 10.0
CHI = MU_R - 1.0
H0 = 1.0e3
G_QUAD = 2.0e5


@dataclass(frozen=True)
class HarmonicCase:
    name: str
    degree_l: int
    amplitude: float
    demag_exact: float
    h_mode: Callable[[np.ndarray], np.ndarray]
    exact_magnetization: Callable[[np.ndarray], np.ndarray]
    exact_scattered_h: Callable[[np.ndarray], np.ndarray]
    observation_points: np.ndarray


def _csr(mat) -> sp.csr_matrix:
    rows, cols, vals = mat.COO()
    return sp.csr_matrix((np.array(vals), (np.array(rows), np.array(cols))), shape=(mat.height, mat.width))


def _dense_base_matrix(mat, space) -> np.ndarray:
    gf = ng.GridFunction(space)
    x = gf.vec.CreateVector()
    y = gf.vec.CreateVector()
    out = np.zeros((space.ndof, space.ndof))
    for j in range(space.ndof):
        x[:] = 0.0
        x[j] = 1.0
        y.data = mat * x
        out[:, j] = y.FV().NumPy()
    return out


def _dipole_h(points: np.ndarray) -> np.ndarray:
    return np.tile([0.0, 0.0, H0], (len(points), 1))


def _dipole_magnetization(points: np.ndarray) -> np.ndarray:
    factor = CHI / (1.0 + CHI / 3.0)
    return factor * _dipole_h(points)


def _dipole_scattered_h(points: np.ndarray) -> np.ndarray:
    coeff = H0 * RADIUS**3 * (MU_R - 1.0) / (MU_R + 2.0)
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    r2 = x * x + y * y + z * z
    r = np.sqrt(r2)
    return np.stack(
        [
            3.0 * coeff * x * z / r**5,
            3.0 * coeff * y * z / r**5,
            -coeff * (r2 - 3.0 * z * z) / r**5,
        ],
        axis=1,
    )


def _quadrupole_h(points: np.ndarray) -> np.ndarray:
    x, y = points[:, 0], points[:, 1]
    return G_QUAD * np.stack([x, -y, np.zeros_like(x)], axis=1)


def _quadrupole_magnetization(points: np.ndarray) -> np.ndarray:
    factor = CHI / (1.0 + CHI * 2.0 / 5.0)
    return factor * _quadrupole_h(points)


def _quadrupole_scattered_h(points: np.ndarray) -> np.ndarray:
    coeff = G_QUAD * RADIUS**5 * (2.0 * (MU_R - 1.0)) / (2.0 * MU_R + 3.0)
    x, y = points[:, 0], points[:, 1]
    r = np.linalg.norm(points, axis=1)
    solid = 0.5 * (x * x - y * y)
    grad_solid = np.stack([x, -y, np.zeros_like(x)], axis=1)
    return -coeff * (grad_solid / r[:, None] ** 5 - 5.0 * solid[:, None] * points / r[:, None] ** 7)


def _cases() -> list[HarmonicCase]:
    return [
        HarmonicCase(
            name="dipole_z",
            degree_l=1,
            amplitude=H0,
            demag_exact=1.0 / 3.0,
            h_mode=_dipole_h,
            exact_magnetization=_dipole_magnetization,
            exact_scattered_h=_dipole_scattered_h,
            observation_points=RADIUS
            * np.array(
                [[0.0, 0.0, 2.0], [2.0, 0.0, 0.0], [1.7, 0.8, 0.9], [-1.6, 0.7, 1.1]],
                float,
            ),
        ),
        HarmonicCase(
            name="quadrupole_x2_y2",
            degree_l=2,
            amplitude=G_QUAD,
            demag_exact=2.0 / 5.0,
            h_mode=_quadrupole_h,
            exact_magnetization=_quadrupole_magnetization,
            exact_scattered_h=_quadrupole_scattered_h,
            observation_points=RADIUS
            * np.array(
                [
                    [2.0, 0.8, 0.6],
                    [-1.7, 1.1, 0.5],
                    [1.4, -1.6, 0.7],
                    [0.9, 1.2, 1.8],
                    [-1.3, -0.9, 1.6],
                ],
                float,
            ),
        ),
    ]


def _mesh_volume_error(mesh: Mesh) -> float:
    exact = 4.0 / 3.0 * math.pi * RADIUS**3
    return float(Integrate(CoefficientFunction(1.0), mesh, VOL)) / exact - 1.0


def _surface_map_and_gram(mesh: Mesh, fes, order: int):
    u = fes.TrialFunction()
    sl2 = ng.SurfaceL2(mesh, order=order)
    sigma_trial, sigma_test = sl2.TnT()
    normal = ng.specialcf.normal(3)

    mass_b = ng.BilinearForm(sl2)
    mass_b += sigma_trial * sigma_test * ng.ds
    mass_b.Assemble()

    trace_b = ng.BilinearForm(trialspace=fes, testspace=sl2)
    trace_b += (u.Trace() * normal) * sigma_test * ng.ds
    trace_b.Assemble()

    single_layer = SingleLayerPotentialOperator(sl2, intorder=max(8, 2 * order + 6))
    V_dense = _dense_base_matrix(single_layer.mat, sl2)
    C = spla.spsolve(sp.csc_matrix(_csr(mass_b.mat)), sp.csc_matrix(_csr(trace_b.mat))).tocsr()
    return C, V_dense


def _mass_matrix(fes) -> sp.csr_matrix:
    u, v = fes.TnT()
    mass = ng.BilinearForm(fes)
    mass += u * v * ng.dx
    mass.Assemble()
    return _csr(mass.mat)


def _project_mode(mesh: Mesh, order: int, mode_fn: Callable[[np.ndarray], np.ndarray]):
    with TaskManager():
        gf = _project_pointdata_to_hdiv(mesh, mode_fn, order, nq=max(4, 2 * order + 4))
    return gf


def _relative_l2_error(mesh: Mesh, gf, exact_fn: Callable[[np.ndarray], np.ndarray], order: int) -> float:
    num = 0.0
    den = 0.0
    intorder = max(8, 2 * order + 8)
    for i in range(mesh.GetNE(ng.VOL)):
        ei = ng.ElementId(ng.VOL, i)
        trafo = mesh.GetTrafo(ei)
        for ip in ng.IntegrationRule(mesh[ei].type, intorder):
            mip = trafo(ip)
            xq = np.array([[mip.point[0], mip.point[1], mip.point[2]]])
            exact = exact_fn(xq)[0]
            got = np.array([float(gf[k](mip)) for k in range(3)])
            weight = ip.weight * mip.measure
            num += weight * float(np.dot(got - exact, got - exact))
            den += weight * float(np.dot(exact, exact))
    return math.sqrt(num / den)


def _field_relative_error(mesh: Mesh, gf, case: HarmonicCase, order: int) -> float:
    got = reconstruct_field_polynomial(
        mesh,
        gf,
        case.observation_points,
        quad=max(4, order + 2),
        quantity="h",
    )
    exact = case.exact_scattered_h(case.observation_points)
    return float(np.linalg.norm(got - exact) / np.linalg.norm(exact))


def _solve_case(mesh: Mesh, order: int, case: HarmonicCase) -> dict:
    fes = ng.HDiv(mesh, order=order)
    with TaskManager():
        mass = _mass_matrix(fes)
        C, V_dense = _surface_map_and_gram(mesh, fes, order)
    gf_h = _project_mode(mesh, order, case.h_mode)
    h_vec = np.asarray(gf_h.vec)
    sigma = C @ h_vec
    denom = float(h_vec @ (mass @ h_vec))
    demag = float(sigma @ (V_dense @ sigma) / denom)
    response = CHI / (1.0 + CHI * demag)

    gf_m = ng.GridFunction(fes)
    gf_m.vec.FV().NumPy()[:] = response * h_vec
    m_l2_rel = _relative_l2_error(mesh, gf_m, case.exact_magnetization, order)
    field_rel = _field_relative_error(mesh, gf_m, case, order)
    return {
        "degree_l": case.degree_l,
        "hdiv_ndof": int(fes.ndof),
        "surface_l2_ndof": int(V_dense.shape[0]),
        "demag": demag,
        "demag_exact": case.demag_exact,
        "demag_rel_err": demag / case.demag_exact - 1.0,
        "response": response,
        "response_exact": CHI / (1.0 + CHI * case.demag_exact),
        "response_rel_err": response / (CHI / (1.0 + CHI * case.demag_exact)) - 1.0,
        "magnetization_l2_rel": m_l2_rel,
        "external_field_rel": field_rel,
    }


def run(mesh_dir: str = DEFAULT_MESH_DIR) -> dict:
    cases = _cases()
    coupled = []
    for order in (1, 2, 3):
        mesh_path = os.path.join(mesh_dir, f"hexsph_o{order}.vol")
        mesh = Mesh(mesh_path)
        coupled.append(
            {
                "mesh_order": order,
                "hdiv_order": order,
                "mesh_path": os.path.relpath(mesh_path, REPO).replace("\\", "/"),
                "n_hex": int(mesh.GetNE(ng.VOL)),
                "volume_error_rel": _mesh_volume_error(mesh),
                "cases": {case.name: _solve_case(mesh, order, case) for case in cases},
            }
        )
    return {
        "description": "Harmonic-subspace self-consistent HDiv-VIM solve on Cubit curved hex sphere.",
        "radius_m": RADIUS,
        "mu_r": MU_R,
        "chi": CHI,
        "mesh_source": os.path.relpath(mesh_dir, REPO).replace("\\", "/"),
        "analytics": {
            "sphere_demag_l": "D_l = l/(2*l+1)",
            "response": "M_l = chi/(1+chi*D_l) * H_ext,l",
        },
        "coupled_geometry_and_hdiv_order": coupled,
    }


def _validate(data: dict) -> None:
    recs = data["coupled_geometry_and_hdiv_order"]
    for name in ("dipole_z", "quadrupole_x2_y2"):
        field_errs = [rec["cases"][name]["external_field_rel"] for rec in recs]
        assert field_errs[1] < 0.02 * field_errs[0], f"{name}: order-2 field must collapse from order-1"
        assert field_errs[2] < field_errs[1], f"{name}: order-3 field must improve over order-2"
        assert field_errs[2] < 2e-3, f"{name}: order-3 external field should be below 0.2%"

    dipole_demag_errs = [abs(rec["cases"]["dipole_z"]["demag_rel_err"]) for rec in recs]
    assert max(dipole_demag_errs) < 5e-4, "dipole demag should be accurate for every curved-hex order"

    quad_demag_errs = [abs(rec["cases"]["quadrupole_x2_y2"]["demag_rel_err"]) for rec in recs]
    assert quad_demag_errs[1] < quad_demag_errs[0], "quadrupole order-2 demag must improve over order-1"
    assert quad_demag_errs[2] < quad_demag_errs[1], "quadrupole order-3 demag must improve over order-2"
    assert quad_demag_errs[2] < 3e-3, "quadrupole order-3 demag should be below 0.3%"


def _print_summary(data: dict) -> None:
    print("Cubit curved-hex sphere harmonic HDiv-VIM solve")
    print(f"  mesh source: {data['mesh_source']}")
    print(f"  radius={data['radius_m']} m  mu_r={data['mu_r']}")
    for rec in data["coupled_geometry_and_hdiv_order"]:
        print(f"\n  order {rec['hdiv_order']}  hex={rec['n_hex']}  volume_err={100*rec['volume_error_rel']:+.3f}%")
        for name, case in rec["cases"].items():
            print(
                f"    {name:18s} ndof={case['hdiv_ndof']:5d}  "
                f"D={case['demag']:.6f} ({100*case['demag_rel_err']:+.3f}%)  "
                f"M_L2={case['magnetization_l2_rel']:.3e}  "
                f"Hext_rel={case['external_field_rel']:.3e}"
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mesh-dir", default=DEFAULT_MESH_DIR, help="directory containing hexsph_o{1,2,3}.vol")
    parser.add_argument("--out", default=OUT_JSON, help="output JSON path")
    args = parser.parse_args()

    ng.SetNumThreads(4)
    data = run(args.mesh_dir)
    _validate(data)
    _print_summary(data)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
        f.write("\n")
    print(f"\nsaved {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
