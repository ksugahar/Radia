"""HDiv p-convergence on a Cubit curved-hex sphere.

This example uses the committed Cubit-generated sphere meshes from
``docs/cubit_mesh_export/hex_sphere_highorder``:

    volume 1 scheme sphere
    block 1 add hex all
    export netgen "hexsph_oN.vol" order N overwrite

The benchmark is deliberately analytic.  For a linear permeable sphere in an
applied solid-harmonic field of degree ``l``, the field inside the sphere is

    H_in = ((2*l + 1) / (l*mu_r + l + 1)) H_ext.

We project that exact internal magnetization into NGSolve ``HDiv(order=p)`` on
the Cubit hex mesh, reconstruct the external magnetic field from the polynomial
volume/surface charges, and compare against the exact scattered multipole field.

Two readings are reported:

* coupled mesh/order: load Cubit's order-1/2/3 curved hex .vol and use the same
  HDiv order.  This is the public "curved Cubit hex + HDiv p" evidence.
* fixed geometry: keep Cubit's order-3 curved hex mesh and vary only HDiv order.
  This separates basis p-convergence from geometry curving.

This is a projection + field-reconstruction benchmark, not a self-consistent
high-order hex VIM solve.  The high-order DemagOperator charge Gram is currently
tet/tri based; using it blindly on hex/quad high-order geometry would be the
wrong proof.  This example tests the Cubit curved-hex geometry, NGSolve HDiv
basis, and element-type-agnostic polynomial charge field path directly.

Run:
    python hdiv_cubit_hex_sphere_pconv.py
"""

from __future__ import annotations

import argparse
import json
import math
import os
from dataclasses import dataclass
from typing import Callable

import numpy as np
import ngsolve as ng
from ngsolve import CoefficientFunction, Integrate, Mesh, TaskManager, VOL

from radia.vim._field import _project_pointdata_to_hdiv, reconstruct_field_polynomial


HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DEFAULT_MESH_DIR = os.path.join(REPO, "examples", "cubit_mesh_export", "hex_sphere_highorder")
OUT_JSON = os.path.join(HERE, "hdiv_cubit_hex_sphere_pconv.json")

RADIUS = 0.05
MU_R = 10.0
CHI = MU_R - 1.0
H0 = 1.0e3
G_QUAD = 2.0e5


@dataclass(frozen=True)
class MultipoleCase:
    name: str
    degree_l: int
    amplitude: float
    exact_magnetization: Callable[[np.ndarray], np.ndarray]
    exact_scattered_h: Callable[[np.ndarray], np.ndarray]
    observation_points: np.ndarray


def _dipole_magnetization(points: np.ndarray) -> np.ndarray:
    coeff = 3.0 / (MU_R + 2.0)
    m = np.array([0.0, 0.0, CHI * coeff * H0])
    return np.tile(m, (len(points), 1))


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


def _quadrupole_magnetization(points: np.ndarray) -> np.ndarray:
    coeff = 5.0 / (2.0 * MU_R + 3.0)
    x, y = points[:, 0], points[:, 1]
    return CHI * coeff * G_QUAD * np.stack([x, -y, np.zeros_like(x)], axis=1)


def _quadrupole_scattered_h(points: np.ndarray) -> np.ndarray:
    coeff = G_QUAD * RADIUS**5 * (2.0 * (MU_R - 1.0)) / (2.0 * MU_R + 3.0)
    x, y = points[:, 0], points[:, 1]
    r = np.linalg.norm(points, axis=1)
    solid = 0.5 * (x * x - y * y)
    grad_solid = np.stack([x, -y, np.zeros_like(x)], axis=1)
    grad_phi = coeff * (grad_solid / r[:, None] ** 5 - 5.0 * solid[:, None] * points / r[:, None] ** 7)
    return -grad_phi


def _cases() -> list[MultipoleCase]:
    return [
        MultipoleCase(
            name="dipole_z",
            degree_l=1,
            amplitude=H0,
            exact_magnetization=_dipole_magnetization,
            exact_scattered_h=_dipole_scattered_h,
            observation_points=RADIUS
            * np.array(
                [
                    [0.0, 0.0, 2.0],
                    [2.0, 0.0, 0.0],
                    [1.7, 0.8, 0.9],
                    [-1.6, 0.7, 1.1],
                ],
                float,
            ),
        ),
        MultipoleCase(
            name="quadrupole_x2_y2",
            degree_l=2,
            amplitude=G_QUAD,
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
    got = float(Integrate(CoefficientFunction(1.0), mesh, VOL))
    return got / exact - 1.0


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


def _field_relative_error(mesh: Mesh, gf, case: MultipoleCase, order: int) -> float:
    quad = max(4, order + 2)
    got = reconstruct_field_polynomial(mesh, gf, case.observation_points, quad=quad, quantity="h")
    exact = case.exact_scattered_h(case.observation_points)
    return float(np.linalg.norm(got - exact) / np.linalg.norm(exact))


def _run_one(mesh: Mesh, hdiv_order: int, cases: list[MultipoleCase]) -> dict:
    out = {}
    for case in cases:
        with TaskManager():
            gf = _project_pointdata_to_hdiv(
                mesh,
                case.exact_magnetization,
                hdiv_order,
                nq=max(4, 2 * hdiv_order + 4),
            )
            l2_rel = _relative_l2_error(mesh, gf, case.exact_magnetization, hdiv_order)
            field_rel = _field_relative_error(mesh, gf, case, hdiv_order)
        out[case.name] = {
            "degree_l": case.degree_l,
            "amplitude": case.amplitude,
            "hdiv_ndof": int(gf.space.ndof),
            "magnetization_l2_rel": l2_rel,
            "external_field_rel": field_rel,
        }
    return out


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
                "cases": _run_one(mesh, order, cases),
            }
        )

    fixed_mesh_path = os.path.join(mesh_dir, "hexsph_o3.vol")
    fixed_mesh = Mesh(fixed_mesh_path)
    fixed = []
    for hdiv_order in (0, 1, 2, 3):
        fixed.append(
            {
                "mesh_order": 3,
                "hdiv_order": hdiv_order,
                "mesh_path": os.path.relpath(fixed_mesh_path, REPO).replace("\\", "/"),
                "n_hex": int(fixed_mesh.GetNE(ng.VOL)),
                "volume_error_rel": _mesh_volume_error(fixed_mesh),
                "cases": _run_one(fixed_mesh, hdiv_order, cases),
            }
        )

    return {
        "description": "HDiv p-convergence on Cubit-generated curved hex sphere; analytic dipole and quadrupole sphere response.",
        "radius_m": RADIUS,
        "mu_r": MU_R,
        "mesh_source": os.path.relpath(mesh_dir, REPO).replace("\\", "/"),
        "cubit_recipe": [
            "create sphere radius 0.05",
            "volume 1 size 0.012",
            "volume 1 scheme sphere",
            "mesh volume 1",
            "block 1 add hex all",
            "export netgen hexsph_oN.vol order N overwrite",
        ],
        "analytics": {
            "response_factor": "(2*l+1)/(l*mu_r+l+1)",
            "dipole_z": "H_ext=(0,0,H0), phi_s=H0*a^3*(mu_r-1)/(mu_r+2)*z/r^3",
            "quadrupole_x2_y2": "H_ext=G*(x,-y,0), phi_s=G*a^5*2*(mu_r-1)/(2*mu_r+3)*(x^2-y^2)/(2*r^5)",
        },
        "coupled_geometry_and_hdiv_order": coupled,
        "fixed_order3_geometry_hdiv_p": fixed,
    }


def _print_summary(data: dict) -> None:
    print("Cubit curved-hex sphere HDiv p-convergence")
    print(f"  mesh source: {data['mesh_source']}")
    print(f"  radius={data['radius_m']} m  mu_r={data['mu_r']}")
    print("\nCoupled Cubit mesh order and HDiv order:")
    for rec in data["coupled_geometry_and_hdiv_order"]:
        print(f"  order {rec['hdiv_order']}  hex={rec['n_hex']}  volume_err={100*rec['volume_error_rel']:+.3f}%")
        for name, case in rec["cases"].items():
            print(
                f"    {name:18s} ndof={case['hdiv_ndof']:5d}  "
                f"M_L2={case['magnetization_l2_rel']:.3e}  "
                f"Hext_rel={case['external_field_rel']:.3e}"
            )
    print("\nFixed Cubit order-3 curved geometry, vary HDiv order only:")
    for rec in data["fixed_order3_geometry_hdiv_p"]:
        print(f"  hdiv p={rec['hdiv_order']}  hex={rec['n_hex']}  volume_err={100*rec['volume_error_rel']:+.3f}%")
        for name, case in rec["cases"].items():
            print(
                f"    {name:18s} ndof={case['hdiv_ndof']:5d}  "
                f"M_L2={case['magnetization_l2_rel']:.3e}  "
                f"Hext_rel={case['external_field_rel']:.3e}"
            )


def _validate(data: dict) -> None:
    coupled = data["coupled_geometry_and_hdiv_order"]
    for name in ("dipole_z", "quadrupole_x2_y2"):
        e1 = coupled[0]["cases"][name]["external_field_rel"]
        e2 = coupled[1]["cases"][name]["external_field_rel"]
        e3 = coupled[2]["cases"][name]["external_field_rel"]
        assert e2 < 0.05 * e1, f"{name}: curved order-2 field error did not collapse vs order-1"
        assert e3 < e2, f"{name}: coupled order-3 field error must improve over order-2"

    fixed = data["fixed_order3_geometry_hdiv_p"]
    for name in ("dipole_z", "quadrupole_x2_y2"):
        errs = [rec["cases"][name]["external_field_rel"] for rec in fixed]
        assert errs[1] < errs[0], f"{name}: HDiv p=1 must improve over p=0 on fixed curved geometry"
        assert errs[2] < errs[1], f"{name}: HDiv p=2 must improve over p=1 on fixed curved geometry"
        if name == "dipole_z":
            # The dipole case has constant internal M; after p=2 it is dominated by the fixed curved
            # geometry / quadrature floor.  The quadrupole carries the strict basis-order staircase.
            assert errs[2] < 2.0e-3 and errs[3] < 2.0e-3, \
                f"{name}: HDiv p>=2 should stay at the fixed-geometry field-error floor"
        else:
            assert errs[3] < errs[2], f"{name}: HDiv p=3 must improve over p=2 on fixed curved geometry"


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
