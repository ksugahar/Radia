#!/usr/bin/env python
"""Cross-element validation: the SAME physical problem solved with hex, wedge, AND
tet HDiv-VIM meshes must give the SAME volume-averaged magnetization.

This is the concrete demonstration of the manuscript's "任意形状メッシュに共通"
(element-independent) claim: the build techniques (block ledger, adaptive near/far
quadrature, one-sided symmetry) carry across element families, and correctness does
not depend on the element type.  The 1 m cube in a uniform +z field, linear
mu_r = 1000, is meshed three completely different ways --

  * structured pure-HEX lattice (Q2 isoparametric RT1),
  * structured pure-WEDGE / prism lattice (3 sub-tets per prism),
  * unstructured pure-TET (Netgen OCC),

-- and the three volume-averaged M_z values must agree.  They cannot agree by
construction (different DOF counts, different charge supports, different quadrature),
so agreement is real cross-validation of the whole assembly + solve path.

NOTE the demag is NON-analytic here: a high-mu cube magnetizes NON-uniformly, so the
effective demag factor is ~0.27 (M_avg ~ 726 kA/m), NOT the uniform-body 1/3.  The
analytic anchor for this method lives in the SEPARATE sphere check (external field vs
the exact point dipole, docs/hdiv_vim); THIS test locks element-INDEPENDENCE.

Correctness only (numerical agreement), so it runs on LAB per the Benchmark Policy;
timing lives in the mdx scaling study.
"""
import json
import platform
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("ngsolve")

import ngsolve as ng                                    # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh         # noqa: E402

import radia.vim as vim                                 # noqa: E402

HERE = Path(__file__).resolve().parent
CUBE = 1.0
H0 = 200.0e3
MU_R = 1000.0


def _map(x, y, z):
    return (CUBE * (x - 0.5), CUBE * (y - 0.5), CUBE * (z - 0.5))


def _hex(n):
    return MakeStructured3DMesh(hexes=True, nx=n, ny=n, nz=n, mapping=_map)


def _wedge(n):
    return MakeStructured3DMesh(prism=True, nx=n, ny=n, nz=n, mapping=_map)


def _tet(maxh):
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(-0.5 * CUBE, -0.5 * CUBE, -0.5 * CUBE),
                          Pnt(0.5 * CUBE, 0.5 * CUBE, 0.5 * CUBE)))
    return ng.Mesh(geo.GenerateMesh(maxh=maxh))


# Three element families at COMPARABLE resolution (matched to ~40-65k face DOFs).
CASES = [
    ("hex", lambda: _hex(12)),
    ("wedge", lambda: _wedge(12)),
    ("tet", lambda: _tet(0.09)),
]


def _solve_mavg(mesh_builder):
    with ng.TaskManager():
        mesh = mesh_builder()
        r = vim.Solve(mesh, MU_R, ng.CoefficientFunction((0.0, 0.0, H0)), order=1)
    return float(r["M_avg"][2]), int(r["ndof"]), int(mesh.GetNE(ng.VOL))


def test_cross_element_cube_mavg_agrees():
    """hex / wedge / tet on the same cube agree in volume-averaged M_z to < 0.3%."""
    results = {}
    for name, builder in CASES:
        mz, ndof, n_el = _solve_mavg(builder)
        results[name] = dict(M_avg_z=mz, ndof=ndof, n_el=n_el)

    mz = np.array([results[n]["M_avg_z"] for n, _ in CASES])
    mean = float(np.mean(mz))
    spread = float((mz.max() - mz.min()) / mean)

    # Physically sane magnetization (high-mu cube, ~726 kA/m; NOT the uniform 1/3).
    assert 6.0e5 < mean < 8.0e5, "cube M_avg_z out of physical range: %.0f A/m" % mean
    # Element-independence: the three discretizations agree to well under 0.3%.
    assert spread < 3.0e-3, (
        "hex/wedge/tet disagree on the cube M_avg_z by %.2e (hex=%.0f wedge=%.0f tet=%.0f) "
        "-- element-independence broken" % (spread, mz[0], mz[1], mz[2]))

    payload = dict(
        description="HDiv-VIM cross-element cube validation (hex/wedge/tet, same 1 m cube, mu_r=1000, +z 200 kA/m)",
        timestamp=datetime.now().isoformat(),
        hostname=platform.node(),
        H_ext_A_per_m=H0, mu_r=MU_R, cube_size_m=CUBE,
        results={n: results[n] for n, _ in CASES},
        M_avg_z_mean_A_per_m=mean,
        mutual_spread_rel=spread,
    )
    (HERE / "cross_element_cube.json").write_text(json.dumps(payload, indent=1))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v", "-x", "-s"]))
