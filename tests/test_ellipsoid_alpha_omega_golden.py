"""Golden test: full-frequency axial eddy polarizability alpha_c(omega) of a
conducting spheroid (axisymmetric FEM), validated against the analytic sphere
and the analytic demag-tensor HF anchor.

Locks docs/maglev/demos/ellipsoid/ellipsoid_alpha_omega_axisym.py:
  - the FEM alpha_c(omega) of a SPHERE reproduces the analytic 4 pi a^3 G(x)
    (the uniform-field BC + moment extractor are correct),
  - the high-frequency limit of a spheroid reaches the analytic perfect-
    conductor value -V/(1-N_c) from the demag tensor.

Runs a few real axisymmetric NGSolve solves (~30-60 s); skipped cleanly if
ngsolve / netgen are unavailable.
"""
import os
import sys

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

_HERE = os.path.dirname(os.path.abspath(__file__))
# CLN was absorbed into radia.maglev ->
# docs/maglev/demos/{sphere,ellipsoid}. Add both so sphere + ellipsoid modules resolve.
_LEV = os.path.join(_HERE, "..", "docs", "maglev", "demos")
sys.path.insert(0, os.path.join(_LEV, "sphere"))
sys.path.insert(0, os.path.join(_LEV, "ellipsoid"))

import ellipsoid_alpha_omega_axisym as M  # noqa: E402
import ellipsoid_alpha_tensor as ET  # noqa: E402


@pytest.fixture(scope="module")
def sphere_mesh():
    return M.build_mesh(M.A_SPH, M.A_SPH, box=0.10, maxh=0.012)


def test_fem_sphere_reproduces_analytic_G(sphere_mesh):
    # FEM alpha_c(omega) must match the analytic Landau sphere 4 pi a^3 G(x)
    for f in (1e4, 5e4):
        w = 2 * np.pi * f
        afem = M.alpha_c(sphere_mesh, w)
        aan = M.sphere_analytic(w)
        rel = abs(afem - aan) / abs(aan)
        assert rel < 0.05, f"f={f}: FEM {afem:.3e} vs analytic {aan:.3e} ({rel*100:.1f}%)"


def test_fem_sphere_real_part_is_lift_and_negative(sphere_mesh):
    # Re[alpha_c] < 0 (diamagnetic / flux-excluding) -> lift
    a = M.alpha_c(sphere_mesh, 2 * np.pi * 5e4)
    assert a.real < 0.0


def test_hf_anchor_oblate_matches_demag_tensor():
    # an oblate spheroid (cheap: c < a) -> alpha_c(HF) ~ -V/(1-N_c)
    ae, c = M.A_SPH, 0.5 * M.A_SPH
    N = ET.demag_tensor((ae, ae, c))
    V = 4 / 3 * np.pi * ae * ae * c
    kappa_hf = -V / (1.0 - N[2])
    m = M.build_mesh(ae, c, box=0.10, maxh=0.012)
    afem = M.alpha_c(m, 2 * np.pi * 5e5)
    rel = abs(afem.real - kappa_hf) / abs(kappa_hf)
    assert rel < 0.08, (
        f"oblate HF alpha_c {afem.real:.3e} vs -V/(1-N_c) {kappa_hf:.3e} ({rel*100:.1f}%)")


def test_shape_changes_per_volume_response():
    # the per-volume HF response 1/(1-N_c) is larger for the oblate than the
    # sphere -> shape (not material) anisotropy
    ae = M.A_SPH
    N_sph = ET.demag_tensor((ae, ae, ae))[2]
    N_obl = ET.demag_tensor((ae, ae, 0.5 * ae))[2]
    assert 1 / (1 - N_obl) > 1 / (1 - N_sph)
