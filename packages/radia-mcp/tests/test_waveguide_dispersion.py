r"""Waveguide dispersion -- propagation above cutoff + evanescence below (#65) -- test.

beta=(2pi/c)sqrt(f^2-fc^2); lambda_g=lambda0/sqrt(1-(fc/f)^2); v_p=c/sqrt(1-(fc/f)^2);
v_g=c sqrt(1-(fc/f)^2); v_p*v_g=c^2; below cutoff alpha=(2pi/c)sqrt(fc^2-f^2). Closed-form
identities + a numerical d omega/d beta check + an FE-anchored cutoff (the propagation sequel to
the #53/#62 cutoff eigenvalues)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.waveguide import (rectangular_waveguide_cutoff, cutoff_frequency,
                                               waveguide_dispersion, guide_wavelength,
                                               waveguide_evanescent_attenuation, C0)


def test_dispersion_identities():
    fc = rectangular_waveguide_cutoff(0.02286, 0.01016, 1, 0)     # WR-90 TE10
    for f in (8e9, 10e9, 12e9, 15e9):
        d = waveguide_dispersion(f, fc)
        # the reciprocal velocity identity v_p * v_g = c^2 (exact)
        assert math.isclose(d["v_phase"] * d["v_group"], C0 ** 2, rel_tol=1e-12)
        # phase super-luminal, group sub-luminal
        assert d["v_phase"] > C0 > d["v_group"] > 0.0
        # guide-wavelength relation 1/lg^2 = 1/l0^2 - 1/lc^2 (exact) and lg > l0
        lam0, lam_c = C0 / f, C0 / fc
        assert math.isclose(1.0 / d["lambda_g"] ** 2, 1.0 / lam0 ** 2 - 1.0 / lam_c ** 2, rel_tol=1e-12)
        assert d["lambda_g"] > lam0
        assert math.isclose(guide_wavelength(f, fc), d["lambda_g"], rel_tol=1e-12)
        # group velocity == d omega/d beta from omega(beta) = c sqrt(beta^2 + kc^2) (central diff)
        kc = 2.0 * math.pi * fc / C0
        db = d["beta"] * 1e-6
        w = lambda bb: C0 * math.sqrt(bb * bb + kc * kc)
        vg_num = (w(d["beta"] + db) - w(d["beta"] - db)) / (2.0 * db)
        assert math.isclose(vg_num, d["v_group"], rel_tol=1e-6)


def test_cutoff_and_high_frequency_limits():
    fc = 6.5571e9
    # propagation branch raises at/below cutoff; evanescent branch raises at/above it
    with pytest.raises(ValueError):
        waveguide_dispersion(fc, fc)
    with pytest.raises(ValueError):
        waveguide_dispersion(0.9 * fc, fc)
    with pytest.raises(ValueError):
        waveguide_evanescent_attenuation(fc, fc)
    # just above cutoff: v_g -> 0, v_p -> inf, beta -> 0
    d = waveguide_dispersion(1.0001 * fc, fc)
    assert d["v_group"] < 0.02 * C0 and d["v_phase"] > 50.0 * C0
    assert waveguide_dispersion(1.000001 * fc, fc)["beta"] < 1.0
    # high frequency: both velocities -> c
    d = waveguide_dispersion(100.0 * fc, fc)
    assert math.isclose(d["v_phase"], C0, rel_tol=1e-3)
    assert math.isclose(d["v_group"], C0, rel_tol=1e-3)


def test_evanescent_limits():
    fc = 6.5571e9
    kc = 2.0 * math.pi * fc / C0
    # alpha -> k_c as f -> 0, and alpha -> 0 as f -> fc
    assert math.isclose(waveguide_evanescent_attenuation(1e-3 * fc, fc), kc, rel_tol=1e-3)
    assert waveguide_evanescent_attenuation(0.9999 * fc, fc) < 0.02 * kc
    assert waveguide_evanescent_attenuation(0.999999 * fc, fc) < 1.0     # continuity at cutoff


@pytest.mark.xval
def test_dispersion_fe_anchor():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, TaskManager
    from netgen.occ import OCCGeometry, WorkPlane
    from radia_mcp.radia_ngsolve.waveguide import helmholtz_cutoff_wavenumbers_2d

    a, b = 0.02286, 0.01016
    rect = WorkPlane().Rectangle(a, b).Face(); rect.edges.name = "wall"
    mesh = Mesh(OCCGeometry(rect, dim=2).GenerateMesh(maxh=min(a, b) / 12))
    with TaskManager():
        kc_fe = helmholtz_cutoff_wavenumbers_2d(mesh, 1, bc="neumann")[0]

    fc_fe = cutoff_frequency(kc_fe)
    fc_an = rectangular_waveguide_cutoff(a, b, 1, 0)
    assert math.isclose(fc_fe, fc_an, rel_tol=5e-3)                 # FE cutoff matches analytic
    # the dispersion built on the FE cutoff agrees with the analytic-cutoff dispersion
    d_fe = waveguide_dispersion(10e9, fc_fe)
    d_an = waveguide_dispersion(10e9, fc_an)
    assert math.isclose(d_fe["v_group"], d_an["v_group"], rel_tol=1e-2)
    assert math.isclose(d_fe["lambda_g"], d_an["lambda_g"], rel_tol=1e-2)


if __name__ == "__main__":
    test_dispersion_identities()
    test_cutoff_and_high_frequency_limits()
    test_evanescent_limits()
    test_dispersion_fe_anchor()
    print("[OK] waveguide dispersion: v_p*v_g=c^2, lambda_g relation, v_g==d omega/d beta, FE-anchored fc.")
