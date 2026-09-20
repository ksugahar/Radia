"""Finite rectangular-section coil: analytic field, not a filament limit."""
import numpy as np
import pytest
import radia as rad
import ngsolve as ng
from netgen.occ import Box, OCCGeometry, Pnt
from radia.coil_builder import CoilBuilder


@pytest.mark.parametrize("scale", [0.1, 1.0, 10.0])
@pytest.mark.parametrize("current", [1000.0, -1000.0])
@pytest.mark.xfail(strict=True, reason="unmet premise: a systematic 6.7787e-04 relative error remains at scale=1.0 with FldLenRndSw(off), identical for arc_max_segment_length default, /4, /16 and /64, so it is not discretisation; it scales as 1/scale (6.78e-3, 6.78e-4, 6.80e-5 at scale 0.1, 1.0, 10.0), the signature of an absolute length effect that FldLenRndSw(off) reduces by only about 4e-8. Drop this marker together with whatever still perturbs the finite-section ring.")
def test_finite_section_axis_field_without_length_perturbation(scale, current):
    r1, r2, half = scale * np.array([.010, .012, .001])
    positions = scale * np.array([-.04, -.02, -.01, -.005, 0, .005, .01, .02, .04])
    density = current / ((r2-r1)*2*half)

    def primitive(u):
        return 0.0 if u == 0 else u*(np.arcsinh(r2/abs(u))-np.arcsinh(r1/abs(u)))

    expected = np.array([
        [0., -2*np.pi*1e-7*density*(primitive(y+half)-primitive(y-half)), 0.]
        for y in positions])
    mesh = ng.Mesh(OCCGeometry(Box(Pnt(*([-scale*.06]*3)),
                                  Pnt(*([scale*.06]*3)))).GenerateMesh(maxh=scale*.08))
    # Both construction and evaluation must use the same precision policy.
    rad.FldLenRndSw("off")
    try:
        coil = (CoilBuilder(current).set_start([(r1+r2)/2,0,0],
                orientation=[[1,0,0],[0,0,1],[0,-1,0]])
                .set_cross_section(r2-r1,2*half).add_arc((r1+r2)/2,360))
        field = rad.RadiaField(rad.ObjCnt(coil.to_radia()), "b")
        actual = np.array([field(mesh(0,float(y),0)) for y in positions])
        error = np.linalg.norm(actual-expected,axis=1)/np.linalg.norm(expected,axis=1)
        assert np.max(error) < 2e-12
    finally:
        rad.FldLenRndSw("on")
