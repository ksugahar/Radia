"""Guards: the curved device charge map reproduces NGSolve's default curved charge rules.

Assumes: NGSolve only, no GPU. The CPU curved route (``_charge_basis_curved``)
assembles its L2/SurfaceL2 masses and its moment forms with NGSolve's default
integration orders on a ``Curve(2)`` mesh. The device stages explicit
``IntegrationRule`` orders (``hdiv_charge_staging.hpp``); the masses are not
integrated exactly on curved cells, so the device must use the very rule
NGSolve picks. Device-vs-CPU B agreement is a separate GPU run.
"""
from pathlib import Path
import re

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")

ROOT = Path(__file__).resolve().parents[1]
STAGING = ROOT / "src/core/gpu/hdiv/hdiv_charge_staging.hpp"


def _orders(p):
    """Python mirror of the staging helpers, keyed by their C++ names."""
    pv = max(p - 1, 0)
    return dict(curvedCellMassOrder=2*pv, curvedFaceMassOrder=2*p,
                curvedCellMomentOrder=max(p + pv - 1, 0), curvedFaceMomentOrder=2*p + 1)


@pytest.fixture(scope="module")
def curved_ball():
    from netgen.occ import OCCGeometry, Pnt, Sphere
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=.6))
    mesh.Curve(2)
    assert mesh.GetCurveOrder() == 2
    return mesh


def _values(form):
    form.Assemble()
    rows, cols, vals = form.mat.COO()
    order = np.lexsort((np.asarray(cols), np.asarray(rows)))
    return np.asarray(vals)[order]


def test_staging_helpers_match_the_mirrored_orders():
    source = STAGING.read_text(encoding="utf-8")
    for p in (1, 2):
        for name, value in _orders(p).items():
            body = re.search(r"inline int %s\(int order\)\{return (.*?);\}" % name, source)
            assert body, name
            expression = body.group(1).replace("std::max", "max")
            actual = eval(expression, {"max": max, "order": p})  # noqa: S307 - pinned helper body
            assert actual == value, (name, p, actual, value)


@pytest.mark.parametrize("p", (1, 2))
def test_default_curved_charge_rules_are_the_staged_rules(curved_ball, p):
    mesh = curved_ball
    pv = max(p - 1, 0)
    orders = _orders(p)
    nn = ng.specialcf.normal(3)
    fes = ng.HDiv(mesh, order=p)
    L2v, L2b = ng.L2(mesh, order=pv), ng.SurfaceL2(mesh, order=p)
    tet = lambda order: {ng.TET: ng.IntegrationRule(ng.TET, order)}
    trig = lambda order: {ng.TRIG: ng.IntegrationRule(ng.TRIG, order)}

    def masses(dx, ds):
        mv = ng.BilinearForm(L2v); mv += L2v.TrialFunction()*L2v.TestFunction()*dx
        mb = ng.BilinearForm(L2b); mb += L2b.TrialFunction()*L2b.TestFunction()*ds
        return _values(mv), _values(mb)

    def moments(dx, ds):
        u = fes.TrialFunction()
        bv = ng.BilinearForm(trialspace=fes, testspace=L2v); bv += (-ng.div(u))*L2v.TestFunction()*dx
        bb = ng.BilinearForm(trialspace=fes, testspace=L2b); bb += (u.Trace()*nn)*L2b.TestFunction()*ds
        return _values(bv), _values(bb)

    default = masses(ng.dx, ng.ds) + moments(ng.dx, ng.ds)
    staged = masses(ng.dx(intrules=tet(orders["curvedCellMassOrder"])),
                    ng.ds(intrules=trig(orders["curvedFaceMassOrder"]))) + \
        moments(ng.dx(intrules=tet(orders["curvedCellMomentOrder"])),
                ng.ds(intrules=trig(orders["curvedFaceMomentOrder"])))
    for name, actual, expected in zip(("cell mass", "face mass", "cell moment", "face moment"), staged, default):
        assert actual.shape == expected.shape and np.array_equal(actual, expected), \
            (name, p, float(np.abs(actual - expected).max()) if actual.shape == expected.shape else actual.shape)
    # The face mass is under-integrated on curved faces, so its exact rule is load-bearing.
    exact_face_mass = masses(ng.dx, ng.ds(intrules=trig(2*p + 8)))[1]
    assert np.abs(exact_face_mass - default[1]).max() > 1e-10*np.abs(default[1]).max(), p
