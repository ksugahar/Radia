"""An axisymmetric mesh has two measures, and code must say which it meant.

``ngsolve.Integrate`` on a meridian mesh returns the plane integral
``int f dr dz``.  That is right for a current or a flux through the section
and wrong for anything living in the revolved solid, which needs
``int f 2 pi r dr dz``.  Both are one identical-looking call, the error is a
bias rather than noise, and it is smallest when the section sits far from the
axis -- that is, when the model looks most trustworthy.  This repository has
paid for that more than once.

So these tests do two things: check that the volume measure is arithmetically
what Pappus says it is, and refuse a bare ``Integrate`` in axisymmetric code,
so the next author has to name the measure rather than remember it.
"""
from __future__ import annotations

import ast
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

# Axisymmetric code: every integral here is over a meridian mesh.
GOVERNED = (
    ROOT / "src/radia/axisym_measure.py",
    ROOT / "src/radia/eddy_axisym_ring.py",
    ROOT / "validation_test/induction_heating/aphi_beak_section_sweep.py",
)
SANCTIONED = {"axi_volume_integral", "axi_section_integral",
              "axi_volume_beyond"}


def _integrate_calls(path):
    """Every ``Integrate``-like call in ``path``, as (line, qualified name)."""
    tree = ast.parse(path.read_text(encoding="utf-8"), str(path))
    found = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if isinstance(func, ast.Attribute):
            name = func.attr
        elif isinstance(func, ast.Name):
            name = func.id
        else:
            continue
        if name == "Integrate":
            found.append((node.lineno, name))
    return found


@pytest.mark.parametrize("path", GOVERNED, ids=lambda p: p.name)
def test_axisymmetric_code_names_its_measure(path):
    """No bare ``Integrate`` outside the module that defines the measures."""
    assert path.exists(), f"{path} moved; update the governed list"
    calls = _integrate_calls(path)
    if path.name == "axisym_measure.py":
        # The one place a bare Integrate is correct: it is what defines both
        # measures, and both are defined there.
        assert calls, "axisym_measure must actually call Integrate"
        source = path.read_text(encoding="utf-8")
        for helper in SANCTIONED:
            assert f"def {helper}(" in source, f"{helper} disappeared"
        return
    assert not calls, (
        f"{path.name} calls Integrate directly at line(s) "
        f"{[line for line, _ in calls]}.  On a meridian mesh that is the "
        f"plane measure.  Use radia.axisym_measure.axi_volume_integral for "
        f"anything in the revolved solid (power, energy, mass, volume) or "
        f"axi_section_integral for anything through the cross-section "
        f"(current, flux), so the choice is visible in the source.")


def test_volume_measure_is_pappus():
    """``axi_volume_integral(1)`` over a rectangle is the torus it revolves to.

    Pappus: the volume is the section area times the circumference travelled
    by its centroid.  A plain meridian integral would return the area instead,
    which is the mistake this module exists to stop.
    """
    ng = pytest.importorskip("ngsolve")
    from netgen.occ import MoveTo, OCCGeometry

    from radia.axisym_measure import axi_section_integral, axi_volume_integral

    r0, width, height = 0.4, 0.008, 0.004
    face = MoveTo(r0, -0.5 * height).Rectangle(width, height).Face()
    face.faces.name = "ring"
    mesh = ng.Mesh(OCCGeometry(face, dim=2).GenerateMesh(maxh=0.001))

    area = float(axi_section_integral(ng.CF(1.0), mesh).real)
    volume = float(axi_volume_integral(ng.CF(1.0), mesh).real)
    pappus = width * height * 2.0 * math.pi * (r0 + 0.5 * width)

    assert area == pytest.approx(width * height, rel=1e-10)
    assert volume == pytest.approx(pappus, rel=1e-10)
    # And the two are genuinely different, by the factor that was being lost.
    assert volume / area == pytest.approx(
        2.0 * math.pi * (r0 + 0.5 * width), rel=1e-10)


def test_volume_measure_weights_the_outer_part_more():
    """The half at larger radius carries more volume than the inner half.

    This is the asymmetry a plain meridian integral erases, and it is the one
    that biased the beak-versus-body loss split: the beak sits at the larger
    radius.  The two halves are separate regions rather than one region cut by
    ``IfPos``, so this measures the measure and not the quadrature's view of a
    discontinuity.
    """
    ng = pytest.importorskip("ngsolve")
    from netgen.occ import Glue, MoveTo, OCCGeometry

    from radia.axisym_measure import axi_section_integral, axi_volume_integral

    r0, width, height = 0.4, 0.008, 0.004
    half = 0.5 * width
    inner = MoveTo(r0, -0.5 * height).Rectangle(half, height).Face()
    inner.faces.name = "inner"
    outer = MoveTo(r0 + half, -0.5 * height).Rectangle(half, height).Face()
    outer.faces.name = "outer"
    mesh = ng.Mesh(OCCGeometry(Glue([inner, outer]),
                               dim=2).GenerateMesh(maxh=0.001))

    one = ng.CF(1.0)
    v_in = float(axi_volume_integral(
        one, mesh, definedon=mesh.Materials("inner")).real)
    v_out = float(axi_volume_integral(
        one, mesh, definedon=mesh.Materials("outer")).real)
    a_in = float(axi_section_integral(
        one, mesh, definedon=mesh.Materials("inner")).real)
    a_out = float(axi_section_integral(
        one, mesh, definedon=mesh.Materials("outer")).real)

    # The meridian measure cannot tell the two halves apart at all.
    assert a_in == pytest.approx(a_out, rel=1e-10)
    # The volume measure gives each half the circumference of its own centroid.
    assert v_in == pytest.approx(
        half * height * 2.0 * math.pi * (r0 + 0.25 * width), rel=1e-10)
    assert v_out == pytest.approx(
        half * height * 2.0 * math.pi * (r0 + 0.75 * width), rel=1e-10)
    assert v_out > v_in
