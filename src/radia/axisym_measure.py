"""Integration measures for axisymmetric meridian meshes.

An axisymmetric model is solved on a two-dimensional ``(r, z)`` mesh, and that
mesh carries two different integration measures that look identical in code:

* the **meridian** measure ``dr dz``, which is what a current or a flux
  *through the cross-section* is integrated in, and
* the **volume** measure ``2 pi r dr dz``, which is what any quantity living
  in the revolved solid -- energy, dissipation, mass, volume -- is integrated
  in.

``ngsolve.Integrate`` on such a mesh gives the meridian measure, silently.  A
caller who wants dissipation and writes the obvious call gets an answer that
is wrong by a factor varying across the section, which is a bias rather than
noise and does not announce itself: it is small when the section is far from
the axis, which is exactly when the model looks most trustworthy.  That has
now cost this repository the same defect more than once.

So there is no correct bare ``Integrate`` on an axisymmetric mesh here.  Use
:func:`axi_volume_integral` or :func:`axi_section_integral`, whose names state
which measure was meant, and let ``tests/test_axisym_measure_contract.py``
keep it that way.
"""
from __future__ import annotations

import numpy as np


def axi_volume_integral(cf, mesh, *, definedon=None, order=10):
    """Integrate ``cf`` over the revolved solid: ``int cf 2 pi r dr dz``.

    This is the measure for every quantity that lives in the three-dimensional
    body -- dissipated power, stored energy, volume, mass, heat content.

    Args:
        cf: CoefficientFunction to integrate.
        mesh: two-dimensional ``(r, z)`` mesh with ``r`` the first coordinate.
        definedon: optional region, e.g. ``mesh.Materials("conductor")``.
        order: integration order.

    Returns:
        complex or float, as ``cf`` dictates.
    """
    import ngsolve as ng

    return ng.Integrate(cf * 2.0 * np.pi * ng.x, mesh, definedon=definedon,
                        order=order)


def axi_section_integral(cf, mesh, *, definedon=None, order=10):
    """Integrate ``cf`` over the meridian section: ``int cf dr dz``.

    This is the measure for a quantity carried *through* the cross-section --
    a terminal current, a flux through the meridian, a section area.  It is
    the plain two-dimensional integral, named so that choosing it is visible.
    """
    import ngsolve as ng

    return ng.Integrate(cf, mesh, definedon=definedon, order=order)


def axi_volume_beyond(cf, mesh, r_m, *, definedon=None, order=10):
    """:func:`axi_volume_integral` restricted to radii beyond ``r_m``."""
    import ngsolve as ng

    return axi_volume_integral(cf * ng.IfPos(ng.x - float(r_m), 1.0, 0.0),
                               mesh, definedon=definedon, order=order)
