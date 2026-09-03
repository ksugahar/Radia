"""CAD-shape dispatch for PEEC internal impedance models."""

import math

import pytest

from radia._b3d_shim import GeomType
from radia.coil_from_cad import _rectangular_face_dimensions
from radia.panels import calc_inductance


def test_strict_rectangle_selects_dowell_and_keeps_real_perimeter():
    topology = {
        "cross_section_kind": "rect",
        "cross_section_width_m_mean": 8.0e-3,
        "cross_section_thickness_m_mean": 1.0e-3,
        "cross_section_perimeter_m_mean": 18.0e-3,
    }
    model = calc_inductance._peec_cross_section_model(
        topology, 5.8e7, 2.0 * math.pi * 100.0e3, 16)
    assert model["name"] == "rectangular-dowell"
    assert model["perimeter_m"] == pytest.approx(18.0e-3)
    assert model["width_m"] == pytest.approx(8.0e-3)
    assert model["thickness_m"] == pytest.approx(1.0e-3)


def test_physical_circle_and_unknown_area_remain_distinguishable():
    circle = calc_inductance._peec_cross_section_model(
        {
            "cross_section_kind": "circle",
            "cross_section_radius_m_mean": 1.0e-3,
        },
        5.8e7, 2.0 * math.pi * 100.0e3, 16,
    )
    unknown = calc_inductance._peec_cross_section_model(
        {
            "cross_section_kind": "unknown",
            "cross_section_area_m2_mean": math.pi * 1.0e-6,
            "cross_section_perimeter_m_mean": 10.0e-3,
        },
        5.8e7, 2.0 * math.pi * 100.0e3, 16,
    )
    assert circle["name"] == "round-bessel"
    assert unknown["name"] == "equivalent-round-bessel"
    assert unknown["perimeter_m"] == pytest.approx(10.0e-3)


class _Edge:
    geom_type = GeomType.LINE

    def __init__(self, length):
        self.length = length


class _Wire:
    def __init__(self, lengths):
        self._edges = [_Edge(length) for length in lengths]

    def edges(self):
        return self._edges


class _Face:
    def __init__(self, lengths, area):
        self._wire = _Wire(lengths)
        self.area = area

    def outer_wire(self):
        return self._wire


def test_face_classifier_rejects_a_nonrectangular_four_edge_polygon():
    rectangle = _Face([8.0, 1.0, 8.0, 1.0], area=8.0)
    rhombus = _Face([8.0, 1.0, 8.0, 1.0], area=6.0)
    assert _rectangular_face_dimensions(rectangle) == pytest.approx((8.0, 1.0))
    assert _rectangular_face_dimensions(rhombus) is None
