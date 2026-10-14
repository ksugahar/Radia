"""Reusable TEAM21a geometry and opposed-coil source, from the public 2009 TEAM specification.

Geometry/source validation only: no eddy solve or benchmark acceptance is claimed.
Source: https://www.compumag.org/wp/wp-content/uploads/2018/06/problem-21-family-v-2009.pdf

Recovered from the dated 2026-07-18 fixture; obsolete solver paths are omitted.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

import radia


@dataclass(frozen=True)
class Team21ASpec:
    frequency_hz: float = 50.0
    coil_turns: int = 300
    coil_current_rms_A: float = 10.0
    plate_width_m: float = 0.360
    plate_length_m: float = 0.820
    plate_thickness_m: float = 0.010
    slit_length_m: float = 0.660
    slit_width_m: float = 0.010
    conductivity_S_per_m: float = 1.3889e6
    relative_permeability: float = 1.0
    coil_outer_size_m: float = 0.270
    coil_inner_size_m: float = 0.200
    coil_outer_corner_radius_m: float = 0.045
    coil_inner_corner_radius_m: float = 0.010
    coil_axial_length_m: float = 0.217
    coil_axial_gap_m: float = 0.024
    coil_to_plate_gap_m: float = 0.012

    @property
    def current_peak_A(self) -> float:
        return math.sqrt(2.0) * self.coil_current_rms_A

    @property
    def skin_depth_m(self) -> float:
        omega = 2.0 * math.pi * self.frequency_hz
        return math.sqrt(
            2.0
            / (
                omega
                * (4.0e-7 * math.pi)
                * self.relative_permeability
                * self.conductivity_S_per_m
            )
        )


SPEC = Team21ASpec()
def slit_centers_m(model: int) -> tuple[float, ...]:
    """Return the equally spaced slit centerlines from Figure A1-2."""

    centers = {
        0: (),
        1: (0.0,),
        2: (-0.060, 0.060),
        3: (-0.090, 0.0, 0.090),
    }
    try:
        return centers[int(model)]
    except KeyError as exc:
        raise ValueError("TEAM 21a model must be 0, 1, 2, or 3") from exc


def _rounded_loop_points(
    size: float,
    radius: float,
    *,
    center_x: float,
    center_y: float,
    arc_segments: int,
) -> list[list[float]]:
    half = 0.5 * size
    corners = (
        (center_x + half - radius, center_y + half - radius, 0.0, 0.5 * math.pi),
        (
            center_x - half + radius,
            center_y + half - radius,
            0.5 * math.pi,
            math.pi,
        ),
        (
            center_x - half + radius,
            center_y - half + radius,
            math.pi,
            1.5 * math.pi,
        ),
        (
            center_x + half - radius,
            center_y - half + radius,
            1.5 * math.pi,
            2.0 * math.pi,
        ),
    )
    points: list[list[float]] = []
    for cx, cy, angle0, angle1 in corners:
        angles = np.linspace(angle0, angle1, arc_segments + 1)
        arc = [
            [cx + radius * math.cos(angle), cy + radius * math.sin(angle)]
            for angle in angles
        ]
        points.extend(arc)
    points.append(points[0])
    return points


def build_mother_source(
    *,
    radial_order: int = 3,
    axial_order: int = 5,
    arc_segments: int = 5,
):
    """Build the two distributed TEAM coils and their separate source handles."""

    radia.UtiDelAll()
    radial_nodes, radial_weights = np.polynomial.legendre.leggauss(radial_order)
    axial_nodes, axial_weights = np.polynomial.legendre.leggauss(axial_order)
    radial_weights *= 0.5
    axial_weights *= 0.5
    outer_near_x = 0.5 * SPEC.plate_thickness_m + SPEC.coil_to_plate_gap_m
    loop_center_x = outer_near_x + 0.5 * SPEC.coil_outer_size_m
    axial_center = 0.5 * (SPEC.coil_axial_gap_m + SPEC.coil_axial_length_m)
    ampere_turns_peak = SPEC.coil_turns * SPEC.current_peak_A
    coil_containers = []
    for sign, z_center in ((1.0, -axial_center), (-1.0, axial_center)):
        filaments = []
        for radial_node, radial_weight in zip(radial_nodes, radial_weights):
            fraction = 0.5 * (radial_node + 1.0)
            size = SPEC.coil_inner_size_m + (
                SPEC.coil_outer_size_m - SPEC.coil_inner_size_m
            ) * fraction
            radius = SPEC.coil_inner_corner_radius_m + (
                SPEC.coil_outer_corner_radius_m
                - SPEC.coil_inner_corner_radius_m
            ) * fraction
            loop_xy = _rounded_loop_points(
                size,
                radius,
                center_x=loop_center_x,
                center_y=0.0,
                arc_segments=arc_segments,
            )
            for axial_node, axial_weight in zip(axial_nodes, axial_weights):
                z_value = z_center + 0.5 * SPEC.coil_axial_length_m * axial_node
                points = [[x_value, y_value, z_value] for x_value, y_value in loop_xy]
                filaments.append(
                    radia.ObjFlmCur(
                        points,
                        sign
                        * ampere_turns_peak
                        * float(radial_weight)
                        * float(axial_weight),
                    )
                )
        coil_containers.append(radia.ObjCnt(filaments))
    source = radia.ObjCnt(coil_containers)
    return source, tuple(coil_containers)


def _conductor_shape(model: int):
    import netgen.occ as occ

    plate = occ.Box(
        occ.Pnt(
            -0.5 * SPEC.plate_thickness_m,
            -0.5 * SPEC.plate_width_m,
            -0.5 * SPEC.plate_length_m,
        ),
        occ.Pnt(
            0.5 * SPEC.plate_thickness_m,
            0.5 * SPEC.plate_width_m,
            0.5 * SPEC.plate_length_m,
        ),
    )
    for center in slit_centers_m(model):
        cutter = occ.Box(
            occ.Pnt(
                -0.6 * SPEC.plate_thickness_m,
                center - 0.5 * SPEC.slit_width_m,
                -0.5 * SPEC.slit_length_m,
            ),
            occ.Pnt(
                0.6 * SPEC.plate_thickness_m,
                center + 0.5 * SPEC.slit_width_m,
                0.5 * SPEC.slit_length_m,
            ),
        )
        plate = plate - cutter
    plate.mat("cond")
    for face in plate.faces:
        face.name = "skin"
    return plate


def build_conductor_mesh(model: int, maxh: float):
    import netgen.occ as occ
    import ngsolve as ng

    return ng.Mesh(
        occ.OCCGeometry(_conductor_shape(model)).GenerateMesh(
            maxh=float(maxh),
            grading=0.35,
        )
    )



