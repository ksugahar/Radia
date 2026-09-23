"""Partial-element geometry shared by the fin PEEC gate-2 measurements.

Two frames appear here.  The kernel's frame is what
`MutualInductanceRectBar` constructs when it is reached at all: it takes the
segment axis and a global fallback axis, because `PEECSegment` carries no
cross-section orientation.  The surface frame is the physical one for a band
on a folded sheet -- width along the perimeter tangent, height along the
outward normal.

`filament_mutual` is Grover's closed form for two equal, aligned, parallel
filaments.  It is not an approximation of the kernel: the built L matrix
reproduces it to machine precision for the fin's segments, because those
segments are unsubdivided and therefore never reach the averaged kernel.
"""
from __future__ import annotations

import math

import numpy as np

MU0 = 4e-7 * math.pi
GAUSS3_PTS = np.array([-math.sqrt(3.0 / 5.0), 0.0, math.sqrt(3.0 / 5.0)])
GAUSS3_WTS = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])


def kernel_frame(direction):
    """The frame `MutualInductanceRectBar` builds, reproduced exactly."""
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d)
    ref = (np.array([1.0, 0.0, 0.0]) if abs(d[0]) < 0.9
           else np.array([0.0, 1.0, 0.0]))
    e_w = np.cross(ref, d)
    e_w /= np.linalg.norm(e_w)
    e_h = np.cross(d, e_w)
    return e_w, e_h


def surface_frame(xy_prev, xy, xy_next):
    """Width along the perimeter tangent, height along the outward normal."""
    tangent = np.asarray(xy_next, dtype=float) - np.asarray(xy_prev, dtype=float)
    tangent = tangent / np.linalg.norm(tangent)
    normal = np.array([tangent[1], -tangent[0]])
    if np.dot(normal, np.asarray(xy, dtype=float)) < 0:
        normal = -normal
    return (np.array([tangent[0], tangent[1], 0.0]),
            np.array([normal[0], normal[1], 0.0]))


def filament_mutual(length, offset):
    """Grover's equal, aligned, parallel filaments."""
    if offset <= 0.0:
        raise ValueError("coincident filaments have no finite mutual")
    ratio = length / offset
    return (MU0 * length / (2.0 * math.pi)) * (
        math.log(ratio + math.sqrt(1.0 + ratio * ratio))
        - math.sqrt(1.0 + 1.0 / (ratio * ratio)) + 1.0 / ratio)


def bar_mutual(center_i, center_j, direction, length, w_i, h_i, w_j, h_j,
               frame_i, frame_j):
    """Volume-averaged mutual, 3x3 Gauss per cross-section.

    This is the scheme `MutualInductanceRectBar` uses, evaluated here with
    whichever frames the caller supplies so the orientation can be varied.
    """
    e_wi, e_hi = frame_i
    e_wj, e_hj = frame_j
    axis = np.asarray(direction, dtype=float)
    axis = axis / np.linalg.norm(axis)
    total = 0.0
    weight = 0.0
    for a, wa in zip(GAUSS3_PTS, GAUSS3_WTS):
        for b, wb in zip(GAUSS3_PTS, GAUSS3_WTS):
            p_i = (np.asarray(center_i, dtype=float)
                   + 0.5 * w_i * a * e_wi + 0.5 * h_i * b * e_hi)
            for c, wc in zip(GAUSS3_PTS, GAUSS3_WTS):
                for d, wd in zip(GAUSS3_PTS, GAUSS3_WTS):
                    p_j = (np.asarray(center_j, dtype=float)
                           + 0.5 * w_j * c * e_wj + 0.5 * h_j * d * e_hj)
                    offset_vec = p_j - p_i
                    perp = offset_vec - np.dot(offset_vec, axis) * axis
                    offset = float(np.linalg.norm(perp))
                    if offset < 1e-12:
                        continue
                    ww = wa * wb * wc * wd
                    total += ww * filament_mutual(length, offset)
                    weight += ww
    return total / weight if weight else 0.0
