"""
Biot-Savart field computation via ngsolve.bem BiotSavartCF.

Provides H and B field from wire currents as NGSolve CoefficientFunctions.
All implementations use ngsolve.bem's multipole-accelerated BiotSavartCF.

Usage:
    from radia.biot_savart import biot_savart_wire, biot_savart_loop

    # Wire segments: list of (p1, p2) tuples
    H_cf = biot_savart_wire(segments, current=1.0)
    B_cf = MU0 * H_cf

    # Circular loop
    H_cf = biot_savart_loop(center, radius, normal, current=1.0)

    # Evaluate on any mesh
    B_at_point = (MU0 * H_cf)(mesh(x, y, z))

Note: BiotSavartCF returns H field (A/m), not B field (T).
      Multiply by mu_0 to get B.
      Uses multipole expansion: accurate for observation points outside
      the source region (r_obs > r_source from expansion center).
"""

import math
import numpy as np
from ngsolve.bem import BiotSavartCF
from ngsolve.bla import Vec3D


MU0 = 4e-7 * math.pi


def biot_savart_wire(segments, current=1.0, center=(0, 0, 0), rad=None,
                     order=15, n_quad=10):
    """Create H-field CoefficientFunction from wire current segments.

    Args:
        segments: list of ((x1,y1,z1), (x2,y2,z2)) wire segment endpoints [m]
        current: total current [A] (or complex for AC)
        center: expansion center [m] (should be near the wire centroid)
        rad: expansion radius [m] (None = auto, 1.5x max distance from center)
        order: multipole expansion order (higher = more accurate near source)
        n_quad: quadrature points per segment for AddCurrent

    Returns:
        CoefficientFunction (dim=3) giving H field [A/m]
    """
    segments = list(segments)
    if not segments:
        raise ValueError("No wire segments provided")

    center = np.asarray(center, dtype=float)

    # Auto-detect expansion radius
    if rad is None:
        all_pts = np.array([p for seg in segments for p in seg])
        dists = np.linalg.norm(all_pts - center[np.newaxis, :], axis=1)
        rad = float(max(np.max(dists) * 1.5, 1e-10))

    bs = BiotSavartCF(order=order, kappa=1e-10,
                      center=Vec3D(*center), rad=rad)

    for p1, p2 in segments:
        bs.AddCurrent(Vec3D(*p1), Vec3D(*p2), complex(current), n_quad)

    return bs


def biot_savart_loop(center=(0, 0, 0), radius=0.030, normal=(0, 0, 1),
                     current=1.0, n_segments=100, order=15, n_quad=10):
    """Create H-field CoefficientFunction from a circular current loop.

    Args:
        center: loop center [m]
        radius: loop radius [m]
        normal: loop normal direction (unit vector)
        current: total current [A]
        n_segments: number of straight segments approximating the circle
        order: multipole expansion order
        n_quad: quadrature points per segment

    Returns:
        CoefficientFunction (dim=3) giving H field [A/m]
    """
    center = np.asarray(center, dtype=float)
    normal = np.asarray(normal, dtype=float)
    normal = normal / np.linalg.norm(normal)

    # Build orthonormal frame
    if abs(normal[2]) < 0.9:
        t1 = np.cross(normal, [0, 0, 1])
    else:
        t1 = np.cross(normal, [1, 0, 0])
    t1 = t1 / np.linalg.norm(t1)
    t2 = np.cross(normal, t1)

    # Generate segments
    segments = []
    for i in range(n_segments):
        theta1 = 2 * math.pi * i / n_segments
        theta2 = 2 * math.pi * (i + 1) / n_segments
        p1 = center + radius * (math.cos(theta1) * t1 + math.sin(theta1) * t2)
        p2 = center + radius * (math.cos(theta2) * t1 + math.sin(theta2) * t2)
        segments.append((tuple(p1), tuple(p2)))

    return biot_savart_wire(segments, current=current, center=tuple(center),
                            order=order, n_quad=n_quad)


def biot_savart_coilbuilder(coil_path, current=1.0, order=15, n_quad=10):
    """Create H-field CoefficientFunction from CoilBuilder path.

    Args:
        coil_path: (N, 3) array of coil center line vertices [m]
        current: total current [A]
        order: multipole expansion order
        n_quad: quadrature points per segment

    Returns:
        CoefficientFunction (dim=3) giving H field [A/m]
    """
    path = np.asarray(coil_path, dtype=float)
    if path.ndim != 2 or path.shape[1] != 3:
        raise ValueError(f"coil_path must be (N, 3), got {path.shape}")

    segments = [(tuple(path[i]), tuple(path[i + 1]))
                for i in range(len(path) - 1)]
    center = tuple(path.mean(axis=0))

    return biot_savart_wire(segments, current=current, center=center,
                            order=order, n_quad=n_quad)
