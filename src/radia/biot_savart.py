"""
Biot-Savart field computation.

Two APIs:
  1. Point evaluation (PEEC coupling, ObjBckg callbacks):
     h_filament(p1, p2, obs) -> H [A/m] at a single point (analytical, exact)

  2. CoefficientFunction (FEM source, IH coil field):
     biot_savart_wire(segments) -> CF(dim=3) H field (BiotSavartCF multipole)
     biot_savart_loop(center, radius) -> CF
     biot_savart_coilbuilder(path) -> CF

BiotSavartCF returns H field (A/m). Multiply by mu_0 for B (Tesla).
"""

import math
import numpy as np


MU0 = 4e-7 * math.pi
INV_4PI = 1.0 / (4.0 * math.pi)


# ============================================================
# Point evaluation (analytical, no mesh required)
# ============================================================

def h_filament(p1, p2, obs_point, current=1.0):
    """H-field at obs_point from a finite current filament (p1 -> p2).

    Exact analytical formula:
        H = I/(4*pi*d) * (cos(alpha1) - cos(alpha2)) * e_perp

    Args:
        p1, p2: filament endpoints [x, y, z] in meters
        obs_point: observation point [x, y, z] in meters
        current: current [A] (default 1.0)

    Returns:
        H-field vector (3,) in A/m
    """
    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)
    r = np.asarray(obs_point, dtype=float)

    dl = p2 - p1
    L = np.linalg.norm(dl)
    if L < 1e-30:
        return np.zeros(3)
    e_l = dl / L

    r1 = r - p1
    cross = np.cross(e_l, r1)
    d = np.linalg.norm(cross)
    if d < 1e-30:
        return np.zeros(3)

    e_perp = cross / d
    r1_mag = np.linalg.norm(r1)
    r2_mag = np.linalg.norm(r - p2)
    if r1_mag < 1e-30 or r2_mag < 1e-30:
        return np.zeros(3)

    cos_a1 = np.dot(e_l, r1) / r1_mag
    cos_a2 = np.dot(e_l, r - p2) / r2_mag

    return current * INV_4PI / d * (cos_a1 - cos_a2) * e_perp


def h_segments(segments, obs_point, current=1.0):
    """H-field at obs_point from multiple wire segments.

    Args:
        segments: list of (p1, p2) endpoint tuples
        obs_point: observation point [x, y, z]
        current: current [A]

    Returns:
        H-field vector (3,) in A/m
    """
    H = np.zeros(3)
    for p1, p2 in segments:
        H += h_filament(p1, p2, obs_point, current)
    return H


def a_filament(p1, p2, obs_point, current=1.0):
    """Vector potential A at obs_point from a finite current filament.

    A = mu0/(4*pi) * I * direction * [arsinh((L-t)/d) + arsinh(t/d)]

    Args:
        p1, p2: filament endpoints [x, y, z] in meters
        obs_point: observation point [x, y, z] in meters
        current: current [A]

    Returns:
        A vector (3,) in T*m
    """
    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)
    r = np.asarray(obs_point, dtype=float)

    dl = p2 - p1
    L = np.linalg.norm(dl)
    if L < 1e-30:
        return np.zeros(3)
    e_l = dl / L

    r1 = r - p1
    t = np.dot(r1, e_l)
    perp = r1 - t * e_l
    d = np.linalg.norm(perp)
    if d < 1e-15:
        d = 1e-15

    integral = np.arcsinh((L - t) / d) + np.arcsinh(t / d)
    return current * MU0 * INV_4PI * integral * e_l


# ============================================================
# CoefficientFunction API (BiotSavartCF multipole, needs mesh)
# ============================================================


def biot_savart_wire(segments, current=1.0, center=(0, 0, 0), rad=None,
                     order=15, n_quad=10):
    """Create H-field CoefficientFunction from wire current segments.

    Uses ngsolve.bem BiotSavartCF (multipole-accelerated).

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
    from ngsolve.bem import BiotSavartCF
    from ngsolve.bla import Vec3D

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
