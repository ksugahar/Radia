"""Small geometry checks for curved-boundary FEM meshes.

The helpers here quantify the first-order error made when a circular boundary
is represented by straight chords.  They are useful as fast, license-free gates
before spending heavier mesh-generation time on high-order curved exports.
"""

import math

TAU = 2.0 * math.pi


def circular_arc_chord_area_deficit(radius, angle):
    """Area between a circular arc and its straight chord.

    Parameters
    ----------
    radius : float
        Circle radius.
    angle : float
        Arc angle in radians.  The practical non-self-intersecting mesh segment
        range is ``0 < angle <= pi``.

    Returns
    -------
    float
        Sector area minus chord triangle area:
        ``0.5 * radius**2 * (angle - sin(angle))``.
    """
    r = float(radius)
    theta = float(angle)
    if r <= 0.0:
        raise ValueError("radius must be positive")
    if theta <= 0.0 or theta > math.pi:
        raise ValueError("angle must be in (0, pi]")
    return 0.5 * r * r * (theta - math.sin(theta))


def circle_polyline_area_metrics(radius, segment_angles):
    """Area metrics for a circle sector represented by chord segments.

    ``segment_angles`` are positive arc angles in radians and must sum to at
    most one full circle.  For a closed full-circle boundary use angles that sum
    to ``2*pi``.  The returned ``faceted_area`` is the sum of signed chord
    triangles, and ``area_deficit`` is the missing curved cap area.
    """
    r = float(radius)
    if r <= 0.0:
        raise ValueError("radius must be positive")
    angles = [float(a) for a in segment_angles]
    if not angles:
        raise ValueError("segment_angles must be non-empty")
    if any(a <= 0.0 or a > math.pi for a in angles):
        raise ValueError("each segment angle must be in (0, pi]")
    total_angle = sum(angles)
    if total_angle > TAU * (1.0 + 1e-12):
        raise ValueError("segment angles must not exceed a full circle")

    exact_area = 0.5 * r * r * total_angle
    faceted_area = 0.5 * r * r * sum(math.sin(a) for a in angles)
    deficit = exact_area - faceted_area
    return {
        "radius": r,
        "segment_count": len(angles),
        "total_angle": total_angle,
        "exact_area": exact_area,
        "faceted_area": faceted_area,
        "area_deficit": deficit,
        "relative_deficit": deficit / exact_area,
    }


def regular_polygon_circle_area_metrics(radius, n_sides):
    """Area metrics for an inscribed regular ``n_sides`` polygon."""
    n = int(n_sides)
    if n < 3:
        raise ValueError("n_sides must be at least 3")
    theta = TAU / n
    metrics = circle_polyline_area_metrics(radius, [theta] * n)
    metrics["n_sides"] = n
    return metrics


def regular_polygon_circle_area(radius, n_sides):
    """Area of an inscribed regular polygon approximating a circle."""
    return regular_polygon_circle_area_metrics(radius, n_sides)["faceted_area"]


def required_regular_polygon_sides_for_area_error(max_relative_deficit):
    """Smallest regular polygon side count whose circle area deficit is below a target."""
    target = float(max_relative_deficit)
    if target <= 0.0 or target >= 1.0:
        raise ValueError("max_relative_deficit must be in (0, 1)")

    lo = 3
    hi = 3
    while regular_polygon_circle_area_metrics(1.0, hi)["relative_deficit"] > target:
        hi *= 2

    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if regular_polygon_circle_area_metrics(1.0, mid)["relative_deficit"] <= target:
            hi = mid
        else:
            lo = mid
    return hi
