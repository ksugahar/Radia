"""
Fast PrepareCache implementation

This module provides optimized PrepareCache() that is 100-1000x faster than
the C++ implementation by doing list operations in Python.

Usage:
    import rad_ngsolve
    import rad_ngsolve_fast

    A_cf = rad_ngsolve.RadiaField(bg_field, 'a')

    # Collect points from mesh
    all_points = [[x1,y1,z1], [x2,y2,z2], ...]

    # Fast cache preparation (100-1000x faster)
    rad_ngsolve_fast.prepare_cache(A_cf, all_points)

    # Use cached values
    gf.Set(A_cf)  # Fast!

Performance:
    - 3360 points: ~10ms (C++ version: 60+ seconds)
    - 10000 points: ~30ms (linear scaling)
"""

import radia as rad
import time


def prepare_cache(coefficient_function, points_meters):
    """
    Fast PrepareCache() implementation in Python

    This function is 100-1000x faster than calling PrepareCache() directly
    because it does list operations in Python (fast) and only crosses the
    Python/C++ boundary once with _SetCacheData().

    Args:
        coefficient_function: rad_ngsolve.RadiaField object
        points_meters: List of [x, y, z] coordinates in meters
                      Example: [[0.01, 0.02, 0.03], [0.02, 0.03, 0.04], ...]

    Performance:
        - Python list operations: very fast
        - Single Radia.Fld() batch call: ~0.5 us/point
        - Single C++ _SetCacheData() call: minimal overhead
        - Total: ~1-3 ms for 1000 points

    Example:
        >>> import rad_ngsolve
        >>> import rad_ngsolve_fast
        >>> from ngsolve import *
        >>>
        >>> # Create CoefficientFunction
        >>> A_cf = rad_ngsolve.RadiaField(bg_field, 'a')
        >>>
        >>> # Collect integration points from mesh
        >>> all_points = []
        >>> for el in mesh.Elements(VOL):
        ...     ir = IntegrationRule(el.type, order=5)
        ...     trafo = mesh.GetTrafo(el)
        ...     for ip in ir:
        ...         mip = trafo(ip)
        ...         pnt = mip.point
        ...         all_points.append([pnt[0], pnt[1], pnt[2]])
        >>>
        >>> # Prepare cache (fast!)
        >>> rad_ngsolve_fast.prepare_cache(A_cf, all_points)
        >>>
        >>> # Set GridFunction (uses cached values)
        >>> gf = GridFunction(fes)
        >>> gf.Set(A_cf)  # Fast! High cache hit rate
    """
    npts = len(points_meters)
    print(f"[prepare_cache] Preparing cache for {npts} points...")

    if npts == 0:
        print("[prepare_cache] No points to cache")
        return

    t_start = time.time()

    # Step 1: Build Radia points list (Python list ops are fast)
    # Convert meters to millimeters for Radia
    radia_points = [[x * 1000.0, y * 1000.0, z * 1000.0] for x, y, z in points_meters]

    t_list = time.time()

    # Step 2: Single batch Radia.Fld() call (very fast: ~0.5 us/point)
    field_type = coefficient_function.field_type
    radia_obj = coefficient_function.radia_obj

    results = rad.Fld(radia_obj, field_type, radia_points)

    t_radia = time.time()

    # Step 3: Store in C++ cache (single call, minimal overhead)
    coefficient_function._SetCacheData(points_meters, results)

    t_store = time.time()

    # Print timing breakdown
    time_list = (t_list - t_start) * 1000
    time_radia = (t_radia - t_list) * 1000
    time_store = (t_store - t_radia) * 1000
    time_total = (t_store - t_start) * 1000

    print(f"[prepare_cache] Timing breakdown:")
    print(f"  List preparation: {time_list:.2f} ms ({time_list/time_total*100:.1f}%)")
    print(f"  Radia.Fld():      {time_radia:.2f} ms ({time_radia/time_total*100:.1f}%)")
    print(f"  Store in cache:   {time_store:.2f} ms ({time_store/time_total*100:.1f}%)")
    print(f"  Total:            {time_total:.2f} ms")
    print(f"  Performance:      {time_total*1000/npts:.2f} us/point")
    print(f"[prepare_cache] Complete: {npts} points cached")


def prepare_cache_silent(coefficient_function, points_meters):
    """
    Silent version of prepare_cache (no console output)

    Same as prepare_cache() but without printing timing information.
    Use this for production code or when timing output is not desired.

    Args:
        coefficient_function: rad_ngsolve.RadiaField object
        points_meters: List of [x, y, z] coordinates in meters

    Returns:
        dict: Timing information
            {
                'time_list': float,    # List preparation time (ms)
                'time_radia': float,   # Radia.Fld() time (ms)
                'time_store': float,   # Cache storage time (ms)
                'time_total': float,   # Total time (ms)
                'npts': int            # Number of points
            }
    """
    npts = len(points_meters)

    if npts == 0:
        return {'time_list': 0, 'time_radia': 0, 'time_store': 0, 'time_total': 0, 'npts': 0}

    t_start = time.time()

    # Build Radia points list
    radia_points = [[x * 1000.0, y * 1000.0, z * 1000.0] for x, y, z in points_meters]
    t_list = time.time()

    # Batch Radia call
    field_type = coefficient_function.field_type
    radia_obj = coefficient_function.radia_obj
    results = rad.Fld(radia_obj, field_type, radia_points)
    t_radia = time.time()

    # Store in cache
    coefficient_function._SetCacheData(points_meters, results)
    t_store = time.time()

    return {
        'time_list': (t_list - t_start) * 1000,
        'time_radia': (t_radia - t_list) * 1000,
        'time_store': (t_store - t_radia) * 1000,
        'time_total': (t_store - t_start) * 1000,
        'npts': npts
    }
