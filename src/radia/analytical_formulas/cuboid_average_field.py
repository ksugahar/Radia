"""Average magnetic flux density over a target rectangular box.

Part 6 §7 of the Wakao-Igarashi-Fujiwara-Kameari series. Originally
derived for micromagnetics: average ``B`` in one cubic lattice cell
due to the magnetisation of another lattice cell (same size).

Implementation status (2026-05-01)
----------------------------------
The **closed-form** path is now the default. The 3-fold antiderivative
primitives ``G1`` (diagonal demag) and ``G2`` (off-diagonal demag) were
derived by direct sympy integration of ``-d^2/dX^2 (1/R)`` and
``-d^2/dXdY (1/R)``; the resulting 64-corner sum on the displacement
box gives the spatial-averaged demagnetisation tensor exactly (no
quadrature error). Verified against the previous Gauss-Legendre
numerical path at ``rel_err ~3e-14`` over 8 source/target geometries.

The previously deferred Wakao Part 6 §7 OCR ambiguity is resolved
by the sympy derivation; Stafl §3.4 (initially named as the "parent
reference") is in fact a 2-D rectangular conductor and is unrelated
to the 3-D cuboid magnetisation problem. See
``memory/project_phase_alpha_g1_g2_derived.md`` for the derivation
record.

The Gauss-Legendre numerical path is preserved as
:func:`average_B_in_box_numerical` for diagnostic / cross-validation
use.

Closed-form summary
-------------------

For uniform magnetisation **M** in the source cuboid, the spatial
average of **B** over a target cuboid is

::

    <B>_T = mu_0 * (A_T @ M  +  M * V_overlap / V_T)

where ``A_T`` is the 3x3 average-demag-tensor obtained from a 64-corner
sum of the antiderivatives ``G1, G2``, and ``V_overlap`` is the volume
of the source-target intersection (zero for non-overlapping cuboids).

Useful for
----------
* Micromagnetics on a regular cubic lattice: cell-cell magnetic
  interaction.
* Mesh-to-mesh transfer between two rectangular grids.
* FEM-magnetic coupling: cell-averaged B for a rectangular FEM element.

Numerical caveat
----------------
The closed-form 64-corner sum suffers from catastrophic cancellation
when ``V_T`` is several orders of magnitude smaller than the source
box: the corner-primitive values are ``O(1)`` and the alternating
sum needs to evaluate to a much smaller number, with ULP roundoff
dominating. Keep target and source box dimensions of similar order;
for true *point* field evaluation (not box average), call
:meth:`radia.analytical_magnet.CuboidMagnet.get_B` directly.

References
----------
Wakao S., Fujiwara K., Tokumasu T., Kameari A.,
  "Useful Formulas of Analytical Integration in Electromagnetic Field
  Computations (Part 6)", IEE Japan SA-05-15 / RM-05-15 (2005),
  §7 eq 53-56.

Newell A. J., Williams W., Dunlop D. J.,
  "A generalization of the demagnetizing tensor for nonuniform
  magnetization", J. Geophys. Res. 98 (1993), pp. 9551-9555.

Yang Z. J. et al., "Potential and force between a magnet and a bulk
  Y1Ba2Cu3O7 superconductor studied by a mechanical pendulum",
  Supercond. Sci. Technol. 3 (1990), 591. -- The point B(x_p) formula
  used inside :class:`radia.analytical_magnet.CuboidMagnet`.
"""

from __future__ import annotations

import math

import numpy as np


MU_0 = 4.0e-7 * math.pi
_EPS = 1.0e-15  # cutoff for X=0 / Y=0 / Z=0 limits in atan/asinh terms


# ---------------------------------------------------------------------------
# Closed-form 3-fold antiderivative primitives (Phase α derivation)
# ---------------------------------------------------------------------------
def _G1(X, Y, Z):
    """Antiderivative of ``g1 = atan(YZ/(X*R))`` -- diagonal-demag corner primitive.

    Used for ``<N_xx>_T`` 64-corner sum. For ``<N_yy>_T`` call ``_G1(Y, X, Z)``;
    for ``<N_zz>_T`` call ``_G1(Z, X, Y)``.
    """
    R = math.sqrt(X * X + Y * Y + Z * Z)
    # XYZ * atan(YZ/(X*R)): 0 in the X -> 0 limit (XYZ vanishes faster than atan diverges)
    if abs(X) < _EPS:
        atan_term = 0.0
    else:
        atan_term = X * Y * Z * math.atan(Y * Z / (X * R))
    if (X * X + Z * Z) > 0:
        asinh_y = (Y / 2.0) * (X * X - Z * Z) * math.asinh(Y / math.sqrt(X * X + Z * Z))
    else:
        asinh_y = 0.0
    if (X * X + Y * Y) > 0:
        asinh_z = (Z / 2.0) * (X * X - Y * Y) * math.asinh(Z / math.sqrt(X * X + Y * Y))
    else:
        asinh_z = 0.0
    R_term = -R * (2.0 * X * X - Y * Y - Z * Z) / 6.0
    return atan_term + asinh_y + asinh_z + R_term


def _G2(X, Y, Z):
    """Antiderivative of ``g2 = -asinh(Z/sqrt(X^2+Y^2))`` -- off-diagonal corner primitive.

    Used for ``<N_xy>_T`` 64-corner sum. For ``<N_xz>_T`` call ``_G2(X, Z, Y)``;
    for ``<N_yz>_T`` call ``_G2(Y, Z, X)``.
    """
    R = math.sqrt(X * X + Y * Y + Z * Z)
    if (X * X + Y * Y) > 0:
        t1 = -X * Y * Z * math.asinh(Z / math.sqrt(X * X + Y * Y))
    else:
        t1 = 0.0
    if (Y * Y + Z * Z) > 0:
        t2 = (Y * (Y * Y - 3.0 * Z * Z) / 6.0) * math.asinh(X / math.sqrt(Y * Y + Z * Z))
    else:
        t2 = 0.0
    if (X * X + Z * Z) > 0:
        t3 = (X * (X * X - 3.0 * Z * Z) / 6.0) * math.asinh(Y / math.sqrt(X * X + Z * Z))
    else:
        t3 = 0.0
    if abs(X) < _EPS:
        atan_yz = 0.0
    else:
        atan_yz = (X * X * Z / 2.0) * math.atan(Y * Z / (X * R))
    if abs(Y) < _EPS:
        atan_xz = 0.0
    else:
        atan_xz = (Y * Y * Z / 2.0) * math.atan(X * Z / (Y * R))
    if abs(Z) < _EPS:
        atan_xy = 0.0
    else:
        atan_xy = (Z * Z * Z / 6.0) * math.atan(X * Y / (Z * R))
    R_term = X * Y * R / 3.0
    return t1 + t2 + t3 + atan_yz + atan_xz + atan_xy + R_term


def _G_ij(i, j, X, Y, Z):
    """Average-demag corner primitive for tensor element (i, j) at displacement (X, Y, Z)."""
    if i == j:
        if i == 0:
            return _G1(X, Y, Z)
        if i == 1:
            return _G1(Y, X, Z)
        return _G1(Z, X, Y)
    if (i, j) in ((0, 1), (1, 0)):
        return _G2(X, Y, Z)
    if (i, j) in ((0, 2), (2, 0)):
        return _G2(X, Z, Y)
    return _G2(Y, Z, X)


def _corner_signs():
    """Signed corners of an axis-aligned cuboid for the 8-term inclusion-exclusion sum."""
    out = []
    for ix, iy, iz in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0),
                       (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)):
        out.append((ix, iy, iz, (-1) ** (ix + iy + iz)))
    return out


try:
    from radia._radia_pybind import (
        _average_B_in_box as _cpp_average_B_in_box,
        _average_demag_tensor as _cpp_average_demag_tensor,
    )
    _HAVE_CPP = True
except ImportError:  # pragma: no cover - C++ extension absent in pure-Python builds
    _HAVE_CPP = False


def _average_demag_tensor_python(source_min, source_max, target_min, target_max):
    """Pure-Python reference implementation of the 64-corner sum.

    Used as the agreement-check baseline against the C++ kernel and as
    a fallback when the C++ extension is unavailable.
    """
    smin = np.asarray(source_min, dtype=float)
    smax = np.asarray(source_max, dtype=float)
    tmin = np.asarray(target_min, dtype=float)
    tmax = np.asarray(target_max, dtype=float)
    V_T = float(np.prod(tmax - tmin))

    src_extents = (
        (smin[0], smax[0]),
        (smin[1], smax[1]),
        (smin[2], smax[2]),
    )
    tgt_extents = (
        (tmin[0], tmax[0]),
        (tmin[1], tmax[1]),
        (tmin[2], tmax[2]),
    )
    corners = _corner_signs()
    src_corners = [
        (src_extents[0][ix], src_extents[1][iy], src_extents[2][iz], sign)
        for (ix, iy, iz, sign) in corners
    ]
    tgt_corners = [
        (tgt_extents[0][ix], tgt_extents[1][iy], tgt_extents[2][iz], sign)
        for (ix, iy, iz, sign) in corners
    ]

    A = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            total = 0.0
            for sx, sy, sz, ss in src_corners:
                for tx, ty, tz, ts in tgt_corners:
                    X = tx - sx
                    Y = ty - sy
                    Z = tz - sz
                    total += ss * ts * _G_ij(i, j, X, Y, Z)
            A[i, j] = total / (4.0 * math.pi * V_T)
    return A


def _overlap_volume(source_min, source_max, target_min, target_max):
    """Axis-aligned-cuboid intersection volume."""
    smin = np.asarray(source_min, dtype=float)
    smax = np.asarray(source_max, dtype=float)
    tmin = np.asarray(target_min, dtype=float)
    tmax = np.asarray(target_max, dtype=float)
    overlap = np.maximum(0.0, np.minimum(smax, tmax) - np.maximum(smin, tmin))
    return float(np.prod(overlap))


def average_demag_tensor(source_min, source_max, target_min, target_max):
    """3x3 average demag tensor ``A`` such that ``<B>_T = mu_0 (A @ M + M * V_overlap/V_T)``.

    Closed-form 64-corner sum of the Phase-alpha antiderivatives ``G1``, ``G2``.
    Calls the C++ kernel when available (typically 100x+ faster than the
    Python reference path).

    Parameters
    ----------
    source_min, source_max : (3,) array_like
        ``(x, y, z)`` lower / upper corners of the source box [m].
    target_min, target_max : (3,) array_like
        Lower / upper corners of the target box [m].

    Returns
    -------
    A : ndarray, shape (3, 3)
        Average demag tensor (dimensionless).
    """
    if _HAVE_CPP:
        smin = np.ascontiguousarray(source_min, dtype=float)
        smax = np.ascontiguousarray(source_max, dtype=float)
        tmin = np.ascontiguousarray(target_min, dtype=float)
        tmax = np.ascontiguousarray(target_max, dtype=float)
        return _cpp_average_demag_tensor(smin, smax, tmin, tmax)
    return _average_demag_tensor_python(source_min, source_max, target_min, target_max)


def average_B_in_box(
    M,
    source_min,
    source_max,
    target_min,
    target_max,
    *,
    method: str = "closed",
    n_quad: int = 8,
):
    """Average ``B`` over a target box from uniform ``M`` in a source box.

    Parameters
    ----------
    M : (3,) array_like
        Magnetisation in the source box [A/m] (Radia convention;
        ``J = mu_0 M`` Tesla).
    source_min, source_max : (3,) array_like
        ``(x, y, z)`` lower / upper corners of the source box [m].
    target_min, target_max : (3,) array_like
        Lower / upper corners of the target box [m].
    method : {"closed", "numerical"}
        ``"closed"`` (default) uses the Phase-alpha 64-corner-sum closed
        form (no quadrature error; calls C++ kernel when available).
        ``"numerical"`` uses the legacy tensor-product Gauss-Legendre
        quadrature kept for cross-checks.
    n_quad : int
        Gauss-Legendre points per axis -- only used when
        ``method == "numerical"``.

    Returns
    -------
    B_avg : ndarray, shape (3,)
        Average ``(<B_x>, <B_y>, <B_z>)`` over the target box [Tesla].
    """
    if method == "closed":
        if _HAVE_CPP:
            M_arr = np.ascontiguousarray(M, dtype=float)
            smin = np.ascontiguousarray(source_min, dtype=float)
            smax = np.ascontiguousarray(source_max, dtype=float)
            tmin = np.ascontiguousarray(target_min, dtype=float)
            tmax = np.ascontiguousarray(target_max, dtype=float)
            return _cpp_average_B_in_box(M_arr, smin, smax, tmin, tmax)
        A = average_demag_tensor(source_min, source_max, target_min, target_max)
        M_arr = np.asarray(M, dtype=float)
        V_T = float(np.prod(np.asarray(target_max) - np.asarray(target_min)))
        V_overlap = _overlap_volume(source_min, source_max, target_min, target_max)
        return MU_0 * (A @ M_arr) + MU_0 * M_arr * (V_overlap / V_T)
    if method == "numerical":
        return average_B_in_box_numerical(M, source_min, source_max, target_min, target_max, n_quad=n_quad)
    raise ValueError(f"unknown method {method!r}; expected 'closed' or 'numerical'")


def average_B_in_box_numerical(
    M,
    source_min,
    source_max,
    target_min,
    target_max,
    n_quad: int = 8,
):
    """Tensor-product Gauss-Legendre numerical integration of ``B(r)`` over target box.

    Kept for diagnostic / cross-validation use; production calls go through
    :func:`average_B_in_box` (closed form).
    """
    from radia.analytical_magnet import CuboidMagnet
    from .gauss_legendre import nodes_weights

    src_min = np.asarray(source_min, dtype=float)
    src_max = np.asarray(source_max, dtype=float)
    tgt_min = np.asarray(target_min, dtype=float)
    tgt_max = np.asarray(target_max, dtype=float)
    src_centre = 0.5 * (src_min + src_max)
    src_size = src_max - src_min
    tgt_centre = 0.5 * (tgt_min + tgt_max)
    tgt_size = tgt_max - tgt_min

    magnet = CuboidMagnet(
        center=src_centre.tolist(),
        dimensions=src_size.tolist(),
        magnetization=list(M),
        units="m",
    )

    x_n, w_x = nodes_weights(min(n_quad, 24))
    y_n, w_y = nodes_weights(min(n_quad, 24))
    z_n, w_z = nodes_weights(min(n_quad, 24))

    half = 0.5 * tgt_size
    xp = half[0] * x_n + tgt_centre[0]
    yp = half[1] * y_n + tgt_centre[1]
    zp = half[2] * z_n + tgt_centre[2]

    B_acc = np.zeros(3)
    W_acc = 0.0
    for i, x in enumerate(xp):
        for j, y in enumerate(yp):
            for k, z in enumerate(zp):
                w = w_x[i] * w_y[j] * w_z[k]
                B = np.array(magnet.get_B([x, y, z]))
                B_acc += w * B
                W_acc += w
    return B_acc / W_acc
