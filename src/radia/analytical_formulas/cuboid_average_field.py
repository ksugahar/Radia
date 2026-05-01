"""Average magnetic flux density over a target rectangular box.

Part 6 §7 of the Wakao-Igarashi-Fujiwara-Kameari series. Originally
derived for micromagnetics: average ``B`` in one cubic lattice cell
due to the magnetisation of another lattice cell (same size).

Implementation status (2026-05-01)
----------------------------------
This module currently provides the **numerical-integration** path
only. Wrapping :func:`radia.analytical_magnet.CuboidMagnet.get_B`
with a tensor-product Gauss-Legendre rule over the target box gives
the average to spectral precision and is unambiguous.

The PDF Part 6 §7 closed-form (eq 53-56, with antiderivatives F1, F2)
was attempted but the OCR-extracted F2 formula could not be reconciled
with target-box-symmetry constraints (B_y, B_z components of the
average leak when only M_x is non-zero, in violation of the source/
target y- and z-mirror symmetries). The correct closed-form
references back to Newell, Williams & Dunlop, J. Geophys. Res. 98
(1993), pp. 9551-9555 -- a clean re-derivation against that paper is
deferred until a high-fidelity copy of Part 6 is available.

A C++ kernel (``src/core/rad_average_field.cpp``) is similarly
deferred: the Python numerical-integration path is fast enough for
all current Radia call sites (FEM-MMM coupling, micromagnetics
prototypes), and a closed-form C++ version would be valuable only if
the call site becomes an inner-loop bottleneck.

Useful for
----------
* Micromagnetics on a regular cubic lattice: cell-cell magnetic
  interaction.
* Mesh-to-mesh transfer between two rectangular grids.
* FEM-MMM coupling: cell-averaged B for a rectangular FEM element.

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


def average_B_in_box(
    M,
    source_min,
    source_max,
    target_min,
    target_max,
    n_quad: int = 8,
):
    """Average ``B`` over a target box from uniform ``M`` in a source box.

    Numerical-integration path: tensor-product Gauss-Legendre quadrature
    of the analytical point-field formula over the target box.

    Source and target boxes do **not** need to share dimensions;
    overlapping boxes are also supported (the point-field formula
    handles the source interior correctly).

    Parameters
    ----------
    M : (3,) array_like
        Magnetisation in the source box [A/m] (Radia convention;
        ``J = mu_0 M`` Tesla).
    source_min, source_max : (3,) array_like
        ``(x, y, z)`` lower / upper corners of the source box [m].
    target_min, target_max : (3,) array_like
        Lower / upper corners of the target box [m].
    n_quad : int
        Gauss-Legendre points per axis for the area integral
        (default 8; gives ~1e-9 relative error for boxes well clear
        of each other; increase for nearby / overlapping geometries).

    Returns
    -------
    B_avg : ndarray, shape (3,)
        Average ``(<B_x>, <B_y>, <B_z>)`` over the target box [Tesla].
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

    # CuboidMagnet uses mm coordinates by default; pass meters explicitly.
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

    # Sum: B_sum = sum_ijk w_x[i] w_y[j] w_z[k] * B_at(xp[i], yp[j], zp[k]).
    # Volume-weighted normalisation cancels the half**3 prefactor against
    # the V_target factor.
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
