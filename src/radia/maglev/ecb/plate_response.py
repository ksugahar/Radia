"""plate_response.py -- alpha(s) frequency sweep + drag/lift crossover.

Helper utilities for ECB plate design: sweep alpha(s) over a frequency band,
identify the drag peak (where -Im[alpha/V] is max) and the lift crossover
(where Re[alpha/V] = 0.5).
"""
from __future__ import annotations

import math
import numpy as np

from ..mixed_galerkin.alpha import Y_mixed, alpha_from_Y


def sweep_alpha_response(lam, tau, g_n, K_SIBC, c1, V, sigma,
                           f_min=1e-1, f_max=1e8, n_pts=73):
    """Sweep alpha(s) over [f_min, f_max] (log-spaced).

    Returns
    -------
    f_arr : ndarray
    alpha_re_over_V : ndarray
    alpha_im_over_V : ndarray
    """
    f_arr = np.logspace(math.log10(f_min), math.log10(f_max), n_pts)
    alpha_re = np.zeros_like(f_arr)
    alpha_im = np.zeros_like(f_arr)
    for i, f in enumerate(f_arr):
        s = 1j * 2 * math.pi * f
        Y = Y_mixed(s, lam, tau, g_n, K_SIBC, c1)
        a = alpha_from_Y(Y, V, sigma)
        alpha_re[i] = a.real / V
        alpha_im[i] = a.imag / V
    return f_arr, alpha_re, alpha_im


def find_drag_peak(f_arr, alpha_im):
    """Find frequency where -Im[alpha/V] is maximum (drag dissipation peak).

    Returns (f_peak_Hz, neg_im_at_peak).  Beware: at very low f the bulk
    Foster truncation can produce artifacts -- check the response curve
    shape.
    """
    i = int(np.argmax(-alpha_im))
    return float(f_arr[i]), float(-alpha_im[i])


def find_lift_crossover(f_arr, alpha_re, threshold=0.5):
    """Find smallest f where Re[alpha/V] >= threshold (lift develops).

    Returns (f_crossover_Hz, None if not crossed).
    """
    for i, x in enumerate(alpha_re):
        if x >= threshold:
            return float(f_arr[i])
    return None
