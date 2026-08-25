"""Proximity-aware iterative PEEC for perimeter filaments.

Approach
========
The PEEC pipeline uses a shape-specific isolated-conductor model for
``Zs_fil`` (round-wire Bessel or rectangular Dowell).  That captures
self-skin but misses proximity (eddy currents driven by neighbour-turn
H field).

This module wraps ``solve_loop_bundle`` with an outer iteration:

1. Solve bundle with current ``Zs_fil`` → ``I_fil``
2. From ``I_fil``, evaluate the total tangential H at each
   filament's surface position (Biot-Savart from all segments).
3. The Leontovich SIBC dissipation at filament k's local arc s_k
   is ``½ Re(Z_s) × |H_t_k|² × s_k × L_k``.  Set
   ``Re(Zs_fil_k) = (½ Re(Z_s) × |H_t_k|² × s_k × L_k) / |I_k|²``
   so the bundle re-solve dissipates the right amount per filament.
4. Repeat to convergence.

Implementation notes (debugged 2026-05-20)
==========================================
* Biot-Savart angular sign convention: ``l·r1/|r1| - l·r2/|r2|``
  (this is the cos α₁ − cos α₂ form, NOT the swapped one).  Confirmed
  against right-hand rule on a 2-m wire at perp 1 m.
* ``r_min`` clamp must be **very small** (1e-12).  An earlier draft
  used 0.1 × wire_radius which over-clamped close-range contributions
  and made the wire's own surface |H| 360× too small (0.14 vs 50 A/m
  Ampere's-law expectation).
* Self-exclusion: only the eval segment's OWN contribution is zeroed
  (avoids the Biot-Savart singularity at the segment midpoint).
  Adjacent segments of the same filament ARE included — they
  contribute the bulk of |H_t_self| at the eval surface point.
* For divergence prevention: under-relax ``Zs_fil`` updates with
  ``relax`` (default 0.3) and clamp tiny |I_k| to a floor of
  ``I_port / (10 × K)`` so the ``Zs_fil ∝ 1/|I_k|²`` update doesn't
  blow up when a filament's current shrinks.

Validation on the 3-turn pancake (Cu, 150 kHz, n_peri=16):
* Cylinder_ac_impedance self-skin asymptote: 3.97 mΩ.
* Bundle with isolated-wire Bessel ``Zs_fil``: ``R = 3.68 mΩ``.
* This module with proximity iteration: see ``R_coil`` field.
* Empirical target: ~15 mΩ on the kubota 3turncoil LCR hi-tester.
"""
from __future__ import annotations

import math
from typing import Optional

import numpy as np


# ----------------------------------------------------------------------
# Biot-Savart kernel (validated against 2-wire analytical 2026-05-20)
# ----------------------------------------------------------------------

def _biot_savart_H(eval_pts: np.ndarray,
                     starts: np.ndarray, ends: np.ndarray,
                     currents: np.ndarray,
                     self_index: Optional[np.ndarray] = None,
                     r_min: float = 1.0e-12) -> np.ndarray:
    """Vectorised H field at ``eval_pts`` from straight current segments.

    Args:
        eval_pts: (M, 3) where H is to be evaluated.
        starts, ends: (N, 3) source segment endpoints.
        currents: (N,) complex current along each segment.
        self_index: optional (M,) integer mapping eval m → segment to
            EXCLUDE.  Use -1 (or any negative) for no exclusion.
        r_min: regularisation floor on |r| (m).  Must be SMALL.

    Returns:
        H: (M, 3) complex H field.
    """
    M = eval_pts.shape[0]
    N = starts.shape[0]
    r1 = eval_pts[:, None, :] - starts[None, :, :]           # (M, N, 3)
    r2 = eval_pts[:, None, :] - ends[None, :, :]
    r1_mag = np.linalg.norm(r1, axis=2)                       # (M, N)
    r2_mag = np.linalg.norm(r2, axis=2)
    r1_mag = np.where(r1_mag > r_min, r1_mag, r_min)
    r2_mag = np.where(r2_mag > r_min, r2_mag, r_min)
    seg_vec = ends - starts
    l = seg_vec[None, :, :]
    l_cross_r1 = np.cross(l, r1, axis=2)                      # (M, N, 3)
    l_cross_r1_mag2 = np.sum(l_cross_r1 * l_cross_r1, axis=2)
    l_cross_r1_mag2 = np.where(l_cross_r1_mag2 > r_min * r_min,
                                  l_cross_r1_mag2, r_min * r_min)
    l_dot_r1 = np.sum(l * r1, axis=2)
    l_dot_r2 = np.sum(l * r2, axis=2)
    # Correct sign convention (cos α₁ - cos α₂):
    angular = l_dot_r1 / r1_mag - l_dot_r2 / r2_mag
    factor = currents[None, :] * angular / l_cross_r1_mag2

    if self_index is not None:
        rows = np.arange(M)
        valid = self_index >= 0
        if valid.any():
            factor[rows[valid], self_index[valid]] = 0.0

    H = (1.0 / (4.0 * math.pi)) * (factor[:, :, None] * l_cross_r1).sum(axis=1)
    return H


def _flatten_paths(filament_paths):
    starts, ends, fil_of_seg = [], [], []
    for k, path_k in enumerate(filament_paths):
        for (p1, p2) in path_k:
            starts.append(np.asarray(p1, dtype=float))
            ends.append(np.asarray(p2, dtype=float))
            fil_of_seg.append(k)
    return (np.asarray(starts), np.asarray(ends),
            np.asarray(fil_of_seg, dtype=int))


def solve_proximity_iterative(
    R_f: np.ndarray,
    L_f: np.ndarray,
    filament_paths,
    frequency: float,
    sigma: float,
    wire_radius_m: float,
    n_peri: int,
    I_port: complex = 1.0,
    max_iter: int = 30,
    tol: float = 1.0e-3,
    relax: float = 0.3,
    verbose: bool = False,
    self_impedance_per_m: complex | None = None,
    dc_resistance_per_m: float | None = None,
    perimeter_m: float | None = None,
    internal_impedance_model: str | None = None,
) -> dict:
    """Outer iteration capturing proximity in perim PEEC bundle.

    Returns dict with R_coil, L_coil, history etc.  See module
    docstring for the algorithm and validation notes.
    """
    from radia.peec_bundle import solve_loop_bundle

    omega = 2.0 * math.pi * frequency
    if omega <= 0:
        I_fil, V_port = solve_loop_bundle(R_f, L_f, frequency, I_port=I_port)
        return {
            "I_fil": I_fil, "V_port": V_port,
            "L_coil": 0.0,
            "R_coil": float(V_port.real / abs(I_port) ** 2) if I_port != 0 else 0.0,
            "Zs_fil": np.zeros(R_f.shape[0], dtype=complex),
            "n_iter": 0, "converged": True, "history": [],
        }

    if n_peri <= 0:
        raise ValueError(f"n_peri must be > 0, got {n_peri}")
    K = R_f.shape[0]
    if len(filament_paths) != K:
        raise ValueError(
            "filament_paths must contain one path per PEEC filament: "
            f"got {len(filament_paths)} paths for {K} matrix rows")
    if K % n_peri != 0:
        raise ValueError(
            "the total PEEC filament count must be an integer multiple "
            "of n_peri: "
            f"got n_peri={n_peri}, K={K}")
    MU_0 = 4.0e-7 * math.pi
    delta = math.sqrt(2.0 / (omega * MU_0 * sigma))
    Re_Zs = 1.0 / (sigma * delta)                       # surface resistance [Ω/sq]
    # Build flattened segment arrays
    starts, ends, fil_of_seg = _flatten_paths(filament_paths)

    # Per-filament axial length
    L_fil = np.array([
        sum(np.linalg.norm(np.asarray(p2) - np.asarray(p1))
            for (p1, p2) in path_k)
        for path_k in filament_paths
    ])
    # Per-filament arc length on the actual conductor perimeter.  Round
    # wire keeps the historical 2*pi*a path; rectangular CAD supplies
    # 2*(w+t), so no equivalent-circle perimeter enters proximity loss.
    if perimeter_m is None:
        if wire_radius_m <= 0.0:
            raise ValueError(
                "perimeter_m is required when wire_radius_m is unavailable")
        perimeter = 2.0 * math.pi * wire_radius_m
    else:
        perimeter = float(perimeter_m)
    if perimeter <= 0.0:
        raise ValueError(f"perimeter_m must be > 0, got {perimeter}")
    s_k = perimeter / float(n_peri)

    # Eval points: one per filament = the filament's geometric midpoint
    # (segment midpoints averaged).  For each k, the eval m=k.
    # Compute per-filament midpoint:
    mid_per_fil = np.zeros((K, 3), dtype=float)
    seg_mid_per_fil_count = np.zeros(K, dtype=int)
    for n in range(starts.shape[0]):
        k = fil_of_seg[n]
        mid_per_fil[k] += 0.5 * (starts[n] + ends[n])
        seg_mid_per_fil_count[k] += 1
    mid_per_fil = mid_per_fil / seg_mid_per_fil_count[:, None]

    # Initialize from the selected isolated-conductor model.  Explicit
    # values carry the rectangular Dowell path; omitted values preserve
    # the public round-wire Bessel API.
    if self_impedance_per_m is None or dc_resistance_per_m is None:
        if wire_radius_m <= 0.0:
            raise ValueError(
                "wire_radius_m is required for the round-wire Bessel model")
        from radia.analytical_formulas.conductor_impedance import (
            cylinder_ac_impedance, cylinder_dc_resistance,
        )
        Z_ac = cylinder_ac_impedance(wire_radius_m, sigma, omega)
        R_dc_per_m = cylinder_dc_resistance(wire_radius_m, sigma)
        resolved_model = internal_impedance_model or "round-bessel"
    else:
        Z_ac = complex(self_impedance_per_m)
        R_dc_per_m = float(dc_resistance_per_m)
        if R_dc_per_m <= 0.0:
            raise ValueError(
                f"dc_resistance_per_m must be > 0, got {R_dc_per_m}")
        resolved_model = internal_impedance_model or "explicit"
    dZ_per_m = complex(Z_ac) - complex(R_dc_per_m, 0.0)
    self_delta = n_peri * dZ_per_m * L_fil
    Zs_fil = self_delta.astype(complex)

    history = []
    converged = False
    I_floor = abs(I_port) / (10.0 * K)   # protect against |I_k| → 0

    for it in range(max_iter):
        I_fil, V_port = solve_loop_bundle(
            R_f, L_f, frequency, I_port=I_port, Zs_fil=Zs_fil)
        R_now = float(V_port.real / abs(I_port) ** 2)

        # Compute H_t at each filament's "representative" eval point.
        # For proximity-aware Zs, we want |H_t_total|² INCLUDING all
        # contributions (self + neighbour).  The Bessel "self" part is
        # baked into the initial Zs_fil (cylinder_ac_impedance); to
        # avoid double-counting, we evaluate only the DIFFERENCE
        # between the bundle's actual surface H and the isolated-wire
        # symmetric value.  Simpler: just compute total |H_t|² and
        # SET R from it (replace, not add).  See solve_proximity_total.
        #
        # For now use the additive iteration: evaluate H from cross-
        # filament contributions only (this is the proximity term).
        H_cross = _biot_savart_H(mid_per_fil, starts, ends,
                                  I_fil[fil_of_seg])
        # Subtract own-filament's contribution (self-skin already in
        # the Zs_fil_self bake-in).  Build a mask per eval point.
        H_self = np.zeros((K, 3), dtype=complex)
        for k in range(K):
            sel = (fil_of_seg == k)
            H_k = _biot_savart_H(
                mid_per_fil[k:k+1], starts[sel], ends[sel],
                I_fil[fil_of_seg][sel])
            H_self[k] = H_k[0]
        H_t_proximity = H_cross - H_self    # neighbour contribution only
        H_mag2_prox = np.real(
            np.sum(H_t_proximity * np.conj(H_t_proximity), axis=1)
        )

        # Per-filament proximity dissipation (Leontovich):
        # P_k_proximity = ½ Re(Z_s) × |H_t_prox|² × s_k × L_k
        # Set ΔZs_fil_k so |I_k|² × Re(ΔZs_fil_k) = P_k_proximity:
        I_mag2 = np.real(I_fil * np.conj(I_fil))
        I_mag2 = np.maximum(I_mag2, I_floor * I_floor)    # floor
        prox_real = 0.5 * Re_Zs * H_mag2_prox * s_k * L_fil / I_mag2
        # The reactive component scales similarly (Re(Z_s) = Im(Z_s)
        # for Leontovich, so use complex (1+j)/(σδ)):
        prox_dZ = prox_real * (1.0 + 1.0j)

        # New target Zs_fil = initial self + proximity correction
        Zs_new = self_delta + prox_dZ
        # Under-relax
        diff_norm = np.linalg.norm(Zs_new - Zs_fil) / max(
            np.linalg.norm(Zs_fil), 1e-30)
        history.append((it, R_now, float(diff_norm),
                         float(np.linalg.norm(H_mag2_prox.max() ** 0.5))))
        if verbose:
            print(f"  iter {it:2d}: R={R_now*1e3:.4f} mΩ, "
                  f"||ΔZ||/||Z||={diff_norm:.4f}, "
                  f"max|H_prox|={H_mag2_prox.max()**0.5:.2f} A/m")
        if diff_norm < tol:
            Zs_fil = Zs_new
            converged = True
            break
        Zs_fil = (1.0 - relax) * Zs_fil + relax * Zs_new

    # Final consistent solve
    I_fil, V_port = solve_loop_bundle(
        R_f, L_f, frequency, I_port=I_port, Zs_fil=Zs_fil)
    Z_coil = V_port / I_port
    L_coil = float(Z_coil.imag / omega)
    R_coil = float(Z_coil.real)

    return {
        "I_fil": I_fil, "V_port": V_port,
        "L_coil": L_coil, "R_coil": R_coil,
        "Zs_fil": Zs_fil,
        "n_iter": len(history),
        "converged": converged,
        "history": history,
        "internal_impedance_model": resolved_model,
        "perimeter_m": perimeter,
    }
