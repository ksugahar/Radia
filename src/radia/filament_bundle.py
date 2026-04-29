"""Radia filament bundle with COMPLEX per-filament currents.

PEEC AC analysis produces a bundle of K filaments (closed polygonal
loops of N straight segments each) with COMPLEX per-filament currents.
This module provides a Radia-style API to query SI-correct AC fields:

    >>> from radia.filament_bundle import FilamentBundleAC
    >>> bundle = FilamentBundleAC(filament_paths, currents)
    >>> H_complex = bundle.fld('h', obs_points)   # (N, 3) complex
    >>> B_complex = bundle.fld('b', obs_points)   # (N, 3) complex (Tesla)
    >>> phi_complex = bundle.fld('phi', obs_points)  # (N,) complex

Internally backed by the C++ Biot-Savart kernel
``radia._radia_pybind._HFromSegmentsComplex`` (Phase 1.12) and the
SIBC path-integral phi via ``compute_phi_inc_from_filaments`` (also
C++-accelerated since Phase 1.12).

Note: as of 2026-04-30 the core ``radTFlmLinCur`` (used by
``rad.ObjFlmCur`` + ``rad.Fld``) is also SI-correct (legacy
``ConI = 1.E-04 * I`` was replaced with proper ``INV_FOUR_PI``
for B/H and ``MU_0_OVER_FOUR_PI`` for A in
``src/core/rad_filament.cpp``).  This module continues to use
``_HFromSegmentsComplex`` directly because it is faster for large
AC bundles (one batched C++ call vs many ``rad.ObjFlmCur`` element
allocations + ``rad.Fld`` dispatches).
"""
from __future__ import annotations
import numpy as np

MU_0 = 4.0 * np.pi * 1e-7    # T·m/A


def _cpp_h_segments_complex():
    try:
        from radia import _radia_pybind as _rpb
        return _rpb._HFromSegmentsComplex
    except (ImportError, AttributeError):
        return None


def _cpp_a_segments_complex():
    try:
        from radia import _radia_pybind as _rpb
        return _rpb._AFromSegmentsComplex
    except (ImportError, AttributeError):
        return None


def _flatten_filament_paths(filament_paths, currents):
    """Convert (filament_paths, currents) into (segs_array, I_re_array,
    I_im_array) flat numpy arrays suitable for the C++ kernel.
    """
    flat_segs = []
    flat_I = []
    for fil_segs, Ik in zip(filament_paths, currents):
        Ik_c = complex(Ik)
        for (p1, p2) in fil_segs:
            flat_segs.append([list(p1), list(p2)])
            flat_I.append(Ik_c)
    if not flat_segs:
        return None, None, None
    segs = np.ascontiguousarray(flat_segs, dtype=np.float64)
    I_arr = np.asarray(flat_I, dtype=complex)
    I_re = np.ascontiguousarray(I_arr.real)
    I_im = np.ascontiguousarray(I_arr.imag)
    return segs, I_re, I_im


class FilamentBundleAC:
    """A bundle of straight-segment filaments with complex per-filament
    currents.  SI-correct field queries via our C++ Biot-Savart kernel.

    Args:
        filament_paths: list of K filament segment lists.  Each filament
            segment list is [(p1, p2), (p2, p3), ...] tracing a polygonal
            current path (closed loop for PEEC AC bundle).
        currents: length-K array of per-filament complex currents [A].
    """

    def __init__(self, filament_paths, currents):
        currents = np.asarray(currents, dtype=complex)
        if len(filament_paths) != len(currents):
            raise ValueError(
                f"filament_paths has {len(filament_paths)} entries but "
                f"currents has {len(currents)} -- they must match.")
        self._segs, self._I_re, self._I_im = _flatten_filament_paths(
            filament_paths, currents)
        if self._segs is None:
            raise ValueError("FilamentBundleAC: no segments in any filament")
        self._n_filaments = len(filament_paths)
        self._n_segments = len(self._segs)
        # Cache the input for path-integral phi
        self._filament_paths = filament_paths
        self._currents = currents

    def fld(self, kind, obs_points, n_threads=0):
        """Evaluate the requested field at observation points, returning
        complex result (SI units).

        Args:
            kind: one of 'h' (A/m), 'b' (Tesla), 'a' (T*m), 'phi' (A).
                  Vector fields return shape (N_obs, 3); scalar 'phi'
                  returns shape (N_obs,).  Per Radia convention.
            obs_points: (N_obs, 3) observation points [m].
            n_threads: 0 = TaskManager default; > 0 = request that many.

        Returns:
            complex ndarray, dtype=complex128.
        """
        kind_low = kind.lower()
        obs = np.asarray(obs_points, dtype=float)
        if obs.ndim == 1:
            obs = obs.reshape(1, 3)
        obs_c = np.ascontiguousarray(obs)

        if kind_low == 'h':
            return self._compute_H(obs_c, n_threads)
        if kind_low == 'b':
            H = self._compute_H(obs_c, n_threads)
            return MU_0 * H
        if kind_low == 'a':
            return self._compute_A(obs_c, n_threads)
        if kind_low == 'phi':
            return self._compute_phi(obs_c)
        raise ValueError(f"unknown field kind {kind!r}, expected one of "
                          f"'h', 'b', 'a', 'phi'")

    def _compute_H(self, obs, n_threads):
        cpp = _cpp_h_segments_complex()
        if cpp is None:
            raise RuntimeError(
                "C++ _HFromSegmentsComplex not available; rebuild "
                "radia with Phase 1.12 (commit b1213a7b).")
        H_re, H_im = cpp(self._segs, obs, self._I_re, self._I_im, n_threads)
        return H_re + 1j * H_im

    def _compute_A(self, obs, n_threads):
        cpp = _cpp_a_segments_complex()
        if cpp is None:
            raise RuntimeError(
                "C++ _AFromSegmentsComplex not available; rebuild "
                "radia with the 2026-04-30 vector-potential kernel.")
        A_re, A_im = cpp(self._segs, obs, self._I_re, self._I_im, n_threads)
        return A_re + 1j * A_im

    def _compute_phi(self, obs):
        """phi from path integral phi(P) = -int H.dl from far to P.
        Uses the Phase 1.12 C++-accelerated compute_phi_inc_from_filaments
        helper (same gauge: phi -> 0 at z_far).
        """
        # Local import to keep this module Radia-only
        try:
            import sys
            sys.path.insert(0,
                "S:/Radia/01_GitHub/examples/induction_heating/bem_reference")
            from bem_sibc_solver import compute_phi_inc_from_filaments
        except ImportError as e:
            raise RuntimeError(
                "compute_phi_inc_from_filaments unavailable: " + str(e))
        return compute_phi_inc_from_filaments(
            obs, self._filament_paths, self._currents)

    def __len__(self):
        return self._n_filaments

    def __repr__(self):
        return (f"FilamentBundleAC(n_filaments={self._n_filaments}, "
                f"n_segments={self._n_segments})")
