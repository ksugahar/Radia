"""Workpiece surface bundle for SIBC-induced fields and Telegen reciprocity ΔL.

The SIBC scalar BIE solve (Phase 1.7+) produces the magnetic scalar potential
phi at workpiece surface vertices.  The induced surface current density is

    H_t = -grad_s(phi)
    J_s = n x H_t = -n x grad_s(phi)

For P1 H1 phi on a flat triangle, grad_s(phi) is constant per triangle, hence
J_s is constant per triangle.

This module provides:

* :func:`phi_to_Js` -- convert per-vertex complex phi to per-triangle complex J_s.
* :class:`WorkpieceSurfaceBundle` -- Radia-style ``.fld('b'|'a', obs)`` access.
* :func:`delta_L_telegen` -- ΔL from coil reciprocity (uses PEEC A_inc on
  workpiece centroids; returns complex ΔL [Henry]).

All field queries use the C++ kernels in ``rad_biot_savart_surface.{h,cpp}``
for performance (TaskManager parallel over obs points).
"""
from __future__ import annotations
import numpy as np

MU_0 = 4.0 * np.pi * 1e-7


def _cpp_b_triangles_complex():
    try:
        from radia import _radia_pybind as _rpb
        return _rpb._BFromTrianglesComplex
    except (ImportError, AttributeError):
        return None


def _cpp_a_triangles_complex():
    try:
        from radia import _radia_pybind as _rpb
        return _rpb._AFromTrianglesComplex
    except (ImportError, AttributeError):
        return None


def _cpp_a_segments_complex():
    try:
        from radia import _radia_pybind as _rpb
        return _rpb._AFromSegmentsComplex
    except (ImportError, AttributeError):
        return None


def _cpp_h_segments_complex():
    try:
        from radia import _radia_pybind as _rpb
        return _rpb._HFromSegmentsComplex
    except (ImportError, AttributeError):
        return None


def phi_to_Js(vertices, triangles, phi_v, normals=None):
    """Convert per-vertex complex phi to per-triangle complex J_s.

    For a P1 hat-function phi on a flat triangle T with vertices v1, v2, v3
    and outward normal n_T, the surface gradient is constant:

        grad_s phi = (1 / (2 A)) * sum_i phi_i * (n_T x e_i)

    where e_i = v_{i+2} - v_{i+1} (cyclic, opposite vertex i).  Then

        J_s = -n x grad_s phi

    Args:
        vertices: (N_v, 3) float vertex coordinates [m]
        triangles: (N_t, 3) int vertex indices (0-based, CCW outward)
        phi_v: (N_v,) complex per-vertex magnetic scalar potential [A]
        normals: optional (N_t, 3) float outward unit normals.  If None,
            computed from CCW vertex ordering.

    Returns:
        centroids: (N_t, 3) float triangle centroids [m]
        areas:     (N_t,)   float triangle areas [m^2]
        J_s:       (N_t, 3) complex per-triangle surface current density [A/m]
    """
    V = np.asarray(vertices, dtype=float)
    T = np.asarray(triangles, dtype=int)
    phi = np.asarray(phi_v, dtype=complex)

    v1 = V[T[:, 0]]
    v2 = V[T[:, 1]]
    v3 = V[T[:, 2]]
    e12 = v2 - v1
    e13 = v3 - v1
    cr = np.cross(e12, e13)             # (N_t, 3)
    cr_mag = np.linalg.norm(cr, axis=1)
    valid = cr_mag > 0
    areas = 0.5 * cr_mag

    if normals is None:
        n = np.zeros_like(cr)
        n[valid] = cr[valid] / cr_mag[valid, None]
    else:
        n = np.asarray(normals, dtype=float)

    centroids = (v1 + v2 + v3) / 3.0

    # Edges opposite each vertex (cyclic):
    #   e_opp[0] = v3 - v2  (opposite v1)
    #   e_opp[1] = v1 - v3  (opposite v2)
    #   e_opp[2] = v2 - v1  (opposite v3)
    e_opp = np.stack([v3 - v2, v1 - v3, v2 - v1], axis=1)   # (N_t, 3, 3)

    phi1 = phi[T[:, 0]][:, None]         # (N_t, 1)
    phi2 = phi[T[:, 1]][:, None]
    phi3 = phi[T[:, 2]][:, None]

    # n x e_opp_i for each triangle
    n_x_e0 = np.cross(n, e_opp[:, 0])    # (N_t, 3)
    n_x_e1 = np.cross(n, e_opp[:, 1])
    n_x_e2 = np.cross(n, e_opp[:, 2])

    # grad_s phi = (1/(2 A)) * sum_i phi_i * (n x e_opp_i)
    inv_2A = np.zeros(len(T), dtype=float)
    inv_2A[valid] = 1.0 / (2.0 * areas[valid])
    grad_s_phi = inv_2A[:, None] * (
        phi1 * n_x_e0 + phi2 * n_x_e1 + phi3 * n_x_e2
    )                                     # (N_t, 3) complex

    # J_s = -n x grad_s phi
    J_s = -np.cross(n, grad_s_phi)
    return centroids, areas, J_s


class WorkpieceSurfaceBundle:
    """Radia-style bundle for the SIBC-induced surface current on a
    triangulated workpiece surface.

    Args:
        vertices:  (N_v, 3) float vertex coordinates [m]
        triangles: (N_t, 3) int vertex indices (0-based)
        J_s:       (N_t, 3) complex per-triangle surface current density [A/m]
                   (already computed; e.g., via :func:`phi_to_Js`)

    Use :meth:`fld` to evaluate B [T] or A [T*m] at observation points.
    """

    def __init__(self, vertices, triangles, J_s):
        V = np.asarray(vertices, dtype=float)
        T = np.asarray(triangles, dtype=int)
        J = np.asarray(J_s, dtype=complex)
        if T.ndim != 2 or T.shape[1] != 3:
            raise ValueError("triangles must be (N_t, 3)")
        if J.shape != T.shape:
            raise ValueError(
                f"J_s must be (N_t, 3) matching triangles; got {J.shape} "
                f"vs {T.shape}.")
        # Materialize the (N_t, 3, 3) per-triangle vertex array for the
        # C++ kernel (it expects vertices laid out per-triangle, not
        # indexed via a separate vertex table).
        self._verts_tri = np.ascontiguousarray(V[T], dtype=float)   # (N_t, 3, 3)
        self._J_re = np.ascontiguousarray(J.real, dtype=float)
        self._J_im = np.ascontiguousarray(J.imag, dtype=float)
        self._n_t = len(T)

    def __len__(self):
        return self._n_t

    def __repr__(self):
        return f"WorkpieceSurfaceBundle(n_triangles={self._n_t})"

    def fld(self, kind, obs, n_threads=0):
        """Evaluate field at observation points.

        Args:
            kind: 'b' (Tesla) or 'a' (T*m).  'h' = b/mu_0 also accepted.
            obs:  (N, 3) observation points [m].
            n_threads: 0 = TaskManager default.

        Returns:
            complex ndarray, shape (N, 3).
        """
        kind_low = kind.lower()
        obs_arr = np.asarray(obs, dtype=float)
        if obs_arr.ndim == 1:
            obs_arr = obs_arr.reshape(1, 3)
        obs_c = np.ascontiguousarray(obs_arr)

        if kind_low == 'b':
            cpp = _cpp_b_triangles_complex()
            if cpp is None:
                raise RuntimeError(
                    "C++ _BFromTrianglesComplex unavailable; rebuild radia.")
            B_re, B_im = cpp(self._verts_tri, self._J_re, self._J_im,
                              obs_c, n_threads)
            return B_re + 1j * B_im
        if kind_low == 'a':
            cpp = _cpp_a_triangles_complex()
            if cpp is None:
                raise RuntimeError(
                    "C++ _AFromTrianglesComplex unavailable; rebuild radia.")
            A_re, A_im = cpp(self._verts_tri, self._J_re, self._J_im,
                              obs_c, n_threads)
            return A_re + 1j * A_im
        if kind_low == 'h':
            B = self.fld('b', obs_c, n_threads)
            return B / MU_0
        raise ValueError(f"unknown field kind {kind!r}; expected 'b', 'a', 'h'")


def delta_L_telegen(filament_paths, currents,
                     centroids, areas, J_s,
                     I_port,
                     n_threads=0):
    """Compute Δ L_coil due to workpiece via Telegen reciprocity.

    The mutual flux change in the coil (driven with port current I_port) due
    to the workpiece's induced surface current J_s is, via Telegen reciprocity:

        ΔΦ_coil = (1/I_port) * integral_S J_s . A_inc dS

    where A_inc is the vector potential generated by the SAME coil filament
    bundle (currents = the per-filament currents I_fil obtained from the
    PEEC forward solve at port current I_port) evaluated AT the workpiece
    surface, and J_s is the SIBC-induced surface current.

    The change in coil port inductance is

        Δ L_port = ΔΦ_coil / I_port
                 = (1/I_port^2) * integral_S J_s . A_inc dS

    For a non-magnetic conductive workpiece (Lenz's law), Re(Δ L_port) < 0;
    Im(Δ L_port) corresponds to a resistance addition Δ R = -ω · Im(Δ L_port).

    Args:
        filament_paths: list of K filament polylines (each [(p1, p2), ...])
        currents: (K,) complex per-filament currents I_fil from the PEEC
                  forward solve at I_port [A].  These define A_inc on S.
        centroids: (N_t, 3) workpiece triangle centroids [m]
        areas: (N_t,) workpiece triangle areas [m^2]
        J_s: (N_t, 3) complex workpiece per-triangle J_s [A/m]
        I_port: scalar (real or complex) coil port current [A].  This is the
                normalisation that makes Δ L_port a port-level inductance.
                Use the SAME drive amplitude that produced ``currents`` and
                ``J_s`` (typically the user-supplied ``--current`` flag).
        n_threads: thread count for A_inc evaluation (0 = TaskManager default)

    Returns:
        delta_L: complex Δ L_port [Henry]; Re(ΔL) = inductance change,
                 Im(ΔL) → resistance addition via Δ R = -ω · Im(ΔL).

    Implementation: 1-point centroid quadrature for ∫ J_s · A_inc dS.  J_s
    is constant per triangle and A_inc varies smoothly on the wp scale, so
    centroid quadrature is sufficient (this is a single area integral, NOT
    a 4D BEM kernel).
    """
    cpp = _cpp_a_segments_complex()
    if cpp is None:
        raise RuntimeError(
            "C++ _AFromSegmentsComplex unavailable; rebuild radia.")

    # Flatten filament_paths to a single (N_seg, 2, 3) segment array with
    # broadcast complex per-segment currents.
    flat_segs = []
    flat_I = []
    currents_arr = np.asarray(currents, dtype=complex)
    for fil_segs, Ik in zip(filament_paths, currents_arr):
        Ik_c = complex(Ik)
        for (p1, p2) in fil_segs:
            flat_segs.append([list(p1), list(p2)])
            flat_I.append(Ik_c)
    if not flat_segs:
        return 0.0 + 0.0j
    segs = np.ascontiguousarray(flat_segs, dtype=float)
    I_arr = np.asarray(flat_I, dtype=complex)
    I_re = np.ascontiguousarray(I_arr.real)
    I_im = np.ascontiguousarray(I_arr.imag)

    obs = np.ascontiguousarray(centroids, dtype=float)
    A_re, A_im = cpp(segs, obs, I_re, I_im, n_threads)
    A_inc = A_re + 1j * A_im                          # (N_t, 3) complex

    J = np.asarray(J_s, dtype=complex)
    a = np.asarray(areas, dtype=float)

    # integral = sum_T  (J_s_T . A_inc_T) * area_T
    integ = np.sum(np.einsum('tc,tc->t', J, A_inc) * a)

    Ip = complex(I_port)
    if abs(Ip) < 1e-30:
        return 0.0 + 0.0j
    delta_L = integ / (Ip * Ip)
    return delta_L


def delta_L_telegen_phiB(filament_paths, currents,
                          centroids, areas, normals, phi_avg,
                          I_port, n_threads=0):
    """Gauge-invariant Δ L_port via the ∫_S φ (n · B_inc) dS form.

    Mathematically equivalent to :func:`delta_L_telegen` in the
    continuum (J_s · A_inc form), but stable at the discrete level.

    Continuum identity (closed surface S, J_s = -n × ∇_s φ):

        ∫_S J_s · A dS  =  ∫_S ∇_s φ · (n × A) dS
                        =  -∫_S φ · ∇_s · (n × A) dS    [IBP, closed S]
                        =  ∫_S φ · (n · curl A) dS      [∇_s · (n × A) = -n · curl A]
                        =  ∫_S φ · (n · B) dS

    The right-hand side uses B = curl A directly, hence it is gauge-
    invariant under A → A + ∇χ.  The left-hand side is gauge-invariant
    in continuum because the J_s = -n × ∇_s φ field is exactly
    surface-divergence-free, but in a discrete H1 P1 setting J_s_h has
    element-edge jumps that act as a weak surface divergence.  This
    spurious divergence couples to the gauge of A and contaminates
    the imaginary part of Δ L (observed empirically: ~100x off vs the
    energy-balance Δ R = 2 P_wp / I^2 prediction).

    Args:
        filament_paths: K filament polylines (each [(p1, p2), ...])
        currents: (K,) complex per-filament currents I_fil [A] (same as
                  for delta_L_telegen)
        centroids: (N_t, 3) workpiece triangle centroids [m]
        areas: (N_t,) workpiece triangle areas [m^2]
        normals: (N_t, 3) per-triangle outward unit normals.  MUST match
                 the orientation convention used to extract J_s (i.e. the
                 same NGSolve specialcf.normal that fed the J_s = -n × ∇φ
                 expression in the legacy form).
        phi_avg: (N_t,) complex per-triangle area-averaged scalar
                 potential, phi_avg_T = (1/area_T) ∫_T φ_h dT [A].  For
                 P1 phi_h on flat triangles this equals (φ_1+φ_2+φ_3)/3.
        I_port: scalar coil port current [A] (same as delta_L_telegen).
        n_threads: TaskManager thread count for B_inc evaluation.

    Returns:
        delta_L: complex Δ L_port [Henry].

    Implementation: 1-point centroid quadrature for B_inc (smooth on the
    workpiece scale), exact integration of phi_h (P1 → φ_avg · area).
    """
    cpp = _cpp_h_segments_complex()
    if cpp is None:
        raise RuntimeError(
            "C++ _HFromSegmentsComplex unavailable; rebuild radia.")

    # Flatten filament paths to (N_seg, 2, 3) + complex per-segment I.
    flat_segs = []
    flat_I = []
    currents_arr = np.asarray(currents, dtype=complex)
    for fil_segs, Ik in zip(filament_paths, currents_arr):
        Ik_c = complex(Ik)
        for (p1, p2) in fil_segs:
            flat_segs.append([list(p1), list(p2)])
            flat_I.append(Ik_c)
    if not flat_segs:
        return 0.0 + 0.0j
    segs = np.ascontiguousarray(flat_segs, dtype=float)
    I_arr = np.asarray(flat_I, dtype=complex)
    I_re = np.ascontiguousarray(I_arr.real)
    I_im = np.ascontiguousarray(I_arr.imag)

    # B_inc = mu_0 * H_inc at triangle centroids.  H is in A/m;
    # multiplying by mu_0 yields B in Tesla.
    obs = np.ascontiguousarray(centroids, dtype=float)
    H_re, H_im = cpp(segs, obs, I_re, I_im, n_threads)
    B_inc = MU_0 * (H_re + 1j * H_im)            # (N_t, 3) complex

    n = np.asarray(normals, dtype=float)
    n_dot_B = np.sum(n * B_inc, axis=1)          # (N_t,) complex

    a = np.asarray(areas, dtype=float)
    phi_a = np.asarray(phi_avg, dtype=complex)

    # ∫_S phi (n · B_inc) dS  ≈  Σ_T phi_avg_T · (n_T · B_T) · area_T
    integ = np.sum(phi_a * n_dot_B * a)

    Ip = complex(I_port)
    if abs(Ip) < 1e-30:
        return 0.0 + 0.0j
    return integ / (Ip * Ip)
