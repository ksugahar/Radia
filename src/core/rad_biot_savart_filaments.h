/*-------------------------------------------------------------------------
*
* File name:      rad_biot_savart_filaments.h
*
* Project:        RADIA
*
* Description:    Fast finite-segment Biot-Savart H field from a collection
*                 of straight wire segments at arbitrary observation points.
*                 Supports COMPLEX per-segment currents (for AC PEEC bundles
*                 with skin-effect / proximity-effect amplitude per filament).
*
* First release:  2026-04-30
*
* Used by:        compute_phi_inc_from_filaments() in bem_sibc_solver.py
*                 to evaluate the incident magnetic scalar potential at
*                 workpiece surface vertices via the path-integral
*                 phi(P) = - int_{far -> P} H . dl.
*
-------------------------------------------------------------------------*/

#ifndef __RAD_BIOT_SAVART_FILAMENTS_H
#define __RAD_BIOT_SAVART_FILAMENTS_H

#include <cstdint>

namespace radia { namespace bs {

/**
 * Finite-segment Biot-Savart H at observation points with complex per-segment
 * currents.  Real and imaginary parts of H are returned separately.
 *
 * Kernel: for each segment k (endpoints p1_k, p2_k) and obs point j,
 *   e_l = (p2 - p1) / |p2 - p1|
 *   r1 = obs - p1,  r2 = obs - p2
 *   d = |e_l x r1|
 *   cos_a1 = (r1 . e_l) / |r1|
 *   cos_a2 = (r2 . e_l) / |r2|
 *   H_k(obs) = I_k / (4 pi d) * (cos_a1 - cos_a2) * (e_l x r1) / d
 * Sum over k.  Skips degenerate segments and on-segment obs points.
 *
 * Complex I = I_re + j*I_im.  H = H_re + j*H_im, distributed over the
 * sum component-wise.
 *
 * Parallelization: outer loop over obs points (lock-free, independent).
 *
 * @param segments   (n_seg, 2, 3) row-major: [k][end][xyz]; double
 * @param n_seg      number of segments
 * @param obs        (n_obs, 3) row-major; double
 * @param n_obs      number of obs points
 * @param I_re,I_im  (n_seg,) real and imaginary parts of segment currents
 * @param H_re_out   (n_obs, 3) row-major: real part of H output (overwritten)
 * @param H_im_out   (n_obs, 3) row-major: imag part of H output (overwritten)
 * @param n_threads  0 = TaskManager default; > 0 = request that many
 */
void HFromSegmentsComplex(
    const double* segments, int n_seg,
    const double* obs,      int n_obs,
    const double* I_re,
    const double* I_im,
    double* H_re_out,
    double* H_im_out,
    int n_threads);

}}  // namespace radia::bs

#endif  // __RAD_BIOT_SAVART_FILAMENTS_H
