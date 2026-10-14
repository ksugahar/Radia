/*-------------------------------------------------------------------------
*
* File name:      rad_average_field.h
*
* Project:        RADIA
*
* Description:    Closed-form spatial average of B over a target rectangular
*                 box from a uniform-magnetisation source rectangular box.
*                 Implements the Phase-alpha (2026-05-01) sympy-derived
*                 antiderivatives G1, G2 with a 64-corner inclusion-
*                 exclusion sum.
*
* First release:  2026-05-01
*
* References:
*   - Wakao S., Fujiwara K., Tokumasu T., Kameari A., "Useful Formulas of
*     Analytical Integration in Electromagnetic Field Computations (Part 6)",
*     IEE Japan SA-05-15 / RM-05-15 (2005), section 7 (cuboid average).
*   - Newell A. J., Williams W., Dunlop D. J., "A generalization of the
*     demagnetizing tensor for nonuniform magnetization", J. Geophys. Res.
*     98 (1993), pp. 9551-9555.
*   - memory/project_phase_alpha_g1_g2_derived.md (Radia Phase-alpha record).
*
-------------------------------------------------------------------------*/

#ifndef __RAD_AVERAGE_FIELD_H
#define __RAD_AVERAGE_FIELD_H

namespace radia {
namespace average_field {

/**
 * Compute the 3x3 average demag tensor A_T such that
 *
 *     <B>_T  =  mu_0 * (A_T @ M  +  M * V_overlap / V_T)
 *
 * Closed-form 64-corner sum of the 3-fold antiderivative primitives
 * G1 (diagonal), G2 (off-diagonal). Both source and target are
 * axis-aligned cuboids.
 *
 * Args:
 *   src_min, src_max  3-vector lower/upper corners of the source box [m]
 *   tgt_min, tgt_max  3-vector lower/upper corners of the target box [m]
 *
 * Output:
 *   A_out             pre-allocated double[9] array, row-major [i*3+j].
 *                     Filled with the dimensionless demag tensor.
 */
void AverageDemagTensor(
    const double src_min[3],
    const double src_max[3],
    const double tgt_min[3],
    const double tgt_max[3],
    double A_out[9]);

/**
 * Compute the spatial average of B over the target box from uniform
 * magnetisation in the source box (Tesla).
 *
 *     <B>_T  =  mu_0 * (A_T @ M  +  M * V_overlap / V_T)
 *
 * For non-overlapping source-target, V_overlap = 0; this gives the
 * spatial-averaged "external" field.
 *
 * Args:
 *   M                 magnetisation [A/m]  (3-vector)
 *   src_min, src_max  source box corners [m]
 *   tgt_min, tgt_max  target box corners [m]
 *
 * Output:
 *   B_out             pre-allocated double[3] array. Filled with
 *                     <B_x>, <B_y>, <B_z> in Tesla.
 */
void AverageBInBox(
    const double M[3],
    const double src_min[3],
    const double src_max[3],
    const double tgt_min[3],
    const double tgt_max[3],
    double B_out[3]);

}}  // namespace radia::average_field

#endif  // __RAD_AVERAGE_FIELD_H
