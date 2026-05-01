/*-------------------------------------------------------------------------
*
* File name:      rad_average_field.cpp
*
* Project:        RADIA
*
* Description:    Implementation of the closed-form average-B-in-box
*                 kernel (Phase alpha, 2026-05-01).
*
*                 g1(X,Y,Z)  =  atan(YZ / (X * R))
*                 g2(X,Y,Z)  = -asinh(Z / sqrt(X^2 + Y^2))
*
*                 G1 = X*Y*Z*atan(YZ/(X*R))
*                    + Y*(X^2 - Z^2)/2 * asinh(Y/sqrt(X^2+Z^2))
*                    + Z*(X^2 - Y^2)/2 * asinh(Z/sqrt(X^2+Y^2))
*                    - R*(2X^2 - Y^2 - Z^2) / 6
*
*                 G2 = -X*Y*Z*asinh(Z/sqrt(X^2+Y^2))
*                    + Y*(Y^2 - 3Z^2)/6 * asinh(X/sqrt(Y^2+Z^2))
*                    + X*(X^2 - 3Z^2)/6 * asinh(Y/sqrt(X^2+Z^2))
*                    + X^2*Z/2 * atan(YZ/(X*R))
*                    + Y^2*Z/2 * atan(XZ/(Y*R))
*                    + Z^3/6  * atan(XY/(Z*R))
*                    + X*Y*R / 3
*
*                 Permutation for tensor element (i,j):
*                   xx -> G1(X,Y,Z),  yy -> G1(Y,X,Z),  zz -> G1(Z,X,Y)
*                   xy -> G2(X,Y,Z),  xz -> G2(X,Z,Y),  yz -> G2(Y,Z,X)
*
*                 64-corner sum (8 source x 8 target) gives <A_ij>_T.
*
-------------------------------------------------------------------------*/

#include "rad_average_field.h"

#include <cmath>
#include <algorithm>

namespace radia {
namespace average_field {

// ----- physical constants -------------------------------------------
namespace {
constexpr double MU_0 = 4.0e-7 * 3.141592653589793238462643383279502884;
constexpr double INV_FOUR_PI = 1.0 / (4.0 * 3.141592653589793238462643383279502884);
constexpr double EPS_X = 1.0e-15;  // X = 0 cutoff for X*Y*Z*atan / X^2*atan terms
}

// ----- corner-primitive G1: diagonal-demag (Y, Z symmetric) ---------
static double G1(double X, double Y, double Z)
{
    const double X2 = X * X;
    const double Y2 = Y * Y;
    const double Z2 = Z * Z;
    const double R2 = X2 + Y2 + Z2;
    const double R = std::sqrt(R2);

    // X*Y*Z*atan(YZ/(X*R)): 0 in the X -> 0 limit
    double atan_term = 0.0;
    if (std::fabs(X) > EPS_X)
        atan_term = X * Y * Z * std::atan((Y * Z) / (X * R));

    // (Y/2) * (X^2 - Z^2) * asinh(Y / sqrt(X^2 + Z^2))
    double asinh_y = 0.0;
    if (X2 + Z2 > 0.0)
        asinh_y = 0.5 * Y * (X2 - Z2) * std::asinh(Y / std::sqrt(X2 + Z2));

    // (Z/2) * (X^2 - Y^2) * asinh(Z / sqrt(X^2 + Y^2))
    double asinh_z = 0.0;
    if (X2 + Y2 > 0.0)
        asinh_z = 0.5 * Z * (X2 - Y2) * std::asinh(Z / std::sqrt(X2 + Y2));

    const double R_term = -R * (2.0 * X2 - Y2 - Z2) / 6.0;

    return atan_term + asinh_y + asinh_z + R_term;
}

// ----- corner-primitive G2: off-diagonal (X, Y symmetric) -----------
static double G2(double X, double Y, double Z)
{
    const double X2 = X * X;
    const double Y2 = Y * Y;
    const double Z2 = Z * Z;
    const double R2 = X2 + Y2 + Z2;
    const double R = std::sqrt(R2);

    double t1 = 0.0;
    if (X2 + Y2 > 0.0)
        t1 = -X * Y * Z * std::asinh(Z / std::sqrt(X2 + Y2));

    double t2 = 0.0;
    if (Y2 + Z2 > 0.0)
        t2 = (Y * (Y2 - 3.0 * Z2) / 6.0) * std::asinh(X / std::sqrt(Y2 + Z2));

    double t3 = 0.0;
    if (X2 + Z2 > 0.0)
        t3 = (X * (X2 - 3.0 * Z2) / 6.0) * std::asinh(Y / std::sqrt(X2 + Z2));

    double atan_yz = 0.0;
    if (std::fabs(X) > EPS_X)
        atan_yz = 0.5 * X2 * Z * std::atan((Y * Z) / (X * R));

    double atan_xz = 0.0;
    if (std::fabs(Y) > EPS_X)
        atan_xz = 0.5 * Y2 * Z * std::atan((X * Z) / (Y * R));

    double atan_xy = 0.0;
    if (std::fabs(Z) > EPS_X)
        atan_xy = (Z2 * Z / 6.0) * std::atan((X * Y) / (Z * R));

    const double R_term = X * Y * R / 3.0;

    return t1 + t2 + t3 + atan_yz + atan_xz + atan_xy + R_term;
}

// ----- corner primitive for tensor element (i, j) -------------------
static double G_ij(int i, int j, double X, double Y, double Z)
{
    if (i == j) {
        if (i == 0) return G1(X, Y, Z);
        if (i == 1) return G1(Y, X, Z);
        return G1(Z, X, Y);
    }
    // off-diagonal symmetric: <A_ij> == <A_ji>
    if ((i == 0 && j == 1) || (i == 1 && j == 0)) return G2(X, Y, Z);
    if ((i == 0 && j == 2) || (i == 2 && j == 0)) return G2(X, Z, Y);
    return G2(Y, Z, X);
}

// ----- public: 3x3 average-demag tensor (closed-form 64-corner sum) -
void AverageDemagTensor(
    const double src_min[3],
    const double src_max[3],
    const double tgt_min[3],
    const double tgt_max[3],
    double A_out[9])
{
    const double V_T = (tgt_max[0] - tgt_min[0])
                     * (tgt_max[1] - tgt_min[1])
                     * (tgt_max[2] - tgt_min[2]);
    const double inv_V_T_4pi = INV_FOUR_PI / V_T;

    // Pre-build the 8 signed corners of source and target boxes.
    double src_corners[8][3];
    double src_signs[8];
    double tgt_corners[8][3];
    double tgt_signs[8];

    int idx = 0;
    for (int ix = 0; ix < 2; ++ix)
    for (int iy = 0; iy < 2; ++iy)
    for (int iz = 0; iz < 2; ++iz, ++idx) {
        src_corners[idx][0] = (ix == 0) ? src_min[0] : src_max[0];
        src_corners[idx][1] = (iy == 0) ? src_min[1] : src_max[1];
        src_corners[idx][2] = (iz == 0) ? src_min[2] : src_max[2];
        src_signs[idx] = ((ix + iy + iz) % 2 == 0) ? -1.0 : +1.0;

        tgt_corners[idx][0] = (ix == 0) ? tgt_min[0] : tgt_max[0];
        tgt_corners[idx][1] = (iy == 0) ? tgt_min[1] : tgt_max[1];
        tgt_corners[idx][2] = (iz == 0) ? tgt_min[2] : tgt_max[2];
        tgt_signs[idx] = ((ix + iy + iz) % 2 == 0) ? -1.0 : +1.0;
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double total = 0.0;
            for (int s = 0; s < 8; ++s) {
                const double sx = src_corners[s][0];
                const double sy = src_corners[s][1];
                const double sz = src_corners[s][2];
                const double ss = src_signs[s];
                for (int t = 0; t < 8; ++t) {
                    const double X = tgt_corners[t][0] - sx;
                    const double Y = tgt_corners[t][1] - sy;
                    const double Z = tgt_corners[t][2] - sz;
                    total += ss * tgt_signs[t] * G_ij(i, j, X, Y, Z);
                }
            }
            A_out[i * 3 + j] = total * inv_V_T_4pi;
        }
    }
}

// ----- public: <B>_T = mu_0 (A M + M V_overlap/V_T) -----------------
void AverageBInBox(
    const double M[3],
    const double src_min[3],
    const double src_max[3],
    const double tgt_min[3],
    const double tgt_max[3],
    double B_out[3])
{
    double A[9];
    AverageDemagTensor(src_min, src_max, tgt_min, tgt_max, A);

    const double V_T = (tgt_max[0] - tgt_min[0])
                     * (tgt_max[1] - tgt_min[1])
                     * (tgt_max[2] - tgt_min[2]);
    double overlap[3];
    overlap[0] = std::max(0.0, std::min(src_max[0], tgt_max[0]) - std::max(src_min[0], tgt_min[0]));
    overlap[1] = std::max(0.0, std::min(src_max[1], tgt_max[1]) - std::max(src_min[1], tgt_min[1]));
    overlap[2] = std::max(0.0, std::min(src_max[2], tgt_max[2]) - std::max(src_min[2], tgt_min[2]));
    const double V_overlap = overlap[0] * overlap[1] * overlap[2];
    const double frac = V_overlap / V_T;

    for (int i = 0; i < 3; ++i) {
        double Hi = 0.0;
        for (int j = 0; j < 3; ++j)
            Hi += A[i * 3 + j] * M[j];
        B_out[i] = MU_0 * (Hi + M[i] * frac);
    }
}

}}  // namespace radia::average_field
