#include "axi_henrotte_numeric.hpp"

#include <cmath>
#include <stdexcept>

namespace axifem::numeric {
namespace {

constexpr double EPS_AXIS = 1.0e-14;
constexpr double PI = 3.14159265358979323846;

using CoefficientMatrix = std::array<double, 16>;

void ValidateInputs(double ra, double rb, double za, double zb,
                    double mu, double sigma) {
    if (!std::isfinite(ra) || !std::isfinite(rb) ||
        !std::isfinite(za) || !std::isfinite(zb) ||
        !std::isfinite(mu) || !std::isfinite(sigma))
        throw std::invalid_argument("axifem Q1 inputs must be finite");
    if (ra < 0.0 || rb <= ra)
        throw std::invalid_argument("axifem Q1 requires 0 <= ra < rb");
    if (zb <= za)
        throw std::invalid_argument("axifem Q1 requires za < zb");
    if (mu <= 0.0)
        throw std::invalid_argument("axifem Q1 permeability must be positive");
    if (sigma < 0.0)
        throw std::invalid_argument("axifem Q1 conductivity must be nonnegative");
}

// Row k, column i contains coefficient k of nodal basis function i in
// {1, s, z, s*z}, with s=r^2. This is the analytic inverse Vandermonde.
CoefficientMatrix Q1InverseVandermondeUnchecked(double ra, double rb,
                                                double za, double zb) {
    const double sa = ra * ra;
    const double sb = rb * rb;
    const double scale = 1.0 / ((sb - sa) * (zb - za));
    return {
         sb * zb * scale, -sa * zb * scale,  sa * za * scale, -sb * za * scale,
             -zb * scale,       zb * scale,      -za * scale,       za * scale,
             -sb * scale,       sa * scale,      -sa * scale,       sb * scale,
                  scale,            -scale,            scale,            -scale,
    };
}

Matrix4 Q1Stiffness(double ra, double rb, double za, double zb, double mu) {
    const CoefficientMatrix inv_v =
        Q1InverseVandermondeUnchecked(ra, rb, za, zb);
    const double r_nodes[4] = {ra, rb, rb, ra};
    static constexpr double gp[4] = {
        -0.8611363115940526, -0.3399810435848563,
         0.3399810435848563,  0.8611363115940526,
    };
    static constexpr double gw[4] = {
        0.3478548451374538, 0.6521451548625461,
        0.6521451548625461, 0.3478548451374538,
    };

    Matrix4 integral{};
    for (int ir = 0; ir < 4; ++ir) {
        for (int iz = 0; iz < 4; ++iz) {
            const double rq = 0.5 * (ra + rb) + 0.5 * (rb - ra) * gp[ir];
            const double zq = 0.5 * (za + zb) + 0.5 * (zb - za) * gp[iz];
            if (rq <= EPS_AXIS)
                continue;
            const double weight =
                gw[ir] * gw[iz] * 0.25 * (rb - ra) * (zb - za);
            const double s = rq * rq;
            double dpsi_dr[4];
            double dpsi_dz[4];
            for (int i = 0; i < 4; ++i) {
                const double b = inv_v[4 + i];
                const double c = inv_v[8 + i];
                const double d = inv_v[12 + i];
                dpsi_dr[i] = 2.0 * rq * (b + d * zq);
                dpsi_dz[i] = c + d * s;
            }
            const double factor = weight / rq;
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                    integral[4 * i + j] += factor *
                        (dpsi_dr[i] * dpsi_dr[j] +
                         dpsi_dz[i] * dpsi_dz[j]);
        }
    }

    Matrix4 result{};
    const double scale = 2.0 * PI / mu;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            const double value =
                scale * r_nodes[i] * r_nodes[j] * integral[4 * i + j];
            const double transpose =
                scale * r_nodes[j] * r_nodes[i] * integral[4 * j + i];
            result[4 * i + j] = 0.5 * (value + transpose);
        }
    return result;
}

Matrix4 Q1SigmaMass(double ra, double rb, double za, double zb,
                    double sigma) {
    const double sa = ra * ra;
    const double sb = rb * rb;
    const double ds1 = sb - sa;
    const double ds2 = 0.5 * (sb * sb - sa * sa);
    const double iz0 = zb - za;
    const double iz1 = 0.5 * (zb * zb - za * za);
    const double iz2 = (zb * zb * zb - za * za * za) / 3.0;

    CoefficientMatrix coefficient{};
    if (sa < EPS_AXIS) {
        coefficient[5] = iz0 * 0.5 * sb * sb;
        coefficient[7] = iz1 * 0.5 * sb * sb;
        coefficient[13] = coefficient[7];
        coefficient[15] = iz2 * 0.5 * sb * sb;
    } else {
        const double log_ratio = std::log(sb / sa);
        coefficient[0] = iz0 * log_ratio;
        coefficient[1] = coefficient[4] = iz0 * ds1;
        coefficient[2] = coefficient[8] = iz1 * log_ratio;
        coefficient[3] = coefficient[12] = iz1 * ds1;
        coefficient[5] = iz0 * ds2;
        coefficient[6] = coefficient[9] = iz1 * ds1;
        coefficient[7] = coefficient[13] = iz1 * ds2;
        coefficient[10] = iz2 * log_ratio;
        coefficient[11] = coefficient[14] = iz2 * ds1;
        coefficient[15] = iz2 * ds2;
    }
    const double coefficient_scale = sigma / (4.0 * PI);
    for (double& value : coefficient)
        value *= coefficient_scale;

    const CoefficientMatrix inv_v =
        Q1InverseVandermondeUnchecked(ra, rb, za, zb);
    Matrix4 nodal{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k)
                for (int l = 0; l < 4; ++l)
                    nodal[4 * i + j] +=
                        inv_v[4 * k + i] * coefficient[4 * k + l] *
                        inv_v[4 * l + j];

    const double r_nodes[4] = {ra, rb, rb, ra};
    Matrix4 result{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            const double ti = 2.0 * PI * r_nodes[i];
            const double tj = 2.0 * PI * r_nodes[j];
            result[4 * i + j] =
                0.5 * ti * tj * (nodal[4 * i + j] + nodal[4 * j + i]);
        }
    return result;
}

}  // namespace

Matrix4 ComputeQ1InverseVandermonde(double ra, double rb,
                                    double za, double zb) {
    ValidateInputs(ra, rb, za, zb, 1.0, 0.0);
    return Q1InverseVandermondeUnchecked(ra, rb, za, zb);
}

Matrix4 ComputeQ1MagneticStiffness(double ra, double rb, double za,
                                   double zb, double mu) {
    ValidateInputs(ra, rb, za, zb, mu, 0.0);
    return Q1Stiffness(ra, rb, za, zb, mu);
}

Matrix4 ComputeQ1SigmaMass(double ra, double rb, double za, double zb,
                           double sigma) {
    ValidateInputs(ra, rb, za, zb, 1.0, sigma);
    return Q1SigmaMass(ra, rb, za, zb, sigma);
}

Q1MagneticElementMatrices ComputeQ1MagneticElementMatrices(
    double ra, double rb, double za, double zb, double mu, double sigma) {
    ValidateInputs(ra, rb, za, zb, mu, sigma);
    return {Q1Stiffness(ra, rb, za, zb, mu),
            Q1SigmaMass(ra, rb, za, zb, sigma)};
}

}  // namespace axifem::numeric
