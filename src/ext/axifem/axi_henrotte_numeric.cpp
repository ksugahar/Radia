#include "axi_henrotte_numeric.hpp"
#include "q2_henrotte_generated.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace axifem::numeric {
namespace {

constexpr double EPS_AXIS = 1.0e-14;
constexpr double PI = 3.14159265358979323846;

using CoefficientMatrix = std::array<double, 16>;

template <std::size_t N>
using Square = std::array<double, N * N>;

template <std::size_t N>
Square<N> Invert(const Square<N>& matrix, const char* message) {
    std::array<double, N * 2 * N> work{};
    for (std::size_t row = 0; row < N; ++row) {
        for (std::size_t col = 0; col < N; ++col) {
            const double value = matrix[N * row + col];
            work[(2 * N) * row + col] = value;
        }
        work[(2 * N) * row + N + row] = 1.0;
    }
    for (std::size_t pivot_col = 0; pivot_col < N; ++pivot_col) {
        std::size_t pivot_row = pivot_col;
        double pivot_size = std::abs(work[(2 * N) * pivot_row + pivot_col]);
        for (std::size_t row = pivot_col + 1; row < N; ++row) {
            const double candidate = std::abs(work[(2 * N) * row + pivot_col]);
            if (candidate > pivot_size) {
                pivot_size = candidate;
                pivot_row = row;
            }
        }
        // The physical monomial basis contains s^2 and z^2. Millimetre-scale
        // elements therefore have valid pivots far below machine epsilon even
        // though ra < rb and za < zb make this tensor-product Vandermonde
        // analytically nonsingular. Reject only an actual zero/underflow pivot.
        if (pivot_size < std::numeric_limits<double>::min() * N)
            throw std::invalid_argument(message);
        if (pivot_row != pivot_col)
            for (std::size_t col = 0; col < 2 * N; ++col)
                std::swap(work[(2 * N) * pivot_row + col],
                          work[(2 * N) * pivot_col + col]);
        const double pivot = work[(2 * N) * pivot_col + pivot_col];
        for (std::size_t col = 0; col < 2 * N; ++col)
            work[(2 * N) * pivot_col + col] /= pivot;
        for (std::size_t row = 0; row < N; ++row) {
            if (row == pivot_col)
                continue;
            const double factor = work[(2 * N) * row + pivot_col];
            for (std::size_t col = 0; col < 2 * N; ++col)
                work[(2 * N) * row + col] -=
                    factor * work[(2 * N) * pivot_col + col];
        }
    }
    Square<N> inverse{};
    for (std::size_t row = 0; row < N; ++row)
        for (std::size_t col = 0; col < N; ++col)
            inverse[N * row + col] = work[(2 * N) * row + N + col];
    return inverse;
}

void Q2GeneralMonomials(double s, double z, double values[9]) {
    values[0] = 1.0;
    values[1] = s;
    values[2] = s * s;
    values[3] = z;
    values[4] = s * z;
    values[5] = s * s * z;
    values[6] = z * z;
    values[7] = s * z * z;
    values[8] = s * s * z * z;
}

void Q2AxisMonomials(double s, double z, double values[6]) {
    values[0] = s;
    values[1] = s * s;
    values[2] = s * z;
    values[3] = s * s * z;
    values[4] = s * z * z;
    values[5] = s * s * z * z;
}

template <std::size_t N>
Matrix9 Q2MonomialToVDof(
    const Square<N>& monomial,
    const Square<N>& inverse,
    const std::array<int, N>& active,
    const std::array<double, 9>& radius) {
    Square<N> temporary{};
    for (std::size_t row = 0; row < N; ++row)
        for (std::size_t col = 0; col < N; ++col)
            for (std::size_t inner = 0; inner < N; ++inner)
                temporary[N * row + col] +=
                    monomial[N * row + inner] * inverse[N * inner + col];

    Matrix9 result{};
    for (std::size_t local_i = 0; local_i < N; ++local_i) {
        const int i = active[local_i];
        for (std::size_t local_j = 0; local_j < N; ++local_j) {
            const int j = active[local_j];
            double value = 0.0;
            for (std::size_t inner = 0; inner < N; ++inner)
                value += inverse[N * inner + local_i] *
                         temporary[N * inner + local_j];
            result[9 * i + j] =
                (2.0 * PI * radius[i]) * value * (2.0 * PI * radius[j]);
        }
    }
    for (int i = 0; i < 9; ++i)
        for (int j = i + 1; j < 9; ++j) {
            const double average = 0.5 * (result[9 * i + j] + result[9 * j + i]);
            result[9 * i + j] = average;
            result[9 * j + i] = average;
        }
    return result;
}

Q2MagneticElementMatrices Q2Matrices(
    double ra, double rb, double za, double zb, double mu, double sigma) {
    const double sa = ra * ra;
    const double sb = rb * rb;
    const double sm = 0.5 * (sa + sb);
    const double zm = 0.5 * (za + zb);
    const double rm = std::sqrt(sm);
    const std::array<double, 9> radius =
        {ra, rb, rb, ra, rm, rb, rm, ra, rm};
    const std::array<double, 9> s_nodes =
        {sa, sb, sb, sa, sm, sb, sm, sa, sm};
    const std::array<double, 9> z_nodes =
        {za, za, zb, zb, za, zm, zb, zm, zm};
    const bool axis_touching = ra < EPS_AXIS;

    if (!axis_touching) {
        Square<9> vandermonde{};
        for (int row = 0; row < 9; ++row) {
            double monomial[9];
            Q2GeneralMonomials(s_nodes[row], z_nodes[row], monomial);
            for (int col = 0; col < 9; ++col)
                vandermonde[9 * row + col] = monomial[col];
        }
        const auto inverse = Invert<9>(
            vandermonde, "axifem Q2 interior Vandermonde is singular");
        Square<9> stiffness_monomial{};
        Square<9> mass_monomial{};
        q2_henrotte::KPhiGeneral(
            sa, sb, za, zb, mu, mu, stiffness_monomial.data());
        q2_henrotte::MSigmaPhiGeneral(
            sa, sb, za, zb, sigma, mass_monomial.data());
        const std::array<int, 9> active = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        return {
            Q2MonomialToVDof(stiffness_monomial, inverse, active, radius),
            Q2MonomialToVDof(mass_monomial, inverse, active, radius),
            false,
        };
    }

    const std::array<int, 6> active = {1, 2, 4, 5, 6, 8};
    Square<6> vandermonde{};
    for (int row = 0; row < 6; ++row) {
        double monomial[6];
        const int local = active[row];
        Q2AxisMonomials(s_nodes[local], z_nodes[local], monomial);
        for (int col = 0; col < 6; ++col)
            vandermonde[6 * row + col] = monomial[col];
    }
    const auto inverse = Invert<6>(
        vandermonde, "axifem Q2 axis Vandermonde is singular");
    Square<6> stiffness_monomial{};
    Square<6> mass_monomial{};
    q2_henrotte::KPhiAxis(sb, za, zb, mu, mu, stiffness_monomial.data());
    q2_henrotte::MSigmaPhiAxis(sb, za, zb, sigma, mass_monomial.data());
    return {
        Q2MonomialToVDof(stiffness_monomial, inverse, active, radius),
        Q2MonomialToVDof(mass_monomial, inverse, active, radius),
        true,
    };
}

void ValidateInputs(double ra, double rb, double za, double zb,
                    double mu, double sigma) {
    if (!std::isfinite(ra) || !std::isfinite(rb) ||
        !std::isfinite(za) || !std::isfinite(zb) ||
        !std::isfinite(mu) || !std::isfinite(sigma))
        throw std::invalid_argument("axifem element inputs must be finite");
    if (ra < 0.0 || rb <= ra)
        throw std::invalid_argument("axifem element requires 0 <= ra < rb");
    if (zb <= za)
        throw std::invalid_argument("axifem element requires za < zb");
    if (mu <= 0.0)
        throw std::invalid_argument("axifem element permeability must be positive");
    if (sigma < 0.0)
        throw std::invalid_argument("axifem element conductivity must be nonnegative");
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

Q2MagneticElementMatrices ComputeQ2MagneticElementMatrices(
    double ra, double rb, double za, double zb, double mu, double sigma) {
    ValidateInputs(ra, rb, za, zb, mu, sigma);
    return Q2Matrices(ra, rb, za, zb, mu, sigma);
}

}  // namespace axifem::numeric
