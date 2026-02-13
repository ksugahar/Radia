/*
 * Radia PEEC - Surface Impedance Boundary Condition (SIBC)
 *
 * Implementation of frequency-dependent surface impedance for conductors.
 *
 * Author: Radia Development Team
 * Date: 2026-02-13
 */

#define _USE_MATH_DEFINES  // For M_PI on Windows
#include "rad_peec_surface_impedance.h"
#include <cmath>
#include <stdexcept>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace radPEEC {

// ============================================================================
// Public Methods
// ============================================================================

std::complex<double> SurfaceImpedance::ComputeRectangular(
    double frequency,
    double sigma,
    double mu_r,
    double width,
    double height)
{
    if (frequency <= 0.0) {
        throw std::invalid_argument("Frequency must be positive");
    }
    if (sigma <= 0.0) {
        throw std::invalid_argument("Conductivity must be positive");
    }
    if (width <= 0.0) {
        throw std::invalid_argument("Width must be positive");
    }

    // Use square cross-section if height not specified
    if (height <= 0.0) {
        height = width;
    }

    // Skin depth
    double delta = SkinDepth(frequency, sigma, mu_r);

    // Surface impedance (classical SIBC formula)
    // Z_s = (1 + j) / (σ δ w)
    // Valid when δ << w, h
    double omega = 2.0 * M_PI * frequency;
    double mu = mu_r * MU_0;

    // Alternative formula: Z_s = √(jωμ/σ) / w
    // This gives the same result for rectangular conductors
    std::complex<double> j(0.0, 1.0);
    std::complex<double> sqrt_jwmu_over_sigma = std::sqrt(j * omega * mu / sigma);

    // For rectangular cross-section, use perimeter-weighted formula
    double perimeter = 2.0 * (width + height);
    double area = width * height;

    // Surface resistance and reactance per unit length
    // R_s = Re{√(jωμ/σ)} * (perimeter / area)
    // X_s = Im{√(jωμ/σ)} * (perimeter / area)
    std::complex<double> Z_s = sqrt_jwmu_over_sigma * (perimeter / area);

    return Z_s;
}

std::complex<double> SurfaceImpedance::ComputeCircular(
    double frequency,
    double sigma,
    double mu_r,
    double radius)
{
    if (frequency <= 0.0) {
        throw std::invalid_argument("Frequency must be positive");
    }
    if (sigma <= 0.0) {
        throw std::invalid_argument("Conductivity must be positive");
    }
    if (radius <= 0.0) {
        throw std::invalid_argument("Radius must be positive");
    }

    double omega = 2.0 * M_PI * frequency;
    double mu = mu_r * MU_0;

    // For circular conductor, use Bessel function correction
    // k = √(jωμσ)
    std::complex<double> j(0.0, 1.0);
    std::complex<double> k_squared = j * omega * mu * sigma;
    double k_abs = std::sqrt(std::abs(k_squared));
    double kr = k_abs * radius;

    // Bessel functions and derivatives
    double ber_kr = BesselBer(kr);
    double bei_kr = BesselBei(kr);
    double ber_prime_kr = BesselBerPrime(kr);
    double bei_prime_kr = BesselBeiPrime(kr);

    // Surface impedance (per unit length)
    // Z_s = (jωμ / 2πr) * [ber(kr) + j·bei(kr)] / [ber'(kr) + j·bei'(kr)]
    std::complex<double> numerator(ber_kr, bei_kr);
    std::complex<double> denominator(ber_prime_kr, bei_prime_kr);
    std::complex<double> Z_s = (j * omega * mu / (2.0 * M_PI * radius)) * (numerator / denominator);

    return Z_s;
}

double SurfaceImpedance::ComputeDowellFactor(
    double frequency,
    double sigma,
    double mu_r,
    double thickness,
    int num_layers)
{
    if (frequency <= 0.0) {
        return 1.0;  // DC: F_R = 1
    }
    if (sigma <= 0.0) {
        throw std::invalid_argument("Conductivity must be positive");
    }
    if (thickness <= 0.0) {
        throw std::invalid_argument("Thickness must be positive");
    }
    if (num_layers < 1) {
        throw std::invalid_argument("Number of layers must be >= 1");
    }

    // Skin depth
    double delta = SkinDepth(frequency, sigma, mu_r);

    // ξ = d/δ (conductor thickness / skin depth)
    double xi = thickness / delta;

    // Dowell's formula for AC resistance factor
    // F_R = R_ac / R_dc

    // Term 1: Skin effect (always present)
    double sinh_2xi = std::sinh(2.0 * xi);
    double sin_2xi = std::sin(2.0 * xi);
    double cosh_2xi = std::cosh(2.0 * xi);
    double cos_2xi = std::cos(2.0 * xi);

    double denominator_skin = cosh_2xi - cos_2xi;
    if (std::abs(denominator_skin) < 1e-12) {
        // Avoid division by zero (should not happen for physical parameters)
        return 1.0;
    }

    double F_skin = xi * (sinh_2xi + sin_2xi) / denominator_skin;

    // Term 2: Proximity effect (only for multi-layer, m > 1)
    double F_proximity = 0.0;
    if (num_layers > 1) {
        double sinh_xi = std::sinh(xi);
        double sin_xi = std::sin(xi);
        double cosh_xi = std::cosh(xi);
        double cos_xi = std::cos(xi);

        double denominator_prox = cosh_xi + cos_xi;
        if (std::abs(denominator_prox) > 1e-12) {
            double m = static_cast<double>(num_layers);
            F_proximity = (2.0 / 3.0) * (m * m - 1.0) * xi * (sinh_xi - sin_xi) / denominator_prox;
        }
    }

    double F_R = F_skin + F_proximity;

    return F_R;
}

double SurfaceImpedance::SkinDepth(
    double frequency,
    double sigma,
    double mu_r)
{
    if (frequency <= 0.0) {
        return std::numeric_limits<double>::infinity();  // DC: infinite skin depth
    }
    if (sigma <= 0.0) {
        throw std::invalid_argument("Conductivity must be positive");
    }

    double omega = 2.0 * M_PI * frequency;
    double mu = mu_r * MU_0;

    // δ = √(2/(ωμσ))
    double delta = std::sqrt(2.0 / (omega * mu * sigma));

    return delta;
}

bool SurfaceImpedance::IsSIBCValid(
    double frequency,
    double sigma,
    double mu_r,
    double characteristic_size)
{
    double delta = SkinDepth(frequency, sigma, mu_r);

    // SIBC is valid when skin depth is much smaller than conductor size
    // Use δ < 0.1 * L as criterion
    return (delta < 0.1 * characteristic_size);
}

double SurfaceImpedance::DCResistancePerLength(
    double sigma,
    double cross_section_area)
{
    if (sigma <= 0.0) {
        throw std::invalid_argument("Conductivity must be positive");
    }
    if (cross_section_area <= 0.0) {
        throw std::invalid_argument("Cross-section area must be positive");
    }

    // R_dc = 1 / (σ * A) [Ohm/m]
    return 1.0 / (sigma * cross_section_area);
}

// ============================================================================
// Private Methods - Bessel Functions
// ============================================================================

double SurfaceImpedance::BesselBer(double x)
{
    // Polynomial approximation for ber(x)
    // Valid for 0 <= x <= 8
    if (x < 0.0) x = -x;  // ber is even function

    if (x <= 8.0) {
        double t = x / 8.0;
        double t2 = t * t;
        double t4 = t2 * t2;
        double t6 = t4 * t2;
        double t8 = t4 * t4;

        // Polynomial coefficients (from Abramowitz & Stegun)
        return 1.0 - 64.0 * t4 + 113.77777774 * t8
               - 32.36345652 * t8 * t4 + 2.64191397 * t8 * t8
               - 0.08349609 * t8 * t8 * t4;
    } else {
        // Asymptotic expansion for large x
        double exp_term = std::exp(x / std::sqrt(2.0));
        double cos_term = std::cos(x / std::sqrt(2.0) - M_PI / 8.0);
        return exp_term * cos_term / std::sqrt(2.0 * M_PI * x);
    }
}

double SurfaceImpedance::BesselBei(double x)
{
    // Polynomial approximation for bei(x)
    // Valid for 0 <= x <= 8
    if (x < 0.0) x = -x;  // bei(-x) = -bei(x), but we handle absolute value

    if (x <= 8.0) {
        double t = x / 8.0;
        double t2 = t * t;
        double t4 = t2 * t2;
        double t6 = t4 * t2;

        // Polynomial coefficients
        return 16.0 * t4 - 113.77777774 * t4 * t4
               + 72.81777742 * t4 * t4 * t2 - 10.56765779 * t4 * t4 * t4
               + 0.52185615 * t4 * t4 * t4 * t2;
    } else {
        // Asymptotic expansion for large x
        double exp_term = std::exp(x / std::sqrt(2.0));
        double sin_term = std::sin(x / std::sqrt(2.0) - M_PI / 8.0);
        return exp_term * sin_term / std::sqrt(2.0 * M_PI * x);
    }
}

double SurfaceImpedance::BesselBerPrime(double x)
{
    // Derivative of ber(x)
    // ber'(x) = -bei(x) - (x/2) * [ber_1(x) + bei_1(x)]
    // Simplified approximation using numerical derivative
    double h = 1e-6;
    return (BesselBer(x + h) - BesselBer(x - h)) / (2.0 * h);
}

double SurfaceImpedance::BesselBeiPrime(double x)
{
    // Derivative of bei(x)
    // bei'(x) = ber(x) - (x/2) * [ber_1(x) - bei_1(x)]
    // Simplified approximation using numerical derivative
    double h = 1e-6;
    return (BesselBei(x + h) - BesselBei(x - h)) / (2.0 * h);
}

} // namespace radPEEC
