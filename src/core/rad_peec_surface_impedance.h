/*
 * Radia PEEC - Surface Impedance Boundary Condition (SIBC)
 *
 * Implements skin effect and proximity effect via surface impedance.
 *
 * Reference:
 *   R. F. Harrington, "Field Computation by Moment Methods", 1968
 *   C. R. Paul, "Inductance: Loop and Partial", 2010
 *
 * Author: Radia Development Team
 * Date: 2026-02-13
 */

#ifndef RAD_PEEC_SURFACE_IMPEDANCE_H
#define RAD_PEEC_SURFACE_IMPEDANCE_H

#define _USE_MATH_DEFINES  // For M_PI on Windows
#include <complex>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace radPEEC {

/**
 * @brief Surface Impedance Boundary Condition (SIBC) calculator
 *
 * Computes frequency-dependent surface impedance for conductors,
 * accounting for skin effect and proximity effect.
 */
class SurfaceImpedance {
public:
    /**
     * @brief Compute surface impedance for rectangular cross-section conductor
     *
     * Uses classical SIBC formula for high-frequency skin effect.
     * Valid when skin depth << conductor dimensions.
     *
     * @param frequency Frequency [Hz]
     * @param sigma Conductivity [S/m]
     * @param mu_r Relative permeability (default 1.0 for copper/aluminum)
     * @param width Conductor width [m]
     * @param height Conductor height [m] (optional, defaults to width)
     * @return Complex surface impedance Z_s = R_s + j*X_s [Ohm]
     */
    static std::complex<double> ComputeRectangular(
        double frequency,
        double sigma,
        double mu_r,
        double width,
        double height = -1.0  // If negative, use width (square cross-section)
    );

    /**
     * @brief Compute surface impedance for circular cross-section conductor
     *
     * Uses Bessel function approximation for round wires.
     * Valid for all frequency ranges (DC to high frequency).
     *
     * @param frequency Frequency [Hz]
     * @param sigma Conductivity [S/m]
     * @param mu_r Relative permeability
     * @param radius Wire radius [m]
     * @return Complex surface impedance Z_s [Ohm/m] (per unit length)
     */
    static std::complex<double> ComputeCircular(
        double frequency,
        double sigma,
        double mu_r,
        double radius
    );

    /**
     * @brief Compute AC resistance factor using Dowell's formula
     *
     * Dowell's formula accurately computes AC resistance for rectangular
     * conductors including skin effect (and optionally proximity effect).
     *
     * For single layer (m=1), proximity effect is zero:
     *   F_R = ξ * [sinh(2ξ) + sin(2ξ)] / [cosh(2ξ) - cos(2ξ)]
     *
     * where ξ = d/δ (conductor thickness / skin depth)
     *
     * Reference: P.L. Dowell, "Effects of eddy currents in transformer windings,"
     *            Proceedings of the IEE, Vol. 113, No. 8, pp. 1387-1394, 1966.
     *
     * @param frequency Frequency [Hz]
     * @param sigma Conductivity [S/m]
     * @param mu_r Relative permeability
     * @param thickness Conductor thickness (height in cross-section) [m]
     * @param num_layers Number of layers (default=1 for single conductor)
     * @return AC resistance factor F_R = R_ac / R_dc
     */
    static double ComputeDowellFactor(
        double frequency,
        double sigma,
        double mu_r,
        double thickness,
        int num_layers = 1
    );

    /**
     * @brief Compute skin depth
     *
     * δ = √(2/(ωμσ))
     *
     * @param frequency Frequency [Hz]
     * @param sigma Conductivity [S/m]
     * @param mu_r Relative permeability
     * @return Skin depth [m]
     */
    static double SkinDepth(
        double frequency,
        double sigma,
        double mu_r = 1.0
    );

    /**
     * @brief Check if SIBC is valid (skin depth << conductor size)
     *
     * SIBC approximation is valid when δ < 0.1 * min(width, height)
     *
     * @param frequency Frequency [Hz]
     * @param sigma Conductivity [S/m]
     * @param mu_r Relative permeability
     * @param characteristic_size Minimum conductor dimension [m]
     * @return true if SIBC is valid, false if volumetric meshing is needed
     */
    static bool IsSIBCValid(
        double frequency,
        double sigma,
        double mu_r,
        double characteristic_size
    );

    /**
     * @brief Compute DC resistance per unit length
     *
     * R_dc = 1 / (σ * A)
     *
     * @param sigma Conductivity [S/m]
     * @param cross_section_area Cross-section area [m^2]
     * @return DC resistance per unit length [Ohm/m]
     */
    static double DCResistancePerLength(
        double sigma,
        double cross_section_area
    );

private:
    static constexpr double MU_0 = 4.0 * M_PI * 1.0e-7;  // H/m

    /**
     * @brief Bessel function ber(x) approximation
     *
     * Used for circular conductor impedance calculation.
     */
    static double BesselBer(double x);

    /**
     * @brief Bessel function bei(x) approximation
     */
    static double BesselBei(double x);

    /**
     * @brief Bessel function derivative ber'(x)
     */
    static double BesselBerPrime(double x);

    /**
     * @brief Bessel function derivative bei'(x)
     */
    static double BesselBeiPrime(double x);
};

} // namespace radPEEC

#endif // RAD_PEEC_SURFACE_IMPEDANCE_H
