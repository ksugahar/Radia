/*
 * rad_green_fullwave.h
 *
 * Full-wave and quasi-static Green's functions for integral equation formulation
 * Supports DC, MQS, EMQS, and Full-wave regimes
 *
 * References:
 * [1] Z. Zhu et al., "Algorithms in FastImp", IEEE TCAD, 2005
 * [2] S. Bilicz et al., "Nonlocal SIBC", ISEM 2023
 * [3] W. C. Gibson, "The Method of Moments in Electromagnetics", 2008
 *
 * Part of Radia project
 */

#ifndef RAD_GREEN_FULLWAVE_H
#define RAD_GREEN_FULLWAVE_H

#include "gmvect.h"
#include "rad_constants.h"
#include <complex>
#include <cmath>

namespace radia {

// Physical constants (using RadConst from rad_constants.h)
constexpr double CONST_C0 = 299792458.0;                       // Speed of light [m/s]
constexpr double CONST_MU0 = RadConst::MU_0;                   // Permeability of free space [H/m]
constexpr double CONST_EPS0 = 8.854187817e-12;                 // Permittivity of free space [F/m]
constexpr double CONST_ETA0 = 376.730313668;                   // Intrinsic impedance of free space [Ohm]
constexpr double CONST_INV_FOUR_PI = RadConst::INV_FOUR_PI;    // 1/(4*pi)

/**
 * @brief Green's function calculator for electromagnetic analysis
 *
 * Computes scalar and vector Green's functions for various formulations:
 * - DC/MQS: G(r) = 1/(4*pi*r) - quasi-static
 * - Full-wave: G(r) = exp(-j*k*r)/(4*pi*r) - including wave propagation
 *
 * Key quantities:
 * - g: Scalar Green's function
 * - grad_g: Gradient of Green's function
 * - A: Vector potential contribution
 * - Phi: Scalar potential contribution
 */
class radTGreenFunction {
public:
    using Complex = std::complex<double>;

    radTGreenFunction();
    ~radTGreenFunction();

    /**
     * @brief Set operating frequency
     * @param frequency Frequency in Hz (0 for DC/quasi-static)
     */
    void SetFrequency(double frequency);

    /**
     * @brief Set material properties
     * @param eps_r Relative permittivity
     * @param mu_r Relative permeability
     * @param sigma Conductivity [S/m]
     */
    void SetMaterial(double eps_r, double mu_r, double sigma);

    /**
     * @brief Get wave number k
     */
    Complex GetWaveNumber() const { return k_; }

    /**
     * @brief Get skin depth
     */
    double GetSkinDepth() const { return skinDepth_; }

    // ========== Scalar Green's Function ==========

    /**
     * @brief Scalar Green's function G(r)
     * DC/MQS: 1/(4*pi*r)
     * Full-wave: exp(-j*k*r)/(4*pi*r)
     * @param r Distance
     * @return G(r)
     */
    Complex ScalarGreen(double r) const;

    /**
     * @brief Gradient of scalar Green's function dG/dr
     * DC/MQS: -1/(4*pi*r^2)
     * Full-wave: -exp(-j*k*r)*(1 + j*k*r)/(4*pi*r^2)
     * @param r Distance
     * @return dG/dr
     */
    Complex ScalarGreenGrad(double r) const;

    /**
     * @brief Full gradient vector grad(G)
     * Returns grad(G) = (dG/dr) * r_hat
     * @param r_vec Vector from source to observation point
     * @return grad(G) as 3-component complex vector
     */
    void ScalarGreenGradVec(const TVector3d& r_vec,
                            Complex& gx, Complex& gy, Complex& gz) const;

    // ========== Vector Potential Kernel ==========

    /**
     * @brief Vector potential kernel A = mu * G * J
     * For MQS: A = (mu/4pi) * J / r
     * For Full-wave: A = (mu/4pi) * J * exp(-jkr) / r
     * @param r Distance
     * @return Kernel value (multiply by J to get A contribution)
     */
    Complex VectorPotentialKernel(double r) const;

    // ========== Scalar Potential Kernel ==========

    /**
     * @brief Scalar potential kernel Phi = (1/eps) * G * rho
     * For MQS: Phi = (1/4pi*eps) * rho / r
     * For Full-wave: Phi = (1/4pi*eps) * rho * exp(-jkr) / r
     * @param r Distance
     * @return Kernel value (multiply by rho to get Phi contribution)
     */
    Complex ScalarPotentialKernel(double r) const;

    // ========== Electric Field Kernels ==========

    /**
     * @brief Electric field from vector potential E_A = -j*omega*A
     * @param r Distance
     * @return Kernel for E_A (multiply by J)
     */
    Complex EFieldFromAPotential(double r) const;

    /**
     * @brief Electric field from scalar potential E_Phi = -grad(Phi)
     * @param r Distance
     * @param dr_dx Derivative of r with respect to x (= x/r)
     * @return Kernel for E_Phi component (multiply by rho)
     */
    Complex EFieldFromPhiPotential(double r, double dr_dx) const;

    // ========== Magnetic Field Kernels ==========

    /**
     * @brief Magnetic field from vector potential B = curl(A)
     * Uses curl(G*J) = grad(G) x J
     * @param r_vec Vector from source to observation point
     * @param J_vec Current density vector
     * @param Bx, By, Bz Output B field components
     */
    void BFieldFromCurrent(const TVector3d& r_vec,
                           const TVector3d& J_vec,
                           Complex& Bx, Complex& By, Complex& Bz) const;

    // ========== Dyadic Green's Function ==========

    /**
     * @brief Dyadic Green's function for EFIE
     * G_bar = (I + grad grad / k^2) * g
     * For PEC surface integral equation
     * @param r_vec Vector from source to observation
     * @param G_bar 3x3 dyadic output
     */
    void DyadicGreen(const TVector3d& r_vec,
                     Complex G_bar[3][3]) const;

    // ========== Surface Integral Kernels ==========

    /**
     * @brief EFIE kernel for surface current K
     * L(K) = j*omega*mu * integral{ G * K } dS'
     *      + (1/j*omega*eps) * integral{ grad(G) * div_s(K) } dS'
     */
    struct EFIEKernel {
        Complex L_self;      // Self-term (diagonal)
        Complex L_mutual;    // Mutual term (off-diagonal)
    };

    /**
     * @brief Compute EFIE kernel for panel pair
     * @param obs_center Observation panel center
     * @param src_center Source panel center
     * @param src_area Source panel area
     * @return EFIE kernel components
     */
    EFIEKernel ComputeEFIEKernel(const TVector3d& obs_center,
                                  const TVector3d& src_center,
                                  double src_area) const;

    // ========== Impedance Extraction Kernels ==========

    /**
     * @brief Inductance extraction kernel (imaginary part of Z)
     * L = Im(Z) / omega
     */
    Complex InductanceKernel(double r) const;

    /**
     * @brief Resistance extraction kernel (real part of Z)
     * R = Re(Z)
     */
    Complex ResistanceKernel(double r, double sigma) const;

    /**
     * @brief Capacitance extraction kernel
     * C = -1 / (omega * Im(1/Z))
     */
    Complex CapacitanceKernel(double r) const;

    // ========== Formulation Selection ==========

    /**
     * @brief Check if quasi-static approximation is valid
     * Valid when k*r << 1 (dimension << wavelength)
     * @param max_dimension Maximum problem dimension
     * @return true if quasi-static is valid
     */
    bool IsQuasiStatic(double max_dimension) const;

    /**
     * @brief Get wavelength
     */
    double GetWavelength() const;

    /**
     * @brief Get electrical size (dimension / wavelength)
     */
    double GetElectricalSize(double dimension) const;

private:
    double frequency_;     // Operating frequency [Hz]
    double omega_;         // Angular frequency [rad/s]
    double eps_r_;         // Relative permittivity
    double mu_r_;          // Relative permeability
    double sigma_;         // Conductivity [S/m]

    Complex k_;            // Wave number
    double skinDepth_;     // Skin depth [m]

    // Cached values
    double invFourPi_;     // 1/(4*pi)
    Complex jOmegaMu_;     // j*omega*mu
    Complex invJOmegaEps_; // 1/(j*omega*eps)

    // Update derived quantities
    void UpdateDerivedQuantities();
};

/**
 * @brief Panel-to-panel interaction calculator
 *
 * Computes interaction integrals between surface panels for MoM.
 * Handles:
 * - Near-field: Full numerical integration
 * - Far-field: Centroid approximation
 * - Self-term: Analytical formulas
 */
class radTPanelInteraction {
public:
    using Complex = std::complex<double>;

    radTPanelInteraction();
    ~radTPanelInteraction();

    /**
     * @brief Set Green's function calculator
     */
    void SetGreenFunction(radTGreenFunction* green) { green_ = green; }

    /**
     * @brief Compute interaction between two triangular panels
     * @param obs_vertices Observation panel vertices (3)
     * @param src_vertices Source panel vertices (3)
     * @return Interaction integral
     */
    Complex TriangleToTriangle(const TVector3d obs_vertices[3],
                                const TVector3d src_vertices[3]) const;

    /**
     * @brief Compute interaction between two quadrilateral panels
     * @param obs_vertices Observation panel vertices (4)
     * @param src_vertices Source panel vertices (4)
     * @return Interaction integral
     */
    Complex QuadToQuad(const TVector3d obs_vertices[4],
                       const TVector3d src_vertices[4]) const;

    /**
     * @brief Compute self-interaction of triangular panel
     * Uses analytical formula for singularity handling
     * @param vertices Triangle vertices (3)
     * @return Self-interaction integral
     */
    Complex TriangleSelf(const TVector3d vertices[3]) const;

    /**
     * @brief Compute self-interaction of quadrilateral panel
     * @param vertices Quad vertices (4)
     * @return Self-interaction integral
     */
    Complex QuadSelf(const TVector3d vertices[4]) const;

    /**
     * @brief Set integration order for numerical quadrature
     * @param order Gaussian quadrature order (1-7)
     */
    void SetIntegrationOrder(int order) { integOrder_ = order; }

private:
    radTGreenFunction* green_;
    int integOrder_;

    // Gaussian quadrature points and weights for triangle
    static const int MAX_GAUSS_PTS = 7;
    double triangleGaussPts_[MAX_GAUSS_PTS][3];   // (xi, eta, weight)
    double quadGaussPts_[MAX_GAUSS_PTS][3];       // (xi, eta, weight)

    void InitGaussPoints();

    // Helper for panel centroid and area
    void ComputeTriangleGeometry(const TVector3d vertices[3],
                                  TVector3d& centroid,
                                  TVector3d& normal,
                                  double& area) const;

    void ComputeQuadGeometry(const TVector3d vertices[4],
                              TVector3d& centroid,
                              TVector3d& normal,
                              double& area) const;
};

// ========== Inline implementations ==========

inline radTGreenFunction::Complex radTGreenFunction::ScalarGreen(double r) const {
    if (r < 1e-15) return Complex(0.0, 0.0);  // Singularity protection

    double gr = invFourPi_ / r;

    if (std::abs(k_) < 1e-15) {
        // DC/MQS: 1/(4*pi*r)
        return Complex(gr, 0.0);
    } else {
        // Full-wave: exp(-jkr)/(4*pi*r)
        Complex jkr = Complex(0.0, 1.0) * k_ * r;
        return std::exp(-jkr) * gr;
    }
}

inline radTGreenFunction::Complex radTGreenFunction::ScalarGreenGrad(double r) const {
    if (r < 1e-15) return Complex(0.0, 0.0);  // Singularity protection

    double r2 = r * r;
    double gr2 = -invFourPi_ / r2;

    if (std::abs(k_) < 1e-15) {
        // DC/MQS: -1/(4*pi*r^2)
        return Complex(gr2, 0.0);
    } else {
        // Full-wave: -exp(-jkr)*(1 + jkr)/(4*pi*r^2)
        Complex jkr = Complex(0.0, 1.0) * k_ * r;
        return std::exp(-jkr) * (1.0 + jkr) * gr2;
    }
}

inline radTGreenFunction::Complex radTGreenFunction::VectorPotentialKernel(double r) const {
    // A kernel = mu * G(r)
    return CONST_MU0 * mu_r_ * ScalarGreen(r);
}

inline radTGreenFunction::Complex radTGreenFunction::ScalarPotentialKernel(double r) const {
    // Phi kernel = G(r) / eps
    return ScalarGreen(r) / (CONST_EPS0 * eps_r_);
}

inline double radTGreenFunction::GetWavelength() const {
    if (frequency_ < 1e-10) return 1e30;  // Effectively infinite for DC
    return CONST_C0 / frequency_;
}

inline double radTGreenFunction::GetElectricalSize(double dimension) const {
    return dimension / GetWavelength();
}

inline bool radTGreenFunction::IsQuasiStatic(double max_dimension) const {
    // Quasi-static valid when electrical size < 0.1
    return GetElectricalSize(max_dimension) < 0.1;
}

} // namespace radia

#endif // RAD_GREEN_FULLWAVE_H
