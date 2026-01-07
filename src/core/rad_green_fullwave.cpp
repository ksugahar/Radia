/*
 * rad_green_fullwave.cpp
 *
 * Full-wave and quasi-static Green's functions for integral equation formulation
 *
 * References:
 * [1] Z. Zhu et al., "Algorithms in FastImp", IEEE TCAD, 2005
 * [2] S. Bilicz et al., "Nonlocal SIBC", ISEM 2023
 * [3] W. C. Gibson, "The Method of Moments in Electromagnetics", 2008
 *
 * Part of Radia project
 */

#include "rad_green_fullwave.h"
#include <algorithm>

namespace radia {

// ============================================================================
// radTGreenFunction implementation
// ============================================================================

radTGreenFunction::radTGreenFunction()
    : frequency_(0.0)
    , omega_(0.0)
    , eps_r_(1.0)
    , mu_r_(1.0)
    , sigma_(0.0)
    , k_(0.0, 0.0)
    , skinDepth_(1e30)
    , invFourPi_(RadConst::INV_FOUR_PI)
    , jOmegaMu_(0.0, 0.0)
    , invJOmegaEps_(0.0, 0.0)
{
}

radTGreenFunction::~radTGreenFunction() {
}

void radTGreenFunction::SetFrequency(double frequency) {
    frequency_ = frequency;
    omega_ = 2.0 * RadConst::PI * frequency;
    UpdateDerivedQuantities();
}

void radTGreenFunction::SetMaterial(double eps_r, double mu_r, double sigma) {
    eps_r_ = eps_r;
    mu_r_ = mu_r;
    sigma_ = sigma;
    UpdateDerivedQuantities();
}

void radTGreenFunction::UpdateDerivedQuantities() {
    // Wave number: k = omega * sqrt(mu * eps_c)
    // where eps_c = eps - j*sigma/omega (complex permittivity)

    if (frequency_ < 1e-10) {
        // DC case
        k_ = Complex(0.0, 0.0);
        skinDepth_ = 1e30;  // Infinite for DC
        jOmegaMu_ = Complex(0.0, 0.0);
        invJOmegaEps_ = Complex(0.0, 0.0);
        return;
    }

    double mu = CONST_MU0 * mu_r_;
    double eps = CONST_EPS0 * eps_r_;

    // Complex permittivity including conductivity
    Complex eps_c = Complex(eps, -sigma_ / omega_);

    // Wave number
    Complex k_squared = Complex(0.0, 1.0) * omega_ * mu * (Complex(0.0, 1.0) * omega_ * eps_c);
    // k^2 = -omega^2 * mu * eps_c = omega^2 * mu * eps - j*omega*mu*sigma

    k_squared = omega_ * omega_ * mu * eps - Complex(0.0, omega_ * mu * sigma_);
    k_ = std::sqrt(k_squared);

    // Ensure k has positive imaginary part (lossy medium, decaying wave)
    if (k_.imag() < 0) {
        k_ = -k_;
    }

    // Skin depth: delta = sqrt(2 / (omega * mu * sigma))
    if (sigma_ > 1e-10) {
        skinDepth_ = std::sqrt(2.0 / (omega_ * mu * sigma_));
    } else {
        skinDepth_ = 1e30;  // Effectively infinite for non-conductors
    }

    // Cached quantities
    jOmegaMu_ = Complex(0.0, omega_ * mu);
    invJOmegaEps_ = 1.0 / (Complex(0.0, omega_) * eps_c);
}

void radTGreenFunction::ScalarGreenGradVec(const TVector3d& r_vec,
                                            Complex& gx, Complex& gy, Complex& gz) const {
    double r = std::sqrt(r_vec.x * r_vec.x + r_vec.y * r_vec.y + r_vec.z * r_vec.z);

    if (r < 1e-15) {
        gx = gy = gz = Complex(0.0, 0.0);
        return;
    }

    Complex dGdr = ScalarGreenGrad(r);

    // grad(G) = dG/dr * r_hat = dG/dr * r_vec / r
    double invR = 1.0 / r;
    gx = dGdr * r_vec.x * invR;
    gy = dGdr * r_vec.y * invR;
    gz = dGdr * r_vec.z * invR;
}

radTGreenFunction::Complex radTGreenFunction::EFieldFromAPotential(double r) const {
    // E_A = -j*omega*A = -j*omega*mu*G*J
    return -jOmegaMu_ * ScalarGreen(r);
}

radTGreenFunction::Complex radTGreenFunction::EFieldFromPhiPotential(double r, double dr_dx) const {
    // E_Phi = -grad(Phi) = -grad(G/eps * rho) = -(dG/dr)/eps * dr/dx * rho
    return -ScalarGreenGrad(r) / (CONST_EPS0 * eps_r_) * dr_dx;
}

void radTGreenFunction::BFieldFromCurrent(const TVector3d& r_vec,
                                           const TVector3d& J_vec,
                                           Complex& Bx, Complex& By, Complex& Bz) const {
    // B = curl(A) = curl(mu*G*J) = mu * grad(G) x J
    Complex gradGx, gradGy, gradGz;
    ScalarGreenGradVec(r_vec, gradGx, gradGy, gradGz);

    double mu = CONST_MU0 * mu_r_;

    // Cross product: grad(G) x J
    Bx = mu * (gradGy * J_vec.z - gradGz * J_vec.y);
    By = mu * (gradGz * J_vec.x - gradGx * J_vec.z);
    Bz = mu * (gradGx * J_vec.y - gradGy * J_vec.x);
}

void radTGreenFunction::DyadicGreen(const TVector3d& r_vec,
                                     Complex G_bar[3][3]) const {
    double r = std::sqrt(r_vec.x * r_vec.x + r_vec.y * r_vec.y + r_vec.z * r_vec.z);

    if (r < 1e-15) {
        // Return identity * g(0) for self-term
        Complex g0 = ScalarGreen(1e-10);  // Regularized
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                G_bar[i][j] = (i == j) ? g0 : Complex(0.0, 0.0);
            }
        }
        return;
    }

    Complex g = ScalarGreen(r);
    Complex dGdr = ScalarGreenGrad(r);

    // Unit vector
    double invR = 1.0 / r;
    double rx = r_vec.x * invR;
    double ry = r_vec.y * invR;
    double rz = r_vec.z * invR;

    // Second derivative d2G/dr2
    Complex d2Gdr2;
    if (std::abs(k_) < 1e-15) {
        // DC: d2G/dr2 = 2/(4*pi*r^3)
        d2Gdr2 = Complex(2.0 * invFourPi_ / (r * r * r), 0.0);
    } else {
        // Full-wave: d2G/dr2 = exp(-jkr) * (2 + 2jkr - k^2*r^2) / (4*pi*r^3)
        Complex jkr = Complex(0.0, 1.0) * k_ * r;
        Complex expjkr = std::exp(-jkr);
        d2Gdr2 = expjkr * (2.0 + 2.0 * jkr - k_ * k_ * r * r) * invFourPi_ / (r * r * r);
    }

    // Dyadic Green's function: G_bar = (I + grad grad / k^2) * g
    // For Full-wave with k != 0:
    // G_bar_ij = g * delta_ij + (1/k^2) * (d2G/dr2 - dG/dr/r) * r_i * r_j
    //          + (1/k^2) * dG/dr/r * delta_ij

    if (std::abs(k_) < 1e-15) {
        // MQS approximation: just use scalar Green's function
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                G_bar[i][j] = (i == j) ? g : Complex(0.0, 0.0);
            }
        }
    } else {
        // Full-wave dyadic
        Complex k2inv = 1.0 / (k_ * k_);
        Complex dGdr_over_r = dGdr * invR;

        double r_components[3] = {rx, ry, rz};

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                Complex tensor_term = (d2Gdr2 - dGdr_over_r) * r_components[i] * r_components[j];
                Complex identity_term = (i == j) ? (g + k2inv * dGdr_over_r) : Complex(0.0, 0.0);
                G_bar[i][j] = identity_term + k2inv * tensor_term;
            }
        }
    }
}

radTGreenFunction::EFIEKernel radTGreenFunction::ComputeEFIEKernel(
    const TVector3d& obs_center,
    const TVector3d& src_center,
    double src_area) const {

    EFIEKernel kernel;

    TVector3d r_vec;
    r_vec.x = obs_center.x - src_center.x;
    r_vec.y = obs_center.y - src_center.y;
    r_vec.z = obs_center.z - src_center.z;

    double r = std::sqrt(r_vec.x * r_vec.x + r_vec.y * r_vec.y + r_vec.z * r_vec.z);

    if (r < 1e-15) {
        // Self-term: use analytical approximation
        // For circular panel of area A: self-term ~ A / (4*pi) * (ln(A/pi)/2 - 1)
        double effective_radius = std::sqrt(src_area / RadConst::PI);
        kernel.L_self = jOmegaMu_ * invFourPi_ * src_area *
                        (0.5 * std::log(src_area / RadConst::PI) - 1.0);
        kernel.L_mutual = Complex(0.0, 0.0);
    } else {
        // Mutual term
        Complex g = ScalarGreen(r);
        kernel.L_self = Complex(0.0, 0.0);
        kernel.L_mutual = jOmegaMu_ * g * src_area;
    }

    return kernel;
}

radTGreenFunction::Complex radTGreenFunction::InductanceKernel(double r) const {
    // L kernel = mu * G(r)
    return VectorPotentialKernel(r);
}

radTGreenFunction::Complex radTGreenFunction::ResistanceKernel(double r, double sigma) const {
    // Skin effect resistance: Rs = sqrt(omega * mu / (2 * sigma))
    if (sigma < 1e-10) return Complex(0.0, 0.0);

    double mu = CONST_MU0 * mu_r_;
    double Rs = std::sqrt(omega_ * mu / (2.0 * sigma));

    // Return resistance contribution (local surface impedance)
    return Complex(Rs, Rs);  // (1+j) * Rs / sqrt(2) = sqrt(j*omega*mu/sigma)
}

radTGreenFunction::Complex radTGreenFunction::CapacitanceKernel(double r) const {
    // C kernel = eps * G(r)
    return CONST_EPS0 * eps_r_ * ScalarGreen(r);
}

// ============================================================================
// radTPanelInteraction implementation
// ============================================================================

radTPanelInteraction::radTPanelInteraction()
    : green_(nullptr)
    , integOrder_(4)
{
    InitGaussPoints();
}

radTPanelInteraction::~radTPanelInteraction() {
}

void radTPanelInteraction::InitGaussPoints() {
    // Gaussian quadrature points for unit triangle (0 <= xi, eta; xi + eta <= 1)
    // Order 4 rule (6 points)
    // Reference: Cowper, "Gaussian quadrature formulas for triangles", 1973

    // Weights and points for 6-point rule
    double a = 0.445948490915965;
    double b = 0.091576213509771;
    double w1 = 0.111690794839005;
    double w2 = 0.054975871827661;

    triangleGaussPts_[0][0] = a;     triangleGaussPts_[0][1] = a;     triangleGaussPts_[0][2] = w1;
    triangleGaussPts_[1][0] = 1-2*a; triangleGaussPts_[1][1] = a;     triangleGaussPts_[1][2] = w1;
    triangleGaussPts_[2][0] = a;     triangleGaussPts_[2][1] = 1-2*a; triangleGaussPts_[2][2] = w1;
    triangleGaussPts_[3][0] = b;     triangleGaussPts_[3][1] = b;     triangleGaussPts_[3][2] = w2;
    triangleGaussPts_[4][0] = 1-2*b; triangleGaussPts_[4][1] = b;     triangleGaussPts_[4][2] = w2;
    triangleGaussPts_[5][0] = b;     triangleGaussPts_[5][1] = 1-2*b; triangleGaussPts_[5][2] = w2;

    // For quads: tensor product Gauss-Legendre on [-1,1] x [-1,1]
    // 2x2 rule (4 points)
    double gp = 1.0 / std::sqrt(3.0);
    quadGaussPts_[0][0] = -gp; quadGaussPts_[0][1] = -gp; quadGaussPts_[0][2] = 1.0;
    quadGaussPts_[1][0] =  gp; quadGaussPts_[1][1] = -gp; quadGaussPts_[1][2] = 1.0;
    quadGaussPts_[2][0] = -gp; quadGaussPts_[2][1] =  gp; quadGaussPts_[2][2] = 1.0;
    quadGaussPts_[3][0] =  gp; quadGaussPts_[3][1] =  gp; quadGaussPts_[3][2] = 1.0;
}

void radTPanelInteraction::ComputeTriangleGeometry(const TVector3d vertices[3],
                                                    TVector3d& centroid,
                                                    TVector3d& normal,
                                                    double& area) const {
    // Centroid
    centroid.x = (vertices[0].x + vertices[1].x + vertices[2].x) / 3.0;
    centroid.y = (vertices[0].y + vertices[1].y + vertices[2].y) / 3.0;
    centroid.z = (vertices[0].z + vertices[1].z + vertices[2].z) / 3.0;

    // Edges
    TVector3d e1, e2;
    e1.x = vertices[1].x - vertices[0].x;
    e1.y = vertices[1].y - vertices[0].y;
    e1.z = vertices[1].z - vertices[0].z;

    e2.x = vertices[2].x - vertices[0].x;
    e2.y = vertices[2].y - vertices[0].y;
    e2.z = vertices[2].z - vertices[0].z;

    // Normal = e1 x e2
    normal.x = e1.y * e2.z - e1.z * e2.y;
    normal.y = e1.z * e2.x - e1.x * e2.z;
    normal.z = e1.x * e2.y - e1.y * e2.x;

    double norm = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    area = 0.5 * norm;

    if (norm > 1e-15) {
        normal.x /= norm;
        normal.y /= norm;
        normal.z /= norm;
    }
}

void radTPanelInteraction::ComputeQuadGeometry(const TVector3d vertices[4],
                                                TVector3d& centroid,
                                                TVector3d& normal,
                                                double& area) const {
    // Centroid
    centroid.x = (vertices[0].x + vertices[1].x + vertices[2].x + vertices[3].x) / 4.0;
    centroid.y = (vertices[0].y + vertices[1].y + vertices[2].y + vertices[3].y) / 4.0;
    centroid.z = (vertices[0].z + vertices[1].z + vertices[2].z + vertices[3].z) / 4.0;

    // Diagonals for normal computation
    TVector3d d1, d2;
    d1.x = vertices[2].x - vertices[0].x;
    d1.y = vertices[2].y - vertices[0].y;
    d1.z = vertices[2].z - vertices[0].z;

    d2.x = vertices[3].x - vertices[1].x;
    d2.y = vertices[3].y - vertices[1].y;
    d2.z = vertices[3].z - vertices[1].z;

    // Normal = d1 x d2
    normal.x = d1.y * d2.z - d1.z * d2.y;
    normal.y = d1.z * d2.x - d1.x * d2.z;
    normal.z = d1.x * d2.y - d1.y * d2.x;

    double norm = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    area = 0.5 * norm;  // Area = 0.5 * |d1 x d2| for planar quad

    if (norm > 1e-15) {
        normal.x /= norm;
        normal.y /= norm;
        normal.z /= norm;
    }
}

radTPanelInteraction::Complex radTPanelInteraction::TriangleToTriangle(
    const TVector3d obs_vertices[3],
    const TVector3d src_vertices[3]) const {

    if (!green_) return Complex(0.0, 0.0);

    TVector3d obs_centroid, obs_normal;
    double obs_area;
    ComputeTriangleGeometry(obs_vertices, obs_centroid, obs_normal, obs_area);

    TVector3d src_centroid, src_normal;
    double src_area;
    ComputeTriangleGeometry(src_vertices, src_centroid, src_normal, src_area);

    // Check distance for near/far field decision
    TVector3d r_vec;
    r_vec.x = obs_centroid.x - src_centroid.x;
    r_vec.y = obs_centroid.y - src_centroid.y;
    r_vec.z = obs_centroid.z - src_centroid.z;
    double r = std::sqrt(r_vec.x * r_vec.x + r_vec.y * r_vec.y + r_vec.z * r_vec.z);

    double char_size = std::sqrt(std::max(obs_area, src_area));

    if (r > 3.0 * char_size) {
        // Far-field: centroid approximation
        return green_->ScalarGreen(r) * src_area;
    }

    // Near-field: numerical integration
    // Source: edges for parametric coordinates
    TVector3d e1_src, e2_src;
    e1_src.x = src_vertices[1].x - src_vertices[0].x;
    e1_src.y = src_vertices[1].y - src_vertices[0].y;
    e1_src.z = src_vertices[1].z - src_vertices[0].z;

    e2_src.x = src_vertices[2].x - src_vertices[0].x;
    e2_src.y = src_vertices[2].y - src_vertices[0].y;
    e2_src.z = src_vertices[2].z - src_vertices[0].z;

    Complex result(0.0, 0.0);
    int npts = 6;  // Use 6-point rule

    for (int i = 0; i < npts; i++) {
        double xi = triangleGaussPts_[i][0];
        double eta = triangleGaussPts_[i][1];
        double w = triangleGaussPts_[i][2];

        // Source point
        TVector3d src_pt;
        src_pt.x = src_vertices[0].x + xi * e1_src.x + eta * e2_src.x;
        src_pt.y = src_vertices[0].y + xi * e1_src.y + eta * e2_src.y;
        src_pt.z = src_vertices[0].z + xi * e1_src.z + eta * e2_src.z;

        // Distance from observation centroid
        r_vec.x = obs_centroid.x - src_pt.x;
        r_vec.y = obs_centroid.y - src_pt.y;
        r_vec.z = obs_centroid.z - src_pt.z;
        r = std::sqrt(r_vec.x * r_vec.x + r_vec.y * r_vec.y + r_vec.z * r_vec.z);

        result += w * green_->ScalarGreen(r);
    }

    // Jacobian for triangle: 2 * area
    return result * 2.0 * src_area;
}

radTPanelInteraction::Complex radTPanelInteraction::QuadToQuad(
    const TVector3d obs_vertices[4],
    const TVector3d src_vertices[4]) const {

    if (!green_) return Complex(0.0, 0.0);

    TVector3d obs_centroid, obs_normal;
    double obs_area;
    ComputeQuadGeometry(obs_vertices, obs_centroid, obs_normal, obs_area);

    TVector3d src_centroid, src_normal;
    double src_area;
    ComputeQuadGeometry(src_vertices, src_centroid, src_normal, src_area);

    // Check distance
    TVector3d r_vec;
    r_vec.x = obs_centroid.x - src_centroid.x;
    r_vec.y = obs_centroid.y - src_centroid.y;
    r_vec.z = obs_centroid.z - src_centroid.z;
    double r = std::sqrt(r_vec.x * r_vec.x + r_vec.y * r_vec.y + r_vec.z * r_vec.z);

    double char_size = std::sqrt(std::max(obs_area, src_area));

    if (r > 3.0 * char_size) {
        // Far-field
        return green_->ScalarGreen(r) * src_area;
    }

    // Near-field: 2x2 Gauss integration
    Complex result(0.0, 0.0);

    for (int i = 0; i < 4; i++) {
        double xi = quadGaussPts_[i][0];
        double eta = quadGaussPts_[i][1];
        double w = quadGaussPts_[i][2];

        // Bilinear interpolation for quad
        double N0 = 0.25 * (1 - xi) * (1 - eta);
        double N1 = 0.25 * (1 + xi) * (1 - eta);
        double N2 = 0.25 * (1 + xi) * (1 + eta);
        double N3 = 0.25 * (1 - xi) * (1 + eta);

        TVector3d src_pt;
        src_pt.x = N0 * src_vertices[0].x + N1 * src_vertices[1].x +
                   N2 * src_vertices[2].x + N3 * src_vertices[3].x;
        src_pt.y = N0 * src_vertices[0].y + N1 * src_vertices[1].y +
                   N2 * src_vertices[2].y + N3 * src_vertices[3].y;
        src_pt.z = N0 * src_vertices[0].z + N1 * src_vertices[1].z +
                   N2 * src_vertices[2].z + N3 * src_vertices[3].z;

        // Distance
        r_vec.x = obs_centroid.x - src_pt.x;
        r_vec.y = obs_centroid.y - src_pt.y;
        r_vec.z = obs_centroid.z - src_pt.z;
        r = std::sqrt(r_vec.x * r_vec.x + r_vec.y * r_vec.y + r_vec.z * r_vec.z);

        // Jacobian: approximate as constant
        result += w * green_->ScalarGreen(r);
    }

    // Jacobian for bilinear quad: area/4 per Gauss point
    return result * src_area;
}

radTPanelInteraction::Complex radTPanelInteraction::TriangleSelf(
    const TVector3d vertices[3]) const {

    if (!green_) return Complex(0.0, 0.0);

    // Self-term for triangle: analytical formula
    // Reference: Wilton et al., "Potential Integrals for Uniform and Linear
    // Source Distributions on Polygonal and Polyhedral Domains", IEEE TAP, 1984

    TVector3d centroid, normal;
    double area;
    ComputeTriangleGeometry(vertices, centroid, normal, area);

    // For MQS (1/r kernel), self-integral of triangle:
    // I_self = (1/4pi) * integral_S dS'/|r-r'|
    //        ≈ sqrt(A) / (2*pi) * ln(sqrt(A)) for characteristic dimension
    // More accurate: use average distance from centroid

    // Effective radius
    double R_eff = std::sqrt(area / RadConst::PI);

    // Self-term approximation (regularized)
    // For circular panel: I = R_eff * (2*ln(2) - 1) / (4*pi)
    double I_real = R_eff * (2.0 * std::log(2.0) - 1.0) / (4.0 * RadConst::PI);

    Complex g_self(I_real, 0.0);

    // For full-wave, add phase correction
    if (std::abs(green_->GetWaveNumber()) > 1e-15) {
        Complex k = green_->GetWaveNumber();
        // First-order correction
        g_self -= Complex(0.0, 1.0) * k * area / (4.0 * RadConst::PI);
    }

    return g_self;
}

radTPanelInteraction::Complex radTPanelInteraction::QuadSelf(
    const TVector3d vertices[4]) const {

    if (!green_) return Complex(0.0, 0.0);

    TVector3d centroid, normal;
    double area;
    ComputeQuadGeometry(vertices, centroid, normal, area);

    // Similar to triangle self-term
    double R_eff = std::sqrt(area / RadConst::PI);
    double I_real = R_eff * (2.0 * std::log(2.0) - 1.0) / (4.0 * RadConst::PI);

    Complex g_self(I_real, 0.0);

    if (std::abs(green_->GetWaveNumber()) > 1e-15) {
        Complex k = green_->GetWaveNumber();
        g_self -= Complex(0.0, 1.0) * k * area / (4.0 * RadConst::PI);
    }

    return g_self;
}

} // namespace radia
