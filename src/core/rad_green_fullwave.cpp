/*
 * rad_green_fullwave.cpp
 *
 * Quasi-static Green's functions for integral equation formulation
 * Uses Laplace kernel (1/R) only - Helmholtz kernel removed.
 *
 * References:
 * [1] Z. Zhu et al., "Algorithms in FastImp: A Fast and Accurate
 *     Impedance Extraction Program", IEEE TCAD, Vol. 24, No. 7, 2005
 * [2] W. C. Gibson, "The Method of Moments in Electromagnetics",
 *     Chapman & Hall/CRC, 2008
 * [3] D. R. Wilton, S. M. Rao, A. W. Glisson, D. H. Schaubert,
 *     O. M. Al-Bundak, C. M. Butler, "Potential Integrals for Uniform
 *     and Linear Source Distributions on Polygonal and Polyhedral
 *     Domains", IEEE Trans. Antennas Propag., Vol. 32, No. 3, 1984
 *     (Analytical 1/R integrals over triangular panels)
 * [4] R. D. Graglia, "On the Numerical Integration of the Linear
 *     Shape Functions Times the 3-D Green's Function or its Gradient
 *     on a Plane Triangle", IEEE Trans. AP, Vol. 41, No. 10, 1993
 *
 * Part of Radia project
 */

#include "rad_green_fullwave.h"
#include <algorithm>

namespace radia {

// ============================================================================
// Wilton Analytical 1/R Integral over Triangle (Ref [3])
// ============================================================================
//
// Computes: I = integral_triangle 1/|r - r'| dS'
//
// Using the Wilton et al. formula, the integral is decomposed into
// edge contributions involving logarithms and arctangents:
//
// I = sum_edges [ P0 * ln((R+ + l+)/(R- + l-))
//               - |h| * (atan2(P0*l+, R0^2 + |h|*R+) - atan2(P0*l-, R0^2 + |h|*R-)) ]
//
// where for each edge from vertex A to vertex B:
//   P0 = perpendicular distance from observation point projection to edge
//   R0 = sqrt(P0^2 + h^2), h = height of obs point above triangle plane
//   l+, l- = projections along edge direction from A and B
//   R+, R- = 3D distances from obs point to vertices A and B
//
// This analytical formula is exact for the 1/R (Laplace) kernel and avoids
// numerical integration, providing both accuracy and efficiency.
// ============================================================================

static double WiltonTriangleIntegral1OverR(
    const TVector3d& obs,
    const TVector3d& v0, const TVector3d& v1, const TVector3d& v2)
{
    const double EPS = 1.0e-14;

    // Compute triangle normal and area
    TVector3d e1 = {v1.x - v0.x, v1.y - v0.y, v1.z - v0.z};
    TVector3d e2 = {v2.x - v0.x, v2.y - v0.y, v2.z - v0.z};
    TVector3d normal = {
        e1.y * e2.z - e1.z * e2.y,
        e1.z * e2.x - e1.x * e2.z,
        e1.x * e2.y - e1.y * e2.x
    };
    double norm_mag = std::sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);

    if (norm_mag < EPS) return 0.0;  // Degenerate triangle

    // Unit normal
    normal.x /= norm_mag;
    normal.y /= norm_mag;
    normal.z /= norm_mag;

    // Vector from v0 to observation point
    TVector3d r0 = {obs.x - v0.x, obs.y - v0.y, obs.z - v0.z};

    // Height h: signed distance from obs point to triangle plane
    double h = r0.x * normal.x + r0.y * normal.y + r0.z * normal.z;
    double absH = std::abs(h);
    double h2 = h * h;

    // Projection of observation point onto triangle plane
    TVector3d obs_proj = {
        obs.x - h * normal.x,
        obs.y - h * normal.y,
        obs.z - h * normal.z
    };

    // Store vertices in array for edge loop
    const TVector3d* verts[3] = {&v0, &v1, &v2};

    // Compute distances from obs to each vertex
    double R[3];
    for (int i = 0; i < 3; i++) {
        double dx = obs.x - verts[i]->x;
        double dy = obs.y - verts[i]->y;
        double dz = obs.z - verts[i]->z;
        R[i] = std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    // Accumulate edge contributions
    double I_total = 0.0;

    for (int i = 0; i < 3; i++) {
        int j = (i + 1) % 3;

        // Edge from vertex i to vertex j
        TVector3d edge = {
            verts[j]->x - verts[i]->x,
            verts[j]->y - verts[i]->y,
            verts[j]->z - verts[i]->z
        };
        double L = std::sqrt(edge.x*edge.x + edge.y*edge.y + edge.z*edge.z);
        if (L < EPS) continue;

        // Unit tangent along edge
        double tx = edge.x / L;
        double ty = edge.y / L;
        double tz = edge.z / L;

        // Inward unit normal to edge in triangle plane: m = n x t
        double mx = normal.y * tz - normal.z * ty;
        double my = normal.z * tx - normal.x * tz;
        double mz = normal.x * ty - normal.y * tx;

        // P0: perpendicular distance from obs_proj to edge line (in plane)
        double dx_proj = obs_proj.x - verts[i]->x;
        double dy_proj = obs_proj.y - verts[i]->y;
        double dz_proj = obs_proj.z - verts[i]->z;
        double P0 = dx_proj * mx + dy_proj * my + dz_proj * mz;

        // R0_edge = sqrt(P0^2 + h^2): perpendicular distance in 3D
        double R0_sq = P0 * P0 + h2;

        // l+ = projection of (obs - vertex_i) onto edge direction
        double dx_i = obs.x - verts[i]->x;
        double dy_i = obs.y - verts[i]->y;
        double dz_i = obs.z - verts[i]->z;
        double l_plus = dx_i * tx + dy_i * ty + dz_i * tz;

        // l- = projection of (obs - vertex_j) onto edge direction
        double dx_j = obs.x - verts[j]->x;
        double dy_j = obs.y - verts[j]->y;
        double dz_j = obs.z - verts[j]->z;
        double l_minus = dx_j * tx + dy_j * ty + dz_j * tz;

        double R_plus = R[i];
        double R_minus = R[j];

        // Log term: P0 * ln((R+ + l+) / (R- + l-))
        double log_term = 0.0;
        double num = R_plus + l_plus;
        double den = R_minus + l_minus;
        if (std::abs(num) > EPS && std::abs(den) > EPS && std::abs(num / den) > EPS) {
            log_term = P0 * std::log(std::abs(num / den));
        }

        // Arctan term: |h| * [atan2(P0*l+, R0^2 + |h|*R+) - atan2(P0*l-, R0^2 + |h|*R-)]
        double atan_term = 0.0;
        if (absH > EPS) {
            double atan_plus = std::atan2(P0 * l_plus, R0_sq + absH * R_plus);
            double atan_minus = std::atan2(P0 * l_minus, R0_sq + absH * R_minus);
            atan_term = absH * (atan_plus - atan_minus);
        }

        // Edge contribution (note: sign convention for inward normal)
        I_total += log_term - atan_term;
    }

    return I_total;
}

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

    // Laplace kernel only: use scalar Green's function times identity
    // (Full-wave dyadic removed)
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            G_bar[i][j] = (i == j) ? g : Complex(0.0, 0.0);
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

    // Near-field: Use Wilton analytical formula (Ref [3])
    // The Wilton formula computes the integral of 1/R over the source triangle
    // exactly using edge-based logarithm and arctangent terms.
    // This is more accurate and often faster than Gauss quadrature.
    double I_wilton = WiltonTriangleIntegral1OverR(
        obs_centroid, src_vertices[0], src_vertices[1], src_vertices[2]);

    // For Laplace kernel: G(r) = 1/(4*pi*r)
    // The Wilton formula gives integral of 1/r, so multiply by 1/(4*pi)
    return Complex(I_wilton * RadConst::INV_FOUR_PI, 0.0);
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
        // Far-field: centroid approximation
        return green_->ScalarGreen(r) * src_area;
    }

    // Near-field: Use Wilton analytical formula (Ref [3])
    // Split source quad into two triangles: (v0, v1, v2) and (v0, v2, v3)
    double I_tri1 = WiltonTriangleIntegral1OverR(
        obs_centroid, src_vertices[0], src_vertices[1], src_vertices[2]);
    double I_tri2 = WiltonTriangleIntegral1OverR(
        obs_centroid, src_vertices[0], src_vertices[2], src_vertices[3]);

    double I_total = I_tri1 + I_tri2;

    // For Laplace kernel: G(r) = 1/(4*pi*r)
    return Complex(I_total * RadConst::INV_FOUR_PI, 0.0);
}

radTPanelInteraction::Complex radTPanelInteraction::TriangleSelf(
    const TVector3d vertices[3]) const {

    if (!green_) return Complex(0.0, 0.0);

    // Self-term for triangle: Wilton analytical formula
    // Reference: D.R. Wilton, S.M. Rao, A.W. Glisson, D.H. Schaubert,
    //   O.M. Al-Bundak, C.M. Butler, "Potential Integrals for Uniform and Linear
    //   Source Distributions on Polygonal and Polyhedral Domains",
    //   IEEE Trans. Antennas Propag., Vol. 32, No. 3, March 1984
    //
    // For self-term, observation point is at centroid (h = 0).
    // The Wilton formula simplifies: arctan terms vanish when h -> 0,
    // leaving only the logarithm edge contributions:
    //   I = sum_edges [ P0 * ln((R+ + l+)/(R- + l-)) ]

    TVector3d centroid, normal;
    double area;
    ComputeTriangleGeometry(vertices, centroid, normal, area);

    // Use Wilton analytical formula with observation point at centroid
    // This gives exact result for the 1/R integral over the triangle
    double I_wilton = WiltonTriangleIntegral1OverR(
        centroid, vertices[0], vertices[1], vertices[2]);

    // For Laplace kernel: G(r) = 1/(4*pi*r)
    return Complex(I_wilton * RadConst::INV_FOUR_PI, 0.0);
}

radTPanelInteraction::Complex radTPanelInteraction::QuadSelf(
    const TVector3d vertices[4]) const {

    if (!green_) return Complex(0.0, 0.0);

    // Self-term for quad: Split into two triangles and use Wilton formula
    // Reference: D.R. Wilton et al., IEEE TAP, Vol. 32, No. 3, 1984
    //
    // Quad vertices: v0, v1, v2, v3 (counterclockwise)
    // Split into triangles: (v0, v1, v2) and (v0, v2, v3)

    TVector3d centroid, normal;
    double area;
    ComputeQuadGeometry(vertices, centroid, normal, area);

    // Use Wilton formula for each sub-triangle with obs at quad centroid
    double I_tri1 = WiltonTriangleIntegral1OverR(
        centroid, vertices[0], vertices[1], vertices[2]);
    double I_tri2 = WiltonTriangleIntegral1OverR(
        centroid, vertices[0], vertices[2], vertices[3]);

    double I_total = I_tri1 + I_tri2;

    // For Laplace kernel: G(r) = 1/(4*pi*r)
    return Complex(I_total * RadConst::INV_FOUR_PI, 0.0);
}

} // namespace radia
