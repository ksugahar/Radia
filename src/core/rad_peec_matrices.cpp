/*
 * rad_peec_matrices.cpp
 *
 * PEEC Matrix Construction Implementation
 *
 * Part of Radia project
 */

#include "rad_peec_matrices.h"
#include <cmath>
#include <algorithm>
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace radia {

// ============================================================================
// PEECMatrixBuilder
// ============================================================================

PEECMatrixBuilder::PEECMatrixBuilder() {}

PEECMatrixBuilder::~PEECMatrixBuilder() {}

void PEECMatrixBuilder::AddSegment(const PEECSegment& segment) {
    segments_.push_back(segment);
}

void PEECMatrixBuilder::AddSegments(const std::vector<PEECSegment>& segments) {
    segments_.insert(segments_.end(), segments.begin(), segments.end());
}

void PEECMatrixBuilder::AddNode(const PEECNode& node) {
    nodes_.push_back(node);
}

void PEECMatrixBuilder::AddNodes(const std::vector<PEECNode>& nodes) {
    nodes_.insert(nodes_.end(), nodes.begin(), nodes.end());
}

void PEECMatrixBuilder::Clear() {
    segments_.clear();
    nodes_.clear();
}

void PEECMatrixBuilder::AutoGenerateNodes() {
    nodes_.clear();

    // Use set to track unique positions (rounded to avoid floating point issues)
    struct Vec3Compare {
        bool operator()(const TVector3d& a, const TVector3d& b) const {
            const double tol = 1e-10;
            if (std::abs(a.x - b.x) > tol) return a.x < b.x;
            if (std::abs(a.y - b.y) > tol) return a.y < b.y;
            return a.z < b.z - tol;
        }
    };

    std::set<TVector3d, Vec3Compare> uniquePositions;
    std::vector<std::pair<TVector3d, double>> positionAreaPairs;

    for (const auto& seg : segments_) {
        double halfLen = seg.length / 2.0;
        TVector3d p1, p2;

        p1.x = seg.center.x - halfLen * seg.direction.x;
        p1.y = seg.center.y - halfLen * seg.direction.y;
        p1.z = seg.center.z - halfLen * seg.direction.z;

        p2.x = seg.center.x + halfLen * seg.direction.x;
        p2.y = seg.center.y + halfLen * seg.direction.y;
        p2.z = seg.center.z + halfLen * seg.direction.z;

        double area = seg.area();

        if (uniquePositions.find(p1) == uniquePositions.end()) {
            uniquePositions.insert(p1);
            positionAreaPairs.push_back({p1, area});
        }
        if (uniquePositions.find(p2) == uniquePositions.end()) {
            uniquePositions.insert(p2);
            positionAreaPairs.push_back({p2, area});
        }
    }

    // Create nodes
    for (const auto& pair : positionAreaPairs) {
        nodes_.push_back(PEECNode(pair.first, pair.second));
    }
}

PEECMatrices PEECMatrixBuilder::Build(bool includeStar) {
    PEECMatrices matrices;

    int n_loop = static_cast<int>(segments_.size());
    matrices.n_loop = n_loop;

    // Allocate and compute L
    matrices.L.resize(n_loop * n_loop, 0.0);
    ComputeL(matrices);

    // Compute R
    matrices.R.resize(n_loop, 0.0);
    ComputeR(matrices);

    if (includeStar) {
        // Auto-generate nodes if not provided
        if (nodes_.empty()) {
            AutoGenerateNodes();
        }

        int n_star = static_cast<int>(nodes_.size());
        matrices.n_star = n_star;

        // Allocate and compute P
        matrices.P.resize(n_star * n_star, 0.0);
        ComputeP(matrices);

        // Allocate and compute M_LS
        matrices.M_LS.resize(n_loop * n_star, 0.0);
        ComputeM_LS(matrices);
    } else {
        matrices.n_star = 0;
    }

    return matrices;
}

void PEECMatrixBuilder::ComputeL(PEECMatrices& matrices) {
    int n = static_cast<int>(segments_.size());

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < n; ++i) {
        // Self-inductance
        matrices.L[i * n + i] = SelfInductance(segments_[i]);

        // Mutual inductance (upper triangle)
        for (int j = i + 1; j < n; ++j) {
            double Lij = MutualInductance(segments_[i], segments_[j]);
            matrices.L[i * n + j] = Lij;
            matrices.L[j * n + i] = Lij;  // Symmetric
        }
    }
}

void PEECMatrixBuilder::ComputeP(PEECMatrices& matrices) {
    int n = static_cast<int>(nodes_.size());

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < n; ++i) {
        // Self-potential
        matrices.P[i * n + i] = SelfPotential(nodes_[i]);

        // Mutual potential (upper triangle)
        for (int j = i + 1; j < n; ++j) {
            double Pij = MutualPotential(nodes_[i], nodes_[j]);
            matrices.P[i * n + j] = Pij;
            matrices.P[j * n + i] = Pij;  // Symmetric
        }
    }
}

void PEECMatrixBuilder::ComputeR(PEECMatrices& matrices) {
    int n = static_cast<int>(segments_.size());

    for (int i = 0; i < n; ++i) {
        matrices.R[i] = segments_[i].resistance();
    }
}

void PEECMatrixBuilder::ComputeM_LS(PEECMatrices& matrices) {
    int n_loop = static_cast<int>(segments_.size());
    int n_star = static_cast<int>(nodes_.size());

    double coeff = PEEC_MU_0 * PEEC_INV_FOUR_PI;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < n_loop; ++i) {
        const PEECSegment& seg = segments_[i];
        double l_i = seg.length;

        for (int j = 0; j < n_star; ++j) {
            const PEECNode& node = nodes_[j];

            // Distance from segment center to node
            double dx = seg.center.x - node.position.x;
            double dy = seg.center.y - node.position.y;
            double dz = seg.center.z - node.position.z;
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (r > 1e-15) {
                matrices.M_LS[i * n_star + j] = coeff * l_i / r;
            }
        }
    }
}

double PEECMatrixBuilder::SelfInductance(const PEECSegment& seg) const {
    // GMD approximation for self-inductance
    // L = (mu_0 / 2*pi) * l * (ln(2*l/GMD) - 1)
    // GMD for rectangular cross-section: GMD ~ 0.2235 * (w + h)

    double l = seg.length;
    double gmd = 0.2235 * (seg.width + seg.height);

    if (gmd < 1e-15) gmd = 1e-6;  // Minimum GMD

    if (l > gmd) {
        return (PEEC_MU_0 / (2.0 * RadConst::PI)) * l * (std::log(2.0 * l / gmd) - 1.0);
    } else {
        // Short segment approximation
        return (PEEC_MU_0 / (2.0 * RadConst::PI)) * l * 0.5;
    }
}

double PEECMatrixBuilder::MutualInductance(const PEECSegment& seg_i,
                                            const PEECSegment& seg_j) const {
    // Neumann formula approximation (point matching)
    // L_ij = (mu_0 / 4*pi) * (d_i . d_j) * l_i * l_j / r_ij

    double dx = seg_i.center.x - seg_j.center.x;
    double dy = seg_i.center.y - seg_j.center.y;
    double dz = seg_i.center.z - seg_j.center.z;
    double r = std::sqrt(dx*dx + dy*dy + dz*dz);

    if (r < 1e-15) return 0.0;

    // Direction dot product
    double dot = seg_i.direction.x * seg_j.direction.x +
                 seg_i.direction.y * seg_j.direction.y +
                 seg_i.direction.z * seg_j.direction.z;

    return (PEEC_MU_0 * PEEC_INV_FOUR_PI) * dot * seg_i.length * seg_j.length / r;
}

double PEECMatrixBuilder::SelfPotential(const PEECNode& node) const {
    // Self-potential using disk approximation
    // P_ii = 1 / (4 * eps_0 * sqrt(pi * A))
    // where A is the associated area

    double coeff = 1.0 / (4.0 * RadConst::PI * PEEC_EPS_0);

    if (node.area > 0) {
        double equiv_radius = std::sqrt(node.area / RadConst::PI);
        return coeff / equiv_radius;
    } else {
        // Default small radius
        return coeff / 1e-6;
    }
}

double PEECMatrixBuilder::MutualPotential(const PEECNode& node_i,
                                           const PEECNode& node_j) const {
    // Mutual potential coefficient
    // P_ij = 1 / (4*pi*eps_0 * r_ij)

    double dx = node_i.position.x - node_j.position.x;
    double dy = node_i.position.y - node_j.position.y;
    double dz = node_i.position.z - node_j.position.z;
    double r = std::sqrt(dx*dx + dy*dy + dz*dz);

    if (r < 1e-15) return 0.0;

    return 1.0 / (4.0 * RadConst::PI * PEEC_EPS_0 * r);
}

// ============================================================================
// PEECSolver
// ============================================================================

PEECSolver::PEECSolver()
    : frequency_(0), omega_(0), hasSurfaceImpedance_(false) {}

PEECSolver::~PEECSolver() {}

void PEECSolver::SetMatrices(const PEECMatrices& matrices) {
    matrices_ = matrices;
}

void PEECSolver::SetFrequency(double freq_hz) {
    frequency_ = freq_hz;
    omega_ = 2.0 * RadConst::PI * freq_hz;
}

void PEECSolver::SetSurfaceImpedance(const std::vector<std::complex<double>>& Zs_diag) {
    Zs_ = Zs_diag;
    hasSurfaceImpedance_ = !Zs_.empty();
}

void PEECSolver::BuildImpedanceMatrix(std::vector<std::complex<double>>& Z) {
    int n_loop = matrices_.n_loop;
    int n_star = matrices_.n_star;
    int n_total = n_loop + n_star;

    Z.resize(n_total * n_total, std::complex<double>(0, 0));

    // Z_LL = R + jw*L
    for (int i = 0; i < n_loop; ++i) {
        for (int j = 0; j < n_loop; ++j) {
            double Rij = (i == j) ? matrices_.R[i] : 0.0;
            double Lij = matrices_.L_at(i, j);
            Z[i * n_total + j] = std::complex<double>(Rij, omega_ * Lij);
        }

        // Add surface impedance (diagonal)
        if (hasSurfaceImpedance_ && i < static_cast<int>(Zs_.size())) {
            Z[i * n_total + i] += Zs_[i];
        }
    }

    if (n_star > 0) {
        // Z_SS = P / (jw)
        std::complex<double> inv_jw;
        if (std::abs(omega_) > 1e-15) {
            inv_jw = std::complex<double>(0, -1.0 / omega_);
        } else {
            // DC: very large impedance (open circuit)
            inv_jw = std::complex<double>(1e15, 0);
        }

        for (int i = 0; i < n_star; ++i) {
            for (int j = 0; j < n_star; ++j) {
                double Pij = matrices_.P_at(i, j);
                Z[(n_loop + i) * n_total + (n_loop + j)] = Pij * inv_jw;
            }
        }

        // Z_LS = jw * M_LS
        std::complex<double> jw(0, omega_);
        for (int i = 0; i < n_loop; ++i) {
            for (int j = 0; j < n_star; ++j) {
                double Mij = matrices_.M_LS_at(i, j);
                std::complex<double> Zij = jw * Mij;
                Z[i * n_total + (n_loop + j)] = Zij;
                Z[(n_loop + j) * n_total + i] = Zij;  // Z_SL = Z_LS^T
            }
        }
    }
}

std::complex<double> PEECSolver::ComputePortImpedance(const std::vector<double>& portVector) {
    int n_loop = matrices_.n_loop;
    int n_star = matrices_.n_star;
    int n_total = n_loop + n_star;

    // Build impedance matrix
    std::vector<std::complex<double>> Z;
    BuildImpedanceMatrix(Z);

    // Create voltage vector
    std::vector<std::complex<double>> V(n_total, std::complex<double>(0, 0));
    for (int i = 0; i < n_loop && i < static_cast<int>(portVector.size()); ++i) {
        V[i] = std::complex<double>(portVector[i], 0);
    }

    // Solve Z * I = V using LU decomposition (simple implementation)
    // For production, use LAPACK
    std::vector<std::complex<double>> I(n_total);
    Solve(V, I);

    // Port impedance: Z_port = v^T * I_loop
    std::complex<double> Z_port(0, 0);
    for (int i = 0; i < n_loop && i < static_cast<int>(portVector.size()); ++i) {
        Z_port += portVector[i] * I[i];
    }

    return Z_port;
}

void PEECSolver::Solve(const std::vector<std::complex<double>>& V,
                        std::vector<std::complex<double>>& I) {
    int n = static_cast<int>(V.size());
    I.resize(n);

    // Build impedance matrix
    std::vector<std::complex<double>> Z;
    BuildImpedanceMatrix(Z);

    // Simple Gaussian elimination with partial pivoting
    // For production, use LAPACK zgesv

    std::vector<std::complex<double>> A = Z;  // Copy for modification
    std::vector<std::complex<double>> b = V;

    // Forward elimination
    for (int k = 0; k < n - 1; ++k) {
        // Find pivot
        int maxIdx = k;
        double maxVal = std::abs(A[k * n + k]);
        for (int i = k + 1; i < n; ++i) {
            double val = std::abs(A[i * n + k]);
            if (val > maxVal) {
                maxVal = val;
                maxIdx = i;
            }
        }

        // Swap rows
        if (maxIdx != k) {
            for (int j = 0; j < n; ++j) {
                std::swap(A[k * n + j], A[maxIdx * n + j]);
            }
            std::swap(b[k], b[maxIdx]);
        }

        // Eliminate
        std::complex<double> pivot = A[k * n + k];
        if (std::abs(pivot) < 1e-30) continue;

        for (int i = k + 1; i < n; ++i) {
            std::complex<double> factor = A[i * n + k] / pivot;
            for (int j = k; j < n; ++j) {
                A[i * n + j] -= factor * A[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }

    // Back substitution
    for (int i = n - 1; i >= 0; --i) {
        std::complex<double> sum = b[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= A[i * n + j] * I[j];
        }
        std::complex<double> diag = A[i * n + i];
        I[i] = (std::abs(diag) > 1e-30) ? sum / diag : std::complex<double>(0, 0);
    }
}

// ============================================================================
// Utility Functions
// ============================================================================

std::vector<PEECSegment> CreateWireSegments(
    const TVector3d& start,
    const TVector3d& end,
    double width,
    double height,
    int n_segments,
    double sigma)
{
    std::vector<PEECSegment> segments;

    TVector3d dir;
    dir.x = end.x - start.x;
    dir.y = end.y - start.y;
    dir.z = end.z - start.z;

    double totalLen = std::sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
    if (totalLen < 1e-15 || n_segments < 1) return segments;

    // Normalize direction
    dir.x /= totalLen;
    dir.y /= totalLen;
    dir.z /= totalLen;

    double segLen = totalLen / n_segments;

    for (int i = 0; i < n_segments; ++i) {
        double t = (i + 0.5) / n_segments;
        TVector3d center;
        center.x = start.x + t * (end.x - start.x);
        center.y = start.y + t * (end.y - start.y);
        center.z = start.z + t * (end.z - start.z);

        segments.push_back(PEECSegment(center, dir, segLen, width, height, sigma));
    }

    return segments;
}

std::vector<PEECSegment> CreateLoopSegments(
    const TVector3d& center,
    double radius,
    const TVector3d& normal,
    double width,
    double height,
    int n_segments,
    double sigma)
{
    std::vector<PEECSegment> segments;

    if (radius < 1e-15 || n_segments < 3) return segments;

    // Create local coordinate system
    // normal is z-axis, find x and y axes
    TVector3d n = normal;
    double nLen = std::sqrt(n.x*n.x + n.y*n.y + n.z*n.z);
    if (nLen < 1e-15) return segments;
    n.x /= nLen; n.y /= nLen; n.z /= nLen;

    // Find a vector not parallel to normal
    TVector3d temp;
    if (std::abs(n.x) < 0.9) {
        temp = TVector3d(1, 0, 0);
    } else {
        temp = TVector3d(0, 1, 0);
    }

    // x-axis = temp x n (normalized)
    TVector3d ex;
    ex.x = temp.y * n.z - temp.z * n.y;
    ex.y = temp.z * n.x - temp.x * n.z;
    ex.z = temp.x * n.y - temp.y * n.x;
    double exLen = std::sqrt(ex.x*ex.x + ex.y*ex.y + ex.z*ex.z);
    ex.x /= exLen; ex.y /= exLen; ex.z /= exLen;

    // y-axis = n x ex
    TVector3d ey;
    ey.x = n.y * ex.z - n.z * ex.y;
    ey.y = n.z * ex.x - n.x * ex.z;
    ey.z = n.x * ex.y - n.y * ex.x;

    double dtheta = 2.0 * RadConst::PI / n_segments;
    double arcLen = radius * dtheta;

    for (int i = 0; i < n_segments; ++i) {
        double theta = (i + 0.5) * dtheta;

        // Position on circle
        double cosT = std::cos(theta);
        double sinT = std::sin(theta);

        TVector3d pos;
        pos.x = center.x + radius * (cosT * ex.x + sinT * ey.x);
        pos.y = center.y + radius * (cosT * ex.y + sinT * ey.y);
        pos.z = center.z + radius * (cosT * ex.z + sinT * ey.z);

        // Tangent direction (perpendicular to radial)
        TVector3d dir;
        dir.x = -sinT * ex.x + cosT * ey.x;
        dir.y = -sinT * ex.y + cosT * ey.y;
        dir.z = -sinT * ex.z + cosT * ey.z;

        segments.push_back(PEECSegment(pos, dir, arcLen, width, height, sigma));
    }

    return segments;
}

} // namespace radia
