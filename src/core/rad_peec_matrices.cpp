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
// PEECPanel - Geometry Computation
// ============================================================================

void PEECPanel::ComputeGeometry() {
    if (vertices.size() < 3) {
        area = 0;
        return;
    }

    // Compute centroid
    center = TVector3d(0, 0, 0);
    for (const auto& v : vertices) {
        center.x += v.x;
        center.y += v.y;
        center.z += v.z;
    }
    center.x /= vertices.size();
    center.y /= vertices.size();
    center.z /= vertices.size();

    if (type == Triangle) {
        // Triangle: Area = 0.5 * |v1-v0 × v2-v0|
        TVector3d v0 = vertices[0];
        TVector3d v1 = vertices[1];
        TVector3d v2 = vertices[2];

        TVector3d edge1 = v1 - v0;
        TVector3d edge2 = v2 - v0;

        // Cross product
        normal.x = edge1.y * edge2.z - edge1.z * edge2.y;
        normal.y = edge1.z * edge2.x - edge1.x * edge2.z;
        normal.z = edge1.x * edge2.y - edge1.y * edge2.x;

        area = 0.5 * std::sqrt(normal.x * normal.x +
                              normal.y * normal.y +
                              normal.z * normal.z);

        // Normalize normal vector
        double norm_len = std::sqrt(normal.x * normal.x +
                                   normal.y * normal.y +
                                   normal.z * normal.z);
        if (norm_len > 0) {
            normal.x /= norm_len;
            normal.y /= norm_len;
            normal.z /= norm_len;
        }

    } else {  // Quadrilateral
        // Split into two triangles and sum areas
        TVector3d v0 = vertices[0];
        TVector3d v1 = vertices[1];
        TVector3d v2 = vertices[2];
        TVector3d v3 = vertices[3];

        // Triangle 1: v0-v1-v2
        TVector3d edge1 = v1 - v0;
        TVector3d edge2 = v2 - v0;
        TVector3d cross1;
        cross1.x = edge1.y * edge2.z - edge1.z * edge2.y;
        cross1.y = edge1.z * edge2.x - edge1.x * edge2.z;
        cross1.z = edge1.x * edge2.y - edge1.y * edge2.x;
        double area1 = 0.5 * std::sqrt(cross1.x * cross1.x +
                                      cross1.y * cross1.y +
                                      cross1.z * cross1.z);

        // Triangle 2: v0-v2-v3
        TVector3d edge3 = v3 - v0;
        TVector3d cross2;
        cross2.x = edge2.y * edge3.z - edge2.z * edge3.y;
        cross2.y = edge2.z * edge3.x - edge2.x * edge3.z;
        cross2.z = edge2.x * edge3.y - edge2.y * edge3.x;
        double area2 = 0.5 * std::sqrt(cross2.x * cross2.x +
                                      cross2.y * cross2.y +
                                      cross2.z * cross2.z);

        area = area1 + area2;

        // Normal: average of two triangle normals
        normal.x = cross1.x + cross2.x;
        normal.y = cross1.y + cross2.y;
        normal.z = cross1.z + cross2.z;

        double norm_len = std::sqrt(normal.x * normal.x +
                                   normal.y * normal.y +
                                   normal.z * normal.z);
        if (norm_len > 0) {
            normal.x /= norm_len;
            normal.y /= norm_len;
            normal.z /= norm_len;
        }
    }
}

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

void PEECMatrixBuilder::AddPanel(const PEECPanel& panel) {
    panels_.push_back(panel);
}

void PEECMatrixBuilder::AddPanels(const std::vector<PEECPanel>& panels) {
    panels_.insert(panels_.end(), panels.begin(), panels.end());
}

void PEECMatrixBuilder::Clear() {
    segments_.clear();
    nodes_.clear();
    panels_.clear();
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
        // Determine number of Star elements (panels or nodes)
        int n_star = 0;
        if (!panels_.empty()) {
            // Use panels for Star elements (true 2D surface integration)
            n_star = static_cast<int>(panels_.size());
        } else {
            // Auto-generate nodes if not provided (point approximation)
            if (nodes_.empty()) {
                AutoGenerateNodes();
            }
            n_star = static_cast<int>(nodes_.size());
        }
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
    // Use panels if available, otherwise use nodes (point approximation)
    if (!panels_.empty()) {
        // Panel-based P matrix (analytical surface integration)
        int n = static_cast<int>(panels_.size());

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (int i = 0; i < n; ++i) {
            const PEECPanel& panel_i = panels_[i];

            // Self-potential (analytical integration)
            if (panel_i.type == PEECPanel::Triangle) {
                matrices.P[i * n + i] = SelfPotentialPanelTriangle(panel_i);
            } else if (panel_i.type == PEECPanel::Quadrilateral) {
                matrices.P[i * n + i] = SelfPotentialPanelQuad(panel_i);
            } else {
                // Fallback to zero for unknown panel types
                matrices.P[i * n + i] = 0.0;
            }

            // Mutual potential (upper triangle)
            for (int j = i + 1; j < n; ++j) {
                const PEECPanel& panel_j = panels_[j];
                double Pij = MutualPotentialPanelTriangle(panel_i, panel_j);
                matrices.P[i * n + j] = Pij;
                matrices.P[j * n + i] = Pij;  // Symmetric
            }
        }
    } else {
        // Node-based P matrix (point approximation - legacy)
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
}

void PEECMatrixBuilder::ComputeR(PEECMatrices& matrices) {
    // DC resistance only. Frequency-dependent Z_s is computed in Python
    // (scipy.special.jv for Bessel SIBC, Dowell formula, or ESIM).
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
    // Dispatch to rectangular or circular formula based on cross-section type
    if (seg.cross_section_type == CrossSectionType::CIRCULAR) {
        return SelfInductanceCircular(seg);
    } else {
        return SelfInductanceRectangular(seg);
    }
}

double PEECMatrixBuilder::SelfInductanceRectangular(const PEECSegment& seg) const {
    // Grover formula for rectangular cross-section (EXACT, no GMD approximation)
    // Reference: F. W. Grover, "Inductance Calculations", Dover, 1946
    //
    // L = (mu_0/2pi) * l * [ln(2*l/sqrt(w^2+h^2)) + 0.25 + (w^2+h^2)/(12*l^2)]
    //
    // This is the EXACT formula for a straight rectangular conductor segment
    // NO conversion to circular cross-section (FastImp approach)

    double l = seg.length;
    double w = seg.width;
    double h = seg.height;

    // Rectangular cross-section diagonal
    double d_rect = std::sqrt(w*w + h*h);

    // Minimum cross-section check
    if (d_rect < 1e-15) d_rect = 1e-6;

    if (l > d_rect) {
        // Grover formula (exact for rectangular cross-section)
        double term1 = std::log(2.0 * l / d_rect);
        double term2 = 0.25;
        double term3 = (w*w + h*h) / (12.0 * l*l);

        return (PEEC_MU_0 / (2.0 * RadConst::PI)) * l * (term1 + term2 + term3);
    } else {
        // Short segment approximation (l << cross-section)
        return (PEEC_MU_0 / (2.0 * RadConst::PI)) * l * 0.5;
    }
}

double PEECMatrixBuilder::SelfInductanceCircular(const PEECSegment& seg) const {
    // Grover formula for circular cross-section
    // Reference: F. W. Grover, "Inductance Calculations", Dover, 1946
    //
    // L = (mu_0/2pi) * l * [ln(2*l/r) - 0.75]
    //
    // Where:
    //   l = conductor length [m]
    //   r = wire radius [m]
    //   -0.75 = internal inductance constant for circular cross-section
    //
    // Note: This differs from rectangular formula by the constant term:
    //   Rectangular: +0.25
    //   Circular: -0.75
    //   Difference: 1.0 (accounts for ~17% difference in internal inductance)

    double l = seg.length;

    // Extract radius from cross-sectional area
    // For circular: area = pi * r^2, so r = sqrt(area / pi)
    double area = seg.area();
    double r = std::sqrt(area / RadConst::PI);

    // Minimum cross-section check
    if (r < 1e-15) r = 1e-6;

    if (l > 2.0 * r) {
        // Grover formula for circular cross-section
        double term1 = std::log(2.0 * l / r);
        double term2 = -0.75;  // Circular cross-section constant

        return (PEEC_MU_0 / (2.0 * RadConst::PI)) * l * (term1 + term2);
    } else {
        // Short segment approximation (l << diameter)
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
// Panel Analytical Integration (Wilton et al., 1984)
// ============================================================================

double PEECMatrixBuilder::SelfPotentialPanelTriangle(const PEECPanel& panel) const {
    // Analytical self-potential for triangular panel using Wilton formula
    // Reference: Wilton et al., IEEE TAP, vol. 32, no. 3, pp. 276-281, 1984
    //
    // Formula: P_self = (1/4πε₀) * Σ_edges [analytical_edge_integral]

    if (panel.type != PEECPanel::Triangle || panel.vertices.size() != 3) {
        return 0.0;  // Fallback to node approximation
    }

    const TVector3d& v0 = panel.vertices[0];
    const TVector3d& v1 = panel.vertices[1];
    const TVector3d& v2 = panel.vertices[2];

    double sum = 0.0;

    // Edge integration over 3 edges
    for (int edge_idx = 0; edge_idx < 3; ++edge_idx) {
        // Edge endpoints
        const TVector3d& p0 = panel.vertices[edge_idx];
        const TVector3d& p1 = panel.vertices[(edge_idx + 1) % 3];
        const TVector3d& p2 = panel.vertices[(edge_idx + 2) % 3];  // Opposite vertex

        // Edge vector
        TVector3d edge;
        edge.x = p1.x - p0.x;
        edge.y = p1.y - p0.y;
        edge.z = p1.z - p0.z;
        double l_edge = std::sqrt(edge.x * edge.x + edge.y * edge.y + edge.z * edge.z);

        if (l_edge < 1e-15) continue;  // Degenerate edge

        // Vector from p0 to opposite vertex
        TVector3d r0;
        r0.x = p2.x - p0.x;
        r0.y = p2.y - p0.y;
        r0.z = p2.z - p0.z;

        // Vector from p1 to opposite vertex
        TVector3d r1;
        r1.x = p2.x - p1.x;
        r1.y = p2.y - p1.y;
        r1.z = p2.z - p1.z;

        // Distances
        double R0 = std::sqrt(r0.x * r0.x + r0.y * r0.y + r0.z * r0.z);
        double R1 = std::sqrt(r1.x * r1.x + r1.y * r1.y + r1.z * r1.z);

        if (R0 < 1e-15 || R1 < 1e-15) continue;

        // Cross product: edge × r0 (to get height vector)
        TVector3d cross;
        cross.x = edge.y * r0.z - edge.z * r0.y;
        cross.y = edge.z * r0.x - edge.x * r0.z;
        cross.z = edge.x * r0.y - edge.y * r0.x;
        double h = std::sqrt(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z) / l_edge;

        if (h < 1e-15) continue;  // Degenerate triangle

        // Wilton analytical formula for edge contribution
        double arg = (R0 + R1 + l_edge) / (R0 + R1 - l_edge);
        if (arg > 0 && arg < 1e15) {  // Valid logarithm argument
            double ln_term = std::log(arg);
            sum += l_edge * ln_term;
        }
    }

    return sum / (4.0 * RadConst::PI * PEEC_EPS_0);
}

double PEECMatrixBuilder::SelfPotentialPanelQuad(const PEECPanel& panel) const {
    // Quadrilateral self-potential by splitting into 2 triangles
    //
    // Split quad (v0, v1, v2, v3) into:
    //   Triangle 1: (v0, v1, v2)
    //   Triangle 2: (v0, v2, v3)
    //
    // Sum contributions from both triangles

    if (panel.type != PEECPanel::Quadrilateral || panel.vertices.size() != 4) {
        return 0.0;  // Invalid quad
    }

    // Create temporary triangle panels
    PEECPanel tri1, tri2;
    tri1.type = PEECPanel::Triangle;
    tri2.type = PEECPanel::Triangle;

    // Triangle 1: v0-v1-v2
    tri1.vertices.push_back(panel.vertices[0]);
    tri1.vertices.push_back(panel.vertices[1]);
    tri1.vertices.push_back(panel.vertices[2]);
    tri1.ComputeGeometry();

    // Triangle 2: v0-v2-v3
    tri2.vertices.push_back(panel.vertices[0]);
    tri2.vertices.push_back(panel.vertices[2]);
    tri2.vertices.push_back(panel.vertices[3]);
    tri2.ComputeGeometry();

    // Compute self-potential for each triangle and sum
    double P_tri1 = SelfPotentialPanelTriangle(tri1);
    double P_tri2 = SelfPotentialPanelTriangle(tri2);

    // Average the two contributions (both triangles are part of the same quad)
    // Actually, we need to be more careful here. The self-potential of a quad
    // is NOT simply the sum of the triangle self-potentials, because the triangles
    // are not independent - they share edges.
    //
    // For now, we approximate by averaging. This is a known limitation.
    // A more accurate approach would require proper quad integration formulas.

    return 0.5 * (P_tri1 + P_tri2);
}

double PEECMatrixBuilder::MutualPotentialPanelTriangle(const PEECPanel& panel_i,
                                                        const PEECPanel& panel_j) const {
    // Mutual potential between two triangular panels
    //
    // Strategy:
    //   - Far-field (distance > 3 * panel_size): Centroid approximation
    //   - Near-field (distance < 3 * panel_size): TODO: Hess-Smith edge integration

    if (panel_i.type != PEECPanel::Triangle || panel_j.type != PEECPanel::Triangle) {
        return 0.0;  // Fallback
    }

    // Distance between panel centroids
    double dx = panel_i.center.x - panel_j.center.x;
    double dy = panel_i.center.y - panel_j.center.y;
    double dz = panel_i.center.z - panel_j.center.z;
    double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

    // Characteristic panel size (sqrt of area)
    double char_size_i = std::sqrt(panel_i.area);
    double char_size_j = std::sqrt(panel_j.area);
    double char_size = std::max(char_size_i, char_size_j);

    if (dist > 3.0 * char_size) {
        // Far-field: Centroid approximation (monopole-monopole)
        // P_ij ≈ (Area_i * Area_j) / (4πε₀ * r)
        return (panel_i.area * panel_j.area) / (4.0 * RadConst::PI * PEEC_EPS_0 * dist);
    } else {
        // Near-field: Use Gauss quadrature integration
        // 3-point Gauss quadrature for triangles (barycentric coordinates)
        //
        // Reference: Gauss integration on triangular domains

        // 3-point Gauss rule for triangles
        const double w = 1.0 / 6.0;  // Each point has equal weight (total = 1/2 for triangle)

        // Barycentric coordinates for 3-point rule
        const double xi[3][3] = {
            {0.5, 0.5, 0.0},   // Midpoint of edge 0-1
            {0.0, 0.5, 0.5},   // Midpoint of edge 1-2
            {0.5, 0.0, 0.5}    // Midpoint of edge 2-0
        };

        double sum = 0.0;

        // Get vertices for both panels
        const TVector3d& v0_i = panel_i.vertices[0];
        const TVector3d& v1_i = panel_i.vertices[1];
        const TVector3d& v2_i = panel_i.vertices[2];

        const TVector3d& v0_j = panel_j.vertices[0];
        const TVector3d& v1_j = panel_j.vertices[1];
        const TVector3d& v2_j = panel_j.vertices[2];

        // Double loop over quadrature points
        for (int qi = 0; qi < 3; ++qi) {
            // Compute physical point on panel_i using barycentric coordinates
            TVector3d point_i;
            point_i.x = xi[qi][0] * v0_i.x + xi[qi][1] * v1_i.x + xi[qi][2] * v2_i.x;
            point_i.y = xi[qi][0] * v0_i.y + xi[qi][1] * v1_i.y + xi[qi][2] * v2_i.y;
            point_i.z = xi[qi][0] * v0_i.z + xi[qi][1] * v1_i.z + xi[qi][2] * v2_i.z;

            for (int qj = 0; qj < 3; ++qj) {
                // Compute physical point on panel_j
                TVector3d point_j;
                point_j.x = xi[qj][0] * v0_j.x + xi[qj][1] * v1_j.x + xi[qj][2] * v2_j.x;
                point_j.y = xi[qj][0] * v0_j.y + xi[qj][1] * v1_j.y + xi[qj][2] * v2_j.y;
                point_j.z = xi[qj][0] * v0_j.z + xi[qj][1] * v1_j.z + xi[qj][2] * v2_j.z;

                // Distance between quadrature points
                double dx_q = point_i.x - point_j.x;
                double dy_q = point_i.y - point_j.y;
                double dz_q = point_i.z - point_j.z;
                double R = std::sqrt(dx_q*dx_q + dy_q*dy_q + dz_q*dz_q);

                if (R > 1e-10) {  // Avoid singularity
                    sum += w * w / R;
                }
            }
        }

        // Multiply by panel areas (Jacobian of the barycentric transformation)
        return (panel_i.area * panel_j.area * sum) / (4.0 * RadConst::PI * PEEC_EPS_0);
    }
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
