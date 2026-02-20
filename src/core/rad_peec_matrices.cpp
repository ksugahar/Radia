/*
 * rad_peec_matrices.cpp
 *
 * PEEC Matrix Construction Implementation
 *
 * Part of Radia project
 */

#include "rad_peec_matrices.h"
#include "rad_bicgstab.h"
#include <cmath>
#include <algorithm>
#include <set>
#include <queue>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef HAVE_LAPACK
#include "mkl_cblas.h"
#include "mkl_lapack.h"
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

PEECMatrixBuilder::PEECMatrixBuilder() : nextPortId_(0) {}

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
    ports_.clear();
    nextPortId_ = 0;
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

// ============================================================================
// Topology-aware methods (node-segment model)
// ============================================================================

int PEECMatrixBuilder::AddNodeAt(const TVector3d& position, double area) {
    int id = static_cast<int>(nodes_.size());
    nodes_.push_back(PEECNode(position, area));
    return id;
}

void PEECMatrixBuilder::AddConnectedSegment(int node_from, int node_to,
                                             double width, double height,
                                             double sigma,
                                             CrossSectionType type,
                                             int nwinc, int nhinc) {
    int n_nodes = static_cast<int>(nodes_.size());
    if (node_from < 0 || node_from >= n_nodes ||
        node_to < 0 || node_to >= n_nodes) {
        return;  // Invalid node IDs
    }

    const TVector3d& p1 = nodes_[node_from].position;
    const TVector3d& p2 = nodes_[node_to].position;

    // Compute center
    TVector3d center;
    center.x = 0.5 * (p1.x + p2.x);
    center.y = 0.5 * (p1.y + p2.y);
    center.z = 0.5 * (p1.z + p2.z);

    // Compute direction and length
    TVector3d dir;
    dir.x = p2.x - p1.x;
    dir.y = p2.y - p1.y;
    dir.z = p2.z - p1.z;
    double len = std::sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);

    if (len < 1e-15) return;  // Zero-length segment

    dir.x /= len;
    dir.y /= len;
    dir.z /= len;

    PEECSegment seg(center, dir, len, width, height, sigma, type);
    seg.node_from = node_from;
    seg.node_to = node_to;
    seg.nwinc = (nwinc > 0) ? nwinc : 1;
    seg.nhinc = (nhinc > 0) ? nhinc : 1;
    seg.parent_segment = -1;
    segments_.push_back(seg);
}

void PEECMatrixBuilder::ExpandFilaments() {
    // Expand segments with nwinc*nhinc > 1 into sub-filaments
    // Sub-filaments are parallel (same node_from, node_to) with smaller cross-sections
    // and offset centers in the local cross-section coordinate system.

    std::vector<PEECSegment> expanded;

    for (int s = 0; s < static_cast<int>(segments_.size()); ++s) {
        const PEECSegment& seg = segments_[s];

        int nw = seg.nwinc;
        int nh = seg.nhinc;

        if (nw <= 1 && nh <= 1) {
            // No subdivision needed - keep original
            expanded.push_back(seg);
            continue;
        }

        // Build local coordinate system perpendicular to segment direction
        // dir = segment direction (unit vector)
        // e_w = width direction (perpendicular to dir)
        // e_h = height direction (perpendicular to dir and e_w)
        TVector3d dir = seg.direction;

        // Choose a reference vector not parallel to dir
        TVector3d ref;
        if (std::abs(dir.x) < 0.9) {
            ref = TVector3d(1, 0, 0);
        } else {
            ref = TVector3d(0, 1, 0);
        }

        // e_w = ref x dir (normalized)
        TVector3d e_w;
        e_w.x = ref.y * dir.z - ref.z * dir.y;
        e_w.y = ref.z * dir.x - ref.x * dir.z;
        e_w.z = ref.x * dir.y - ref.y * dir.x;
        double e_w_len = std::sqrt(e_w.x*e_w.x + e_w.y*e_w.y + e_w.z*e_w.z);
        if (e_w_len < 1e-15) continue;
        e_w.x /= e_w_len;
        e_w.y /= e_w_len;
        e_w.z /= e_w_len;

        // e_h = dir x e_w (already unit vector)
        TVector3d e_h;
        e_h.x = dir.y * e_w.z - dir.z * e_w.y;
        e_h.y = dir.z * e_w.x - dir.x * e_w.z;
        e_h.z = dir.x * e_w.y - dir.y * e_w.x;

        // Sub-filament dimensions
        double sub_w = seg.width / nw;
        double sub_h = seg.height / nh;

        // Create sub-filaments
        for (int iw = 0; iw < nw; ++iw) {
            for (int ih = 0; ih < nh; ++ih) {
                // Offset from parent center in local coordinates
                // Center of sub-filament (iw, ih) relative to parent center
                double offset_w = sub_w * (iw - (nw - 1) * 0.5);
                double offset_h = sub_h * (ih - (nh - 1) * 0.5);

                TVector3d sub_center;
                sub_center.x = seg.center.x + offset_w * e_w.x + offset_h * e_h.x;
                sub_center.y = seg.center.y + offset_w * e_w.y + offset_h * e_h.y;
                sub_center.z = seg.center.z + offset_w * e_w.z + offset_h * e_h.z;

                PEECSegment sub(sub_center, seg.direction, seg.length,
                                sub_w, sub_h, seg.sigma, seg.cross_section_type);
                sub.node_from = seg.node_from;
                sub.node_to = seg.node_to;
                sub.nwinc = 1;
                sub.nhinc = 1;
                sub.parent_segment = s;  // Index in original segment list

                expanded.push_back(sub);
            }
        }
    }

    segments_ = std::move(expanded);
}

int PEECMatrixBuilder::AddPort(int node_positive, int node_negative) {
    int id = nextPortId_++;
    ports_.push_back(PEECPort(node_positive, node_negative, id));
    return id;
}

void PEECMatrixBuilder::BuildIncidenceMatrix(PEECMatrices& matrices) {
    // Identify junction nodes: internal nodes that are NOT port terminals
    // A junction node has KCL constraint: sum of currents = 0
    int n_nodes = static_cast<int>(nodes_.size());
    int n_filaments = static_cast<int>(segments_.size());

    // Mark port terminal nodes
    std::vector<bool> is_port_terminal(n_nodes, false);
    for (const auto& port : ports_) {
        if (port.node_positive >= 0 && port.node_positive < n_nodes)
            is_port_terminal[port.node_positive] = true;
        if (port.node_negative >= 0 && port.node_negative < n_nodes)
            is_port_terminal[port.node_negative] = true;
    }

    // Count connections per node
    std::vector<int> connection_count(n_nodes, 0);
    for (const auto& seg : segments_) {
        if (seg.node_from >= 0) connection_count[seg.node_from]++;
        if (seg.node_to >= 0) connection_count[seg.node_to]++;
    }

    // Identify junction nodes: connected to 2+ segments and NOT a port terminal
    std::vector<int> node_to_junction(n_nodes, -1);  // maps node ID -> junction index
    int n_junction = 0;
    for (int i = 0; i < n_nodes; ++i) {
        if (!is_port_terminal[i] && connection_count[i] >= 2) {
            node_to_junction[i] = n_junction++;
        }
    }

    matrices.n_junction = n_junction;

    if (n_junction == 0) {
        // No junctions - all segments in series or simple topology
        matrices.incidence_indptr.assign(1, 0);
        matrices.incidence_indices.clear();
        matrices.incidence_data.clear();
        matrices.ports = ports_;
        return;
    }

    // Build CSR incidence matrix A (n_junction x n_filaments)
    // First pass: count non-zeros per row
    std::vector<int> row_nnz(n_junction, 0);
    for (int f = 0; f < n_filaments; ++f) {
        const auto& seg = segments_[f];
        if (seg.node_from >= 0 && node_to_junction[seg.node_from] >= 0)
            row_nnz[node_to_junction[seg.node_from]]++;
        if (seg.node_to >= 0 && node_to_junction[seg.node_to] >= 0)
            row_nnz[node_to_junction[seg.node_to]]++;
    }

    // Build row pointers
    matrices.incidence_indptr.resize(n_junction + 1);
    matrices.incidence_indptr[0] = 0;
    for (int i = 0; i < n_junction; ++i) {
        matrices.incidence_indptr[i + 1] = matrices.incidence_indptr[i] + row_nnz[i];
    }

    int total_nnz = matrices.incidence_indptr[n_junction];
    matrices.incidence_indices.resize(total_nnz);
    matrices.incidence_data.resize(total_nnz);

    // Second pass: fill entries
    std::vector<int> row_pos(n_junction, 0);  // current insert position per row
    for (int i = 0; i < n_junction; ++i) {
        row_pos[i] = matrices.incidence_indptr[i];
    }

    for (int f = 0; f < n_filaments; ++f) {
        const auto& seg = segments_[f];

        // Filament leaves node_from: A[junction, filament] = +1
        if (seg.node_from >= 0) {
            int j = node_to_junction[seg.node_from];
            if (j >= 0) {
                int pos = row_pos[j]++;
                matrices.incidence_indices[pos] = f;
                matrices.incidence_data[pos] = +1.0;
            }
        }

        // Filament enters node_to: A[junction, filament] = -1
        if (seg.node_to >= 0) {
            int j = node_to_junction[seg.node_to];
            if (j >= 0) {
                int pos = row_pos[j]++;
                matrices.incidence_indices[pos] = f;
                matrices.incidence_data[pos] = -1.0;
            }
        }
    }

    // Copy port definitions
    matrices.ports = ports_;
}

PEECMatrices PEECMatrixBuilder::Build(bool includeStar) {
    // Expand multi-filament segments before computing matrices
    ExpandFilaments();

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
            // In topology mode, nodes are already defined by AddNodeAt
            // In legacy mode, auto-generate from segment endpoints
            bool has_topology = false;
            for (const auto& seg : segments_) {
                if (seg.node_from >= 0) { has_topology = true; break; }
            }
            if (!has_topology && nodes_.empty()) {
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

    // Build incidence matrix if topology is present
    bool has_topology = false;
    for (const auto& seg : segments_) {
        if (seg.node_from >= 0) { has_topology = true; break; }
    }
    if (has_topology) {
        BuildIncidenceMatrix(matrices);
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
    // Full Rosa/Grover formula for rectangular cross-section
    //
    // Reference:
    //   Rosa, NBS Scientific Paper 169 (1908)
    //   Grover, "Inductance Calculations", Dover, 1946
    //   Ruehli, IBM J. Res. Dev., 16(5), 470-481, 1972
    //
    // Exact double-volume integral of 1/R for rectangular conductor.
    // L = mu_0 * l * Y(l, w, h)
    // where Y contains asinh, atan, and rational correction terms.

    double l = seg.length;
    double w = seg.width;
    double h = seg.height;

    // Minimum dimension check
    if (w < 1e-15) w = 1e-10;
    if (h < 1e-15) h = 1e-10;
    if (l < 1e-15) return 0.0;

    // Normalize by length
    double wn = w / l;
    double tn = h / l;

    // Derived quantities
    double aw = std::sqrt(wn * wn + 1.0);
    double at = std::sqrt(tn * tn + 1.0);
    double r  = std::sqrt(wn * wn + tn * tn);
    double ar = std::sqrt(wn * wn + tn * tn + 1.0);

    double z = 0.0;

    // Base terms (1/4 coefficient)
    z += 0.25 * (
        (1.0 / wn) * std::asinh(wn / at) +
        (1.0 / tn) * std::asinh(tn / aw) +
        std::asinh(1.0 / r)
    );

    // Correction terms (1/24 coefficient) - 6 additional terms
    z += (1.0 / 24.0) * (
        (tn * tn / wn) * std::asinh(wn / (tn * at * (r + ar))) +
        (wn * wn / tn) * std::asinh(tn / (wn * aw * (r + ar))) +
        (tn * tn / (wn * wn)) * std::asinh(wn * wn / (tn * r * (at + ar))) +
        (wn * wn / (tn * tn)) * std::asinh(tn * tn / (wn * r * (aw + ar))) +
        (1.0 / (wn * tn * tn)) * std::asinh(wn * tn * tn / (at * (aw + ar))) +
        (1.0 / (tn * wn * wn)) * std::asinh(tn * wn * wn / (aw * (at + ar)))
    );

    // Arctangent corrections (1/6 coefficient)
    z -= (1.0 / 6.0) * (
        (1.0 / (wn * tn)) * std::atan(wn * tn / ar) +
        (tn / wn) * std::atan(wn / (tn * ar)) +
        (wn / tn) * std::atan(tn / (wn * ar))
    );

    // Rational corrections (1/60 coefficient)
    z -= (1.0 / 60.0) * (
        ((ar + r + tn + at) * tn * tn) /
            ((ar + r) * (r + tn) * (tn + at) * (at + ar)) +
        ((ar + r + wn + aw) * wn * wn) /
            ((ar + r) * (r + wn) * (wn + aw) * (aw + ar)) +
        (ar + aw + 1.0 + at) /
            ((ar + aw) * (aw + 1.0) * (1.0 + at) * (at + ar))
    );

    // Inverse distance corrections (1/20 coefficient)
    z -= (1.0 / 20.0) * (
        1.0 / (r + ar) + 1.0 / (aw + ar) + 1.0 / (at + ar)
    );

    z *= (2.0 / RadConst::PI);
    z *= l;

    return z * PEEC_MU_0;
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
    // Neumann integral for mutual inductance between two straight filaments
    //
    // M_ij = (mu_0 / 4*pi) * (d_i . d_j) * integral integral dt ds / |r_i(t) - r_j(s)|
    //
    // For parallel filaments (|d_i . d_j| > 0.999): analytical Rosa/Grover formula
    // For general case: 4-point Gauss quadrature

    // Direction dot product
    double dot = seg_i.direction.x * seg_j.direction.x +
                 seg_i.direction.y * seg_j.direction.y +
                 seg_i.direction.z * seg_j.direction.z;

    if (std::abs(dot) < 1e-10) return 0.0;  // Perpendicular filaments

    // Center-to-center vector
    double rx = seg_j.center.x - seg_i.center.x;
    double ry = seg_j.center.y - seg_i.center.y;
    double rz = seg_j.center.z - seg_i.center.z;

    double l_i = seg_i.length;
    double l_j = seg_j.length;

    // Check if filaments are nearly parallel
    if (std::abs(std::abs(dot) - 1.0) < 1e-3) {
        // Parallel filaments: use analytical Neumann/Rosa/Grover formula
        //
        // Reference: F. W. Grover, "Inductance Calculations", Dover, 1946
        //
        // M = (mu_0/(4*pi)) * [F(alpha, d) + F(beta, d) - F(gamma, d) - F(delta, d)]
        //
        // where F(x, d) = x * arsinh(x/d) - sqrt(x^2 + d^2)
        //   (F is an even function of x)
        //
        // alpha = (l_i + l_j)/2 + p,  beta = (l_i + l_j)/2 - p
        // gamma = (l_i - l_j)/2 + p,  delta = (l_j - l_i)/2 + p = -gamma + ... hmm
        //
        // Using exact derivation from double integration:
        //   b1 = l_i/2, a1 = -l_i/2 (filament i limits)
        //   b2 = p + l_j/2, a2 = p - l_j/2 (filament j limits shifted by axial offset p)
        //   M = (mu_0/(4*pi)) * [F(b1-a2) + F(a1-b2) - F(b1-b2) - F(a1-a2)]

        // Axial offset: projection of r onto direction
        double p = rx * seg_i.direction.x + ry * seg_i.direction.y + rz * seg_i.direction.z;
        if (dot < 0) p = -p;  // Account for anti-parallel

        // Perpendicular distance
        // d_perp = |r - p*d_i|
        double px = p * seg_i.direction.x;
        double py = p * seg_i.direction.y;
        double pz = p * seg_i.direction.z;
        double dpx = rx - px;
        double dpy = ry - py;
        double dpz = rz - pz;
        double d_perp = std::sqrt(dpx*dpx + dpy*dpy + dpz*dpz);

        if (d_perp < 1e-15) {
            // Collinear filaments - use Gauss quadrature (falls through below)
        } else {
            // F(x, d) = x * arsinh(x/d) - sqrt(x^2 + d^2)
            // F is even in x, so F(-x, d) = F(x, d)
            // F(x, d) = x * arsinh(x/d) - sqrt(x^2 + d^2)
            // Note: x + sqrt(x^2+d^2) > 0 for all x when d > 0, so log is safe.
            // Do NOT use abs(x) here - it breaks the formula for negative x.
            auto F = [](double x, double d) -> double {
                double x2d2 = x*x + d*d;
                return x * std::log((x + std::sqrt(x2d2)) / d) - std::sqrt(x2d2);
            };

            // Filament i: from -l_i/2 to +l_i/2
            // Filament j: from p - l_j/2 to p + l_j/2
            double b1 = l_i / 2.0;
            double a1 = -l_i / 2.0;
            double b2 = p + l_j / 2.0;
            double a2 = p - l_j / 2.0;

            double M = (PEEC_MU_0 * PEEC_INV_FOUR_PI) *
                       (F(b1 - a2, d_perp) + F(a1 - b2, d_perp) -
                        F(b1 - b2, d_perp) - F(a1 - a2, d_perp));

            return dot * M;  // Preserve sign for anti-parallel filaments
        }
    }

    // Check if segments are close (use fourfil for near-field)
    double r_dist = std::sqrt(rx*rx + ry*ry + rz*rz);
    double max_cs = std::max({seg_i.width, seg_i.height, seg_j.width, seg_j.height});
    if (r_dist < 3.0 * max_cs && max_cs > 1e-15) {
        // Near-field: use recursive fourfil subdivision
        return MutualInductanceFourfil(seg_i, seg_j, 2);
    }

    // General case: 8-point Gauss-Legendre quadrature
    static const double gp[] = {
        -0.9602898564975363, -0.7966664774136267,
        -0.5255324099163290, -0.1834346424956498,
         0.1834346424956498,  0.5255324099163290,
         0.7966664774136267,  0.9602898564975363
    };
    static const double gw[] = {
         0.1012285362903763,  0.2223810344533745,
         0.3137066458778873,  0.3626837833783620,
         0.3626837833783620,  0.3137066458778873,
         0.2223810344533745,  0.1012285362903763
    };
    static const int ng = 8;

    double sum = 0.0;
    for (int ki = 0; ki < ng; ++ki) {
        double ti = gp[ki] * (l_i / 2.0);
        double xi = seg_i.center.x + ti * seg_i.direction.x;
        double yi = seg_i.center.y + ti * seg_i.direction.y;
        double zi = seg_i.center.z + ti * seg_i.direction.z;

        for (int kj = 0; kj < ng; ++kj) {
            double tj = gp[kj] * (l_j / 2.0);
            double xj = seg_j.center.x + tj * seg_j.direction.x;
            double yj = seg_j.center.y + tj * seg_j.direction.y;
            double zj = seg_j.center.z + tj * seg_j.direction.z;

            double ddx = xi - xj;
            double ddy = yi - yj;
            double ddz = zi - zj;
            double dist = std::sqrt(ddx*ddx + ddy*ddy + ddz*ddz);

            if (dist > 1e-15) {
                sum += gw[ki] * gw[kj] / dist;
            }
        }
    }

    // Scale: integral was over [-1,1]x[-1,1], actual limits are [-l/2,l/2]
    // Jacobian: (l_i/2) * (l_j/2)
    return (PEEC_MU_0 * PEEC_INV_FOUR_PI) * dot * sum * (l_i / 2.0) * (l_j / 2.0);
}

double PEECMatrixBuilder::MutualInductanceFourfil(const PEECSegment& seg_i,
                                                   const PEECSegment& seg_j,
                                                   int depth) const {
    // Recursive fourfil subdivision for near-field mutual inductance.
    // Each filament is split into 2x2 cross-section sub-filaments.
    // M = average of all 4 sub-filament pair mutual inductances.
    //
    // Reference: A. E. Ruehli, "Inductance Calculations in a Complex
    //   Integrated Circuit Environment", IBM J. Res. Dev., 1972.
    // Also: FastMaxwell src/calcaoneoverr.h fourfil() approach.

    if (depth <= 0) {
        // Base case: compute mutual inductance via 8-point Gauss quadrature
        double dot = seg_i.direction.x * seg_j.direction.x +
                     seg_i.direction.y * seg_j.direction.y +
                     seg_i.direction.z * seg_j.direction.z;
        if (std::abs(dot) < 1e-10) return 0.0;

        double l_i = seg_i.length;
        double l_j = seg_j.length;

        static const double gp[] = {
            -0.9602898564975363, -0.7966664774136267,
            -0.5255324099163290, -0.1834346424956498,
             0.1834346424956498,  0.5255324099163290,
             0.7966664774136267,  0.9602898564975363
        };
        static const double gw[] = {
             0.1012285362903763,  0.2223810344533745,
             0.3137066458778873,  0.3626837833783620,
             0.3626837833783620,  0.3137066458778873,
             0.2223810344533745,  0.1012285362903763
        };

        double sum = 0.0;
        for (int ki = 0; ki < 8; ++ki) {
            double ti = gp[ki] * (l_i / 2.0);
            double xi = seg_i.center.x + ti * seg_i.direction.x;
            double yi = seg_i.center.y + ti * seg_i.direction.y;
            double zi = seg_i.center.z + ti * seg_i.direction.z;
            for (int kj = 0; kj < 8; ++kj) {
                double tj = gp[kj] * (l_j / 2.0);
                double xj = seg_j.center.x + tj * seg_j.direction.x;
                double yj = seg_j.center.y + tj * seg_j.direction.y;
                double zj = seg_j.center.z + tj * seg_j.direction.z;
                double ddx = xi - xj; double ddy = yi - yj; double ddz = zi - zj;
                double dist = std::sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
                if (dist > 1e-15) sum += gw[ki] * gw[kj] / dist;
            }
        }
        return (PEEC_MU_0 * PEEC_INV_FOUR_PI) * dot * sum * (l_i / 2.0) * (l_j / 2.0);
    }

    // Build local coordinate system for cross-section subdivision.
    // We need two vectors perpendicular to the filament direction.
    // Use the segment's width/height axes (approximate from direction).

    auto make_perp_axes = [](const TVector3d& dir, TVector3d& u, TVector3d& v) {
        // Find a vector not parallel to dir
        TVector3d ref;
        if (std::fabs(dir.x) < 0.9) { ref.x = 1; ref.y = 0; ref.z = 0; }
        else { ref.x = 0; ref.y = 1; ref.z = 0; }
        // u = dir × ref (normalized)
        u.x = dir.y * ref.z - dir.z * ref.y;
        u.y = dir.z * ref.x - dir.x * ref.z;
        u.z = dir.x * ref.y - dir.y * ref.x;
        double norm_u = std::sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
        u.x /= norm_u; u.y /= norm_u; u.z /= norm_u;
        // v = dir × u
        v.x = dir.y * u.z - dir.z * u.y;
        v.y = dir.z * u.x - dir.x * u.z;
        v.z = dir.x * u.y - dir.y * u.x;
    };

    // Create 4 sub-filaments for seg_j (2x2 subdivision of cross-section)
    TVector3d u_j, v_j;
    make_perp_axes(seg_j.direction, u_j, v_j);

    double half_w = seg_j.width / 4.0;   // offset = w/4 from center
    double half_h = seg_j.height / 4.0;

    double M_sum = 0.0;
    for (int iw = -1; iw <= 1; iw += 2) {
        for (int ih = -1; ih <= 1; ih += 2) {
            PEECSegment sub_j = seg_j;
            sub_j.width = seg_j.width / 2.0;
            sub_j.height = seg_j.height / 2.0;
            sub_j.center.x = seg_j.center.x + iw * half_w * u_j.x + ih * half_h * v_j.x;
            sub_j.center.y = seg_j.center.y + iw * half_w * u_j.y + ih * half_h * v_j.y;
            sub_j.center.z = seg_j.center.z + iw * half_w * u_j.z + ih * half_h * v_j.z;

            M_sum += MutualInductanceFourfil(seg_i, sub_j, depth - 1);
        }
    }

    return M_sum / 4.0;  // Average over 4 sub-filaments
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
// Panel Analytical Integration - Hess-Smith Edge Integration
// Reference: Arcioni, Bressan, Perregrini, IEEE MTT, vol. 45, 1997
// ============================================================================

double PEECMatrixBuilder::HessSmithPotential(const std::vector<TVector3d>& vertices,
                                              const TVector3d& normal,
                                              const TVector3d& obs_point) const {
    // Analytical integration of 1/R over a flat polygon, evaluated at obs_point.
    // Returns integral = integral_S (1/|r - r'|) dS'
    //
    // Uses per-edge analytical primitives (log for SLP, atan for solid angle).
    // Works for triangles and quads (any convex polygon).

    const int nv = static_cast<int>(vertices.size());
    if (nv < 3) return 0.0;

    // Height of observation point above panel plane
    double hx = obs_point.x - vertices[0].x;
    double hy = obs_point.y - vertices[0].y;
    double hz = obs_point.z - vertices[0].z;
    double h = hx * normal.x + hy * normal.y + hz * normal.z;
    double abs_h = std::fabs(h);

    double slp_sum = 0.0;  // Single layer potential (log terms)
    double dlp_sum = 0.0;  // Solid angle (atan terms)

    for (int i = 0; i < nv; ++i) {
        int j = (i + 1) % nv;

        // Vector from obs_point to edge vertices
        double r0x = vertices[i].x - obs_point.x;
        double r0y = vertices[i].y - obs_point.y;
        double r0z = vertices[i].z - obs_point.z;
        double r1x = vertices[j].x - obs_point.x;
        double r1y = vertices[j].y - obs_point.y;
        double r1z = vertices[j].z - obs_point.z;

        double R0 = std::sqrt(r0x*r0x + r0y*r0y + r0z*r0z);
        double R1 = std::sqrt(r1x*r1x + r1y*r1y + r1z*r1z);

        if (R0 < 1e-15 || R1 < 1e-15) continue;

        // Edge vector and length
        double ex = vertices[j].x - vertices[i].x;
        double ey = vertices[j].y - vertices[i].y;
        double ez = vertices[j].z - vertices[i].z;
        double l_edge = std::sqrt(ex*ex + ey*ey + ez*ez);
        if (l_edge < 1e-15) continue;

        // Unit edge tangent
        double tx = ex / l_edge;
        double ty = ey / l_edge;
        double tz = ez / l_edge;

        // Outward edge normal in panel plane: m = n × t
        double mx = normal.y * tz - normal.z * ty;
        double my = normal.z * tx - normal.x * tz;
        double mz = normal.x * ty - normal.y * tx;

        // Projections onto edge coordinate system
        double d = r0x * mx + r0y * my + r0z * mz;  // perpendicular distance to edge line
        double s0 = r0x * tx + r0y * ty + r0z * tz;  // projection along edge from vertex i
        double s1 = r1x * tx + r1y * ty + r1z * tz;  // projection along edge from vertex j

        // SLP term: log((R0 + s0) / (R1 + s1)) * |d| or equivalent stable form
        // Using Newman formula: log((R0 - s0 + R1 + s1) / (R0 + s0 + R1 - s1)) * d
        // which is equivalent but more numerically stable
        double num = R0 + R1 + l_edge;
        double den = R0 + R1 - l_edge;
        if (den > 1e-15 && num > 1e-15) {
            double log_term = std::log(num / den);
            slp_sum += d * log_term;
        }

        // DLP (solid angle) term - only if h != 0
        if (abs_h > 1e-15) {
            // Solid angle contribution from this edge
            // atan2(h * l_edge * d, R0*R1*(R0*R1 + r0·r1))
            // where r0·r1 = dot product of r0, r1
            double dot01 = r0x*r1x + r0y*r1y + r0z*r1z;
            double denom = R0 * R1 + dot01;
            if (std::fabs(denom) > 1e-30) {
                // Cross product component along normal: (r0 × r1) · n
                double cx = r0y*r1z - r0z*r1y;
                double cy = r0z*r1x - r0x*r1z;
                double cz = r0x*r1y - r0y*r1x;
                double cross_n = cx * normal.x + cy * normal.y + cz * normal.z;
                dlp_sum += std::atan2(h * cross_n, R0 * R1 + dot01 * 1.0);
            }
        }
    }

    // Result: integral_S 1/R dS = slp_sum - |h| * dlp_sum  (for h > 0)
    //                            = slp_sum + |h| * dlp_sum  (for h < 0)
    // General: slp_sum - h * sign(h) * |dlp_sum| ... actually:
    // The formula is: integral = slp_sum - h * dlp_sum
    // where dlp_sum already has the correct sign from atan2
    double result = slp_sum - h * dlp_sum;

    return result;
}

double PEECMatrixBuilder::SelfPotentialPanelTriangle(const PEECPanel& panel) const {
    // Self-potential for triangular panel using Hess-Smith analytical integration.
    // Evaluates integral_S integral_S' 1/|r-r'| dS dS' / area
    // by 7-point Gauss quadrature on test panel + analytical on source.

    if (panel.type != PEECPanel::Triangle || panel.vertices.size() != 3) {
        return 0.0;
    }

    // 7-point Gauss quadrature for triangles (Strang & Fix)
    const double q7_w[7] = {
        0.225 / 2.0,
        0.13239415, 0.13239415, 0.13239415,
        0.12593918, 0.12593918, 0.12593918
    };
    const double q7_bc[7][3] = {
        {1.0/3.0, 1.0/3.0, 1.0/3.0},
        {0.05971587, 0.47014206, 0.47014206},
        {0.47014206, 0.05971587, 0.47014206},
        {0.47014206, 0.47014206, 0.05971587},
        {0.79742699, 0.10128651, 0.10128651},
        {0.10128651, 0.79742699, 0.10128651},
        {0.10128651, 0.10128651, 0.79742699}
    };

    const TVector3d& v0 = panel.vertices[0];
    const TVector3d& v1 = panel.vertices[1];
    const TVector3d& v2 = panel.vertices[2];

    // Small offset along normal for self-term (avoid singularity)
    double char_len = std::sqrt(panel.area);
    double offset = char_len * 0.01;  // 1% of characteristic length

    double sum = 0.0;
    for (int q = 0; q < 7; ++q) {
        TVector3d obs;
        obs.x = q7_bc[q][0] * v0.x + q7_bc[q][1] * v1.x + q7_bc[q][2] * v2.x;
        obs.y = q7_bc[q][0] * v0.y + q7_bc[q][1] * v1.y + q7_bc[q][2] * v2.y;
        obs.z = q7_bc[q][0] * v0.z + q7_bc[q][1] * v1.z + q7_bc[q][2] * v2.z;

        // Offset observation point slightly along normal
        obs.x += offset * panel.normal.x;
        obs.y += offset * panel.normal.y;
        obs.z += offset * panel.normal.z;

        double val = HessSmithPotential(panel.vertices, panel.normal, obs);
        sum += q7_w[q] * val;
    }

    // Jacobian: 2*A for triangle quadrature
    return (2.0 * panel.area * sum) / (4.0 * RadConst::PI * PEEC_EPS_0);
}

double PEECMatrixBuilder::SelfPotentialPanelQuad(const PEECPanel& panel) const {
    // Self-potential for quad panel using Hess-Smith analytical integration.
    // 4-point Gauss quadrature on test quad + HessSmith on source quad.

    if (panel.type != PEECPanel::Quadrilateral || panel.vertices.size() != 4) {
        return 0.0;
    }

    // 2x2 Gauss-Legendre on quad mapped via bilinear shape functions
    const double gp = 1.0 / std::sqrt(3.0);
    const double quad_pts[4][2] = {
        {-gp, -gp}, {gp, -gp}, {gp, gp}, {-gp, gp}
    };

    const TVector3d& v0 = panel.vertices[0];
    const TVector3d& v1 = panel.vertices[1];
    const TVector3d& v2 = panel.vertices[2];
    const TVector3d& v3 = panel.vertices[3];

    double char_len = std::sqrt(panel.area);
    double offset = char_len * 0.01;

    double sum = 0.0;
    for (int q = 0; q < 4; ++q) {
        double xi = quad_pts[q][0];
        double eta = quad_pts[q][1];

        // Bilinear shape functions
        double N0 = 0.25 * (1 - xi) * (1 - eta);
        double N1 = 0.25 * (1 + xi) * (1 - eta);
        double N2 = 0.25 * (1 + xi) * (1 + eta);
        double N3 = 0.25 * (1 - xi) * (1 + eta);

        TVector3d obs;
        obs.x = N0*v0.x + N1*v1.x + N2*v2.x + N3*v3.x + offset*panel.normal.x;
        obs.y = N0*v0.y + N1*v1.y + N2*v2.y + N3*v3.y + offset*panel.normal.y;
        obs.z = N0*v0.z + N1*v1.z + N2*v2.z + N3*v3.z + offset*panel.normal.z;

        // Jacobian for bilinear mapping
        double dxdxi  = 0.25*(-(1-eta)*v0.x + (1-eta)*v1.x + (1+eta)*v2.x - (1+eta)*v3.x);
        double dydxi  = 0.25*(-(1-eta)*v0.y + (1-eta)*v1.y + (1+eta)*v2.y - (1+eta)*v3.y);
        double dzdxi  = 0.25*(-(1-eta)*v0.z + (1-eta)*v1.z + (1+eta)*v2.z - (1+eta)*v3.z);
        double dxdeta = 0.25*(-(1-xi)*v0.x - (1+xi)*v1.x + (1+xi)*v2.x + (1-xi)*v3.x);
        double dydeta = 0.25*(-(1-xi)*v0.y - (1+xi)*v1.y + (1+xi)*v2.y + (1-xi)*v3.y);
        double dzdeta = 0.25*(-(1-xi)*v0.z - (1+xi)*v1.z + (1+xi)*v2.z + (1-xi)*v3.z);

        // Cross product for Jacobian magnitude
        double cx = dydxi*dzdeta - dzdxi*dydeta;
        double cy = dzdxi*dxdeta - dxdxi*dzdeta;
        double cz = dxdxi*dydeta - dydxi*dxdeta;
        double jac = std::sqrt(cx*cx + cy*cy + cz*cz);

        double val = HessSmithPotential(panel.vertices, panel.normal, obs);
        sum += val * jac;  // weight = 1.0 for 2x2 Gauss
    }

    return sum / (4.0 * RadConst::PI * PEEC_EPS_0);
}

double PEECMatrixBuilder::MutualPotentialPanelTriangle(const PEECPanel& panel_i,
                                                        const PEECPanel& panel_j) const {
    // Mutual potential between two triangular panels.
    // Uses 7-point Gauss on test panel_i + HessSmith analytical on source panel_j.
    //
    // Far-field (distance > 5 * panel_size): centroid approximation.
    // Near/mid-field: Gauss + HessSmith (no singularity issues).

    if (panel_i.type != PEECPanel::Triangle || panel_j.type != PEECPanel::Triangle) {
        return 0.0;
    }

    double dx = panel_i.center.x - panel_j.center.x;
    double dy = panel_i.center.y - panel_j.center.y;
    double dz = panel_i.center.z - panel_j.center.z;
    double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

    double char_size = std::max(std::sqrt(panel_i.area), std::sqrt(panel_j.area));

    if (dist > 5.0 * char_size) {
        // Far-field: centroid approximation
        return (panel_i.area * panel_j.area) / (4.0 * RadConst::PI * PEEC_EPS_0 * dist);
    }

    // 7-point Gauss quadrature on test panel_i (Strang & Fix)
    const double q7_w[7] = {
        0.225 / 2.0,
        0.13239415, 0.13239415, 0.13239415,
        0.12593918, 0.12593918, 0.12593918
    };
    const double q7_bc[7][3] = {
        {1.0/3.0, 1.0/3.0, 1.0/3.0},
        {0.05971587, 0.47014206, 0.47014206},
        {0.47014206, 0.05971587, 0.47014206},
        {0.47014206, 0.47014206, 0.05971587},
        {0.79742699, 0.10128651, 0.10128651},
        {0.10128651, 0.79742699, 0.10128651},
        {0.10128651, 0.10128651, 0.79742699}
    };

    const TVector3d& v0 = panel_i.vertices[0];
    const TVector3d& v1 = panel_i.vertices[1];
    const TVector3d& v2 = panel_i.vertices[2];

    double sum = 0.0;
    for (int q = 0; q < 7; ++q) {
        TVector3d obs;
        obs.x = q7_bc[q][0]*v0.x + q7_bc[q][1]*v1.x + q7_bc[q][2]*v2.x;
        obs.y = q7_bc[q][0]*v0.y + q7_bc[q][1]*v1.y + q7_bc[q][2]*v2.y;
        obs.z = q7_bc[q][0]*v0.z + q7_bc[q][1]*v1.z + q7_bc[q][2]*v2.z;

        double val = HessSmithPotential(panel_j.vertices, panel_j.normal, obs);
        sum += q7_w[q] * val;
    }

    // Jacobian: 2*A_i for triangle quadrature on test panel
    return (2.0 * panel_i.area * sum) / (4.0 * RadConst::PI * PEEC_EPS_0);
}

// ============================================================================
// PEECSolver
//
// MNA multi-port solver using LAPACK zgetrf_/zgetrs_ (shared MKL with MSC).
// Same LAPACK include pattern as rad_relaxation_methods.cpp (dgesv_).
// ============================================================================

PEECSolver::PEECSolver()
    : frequency_(0), omega_(0), hasSurfaceImpedance_(false),
      n_nodes_(0), hasTopology_(false),
      solverMethod_(0), bicgstab_tol_(1e-10), bicgstab_max_iter_(1000) {}

PEECSolver::~PEECSolver() {}

void PEECSolver::SetMatrices(const PEECMatrices& matrices) {
    matrices_ = matrices;
}

void PEECSolver::SetSegmentNodes(const std::vector<std::pair<int,int>>& seg_nodes, int n_nodes) {
    segment_nodes_ = seg_nodes;
    n_nodes_ = n_nodes;
    hasTopology_ = true;
}

void PEECSolver::SetPorts(const std::vector<PEECPort>& ports) {
    ports_ = ports;
}

// ---- Legacy API (upgraded to LAPACK zgesv_) ----

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

    std::vector<std::complex<double>> Z;
    BuildImpedanceMatrix(Z);

    std::vector<std::complex<double>> V(n_total, std::complex<double>(0, 0));
    for (int i = 0; i < n_loop && i < static_cast<int>(portVector.size()); ++i) {
        V[i] = std::complex<double>(portVector[i], 0);
    }

    std::vector<std::complex<double>> I(n_total);
    Solve(V, I);

    std::complex<double> Z_port(0, 0);
    for (int i = 0; i < n_loop && i < static_cast<int>(portVector.size()); ++i) {
        Z_port += portVector[i] * I[i];
    }

    return Z_port;
}

void PEECSolver::Solve(const std::vector<std::complex<double>>& V,
                        std::vector<std::complex<double>>& I) {
    int n = static_cast<int>(V.size());
    I = V;  // RHS will be overwritten with solution

    std::vector<std::complex<double>> Z;
    BuildImpedanceMatrix(Z);

#ifdef HAVE_LAPACK
    // Use LAPACK zgesv_ (complex LU) - same pattern as MSC dgesv_ in
    // rad_relaxation_methods.cpp:1522-1532

    // Transpose row-major -> column-major for LAPACK (same approach as MSC)
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            std::swap(Z[i * n + j], Z[j * n + i]);
        }
    }

    int ln = n, nrhs = 1, info = 0;
    std::vector<int> ipiv(n);
    zgesv_(&ln, &nrhs, reinterpret_cast<MKL_Complex16*>(Z.data()), &ln,
           ipiv.data(), reinterpret_cast<MKL_Complex16*>(I.data()), &ln, &info);
#else
    // Fallback: naive Gaussian elimination with partial pivoting
    std::vector<std::complex<double>> A = Z;
    std::vector<std::complex<double>> b = V;

    for (int k = 0; k < n - 1; ++k) {
        int maxIdx = k;
        double maxVal = std::abs(A[k * n + k]);
        for (int i = k + 1; i < n; ++i) {
            double val = std::abs(A[i * n + k]);
            if (val > maxVal) { maxVal = val; maxIdx = i; }
        }
        if (maxIdx != k) {
            for (int j = 0; j < n; ++j) std::swap(A[k * n + j], A[maxIdx * n + j]);
            std::swap(b[k], b[maxIdx]);
        }
        std::complex<double> pivot = A[k * n + k];
        if (std::abs(pivot) < 1e-30) continue;
        for (int i = k + 1; i < n; ++i) {
            std::complex<double> factor = A[i * n + k] / pivot;
            for (int j = k; j < n; ++j) A[i * n + j] -= factor * A[k * n + j];
            b[i] -= factor * b[k];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        std::complex<double> sum = b[i];
        for (int j = i + 1; j < n; ++j) sum -= A[i * n + j] * I[j];
        std::complex<double> diag = A[i * n + i];
        I[i] = (std::abs(diag) > 1e-30) ? sum / diag : std::complex<double>(0, 0);
    }
#endif
}

// ---- MNA Internal Methods ----

void PEECSolver::BuildFullIncidenceMatrix(std::vector<double>& A_full) {
    int n_loop = matrices_.n_loop;
    A_full.assign(n_nodes_ * n_loop, 0.0);

    for (int f = 0; f < n_loop && f < static_cast<int>(segment_nodes_.size()); ++f) {
        int nf = segment_nodes_[f].first;   // node_from
        int nt = segment_nodes_[f].second;  // node_to

        if (nf >= 0 && nf < n_nodes_) {
            A_full[nf * n_loop + f] = +1.0;  // filament leaves node_from
        }
        if (nt >= 0 && nt < n_nodes_) {
            A_full[nt * n_loop + f] = -1.0;  // filament enters node_to
        }
    }
}

void PEECSolver::FindConnectedComponents(std::vector<std::vector<int>>& components) {
    components.clear();

    // Build adjacency list from segment_nodes
    std::vector<std::set<int>> adj(n_nodes_);
    for (const auto& seg : segment_nodes_) {
        int nf = seg.first;
        int nt = seg.second;
        if (nf >= 0 && nf < n_nodes_ && nt >= 0 && nt < n_nodes_) {
            adj[nf].insert(nt);
            adj[nt].insert(nf);
        }
    }

    // BFS to find connected components
    std::vector<bool> visited(n_nodes_, false);
    for (int start = 0; start < n_nodes_; ++start) {
        if (visited[start]) continue;

        std::vector<int> comp;
        std::queue<int> q;
        q.push(start);
        visited[start] = true;

        while (!q.empty()) {
            int node = q.front();
            q.pop();
            comp.push_back(node);

            for (int neighbor : adj[node]) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    q.push(neighbor);
                }
            }
        }
        components.push_back(std::move(comp));
    }
}

void PEECSolver::SelectGroundNodes(const std::vector<std::vector<int>>& components,
                                    std::vector<int>& ground_nodes) {
    ground_nodes.clear();
    std::set<int> comp_set;

    for (const auto& comp : components) {
        comp_set.clear();
        for (int n : comp) comp_set.insert(n);

        // Prefer negative terminal of a port in this component
        bool grounded = false;
        for (const auto& port : ports_) {
            if (comp_set.count(port.node_negative) && port.node_negative != port.node_positive) {
                ground_nodes.push_back(port.node_negative);
                grounded = true;
                break;
            }
        }
        if (!grounded) {
            // Ground the smallest-numbered node
            ground_nodes.push_back(*std::min_element(comp.begin(), comp.end()));
        }
    }
}

void PEECSolver::BuildZBranch(double freq,
                               const std::complex<double>* Zs, int n_Zs,
                               std::vector<std::complex<double>>& Z_branch) {
    int n_loop = matrices_.n_loop;
    double omega = 2.0 * RadConst::PI * freq;

    Z_branch.assign(n_loop * n_loop, std::complex<double>(0, 0));

    // Z_branch = diag(R) + jw*L
    for (int i = 0; i < n_loop; ++i) {
        for (int j = 0; j < n_loop; ++j) {
            Z_branch[i * n_loop + j] = std::complex<double>(0, omega * matrices_.L_at(i, j));
        }
        Z_branch[i * n_loop + i] += std::complex<double>(matrices_.R[i], 0);
    }

    // Add surface impedance (diagonal)
    if (Zs != nullptr) {
        int n = std::min(n_Zs, n_loop);
        for (int i = 0; i < n; ++i) {
            Z_branch[i * n_loop + i] += Zs[i];
        }
    }
}

void PEECSolver::BuildZEff(double freq,
                             const std::complex<double>* Zs, int n_Zs,
                             std::vector<std::complex<double>>& Z_eff) {
    int n_loop = matrices_.n_loop;
    int n_star = matrices_.n_star;
    double omega = 2.0 * RadConst::PI * freq;

    BuildZBranch(freq, Zs, n_Zs, Z_eff);

    // If no panels or DC, Z_eff = Z_branch
    if (n_star <= 0 || matrices_.P.empty() || matrices_.M_LS.empty() || freq < 1e-10) {
        return;
    }

#ifdef HAVE_LAPACK
    // Schur complement: Z_eff = Z_LL - Z_LS * Z_SS^{-1} * Z_SL
    // Z_SS = P / (jw), Z_LS = jw * M_LS, Z_SL = Z_LS^T

    std::complex<double> jw(0, omega);
    std::complex<double> inv_jw(0, -1.0 / omega);

    // Build Z_SS (column-major for LAPACK inversion)
    std::vector<std::complex<double>> Z_SS(n_star * n_star);
    for (int i = 0; i < n_star; ++i) {
        for (int j = 0; j < n_star; ++j) {
            Z_SS[j * n_star + i] = matrices_.P_at(i, j) * inv_jw;  // column-major
        }
    }

    // Build Z_LS (row-major)
    std::vector<std::complex<double>> Z_LS(n_loop * n_star);
    for (int i = 0; i < n_loop; ++i) {
        for (int j = 0; j < n_star; ++j) {
            Z_LS[i * n_star + j] = jw * std::complex<double>(matrices_.M_LS_at(i, j), 0);
        }
    }

    // Invert Z_SS using zgesv_ with identity RHS
    // Solve Z_SS * X = I -> X = Z_SS^{-1}
    std::vector<std::complex<double>> Z_SS_inv(n_star * n_star, std::complex<double>(0, 0));
    for (int i = 0; i < n_star; ++i) Z_SS_inv[i * n_star + i] = std::complex<double>(1, 0);

    int ln_star = n_star, nrhs = n_star, info = 0;
    std::vector<int> ipiv(n_star);
    zgesv_(&ln_star, &nrhs, reinterpret_cast<MKL_Complex16*>(Z_SS.data()), &ln_star,
           ipiv.data(), reinterpret_cast<MKL_Complex16*>(Z_SS_inv.data()), &ln_star, &info);
    // Z_SS_inv is now Z_SS^{-1} in column-major

    // Transpose Z_SS_inv back to row-major
    std::vector<std::complex<double>> Z_SS_inv_row(n_star * n_star);
    for (int i = 0; i < n_star; ++i) {
        for (int j = 0; j < n_star; ++j) {
            Z_SS_inv_row[i * n_star + j] = Z_SS_inv[j * n_star + i];
        }
    }

    // temp = Z_LS * Z_SS_inv (n_loop x n_star)
    std::complex<double> alpha(1, 0), beta(0, 0);
    std::vector<std::complex<double>> temp(n_loop * n_star);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n_loop, n_star, n_star,
                &alpha, Z_LS.data(), n_star,
                Z_SS_inv_row.data(), n_star,
                &beta, temp.data(), n_star);

    // Z_eff -= temp * Z_SL = temp * Z_LS^T
    std::complex<double> neg_alpha(-1, 0);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                n_loop, n_loop, n_star,
                &neg_alpha, temp.data(), n_star,
                Z_LS.data(), n_star,
                &alpha, Z_eff.data(), n_loop);  // alpha=1 for additive
#else
    // Fallback: manual Schur complement (no BLAS)
    std::complex<double> jw(0, omega);
    std::complex<double> inv_jw(0, -1.0 / omega);

    // This path is slow but correct; HAVE_LAPACK should always be defined
    // for production builds.
#endif
}

void PEECSolver::MNASolveMultiPort(const std::vector<std::complex<double>>& Z_eff,
                                    std::complex<double>* Z_out, int n_ports_out) {
    int n_loop = matrices_.n_loop;
    int n_ports = static_cast<int>(ports_.size());
    if (n_ports == 0 || n_ports_out == 0) return;

    // Zero output
    for (int i = 0; i < n_ports * n_ports; ++i) Z_out[i] = std::complex<double>(0, 0);

#ifdef HAVE_LAPACK
    int info = 0;

    // Step 1: Invert Z_eff -> Y_branch
    // Method 0: LU via LAPACK zgesv_ (default)
    // Method 1: BiCGSTAB iterative solver (column-by-column)
    std::vector<std::complex<double>> Y_branch(n_loop * n_loop, std::complex<double>(0, 0));

    if (solverMethod_ == 1) {
        // BiCGSTAB: solve Z_eff * Y[:,j] = I[:,j] for each column
        bicgstab::DenseInvert<std::complex<double>>(
            n_loop, Z_eff.data(), Y_branch.data(),
            bicgstab_tol_, bicgstab_max_iter_);
    } else {
        // LU: zgesv_ with identity RHS
        std::vector<std::complex<double>> Z_copy(Z_eff);
        for (int i = 0; i < n_loop; ++i) Y_branch[i * n_loop + i] = std::complex<double>(1, 0);

        // Transpose to column-major for LAPACK
        for (int i = 0; i < n_loop; ++i) {
            for (int j = i + 1; j < n_loop; ++j) {
                std::swap(Z_copy[i * n_loop + j], Z_copy[j * n_loop + i]);
            }
        }
        for (int i = 0; i < n_loop; ++i) {
            for (int j = i + 1; j < n_loop; ++j) {
                std::swap(Y_branch[i * n_loop + j], Y_branch[j * n_loop + i]);
            }
        }

        int ln = n_loop, nrhs = n_loop, info = 0;
        std::vector<int> ipiv(n_loop);
        zgesv_(&ln, &nrhs, reinterpret_cast<MKL_Complex16*>(Z_copy.data()), &ln,
               ipiv.data(), reinterpret_cast<MKL_Complex16*>(Y_branch.data()), &ln, &info);

        if (info != 0) return;  // Singular matrix

        // Transpose back to row-major
        for (int i = 0; i < n_loop; ++i) {
            for (int j = i + 1; j < n_loop; ++j) {
                std::swap(Y_branch[i * n_loop + j], Y_branch[j * n_loop + i]);
            }
        }
    }

    // Step 2: Build A_full incidence matrix (n_nodes x n_loop, real)
    std::vector<double> A_full;
    BuildFullIncidenceMatrix(A_full);

    // Convert A_full to complex for BLAS zgemm
    std::vector<std::complex<double>> A_cmplx(n_nodes_ * n_loop);
    for (int i = 0; i < n_nodes_ * n_loop; ++i) {
        A_cmplx[i] = std::complex<double>(A_full[i], 0);
    }

    // Step 3: Y_node = A * Y_branch * A^T  (n_nodes x n_nodes)
    // temp = Y_branch * A^T  (n_loop x n_nodes)
    std::complex<double> alpha(1, 0), beta_zero(0, 0);
    std::vector<std::complex<double>> temp(n_loop * n_nodes_);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                n_loop, n_nodes_, n_loop,
                &alpha, Y_branch.data(), n_loop,
                A_cmplx.data(), n_loop,
                &beta_zero, temp.data(), n_nodes_);

    // Y_node = A * temp  (n_nodes x n_nodes)
    std::vector<std::complex<double>> Y_node(n_nodes_ * n_nodes_);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n_nodes_, n_nodes_, n_loop,
                &alpha, A_cmplx.data(), n_loop,
                temp.data(), n_nodes_,
                &beta_zero, Y_node.data(), n_nodes_);

    // Step 4: Find connected components and select ground nodes
    std::vector<std::vector<int>> components;
    FindConnectedComponents(components);

    std::vector<int> ground_nodes;
    SelectGroundNodes(components, ground_nodes);

    std::set<int> ground_set(ground_nodes.begin(), ground_nodes.end());
    int n_ground = static_cast<int>(ground_set.size());
    int n_reduced = n_nodes_ - n_ground;

    if (n_reduced <= 0) return;

    // Step 5: Build node ordering (non-ground first, ground last)
    std::vector<int> non_ground;
    non_ground.reserve(n_reduced);
    for (int i = 0; i < n_nodes_; ++i) {
        if (!ground_set.count(i)) non_ground.push_back(i);
    }

    std::vector<int> node_order = non_ground;
    for (int g : ground_nodes) node_order.push_back(g);

    // Inverse mapping: original node -> permuted index
    std::vector<int> inv_order(n_nodes_, -1);
    for (int i = 0; i < static_cast<int>(node_order.size()); ++i) {
        inv_order[node_order[i]] = i;
    }

    // Step 6: Extract Y_reduced (top-left n_reduced x n_reduced of permuted Y_node)
    // Store in column-major for LAPACK
    std::vector<std::complex<double>> Y_reduced(n_reduced * n_reduced);
    for (int i = 0; i < n_reduced; ++i) {
        for (int j = 0; j < n_reduced; ++j) {
            int orig_i = node_order[i];
            int orig_j = node_order[j];
            Y_reduced[j * n_reduced + i] = Y_node[orig_i * n_nodes_ + orig_j];  // column-major
        }
    }

    // Step 7: LU factorize Y_reduced (factor once, reuse for each port)
    int ln_red = n_reduced;
    std::vector<int> ipiv_red(n_reduced);
    zgetrf_(&ln_red, &ln_red, reinterpret_cast<MKL_Complex16*>(Y_reduced.data()),
            &ln_red, ipiv_red.data(), &info);

    if (info != 0) return;  // Singular matrix

    // Step 8: Solve for each port excitation
    for (int j = 0; j < n_ports; ++j) {
        int node_pos = ports_[j].node_positive;
        int node_neg = ports_[j].node_negative;

        // Current injection: +1A at positive, -1A at negative
        std::vector<std::complex<double>> I_ext(n_reduced, std::complex<double>(0, 0));

        if (!ground_set.count(node_pos)) {
            int idx = inv_order[node_pos];
            if (idx >= 0 && idx < n_reduced) I_ext[idx] += std::complex<double>(1, 0);
        }
        if (!ground_set.count(node_neg)) {
            int idx = inv_order[node_neg];
            if (idx >= 0 && idx < n_reduced) I_ext[idx] -= std::complex<double>(1, 0);
        }

        // Solve using pre-factored LU
        char trans = 'N';
        int nrhs_one = 1;
        zgetrs_(&trans, &ln_red, &nrhs_one,
                reinterpret_cast<MKL_Complex16*>(Y_reduced.data()), &ln_red,
                ipiv_red.data(), reinterpret_cast<MKL_Complex16*>(I_ext.data()),
                &ln_red, &info);

        // Map back to full voltage vector (ground nodes = 0V)
        std::vector<std::complex<double>> V_full(n_nodes_, std::complex<double>(0, 0));
        for (int idx = 0; idx < n_reduced; ++idx) {
            V_full[node_order[idx]] = I_ext[idx];
        }

        // Extract Z_mat: Z[i,j] = V_port_i when I_port_j = 1A
        for (int i = 0; i < n_ports; ++i) {
            Z_out[i * n_ports + j] = V_full[ports_[i].node_positive]
                                   - V_full[ports_[i].node_negative];
        }
    }

#else
    // No LAPACK: fallback not implemented for MNA multi-port
    // (requires topology support which is only meaningful with LAPACK)
#endif
}

// ---- MNA Public API ----

void PEECSolver::ComputeZMatrix(double freq,
                                 const std::complex<double>* Zs, int n_Zs,
                                 std::complex<double>* Z_out, int n_ports_out) {
    if (!hasTopology_) return;

    std::vector<std::complex<double>> Z_eff;
    BuildZEff(freq, Zs, n_Zs, Z_eff);
    MNASolveMultiPort(Z_eff, Z_out, n_ports_out);
}

void PEECSolver::ComputeCouplingCoefficient(double freq,
                                             const std::complex<double>* Zs, int n_Zs,
                                             double* k_out, double* L_out,
                                             int n_ports_out) {
    int n_ports = static_cast<int>(ports_.size());
    if (n_ports < 2) return;

    std::vector<std::complex<double>> Z_mat(n_ports * n_ports);
    ComputeZMatrix(freq, Zs, n_Zs, Z_mat.data(), n_ports_out);

    double omega = 2.0 * RadConst::PI * freq;
    if (omega < 1e-10) return;

    // L_ij = Im(Z_ij) / omega
    for (int i = 0; i < n_ports; ++i) {
        for (int j = 0; j < n_ports; ++j) {
            L_out[i * n_ports + j] = Z_mat[i * n_ports + j].imag() / omega;
        }
    }

    // k_ij = L_ij / sqrt(L_ii * L_jj)
    for (int i = 0; i < n_ports; ++i) {
        for (int j = 0; j < n_ports; ++j) {
            if (L_out[i * n_ports + i] > 0 && L_out[j * n_ports + j] > 0) {
                k_out[i * n_ports + j] = L_out[i * n_ports + j]
                    / std::sqrt(L_out[i * n_ports + i] * L_out[j * n_ports + j]);
            } else {
                k_out[i * n_ports + j] = 0.0;
            }
        }
    }
}

void PEECSolver::FrequencySweep(const double* freqs, int n_freq,
                                 const std::complex<double>* Zs_all,
                                 std::complex<double>* Z_out, int n_ports_out) {
    int n_ports = static_cast<int>(ports_.size());
    int n_loop = matrices_.n_loop;

    for (int f = 0; f < n_freq; ++f) {
        const std::complex<double>* Zs_f = (Zs_all != nullptr) ? &Zs_all[f * n_loop] : nullptr;
        int n_Zs = (Zs_all != nullptr) ? n_loop : 0;
        ComputeZMatrix(freqs[f], Zs_f, n_Zs,
                       &Z_out[f * n_ports * n_ports], n_ports_out);
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
