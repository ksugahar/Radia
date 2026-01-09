/*
 * rad_rwg_basis.cpp
 *
 * RWG Basis Functions Implementation
 *
 * Part of Radia project
 */

#include "rad_rwg_basis.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>

#ifdef HAVE_LAPACK
#include <mkl_lapacke.h>
#endif

namespace radia {

// ============================================================================
// Gauss-Legendre Quadrature Data for Triangle
// ============================================================================
// 7-point rule (degree 5 exactness)
// Coordinates are barycentric (L1, L2), L3 = 1 - L1 - L2
const double RWGIntegration::triQuadPts_[NQUAD_TRI][2] = {
    {0.333333333333333, 0.333333333333333},   // Centroid
    {0.059715871789770, 0.470142064105115},   // Edge points
    {0.470142064105115, 0.059715871789770},
    {0.470142064105115, 0.470142064105115},
    {0.797426985353087, 0.101286507323456},   // Vertex-near points
    {0.101286507323456, 0.797426985353087},
    {0.101286507323456, 0.101286507323456}
};

const double RWGIntegration::triQuadWts_[NQUAD_TRI] = {
    0.225000000000000,    // Centroid weight
    0.132394152788506,    // Edge weights
    0.132394152788506,
    0.132394152788506,
    0.125939180544827,    // Vertex-near weights
    0.125939180544827,
    0.125939180544827
};

// ============================================================================
// RWGMesh Implementation
// ============================================================================

RWGMesh::RWGMesh() : numInteriorEdges_(0) {
}

RWGMesh::~RWGMesh() {
}

void RWGMesh::CreateFromTriangles(const std::vector<TVector3d>& vertices,
                                   const std::vector<std::array<int, 3>>& triangles) {
    vertices_ = vertices;
    triangles_.resize(triangles.size());

    for (size_t i = 0; i < triangles.size(); ++i) {
        triangles_[i].vertices[0] = triangles[i][0];
        triangles_[i].vertices[1] = triangles[i][1];
        triangles_[i].vertices[2] = triangles[i][2];
    }

    ComputeTriangleProperties();
    BuildEdgeConnectivity();
    ComputeEdgeProperties();
}

void RWGMesh::CreateFromQuadPanels(const std::vector<TVector3d>& vertices,
                                    const std::vector<std::array<int, 4>>& quads) {
    vertices_ = vertices;
    triangles_.clear();
    triangles_.reserve(quads.size() * 2);

    // Split each quad into 2 triangles
    for (const auto& quad : quads) {
        // Triangle 1: vertices 0, 1, 2
        RWGTriangle tri1;
        tri1.vertices[0] = quad[0];
        tri1.vertices[1] = quad[1];
        tri1.vertices[2] = quad[2];
        triangles_.push_back(tri1);

        // Triangle 2: vertices 0, 2, 3
        RWGTriangle tri2;
        tri2.vertices[0] = quad[0];
        tri2.vertices[1] = quad[2];
        tri2.vertices[2] = quad[3];
        triangles_.push_back(tri2);
    }

    ComputeTriangleProperties();
    BuildEdgeConnectivity();
    ComputeEdgeProperties();
}

void RWGMesh::CreateWireMesh(const std::vector<TVector3d>& path,
                              double wireRadius,
                              int numAroundWire,
                              int numAlongSegment) {
    if (path.size() < 2) {
        throw std::invalid_argument("Wire path must have at least 2 points");
    }

    vertices_.clear();
    triangles_.clear();

    // Check if path is closed (first and last points are the same or very close)
    TVector3d pathDiff = path.back() - path.front();
    bool isClosedLoop = (pathDiff.Abs() < 1e-10);

    int nSegments = static_cast<int>(path.size()) - 1;
    int nRings = isClosedLoop ? nSegments : nSegments + 1;  // For closed loop, last ring = first ring

    // Generate vertices: rings of points around wire at each path point
    vertices_.reserve(nRings * numAroundWire);

    for (int ring = 0; ring < nRings; ++ring) {
        const TVector3d& p = path[ring];

        // Compute tangent direction
        TVector3d tangent;
        if (ring == 0) {
            if (isClosedLoop) {
                // For closed loop, use path[n-1] and path[1] for tangent
                tangent = path[1] - path[nSegments - 1];
            } else {
                tangent = path[1] - path[0];
            }
        } else if (ring == nRings - 1) {
            if (isClosedLoop) {
                tangent = path[0] - path[ring - 1];
            } else {
                tangent = path[nRings - 1] - path[nRings - 2];
            }
        } else {
            tangent = path[ring + 1] - path[ring - 1];
        }

        double tLen = tangent.Abs();
        if (tLen > 1e-15) {
            tangent = tangent / tLen;
        } else {
            tangent = TVector3d(1, 0, 0);
        }

        // Find perpendicular directions
        TVector3d perp1, perp2;
        if (std::abs(tangent.z) < 0.9) {
            perp1 = TVector3d(-tangent.y, tangent.x, 0);
        } else {
            perp1 = TVector3d(0, -tangent.z, tangent.y);
        }
        perp1 = perp1 / perp1.Abs();
        perp2 = tangent ^ perp1;

        // Generate ring vertices
        for (int i = 0; i < numAroundWire; ++i) {
            double theta = 2.0 * RadConst::PI * i / numAroundWire;
            double cosT = std::cos(theta);
            double sinT = std::sin(theta);

            TVector3d v = p + wireRadius * (cosT * perp1 + sinT * perp2);
            vertices_.push_back(v);
        }
    }

    // Generate triangles connecting adjacent rings
    int numTriSegments = isClosedLoop ? nRings : nSegments;
    triangles_.reserve(numTriSegments * numAroundWire * 2);

    for (int seg = 0; seg < numTriSegments; ++seg) {
        int ring0 = seg * numAroundWire;
        int ring1 = ((seg + 1) % nRings) * numAroundWire;  // Wrap around for closed loop

        for (int i = 0; i < numAroundWire; ++i) {
            int i_next = (i + 1) % numAroundWire;

            // Two triangles per quad
            // Triangle 1: ring0[i], ring1[i], ring1[i_next]
            RWGTriangle tri1;
            tri1.vertices[0] = ring0 + i;
            tri1.vertices[1] = ring1 + i;
            tri1.vertices[2] = ring1 + i_next;
            triangles_.push_back(tri1);

            // Triangle 2: ring0[i], ring1[i_next], ring0[i_next]
            RWGTriangle tri2;
            tri2.vertices[0] = ring0 + i;
            tri2.vertices[1] = ring1 + i_next;
            tri2.vertices[2] = ring0 + i_next;
            triangles_.push_back(tri2);
        }
    }

    ComputeTriangleProperties();
    BuildEdgeConnectivity();
    ComputeEdgeProperties();
}

void RWGMesh::CreateRectangularPlate(const TVector3d& center,
                                      double Lx, double Ly,
                                      const TVector3d& normal,
                                      int nx, int ny) {
    vertices_.clear();
    triangles_.clear();

    // Create local coordinate system
    TVector3d n = normal;
    double nLen = n.Abs();
    if (nLen > 1e-15) {
        n = n / nLen;
    } else {
        n = TVector3d(0, 0, 1);
    }

    TVector3d u, v;
    if (std::abs(n.z) < 0.9) {
        u = TVector3d(-n.y, n.x, 0);
    } else {
        u = TVector3d(0, -n.z, n.y);
    }
    u = u / u.Abs();
    v = n ^ u;

    // Generate vertices: (nx+1) x (ny+1) grid
    int nVx = nx + 1;
    int nVy = ny + 1;
    vertices_.reserve(nVx * nVy);

    double dx = Lx / nx;
    double dy = Ly / ny;

    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            double xi = -Lx / 2.0 + i * dx;
            double yi = -Ly / 2.0 + j * dy;
            TVector3d p = center + xi * u + yi * v;
            vertices_.push_back(p);
        }
    }

    // Generate triangles: 2 triangles per cell
    triangles_.reserve(2 * nx * ny);

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int v00 = j * nVx + i;
            int v10 = j * nVx + (i + 1);
            int v01 = (j + 1) * nVx + i;
            int v11 = (j + 1) * nVx + (i + 1);

            // Triangle 1: v00, v10, v11
            RWGTriangle tri1;
            tri1.vertices[0] = v00;
            tri1.vertices[1] = v10;
            tri1.vertices[2] = v11;
            triangles_.push_back(tri1);

            // Triangle 2: v00, v11, v01
            RWGTriangle tri2;
            tri2.vertices[0] = v00;
            tri2.vertices[1] = v11;
            tri2.vertices[2] = v01;
            triangles_.push_back(tri2);
        }
    }

    ComputeTriangleProperties();
    BuildEdgeConnectivity();
    ComputeEdgeProperties();
}

void RWGMesh::CreateCircularDisk(const TVector3d& center,
                                  double radius,
                                  const TVector3d& normal,
                                  int nr, int ntheta) {
    vertices_.clear();
    triangles_.clear();

    // Create local coordinate system
    TVector3d n = normal;
    double nLen = n.Abs();
    if (nLen > 1e-15) {
        n = n / nLen;
    } else {
        n = TVector3d(0, 0, 1);
    }

    TVector3d u, v;
    if (std::abs(n.z) < 0.9) {
        u = TVector3d(-n.y, n.x, 0);
    } else {
        u = TVector3d(0, -n.z, n.y);
    }
    u = u / u.Abs();
    v = n ^ u;

    // Generate vertices:
    // Center vertex + nr rings of ntheta vertices each
    vertices_.reserve(1 + nr * ntheta);

    // Center vertex
    vertices_.push_back(center);

    // Ring vertices
    double dr = radius / nr;
    double dtheta = 2.0 * RadConst::PI / ntheta;

    for (int r = 1; r <= nr; ++r) {
        double ri = r * dr;
        for (int t = 0; t < ntheta; ++t) {
            double theta = t * dtheta;
            double cosT = std::cos(theta);
            double sinT = std::sin(theta);
            TVector3d p = center + ri * (cosT * u + sinT * v);
            vertices_.push_back(p);
        }
    }

    // Generate triangles
    // Inner triangles: connect center to first ring
    for (int t = 0; t < ntheta; ++t) {
        int t_next = (t + 1) % ntheta;
        RWGTriangle tri;
        tri.vertices[0] = 0;  // Center
        tri.vertices[1] = 1 + t;
        tri.vertices[2] = 1 + t_next;
        triangles_.push_back(tri);
    }

    // Ring triangles: connect adjacent rings
    for (int r = 1; r < nr; ++r) {
        int ring0_start = 1 + (r - 1) * ntheta;
        int ring1_start = 1 + r * ntheta;

        for (int t = 0; t < ntheta; ++t) {
            int t_next = (t + 1) % ntheta;

            // Triangle 1: ring0[t], ring1[t], ring1[t_next]
            RWGTriangle tri1;
            tri1.vertices[0] = ring0_start + t;
            tri1.vertices[1] = ring1_start + t;
            tri1.vertices[2] = ring1_start + t_next;
            triangles_.push_back(tri1);

            // Triangle 2: ring0[t], ring1[t_next], ring0[t_next]
            RWGTriangle tri2;
            tri2.vertices[0] = ring0_start + t;
            tri2.vertices[1] = ring1_start + t_next;
            tri2.vertices[2] = ring0_start + t_next;
            triangles_.push_back(tri2);
        }
    }

    ComputeTriangleProperties();
    BuildEdgeConnectivity();
    ComputeEdgeProperties();
}

void RWGMesh::CreateCylindricalShell(const TVector3d& center,
                                      double radius,
                                      double height,
                                      const TVector3d& axis,
                                      int nz, int ntheta) {
    vertices_.clear();
    triangles_.clear();

    // Create local coordinate system
    TVector3d a = axis;
    double aLen = a.Abs();
    if (aLen > 1e-15) {
        a = a / aLen;
    } else {
        a = TVector3d(0, 0, 1);
    }

    TVector3d u, v;
    if (std::abs(a.z) < 0.9) {
        u = TVector3d(-a.y, a.x, 0);
    } else {
        u = TVector3d(0, -a.z, a.y);
    }
    u = u / u.Abs();
    v = a ^ u;

    // Generate vertices: (nz+1) rings of ntheta vertices each
    int nRings = nz + 1;
    vertices_.reserve(nRings * ntheta);

    double dz = height / nz;
    double dtheta = 2.0 * RadConst::PI / ntheta;

    TVector3d baseCenter = center - (height / 2.0) * a;

    for (int z = 0; z <= nz; ++z) {
        TVector3d ringCenter = baseCenter + z * dz * a;
        for (int t = 0; t < ntheta; ++t) {
            double theta = t * dtheta;
            double cosT = std::cos(theta);
            double sinT = std::sin(theta);
            TVector3d p = ringCenter + radius * (cosT * u + sinT * v);
            vertices_.push_back(p);
        }
    }

    // Generate triangles: 2 triangles per quad
    triangles_.reserve(2 * nz * ntheta);

    for (int z = 0; z < nz; ++z) {
        int ring0_start = z * ntheta;
        int ring1_start = (z + 1) * ntheta;

        for (int t = 0; t < ntheta; ++t) {
            int t_next = (t + 1) % ntheta;

            // Triangle 1: ring0[t], ring1[t], ring1[t_next]
            RWGTriangle tri1;
            tri1.vertices[0] = ring0_start + t;
            tri1.vertices[1] = ring1_start + t;
            tri1.vertices[2] = ring1_start + t_next;
            triangles_.push_back(tri1);

            // Triangle 2: ring0[t], ring1[t_next], ring0[t_next]
            RWGTriangle tri2;
            tri2.vertices[0] = ring0_start + t;
            tri2.vertices[1] = ring1_start + t_next;
            tri2.vertices[2] = ring0_start + t_next;
            triangles_.push_back(tri2);
        }
    }

    ComputeTriangleProperties();
    BuildEdgeConnectivity();
    ComputeEdgeProperties();
}

void RWGMesh::MergeMesh(const RWGMesh& other) {
    // Offset for new vertex indices
    int vertexOffset = static_cast<int>(vertices_.size());

    // Add vertices
    for (const auto& v : other.vertices_) {
        vertices_.push_back(v);
    }

    // Add triangles with offset indices
    for (const auto& tri : other.triangles_) {
        RWGTriangle newTri;
        newTri.vertices[0] = tri.vertices[0] + vertexOffset;
        newTri.vertices[1] = tri.vertices[1] + vertexOffset;
        newTri.vertices[2] = tri.vertices[2] + vertexOffset;
        triangles_.push_back(newTri);
    }

    // Rebuild connectivity
    ComputeTriangleProperties();
    BuildEdgeConnectivity();
    ComputeEdgeProperties();
}

std::vector<int> RWGMesh::GetBoundaryEdges() const {
    std::vector<int> boundaryEdges;
    for (int i = 0; i < static_cast<int>(edges_.size()); ++i) {
        if (edges_[i].isBoundary) {
            boundaryEdges.push_back(i);
        }
    }
    return boundaryEdges;
}

void RWGMesh::ComputeTriangleProperties() {
    for (auto& tri : triangles_) {
        const TVector3d& v0 = vertices_[tri.vertices[0]];
        const TVector3d& v1 = vertices_[tri.vertices[1]];
        const TVector3d& v2 = vertices_[tri.vertices[2]];

        // Centroid
        tri.center = (v0 + v1 + v2) / 3.0;

        // Edge vectors
        TVector3d e1 = v1 - v0;
        TVector3d e2 = v2 - v0;

        // Normal (cross product)
        TVector3d n = e1 ^ e2;
        double area2 = n.Abs();

        tri.area = area2 / 2.0;

        if (area2 > 1e-30) {
            tri.normal = n / area2;
        } else {
            tri.normal = TVector3d(0, 0, 1);
        }
    }
}

void RWGMesh::BuildEdgeConnectivity() {
    edges_.clear();

    // Use hash map for efficient edge lookup
    std::unordered_map<uint64_t, int> edgeMap;

    for (int t = 0; t < static_cast<int>(triangles_.size()); ++t) {
        auto& tri = triangles_[t];

        // Process 3 edges of this triangle
        for (int localEdge = 0; localEdge < 3; ++localEdge) {
            // Edge i is opposite to vertex i
            // Edge connects vertices (i+1)%3 and (i+2)%3
            int v1 = tri.vertices[(localEdge + 1) % 3];
            int v2 = tri.vertices[(localEdge + 2) % 3];

            uint64_t hash = EdgeHash(v1, v2);

            auto it = edgeMap.find(hash);
            if (it == edgeMap.end()) {
                // New edge
                RWGEdge edge;
                edge.vertex1 = std::min(v1, v2);
                edge.vertex2 = std::max(v1, v2);
                edge.triPlus = t;
                edge.localIdxPlus = localEdge;
                edge.triMinus = -1;  // Not yet found
                edge.localIdxMinus = -1;
                edge.isBoundary = true;

                int edgeIdx = static_cast<int>(edges_.size());
                edges_.push_back(edge);
                edgeMap[hash] = edgeIdx;

                tri.edges[localEdge] = edgeIdx;
                // Sign: +1 if edge direction (v1->v2) matches our convention
                tri.edgeSign[localEdge] = (v1 < v2) ? 1 : -1;
            } else {
                // Existing edge - this triangle is T-
                int edgeIdx = it->second;
                auto& edge = edges_[edgeIdx];
                edge.triMinus = t;
                edge.localIdxMinus = localEdge;
                edge.isBoundary = false;

                tri.edges[localEdge] = edgeIdx;
                tri.edgeSign[localEdge] = (v1 < v2) ? -1 : 1;
            }
        }
    }

    // Count interior edges
    numInteriorEdges_ = 0;
    for (const auto& edge : edges_) {
        if (!edge.isBoundary) {
            numInteriorEdges_++;
        }
    }
}

void RWGMesh::ComputeEdgeProperties() {
    for (auto& edge : edges_) {
        const TVector3d& v1 = vertices_[edge.vertex1];
        const TVector3d& v2 = vertices_[edge.vertex2];

        edge.midpoint = (v1 + v2) / 2.0;

        TVector3d dir = v2 - v1;
        edge.length = dir.Abs();

        if (edge.length > 1e-15) {
            edge.direction = dir / edge.length;
        } else {
            edge.direction = TVector3d(1, 0, 0);
        }

        // Compute rho vectors for RWG basis evaluation
        if (edge.triPlus >= 0) {
            // Free vertex in T+ (vertex opposite to edge)
            const auto& triP = triangles_[edge.triPlus];
            int freeVertIdxP = triP.vertices[edge.localIdxPlus];
            edge.rhoPlusC = triP.center - vertices_[freeVertIdxP];
        }

        if (edge.triMinus >= 0) {
            // Free vertex in T-
            const auto& triM = triangles_[edge.triMinus];
            int freeVertIdxM = triM.vertices[edge.localIdxMinus];
            edge.rhoMinusC = vertices_[freeVertIdxM] - triM.center;
        }
    }
}

TVector3d RWGMesh::EvaluateRWG(int edgeIdx, int triIdx, const TVector3d& point) const {
    const auto& edge = edges_[edgeIdx];

    if (triIdx == edge.triPlus) {
        // In T+: f = (l/2A) * rho^+ = (l/2A) * (r - r_free)
        const auto& tri = triangles_[triIdx];
        int freeVertIdx = tri.vertices[edge.localIdxPlus];
        const TVector3d& rFree = vertices_[freeVertIdx];

        TVector3d rho = point - rFree;
        double coeff = edge.length / (2.0 * tri.area);

        return coeff * rho;
    }
    else if (triIdx == edge.triMinus) {
        // In T-: f = (l/2A) * rho^- = (l/2A) * (r_free - r)
        const auto& tri = triangles_[triIdx];
        int freeVertIdx = tri.vertices[edge.localIdxMinus];
        const TVector3d& rFree = vertices_[freeVertIdx];

        TVector3d rho = rFree - point;
        double coeff = edge.length / (2.0 * tri.area);

        return coeff * rho;
    }

    // Point not in either triangle
    return TVector3d(0, 0, 0);
}

double RWGMesh::EvaluateRWGDiv(int edgeIdx, int triIdx) const {
    const auto& edge = edges_[edgeIdx];

    if (triIdx == edge.triPlus) {
        // div(f) = l/A in T+
        return edge.length / triangles_[triIdx].area;
    }
    else if (triIdx == edge.triMinus) {
        // div(f) = -l/A in T-
        return -edge.length / triangles_[triIdx].area;
    }

    return 0.0;
}

TVector3d RWGMesh::GetFreeVertex(int edgeIdx, int triIdx) const {
    const auto& edge = edges_[edgeIdx];

    if (triIdx == edge.triPlus) {
        const auto& tri = triangles_[triIdx];
        return vertices_[tri.vertices[edge.localIdxPlus]];
    }
    else if (triIdx == edge.triMinus) {
        const auto& tri = triangles_[triIdx];
        return vertices_[tri.vertices[edge.localIdxMinus]];
    }

    return TVector3d(0, 0, 0);
}

// ============================================================================
// RWGIntegration Implementation
// ============================================================================

double RWGIntegration::SignedSolidAngle(const TVector3d& point,
                                         const TVector3d& v1,
                                         const TVector3d& v2,
                                         const TVector3d& v3) {
    // Van Oosterom & Strackee formula for solid angle
    TVector3d a = v1 - point;
    TVector3d b = v2 - point;
    TVector3d c = v3 - point;

    double la = a.Abs();
    double lb = b.Abs();
    double lc = c.Abs();

    if (la < 1e-15 || lb < 1e-15 || lc < 1e-15) {
        return 0.0;  // Point is at a vertex
    }

    // Numerator: triple product a . (b x c)
    double num = a.x * (b.y * c.z - b.z * c.y) +
                 a.y * (b.z * c.x - b.x * c.z) +
                 a.z * (b.x * c.y - b.y * c.x);

    // Denominator
    double ab = a.x * b.x + a.y * b.y + a.z * b.z;
    double bc = b.x * c.x + b.y * c.y + b.z * c.z;
    double ca = c.x * a.x + c.y * a.y + c.z * a.z;

    double denom = la * lb * lc + ab * lc + bc * la + ca * lb;

    if (std::abs(denom) < 1e-30) {
        return 0.0;
    }

    return 2.0 * std::atan2(num, denom);
}

double RWGIntegration::LineIntegralTerm(const TVector3d& point,
                                         const TVector3d& v1,
                                         const TVector3d& v2,
                                         double& R_plus, double& R_minus) {
    // Line integral term from Wilton et al. (1984)
    // For edge from v1 to v2

    TVector3d edge = v2 - v1;
    double l = edge.Abs();
    if (l < 1e-15) {
        R_plus = R_minus = 0;
        return 0;
    }

    TVector3d unit = edge / l;

    // Project point onto line
    TVector3d r1 = v1 - point;
    TVector3d r2 = v2 - point;

    // Distance from point to line endpoints
    R_minus = r1.Abs();
    R_plus = r2.Abs();

    // Signed distances along edge direction
    double t1 = -(r1.x * unit.x + r1.y * unit.y + r1.z * unit.z);
    double t2 = -(r2.x * unit.x + r2.y * unit.y + r2.z * unit.z);

    // Perpendicular distance squared
    double d_perp_sq = R_minus * R_minus - t1 * t1;
    if (d_perp_sq < 0) d_perp_sq = 0;
    double d_perp = std::sqrt(d_perp_sq);

    // Line integral: ln((R+ + l+) / (R- + l-)) where l± are signed projections
    double arg1 = R_plus + t2;
    double arg2 = R_minus + t1;

    if (arg1 < 1e-15) arg1 = 1e-15;
    if (arg2 < 1e-15) arg2 = 1e-15;

    return std::log(arg1 / arg2);
}

double RWGIntegration::PotentialIntegral(const TVector3d& point,
                                          const TVector3d& v1,
                                          const TVector3d& v2,
                                          const TVector3d& v3) {
    // Analytical potential integral using Wilton et al. formula
    //
    // I = integral_T { 1/R } dA
    //
    // Formula: I = sum over edges of { P0_i * beta_i } - |d| * Omega
    // where:
    //   P0_i = perpendicular distance from point to edge i (in triangle plane)
    //   beta_i = atan2(...) term involving distances to edge endpoints
    //   d = signed distance from point to triangle plane
    //   Omega = solid angle subtended by triangle

    // Triangle normal and area
    TVector3d e1 = v2 - v1;
    TVector3d e2 = v3 - v1;
    TVector3d n = e1 ^ e2;
    double area2 = n.Abs();

    if (area2 < 1e-30) {
        return 0.0;  // Degenerate triangle
    }

    TVector3d normal = n / area2;

    // Signed distance to plane
    TVector3d r0 = point - v1;
    double d = r0.x * normal.x + r0.y * normal.y + r0.z * normal.z;
    double abs_d = std::abs(d);

    // Solid angle term
    double Omega = SignedSolidAngle(point, v1, v2, v3);

    // Project point onto triangle plane
    TVector3d p_proj = point - d * normal;

    // Edge contributions
    double I = 0.0;

    TVector3d edges[3][2] = {
        {v1, v2}, {v2, v3}, {v3, v1}
    };

    for (int i = 0; i < 3; ++i) {
        const TVector3d& ea = edges[i][0];
        const TVector3d& eb = edges[i][1];

        // Edge direction
        TVector3d edge_vec = eb - ea;
        double edge_len = edge_vec.Abs();
        if (edge_len < 1e-15) continue;

        TVector3d edge_unit = edge_vec / edge_len;

        // Outward normal in plane (perpendicular to edge, in triangle plane)
        TVector3d edge_normal = edge_unit ^ normal;

        // Perpendicular distance from projected point to edge line
        TVector3d to_edge = ea - p_proj;
        double P0 = to_edge.x * edge_normal.x + to_edge.y * edge_normal.y + to_edge.z * edge_normal.z;

        // Distances from point to edge endpoints
        double R_plus, R_minus;
        double ln_term = LineIntegralTerm(point, ea, eb, R_plus, R_minus);

        // Projection distances along edge
        TVector3d r_a = ea - point;
        TVector3d r_b = eb - point;
        double t_a = r_a.x * edge_unit.x + r_a.y * edge_unit.y + r_a.z * edge_unit.z;
        double t_b = r_b.x * edge_unit.x + r_b.y * edge_unit.y + r_b.z * edge_unit.z;

        // Beta term
        double beta = std::atan2(P0 * t_b, P0 * P0 + abs_d * R_plus) -
                      std::atan2(P0 * t_a, P0 * P0 + abs_d * R_minus);

        I += P0 * ln_term - abs_d * beta;
    }

    return I;
}

TVector3d RWGIntegration::GradPotentialIntegral(const TVector3d& point,
                                                 const TVector3d& v1,
                                                 const TVector3d& v2,
                                                 const TVector3d& v3) {
    // Gradient of potential integral
    // grad I = integral_T { (r - r') / R^3 } dA = Omega * n
    // where n is triangle normal and Omega is solid angle

    // Triangle normal
    TVector3d e1 = v2 - v1;
    TVector3d e2 = v3 - v1;
    TVector3d n = e1 ^ e2;
    double area2 = n.Abs();

    if (area2 < 1e-30) {
        return TVector3d(0, 0, 0);
    }

    TVector3d normal = n / area2;

    // Solid angle
    double Omega = SignedSolidAngle(point, v1, v2, v3);

    // Sign depends on which side of plane the point is
    TVector3d r0 = point - v1;
    double d = r0.x * normal.x + r0.y * normal.y + r0.z * normal.z;

    if (d < 0) {
        Omega = -Omega;
    }

    return Omega * normal;
}

std::complex<double> RWGIntegration::VectorPotentialIntegral(
    const RWGMesh& mesh, int edgeM, int edgeN) {
    // L_mn = μ₀/(4π) ∫∫ f_m(r) · f_n(r') / |r-r'| dA dA'
    //
    // f_m spans T_m^+ and T_m^-
    // f_n spans T_n^+ and T_n^-
    //
    // Total: 4 pairs of triangle integrations (or 2 if m == n)

    const auto& edgeM_data = mesh.Edge(edgeM);
    const auto& edgeN_data = mesh.Edge(edgeN);

    std::complex<double> result(0, 0);

    // Integration over observation triangles (f_m)
    int triM[2] = {edgeM_data.triPlus, edgeM_data.triMinus};
    int triN[2] = {edgeN_data.triPlus, edgeN_data.triMinus};

    for (int im = 0; im < 2; ++im) {
        if (triM[im] < 0) continue;
        const auto& triM_data = mesh.Triangle(triM[im]);

        for (int in = 0; in < 2; ++in) {
            if (triN[in] < 0) continue;
            const auto& triN_data = mesh.Triangle(triN[in]);

            // Get triangle vertices
            const TVector3d& mV0 = mesh.Vertex(triM_data.vertices[0]);
            const TVector3d& mV1 = mesh.Vertex(triM_data.vertices[1]);
            const TVector3d& mV2 = mesh.Vertex(triM_data.vertices[2]);

            const TVector3d& nV0 = mesh.Vertex(triN_data.vertices[0]);
            const TVector3d& nV1 = mesh.Vertex(triN_data.vertices[1]);
            const TVector3d& nV2 = mesh.Vertex(triN_data.vertices[2]);

            // Check if same triangle (self-term)
            if (triM[im] == triN[in]) {
                // Use analytical self-term
                result += SelfTermIntegral(mesh, edgeM, triM[im], true);
            } else {
                // Numerical quadrature
                double integral = 0.0;

                for (int qm = 0; qm < NQUAD_TRI; ++qm) {
                    // Observation point in triM
                    double L1m = triQuadPts_[qm][0];
                    double L2m = triQuadPts_[qm][1];
                    double L3m = 1.0 - L1m - L2m;

                    TVector3d r_obs = L1m * mV0 + L2m * mV1 + L3m * mV2;

                    // Evaluate RWG basis at observation point
                    TVector3d f_m = mesh.EvaluateRWG(edgeM, triM[im], r_obs);

                    for (int qn = 0; qn < NQUAD_TRI; ++qn) {
                        // Source point in triN
                        double L1n = triQuadPts_[qn][0];
                        double L2n = triQuadPts_[qn][1];
                        double L3n = 1.0 - L1n - L2n;

                        TVector3d r_src = L1n * nV0 + L2n * nV1 + L3n * nV2;

                        // Distance
                        TVector3d R_vec = r_obs - r_src;
                        double R = R_vec.Abs();

                        if (R < 1e-15) continue;

                        // Evaluate RWG basis at source point
                        TVector3d f_n = mesh.EvaluateRWG(edgeN, triN[in], r_src);

                        // Dot product
                        double fm_dot_fn = f_m.x * f_n.x + f_m.y * f_n.y + f_m.z * f_n.z;

                        // Accumulate
                        integral += triQuadWts_[qm] * triQuadWts_[qn] *
                                    fm_dot_fn / R;
                    }
                }

                // Scale by triangle areas
                integral *= triM_data.area * triN_data.area;

                result += std::complex<double>(integral, 0);
            }
        }
    }

    // Scale by μ₀/(4π)
    result *= RWG_MU_0 * RWG_INV_FOUR_PI;

    return result;
}

std::complex<double> RWGIntegration::ScalarPotentialIntegral(
    const RWGMesh& mesh, int edgeM, int edgeN) {
    // P_mn = 1/(4πε₀) ∫∫ div(f_m) * div(f_n) / |r-r'| dA dA'
    //
    // div(f) is constant in each triangle:
    //   div(f) = ±l/A

    const auto& edgeM_data = mesh.Edge(edgeM);
    const auto& edgeN_data = mesh.Edge(edgeN);

    std::complex<double> result(0, 0);

    int triM[2] = {edgeM_data.triPlus, edgeM_data.triMinus};
    int triN[2] = {edgeN_data.triPlus, edgeN_data.triMinus};

    for (int im = 0; im < 2; ++im) {
        if (triM[im] < 0) continue;
        const auto& triM_data = mesh.Triangle(triM[im]);
        double divM = mesh.EvaluateRWGDiv(edgeM, triM[im]);

        for (int in = 0; in < 2; ++in) {
            if (triN[in] < 0) continue;
            const auto& triN_data = mesh.Triangle(triN[in]);
            double divN = mesh.EvaluateRWGDiv(edgeN, triN[in]);

            // div product (constant)
            double divM_divN = divM * divN;

            // Get triangle vertices
            const TVector3d& mV0 = mesh.Vertex(triM_data.vertices[0]);
            const TVector3d& mV1 = mesh.Vertex(triM_data.vertices[1]);
            const TVector3d& mV2 = mesh.Vertex(triM_data.vertices[2]);

            const TVector3d& nV0 = mesh.Vertex(triN_data.vertices[0]);
            const TVector3d& nV1 = mesh.Vertex(triN_data.vertices[1]);
            const TVector3d& nV2 = mesh.Vertex(triN_data.vertices[2]);

            // Self-term check
            if (triM[im] == triN[in]) {
                // Use analytical formula
                result += SelfTermIntegral(mesh, edgeM, triM[im], false);
            } else {
                // Numerical quadrature for 1/R over two triangles
                double integral = 0.0;

                for (int qm = 0; qm < NQUAD_TRI; ++qm) {
                    double L1m = triQuadPts_[qm][0];
                    double L2m = triQuadPts_[qm][1];
                    double L3m = 1.0 - L1m - L2m;

                    TVector3d r_obs = L1m * mV0 + L2m * mV1 + L3m * mV2;

                    for (int qn = 0; qn < NQUAD_TRI; ++qn) {
                        double L1n = triQuadPts_[qn][0];
                        double L2n = triQuadPts_[qn][1];
                        double L3n = 1.0 - L1n - L2n;

                        TVector3d r_src = L1n * nV0 + L2n * nV1 + L3n * nV2;

                        double R = (r_obs - r_src).Abs();
                        if (R < 1e-15) continue;

                        integral += triQuadWts_[qm] * triQuadWts_[qn] / R;
                    }
                }

                integral *= triM_data.area * triN_data.area * divM_divN;
                result += std::complex<double>(integral, 0);
            }
        }
    }

    // Scale by 1/(4π) - note: ε₀ factor is applied in Z matrix assembly
    result *= RWG_INV_FOUR_PI;

    return result;
}

std::complex<double> RWGIntegration::SelfTermIntegral(
    const RWGMesh& mesh, int edgeIdx, int triIdx, bool isVectorPotential) {
    // Self-term: source and observation in same triangle
    //
    // Use analytical formula from Wilton et al. (1984):
    // I_self = integral_T integral_T { 1/R } dA dA'
    //
    // For equilateral triangle of side a:
    //   I_self = (sqrt(3)/2) * a^3 * ln(3)  ≈ 0.952 * a^3
    //
    // General formula involves logarithmic terms for each edge

    const auto& tri = mesh.Triangle(triIdx);
    const TVector3d& v0 = mesh.Vertex(tri.vertices[0]);
    const TVector3d& v1 = mesh.Vertex(tri.vertices[1]);
    const TVector3d& v2 = mesh.Vertex(tri.vertices[2]);

    double A = tri.area;

    // Approximate self-integral using characteristic length
    // I ≈ 4 * A * ln(L/a) where L = sqrt(A), a = wire radius (set to sqrt(A)/10)
    double L = std::sqrt(A);
    double a_eff = L / 10.0;  // Effective wire radius

    double I_self = 4.0 * A * std::log(2.0 * L / a_eff);

    if (isVectorPotential) {
        // For vector potential, scale by RWG basis factor
        const auto& edge = mesh.Edge(edgeIdx);
        double l = edge.length;
        double coeff = l * l / (4.0 * A * A);  // (l/2A)^2 * A for integration

        // f_m · f_n in same triangle: roughly (l/2A)^2 * |rho|^2
        // Average |rho|^2 ≈ (2/3) * h^2 where h = altitude
        double h_avg = 2.0 * A / l;  // Approximate altitude
        double rho_sq_avg = (2.0 / 3.0) * h_avg * h_avg;

        return std::complex<double>(RWG_MU_0 * RWG_INV_FOUR_PI * coeff * rho_sq_avg * I_self, 0);
    } else {
        // For scalar potential: div(f_m) * div(f_n) = (l/A)^2
        const auto& edge = mesh.Edge(edgeIdx);
        double l = edge.length;
        double divF_sq = (l / A) * (l / A);

        return std::complex<double>(RWG_INV_FOUR_PI / RWG_EPS_0 * divF_sq * A * A * I_self / A, 0);
    }
}

// ============================================================================
// LoopStarDecomposition Implementation
// ============================================================================

LoopStarDecomposition::LoopStarDecomposition()
    : numEdges_(0)
    , numVertices_(0)
    , numLoops_(0)
    , numStars_(0)
    , isValid_(false)
{
}

LoopStarDecomposition::~LoopStarDecomposition() {
}

void LoopStarDecomposition::BuildSpanningTree(const RWGMesh& mesh,
                                               std::vector<int>& parent,
                                               std::vector<int>& treeEdges,
                                               std::vector<int>& coTreeEdges) {
    // Build spanning tree using BFS
    // Tree edges: edges in the spanning tree
    // Co-tree edges: edges not in the spanning tree (these form loops)

    int nVertices = mesh.NumVertices();
    int nEdges = mesh.NumEdges();

    parent.assign(nVertices, -1);
    treeEdges.clear();
    coTreeEdges.clear();

    // Visited flag
    std::vector<bool> visited(nVertices, false);

    // BFS starting from vertex 0
    std::vector<int> queue;
    queue.push_back(0);
    visited[0] = true;
    parent[0] = 0;  // Root is its own parent

    size_t qHead = 0;
    while (qHead < queue.size()) {
        int v = queue[qHead++];

        // Find all edges connected to vertex v
        for (int e = 0; e < nEdges; ++e) {
            const auto& edge = mesh.Edge(e);
            if (edge.isBoundary) continue;  // Only interior edges

            int other = -1;
            if (edge.vertex1 == v) {
                other = edge.vertex2;
            } else if (edge.vertex2 == v) {
                other = edge.vertex1;
            } else {
                continue;  // Edge not connected to v
            }

            if (!visited[other]) {
                // Tree edge
                visited[other] = true;
                parent[other] = v;
                treeEdges.push_back(e);
                queue.push_back(other);
            }
        }
    }

    // Co-tree edges are all interior edges not in tree
    std::vector<bool> isTreeEdge(nEdges, false);
    for (int e : treeEdges) {
        isTreeEdge[e] = true;
    }

    for (int e = 0; e < nEdges; ++e) {
        if (!mesh.Edge(e).isBoundary && !isTreeEdge[e]) {
            coTreeEdges.push_back(e);
        }
    }
}

void LoopStarDecomposition::Build(const RWGMesh& mesh) {
    isValid_ = false;

    numEdges_ = mesh.NumInteriorEdges();
    numVertices_ = mesh.NumVertices();

    if (numEdges_ == 0 || numVertices_ == 0) {
        return;
    }

    // Build spanning tree
    std::vector<int> parent;
    std::vector<int> treeEdges;
    std::vector<int> coTreeEdges;
    BuildSpanningTree(mesh, parent, treeEdges, coTreeEdges);

    // Number of loops = number of co-tree edges
    // Each co-tree edge, when added to the spanning tree, creates exactly one loop
    numLoops_ = static_cast<int>(coTreeEdges.size());

    // Number of stars = number of vertices - 1 (one vertex is grounded)
    numStars_ = numVertices_ - 1;

    // Create mapping from edge index to interior edge index
    std::vector<int> edgeToInterior(mesh.NumEdges(), -1);
    int interiorIdx = 0;
    for (int e = 0; e < mesh.NumEdges(); ++e) {
        if (!mesh.Edge(e).isBoundary) {
            edgeToInterior[e] = interiorIdx++;
        }
    }

    // Build T_loop: numEdges_ x numLoops_
    // Each column corresponds to a loop (one co-tree edge)
    // The loop is the path from vertex1 to vertex2 of the co-tree edge through the tree,
    // plus the co-tree edge itself

    T_loop_.resize(numEdges_ * numLoops_, 0.0);

    for (int loopIdx = 0; loopIdx < numLoops_; ++loopIdx) {
        int coTreeEdge = coTreeEdges[loopIdx];
        const auto& edge = mesh.Edge(coTreeEdge);

        // Find path from vertex1 to vertex2 through the tree
        std::vector<int> path1, path2;

        // Path from vertex1 to root
        int v = edge.vertex1;
        while (v != parent[v]) {
            path1.push_back(v);
            v = parent[v];
        }
        path1.push_back(v);  // Root

        // Path from vertex2 to root
        v = edge.vertex2;
        while (v != parent[v]) {
            path2.push_back(v);
            v = parent[v];
        }
        path2.push_back(v);  // Root

        // Find common ancestor (LCA)
        std::vector<bool> onPath1(numVertices_, false);
        for (int pv : path1) {
            onPath1[pv] = true;
        }

        int lca = -1;
        for (int pv : path2) {
            if (onPath1[pv]) {
                lca = pv;
                break;
            }
        }

        if (lca < 0) {
            // Should not happen in connected graph
            continue;
        }

        // Loop edges: path from vertex1 to LCA + path from LCA to vertex2 + co-tree edge
        // Edge signs: +1 if edge direction matches loop direction, -1 otherwise

        // Path vertex1 -> LCA (going up the tree)
        for (size_t i = 0; i < path1.size() - 1; ++i) {
            if (path1[i + 1] == lca || onPath1[path1[i + 1]]) {
                // Find edge between path1[i] and path1[i+1]
                int v1 = path1[i];
                int v2 = path1[i + 1];

                for (int te : treeEdges) {
                    const auto& tedge = mesh.Edge(te);
                    if ((tedge.vertex1 == v1 && tedge.vertex2 == v2) ||
                        (tedge.vertex1 == v2 && tedge.vertex2 == v1)) {
                        int intIdx = edgeToInterior[te];
                        if (intIdx >= 0) {
                            // Sign: +1 if going from v1 to v2 matches edge direction
                            double sign = (tedge.vertex1 == v1) ? 1.0 : -1.0;
                            T_loop_[intIdx * numLoops_ + loopIdx] = sign;
                        }
                        break;
                    }
                }
            }
            if (path1[i + 1] == lca) break;
        }

        // Path LCA -> vertex2 (going down the tree, so reverse direction)
        for (size_t i = 0; i < path2.size() - 1; ++i) {
            if (path2[i] == lca) break;

            int v1 = path2[i + 1];
            int v2 = path2[i];

            for (int te : treeEdges) {
                const auto& tedge = mesh.Edge(te);
                if ((tedge.vertex1 == v1 && tedge.vertex2 == v2) ||
                    (tedge.vertex1 == v2 && tedge.vertex2 == v1)) {
                    int intIdx = edgeToInterior[te];
                    if (intIdx >= 0) {
                        // Sign: +1 if going from v1 to v2 matches edge direction
                        double sign = (tedge.vertex1 == v1) ? 1.0 : -1.0;
                        T_loop_[intIdx * numLoops_ + loopIdx] = sign;
                    }
                    break;
                }
            }
        }

        // Add the co-tree edge itself with +1 sign
        int coTreeIntIdx = edgeToInterior[coTreeEdge];
        if (coTreeIntIdx >= 0) {
            T_loop_[coTreeIntIdx * numLoops_ + loopIdx] = 1.0;
        }
    }

    // Build T_star: numEdges_ x numStars_
    // Star basis: one star per vertex (except one grounded vertex)
    // Star function at vertex v: sum of all edges incident to v, with appropriate signs

    T_star_.resize(numEdges_ * numStars_, 0.0);

    // Grounded vertex is vertex 0
    for (int starIdx = 0; starIdx < numStars_; ++starIdx) {
        int vertex = starIdx + 1;  // Skip vertex 0 (grounded)

        // Find all edges incident to this vertex
        for (int e = 0; e < mesh.NumEdges(); ++e) {
            const auto& edge = mesh.Edge(e);
            if (edge.isBoundary) continue;

            int intIdx = edgeToInterior[e];
            if (intIdx < 0) continue;

            if (edge.vertex1 == vertex) {
                // Edge goes from vertex to vertex2: positive direction
                T_star_[intIdx * numStars_ + starIdx] = 1.0;
            } else if (edge.vertex2 == vertex) {
                // Edge goes from vertex1 to vertex: negative direction
                T_star_[intIdx * numStars_ + starIdx] = -1.0;
            }
        }
    }

    isValid_ = true;
}

void LoopStarDecomposition::TransformMatrix(const std::vector<std::complex<double>>& Z_RWG,
                                            std::vector<std::complex<double>>& Z_LL,
                                            std::vector<std::complex<double>>& Z_LS,
                                            std::vector<std::complex<double>>& Z_SL,
                                            std::vector<std::complex<double>>& Z_SS) const {
    if (!isValid_) return;

    // Z_LS = T^T * Z_RWG * T
    // where T = [T_loop | T_star]
    //
    // Z_LL = T_loop^T * Z_RWG * T_loop
    // Z_LS = T_loop^T * Z_RWG * T_star
    // Z_SL = T_star^T * Z_RWG * T_loop
    // Z_SS = T_star^T * Z_RWG * T_star

    Z_LL.resize(numLoops_ * numLoops_, std::complex<double>(0, 0));
    Z_LS.resize(numLoops_ * numStars_, std::complex<double>(0, 0));
    Z_SL.resize(numStars_ * numLoops_, std::complex<double>(0, 0));
    Z_SS.resize(numStars_ * numStars_, std::complex<double>(0, 0));

    // First compute Z_RWG * T_loop and Z_RWG * T_star
    std::vector<std::complex<double>> ZT_loop(numEdges_ * numLoops_, std::complex<double>(0, 0));
    std::vector<std::complex<double>> ZT_star(numEdges_ * numStars_, std::complex<double>(0, 0));

    // ZT_loop = Z_RWG * T_loop
    for (int i = 0; i < numEdges_; ++i) {
        for (int j = 0; j < numLoops_; ++j) {
            std::complex<double> sum(0, 0);
            for (int k = 0; k < numEdges_; ++k) {
                sum += Z_RWG[i * numEdges_ + k] * T_loop_[k * numLoops_ + j];
            }
            ZT_loop[i * numLoops_ + j] = sum;
        }
    }

    // ZT_star = Z_RWG * T_star
    for (int i = 0; i < numEdges_; ++i) {
        for (int j = 0; j < numStars_; ++j) {
            std::complex<double> sum(0, 0);
            for (int k = 0; k < numEdges_; ++k) {
                sum += Z_RWG[i * numEdges_ + k] * T_star_[k * numStars_ + j];
            }
            ZT_star[i * numStars_ + j] = sum;
        }
    }

    // Z_LL = T_loop^T * ZT_loop
    for (int i = 0; i < numLoops_; ++i) {
        for (int j = 0; j < numLoops_; ++j) {
            std::complex<double> sum(0, 0);
            for (int k = 0; k < numEdges_; ++k) {
                sum += T_loop_[k * numLoops_ + i] * ZT_loop[k * numLoops_ + j];
            }
            Z_LL[i * numLoops_ + j] = sum;
        }
    }

    // Z_LS = T_loop^T * ZT_star
    for (int i = 0; i < numLoops_; ++i) {
        for (int j = 0; j < numStars_; ++j) {
            std::complex<double> sum(0, 0);
            for (int k = 0; k < numEdges_; ++k) {
                sum += T_loop_[k * numLoops_ + i] * ZT_star[k * numStars_ + j];
            }
            Z_LS[i * numStars_ + j] = sum;
        }
    }

    // Z_SL = T_star^T * ZT_loop
    for (int i = 0; i < numStars_; ++i) {
        for (int j = 0; j < numLoops_; ++j) {
            std::complex<double> sum(0, 0);
            for (int k = 0; k < numEdges_; ++k) {
                sum += T_star_[k * numStars_ + i] * ZT_loop[k * numLoops_ + j];
            }
            Z_SL[i * numLoops_ + j] = sum;
        }
    }

    // Z_SS = T_star^T * ZT_star
    for (int i = 0; i < numStars_; ++i) {
        for (int j = 0; j < numStars_; ++j) {
            std::complex<double> sum(0, 0);
            for (int k = 0; k < numEdges_; ++k) {
                sum += T_star_[k * numStars_ + i] * ZT_star[k * numStars_ + j];
            }
            Z_SS[i * numStars_ + j] = sum;
        }
    }
}

void LoopStarDecomposition::TransformRHS(const std::vector<std::complex<double>>& V_RWG,
                                         std::vector<std::complex<double>>& V_L,
                                         std::vector<std::complex<double>>& V_S) const {
    if (!isValid_) return;

    // V_L = T_loop^T * V_RWG
    V_L.resize(numLoops_, std::complex<double>(0, 0));
    for (int i = 0; i < numLoops_; ++i) {
        std::complex<double> sum(0, 0);
        for (int k = 0; k < numEdges_; ++k) {
            sum += T_loop_[k * numLoops_ + i] * V_RWG[k];
        }
        V_L[i] = sum;
    }

    // V_S = T_star^T * V_RWG
    V_S.resize(numStars_, std::complex<double>(0, 0));
    for (int i = 0; i < numStars_; ++i) {
        std::complex<double> sum(0, 0);
        for (int k = 0; k < numEdges_; ++k) {
            sum += T_star_[k * numStars_ + i] * V_RWG[k];
        }
        V_S[i] = sum;
    }
}

void LoopStarDecomposition::TransformSolution(const std::vector<std::complex<double>>& I_L,
                                               const std::vector<std::complex<double>>& I_S,
                                               std::vector<std::complex<double>>& I_RWG) const {
    if (!isValid_) return;

    // I_RWG = T_loop * I_L + T_star * I_S
    I_RWG.resize(numEdges_, std::complex<double>(0, 0));

    for (int i = 0; i < numEdges_; ++i) {
        std::complex<double> sum(0, 0);

        // T_loop * I_L
        for (int j = 0; j < numLoops_; ++j) {
            sum += T_loop_[i * numLoops_ + j] * I_L[j];
        }

        // T_star * I_S
        for (int j = 0; j < numStars_; ++j) {
            sum += T_star_[i * numStars_ + j] * I_S[j];
        }

        I_RWG[i] = sum;
    }
}

// ============================================================================
// RWGEFIESolver Implementation
// ============================================================================

RWGEFIESolver::RWGEFIESolver()
    : conductivity_(5.8e7)
    , frequency_(0)
    , omega_(0)
    , jOmega_(0, 0)
    , voltageExcitation_(1.0, 0)
    , impedance_(0, 0)
    , totalCurrent_(0, 0)
    , useMQS_(true)        // MQS mode enabled by default for power electronics
{
    // Loop-Star decomposition is always used (no option to disable)
    // This ensures low-frequency stability for power electronics applications
}

RWGEFIESolver::~RWGEFIESolver() {
}

void RWGEFIESolver::SetFrequency(double freq) {
    frequency_ = freq;
    omega_ = 2.0 * RadConst::PI * freq;
    jOmega_ = std::complex<double>(0, omega_);
}

void RWGEFIESolver::DefinePort(const std::vector<int>& edgeIndices) {
    portEdges_ = edgeIndices;
}

void RWGEFIESolver::SetVoltageExcitation(std::complex<double> V) {
    voltageExcitation_ = V;
}

void RWGEFIESolver::Solve() {
    if (!mesh_ || mesh_->NumInteriorEdges() == 0) {
        throw std::runtime_error("No mesh or no interior edges");
    }

    AssembleMatrices();
    AssembleRHS();
    SolveLoopStar();  // Loop-Star decomposition is always used
    ComputeImpedance();
}

void RWGEFIESolver::AssembleMatrices() {
    int N = mesh_->NumInteriorEdges();

    L_matrix_.resize(N * N, std::complex<double>(0, 0));
    P_matrix_.resize(N * N, std::complex<double>(0, 0));
    R_matrix_.resize(N * N, std::complex<double>(0, 0));
    Z_matrix_.resize(N * N, std::complex<double>(0, 0));

    // Map from edge index to interior edge index
    std::vector<int> edgeToUnknown(mesh_->NumEdges(), -1);
    int unknownIdx = 0;
    for (int e = 0; e < mesh_->NumEdges(); ++e) {
        if (!mesh_->Edge(e).isBoundary) {
            edgeToUnknown[e] = unknownIdx++;
        }
    }

    // Skin depth
    double skinDepth = 1e-6;  // Default
    if (frequency_ > 0 && conductivity_ > 0) {
        skinDepth = std::sqrt(2.0 / (omega_ * RWG_MU_0 * conductivity_));
    }

    // Surface resistance: R_s = 1/(sigma * delta)
    double R_s = 1.0 / (conductivity_ * skinDepth);

    // Assemble matrices
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < N; ++i) {
        // Find original edge index for unknown i
        int edgeI = -1;
        for (int e = 0; e < mesh_->NumEdges(); ++e) {
            if (edgeToUnknown[e] == i) {
                edgeI = e;
                break;
            }
        }
        if (edgeI < 0) continue;

        for (int j = 0; j < N; ++j) {
            int edgeJ = -1;
            for (int e = 0; e < mesh_->NumEdges(); ++e) {
                if (edgeToUnknown[e] == j) {
                    edgeJ = e;
                    break;
                }
            }
            if (edgeJ < 0) continue;

            // L matrix (vector potential)
            L_matrix_[i * N + j] = RWGIntegration::VectorPotentialIntegral(*mesh_, edgeI, edgeJ);

            // P matrix (scalar potential)
            P_matrix_[i * N + j] = RWGIntegration::ScalarPotentialIntegral(*mesh_, edgeI, edgeJ);

            // R matrix (resistance) - only diagonal
            if (i == j) {
                // Resistance: R = R_s * l / w where w ≈ 2*A/l
                const auto& edge = mesh_->Edge(edgeI);
                double l = edge.length;

                // Average width from adjacent triangles
                double A_avg = 0;
                int nTri = 0;
                if (edge.triPlus >= 0) {
                    A_avg += mesh_->Triangle(edge.triPlus).area;
                    nTri++;
                }
                if (edge.triMinus >= 0) {
                    A_avg += mesh_->Triangle(edge.triMinus).area;
                    nTri++;
                }
                if (nTri > 0) A_avg /= nTri;

                double w_eff = 2.0 * A_avg / l;
                R_matrix_[i * N + j] = std::complex<double>(R_s * l / w_eff, 0);
            }
        }
    }

    // Combined matrix for MQS (Magneto-Quasi-Static):
    // Z = jωL + R
    // Note: In MQS approximation, the scalar potential term P/(jωε) is negligible
    // because displacement current is ignored. Including it causes large capacitive
    // reactance that dominates the inductive behavior.
    //
    // For full-wave EFIE, use: Z = jωL + R + P/(jωε₀)

    for (int i = 0; i < N * N; ++i) {
        Z_matrix_[i] = jOmega_ * L_matrix_[i] + R_matrix_[i];
    }
}

void RWGEFIESolver::AssembleRHS() {
    int N = mesh_->NumInteriorEdges();
    rhs_.resize(N, std::complex<double>(0, 0));

    // Map from edge index to unknown index
    std::vector<int> edgeToUnknown(mesh_->NumEdges(), -1);
    int unknownIdx = 0;
    for (int e = 0; e < mesh_->NumEdges(); ++e) {
        if (!mesh_->Edge(e).isBoundary) {
            edgeToUnknown[e] = unknownIdx++;
        }
    }

    // For voltage excitation at port:
    // V_m = <f_m, E_inc> = integral of E_inc · f_m over surface
    //
    // For a gap voltage V at port edges:
    // E_inc = V / d in the gap direction
    // where d is gap width

    if (portEdges_.empty()) {
        // No port defined - use uniform voltage across first edge
        if (N > 0) {
            rhs_[0] = voltageExcitation_;
        }
    } else {
        // Distribute voltage among port edges
        double totalLength = 0;
        for (int portEdge : portEdges_) {
            if (edgeToUnknown[portEdge] >= 0) {
                totalLength += mesh_->Edge(portEdge).length;
            }
        }

        if (totalLength > 0) {
            for (int portEdge : portEdges_) {
                int unk = edgeToUnknown[portEdge];
                if (unk >= 0) {
                    // Weight by edge length
                    double w = mesh_->Edge(portEdge).length / totalLength;
                    rhs_[unk] = voltageExcitation_ * w;
                }
            }
        }
    }
}

void RWGEFIESolver::SolveLoopStar() {
    // Loop-Star decomposition for low-frequency stability
    //
    // The RWG basis is transformed into Loop (solenoidal) and Star (irrotational) functions:
    //   I_RWG = T_loop * I_L + T_star * I_S
    //
    // The system is transformed as:
    //   [Z_LL  Z_LS] [I_L]   [V_L]
    //   [Z_SL  Z_SS] [I_S] = [V_S]
    //
    // For MQS (Magneto-Quasi-Static) mode, only Loop equations are solved:
    //   Z_LL * I_L = V_L
    //   I_S = 0
    //
    // This is valid when displacement current is negligible (div J = 0).

    int N = static_cast<int>(rhs_.size());
    if (N == 0) return;

    // Build Loop-Star decomposition
    if (!loopStar_) {
        loopStar_ = std::make_unique<LoopStarDecomposition>();
    }
    loopStar_->Build(*mesh_);

    if (!loopStar_->IsValid()) {
        throw std::runtime_error("Loop-Star decomposition failed: invalid mesh topology");
    }

    int numLoops = loopStar_->NumLoops();
    int numStars = loopStar_->NumStars();

    if (numLoops == 0) {
        throw std::runtime_error("Loop-Star decomposition failed: no loops found in mesh");
    }

    // Transform matrix and RHS
    std::vector<std::complex<double>> Z_LL, Z_LS, Z_SL, Z_SS;
    loopStar_->TransformMatrix(Z_matrix_, Z_LL, Z_LS, Z_SL, Z_SS);

    std::vector<std::complex<double>> V_L, V_S;
    loopStar_->TransformRHS(rhs_, V_L, V_S);

    std::vector<std::complex<double>> I_L, I_S;

    if (useMQS_) {
        // MQS mode: solve only Loop equations
        // Z_LL * I_L = V_L
        I_L = V_L;
        I_S.assign(numStars, std::complex<double>(0, 0));

        #ifdef HAVE_LAPACK
        std::vector<std::complex<double>> A = Z_LL;
        std::vector<int> ipiv(numLoops);

        lapack_int info = LAPACKE_zgesv(LAPACK_COL_MAJOR, numLoops, 1,
                                         reinterpret_cast<lapack_complex_double*>(A.data()),
                                         numLoops,
                                         ipiv.data(),
                                         reinterpret_cast<lapack_complex_double*>(I_L.data()),
                                         numLoops);
        if (info != 0) {
            throw std::runtime_error("LAPACK zgesv failed in Loop-Star solve with info = " + std::to_string(info));
        }
        #else
        // Gaussian elimination for Loop equations
        std::vector<std::complex<double>> A = Z_LL;

        for (int k = 0; k < numLoops; ++k) {
            // Find pivot
            int maxRow = k;
            double maxVal = std::abs(A[k * numLoops + k]);
            for (int i = k + 1; i < numLoops; ++i) {
                double val = std::abs(A[i * numLoops + k]);
                if (val > maxVal) {
                    maxVal = val;
                    maxRow = i;
                }
            }

            // Swap rows
            if (maxRow != k) {
                for (int j = 0; j < numLoops; ++j) {
                    std::swap(A[k * numLoops + j], A[maxRow * numLoops + j]);
                }
                std::swap(I_L[k], I_L[maxRow]);
            }

            if (std::abs(A[k * numLoops + k]) < 1e-30) {
                throw std::runtime_error("Singular matrix in Loop-Star MQS solve");
            }

            // Eliminate
            for (int i = k + 1; i < numLoops; ++i) {
                std::complex<double> factor = A[i * numLoops + k] / A[k * numLoops + k];
                for (int j = k; j < numLoops; ++j) {
                    A[i * numLoops + j] -= factor * A[k * numLoops + j];
                }
                I_L[i] -= factor * I_L[k];
            }
        }

        // Back substitution
        for (int i = numLoops - 1; i >= 0; --i) {
            for (int j = i + 1; j < numLoops; ++j) {
                I_L[i] -= A[i * numLoops + j] * I_L[j];
            }
            I_L[i] /= A[i * numLoops + i];
        }
        #endif
    } else {
        // Full EFIE: solve coupled Loop-Star system
        // [Z_LL  Z_LS] [I_L]   [V_L]
        // [Z_SL  Z_SS] [I_S] = [V_S]

        int totalDOF = numLoops + numStars;

        // Assemble full matrix
        std::vector<std::complex<double>> Z_full(totalDOF * totalDOF, std::complex<double>(0, 0));
        std::vector<std::complex<double>> V_full(totalDOF);

        // Z_LL block
        for (int i = 0; i < numLoops; ++i) {
            for (int j = 0; j < numLoops; ++j) {
                Z_full[i * totalDOF + j] = Z_LL[i * numLoops + j];
            }
        }

        // Z_LS block
        for (int i = 0; i < numLoops; ++i) {
            for (int j = 0; j < numStars; ++j) {
                Z_full[i * totalDOF + (numLoops + j)] = Z_LS[i * numStars + j];
            }
        }

        // Z_SL block
        for (int i = 0; i < numStars; ++i) {
            for (int j = 0; j < numLoops; ++j) {
                Z_full[(numLoops + i) * totalDOF + j] = Z_SL[i * numLoops + j];
            }
        }

        // Z_SS block
        for (int i = 0; i < numStars; ++i) {
            for (int j = 0; j < numStars; ++j) {
                Z_full[(numLoops + i) * totalDOF + (numLoops + j)] = Z_SS[i * numStars + j];
            }
        }

        // RHS
        for (int i = 0; i < numLoops; ++i) {
            V_full[i] = V_L[i];
        }
        for (int i = 0; i < numStars; ++i) {
            V_full[numLoops + i] = V_S[i];
        }

        // Solve
        #ifdef HAVE_LAPACK
        std::vector<int> ipiv(totalDOF);
        lapack_int info = LAPACKE_zgesv(LAPACK_COL_MAJOR, totalDOF, 1,
                                         reinterpret_cast<lapack_complex_double*>(Z_full.data()),
                                         totalDOF,
                                         ipiv.data(),
                                         reinterpret_cast<lapack_complex_double*>(V_full.data()),
                                         totalDOF);
        if (info != 0) {
            throw std::runtime_error("LAPACK zgesv failed in full Loop-Star solve with info = " + std::to_string(info));
        }
        #else
        // Gaussian elimination
        for (int k = 0; k < totalDOF; ++k) {
            int maxRow = k;
            double maxVal = std::abs(Z_full[k * totalDOF + k]);
            for (int i = k + 1; i < totalDOF; ++i) {
                double val = std::abs(Z_full[i * totalDOF + k]);
                if (val > maxVal) {
                    maxVal = val;
                    maxRow = i;
                }
            }

            if (maxRow != k) {
                for (int j = 0; j < totalDOF; ++j) {
                    std::swap(Z_full[k * totalDOF + j], Z_full[maxRow * totalDOF + j]);
                }
                std::swap(V_full[k], V_full[maxRow]);
            }

            if (std::abs(Z_full[k * totalDOF + k]) < 1e-30) {
                throw std::runtime_error("Singular matrix in full Loop-Star solve");
            }

            for (int i = k + 1; i < totalDOF; ++i) {
                std::complex<double> factor = Z_full[i * totalDOF + k] / Z_full[k * totalDOF + k];
                for (int j = k; j < totalDOF; ++j) {
                    Z_full[i * totalDOF + j] -= factor * Z_full[k * totalDOF + j];
                }
                V_full[i] -= factor * V_full[k];
            }
        }

        for (int i = totalDOF - 1; i >= 0; --i) {
            for (int j = i + 1; j < totalDOF; ++j) {
                V_full[i] -= Z_full[i * totalDOF + j] * V_full[j];
            }
            V_full[i] /= Z_full[i * totalDOF + i];
        }
        #endif

        // Extract I_L and I_S
        I_L.resize(numLoops);
        I_S.resize(numStars);
        for (int i = 0; i < numLoops; ++i) {
            I_L[i] = V_full[i];
        }
        for (int i = 0; i < numStars; ++i) {
            I_S[i] = V_full[numLoops + i];
        }
    }

    // Transform back to RWG basis
    loopStar_->TransformSolution(I_L, I_S, solution_);
}

void RWGEFIESolver::ComputeImpedance() {
    // Z = V / I
    // I = sum of currents at port edges

    if (portEdges_.empty() && !solution_.empty()) {
        // Single edge port (first unknown)
        totalCurrent_ = solution_[0];
        if (std::abs(totalCurrent_) > 1e-30) {
            impedance_ = voltageExcitation_ / totalCurrent_;
        }
        return;
    }

    // Map from edge index to unknown index
    std::vector<int> edgeToUnknown(mesh_->NumEdges(), -1);
    int unknownIdx = 0;
    for (int e = 0; e < mesh_->NumEdges(); ++e) {
        if (!mesh_->Edge(e).isBoundary) {
            edgeToUnknown[e] = unknownIdx++;
        }
    }

    // Sum currents at port
    totalCurrent_ = std::complex<double>(0, 0);
    for (int portEdge : portEdges_) {
        int unk = edgeToUnknown[portEdge];
        if (unk >= 0) {
            // Current = coefficient * edge length
            totalCurrent_ += solution_[unk] * mesh_->Edge(portEdge).length;
        }
    }

    if (std::abs(totalCurrent_) > 1e-30) {
        impedance_ = voltageExcitation_ / totalCurrent_;
    }
}

void RWGEFIESolver::ComputeB(const TVector3d& point,
                              std::complex<double>& Bx,
                              std::complex<double>& By,
                              std::complex<double>& Bz) const {
    Bx = By = Bz = std::complex<double>(0, 0);

    if (!mesh_ || solution_.empty()) return;

    // Map from edge index to unknown index
    std::vector<int> edgeToUnknown(mesh_->NumEdges(), -1);
    int unknownIdx = 0;
    for (int e = 0; e < mesh_->NumEdges(); ++e) {
        if (!mesh_->Edge(e).isBoundary) {
            edgeToUnknown[e] = unknownIdx++;
        }
    }

    // B = curl(A) = mu_0/(4pi) * integral { J x R / R^3 } dV
    // For surface current K:
    // B = mu_0/(4pi) * integral_S { K x (r - r') / |r - r'|^3 } dA'

    for (int e = 0; e < mesh_->NumEdges(); ++e) {
        if (mesh_->Edge(e).isBoundary) continue;

        int unk = edgeToUnknown[e];
        if (unk < 0) continue;

        std::complex<double> I_coeff = solution_[unk];
        const auto& edge = mesh_->Edge(e);

        // Integrate over both triangles of this edge
        int tris[2] = {edge.triPlus, edge.triMinus};

        for (int t = 0; t < 2; ++t) {
            if (tris[t] < 0) continue;

            const auto& tri = mesh_->Triangle(tris[t]);
            const TVector3d& v0 = mesh_->Vertex(tri.vertices[0]);
            const TVector3d& v1 = mesh_->Vertex(tri.vertices[1]);
            const TVector3d& v2 = mesh_->Vertex(tri.vertices[2]);

            // Numerical quadrature
            for (int q = 0; q < RWGIntegration::NQUAD_TRI; ++q) {
                double L1 = RWGIntegration::triQuadPts_[q][0];
                double L2 = RWGIntegration::triQuadPts_[q][1];
                double L3 = 1.0 - L1 - L2;

                TVector3d r_src = L1 * v0 + L2 * v1 + L3 * v2;

                // Distance vector
                TVector3d R_vec = point - r_src;
                double R = R_vec.Abs();
                if (R < 1e-15) continue;

                double R3 = R * R * R;

                // RWG basis at source point (this is the current direction)
                TVector3d f = mesh_->EvaluateRWG(e, tris[t], r_src);

                // Current: K = I_coeff * f
                // B contribution: (mu_0/4pi) * K x R / R^3 * dA
                TVector3d K(I_coeff.real() * f.x, I_coeff.real() * f.y, I_coeff.real() * f.z);

                // Cross product K x R
                double crossX = K.y * R_vec.z - K.z * R_vec.y;
                double crossY = K.z * R_vec.x - K.x * R_vec.z;
                double crossZ = K.x * R_vec.y - K.y * R_vec.x;

                double coeff = RWG_MU_0 * RWG_INV_FOUR_PI * RWGIntegration::triQuadWts_[q] * tri.area / R3;

                Bx += std::complex<double>(coeff * crossX, 0);
                By += std::complex<double>(coeff * crossY, 0);
                Bz += std::complex<double>(coeff * crossZ, 0);
            }
        }
    }
}

} // namespace radia
