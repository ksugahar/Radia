/*
 * rad_conductor.cpp
 *
 * Conductor element implementation
 *
 * Part of Radia project
 */

#include "rad_conductor.h"
#include "rad_constants.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

#ifdef HAVE_LAPACK
#include <mkl_lapacke.h>
#endif

namespace radia {

// ============================================================================
// radTConductor implementation
// ============================================================================

radTConductor::radTConductor()
    : conductivity_(5.8e7)  // Default: copper
    , frequency_(0)
    , formulation_(ConductorFormulation::DC)
    , portImpedance_(0, 0)
    , excitationType_(0)
    , excitationValue_(0, 0)
    , totalCurrent_(0, 0)
{
}

radTConductor::~radTConductor() {
}

void radTConductor::SetConductivity(double conductivity) {
    if (conductivity <= 0) {
        throw std::invalid_argument("Conductivity must be positive");
    }
    conductivity_ = conductivity;
}

void radTConductor::SetFrequency(double frequency) {
    frequency_ = frequency;

    // Auto-select formulation based on frequency
    if (frequency <= 0) {
        formulation_ = ConductorFormulation::DC;
    } else {
        // Keep current formulation if already set for AC
        if (formulation_ == ConductorFormulation::DC) {
            formulation_ = ConductorFormulation::MQS;
        }
    }
}

void radTConductor::CreateFromRecBlock(const TVector3d& center,
                                        const TVector3d& dimensions,
                                        int numPanelsPerFace) {
    panels_.clear();

    double hx = dimensions.x / 2.0;
    double hy = dimensions.y / 2.0;
    double hz = dimensions.z / 2.0;

    // 8 corners of the block
    TVector3d corners[8] = {
        TVector3d(center.x - hx, center.y - hy, center.z - hz),
        TVector3d(center.x + hx, center.y - hy, center.z - hz),
        TVector3d(center.x + hx, center.y + hy, center.z - hz),
        TVector3d(center.x - hx, center.y + hy, center.z - hz),
        TVector3d(center.x - hx, center.y - hy, center.z + hz),
        TVector3d(center.x + hx, center.y - hy, center.z + hz),
        TVector3d(center.x + hx, center.y + hy, center.z + hz),
        TVector3d(center.x - hx, center.y + hy, center.z + hz)
    };

    int n = numPanelsPerFace;

    // Face -Z (bottom)
    GenerateRectangularFacePanels(corners[0],
                                  TVector3d(dimensions.x, 0, 0),
                                  TVector3d(0, dimensions.y, 0),
                                  n, n);

    // Face +Z (top)
    GenerateRectangularFacePanels(corners[4],
                                  TVector3d(dimensions.x, 0, 0),
                                  TVector3d(0, dimensions.y, 0),
                                  n, n);

    // Face -Y
    GenerateRectangularFacePanels(corners[0],
                                  TVector3d(dimensions.x, 0, 0),
                                  TVector3d(0, 0, dimensions.z),
                                  n, n);

    // Face +Y
    GenerateRectangularFacePanels(corners[3],
                                  TVector3d(dimensions.x, 0, 0),
                                  TVector3d(0, 0, dimensions.z),
                                  n, n);

    // Face -X
    GenerateRectangularFacePanels(corners[0],
                                  TVector3d(0, dimensions.y, 0),
                                  TVector3d(0, 0, dimensions.z),
                                  n, n);

    // Face +X
    GenerateRectangularFacePanels(corners[1],
                                  TVector3d(0, dimensions.y, 0),
                                  TVector3d(0, 0, dimensions.z),
                                  n, n);

    ComputePanelNormals();

    // Initialize solution vectors
    int nPanels = static_cast<int>(panels_.size());
    surfaceCurrent_.resize(3 * nPanels, std::complex<double>(0, 0));
    surfaceCharge_.resize(nPanels, std::complex<double>(0, 0));
}

void radTConductor::CreateFromHexahedron(const std::vector<TVector3d>& vertices,
                                          int numPanelsPerFace) {
    if (vertices.size() != 8) {
        throw std::invalid_argument("Hexahedron requires exactly 8 vertices");
    }

    panels_.clear();

    int n = numPanelsPerFace;

    // Hexahedron face vertex indices (same as Radia convention)
    // Face 0: bottom (0,1,2,3)
    // Face 1: top (4,5,6,7)
    // Face 2: front (0,1,5,4)
    // Face 3: back (2,3,7,6)
    // Face 4: left (0,3,7,4)
    // Face 5: right (1,2,6,5)

    int faceIndices[6][4] = {
        {0, 1, 2, 3},  // bottom
        {4, 5, 6, 7},  // top
        {0, 1, 5, 4},  // front
        {2, 3, 7, 6},  // back
        {0, 3, 7, 4},  // left
        {1, 2, 6, 5}   // right
    };

    for (int f = 0; f < 6; ++f) {
        const TVector3d& v0 = vertices[faceIndices[f][0]];
        const TVector3d& v1 = vertices[faceIndices[f][1]];
        const TVector3d& v2 = vertices[faceIndices[f][2]];
        const TVector3d& v3 = vertices[faceIndices[f][3]];

        // Edge vectors
        TVector3d edge1(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z);
        TVector3d edge2(v3.x - v0.x, v3.y - v0.y, v3.z - v0.z);

        GenerateRectangularFacePanels(v0, edge1, edge2, n, n);
    }

    ComputePanelNormals();

    int nPanels = static_cast<int>(panels_.size());
    surfaceCurrent_.resize(3 * nPanels, std::complex<double>(0, 0));
    surfaceCharge_.resize(nPanels, std::complex<double>(0, 0));
}

void radTConductor::CreateFromRadiaObject(radTg3d* radiaObj, int numPanelsPerFace) {
    // TODO: Extract geometry from Radia object
    // This will interface with existing Radia geometry classes
    throw std::runtime_error("CreateFromRadiaObject not yet implemented");
}

void radTConductor::CreateFromSurfaceMesh(const std::vector<TVector3d>& vertices,
                                           const std::vector<std::array<int, 3>>& triangles) {
    panels_.clear();
    panels_.reserve(triangles.size());

    for (const auto& tri : triangles) {
        SurfacePanel panel;
        panel.type = SurfacePanel::Triangle;

        panel.vertices.resize(3);
        panel.vertices[0] = vertices[tri[0]];
        panel.vertices[1] = vertices[tri[1]];
        panel.vertices[2] = vertices[tri[2]];

        // Compute center
        panel.center.x = (panel.vertices[0].x + panel.vertices[1].x + panel.vertices[2].x) / 3.0;
        panel.center.y = (panel.vertices[0].y + panel.vertices[1].y + panel.vertices[2].y) / 3.0;
        panel.center.z = (panel.vertices[0].z + panel.vertices[1].z + panel.vertices[2].z) / 3.0;

        // Compute normal and area
        TVector3d e1(panel.vertices[1].x - panel.vertices[0].x,
                     panel.vertices[1].y - panel.vertices[0].y,
                     panel.vertices[1].z - panel.vertices[0].z);
        TVector3d e2(panel.vertices[2].x - panel.vertices[0].x,
                     panel.vertices[2].y - panel.vertices[0].y,
                     panel.vertices[2].z - panel.vertices[0].z);

        // Cross product
        panel.normal.x = e1.y * e2.z - e1.z * e2.y;
        panel.normal.y = e1.z * e2.x - e1.x * e2.z;
        panel.normal.z = e1.x * e2.y - e1.y * e2.x;

        double len = std::sqrt(panel.normal.x * panel.normal.x +
                               panel.normal.y * panel.normal.y +
                               panel.normal.z * panel.normal.z);

        panel.area = len / 2.0;

        if (len > 1e-15) {
            panel.normal.x /= len;
            panel.normal.y /= len;
            panel.normal.z /= len;
        }

        panels_.push_back(panel);
    }

    int nPanels = static_cast<int>(panels_.size());
    surfaceCurrent_.resize(3 * nPanels, std::complex<double>(0, 0));
    surfaceCharge_.resize(nPanels, std::complex<double>(0, 0));
}

void radTConductor::CreateFromPanels(const std::vector<SurfacePanel>& panels) {
    panels_ = panels;

    int nPanels = static_cast<int>(panels_.size());
    surfaceCurrent_.resize(3 * nPanels, std::complex<double>(0, 0));
    surfaceCharge_.resize(nPanels, std::complex<double>(0, 0));
}

void radTConductor::CreateWire(const std::vector<TVector3d>& path,
                                const std::string& crossSection,
                                double width,
                                double height,
                                int numPanelsAround,
                                int numPanelsAlongLength) {
    panels_.clear();

    if (path.size() < 2) {
        throw std::invalid_argument("Wire path must have at least 2 points");
    }

    bool isCircular = (crossSection == "circular");
    if (height <= 0 && !isCircular) {
        height = width;  // Square cross-section
    }

    // Auto-compute panels along length
    if (numPanelsAlongLength <= 0) {
        double totalLength = 0;
        for (size_t i = 1; i < path.size(); ++i) {
            double dx = path[i].x - path[i-1].x;
            double dy = path[i].y - path[i-1].y;
            double dz = path[i].z - path[i-1].z;
            totalLength += std::sqrt(dx*dx + dy*dy + dz*dz);
        }
        double avgPanelSize = (isCircular ? RadConst::PI * width : 2*(width + height)) / numPanelsAround;
        numPanelsAlongLength = std::max(1, static_cast<int>(totalLength / avgPanelSize));
    }

    // Generate panels along wire
    for (size_t seg = 0; seg < path.size() - 1; ++seg) {
        const TVector3d& p0 = path[seg];
        const TVector3d& p1 = path[seg + 1];

        // Tangent direction
        TVector3d tangent;
        tangent.x = p1.x - p0.x;
        tangent.y = p1.y - p0.y;
        tangent.z = p1.z - p0.z;

        double segLen = std::sqrt(tangent.x*tangent.x + tangent.y*tangent.y + tangent.z*tangent.z);
        if (segLen < 1e-15) continue;

        tangent.x /= segLen;
        tangent.y /= segLen;
        tangent.z /= segLen;

        // Find perpendicular directions
        TVector3d perp1, perp2;

        // Choose initial perpendicular
        if (std::abs(tangent.z) < 0.9) {
            perp1.x = -tangent.y;
            perp1.y = tangent.x;
            perp1.z = 0;
        } else {
            perp1.x = 0;
            perp1.y = -tangent.z;
            perp1.z = tangent.y;
        }

        double len1 = std::sqrt(perp1.x*perp1.x + perp1.y*perp1.y + perp1.z*perp1.z);
        perp1.x /= len1;
        perp1.y /= len1;
        perp1.z /= len1;

        // Cross product for second perpendicular
        perp2.x = tangent.y * perp1.z - tangent.z * perp1.y;
        perp2.y = tangent.z * perp1.x - tangent.x * perp1.z;
        perp2.z = tangent.x * perp1.y - tangent.y * perp1.x;

        // Generate panels around circumference
        int nAlong = numPanelsAlongLength / static_cast<int>(path.size() - 1);
        nAlong = std::max(1, nAlong);

        double dl = segLen / nAlong;

        for (int iAlong = 0; iAlong < nAlong; ++iAlong) {
            double t0 = static_cast<double>(iAlong) / nAlong;
            double t1 = static_cast<double>(iAlong + 1) / nAlong;

            TVector3d c0, c1;
            c0.x = p0.x + t0 * (p1.x - p0.x);
            c0.y = p0.y + t0 * (p1.y - p0.y);
            c0.z = p0.z + t0 * (p1.z - p0.z);
            c1.x = p0.x + t1 * (p1.x - p0.x);
            c1.y = p0.y + t1 * (p1.y - p0.y);
            c1.z = p0.z + t1 * (p1.z - p0.z);

            if (isCircular) {
                // Circular cross-section
                double radius = width / 2.0;
                for (int iAround = 0; iAround < numPanelsAround; ++iAround) {
                    double theta0 = 2.0 * RadConst::PI * iAround / numPanelsAround;
                    double theta1 = 2.0 * RadConst::PI * (iAround + 1) / numPanelsAround;

                    SurfacePanel panel;
                    panel.type = SurfacePanel::Quadrilateral;
                    panel.vertices.resize(4);

                    double cos0 = std::cos(theta0), sin0 = std::sin(theta0);
                    double cos1 = std::cos(theta1), sin1 = std::sin(theta1);

                    // 4 corners of panel
                    for (int corner = 0; corner < 4; ++corner) {
                        double c_theta = (corner < 2) ? cos0 : cos1;
                        double s_theta = (corner < 2) ? sin0 : sin1;
                        const TVector3d& c_along = (corner == 0 || corner == 3) ? c0 : c1;

                        panel.vertices[corner].x = c_along.x + radius * (c_theta * perp1.x + s_theta * perp2.x);
                        panel.vertices[corner].y = c_along.y + radius * (c_theta * perp1.y + s_theta * perp2.y);
                        panel.vertices[corner].z = c_along.z + radius * (c_theta * perp1.z + s_theta * perp2.z);
                    }

                    // Compute center, normal, area
                    panel.center.x = 0.25 * (panel.vertices[0].x + panel.vertices[1].x +
                                             panel.vertices[2].x + panel.vertices[3].x);
                    panel.center.y = 0.25 * (panel.vertices[0].y + panel.vertices[1].y +
                                             panel.vertices[2].y + panel.vertices[3].y);
                    panel.center.z = 0.25 * (panel.vertices[0].z + panel.vertices[1].z +
                                             panel.vertices[2].z + panel.vertices[3].z);

                    double thetaC = (theta0 + theta1) / 2.0;
                    panel.normal.x = std::cos(thetaC) * perp1.x + std::sin(thetaC) * perp2.x;
                    panel.normal.y = std::cos(thetaC) * perp1.y + std::sin(thetaC) * perp2.y;
                    panel.normal.z = std::cos(thetaC) * perp1.z + std::sin(thetaC) * perp2.z;

                    double dtheta = (theta1 - theta0);
                    panel.area = radius * dtheta * dl;

                    panels_.push_back(panel);
                }
            } else {
                // Rectangular cross-section: 4 faces (top, bottom, left, right)
                double hw = width / 2.0;   // half width
                double hh = height / 2.0;  // half height

                // Face offsets in local (perp1, perp2) coordinates
                // Face 0: +perp1 direction (right face)
                // Face 1: -perp1 direction (left face)
                // Face 2: +perp2 direction (top face)
                // Face 3: -perp2 direction (bottom face)

                struct FaceInfo {
                    double offset1, offset2;  // offset in perp1, perp2
                    double dir1, dir2;        // normal direction in perp1, perp2
                    double faceWidth;         // face width
                };

                FaceInfo faces[4] = {
                    { hw, 0, 1, 0, height},   // +perp1 (right)
                    {-hw, 0,-1, 0, height},   // -perp1 (left)
                    { 0, hh, 0, 1, width},    // +perp2 (top)
                    { 0,-hh, 0,-1, width}     // -perp2 (bottom)
                };

                for (int face = 0; face < 4; ++face) {
                    const FaceInfo& fi = faces[face];

                    SurfacePanel panel;
                    panel.type = SurfacePanel::Quadrilateral;
                    panel.vertices.resize(4);

                    // Face center offset
                    double offsetX = fi.offset1 * perp1.x + fi.offset2 * perp2.x;
                    double offsetY = fi.offset1 * perp1.y + fi.offset2 * perp2.y;
                    double offsetZ = fi.offset1 * perp1.z + fi.offset2 * perp2.z;

                    // Width direction (perpendicular to normal and tangent)
                    TVector3d widthDir;
                    if (face < 2) {
                        // Left/right faces: width along perp2
                        widthDir = perp2;
                    } else {
                        // Top/bottom faces: width along perp1
                        widthDir = perp1;
                    }

                    double halfFaceWidth = fi.faceWidth / 2.0;

                    // 4 corners: (c0, -halfWidth), (c1, -halfWidth), (c1, +halfWidth), (c0, +halfWidth)
                    panel.vertices[0].x = c0.x + offsetX - halfFaceWidth * widthDir.x;
                    panel.vertices[0].y = c0.y + offsetY - halfFaceWidth * widthDir.y;
                    panel.vertices[0].z = c0.z + offsetZ - halfFaceWidth * widthDir.z;

                    panel.vertices[1].x = c1.x + offsetX - halfFaceWidth * widthDir.x;
                    panel.vertices[1].y = c1.y + offsetY - halfFaceWidth * widthDir.y;
                    panel.vertices[1].z = c1.z + offsetZ - halfFaceWidth * widthDir.z;

                    panel.vertices[2].x = c1.x + offsetX + halfFaceWidth * widthDir.x;
                    panel.vertices[2].y = c1.y + offsetY + halfFaceWidth * widthDir.y;
                    panel.vertices[2].z = c1.z + offsetZ + halfFaceWidth * widthDir.z;

                    panel.vertices[3].x = c0.x + offsetX + halfFaceWidth * widthDir.x;
                    panel.vertices[3].y = c0.y + offsetY + halfFaceWidth * widthDir.y;
                    panel.vertices[3].z = c0.z + offsetZ + halfFaceWidth * widthDir.z;

                    // Center
                    panel.center.x = 0.25 * (panel.vertices[0].x + panel.vertices[1].x +
                                             panel.vertices[2].x + panel.vertices[3].x);
                    panel.center.y = 0.25 * (panel.vertices[0].y + panel.vertices[1].y +
                                             panel.vertices[2].y + panel.vertices[3].y);
                    panel.center.z = 0.25 * (panel.vertices[0].z + panel.vertices[1].z +
                                             panel.vertices[2].z + panel.vertices[3].z);

                    // Normal (outward from wire center)
                    panel.normal.x = fi.dir1 * perp1.x + fi.dir2 * perp2.x;
                    panel.normal.y = fi.dir1 * perp1.y + fi.dir2 * perp2.y;
                    panel.normal.z = fi.dir1 * perp1.z + fi.dir2 * perp2.z;

                    // Area
                    panel.area = fi.faceWidth * dl;

                    panels_.push_back(panel);
                }
            }
        }
    }

    int nPanels = static_cast<int>(panels_.size());
    surfaceCurrent_.resize(3 * nPanels, std::complex<double>(0, 0));
    surfaceCharge_.resize(nPanels, std::complex<double>(0, 0));
}

void radTConductor::CreateLoop(const TVector3d& center,
                                double radius,
                                const TVector3d& normal,
                                const std::string& crossSection,
                                double wireWidth,
                                double wireHeight,
                                int numPanelsAround,
                                int numPanelsLoop) {
    // Generate path along loop
    std::vector<TVector3d> path(numPanelsLoop + 1);

    // Normalize loop normal
    double nLen = std::sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
    TVector3d n;
    n.x = normal.x / nLen;
    n.y = normal.y / nLen;
    n.z = normal.z / nLen;

    // Find two perpendicular directions in loop plane
    TVector3d u, v;
    if (std::abs(n.z) < 0.9) {
        u.x = -n.y;
        u.y = n.x;
        u.z = 0;
    } else {
        u.x = 0;
        u.y = -n.z;
        u.z = n.y;
    }
    double uLen = std::sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
    u.x /= uLen; u.y /= uLen; u.z /= uLen;

    v.x = n.y * u.z - n.z * u.y;
    v.y = n.z * u.x - n.x * u.z;
    v.z = n.x * u.y - n.y * u.x;

    for (int i = 0; i <= numPanelsLoop; ++i) {
        double theta = 2.0 * RadConst::PI * i / numPanelsLoop;
        double cosT = std::cos(theta);
        double sinT = std::sin(theta);

        path[i].x = center.x + radius * (cosT * u.x + sinT * v.x);
        path[i].y = center.y + radius * (cosT * u.y + sinT * v.y);
        path[i].z = center.z + radius * (cosT * u.z + sinT * v.z);
    }

    CreateWire(path, crossSection, wireWidth, wireHeight, numPanelsAround, numPanelsLoop);
}

void radTConductor::CreateSpiral(const TVector3d& center,
                                  double innerRadius,
                                  double outerRadius,
                                  double pitch,
                                  int numTurns,
                                  const TVector3d& axis,
                                  const std::string& crossSection,
                                  double wireWidth,
                                  double wireHeight,
                                  int numPanelsAround) {
    // Generate spiral path
    int pointsPerTurn = 30;
    int totalPoints = numTurns * pointsPerTurn + 1;
    std::vector<TVector3d> path(totalPoints);

    // Normalize axis
    double aLen = std::sqrt(axis.x*axis.x + axis.y*axis.y + axis.z*axis.z);
    TVector3d a;
    a.x = axis.x / aLen;
    a.y = axis.y / aLen;
    a.z = axis.z / aLen;

    // Find perpendicular directions
    TVector3d u, v;
    if (std::abs(a.z) < 0.9) {
        u.x = -a.y;
        u.y = a.x;
        u.z = 0;
    } else {
        u.x = 0;
        u.y = -a.z;
        u.z = a.y;
    }
    double uLen = std::sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
    u.x /= uLen; u.y /= uLen; u.z /= uLen;

    v.x = a.y * u.z - a.z * u.y;
    v.y = a.z * u.x - a.x * u.z;
    v.z = a.x * u.y - a.y * u.x;

    for (int i = 0; i < totalPoints; ++i) {
        double t = static_cast<double>(i) / (totalPoints - 1);
        double theta = 2.0 * RadConst::PI * numTurns * t;
        double r = innerRadius + (outerRadius - innerRadius) * t;
        double z = pitch * numTurns * t;

        double cosT = std::cos(theta);
        double sinT = std::sin(theta);

        path[i].x = center.x + r * (cosT * u.x + sinT * v.x) + z * a.x;
        path[i].y = center.y + r * (cosT * u.y + sinT * v.y) + z * a.y;
        path[i].z = center.z + r * (cosT * u.z + sinT * v.z) + z * a.z;
    }

    CreateWire(path, crossSection, wireWidth, wireHeight, numPanelsAround, totalPoints - 1);
}

std::vector<TVector3d> radTConductor::GetPanelCenters() const {
    std::vector<TVector3d> centers(panels_.size());
    for (size_t i = 0; i < panels_.size(); ++i) {
        centers[i] = panels_[i].center;
    }
    return centers;
}

std::vector<double> radTConductor::GetPanelAreas() const {
    std::vector<double> areas(panels_.size());
    for (size_t i = 0; i < panels_.size(); ++i) {
        areas[i] = panels_[i].area;
    }
    return areas;
}

void radTConductor::GenerateRectangularFacePanels(const TVector3d& corner,
                                                   const TVector3d& edge1,
                                                   const TVector3d& edge2,
                                                   int n1, int n2) {
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            double t1_0 = static_cast<double>(i) / n1;
            double t1_1 = static_cast<double>(i + 1) / n1;
            double t2_0 = static_cast<double>(j) / n2;
            double t2_1 = static_cast<double>(j + 1) / n2;

            SurfacePanel panel;
            panel.type = SurfacePanel::Quadrilateral;
            panel.vertices.resize(4);

            // 4 corners
            panel.vertices[0].x = corner.x + t1_0 * edge1.x + t2_0 * edge2.x;
            panel.vertices[0].y = corner.y + t1_0 * edge1.y + t2_0 * edge2.y;
            panel.vertices[0].z = corner.z + t1_0 * edge1.z + t2_0 * edge2.z;

            panel.vertices[1].x = corner.x + t1_1 * edge1.x + t2_0 * edge2.x;
            panel.vertices[1].y = corner.y + t1_1 * edge1.y + t2_0 * edge2.y;
            panel.vertices[1].z = corner.z + t1_1 * edge1.z + t2_0 * edge2.z;

            panel.vertices[2].x = corner.x + t1_1 * edge1.x + t2_1 * edge2.x;
            panel.vertices[2].y = corner.y + t1_1 * edge1.y + t2_1 * edge2.y;
            panel.vertices[2].z = corner.z + t1_1 * edge1.z + t2_1 * edge2.z;

            panel.vertices[3].x = corner.x + t1_0 * edge1.x + t2_1 * edge2.x;
            panel.vertices[3].y = corner.y + t1_0 * edge1.y + t2_1 * edge2.y;
            panel.vertices[3].z = corner.z + t1_0 * edge1.z + t2_1 * edge2.z;

            // Center
            panel.center.x = 0.25 * (panel.vertices[0].x + panel.vertices[1].x +
                                     panel.vertices[2].x + panel.vertices[3].x);
            panel.center.y = 0.25 * (panel.vertices[0].y + panel.vertices[1].y +
                                     panel.vertices[2].y + panel.vertices[3].y);
            panel.center.z = 0.25 * (panel.vertices[0].z + panel.vertices[1].z +
                                     panel.vertices[2].z + panel.vertices[3].z);

            // Area
            double dx1 = edge1.x / n1;
            double dy1 = edge1.y / n1;
            double dz1 = edge1.z / n1;
            double dx2 = edge2.x / n2;
            double dy2 = edge2.y / n2;
            double dz2 = edge2.z / n2;

            // Cross product magnitude
            double cx = dy1 * dz2 - dz1 * dy2;
            double cy = dz1 * dx2 - dx1 * dz2;
            double cz = dx1 * dy2 - dy1 * dx2;
            panel.area = std::sqrt(cx*cx + cy*cy + cz*cz);

            panels_.push_back(panel);
        }
    }
}

void radTConductor::ComputePanelNormals() {
    for (auto& panel : panels_) {
        if (panel.vertices.size() >= 3) {
            TVector3d e1, e2;
            e1.x = panel.vertices[1].x - panel.vertices[0].x;
            e1.y = panel.vertices[1].y - panel.vertices[0].y;
            e1.z = panel.vertices[1].z - panel.vertices[0].z;

            int idx2 = (panel.vertices.size() > 3) ? 3 : 2;
            e2.x = panel.vertices[idx2].x - panel.vertices[0].x;
            e2.y = panel.vertices[idx2].y - panel.vertices[0].y;
            e2.z = panel.vertices[idx2].z - panel.vertices[0].z;

            panel.normal.x = e1.y * e2.z - e1.z * e2.y;
            panel.normal.y = e1.z * e2.x - e1.x * e2.z;
            panel.normal.z = e1.x * e2.y - e1.y * e2.x;

            double len = std::sqrt(panel.normal.x * panel.normal.x +
                                   panel.normal.y * panel.normal.y +
                                   panel.normal.z * panel.normal.z);
            if (len > 1e-15) {
                panel.normal.x /= len;
                panel.normal.y /= len;
                panel.normal.z /= len;
            }
        }
    }
}

std::complex<double> radTConductor::GreenFunction(double r) const {
    if (r < 1e-15) {
        return std::complex<double>(0, 0);
    }

    switch (formulation_) {
        case ConductorFormulation::DC:
        case ConductorFormulation::MQS:
            // G(r) = 1 / (4*pi*r)
            return std::complex<double>(INV_FOUR_PI / r, 0);

        case ConductorFormulation::EMQS:
        case ConductorFormulation::FullWave: {
            // G(r) = exp(-jkr) / (4*pi*r)
            double omega = 2.0 * RadConst::PI * frequency_;
            double k = omega * std::sqrt(MU_0 * EPS_0);
            std::complex<double> jkr(0, k * r);
            return std::exp(-jkr) * INV_FOUR_PI / r;
        }
    }

    return std::complex<double>(0, 0);
}

void radTConductor::DefinePort(const std::vector<int>& terminal1,
                                const std::vector<int>& terminal2) {
    portTerminal1_ = terminal1;
    portTerminal2_ = terminal2;
}

void radTConductor::SetVoltageExcitation(double V_real, double V_imag) {
    excitationType_ = 1;  // Voltage excitation
    excitationValue_ = std::complex<double>(V_real, V_imag);
}

void radTConductor::SetCurrentExcitation(double I_real, double I_imag) {
    excitationType_ = 2;  // Current excitation
    excitationValue_ = std::complex<double>(I_real, I_imag);
}

void radTConductor::ComputeB(const TVector3d& point,
                              std::complex<double>& Bx,
                              std::complex<double>& By,
                              std::complex<double>& Bz) const {
    Bx = By = Bz = std::complex<double>(0, 0);

    // Biot-Savart or full-wave depending on formulation
    for (size_t i = 0; i < panels_.size(); ++i) {
        const auto& panel = panels_[i];

        double dx = point.x - panel.center.x;
        double dy = point.y - panel.center.y;
        double dz = point.z - panel.center.z;
        double r = std::sqrt(dx*dx + dy*dy + dz*dz);

        if (r < 1e-15) continue;

        // Get surface current at this panel
        std::complex<double> Kx = surfaceCurrent_[3*i];
        std::complex<double> Ky = surfaceCurrent_[3*i + 1];
        std::complex<double> Kz = surfaceCurrent_[3*i + 2];

        // B = mu_0 / (4*pi) * integral{ K x r_hat / r^2 } dA
        // Simplified: assume uniform K over panel
        std::complex<double> G = GreenFunction(r);
        double r3 = r * r * r;

        // K x r contribution
        double rx = dx / r, ry = dy / r, rz = dz / r;

        Bx += MU_0 * G * (Ky * rz - Kz * ry) * panel.area / r;
        By += MU_0 * G * (Kz * rx - Kx * rz) * panel.area / r;
        Bz += MU_0 * G * (Kx * ry - Ky * rx) * panel.area / r;
    }
}

std::complex<double> radTConductor::ComputeB(const TVector3d& point, int component) const {
    std::complex<double> Bx, By, Bz;
    ComputeB(point, Bx, By, Bz);

    switch (component) {
        case 0: return Bx;
        case 1: return By;
        case 2: return Bz;
        default: return std::complex<double>(0, 0);
    }
}

void radTConductor::ComputeE(const TVector3d& point,
                              std::complex<double>& Ex,
                              std::complex<double>& Ey,
                              std::complex<double>& Ez) const {
    Ex = Ey = Ez = std::complex<double>(0, 0);

    // E = -grad(phi) - j*omega*A
    // For now, simplified implementation
    double omega = 2.0 * RadConst::PI * frequency_;

    for (size_t i = 0; i < panels_.size(); ++i) {
        const auto& panel = panels_[i];

        double dx = point.x - panel.center.x;
        double dy = point.y - panel.center.y;
        double dz = point.z - panel.center.z;
        double r = std::sqrt(dx*dx + dy*dy + dz*dz);

        if (r < 1e-15) continue;

        std::complex<double> sigma = surfaceCharge_[i];
        std::complex<double> G = GreenFunction(r);

        // -grad(phi) contribution from charges
        double r3 = r * r * r;
        Ex += sigma * dx / (4.0 * RadConst::PI * EPS_0 * r3) * panel.area;
        Ey += sigma * dy / (4.0 * RadConst::PI * EPS_0 * r3) * panel.area;
        Ez += sigma * dz / (4.0 * RadConst::PI * EPS_0 * r3) * panel.area;
    }
}

void radTConductor::ComputeA(const TVector3d& point,
                              std::complex<double>& Ax,
                              std::complex<double>& Ay,
                              std::complex<double>& Az) const {
    Ax = Ay = Az = std::complex<double>(0, 0);

    // Vector potential: A = mu_0/(4*pi) * integral{ K / |r - r'| } dA'
    // Using Green's function: G(r) = 1/(4*pi*r) for DC/MQS, exp(-jkr)/(4*pi*r) for full-wave
    for (size_t i = 0; i < panels_.size(); ++i) {
        const auto& panel = panels_[i];

        double dx = point.x - panel.center.x;
        double dy = point.y - panel.center.y;
        double dz = point.z - panel.center.z;
        double r = std::sqrt(dx*dx + dy*dy + dz*dz);

        if (r < 1e-15) continue;

        // Get surface current at this panel (3 components)
        std::complex<double> Kx = surfaceCurrent_[3*i];
        std::complex<double> Ky = surfaceCurrent_[3*i + 1];
        std::complex<double> Kz = surfaceCurrent_[3*i + 2];

        // A = mu_0 * G(r) * K * dA
        std::complex<double> G = GreenFunction(r);
        std::complex<double> coeff = MU_0 * G * panel.area;

        Ax += coeff * Kx;
        Ay += coeff * Ky;
        Az += coeff * Kz;
    }
}

std::complex<double> radTConductor::ComputePhi(const TVector3d& point) const {
    std::complex<double> phi(0, 0);

    // Scalar potential: Phi = 1/(4*pi*eps_0) * integral{ sigma / |r - r'| } dA'
    for (size_t i = 0; i < panels_.size(); ++i) {
        const auto& panel = panels_[i];

        double dx = point.x - panel.center.x;
        double dy = point.y - panel.center.y;
        double dz = point.z - panel.center.z;
        double r = std::sqrt(dx*dx + dy*dy + dz*dz);

        if (r < 1e-15) continue;

        std::complex<double> sigma = surfaceCharge_[i];
        std::complex<double> G = GreenFunction(r);

        // Phi = sigma / eps_0 * G(r) * dA
        phi += sigma / EPS_0 * G * panel.area;
    }

    return phi;
}

// ============================================================================
// radTConductorSolver implementation
// ============================================================================

radTConductorSolver::radTConductorSolver()
    : frequency_(0)
    , formulation_(ConductorFormulation::MQS)
    , usePfft_(true)
    , useLoopStar_(false)  // Disabled by default for now
    , directImpedanceMode_(false)
    , numLoops_(0)
    , numStars_(0)
{
}

radTConductorSolver::~radTConductorSolver() {
}

void radTConductorSolver::AddConductor(std::shared_ptr<radTConductor> conductor) {
    conductors_.push_back(conductor);
}

void radTConductorSolver::Clear() {
    conductors_.clear();
    pfft_.reset();
}

void radTConductorSolver::SetFrequency(double frequency) {
    frequency_ = frequency;
    for (auto& cond : conductors_) {
        cond->SetFrequency(frequency);
    }
}

void radTConductorSolver::SetFormulation(ConductorFormulation form) {
    formulation_ = form;
    for (auto& cond : conductors_) {
        cond->SetFormulation(form);
    }
}

void radTConductorSolver::Solve() {
    BuildSystemMatrix();
    SolveLinearSystem();
    ExtractSolution();
}

std::complex<double> radTConductorSolver::SolveImpedance(double frequency) {
    SetFrequency(frequency);
    Solve();

    // Sum port impedances from all conductors
    std::complex<double> Z(0, 0);
    for (const auto& cond : conductors_) {
        Z += cond->GetPortImpedance();
    }
    return Z;
}

std::vector<std::complex<double>> radTConductorSolver::ImpedanceSweep(
    const std::vector<double>& frequencies) {

    std::vector<std::complex<double>> result(frequencies.size());

    for (size_t i = 0; i < frequencies.size(); ++i) {
        result[i] = SolveImpedance(frequencies[i]);
    }

    return result;
}

void radTConductorSolver::BuildSystemMatrix() {
    // ==========================================================================
    // Segment-Based Inductance Calculation (Darwin Approximation)
    // ==========================================================================
    //
    // For wire-type conductors (loops, spirals), we use segment-based unknowns:
    //
    // 1. Divide wire into N_seg segments (N_seg = nPanels / nPanelsAround)
    // 2. Unknown: I_seg (segment current), NOT K_panel (surface current density)
    // 3. Matrix size: N_seg x N_seg (much smaller than N_panels x N_panels)
    //
    // Self inductance of segment (partial self-inductance):
    //   L_self = (μ₀/2π) * l * [ln(2l/GMD) - 1]
    //   where GMD = geometric mean distance of wire cross-section
    //
    // Mutual inductance between segments (Neumann formula):
    //   M_ij = (μ₀/4π) * ∮∮ (dl_i · dl_j) / r_ij
    //   Simplified for short segments: M_ij ≈ (μ₀/4π) * l_i * l_j * cos(θ) / r_ij
    //
    // Matrix equation: (jωL + R) * I = V
    //   where L is N_seg x N_seg inductance matrix
    //         R is diagonal resistance matrix
    //         I is segment current vector
    //         V is voltage excitation
    //
    // Total loop inductance: L_total = Σ_i Σ_j L_ij (all segments connected in series)
    // ==========================================================================

    // Count total panels and determine segment structure
    int totalPanels = 0;
    for (const auto& cond : conductors_) {
        totalPanels += cond->NumPanels();
    }

    if (totalPanels == 0) return;

    // Collect all panels
    std::vector<SurfacePanel> allPanels;
    std::vector<double> allConductivities;
    std::vector<int> panelToConductor;

    for (size_t c = 0; c < conductors_.size(); ++c) {
        const auto& cond = conductors_[c];
        const auto& panels = cond->GetPanels();
        for (const auto& panel : panels) {
            allPanels.push_back(panel);
            allConductivities.push_back(cond->GetConductivity());
            panelToConductor.push_back(static_cast<int>(c));
        }
    }

    // Physical parameters
    double omega = 2.0 * RadConst::PI * frequency_;
    std::complex<double> jOmega(0, omega);

    // ==========================================================================
    // Determine segment structure from panel geometry
    // ==========================================================================
    // For CreateLoop/CreateWire: panels are arranged around wire circumference
    // nPanelsAround = number of panels forming the wire cross-section perimeter
    // nSegments = number of segments along the wire length

    int nPanelsAround = 8;  // Typical value from CreateLoop/CreateWire
    int nSegments = totalPanels / nPanelsAround;
    if (nSegments < 1) nSegments = 1;

    // Compute geometric properties from panels
    double totalArea = 0;
    for (const auto& p : allPanels) totalArea += p.area;
    double avgPanelArea = totalArea / totalPanels;

    // Estimate panel dimensions from geometry
    // Panels are quadrilaterals: we need to distinguish circumferential vs axial dimension
    // Method: compute average edge lengths to detect aspect ratio

    double avgCircumEdge = 0;  // Edge in circumferential direction (around wire)
    double avgAxialEdge = 0;   // Edge in axial direction (along wire)
    int edgeCount = 0;

    for (const auto& p : allPanels) {
        if (p.vertices.size() >= 4) {
            // Compute 4 edge lengths
            double edges[4];
            for (int e = 0; e < 4; ++e) {
                int e1 = (e + 1) % 4;
                double dx = p.vertices[e1].x - p.vertices[e].x;
                double dy = p.vertices[e1].y - p.vertices[e].y;
                double dz = p.vertices[e1].z - p.vertices[e].z;
                edges[e] = std::sqrt(dx*dx + dy*dy + dz*dz);
            }
            // Opposite edges should be similar; shorter pair is circumferential
            double avgEdge01 = (edges[0] + edges[2]) / 2.0;  // edges 0 and 2
            double avgEdge12 = (edges[1] + edges[3]) / 2.0;  // edges 1 and 3
            if (avgEdge01 < avgEdge12) {
                avgCircumEdge += avgEdge01;
                avgAxialEdge += avgEdge12;
            } else {
                avgCircumEdge += avgEdge12;
                avgAxialEdge += avgEdge01;
            }
            edgeCount++;
        }
    }

    if (edgeCount > 0) {
        avgCircumEdge /= edgeCount;
        avgAxialEdge /= edgeCount;
    }

    // Wire perimeter from circumferential panel edge
    double wirePerimeter = nPanelsAround * avgCircumEdge;

    // Wire cross-section area (circular wire)
    double wireArea = wirePerimeter * wirePerimeter / (4.0 * RadConst::PI);

    // Wire radius (equivalent circular)
    double wireRadius = wirePerimeter / (2.0 * RadConst::PI);

    // GMD for circular wire: GMD = r * exp(-1/4) ≈ 0.7788 * r
    double GMD = 0.7788 * wireRadius;

    // Total wire length from axial panel edge
    double wireLength = nSegments * avgAxialEdge;

    // Segment length
    double segLength = avgAxialEdge;

    // ==========================================================================
    // Build segment data structures
    // ==========================================================================
    struct Segment {
        TVector3d center;       // Segment center point
        TVector3d direction;    // Unit tangent direction
        double length;          // Segment length
        double conductivity;    // Material conductivity
    };

    std::vector<Segment> segments(nSegments);

    // First pass: compute all segment centers
    for (int seg = 0; seg < nSegments; ++seg) {
        TVector3d center(0, 0, 0);
        int startPanel = seg * nPanelsAround;
        int endPanel = std::min(startPanel + nPanelsAround, totalPanels);
        int nInRing = endPanel - startPanel;

        for (int p = startPanel; p < endPanel; ++p) {
            center.x += allPanels[p].center.x;
            center.y += allPanels[p].center.y;
            center.z += allPanels[p].center.z;
        }
        if (nInRing > 0) {
            center.x /= nInRing;
            center.y /= nInRing;
            center.z /= nInRing;
        }

        segments[seg].center = center;
        segments[seg].length = segLength;
        segments[seg].conductivity = allConductivities[startPanel];
    }

    // Second pass: compute segment directions from centers
    for (int seg = 0; seg < nSegments; ++seg) {
        TVector3d dir(0, 0, 0);

        if (nSegments == 1) {
            // Single segment: use default direction
            dir = TVector3d(1, 0, 0);
        } else if (seg == 0) {
            // First segment: use direction to next
            dir.x = segments[1].center.x - segments[0].center.x;
            dir.y = segments[1].center.y - segments[0].center.y;
            dir.z = segments[1].center.z - segments[0].center.z;
        } else if (seg == nSegments - 1) {
            // Last segment: use direction from previous
            dir.x = segments[seg].center.x - segments[seg-1].center.x;
            dir.y = segments[seg].center.y - segments[seg-1].center.y;
            dir.z = segments[seg].center.z - segments[seg-1].center.z;
        } else {
            // Middle segment: use direction from prev to next (central difference)
            dir.x = segments[seg+1].center.x - segments[seg-1].center.x;
            dir.y = segments[seg+1].center.y - segments[seg-1].center.y;
            dir.z = segments[seg+1].center.z - segments[seg-1].center.z;
        }

        // Normalize
        double len = std::sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
        if (len > 1e-12) {
            dir.x /= len;
            dir.y /= len;
            dir.z /= len;
        } else {
            dir = TVector3d(1, 0, 0);  // Fallback
        }

        segments[seg].direction = dir;
    }

    // ==========================================================================
    // Build N_seg x N_seg inductance matrix
    // ==========================================================================
    int N = nSegments;  // Unknown: segment currents I_seg

    systemMatrix_.resize(N * N, std::complex<double>(0, 0));
    rhs_.resize(N, std::complex<double>(0, 0));
    solution_.resize(N, std::complex<double>(0, 0));

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < N; ++i) {
        const Segment& segI = segments[i];
        double l_i = segI.length;

        // Resistance of segment: R = l / (sigma * A)
        double R_seg = 0;
        if (segI.conductivity > 0) {
            if (frequency_ > 0) {
                // AC resistance with skin effect
                double skinDepth = std::sqrt(2.0 / (omega * MU_0 * segI.conductivity));
                // Effective area for skin effect (approximate)
                double effectiveArea = wirePerimeter * skinDepth;
                if (effectiveArea > wireArea) effectiveArea = wireArea;  // DC limit
                R_seg = l_i / (segI.conductivity * effectiveArea);
            } else {
                // DC resistance
                R_seg = l_i / (segI.conductivity * wireArea);
            }
        }

        for (int j = 0; j < N; ++j) {
            const Segment& segJ = segments[j];
            double l_j = segJ.length;

            std::complex<double> L_ij;

            if (i == j) {
                // Self-inductance of segment (partial self-inductance)
                // L_self = (μ₀/2π) * l * [ln(2l/GMD) - 1]
                double ln_term = std::log(2.0 * l_i / GMD);
                if (ln_term < 0.25) ln_term = 0.25;  // Prevent negative for short segments
                double L_self = (MU_0 / (2.0 * RadConst::PI)) * l_i * (ln_term - 1.0);
                if (L_self < 0) L_self = (MU_0 / (2.0 * RadConst::PI)) * l_i * 0.25;  // Minimum
                L_ij = std::complex<double>(L_self, 0);
            } else {
                // Mutual inductance between segments (Neumann formula)
                // M_ij = (μ₀/4π) * l_i * l_j * cos(θ_ij) / r_ij

                // Distance between segment centers
                double dx = segI.center.x - segJ.center.x;
                double dy = segI.center.y - segJ.center.y;
                double dz = segI.center.z - segJ.center.z;
                double r = std::sqrt(dx*dx + dy*dy + dz*dz);

                // Prevent division by zero
                if (r < 0.001 * l_i) r = 0.001 * l_i;

                // cos(θ) = dot product of unit tangent vectors
                double cos_theta = segI.direction.x * segJ.direction.x +
                                   segI.direction.y * segJ.direction.y +
                                   segI.direction.z * segJ.direction.z;

                // Mutual inductance
                // M = (μ₀/4π) * l_i * l_j * cos(θ) / r
                double M_ij = (MU_0 * INV_FOUR_PI) * l_i * l_j * cos_theta / r;

                // For full-wave, add phase factor exp(-jkr)
                if (formulation_ == ConductorFormulation::FullWave ||
                    formulation_ == ConductorFormulation::EMQS) {
                    double k = omega * std::sqrt(MU_0 * EPS_0);
                    std::complex<double> phase = std::exp(std::complex<double>(0, -k * r));
                    L_ij = M_ij * phase;
                } else {
                    L_ij = std::complex<double>(M_ij, 0);
                }
            }

            // Matrix element: Z_ij = jωL_ij + R*δ_ij
            std::complex<double> Z_ij = jOmega * L_ij;
            if (i == j) {
                Z_ij += std::complex<double>(R_seg, 0);
            }

            systemMatrix_[i * N + j] = Z_ij;
        }
    }

    // ==========================================================================
    // For series-connected loop coil: Direct impedance calculation
    // ==========================================================================
    //
    // For a loop coil, all segments are in series carrying the same current I.
    // The total impedance is:
    //
    //   Z_total = Σ_i Σ_j Z_ij = Σ_i Σ_j (jωL_ij + R_i*δ_ij)
    //
    // This is equivalent to solving the scalar equation:
    //   Z_total * I = V
    //
    // Rather than solving an N x N system (which would give different I_i),
    // we directly compute Z_total by summing all matrix elements.
    //
    // Note: This approach is correct for single-path conductors (loops, spirals).
    // For complex topologies with parallel branches, the full matrix solve is needed.

    // Compute total impedance by summing all matrix elements
    std::complex<double> Z_total(0, 0);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            Z_total += systemMatrix_[i * N + j];
        }
    }

    // Store in solution vector for ComputePortImpedance to retrieve
    // We use a special format: solution_[0] = Z_total (impedance directly)
    // and set a flag to indicate direct impedance mode

    // For voltage excitation, compute current: I = V / Z_total
    std::complex<double> V_applied(1.0, 0);  // Default 1V
    for (size_t c = 0; c < conductors_.size(); ++c) {
        auto& cond = conductors_[c];
        if (cond->GetExcitationType() == 1) {
            V_applied = cond->GetExcitationValue();
        }
    }

    // Store impedance and current in solution
    // solution_[0] = current I = V / Z_total
    if (std::abs(Z_total) > 1e-30) {
        std::complex<double> I = V_applied / Z_total;
        for (int i = 0; i < N; ++i) {
            solution_[i] = I;  // All segments have same current
        }
    }

    // RHS is not used for direct solve, but set for consistency
    for (int i = 0; i < N; ++i) {
        rhs_[i] = V_applied / static_cast<double>(N);
    }

    // Set flag to skip SolveLinearSystem (solution already computed)
    directImpedanceMode_ = true;
}

void radTConductorSolver::SolveLinearSystem() {
    // Skip if using direct impedance mode (solution already computed in BuildSystemMatrix)
    if (directImpedanceMode_) {
        return;
    }

    int nDOF = static_cast<int>(rhs_.size());
    if (nDOF == 0) return;

    // Copy for LAPACK (column-major)
    std::vector<std::complex<double>> A = systemMatrix_;
    solution_ = rhs_;

    // LAPACK zgesv: solve A*x = b
    std::vector<int> ipiv(nDOF);

    // Use MKL LAPACK
    #ifdef HAVE_LAPACK
    lapack_int info = LAPACKE_zgesv(LAPACK_COL_MAJOR, nDOF, 1,
                                     reinterpret_cast<lapack_complex_double*>(A.data()),
                                     nDOF,
                                     ipiv.data(),
                                     reinterpret_cast<lapack_complex_double*>(solution_.data()),
                                     nDOF);
    if (info != 0) {
        throw std::runtime_error("LAPACK zgesv failed with info = " + std::to_string(info));
    }
    #else
    // Fallback: simple Gaussian elimination (for testing)
    // Note: This should be replaced with proper LAPACK in production
    for (int k = 0; k < nDOF; ++k) {
        // Find pivot
        int maxRow = k;
        double maxVal = std::abs(A[k * nDOF + k]);
        for (int i = k + 1; i < nDOF; ++i) {
            double val = std::abs(A[i * nDOF + k]);
            if (val > maxVal) {
                maxVal = val;
                maxRow = i;
            }
        }

        // Swap rows
        if (maxRow != k) {
            for (int j = 0; j < nDOF; ++j) {
                std::swap(A[k * nDOF + j], A[maxRow * nDOF + j]);
            }
            std::swap(solution_[k], solution_[maxRow]);
        }

        // Eliminate
        for (int i = k + 1; i < nDOF; ++i) {
            std::complex<double> factor = A[i * nDOF + k] / A[k * nDOF + k];
            for (int j = k; j < nDOF; ++j) {
                A[i * nDOF + j] -= factor * A[k * nDOF + j];
            }
            solution_[i] -= factor * solution_[k];
        }
    }

    // Back substitution
    for (int i = nDOF - 1; i >= 0; --i) {
        for (int j = i + 1; j < nDOF; ++j) {
            solution_[i] -= A[i * nDOF + j] * solution_[j];
        }
        solution_[i] /= A[i * nDOF + i];
    }
    #endif
}

void radTConductorSolver::ExtractSolution() {
    // Darwin approximation: only K unknowns (no charge σ)
    // Solution ordering: [K_0, K_1, ..., K_N]

    int totalPanels = 0;
    for (const auto& cond : conductors_) {
        totalPanels += cond->NumPanels();
    }

    if (totalPanels == 0) return;

    int panelOffset = 0;

    for (auto& cond : conductors_) {
        int nPanels = cond->NumPanels();
        auto& K = cond->SurfaceCurrent();
        auto& sigma = cond->SurfaceCharge();

        const auto& panels = cond->GetPanels();

        // Resize storage if needed
        if (static_cast<int>(K.size()) != 3 * nPanels) {
            K.resize(3 * nPanels);
        }
        if (static_cast<int>(sigma.size()) != nPanels) {
            sigma.resize(nPanels);
        }

        for (int i = 0; i < nPanels; ++i) {
            int globalIdx = panelOffset + i;

            // Get scalar current magnitude from solution
            std::complex<double> K_mag = solution_[globalIdx];

            // Convert to vector current (along tangent direction)
            // For wire-like structures, current flows along the wire
            // Use panel normal to determine tangent
            const auto& panel = panels[i];

            // Simple model: current flows perpendicular to normal
            // Project onto tangent plane: K_vec = K_mag * (I - n⊗n) · e_x
            // This gives current flowing in x-direction projected onto surface
            K[3*i] = K_mag * (1.0 - panel.normal.x * panel.normal.x);
            K[3*i + 1] = K_mag * (-panel.normal.x * panel.normal.y);
            K[3*i + 2] = K_mag * (-panel.normal.x * panel.normal.z);

            // Darwin approximation: σ not directly computed
            // Set to zero (can be estimated from div(K) if needed)
            sigma[i] = std::complex<double>(0, 0);
        }

        panelOffset += nPanels;
    }

    // Compute port impedance
    ComputePortImpedance();
}

void radTConductorSolver::ComputePortImpedance() {
    // ==========================================================================
    // Compute impedance from segment-based solution
    // ==========================================================================
    //
    // For segment-based solution:
    //   solution_[i] = I_seg (segment current in Amperes)
    //   Matrix equation: [Z_mat] I = V
    //   where Z_mat = jωL + R (N_seg x N_seg)
    //
    // Impedance calculation:
    //   Z = V_total / I (for series-connected segments, I is same for all)
    //
    // Since segments are in series for a single coil:
    //   - Total voltage: V_total = applied voltage
    //   - Current: I = I_seg (same in all segments)
    //   - Impedance: Z = V / I
    //
    // Verification via energy method:
    //   Complex power: P = V · I*
    //   Impedance: Z = V² / P* = V / I

    int nSegments = static_cast<int>(solution_.size());

    for (auto& cond : conductors_) {
        int nPanels = cond->NumPanels();
        if (nPanels == 0) continue;

        // Get excitation voltage
        std::complex<double> V_applied(1.0, 0);  // Default 1V
        if (cond->GetExcitationType() == 1) {
            V_applied = cond->GetExcitationValue();
        }

        // For series-connected segments, current should be same in all segments
        // Use average of segment currents (should all be equal for proper solution)
        std::complex<double> I_avg(0, 0);
        if (nSegments > 0) {
            for (int i = 0; i < nSegments; ++i) {
                I_avg += solution_[i];
            }
            I_avg /= static_cast<double>(nSegments);
        }

        // Impedance: Z = V / I
        std::complex<double> Z(0, 0);
        if (std::abs(I_avg) > 1e-30) {
            Z = V_applied / I_avg;
        }

        // Store results
        cond->SetPortImpedance(Z);
        cond->SetTotalCurrent(I_avg);
    }
}

// ============================================================================
// Loop-Star Decomposition Implementation
// ============================================================================
//
// Loop-Star decomposition separates the current into solenoidal (loop) and
// non-solenoidal (star) components:
//   J = J_loop + J_star
//
// where:
//   div(J_loop) = 0     (solenoidal, no charge contribution)
//   curl(J_star) = 0    (irrotational, contributes to charge)
//
// This decomposition improves low-frequency stability because:
//   - Loop basis: Only contributes to L matrix (inductance)
//   - Star basis: Only contributes to P matrix (capacitance/charge)
//
// For wire-type conductors (loops, spirals), the current is mostly solenoidal,
// so we use a simplified approach based on the wire path connectivity.

void radTConductorSolver::BuildEdgeConnectivity() {
    // Build edge connectivity for all panels
    // This identifies shared edges between adjacent panels

    edges_.clear();

    // Collect all panels
    std::vector<SurfacePanel> allPanels;
    for (const auto& cond : conductors_) {
        const auto& panels = cond->GetPanels();
        for (const auto& panel : panels) {
            allPanels.push_back(panel);
        }
    }

    int N = static_cast<int>(allPanels.size());
    if (N == 0) return;

    // Tolerance for vertex matching
    double tol = 1e-10;

    // For each panel, extract edges and find matching edges in other panels
    for (int i = 0; i < N; ++i) {
        const auto& panelI = allPanels[i];
        int nVertI = static_cast<int>(panelI.vertices.size());

        for (int ei = 0; ei < nVertI; ++ei) {
            // Edge from vertex ei to vertex (ei+1) mod nVertI
            TVector3d v1 = panelI.vertices[ei];
            TVector3d v2 = panelI.vertices[(ei + 1) % nVertI];

            // Check if this edge already exists in edges_ list
            bool found = false;
            for (auto& edge : edges_) {
                if (edge.panel1 == i && edge.localIdx1 == ei) {
                    found = true;
                    break;
                }
                if (edge.panel2 == i && edge.localIdx2 == ei) {
                    found = true;
                    break;
                }
            }

            if (found) continue;

            // Search for matching edge in other panels
            int matchPanel = -1;
            int matchEdge = -1;

            for (int j = i + 1; j < N; ++j) {
                const auto& panelJ = allPanels[j];
                int nVertJ = static_cast<int>(panelJ.vertices.size());

                for (int ej = 0; ej < nVertJ; ++ej) {
                    TVector3d u1 = panelJ.vertices[ej];
                    TVector3d u2 = panelJ.vertices[(ej + 1) % nVertJ];

                    // Check if edges match (same vertices, possibly reversed)
                    TVector3d d1a = v1 - u1;
                    TVector3d d1b = v2 - u2;
                    TVector3d d2a = v1 - u2;
                    TVector3d d2b = v2 - u1;
                    double dist1 = d1a.Abs() + d1b.Abs();
                    double dist2 = d2a.Abs() + d2b.Abs();

                    if (dist1 < tol || dist2 < tol) {
                        matchPanel = j;
                        matchEdge = ej;
                        break;
                    }
                }

                if (matchPanel >= 0) break;
            }

            // Create edge entry
            Edge edge;
            edge.panel1 = i;
            edge.localIdx1 = ei;
            edge.panel2 = matchPanel;
            edge.localIdx2 = matchEdge;
            edge.midpoint = TVector3d(
                0.5 * (v1.x + v2.x),
                0.5 * (v1.y + v2.y),
                0.5 * (v1.z + v2.z)
            );
            TVector3d dir = v2 - v1;
            edge.length = dir.Abs();
            if (edge.length > 1e-15) {
                edge.direction = TVector3d(
                    dir.x / edge.length,
                    dir.y / edge.length,
                    dir.z / edge.length
                );
            } else {
                edge.direction = TVector3d(1, 0, 0);
            }

            edges_.push_back(edge);
        }
    }
}

void radTConductorSolver::BuildLoopStarBasis() {
    // Build Loop and Star basis functions
    //
    // For a mesh with N panels and E edges:
    //   - Number of internal edges: E_int (edges shared by 2 panels)
    //   - Number of boundary edges: E_bnd (edges on mesh boundary)
    //   - Number of vertices: V
    //
    // By Euler's formula for planar graphs:
    //   V - E + F = 2  (for simply connected surface)
    //
    // Loop basis dimension: N_loop = E_int - N + 1 (independent loops)
    // Star basis dimension: N_star = N - 1 (tree edges connecting panels)
    //
    // For wire-type conductors (closed loop):
    //   - One global loop (current flowing around the wire)
    //   - N-1 "star" basis (local current redistribution)

    // Count internal and boundary edges
    int numInternal = 0;
    int numBoundary = 0;

    for (const auto& edge : edges_) {
        if (edge.panel2 >= 0) {
            numInternal++;
        } else {
            numBoundary++;
        }
    }

    int N = 0;
    for (const auto& cond : conductors_) {
        N += cond->NumPanels();
    }

    // For simply connected surface:
    // N_loop = E_int - N + 1
    // N_star = N - 1

    numLoops_ = numInternal - N + 1;
    if (numLoops_ < 0) numLoops_ = 0;
    numStars_ = N - 1;
    if (numStars_ < 0) numStars_ = 0;

    // Simplified approach for wire-type conductors:
    // Assume the primary current mode is a single loop around the wire
    // Additional modes are local perturbations

    // Initialize basis matrices
    loopBasis_.clear();
    starBasis_.clear();

    // For now, use identity basis (no decomposition)
    // This maintains compatibility with the current solver
    // Full Loop-Star implementation requires:
    // 1. Graph-based loop detection (cycle finding)
    // 2. Spanning tree construction for star basis

    // Placeholder: Use panel-based basis
    loopBasis_.resize(N, std::vector<double>(N, 0.0));
    starBasis_.resize(N, std::vector<double>(N, 0.0));

    for (int i = 0; i < N; ++i) {
        loopBasis_[i][i] = 1.0;  // Identity for now
        starBasis_[i][i] = 1.0;
    }
}

void radTConductorSolver::BuildSystemMatrixLoopStar() {
    // Build system matrix with Loop-Star decomposition
    //
    // The standard EFIE:
    //   [jωL + R + (1/jω)P] J = E_inc
    //
    // With Loop-Star decomposition:
    //   J = T_L * J_L + T_S * J_S
    //
    // where T_L and T_S are transformation matrices
    //
    // The decomposed system:
    //   [T_L^T * (jωL + R) * T_L    T_L^T * (jωL + R) * T_S  ] [J_L]   [T_L^T * E_inc]
    //   [T_S^T * (jωL + R) * T_L    T_S^T * (jωL + R + P/jω) * T_S] [J_S] = [T_S^T * E_inc]
    //
    // Key insight: P matrix only affects Star-Star block
    // This improves conditioning at low frequency

    // For now, fall back to standard Darwin approximation
    // Full Loop-Star requires proper T_L and T_S matrices

    // Build edge connectivity if not done
    if (edges_.empty()) {
        BuildEdgeConnectivity();
    }

    // Build Loop-Star basis
    BuildLoopStarBasis();

    // Fall back to standard method until full implementation
    // (The Loop-Star basis construction is complex and requires
    // graph algorithms for cycle detection)

    // For wire-type conductors at low frequency, the Darwin approximation
    // with proper self-term handling should be sufficient

    // Call standard BuildSystemMatrix
    BuildSystemMatrix();
}

} // namespace radia
