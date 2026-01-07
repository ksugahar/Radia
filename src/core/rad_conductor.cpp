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

// ============================================================================
// radTConductorSolver implementation
// ============================================================================

radTConductorSolver::radTConductorSolver()
    : frequency_(0)
    , formulation_(ConductorFormulation::MQS)
    , usePfft_(true)
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
    // Count total unknowns
    int totalPanels = 0;
    for (const auto& cond : conductors_) {
        totalPanels += cond->NumPanels();
    }

    if (totalPanels == 0) return;

    // FastImp formulation: EFIE + continuity equation
    // Unknowns per panel: scalar surface current magnitude + charge
    // For each panel i:
    //   - K_i: tangential surface current density [A/m]
    //   - sigma_i: surface charge density [C/m^2]
    //
    // Equations:
    // 1. EFIE on conductor surface: E_tan = 0 (PEC) or E_tan = Z_s * K (impedance BC)
    //    n x (E_inc + E_scat) = n x (Z_s * K)
    //    where E_scat = -j*omega*A - grad(Phi)
    //
    // 2. Charge conservation: div_s(K) + j*omega*sigma = 0
    //
    // Matrix structure (simplified scalar formulation):
    // [L + R   P  ] [K    ]   [V_inc]
    // [D      C  ] [sigma] = [0    ]
    //
    // L: inductance matrix from A potential
    // R: resistance matrix (skin effect)
    // P: gradient of scalar potential
    // D: surface divergence
    // C: capacitance matrix from charge

    int nDOF = 2 * totalPanels;  // K and sigma for each panel

    systemMatrix_.resize(nDOF * nDOF, std::complex<double>(0, 0));
    rhs_.resize(nDOF, std::complex<double>(0, 0));
    solution_.resize(nDOF, std::complex<double>(0, 0));

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
    std::complex<double> jOmegaMu = jOmega * MU_0;
    std::complex<double> invJOmegaEps = (omega > 1e-10) ?
        1.0 / (jOmega * EPS_0) : std::complex<double>(0, 0);

    // Build matrix using pFFT if enabled and beneficial
    bool usePfftAccel = usePfft_ && totalPanels > 100;

    if (usePfftAccel) {
        std::vector<TVector3d> allCenters;
        std::vector<double> allAreas;

        for (const auto& panel : allPanels) {
            allCenters.push_back(panel.center);
            allAreas.push_back(panel.area);
        }

        pfft_ = std::make_unique<radTPfft>();
        pfft_->Initialize(allCenters, allAreas);

        if (formulation_ == ConductorFormulation::FullWave ||
            formulation_ == ConductorFormulation::EMQS) {
            pfft_->SetupKernelFullWave(frequency_, EPS_0, MU_0);
        } else {
            pfft_->SetupKernelMQS();
        }
    }

    // Fill system matrix
    // Row ordering: [K_0, sigma_0, K_1, sigma_1, ...]
    // Or block form: [K_0..K_N, sigma_0..sigma_N]

    // Using block form for clarity
    // Block (0,0): L + R (N x N) - inductance + resistance
    // Block (0,1): P (N x N) - scalar potential gradient
    // Block (1,0): D (N x N) - surface divergence
    // Block (1,1): C (N x N) - capacitance

    int N = totalPanels;

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < N; ++i) {
        const SurfacePanel& panelI = allPanels[i];
        double sigmaI = allConductivities[i];

        // Surface impedance (skin effect)
        double Rs = 0;
        if (frequency_ > 0 && sigmaI > 0) {
            double skinDepth = std::sqrt(2.0 / (omega * MU_0 * sigmaI));
            Rs = 1.0 / (sigmaI * skinDepth);  // DC + skin effect
        }

        for (int j = 0; j < N; ++j) {
            const SurfacePanel& panelJ = allPanels[j];

            // Distance between panel centers
            double dx = panelI.center.x - panelJ.center.x;
            double dy = panelI.center.y - panelJ.center.y;
            double dz = panelI.center.z - panelJ.center.z;
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            // Green's function
            std::complex<double> G;
            std::complex<double> dGdr;

            if (r < 1e-12) {
                // Self-term: use analytical approximation
                double R_eff = std::sqrt(panelJ.area / RadConst::PI);
                G = std::complex<double>(R_eff * (2.0 * std::log(2.0) - 1.0) * INV_FOUR_PI, 0);
                dGdr = std::complex<double>(0, 0);
            } else {
                if (formulation_ == ConductorFormulation::FullWave ||
                    formulation_ == ConductorFormulation::EMQS) {
                    // Full-wave Green's function
                    double k = omega * std::sqrt(MU_0 * EPS_0);
                    std::complex<double> jkr(0, k * r);
                    std::complex<double> expjkr = std::exp(-jkr);
                    G = expjkr * INV_FOUR_PI / r;
                    dGdr = -expjkr * (1.0 + jkr) * INV_FOUR_PI / (r * r);
                } else {
                    // MQS Green's function
                    G = std::complex<double>(INV_FOUR_PI / r, 0);
                    dGdr = std::complex<double>(-INV_FOUR_PI / (r * r), 0);
                }
            }

            // L matrix: inductance from vector potential
            // L_ij = mu * integral{ G(r,r') } dS_j
            std::complex<double> L_ij = MU_0 * G * panelJ.area;

            // R matrix: resistance (diagonal only for simplified model)
            std::complex<double> R_ij = (i == j) ?
                std::complex<double>(Rs * panelI.area, 0) : std::complex<double>(0, 0);

            // Block (0,0): jOmega*L + R
            int idx00 = i * (2*N) + j;
            systemMatrix_[idx00] = jOmega * L_ij + R_ij;

            // P matrix: scalar potential gradient
            // P_ij = (1/eps) * integral{ grad(G) . n } dS_j
            // Simplified: use centroid approximation
            double n_dot_r = 0;
            if (r > 1e-12) {
                n_dot_r = (panelI.normal.x * dx + panelI.normal.y * dy +
                           panelI.normal.z * dz) / r;
            }
            std::complex<double> P_ij = dGdr * n_dot_r * panelJ.area / EPS_0;

            // Block (0,1): scalar potential contribution
            int idx01 = i * (2*N) + (N + j);
            systemMatrix_[idx01] = -P_ij;  // E = -grad(Phi)

            // D matrix: surface divergence (connectivity)
            // div_s(K) relates currents between adjacent panels
            // Simplified: use area ratio for neighboring panels
            std::complex<double> D_ij;
            if (i == j) {
                D_ij = std::complex<double>(1.0 / panelI.area, 0);
            } else {
                // Off-diagonal: check if panels are neighbors
                // For now, use simple geometric proximity
                double charSize = std::sqrt(panelI.area);
                if (r < 2.0 * charSize) {
                    D_ij = std::complex<double>(-1.0 / (2.0 * N * panelJ.area), 0);
                } else {
                    D_ij = std::complex<double>(0, 0);
                }
            }

            // Block (1,0): divergence operator
            int idx10 = (N + i) * (2*N) + j;
            systemMatrix_[idx10] = D_ij;

            // C matrix: capacitance from charge
            // C_ij = integral{ G(r,r') } dS_j
            std::complex<double> C_ij = G * panelJ.area;

            // Block (1,1): jOmega * C (from continuity: div_s(K) + jOmega*sigma = 0)
            int idx11 = (N + i) * (2*N) + (N + j);
            systemMatrix_[idx11] = jOmega * C_ij;
        }
    }

    // Setup RHS from port excitation
    // Support both voltage and current excitation modes
    int panelOffset = 0;
    for (size_t c = 0; c < conductors_.size(); ++c) {
        auto& cond = conductors_[c];
        int nPanels = cond->NumPanels();

        if (nPanels == 0) continue;

        int excType = cond->GetExcitationType();
        std::complex<double> excValue = cond->GetExcitationValue();

        if (excType == 1) {
            // Voltage excitation: apply voltage across conductor
            // V_applied appears in the EFIE equation
            // For a wire conductor, voltage is applied at the port

            // Distribute voltage source across first few panels (port terminals)
            int numPortPanels = std::min(nPanels, 4);  // Port size
            std::complex<double> V_per_panel = excValue / static_cast<double>(numPortPanels);

            for (int i = 0; i < numPortPanels; ++i) {
                rhs_[panelOffset + i] = V_per_panel;
            }

        } else if (excType == 2) {
            // Current excitation: impose total current constraint
            // This modifies the system to enforce I_total = I_specified

            // For current-driven analysis, we add a constraint equation
            // Sum of currents across a cross-section equals I_specified

            // Get total cross-sectional area for normalization
            double totalArea = 0;
            const auto& panels = cond->GetPanels();
            for (int i = 0; i < nPanels; ++i) {
                totalArea += panels[i].area;
            }

            // Distribute current source: J = I / A
            // The current density K [A/m] on each panel contributes to total current
            // For a wire, I_total = integral(K * dl) around circumference
            // Simplified: distribute uniformly

            std::complex<double> K_avg = excValue / std::sqrt(totalArea);  // Approximate uniform distribution

            // Set as Dirichlet-like constraint on surface current
            for (int i = 0; i < nPanels; ++i) {
                // Modify matrix to enforce current level
                // Add penalty term to enforce average current
                int idx = panelOffset + i;
                systemMatrix_[idx * (2*N) + idx] += std::complex<double>(1e6, 0);  // Strong penalty
                rhs_[idx] = K_avg * std::complex<double>(1e6, 0);
            }

            // Store excitation current for later retrieval
            cond->SetTotalCurrent(excValue);

        } else {
            // No excitation: apply default 1V for impedance calculation
            if (nPanels > 0) {
                rhs_[panelOffset] = std::complex<double>(1.0, 0);
            }
        }

        panelOffset += nPanels;
    }
}

void radTConductorSolver::SolveLinearSystem() {
    int nDOF = static_cast<int>(rhs_.size());
    if (nDOF == 0) return;

    // Copy for LAPACK (column-major)
    std::vector<std::complex<double>> A = systemMatrix_;
    solution_ = rhs_;

    // LAPACK zgesv: solve A*x = b
    std::vector<int> ipiv(nDOF);

    // Use MKL LAPACK
    #ifdef USE_MKL
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
    // Extract K and sigma from solution vector
    // Solution ordering: [K_0..K_N, sigma_0..sigma_N]

    int totalPanels = 0;
    for (const auto& cond : conductors_) {
        totalPanels += cond->NumPanels();
    }

    if (totalPanels == 0) return;

    int N = totalPanels;
    int panelOffset = 0;

    for (auto& cond : conductors_) {
        int nPanels = cond->NumPanels();
        auto& K = cond->SurfaceCurrent();
        auto& sigma = cond->SurfaceCharge();

        const auto& panels = cond->GetPanels();

        for (int i = 0; i < nPanels; ++i) {
            int globalIdx = panelOffset + i;

            // Get scalar current magnitude
            std::complex<double> K_mag = solution_[globalIdx];

            // Convert to vector current (along tangent direction)
            // For wire-like structures, current flows along the wire
            // Use panel normal to determine tangent
            const auto& panel = panels[i];

            // Simple model: current flows perpendicular to normal
            // For more accuracy, track wire direction
            K[3*i] = K_mag * (1.0 - panel.normal.x * panel.normal.x);
            K[3*i + 1] = K_mag * (-panel.normal.x * panel.normal.y);
            K[3*i + 2] = K_mag * (-panel.normal.x * panel.normal.z);

            // Get charge density
            sigma[i] = solution_[N + globalIdx];
        }

        panelOffset += nPanels;
    }

    // Compute port impedance
    // Z = V / I, where V = 1V (applied), I = integral of K over port
    ComputePortImpedance();
}

void radTConductorSolver::ComputePortImpedance() {
    // Compute impedance from port terminals
    for (auto& cond : conductors_) {
        // Sum current at port terminals
        std::complex<double> I_total(0, 0);
        const auto& K = cond->SurfaceCurrent();
        const auto& panels = cond->GetPanels();

        // Simple model: sum all currents (for single port)
        for (size_t i = 0; i < panels.size(); ++i) {
            double Kx = std::abs(K[3*i]);
            double Ky = std::abs(K[3*i + 1]);
            double Kz = std::abs(K[3*i + 2]);
            double K_mag = std::sqrt(Kx*Kx + Ky*Ky + Kz*Kz);
            I_total += std::complex<double>(K_mag * panels[i].area, 0);
        }

        // Z = V / I (V = 1V applied)
        if (std::abs(I_total) > 1e-15) {
            // Store in conductor's port impedance
            // Note: This is simplified; proper port extraction needs more work
        }
    }
}

} // namespace radia
