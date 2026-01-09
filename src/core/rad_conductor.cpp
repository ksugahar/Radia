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

    // Store wire centerline path for Biot-Savart field computation
    wireCenterPath_ = path;

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

    if (panels_.empty()) return;

    // If totalCurrent_ is set (from Solve), use filament Biot-Savart
    // Use stored wire centerline path for accurate field computation
    if (std::abs(totalCurrent_) > 1e-30) {
        const double coeff = MU_0 / (4.0 * RadConst::PI);

        // Use wireCenterPath_ if available (from CreateWire/CreateLoop)
        if (!wireCenterPath_.empty() && wireCenterPath_.size() > 1) {
            // Biot-Savart: B = mu_0/(4*pi) * I * integral{ dl x (r-r') / |r-r'|^3 }
            for (size_t i = 0; i < wireCenterPath_.size() - 1; ++i) {
                // Segment from point i to i+1
                const TVector3d& p0 = wireCenterPath_[i];
                const TVector3d& p1 = wireCenterPath_[i + 1];

                // dl = segment vector
                TVector3d dl = p1 - p0;

                // Midpoint of segment
                TVector3d midPt;
                midPt.x = (p0.x + p1.x) * 0.5;
                midPt.y = (p0.y + p1.y) * 0.5;
                midPt.z = (p0.z + p1.z) * 0.5;

                // Vector from segment midpoint to observation point
                double dx = point.x - midPt.x;
                double dy = point.y - midPt.y;
                double dz = point.z - midPt.z;
                double r = std::sqrt(dx*dx + dy*dy + dz*dz);

                if (r < 1e-12) continue;

                double r3 = r * r * r;

                // dl x (r - r') where (r-r') points from source to obs
                double crossX = dl.y * dz - dl.z * dy;
                double crossY = dl.z * dx - dl.x * dz;
                double crossZ = dl.x * dy - dl.y * dx;

                Bx += totalCurrent_ * coeff * crossX / r3;
                By += totalCurrent_ * coeff * crossY / r3;
                Bz += totalCurrent_ * coeff * crossZ / r3;
            }
        } else {
            // Fallback: use panel centers (less accurate, may have overcounting)
            // This path is used for non-wire conductors (blocks, etc.)
            for (size_t i = 0; i < panels_.size(); ++i) {
                TVector3d dl;
                if (i == 0 && panels_.size() > 1) {
                    dl = panels_[1].center - panels_[0].center;
                } else if (i == panels_.size() - 1 && panels_.size() > 1) {
                    dl = panels_[i].center - panels_[i-1].center;
                } else if (panels_.size() > 2) {
                    dl = panels_[i+1].center - panels_[i-1].center;
                    dl = dl * 0.5;
                } else {
                    continue;
                }

                double dx = point.x - panels_[i].center.x;
                double dy = point.y - panels_[i].center.y;
                double dz = point.z - panels_[i].center.z;
                double r = std::sqrt(dx*dx + dy*dy + dz*dz);

                if (r < 1e-12) continue;

                double r3 = r * r * r;

                double crossX = dl.y * dz - dl.z * dy;
                double crossY = dl.z * dx - dl.x * dz;
                double crossZ = dl.x * dy - dl.y * dx;

                Bx += totalCurrent_ * coeff * crossX / r3;
                By += totalCurrent_ * coeff * crossY / r3;
                Bz += totalCurrent_ * coeff * crossZ / r3;
            }
        }
    } else {
        // Fallback: use surface current if available (for RWG-EFIE solution)
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
            std::complex<double> G = GreenFunction(r);

            // K x r contribution
            double rx = dx / r, ry = dy / r, rz = dz / r;

            Bx += MU_0 * G * (Ky * rz - Kz * ry) * panel.area / r;
            By += MU_0 * G * (Kz * rx - Kx * rz) * panel.area / r;
            Bz += MU_0 * G * (Kx * ry - Ky * rx) * panel.area / r;
        }
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

    // Vector potential: A = mu_0/(4*pi) * integral{ I * dl / |r - r'| }
    // For wire conductors, use filament approximation along wire centerline

    if (std::abs(totalCurrent_) > 1e-30 && !wireCenterPath_.empty() && wireCenterPath_.size() > 1) {
        // Use wire centerline path for filament vector potential
        // A = mu_0/(4*pi) * I * integral{ dl / |r - r'| }
        const double coeff = MU_0 / (4.0 * RadConst::PI);

        for (size_t i = 0; i < wireCenterPath_.size() - 1; ++i) {
            const TVector3d& p0 = wireCenterPath_[i];
            const TVector3d& p1 = wireCenterPath_[i + 1];

            // dl = segment vector
            TVector3d dl = p1 - p0;

            // Midpoint of segment
            TVector3d midPt;
            midPt.x = (p0.x + p1.x) * 0.5;
            midPt.y = (p0.y + p1.y) * 0.5;
            midPt.z = (p0.z + p1.z) * 0.5;

            // Vector from segment midpoint to observation point
            double dx = point.x - midPt.x;
            double dy = point.y - midPt.y;
            double dz = point.z - midPt.z;
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (r < 1e-12) continue;

            // A = mu_0/(4*pi) * I * dl / r
            Ax += totalCurrent_ * coeff * dl.x / r;
            Ay += totalCurrent_ * coeff * dl.y / r;
            Az += totalCurrent_ * coeff * dl.z / r;
        }
    } else {
        // Fallback: use surface current if available (for RWG-EFIE solution)
        // A = mu_0/(4*pi) * integral{ K / |r - r'| } dA'
        for (size_t i = 0; i < panels_.size(); ++i) {
            const auto& panel = panels_[i];

            double dx = point.x - panel.center.x;
            double dy = point.y - panel.center.y;
            double dz = point.z - panel.center.z;
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (r < 1e-15) continue;

            // Get surface current at this panel (3 components)
            if (surfaceCurrent_.size() <= 3*i + 2) continue;
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
    , directImpedanceMode_(false)
    , portImpedance_(0, 0)
    , totalCurrent_(0, 0)
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
    rwgMesh_.reset();
    rwgSolver_.reset();
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
    if (conductors_.empty()) return;

    // For wire conductors: use analytical formula for self-inductance
    // For general conductors: use simple panel-based approximation

    auto& cond = conductors_[0];
    double sigma = cond->GetConductivity();
    double omega = 2.0 * RadConst::PI * frequency_;

    // Skin depth
    double skinDepth = 1e-3;  // default 1mm
    if (omega > 0 && sigma > 0) {
        skinDepth = std::sqrt(2.0 / (omega * MU_0 * sigma));
    }

    // Surface resistance
    double Rs = 1.0 / (sigma * skinDepth);

    const auto& panels = cond->GetPanels();
    int N = static_cast<int>(panels.size());

    if (N == 0) return;

    // Estimate wire/coil geometry from panels
    // Compute bounding box and total path length
    TVector3d minPt = panels[0].center;
    TVector3d maxPt = panels[0].center;

    for (const auto& p : panels) {
        minPt.x = std::min(minPt.x, p.center.x);
        minPt.y = std::min(minPt.y, p.center.y);
        minPt.z = std::min(minPt.z, p.center.z);
        maxPt.x = std::max(maxPt.x, p.center.x);
        maxPt.y = std::max(maxPt.y, p.center.y);
        maxPt.z = std::max(maxPt.z, p.center.z);
    }

    TVector3d size = maxPt - minPt;
    double Lx = size.x, Ly = size.y, Lz = size.z;

    // Total surface area
    double totalArea = 0;
    for (const auto& p : panels) {
        totalArea += p.area;
    }

    // Estimate wire radius from panel area / perimeter
    // For a cylinder: Area = 2*pi*r*L, so r ~ Area/(2*pi*L)
    double pathLength = 0;
    for (size_t i = 1; i < panels.size(); ++i) {
        TVector3d d = panels[i].center - panels[i-1].center;
        pathLength += d.Abs();
    }
    if (pathLength < 1e-10) {
        pathLength = std::max({Lx, Ly, Lz});
    }

    double wireRadius = totalArea / (2.0 * RadConst::PI * pathLength);
    if (wireRadius < 1e-10) wireRadius = 1e-3;  // default 1mm

    // Estimate loop radius from bounding box
    double loopRadius = std::max(Lx, Ly) / 2.0;
    if (loopRadius < 1e-10) loopRadius = 0.05;  // default 50mm

    // Calculate self-inductance using analytical formula
    // For circular loop: L = μ₀ * R * (ln(8R/a) - 2 + Li/4 + a²/(8R²))
    // where Li = 0 for uniform current, 1/4 for skin effect
    double L_self;

    if (Lz < wireRadius) {
        // Flat coil - use circular loop formula
        L_self = MU_0 * loopRadius * (std::log(8.0 * loopRadius / wireRadius) - 2.0);
    } else if (Lz > 4.0 * loopRadius) {
        // Solenoid - use solenoid formula
        // L = μ₀ * n² * V where n = turns/length, V = π*R²*h
        double nTurns = pathLength / (2.0 * RadConst::PI * loopRadius);
        double height = Lz;
        double n = nTurns / height;
        double volume = RadConst::PI * loopRadius * loopRadius * height;
        L_self = MU_0 * n * n * volume;
    } else {
        // General coil - use circular loop as approximation
        L_self = MU_0 * loopRadius * (std::log(8.0 * loopRadius / wireRadius) - 2.0);
    }

    // DC resistance
    double wireCrossSection = RadConst::PI * wireRadius * wireRadius;
    double R_dc = pathLength / (sigma * wireCrossSection);

    // AC resistance (skin effect)
    double R_ac;
    if (skinDepth < wireRadius) {
        // Skin effect dominant: R ~ R_dc * r / (2*delta)
        R_ac = R_dc * wireRadius / (2.0 * skinDepth);
    } else {
        R_ac = R_dc;
    }

    // Impedance: Z = R + jωL
    std::complex<double> jOmega(0, omega);
    portImpedance_ = std::complex<double>(R_ac, 0) + jOmega * L_self;

    // Current from excitation
    int excType = cond->GetExcitationType();
    std::complex<double> excValue = cond->GetExcitationValue();

    if (excType == 2) {
        // Current excitation: use current directly
        totalCurrent_ = excValue;
    } else if (excType == 1) {
        // Voltage excitation: I = V / Z
        if (std::abs(portImpedance_) > 1e-30) {
            totalCurrent_ = excValue / portImpedance_;
        } else {
            totalCurrent_ = std::complex<double>(0, 0);
        }
    } else {
        // No excitation: use default 1V
        if (std::abs(portImpedance_) > 1e-30) {
            totalCurrent_ = std::complex<double>(1.0, 0) / portImpedance_;
        } else {
            totalCurrent_ = std::complex<double>(0, 0);
        }
    }

    // Update conductor
    cond->SetPortImpedance(portImpedance_);
    cond->SetTotalCurrent(totalCurrent_);
}

std::complex<double> radTConductorSolver::SolveImpedance(double frequency) {
    SetFrequency(frequency);
    Solve();
    return portImpedance_;
}

std::vector<std::complex<double>> radTConductorSolver::ImpedanceSweep(
    const std::vector<double>& frequencies) {

    std::vector<std::complex<double>> result(frequencies.size());

    for (size_t i = 0; i < frequencies.size(); ++i) {
        result[i] = SolveImpedance(frequencies[i]);
    }

    return result;
}

// ============================================================================
// RWG mesh building and field computation
// ============================================================================

void radTConductorSolver::BuildRWGMesh() {
    // Collect all vertices and triangles from conductor panels
    // Use vertex merging to share vertices between adjacent panels
    std::vector<TVector3d> vertices;
    std::vector<std::array<int, 3>> triangles;

    // Tolerance for vertex matching
    const double tol = 1e-10;

    // Lambda to find or add vertex
    auto findOrAddVertex = [&](const TVector3d& v) -> int {
        for (size_t i = 0; i < vertices.size(); ++i) {
            TVector3d diff = vertices[i] - v;
            if (diff.Abs() < tol) {
                return static_cast<int>(i);
            }
        }
        vertices.push_back(v);
        return static_cast<int>(vertices.size() - 1);
    };

    for (const auto& cond : conductors_) {
        const auto& panels = cond->GetPanels();

        for (const auto& panel : panels) {
            // Get vertex indices (with merging)
            std::vector<int> vIdx;
            for (const auto& v : panel.vertices) {
                vIdx.push_back(findOrAddVertex(v));
            }

            // Create triangles from panel
            if (panel.type == SurfacePanel::Triangle && vIdx.size() >= 3) {
                // Single triangle
                triangles.push_back({vIdx[0], vIdx[1], vIdx[2]});
            } else if (vIdx.size() >= 4) {
                // Quad -> 2 triangles
                triangles.push_back({vIdx[0], vIdx[1], vIdx[2]});
                triangles.push_back({vIdx[0], vIdx[2], vIdx[3]});
            }
        }
    }

    // Create RWG mesh
    rwgMesh_ = std::make_shared<RWGMesh>();
    rwgMesh_->CreateFromTriangles(vertices, triangles);
}

void radTConductorSolver::ComputeB(const TVector3d& point,
                                   std::complex<double>& Bx,
                                   std::complex<double>& By,
                                   std::complex<double>& Bz) const {
    Bx = By = Bz = std::complex<double>(0, 0);

    if (!rwgSolver_) return;

    rwgSolver_->ComputeB(point, Bx, By, Bz);
}

void radTConductorSolver::BuildSystemMatrix() {
    // Not used - BuildRWGMesh creates mesh instead
}

void radTConductorSolver::SolveLinearSystem() {
    // Not used - RWGEFIESolver handles linear solve
}

void radTConductorSolver::ExtractSolution() {
    // Not used - RWGEFIESolver extracts its own solution
}

void radTConductorSolver::ComputePortImpedance() {
    // Not used - RWGEFIESolver computes impedance
}

} // namespace radia
