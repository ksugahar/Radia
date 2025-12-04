/*-------------------------------------------------------------------------
*
* File name:      rad_pyramid_msc.h
*
* Project:        RADIA
*
* Description:    MMM with Magnetic Surface Charge (MSC) method
*                 Reference: H. Yano, K. Sugahara, "Magnetic Moment Method
*                 with the Idea of Magnetic Surface Charge Method"
*                 J. Magn. Soc. Jpn., 2023
*
* Author(s):      ELF Corporation / Claude Code
*
* First release:  2025
*
* Copyright (C):  2025
*
-------------------------------------------------------------------------*/

#ifndef __RAD_PYRAMID_MSC_H
#define __RAD_PYRAMID_MSC_H

#include "rad_geometry_3d.h"
#include "rad_planar_2d.h"
#include <vector>
#include <array>

//-------------------------------------------------------------------------
// Forward declarations
//-------------------------------------------------------------------------

class radTPolyhedron;

//-------------------------------------------------------------------------
// Pyramid structure for MSC method
// Each hexahedron is subdivided into 6 pyramids
//-------------------------------------------------------------------------

struct radTPyramidFace {
    TVector3d center;           // Face center point
    TVector3d normal;           // Outward normal vector
    double area;                // Face area
    std::array<TVector3d, 4> vertices;  // Face vertices (up to 4 for quad faces)
    int num_vertices;           // Number of vertices (3 for triangle, 4 for quad)

    radTPyramidFace() : area(0.0), num_vertices(0) {
        center.x = center.y = center.z = 0.0;
        normal.x = normal.y = normal.z = 0.0;
    }
};

//-------------------------------------------------------------------------
// MSC element data structure
// Stores pyramid decomposition and surface charge data
//-------------------------------------------------------------------------

struct radTMSCElementData {
    // Geometry
    TVector3d element_center;       // Element center (common apex for all pyramids)
    std::array<radTPyramidFace, 6> faces;  // 6 faces of hexahedron
    int num_faces;                  // Number of faces (6 for hex, 4 for tet)

    // Field evaluation points
    // H is evaluated at midpoint between face center and element center
    std::array<TVector3d, 6> eval_points;

    // Surface charge densities (unknowns)
    std::array<double, 6> sigma;    // sigma_n for each face

    // Point charges at element center (derived from sigma)
    // m_n = -sigma_n * area_n (negative to satisfy Gauss's law)
    std::array<double, 6> point_charges;

    radTMSCElementData() : num_faces(0) {
        element_center.x = element_center.y = element_center.z = 0.0;
        for(int i = 0; i < 6; i++) {
            sigma[i] = 0.0;
            point_charges[i] = 0.0;
        }
    }

    // Compute point charges from sigma values
    void UpdatePointCharges() {
        for(int i = 0; i < num_faces; i++) {
            point_charges[i] = -sigma[i] * faces[i].area;
        }
    }

    // Total point charge at element center (should be zero for Gauss's law)
    double TotalPointCharge() const {
        double total = 0.0;
        for(int i = 0; i < num_faces; i++) {
            total += point_charges[i];
        }
        return total;
    }
};

//-------------------------------------------------------------------------
// MSC Method class
// Handles the MMM with MSC formulation
//-------------------------------------------------------------------------

class radTMSCMethod {

public:
    radTMSCMethod();
    ~radTMSCMethod();

    // Initialize MSC data from polyhedron elements
    int InitializeFromPolyhedrons(const std::vector<radTPolyhedron*>& polyhedrons);

    // Setup interaction matrix for MSC formulation
    // Matrix size: 6N x 6N (6 DOF per element)
    int SetupInteractionMatrix();

    // Compute magnetic field H at point r from surface charge sigma on a face
    // H = (1/4pi) * integral of sigma * (r' - r) / |r' - r|^3 dS'
    TVector3d FieldFromSurfaceCharge(const TVector3d& obs_point,
                                      const radTPyramidFace& face,
                                      double sigma) const;

    // Compute magnetic field H at point r from point charge m at position p
    // H = (1/4pi) * m * (r - p) / |r - p|^3
    TVector3d FieldFromPointCharge(const TVector3d& obs_point,
                                    const TVector3d& charge_pos,
                                    double charge) const;

    // Compute total field at evaluation point from all sources
    TVector3d ComputeFieldAtPoint(const TVector3d& obs_point, int exclude_elem = -1) const;

    // Solve the MSC system
    // B_n = mu_0 * (H_ext + H_ind)_n = mu_0 * mu_r * H_n
    // This gives: sigma_n as function of H_ext
    int Solve(double precision, int max_iter);

    // Get number of elements
    int GetNumElements() const { return static_cast<int>(m_elements.size()); }

    // Get number of DOF (6 per element)
    int GetNumDOF() const { return GetNumElements() * 6; }

    // Access element data
    const radTMSCElementData& GetElement(int idx) const { return m_elements[idx]; }
    radTMSCElementData& GetElement(int idx) { return m_elements[idx]; }

private:
    std::vector<radTMSCElementData> m_elements;

    // Interaction matrix: A[i*6+f1][j*6+f2] = influence of sigma_f2 of element j on H at face f1 of element i
    std::vector<std::vector<double>> m_interaction_matrix;

    // External field at each evaluation point
    std::vector<double> m_external_field;  // Size: 6N (scalar H dot n for each face)

    // Relative permeability for each element
    std::vector<double> m_mu_r;

    // Helper functions
    void ComputeFaceGeometry(radTPyramidFace& face);
    void ComputeEvaluationPoints(radTMSCElementData& elem);

    // Analytic field from uniformly charged triangle
    TVector3d FieldFromChargedTriangle(const TVector3d& obs,
                                        const TVector3d& v0,
                                        const TVector3d& v1,
                                        const TVector3d& v2,
                                        double sigma) const;

    // Analytic field from uniformly charged quadrilateral (split into 2 triangles)
    TVector3d FieldFromChargedQuad(const TVector3d& obs,
                                    const std::array<TVector3d, 4>& vertices,
                                    double sigma) const;
};

//-------------------------------------------------------------------------
// Utility functions for hexahedron to pyramid decomposition
//-------------------------------------------------------------------------

namespace RadMSCUtils {

    // Get face vertices of hexahedron in standard ordering
    // face_idx: 0=bottom, 1=top, 2=front, 3=back, 4=left, 5=right
    void GetHexahedronFaceVertices(const std::array<TVector3d, 8>& hex_vertices,
                                   int face_idx,
                                   std::array<TVector3d, 4>& face_vertices);

    // Compute outward normal for a face
    TVector3d ComputeFaceNormal(const std::array<TVector3d, 4>& vertices, int num_vertices);

    // Compute face area
    double ComputeFaceArea(const std::array<TVector3d, 4>& vertices, int num_vertices);

    // Compute face center
    TVector3d ComputeFaceCenter(const std::array<TVector3d, 4>& vertices, int num_vertices);

    // Check if polyhedron is hexahedron (6 faces)
    bool IsHexahedron(const radTPolyhedron* poly);

    // Extract vertices from polyhedron
    int ExtractVertices(const radTPolyhedron* poly, std::vector<TVector3d>& vertices);
}

//-------------------------------------------------------------------------

#endif // __RAD_PYRAMID_MSC_H
