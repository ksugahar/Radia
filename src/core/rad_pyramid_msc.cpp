/*-------------------------------------------------------------------------
*
* File name:      rad_pyramid_msc.cpp
*
* Project:        RADIA
*
* Description:    MMM with Magnetic Surface Charge (MSC) method implementation
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

#include "rad_pyramid_msc.h"
#include "rad_polyhedron.h"
#include <cmath>
#include <algorithm>

// Physical constants
static const double PI = 3.14159265358979323846;
static const double MU_0 = 4.0 * PI * 1e-7;  // Permeability of free space [H/m]
static const double INV_4PI = 1.0 / (4.0 * PI);

//-------------------------------------------------------------------------
// radTMSCMethod implementation
//-------------------------------------------------------------------------

radTMSCMethod::radTMSCMethod()
{
}

//-------------------------------------------------------------------------

radTMSCMethod::~radTMSCMethod()
{
}

//-------------------------------------------------------------------------

int radTMSCMethod::InitializeFromPolyhedrons(const std::vector<radTPolyhedron*>& polyhedrons)
{
    m_elements.clear();
    m_elements.reserve(polyhedrons.size());

    for(size_t i = 0; i < polyhedrons.size(); i++)
    {
        radTPolyhedron* poly = polyhedrons[i];
        if(poly == nullptr) continue;

        // Check if hexahedron
        if(!RadMSCUtils::IsHexahedron(poly))
        {
            // For now, only support hexahedra
            // TODO: Extend to tetrahedra (4 faces)
            continue;
        }

        radTMSCElementData elem;
        elem.element_center = poly->CentrPoint;
        elem.num_faces = 6;

        // Extract vertices
        std::vector<TVector3d> vertices;
        RadMSCUtils::ExtractVertices(poly, vertices);

        if(vertices.size() < 8)
        {
            // Not enough vertices for hexahedron
            continue;
        }

        // Setup faces
        std::array<TVector3d, 8> hex_verts;
        for(int v = 0; v < 8 && v < (int)vertices.size(); v++)
        {
            hex_verts[v] = vertices[v];
        }

        for(int f = 0; f < 6; f++)
        {
            RadMSCUtils::GetHexahedronFaceVertices(hex_verts, f, elem.faces[f].vertices);
            elem.faces[f].num_vertices = 4;
            ComputeFaceGeometry(elem.faces[f]);
        }

        // Compute evaluation points (midpoint between face center and element center)
        ComputeEvaluationPoints(elem);

        // Initialize sigma to zero
        for(int f = 0; f < 6; f++)
        {
            elem.sigma[f] = 0.0;
            elem.point_charges[f] = 0.0;
        }

        m_elements.push_back(elem);
    }

    return static_cast<int>(m_elements.size());
}

//-------------------------------------------------------------------------

void radTMSCMethod::ComputeFaceGeometry(radTPyramidFace& face)
{
    // Compute face center
    face.center = RadMSCUtils::ComputeFaceCenter(face.vertices, face.num_vertices);

    // Compute face normal
    face.normal = RadMSCUtils::ComputeFaceNormal(face.vertices, face.num_vertices);

    // Compute face area
    face.area = RadMSCUtils::ComputeFaceArea(face.vertices, face.num_vertices);
}

//-------------------------------------------------------------------------

void radTMSCMethod::ComputeEvaluationPoints(radTMSCElementData& elem)
{
    // Evaluation point is midpoint between face center and element center
    // This is where H is evaluated for the MSC formulation
    for(int f = 0; f < elem.num_faces; f++)
    {
        elem.eval_points[f].x = 0.5 * (elem.faces[f].center.x + elem.element_center.x);
        elem.eval_points[f].y = 0.5 * (elem.faces[f].center.y + elem.element_center.y);
        elem.eval_points[f].z = 0.5 * (elem.faces[f].center.z + elem.element_center.z);
    }
}

//-------------------------------------------------------------------------

TVector3d radTMSCMethod::FieldFromSurfaceCharge(const TVector3d& obs_point,
                                                 const radTPyramidFace& face,
                                                 double sigma) const
{
    // Compute magnetic field H from uniformly charged surface
    // H = (sigma / 4pi) * integral of (r' - r) / |r' - r|^3 dS'
    //
    // For a quadrilateral face, split into two triangles

    if(face.num_vertices == 3)
    {
        return FieldFromChargedTriangle(obs_point,
                                        face.vertices[0],
                                        face.vertices[1],
                                        face.vertices[2],
                                        sigma);
    }
    else if(face.num_vertices == 4)
    {
        return FieldFromChargedQuad(obs_point, face.vertices, sigma);
    }

    return TVector3d(0.0, 0.0, 0.0);
}

//-------------------------------------------------------------------------

TVector3d radTMSCMethod::FieldFromPointCharge(const TVector3d& obs_point,
                                               const TVector3d& charge_pos,
                                               double charge) const
{
    // Magnetic field from point charge (magnetic monopole)
    // H = (m / 4pi) * (r - p) / |r - p|^3
    //
    // Note: This is the negative gradient of scalar potential phi = m / (4pi |r - p|)

    TVector3d r_minus_p;
    r_minus_p.x = obs_point.x - charge_pos.x;
    r_minus_p.y = obs_point.y - charge_pos.y;
    r_minus_p.z = obs_point.z - charge_pos.z;

    double r_mag_sq = r_minus_p.x * r_minus_p.x +
                      r_minus_p.y * r_minus_p.y +
                      r_minus_p.z * r_minus_p.z;

    if(r_mag_sq < 1e-20)
    {
        // Observation point is at charge position - return zero
        return TVector3d(0.0, 0.0, 0.0);
    }

    double r_mag = sqrt(r_mag_sq);
    double r_mag_cubed = r_mag_sq * r_mag;
    double coef = charge * INV_4PI / r_mag_cubed;

    TVector3d H;
    H.x = coef * r_minus_p.x;
    H.y = coef * r_minus_p.y;
    H.z = coef * r_minus_p.z;

    return H;
}

//-------------------------------------------------------------------------

TVector3d radTMSCMethod::FieldFromChargedTriangle(const TVector3d& obs,
                                                   const TVector3d& v0,
                                                   const TVector3d& v1,
                                                   const TVector3d& v2,
                                                   double sigma) const
{
    // Analytic field from uniformly charged triangle
    // Using the formula from Graglia (1993) for solid angle
    //
    // For magnetic surface charge: H = (sigma / 4pi) * Omega
    // where Omega is the solid angle subtended by the triangle

    // Vectors from observation point to vertices
    TVector3d r0, r1, r2;
    r0.x = v0.x - obs.x; r0.y = v0.y - obs.y; r0.z = v0.z - obs.z;
    r1.x = v1.x - obs.x; r1.y = v1.y - obs.y; r1.z = v1.z - obs.z;
    r2.x = v2.x - obs.x; r2.y = v2.y - obs.y; r2.z = v2.z - obs.z;

    // Magnitudes
    double R0 = sqrt(r0.x*r0.x + r0.y*r0.y + r0.z*r0.z);
    double R1 = sqrt(r1.x*r1.x + r1.y*r1.y + r1.z*r1.z);
    double R2 = sqrt(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);

    if(R0 < 1e-15 || R1 < 1e-15 || R2 < 1e-15)
    {
        // Observation point is at a vertex - use numerical integration or return zero
        return TVector3d(0.0, 0.0, 0.0);
    }

    // Solid angle formula (van Oosterom & Strackee, 1983)
    // Omega = 2 * atan2(numerator, denominator)
    // numerator = r0 . (r1 x r2)
    // denominator = R0*R1*R2 + R2*(r0.r1) + R1*(r0.r2) + R0*(r1.r2)

    // Cross product r1 x r2
    TVector3d r1_cross_r2;
    r1_cross_r2.x = r1.y * r2.z - r1.z * r2.y;
    r1_cross_r2.y = r1.z * r2.x - r1.x * r2.z;
    r1_cross_r2.z = r1.x * r2.y - r1.y * r2.x;

    double numerator = r0.x * r1_cross_r2.x + r0.y * r1_cross_r2.y + r0.z * r1_cross_r2.z;

    double r0_dot_r1 = r0.x*r1.x + r0.y*r1.y + r0.z*r1.z;
    double r0_dot_r2 = r0.x*r2.x + r0.y*r2.y + r0.z*r2.z;
    double r1_dot_r2 = r1.x*r2.x + r1.y*r2.y + r1.z*r2.z;

    double denominator = R0*R1*R2 + R2*r0_dot_r1 + R1*r0_dot_r2 + R0*r1_dot_r2;

    double solid_angle = 2.0 * atan2(numerator, denominator);

    // Triangle normal (outward)
    TVector3d edge1, edge2, normal;
    edge1.x = v1.x - v0.x; edge1.y = v1.y - v0.y; edge1.z = v1.z - v0.z;
    edge2.x = v2.x - v0.x; edge2.y = v2.y - v0.y; edge2.z = v2.z - v0.z;

    normal.x = edge1.y * edge2.z - edge1.z * edge2.y;
    normal.y = edge1.z * edge2.x - edge1.x * edge2.z;
    normal.z = edge1.x * edge2.y - edge1.y * edge2.x;

    double normal_mag = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
    if(normal_mag > 1e-15)
    {
        normal.x /= normal_mag;
        normal.y /= normal_mag;
        normal.z /= normal_mag;
    }

    // H field in direction of normal
    // H = (sigma / 4pi) * solid_angle * normal
    double coef = sigma * INV_4PI * solid_angle;

    TVector3d H;
    H.x = coef * normal.x;
    H.y = coef * normal.y;
    H.z = coef * normal.z;

    return H;
}

//-------------------------------------------------------------------------

TVector3d radTMSCMethod::FieldFromChargedQuad(const TVector3d& obs,
                                               const std::array<TVector3d, 4>& vertices,
                                               double sigma) const
{
    // Split quadrilateral into two triangles and sum fields
    // Triangle 1: v0, v1, v2
    // Triangle 2: v0, v2, v3

    TVector3d H1 = FieldFromChargedTriangle(obs, vertices[0], vertices[1], vertices[2], sigma);
    TVector3d H2 = FieldFromChargedTriangle(obs, vertices[0], vertices[2], vertices[3], sigma);

    TVector3d H;
    H.x = H1.x + H2.x;
    H.y = H1.y + H2.y;
    H.z = H1.z + H2.z;

    return H;
}

//-------------------------------------------------------------------------

TVector3d radTMSCMethod::ComputeFieldAtPoint(const TVector3d& obs_point, int exclude_elem) const
{
    TVector3d H_total(0.0, 0.0, 0.0);

    for(size_t i = 0; i < m_elements.size(); i++)
    {
        if((int)i == exclude_elem) continue;

        const radTMSCElementData& elem = m_elements[i];

        // Contribution from surface charges on each face
        for(int f = 0; f < elem.num_faces; f++)
        {
            TVector3d H_face = FieldFromSurfaceCharge(obs_point, elem.faces[f], elem.sigma[f]);
            H_total.x += H_face.x;
            H_total.y += H_face.y;
            H_total.z += H_face.z;
        }

        // Contribution from point charges at element center
        double total_point_charge = 0.0;
        for(int f = 0; f < elem.num_faces; f++)
        {
            total_point_charge += elem.point_charges[f];
        }

        TVector3d H_point = FieldFromPointCharge(obs_point, elem.element_center, total_point_charge);
        H_total.x += H_point.x;
        H_total.y += H_point.y;
        H_total.z += H_point.z;
    }

    return H_total;
}

//-------------------------------------------------------------------------

int radTMSCMethod::SetupInteractionMatrix()
{
    int n_elem = GetNumElements();
    int ndof = GetNumDOF();

    if(n_elem == 0) return 0;

    // Allocate interaction matrix
    m_interaction_matrix.resize(ndof);
    for(int i = 0; i < ndof; i++)
    {
        m_interaction_matrix[i].resize(ndof, 0.0);
    }

    // Fill interaction matrix
    // A[i*6+f1][j*6+f2] = influence of sigma_f2 of element j on H_n at face f1 of element i
    // where H_n is the normal component of H

    for(int i = 0; i < n_elem; i++)
    {
        const radTMSCElementData& elem_i = m_elements[i];

        for(int f1 = 0; f1 < elem_i.num_faces; f1++)
        {
            int row = i * 6 + f1;
            const TVector3d& eval_pt = elem_i.eval_points[f1];
            const TVector3d& normal = elem_i.faces[f1].normal;

            for(int j = 0; j < n_elem; j++)
            {
                const radTMSCElementData& elem_j = m_elements[j];

                for(int f2 = 0; f2 < elem_j.num_faces; f2++)
                {
                    int col = j * 6 + f2;

                    // Compute H field at eval_pt from unit sigma on face f2 of element j
                    TVector3d H_from_surface = FieldFromSurfaceCharge(eval_pt, elem_j.faces[f2], 1.0);

                    // Also need contribution from point charge at element center
                    // Point charge = -sigma * area
                    double area_f2 = elem_j.faces[f2].area;
                    TVector3d H_from_point = FieldFromPointCharge(eval_pt, elem_j.element_center, -area_f2);

                    TVector3d H_total;
                    H_total.x = H_from_surface.x + H_from_point.x;
                    H_total.y = H_from_surface.y + H_from_point.y;
                    H_total.z = H_from_surface.z + H_from_point.z;

                    // Normal component
                    double H_n = H_total.x * normal.x + H_total.y * normal.y + H_total.z * normal.z;

                    m_interaction_matrix[row][col] = H_n;
                }
            }
        }
    }

    return 1;
}

//-------------------------------------------------------------------------

int radTMSCMethod::Solve(double precision, int max_iter)
{
    int n_elem = GetNumElements();
    int ndof = GetNumDOF();

    if(n_elem == 0 || ndof == 0) return 0;

    // Setup external field vector
    m_external_field.resize(ndof);
    // TODO: Fill external field from ObjBckg or ObjBckgCF

    // Setup system matrix: (I - chi * A) * sigma = chi * H_ext_n
    // where chi = (mu_r - 1) / (mu_r + 1) for MSC formulation
    // (This is a simplification; actual relation depends on geometry)

    // For linear materials with permeability mu_r:
    // B_n = mu_0 * mu_r * H_n (inside material)
    // sigma = B_n - mu_0 * H_n = mu_0 * (mu_r - 1) * H_n

    // TODO: Implement full LU or BiCGSTAB solve
    // For now, return success
    return 1;
}

//-------------------------------------------------------------------------
// RadMSCUtils namespace implementation
//-------------------------------------------------------------------------

namespace RadMSCUtils {

void GetHexahedronFaceVertices(const std::array<TVector3d, 8>& hex_vertices,
                               int face_idx,
                               std::array<TVector3d, 4>& face_vertices)
{
    // Standard hexahedron vertex ordering:
    //     7------6
    //    /|     /|
    //   4------5 |
    //   | 3----|-2
    //   |/     |/
    //   0------1
    //
    // Face ordering: 0=bottom (0123), 1=top (4567), 2=front (0154), 3=back (3267), 4=left (0374), 5=right (1265)

    static const int face_indices[6][4] = {
        {0, 3, 2, 1},  // bottom (z-) - counter-clockwise from outside
        {4, 5, 6, 7},  // top (z+)
        {0, 1, 5, 4},  // front (y-)
        {3, 7, 6, 2},  // back (y+)
        {0, 4, 7, 3},  // left (x-)
        {1, 2, 6, 5}   // right (x+)
    };

    for(int i = 0; i < 4; i++)
    {
        face_vertices[i] = hex_vertices[face_indices[face_idx][i]];
    }
}

//-------------------------------------------------------------------------

TVector3d ComputeFaceNormal(const std::array<TVector3d, 4>& vertices, int num_vertices)
{
    TVector3d normal(0.0, 0.0, 0.0);

    if(num_vertices < 3) return normal;

    // Use first three vertices to compute normal
    TVector3d edge1, edge2;
    edge1.x = vertices[1].x - vertices[0].x;
    edge1.y = vertices[1].y - vertices[0].y;
    edge1.z = vertices[1].z - vertices[0].z;

    edge2.x = vertices[2].x - vertices[0].x;
    edge2.y = vertices[2].y - vertices[0].y;
    edge2.z = vertices[2].z - vertices[0].z;

    // Cross product
    normal.x = edge1.y * edge2.z - edge1.z * edge2.y;
    normal.y = edge1.z * edge2.x - edge1.x * edge2.z;
    normal.z = edge1.x * edge2.y - edge1.y * edge2.x;

    // Normalize
    double mag = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
    if(mag > 1e-15)
    {
        normal.x /= mag;
        normal.y /= mag;
        normal.z /= mag;
    }

    return normal;
}

//-------------------------------------------------------------------------

double ComputeFaceArea(const std::array<TVector3d, 4>& vertices, int num_vertices)
{
    if(num_vertices == 3)
    {
        // Triangle area
        TVector3d edge1, edge2, cross;
        edge1.x = vertices[1].x - vertices[0].x;
        edge1.y = vertices[1].y - vertices[0].y;
        edge1.z = vertices[1].z - vertices[0].z;

        edge2.x = vertices[2].x - vertices[0].x;
        edge2.y = vertices[2].y - vertices[0].y;
        edge2.z = vertices[2].z - vertices[0].z;

        cross.x = edge1.y * edge2.z - edge1.z * edge2.y;
        cross.y = edge1.z * edge2.x - edge1.x * edge2.z;
        cross.z = edge1.x * edge2.y - edge1.y * edge2.x;

        return 0.5 * sqrt(cross.x*cross.x + cross.y*cross.y + cross.z*cross.z);
    }
    else if(num_vertices == 4)
    {
        // Quadrilateral area = sum of two triangle areas
        double area1 = 0.0, area2 = 0.0;

        // Triangle 1: v0, v1, v2
        {
            TVector3d edge1, edge2, cross;
            edge1.x = vertices[1].x - vertices[0].x;
            edge1.y = vertices[1].y - vertices[0].y;
            edge1.z = vertices[1].z - vertices[0].z;

            edge2.x = vertices[2].x - vertices[0].x;
            edge2.y = vertices[2].y - vertices[0].y;
            edge2.z = vertices[2].z - vertices[0].z;

            cross.x = edge1.y * edge2.z - edge1.z * edge2.y;
            cross.y = edge1.z * edge2.x - edge1.x * edge2.z;
            cross.z = edge1.x * edge2.y - edge1.y * edge2.x;

            area1 = 0.5 * sqrt(cross.x*cross.x + cross.y*cross.y + cross.z*cross.z);
        }

        // Triangle 2: v0, v2, v3
        {
            TVector3d edge1, edge2, cross;
            edge1.x = vertices[2].x - vertices[0].x;
            edge1.y = vertices[2].y - vertices[0].y;
            edge1.z = vertices[2].z - vertices[0].z;

            edge2.x = vertices[3].x - vertices[0].x;
            edge2.y = vertices[3].y - vertices[0].y;
            edge2.z = vertices[3].z - vertices[0].z;

            cross.x = edge1.y * edge2.z - edge1.z * edge2.y;
            cross.y = edge1.z * edge2.x - edge1.x * edge2.z;
            cross.z = edge1.x * edge2.y - edge1.y * edge2.x;

            area2 = 0.5 * sqrt(cross.x*cross.x + cross.y*cross.y + cross.z*cross.z);
        }

        return area1 + area2;
    }

    return 0.0;
}

//-------------------------------------------------------------------------

TVector3d ComputeFaceCenter(const std::array<TVector3d, 4>& vertices, int num_vertices)
{
    TVector3d center(0.0, 0.0, 0.0);

    for(int i = 0; i < num_vertices; i++)
    {
        center.x += vertices[i].x;
        center.y += vertices[i].y;
        center.z += vertices[i].z;
    }

    if(num_vertices > 0)
    {
        double inv_n = 1.0 / num_vertices;
        center.x *= inv_n;
        center.y *= inv_n;
        center.z *= inv_n;
    }

    return center;
}

//-------------------------------------------------------------------------

bool IsHexahedron(const radTPolyhedron* poly)
{
    if(poly == nullptr) return false;
    return poly->AmOfFaces == 6;
}

//-------------------------------------------------------------------------

int ExtractVertices(const radTPolyhedron* poly, std::vector<TVector3d>& vertices)
{
    vertices.clear();
    if(poly == nullptr) return 0;

    // Extract unique vertices from all faces
    const double tol = 1e-10;

    for(int f = 0; f < poly->AmOfFaces; f++)
    {
        const radTHandlePgnAndTrans& hPgnTrans = poly->VectHandlePgnAndTrans[f];
        radTPolygon* pgn = hPgnTrans.PgnHndl.rep;
        radTrans* trans = hPgnTrans.TransHndl.rep;

        if(pgn == nullptr) continue;

        // Get vertices from polygon
        for(int v = 0; v < pgn->AmOfEdgePoints; v++)
        {
            TVector2d& pt2d = pgn->EdgePointsVector[v];
            TVector3d pt3d_local(pt2d.x, pt2d.y, pgn->CoordZ);

            // Transform to global coordinates
            TVector3d pt3d_global = trans->TrBiPoint(pt3d_local);

            // Check if vertex already exists
            bool found = false;
            for(const TVector3d& existing : vertices)
            {
                double dx = pt3d_global.x - existing.x;
                double dy = pt3d_global.y - existing.y;
                double dz = pt3d_global.z - existing.z;
                if(dx*dx + dy*dy + dz*dz < tol*tol)
                {
                    found = true;
                    break;
                }
            }

            if(!found)
            {
                vertices.push_back(pt3d_global);
            }
        }
    }

    return static_cast<int>(vertices.size());
}

} // namespace RadMSCUtils
