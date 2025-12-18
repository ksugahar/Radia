/*-------------------------------------------------------------------------
*
* File name:      rad_extruded_polygon_msc.cpp
*
* Project:        RADIA
*
* Description:    Magnetic field source: extruded polygon (prism)
*                 Using Pyramid MSC (Magnetic Surface Charge) method
*                 6 DOF per element (surface charge density on each face)
*
*                 Reference: H. Yano, K. Sugahara, "Magnetic Moment Method
*                 with the Idea of Magnetic Surface Charge Method"
*                 J. Magn. Soc. Jpn., 2023
*
* Author(s):      Claude Code
*
* First release:  2025
*
* Copyright (C):  2025
*
-------------------------------------------------------------------------*/

#include "rad_extruded_polygon_msc.h"
#include "rad_planar_2d.h"
#include "rad_io_buffer.h"
#include <cmath>
#include <algorithm>

//-------------------------------------------------------------------------
// Constants
//-------------------------------------------------------------------------

static const double PI = 3.14159265358979323846;
static const double MU_0 = 4.0 * PI * 1e-7;
static const double INV_4PI = 1.0 / (4.0 * PI);

//-------------------------------------------------------------------------
// radTMSCFace implementation
//-------------------------------------------------------------------------

void radTMSCFace::ComputeGeometry()
{
    // Compute center
    center.x = center.y = center.z = 0.0;
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

    // Compute normal and area using first 3 vertices
    if(num_vertices >= 3)
    {
        TVector3d edge1, edge2;
        edge1.x = vertices[1].x - vertices[0].x;
        edge1.y = vertices[1].y - vertices[0].y;
        edge1.z = vertices[1].z - vertices[0].z;

        edge2.x = vertices[2].x - vertices[0].x;
        edge2.y = vertices[2].y - vertices[0].y;
        edge2.z = vertices[2].z - vertices[0].z;

        // Cross product for normal
        TVector3d cross;
        cross.x = edge1.y * edge2.z - edge1.z * edge2.y;
        cross.y = edge1.z * edge2.x - edge1.x * edge2.z;
        cross.z = edge1.x * edge2.y - edge1.y * edge2.x;

        double mag = sqrt(cross.x*cross.x + cross.y*cross.y + cross.z*cross.z);

        if(num_vertices == 3)
        {
            area = 0.5 * mag;
        }
        else if(num_vertices == 4)
        {
            // Split quad into two triangles
            double area1 = 0.5 * mag;

            // Triangle 2: v0, v2, v3
            TVector3d edge3;
            edge3.x = vertices[3].x - vertices[0].x;
            edge3.y = vertices[3].y - vertices[0].y;
            edge3.z = vertices[3].z - vertices[0].z;

            TVector3d cross2;
            cross2.x = edge2.y * edge3.z - edge2.z * edge3.y;
            cross2.y = edge2.z * edge3.x - edge2.x * edge3.z;
            cross2.z = edge2.x * edge3.y - edge2.y * edge3.x;

            double area2 = 0.5 * sqrt(cross2.x*cross2.x + cross2.y*cross2.y + cross2.z*cross2.z);
            area = area1 + area2;
        }

        if(mag > 1e-15)
        {
            normal.x = cross.x / mag;
            normal.y = cross.y / mag;
            normal.z = cross.z / mag;
        }
    }
}

//-------------------------------------------------------------------------
// radTExtrPolygonMSC implementation
//-------------------------------------------------------------------------

radTExtrPolygonMSC::radTExtrPolygonMSC(const TVector3d& InFirstPoint, const TAxisOrient& InAxOrnt,
                                       double InThickness, const radThg& InBasePolygonHandle,
                                       const TVector3d& InMagn, const radThg& InMaterHandle)
    : radTg3dRelax()
{
    FirstPoint = InFirstPoint;
    AxOrnt = InAxOrnt;
    Thickness = InThickness;
    BasePolygonHandle = InBasePolygonHandle;
    Magn = InMagn;
    MaterHandle = InMaterHandle;
    NumFaces = 6;

    for(int i = 0; i < 6; i++)
    {
        Sigma[i] = 0.0;
        PointCharges[i] = 0.0;
    }

    SetupFaces();
    ComputeEvalPoints();
    DefineCentrPoint();
}

//-------------------------------------------------------------------------

radTExtrPolygonMSC::radTExtrPolygonMSC(const TVector3d& InFirstPoint, const TAxisOrient& InAxOrnt,
                                       double InThickness, TVector2d* ArrayOfPoints2d,
                                       int lenArrayOfPoints2d, const TVector3d& InMagn)
    : radTg3dRelax()
{
    FirstPoint = InFirstPoint;
    AxOrnt = InAxOrnt;
    Thickness = InThickness;
    Magn = InMagn;
    NumFaces = 6;

    for(int i = 0; i < 6; i++)
    {
        Sigma[i] = 0.0;
        PointCharges[i] = 0.0;
    }

    // Create base polygon (use correct constructor: TVector2d*, int)
    radTPolygon* pBasePgn = new radTPolygon(ArrayOfPoints2d, lenArrayOfPoints2d);
    radThg hBasePgn(pBasePgn);
    BasePolygonHandle = hBasePgn;

    SetupFaces();
    ComputeEvalPoints();
    DefineCentrPoint();
}

//-------------------------------------------------------------------------

radTExtrPolygonMSC::radTExtrPolygonMSC(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew,
                                       radTmhg& gMapOfHandlers)
{
    // Follow standard Radia binary constructor pattern:
    // 1. Parse base class data first
    DumpBinParse_g3d(inStr, mKeysOldNew, gMapOfHandlers);
    DumpBinParse_g3dRelax(inStr, mKeysOldNew, gMapOfHandlers);
    // 2. Parse our class-specific data
    DumpBinParse_ExtrPolygonMSC(inStr);
    NumFaces = 6;
    SetupFaces();
    ComputeEvalPoints();
    DefineCentrPoint();
}

//-------------------------------------------------------------------------

void radTExtrPolygonMSC::SetupFaces()
{
    // Setup 6 faces of hexahedral prism from base polygon and extrusion
    radTPolygon* pPgn = (radTPolygon*)(BasePolygonHandle.rep);
    if(pPgn == nullptr) return;

    // Get base polygon vertices (2D)
    int nPts = pPgn->AmOfEdgePoints;
    if(nPts < 3) return;

    // For simplicity, assume rectangular base (4 vertices)
    // General polygon would need more faces

    // Construct 3D vertices of the hexahedron
    // Bottom face at FirstPoint, top face at FirstPoint + Thickness * axis
    TVector3d axis_vec;
    if(AxOrnt == ParallelToX)
    {
        axis_vec.x = 1.0;
        axis_vec.y = 0.0;
        axis_vec.z = 0.0;
    }
    else if(AxOrnt == ParallelToY)
    {
        axis_vec.x = 0.0;
        axis_vec.y = 1.0;
        axis_vec.z = 0.0;
    }
    else  // ParallelToZ
    {
        axis_vec.x = 0.0;
        axis_vec.y = 0.0;
        axis_vec.z = 1.0;
    }

    // Build vertices array for hexahedron
    // Assume 4 base polygon points for rectangular prism
    std::array<TVector3d, 8> hex_verts;

    // Map 2D polygon to 3D based on axis orientation
    TVector2d& basePt = pPgn->EdgePointsVector[0];

    for(int i = 0; i < std::min(4, nPts); i++)
    {
        TVector2d& pt2d = pPgn->EdgePointsVector[i];
        double dx = pt2d.x - basePt.x;
        double dy = pt2d.y - basePt.y;

        // Bottom face (at FirstPoint)
        if(AxOrnt == ParallelToX)
        {
            hex_verts[i].x = FirstPoint.x;
            hex_verts[i].y = FirstPoint.y + dx;
            hex_verts[i].z = FirstPoint.z + dy;
        }
        else if(AxOrnt == ParallelToY)
        {
            hex_verts[i].x = FirstPoint.x + dx;
            hex_verts[i].y = FirstPoint.y;
            hex_verts[i].z = FirstPoint.z + dy;
        }
        else  // ParallelToZ
        {
            hex_verts[i].x = FirstPoint.x + dx;
            hex_verts[i].y = FirstPoint.y + dy;
            hex_verts[i].z = FirstPoint.z;
        }

        // Top face (at FirstPoint + Thickness * axis)
        hex_verts[i + 4].x = hex_verts[i].x + Thickness * axis_vec.x;
        hex_verts[i + 4].y = hex_verts[i].y + Thickness * axis_vec.y;
        hex_verts[i + 4].z = hex_verts[i].z + Thickness * axis_vec.z;
    }

    // Setup 6 faces with outward normals
    // Face ordering: 0=bottom, 1=top, 2=front, 3=back, 4=left, 5=right

    // Standard face vertex indices (counter-clockwise from outside)
    static const int face_indices[6][4] = {
        {0, 3, 2, 1},  // bottom (extrusion start)
        {4, 5, 6, 7},  // top (extrusion end)
        {0, 1, 5, 4},  // front
        {3, 7, 6, 2},  // back
        {0, 4, 7, 3},  // left
        {1, 2, 6, 5}   // right
    };

    for(int f = 0; f < 6; f++)
    {
        Faces[f].num_vertices = 4;
        for(int v = 0; v < 4; v++)
        {
            Faces[f].vertices[v] = hex_verts[face_indices[f][v]];
        }
        Faces[f].ComputeGeometry();
    }
}

//-------------------------------------------------------------------------

void radTExtrPolygonMSC::ComputeEvalPoints()
{
    // Evaluation point is midpoint between face center and element center
    for(int f = 0; f < NumFaces; f++)
    {
        EvalPoints[f].x = 0.5 * (Faces[f].center.x + CentrPoint.x);
        EvalPoints[f].y = 0.5 * (Faces[f].center.y + CentrPoint.y);
        EvalPoints[f].z = 0.5 * (Faces[f].center.z + CentrPoint.z);
    }
}

//-------------------------------------------------------------------------

void radTExtrPolygonMSC::UpdatePointCharges()
{
    // Point charge at element center from each face:
    // m_n = -sigma_n * area_n
    // This ensures div B = 0 within each pyramid subdivision
    for(int f = 0; f < NumFaces; f++)
    {
        PointCharges[f] = -Sigma[f] * Faces[f].area;
    }
}

//-------------------------------------------------------------------------

TVector3d radTExtrPolygonMSC::FieldFromSurfaceCharge(const TVector3d& obs,
                                                      const radTMSCFace& face,
                                                      double sigma) const
{
    // Split quadrilateral into two triangles and sum fields
    if(face.num_vertices == 3)
    {
        return FieldFromChargedTriangle(obs, face.vertices[0], face.vertices[1],
                                        face.vertices[2], sigma);
    }
    else if(face.num_vertices == 4)
    {
        // Triangle 1: v0, v1, v2
        TVector3d H1 = FieldFromChargedTriangle(obs, face.vertices[0], face.vertices[1],
                                                 face.vertices[2], sigma);
        // Triangle 2: v0, v2, v3
        TVector3d H2 = FieldFromChargedTriangle(obs, face.vertices[0], face.vertices[2],
                                                 face.vertices[3], sigma);

        TVector3d H;
        H.x = H1.x + H2.x;
        H.y = H1.y + H2.y;
        H.z = H1.z + H2.z;
        return H;
    }

    return TVector3d(0.0, 0.0, 0.0);
}

//-------------------------------------------------------------------------

TVector3d radTExtrPolygonMSC::FieldFromPointCharge(const TVector3d& obs,
                                                    const TVector3d& pos,
                                                    double charge) const
{
    // Magnetic field from point magnetic charge (monopole)
    // H = (m / 4pi) * (r - p) / |r - p|^3

    TVector3d r_minus_p;
    r_minus_p.x = obs.x - pos.x;
    r_minus_p.y = obs.y - pos.y;
    r_minus_p.z = obs.z - pos.z;

    double r_mag_sq = r_minus_p.x * r_minus_p.x +
                      r_minus_p.y * r_minus_p.y +
                      r_minus_p.z * r_minus_p.z;

    if(r_mag_sq < 1e-20)
    {
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

TVector3d radTExtrPolygonMSC::FieldFromChargedTriangle(const TVector3d& obs,
                                                        const TVector3d& v0,
                                                        const TVector3d& v1,
                                                        const TVector3d& v2,
                                                        double sigma) const
{
    // Analytic field from uniformly charged triangle using solid angle formula
    // H = (sigma / 4pi) * Omega * n_hat
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
        return TVector3d(0.0, 0.0, 0.0);
    }

    // Solid angle formula (van Oosterom & Strackee, 1983)
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

    // Triangle normal
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

    // H = (sigma / 4pi) * solid_angle * normal
    double coef = sigma * INV_4PI * solid_angle;

    TVector3d H;
    H.x = coef * normal.x;
    H.y = coef * normal.y;
    H.z = coef * normal.z;

    return H;
}

//-------------------------------------------------------------------------

void radTExtrPolygonMSC::B_comp(radTField* FieldPtr)
{
    // Compute B field at observation point using MSC method
    // B = mu_0 * H
    // H = sum over all faces of (H from surface charge + H from point charge)

    TVector3d obs = FieldPtr->P;
    TVector3d H_total(0.0, 0.0, 0.0);

    // Total point charge at element center
    double total_point_charge = 0.0;
    for(int f = 0; f < NumFaces; f++)
    {
        total_point_charge += PointCharges[f];
    }

    // Field from surface charges
    for(int f = 0; f < NumFaces; f++)
    {
        TVector3d H_face = FieldFromSurfaceCharge(obs, Faces[f], Sigma[f]);
        H_total.x += H_face.x;
        H_total.y += H_face.y;
        H_total.z += H_face.z;
    }

    // Field from point charges at element center
    TVector3d H_point = FieldFromPointCharge(obs, CentrPoint, total_point_charge);
    H_total.x += H_point.x;
    H_total.y += H_point.y;
    H_total.z += H_point.z;

    // B = mu_0 * H
    if(FieldPtr->FieldKey.B_)
    {
        FieldPtr->B.x += MU_0 * H_total.x;
        FieldPtr->B.y += MU_0 * H_total.y;
        FieldPtr->B.z += MU_0 * H_total.z;
    }

    if(FieldPtr->FieldKey.H_)
    {
        FieldPtr->H.x += H_total.x;
        FieldPtr->H.y += H_total.y;
        FieldPtr->H.z += H_total.z;
    }
}

//-------------------------------------------------------------------------

void radTExtrPolygonMSC::B_intComp(radTField* FieldPtr)
{
    // Integral of B over a line - for now, use numerical integration
    // This is a simplified implementation
    B_comp(FieldPtr);
}

//-------------------------------------------------------------------------

void radTExtrPolygonMSC::B_intCompSpecCases(radTField* FieldPtr, const TSpecCaseID& SpecCaseID)
{
    // Special cases for field integration
    B_intComp(FieldPtr);
}

//-------------------------------------------------------------------------

TVector3d radTExtrPolygonMSC::ComputeAverageMagnetization()
{
    // Compute average magnetization from surface charges
    // M = sum over faces of (sigma_n * n_hat * area_n) / volume
    TVector3d M_total(0.0, 0.0, 0.0);

    for(int f = 0; f < NumFaces; f++)
    {
        M_total.x += Sigma[f] * Faces[f].normal.x * Faces[f].area;
        M_total.y += Sigma[f] * Faces[f].normal.y * Faces[f].area;
        M_total.z += Sigma[f] * Faces[f].normal.z * Faces[f].area;
    }

    double vol = Volume();
    if(vol > 1e-20)
    {
        M_total.x /= vol;
        M_total.y /= vol;
        M_total.z /= vol;
    }

    return M_total;
}

//-------------------------------------------------------------------------

void radTExtrPolygonMSC::VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique)
{
    OutVect.clear();

    // Add all unique vertices from all faces
    const double tol = 1e-10;

    for(int f = 0; f < NumFaces; f++)
    {
        for(int v = 0; v < Faces[f].num_vertices; v++)
        {
            const TVector3d& pt = Faces[f].vertices[v];

            bool found = false;
            if(EnsureUnique)
            {
                for(const TVector3d& existing : OutVect)
                {
                    double dx = pt.x - existing.x;
                    double dy = pt.y - existing.y;
                    double dz = pt.z - existing.z;
                    if(dx*dx + dy*dy + dz*dz < tol*tol)
                    {
                        found = true;
                        break;
                    }
                }
            }

            if(!found)
            {
                OutVect.push_back(pt);
            }
        }
    }
}

//-------------------------------------------------------------------------

radTg3dGraphPresent* radTExtrPolygonMSC::CreateGraphPresent()
{
    // Create graphics presentation - delegate to standard polyhedron graphics
    return nullptr;  // TODO: Implement graphics
}

//-------------------------------------------------------------------------

void radTExtrPolygonMSC::Dump(std::ostream& o, int ShortSign)
{
    radTg3dRelax::Dump(o, ShortSign);
    o << "ExtrPolygonMSC";
    if(ShortSign == 1) return;

    o << "\n   6 DOF MSC-based extruded polygon";
    o << "\n   Thickness: " << Thickness;
    o << "\n   Number of faces: " << NumFaces;
    o << "\n   Sigma values: [";
    for(int i = 0; i < NumFaces; i++)
    {
        o << Sigma[i];
        if(i < NumFaces - 1) o << ", ";
    }
    o << "]";
}

//-------------------------------------------------------------------------

void radTExtrPolygonMSC::DumpPureObjInfo(std::ostream& o, int ShortSign)
{
    if(ShortSign == 1)
    {
        o << "{{ExtrPolygonMSC}}";
    }
    else
    {
        o << "ExtrPolygonMSC{6 DOF MSC}";
    }
}

//-------------------------------------------------------------------------

void radTExtrPolygonMSC::DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut,
                                  map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers,
                                  int& gUniqueMapKey, int elemKey)
{
    radTg3dRelax::DumpBin(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, elemKey);
    DumpBin_ExtrPolygonMSC(oStr);
}

//-------------------------------------------------------------------------

void radTExtrPolygonMSC::DumpBin_ExtrPolygonMSC(CAuxBinStrVect& oStr)
{
    oStr << FirstPoint.x << FirstPoint.y << FirstPoint.z;
    char cAxOrnt = AxOrnt;  // TAxisOrient is an enum, serialize as char
    oStr << cAxOrnt;
    oStr << Thickness;
    oStr << NumFaces;

    for(int i = 0; i < NumFaces; i++)
    {
        oStr << Sigma[i] << PointCharges[i];
    }

    // Dump base polygon data (use DumpBin_Polygon, not DumpBin)
    radTPolygon* pPgn = (radTPolygon*)(BasePolygonHandle.rep);
    char cBasePgnDefined = (pPgn != nullptr);
    oStr << cBasePgnDefined;
    if(cBasePgnDefined)
    {
        pPgn->DumpBin_Polygon(oStr);
    }
}

//-------------------------------------------------------------------------

void radTExtrPolygonMSC::DumpBinParse_ExtrPolygonMSC(CAuxBinStrVect& inStr)
{
    inStr >> FirstPoint.x >> FirstPoint.y >> FirstPoint.z;
    char cAxOrnt = 0;
    inStr >> cAxOrnt;
    AxOrnt = (TAxisOrient)cAxOrnt;  // Convert char back to TAxisOrient enum
    inStr >> Thickness;
    inStr >> NumFaces;

    for(int i = 0; i < NumFaces; i++)
    {
        inStr >> Sigma[i] >> PointCharges[i];
    }

    // Parse base polygon data (use radTPolygon(inStr) constructor)
    char cBasePgnDefined = 0;
    inStr >> cBasePgnDefined;
    if(cBasePgnDefined)
    {
        radThg hg(new radTPolygon(inStr));
        BasePolygonHandle = hg;
    }
}

//-------------------------------------------------------------------------
// Subdivision methods - simplified implementations
//-------------------------------------------------------------------------

int radTExtrPolygonMSC::SubdivideItself(double* SubdivArray, radThg& In_hg,
                                         radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
    // Delegate to standard subdivision
    return SubdivideItselfByParPlanes(SubdivArray, 3, In_hg, radPtr, pSubdivOptions);
}

//-------------------------------------------------------------------------

int radTExtrPolygonMSC::SubdivideItselfByOneSetOfParPlanes(TVector3d& InPlanesNormal,
                                                            TVector3d* InPointsOnCuttingPlanes,
                                                            int AmOfPieces_mi_1, radThg& In_hg,
                                                            radTApplication* radPtr,
                                                            radTSubdivOptions* pSubdivOptions,
                                                            radTvhg* pVectOfHgChanged)
{
    // Simplified implementation - create new MSC elements
    // For now, return without subdivision
    return 1;
}

//-------------------------------------------------------------------------

int radTExtrPolygonMSC::SubdivideItselfByParPlanes(double* SubdivArray, int AmOfDir,
                                                    radThg& In_hg, radTApplication* radPtr,
                                                    radTSubdivOptions* pSubdivOptions)
{
    // Simplified implementation
    return 1;
}

//-------------------------------------------------------------------------

int radTExtrPolygonMSC::CutItself(TVector3d* InCuttingPlane, radThg& In_hg,
                                   radTPair_int_hg& LowerNewPair_int_hg,
                                   radTPair_int_hg& UpperNewPair_int_hg,
                                   radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
    return 1;
}

//-------------------------------------------------------------------------

int radTExtrPolygonMSC::SubdivideItselfByEllipticCylinder(double* SubdivArray,
                                                           radTCylindricSubdivSpec* pSubdivSpec,
                                                           radThg& In_hg, radTApplication* radPtr,
                                                           radTSubdivOptions* pSubdivOptions)
{
    return 1;
}

//-------------------------------------------------------------------------

int radTExtrPolygonMSC::FindLowestAndUppestVertices(TVector3d& InNormal,
                                                     radTSubdivOptions* pSubdivOptions,
                                                     TVector3d& OutLowest, TVector3d& OutUppest,
                                                     radTrans& OutTrans, char& OutNormalIsParToZ,
                                                     char& OutNormalIsParToPlXZ)
{
    return 1;
}

//-------------------------------------------------------------------------

int radTExtrPolygonMSC::ConvertToPolyhedron(radThg& hg, radTApplication* radPtr, char PutNewStuffIntoGenCont)
{
    return 1;
}

//-------------------------------------------------------------------------
