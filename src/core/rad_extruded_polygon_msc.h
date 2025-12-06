/*-------------------------------------------------------------------------
*
* File name:      rad_extruded_polygon_msc.h
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
* Author(s):      Claude Code / ELF Corporation
*
* First release:  2025
*
* Copyright (C):  2025
*
-------------------------------------------------------------------------*/

#ifndef __RAD_EXTRUDED_POLYGON_MSC_H
#define __RAD_EXTRUDED_POLYGON_MSC_H

#include "rad_geometry_3d.h"
#include "rad_planar_2d.h"
#include "rad_transform_def.h"
#include <array>
#include <vector>

//-------------------------------------------------------------------------
// MSC Face data structure
//-------------------------------------------------------------------------

struct radTMSCFace {
    std::array<TVector3d, 4> vertices;  // Face vertices (4 for quad, 3 for triangle)
    int num_vertices;                    // 3 or 4
    TVector3d center;                    // Face center
    TVector3d normal;                    // Outward normal
    double area;                         // Face area

    radTMSCFace() : num_vertices(0), area(0.0) {
        center.x = center.y = center.z = 0.0;
        normal.x = normal.y = normal.z = 0.0;
    }

    void ComputeGeometry();  // Compute center, normal, area from vertices
};

//-------------------------------------------------------------------------
// radTExtrPolygonMSC - 6 DOF MSC-based extruded polygon
//-------------------------------------------------------------------------

class radTExtrPolygonMSC : public radTg3dRelax {
public:
    // Geometry
    TVector3d FirstPoint;
    TAxisOrient AxOrnt;
    radThg BasePolygonHandle;
    double Thickness;

    // MSC data
    std::array<radTMSCFace, 6> Faces;     // 6 faces of hexahedron
    std::array<TVector3d, 6> EvalPoints;  // Field evaluation points (midpoint)
    std::array<double, 6> Sigma;          // Surface charge densities (unknowns)
    std::array<double, 6> PointCharges;   // Point charges at center (= -sigma * area)
    int NumFaces;

    // Constructors
    radTExtrPolygonMSC(const TVector3d& InFirstPoint, const TAxisOrient& InAxOrnt, double InThickness,
                       const radThg& InBasePolygonHandle, const TVector3d& InMagn,
                       const radThg& InMaterHandle);

    radTExtrPolygonMSC(const TVector3d& InFirstPoint, const TAxisOrient& InAxOrnt, double InThickness,
                       TVector2d* ArrayOfPoints2d, int lenArrayOfPoints2d, const TVector3d& InMagn);

    radTExtrPolygonMSC(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers);

    radTExtrPolygonMSC() : radTg3dRelax(), NumFaces(0) {
        for(int i = 0; i < 6; i++) {
            Sigma[i] = 0.0;
            PointCharges[i] = 0.0;
        }
    }

    // Type identification
    int Type_g3dRelax() { return 2; }  // Same as original ExtrPolygon

    // Field computation - MSC method
    void B_comp(radTField*);
    void B_intComp(radTField*);
    void B_intCompSpecCases(radTField*, const TSpecCaseID&);

    // MSC-specific methods
    TVector3d FieldFromSurfaceCharge(const TVector3d& obs, const radTMSCFace& face, double sigma) const;
    TVector3d FieldFromPointCharge(const TVector3d& obs, const TVector3d& pos, double charge) const;
    TVector3d FieldFromChargedTriangle(const TVector3d& obs, const TVector3d& v0,
                                       const TVector3d& v1, const TVector3d& v2, double sigma) const;

    // Geometry setup
    void SetupFaces();
    void ComputeEvalPoints();
    void UpdatePointCharges();

    // Dump/IO
    void Dump(std::ostream&, int ShortSign = 0);
    void DumpPureObjInfo(std::ostream&, int);
    void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut,
                 map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey);
    void DumpBin_ExtrPolygonMSC(CAuxBinStrVect& oStr);
    void DumpBinParse_ExtrPolygonMSC(CAuxBinStrVect& inStr);

    // Graphics
    radTg3dGraphPresent* CreateGraphPresent();

    // Geometry queries
    double Volume() { return Thickness * (((radTPolygon*)(BasePolygonHandle.rep))->Area()); }
    void VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique);

    // Subdivision
    int SubdivideItself(double*, radThg&, radTApplication*, radTSubdivOptions*);
    int SubdivideItselfByOneSetOfParPlanes(TVector3d& InPlanesNormal, TVector3d* InPointsOnCuttingPlanes,
                                            int AmOfPieces_mi_1, radThg& In_hg, radTApplication* radPtr,
                                            radTSubdivOptions* pSubdivOptions, radTvhg* pVectOfHgChanged);
    int SubdivideItselfByParPlanes(double* SubdivArray, int AmOfDir, radThg& In_hg,
                                    radTApplication* radPtr, radTSubdivOptions* pSubdivOptions);
    int CutItself(TVector3d* InCuttingPlane, radThg& In_hg, radTPair_int_hg& LowerNewPair_int_hg,
                  radTPair_int_hg& UpperNewPair_int_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions);
    int SubdivideItselfByEllipticCylinder(double* SubdivArray, radTCylindricSubdivSpec* pSubdivSpec,
                                           radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions);

    int FindLowestAndUppestVertices(TVector3d&, radTSubdivOptions*, TVector3d&, TVector3d&, radTrans&, char&, char&);
    int ConvertToPolyhedron(radThg&, radTApplication*, char);

    int DuplicateItself(radThg& hg, radTApplication*, char) {
        return FinishDuplication(new radTExtrPolygonMSC(*this), hg);
    }

    void DefineCentrPoint() {
        TVector2d& CentrPo2d = ((radTPolygon*)(BasePolygonHandle.rep))->CentrPoint;
        TVector2d& BaseFirstPo2d = ((radTPolygon*)(BasePolygonHandle.rep))->EdgePointsVector[0];
        CentrPoint.x = FirstPoint.x + 0.5 * Thickness;
        CentrPoint.y = FirstPoint.y + CentrPo2d.x - BaseFirstPo2d.x;
        CentrPoint.z = FirstPoint.z + CentrPo2d.y - BaseFirstPo2d.y;
    }

    // DOF for MSC: 6 surface charge densities
    int NumberOfDegOfFreedom() { return 6; }

    int SizeOfThis() {
        int GenSize = sizeof(radTExtrPolygonMSC);
        GenSize += sizeof(radTPolygon);
        GenSize += ((radTPolygon*)(BasePolygonHandle.rep))->AmOfEdgePoints * sizeof(TVector2d);
        return GenSize;
    }

    // Get/Set sigma values (for solver interface)
    void GetSigma(double* sigma_out) const {
        for(int i = 0; i < NumFaces; i++) sigma_out[i] = Sigma[i];
    }

    void SetSigma(const double* sigma_in) {
        for(int i = 0; i < NumFaces; i++) Sigma[i] = sigma_in[i];
        UpdatePointCharges();
    }

    // Convert sigma to magnetization (for compatibility with existing code)
    // M = average magnetization computed from surface charges
    TVector3d ComputeAverageMagnetization();
};

//-------------------------------------------------------------------------

#endif // __RAD_EXTRUDED_POLYGON_MSC_H
