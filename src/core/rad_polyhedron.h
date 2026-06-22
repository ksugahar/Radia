/*-------------------------------------------------------------------------
*
* File name:      rad_polyhedron.h
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 polyhedron with constant magnetization
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#ifndef __RAD_POLYHEDRON_H
#define __RAD_POLYHEDRON_H

#include "rad_geometry_3d.h"
#include "rad_planar_2d.h"
#include "rad_transform_def.h"
#include "rad_convergence.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

extern radTConvergRepair& radCR;

//-------------------------------------------------------------------------
// Opt-in "improved yano-type" surface-charge (MSC) kernel (pyramid-cloud + pyramid-centroid eval).
// Set via rad.SolverConfig(yano_pyramid_cloud=True).  Default false == the EIEM2 single-point kernel
// (bit-identical to the historical default).  See radTPolyhedron::MscEvalPoint / MscCompensationField.
extern bool g_yano_pyramid_cloud;
extern bool g_yano_no_center_charge;   // research: drop the element-center cancellation charge (raw collocation)
extern double g_yano_eval_alpha;   // research: override EIEM2 collocation-point alpha (-1 = default 0.5)
extern bool g_yano_moment;   // upgrade: solve the parameter-free MOMENT system (BuildMomentSystemCore) instead of EIEM2 collocation (hex-only, opt-in)

//-------------------------------------------------------------------------

struct radTHandlePgnAndTrans {
	radTHandle<radTPolygon> PgnHndl;
	radTHandle<radTrans> TransHndl;
	bool FaceIsInternalAfterCut;

	radTHandlePgnAndTrans(radTHandle<radTPolygon>& inPgnHndl, radTHandle<radTrans>& inTransHndl) 
	{
		PgnHndl = inPgnHndl; TransHndl = inTransHndl;
	}
	radTHandlePgnAndTrans() { FaceIsInternalAfterCut = false;}
	inline friend int operator <(const radTHandlePgnAndTrans&, const radTHandlePgnAndTrans&);
	inline friend int operator ==(const radTHandlePgnAndTrans&, const radTHandlePgnAndTrans&);
};

//-------------------------------------------------------------------------

inline int operator <(const radTHandlePgnAndTrans&, const radTHandlePgnAndTrans&) { return 1;}

//-------------------------------------------------------------------------

inline int operator ==(const radTHandlePgnAndTrans& h1, const radTHandlePgnAndTrans& h2)
{
	return (h1.PgnHndl.rep == h2.PgnHndl.rep) && (h1.TransHndl.rep == h2.TransHndl.rep);
}

//-------------------------------------------------------------------------

#ifdef __GNUC__
using radTVectHandlePgnAndTrans = vector<radTHandlePgnAndTrans>;
using radTVectOfPtrToVect3d = vector<std::array<TVector3d, 2>>;
using radTListOfVector3d = list<TVector3d>;
using radTVectVect3d = vector<TVector3d>;
#else
using radTVectHandlePgnAndTrans = vector<radTHandlePgnAndTrans, allocator<radTHandlePgnAndTrans> >;
using radTVectOfPtrToVect3d = vector<std::array<TVector3d, 2>, allocator<std::array<TVector3d, 2>> >;
using radTListOfVector3d = list<TVector3d, allocator<TVector3d> >;
using radTVectVect3d = vector<TVector3d, allocator<TVector3d> >;
#endif

#ifdef __MWERKS__
/*
null_template
struct iterator_traits <TVector3d*> {
	 typedef ptrdiff_t difference_type;
	 typedef TVector3d value_type;
	 typedef TVector3d* pointer;
	 typedef TVector3d& reference;
	 typedef random_access_iterator_tag iterator_category;
};
null_template
struct iterator_traits <radTHandlePgnAndTrans*> {
	 typedef ptrdiff_t difference_type;
	 typedef radTHandlePgnAndTrans value_type;
	 typedef radTHandlePgnAndTrans* pointer;
	 typedef radTHandlePgnAndTrans& reference;
	 typedef random_access_iterator_tag iterator_category;
};
null_template
struct iterator_traits <TVector3d**> {
	 typedef ptrdiff_t difference_type;
	 typedef TVector3d* value_type;
	 typedef TVector3d** pointer;
	 typedef TVector3d*& reference;
	 typedef random_access_iterator_tag iterator_category;
};
*/
#endif

//-------------------------------------------------------------------------

class radTPolyhedron : public radTg3dRelax {

	//const TMatrix3d* pJ_LinCoef;
	TMatrix3d* pJ_LinCoef;
	char mLinTreat; //0- treat as relative

public:
	radTVectHandlePgnAndTrans VectHandlePgnAndTrans;
	int AmOfFaces;

	TVector3d J; //to move to base?
	bool J_IsNotZero;

	short SomethingIsWrong;
	radTPairOfDouble AuxPairOfDouble; // Used for cylindrical subdivision

	// MSC (Magnetic Surface Charge) support for hexahedra and wedges
	// For 5+ face elements, we use surface charge density (sigma) on each face
	// instead of magnetization vector (Mx, My, Mz)
	// Hexahedra: 6 faces -> 6 DOF, Wedges: 5 faces -> 5 DOF
	double Sigma[6];           // Surface charge densities for each face (max 6 DOF)
	double FaceArea[6];        // Face areas (max 6)
	TVector3d FaceNormal[6];   // Face normals (outward, max 6)
	TVector3d FaceCenter[6];   // Face centers (max 6)
	bool Use6DOF_MSC;          // Flag: true if element uses per-face MSC (5+ faces)
	double CurrentChi;         // Chi used for current solve (for H = M/chi update)

	radTPolyhedron(TVector3d* ArrayOfPoints, int lenArrayOfPoints, int** ArrayOfFaces, int* ArrayOfLengths, int lenArrayOfFaces, const TVector3d& InMagn)
		: radTg3dRelax(InMagn)
	{
		//Magn = InMagn; AmOfFaces = lenArrayOfFaces; SomethingIsWrong = 0;
		AmOfFaces = lenArrayOfFaces; SomethingIsWrong = 0;
		pJ_LinCoef = 0; mLinTreat = 0;
		J_IsNotZero = false;

		// Initialize 6 DOF MSC data for hexahedra
		Use6DOF_MSC = (AmOfFaces >= 5);  // Wedges (5) and hexahedra (6)
		CurrentChi = 1.0;  // Default chi
		for(int i = 0; i < 6; i++) {
			Sigma[i] = 0.0;
			FaceArea[i] = 0.0;
			FaceNormal[i].x = FaceNormal[i].y = FaceNormal[i].z = 0.0;
			FaceCenter[i].x = FaceCenter[i].y = FaceCenter[i].z = 0.0;
		}

		//DefineCentrPoint(ArrayOfPoints, lenArrayOfPoints); //OC090908
		ShiftFacesNumeration(ArrayOfFaces, ArrayOfLengths);
		FillInVectHandlePgnAndTrans(ArrayOfPoints, lenArrayOfPoints, ArrayOfFaces, ArrayOfLengths);
		if(SomethingIsWrong) return;
		DefineCentrPoint(ArrayOfPoints, lenArrayOfPoints);

		// Setup face geometry for 6 DOF MSC
		if(Use6DOF_MSC) SetupFaceGeometry();
	}
	radTPolyhedron(TVector3d* ArrayOfPoints, int lenArrayOfPoints, int** ArrayOfFaces, int* ArrayOfLengths, int lenArrayOfFaces,
		const TVector3d& InMagn, TMatrix3d& InM_LinCoef, TVector3d& InJ, TMatrix3d& InJ_LinCoef, char LinTreat)
		: radTg3dRelax(InMagn, InM_LinCoef)
	{
		AmOfFaces = lenArrayOfFaces; SomethingIsWrong = 0;

		// Initialize 6 DOF MSC data for hexahedra
		Use6DOF_MSC = (AmOfFaces >= 5);  // Wedges (5) and hexahedra (6)
		CurrentChi = 1.0;  // Default chi
		for(int i = 0; i < 6; i++) {
			Sigma[i] = 0.0;
			FaceArea[i] = 0.0;
			FaceNormal[i].x = FaceNormal[i].y = FaceNormal[i].z = 0.0;
			FaceCenter[i].x = FaceCenter[i].y = FaceCenter[i].z = 0.0;
		}

		//DefineCentrPoint(ArrayOfPoints, lenArrayOfPoints); //OC090908
		ShiftFacesNumeration(ArrayOfFaces, ArrayOfLengths);
		FillInVectHandlePgnAndTrans(ArrayOfPoints, lenArrayOfPoints, ArrayOfFaces, ArrayOfLengths);
		if(SomethingIsWrong) return;
		DefineCentrPoint(ArrayOfPoints, lenArrayOfPoints);

		// Setup face geometry for 6 DOF MSC
		if(Use6DOF_MSC) SetupFaceGeometry();

		J = InJ;
		bool J_LinCoefIsNotZero = !InJ_LinCoef.isZero();
		J_IsNotZero = (!InJ.isZero()) || J_LinCoefIsNotZero; //??
		if(J_LinCoefIsNotZero) pJ_LinCoef = new TMatrix3d(InJ_LinCoef);
		else pJ_LinCoef = 0;

		mLinTreat = LinTreat;
	}
	radTPolyhedron(TVector3d** ArrayOfFaces, int* ArrayOfLengths, int lenArrayOfFaces, const TVector3d& InMagn)
		: radTg3dRelax(InMagn)
	{
		//Magn = InMagn; AmOfFaces = lenArrayOfFaces; SomethingIsWrong = 0;
		AmOfFaces = lenArrayOfFaces; SomethingIsWrong = 0;
		pJ_LinCoef = 0; mLinTreat = 0;
		J_IsNotZero = false;

		// Initialize 6 DOF MSC data for hexahedra
		Use6DOF_MSC = (AmOfFaces >= 5);  // Wedges (5) and hexahedra (6)
		for(int i = 0; i < 6; i++) {
			Sigma[i] = 0.0;
			FaceArea[i] = 0.0;
			FaceNormal[i].x = FaceNormal[i].y = FaceNormal[i].z = 0.0;
			FaceCenter[i].x = FaceCenter[i].y = FaceCenter[i].z = 0.0;
		}

		TVector3d* OutArrayOfPoints;
		int lenArrayOfPoints;
		int** OutArrayOfFaces;
		MakeNormalPresentation(ArrayOfFaces, ArrayOfLengths, OutArrayOfPoints, lenArrayOfPoints, OutArrayOfFaces);
		if(SomethingIsWrong) { DeleteInputArrays(OutArrayOfPoints, OutArrayOfFaces); return;}

		//DefineCentrPoint(OutArrayOfPoints, lenArrayOfPoints); //OC090908
		FillInVectHandlePgnAndTrans(OutArrayOfPoints, lenArrayOfPoints, OutArrayOfFaces, ArrayOfLengths);
		if(SomethingIsWrong) { DeleteInputArrays(OutArrayOfPoints, OutArrayOfFaces); return;}
		DefineCentrPoint(OutArrayOfPoints, lenArrayOfPoints);
		DeleteInputArrays(OutArrayOfPoints, OutArrayOfFaces);

		// Setup face geometry for 6 DOF MSC
		if(Use6DOF_MSC) SetupFaceGeometry();
	}
	radTPolyhedron(const radTVectHandlePgnAndTrans& InVectHandlePgnAndTrans,
		const TVector3d* pInMagn, TMatrix3d* pInM_LinCoef, const radThg& InMatHandle,
		const TVector3d* pInJ, TMatrix3d* pInJ_LinCoef, char LinTreat, const TVector3d* pPrevLinRefP) //used at cutting / subdivision
	//radTPolyhedron(const radTVectHandlePgnAndTrans& InVectHandlePgnAndTrans, const TVector3d& InMagn, const radThg& InMatHandle)
	//radTPolyhedron(const radTVectHandlePgnAndTrans& InVectHandlePgnAndTrans, const TVector3d& InCentrPoint, const TVector3d& InMagn, const radThg& InMatHandle)
		: radTg3dRelax(pInMagn, pInM_LinCoef, InMatHandle)
	{
		//pJ_LinCoef = 0; mLinTreat = 0;
		//J_IsNotZero = false;

		AmOfFaces = (int)InVectHandlePgnAndTrans.size();
		for(int i=0; i<AmOfFaces; i++) VectHandlePgnAndTrans.push_back(InVectHandlePgnAndTrans[i]);
		SomethingIsWrong = 0;
		DefineCentrPoint();

		// Initialize 6 DOF MSC data for hexahedra
		Use6DOF_MSC = (AmOfFaces >= 5);  // Wedges (5) and hexahedra (6)
		for(int i = 0; i < 6; i++) {
			Sigma[i] = 0.0;
			FaceArea[i] = 0.0;
			FaceNormal[i].x = FaceNormal[i].y = FaceNormal[i].z = 0.0;
			FaceCenter[i].x = FaceCenter[i].y = FaceCenter[i].z = 0.0;
		}

		// Setup face geometry for 6 DOF MSC
		if(Use6DOF_MSC) SetupFaceGeometry();

		J_IsNotZero = false;
		J.Zero();
		pJ_LinCoef = 0;
		if(pInJ != 0)
		{
			J = *pInJ;
			if(!J.isZero()) J_IsNotZero = true;
		}
		bool J_LinCoef_AreNotZero = false;
		if(pInJ_LinCoef != 0)
		{
			if(!pInJ_LinCoef->isZero())
			{
				pJ_LinCoef = new TMatrix3d(*pInJ_LinCoef);
				J_LinCoef_AreNotZero = true;
				J_IsNotZero = true;
			}
		}
		mLinTreat = LinTreat;

		if((LinTreat == 0) && (pPrevLinRefP != 0) && ((pM_LinCoef != 0) || J_LinCoef_AreNotZero)) //Linear dependence is Relative w.r. to object Center
		{//Do this correction only after M and/or J are defined!
		 //Consider moving this to base class(es)
			TVector3d dCenP = CentrPoint - (*pPrevLinRefP);
			if(pM_LinCoef != 0)
			{
				if(!pM_LinCoef->isZero()) Magn += ((*pM_LinCoef)*dCenP);
			}
			if(J_LinCoef_AreNotZero) J += ((*pJ_LinCoef)*dCenP);
		}
	}
	radTPolyhedron(const radTHandlePgnAndTrans& inHandleBasePgnAndTrf1, const radTHandlePgnAndTrans& inHandleBasePgnAndTrf2, double avgCur=0, double* arMagComp=0)
	{//tries to generate a convex polyhedron from two base face polygons
		SomethingIsWrong = 0;
		mLinTreat = 0; J_IsNotZero = false;
		J.x = J.y = J.z = 0.; pJ_LinCoef = 0;
		Magn.x = Magn.y = Magn.z = 0.; pM_LinCoef = 0;

		// Initialize 6 DOF MSC data (will be updated after AmOfFaces is known)
		Use6DOF_MSC = false;
		for(int i = 0; i < 6; i++) {
			Sigma[i] = 0.0;
			FaceArea[i] = 0.0;
			FaceNormal[i].x = FaceNormal[i].y = FaceNormal[i].z = 0.0;
			FaceCenter[i].x = FaceCenter[i].y = FaceCenter[i].z = 0.0;
		}

		AttemptToCreateConvexPolyhedronFromTwoBaseFaces(inHandleBasePgnAndTrf1, inHandleBasePgnAndTrf2);
		if(SomethingIsWrong) return;

		// Now that AmOfFaces is set, update 6 DOF flag for hexahedra
		Use6DOF_MSC = (AmOfFaces >= 5);  // Wedges (5) and hexahedra (6)
		if(Use6DOF_MSC) SetupFaceGeometry();

		if(avgCur != 0)
		{
			SetCurrentDensityForConstCurrent(avgCur, 0, 1); //may set J and pJ_LinCoef
		}
		if(arMagComp != 0)
		{
			Magn.x = *arMagComp;
			Magn.y = *(arMagComp + 1);
			Magn.z = *(arMagComp + 2);
		}
	}
	radTPolyhedron(radTPolyhedron& aPlhdr) : radTg3dRelax(aPlhdr)
	{
		*this = aPlhdr;
		if(aPlhdr.pJ_LinCoef != 0) pJ_LinCoef = new TMatrix3d(*(aPlhdr.pJ_LinCoef));
	}
	// radTPolyhedron(CAuxBinStrVect&, ...) REMOVED (Phase B2c, 2026-04-15)
	radTPolyhedron() : radTg3dRelax()
	{
		pJ_LinCoef = 0; mLinTreat = 0;
		J_IsNotZero = false;
		SomethingIsWrong = 0;

		// Initialize 6 DOF MSC data
		Use6DOF_MSC = false;
		AmOfFaces = 0;
		for(int i = 0; i < 6; i++) {
			Sigma[i] = 0.0;
			FaceArea[i] = 0.0;
			FaceNormal[i].x = FaceNormal[i].y = FaceNormal[i].z = 0.0;
			FaceCenter[i].x = FaceCenter[i].y = FaceCenter[i].z = 0.0;
		}
	}
	~radTPolyhedron() 
	{
		if(pJ_LinCoef != 0) delete pJ_LinCoef;
	}

	int Type_g3dRelax() { return 5;}
	// DOF: AmOfFaces for MSC elements (sigma per face), 3 for tetrahedra (Mx, My, Mz)
	// Hexahedra: 6 DOF, Wedges: 5 DOF, Tetrahedra: 3 DOF
	// Returns 0 if no material is applied (same behavior as radTRecMag)
	int NumberOfDegOfFreedom() { return (MaterHandle.rep == 0) ? 0 : (Use6DOF_MSC ? AmOfFaces : 3); }

	void FillInVectHandlePgnAndTrans(TVector3d*, int, int**, int*);
	void MakeNormalPresentation(TVector3d**, int*, TVector3d*&, int&, int**&);
	int CheckIfFacePolygonsArePlanar(TVector3d*, int**, int*, TVector3d*);
	int DetermineActualFacesNormals(TVector3d*, int, int**, int*, TVector3d*);
	int FillInTransAndFacesInLocFrames(TVector3d*, int**, int*, TVector3d*);

	//void B_comp(radTField*);
	//void B_intComp(radTField*);

	void B_comp_frM(radTField*);
	void B_comp_frJ(radTField*);
	void B_intComp_frM(radTField*);
	void B_intComp_frJ(radTField*);

	// Element type detection (MSC method support)
	bool IsTetrahedron() const { return AmOfFaces == 4; }
	bool IsWedge() const { return AmOfFaces == 5; }
	bool IsHexahedron() const { return AmOfFaces == 6; }
	bool IsMSCElement() const { return AmOfFaces >= 5; }

	// MSC (Magnetic Surface Charge) methods for supported element types
	void B_comp_tetrahedron_analytical(radTField*);
	void B_comp_wedge_analytical(radTField*);  // 3DOF wedge/prism (5 faces: 2 tri + 3 quad)
	void B_comp_wedge_MSC(radTField*);         // 5DOF wedge MSC (sigma per face)
	void B_comp_hexahedron_MSC(radTField*);

	// 6 DOF MSC field computation for hexahedra
	TVector3d FieldFromChargedTriangle(const TVector3d& obs, const TVector3d& v0,
	                                    const TVector3d& v1, const TVector3d& v2, double sigma) const;
	// Version with explicit normal (for IMA boundary faces where computed normal is wrong)
	TVector3d FieldFromChargedTriangleWithNormal(const TVector3d& obs, const TVector3d& v0,
	                                              const TVector3d& v1, const TVector3d& v2,
	                                              double sigma, const TVector3d& explicitNormal) const;
	TVector3d FieldFromFace(const TVector3d& obs, int faceIdx, double sigma) const;  // Generalized: tri or quad
	TVector3d FieldFromQuadFace(const TVector3d& obs, int faceIdx, double sigma) const;
	TVector3d FieldFromQuadFaceMirrored(const TVector3d& obs, const TVector3d& V0,
	                                     const TVector3d& V1, const TVector3d& V2,
	                                     const TVector3d& V3, double sigma,
	                                     bool flipNormal, const TVector3d& mirrorCenter) const;
	// Version with explicit normals for IMA boundary faces
	TVector3d FieldFromQuadFaceMirroredWithNormals(const TVector3d& obs, const TVector3d& V0,
	                                                const TVector3d& V1, const TVector3d& V2,
	                                                const TVector3d& V3, double sigma,
	                                                const TVector3d& tri1Normal,
	                                                const TVector3d& tri2Normal) const;
	TVector3d FieldFromPointCharge(const TVector3d& obs, double charge) const;
	// Generalized mirrored face field: dispatches to tri or quad based on numVerts
	TVector3d FieldFromFaceMirrored(const TVector3d& obs,
	                                 const TVector3d* mirrorVerts, int numVerts,
	                                 double sigma, bool flipNormal,
	                                 const TVector3d& mirrorCenter) const;

	// --- "improved yano-type" MSC kernel (opt-in via g_yano_pyramid_cloud) ---
	// MscEvalPoint: the per-face collocation evaluation point.
	//   default (flag off): EIEM2 midpoint 0.5*(FaceCenter[i] + CentrPoint);
	//   flag on: pyramid centroid 0.75*MscFaceAreaCentroid(i) + 0.25*MscVolumeCentroid().
	// MscCompensationField: the source-face charge-neutralising compensation field (no 1/4pi).
	//   default (flag off): single point charge -FaceArea[i] at the element center (FieldFromPointCharge);
	//   flag on: element-common cloud (all (volume-centroid, edge) partition triangles, 3-pt quad,
	//            normalised), total charge -FaceArea[i] -> per-DOF charge-neutral, loop-source-null.
	// Flag-off branches return EXACTLY the historical inline expressions (bit-identical goldens).
	TVector3d FieldFromPointChargeAt(const TVector3d& obs, const TVector3d& src, double charge) const;
	TVector3d MscFaceAreaCentroid(int faceIdx) const;
	TVector3d MscVolumeCentroid() const;
	TVector3d MscEvalPoint(int faceIdx) const;
	TVector3d MscCompensationField(const TVector3d& obs, int faceIdx) const;
	// Element-common compensation cloud (the loop-source-null pyramid cloud): fills pts[] (global) and
	// wts[] (normalised, sum=1) with up to `cap` (volume-centroid, edge) 3-pt-quadrature points; returns
	// the count.  Single source for both MscCompensationField (inline path) and the precomputed
	// Compute6x6BlockFast fast path, so the two agree by construction.  For a hexahedron count==72.
	int MscCompensationCloud(TVector3d* pts, double* wts, int cap) const;

	// 6 DOF MSC setup for hexahedra
	// IMPORTANT: This relies on Netgen face winding convention for correct normal direction.
	// No inside/outside check is performed - the normal is computed mechanically from
	// the polygon's local coordinate system which was set up from vertex winding order.
	// Face ordering (Netgen convention): 0=z-, 1=x+, 2=y-, 3=x-, 4=y+, 5=z+
	void SetupFaceGeometry()
	{
		// Compute face normals, areas, and centers for MSC elements
		// Supports wedges (5 faces) and hexahedra (6 faces)
		if(AmOfFaces < 5 || AmOfFaces > 6) return;

		for(int i = 0; i < AmOfFaces; i++)
		{
			radTHandlePgnAndTrans& hPgnTrans = VectHandlePgnAndTrans[i];
			radTPolygon* pPgn = hPgnTrans.PgnHndl.rep;
			radTrans* pTrans = hPgnTrans.TransHndl.rep;

			// Get face center in local frame, then transform to global
			TVector2d& locCP = pPgn->CentrPoint;
			TVector3d localCenter(locCP.x, locCP.y, pPgn->CoordZ);
			FaceCenter[i] = pTrans->TrBiPoint(localCenter);

			// Face normal is the Z-axis of the local frame transformed to global
			// The local frame was set up from the face vertex winding order,
			// so the Z-axis (0,0,1) transformed to global gives the outward normal.
			// This is the Netgen convention - no inside/outside check needed.
			TVector3d localNormal(0.0, 0.0, 1.0);
			FaceNormal[i] = pTrans->TrBiPoint(localNormal) - pTrans->TrBiPoint(TVector3d(0, 0, 0));

			// Normalize the normal vector
			double normLen = sqrt(FaceNormal[i].x * FaceNormal[i].x +
			                      FaceNormal[i].y * FaceNormal[i].y +
			                      FaceNormal[i].z * FaceNormal[i].z);
			if(normLen > 1e-15)
			{
				FaceNormal[i].x /= normLen;
				FaceNormal[i].y /= normLen;
				FaceNormal[i].z /= normLen;
			}

			// NOTE: Previous implementation had inside/outside check here:
			// if(dotProd < 0.0) { flip normal }
			// This was removed to follow Netgen convention strictly.
			// The face winding order determines normal direction mechanically.

			// Compute face area from polygon
			FaceArea[i] = pPgn->Area();
		}
	}

	// CutItself / FindIntersectionWithFace / SetUpUpperAndLowerPolygon / DetermineNewFaceAndTrans
	// FillInNewHandlePgnAndTransFrom3d / CheckIfTwoPointAlreadyMapped / SubdivideItself* /
	// KsFromSizeToNumb / DeterminePointsOnCuttingPlanes / FindLowestAndUppestVertices /
	// IntrsctOfTwoLines / EstimateSize / FindEdgePoints* / FindLocalEllipticCoord /
	// SubdivideByEllipses / SubdivideItselfOverAzimuth / SubdivideItselfByEllipticCylinder
	// REMOVED (Phase C, 2026-04-16)
	void DefineRelAndAbsTol(double*);

	int CheckForSpecialShapes(radTVectHandlePgnAndTrans&, radThg&, double*);

	double Volume();
	void VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique);

	void FindTypicalSize(TVector3d*, int, double&);

	// Dump / DumpPureObjInfo / DumpBin / DumpBin_Polyhedron / DumpBinParse_Polyhedron REMOVED (Phase B2b/B2c, 2026-04-15)

	void Push_backCenterPointAndField(radTFieldKey* pFieldKey, radTVectPairOfVect3d* pVectPairOfVect3d, radTrans* pBaseTrans, radTg3d* g3dSrcPtr, radTApplication* pAppl);
	
	void AttemptToCreateConvexPolyhedronFromTwoBaseFaces(const radTHandlePgnAndTrans& inHandleBasePgnAndTrf1, const radTHandlePgnAndTrans& inHandleBasePgnAndTrf2);
	bool CheckIfAllPolygonVerticesAreOnOneSideOfPlane(const radTHandlePgnAndTrans& hPgnAndTrf, const TVector3d& vPoint, const TVector3d& vNorm, double AbsTol);
	void CollectAndMapUniquePolygonPoints(const radTHandlePgnAndTrans& hPgnAndTrf, vector<TVector3d>& vectPoints, vector<int>& vectInd, double AbsTol);
	void GenerateSideFacesContainingSegmentsOfBaseFace(const vector<TVector3d>& vectVertexPoints, vector<vector<int> >& vectIndAllFaces, int indBaseFace, double* arTol);
	void ReorderPointsToEnsureNonSelfIntersectingPolygon(const vector<TVector3d>& vectPoints, vector<int>& vectIndPgnPoints, double* arTol);
	void SetCurrentDensityForConstCurrent(double avgCur, int indBaseFace1, int indBaseFace2);

	int ConvertToPolyhedron(radThg&, radTApplication*, char) { return 1;} // 1 is essential

	// Override B_genComp for field computation
	void B_genComp(radTField* pField) override;

	void B_comp(radTField* pField)
	{
		bool M_IsNotZero = !Magn.isZero();
		if(M_IsNotZero || (pField->FieldKey.PreRelax_)) B_comp_frM(pField);
		if(J_IsNotZero) B_comp_frJ(pField);
	}
	void B_intComp(radTField* pField)
	{
		bool M_IsNotZero = !Magn.isZero();
		if(M_IsNotZero) B_intComp_frM(pField);
		if(J_IsNotZero) B_intComp_frJ(pField);
	}

	int DuplicateItself(radThg& hg, radTApplication*, char)
	{
		return FinishDuplication(new radTPolyhedron(*this), hg);
	}
	int SizeOfThis()
	{
		int GenSize = sizeof(radTPolyhedron);
		int BufSize = sizeof(radTrans);
		GenSize += AmOfFaces*BufSize;
		for(int i=0; i<AmOfFaces; i++) GenSize += (VectHandlePgnAndTrans[i].PgnHndl.rep)->SizeOfThis();
		return GenSize;
	}
	void DefineCentrPoint(TVector3d* ArrayOfPoints, int AmOfPoints)
	{// Modify later if necessary
		TVector3d Sum(0.,0.,0.);
		for(int k=0; k<AmOfPoints; k++) Sum = Sum + ArrayOfPoints[k];
		double Buf = 1./AmOfPoints;
		CentrPoint = TVector3d(radCR.Double(Buf*Sum.x), radCR.Double(Buf*Sum.y), radCR.Double(Buf*Sum.z));
	}
	void DefineCentrPoint()
	{// This algorithm differs from the above (and may give different results)
		TVector3d Sum(0.,0.,0.);
		for(int k=0; k<AmOfFaces; k++)
		{
			radTHandlePgnAndTrans& HandlePgnAndTrans = VectHandlePgnAndTrans[k];
			radTPolygon* PgnPtr = HandlePgnAndTrans.PgnHndl.rep;
			TVector2d& LocFaceCP = PgnPtr->CentrPoint;
			TVector3d FaceCP(LocFaceCP.x, LocFaceCP.y, PgnPtr->CoordZ);
			Sum = Sum + HandlePgnAndTrans.TransHndl.rep->TrBiPoint(FaceCP);
		}
		double Buf = 1./double(AmOfFaces);
		CentrPoint = TVector3d(radCR.Double(Buf*Sum.x), radCR.Double(Buf*Sum.y), radCR.Double(Buf*Sum.z));
	}

	//void ReCalcCentrPointFromPgnAndTrans(const radTVectHandlePgnAndTrans& vHandlePgnAndTrans, const TVector3d& oldCenPoint, TVector3d& newCenPoint) //OC090908
	//{// This algorithm differs from the above (and may give different results)
	//	TVector3d Sum(0.,0.,0.);
	//	int locAmOfFaces = (int)vHandlePgnAndTrans.size();
	//	if(locAmOfFaces <= 0) return;
	//	for(int k=0; k<locAmOfFaces; k++)
	//	{
	//		const radTHandlePgnAndTrans& hPgnAndTrans = vHandlePgnAndTrans[k];
	//		radTPolygon* pPgn = hPgnAndTrans.PgnHndl.rep;
	//		TVector2d& LocFaceCP = pPgn->CentrPoint;
	//		TVector3d FaceCP(LocFaceCP.x, LocFaceCP.y, pPgn->CoordZ);
	//		Sum = Sum + hPgnAndTrans.TransHndl.rep->TrBiPoint(FaceCP);
	//	}
	//	double Buf = 1./double(locAmOfFaces);
	//	//CentrPoint = TVector3d(radCR.Double(Buf*Sum.x), radCR.Double(Buf*Sum.y), radCR.Double(Buf*Sum.z));
	//	newCenPoint.x = radCR.Double(Buf*Sum.x); //??
	//	newCenPoint.y = radCR.Double(Buf*Sum.y); 
	//	newCenPoint.z = radCR.Double(Buf*Sum.z);
	//	newCenPoint += oldCenPoint;
	//}

	//void CorrectFacePolygonsForNewCenPoint(radTVectHandlePgnAndTrans& vHandlePgnAndTrans, const TVector3d& difCenPoints) //OC090908
	//{//difCenPoints = Old - New
	//	int locAmOfFaces = (int)vHandlePgnAndTrans.size();
	//	if(locAmOfFaces <= 0) return;
	//	for(int k=0; k<locAmOfFaces; k++)
	//	{
	//		radTHandlePgnAndTrans& hPgnAndTrans = vHandlePgnAndTrans[k];
	//		radTPolygon* pPgn = hPgnAndTrans.PgnHndl.rep;
	//		radTrans* pTrans = hPgnAndTrans.TransHndl.rep;
	//		TVector3d addVertexPoint = pTrans->TrBiPoint(difCenPoints);
	//		pPgn->CoordZ += addVertexPoint.z;
	//		TVector2d addVertexPoint2d(addVertexPoint.x, addVertexPoint.y);
	//		pPgn->CentrPoint += addVertexPoint2d;
	//		radTVect2dVect& curVectEdgePoints = pPgn->EdgePointsVector;
	//		for(int j=0; j<pPgn->AmOfEdgePoints; j++)
	//		{
	//			curVectEdgePoints[j] += addVertexPoint2d;
	//		}
	//	}
	//}

	void ShiftFacesNumeration(int** ArrayOfFaces, int* ArrayOfLengths)
	{
		for(int i=0; i<AmOfFaces; i++)
		{
			int* CurrentFace = ArrayOfFaces[i];
			for(int j=0; j<ArrayOfLengths[i]; j++) (CurrentFace[j])--;
		}
	}
	void DefineNormalVia3Points(const TVector3d& P1, const TVector3d& P2, const TVector3d& P3, TVector3d& Normal)
	{
		TVector3d R1 = P2 - P1, R2 = P3 - P1;
		Normal.x = R1.y*R2.z - R2.y*R1.z;
		Normal.y = R2.x*R1.z - R1.x*R2.z;
		Normal.z = R1.x*R2.y - R2.x*R1.y;
	}
	double Vect3dNorm(const TVector3d& V)
	{
		double AbsX = Abs(V.x), AbsY = Abs(V.y), AbsZ = Abs(V.z);
		double MaxXY = (AbsX>AbsY)? AbsX : AbsY;
		return (MaxXY>AbsZ)? MaxXY : AbsZ;
	}
	int NextCircularNumber(int CurrentNo, int Total)
	{
		return (CurrentNo == Total-1)? 0 : CurrentNo + 1;
	}
	void ReverseArrayOfInt(int* ArrayOfInt, int lenArrayOfInt)
	{
		int *DirPtr = ArrayOfInt, *RevPtr = &(ArrayOfInt[lenArrayOfInt-1]);
		for(int i=0; i < (lenArrayOfInt >> 1); i++)
		{
			int Buf = *DirPtr; *(DirPtr++) = *RevPtr; *(RevPtr--) = Buf;
		}
	}
	void ReverseArrayOfVect3dPtr(TVector3d** ArrayOfVect3dPtr, int lenArrayOfVect3dPtr)
	{
		TVector3d** DirPtr = ArrayOfVect3dPtr;
		TVector3d** RevPtr = &(ArrayOfVect3dPtr[lenArrayOfVect3dPtr-1]);
		for(int i=0; i < (lenArrayOfVect3dPtr >> 1); i++)
		{
			TVector3d* Buf = *DirPtr; *(DirPtr++) = *RevPtr; *(RevPtr--) = Buf;
		}
	}
	int CheckIfJunctionIsConvex(const TVector3d& N1, const TVector3d& JointSegm, const TVector3d& N2)
	{
		TVector3d N1_vect_by_N2(N1.y*N2.z - N2.y*N1.z, N2.x*N1.z - N1.x*N2.z, N1.x*N2.y - N2.x*N1.y);
		return (N1_vect_by_N2*JointSegm >= 0.)? 1 : 0; // Or ">"?
	}
	void DeleteInputArrays(TVector3d* ArrayOfPoints, int** ArrayOfFaces, int* ArrayOfLengths =nullptr)
	{
		if(ArrayOfFaces != nullptr)
		{
			for(int i=0; i<AmOfFaces; i++) delete[] (ArrayOfFaces[i]);
			delete[] ArrayOfFaces;
		}
		if(ArrayOfLengths != nullptr) delete[] ArrayOfLengths;
		if(ArrayOfPoints != nullptr) delete[] ArrayOfPoints;

		if(pJ_LinCoef != 0) delete pJ_LinCoef;
		pJ_LinCoef = 0;
	}
	void DeleteAuxInputArrays(TVector3d** ArrayOfFaces)
	{
		if(ArrayOfFaces != nullptr)
		{
			for(int i=0; i<AmOfFaces; i++) delete[] (ArrayOfFaces[i]);
			delete[] ArrayOfFaces;
		}
	}
	void DeleteAuxInputArrays(short** ArrayOfFaces)
	{
		if(ArrayOfFaces != nullptr)
		{
			for(int i=0; i<AmOfFaces; i++) delete[] (ArrayOfFaces[i]);
			delete[] ArrayOfFaces;
		}
	}
	int ItemIsNotFullyInternalAfterCut()
	{
		for(int k=0; k<AmOfFaces; k++)
			if(!VectHandlePgnAndTrans[k].FaceIsInternalAfterCut) return 1;
		return 0;
	}
	int CreateNewEntity(radTVectHandlePgnAndTrans& vHandlePgnAndTrans, radThg& hg, short RecognizeRecMagsInPolyhedrons, double* RelAbsTol)
	{
		short CreateA_Polyhedron = 1;
		if(RecognizeRecMagsInPolyhedrons)
			if(CheckForSpecialShapes(vHandlePgnAndTrans, hg, RelAbsTol)) CreateA_Polyhedron = 0;
		if(CreateA_Polyhedron)
		{
			radTSend Send;
			radTPolyhedron* PolyhedronPtr = new radTPolyhedron(vHandlePgnAndTrans, &Magn, pM_LinCoef, MaterHandle, &J, pJ_LinCoef, mLinTreat, &CentrPoint);
			//radTPolyhedron* PolyhedronPtr = new radTPolyhedron(vHandlePgnAndTrans, Magn, MaterHandle);

			if(PolyhedronPtr == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}
			hg = radThg(PolyhedronPtr);
		}
		return 1;
	}
	int FindTwoOrtogonalVectors(TVector3d& InV, TVector3d* TwoVect)
	{
		TVector3d V1(0.,0.,0.);
		if(!((InV.x==0.) && (InV.y==0.))) { V1.x = -InV.y; V1.y = InV.x;}
		else if(!((InV.x==0.) && (InV.z==0.))) { V1.x = -InV.z; V1.z = InV.x;}
		else if(!((InV.y==0.) && (InV.z==0.))) { V1.y = -InV.z; V1.z = InV.y;}
		else return 0;
		TVector3d V2 = InV^V1;
		*TwoVect = V1; *(TwoVect+1) = V2; 
		return 1;
	}
	char CheckIfOnlyNeighbouringEdgePointsTrapped(int* IntersectingBoundsNos, int AmOfEdgePo) //OC291003
	{
		if((AmOfEdgePo <= 0) || (IntersectingBoundsNos == 0)) return 0;
		int FirstNo = IntersectingBoundsNos[0], SecondNo = IntersectingBoundsNos[1];

		if((FirstNo == SecondNo) || (SecondNo == NextCircularNumber(FirstNo, AmOfEdgePo)) || (FirstNo == NextCircularNumber(SecondNo, AmOfEdgePo))) return 1;
		else return 0;
	}
	int ScaleCurrent(double scaleCoef) //OC250713 //virtual in g3d
	{//note: if(scaleCoef == 0) this still doesn't change J_IsNotZero
		if(J_IsNotZero) 
		{
			J *= scaleCoef; 

			if(pJ_LinCoef != 0)
			{
				*pJ_LinCoef *= scaleCoef; //??
			}
			return 1;
		}
		else return 0;
	}
};

//-------------------------------------------------------------------------

#endif
