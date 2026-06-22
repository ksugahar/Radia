/*-------------------------------------------------------------------------
*
* File name:      radvlpgn.cpp
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

#include "rad_serialization.h"
#include "rad_polyhedron.h"
#include <array>
#include <cstdlib>    // std::getenv
#include <string>     // std::string
#include <stdexcept>  // std::runtime_error
// #include "rad_subdivided_polyhedron.h" REMOVED (Phase C, 2026-04-16)
#include "rad_geometry_3d_aux.h"
#include "rad_application.h"
#include "auxparse.h"
#include "rad_poly_analytical.h"
#include "radentry.h"  // For RadSolverGetTetraMethod()
#include "rad_point_classify.h"  // For solid angle point-in-polyhedron test
#include "rad_constants.h"  // Unified mathematical/physical constants
#include "rad_interaction.h"  // For RadIMAFieldContext (IMA field computation)

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTPolyhedron::FillInVectHandlePgnAndTrans(TVector3d* ArrayOfPoints, int lenArrayOfPoints, int** ArrayOfFaces, int* ArrayOfLengths)
{
	radTSend Send;

	std::vector<TVector3d> vArrayOfFacesNormals(AmOfFaces);
	TVector3d* ArrayOfFacesNormals = vArrayOfFacesNormals.data();

	if(!CheckIfFacePolygonsArePlanar(ArrayOfPoints, ArrayOfFaces, ArrayOfLengths, ArrayOfFacesNormals)) return;
	if(!DetermineActualFacesNormals(ArrayOfPoints, lenArrayOfPoints, ArrayOfFaces, ArrayOfLengths, ArrayOfFacesNormals)) return;
	if(!FillInTransAndFacesInLocFrames(ArrayOfPoints, ArrayOfFaces, ArrayOfLengths, ArrayOfFacesNormals)) return;
	// RAII: automatic cleanup
}

//-------------------------------------------------------------------------

void radTPolyhedron::MakeNormalPresentation(TVector3d** ArrayOfFaces, int* ArrayOfLengths, TVector3d*& OutArrayOfPoints, int& lenArrayOfPoints, int**& OutArrayOfFaces)
{
	radTVectVect3d VectOfPoints;
	radTSend Send;

	OutArrayOfFaces = 0;
	OutArrayOfFaces = new int*[AmOfFaces];
	if(OutArrayOfFaces == 0) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error900"); return;}

	int PointsCount = 0;
	for(int i=0; i<AmOfFaces; i++)
	{
		TVector3d* CurrentFace = ArrayOfFaces[i];

		int* CurrentFaceInt = new int[ArrayOfLengths[i]];
		if(CurrentFaceInt == 0) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error900"); return;}
		OutArrayOfFaces[i] = CurrentFaceInt;

		for(int k=0; k<ArrayOfLengths[i]; k++)
		{
			TVector3d& CurrentPoint = CurrentFace[k];
			short PointAlreadyThere = 0;
			int LocCount = -1;
			for(radTVectVect3d::iterator Iter = VectOfPoints.begin(); Iter != VectOfPoints.end(); ++Iter)
			{
				++LocCount;
				TVector3d& PointInCont = *Iter;
				if((PointInCont.x == CurrentPoint.x) && (PointInCont.y == CurrentPoint.y) && (PointInCont.z == CurrentPoint.z))
				{
					PointAlreadyThere = 1; break;
				}
			}
			if(PointAlreadyThere) CurrentFaceInt[k] = LocCount;
			else
			{
				VectOfPoints.push_back(CurrentPoint);
				CurrentFaceInt[k] = PointsCount;
				++PointsCount;
			}
		}
	}
	lenArrayOfPoints = PointsCount;
	OutArrayOfPoints = new TVector3d[lenArrayOfPoints];
	if(OutArrayOfPoints == 0) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error900"); return;}
	for(int j=0; j<lenArrayOfPoints; j++)
	{
		OutArrayOfPoints[j] = VectOfPoints[j];
	}
}

//-------------------------------------------------------------------------

int radTPolyhedron::CheckIfFacePolygonsArePlanar(TVector3d* ArrayOfPoints, int** ArrayOfFaces, int* ArrayOfLengths, TVector3d* ArrayOfFacesNormals)
{
	radTSend Send;
	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	const double RelLengthTol = RelAbsTol[0];

	for(int i=0; i<AmOfFaces; i++)
	{
		int* CurrentFace = ArrayOfFaces[i];

		TVector3d& P0 = ArrayOfPoints[CurrentFace[0]];
		TVector3d& P1 = ArrayOfPoints[CurrentFace[1]];
		TVector3d Normal;

		double BufAbsLenTol0 = RelLengthTol*Vect3dNorm(P0);
		double BufAbsLenTol1 = RelLengthTol*Vect3dNorm(P1);
		double AbsLenTol = (BufAbsLenTol1>BufAbsLenTol0)? BufAbsLenTol1 : BufAbsLenTol0;

		int CurrentLength = ArrayOfLengths[i];
		if(CurrentLength < 3) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error048"); return 0;}
		
		int j;
		short NormalDetermined = 0;
		for(j=2; j<CurrentLength; j++)
		{
			TVector3d& P2 = ArrayOfPoints[CurrentFace[j]];

			BufAbsLenTol0 = RelLengthTol*Vect3dNorm(P2);
			if(BufAbsLenTol0>AbsLenTol) AbsLenTol = BufAbsLenTol0;

			double AbsLenTol_e2 = AbsLenTol*AbsLenTol;

			DefineNormalVia3Points(P0, P1, P2, Normal);
			if(Abs(Normal.x)>AbsLenTol_e2 || Abs(Normal.y)>AbsLenTol_e2 || Abs(Normal.z)>AbsLenTol_e2)
			{
				NormalDetermined = 1; break;
			}
		}
		if(!NormalDetermined) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error049"); return 0;}

		for(int k=j+1; k<CurrentLength; k++)
		{
			TVector3d& aPoint = ArrayOfPoints[CurrentFace[k]];
			TVector3d anR = aPoint - P0;

			double AbsBufVal = Abs(Normal*anR);
			double CompareValue = Vect3dNorm(Normal)*AbsLenTol;

			if(AbsBufVal > CompareValue) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error047"); return 0;}
		}
		ArrayOfFacesNormals[i] = Normal;
	}
	return 1;
}

//-------------------------------------------------------------------------

void radTPolyhedron::FindTypicalSize(TVector3d* ArrayOfPoints, int AmOfPoints, double& TypicalSize)
{
	double MaxX = 0., MaxY = 0., MaxZ = 0.;
	TVector3d& FirstPo = *ArrayOfPoints;
	for(int i=1; i<AmOfPoints; i++)
	{
		TVector3d& aPoint = ArrayOfPoints[i];
		double rX = fabs(aPoint.x - FirstPo.x), rY = fabs(aPoint.y - FirstPo.y), rZ = fabs(aPoint.z - FirstPo.z);
		if(rX > MaxX) MaxX = rX;
		if(rY > MaxY) MaxY = rY;
		if(rZ > MaxZ) MaxZ = rZ;
	}
	TypicalSize = (MaxX > MaxY)? MaxX : MaxY;
	if(TypicalSize < MaxZ) TypicalSize = MaxZ;
}

//-------------------------------------------------------------------------

int radTPolyhedron::DetermineActualFacesNormals(TVector3d* ArrayOfPoints, int AmOfPoints, int** ArrayOfFaces, int* ArrayOfLengths, TVector3d* ArrayOfFacesNormals)
{
	radTSend Send;
	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	const double RelLengthTol = RelAbsTol[0];
	double TypicalSize;
	FindTypicalSize(ArrayOfPoints, AmOfPoints, TypicalSize);
	double AbsLengthTol = RelLengthTol*TypicalSize;

	short GoodFaceFound = 0, NormalJustReversed = 0;

// Looking for the first good face
	int i;
	for(i=0; i<AmOfFaces; i++)
	{
		int* CurrentFace = ArrayOfFaces[i];
		TVector3d& Normal = ArrayOfFacesNormals[i];

		double InvNorm = 1./sqrt(Normal.x*Normal.x + Normal.y*Normal.y + Normal.z*Normal.z);
		Normal = InvNorm*Normal;

		TVector3d& P0 = ArrayOfPoints[CurrentFace[0]];

		short AtLeastOnePointIsOK = 0, SomePointsAreNotOK = 0;
		NormalJustReversed = 0;

		int j;
		for(j=0; j<AmOfPoints; j++)
		{
			TVector3d anR = P0 - ArrayOfPoints[j];
			double ScalProd = Normal*anR;

			if(Abs(ScalProd) > AbsLengthTol)
			{
				if(ScalProd < -AbsLengthTol) 
				{
					if(!AtLeastOnePointIsOK)
					{
						Normal = (-1.)*Normal; NormalJustReversed = 1; 
					}
					else SomePointsAreNotOK = 1;
					break;
				}
				else AtLeastOnePointIsOK = 1;
			}
		}
		if(NormalJustReversed)
			for(j=0; j<AmOfPoints; j++)
			{
				TVector3d anR = P0 - ArrayOfPoints[j];
				double ScalProd = Normal*anR;
				if(Normal*anR < -AbsLengthTol) { SomePointsAreNotOK = 1; break;}
			}
		if(!SomePointsAreNotOK) { GoodFaceFound = 1; break;}
	}
	if(!GoodFaceFound) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error050"); return 0;}
	int NoOfFirstGoodFace = i;

// Determining all the rest normals
	std::vector<std::vector<short>> vSegmentPassed(AmOfFaces);
	std::vector<short*> vSegmentPassedPtrs(AmOfFaces);
	for(i=0; i<AmOfFaces; i++)
	{
		int CurrentLength = ArrayOfLengths[i];
		vSegmentPassed[i].resize(CurrentLength, 0);
		vSegmentPassedPtrs[i] = vSegmentPassed[i].data();
	}
	short** SegmentPassed = vSegmentPassedPtrs.data();
	int ii = NoOfFirstGoodFace;

	std::vector<char> vGenFacesPassed(AmOfFaces, 0);
	char* GenFacesPassed = vGenFacesPassed.data();
	std::vector<char> vPossibleNextFaces(AmOfFaces);
	char* PossibleNextFaces = vPossibleNextFaces.data();

	for(i=0; i<AmOfFaces; i++)
	{
		for(int p=0; p<AmOfFaces; p++) PossibleNextFaces[p] = 0;

		int CurrentLength = ArrayOfLengths[ii];
		int* CurrentFace = ArrayOfFaces[ii];

		TVector3d& FirstNormal = ArrayOfFacesNormals[ii];
		if(NormalJustReversed)
		{ 
			ReverseArrayOfInt(CurrentFace, CurrentLength); NormalJustReversed = 0;
		}

		short* SegmentPassedOnCurrentFace = SegmentPassed[ii];
		for(int j=0; j<CurrentLength; j++)
		{
			int NoOfSegmStPo = CurrentFace[j]; 
			int NoOfSegmFiPo = CurrentFace[NextCircularNumber(j, CurrentLength)];
			
			if(!(SegmentPassedOnCurrentFace[j]))
			{
				short ThisSegmentIsSomeWhereElse = 0;
				for(int iii=0; iii<AmOfFaces; iii++)
				{
					if(iii != ii)
					{
						int LocCurrentLength = ArrayOfLengths[iii];
						int* LocCurrentFace = ArrayOfFaces[iii];

						int jjSt = -1, jjFi = -1;
						for(int jj=0; jj<LocCurrentLength; jj++)
						{
							int CurrentNo = LocCurrentFace[jj];
							if(CurrentNo == NoOfSegmStPo)
							{
								jjSt = jj; if(jjFi > -1) break;
							}
							if(CurrentNo == NoOfSegmFiPo)
							{
								jjFi = jj; if(jjSt > -1) break;
							}
						}
						if((jjSt>-1) && (jjFi>-1))
						{
							ThisSegmentIsSomeWhereElse = 1;
							short& ThisSegmentAlreadyPassed = (SegmentPassed[iii])[jjFi];
							if(ThisSegmentAlreadyPassed) break;

							PossibleNextFaces[iii] = 1;

							if(NextCircularNumber(jjFi, LocCurrentLength) != jjSt)
							{
								ReverseArrayOfInt(LocCurrentFace, LocCurrentLength);
								TVector3d& Normal = ArrayOfFacesNormals[iii];
								Normal = (-1.)*Normal;
							}
							TVector3d SegmVect = ArrayOfPoints[NoOfSegmFiPo] - ArrayOfPoints[NoOfSegmStPo];
							if(!CheckIfJunctionIsConvex(FirstNormal, SegmVect, ArrayOfFacesNormals[iii]))
							{// Modify this if non-convex volumes are supported or treated specially
								SomethingIsWrong = 1;
								Send.ErrorMessage("Radia::Error106");
								// RAII: vSegmentPassed will be automatically cleaned up
								return 0;
							}
							ThisSegmentAlreadyPassed = 1; break;
						}
					}
				}
				if(!ThisSegmentIsSomeWhereElse) // Modify this if unclosed volumes are not supported
				{
					Send.WarningMessage("Radia::Warning013"); 
				}
			}
		}
		GenFacesPassed[ii] = 1;

		//int NextFaceNo;
		int NextFaceNo = -1; //OC 100902
		for(int pp=0; pp<AmOfFaces; pp++)
		{
			if((PossibleNextFaces[pp] == 1) && (!GenFacesPassed[pp])) { NextFaceNo = pp; break;}
		}

		if(NextFaceNo >= 0)//OC 100902
		{
			ii = NextFaceNo;
		}
	}
	// RAII: vSegmentPassed, vGenFacesPassed, and vPossibleNextFaces will be automatically cleaned up
	return 1;
}

//-------------------------------------------------------------------------

int radTPolyhedron::FillInTransAndFacesInLocFrames(TVector3d* ArrayOfPoints, int** ArrayOfFaces, int* ArrayOfLengths, TVector3d* ArrayOfFacesNormals)
{
	radTSend Send;
	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	const double RelLengthTol = RelAbsTol[0];

	for(int i=0; i<AmOfFaces; i++)
	{
		TVector3d& N = ArrayOfFacesNormals[i];
		//double AbsLengthTol = RelLengthTol*Vect3dNorm(N); //OC
		double InvLen = 1./sqrt(N.x*N.x + N.y*N.y + N.z*N.z);
		N = InvLen*N;

		TVector3d St1, St2, St3;
		//if(Abs(N.z+1.) > AbsLengthTol) //OC
		if(Abs(N.z+1.) > RelLengthTol) //OC
		{
			double InvNzp1 = 1./(N.z + 1.);
			St1 = TVector3d(N.y*N.y*InvNzp1 + N.z, -N.x*N.y*InvNzp1, -N.x);
			St2 = TVector3d(St1.y, N.x*N.x*InvNzp1 + N.z, -N.y);
			St3 = TVector3d(-St1.z, -St2.z, N.z);
		}
		else
		{
			St1 = TVector3d(1., 0., 0.);
			St2 = TVector3d(0., -1., 0.);
			St3 = TVector3d(0., 0., -1.);
		}
		TMatrix3d R1(St1, St2, St3);

		int k;
		int* CurrentFace = ArrayOfFaces[i];
		int CurrentLength = ArrayOfLengths[i];

		//TVector3d EdgeVect; //OC
		//for(k=0; k<CurrentLength; k++)
		//{
		//	EdgeVect = ArrayOfPoints[CurrentFace[NextCircularNumber(k, CurrentLength)]] - ArrayOfPoints[CurrentFace[k]];
		//	double EdgeVectLength = sqrt(EdgeVect.x*EdgeVect.x + EdgeVect.y*EdgeVect.y + EdgeVect.z*EdgeVect.z);
		//	if(EdgeVectLength > AbsLengthTol)
		//	{
		//		EdgeVect = (1./EdgeVectLength)*EdgeVect; break;
		//	}
		//}
		//TVector3d V = R1*EdgeVect;

// Uncomment the following to speed-up a bit (however, loss of precision in FieldInt computation was encountered due to this)
		//St1.x = V.y; St1.y = -V.x; St1.z = 0.;
		//St2.x = V.x; St2.y = V.y; St2.z = 0.;
		//St3.x = 0.; St3.y = 0.; St3.z = 1.;
		//TMatrix3d R2(St1, St2, St3);
		//TMatrix3d R = R2*R1;

		TMatrix3d R = R1;
		TVector3d Zero(0.,0.,0.);

		radTrans* RotationPtr = new radTrans(R, Zero, 1., 1., 2);
		if(RotationPtr == 0) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error900"); return 0;}

		radTVect2dVect Vect2dVect;
		double LocCoordZ;
		for(k=0; k<CurrentLength; k++)
		{
			TVector3d P = ArrayOfPoints[CurrentFace[k]];
			//TVector3d P = ArrayOfPoints[CurrentFace[k]] - CentrPoint; //OC090908 assuming that CentrPoint was already defined
			//ATTENTION: this modification requires updates in all Field computation and Subdivision routines !!! //OC090908

			TVector3d P_loc = RotationPtr->TrPoint(P);
			TVector2d P2d(P_loc.x, P_loc.y);
			Vect2dVect.push_back(P2d);
			LocCoordZ = P_loc.z;
		}

		TVector3d LocMagn = RotationPtr->TrVectField(Magn);
		radTPolygon* FacePgnPtr = new radTPolygon(Vect2dVect, LocCoordZ, LocMagn);
		if(FacePgnPtr == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}
		radTHandle<radTPolygon> HandlePgn(FacePgnPtr);

		RotationPtr->Invert();
		radTHandle<radTrans> HandleRotat(RotationPtr);

		radTHandlePgnAndTrans HandlePgnAndTrans;
		HandlePgnAndTrans.PgnHndl = HandlePgn;
		HandlePgnAndTrans.TransHndl = HandleRotat;
		VectHandlePgnAndTrans.push_back(HandlePgnAndTrans);
	}
	return 1;
}

//-------------------------------------------------------------------------
// Tetrahedral mesh support (MSC - Magnetic Surface Charge method)
//-------------------------------------------------------------------------

void radTPolyhedron::B_comp_tetrahedron_analytical(radTField* FieldPtr)
{
	// =========================================================================
	// GLOBAL COORDINATE METHOD for tetrahedral elements
	// =========================================================================
	// This method computes fields using GLOBAL coordinates only.
	// No local coordinate transformation is used, avoiding transformation errors.
	//
	// For each triangular face:
	//   1. Get vertices in GLOBAL coordinates
	//   2. Compute face normal in GLOBAL coordinates
	//   3. Compute H field contribution using RadFieldFromTriangleFaceGlobal
	//   4. Sum contributions from all 4 faces
	//
	// For PreRelax mode (interaction matrix construction):
	//   - Diagonal (self-interaction): Set to 0.5
	//   - Off-diagonal: Compute dH/dM directly in global coordinates

	radTFieldKey& FldKey = FieldPtr->FieldKey;
	TVector3d& obsPoint = FieldPtr->P;

	// Check for self-interaction (observation point at centroid)
	bool isSelfInteraction = false;
	{
		TVector3d diff = obsPoint - CentrPoint;
		double distSq = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
		// Estimate element size from face centers
		double elemSizeSq = 0.0;
		for(int i = 0; i < AmOfFaces; i++)
		{
			// Use const reference to avoid radTHandle copy (important for OpenMP thread safety)
			const radTHandlePgnAndTrans& hpt = VectHandlePgnAndTrans[i];
			radTPolygon* pgn = hpt.PgnHndl.rep;
			radTrans* tr = hpt.TransHndl.rep;
			TVector3d faceCtr = tr->TrPoint(TVector3d(pgn->CentrPoint.x, pgn->CentrPoint.y, pgn->CoordZ));
			TVector3d d = faceCtr - CentrPoint;
			double dSq = d.x*d.x + d.y*d.y + d.z*d.z;
			if(dSq > elemSizeSq) elemSizeSq = dSq;
		}
		double tol = elemSizeSq * 1.0e-10;
		if(tol < 1.0e-20) tol = 1.0e-20;
		isSelfInteraction = (distSq < tol);
	}

	// NOTE: Self-interaction is computed by the formula below, not hardcoded.
	// The geometry of the tetrahedron determines the actual demagnetization factors.
	// Hardcoding 0.5 was incorrect - the formula gives correct values like ~0.1 for
	// the diagonal elements (depends on tetrahedron shape).

	// Get face vertices in GLOBAL coordinates
	// For a tetrahedron, we have 4 triangular faces
	// Use fixed-size array to avoid heap allocation (important for OpenMP thread safety)
	std::array<std::array<TVector3d, 3>, 4> faceVertices;  // Max 4 faces for tetrahedron
	int nFaces = (AmOfFaces <= 4) ? AmOfFaces : 4;
	for(int i = 0; i < nFaces; i++)
	{
		// Use const reference to avoid radTHandle copy (which modifies reference count)
		// This is critical for OpenMP thread safety
		const radTHandlePgnAndTrans& hpt = VectHandlePgnAndTrans[i];
		radTPolygon* pgn = hpt.PgnHndl.rep;
		radTrans* tr = hpt.TransHndl.rep;

		// Get 2D vertices from polygon and transform to global 3D
		const radTVect2dVect& verts2d = pgn->EdgePointsVector;
		if(verts2d.size() < 3) continue;

		// Transform each vertex: local 2D (x, y, CoordZ) -> global 3D
		faceVertices[i][0] = tr->TrPoint(TVector3d(verts2d[0].x, verts2d[0].y, pgn->CoordZ));
		faceVertices[i][1] = tr->TrPoint(TVector3d(verts2d[1].x, verts2d[1].y, pgn->CoordZ));
		faceVertices[i][2] = tr->TrPoint(TVector3d(verts2d[2].x, verts2d[2].y, pgn->CoordZ));
	}

	if(FldKey.PreRelax_)
	{
		// =====================================================================
		// PreRelax mode: Compute interaction matrix coefficients
		// =====================================================================
		// We need to compute dH/dM for off-diagonal elements.
		//
		// For each unit magnetization direction (Mx=1, My=1, Mz=1), compute
		// the H field contribution and store in the appropriate row.
		//
		// Store as: B.row = dH/dMx, H.row = dH/dMy, A.row = dH/dMz

		TVector3d unit_Mx(1., 0., 0.);
		TVector3d unit_My(0., 1., 0.);
		TVector3d unit_Mz(0., 0., 1.);

		TVector3d H_from_Mx(0., 0., 0.);
		TVector3d H_from_My(0., 0., 0.);
		TVector3d H_from_Mz(0., 0., 0.);

		// Sum contributions from all faces
		// Note: CentrPoint is the element centroid, used for outward normal check
		for(int i = 0; i < nFaces; i++)
		{
			const TVector3d& V0 = faceVertices[i][0];
			const TVector3d& V1 = faceVertices[i][1];
			const TVector3d& V2 = faceVertices[i][2];

			H_from_Mx += RadFieldFromTriangleFaceGlobal(V0, V1, V2, unit_Mx, obsPoint, CentrPoint);
			H_from_My += RadFieldFromTriangleFaceGlobal(V0, V1, V2, unit_My, obsPoint, CentrPoint);
			H_from_Mz += RadFieldFromTriangleFaceGlobal(V0, V1, V2, unit_Mz, obsPoint, CentrPoint);
		}

		// Store in the matrix format expected by Radia:
		// B = row for dH/dMx (Hx, Hy, Hz when Mx=1)
		// H = row for dH/dMy (Hx, Hy, Hz when My=1)
		// A = row for dH/dMz (Hx, Hy, Hz when Mz=1)
		// NOTE: Use += to match polygon PreRelax sign convention (which stores N_physical directly)
		// The polygon's B_comp stores H_field without negation in PreRelax mode.
		FieldPtr->B.x += H_from_Mx.x;
		FieldPtr->B.y += H_from_Mx.y;
		FieldPtr->B.z += H_from_Mx.z;
		FieldPtr->H.x += H_from_My.x;
		FieldPtr->H.y += H_from_My.y;
		FieldPtr->H.z += H_from_My.z;
		FieldPtr->A.x += H_from_Mz.x;
		FieldPtr->A.y += H_from_Mz.y;
		FieldPtr->A.z += H_from_Mz.z;
	}
	else
	{
		// =====================================================================
		// Normal mode: Compute actual H field from current magnetization
		// =====================================================================
		TVector3d H_total(0., 0., 0.);

		// Note: CentrPoint is the element centroid, used for outward normal check
		for(int i = 0; i < nFaces; i++)
		{
			const TVector3d& V0 = faceVertices[i][0];
			const TVector3d& V1 = faceVertices[i][1];
			const TVector3d& V2 = faceVertices[i][2];

			H_total += RadFieldFromTriangleFaceGlobal(V0, V1, V2, Magn, obsPoint, CentrPoint);
		}

		if(FldKey.H_) FieldPtr->H += H_total;

		// For B field, check if observation point is inside the tetrahedron
		// If inside, B = mu_0 * (H + M), so we add M to B
		if(FldKey.B_)
		{
			FieldPtr->B += H_total;

			// Extract 4 unique vertices from face vertices for inside test
			// Each face shares vertices with other faces, so we collect unique ones
			TVector3d tetraVerts[4];
			tetraVerts[0] = faceVertices[0][0];
			tetraVerts[1] = faceVertices[0][1];
			tetraVerts[2] = faceVertices[0][2];
			// The 4th vertex is on face 1, 2, or 3 but not on face 0
			// Face 1 shares 2 vertices with face 0, so one vertex is new
			for(int j = 0; j < 3; j++) {
				const TVector3d& v = faceVertices[1][j];
				bool isNew = true;
				for(int k = 0; k < 3; k++) {
					TVector3d diff;
					diff.x = v.x - tetraVerts[k].x;
					diff.y = v.y - tetraVerts[k].y;
					diff.z = v.z - tetraVerts[k].z;
					if(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z < 1e-20) {
						isNew = false;
						break;
					}
				}
				if(isNew) {
					tetraVerts[3] = v;
					break;
				}
			}

			// Check if observation point is inside the tetrahedron
			bool pointInside = RadPointClassify::PointInTetrahedronSolidAngle(
				obsPoint, tetraVerts[0], tetraVerts[1], tetraVerts[2], tetraVerts[3]);

			if(pointInside) {
				// B = mu_0 * (H + M), H_total is H, so add M to B
				FieldPtr->B += Magn;
			}
		}

		// M field (magnetization at observation point)
		if(FldKey.M_)
		{
			TVector3d tetraVerts[4];
			tetraVerts[0] = faceVertices[0][0];
			tetraVerts[1] = faceVertices[0][1];
			tetraVerts[2] = faceVertices[0][2];
			for(int j = 0; j < 3; j++) {
				const TVector3d& v = faceVertices[1][j];
				bool isNew = true;
				for(int k = 0; k < 3; k++) {
					TVector3d diff;
					diff.x = v.x - tetraVerts[k].x;
					diff.y = v.y - tetraVerts[k].y;
					diff.z = v.z - tetraVerts[k].z;
					if(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z < 1e-20) {
						isNew = false;
						break;
					}
				}
				if(isNew) {
					tetraVerts[3] = v;
					break;
				}
			}

			bool pointInside = RadPointClassify::PointInTetrahedronSolidAngle(
				obsPoint, tetraVerts[0], tetraVerts[1], tetraVerts[2], tetraVerts[3]);

			if(pointInside) {
				FieldPtr->M += Magn;
			}
		}

		// =====================================================================
		// Scalar potential (phi) and vector potential (A) computation
		// =====================================================================
		// Uses FACE INTEGRATION (not dipole approximation)
		//
		// Vector potential: A(r) = (mu_0/4pi) * integral_S (M x n) / |r-r'| dS
		//                   B = curl(A)
		//
		// Each triangular face contributes to the total vector potential.
		//
		// NOTE: On symmetry axes, the face-based integration may give A=0 due to
		// symmetric cancellation. This differs from ObjRecMag which uses an
		// 8-corner BufVect formula that doesn't cancel. This is mathematically
		// correct behavior for the face-based approach.
		// =====================================================================
		if(FldKey.A_)
		{
			TVector3d A_total(0., 0., 0.);

			// Sum contributions from all triangular faces
			for(int i = 0; i < nFaces; i++)
			{
				const TVector3d& V0 = faceVertices[i][0];
				const TVector3d& V1 = faceVertices[i][1];
				const TVector3d& V2 = faceVertices[i][2];

				A_total += RadVectorPotentialFromTriangleFaceGlobal(
					V0, V1, V2, Magn, obsPoint, CentrPoint);
			}

			FieldPtr->A += A_total;
		}

		// Scalar potential (phi) using face-based integration (no dipole approximation)
		// Phi = (1/4pi) * M dot BufVect, where BufVect = sum(n * integral(1/|r-r'|) dS)
		if(FldKey.Phi_)
		{
			double Phi_total = 0.0;

			// Sum contributions from all triangular faces
			for(int i = 0; i < nFaces; i++)
			{
				const TVector3d& V0 = faceVertices[i][0];
				const TVector3d& V1 = faceVertices[i][1];
				const TVector3d& V2 = faceVertices[i][2];

				Phi_total += RadScalarPotentialFromTriangleFaceGlobal(
					V0, V1, V2, Magn, obsPoint, CentrPoint);
			}

			FieldPtr->Phi += Phi_total;
		}

		// =====================================================================
		// IMA (Image Method) Field Contributions for tetrahedral elements
		// Same pattern as B_comp_hexahedron_MSC IMA permanent magnet path
		// Magnetization is a pseudo-vector: mirror flips component along axis
		// =====================================================================
		if(RadIMAFieldContext::IsActive() && !FldKey.PreRelax_)
		{
			int imaSym = RadIMAFieldContext::GetSymmetry();
			int signX = RadIMAFieldContext::GetSignX();
			int signY = RadIMAFieldContext::GetSignY();
			int signZ = RadIMAFieldContext::GetSignZ();

			// Helper lambda to compute field from mirrored tetrahedral geometry
			auto computeMirroredField = [&](int mirrorAxis, int sign) -> TVector3d {
				TVector3d H_mirror(0., 0., 0.);

				// Create mirrored face vertices
				std::array<std::array<TVector3d, 3>, 4> mirrorVerts;
				for(int i = 0; i < nFaces; i++)
				{
					for(int j = 0; j < 3; j++)
					{
						mirrorVerts[i][j] = faceVertices[i][j];
						if(mirrorAxis & radTInteraction::IMA_X) mirrorVerts[i][j].x = -mirrorVerts[i][j].x;
						if(mirrorAxis & radTInteraction::IMA_Y) mirrorVerts[i][j].y = -mirrorVerts[i][j].y;
						if(mirrorAxis & radTInteraction::IMA_Z) mirrorVerts[i][j].z = -mirrorVerts[i][j].z;
					}
				}

				// Count mirror axes: odd number requires winding reversal
				int numAxes = 0;
				if(mirrorAxis & radTInteraction::IMA_X) numAxes++;
				if(mirrorAxis & radTInteraction::IMA_Y) numAxes++;
				if(mirrorAxis & radTInteraction::IMA_Z) numAxes++;
				bool reverseWinding = (numAxes % 2 == 1);

				// Mirrored center point
				TVector3d mirrorCenter = CentrPoint;
				if(mirrorAxis & radTInteraction::IMA_X) mirrorCenter.x = -mirrorCenter.x;
				if(mirrorAxis & radTInteraction::IMA_Y) mirrorCenter.y = -mirrorCenter.y;
				if(mirrorAxis & radTInteraction::IMA_Z) mirrorCenter.z = -mirrorCenter.z;

				// Mirror magnetization: pseudo-vector flip + BC sign
				TVector3d mirrorMagn = Magn;
				if(mirrorAxis & radTInteraction::IMA_X) mirrorMagn.x = -mirrorMagn.x;
				if(mirrorAxis & radTInteraction::IMA_Y) mirrorMagn.y = -mirrorMagn.y;
				if(mirrorAxis & radTInteraction::IMA_Z) mirrorMagn.z = -mirrorMagn.z;
				mirrorMagn = mirrorMagn * (double)sign;

				// Compute field from mirrored triangular faces
				for(int i = 0; i < nFaces; i++)
				{
					const TVector3d& MV0 = mirrorVerts[i][0];
					const TVector3d& MV1 = mirrorVerts[i][1];
					const TVector3d& MV2 = mirrorVerts[i][2];

					if(reverseWinding)
						H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV2, MV1, mirrorMagn, obsPoint, mirrorCenter);
					else
						H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV1, MV2, mirrorMagn, obsPoint, mirrorCenter);
				}

				return H_mirror;
			};

			TVector3d H_ima(0., 0., 0.);

			// Single axis contributions
			if(imaSym & radTInteraction::IMA_X)
				H_ima += computeMirroredField(radTInteraction::IMA_X, signX);
			if(imaSym & radTInteraction::IMA_Y)
				H_ima += computeMirroredField(radTInteraction::IMA_Y, signY);
			if(imaSym & radTInteraction::IMA_Z)
				H_ima += computeMirroredField(radTInteraction::IMA_Z, signZ);

			// Dual axis contributions
			if((imaSym & radTInteraction::IMA_X) && (imaSym & radTInteraction::IMA_Y))
				H_ima += computeMirroredField(radTInteraction::IMA_XY, signX * signY);
			if((imaSym & radTInteraction::IMA_X) && (imaSym & radTInteraction::IMA_Z))
				H_ima += computeMirroredField(radTInteraction::IMA_XZ, signX * signZ);
			if((imaSym & radTInteraction::IMA_Y) && (imaSym & radTInteraction::IMA_Z))
				H_ima += computeMirroredField(radTInteraction::IMA_YZ, signY * signZ);

			// Triple axis contribution
			if((imaSym & radTInteraction::IMA_X) && (imaSym & radTInteraction::IMA_Y) && (imaSym & radTInteraction::IMA_Z))
				H_ima += computeMirroredField(radTInteraction::IMA_XYZ, signX * signY * signZ);

			// Add IMA contributions
			if(FldKey.H_) FieldPtr->H += H_ima;
			if(FldKey.B_) FieldPtr->B += H_ima;
		}
	}
}

//-------------------------------------------------------------------------

void radTPolyhedron::B_comp_wedge_analytical(radTField* FieldPtr)
{
	// =========================================================================
	// GLOBAL COORDINATE METHOD for wedge/prism elements (5 faces: 2 tri + 3 quad)
	// =========================================================================
	// This method computes fields using GLOBAL coordinates only.
	// Wedges are 3DOF elements (same as tetrahedra) so they use the same
	// magnetization-based approach.
	//
	// Face types:
	//   - 2 triangular end faces (3 vertices each)
	//   - 3 quadrilateral side faces (4 vertices each, split into 2 triangles)
	//
	// For PreRelax mode (interaction matrix construction):
	//   - Compute dH/dM for all 3 unit magnetization directions
	//
	// For Normal mode:
	//   - Compute H field from current magnetization

	radTFieldKey& FldKey = FieldPtr->FieldKey;
	TVector3d& obsPoint = FieldPtr->P;

	// Get face vertices in GLOBAL coordinates
	// For a wedge, we have 5 faces with variable number of vertices
	// Use fixed-size array to avoid heap allocation (better for OpenMP)
	std::array<int, 8> faceNumVerts;          // Number of vertices per face
	std::array<std::array<TVector3d, 4>, 8> faceVertices;  // Max 4 vertices per face
	int nFaces = (AmOfFaces <= 8) ? AmOfFaces : 8;

	for(int i = 0; i < nFaces; i++)
	{
		// Use const reference to avoid radTHandle copy (which modifies reference count)
		// This is critical for OpenMP thread safety
		const radTHandlePgnAndTrans& hpt = VectHandlePgnAndTrans[i];
		radTPolygon* pgn = hpt.PgnHndl.rep;
		radTrans* tr = hpt.TransHndl.rep;

		// Get 2D vertices from polygon and transform to global 3D
		const radTVect2dVect& verts2d = pgn->EdgePointsVector;
		int nv = (int)verts2d.size();
		if(nv > 4) nv = 4;  // Limit to 4 vertices
		faceNumVerts[i] = nv;

		for(int j = 0; j < nv; j++)
		{
			faceVertices[i][j] = tr->TrPoint(TVector3d(verts2d[j].x, verts2d[j].y, pgn->CoordZ));
		}
	}

	if(FldKey.PreRelax_)
	{
		// =====================================================================
		// PreRelax mode: Compute interaction matrix coefficients
		// =====================================================================
		// The interaction matrix code calls B_comp with the element's Magn set to
		// a unit vector (1,0,0), (0,1,0), or (0,0,1), and expects the H field
		// to be returned in FieldPtr->H.
		//
		// We compute H from the current Magn and also store the result in B
		// for compatibility with both old and new matrix formats.

		TVector3d H_total(0., 0., 0.);

		for(int i = 0; i < nFaces; i++)
		{
			int nv = faceNumVerts[i];
			if(nv == 3)
			{
				// Triangular face: use directly
				const TVector3d& V0 = faceVertices[i][0];
				const TVector3d& V1 = faceVertices[i][1];
				const TVector3d& V2 = faceVertices[i][2];

				H_total += RadFieldFromTriangleFaceGlobal(V0, V1, V2, Magn, obsPoint, CentrPoint);
			}
			else if(nv == 4)
			{
				// Quadrilateral face: split into 2 triangles (V0,V1,V2) and (V0,V2,V3)
				const TVector3d& V0 = faceVertices[i][0];
				const TVector3d& V1 = faceVertices[i][1];
				const TVector3d& V2 = faceVertices[i][2];
				const TVector3d& V3 = faceVertices[i][3];

				// Triangle 1: V0, V1, V2
				H_total += RadFieldFromTriangleFaceGlobal(V0, V1, V2, Magn, obsPoint, CentrPoint);

				// Triangle 2: V0, V2, V3
				H_total += RadFieldFromTriangleFaceGlobal(V0, V2, V3, Magn, obsPoint, CentrPoint);
			}
		}

		// Store in both H (for new interaction matrix code) and B (for old matrix format)
		FieldPtr->H += H_total;
		FieldPtr->B += H_total;
	}
	else
	{
		// =====================================================================
		// Normal mode: Compute actual H field from current magnetization
		// =====================================================================
		TVector3d H_total(0., 0., 0.);

		for(int i = 0; i < nFaces; i++)
		{
			int nv = faceNumVerts[i];
			if(nv == 3)
			{
				// Triangular face
				const TVector3d& V0 = faceVertices[i][0];
				const TVector3d& V1 = faceVertices[i][1];
				const TVector3d& V2 = faceVertices[i][2];
				H_total += RadFieldFromTriangleFaceGlobal(V0, V1, V2, Magn, obsPoint, CentrPoint);
			}
			else if(nv == 4)
			{
				// Quadrilateral face: split into 2 triangles
				const TVector3d& V0 = faceVertices[i][0];
				const TVector3d& V1 = faceVertices[i][1];
				const TVector3d& V2 = faceVertices[i][2];
				const TVector3d& V3 = faceVertices[i][3];

				H_total += RadFieldFromTriangleFaceGlobal(V0, V1, V2, Magn, obsPoint, CentrPoint);
				H_total += RadFieldFromTriangleFaceGlobal(V0, V2, V3, Magn, obsPoint, CentrPoint);
			}
		}

		if(FldKey.H_) FieldPtr->H += H_total;

		// For B field, check if observation point is inside the wedge
		if(FldKey.B_)
		{
			FieldPtr->B += H_total;

			// Use solid angle method to check if point is inside
			// Sum solid angles of all faces - should be 4*pi if inside
			double totalSolidAngle = 0.0;
			for(int i = 0; i < nFaces; i++)
			{
				int nv = faceNumVerts[i];
				if(nv == 3)
				{
					const TVector3d& V0 = faceVertices[i][0];
					const TVector3d& V1 = faceVertices[i][1];
					const TVector3d& V2 = faceVertices[i][2];
					totalSolidAngle += RadPointClassify::ComputeTriangleSolidAngle(obsPoint, V0, V1, V2);
				}
				else if(nv == 4)
				{
					const TVector3d& V0 = faceVertices[i][0];
					const TVector3d& V1 = faceVertices[i][1];
					const TVector3d& V2 = faceVertices[i][2];
					const TVector3d& V3 = faceVertices[i][3];
					totalSolidAngle += RadPointClassify::ComputeTriangleSolidAngle(obsPoint, V0, V1, V2);
					totalSolidAngle += RadPointClassify::ComputeTriangleSolidAngle(obsPoint, V0, V2, V3);
				}
			}

			const double FOUR_PI = 4.0 * RadConst::PI;
			if(std::abs(totalSolidAngle - FOUR_PI) < 0.1)  // Inside check with tolerance
			{
				// B = mu_0 * (H + M), H_total is H, so add M to B
				FieldPtr->B += Magn;
			}
		}

		// M field (magnetization at observation point)
		if(FldKey.M_)
		{
			// Same solid angle check as B field
			double totalSolidAngle = 0.0;
			for(int i = 0; i < nFaces; i++)
			{
				int nv = faceNumVerts[i];
				if(nv == 3)
				{
					const TVector3d& V0 = faceVertices[i][0];
					const TVector3d& V1 = faceVertices[i][1];
					const TVector3d& V2 = faceVertices[i][2];
					totalSolidAngle += RadPointClassify::ComputeTriangleSolidAngle(obsPoint, V0, V1, V2);
				}
				else if(nv == 4)
				{
					const TVector3d& V0 = faceVertices[i][0];
					const TVector3d& V1 = faceVertices[i][1];
					const TVector3d& V2 = faceVertices[i][2];
					const TVector3d& V3 = faceVertices[i][3];
					totalSolidAngle += RadPointClassify::ComputeTriangleSolidAngle(obsPoint, V0, V1, V2);
					totalSolidAngle += RadPointClassify::ComputeTriangleSolidAngle(obsPoint, V0, V2, V3);
				}
			}

			const double FOUR_PI = 4.0 * RadConst::PI;
			if(std::abs(totalSolidAngle - FOUR_PI) < 0.1)
			{
				FieldPtr->M += Magn;
			}
		}

		// =====================================================================
		// IMA (Image Method) Field Contributions for wedge elements (analytical 3DOF)
		// Same pattern as tetrahedral IMA, but handles both tri and quad faces
		// =====================================================================
		if(RadIMAFieldContext::IsActive() && !FldKey.PreRelax_)
		{
			int imaSym = RadIMAFieldContext::GetSymmetry();
			int signX = RadIMAFieldContext::GetSignX();
			int signY = RadIMAFieldContext::GetSignY();
			int signZ = RadIMAFieldContext::GetSignZ();

			auto computeMirroredField = [&](int mirrorAxis, int sign) -> TVector3d {
				TVector3d H_mirror(0., 0., 0.);

				// Create mirrored face vertices
				std::array<std::array<TVector3d, 4>, 8> mirrorVerts;
				std::array<int, 8> mirrorNumVerts;
				for(int i = 0; i < nFaces; i++)
				{
					mirrorNumVerts[i] = faceNumVerts[i];
					int nv = faceNumVerts[i];
					for(int j = 0; j < nv; j++)
					{
						mirrorVerts[i][j] = faceVertices[i][j];
						if(mirrorAxis & radTInteraction::IMA_X) mirrorVerts[i][j].x = -mirrorVerts[i][j].x;
						if(mirrorAxis & radTInteraction::IMA_Y) mirrorVerts[i][j].y = -mirrorVerts[i][j].y;
						if(mirrorAxis & radTInteraction::IMA_Z) mirrorVerts[i][j].z = -mirrorVerts[i][j].z;
					}
				}

				// Count mirror axes: odd number requires winding reversal
				int numAxes = 0;
				if(mirrorAxis & radTInteraction::IMA_X) numAxes++;
				if(mirrorAxis & radTInteraction::IMA_Y) numAxes++;
				if(mirrorAxis & radTInteraction::IMA_Z) numAxes++;
				bool reverseWinding = (numAxes % 2 == 1);

				// Mirrored center point
				TVector3d mirrorCenter = CentrPoint;
				if(mirrorAxis & radTInteraction::IMA_X) mirrorCenter.x = -mirrorCenter.x;
				if(mirrorAxis & radTInteraction::IMA_Y) mirrorCenter.y = -mirrorCenter.y;
				if(mirrorAxis & radTInteraction::IMA_Z) mirrorCenter.z = -mirrorCenter.z;

				// Mirror magnetization: pseudo-vector flip + BC sign
				TVector3d mirrorMagn = Magn;
				if(mirrorAxis & radTInteraction::IMA_X) mirrorMagn.x = -mirrorMagn.x;
				if(mirrorAxis & radTInteraction::IMA_Y) mirrorMagn.y = -mirrorMagn.y;
				if(mirrorAxis & radTInteraction::IMA_Z) mirrorMagn.z = -mirrorMagn.z;
				mirrorMagn = mirrorMagn * (double)sign;

				// Compute field from mirrored faces
				for(int i = 0; i < nFaces; i++)
				{
					int nv = mirrorNumVerts[i];
					if(nv == 3)
					{
						const TVector3d& MV0 = mirrorVerts[i][0];
						const TVector3d& MV1 = mirrorVerts[i][1];
						const TVector3d& MV2 = mirrorVerts[i][2];
						if(reverseWinding)
							H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV2, MV1, mirrorMagn, obsPoint, mirrorCenter);
						else
							H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV1, MV2, mirrorMagn, obsPoint, mirrorCenter);
					}
					else if(nv == 4)
					{
						const TVector3d& MV0 = mirrorVerts[i][0];
						const TVector3d& MV1 = mirrorVerts[i][1];
						const TVector3d& MV2 = mirrorVerts[i][2];
						const TVector3d& MV3 = mirrorVerts[i][3];
						if(reverseWinding)
						{
							H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV3, MV2, mirrorMagn, obsPoint, mirrorCenter);
							H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV2, MV1, mirrorMagn, obsPoint, mirrorCenter);
						}
						else
						{
							H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV1, MV2, mirrorMagn, obsPoint, mirrorCenter);
							H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV2, MV3, mirrorMagn, obsPoint, mirrorCenter);
						}
					}
				}

				return H_mirror;
			};

			TVector3d H_ima(0., 0., 0.);

			// Single axis contributions
			if(imaSym & radTInteraction::IMA_X)
				H_ima += computeMirroredField(radTInteraction::IMA_X, signX);
			if(imaSym & radTInteraction::IMA_Y)
				H_ima += computeMirroredField(radTInteraction::IMA_Y, signY);
			if(imaSym & radTInteraction::IMA_Z)
				H_ima += computeMirroredField(radTInteraction::IMA_Z, signZ);

			// Dual axis contributions
			if((imaSym & radTInteraction::IMA_X) && (imaSym & radTInteraction::IMA_Y))
				H_ima += computeMirroredField(radTInteraction::IMA_XY, signX * signY);
			if((imaSym & radTInteraction::IMA_X) && (imaSym & radTInteraction::IMA_Z))
				H_ima += computeMirroredField(radTInteraction::IMA_XZ, signX * signZ);
			if((imaSym & radTInteraction::IMA_Y) && (imaSym & radTInteraction::IMA_Z))
				H_ima += computeMirroredField(radTInteraction::IMA_YZ, signY * signZ);

			// Triple axis contribution
			if((imaSym & radTInteraction::IMA_X) && (imaSym & radTInteraction::IMA_Y) && (imaSym & radTInteraction::IMA_Z))
				H_ima += computeMirroredField(radTInteraction::IMA_XYZ, signX * signY * signZ);

			// Add IMA contributions
			if(FldKey.H_) FieldPtr->H += H_ima;
			if(FldKey.B_) FieldPtr->B += H_ima;
		}
	}
}

//-------------------------------------------------------------------------

void radTPolyhedron::B_comp_hexahedron_MSC(radTField* FieldPtr)
{
	// =========================================================================
	// GLOBAL COORDINATE METHOD for hexahedral elements (6 quadrilateral faces)
	// =========================================================================
	// This method computes fields using GLOBAL coordinates only.
	// Each quadrilateral face is split into 2 triangles for computation.
	//
	// For each quadrilateral face with vertices [V0, V1, V2, V3]:
	//   - Triangle 1: V0, V1, V2
	//   - Triangle 2: V0, V2, V3
	// This is the standard diagonal split for quadrilateral faces.
	//
	// For PreRelax mode (interaction matrix construction):
	//   - Compute dH/dM for all 3 unit magnetization directions
	//
	// For Normal mode:
	//   - Compute H field from current magnetization

	radTFieldKey& FldKey = FieldPtr->FieldKey;
	TVector3d& obsPoint = FieldPtr->P;

	// Get face vertices in GLOBAL coordinates
	// For a hexahedron, we have 6 quadrilateral faces (4 vertices each)
	// Use fixed-size array to avoid heap allocation (better for OpenMP)
	std::array<std::array<TVector3d, 4>, 8> faceVertices;  // Max 8 faces
	int nFaces = (AmOfFaces <= 8) ? AmOfFaces : 8;
	for(int i = 0; i < nFaces; i++)
	{
		// Use const reference to avoid radTHandle copy (which modifies reference count)
		// This is critical for OpenMP thread safety
		const radTHandlePgnAndTrans& hpt = VectHandlePgnAndTrans[i];
		radTPolygon* pgn = hpt.PgnHndl.rep;
		radTrans* tr = hpt.TransHndl.rep;

		// Get 2D vertices from polygon and transform to global 3D
		const radTVect2dVect& verts2d = pgn->EdgePointsVector;
		if(verts2d.size() < 4) continue;  // Need 4 vertices for quad face

		// Transform each vertex: local 2D (x, y, CoordZ) -> global 3D
		faceVertices[i][0] = tr->TrPoint(TVector3d(verts2d[0].x, verts2d[0].y, pgn->CoordZ));
		faceVertices[i][1] = tr->TrPoint(TVector3d(verts2d[1].x, verts2d[1].y, pgn->CoordZ));
		faceVertices[i][2] = tr->TrPoint(TVector3d(verts2d[2].x, verts2d[2].y, pgn->CoordZ));
		faceVertices[i][3] = tr->TrPoint(TVector3d(verts2d[3].x, verts2d[3].y, pgn->CoordZ));
	}

	if(FldKey.PreRelax_)
	{
		// =====================================================================
		// PreRelax mode: Compute interaction matrix coefficients
		// =====================================================================
		TVector3d unit_Mx(1., 0., 0.);
		TVector3d unit_My(0., 1., 0.);
		TVector3d unit_Mz(0., 0., 1.);

		TVector3d H_from_Mx(0., 0., 0.);
		TVector3d H_from_My(0., 0., 0.);
		TVector3d H_from_Mz(0., 0., 0.);

		// Sum contributions from all faces (split each quad into 2 triangles)
		for(int i = 0; i < nFaces; i++)
		{
			const TVector3d& V0 = faceVertices[i][0];
			const TVector3d& V1 = faceVertices[i][1];
			const TVector3d& V2 = faceVertices[i][2];
			const TVector3d& V3 = faceVertices[i][3];

			// Triangle 1: V0, V1, V2
			H_from_Mx += RadFieldFromTriangleFaceGlobal(V0, V1, V2, unit_Mx, obsPoint, CentrPoint);
			H_from_My += RadFieldFromTriangleFaceGlobal(V0, V1, V2, unit_My, obsPoint, CentrPoint);
			H_from_Mz += RadFieldFromTriangleFaceGlobal(V0, V1, V2, unit_Mz, obsPoint, CentrPoint);

			// Triangle 2: V0, V2, V3
			H_from_Mx += RadFieldFromTriangleFaceGlobal(V0, V2, V3, unit_Mx, obsPoint, CentrPoint);
			H_from_My += RadFieldFromTriangleFaceGlobal(V0, V2, V3, unit_My, obsPoint, CentrPoint);
			H_from_Mz += RadFieldFromTriangleFaceGlobal(V0, V2, V3, unit_Mz, obsPoint, CentrPoint);
		}

		// Store in the matrix format expected by Radia
		FieldPtr->B.x += H_from_Mx.x;
		FieldPtr->B.y += H_from_Mx.y;
		FieldPtr->B.z += H_from_Mx.z;
		FieldPtr->H.x += H_from_My.x;
		FieldPtr->H.y += H_from_My.y;
		FieldPtr->H.z += H_from_My.z;
		FieldPtr->A.x += H_from_Mz.x;
		FieldPtr->A.y += H_from_Mz.y;
		FieldPtr->A.z += H_from_Mz.z;
	}
	else
	{
		// =====================================================================
		// Normal mode: Compute actual H field
		// =====================================================================
		// For permanent magnets (no material, fixed M): compute sigma = M · n
		// For soft materials (after Solve): use solved Sigma[i] values
		// =====================================================================
		TVector3d H_total(0., 0., 0.);

		// Check if Sigma values have been set (by Solve or directly)
		// If all Sigma = 0 and Magn != 0, this is a permanent magnet - compute sigma from M
		bool sigmaIsZero = true;
		for(int i = 0; i < nFaces && sigmaIsZero; i++) {
			if(Sigma[i] != 0.0) sigmaIsZero = false;
		}
		bool magnIsNotZero = (Magn.x != 0.0 || Magn.y != 0.0 || Magn.z != 0.0);

		if(sigmaIsZero && magnIsNotZero)
		{
			// Permanent magnet: compute sigma = M · n for each face
			// Use triangle-based field computation (like tetrahedron)
			for(int i = 0; i < nFaces; i++)
			{
				const TVector3d& V0 = faceVertices[i][0];
				const TVector3d& V1 = faceVertices[i][1];
				const TVector3d& V2 = faceVertices[i][2];
				const TVector3d& V3 = faceVertices[i][3];

				// Split quad into 2 triangles and use triangle field formula
				// Triangle 1: V0, V1, V2
				H_total += RadFieldFromTriangleFaceGlobal(V0, V1, V2, Magn, obsPoint, CentrPoint);
				// Triangle 2: V0, V2, V3
				H_total += RadFieldFromTriangleFaceGlobal(V0, V2, V3, Magn, obsPoint, CentrPoint);
			}
		}
		else
		{
			// Soft material (solved): use Sigma values directly
			// The solver uses A = -K/(4pi) - 1/chi * I (ELF-compatible), giving sigma with correct sign
			// Field formula: H = 2 * (sigma / 4pi) * solid_angle_integral
			//
			// Factor of 2 explanation (verified against ELF full model):
			// The MSC sigma represents sum of charges on both sides of each face.
			// ELF uses this convention for efficiency in the matrix computation.
			// When computing field, we need 2x to account for both charge sheets.
			for(int i = 0; i < nFaces; i++)
			{
				TVector3d H_face = FieldFromQuadFace(obsPoint, i, Sigma[i]);
				H_total.x += H_face.x;
				H_total.y += H_face.y;
				H_total.z += H_face.z;
			}

			// Apply 1/(4*pi) factor - consistent with matrix construction (ELF uses same factor for both)
			H_total.x *= RadConst::INV_FOUR_PI;
			H_total.y *= RadConst::INV_FOUR_PI;
			H_total.z *= RadConst::INV_FOUR_PI;
		}

		if(FldKey.H_) FieldPtr->H += H_total;

		// For B field, check if observation point is inside the hexahedron
		// If inside, B = mu_0 * (H + M), so we add M to B
		if(FldKey.B_)
		{
			FieldPtr->B += H_total;

			// Check if observation point is inside the hexahedron using solid angle
			// Sum solid angles from all face triangles (each quad split into 2 triangles)
			// Inside: |total_solid_angle| = 4*pi
			// Outside: |total_solid_angle| = 0
			const double PI = 3.14159265358979323846;
			const double FOUR_PI = 4.0 * PI;
			const double tolerance = 0.1;

			double total_solid_angle = 0.0;
			for(int i = 0; i < nFaces; i++)
			{
				const TVector3d& V0 = faceVertices[i][0];
				const TVector3d& V1 = faceVertices[i][1];
				const TVector3d& V2 = faceVertices[i][2];
				const TVector3d& V3 = faceVertices[i][3];

				// Triangle 1: V0, V1, V2
				total_solid_angle += RadPointClassify::ComputeTriangleSolidAngle(obsPoint, V0, V1, V2);
				// Triangle 2: V0, V2, V3
				total_solid_angle += RadPointClassify::ComputeTriangleSolidAngle(obsPoint, V0, V2, V3);
			}

			// Check for both winding conventions (inward or outward normals)
			bool pointInside = std::fabs(std::fabs(total_solid_angle) - FOUR_PI) < FOUR_PI * tolerance;
			if(pointInside)
			{
				// B = mu_0 * (H + M), H_total is H, so add M to B
				FieldPtr->B += Magn;
			}
		}

		// M field (magnetization at observation point)
		if(FldKey.M_)
		{
			// Check if observation point is inside the hexahedron using solid angle
			const double PI = 3.14159265358979323846;
			const double FOUR_PI = 4.0 * PI;
			const double tolerance = 0.1;

			double total_solid_angle = 0.0;
			for(int i = 0; i < nFaces; i++)
			{
				const TVector3d& V0 = faceVertices[i][0];
				const TVector3d& V1 = faceVertices[i][1];
				const TVector3d& V2 = faceVertices[i][2];
				const TVector3d& V3 = faceVertices[i][3];

				total_solid_angle += RadPointClassify::ComputeTriangleSolidAngle(obsPoint, V0, V1, V2);
				total_solid_angle += RadPointClassify::ComputeTriangleSolidAngle(obsPoint, V0, V2, V3);
			}

			bool pointInside = std::fabs(std::fabs(total_solid_angle) - FOUR_PI) < FOUR_PI * tolerance;
			if(pointInside)
			{
				FieldPtr->M += Magn;
			}
		}

		// =====================================================================
		// Scalar potential (phi) and vector potential (A) computation
		// =====================================================================
		// Uses FACE INTEGRATION (not dipole approximation)
		//
		// Vector potential: A(r) = (mu_0/4pi) * integral_S (M x n) / |r-r'| dS
		//                   B = curl(A)
		//
		// Each quadrilateral face is split into 2 triangles for integration.
		//
		// NOTE: On symmetry axes, the face-based integration may give A=0 due to
		// symmetric cancellation. This differs from ObjRecMag which uses an
		// 8-corner BufVect formula that doesn't cancel. For rectangular blocks,
		// use ObjRecMag instead for consistent A field on symmetry axes.
		// =====================================================================
		if(FldKey.A_)
		{
			TVector3d A_total(0., 0., 0.);

			for(int i = 0; i < nFaces; i++)
			{
				const TVector3d& V0 = faceVertices[i][0];
				const TVector3d& V1 = faceVertices[i][1];
				const TVector3d& V2 = faceVertices[i][2];
				const TVector3d& V3 = faceVertices[i][3];

				// Triangle 1: V0, V1, V2
				A_total += RadVectorPotentialFromTriangleFaceGlobal(
					V0, V1, V2, Magn, obsPoint, CentrPoint);
				// Triangle 2: V0, V2, V3
				A_total += RadVectorPotentialFromTriangleFaceGlobal(
					V0, V2, V3, Magn, obsPoint, CentrPoint);
			}

			FieldPtr->A += A_total;
		}

		// Scalar potential (phi) using face-based integration (no dipole approximation)
		// Phi = (1/4pi) * M dot BufVect, where BufVect = sum(n * integral(1/|r-r'|) dS)
		// Each quad face is split into 2 triangles for integration.
		if(FldKey.Phi_)
		{
			double Phi_total = 0.0;

			for(int i = 0; i < nFaces; i++)
			{
				const TVector3d& V0 = faceVertices[i][0];
				const TVector3d& V1 = faceVertices[i][1];
				const TVector3d& V2 = faceVertices[i][2];
				const TVector3d& V3 = faceVertices[i][3];

				// Triangle 1: V0, V1, V2
				Phi_total += RadScalarPotentialFromTriangleFaceGlobal(
					V0, V1, V2, Magn, obsPoint, CentrPoint);
				// Triangle 2: V0, V2, V3
				Phi_total += RadScalarPotentialFromTriangleFaceGlobal(
					V0, V2, V3, Magn, obsPoint, CentrPoint);
			}

			FieldPtr->Phi += Phi_total;
		}

		// =====================================================================
		// IMA (Image Method) Field Contributions
		// =====================================================================
		// When IMA is active, we need to add contributions from virtual mirror
		// elements. Each mirror element has:
		// 1. Mirrored geometry (negate x, y, or z coordinates)
		// 2. Permuted DOF values (face sigma values swapped)
		// 3. Sign based on symmetry type (+1 symmetric, -1 antisymmetric)
		// =====================================================================
		if(RadIMAFieldContext::IsActive() && !FldKey.PreRelax_)
		{
			int imaSym = RadIMAFieldContext::GetSymmetry();
			int signX = RadIMAFieldContext::GetSignX();
			int signY = RadIMAFieldContext::GetSignY();
			int signZ = RadIMAFieldContext::GetSignZ();

			// Helper lambda to compute field from mirrored geometry
			// The 'sign' parameter (+1 or -1) indicates the symmetry type:
			// +1 = symmetric BC (magnetization preserved)
			// -1 = antisymmetric BC (magnetization negated)
			auto computeMirroredField = [&](int mirrorAxis, int sign) -> TVector3d {
				TVector3d H_mirror(0., 0., 0.);

				// Create mirrored face vertices
				std::array<std::array<TVector3d, 4>, 8> mirrorVerts;
				for(int i = 0; i < nFaces; i++)
				{
					for(int j = 0; j < 4; j++)
					{
						mirrorVerts[i][j] = faceVertices[i][j];
						if(mirrorAxis & radTInteraction::IMA_X) mirrorVerts[i][j].x = -mirrorVerts[i][j].x;
						if(mirrorAxis & radTInteraction::IMA_Y) mirrorVerts[i][j].y = -mirrorVerts[i][j].y;
						if(mirrorAxis & radTInteraction::IMA_Z) mirrorVerts[i][j].z = -mirrorVerts[i][j].z;
					}
				}

				// Count number of mirror axes (odd number = need winding reversal)
				int numAxes = 0;
				if(mirrorAxis & radTInteraction::IMA_X) numAxes++;
				if(mirrorAxis & radTInteraction::IMA_Y) numAxes++;
				if(mirrorAxis & radTInteraction::IMA_Z) numAxes++;
				bool reverseWinding = (numAxes % 2 == 1);  // Odd number of reflections

				// Compute mirrored center point
				TVector3d mirrorCenter = CentrPoint;
				if(mirrorAxis & radTInteraction::IMA_X) mirrorCenter.x = -mirrorCenter.x;
				if(mirrorAxis & radTInteraction::IMA_Y) mirrorCenter.y = -mirrorCenter.y;
				if(mirrorAxis & radTInteraction::IMA_Z) mirrorCenter.z = -mirrorCenter.z;

				// Compute field from mirrored geometry
				if(sigmaIsZero && magnIsNotZero)
				{
					// Permanent magnet: use mirrored magnetization with sign
					TVector3d mirrorMagn = Magn;
					// For permanent magnets, the magnetization is a pseudo-vector
					// Mirror transformation: M -> M - 2*(M.n)*n for reflection across plane with normal n
					// For axis-aligned planes, this flips the perpendicular component
					if(mirrorAxis & radTInteraction::IMA_X) mirrorMagn.x = -mirrorMagn.x;
					if(mirrorAxis & radTInteraction::IMA_Y) mirrorMagn.y = -mirrorMagn.y;
					if(mirrorAxis & radTInteraction::IMA_Z) mirrorMagn.z = -mirrorMagn.z;
					// Apply BC sign: antisymmetric BC negates the entire magnetization
					mirrorMagn = mirrorMagn * (double)sign;

					for(int i = 0; i < nFaces; i++)
					{
						// Get vertex positions for this face
						const TVector3d& MV0 = mirrorVerts[i][0];
						const TVector3d& MV1 = mirrorVerts[i][1];
						const TVector3d& MV2 = mirrorVerts[i][2];
						const TVector3d& MV3 = mirrorVerts[i][3];

						if(reverseWinding)
						{
							// Reverse winding: use V0, V3, V2, V1 order
							// Triangle 1: (V0, V3, V2), Triangle 2: (V0, V2, V1)
							H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV3, MV2, mirrorMagn, obsPoint, mirrorCenter);
							H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV2, MV1, mirrorMagn, obsPoint, mirrorCenter);
						}
						else
						{
							// Keep original winding
							H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV1, MV2, mirrorMagn, obsPoint, mirrorCenter);
							H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV2, MV3, mirrorMagn, obsPoint, mirrorCenter);
						}
					}
				}
				else
				{
					// Soft material: compute mirror contributions for field
					for(int i = 0; i < nFaces; i++)
					{
						const TVector3d& MV0 = mirrorVerts[i][0];
						const TVector3d& MV1 = mirrorVerts[i][1];
						const TVector3d& MV2 = mirrorVerts[i][2];
						const TVector3d& MV3 = mirrorVerts[i][3];

						double mirrorSigma = sign * Sigma[i];
						TVector3d H_face = FieldFromQuadFaceMirrored(obsPoint, MV0, MV1, MV2, MV3, mirrorSigma, reverseWinding, mirrorCenter);

						H_mirror.x += H_face.x;
						H_mirror.y += H_face.y;
						H_mirror.z += H_face.z;
					}

					// Apply 1/(4*pi) factor - same as main field computation
					H_mirror.x *= RadConst::INV_FOUR_PI;
					H_mirror.y *= RadConst::INV_FOUR_PI;
					H_mirror.z *= RadConst::INV_FOUR_PI;
				}

				return H_mirror;  // No additional sign multiplication
			};

			TVector3d H_ima(0., 0., 0.);

			// Single axis contributions
			if(imaSym & radTInteraction::IMA_X)
				H_ima += computeMirroredField(radTInteraction::IMA_X, signX);
			if(imaSym & radTInteraction::IMA_Y)
				H_ima += computeMirroredField(radTInteraction::IMA_Y, signY);
			if(imaSym & radTInteraction::IMA_Z)
				H_ima += computeMirroredField(radTInteraction::IMA_Z, signZ);

			// Dual axis contributions
			if((imaSym & radTInteraction::IMA_X) && (imaSym & radTInteraction::IMA_Y))
				H_ima += computeMirroredField(radTInteraction::IMA_XY, signX * signY);
			if((imaSym & radTInteraction::IMA_X) && (imaSym & radTInteraction::IMA_Z))
				H_ima += computeMirroredField(radTInteraction::IMA_XZ, signX * signZ);
			if((imaSym & radTInteraction::IMA_Y) && (imaSym & radTInteraction::IMA_Z))
				H_ima += computeMirroredField(radTInteraction::IMA_YZ, signY * signZ);

			// Triple axis contribution
			if((imaSym & radTInteraction::IMA_X) && (imaSym & radTInteraction::IMA_Y) && (imaSym & radTInteraction::IMA_Z))
				H_ima += computeMirroredField(radTInteraction::IMA_XYZ, signX * signY * signZ);

			// Add IMA contributions
			if(FldKey.H_) FieldPtr->H += H_ima;
			if(FldKey.B_) FieldPtr->B += H_ima;
		}
	}
}

//-------------------------------------------------------------------------
// B_comp_wedge_MSC: 5-DOF MSC method for wedge elements (2 tri + 3 quad faces)
// Each face has one surface charge sigma as DOF (total 5 DOF)
//-------------------------------------------------------------------------
void radTPolyhedron::B_comp_wedge_MSC(radTField* FieldPtr)
{
	radTFieldKey& FldKey = FieldPtr->FieldKey;
	TVector3d& obsPoint = FieldPtr->P;

	// Get face vertices in GLOBAL coordinates
	// Wedges have 5 faces with variable vertex count (3 or 4)
	std::array<int, 8> faceNumVerts;
	std::array<std::array<TVector3d, 4>, 8> faceVertices;
	int nFaces = (AmOfFaces <= 8) ? AmOfFaces : 8;

	for(int i = 0; i < nFaces; i++)
	{
		const radTHandlePgnAndTrans& hpt = VectHandlePgnAndTrans[i];
		radTPolygon* pgn = hpt.PgnHndl.rep;
		radTrans* tr = hpt.TransHndl.rep;

		const radTVect2dVect& verts2d = pgn->EdgePointsVector;
		int nv = (int)verts2d.size();
		if(nv > 4) nv = 4;
		faceNumVerts[i] = nv;

		for(int j = 0; j < nv; j++)
		{
			faceVertices[i][j] = tr->TrPoint(TVector3d(verts2d[j].x, verts2d[j].y, pgn->CoordZ));
		}
	}

	if(FldKey.PreRelax_)
	{
		// PreRelax mode: Compute N-matrix (dH/dM) for cross-type interaction blocks
		// Same computation as B_comp_wedge_analytical PreRelax
		TVector3d H_total(0., 0., 0.);

		for(int i = 0; i < nFaces; i++)
		{
			int nv = faceNumVerts[i];
			if(nv == 3)
			{
				const TVector3d& V0 = faceVertices[i][0];
				const TVector3d& V1 = faceVertices[i][1];
				const TVector3d& V2 = faceVertices[i][2];
				H_total += RadFieldFromTriangleFaceGlobal(V0, V1, V2, Magn, obsPoint, CentrPoint);
			}
			else if(nv == 4)
			{
				const TVector3d& V0 = faceVertices[i][0];
				const TVector3d& V1 = faceVertices[i][1];
				const TVector3d& V2 = faceVertices[i][2];
				const TVector3d& V3 = faceVertices[i][3];

				H_total += RadFieldFromTriangleFaceGlobal(V0, V1, V2, Magn, obsPoint, CentrPoint);
				H_total += RadFieldFromTriangleFaceGlobal(V0, V2, V3, Magn, obsPoint, CentrPoint);
			}
		}

		FieldPtr->H += H_total;
		FieldPtr->B += H_total;
	}
	else
	{
		// Normal mode: Compute field from sigma values or magnetization
		TVector3d H_total(0., 0., 0.);

		// Check if Sigma values have been set (by Solve or directly)
		bool sigmaIsZero = true;
		for(int i = 0; i < nFaces && sigmaIsZero; i++) {
			if(Sigma[i] != 0.0) sigmaIsZero = false;
		}
		bool magnIsNotZero = (Magn.x != 0.0 || Magn.y != 0.0 || Magn.z != 0.0);

		if(sigmaIsZero && magnIsNotZero)
		{
			// Permanent magnet: compute field from M using triangle faces
			for(int i = 0; i < nFaces; i++)
			{
				int nv = faceNumVerts[i];
				if(nv == 3)
				{
					const TVector3d& V0 = faceVertices[i][0];
					const TVector3d& V1 = faceVertices[i][1];
					const TVector3d& V2 = faceVertices[i][2];
					H_total += RadFieldFromTriangleFaceGlobal(V0, V1, V2, Magn, obsPoint, CentrPoint);
				}
				else if(nv == 4)
				{
					const TVector3d& V0 = faceVertices[i][0];
					const TVector3d& V1 = faceVertices[i][1];
					const TVector3d& V2 = faceVertices[i][2];
					const TVector3d& V3 = faceVertices[i][3];
					H_total += RadFieldFromTriangleFaceGlobal(V0, V1, V2, Magn, obsPoint, CentrPoint);
					H_total += RadFieldFromTriangleFaceGlobal(V0, V2, V3, Magn, obsPoint, CentrPoint);
				}
			}
		}
		else
		{
			// Soft material (solved): use Sigma values via FieldFromFace
			for(int i = 0; i < nFaces; i++)
			{
				TVector3d H_face = FieldFromFace(obsPoint, i, Sigma[i]);
				H_total.x += H_face.x;
				H_total.y += H_face.y;
				H_total.z += H_face.z;
			}

			// Apply 1/(4*pi) factor - consistent with matrix construction
			H_total.x *= RadConst::INV_FOUR_PI;
			H_total.y *= RadConst::INV_FOUR_PI;
			H_total.z *= RadConst::INV_FOUR_PI;
		}

		if(FldKey.H_) FieldPtr->H += H_total;

		// B field with inside check using solid angle
		if(FldKey.B_)
		{
			FieldPtr->B += H_total;

			double totalSolidAngle = 0.0;
			for(int i = 0; i < nFaces; i++)
			{
				int nv = faceNumVerts[i];
				if(nv == 3)
				{
					totalSolidAngle += RadPointClassify::ComputeTriangleSolidAngle(
						obsPoint, faceVertices[i][0], faceVertices[i][1], faceVertices[i][2]);
				}
				else if(nv == 4)
				{
					totalSolidAngle += RadPointClassify::ComputeTriangleSolidAngle(
						obsPoint, faceVertices[i][0], faceVertices[i][1], faceVertices[i][2]);
					totalSolidAngle += RadPointClassify::ComputeTriangleSolidAngle(
						obsPoint, faceVertices[i][0], faceVertices[i][2], faceVertices[i][3]);
				}
			}

			const double FOUR_PI = 4.0 * RadConst::PI;
			if(std::abs(totalSolidAngle - FOUR_PI) < 0.1)
			{
				FieldPtr->B += Magn;
			}
		}

		// M field
		if(FldKey.M_)
		{
			double totalSolidAngle = 0.0;
			for(int i = 0; i < nFaces; i++)
			{
				int nv = faceNumVerts[i];
				if(nv == 3)
				{
					totalSolidAngle += RadPointClassify::ComputeTriangleSolidAngle(
						obsPoint, faceVertices[i][0], faceVertices[i][1], faceVertices[i][2]);
				}
				else if(nv == 4)
				{
					totalSolidAngle += RadPointClassify::ComputeTriangleSolidAngle(
						obsPoint, faceVertices[i][0], faceVertices[i][1], faceVertices[i][2]);
					totalSolidAngle += RadPointClassify::ComputeTriangleSolidAngle(
						obsPoint, faceVertices[i][0], faceVertices[i][2], faceVertices[i][3]);
				}
			}

			const double FOUR_PI = 4.0 * RadConst::PI;
			if(std::abs(totalSolidAngle - FOUR_PI) < 0.1)
			{
				FieldPtr->M += Magn;
			}
		}

		// A field (vector potential)
		if(FldKey.A_)
		{
			TVector3d A_total(0., 0., 0.);
			for(int i = 0; i < nFaces; i++)
			{
				int nv = faceNumVerts[i];
				if(nv == 3)
				{
					const TVector3d& V0 = faceVertices[i][0];
					const TVector3d& V1 = faceVertices[i][1];
					const TVector3d& V2 = faceVertices[i][2];
					A_total += RadVectorPotentialFromTriangleFaceGlobal(
						V0, V1, V2, Magn, obsPoint, CentrPoint);
				}
				else if(nv == 4)
				{
					const TVector3d& V0 = faceVertices[i][0];
					const TVector3d& V1 = faceVertices[i][1];
					const TVector3d& V2 = faceVertices[i][2];
					const TVector3d& V3 = faceVertices[i][3];
					A_total += RadVectorPotentialFromTriangleFaceGlobal(
						V0, V1, V2, Magn, obsPoint, CentrPoint);
					A_total += RadVectorPotentialFromTriangleFaceGlobal(
						V0, V2, V3, Magn, obsPoint, CentrPoint);
				}
			}
			FieldPtr->A += A_total;
		}

		// Scalar potential (phi)
		if(FldKey.Phi_)
		{
			double Phi_total = 0.0;
			for(int i = 0; i < nFaces; i++)
			{
				int nv = faceNumVerts[i];
				if(nv == 3)
				{
					const TVector3d& V0 = faceVertices[i][0];
					const TVector3d& V1 = faceVertices[i][1];
					const TVector3d& V2 = faceVertices[i][2];
					Phi_total += RadScalarPotentialFromTriangleFaceGlobal(
						V0, V1, V2, Magn, obsPoint, CentrPoint);
				}
				else if(nv == 4)
				{
					const TVector3d& V0 = faceVertices[i][0];
					const TVector3d& V1 = faceVertices[i][1];
					const TVector3d& V2 = faceVertices[i][2];
					const TVector3d& V3 = faceVertices[i][3];
					Phi_total += RadScalarPotentialFromTriangleFaceGlobal(
						V0, V1, V2, Magn, obsPoint, CentrPoint);
					Phi_total += RadScalarPotentialFromTriangleFaceGlobal(
						V0, V2, V3, Magn, obsPoint, CentrPoint);
				}
			}
			FieldPtr->Phi += Phi_total;
		}

		// =====================================================================
		// IMA (Image Method) Field Contributions for wedge elements
		// Same pattern as B_comp_hexahedron_MSC IMA, adapted for variable face vertex count
		// =====================================================================
		if(RadIMAFieldContext::IsActive() && !FldKey.PreRelax_)
		{
			int imaSym = RadIMAFieldContext::GetSymmetry();
			int signX = RadIMAFieldContext::GetSignX();
			int signY = RadIMAFieldContext::GetSignY();
			int signZ = RadIMAFieldContext::GetSignZ();

			auto computeMirroredField = [&](int mirrorAxis, int sign) -> TVector3d {
				TVector3d H_mirror(0., 0., 0.);

				// Create mirrored face vertices (handle variable numVerts per face)
				std::array<std::array<TVector3d, 4>, 8> mirrorVerts;
				for(int i = 0; i < nFaces; i++)
				{
					for(int j = 0; j < faceNumVerts[i]; j++)
					{
						mirrorVerts[i][j] = faceVertices[i][j];
						if(mirrorAxis & radTInteraction::IMA_X) mirrorVerts[i][j].x = -mirrorVerts[i][j].x;
						if(mirrorAxis & radTInteraction::IMA_Y) mirrorVerts[i][j].y = -mirrorVerts[i][j].y;
						if(mirrorAxis & radTInteraction::IMA_Z) mirrorVerts[i][j].z = -mirrorVerts[i][j].z;
					}
				}

				// Count number of mirror axes (odd number = need winding reversal)
				int numAxes = 0;
				if(mirrorAxis & radTInteraction::IMA_X) numAxes++;
				if(mirrorAxis & radTInteraction::IMA_Y) numAxes++;
				if(mirrorAxis & radTInteraction::IMA_Z) numAxes++;
				bool reverseWinding = (numAxes % 2 == 1);

				// Compute mirrored center point
				TVector3d mirrorCenter = CentrPoint;
				if(mirrorAxis & radTInteraction::IMA_X) mirrorCenter.x = -mirrorCenter.x;
				if(mirrorAxis & radTInteraction::IMA_Y) mirrorCenter.y = -mirrorCenter.y;
				if(mirrorAxis & radTInteraction::IMA_Z) mirrorCenter.z = -mirrorCenter.z;

				if(sigmaIsZero && magnIsNotZero)
				{
					// Permanent magnet path
					TVector3d mirrorMagn = Magn;
					if(mirrorAxis & radTInteraction::IMA_X) mirrorMagn.x = -mirrorMagn.x;
					if(mirrorAxis & radTInteraction::IMA_Y) mirrorMagn.y = -mirrorMagn.y;
					if(mirrorAxis & radTInteraction::IMA_Z) mirrorMagn.z = -mirrorMagn.z;
					mirrorMagn = mirrorMagn * (double)sign;

					for(int i = 0; i < nFaces; i++)
					{
						int nv = faceNumVerts[i];
						const TVector3d& MV0 = mirrorVerts[i][0];
						const TVector3d& MV1 = mirrorVerts[i][1];
						const TVector3d& MV2 = mirrorVerts[i][2];

						if(nv == 3)
						{
							// Triangular face: single call
							if(reverseWinding)
								H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV2, MV1, mirrorMagn, obsPoint, mirrorCenter);
							else
								H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV1, MV2, mirrorMagn, obsPoint, mirrorCenter);
						}
						else if(nv == 4)
						{
							// Quad face: split into 2 triangles
							const TVector3d& MV3 = mirrorVerts[i][3];
							if(reverseWinding)
							{
								H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV3, MV2, mirrorMagn, obsPoint, mirrorCenter);
								H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV2, MV1, mirrorMagn, obsPoint, mirrorCenter);
							}
							else
							{
								H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV1, MV2, mirrorMagn, obsPoint, mirrorCenter);
								H_mirror += RadFieldFromTriangleFaceGlobal(MV0, MV2, MV3, mirrorMagn, obsPoint, mirrorCenter);
							}
						}
					}
				}
				else
				{
					// Soft material path: use FieldFromFaceMirrored
					for(int i = 0; i < nFaces; i++)
					{
						double mirrorSigma = sign * Sigma[i];
						TVector3d H_face = FieldFromFaceMirrored(obsPoint,
						                                          mirrorVerts[i].data(), faceNumVerts[i],
						                                          mirrorSigma, reverseWinding,
						                                          mirrorCenter);
						H_mirror.x += H_face.x;
						H_mirror.y += H_face.y;
						H_mirror.z += H_face.z;
					}

					// Apply 1/(4*pi) factor
					H_mirror.x *= RadConst::INV_FOUR_PI;
					H_mirror.y *= RadConst::INV_FOUR_PI;
					H_mirror.z *= RadConst::INV_FOUR_PI;
				}

				return H_mirror;
			};

			TVector3d H_ima(0., 0., 0.);

			// Single axis contributions
			if(imaSym & radTInteraction::IMA_X)
				H_ima += computeMirroredField(radTInteraction::IMA_X, signX);
			if(imaSym & radTInteraction::IMA_Y)
				H_ima += computeMirroredField(radTInteraction::IMA_Y, signY);
			if(imaSym & radTInteraction::IMA_Z)
				H_ima += computeMirroredField(radTInteraction::IMA_Z, signZ);

			// Dual axis contributions
			if((imaSym & radTInteraction::IMA_X) && (imaSym & radTInteraction::IMA_Y))
				H_ima += computeMirroredField(radTInteraction::IMA_XY, signX * signY);
			if((imaSym & radTInteraction::IMA_X) && (imaSym & radTInteraction::IMA_Z))
				H_ima += computeMirroredField(radTInteraction::IMA_XZ, signX * signZ);
			if((imaSym & radTInteraction::IMA_Y) && (imaSym & radTInteraction::IMA_Z))
				H_ima += computeMirroredField(radTInteraction::IMA_YZ, signY * signZ);

			// Triple axis contribution
			if((imaSym & radTInteraction::IMA_X) && (imaSym & radTInteraction::IMA_Y) && (imaSym & radTInteraction::IMA_Z))
				H_ima += computeMirroredField(radTInteraction::IMA_XYZ, signX * signY * signZ);

			// Add IMA contributions
			if(FldKey.H_) FieldPtr->H += H_ima;
			if(FldKey.B_) FieldPtr->B += H_ima;
		}
	}
}

//-------------------------------------------------------------------------
// FieldFromQuadFaceMirrored: Helper for IMA field computation
// Computes field from a quad face with explicit vertex positions
// Uses mirrorCenter to determine correct outward normal direction
// If flipNormal is true, swaps triangle vertex order to flip normal direction
//-------------------------------------------------------------------------
TVector3d radTPolyhedron::FieldFromQuadFaceMirrored(const TVector3d& obs,
                                                     const TVector3d& V0,
                                                     const TVector3d& V1,
                                                     const TVector3d& V2,
                                                     const TVector3d& V3,
                                                     double sigma,
                                                     bool flipNormal,
                                                     const TVector3d& mirrorCenter) const
{
	// Split quad into 2 triangles, applying sign_factor correction
	// based on whether triangle normal points outward from mirrorCenter.
	// This matches the logic in FieldFromQuadFace().

	TVector3d H_total(0.0, 0.0, 0.0);

	// Triangle split: (V0,V1,V2) and (V0,V2,V3) - or reversed if flipNormal
	TVector3d tri_verts[2][3];
	if(flipNormal)
	{
		// Flip normal by swapping second and third vertex
		tri_verts[0][0] = V0; tri_verts[0][1] = V2; tri_verts[0][2] = V1;
		tri_verts[1][0] = V0; tri_verts[1][1] = V3; tri_verts[1][2] = V2;
	}
	else
	{
		tri_verts[0][0] = V0; tri_verts[0][1] = V1; tri_verts[0][2] = V2;
		tri_verts[1][0] = V0; tri_verts[1][1] = V2; tri_verts[1][2] = V3;
	}

	for(int iTri = 0; iTri < 2; iTri++)
	{
		const TVector3d& T0 = tri_verts[iTri][0];
		const TVector3d& T1 = tri_verts[iTri][1];
		const TVector3d& T2 = tri_verts[iTri][2];

		// Compute triangle normal
		TVector3d edge1, edge2, tri_normal;
		edge1.x = T1.x - T0.x; edge1.y = T1.y - T0.y; edge1.z = T1.z - T0.z;
		edge2.x = T2.x - T0.x; edge2.y = T2.y - T0.y; edge2.z = T2.z - T0.z;

		tri_normal.x = edge1.y * edge2.z - edge1.z * edge2.y;
		tri_normal.y = edge1.z * edge2.x - edge1.x * edge2.z;
		tri_normal.z = edge1.x * edge2.y - edge1.y * edge2.x;

		double norm_len = sqrt(tri_normal.x*tri_normal.x + tri_normal.y*tri_normal.y + tri_normal.z*tri_normal.z);
		if(norm_len < 1e-20) continue;

		// Normalize
		tri_normal.x /= norm_len;
		tri_normal.y /= norm_len;
		tri_normal.z /= norm_len;

		// Compute triangle center
		TVector3d tri_center;
		tri_center.x = (T0.x + T1.x + T2.x) / 3.0;
		tri_center.y = (T0.y + T1.y + T2.y) / 3.0;
		tri_center.z = (T0.z + T1.z + T2.z) / 3.0;

		// Check if normal points outward (away from MIRROR element center)
		TVector3d to_center;
		to_center.x = tri_center.x - mirrorCenter.x;
		to_center.y = tri_center.y - mirrorCenter.y;
		to_center.z = tri_center.z - mirrorCenter.z;

		double dot_prod = tri_normal.x * to_center.x + tri_normal.y * to_center.y + tri_normal.z * to_center.z;

		// Sign factor: +1 if normal points outward, -1 if inward
		double sign_factor = (dot_prod >= 0.0) ? 1.0 : -1.0;

		// Compute field from this triangle with sign-corrected sigma
		TVector3d H_tri = FieldFromChargedTriangle(obs, T0, T1, T2, sigma * sign_factor);

		H_total.x += H_tri.x;
		H_total.y += H_tri.y;
		H_total.z += H_tri.z;
	}

	return H_total;
}

TVector3d radTPolyhedron::FieldFromFaceMirrored(const TVector3d& obs,
                                                 const TVector3d* mirrorVerts, int numVerts,
                                                 double sigma, bool flipNormal,
                                                 const TVector3d& mirrorCenter) const
{
	// Generalized mirrored face field: dispatches to tri or quad based on numVerts
	if(numVerts == 3)
	{
		// Triangular face: apply flipNormal by swapping V1/V2, then check outward normal
		TVector3d T0, T1, T2;
		if(flipNormal)
		{
			T0 = mirrorVerts[0]; T1 = mirrorVerts[2]; T2 = mirrorVerts[1];
		}
		else
		{
			T0 = mirrorVerts[0]; T1 = mirrorVerts[1]; T2 = mirrorVerts[2];
		}

		// Compute triangle normal from cross product
		TVector3d edge1, edge2, tri_normal;
		edge1.x = T1.x - T0.x; edge1.y = T1.y - T0.y; edge1.z = T1.z - T0.z;
		edge2.x = T2.x - T0.x; edge2.y = T2.y - T0.y; edge2.z = T2.z - T0.z;

		tri_normal.x = edge1.y * edge2.z - edge1.z * edge2.y;
		tri_normal.y = edge1.z * edge2.x - edge1.x * edge2.z;
		tri_normal.z = edge1.x * edge2.y - edge1.y * edge2.x;

		double norm_len = sqrt(tri_normal.x*tri_normal.x + tri_normal.y*tri_normal.y + tri_normal.z*tri_normal.z);
		if(norm_len < 1e-20) return TVector3d(0., 0., 0.);

		// Check if normal points outward (away from MIRROR element center)
		TVector3d tri_center;
		tri_center.x = (T0.x + T1.x + T2.x) / 3.0;
		tri_center.y = (T0.y + T1.y + T2.y) / 3.0;
		tri_center.z = (T0.z + T1.z + T2.z) / 3.0;

		TVector3d to_center;
		to_center.x = tri_center.x - mirrorCenter.x;
		to_center.y = tri_center.y - mirrorCenter.y;
		to_center.z = tri_center.z - mirrorCenter.z;

		double dot_prod = tri_normal.x * to_center.x + tri_normal.y * to_center.y + tri_normal.z * to_center.z;
		double sign_factor = (dot_prod >= 0.0) ? 1.0 : -1.0;

		return FieldFromChargedTriangle(obs, T0, T1, T2, sigma * sign_factor);
	}
	else if(numVerts >= 4)
	{
		// Quad face: delegate to FieldFromQuadFaceMirrored (already has sign_factor)
		return FieldFromQuadFaceMirrored(obs, mirrorVerts[0], mirrorVerts[1],
		                                  mirrorVerts[2], mirrorVerts[3],
		                                  sigma, flipNormal, mirrorCenter);
	}

	return TVector3d(0., 0., 0.);
}

//-------------------------------------------------------------------------
// 6 DOF MSC field computation methods for hexahedra
// Note: 1/(4*pi) factor is applied in matrix assembly (rad_interaction.cpp),
// not in these field computation functions. Use RadConst::INV_FOUR_PI.
//-------------------------------------------------------------------------

TVector3d radTPolyhedron::FieldFromChargedTriangle(const TVector3d& obs,
                                                    const TVector3d& v0,
                                                    const TVector3d& v1,
                                                    const TVector3d& v2,
                                                    double sigma) const
{
	// Analytic field from uniformly charged triangle
	// Implements BOTH tangential (log terms) AND normal (atan terms) components
	// Returns field WITHOUT 4pi divisor (4pi is applied in matrix assembly)
	// Matches ELF_MAGIC: m_legacy_surface.f90::calc_field_quad_tria

	const double EPS = 1.0e-20;

	// Build local coordinate system (ELF_MAGIC convention)
	// basis_a = edge (v1-v0) normalized (local X)
	// basis_c = face normal (local Z)
	// basis_b = basis_c x basis_a (local Y)

	TVector3d e1, e2;
	e1.x = v1.x - v0.x; e1.y = v1.y - v0.y; e1.z = v1.z - v0.z;
	e2.x = v2.x - v0.x; e2.y = v2.y - v0.y; e2.z = v2.z - v0.z;

	// Face normal = e1 x e2 (basis_c)
	TVector3d basis_c;
	basis_c.x = e1.y * e2.z - e1.z * e2.y;
	basis_c.y = e1.z * e2.x - e1.x * e2.z;
	basis_c.z = e1.x * e2.y - e1.y * e2.x;

	double cLen = sqrt(basis_c.x*basis_c.x + basis_c.y*basis_c.y + basis_c.z*basis_c.z);
	if(cLen < EPS) return TVector3d(0.0, 0.0, 0.0);
	basis_c.x /= cLen; basis_c.y /= cLen; basis_c.z /= cLen;

	// basis_a = e1 normalized (ELF: tri_verts(:,2) - tri_verts(:,1) = v1-v0 = e1)
	TVector3d basis_a = e1;
	double aLen = sqrt(basis_a.x*basis_a.x + basis_a.y*basis_a.y + basis_a.z*basis_a.z);
	if(aLen < EPS) return TVector3d(0.0, 0.0, 0.0);
	basis_a.x /= aLen; basis_a.y /= aLen; basis_a.z /= aLen;

	// basis_b = basis_c x basis_a (ELF: normalize_cross_product(basis_c, basis_a, basis_b))
	TVector3d basis_b;
	basis_b.x = basis_c.y * basis_a.z - basis_c.z * basis_a.y;
	basis_b.y = basis_c.z * basis_a.x - basis_c.x * basis_a.z;
	basis_b.z = basis_c.x * basis_a.y - basis_c.y * basis_a.x;
	double bLen = sqrt(basis_b.x*basis_b.x + basis_b.y*basis_b.y + basis_b.z*basis_b.z);
	if(bLen < EPS) return TVector3d(0.0, 0.0, 0.0);
	basis_b.x /= bLen; basis_b.y /= bLen; basis_b.z /= bLen;

	TVector3d AA = basis_a;  // Local X
	TVector3d BB = basis_b;  // Local Y

	// Convert vertices to local 2D coordinates (v0 = origin, ELF: face_origin = tri_verts(:,1))
	double xy0_x = 0.0, xy0_y = 0.0;

	double xy1_x = e1.x*AA.x + e1.y*AA.y + e1.z*AA.z;
	double xy1_y = e1.x*BB.x + e1.y*BB.y + e1.z*BB.z;

	double xy2_x = e2.x*AA.x + e2.y*AA.y + e2.z*AA.z;
	double xy2_y = e2.x*BB.x + e2.y*BB.y + e2.z*BB.z;

	// Edge parameters (3 edges for triangle)
	double XY[3][2] = {{xy0_x, xy0_y}, {xy1_x, xy1_y}, {xy2_x, xy2_y}};
	double DS[3], AM[3], SM[3], XD[3], YD[3];
	double EPSG = 0.0;

	for(int j = 0; j < 3; j++)
	{
		int l = (j + 1) % 3;
		double dx = XY[l][0] - XY[j][0];
		double dy = XY[l][1] - XY[j][1];
		if(fabs(dx) < EPS) dx = (dx >= 0) ? EPS : -EPS;

		DS[j] = sqrt(dx*dx + dy*dy);
		AM[j] = dy / dx;
		SM[j] = sqrt(AM[j]*AM[j] + 1.0);
		XD[j] = -dx / DS[j];
		YD[j] =  dy / DS[j];

		if(DS[j] > EPSG) EPSG = DS[j];
	}
	EPSG *= 1.0e-12;

	// Transform observation point to local coordinates
	TVector3d d;
	d.x = obs.x - v0.x;
	d.y = obs.y - v0.y;
	d.z = obs.z - v0.z;

	double EE1 = d.x*AA.x + d.y*AA.y + d.z*AA.z;  // local X
	double EE2 = d.x*BB.x + d.y*BB.y + d.z*BB.z;  // local Y
	double EE3 = d.x*basis_c.x + d.y*basis_c.y + d.z*basis_c.z;  // local Z (height)

	// Distances from observation point to vertices
	double X[3], Y[3], H[3], E[3], R[3];
	for(int j = 0; j < 3; j++)
	{
		X[j] = EE1 - XY[j][0];
		Y[j] = EE2 - XY[j][1];
		H[j] = Y[j] * X[j];
		E[j] = EE3*EE3 + X[j]*X[j];
		R[j] = sqrt(X[j]*X[j] + Y[j]*Y[j] + EE3*EE3);
	}

	double Z = EE3;

	// Edge contributions
	double RM[3], RP[3], RR[3], AL[3];
	for(int j = 0; j < 3; j++)
	{
		int jp1 = (j + 1) % 3;
		RM[j] = R[j] + R[jp1] - DS[j];
		RP[j] = R[j] + R[jp1] + DS[j];
		RR[j] = (RM[j] / RP[j] > EPS) ? (RM[j] / RP[j]) : EPS;
		AL[j] = log(RR[j]);
	}

	// Field components in local frame WITHOUT 4pi divisor
	// (matches ELF_MAGIC convention - 4pi is applied in matrix assembly)
	// Tangential components (log terms): H_tan = sigma * sum of log terms
	double HH1 = sigma * (-YD[0]*AL[0] - YD[1]*AL[1] - YD[2]*AL[2]);
	double HH2 = sigma * (-XD[0]*AL[0] - XD[1]*AL[1] - XD[2]*AL[2]);
	double HH3 = 0.0;

	// Normal component (atan terms) - only if not on surface
	if(fabs(Z) > EPSG)
	{
		double ZR[3];
		for(int j = 0; j < 3; j++)
		{
			ZR[j] = Z * R[j];
		}

		double AT[3], BT[3];
		for(int j = 0; j < 3; j++)
		{
			int jp1 = (j + 1) % 3;
			AT[j] = (AM[j]*E[j] - H[j]) / ZR[j];
			BT[j] = (AM[j]*E[jp1] - H[jp1]) / ZR[jp1];
		}

		HH3 = sigma * (-atan(AT[0]) - atan(AT[1]) - atan(AT[2])
		               +atan(BT[0]) + atan(BT[1]) + atan(BT[2]));
	}

	// Transform back to global coordinates
	TVector3d Hfield;
	Hfield.x = HH1*AA.x + HH2*BB.x + HH3*basis_c.x;
	Hfield.y = HH1*AA.y + HH2*BB.y + HH3*basis_c.y;
	Hfield.z = HH1*AA.z + HH2*BB.z + HH3*basis_c.z;

	return Hfield;
}

//-------------------------------------------------------------------------
// FieldFromChargedTriangleWithNormal: Version with explicit normal for IMA boundary faces
// When vertices don't move under mirroring (boundary faces), the computed normal from
// cross product would be wrong. This version uses an explicit transformed normal.
//-------------------------------------------------------------------------
TVector3d radTPolyhedron::FieldFromChargedTriangleWithNormal(const TVector3d& obs,
                                                              const TVector3d& v0,
                                                              const TVector3d& v1,
                                                              const TVector3d& v2,
                                                              double sigma,
                                                              const TVector3d& explicitNormal) const
{
	// Analytic field from uniformly charged triangle with explicit normal direction
	// Used for IMA boundary faces where the mirrored geometry doesn't give correct normal
	// Returns field WITHOUT 4pi divisor (4pi is applied in matrix assembly)

	const double EPS = 1.0e-20;

	// Use the provided explicit normal as basis_c (local Z axis)
	TVector3d basis_c = explicitNormal;
	double cLen = sqrt(basis_c.x*basis_c.x + basis_c.y*basis_c.y + basis_c.z*basis_c.z);
	if(cLen < EPS) return TVector3d(0.0, 0.0, 0.0);
	basis_c.x /= cLen; basis_c.y /= cLen; basis_c.z /= cLen;

	// Build local coordinate system using explicit normal
	// basis_a = edge (v1-v0) normalized (local X)
	TVector3d e1, e2;
	e1.x = v1.x - v0.x; e1.y = v1.y - v0.y; e1.z = v1.z - v0.z;
	e2.x = v2.x - v0.x; e2.y = v2.y - v0.y; e2.z = v2.z - v0.z;

	TVector3d basis_a = e1;
	double aLen = sqrt(basis_a.x*basis_a.x + basis_a.y*basis_a.y + basis_a.z*basis_a.z);
	if(aLen < EPS) return TVector3d(0.0, 0.0, 0.0);
	basis_a.x /= aLen; basis_a.y /= aLen; basis_a.z /= aLen;

	// basis_b = basis_c x basis_a (local Y)
	TVector3d basis_b;
	basis_b.x = basis_c.y * basis_a.z - basis_c.z * basis_a.y;
	basis_b.y = basis_c.z * basis_a.x - basis_c.x * basis_a.z;
	basis_b.z = basis_c.x * basis_a.y - basis_c.y * basis_a.x;
	double bLen = sqrt(basis_b.x*basis_b.x + basis_b.y*basis_b.y + basis_b.z*basis_b.z);
	if(bLen < EPS) return TVector3d(0.0, 0.0, 0.0);
	basis_b.x /= bLen; basis_b.y /= bLen; basis_b.z /= bLen;

	TVector3d AA = basis_a;  // Local X
	TVector3d BB = basis_b;  // Local Y

	// Convert vertices to local 2D coordinates (v0 = origin)
	double xy0_x = 0.0, xy0_y = 0.0;

	double xy1_x = e1.x*AA.x + e1.y*AA.y + e1.z*AA.z;
	double xy1_y = e1.x*BB.x + e1.y*BB.y + e1.z*BB.z;

	double xy2_x = e2.x*AA.x + e2.y*AA.y + e2.z*AA.z;
	double xy2_y = e2.x*BB.x + e2.y*BB.y + e2.z*BB.z;

	// Edge parameters (3 edges for triangle)
	double XY[3][2] = {{xy0_x, xy0_y}, {xy1_x, xy1_y}, {xy2_x, xy2_y}};
	double DS[3], AM[3], SM[3], XD[3], YD[3];
	double EPSG = 0.0;

	for(int j = 0; j < 3; j++)
	{
		int l = (j + 1) % 3;
		double dx = XY[l][0] - XY[j][0];
		double dy = XY[l][1] - XY[j][1];
		if(fabs(dx) < EPS) dx = (dx >= 0) ? EPS : -EPS;

		DS[j] = sqrt(dx*dx + dy*dy);
		AM[j] = dy / dx;
		SM[j] = sqrt(AM[j]*AM[j] + 1.0);
		XD[j] = -dx / DS[j];
		YD[j] =  dy / DS[j];

		if(DS[j] > EPSG) EPSG = DS[j];
	}
	EPSG *= 1.0e-12;

	// Transform observation point to local coordinates
	TVector3d d;
	d.x = obs.x - v0.x;
	d.y = obs.y - v0.y;
	d.z = obs.z - v0.z;

	double EE1 = d.x*AA.x + d.y*AA.y + d.z*AA.z;  // local X
	double EE2 = d.x*BB.x + d.y*BB.y + d.z*BB.z;  // local Y
	double EE3 = d.x*basis_c.x + d.y*basis_c.y + d.z*basis_c.z;  // local Z (height)

	// Distances from observation point to vertices
	double X[3], Y[3], H[3], E[3], R[3];
	for(int j = 0; j < 3; j++)
	{
		X[j] = EE1 - XY[j][0];
		Y[j] = EE2 - XY[j][1];
		H[j] = Y[j] * X[j];
		E[j] = EE3*EE3 + X[j]*X[j];
		R[j] = sqrt(X[j]*X[j] + Y[j]*Y[j] + EE3*EE3);
	}

	double Z = EE3;

	// Edge contributions
	double RM[3], RP[3], RR[3], AL[3];
	for(int j = 0; j < 3; j++)
	{
		int jp1 = (j + 1) % 3;
		RM[j] = R[j] + R[jp1] - DS[j];
		RP[j] = R[j] + R[jp1] + DS[j];
		RR[j] = (RM[j] / RP[j] > EPS) ? (RM[j] / RP[j]) : EPS;
		AL[j] = log(RR[j]);
	}

	// Field components in local frame WITHOUT 4pi divisor
	double HH1 = sigma * (-YD[0]*AL[0] - YD[1]*AL[1] - YD[2]*AL[2]);
	double HH2 = sigma * (-XD[0]*AL[0] - XD[1]*AL[1] - XD[2]*AL[2]);
	double HH3 = 0.0;

	// Normal component (atan terms) - only if not on surface
	if(fabs(Z) > EPSG)
	{
		double ZR[3];
		for(int j = 0; j < 3; j++)
		{
			ZR[j] = Z * R[j];
		}

		double AT[3], BT[3];
		for(int j = 0; j < 3; j++)
		{
			int jp1 = (j + 1) % 3;
			AT[j] = (AM[j]*E[j] - H[j]) / ZR[j];
			BT[j] = (AM[j]*E[jp1] - H[jp1]) / ZR[jp1];
		}

		HH3 = sigma * (-atan(AT[0]) - atan(AT[1]) - atan(AT[2])
		               +atan(BT[0]) + atan(BT[1]) + atan(BT[2]));
	}

	// Transform back to global coordinates
	TVector3d Hfield;
	Hfield.x = HH1*AA.x + HH2*BB.x + HH3*basis_c.x;
	Hfield.y = HH1*AA.y + HH2*BB.y + HH3*basis_c.y;
	Hfield.z = HH1*AA.z + HH2*BB.z + HH3*basis_c.z;

	return Hfield;
}

//-------------------------------------------------------------------------
// FieldFromQuadFaceMirroredWithNormals: Version with explicit normals for IMA boundary faces
//-------------------------------------------------------------------------
TVector3d radTPolyhedron::FieldFromQuadFaceMirroredWithNormals(const TVector3d& obs,
                                                                const TVector3d& V0,
                                                                const TVector3d& V1,
                                                                const TVector3d& V2,
                                                                const TVector3d& V3,
                                                                double sigma,
                                                                const TVector3d& tri1Normal,
                                                                const TVector3d& tri2Normal) const
{
	// Split quad into 2 triangles with explicit normals
	// Triangle 1: V0, V1, V2
	// Triangle 2: V0, V2, V3
	TVector3d H1 = FieldFromChargedTriangleWithNormal(obs, V0, V1, V2, sigma, tri1Normal);
	TVector3d H2 = FieldFromChargedTriangleWithNormal(obs, V0, V2, V3, sigma, tri2Normal);
	return TVector3d(H1.x + H2.x, H1.y + H2.y, H1.z + H2.z);
}

//-------------------------------------------------------------------------

TVector3d radTPolyhedron::FieldFromQuadFace(const TVector3d& obs, int faceIdx, double sigma) const
{
	// Compute field from a single quadrilateral face with unit surface charge
	// Yano MSC method:
	// - Split quad into 2 triangles
	// - For each triangle, check if normal points outward
	// - Apply sign_factor to ensure outward-pointing normal
	// - Returns field WITHOUT 4pi divisor (applied in matrix assembly)

	// Get face vertices
	const radTHandlePgnAndTrans& hpt = VectHandlePgnAndTrans[faceIdx];
	radTPolygon* pgn = hpt.PgnHndl.rep;
	radTrans* tr = hpt.TransHndl.rep;

	const radTVect2dVect& verts2d = pgn->EdgePointsVector;
	if(verts2d.size() < 4) return TVector3d(0.0, 0.0, 0.0);

	TVector3d V0 = tr->TrPoint(TVector3d(verts2d[0].x, verts2d[0].y, pgn->CoordZ));
	TVector3d V1 = tr->TrPoint(TVector3d(verts2d[1].x, verts2d[1].y, pgn->CoordZ));
	TVector3d V2 = tr->TrPoint(TVector3d(verts2d[2].x, verts2d[2].y, pgn->CoordZ));
	TVector3d V3 = tr->TrPoint(TVector3d(verts2d[3].x, verts2d[3].y, pgn->CoordZ));

	TVector3d H_total(0.0, 0.0, 0.0);

	// Standard triangle split for quadrilateral
	// Triangle 1: V0, V1, V2 (indices 0, 1, 2)
	// Triangle 2: V0, V2, V3 (indices 0, 2, 3)
	TVector3d tri_verts[2][3] = {
		{V0, V1, V2},
		{V0, V2, V3}
	};

	for(int iTri = 0; iTri < 2; iTri++)
	{
		const TVector3d& T0 = tri_verts[iTri][0];
		const TVector3d& T1 = tri_verts[iTri][1];
		const TVector3d& T2 = tri_verts[iTri][2];

		// Compute triangle normal
		TVector3d edge1, edge2, tri_normal;
		edge1.x = T1.x - T0.x; edge1.y = T1.y - T0.y; edge1.z = T1.z - T0.z;
		edge2.x = T2.x - T0.x; edge2.y = T2.y - T0.y; edge2.z = T2.z - T0.z;

		tri_normal.x = edge1.y * edge2.z - edge1.z * edge2.y;
		tri_normal.y = edge1.z * edge2.x - edge1.x * edge2.z;
		tri_normal.z = edge1.x * edge2.y - edge1.y * edge2.x;

		double norm_len = sqrt(tri_normal.x*tri_normal.x + tri_normal.y*tri_normal.y + tri_normal.z*tri_normal.z);
		if(norm_len < 1e-20) continue;

		// Normalize
		tri_normal.x /= norm_len;
		tri_normal.y /= norm_len;
		tri_normal.z /= norm_len;

		// Compute triangle center
		TVector3d tri_center;
		tri_center.x = (T0.x + T1.x + T2.x) / 3.0;
		tri_center.y = (T0.y + T1.y + T2.y) / 3.0;
		tri_center.z = (T0.z + T1.z + T2.z) / 3.0;

		// Check if normal points outward (away from element center)
		TVector3d to_center;
		to_center.x = tri_center.x - CentrPoint.x;
		to_center.y = tri_center.y - CentrPoint.y;
		to_center.z = tri_center.z - CentrPoint.z;

		double dot_prod = tri_normal.x * to_center.x + tri_normal.y * to_center.y + tri_normal.z * to_center.z;

		// Sign factor: +1 if normal points outward, -1 if inward
		double sign_factor = (dot_prod >= 0.0) ? 1.0 : -1.0;

		// Compute field from this triangle with sign-corrected sigma
		TVector3d H_tri = FieldFromChargedTriangle(obs, T0, T1, T2, sigma * sign_factor);

		H_total.x += H_tri.x;
		H_total.y += H_tri.y;
		H_total.z += H_tri.z;
	}

	return H_total;
}

//-------------------------------------------------------------------------
// FieldFromFace: Generalized field from a single face (triangular or quadrilateral)
// Dispatches based on vertex count. Returns field WITHOUT 4pi divisor.
//-------------------------------------------------------------------------
TVector3d radTPolyhedron::FieldFromFace(const TVector3d& obs, int faceIdx, double sigma) const
{
	// Get face polygon info
	const radTHandlePgnAndTrans& hpt = VectHandlePgnAndTrans[faceIdx];
	radTPolygon* pgn = hpt.PgnHndl.rep;
	radTrans* tr = hpt.TransHndl.rep;

	const radTVect2dVect& verts2d = pgn->EdgePointsVector;
	int nv = (int)verts2d.size();

	if(nv >= 4)
	{
		// Quadrilateral face - delegate to existing FieldFromQuadFace
		return FieldFromQuadFace(obs, faceIdx, sigma);
	}
	else if(nv == 3)
	{
		// Triangular face - compute field from single triangle with outward normal check
		TVector3d V0 = tr->TrPoint(TVector3d(verts2d[0].x, verts2d[0].y, pgn->CoordZ));
		TVector3d V1 = tr->TrPoint(TVector3d(verts2d[1].x, verts2d[1].y, pgn->CoordZ));
		TVector3d V2 = tr->TrPoint(TVector3d(verts2d[2].x, verts2d[2].y, pgn->CoordZ));

		// Compute triangle normal from cross product
		TVector3d edge1, edge2, tri_normal;
		edge1.x = V1.x - V0.x; edge1.y = V1.y - V0.y; edge1.z = V1.z - V0.z;
		edge2.x = V2.x - V0.x; edge2.y = V2.y - V0.y; edge2.z = V2.z - V0.z;

		tri_normal.x = edge1.y * edge2.z - edge1.z * edge2.y;
		tri_normal.y = edge1.z * edge2.x - edge1.x * edge2.z;
		tri_normal.z = edge1.x * edge2.y - edge1.y * edge2.x;

		double norm_len = sqrt(tri_normal.x*tri_normal.x + tri_normal.y*tri_normal.y + tri_normal.z*tri_normal.z);
		if(norm_len < 1e-20) return TVector3d(0.0, 0.0, 0.0);

		tri_normal.x /= norm_len;
		tri_normal.y /= norm_len;
		tri_normal.z /= norm_len;

		// Compute triangle center
		TVector3d tri_center;
		tri_center.x = (V0.x + V1.x + V2.x) / 3.0;
		tri_center.y = (V0.y + V1.y + V2.y) / 3.0;
		tri_center.z = (V0.z + V1.z + V2.z) / 3.0;

		// Check if normal points outward (away from element center)
		TVector3d to_center;
		to_center.x = tri_center.x - CentrPoint.x;
		to_center.y = tri_center.y - CentrPoint.y;
		to_center.z = tri_center.z - CentrPoint.z;

		double dot_prod = tri_normal.x * to_center.x + tri_normal.y * to_center.y + tri_normal.z * to_center.z;

		// Sign factor: +1 if normal points outward, -1 if inward
		double sign_factor = (dot_prod >= 0.0) ? 1.0 : -1.0;

		// Compute field from this triangle with sign-corrected sigma
		return FieldFromChargedTriangle(obs, V0, V1, V2, sigma * sign_factor);
	}

	return TVector3d(0.0, 0.0, 0.0);
}

//-------------------------------------------------------------------------

TVector3d radTPolyhedron::FieldFromPointCharge(const TVector3d& obs, double charge) const
{
	// Magnetic field from point magnetic charge (monopole) at element center
	// H = m * (r - p) / |r - p|^3
	// Returns field WITHOUT 4pi divisor (matches ELF_MAGIC convention)
	// 4pi is applied in matrix assembly

	TVector3d r_minus_p;
	r_minus_p.x = obs.x - CentrPoint.x;
	r_minus_p.y = obs.y - CentrPoint.y;
	r_minus_p.z = obs.z - CentrPoint.z;

	double r_mag_sq = r_minus_p.x * r_minus_p.x +
	                  r_minus_p.y * r_minus_p.y +
	                  r_minus_p.z * r_minus_p.z;

	if(r_mag_sq < 1e-20)
	{
		return TVector3d(0.0, 0.0, 0.0);
	}

	double r_mag = sqrt(r_mag_sq);
	double r_mag_cubed = r_mag_sq * r_mag;
	// No 4pi divisor here - matches ELF_MAGIC convention
	double coef = charge / r_mag_cubed;

	TVector3d H;
	H.x = coef * r_minus_p.x;
	H.y = coef * r_minus_p.y;
	H.z = coef * r_minus_p.z;

	return H;
}

//-------------------------------------------------------------------------
// MSC per-face collocation point (MscEvalPoint) for the external-field sampling.
//-------------------------------------------------------------------------

bool g_yano_moment_hacapk = false;   // moment linear step via the HACApK H-matrix + BiCGSTAB (method 2 / scalable storage); set by SolveGen for moment-eligible + method 2 (else dense LU)
// NOTE (Phase 3b-1, 2026-06-22): the g_yano_moment opt-out flag was REMOVED.  moment-yano is now the
// UNCONDITIONAL surface-charge demag (hex 6-DOF + wedge 5-DOF, method 0/1/2); there is no EIEM2 opt-out.
// EIEM2 survives only for mixed tet+MSC until 3b-2 deletes it.




// Per-face MSC collocation point (global frame): midpoint of the face center and the element
// center (alpha = 0.5).  Used by SetupExternFieldArray to sample the external field per MSC face.
TVector3d radTPolyhedron::MscEvalPoint(int faceIdx) const
{
	return TVector3d(0.5*(FaceCenter[faceIdx].x + CentrPoint.x),
	                 0.5*(FaceCenter[faceIdx].y + CentrPoint.y),
	                 0.5*(FaceCenter[faceIdx].z + CentrPoint.z));
}



//-------------------------------------------------------------------------
// B_genComp: Simplified version without TrfMlt support
// TrfMlt has been removed from Radia - use explicit element duplication instead
//-------------------------------------------------------------------------
void radTPolyhedron::B_genComp(radTField* FieldPtr)
{
	radTFieldKey& FieldKey = FieldPtr->FieldKey;

	// Handle special field keys
	if(FieldKey.Ib_ || FieldKey.Ih_)
	{
		B_intComp(FieldPtr);
		return;
	}
	if(FieldKey.Force_)
	{
		IntOverShape(FieldPtr);
		return;
	}

	// Standard B_comp for all cases
	B_comp(FieldPtr);
}

//-------------------------------------------------------------------------

void radTPolyhedron::B_comp_frM(radTField* FieldPtr)
//void radTPolyhedron::B_comp(radTField* FieldPtr)
{
	// =========================================================================
	// Dispatch to specialized MSC methods based on element type
	// =========================================================================
	// Supported element types:
	// - Tetrahedron: 4 triangular faces (AmOfFaces == 4)
	// - Hexahedron: 6 quadrilateral faces (AmOfFaces == 6)
	// =========================================================================

	// For tetrahedral elements, use the analytical MSC method
	// The analytical method uses closed-form surface charge formulas and has been
	// verified to produce identical results to the original Gauss integration method.
	if(IsTetrahedron())
	{
		B_comp_tetrahedron_analytical(FieldPtr);
		return;
	}

	// For hexahedral elements (6 quadrilateral faces), use the MSC method
	if(IsHexahedron())
	{
		B_comp_hexahedron_MSC(FieldPtr);
		return;
	}

	// For wedge elements (5 faces: 2 triangular + 3 quadrilateral)
	// Use 5-DOF MSC when enabled, otherwise fall back to 3-DOF analytical
	if(AmOfFaces == 5)
	{
		if(Use6DOF_MSC)
			B_comp_wedge_MSC(FieldPtr);
		else
			B_comp_wedge_analytical(FieldPtr);
		return;
	}
	else
	{
		// Unsupported polyhedron type - issue warning but continue with generic method
		// throw std::runtime_error("Unsupported polyhedron type: only tetrahedra (4 faces) and hexahedra (6 faces) are supported. AmOfFaces=" + std::to_string(AmOfFaces));
	}

	// NOTE: The following code handles wedge elements (5 faces) and other polyhedra
	// using generic polygon-based field computation (Gauss integration).

	// Use standard polygon-based computation for all polyhedra (including tetrahedra)
	TVector3d Zero(0.,0.,0.);
	
	// Determine if observation point is inside polyhedron using geometric test
	// For each face, compute face normal from vertices and check if point is on inside
	short PointIsInside = 1;
	if(FieldPtr->FieldKey.M_ || FieldPtr->FieldKey.B_)
	{
		TVector3d& P = FieldPtr->P;
		for(int i = 0; i < AmOfFaces && PointIsInside; i++)
		{
			const radTHandlePgnAndTrans& hpt = VectHandlePgnAndTrans[i];
			radTPolygon* pgn = hpt.PgnHndl.rep;
			radTrans* tr = hpt.TransHndl.rep;
			
			// Get first 3 vertices of the polygon in GLOBAL coordinates
			radTVect2dVect& verts2d = pgn->EdgePointsVector;
			if(verts2d.size() < 3) continue;
			
			// Transform vertices from local (2D + CoordZ) to global
			TVector3d v0 = tr->TrPoint(TVector3d(verts2d[0].x, verts2d[0].y, pgn->CoordZ));
			TVector3d v1 = tr->TrPoint(TVector3d(verts2d[1].x, verts2d[1].y, pgn->CoordZ));
			TVector3d v2 = tr->TrPoint(TVector3d(verts2d[2].x, verts2d[2].y, pgn->CoordZ));
			
			// Compute face normal from cross product
			TVector3d edge1 = v1 - v0;
			TVector3d edge2 = v2 - v0;
			TVector3d faceNormal;
			faceNormal.x = edge1.y * edge2.z - edge1.z * edge2.y;
			faceNormal.y = edge1.z * edge2.x - edge1.x * edge2.z;
			faceNormal.z = edge1.x * edge2.y - edge1.y * edge2.x;
			
			// Normalize
			double normLen = sqrt(faceNormal.x*faceNormal.x + faceNormal.y*faceNormal.y + faceNormal.z*faceNormal.z);
			if(normLen < 1e-15) continue;  // Degenerate face
			faceNormal.x /= normLen;
			faceNormal.y /= normLen;
			faceNormal.z /= normLen;
			
			// Compute face center
			TVector3d faceCenter;
			faceCenter.x = faceCenter.y = faceCenter.z = 0;
			for(size_t j = 0; j < verts2d.size(); j++)
			{
				TVector3d vj = tr->TrPoint(TVector3d(verts2d[j].x, verts2d[j].y, pgn->CoordZ));
				faceCenter.x += vj.x;
				faceCenter.y += vj.y;
				faceCenter.z += vj.z;
			}
			faceCenter.x /= verts2d.size();
			faceCenter.y /= verts2d.size();
			faceCenter.z /= verts2d.size();
			
			// Ensure normal points outward (away from polyhedron centroid)
			TVector3d toCenter = CentrPoint - faceCenter;
			double dot_toCenter = toCenter.x*faceNormal.x + toCenter.y*faceNormal.y + toCenter.z*faceNormal.z;
			if(dot_toCenter > 0) {
				// Normal points toward center, flip it
				faceNormal.x = -faceNormal.x;
				faceNormal.y = -faceNormal.y;
				faceNormal.z = -faceNormal.z;
			}
			
			// Check if observation point is on inside of this face
			// Inside means (P - faceCenter) dot normal < 0
			TVector3d toP = P - faceCenter;
			double dot_toP = toP.x*faceNormal.x + toP.y*faceNormal.y + toP.z*faceNormal.z;
			
			// Allow small tolerance for points on the surface
			double faceDist = sqrt(toCenter.x*toCenter.x + toCenter.y*toCenter.y + toCenter.z*toCenter.z);
			double tol = faceDist * 1e-6;  // Use larger tolerance
			if(tol < 1e-10) tol = 1e-10;
			
			if(dot_toP > tol)
			{
				// Point is outside this face
				PointIsInside = 0;
			}
		}
	}

	// =========================================================================
	// For tetrahedra in PreRelax mode, use direct field computation
	// =========================================================================
	// DISABLED: Testing standard polygon path instead
	if(false && FieldPtr->FieldKey.PreRelax_ && IsTetrahedron())
	{
		// Compute the full demagnetization matrix directly
		// by evaluating H for each of the 3 unit magnetization directions
		TVector3d saved_Magn = Magn;
		TVector3d obsP = FieldPtr->P;
		
		// Compute H for Mx=1
		Magn = TVector3d(1., 0., 0.);
		TVector3d H_from_Mx(0., 0., 0.);
		for(int i=0; i<AmOfFaces; i++)
		{
			const radTHandlePgnAndTrans& hpt = VectHandlePgnAndTrans[i];
			radTPolygon* pgn = hpt.PgnHndl.rep;
			radTrans* tr = hpt.TransHndl.rep;
			
			radTFieldKey key;
			key.H_ = 1;
			radTField fld(key, FieldPtr->CompCriterium, tr->TrPoint_inv(obsP), Zero, Zero, Zero, Zero, Zero);
			pgn->Magn = tr->TrVectField_inv(Magn);
			pgn->B_comp(&fld);
			H_from_Mx += tr->TrVectField(fld.H);
		}
		
		// Compute H for My=1
		Magn = TVector3d(0., 1., 0.);
		TVector3d H_from_My(0., 0., 0.);
		for(int i=0; i<AmOfFaces; i++)
		{
			const radTHandlePgnAndTrans& hpt = VectHandlePgnAndTrans[i];
			radTPolygon* pgn = hpt.PgnHndl.rep;
			radTrans* tr = hpt.TransHndl.rep;
			
			radTFieldKey key;
			key.H_ = 1;
			radTField fld(key, FieldPtr->CompCriterium, tr->TrPoint_inv(obsP), Zero, Zero, Zero, Zero, Zero);
			pgn->Magn = tr->TrVectField_inv(Magn);
			pgn->B_comp(&fld);
			H_from_My += tr->TrVectField(fld.H);
		}
		
		// Compute H for Mz=1
		Magn = TVector3d(0., 0., 1.);
		TVector3d H_from_Mz(0., 0., 0.);
		for(int i=0; i<AmOfFaces; i++)
		{
			const radTHandlePgnAndTrans& hpt = VectHandlePgnAndTrans[i];
			radTPolygon* pgn = hpt.PgnHndl.rep;
			radTrans* tr = hpt.TransHndl.rep;
			
			radTFieldKey key;
			key.H_ = 1;
			radTField fld(key, FieldPtr->CompCriterium, tr->TrPoint_inv(obsP), Zero, Zero, Zero, Zero, Zero);
			pgn->Magn = tr->TrVectField_inv(Magn);
			pgn->B_comp(&fld);
			H_from_Mz += tr->TrVectField(fld.H);
		}
		
		// Restore original magnetization
		Magn = saved_Magn;

		// Store negated values to match solver convention A = -N + 1/chi
		// The solver expects N such that H_demag = N * M, but surface charge method
		// computes H = -N_demag * M, so we need to negate here.
		// NOTE: This should match the tetrahedron_analytical code which also uses +=
		// Both tetra and wedge should store the same sign convention.
		FieldPtr->B.x -= H_from_Mx.x;
		FieldPtr->B.y -= H_from_Mx.y;
		FieldPtr->B.z -= H_from_Mx.z;
		FieldPtr->H.x -= H_from_My.x;
		FieldPtr->H.y -= H_from_My.y;
		FieldPtr->H.z -= H_from_My.z;
		FieldPtr->A.x -= H_from_Mz.x;
		FieldPtr->A.y -= H_from_Mz.y;
		FieldPtr->A.z -= H_from_Mz.z;
		return;
	}

	radTFieldKey LocFieldKey = FieldPtr->FieldKey;
	if(LocFieldKey.B_) LocFieldKey.H_ = 1;
	radTField SumLocField(LocFieldKey, FieldPtr->CompCriterium, FieldPtr->P, Zero, Zero, Zero, Zero, Zero);

	for(int i=0; i<AmOfFaces; i++)
	{
		radTHandlePgnAndTrans HandlePgnAndTrans = VectHandlePgnAndTrans[i];

		radTPolygon* PgnPtr = HandlePgnAndTrans.PgnHndl.rep;
		radTrans* TransPtr = HandlePgnAndTrans.TransHndl.rep;

		radTField LocField(LocFieldKey, SumLocField.CompCriterium, Zero, Zero, Zero, Zero, Zero, Zero);
		LocField.P = TransPtr->TrPoint_inv(SumLocField.P);
		TVector3d PrevP = LocField.P;

				if(!LocFieldKey.PreRelax_) PgnPtr->Magn = TransPtr->TrVectField_inv(Magn);
		
		// B_comp computes field in LOCAL coordinates
		// TrField below transforms result to GLOBAL coordinates
		PgnPtr->B_comp(&LocField);

		if(LocField.P != PrevP)
		{
			SumLocField.P = TransPtr->TrPoint(LocField.P); //OC040504 test
		}

		// Note: We no longer use LocField.PointIsInsideFrame for the inside check
		// as it's unreliable. The geometric check above handles this.

		if(!LocFieldKey.PreRelax_) SumLocField += TransPtr->TrField(LocField);
		else
		{
			TMatrix3d Q(LocField.B, LocField.H, LocField.A);
			TransPtr->TrMatrixLeft_inv(Q);
			TransPtr->TrMatrix(Q);
			SumLocField.B += Q.Str0;
			SumLocField.H += Q.Str1;
			SumLocField.A += Q.Str2;
		}
	}

	//FieldPtr->P = OrigP; //OC090908

	radTFieldKey& FldKey = FieldPtr->FieldKey;
	if(FldKey.PreRelax_)
	{
		// Self-term handling for tetrahedral elements
		// When observation point is at element centroid (self-interaction),
		// no additional term needed for direct solvers.
		// The surface integral already computes the correct self-demagnetization.
		if(IsTetrahedron())
		{
			TVector3d diff = FieldPtr->P - CentrPoint;
			double distSq = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
			// Use element size for tolerance
			double elemSizeSq = 0.0;
			for(int ii = 0; ii < AmOfFaces; ii++)
			{
				const radTHandlePgnAndTrans& hpt = VectHandlePgnAndTrans[ii];
				radTPolygon* pgn = hpt.PgnHndl.rep;
				radTrans* tr = hpt.TransHndl.rep;
				TVector3d faceCtr = tr->TrPoint(TVector3d(pgn->CentrPoint.x, pgn->CentrPoint.y, pgn->CoordZ));
				TVector3d d = faceCtr - CentrPoint;
				double dSq = d.x*d.x + d.y*d.y + d.z*d.z;
				if(dSq > elemSizeSq) elemSizeSq = dSq;
			}
			double tol = elemSizeSq * 1.0e-10;
			if(tol < 1.0e-20) tol = 1.0e-20;

			if(distSq < tol)
			{
				// Self-interaction: For direct solver, no additional term needed.
				// The surface integral computes the correct self-demagnetization.
			}
		}
		
		FieldPtr->B += SumLocField.B;
		FieldPtr->H += SumLocField.H;
		FieldPtr->A += SumLocField.A;
		return;
	}
	if(FldKey.H_) FieldPtr->H += SumLocField.H;
	if(FldKey.M_) if(PointIsInside) FieldPtr->M += Magn;
	if(FldKey.B_)
	{
		FieldPtr->B += SumLocField.H;
		if(PointIsInside) FieldPtr->B += Magn;
	}
	if(FldKey.A_) FieldPtr->A += SumLocField.A;
}

//-------------------------------------------------------------------------

void radTPolyhedron::B_comp_frJ(radTField* pField)
{
	TVector3d vEx(1.,0.,0.), vEy(0.,1.,0.), vEz(0.,0.,1.);
	const double Pi = 3.14159265358979;
	const double relPrecSwitchRootDecomp = 1E-08;

	TVector3d &PobsLab = pField->P;
	TMatrix3d QT; //for linear terms of J in local frames of faces
	TVector3d vSumFaces0(0,0,0), vSumFaces1(0,0,0), vSumFaces3(0,0,0); 
	double sumFaces2 = 0;
	TVector2d vQTX, vQTY, vI2, vSegm, vSegmUnit, vSegmExtNorm, vPobsProjToVertex, vProjToVertex1, vProjToVertex2;
	TVector2d vShiftCenPointVertex, PobsProj;
	double qtXZ, qtYZ;
	double I1, partI3;
	double R1, R2, s1_p_R1, s2_p_R2, PiMult1;

	bool PointIsInside = true;

	for(int i=0; i<AmOfFaces; i++)
	{
		radTHandlePgnAndTrans hPgnAndTrans = VectHandlePgnAndTrans[i];
		radTPolygon* pPgn = hPgnAndTrans.PgnHndl.rep;
		radTrans* pTrans = hPgnAndTrans.TransHndl.rep;

		radTVect2dVect &vPgnVertices = pPgn->EdgePointsVector;

		TVector2d &pgnFirstP = vPgnVertices[0];
		double AbsRandX = radCR.AbsRandMagnitude(pgnFirstP.x - pPgn->CentrPoint.x);
		double AbsRandY = radCR.AbsRandMagnitude(pgnFirstP.y - pPgn->CentrPoint.y);
		double AbsRandXY = AbsRandX;
		if(AbsRandXY < AbsRandY) AbsRandXY = AbsRandY;
		double AbsRandZ = radCR.AbsRandMagnitude(pPgn->CoordZ);
		if(AbsRandZ < AbsRandXY) AbsRandZ = AbsRandXY;

		TVector3d Pobs = pTrans->TrPoint_inv(PobsLab);
		double hi = pPgn->CoordZ - Pobs.z; //to check

		// Artificial shift of an observation point
		// if the point is exactly on the border (to avoid "divide by zero" error):
		if(hi == 0.) 
		{
			hi = AbsRandZ;
		}

		if(hi < 0.) PointIsInside = false; //to check

		double hiE2 = hi*hi;
		double abs_hi = ::fabs(hi);

		//TVector2d PobsProj(Pobs.x, Pobs.y);
		PobsProj.x = Pobs.x; PobsProj.y = Pobs.y;
		bool PobsProjIsInside = true;

		int pgnAmOfVertices = pPgn->AmOfEdgePoints;
		int pgnAmOfVertices_mi_1 = pgnAmOfVertices - 1;
		int j2 = 1;
		I1 = 0.;
		partI3 = 0.;
		vI2.x = vI2.y = 0.;
		double sum_hij = 0.;

		for(int j=0; j<pgnAmOfVertices; j++)
		{
			if(j == pgnAmOfVertices_mi_1) j2 = 0;

			TVector2d &vPgnVertex1 = vPgnVertices[j];

			vSegm = vPgnVertices[j2++] - vPgnVertex1; //don't use j2 after this!
			vSegmUnit = vSegm;
			vSegmUnit.Normalize();
			vSegmExtNorm.x = vSegmUnit.y; vSegmExtNorm.y = -vSegmUnit.x;

			vPobsProjToVertex = vPgnVertex1 - PobsProj;
			//vPobsProjToVertex = (vPgnVertex1 - vShiftCenPointVertex) - PobsProj;
			double hij = vPobsProjToVertex*vSegmExtNorm;
			// Artificial shift of an observation point
			// if the point is exactly on the border (to avoid "divide by zero" error):
			if(hij == 0.)
			{
				hij = AbsRandXY;
				//slight displacement of the observation point is also required...?
			}
			double abs_hij = ::fabs(hij);
			double sign_hij = Sign(hij);

			if(hij < 0)
			{
				PobsProjIsInside = false;
			}
			sum_hij += hij;

			vProjToVertex1 = vPobsProjToVertex - (hij*vSegmExtNorm);
			vProjToVertex2 = vProjToVertex1 + vSegm;

			double s1 = vProjToVertex1*vSegmUnit;
			double s2 = vProjToVertex2*vSegmUnit;

			//double scalProdProjToVert = vProjToVertex1*vProjToVertex2;
			//if(scalProdProjToVert >= 0)
			//{
			if(s2 < s1) //??
			{
				double sBuf = s2;
				s2 = s1; s1 = sBuf;
			}
			//}
			//else 
			//{//to check
			//	s1 = -s1;
			//}

			double s1e2 = s1*s1;
			double s2e2 = s2*s2;
			double hiE2_p_hijE2 = hiE2 + hij*hij;

			if((hiE2_p_hijE2 < s1e2*relPrecSwitchRootDecomp) && (s1 < 0.))
			{
				s1_p_R1 = -0.5*hiE2_p_hijE2/s1;
				R1 = -s1 + s1_p_R1;
			}
			else
			{
				R1 = sqrt(hiE2_p_hijE2 + s1e2);
				s1_p_R1 = s1 + R1; 
			}

			if((hiE2_p_hijE2 < s2e2*relPrecSwitchRootDecomp) && (s2 < 0.))
			{
				s2_p_R2 = -0.5*hiE2_p_hijE2/s2;
				R2 = -s2 + s2_p_R2;
			}
			else
			{
				R2 = sqrt(hiE2_p_hijE2 + s2e2);
				s2_p_R2 = s2 + R2; 
			}

			//R1 = sqrt(hiE2_p_hijE2 + s1*s1);
			//R2 = sqrt(hiE2_p_hijE2 + s2*s2);
			//if(s1_p_R1 <= 0) 
			//{
			//	s1_p_R1 = AbsRandXY;
			//}
			//if(s2_p_R2 <= 0) 
			//{
			//	s2_p_R2 = AbsRandXY;
			//}

			double LogDif = log(s2_p_R2/(s1_p_R1)); 
			
			double ArgAtan1 = hi*s1/(hij*R1);
			double ArgAtan2 = hi*s2/(hij*R2);

			double SumAtan1 = atan(TransAtans(ArgAtan2, -ArgAtan1, PiMult1));
			SumAtan1 += Pi*PiMult1;
			double addAtanI1 = hi*SumAtan1 + hij*LogDif;

			I1 += addAtanI1;

			double multI2 = 0.5*(s2*R2 - s1*R1 + hiE2_p_hijE2*LogDif);
			vI2 += multI2*vSegmExtNorm;

			partI3 += hij*multI2;
		}

		double I1_0 = I1;
		if(PobsProjIsInside)
		{
			I1 -= 2*Pi*abs_hi; //!!!
		}

		TVector3d vFaceN = pTrans->TrBiPoint(vEz); //face normal in laboratory frame?
		TVector3d vExLab = pTrans->TrBiPoint(vEx);
		TVector3d vEyLab = pTrans->TrBiPoint(vEy);

		vSumFaces0 += I1*vFaceN;
		
		double hi_I1 = hi*I1; //required for A
		sumFaces2 += hi_I1;

		if(pJ_LinCoef != 0)
		{
			QT = *pJ_LinCoef;
			//pTrans->TrMatrixGeom(QT); //linear terms of J in local frame of the current face - to check!!
			//pTrans->TrMatrixGeomLeft_inv(QT); //linear terms of J in local frame of the current face - to check!!
			pTrans->TrMatrixGeomLeft(QT); //linear terms of J in local frame of the current face - to check!!
			pTrans->TrMatrixGeom_inv(QT); //linear terms of J in local frame of the current face - to check!!

			vQTX.x = QT.Str0.x; vQTX.y = QT.Str0.y;
			vQTY.x = QT.Str1.x; vQTY.y = QT.Str1.y;
			qtXZ = QT.Str0.z;
			qtYZ = QT.Str1.z;

			vSumFaces1 += ((vI2*vQTY)*vExLab) - ((vI2*vQTX)*vEyLab);
			vSumFaces1 += hi_I1*((qtYZ*vExLab) - (qtXZ*vEyLab));

			if(pField->FieldKey.A_)
			{
				double I3 = hi*hi_I1 + partI3;
				vSumFaces3 += I3*vFaceN;
			}
		}
	}

	// Biot-Savart constant: B = (mu_0/4*pi) * integral(J x r / r^3) dV
	// Radia now uses SI units (meters) internally, matching ELF.
	// ConstForJ = mu_0/(4*pi) = 1e-7 H/m
	const double ConstForJ = 1.0e-7;
	TVector3d Jmain = J;
	if((mLinTreat == 0) && (pJ_LinCoef != 0)) //treat as being relative 
	{
		Jmain -= ((*pJ_LinCoef)*CentrPoint);
	}
	if(pJ_LinCoef != 0)
	{
		if(pField->FieldKey.A_ || pField->FieldKey.B_ || pField->FieldKey.H_ || pField->FieldKey.J_)
		{
			Jmain += (*pJ_LinCoef)*PobsLab;
		}
	}

	if(pField->FieldKey.J_) 
	{
		if(PointIsInside) pField->J += Jmain;
	}
	if(pField->FieldKey.A_)
	{
		TVector3d vA(0,0,0);
		if(pJ_LinCoef != 0)
		{
			//Jmain += (*pJ_LinCoef)*PobsLab;
			vA += (*pJ_LinCoef)*((ConstForJ/3.)*vSumFaces3);
		}
		vA += (0.5*ConstForJ*sumFaces2)*Jmain;
		pField->A += vA;
	}
	if(pField->FieldKey.B_ || pField->FieldKey.H_)
	{
		TVector3d vB(0,0,0);
		if(pJ_LinCoef != 0)
		{
			//Jmain += (*pJ_LinCoef)*PobsLab;
			TVector3d vQc(pJ_LinCoef->Str1.z - pJ_LinCoef->Str2.y, 
						  pJ_LinCoef->Str2.x - pJ_LinCoef->Str0.z, 
						  pJ_LinCoef->Str0.y - pJ_LinCoef->Str1.x);
			vB += ConstForJ*(vSumFaces1 - ((0.5*sumFaces2)*vQc));
		}
		vB += ConstForJ*(Jmain^vSumFaces0);
		
		pField->B += vB;
		pField->H += vB;
	}
	//if(pField->FieldKey.J_) 
	//{
	//	if(PointIsInside) 
	//	{
	//		pField->J += J;
	//		if(pJ_LinCoef != 0)
	//		{
	//			pField->J += (mLinTreat == 0)? (*pJ_LinCoef)*(PobsLab - CentrPoint) : (*pJ_LinCoef)*PobsLab;
	//		}
	//	}
	//}
}

//-------------------------------------------------------------------------

void radTPolyhedron::B_intComp_frM(radTField* FieldPtr)
//void radTPolyhedron::B_intComp(radTField* FieldPtr)
{
	if(FieldPtr->FieldKey.FinInt_) { B_intCompFinNum(FieldPtr); return;}

	TVector3d Zero(0.,0.,0.);
	for(int i=0; i<AmOfFaces; i++)
	{
		radTHandlePgnAndTrans HandlePgnAndTrans = VectHandlePgnAndTrans[i];

		radTPolygon* PgnPtr = HandlePgnAndTrans.PgnHndl.rep;
		radTrans* TransPtr = HandlePgnAndTrans.TransHndl.rep;

		radTField LocField(FieldPtr->FieldKey, FieldPtr->CompCriterium, Zero, Zero, Zero, Zero, Zero, Zero);

		PgnPtr->Magn = TransPtr->TrVectField_inv(Magn);

		LocField.P = TransPtr->TrPoint_inv(FieldPtr->P);
		LocField.NextP = TransPtr->TrPoint_inv(FieldPtr->NextP);
		//LocField.P = TransPtr->TrPoint_inv(FieldPtr->P - CentrPoint); //OC090908
		//LocField.NextP = TransPtr->TrPoint_inv(FieldPtr->NextP - CentrPoint); //OC090908

		PgnPtr->B_intComp(&LocField);
		*FieldPtr += TransPtr->TrField(LocField);
	}
}

//-------------------------------------------------------------------------

void radTPolyhedron::B_intComp_frJ(radTField* pField)
{
	if(pField->FieldKey.FinInt_) { B_intCompFinNum(pField); return;} //to check if this still works

	//- find rotation which transforms integration line/vector to Ez;
	//- for each face polygon:
		//* find coordinates of vertex points in the frame where the integration line is along Ez;
		//* calculate face contribution to the total field integral value (in the above frame);
	//- transform the total field integral vector to the laboratory frame

	const double Pi = 3.14159265358979;
	const double inv12 = 1./12.;
	double relTolPerp = 10.*radCR.RelRand; //1.E-11; //to tune
	TVector3d vIntAx = pField->NextP - pField->P, Ez(0,0,1.);
	TVector3d vIntAxUnit = vIntAx; vIntAxUnit.Normalize();

	radTrans trIntAxis2Ez, trAuxRot, trAuxTransl;
	trAuxRot.SetupRotation(CentrPoint, vIntAxUnit, Ez);

	TVector3d vIntAxCenPoint = 0.5*(pField->NextP + pField->P);
	TVector3d vTestIntAxCenPoint = trAuxRot.TrPoint(vIntAxCenPoint);
	TVector3d vTransl(-vTestIntAxCenPoint.x, -vTestIntAxCenPoint.y, 0.);
	//TVector3d vTransl(-vTestIntAxCenPoint.x, -vTestIntAxCenPoint.y, -vTestIntAxCenPoint.z);
	trAuxTransl.SetupTranslation(vTransl);
	TrProduct(&trAuxTransl, &trAuxRot, trIntAxis2Ez);

	TVector3d Jrot = trIntAxis2Ez.TrVectField(J); //to check this !!!
	double qxx, qxy, qxz;
	double qyx, qyy, qyz;
	double qzx, qzy, qzz;
	double qxxmqyy, qxypqyx, twoqxx, twoqzx, twoqzy;

	if(pJ_LinCoef != 0)
	{
		TMatrix3d Qrot = *pJ_LinCoef;
		trIntAxis2Ez.TrMatrixGeomLeft_inv(Qrot); //{ Qrot = Qrot*M_inv;} //to check
		trIntAxis2Ez.TrMatrix(Qrot); //{ Qrot = s*M*Qrot;}

		TVector3d vRefP(0,0,0);

		if(mLinTreat == 0) //treat as being relative 
		{
			vRefP = CentrPoint;
			//TVector3d vCenPtIntFrame = trIntAxis2Ez.TrPoint(CentrPoint);
			//Jrot -= (Qrot*vCenPtIntFrame); //to check
		}
		TVector3d vRefPIntFrame = trIntAxis2Ez.TrPoint(vRefP);
		Jrot -= (Qrot*vRefPIntFrame); //to check

		TVector3d &st0 = Qrot.Str0, &st1 = Qrot.Str1, &st2 = Qrot.Str2;
		qxx = st0.x; qxy = st0.y; qxz = st0.z;
		qyx = st1.x; qyy = st1.y; qyz = st1.z;
		qzx = st2.x; qzy = st2.y; qzz = st2.z;

		qxxmqyy = qxx - qyy; qxypqyx = qxy + qyx;
		twoqxx = 2.*qxx;
		twoqzx = 2.*qzx; twoqzy = 2.*qzy;
	}
	double axe2, aye2, cze2;
	double trecze2qzz, trecze2qxz, tricze2qyz, axe2qzz, axe2qxz, aye2qyz, axqzz, axqyz, ayqzz, ayqxz, ayqyz;
	double axqyzp2qyxp2qxy, qxxpaxqxzmqyy, qxxpaxqxzmqyymayqyz, qxypayqxzpqyxpaxqyz;
	double qzypayqzz, twoqzypayqzz, aytwoqzypayqzz, twoqxypayqxzptwoqyx, twoqxxmqyymayqyz, twoqxxpaxqxzmqyy;
	double twoax, twocz, tricz, tricze2, triczqzz;
	double czqxz, triczqxz, czqyz, triczqyz;
	double qzxpaxqzz, twoaxqzz;

	double jx0 = Jrot.x, jy0 = Jrot.y, jz0 = Jrot.z;
	double quartJz0 = 0.25*jz0;
	double twojx0 = 2.*jx0, trijx0 = 3.*jx0, twojy0 = 2.*jy0, trijy0 = 3.*jy0, twojz0 = 2.*jz0, trijz0 = 3*jz0;

	radTrans trFace2Int;
	TVector3d vR1loc, vR1int, vR2loc, vR2int, vSegm3d, vSegm3dUnit, vSegm3dExtNorm, vResLocIB(0,0,0);
	TVector2d vSegmExtNorm;
	double auxBuf, piMult1 = 0;
	double IBxLoc, IByLoc, IBzLoc;

	for(int i=0; i<AmOfFaces; i++)
	{
		radTHandlePgnAndTrans hPgnAndTrans = VectHandlePgnAndTrans[i];
		radTPolygon* pPgn = hPgnAndTrans.PgnHndl.rep;
		radTrans* pTrans = hPgnAndTrans.TransHndl.rep;
		TrProduct(&trIntAxis2Ez, pTrans, trFace2Int);

		TVector3d vFaceNormIntFrame = trFace2Int.TrBiPoint(Ez);
		//if this vector is perpendicular to Ez, then don't compute contribution from this face
		if(::fabs(vFaceNormIntFrame.z) < relTolPerp) continue; //to check consistency with the subsequent
		
		radTVect2dVect &vPgnVertices = pPgn->EdgePointsVector;
		TVector2d &vert0_2d = vPgnVertices[0];
		TVector3d vR0loc(vert0_2d.x, vert0_2d.y, pPgn->CoordZ);
		TVector3d vR0int = trFace2Int.TrPoint(vR0loc);

		double signZ_FaceNorm = Sign(vFaceNormIntFrame.z);
		double ax = -vFaceNormIntFrame.x/vFaceNormIntFrame.z;
		double ay = -vFaceNormIntFrame.y/vFaceNormIntFrame.z;
		double cz = vR0int.z - ax*vR0int.x - ay*vR0int.y;

		if(pJ_LinCoef != 0)
		{
			axe2 = ax*ax; aye2 = ay*ay; cze2 = cz*cz;
			twoax = 2.*ax;
			twocz = 2.*cz; tricz = 3.*cz; tricze2 = 3.*cze2;
			axe2qzz = axe2*qzz; axe2qxz = axe2*qxz; aye2qyz = aye2*qyz;
			trecze2qzz = tricze2*qzz; trecze2qxz = tricze2*qxz; tricze2qyz = tricze2*qyz;
			triczqzz = tricz*qzz;
			czqxz = cz*qxz; triczqxz = tricz*qxz;
			czqyz = cz*qyz; triczqyz = tricz*qyz; 
			axqzz = ax*qzz; axqyz = ax*qyz;
			ayqzz = ay*qzz; ayqxz = ay*qxz; ayqyz = ay*qyz;

			axqyzp2qyxp2qxy = axqyz + 2.*qyx + 2.*qxy;
			qxxpaxqxzmqyy = qxxmqyy + ax*qxz;
			qxxpaxqxzmqyymayqyz = qxxpaxqxzmqyy - ayqyz;
			qzxpaxqzz = qzx + axqzz; twoaxqzz = 2.*axqzz;
			qzypayqzz = qzy + ayqzz; twoqzypayqzz = 2.*qzy + ayqzz;
			aytwoqzypayqzz = ay*twoqzypayqzz;
			qxypayqxzpqyxpaxqyz = qxypqyx + ayqxz + axqyz;
			twoqxypayqxzptwoqyx = 2.*qxypqyx + ayqxz;
			twoqxxmqyymayqyz = 2.*(qxxmqyy - ayqyz);
			twoqxxpaxqxzmqyy = 2.*qxxpaxqxzmqyy;

		}

		int pgnAmOfVertices = pPgn->AmOfEdgePoints;
		int pgnAmOfVertices_mi_1 = pgnAmOfVertices - 1;
		int j2 = 1;
		vR1loc.z = pPgn->CoordZ;
		vR2loc.z = pPgn->CoordZ;

		double AbsRandX = radCR.AbsRandMagnitude(vert0_2d.x - pPgn->CentrPoint.x);
		double AbsRandY = radCR.AbsRandMagnitude(vert0_2d.y - pPgn->CentrPoint.y);
		double AbsRandXY = AbsRandX;
		if(AbsRandXY < AbsRandY) AbsRandXY = AbsRandY;

		TVector3d vFaceContribIB(0,0,0);
		//bool IntLineCrossesFace = true;

		for(int j=0; j<pgnAmOfVertices; j++)
		{
			if(j == pgnAmOfVertices_mi_1) j2 = 0;

			TVector2d &vPgnVertex1 = vPgnVertices[j];
			TVector2d &vPgnVertex2 = vPgnVertices[j2++]; //don't use j2 after this!

			vR1loc.x = vPgnVertex1.x; vR1loc.y = vPgnVertex1.y;
			vR2loc.x = vPgnVertex2.x; vR2loc.y = vPgnVertex2.y;

			vR1int = trFace2Int.TrPoint(vR1loc);
			vR2int = trFace2Int.TrPoint(vR2loc);
			vSegm3d = vR2int - vR1int;
			vSegm3dUnit = vSegm3d; vSegm3dUnit.Normalize();
			//vSegm3dExtNorm = vSegm3dUnit^vFaceNormIntFrame;
			vSegm3dExtNorm = vSegm3dUnit^(signZ_FaceNorm*Ez); //!!
			
			vSegmExtNorm.x = vSegm3dExtNorm.x; vSegmExtNorm.y = vSegm3dExtNorm.y;
			//vSegmExtNorm.Normalize();
			if(::fabs(vSegmExtNorm.y) < relTolPerp) continue; //to check consistency with the subsequent
			
			double signY_EdgeNorm = Sign(vSegmExtNorm.y);

			double x1 = vR1int.x, x2 = vR2int.x;
			double y1 = vR1int.y, y2 = vR2int.y;
			if(x1 > x2)
			{
				auxBuf = x1; x1 = x2; x2 = auxBuf;
				auxBuf = y1; y1 = y2; y2 = auxBuf;
			}
			else if(x1 == x2)
			{
				x1 -= AbsRandXY; //Artificial shift
			}
			if(x1 == 0.) x1 = -AbsRandXY; //Artificial shift
			if(x2 == 0.) x2 = AbsRandXY;

			//if((x1*vSegmExtNorm.x + y1*vSegmExtNorm.y) < 0.)
			//{
			//	IntLineCrossesFace = false; //to check whether this is true and whether it is necessary
			//}

			double kx = (y2 - y1)/(x2 - x1);
			double by = y1 - kx*x1;

			if(by == 0.) by = AbsRandXY; //?

			double kxe2 = kx*kx;
			double kxe2p1 = kxe2 + 1., kxe2m1 = kxe2 - 1.;
			double y1dx1 = y1/x1, y2dx2 = y2/x2;
			double invkxe2p1 = 1./kxe2p1;
			double invkxe2p1e2 = invkxe2p1*invkxe2p1;
			double x1e2 = x1*x1, x2e2 = x2*x2, y1e2 = y1*y1, y2e2 = y2*y2;
			double x1e2py1e2 = x1e2 + y1e2, x2e2py2e2 = x2e2 + y2e2, x2e2mx1e2 = x2e2 - x1e2, x2mx1 = x2 - x1;
			double x1e2py1e2dx1e2 = x1e2py1e2/x1e2, x2e2py2e2dx2e2 = x2e2py2e2/x2e2, x2e2py2e2dx1e2py1e2 = x2e2py2e2/x1e2py1e2;
			double kxy1px1 = kx*y1 + x1, kxy2px2 = kx*y2 + x2;
			double kxy1px1dby = kxy1px1/by, kxy2px2dby = kxy2px2/by;
			double axx1p2cz = ax*x1 + twocz, axx2p2cz = ax*x2 + twocz;
			double czkxe2p1 = cz*kxe2p1;

			piMult1 = 0;
			double sumAtan1 = atan(TransAtans(kxy2px2dby, -kxy1px1dby, piMult1)); //-ArcTan(kxy1px1dby) + ArcTan(kxy2px2dby)
			sumAtan1 += Pi*piMult1;
			double atan_y1dx1 = atan(y1dx1), atan_y2dx2 = atan(y2dx2);
			double log_x2e2py2e2dx2e2 = log(x2e2py2e2dx2e2), log_x1e2py1e2dx1e2 = log(x1e2py1e2dx1e2), log_x2e2py2e2dx1e2py1e2 = log(x2e2py2e2dx1e2py1e2);
			//consider storing logs and atans to minimize their re-calculation ?

			if(pJ_LinCoef == 0) //without linear term
			{
				IBxLoc = quartJz0*(2*ay*kx*x2e2mx1e2 + 2*by*invkxe2p1*(ay + ax*kx + 2*ay*kxe2)*x2mx1 - 2*by*invkxe2p1e2*(ay*by*kxe2m1 - 2*(czkxe2p1 - ax*by*kx))*sumAtan1 
					+ 2*ay*(x1e2*atan_y1dx1 - x2e2*atan_y2dx2) + axx2p2cz*x2*log_x2e2py2e2dx2e2 - axx1p2cz*x1*log_x1e2py1e2dx1e2 - by*invkxe2p1e2*(ax*by*kxe2m1 - 2*kx*(ay*by + czkxe2p1))*log_x2e2py2e2dx1e2py1e2);
				IByLoc = quartJz0*(-2*by*invkxe2p1*(ax + ay*kx)*x2mx1 - 2*by*invkxe2p1e2*(ax*by*kxe2m1 - 2*kx*(ay*by + czkxe2p1))*sumAtan1 
					+ 2*(axx1p2cz*x1*atan_y1dx1 - axx2p2cz*x2*atan_y2dx2) + ay*(x1e2*log_x1e2py1e2dx1e2 - x2e2*log_x2e2py2e2dx2e2) + by*invkxe2p1e2*(ay*by*kxe2m1 - 2*(czkxe2p1 - ax*by*kx))*log_x2e2py2e2dx1e2py1e2);
				IBzLoc = 0.25*(-2*ay*jx0*kx*x2e2mx1e2 - 2*by*invkxe2p1*(ax*(-jy0 + jx0*kx) + ay*(jx0 - jy0*kx + twojx0*kxe2))*x2mx1 
					+ 2*by*invkxe2p1e2*(ay*by*(-twojy0*kx + jx0*kxe2m1) + ax*by*(twojx0*kx + jy0*kxe2m1) - twocz*(jx0 + jy0*kx)*kxe2p1)*sumAtan1 
					+ 2*x2*(twocz*jy0 + ay*jx0*x2 + ax*jy0*x2)*atan_y2dx2 - 2*x1*(twocz*jy0 + ay*jx0*x1 + ax*jy0*x1)*atan_y1dx1 
					+ x1*(twocz*jx0 + ax*jx0*x1 - ay*jy0*x1)*log_x1e2py1e2dx1e2 - x2*(twocz*jx0 + ax*jx0*x2 - ay*jy0*x2)*log_x2e2py2e2dx2e2 
					- by*invkxe2p1e2*(ax*by*(jx0 + twojy0*kx - jx0*kxe2) + ay*by*(twojx0*kx + jy0*kxe2m1) - twocz*(jy0 - jx0*kx)*kxe2p1)*log_x2e2py2e2dx1e2py1e2);
			}
			else
			{
				double kxe2p1kx = kxe2p1*kx, kxe4 = kxe2*kxe2, invkxe2p1e3 = invkxe2p1e2*invkxe2p1, kxe2m3 = kxe2 - 3.;
				double kxe3 = kxe2*kx, kxe2p1e2 = kxe2p1*kxe2p1, bye2 = by*by;
				double kxe4m1 = kxe4 - 1, trekxe2m1 = 3.*kxe2 - 1, trekxe2p2 = 3.*kxe2 + 2.;
				double kxqxy = kx*qxy, kxe2qzy = kxe2*qzy;
				double x1e3 = x1e2*x1, x2e3 = x2e2*x2;
				double x2e3mx1e3 = x2e3 - x1e3;
				double czkxkxe2p1 = cz*kxe2p1kx, bykxqyy = by*kx*qyy, bykxqxy = by*kxqxy;
				double sixczkxkxe2p1 = 6.*czkxkxe2p1, triczkxkxe2p1 = 3.*czkxkxe2p1, triczkxe2p1 = 3.*czkxe2p1;

				IBxLoc = inv12*(invkxe2p1*(sixczkxkxe2p1*qzy + twoax*by*(kx*qzx + twoqzy + 3*kxe2qzy) + axe2*by*kx*qzz - aye2*by*kx*qzz 
					+ ay*(6*jz0*kxe2p1kx + sixczkxkxe2p1*qzz + by*(4*qzx + 6*kxe2*qzx - twoqzy*kx + 4*axqzz + 6*ax*kxe2*qzz)))*x2e2mx1e2 + 4*kx*(ay*qzxpaxqzz + ax*qzy)*x2e3mx1e3 
					+ (2*by*invkxe2p1e2*(triczkxe2p1*(kx*qzx + qzy + 2*kxe2qzy) - axe2*by*kxe2m1*qzz + aye2*by*kxe2m1*qzz + ax*(trijz0*kxe2p1kx + by*(twoqzx - twoqzx*kxe2 + 4*qzy*kx) + triczkxkxe2p1*qzz) 
					+ ay*(jz0*(3 + 9*kxe2 + 6*kxe4) + 2*by*(2*kx*qzxpaxqzz - qzy + kxe2qzy) + tricz*(1 + 3*kxe2 + 2*kxe4)*qzz)) + aytwoqzypayqzz*(3*bye2 + 3*by*kx*(x1 + x2) + kxe2*(x1e2 + x1*x2 + x2e2)))*x2mx1 
					+ 2*by*invkxe2p1e3*(6*czkxe2p1*(jz0*kxe2p1 - by*kx*qzxpaxqzz) + tricze2*kxe2p1e2*qzz + by*((twoax*by*kx*kxe2m3 - tricz*kxe4m1)*qzy 
					+ ay*(-trijz0*kxe4m1 + 2*by*kx*kxe2m3*qzxpaxqzz - triczqzz*kxe4m1)) + by*(axe2*by*qzz*trekxe2m1 - by*aytwoqzypayqzz*trekxe2m1 - twoax*(trijz0*kxe2p1kx - by*qzx*trekxe2m1)))*sumAtan1 
					+ 2*(qzy*(tricz + twoax*x1) + ay*(trijz0 + triczqzz + twoqzx*x1 + twoaxqzz*x1))*x1e2*atan_y1dx1 - 2*(qzy*(tricz + twoax*x2) + ay*(trijz0 + triczqzz + twoqzx*x2 + twoaxqzz*x2))*x2e2*atan_y2dx2 
					+ x2*(trecze2qzz + tricz*(twojz0 + qzxpaxqzz*x2) + x2*(axe2qzz*x2 - aytwoqzypayqzz*x2 + ax*(trijz0 + twoqzx*x2)))*log_x2e2py2e2dx2e2
					- x1*(trecze2qzz + tricz*(twojz0 + qzxpaxqzz*x1) + x1*(axe2qzz*x1 - aytwoqzypayqzz*x1 + ax*(trijz0 + twoqzx*x1)))*log_x1e2py1e2dx1e2 
					+ by*invkxe2p1e3*(axe2*bye2*kx*kxe2m3*qzz - aye2*bye2*kx*kxe2m3*qzz + triczkxe2p1*(twojz0*kxe2p1kx + by*(qzx - kxe2*qzx + twoqzy*kx) + czkxkxe2p1*qzz) 
					+ 2*ay*by*(trijz0*kxe2p1kx + by*(qzx - 3*kxe2*qzx + 3*qzy*kx - kxe3*qzy) + triczkxkxe2p1*qzz) 
					+ ax*by*(-trijz0*kxe4m1 - triczqzz*kxe4m1 + 2*by*(-3*kx*qzx + kxe3*qzx + qzypayqzz - 3*kxe2*qzypayqzz)))*log_x2e2py2e2dx1e2py1e2);
				IByLoc = inv12*(-(by*invkxe2p1*(axe2qzz + twoax*(qzx + kx*qzypayqzz) + ay*(2*kx*qzx + 4*qzy + 2*ayqzz + 3*kxe2*twoqzypayqzz))*x2e2mx1e2) - 2*kx*aytwoqzypayqzz*x2e3mx1e3 
					+ 2*by*invkxe2p1e2*(-triczkxe2p1*(qzx + kx*qzy) + 2*axe2*by*kx*qzz - 2*aye2*by*kx*qzz - ay*(trijz0*kxe2p1kx + by*(twoqzx - twoqzx*kxe2 + 4*qzy*kx) + triczkxkxe2p1*qzz) 
					- ax*(trijz0*kxe2p1 + triczkxe2p1*qzz - 2*by*(twoqzx*kx - qzy - ayqzz + kxe2*qzypayqzz)))*x2mx1 
					+ 2*by*invkxe2p1e3*(axe2*bye2*kx*kxe2m3*qzz - aye2*bye2*kx*kxe2m3*qzz + triczkxe2p1*(twojz0*kxe2p1kx + by*(qzx - kxe2*qzx + twoqzy*kx) + czkxkxe2p1*qzz) 
					+ 2*ay*by*(trijz0*kxe2p1kx + by*(qzx - 3*qzx*kxe2 + 3*qzy*kx - kxe3*qzy) + triczkxkxe2p1*qzz) 
					+ ax*by*(-trijz0*kxe4m1 - triczqzz*kxe4m1 + 2*by*(-3*qzx*kx + kxe3*qzx + qzypayqzz - 3*kxe2*qzypayqzz)))*sumAtan1 
					+ 2*x1*(trecze2qzz + tricz*(twojz0 + qzxpaxqzz*x1) + x1*(axe2qzz*x1 - aytwoqzypayqzz*x1 + ax*(trijz0 + twoqzx*x1)))*atan_y1dx1 
					- 2*x2*(trecze2qzz + tricz*(twojz0 + qzxpaxqzz*x2) + x2*(axe2qzz*x2 - aytwoqzypayqzz*x2 + ax*(trijz0 + twoqzx*x2)))*atan_y2dx2 
					+ (qzy*(tricz + twoax*x1) + ay*(trijz0 + triczqzz + twoqzx*x1 + twoaxqzz*x1))*x1e2*log_x1e2py1e2dx1e2 
					- (qzy*(tricz + twoax*x2) + ay*(trijz0 + triczqzz + twoqzx*x2 + twoaxqzz*x2))*x2e2*log_x2e2py2e2dx2e2
					+ by*invkxe2p1e3*(-triczkxe2p1*(twojz0*kxe2p1 + by*(-2*kx*qzxpaxqzz + qzy - kxe2qzy)) - tricze2*kxe2p1e2*qzz 
					+ ay*by*(trijz0*kxe4m1 - 2*by*(-3*kx*qzxpaxqzz + kxe3*qzxpaxqzz + qzy - 3*kxe2qzy) + triczqzz*kxe4m1) + aye2*bye2*qzz*trekxe2m1 
					+ ax*by*(6*jz0*kxe2p1kx + by*(6*qzy*kx - twoqzy*kxe3 - axqzz*trekxe2m1 - twoqzx*trekxe2m1)))*log_x2e2py2e2dx1e2py1e2);
				IBzLoc = inv12*(-(invkxe2p1*(sixczkxkxe2p1*qxy + twoax*by*(kx*qxxmqyy + 2*qxy + 3*kxe2*qxy - qyx) + axe2*by*(kx*qxz - qyz) + aye2*by*(kx*qxz - qyz)*trekxe2p2 
					+ ay*(6*jx0*kxe2p1kx + sixczkxkxe2p1*qxz + by*(6*kxe3*qxy + 4*ax*qxz + 6*kxe2*(ax*qxz - qyy) - 4*qyy - 2*kx*(-2*qxy + qyx + axqyz) + twoqxx*trekxe2p2)))*x2e2mx1e2) 
					- kx*(4*ax*qxy + 2*ay*(2*qxxmqyy + kxqxy + twoax*qxz) + aye2*(kx*qxz - 2*qyz))*x2e3mx1e3 
					- by*invkxe2p1e2*(aye2*by*(qxz + 8*kxe2*qxz + 3*kxe4*qxz - 4*kx*qyz) + 2*ay*(-trijy0*kxe2p1kx + jx0*(3 + 9*kxe2 + 6*kxe4) + 4*qxx*by*kx + by*qxy + 8*by*kxe2*qxy + 3*by*kxe4*qxy 
					+ triczqxz + 4*ax*by*kx*qxz + 9*czqxz*kxe2 + 6*czqxz*kxe4 - 2*by*qyx + 2*by*kxe2*qyx - 4*bykxqyy - twoax*by*qyz - triczqyz*kx + twoax*by*kxe2*qyz - triczqyz*kxe3) 
					+ 2*(triczkxe2p1*(kx*qxxmqyy + qxy + 2*kxe2*qxy - qyx) + axe2*by*(qxz - kxe2*qxz + 2*kx*qyz) 
					+ ax*(-trijy0*kxe2p1 + trijx0*kxe2p1kx + twoqxx*by - twoqxx*by*kxe2 + 4*bykxqxy + triczqxz*kx + triczqxz*kxe3 + 4*by*kx*qyx - 2*by*qyy + 2*by*kxe2*qyy - triczqyz - triczqyz*kxe2)))*x2mx1 
					- 2*by*invkxe2p1e3*(ax*by*(-6*jx0*kxe2p1kx - trijy0*kxe4m1 + by*(-3*axqyzp2qyxp2qxy*kx + axqyzp2qyxp2qxy*kxe3 + (-2 + 6*kxe2)*qxx - ax*qxz + 3*kxe2*(ax*qxz - 2*qyy) + 2*qyy)) 
					+ tricze2*kxe2p1e2*(qxz + kx*qyz) + aye2*bye2*(qxz - 3*kxe2*qxz - kx*kxe2m3*qyz) + ay*by*(6*jy0*kxe2p1kx - trijx0*kxe4m1 - 6*qxx*by*kx + twoqxx*by*kxe3 + 2*by*qxy - 6*by*kxe2*qxy 
					+ triczqxz - 6*ax*by*kx*qxz + twoax*by*kxe3*qxz - triczqxz*kxe4 + 2*by*qyx - 6*by*kxe2*qyx + 6*bykxqyy - 2*by*kxe3*qyy + twoax*by*qyz + 6*czqyz*kx - 6*ax*by*kxe2*qyz + 6*czqyz*kxe3) 
					+ triczkxe2p1*(twojx0*kxe2p1 + twojy0*kxe2p1kx + by*(-2*kx*qxxpaxqxzmqyy + qxypqyx + axqyz - kxe2*(qxypqyx + axqyz))))*sumAtan1  
					+ 2*x1*(-tricze2qyz - tricz*(twojy0 + qxypayqxzpqyxpaxqyz*x1) + x1*(aye2qyz*x1 - ax*(trijy0 + axqyzp2qyxp2qxy*x1) - ay*(trijx0 + twoqxxpaxqxzmqyy*x1)))*atan_y1dx1 
					- 2*x2*(-tricze2qyz - tricz*(twojy0 + qxypayqxzpqyxpaxqyz*x2) + x2*(aye2qyz*x2 - ax*(trijy0 + axqyzp2qyxp2qxy*x2) - ay*(trijx0 + twoqxxpaxqxzmqyy*x2)))*atan_y2dx2 
					+ x1*(trecze2qxz + tricz*(twojx0 + qxxpaxqxzmqyymayqyz*x1) + x1*(axe2qxz*x1 - ay*(trijy0 + twoqxypayqxzptwoqyx*x1) + ax*(trijx0 + twoqxxmqyymayqyz*x1)))*log_x1e2py1e2dx1e2 
					- x2*(trecze2qxz + tricz*(twojx0 + qxxpaxqxzmqyymayqyz*x2) + x2*(axe2qxz*x2 - ay*(trijy0 + twoqxypayqxzptwoqyx*x2) + ax*(trijx0 + twoqxxmqyymayqyz*x2)))*log_x2e2py2e2dx2e2
					- by*invkxe2p1e3*(axe2*bye2*(-3*kx*qxz + kxe3*qxz + qyz - 3*kxe2*qyz) - aye2*bye2*(-3*kx*qxz + kxe3*qxz + qyz - 3*kxe2*qyz) 
					+ triczkxe2p1*(-twojy0*kxe2p1 + twojx0*kxe2p1kx + by*qxx - by*kxe2*qxx + 2*bykxqxy + czqxz*kx + czqxz*kxe3 + 2*by*kx*qyx - by*qyy + by*kxe2*qyy - czqyz - czqyz*kxe2) 
					+ ax*by*(6*jy0*kxe2p1kx - trijx0*kxe4m1 - 6*qxx*by*kx + twoqxx*by*kxe3 + 2*by*qxy - 6*by*kxe2*qxy + 2*ay*by*qxz + triczqxz - 6*ay*by*kxe2*qxz - triczqxz*kxe4 + 2*by*qyx - 6*by*kxe2*qyx 
					+ 6*bykxqyy - 2*by*kxe3*qyy + 6*ay*by*kx*qyz + 6*czqyz*kx - 2*ay*by*kxe3*qyz + 6*czqyz*kxe3) + ay*by*(6*jx0*kxe2p1kx + trijy0*kxe4m1 + twoqxx*by - 6*qxx*by*kxe2 
					+ 6*bykxqxy - 2*by*kxe3*qxy + 6*czqxz*kx + 6*czqxz*kxe3 + 6*by*kx*qyx - 2*by*kxe3*qyx - 2*by*qyy + 6*by*kxe2*qyy - triczqyz + triczqyz*kxe4))*log_x2e2py2e2dx1e2py1e2);
			}
			vFaceContribIB.x += signY_EdgeNorm*IBxLoc;
			vFaceContribIB.y += signY_EdgeNorm*IByLoc;
			vFaceContribIB.z += signY_EdgeNorm*IBzLoc;
		}

		vResLocIB += signZ_FaceNorm*vFaceContribIB;
	}
	// Biot-Savart constant: B = (mu_0/4*pi) * integral(J x r / r^3) dV
	// Radia now uses SI units (meters) internally, matching ELF.
	// Note: This uses 2e-7 (= 2 * mu_0/(4*pi)) for specific integral formula
	const double ConstForJ = 2.0e-7;
	vResLocIB *= ConstForJ;
	TVector3d vResIB = trIntAxis2Ez.TrVectField_inv(vResLocIB);
	
	if(pField->FieldKey.Ib_) pField->Ib += vResIB;
	if(pField->FieldKey.Ih_) pField->Ih += vResIB;
}

//-------------------------------------------------------------------------


//-------------------------------------------------------------------------

// radTPolyhedron::Dump / DumpPureObjInfo REMOVED (Phase B2b, 2026-04-15)

//-------------------------------------------------------------------------

// DumpBin / DumpBin_Polyhedron / DumpBinParse_Polyhedron REMOVED (Phase B2c, 2026-04-15)

//-------------------------------------------------------------------------

void radTPolyhedron::DefineRelAndAbsTol(double* RelAbsTol)
{
	double RelZeroToler = 1.E-09;
	//RelZeroToler = 500.*max(RelZeroToler, radCR.RelRand);

	double MaxVal = RelZeroToler;
	//if(radCR.RelRand < radCR.RelRand) MaxVal = radCR.RelRand;
	if(MaxVal < radCR.RelRand) MaxVal = radCR.RelRand;
	//RelZeroToler = 500.*MaxVal;
	RelZeroToler = 100.*MaxVal; //OC291003

	RelAbsTol[0] = RelZeroToler;
	if(VectHandlePgnAndTrans.empty()) return;

	radTPolygon* PgnPtr = VectHandlePgnAndTrans[0].PgnHndl.rep;
	TVector2d& vpBuf = (PgnPtr->EdgePointsVector)[0];
	TVector3d aVertexPoint = TVector3d(vpBuf.x, vpBuf.y, PgnPtr->CoordZ);
	aVertexPoint = VectHandlePgnAndTrans[0].TransHndl.rep->TrPoint(aVertexPoint);
	TVector3d VectToCenter = CentrPoint - aVertexPoint;
	RelAbsTol[1] = RelZeroToler*NormAbs(VectToCenter);
	//RelAbsTol[1] = RelZeroToler*NormAbs(aVertexPoint); //OC090908
}

//-------------------------------------------------------------------------

// radTPolyhedron::CutItself REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::FindIntersectionWithFace REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::IntrsctOfTwoLines REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::CheckIfTwoPointAlreadyMapped REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::SetUpUpperAndLowerPolygon REMOVED (Phase C, 2026-04-16)


//-------------------------------------------------------------------------

// radTPolyhedron::DetermineNewFaceAndTrans REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::FillInNewHandlePgnAndTransFrom3d REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::KsFromSizeToNumb REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::SubdivideItselfByParPlanes REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::SubdivideItselfByOneSetOfParPlanes REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::DeterminePointsOnCuttingPlanes REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::FindLowestAndUppestVertices REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

int radTPolyhedron::CheckForSpecialShapes(radTVectHandlePgnAndTrans& VectHandlePgnAndTrans, radThg& In_hg, double* RelAbsTol)
{// Currently, this distinguishes only parallelepipeds
	int NumberOfFaces = (int)VectHandlePgnAndTrans.size();
	if(NumberOfFaces != 6) return 0;

	double CPx, CPy, CPz, DimX, DimY, DimZ;
	TVector3d aNormal(0.,0.,1.);
	TVector3d AllNormals[6];
	double xSt, xFi, ySt, yFi, zSt, zFi;
	double RelTol = RelAbsTol[0], AbsTol = RelAbsTol[1];
	TVector3d Zero(0.,0.,0.);

	short InternalFacesMap[] = { 0,0,0,0,0,0 };
	short OrientationIsNotGood = 0;
	for(int k=0; k<NumberOfFaces; k++)
	{
		radTHandlePgnAndTrans& HandlePgnAndTrans = VectHandlePgnAndTrans[k];
		radTPolygon* PgnPtr = HandlePgnAndTrans.PgnHndl.rep;
		radTrans* TransPtr = HandlePgnAndTrans.TransHndl.rep;
		if(PgnPtr->AmOfEdgePoints != 4) return 0;

		TVector3d CurrentNormal = TransPtr->TrBiPoint(aNormal);

		double AbsNx = fabs(CurrentNormal.x), AbsNy = fabs(CurrentNormal.y), AbsNz = fabs(CurrentNormal.z);

		if((AbsNy < RelTol) && (AbsNz < RelTol) && (fabs(CurrentNormal.x + 1.) < RelTol))
		{
			TVector2d& AnEdgePo2d = PgnPtr->EdgePointsVector[0];
			TVector3d AnEdgePo3d(AnEdgePo2d.x, AnEdgePo2d.y, PgnPtr->CoordZ);
			AnEdgePo3d = TransPtr->TrPoint(AnEdgePo3d);
			xSt = AnEdgePo3d.x; 
			InternalFacesMap[0] = HandlePgnAndTrans.FaceIsInternalAfterCut;
		}
		else if((AbsNy < RelTol) && (AbsNz < RelTol) && (fabs(CurrentNormal.x - 1.) < RelTol))
		{
			TVector2d& AnEdgePo2d = PgnPtr->EdgePointsVector[0];
			TVector3d AnEdgePo3d(AnEdgePo2d.x, AnEdgePo2d.y, PgnPtr->CoordZ);
			AnEdgePo3d = TransPtr->TrPoint(AnEdgePo3d);
			xFi = AnEdgePo3d.x;
			InternalFacesMap[1] = HandlePgnAndTrans.FaceIsInternalAfterCut;
		}
		else if((AbsNx < RelTol) && (AbsNz < RelTol) && (fabs(CurrentNormal.y + 1.) < RelTol))
		{
			TVector2d& AnEdgePo2d = PgnPtr->EdgePointsVector[0];
			TVector3d AnEdgePo3d(AnEdgePo2d.x, AnEdgePo2d.y, PgnPtr->CoordZ);
			AnEdgePo3d = TransPtr->TrPoint(AnEdgePo3d);
			ySt = AnEdgePo3d.y;
			InternalFacesMap[2] = HandlePgnAndTrans.FaceIsInternalAfterCut;
		}
		else if((AbsNx < RelTol) && (AbsNz < RelTol) && (fabs(CurrentNormal.y - 1.) < RelTol))
		{
			TVector2d& AnEdgePo2d = PgnPtr->EdgePointsVector[0];
			TVector3d AnEdgePo3d(AnEdgePo2d.x, AnEdgePo2d.y, PgnPtr->CoordZ);
			AnEdgePo3d = TransPtr->TrPoint(AnEdgePo3d);
			yFi = AnEdgePo3d.y;
			InternalFacesMap[3] = HandlePgnAndTrans.FaceIsInternalAfterCut;
		}
		else if((AbsNx < RelTol) && (AbsNy < RelTol) && (fabs(CurrentNormal.z + 1.) < RelTol))
		{
			TVector2d& AnEdgePo2d = PgnPtr->EdgePointsVector[0];
			TVector3d AnEdgePo3d(AnEdgePo2d.x, AnEdgePo2d.y, PgnPtr->CoordZ);
			AnEdgePo3d = TransPtr->TrPoint(AnEdgePo3d);
			zSt = AnEdgePo3d.z;
			InternalFacesMap[4] = HandlePgnAndTrans.FaceIsInternalAfterCut;
		}
		else if((AbsNx < RelTol) && (AbsNy < RelTol) && (fabs(CurrentNormal.z - 1.) < RelTol))
		{
			TVector2d& AnEdgePo2d = PgnPtr->EdgePointsVector[0];
			TVector3d AnEdgePo3d(AnEdgePo2d.x, AnEdgePo2d.y, PgnPtr->CoordZ);
			AnEdgePo3d = TransPtr->TrPoint(AnEdgePo3d);
			zFi = AnEdgePo3d.z;
			InternalFacesMap[5] = HandlePgnAndTrans.FaceIsInternalAfterCut;
		}
		else OrientationIsNotGood = 1;
		AllNormals[k] = CurrentNormal;
	}

	radTrans TransForRecMag;
	if(!OrientationIsNotGood)
	{
		CPx = 0.5*(xSt+xFi); DimX = xFi-xSt;
		CPy = 0.5*(ySt+yFi); DimY = yFi-ySt;
		CPz = 0.5*(zSt+zFi); DimZ = zFi-zSt;
	}
	else
	{
		int NoOfFacePerpToFirst=-1;
		TVector3d* TestNormalPtr = AllNormals;
		for(int j=0; j<NumberOfFaces; j++)
		{
			TVector3d* LocPtr = TestNormalPtr;
			for(int jj=j+1; jj<NumberOfFaces; jj++)
			{
				double ScalProd = (*TestNormalPtr)*(*(++LocPtr));
				short NormalsArePerp = (fabs(ScalProd) < RelTol);
				if(j==0) if(NoOfFacePerpToFirst < 0) if(NormalsArePerp) NoOfFacePerpToFirst = jj;

				if(!(NormalsArePerp || (fabs(ScalProd+1.) < RelTol))) return 0;
			}
			TestNormalPtr++;
		}
		if(NoOfFacePerpToFirst < 0) return 0;
		radTHandlePgnAndTrans& FirstHandlePgnAndTrans = VectHandlePgnAndTrans[0];
		radTPolygon* FirstPgnPtr = FirstHandlePgnAndTrans.PgnHndl.rep;
		radTrans* FirstTransPtr = FirstHandlePgnAndTrans.TransHndl.rep;

		radTHandlePgnAndTrans& NextHandlePgnAndTrans = VectHandlePgnAndTrans[NoOfFacePerpToFirst];
		radTPolygon* NextPgnPtr = NextHandlePgnAndTrans.PgnHndl.rep;
		radTrans* NextTransPtr = NextHandlePgnAndTrans.TransHndl.rep;

		TVector2d& p0 = FirstPgnPtr->EdgePointsVector[0];
		TVector2d& p1 = FirstPgnPtr->EdgePointsVector[1];
		TVector2d& p2 = FirstPgnPtr->EdgePointsVector[2];
		TVector2d& p3 = FirstPgnPtr->EdgePointsVector[3];

		TVector2d v0 = p1-p0;
		double InvNorm = 1./sqrt(v0.x*v0.x + v0.y*v0.y);
		v0 = InvNorm*v0;
		TVector3d St0(v0.x, v0.y, 0.);
		TVector3d St1(-v0.y, v0.x, 0.);
		TVector3d St2(0., 0., 1.);
		TMatrix3d M_fr_v0_to_i(St0, St1, St2);
		radTrans R(M_fr_v0_to_i, Zero, 1., 1., 2); // From v0 to i

		double FirstZ = FirstPgnPtr->CoordZ;
		TVector3d P0(p0.x, p0.y, FirstZ), P1(p1.x, p1.y, FirstZ), P2(p2.x, p2.y, FirstZ), P3(p3.x, p3.y, FirstZ);
		P0 = R.TrPoint(P0); P1 = R.TrPoint(P1); P2 = R.TrPoint(P2); P3 = R.TrPoint(P3);

		CPx = 0.5*(P0.x + P1.x), CPy = 0.5*(P1.y + P2.y);
		DimX = fabs(P1.x - P0.x), DimY = fabs(P2.y - P1.y);

		R.Invert();
		TrProduct(FirstTransPtr, &R, TransForRecMag);

		TVector3d OneMorePoint;
		double NextZ = NextPgnPtr->CoordZ;
		for(int kk=0; kk<3; kk++)
		{
			TVector2d& p = NextPgnPtr->EdgePointsVector[kk];
			TVector3d Pt(p.x, p.y, NextZ);
			Pt = NextTransPtr->TrPoint(Pt); Pt = TransForRecMag.TrPoint_inv(Pt);

			if((fabs(Pt.z-P0.z)>AbsTol) && (fabs(Pt.z-P1.z)>AbsTol) && (fabs(Pt.z-P2.z)>AbsTol) && (fabs(Pt.z-P1.z)>AbsTol))
			{
				OneMorePoint = Pt; break;
			}
		}
		CPz = 0.5*(P0.z + OneMorePoint.z);
		DimZ = fabs(P1.z - OneMorePoint.z);

		for(int i=0; i<NumberOfFaces; i++)
		{
			TVector3d NormalInLocFrame = TransForRecMag.TrBiPoint_inv(AllNormals[i]);
			short FaceNoForRecMag;
			if(fabs(NormalInLocFrame.x + 1.) < RelTol) FaceNoForRecMag = 0;
			else if(fabs(NormalInLocFrame.x - 1.) < RelTol) FaceNoForRecMag = 1;
			else if(fabs(NormalInLocFrame.y + 1.) < RelTol) FaceNoForRecMag = 2;
			else if(fabs(NormalInLocFrame.y - 1.) < RelTol) FaceNoForRecMag = 3;
			else if(fabs(NormalInLocFrame.z + 1.) < RelTol) FaceNoForRecMag = 4;
			else if(fabs(NormalInLocFrame.z - 1.) < RelTol) FaceNoForRecMag = 5;
			else return 0;
			InternalFacesMap[FaceNoForRecMag] = VectHandlePgnAndTrans[i].FaceIsInternalAfterCut;
		}
	}

	TVector3d CPoiVect(CPx, CPy, CPz), DimsVect(DimX, DimY, DimZ);
	TVector3d MagnForRecMag = OrientationIsNotGood? TransForRecMag.TrVectField_inv(Magn) : Magn;

	TVector3d J_ForRecMag = Zero;
	short J_IsNotZeroLoc = 0;
	if(J_IsNotZero)
	{
		J_ForRecMag = OrientationIsNotGood? TransForRecMag.TrVectField_inv(J) : J; //??
		J_IsNotZeroLoc = 1;
	}

	//short J_IsNotZero = 0;
	//radTRecMag* NewRecMagPtr = new radTRecMag(CPoiVect, DimsVect, MagnForRecMag, Zero, MaterHandle, J_IsNotZero);
	radTRecMag* NewRecMagPtr = new radTRecMag(CPoiVect, DimsVect, MagnForRecMag, J_ForRecMag, MaterHandle, J_IsNotZeroLoc);
	if(NewRecMagPtr==0) return 0;
	NewRecMagPtr->J_IsNotZero = 0;
	NewRecMagPtr->SetFacesInternalAfterCut(InternalFacesMap);

	if(OrientationIsNotGood)
	{
		radThg hTrans(new radTrans(TransForRecMag));
		NewRecMagPtr->AddTransform(1, hTrans);
		NewRecMagPtr->ConsiderOnlyWithTrans = 1;
	}

	radThg hRecMag(NewRecMagPtr);
	In_hg = hRecMag;
	return 1;
}

//-------------------------------------------------------------------------

double radTPolyhedron::Volume()
{
	const double kMax = 1.E+10;
	double VolSum = 0.;

	radTPolygon* PgnPtr;
	radTrans* TransPtr;
	radTHandlePgnAndTrans* HandlePgnAndTransPtr;

	for(int i=0; i<AmOfFaces; i++)
	{
		HandlePgnAndTransPtr = &(VectHandlePgnAndTrans[i]);
		PgnPtr = HandlePgnAndTransPtr->PgnHndl.rep;
		TransPtr = HandlePgnAndTransPtr->TransHndl.rep;

		TVector2d p1 = PgnPtr->EdgePointsVector[0], p2 = PgnPtr->EdgePointsVector[1], p3 = PgnPtr->EdgePointsVector[PgnPtr->IndexOfGoodThirdPoint()];
		TVector3d P1 = TVector3d(p1.x, p1.y, PgnPtr->CoordZ); P1 = TransPtr->TrPoint(P1);
		TVector3d P2 = TVector3d(p2.x, p2.y, PgnPtr->CoordZ); P2 = TransPtr->TrPoint(P2);
		TVector3d P3 = TVector3d(p3.x, p3.y, PgnPtr->CoordZ); P3 = TransPtr->TrPoint(P3);

		double x1 = P1.x, x2 = P2.x, x3 = P3.x;
		double y1 = P1.y, y2 = P2.y, y3 = P3.y;
		double z1 = P1.z, z2 = P2.z, z3 = P3.z;

		double x2mx1 = x2 - x1, y2my1 = y2 - y1;

		double BufDenom = (x2-x3)*y1 + (x3-x1)*y2 - x2mx1*y3;
		double AbsBufDenom_kMax = fabs(BufDenom)*kMax;
		double aSigBuf = y3*(z1-z2)+y1*(z2-z3)+y2*(z3-z1), aSig;
		double bSigBuf = x3*(z2-z1)+x2*(z1-z3)+x1*(z3-z2), bSig;

		if((fabs(aSigBuf) < AbsBufDenom_kMax) && (fabs(bSigBuf) < AbsBufDenom_kMax))
		{
			aSig = aSigBuf/BufDenom; bSig = bSigBuf/BufDenom;
			double cSig = z1 - aSig*x1 - bSig*y1;

			double xk1 = x1, yk1 = y1, xk2 = x2, yk2 = y2;
			for(int k=1; k<=PgnPtr->AmOfEdgePoints; k++)
			{
				if(k>1)
				{
					p3 = PgnPtr->EdgePointsVector[(k==PgnPtr->AmOfEdgePoints)? 0 : k];
					P3 = TVector3d(p3.x, p3.y, PgnPtr->CoordZ); P3 = TransPtr->TrPoint(P3);

					xk2 = P3.x; yk2 = P3.y;
					x2mx1 = xk2 - xk1; y2my1 = yk2 - yk1;
				}
				if(fabs(y2my1) < fabs(x2mx1)*kMax)
				{
					double aSigK = y2my1/x2mx1, bSigK = yk1 - aSigK*xk1;
					double xk1pxk2 = xk1 + xk2;
					double xk1e2pxk1xk2pxk2e2 = xk1*xk1 + xk1*xk2 + xk2*xk2;

					VolSum += (x2mx1/6.)*(3.*bSigK*(2.*cSig + aSig*xk1pxk2) + aSigK*(3.*cSig*xk1pxk2 + 2.*aSig*xk1e2pxk1xk2pxk2e2) 
							+ bSig*(3.*bSigK*(bSigK + aSigK*xk1pxk2) + aSigK*aSigK*xk1e2pxk1xk2pxk2e2));
				}
				xk1 = xk2; yk1 = yk2;
			}
		}
	}
	return fabs(VolSum);
}

//-------------------------------------------------------------------------

// radTPolyhedron::EstimateSize REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::SubdivideItselfByEllipticCylinder REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::SubdivideByEllipses REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::FindLocalEllipticCoord REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::SubdivideItselfOverAzimuth REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::FindEdgePointsOverEllipseSet0 REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::FindEdgePointsOverEllipseSet REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::FindEdgePointsOverPhiAndAxForCylSubd REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolyhedron::SubdivideItself REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

void radTPolyhedron::VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique)
{
	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	double AbsTol = RelAbsTol[1];
	double AbsTolE2 = AbsTol*AbsTol;

	radTVectorOfVector3d LocPts;

	for(radTVectHandlePgnAndTrans::iterator FaceIter = VectHandlePgnAndTrans.begin(); FaceIter != VectHandlePgnAndTrans.end(); ++FaceIter)
	{
		radTPolygon* PgnPtr = (*FaceIter).PgnHndl.rep;
		radTrans* TransPtr = (*FaceIter).TransHndl.rep;
		double LocZ = PgnPtr->CoordZ;
		for(radTVect2dVect::iterator PointIter = PgnPtr->EdgePointsVector.begin(); PointIter != PgnPtr->EdgePointsVector.end(); ++PointIter)
		{
			TVector2d& p2d = *PointIter;
			TVector3d p3d = TVector3d(p2d.x, p2d.y, LocZ);
			p3d = TransPtr->TrPoint(p3d);
			//p3d = TransPtr->TrPoint(p3d) + CentrPoint; //OC090908

			bool ThisPtAlreadyExists = false;
			if(EnsureUnique)
			{
				for(const auto& pt : LocPts)
				{
					TVector3d dp = pt - p3d;
					if((dp.x*dp.x + dp.y*dp.y + dp.z*dp.z) <= AbsTolE2)
					{
						ThisPtAlreadyExists = true;
						break;
					}
				}
			}
			if(!ThisPtAlreadyExists) LocPts.push_back(p3d);
		}
	}
	int SizeLocPts = (int)LocPts.size();
	if(SizeLocPts == 0) return;

	for(int j=0; j<SizeLocPts; j++) OutVect.push_back(LocPts[j]);

	LocPts.erase(LocPts.begin(), LocPts.end());
}

//-------------------------------------------------------------------------

void radTPolyhedron::Push_backCenterPointAndField(radTFieldKey* pFieldKey, radTVectPairOfVect3d* pVectPairOfVect3d, radTrans* pBaseTrans, radTg3d* g3dSrcPtr, radTApplication* pAppl)
{// Attention: this assumes no more than one transformation with mult. no more than 1 !!!
 //to move to base class?
	radTrans bufTrans, *pTrans=0;
	TVector3d cenPointInLabFr, vZero(0.,0.,0.);
	GetTrfAndCenPointInLabFrame(pBaseTrans, bufTrans, pTrans, cenPointInLabFr);
	radTPairOfVect3d Pair(cenPointInLabFr, vZero);

	if(pFieldKey->M_)
	{
		TVector3d resM = Magn;
		if((pM_LinCoef != 0) && (mLinTreat != 0)) //treat Lin. as absolute
		{
			resM += ((*pM_LinCoef)*CentrPoint); //CentrPoint in Loc. Frame - to check!
		}
		Pair.V2 = (pTrans == 0)? resM : pTrans->TrVectField(resM);
	}
	else if(pFieldKey->J_)
	{
		TVector3d resJ = J;
		if((pJ_LinCoef != 0) && (mLinTreat != 0)) //treat Lin. as absolute
		{
			resJ += ((*pJ_LinCoef)*CentrPoint); //CentrPoint in Loc. Frame - to check!
		}
		Pair.V2 = (pTrans == 0)? resJ : pTrans->TrVectField(resJ);
	}
	else
	{
		radTCompCriterium CompCriterium;
		radTField Field(*pFieldKey, CompCriterium, cenPointInLabFr, vZero, vZero, vZero, vZero, 0.);
		g3dSrcPtr->B_genComp(&Field);
		Pair.V2 = (pFieldKey->B_)? Field.B : ((pFieldKey->H_)? Field.H : ((pFieldKey->A_)? Field.A : vZero));
	}
	pVectPairOfVect3d->push_back(Pair);
}

//-------------------------------------------------------------------------

void radTPolyhedron::AttemptToCreateConvexPolyhedronFromTwoBaseFaces(const radTHandlePgnAndTrans& inHandleBasePgnAndTrf1, const radTHandlePgnAndTrans& inHandleBasePgnAndTrf2)
{//sets SomethingIsWrong = 1 in case of any problem

	radTPolygon *pPgn1 = inHandleBasePgnAndTrf1.PgnHndl.rep;
	radTrans *pTrf1 = inHandleBasePgnAndTrf1.TransHndl.rep;

	radTPolygon *pPgn2 = inHandleBasePgnAndTrf2.PgnHndl.rep;
	radTrans *pTrf2 = inHandleBasePgnAndTrf2.TransHndl.rep;

	double zc1 = pPgn1->CoordZ, zc2 = pPgn2->CoordZ;
	TVector2d vCenPoint1 = pPgn1->CentrPoint, vCenPoint2 = pPgn2->CentrPoint;

	TVector3d vNormZ(0,0,1);
	TVector3d vNorm1 = pTrf1->TrBiPoint(vNormZ), vNorm2 = pTrf2->TrBiPoint(vNormZ);
	TVector3d vCenPoint1_3d(vCenPoint1.x, vCenPoint1.y, zc1), vCenPoint2_3d(vCenPoint2.x, vCenPoint2.y, zc2);
	vCenPoint1_3d = pTrf1->TrPoint(vCenPoint1_3d); 
	vCenPoint2_3d = pTrf2->TrPoint(vCenPoint2_3d);

	double RelZeroToler = 1.E-09;
	RelZeroToler = 10.*((RelZeroToler>radCR.RelRand)? RelZeroToler : radCR.RelRand);
	double AbsTol = RelZeroToler*0.5*(pPgn1->EstimateTypSize() + pPgn2->EstimateTypSize());
	double arTol[] = {RelZeroToler, AbsTol};

	if((!CheckIfAllPolygonVerticesAreOnOneSideOfPlane(inHandleBasePgnAndTrf2, vCenPoint1_3d, vNorm1, AbsTol)) ||
	   (!CheckIfAllPolygonVerticesAreOnOneSideOfPlane(inHandleBasePgnAndTrf1, vCenPoint2_3d, vNorm2, AbsTol)))
	{
		SomethingIsWrong = 1;
		radTSend::ErrorMessage("Radia::Error110"); 
		return;
	}

	//generating mantle (side faces)
	//generate array of unique vertex points and indexes of these points for two bases faces
	vector<TVector3d> vectVertexPoints;
	vector<vector<int> > vectIndAllFaces;
	vector<int> vectIndFacePgn;

	CollectAndMapUniquePolygonPoints(inHandleBasePgnAndTrf1, vectVertexPoints, vectIndFacePgn, AbsTol);
	vectIndAllFaces.push_back(vectIndFacePgn);
	vectIndFacePgn.erase(vectIndFacePgn.begin(), vectIndFacePgn.end());

	CollectAndMapUniquePolygonPoints(inHandleBasePgnAndTrf2, vectVertexPoints, vectIndFacePgn, AbsTol);
	vectIndAllFaces.push_back(vectIndFacePgn);
	vectIndFacePgn.erase(vectIndFacePgn.begin(), vectIndFacePgn.end());

	GenerateSideFacesContainingSegmentsOfBaseFace(vectVertexPoints, vectIndAllFaces, 0, arTol);
	GenerateSideFacesContainingSegmentsOfBaseFace(vectVertexPoints, vectIndAllFaces, 1, arTol);

	int lenArVertexPoints = (int)vectVertexPoints.size();
	TVector3d *arVertexPoints = CAuxParse::Vect2Ar(vectVertexPoints);

	AmOfFaces = (int)vectIndAllFaces.size();
	int *arIndAllFacesLengths = 0;
	int **arIndAllFaces = CAuxParse::Vect2Ar2D(vectIndAllFaces, arIndAllFacesLengths);

	FillInVectHandlePgnAndTrans(arVertexPoints, lenArVertexPoints, arIndAllFaces, arIndAllFacesLengths);
	if(SomethingIsWrong) return;
	DefineCentrPoint(arVertexPoints, lenArVertexPoints);

	if(arVertexPoints != 0) delete[] arVertexPoints;
	if(arIndAllFacesLengths != 0) delete[] arIndAllFacesLengths;
	if(arIndAllFaces != 0)
	{
		for(int i=0; i<AmOfFaces; i++)
		{
			if(arIndAllFaces[i] != 0) delete[] arIndAllFaces[i];
		}
		delete[] arIndAllFaces;
	}
	vectVertexPoints.erase(vectVertexPoints.begin(), vectVertexPoints.end());
	vectIndAllFaces.erase(vectIndAllFaces.begin(), vectIndAllFaces.end());
	vectIndFacePgn.erase(vectIndFacePgn.begin(), vectIndFacePgn.end());
}

//-------------------------------------------------------------------------

void radTPolyhedron::GenerateSideFacesContainingSegmentsOfBaseFace(const vector<TVector3d>& vectVertexPoints, vector<vector<int> >& vectIndAllFaces, int indBaseFace, double* arTol)
{
	int numExistingFaces = (int)vectIndAllFaces.size();
	if((indBaseFace < 0) || (indBaseFace >= numExistingFaces)) return;

	int indOtherBaseFace = (indBaseFace == 0)? 1 : 0;
	if(indOtherBaseFace >= numExistingFaces) return;

	vector<int> vectBaseFace = vectIndAllFaces[indBaseFace];
	int numPointsInBaseFace = (int)vectBaseFace.size();

	vector<int> vectOtherBaseFace = vectIndAllFaces[indOtherBaseFace];
	int numPointsInOtherBaseFace = (int)vectOtherBaseFace.size();

	int totNumVertexPoints = (int)vectVertexPoints.size();
	double AbsTol = arTol[1];

	//loop over segments of base face
	vector<int> vectNewFace, vectNewFaceAux;
	int indP1 = vectBaseFace[0];
	const TVector3d *pP1 = &(vectVertexPoints[indP1]), *pP2;
	TVector3d vPlaneNorm;
	for(int i=0; i<numPointsInBaseFace; i++)
	{
		int i2 = i + 1;
		if(i2 == numPointsInBaseFace) i2 = 0;

		int indP2 = vectBaseFace[i2];
		pP2 = &(vectVertexPoints[indP2]);

		//look for vertex points in other base face, which can define new face
		for(int j=0; j<numPointsInOtherBaseFace; j++)
		{
			int indP = vectOtherBaseFace[j];
			const TVector3d &P = vectVertexPoints[indP];

			//check whether this point doesn't belong to the line passing through *pP1, *pP2
			if(CheckIfThreePointsAreOnOneLine(*pP1, *pP2, P, AbsTol)) continue;

			DefineNormalVia3Points(*pP1, *pP2, P, vPlaneNorm);
			vPlaneNorm.Normalize();

			vectNewFace.erase(vectNewFace.begin(), vectNewFace.end());
			int signPrevScalProd = 0;
			bool newFaceIsGood = true;
			//check whether all other points of the other base face are located on the same side of the test plane as the points of the main base 
			for(int k=0; k<totNumVertexPoints; k++)
			{
				const TVector3d &testP = vectVertexPoints[k];
				if(PracticallyEqual(testP, *pP1, AbsTol) || PracticallyEqual(testP, *pP2, AbsTol) || PracticallyEqual(testP, P, AbsTol)) continue;

				TVector3d dTestP = testP - P;
				double testScalProd = dTestP*vPlaneNorm;
				int signTestScalProd = 0;
				if(testScalProd < -AbsTol) signTestScalProd = -1;
				else if(testScalProd > AbsTol) signTestScalProd = 1;

				//if(signPrevScalProd*testScalProd < 0) { newFaceIsGood = false; break;}
				if(signPrevScalProd*signTestScalProd < 0) //OC140209
				{ 
					newFaceIsGood = false; break;
				}

				if(signTestScalProd != 0) signPrevScalProd = signTestScalProd;
				else
				{
					vectNewFace.push_back(k); //one more vertex point belongs to this face
				}
			}
			if(newFaceIsGood)
			{
				vectNewFace.push_back(indP1);
				vectNewFace.push_back(indP2);
				vectNewFace.push_back(indP);
				//break;

				//check whether the "new" face is not already present in vectIndAllFaces
				if(!CAuxParse::CheckIfVectElemArePresent(vectNewFace, vectIndAllFaces)) //OC150209
				{
					int numVertexPointsInNewFace = (int)vectNewFace.size();
					if(numVertexPointsInNewFace > 3)
					{
						//make sure that vertex indices are listed continuously, without creating self-intersecting polygon
						ReorderPointsToEnsureNonSelfIntersectingPolygon(vectVertexPoints, vectNewFace, arTol);
					}
					vectIndAllFaces.push_back(vectNewFace); //add only if not already present
					vectNewFace.erase(vectNewFace.begin(), vectNewFace.end());
					break;
				}
				else
				{
					vectNewFace.erase(vectNewFace.begin(), vectNewFace.end());
					//continue looking for points to create new face with given segment (*pP1, *pP2)
				}
			}
		}
		
		//check whether the "new" face is not already present in vectIndAllFaces
		//if(!CAuxParse::CheckIfVectElemArePresent(vectNewFace, vectIndAllFaces)) //OC150209
		//{
		//	int numVertexPointsInNewFace = (int)vectNewFace.size();
		//	if(numVertexPointsInNewFace > 3)
		//	{
		//		//make sure that vertex indices are listed continuously, without creating self-intersecting polygon
		//		ReorderPointsToEnsureNonSelfIntersectingPolygon(vectVertexPoints, vectNewFace, arTol);
		//	}
		//	vectIndAllFaces.push_back(vectNewFace); //add only if not already present
		//}
		//vectNewFace.erase(vectNewFace.begin(), vectNewFace.end());

		pP1 = pP2;
		indP1 = indP2;
	}
}

//-------------------------------------------------------------------------

void radTPolyhedron::ReorderPointsToEnsureNonSelfIntersectingPolygon(const vector<TVector3d>& vectPoints, vector<int>& vectIndPgnPoints, double* arTol)
{//may not work for non-convex polygons
	int numInd = (int)vectIndPgnPoints.size();
	if(numInd <= 0) return;

	double RelTol = arTol[0];
	double RelTolE2 = RelTol*RelTol;

	//find face center point
	TVector3d vCenPoint(0,0,0);
	for(int i=0; i<numInd; i++) vCenPoint += vectPoints[vectIndPgnPoints[i]];
	vCenPoint *= (1./numInd);

	TVector3d r0 = vectPoints[vectIndPgnPoints[0]] - vCenPoint;
	r0.Normalize();

	//calculating reference normal to polygon plane
	TVector3d vRefNorm;
	for(int j=1; j<numInd; j++)
	{
		TVector3d rj = vectPoints[vectIndPgnPoints[j]] - vCenPoint;
		rj.Normalize();
		vRefNorm = r0^rj;
		if(vRefNorm.AmpE2() > RelTolE2) break;
	}
	vRefNorm.Normalize();

	//re-ordering points
	list<radTPairIntDouble> listIndAngle;
	radTPairIntDouble firstPair(vectIndPgnPoints[0], 0.);
	listIndAngle.push_back(firstPair);

	for(int k=1; k<numInd; k++)
	{
		int curInd = vectIndPgnPoints[k];
		TVector3d rk = vectPoints[curInd] - vCenPoint;
		rk.Normalize();
		double curAngle = AngleBwUnitVectors(r0, rk, &vRefNorm);
		radTPairIntDouble curPair(curInd, curAngle);
		listIndAngle.push_back(curPair);
	}
	listIndAngle.sort(radTPairIntDouble::less);

	vectIndPgnPoints.erase(vectIndPgnPoints.begin(), vectIndPgnPoints.end());
	for(list<radTPairIntDouble>::const_iterator iter = listIndAngle.begin(); iter != listIndAngle.end(); ++iter)
	{
		vectIndPgnPoints.push_back(iter->mInt);
	}
	listIndAngle.erase(listIndAngle.begin(), listIndAngle.end());
}

//-------------------------------------------------------------------------

void radTPolyhedron::CollectAndMapUniquePolygonPoints(const radTHandlePgnAndTrans& hPgnAndTrf, vector<TVector3d>& vectPoints, vector<int>& vectInd, double AbsTol)
{
	radTPolygon *pPgn = hPgnAndTrf.PgnHndl.rep;
	radTVect2dVect &vPointsPgn = pPgn->EdgePointsVector;
	double zc = pPgn->CoordZ;
	radTrans *pTrf = hPgnAndTrf.TransHndl.rep;
	int numVertexPoints = (int)vPointsPgn.size();
	double AbsTolE2 = AbsTol*AbsTol;

	for(int i=0; i<numVertexPoints; i++)
	{
		TVector2d &vP2d = vPointsPgn[i];
		TVector3d vP(vP2d.x, vP2d.y, zc);
		vP = pTrf->TrPoint(vP);

		int curSizeVectPoints = (int)vectPoints.size();
		int indPointFound = -1;
		for(int j=0; j<curSizeVectPoints; j++)
		{
			TVector3d dP = vP - vectPoints[j];
			if(dP.AmpE2() <= AbsTolE2)
			{
				indPointFound = j; break;
			}
		}
		if(indPointFound < 0)
		{
			vectPoints.push_back(vP);
			indPointFound = curSizeVectPoints;
		}
		vectInd.push_back(indPointFound);
	}
}

//-------------------------------------------------------------------------

bool radTPolyhedron::CheckIfAllPolygonVerticesAreOnOneSideOfPlane(const radTHandlePgnAndTrans& hPgnAndTrf, const TVector3d& vPoint, const TVector3d& vNorm, double AbsTol)
{//vNorm should be unit vector
	radTPolygon *pPgn = hPgnAndTrf.PgnHndl.rep;
	radTVect2dVect &vPointsPgn = pPgn->EdgePointsVector;
	double zc = pPgn->CoordZ;
	radTrans *pTrf = hPgnAndTrf.TransHndl.rep;
	int numVertexPoints = (int)vPointsPgn.size();

	int signPrevScalProd = 0;
	for(int i=0; i<numVertexPoints; i++)
	{
		TVector2d &vP2d = vPointsPgn[i];
		TVector3d vP(vP2d.x, vP2d.y, zc);
		TVector3d vR = pTrf->TrPoint(vP) - vPoint;
		double testScalProd = vR*vNorm;
		int signTestScalProd = 0;
		if(testScalProd < -AbsTol) signTestScalProd = -1;
		else if(testScalProd > AbsTol) signTestScalProd = 1;

		//if(signPrevScalProd*testScalProd < 0) return false;
		if(signPrevScalProd*signTestScalProd < 0) return false; //OC140209
		if(signTestScalProd != 0) signPrevScalProd = signTestScalProd;
	}
	return true;
}

//-------------------------------------------------------------------------

void radTPolyhedron::SetCurrentDensityForConstCurrent(double avgCur, int indBaseFace1, int indBaseFace2)
{//assumes that the member:
 //radTVectHandlePgnAndTrans VectHandlePgnAndTrans;
 //was already set up.
	if(avgCur == 0) return;
	int numFaces = (int)VectHandlePgnAndTrans.size();
	if(numFaces <= 0) { SomethingIsWrong = 1; return;}
	if((indBaseFace1 < 0) || (indBaseFace1 >= numFaces)) { SomethingIsWrong = 1; return;}
	if((indBaseFace2 < 0) || (indBaseFace2 >= numFaces)) { SomethingIsWrong = 1; return;}

	TVector3d Ez(0,0,1);

	radTHandlePgnAndTrans hPgnAndTrans1 = VectHandlePgnAndTrans[indBaseFace1];
	radTPolygon* pPgn1 = hPgnAndTrans1.PgnHndl.rep;
	radTrans* pTrans1 = hPgnAndTrans1.TransHndl.rep;
	TVector3d vNorm1 = (pTrans1->TrBiPoint(Ez));
	TVector3d cenP1(pPgn1->CentrPoint.x, pPgn1->CentrPoint.y, pPgn1->CoordZ);
	cenP1 = pTrans1->TrPoint(cenP1);
	double S1 = pPgn1->Area();

	radTHandlePgnAndTrans hPgnAndTrans2 = VectHandlePgnAndTrans[indBaseFace2];
	radTPolygon* pPgn2 = hPgnAndTrans2.PgnHndl.rep;
	radTrans* pTrans2 = hPgnAndTrans2.TransHndl.rep;
	TVector3d vNorm2 = (pTrans2->TrBiPoint(Ez));
	TVector3d cenP2(pPgn2->CentrPoint.x, pPgn2->CentrPoint.y, pPgn2->CoordZ);
	cenP2 = pTrans2->TrPoint(cenP2);
	double S2 = pPgn2->Area();

	TVector3d vCurAxis = cenP2 - cenP1;
	double height = vCurAxis.Abs();
	if(height == 0)
	{
		SomethingIsWrong = 1;
		radTSend::ErrorMessage("Radia::Error122"); return;
	}
	double invHeight = 1./height;
	vCurAxis *= invHeight;

	double vNorm1_vCurAxis = vNorm1*vCurAxis, vNorm2_vCurAxis = vNorm2*vCurAxis;

	if((S1 <= 0) || (S2 <= 0) || (vNorm1_vCurAxis == 0) || (vNorm2_vCurAxis == 0))
	{
		SomethingIsWrong = 1;
		radTSend::ErrorMessage("Radia::Error122"); return;
	}

	if(vNorm1_vCurAxis < 0) { vNorm1_vCurAxis = -vNorm1_vCurAxis; vNorm1 *= (-1);}
	if(vNorm2_vCurAxis < 0) { vNorm2_vCurAxis = -vNorm2_vCurAxis; vNorm2 *= (-1);}

	double J1 = avgCur/(vNorm1_vCurAxis*S1);
	double J2 = avgCur/(vNorm2_vCurAxis*S2);
	double Js = (J2 - J1)*invHeight;

	J = J1*vCurAxis;

	if(Js != 0.)
	{//Matrix of linear coefficients for J is necessary
		
		J += ((-Js)*(cenP1*vCurAxis))*vCurAxis;
		
		if(pJ_LinCoef == 0) pJ_LinCoef = new TMatrix3d();
		
		double vxvy = Js*vCurAxis.x*vCurAxis.y, vxvz = Js*vCurAxis.x*vCurAxis.z, vyvz = Js*vCurAxis.y*vCurAxis.z;
		pJ_LinCoef->Str0.x = Js*vCurAxis.x*vCurAxis.x; pJ_LinCoef->Str0.y = vxvy; pJ_LinCoef->Str0.z = vxvz;
		pJ_LinCoef->Str1.x = vxvy; pJ_LinCoef->Str1.y = Js*vCurAxis.y*vCurAxis.y; pJ_LinCoef->Str1.z = vyvz;
		pJ_LinCoef->Str2.x = vxvz; pJ_LinCoef->Str2.y = vyvz; pJ_LinCoef->Str2.z = Js*vCurAxis.z*vCurAxis.z;

		J += (*pJ_LinCoef)*CentrPoint;
		J_IsNotZero = true;
		mLinTreat = 0; //treat J and lin. coef. as Relative
	}

	if(!J.isZero()) J_IsNotZero = true;
}

//-------------------------------------------------------------------------
