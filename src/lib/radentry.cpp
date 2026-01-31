
#include "rad_string_long.h"
#include "auxparse.h"
#include <Python.h>

//#ifdef WIN32
//#define WIN32_LEAN_AND_MEAN
////#include <windows.h>
//#endif

#include "radentry.h"
#include "rad_io_buffer.h"

//DEBUG
//#include <mpi.h>
//#endif

extern "C" {

void RecMag( double,double,double, double,double,double, double,double,double );
void ExtrudedPolygonDLL( double, double, double*, int, char, double* );
//void PolyhedronDLL( double*, int, int*, int*, int, double* );
void PolyhedronDLL( double*, int, int*, int*, int, double*, double*, double*, double* );
void MultGenExtrPolygonDLL( double*, int*, double*, int, double* );
void MultGenExtrRectangleDLL( double*, double*, int, double* );
//void MultGenExtrTriangleDLL( double, double, double*, double*, int, char, double*, const char*,const char*,const char* );
void MultGenExtrTriangleDLL( double, double, double*, double*, int, char, double*, const char*,const char*,const char*,const char* ); //OC30072018

void ArcMag( double,double,double, double,double, double,double, double, int, char*, double,double,double );
void ArcPolygon();
void ArcPolygonDLL( double,double, char, double*, int, double,double, int, char, double,double,double );
void CylMag( double,double,double, double, double, int, char*, double,double,double );
void RecCur( double,double,double, double,double,double, double,double,double );
void ArcCur( double,double,double, double,double, double,double, double, int, double, char*, char* );
void RaceTrack( double,double,double, double,double, double,double, double, int, double, char*, char* );
void FlmCurDLL( double*, int, double );
void ScaleCurInObj( int,double );
void BackgroundFieldSource( double,double,double );
void CoefficientFunctionFieldSource( PyObject* );
void Rectngl( double,double,double, double,double );
void Group( int*, long );
void AddToGroup( int, int*, long );
void OutGroupSize( int );
void OutGroupSubObjectKeys( int );
void DuplicateElementG3DOpt( int, const char* );
// CutElementG3DOpt REMOVED (2026-01-14) - Use Cubit/Netgen for mesh operations
void SubdivideElementG3DOpt( int, double*, char, double*, int, const char*, const char*, const char* );
void GeometricalVolume( int );
void GeometricalLimits( int );
void NumberOfDegOfFreedom( int );

void MagnOfObj( int );
void ObjField( int, char* );
void SetObjMagn( int, double,double,double );

void Translation( double,double,double );
void Rotation( double,double,double, double,double,double, double );
void PlaneSym( double,double,double, double,double,double );
void FieldInversion();
void CombineTransformLeft( int, int );
void CombineTransformRight( int, int );
void TransformObject( int, int );
void ApplySymmetry( int, int, int );

void LinearMaterial( double,double, double,double,double );
void LinearMaterial2( double,double, double );
void LinearMaterialIsotropic( double );
void LinearMaterialAnisotropic( double,double, double,double,double );
void PermanentMagnet( double,double, double,double,double );
void MagFixed( double,double,double );
void MagLinear( double,double, double,double,double );
void MagCurve( double*, int, double,double,double );

void NonlinearIsotropMaterial2( double,double, double,double, double,double );
void NonlinearIsotropMaterial3Opt( double**, long );
void NonlinearLaminatedMaterialFrm( double*,double*,double*, double, double* );
void NonlinearLaminatedMaterialTab( double*, int, double, double* );
void NonlinearAnisotropMaterialOpt0( double*, int, double*, int );
void NonlinearAnisotropMaterialOpt1( double**, double** );
void NonlinearAnisotropMaterialOpt2( double**, double );
void NonlinearAnisotropMaterialOpt3( double, double** );
void ApplyMaterial( int, int );
void MvsH( int, char*, double,double,double );

void PreRelax( int, int );
void ShowInteractMatrix(int);
int GetInteractMatrix(int, double*, int*);
void SetRelaxSubInterval(int, int, int, int);
void ShowInteractVector(int, char*);
void ManualRelax( int, int, int, double );
//void AutoRelax( int, double, int, int );
void AutoRelaxOpt( int, double, int, int, const char* );
void UpdateSourcesForRelax( int );
void SolveGen( int, double, int, int, const char* );
void SolveGenNonl( int, double, int, int, int, const char* );
int BuildMatrix( int, const char* );
#ifdef RADIA_USE_HACAPK
void SetHACApKParams( double, int, double );
void GetHACApKStats( double*, int* );
#endif
void GetSolveStats( double*, int* );
void SetBiCGSTABTolerance( double );
double GetBiCGSTABTolerance();
void SetRelaxParam( double );
double GetRelaxParam();
// SetIMASymmetry, BuildIMAMatrix REMOVED (2026-01-31) - Use BuildMatrix(obj, image) instead
void ClassifyPoints( int*, int*, int, double*, int, double );
void ComputeFieldBatch( double*, double*, int, double*, int, int );
void ComputeScalarPotentialBatch( double*, int, double*, int );
void ComputeVectorPotentialBatch( double*, int, double*, int );

void FieldArbitraryPointsArray( long, const char*, double**, long );
void Field( int, char*, double,double,double, double,double,double, int, char*, double );
void FieldEnergy( int, int, int,int,int );
void FieldForce( int, int );
void FieldForceThroughEnergy( int, int, char*, int,int,int );
void FieldTorqueThroughEnergy( int, int, char*, double,double,double, int,int,int );
void FocusingPotential( int, double,double,double, double,double,double, int );
//void FocusingKickPer( int, double,double,double, double,double,double, double,int, double,double,double, double,int,double,int, const char*, int,int,double,double );
void FocusingKickPer( int, double,double,double, double,double,double, double,double, double,double,double, double,int,double,int, const char*, int,int,double,double, const char*, double, const char* ); //OC03112019
void FocusingKickPerFormStrRep( double*,double*,double*,double*,double*, int,int, double, int, const char* );

void ParticleTrajectory( int, double, double,double,double,double, double,double, int );
void FieldInt( int, char*, char*, double,double,double, double,double,double );
void CompCriterium( double, double, double, double, double,double );
void CompPrecisionOpt( const char*, const char*, const char*, const char*, const char*, const char*, const char*, const char* );
void PhysicalUnits();
void PhysicalUnitsSet(const char*);
void RandomizationOnOrOff( char* );
void TolForConvergence( double, double, double );
void ShimSignature( int, char*, double,double,double, double,double,double, double,double,double, int, double,double,double );

int GraphicsForElemVTK( int, const char*, const char*, const char* );

void ApplyDrawAttrToElem( int, double,double,double, double );

void DeleteElement( int );
void DeleteAllElements1();
void InterruptTime( double );
void RadiaVersion();
//void DumpElem( int );
void DumpElemOpt( int*, int, const char* );
void DumpElemParseOpt( const unsigned char*, int );
void GenDump();

void ProcMPI( const char*, double*, long*, long*, long*);
//void ProcMPI( const char* );
}

//-------------------------------------------------------------------------

extern radTIOBuffer ioBuffer;
//radTIOBuffer ioBuffer; //OC, to place back!!!

//-------------------------------------------------------------------------

//#ifdef WIN32
//
//extern HINSTANCE hinstCurrentRadia;
//extern HINSTANCE hinstPreviousRadia;
//extern LPSTR lpszCmdLineRadia;
//extern int nCmdShowRadia;
//
//BOOL APIENTRY DllMain(HANDLE hModule, DWORD ul_reason_for_call, LPVOID lpReserved)
//{
//	hinstCurrentRadia = (HINSTANCE)hModule;
//    return TRUE;
//}
//
//#endif

//-------------------------------------------------------------------------

int (*pgRadYieldExternFunc)() = 0;

//-------------------------------------------------------------------------
// Note: Tetrahedral method selection was removed (2025-12-09).
// Dipole-dipole method was found numerically unstable.
// Surface charge (MSC) method is always used.
//-------------------------------------------------------------------------

int CALL RadUtiYeldFuncSet(int (*pExtFunc)())
{
	if(pExtFunc != 0) 
	{
		pgRadYieldExternFunc = pExtFunc;
	}
	return OK;
}

//-------------------------------------------------------------------------
// Copied from AlpDllEntry.cpp
const char* CALL RadErrGet(int er)
{
	return ioBuffer.GetError(er);
}

int CALL RadErrGetSize(int* siz,int er)
{
	*siz= ioBuffer.GetErrorSize(er);
	return OK;
}

int CALL RadErrGetText(char* t,int er)
{
	//strcpy(t,ioBuffer.GetError(er));
	//OC02102018 (to avoid a need to have a separate call-back function for warning in an interface):
	if(er > 0) strcpy(t,ioBuffer.GetError(er));
	else strcpy(t,ioBuffer.GetWarning(er));
	return OK;
}

//-------------------------------------------------------------------------
// Copied from AlpDllEntry.cpp
const char* CALL RadWarGet(int er)
{
	return ioBuffer.GetWarning(er);
}

int CALL RadWarGetSize(int* siz,int er)
{
	*siz= ioBuffer.GetWarningSize(er);
	return OK;
}

int CALL RadWarGetText(char* t,int er)
{
	strcpy(t,ioBuffer.GetWarning(er));
	return OK;
}

//-------------------------------------------------------------------------

int CALL RadObjRecMag(int* n, double* pP, double* pL, double* pM)
{
	RecMag(pP[0], pP[1], pP[2], pL[0], pL[1], pL[2], pM[0], pM[1], pM[2]);
	
	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjThckPgn(int* n, double xc, double lx, double* pFlatVertices, int NumVertices, char a, double* pM)
{
	ExtrudedPolygonDLL(xc, lx, pFlatVertices, NumVertices, a, pM);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjPolyhdr(int* n, double* pFlatVertices, int NumVertices, int* pFlatFaces, int* pFacesLengths, int NumFaces, double* pM, double* pM_LinCoef, double* pJ, double* pJ_LinCoef)
{
	PolyhedronDLL(pFlatVertices, NumVertices, pFlatFaces, pFacesLengths, NumFaces, pM, pM_LinCoef, pJ, pJ_LinCoef);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjMltExtPgn(int* n, double* pFlatVertices, int* pLayerLengths, double* pAttitudes, int NumLayers, double* pM)
{// pFlatVertices - flat array of 2d points
	MultGenExtrPolygonDLL(pFlatVertices, pLayerLengths, pAttitudes, NumLayers, pM);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjMltExtRtg(int* n, double* pFlatCenPts, double* pFlatRtgSizes, int NumLayers, double* pM)
{
	MultGenExtrRectangleDLL(pFlatCenPts, pFlatRtgSizes, NumLayers, pM);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjMltExtTri(int* n, double xc, double lx, double* pFlatVert, double* pFlatSubd, int nv, char a, double* pM, char* sOpt)
{
	const char *sOpt1=0, *sOpt2=0, *sOpt3=0, *sOpt4=0;
	vector<string> AuxStrings;
	if(sOpt != 0)
	{
		//char *SepStrArr[] = {";", ","};
		//CAuxParse::StringSplit(sOpt, SepStrArr, 2, " ", AuxStrings);
		//OC30072018 
		int lenStrOpt = (int)strlen(sOpt);
		char *sOptLoc = new char[lenStrOpt + 1];
		CAuxParse::StringSymbolsRemove(sOpt, (char*)" ", sOptLoc);
		CAuxParse::StringSplitNested(sOptLoc,";,", AuxStrings);
		delete[] sOptLoc;

		int AmOfTokens = (int)AuxStrings.size();
		if(AmOfTokens > 0) 
		{
			sOpt1 = (AuxStrings[0]).c_str();
			if(AmOfTokens > 1) 
			{
				sOpt2 = (AuxStrings[1]).c_str();
				if(AmOfTokens > 2) 
				{
					sOpt3 = (AuxStrings[2]).c_str();
					if(AmOfTokens > 3) sOpt4 = (AuxStrings[3]).c_str();
				}
			}
		}
	}
	
	MultGenExtrTriangleDLL(xc, lx, pFlatVert, pFlatSubd, nv, a, pM, sOpt1, sOpt2, sOpt3, sOpt4);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

//int CALL RadObjArcMag(int* n, double* P, double* R, double* Phi, double h, int nseg, char a, double* M)
//{
//	ArcMag(P[0], P[1], P[2], R[0], R[1], Phi[0], Phi[1], h, nseg, &a, M[0], M[1], M[2]);
//
//	*n = ioBuffer.OutInt();
//	return ioBuffer.OutErrorStatus();
//}

//-------------------------------------------------------------------------

int CALL RadObjArcPgnMag(int* n, double* P, char a, double* pFlatVert, int nv, double* Phi, int nseg, char sym_no, double* M)
{
	if((P == 0) || (pFlatVert == 0) || (Phi == 0) || (M == 0)) { ioBuffer.StoreErrorMessage("Radia::Error000"); *n = ioBuffer.OutInt(); return ioBuffer.OutErrorStatus();}
	//consider puting this everywhere

	ArcPolygonDLL(P[0], P[1], a, pFlatVert, nv, Phi[0], Phi[1], nseg, sym_no, M[0], M[1], M[2]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjCylMag(int* n, double* P, double r, double h, int nseg, char a, double* M)
{
	CylMag(P[0], P[1], P[2], r, h, nseg, &a, M[0], M[1], M[2]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjRecCur(int* n, double* pP, double* pL, double* pJ)
{
	RecCur(pP[0], pP[1], pP[2], pL[0], pL[1], pL[2], pJ[0], pJ[1], pJ[2]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

//int CALL RadObjArcCur(int* n, double* pP, double* pR, double* pPhi, double h, double j, int nseg, char* cManOrAuto)
int CALL RadObjArcCur(int* n, double* pP, double* pR, double* pPhi, double h, int nseg, char man_auto, char a, double j)
{
	char strManAuto[] = "man\0  ";
	man_auto = (char)toupper(man_auto);
	if(man_auto == 'A') strcpy(strManAuto, "auto\0");

	ArcCur(pP[0], pP[1], pP[2], pR[0], pR[1], pPhi[0], pPhi[1], h, nseg, j, strManAuto, &a);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

//int CALL RadObjRaceTrk(int* n, double* pP, double* pR, double* pL, double h, double j, int nseg, char* cManOrAuto)
int CALL RadObjRaceTrk(int* n, double* pP, double* pR, double* pL, double h, int nseg, char man_auto, char a, double j)
{
	char strManAuto[] = "man\0  ";
	man_auto = (char)toupper(man_auto);
	if(man_auto == 'A') strcpy(strManAuto, "auto\0");

	RaceTrack(pP[0], pP[1], pP[2], pR[0], pR[1], pL[0], pL[1], h, nseg, j, strManAuto, &a);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjFlmCur(int* n, double* pPts, int np, double i)
{
	FlmCurDLL(pPts, np, i);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjScaleCur(int n, double scaleCoef)
{
	ScaleCurInObj(n, scaleCoef);
	ioBuffer.OutInt(); // to clear buffer
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjBckg(int* n, double* pB)
{
	BackgroundFieldSource(pB[0], pB[1], pB[2]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjBckgCF(int* n, PyObject* callback)
{
	CoefficientFunctionFieldSource(callback);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjCnt(int* n, int* pKeys, int NumKeys)
{
	Group(pKeys, NumKeys);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjAddToCnt(int Cnt, int* pKeys, int NumKeys)
{
	AddToGroup(Cnt, pKeys, NumKeys);

	Cnt = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjCntSize(int* n, int Cnt)
{
	OutGroupSize(Cnt);
	//OutGroupSize(Cnt, deep);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjCntStuf(int* pCntIndexes, int Cnt)
{
	OutGroupSubObjectKeys(Cnt);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20], NumDims;
	ioBuffer.OutMultiDimArrayOfInt(pCntIndexes, Dims, NumDims);
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadObjDpl(int* n, int Obj, char* Opt1)
{
	DuplicateElementG3DOpt(Obj, Opt1);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjM(double* pM, int* arMesh, int Obj) //OC21092018
//int CALL RadObjM(int* arMesh, int Obj) //OC15092018
//int CALL RadObjM(double* pM, int Obj)
{
	arMesh[0] = 0; //OC15092018

	MagnOfObj(Obj);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	//ioBuffer.OutMultiDimArrayOfDouble(0, Dims, NumDims); //OC15092018
	ioBuffer.OutMultiDimArrayOfDouble(pM, Dims, NumDims); //OC27092018

	if(arMesh != 0)
	{
		arMesh[0] = NumDims; //OC15092018
		for(int i=0; i<NumDims; i++) arMesh[i+1] = Dims[i];
	}

	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadObjCenFld(double* pM, int* arMesh, int Obj, char type) //OC27092018
//int CALL RadObjCenFld(int* arMesh, int Obj, char type) //OC22092018
//int CALL RadObjCenFld(double* pM, int Obj, char type)
{
	ObjField(Obj, &type);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	//ioBuffer.OutMultiDimArrayOfDouble(0, Dims, NumDims); //OC22092018
	ioBuffer.OutMultiDimArrayOfDouble(pM, Dims, NumDims); //OC27092018

	if(arMesh != 0)
	{
		arMesh[0] = NumDims; //OC22092018
		for(int i=0; i<NumDims; i++) arMesh[i+1] = Dims[i];
	}

	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadObjSetM(int obj, double* pM)
{
	SetObjMagn(obj, pM[0],pM[1],pM[2]);

	ioBuffer.OutInt(); // to clear buffer
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

// RadObjCutMag REMOVED (2026-01-14) - Use Cubit/Netgen for mesh operations

//-------------------------------------------------------------------------

int CALL RadObjGeoVol(double* Vol, int Obj)
{
	GeometricalVolume(Obj);

	*Vol = ioBuffer.OutDouble();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjGeoLim(double* Lim, int Obj)
{
	GeometricalLimits(Obj);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(Lim, Dims, NumDims);
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadObjDegFre(int* Num, int Obj)
{//03102018
	NumberOfDegOfFreedom(Obj);

	*Num = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadTrfCmbL(int* n, int OrigTrf, int trf)
{
	CombineTransformLeft(OrigTrf, trf);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadTrfCmbR(int* n, int OrigTrf, int trf)
{
	CombineTransformRight(OrigTrf, trf);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadTrfInv(int* n)
{
	FieldInversion();

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

// RadTrfMlt REMOVED (2026-01-31) - Use Image symmetry instead
// The shared-DOF approach was fundamentally incompatible with MSC 6DOF hexahedra
// Use: RadSolve(..., image="+x-z") or RadBuildMatrix(obj, image="+x-z")

//-------------------------------------------------------------------------

int CALL RadTrfOrnt(int* n, int obj, int trf)
{
	TransformObject(obj, trf);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

// RadTrfPlSym REMOVED (2026-01-31) - Use Image symmetry instead
// Use: RadSolve(..., image="+x") or RadBuildMatrix(obj, image="+x")

//-------------------------------------------------------------------------

int CALL RadTrfRot(int* n, double* pP, double* pV, double phi)
{
	Rotation(pP[0], pP[1], pP[2], pV[0], pV[1], pV[2], phi);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadTrfTrsl(int* n, double* pV)
{
	Translation(pV[0], pV[1], pV[2]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatApl(int* n, int Obj, int Mat) 
{
	ApplyMaterial(Obj, Mat);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatLin(int* n, double* pKsi, double* pMr, int LenMr)
{
	if(LenMr == 3) LinearMaterial(pKsi[0], pKsi[1], pMr[0] , pMr[1], pMr[2]);
	else if(LenMr == 1) LinearMaterial2(pKsi[0], pKsi[1], pMr[0]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatLinIso(int* n, double ksi)
{
	LinearMaterialIsotropic(ksi);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatLinAniso(int* n, double* pKsi, double* pEasyAxis)
{
	LinearMaterialAnisotropic(pKsi[0], pKsi[1], pEasyAxis[0], pEasyAxis[1], pEasyAxis[2]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatPM(int* n, double Br, double Hc, double* pMagAxis)
{
	PermanentMagnet(Br, Hc, pMagAxis[0], pMagAxis[1], pMagAxis[2]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------
// Fixed magnetization permanent magnet material (Type 100)
// Magnetization does not change with H field (no demagnetization)
int CALL RadMatMagFixed(int* n, double* pMagn)
{
	MagFixed(pMagn[0], pMagn[1], pMagn[2]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------
// Linear demagnetization permanent magnet material (Type 102)
// B = Br + mu_0 * mu_rec * H, where mu_rec = Br / (mu_0 * Hc)
// Currently behaves as fixed magnetization (demagnetization not yet implemented)
int CALL RadMatMagLinear(int* n, double Br, double Hc, double* pMagAxis)
{
	MagLinear(Br, Hc, pMagAxis[0], pMagAxis[1], pMagAxis[2]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------
// User-defined demagnetization curve permanent magnet material (Type 103)
// B-H curve: [[H1, B1], [H2, B2], ...]
// Currently behaves as fixed magnetization (demagnetization not yet implemented)
int CALL RadMatMagCurve(int* n, double* pCurveData, int np, double* pMagAxis)
{
	MagCurve(pCurveData, np, pMagAxis[0], pMagAxis[1], pMagAxis[2]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatMvsH(double* pM, int* pNm, int Obj, char* id, double* pH)
{
	MvsH(Obj, id, pH[0], pH[1], pH[2]);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pM, Dims, NumDims);
	*pNm = Dims[0];
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadMatSatIsoFrm(int* n, double* pKsiMs1, double* pKsiMs2, double* pKsiMs3)
{
	//NonlinearIsotropMaterial2(pKsiMs1[0],pKsiMs1[1], pKsiMs2[0],pKsiMs2[1], pKsiMs3[0],pKsiMs3[1]);
	//OC03102018
	double KsiMs1_0=0, KsiMs1_1=0, KsiMs2_0=0, KsiMs2_1=0, KsiMs3_0=0, KsiMs3_1=0;
	if(pKsiMs1 != 0)
	{
		KsiMs1_0 = *pKsiMs1; KsiMs1_1 = *(pKsiMs1+1);
	}
	if(pKsiMs2 != 0)
	{
		KsiMs2_0 = *pKsiMs2; KsiMs2_1 = *(pKsiMs2+1);
	}
	if(pKsiMs3 != 0)
	{
		KsiMs3_0 = *pKsiMs3; KsiMs3_1 = *(pKsiMs3+1);
	}
	NonlinearIsotropMaterial2(KsiMs1_0,KsiMs1_1, KsiMs2_0,KsiMs2_1, KsiMs3_0,KsiMs3_1);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatSatIsoTab(int* n, double* pFlatMatDef, int AmOfPts)
{
	double **PointsArray = new double*[AmOfPts];
	if(PointsArray == 0) 
	{ 
		ioBuffer.StoreErrorMessage("Radia::Error900"); return ioBuffer.OutErrorStatus();
	}
	double **tPointsArray = PointsArray;
	double *tCoord = pFlatMatDef;
	for(int i=0; i<AmOfPts; i++)
	{
		//*(tPointsArray++) = tCoord;
		//tCoord += 2;

		*tPointsArray = new double[2];
		(*tPointsArray)[0] = *(tCoord++);
		(*tPointsArray)[1] = *(tCoord++);
		tPointsArray++;
	}

	NonlinearIsotropMaterial3Opt(PointsArray, (long)AmOfPts);

	*n = ioBuffer.OutInt();
	if(PointsArray != 0) 
	{
		for(int i=0; i<AmOfPts; i++) if(PointsArray[i] != 0) delete[] (PointsArray[i]);
		delete[] PointsArray;
	}
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatSatLamFrm(int* n, double* pKsiMs1, double* pKsiMs2, double* pKsiMs3, double p, double* N)
{
	NonlinearLaminatedMaterialFrm(pKsiMs1, pKsiMs2, pKsiMs3, p, N);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatSatLamTab(int* n, double* pFlatMatDef, int AmOfMatPts, double p, double* N)
{
	NonlinearLaminatedMaterialTab(pFlatMatDef, AmOfMatPts, p, N);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatSatAniso(int* n, double* pDataPar, int LenDataPar, double* pDataPer, int LenDataPer)
{
	double **ParArray = 0;
	double **PerArray = 0;

	//if(LenDataPar == 11) // {ksi1,ms1,hc1,ksi2,ms2,hc2,ksi3,ms3,hc3,ksi0,hc0}
	//{
	//}
	if(LenDataPar == 8) // {ksi1,ms1,ksi2,ms2,ksi3,ms3,ksi0,hc}
	{
		ParArray = new double*[4];
		double **tParArray = ParArray;
		double *tCoord = pDataPar;
		for(int i=0; i<4; i++)
		{
			*(tParArray++) = tCoord;
			tCoord += 2;
		}
	}

	//if(LenDataPer == 8)
	if(LenDataPer == 7)
	{
		PerArray = new double*[4];
		double **tPerArray = PerArray;
		double *tCoord = pDataPer;
		for(int i=0; i<4; i++)
		{
			*(tPerArray++) = tCoord;
			tCoord += 2;
		}
	}

	if(LenDataPar == 11) NonlinearAnisotropMaterialOpt0(pDataPar, LenDataPar, pDataPer, LenDataPer);
	//else if((LenDataPar == 8) && (LenDataPer == 8)) NonlinearAnisotropMaterialOpt1(ParArray, PerArray);
	else if((LenDataPar == 8) && (LenDataPer == 7)) NonlinearAnisotropMaterialOpt1(ParArray, PerArray);
	else if((LenDataPar == 8) && (LenDataPer == 1)) NonlinearAnisotropMaterialOpt2(ParArray, pDataPer[0]);
	//else if((LenDataPar == 1) && (LenDataPer == 8)) NonlinearAnisotropMaterialOpt3(pDataPar[0], PerArray);
	else if((LenDataPar == 1) && (LenDataPer == 7)) NonlinearAnisotropMaterialOpt3(pDataPar[0], PerArray);
	else
	{
		ioBuffer.StoreErrorMessage("Radia::Error000"); 
		return ioBuffer.OutErrorStatus();
	}

	*n = ioBuffer.OutInt();
	if(ParArray != 0) delete[] ParArray;
	if(PerArray != 0) delete[] PerArray;
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

// RadRlxPre REMOVED (2026-01-31) - Use RadBuildMatrix instead
// int CALL RadRlxPre(int* n, int Obj, int SrcObj) - REMOVED

//-------------------------------------------------------------------------

int CALL RadRlxMan(double* dOut, int* nOut, int Intrc, int Meth, int IterNum, double RlxPar)
{
	ManualRelax(Intrc, Meth, IterNum, RlxPar);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(dOut, Dims, NumDims);
	*nOut = Dims[0];
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadRlxAuto(double* dOut, int* nOut, int Intrc, double Prec, int MaxIter, int Meth, const char* Opt1)
{
	//AutoRelax(Intrc, Prec, MaxIter, Meth);
	AutoRelaxOpt(Intrc, Prec, MaxIter, Meth, Opt1); //OC240408

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(dOut, Dims, NumDims);
	*nOut = Dims[0];
	return ErrStat;
}

//-------------------------------------------------------------------------

EXP int CALL RadRlxUpdSrc(int intrc)
{
	UpdateSourcesForRelax(intrc);

	ioBuffer.OutInt(); // to clean buffer
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

// Forward declarations for conductor field computation (defined in rad_conductor_api.cpp)
bool IsConductorHandle(int handle);
bool ComputeConductorField(double* pB, int* pNb, int cond, const char* fieldId, double* pCoord);

int CALL RadFld(double* pB, int* pNb, int Obj, char* ID, double* pCoord, int Np)
{
	//=========================================================================
	// Check if object is a conductor - dispatch to conductor field computation
	//=========================================================================
	if(IsConductorHandle(Obj))
	{
		// Conductor field computation (AC analysis)
		// For single point, compute directly
		// For multiple points, loop and concatenate results
		if(Np == 1)
		{
			if(ComputeConductorField(pB, pNb, Obj, ID, pCoord))
			{
				return 0;  // Success
			}
			// Fall through to magnetic object handling if conductor computation fails
		}
		else
		{
			// Multiple points: batch computation
			int totalNb = 0;
			double* pOut = pB;

			for(int i = 0; i < Np; i++)
			{
				int nb = 0;
				double tempB[14];  // Max field components
				double* pPt = pCoord + i * 3;

				if(ComputeConductorField(tempB, &nb, Obj, ID, pPt))
				{
					for(int j = 0; j < nb; j++)
					{
						pOut[j] = tempB[j];
					}
					pOut += nb;
					if(i == 0) totalNb = nb;  // Set size from first point
				}
			}

			*pNb = totalNb * Np;
			return 0;  // Success
		}
	}

	//=========================================================================
	// Standard magnetic object field computation
	//=========================================================================
	double **PointsArray = new double*[Np];
	if(PointsArray == 0)
	{
		ioBuffer.StoreErrorMessage("Radia::Error900"); return ioBuffer.OutErrorStatus();
	}
	double **tPointsArray = PointsArray;
	double *tCoord = pCoord;
	for(int i=0; i<Np; i++)
	{
		*(tPointsArray++) = tCoord;
		tCoord += 3;
	}

	FieldArbitraryPointsArray((long)Obj, ID, PointsArray, (long)Np);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0)
	{
		if(PointsArray != 0) delete[] PointsArray;
		return ErrStat;
	}

	int Dims[20];
	int NumDims=0;
	ioBuffer.OutMultiDimArrayOfDouble(pB, Dims, NumDims);

	// Clear internal buffer to prevent memory accumulation
	ioBuffer.EraseDoubleBufferMulti();

	int TotLen = 0; //OC19012020
	if(NumDims > 0)
	{
		TotLen = 1;
		for(int k=0; k<NumDims; k++) TotLen *= Dims[k];
	}
	//int TotLen = 1;
	//for(int k=0; k<NumDims; k++) TotLen *= Dims[k];
	*pNb = TotLen;

	//DEBUG
	//std::cout << "RadFld: Nb=" << *pNb << "\n"; //DEBUG
	//std::cout.flush(); //DEBUG

	//DEBUG
	//int rank = 0;
	//MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//std::cout << "rank=" << rank << " In RadFld before exiting\n"; //DEBUG
	//std::cout.flush(); //DEBUG

	if(PointsArray != 0) delete[] PointsArray;
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldEnr(double* d, int objdst, int objsrc, int* SbdPar)
{
	int LocSbdArr[] = {0,0,0};
	if(SbdPar != 0) { LocSbdArr[0] = SbdPar[0]; LocSbdArr[1] = SbdPar[1]; LocSbdArr[2] = SbdPar[2];}

	FieldEnergy(objdst, objsrc, LocSbdArr[0], LocSbdArr[1], LocSbdArr[2]);

	*d = ioBuffer.OutDouble();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFldEnrFrc(double* pF, int* pNf, int objdst, int objsrc, char* id, int* SbdPar)
{
	int LocSbdArr[] = {0,0,0};
	if(SbdPar != 0) { LocSbdArr[0] = SbdPar[0]; LocSbdArr[1] = SbdPar[1]; LocSbdArr[2] = SbdPar[2];}

	FieldForceThroughEnergy(objdst, objsrc, id, LocSbdArr[0], LocSbdArr[1], LocSbdArr[2]);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pF, Dims, NumDims);
	*pNf = Dims[0];
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldEnrTrq(double* pF, int* pNf, int objdst, int objsrc, char* id, double* pP, int* SbdPar)
{
	int LocSbdArr[] = {0,0,0};
	if(SbdPar != 0) { LocSbdArr[0] = SbdPar[0]; LocSbdArr[1] = SbdPar[1]; LocSbdArr[2] = SbdPar[2];}

	FieldTorqueThroughEnergy(objdst, objsrc, id, pP[0], pP[1], pP[2], LocSbdArr[0], LocSbdArr[1], LocSbdArr[2]);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pF, Dims, NumDims);
	*pNf = Dims[0];
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldFrc(double* pF, int* pNf, int Obj, int Shape)
{
	FieldForce(Obj, Shape);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pF, Dims, NumDims);
	*pNf = Dims[0];
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldFocPot(double* d, int Obj, double* P1, double* P2, int np)
{
	FocusingPotential(Obj, P1[0], P1[1], P1[2], P2[0], P2[1], P2[2], np);

	*d = ioBuffer.OutDouble();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFldFocKickPer(double* pMatr1, double* pMatr2, double* pIntBtrE2, double* pArg1, double* pArg2, int* psize, int obj, double* P1, double* Ns, double per, int nper, int nps, double* Ntr, double r1, int np1, double d1, double r2, int np2, double d2, int nh, char* com, char* unit, double en, char* frm) //OC03112019
//int CALL RadFldFocKickPer(double* pMatr1, double* pMatr2, double* pIntBtrE2, double* pArg1, double* pArg2, int* psize, int obj, double* P1, double* Ns, double per, int nper, int nps, double* Ntr, double r1, int np1, double d1, double r2, int np2, double d2, int nh, char* com)
{
	FocusingKickPer(obj, P1[0], P1[1], P1[2], Ns[0], Ns[1], Ns[2], per, nper, Ntr[0], Ntr[1], Ntr[2], r1, np1, r2, np2, com, nh, nps, d1, d2, unit, en, frm); //OC03112019
	//FocusingKickPer(obj, P1[0], P1[1], P1[2], Ns[0], Ns[1], Ns[2], per, nper, Ntr[0], Ntr[1], Ntr[2], r1, np1, r2, np2, com, nh, nps, d1, d2);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	long NpKick = np1*np2;
	//long LenArr = 2*np1*np2 + np1 + np2 + 1;
	long LenArr = 3*np1*np2 + np1 + np2 + 1;

	double* pAuxBuf = new double[LenArr]; 
	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pAuxBuf, Dims, NumDims);

	double *tAuxBuf = pAuxBuf;

	double *tMatr1 = pMatr1;
	double *tMatr2 = pMatr2;
	double *tIntBtrE2 = pIntBtrE2;

	double *tArg1 = pArg1;
	double *tArg2 = pArg2;

	for(long i=0; i<NpKick; i++) *(tMatr1++) = *(tAuxBuf++);
	for(long j=0; j<NpKick; j++) *(tMatr2++) = *(tAuxBuf++);
	for(long ii=0; ii<NpKick; ii++) *(tIntBtrE2++) = *(tAuxBuf++);

	for(long k=0; k<np1; k++) *(tArg1++) = *(tAuxBuf++);
	for(long m=0; m<np2; m++) *(tArg2++) = *(tAuxBuf++);
	*psize = (int)(*tAuxBuf + 1.E-10);

	if(pAuxBuf != 0) delete[] pAuxBuf;
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldFocKickPerFormStr(char* OutStr, double* pMatr1, double* pMatr2, double* pIntBtrE2, double* pArg1, double* pArg2, int np1, int np2, double per, int nper, char* com)
{
	FocusingKickPerFormStrRep(pMatr1, pMatr2, pIntBtrE2, pArg1, pArg2, np1, np2, per, nper, com);
	
	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	ioBuffer.OutStringClean(OutStr);
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldInt(double* pF, int* pNf, int Obj, char* InfOrFin, char* id, double* P1, double* P2)
{
	FieldInt(Obj, InfOrFin, id, P1[0], P1[1], P1[2], P2[0], P2[1], P2[2]);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pF, Dims, NumDims);

	int TotLen = 1;
	for(int k=0; k<NumDims; k++) TotLen *= Dims[k];
	*pNf = TotLen;
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldLst(double* pF, int* pNf, int Obj, char* id, double* P1, double* P2, int np, char* ArgOrNoArg, double Strt)
{
	Field(Obj, id, P1[0], P1[1], P1[2], P2[0], P2[1], P2[2], np, ArgOrNoArg, Strt);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pF, Dims, NumDims);

	int TotLen = 1;
	for(int k=0; k<NumDims; k++) TotLen *= Dims[k];
	*pNf = TotLen;
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldPtcTrj(double* pF, int* pNf, int Obj, double E, double* pIC, double* pIL, int np)
{
	ParticleTrajectory(Obj, E, pIC[0], pIC[1], pIC[2], pIC[3], pIL[0], pIL[1], np);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pF, Dims, NumDims);

	int TotLen = 1;
	for(int k=0; k<NumDims; k++) TotLen *= Dims[k];
	*pNf = TotLen;
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldShimSig(double* pF, int* pNf, int obj, char* id, double* V, double* P1, double* P2, int np, double* inVi)
{
	double Vi[] = {0,0,0};
	if(inVi != 0) { Vi[0] = inVi[0]; Vi[1] = inVi[1]; Vi[2] = inVi[2];}

	ShimSignature(obj, id, V[0], V[1], V[2], P1[0], P1[1], P1[2], P2[0], P2[1], P2[2], np, Vi[0], Vi[1], Vi[2]);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pF, Dims, NumDims);

	int TotLen = 1;
	for(int k=0; k<NumDims; k++) TotLen *= Dims[k];
	*pNf = TotLen;
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldCmpCrt(int* n, double PrcB, double PrcA, double PrcBInt, double PrcFrc, double PrcTrjCrd, double PrcTrjAng)
{
	CompCriterium(PrcB, PrcA, PrcBInt, PrcFrc, PrcTrjCrd, PrcTrjAng);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFldCmpPrc(int* n, char* Opt)
{
	//const char *Opt1=0, *Opt2=0, *Opt3=0, *Opt4=0, *Opt5=0, *Opt6=0, *Opt7=0, *Opt8=0;
	const char *arOpt[8]; //OC18122019
	for(int i=0; i<8; i++) arOpt[i] = 0;
	vector<string> AuxStrings;
	if(Opt != 0)
	{
		//OC18122019
		int lenStrOpt = (int)strlen(Opt);
		char *sOptLoc = new char[lenStrOpt + 1];
		CAuxParse::StringSymbolsRemove(Opt, (char*)" ", sOptLoc);
		CAuxParse::StringSplitNested(sOptLoc, ";,", AuxStrings);
		delete[] sOptLoc;
		int AmOfOpt = (int)AuxStrings.size();
		for(int j=0; j<AmOfOpt; j++) arOpt[j] = (AuxStrings[j]).c_str();

		//char *SepStrArr[] = {(char*)";", (char*)","}; //OC04082018 (to please GCC 4.9)
		//CAuxParse::StringSplit(Opt, SepStrArr, 2, (char*)" ", AuxStrings);
		//int AmOfTokens = (int)AuxStrings.size();
		//if(AmOfTokens > 0) Opt1 = (AuxStrings[0]).c_str();
		//if(AmOfTokens > 1) Opt2 = (AuxStrings[1]).c_str();
		//if(AmOfTokens > 2) Opt3 = (AuxStrings[2]).c_str();
		//if(AmOfTokens > 3) Opt4 = (AuxStrings[3]).c_str();
		//if(AmOfTokens > 4) Opt5 = (AuxStrings[4]).c_str();
		//if(AmOfTokens > 5) Opt6 = (AuxStrings[5]).c_str();
		//if(AmOfTokens > 6) Opt7 = (AuxStrings[6]).c_str();
		//if(AmOfTokens > 7) Opt8 = (AuxStrings[7]).c_str();
	}

	CompPrecisionOpt(arOpt[0], arOpt[1], arOpt[2], arOpt[3], arOpt[4], arOpt[5], arOpt[6], arOpt[7]); //OC18122019
	//CompPrecisionOpt(Opt1, Opt2, Opt3, Opt4, Opt5, Opt6, Opt7, Opt8);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFldUnits(char* OutStr)
{
	PhysicalUnits();

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	strcpy(OutStr, ioBuffer.OutStringPtr()); //OC27092018
	//strcpy(OutStr, ioBuffer.OutString());
	ioBuffer.EraseStringBuffer();
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldUnitsSize(int* OutSize)
{
	PhysicalUnits();

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	*OutSize = (int)strlen(ioBuffer.OutStringPtr()); //27092018
	//*OutSize = (int)strlen(ioBuffer.OutString());
	ioBuffer.EraseStringBuffer();
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldUnitsSet(const char* UnitStr)
{
	PhysicalUnitsSet(UnitStr);

	int ErrStat = ioBuffer.OutErrorStatus();
	ioBuffer.EraseStringBuffer();
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldFrcShpRtg(int* n, double* pP, double* pW)
{
	Rectngl(pP[0], pP[1], pP[2], pW[0], pW[1]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFldLenRndSw(int* n, char* OnOrOff)
{
	RandomizationOnOrOff(OnOrOff);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFldLenTol(int* n, double AbsVal, double RelVal, double ZeroVal)
{
	TolForConvergence(AbsVal, RelVal, ZeroVal);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------
// RadObjDrwAtr REMOVED (2026-01-14) - Drawing attributes no longer used
// Use VTK export with ParaView for visualization
//-------------------------------------------------------------------------

int CALL RadUtiDel(int* n, int Elem)
{
	DeleteElement(Elem);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadUtiDelAll(int* n)
{
	DeleteAllElements1();

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadUtiDmp(char* OutStr, int* pSize, int* arElem, int nElem, char* AscOrBin) //OC27092018
//int CALL RadUtiDmp(char* OutStr, int* arElem, int nElem, char* AscOrBin)
//int CALL RadUtiDmp(char* OutStr, int Elem)
{
	//DumpElem(Elem);
	DumpElemOpt(arElem, nElem, AscOrBin); //OC230713

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;


	if((strcmp(AscOrBin, "asc\0") == 0) || (strcmp(AscOrBin, "Asc\0") == 0) || (strcmp(AscOrBin, "ASC\0") == 0)) 
	{
		if(pSize != 0) *pSize = (int)strlen(ioBuffer.OutStringPtr()) + 1; //to include terminating '\0' (?)
		//if(pSize != 0) *pSize = (int)strlen(ioBuffer.OutStringPtr()); //OC27092018
		if(OutStr != 0) ioBuffer.OutStringClean(OutStr); 
		//strcpy(OutStr, ioBuffer.OutString());
	}
	else 
	{
		if(pSize != 0) *pSize = (int)ioBuffer.OutByteStringSize(); //27092018
		if(OutStr != 0) ioBuffer.OutByteStringClean(OutStr); //27092018
	}
	//{
	//	long sizeData = ioBuffer.OutByteStringSize();
	//	const char *tData = ioBuffer.OutByteStringPtr(); //27092018
	//	//const char *tData = ioBuffer.OutByteString();
	//	char *tOutStr = OutStr;
	//	for(long i=0; i<sizeData; i++) *(tOutStr++) = *(tData++);
	//}

	//ioBuffer.EraseStringBuffer(); //in any case
	return ErrStat;
}

//-------------------------------------------------------------------------

//int CALL RadUtiDmpRead(char* OutStr, char* AscOrBin)
//{
//	int ErrStat = ioBuffer.OutErrorStatus();
//	if(ErrStat > 0) return ErrStat;
//
//	if((strcmp(AscOrBin, "asc\0") == 0) || (strcmp(AscOrBin, "Asc\0") == 0) || (strcmp(AscOrBin, "ASC\0") == 0))
//		strcpy(OutStr, ioBuffer.OutStringPtr()); //27092018
//		//strcpy(OutStr, ioBuffer.OutString());
//	else
//	{
//		long sizeData = ioBuffer.OutByteStringSize();
//		const char *tData = ioBuffer.OutByteStringPtr(); //27092018
//		//const char *tData = ioBuffer.OutByteString();
//		char *tOutStr = OutStr;
//		for(long i=0; i<sizeData; i++) *(tOutStr++) = *(tData++);
//	}
//
//	ioBuffer.EraseStringBuffer(); //in any case
//	return ErrStat;
//}

//-------------------------------------------------------------------------

////int CALL RadUtiDmpSize(int* OutSize, int Elem)
//int CALL RadUtiDmpSize(int* OutSize, int* arElem, int nElem, char* AscOrBin, bool doEraseBuf)
//{
//	//DumpElem(Elem);
//	DumpElemOpt(arElem, nElem, AscOrBin); //OC230713
//
//	int ErrStat = ioBuffer.OutErrorStatus();
//	if(ErrStat > 0) return ErrStat;
//
//	if((strcmp(AscOrBin, "asc\0") == 0) || (strcmp(AscOrBin, "Asc\0") == 0) || (strcmp(AscOrBin, "ASC\0") == 0))
//		*OutSize = (int)strlen(ioBuffer.OutStringPtr()); //OC27092018
//		//*OutSize = (int)strlen(ioBuffer.OutString());
//	else *OutSize = (int)ioBuffer.OutByteStringSize();
//
//	if(doEraseBuf) ioBuffer.EraseStringBuffer();
//	//leaving buffer not erased may be interesting if immediately after this RadUtiDmp will be called (then the DumpElemOpt(..) call won't need to be repeated)
//	return ErrStat;
//}

//-------------------------------------------------------------------------

int CALL RadUtiDmpPrs(int* arElem, int* pnElem, unsigned char* sBytes, int nBytes)
{//OC01102018
	DumpElemParseOpt(sBytes, nBytes);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	bool resIsList = (bool)sBytes[0];
	if(pnElem != 0)
	{
		*pnElem = 0;
		if(resIsList) 
		{
			int arDims[20], nDims=0;
			if(arElem == 0) ioBuffer.OutMultiDimArrayOfIntDims(arDims, nDims);
			else ioBuffer.OutMultiDimArrayOfInt(arElem, arDims, nDims);
			*pnElem = arDims[0]; //output can only be 1D array in this case
		}
		else 
		{
			*pnElem = 1;
			if(arElem != 0)
			{
				*arElem = ioBuffer.OutInt();
			}
		}
	}
	else if(arElem != 0)
	{
		if(resIsList) 
		{
			int arDims[20], nDims=0;
			ioBuffer.OutMultiDimArrayOfInt(arElem, arDims, nDims);
		}
		else *arElem = ioBuffer.OutInt();
	}
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadUtiIntrptTim(double* d, double t)
{
	InterruptTime(t);

	*d = ioBuffer.OutDouble();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadUtiVer(double* d)
{
	RadiaVersion();

	*d = ioBuffer.OutDouble();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadUtiMPI(int* arPar, char* sOnOff, double* arData, long* pnData, long* pRankFrom, long* pRankTo) //OC19032020
//int CALL RadUtiMPI(int* arPar, char* sOnOff)
{
	ProcMPI(sOnOff, arData, pnData, pRankFrom, pRankTo); //OC19032020
	//ProcMPI(sOnOff);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	//Should we distinguish cases what to extract?
	//if((strcmp(sOnOff, "on") == 0) || (strcmp(sOnOff, "On") == 0) || (strcmp(sOnOff, "ON") == 0) || (strcmp(sOnOff, "off") == 0) || (strcmp(sOnOff, "Off") == 0) || (strcmp(sOnOff, "OFF") == 0))
	if((strcmp(sOnOff, "on") == 0) || (strcmp(sOnOff, "On") == 0) || (strcmp(sOnOff, "ON") == 0) ||
	   (strcmp(sOnOff, "in") == 0) || (strcmp(sOnOff, "In") == 0) || (strcmp(sOnOff, "IN") == 0) ||
	   (strcmp(sOnOff, "off") == 0) || (strcmp(sOnOff, "Off") == 0) || (strcmp(sOnOff, "OFF") == 0))
	{
		int arDims[20], nDims=0;
		ioBuffer.OutMultiDimArrayOfInt(arPar, arDims, nDims);
	}//In other cases, there will be separate extraction after this

	return ErrStat;
}

//-------------------------------------------------------------------------
// "Secondary" functions
//-------------------------------------------------------------------------

// RadTrfZerPara REMOVED (2026-01-31) - Use Image symmetry instead
// RadTrfZerPerp REMOVED (2026-01-31) - Use Image symmetry instead
// These functions used RadTrfMlt internally, which has fundamental issues with MSC 6DOF hexahedra
// Use: RadSolve(..., image="+x-z") or RadBuildMatrix(obj, image="+x-z")

//-------------------------------------------------------------------------

int CALL RadSolve(double* dOut, int* nOut, int obj, double prec, int iter, int meth, const char* image)
{
	SolveGen(obj, prec, iter, meth, image);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(dOut, Dims, NumDims);
	*nOut = Dims[0];
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadSolveNonl(double* dOut, int* nOut, int obj, double prec, int iter, int meth, int nonl_method, const char* image)
{
	SolveGenNonl(obj, prec, iter, meth, nonl_method, image);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(dOut, Dims, NumDims);
	*nOut = Dims[0];
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadBuildMatrix(int* n, int ElemKey, const char* image)
{
	int result = BuildMatrix(ElemKey, image);
	*n = result;
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadUtiDataGet(char* pcData, const char typeData[3], long key) //OC04102018
//int CALL RadUtiDataGet(char* pcData, char typeData[3], long key) //OC27092018
//int CALL RadUtiDataGet(double* pData, long key) //OC15092018
{
	ioBuffer.OutDataClean(pcData, typeData, key);

	int ErrStat = ioBuffer.OutErrorStatus();
	return ErrStat;
}

// RadPreRelax REMOVED (2026-01-31) - Use RadBuildMatrix instead
// int CALL RadPreRelax(int* n, int ElemKey, int SrcElemKey) - REMOVED

//-------------------------------------------------------------------------

int CALL RadGetInteractMatrix(double* pMatrix, int* pDOF, int InteractElemKey)
{
	int result = GetInteractMatrix(InteractElemKey, pMatrix, pDOF);
	if(result == 0) return ioBuffer.OutErrorStatus();
	return 0;
}

//-------------------------------------------------------------------------

int CALL RadSetRelaxSubInterval(int InteractElemKey, int StartNo, int FinNo, int RelaxTogether)
{
	SetRelaxSubInterval(InteractElemKey, StartNo, FinNo, RelaxTogether);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

#ifdef RADIA_USE_HACAPK
int CALL RadSetHACApKParams(int* n, double eps, int leaf_size, double eta)
{
	SetHACApKParams(eps, leaf_size, eta);
	*n = 1;
	return ioBuffer.OutErrorStatus();
}

int CALL RadSetHMatrixEpsilon(int* n, double eps)
{
	// Set only epsilon, keep other params at default
	// This is ELF-compatible: magic.set_hmatrix_epsilon(eps)
	SetHACApKParams(eps, -1, -1.0);  // -1 means keep current value
	*n = 1;
	return ioBuffer.OutErrorStatus();
}

int CALL RadGetHACApKStats(double* dOut, int* nOut)
{
	GetHACApKStats(dOut, nOut);
	return ioBuffer.OutErrorStatus();
}
#endif

//-------------------------------------------------------------------------

int CALL RadGetSolveStats(double* dOut, int* nOut)
{
	GetSolveStats(dOut, nOut);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadSetBiCGSTABTol(int* n, double tol)
{
	SetBiCGSTABTolerance(tol);
	*n = 1;
	return ioBuffer.OutErrorStatus();
}

int CALL RadGetBiCGSTABTol(double* tol)
{
	*tol = GetBiCGSTABTolerance();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadSetRelaxParam(int* n, double relax)
{
	SetRelaxParam(relax);
	*n = 1;
	return ioBuffer.OutErrorStatus();
}

int CALL RadGetRelaxParam(double* relax)
{
	*relax = GetRelaxParam();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------
// Image Symmetry API - REMOVED (2026-01-31)
// Use RadBuildMatrix(obj, image) instead of RadSetIMASymmetry + RadBuildIMAMatrix
//-------------------------------------------------------------------------

// RadSetIMASymmetry REMOVED - use RadBuildMatrix(obj, "+x-z") instead
// RadBuildIMAMatrix REMOVED - use RadBuildMatrix(obj, "+x-z") instead

//-------------------------------------------------------------------------

int CALL RadClassifyPoints(int* classification, int* nearest_elem, int n_points,
                           double* points, int container_handle, double near_threshold)
{
	ClassifyPoints(classification, nearest_elem, n_points, points, container_handle, near_threshold);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFldBatch(double* B_out, double* H_out, int n_points,
                     double* points, int container_handle, int method)
{
	ComputeFieldBatch(B_out, H_out, n_points, points, container_handle, method);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFldPhi(double* phi_out, int n_points, double* points, int container_handle)
{
	ComputeScalarPotentialBatch(phi_out, n_points, points, container_handle);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFldA(double* A_out, int n_points, double* points, int container_handle)
{
	ComputeVectorPotentialBatch(A_out, n_points, points, container_handle);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFldVTS(int container_handle, const char* filename,
                   double x_min, double x_max, int nx,
                   double y_min, double y_max, int ny,
                   double z_min, double z_max, int nz,
                   int include_B, int include_H, double unit_scale)
{
	// Compute total number of points
	int n_points = nx * ny * nz;

	// Allocate arrays for field computation
	double* points = new double[n_points * 3];
	double* B_out = include_B ? new double[n_points * 3] : nullptr;
	double* H_out = include_H ? new double[n_points * 3] : nullptr;

	// Generate grid points in VTS order (i varies fastest, then j, then k)
	double dx = (nx > 1) ? (x_max - x_min) / (nx - 1) : 0.0;
	double dy = (ny > 1) ? (y_max - y_min) / (ny - 1) : 0.0;
	double dz = (nz > 1) ? (z_max - z_min) / (nz - 1) : 0.0;

	int idx = 0;
	for(int k = 0; k < nz; k++) {
		for(int j = 0; j < ny; j++) {
			for(int i = 0; i < nx; i++) {
				points[idx * 3 + 0] = x_min + i * dx;
				points[idx * 3 + 1] = y_min + j * dy;
				points[idx * 3 + 2] = z_min + k * dz;
				idx++;
			}
		}
	}

	// Compute fields using batch API (OpenMP parallelized)
	ComputeFieldBatch(B_out, H_out, n_points, points, container_handle, 0);

	// Write VTS file
	FILE* f = fopen(filename, "w");
	if(f == nullptr) {
		delete[] points;
		if(B_out) delete[] B_out;
		if(H_out) delete[] H_out;
		ioBuffer.StoreErrorMessage("Radia::Error000");
		return ioBuffer.OutErrorStatus();
	}

	// XML header
	fprintf(f, "<?xml version=\"1.0\"?>\n");
	fprintf(f, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf(f, "  <StructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n", nx-1, ny-1, nz-1);
	fprintf(f, "    <Piece Extent=\"0 %d 0 %d 0 %d\">\n", nx-1, ny-1, nz-1);

	// Points (apply unit scale for conversion to meters)
	fprintf(f, "      <Points>\n");
	fprintf(f, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	for(int i = 0; i < n_points; i++) {
		fprintf(f, "          %.15e %.15e %.15e\n",
		        points[i*3+0] * unit_scale,
		        points[i*3+1] * unit_scale,
		        points[i*3+2] * unit_scale);
	}
	fprintf(f, "        </DataArray>\n");
	fprintf(f, "      </Points>\n");

	// Point data
	fprintf(f, "      <PointData>\n");

	// B field
	if(include_B && B_out) {
		fprintf(f, "        <DataArray type=\"Float64\" Name=\"B\" NumberOfComponents=\"3\" format=\"ascii\">\n");
		for(int i = 0; i < n_points; i++) {
			fprintf(f, "          %.15e %.15e %.15e\n", B_out[i*3+0], B_out[i*3+1], B_out[i*3+2]);
		}
		fprintf(f, "        </DataArray>\n");

		// B magnitude
		fprintf(f, "        <DataArray type=\"Float64\" Name=\"B_magnitude\" format=\"ascii\">\n");
		for(int i = 0; i < n_points; i++) {
			double Bx = B_out[i*3+0], By = B_out[i*3+1], Bz = B_out[i*3+2];
			double B_mag = sqrt(Bx*Bx + By*By + Bz*Bz);
			fprintf(f, "          %.15e\n", B_mag);
		}
		fprintf(f, "        </DataArray>\n");
	}

	// H field
	if(include_H && H_out) {
		fprintf(f, "        <DataArray type=\"Float64\" Name=\"H\" NumberOfComponents=\"3\" format=\"ascii\">\n");
		for(int i = 0; i < n_points; i++) {
			fprintf(f, "          %.15e %.15e %.15e\n", H_out[i*3+0], H_out[i*3+1], H_out[i*3+2]);
		}
		fprintf(f, "        </DataArray>\n");

		// H magnitude
		fprintf(f, "        <DataArray type=\"Float64\" Name=\"H_magnitude\" format=\"ascii\">\n");
		for(int i = 0; i < n_points; i++) {
			double Hx = H_out[i*3+0], Hy = H_out[i*3+1], Hz = H_out[i*3+2];
			double H_mag = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
			fprintf(f, "          %.15e\n", H_mag);
		}
		fprintf(f, "        </DataArray>\n");
	}

	fprintf(f, "      </PointData>\n");
	fprintf(f, "    </Piece>\n");
	fprintf(f, "  </StructuredGrid>\n");
	fprintf(f, "</VTKFile>\n");

	fclose(f);

	// Cleanup
	delete[] points;
	if(B_out) delete[] B_out;
	if(H_out) delete[] H_out;

	return ioBuffer.OutErrorStatus();
}

//=========================================================================
// Conductor Analysis API Implementation (FastImp-based)
//=========================================================================

// Forward declarations for conductor functions (implemented in rad_conductor_api.cpp)
void CndRecBlock(double* P, double* L, double sigma, int num_panels);
void CndHexahedron(double* FlatVert, double sigma, int num_panels);
void CndWire(double* FlatPath, int np, char cross_section, double width, double height, double sigma, int num_panels_around);
void CndLoop(double* center, double radius, double* normal, char cross_section, double wire_width, double wire_height, double sigma, int num_panels_around, int num_panels_loop);
void CndSpiral(double* center, double inner_radius, double outer_radius, double pitch, int num_turns, double* axis, char cross_section, double wire_width, double wire_height, double sigma, int num_panels_around);
void CndCnt(int* Conds, int ncond);
void CndCntAdd(int cnt, int* Conds, int ncond);
void CndSetFormulation(int cond, const char* formulation);
void CndSetFrequency(int cond, double frequency);
void CndSetMuR(int cond, double mu_r);
void CndGetSkinDepth(double* delta, int cond);
void CndGetSurfaceImpedance(double* Z_real, double* Z_imag, int cond);
void CndSetVoltage(int cond, double V_real, double V_imag);
void CndSetCurrent(int cond, double I_real, double I_imag);
void CndGetTotalCurrent(double* I_real, double* I_imag, int cond);
void CndSetPfft(int cond, int enable);
void CndDefinePort(int cond, int* terminal1, int n1, int* terminal2, int n2);
void CndDefinePortAuto(int cond);
void CndSolve(int cond);
void CndGetImpedance(double* Z_real, double* Z_imag, int cond);
void CndImpedanceSweep(double* Z_real, double* Z_imag, int cond, double* freqs, int nfreq);
// NOTE: IsConductorHandle and ComputeConductorField declared before RadFld()
void CndGetSurfaceCurrent(double* K_real, double* K_imag, int* npanels, int cond);
void CndGetSurfaceCharge(double* sigma_real, double* sigma_imag, int* npanels, int cond);
void CndNumPanels(int* npanels, int cond);
void CndGetPanelInfo(double* centers, double* areas, int cond);
void CoupledSolve(int cond_cnt, int mag_cnt, double precision, int max_iter);
void MatSIBC(double sigma, double mu_r);
void SIBCSetType(int mat, const char* sibc_type);
void SIBCSetCrossSection(int mat, const char* shape, double* params, int nparams);
void CndFmmSetEnabled(int enabled);
void CndFmmGetEnabled(int* enabled);
void CndFmmSetParameters(int p, int ncrit, int threshold);
void CndFmmGetParameters(int* p, int* ncrit, int* threshold);

//-------------------------------------------------------------------------

int CALL RadCndRecBlock(int* n, double* P, double* L, double sigma, int num_panels)
{
	CndRecBlock(P, L, sigma, num_panels);
	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndHexahedron(int* n, double* FlatVert, double sigma, int num_panels)
{
	CndHexahedron(FlatVert, sigma, num_panels);
	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndWire(int* n, double* FlatPath, int np, char cross_section,
                    double width, double height, double sigma, int num_panels_around)
{
	CndWire(FlatPath, np, cross_section, width, height, sigma, num_panels_around);
	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndLoop(int* n, double* center, double radius, double* normal,
                    char cross_section, double wire_width, double wire_height,
                    double sigma, int num_panels_around, int num_panels_loop)
{
	CndLoop(center, radius, normal, cross_section, wire_width, wire_height,
	        sigma, num_panels_around, num_panels_loop);
	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndSpiral(int* n, double* center, double inner_radius, double outer_radius,
                      double pitch, int num_turns, double* axis,
                      char cross_section, double wire_width, double wire_height,
                      double sigma, int num_panels_around)
{
	CndSpiral(center, inner_radius, outer_radius, pitch, num_turns, axis,
	          cross_section, wire_width, wire_height, sigma, num_panels_around);
	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndCnt(int* n, int* Conds, int ncond)
{
	CndCnt(Conds, ncond);
	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndCntAdd(int cnt, int* Conds, int ncond)
{
	CndCntAdd(cnt, Conds, ncond);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndSetFormulation(int cond, const char* formulation)
{
	CndSetFormulation(cond, formulation);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndSetFrequency(int cond, double frequency)
{
	CndSetFrequency(cond, frequency);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndSetMuR(int cond, double mu_r)
{
	CndSetMuR(cond, mu_r);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndGetSkinDepth(double* delta, int cond)
{
	CndGetSkinDepth(delta, cond);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndGetSurfaceImpedance(double* Z_real, double* Z_imag, int cond)
{
	CndGetSurfaceImpedance(Z_real, Z_imag, cond);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndSetVoltage(int cond, double V_real, double V_imag)
{
	CndSetVoltage(cond, V_real, V_imag);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndSetCurrent(int cond, double I_real, double I_imag)
{
	CndSetCurrent(cond, I_real, I_imag);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndGetTotalCurrent(double* I_real, double* I_imag, int cond)
{
	CndGetTotalCurrent(I_real, I_imag, cond);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndSetPfft(int cond, int enable)
{
	CndSetPfft(cond, enable);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndDefinePort(int cond, int* terminal1, int n1, int* terminal2, int n2)
{
	CndDefinePort(cond, terminal1, n1, terminal2, n2);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndDefinePortAuto(int cond)
{
	CndDefinePortAuto(cond);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndSolve(int cond)
{
	CndSolve(cond);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndGetImpedance(double* Z_real, double* Z_imag, int cond)
{
	CndGetImpedance(Z_real, Z_imag, cond);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndImpedanceSweep(double* Z_real, double* Z_imag, int cond,
                               double* freqs, int nfreq)
{
	CndImpedanceSweep(Z_real, Z_imag, cond, freqs, nfreq);
	return ioBuffer.OutErrorStatus();
}

// NOTE: RadCndFld* functions removed - use unified RadFld() API with conductor detection

//-------------------------------------------------------------------------

int CALL RadCndGetSurfaceCurrent(double* K_real, double* K_imag, int* npanels, int cond)
{
	CndGetSurfaceCurrent(K_real, K_imag, npanels, cond);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndGetSurfaceCharge(double* sigma_real, double* sigma_imag, int* npanels, int cond)
{
	CndGetSurfaceCharge(sigma_real, sigma_imag, npanels, cond);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndNumPanels(int* npanels, int cond)
{
	CndNumPanels(npanels, cond);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndGetPanelInfo(double* centers, double* areas, int cond)
{
	CndGetPanelInfo(centers, areas, cond);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCoupledSolve(int cond_cnt, int mag_cnt, double precision, int max_iter)
{
	CoupledSolve(cond_cnt, mag_cnt, precision, max_iter);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatSIBC(int* mat, double sigma, double mu_r)
{
	MatSIBC(sigma, mu_r);
	*mat = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadSIBCSetType(int mat, const char* sibc_type)
{
	SIBCSetType(mat, sibc_type);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadSIBCSetCrossSection(int mat, const char* shape, double* params, int nparams)
{
	SIBCSetCrossSection(mat, shape, params, nparams);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndFmmSetEnabled(int enabled)
{
	CndFmmSetEnabled(enabled);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndFmmGetEnabled(int* enabled)
{
	CndFmmGetEnabled(enabled);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndFmmSetParameters(int p, int ncrit, int threshold)
{
	CndFmmSetParameters(p, ncrit, threshold);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadCndFmmGetParameters(int* p, int* ncrit, int* threshold)
{
	CndFmmGetParameters(p, ncrit, threshold);
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------
