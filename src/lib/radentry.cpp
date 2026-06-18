
#include "rad_string_long.h"
#include "auxparse.h"
#include <Python.h>

//#ifdef WIN32
//#define WIN32_LEAN_AND_MEAN
////#include <windows.h>
//#endif

#include "radentry.h"
#include "rad_io_buffer.h"
#include "rad_parallel.h"  // RegionTaskManager for auto-starting TaskManager

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
void GeometricalVolume( int );
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
void EnergyHysteresisMaterial( int, double*, double*, double*, int*, double );
void PlayHysteresisMaterial( int, double*, double*, double*, int* );
void ApplyMaterial( int, int );
void MvsH( int, char*, double,double,double );
int MatHysSaveState( int, double*, int* );
int MatHysRestoreState( int, const double*, int );
int MatHysCommitState( int );
void PreRelax( int, int );
void ShowInteractMatrix(int);
int GetInteractMatrix(int, double*, int*);
int HMatrixDensify(int, double*, int*);
double HLUTestOnHACApK(int);
int HLUDebugMaterialize(int, double*, int*, int*);
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
void SetKeepMagnetization( bool );
bool GetKeepMagnetization();
void SetNewtonMethod( bool );
bool GetNewtonMethod();
void SetNewtonDamping( bool, int, double );
void GetNewtonDampingStats( bool*, int*, double* );
void SetBInputNewton( bool );
bool GetBInputNewton();
void SetBInputHantila( bool );
bool GetBInputHantila();
void SetHantilaAlpha( double );
double GetHantilaAlpha();
void SetHantilaRelax( double );
double GetHantilaRelax();
// SetIMASymmetry, BuildIMAMatrix REMOVED (2026-01-31) - Use BuildMatrix(obj, image) instead
void ClassifyPoints( int*, int*, int, double*, int, double );
void ComputeFieldBatch( double*, double*, int, double*, int );
void ComputeScalarPotentialBatch( double*, int, double*, int );
void ComputeVectorPotentialBatch( double*, int, double*, int );

void FieldArbitraryPointsArray( long, const char*, double**, long );
void Field( int, char*, double,double,double, double,double,double, int, char*, double );
void FieldForce( int, int );
// FieldEnergy / FieldForceThroughEnergy / FieldTorqueThroughEnergy REMOVED (Phase C, 2026-04-16)
//void FocusingKickPer( int, double,double,double, double,double,double, double,int, double,double,double, double,int,double,int, const char*, int,int,double,double );
void FocusingKickPer( int, double,double,double, double,double,double, double,double, double,double,double, double,int,double,int, const char*, int,int,double,double, const char*, double, const char* ); //OC03112019
void FocusingKickPerFormStrRep( double*,double*,double*,double*,double*, int,int, double, int, const char* );

void FieldInt( int, char*, char*, double,double,double, double,double,double );
void CompCriterium( double, double, double, double, double,double );
void CompPrecisionOpt( const char*, const char*, const char*, const char*, const char*, const char*, const char*, const char* );
void PhysicalUnits();
void PhysicalUnitsSet(const char*);
void RandomizationOnOrOff( char* );
void TolForConvergence( double, double, double );
void ShimSignature( int, char*, double,double,double, double,double,double, double,double,double, int, double,double,double );

// GraphicsForElemVTK and ApplyDrawAttrToElem removed (Graphics3D code)

void DeleteElement( int );
void DeleteAllElements1();
void InterruptTime( double );
void RadiaVersion();
// DumpElem / DumpElemParseOpt / GenDump REMOVED (Phase B2a, 2026-04-15)

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
	(void)xc;
	(void)lx;
	(void)pFlatVert;
	(void)pFlatSubd;
	(void)nv;
	(void)a;
	(void)pM;
	(void)sOpt;

	if(n != 0) *n = 0;
	ioBuffer.StoreErrorMessage("Radia::Error126");
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

int CALL RadMatEnergyHysteresis(int* n, int K, double* chi,
	double* r_flat, double* f_flat, int* table_sizes, double eps)
{
	EnergyHysteresisMaterial(K, chi, r_flat, f_flat, table_sizes, eps);
	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatPlayHysteresis(int* n, int K, double* eta,
	double* r_flat, double* f_flat, int* table_sizes)
{
	PlayHysteresisMaterial(K, eta, r_flat, f_flat, table_sizes);
	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatHysSaveState(int mat, double* pState, int* pLen)
{
	return MatHysSaveState(mat, pState, pLen);
}

int CALL RadMatHysRestoreState(int mat, const double* pState, int Len)
{
	return MatHysRestoreState(mat, pState, Len);
}

int CALL RadMatHysCommitState(int mat)
{
	return MatHysCommitState(mat);
}

extern int MatHysGetNuRev(int, double*);
int CALL RadMatHysGetNuRev(int mat, double* pNuRev)
{
	return MatHysGetNuRev(mat, pNuRev);
}

extern int MatHysIrreversible(int, double*, double*);
int CALL RadMatHysIrreversible(int mat, double* pB, double* pHirr)
{
	return MatHysIrreversible(mat, pB, pHirr);
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

int CALL RadFld(double* pB, int* pNb, int Obj, char* ID, double* pCoord, int Np)
{
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

// RadFldEnr / RadFldEnrFrc / RadFldEnrTrq (energy-based API) REMOVED (Phase C, 2026-04-16)

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

// RadUtiDmp / RadUtiDmpPrs REMOVED (Phase B1, 2026-04-15) -
// .rad save/load is no longer supported.

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
#ifdef HAVE_LAPACK
	mkl_set_num_threads(1);  // BLAS calls from TaskManager threads must be single-threaded
#endif
	ngcore::RegionTaskManager rtm;  // auto-start TaskManager if not running
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
#ifdef HAVE_LAPACK
	mkl_set_num_threads(1);  // BLAS calls from TaskManager threads must be single-threaded
#endif
	ngcore::RegionTaskManager rtm;  // auto-start TaskManager if not running
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
#ifdef HAVE_LAPACK
	mkl_set_num_threads(1);  // BLAS calls from TaskManager threads must be single-threaded
#endif
	ngcore::RegionTaskManager rtm;  // auto-start TaskManager if not running
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

int CALL RadHMatrixDensify(double* pMatrix, int* pDOF, int InteractElemKey)
{
#ifdef HAVE_LAPACK
	mkl_set_num_threads(1);  // BLAS calls from TaskManager threads must be single-threaded
#endif
	ngcore::RegionTaskManager rtm;  // auto-start TaskManager: BuildHMatrix + MatVec use ngcore::ParallelFor
	int result = HMatrixDensify(InteractElemKey, pMatrix, pDOF);
	if(result == 0) return ioBuffer.OutErrorStatus();
	return 0;
}

//-------------------------------------------------------------------------

double CALL RadHLUTestOnHACApK(int InteractElemKey)
{
#ifdef HAVE_LAPACK
	mkl_set_num_threads(1);
#endif
	ngcore::RegionTaskManager rtm;
	return HLUTestOnHACApK(InteractElemKey);
}

int CALL RadHLUDebugMaterialize(int InteractElemKey, double *A_perm_out, int *lod_out, int *nd_out)
{
#ifdef HAVE_LAPACK
	mkl_set_num_threads(1);
#endif
	ngcore::RegionTaskManager rtm;
	return HLUDebugMaterialize(InteractElemKey, A_perm_out, lod_out, nd_out);
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

int CALL RadSetKeepMagnetization(int* n, int keep)
{
	SetKeepMagnetization(keep != 0);
	*n = 1;
	return ioBuffer.OutErrorStatus();
}

int CALL RadGetKeepMagnetization(int* keep)
{
	*keep = GetKeepMagnetization() ? 1 : 0;
	return ioBuffer.OutErrorStatus();
}

int CALL RadSetNewtonMethod(int* n, int use_newton)
{
	SetNewtonMethod(use_newton != 0);
	*n = 1;
	return ioBuffer.OutErrorStatus();
}

int CALL RadGetNewtonMethod(int* use_newton)
{
	*use_newton = GetNewtonMethod() ? 1 : 0;
	return ioBuffer.OutErrorStatus();
}

int CALL RadSetNewtonDamping(int* n, int enabled, int max_iter, double min_omega)
{
	SetNewtonDamping(enabled != 0, max_iter, min_omega);
	*n = 1;
	return ioBuffer.OutErrorStatus();
}

int CALL RadGetNewtonDampingStats(int* enabled, int* max_iter, double* min_omega)
{
	bool enabled_bool;
	GetNewtonDampingStats(&enabled_bool, max_iter, min_omega);
	*enabled = enabled_bool ? 1 : 0;
	return ioBuffer.OutErrorStatus();
}

int CALL RadSetBInputNewton(int* n, int enabled)
{
	SetBInputNewton(enabled != 0);
	*n = 1;
	return ioBuffer.OutErrorStatus();
}

int CALL RadGetBInputNewton(int* enabled)
{
	*enabled = GetBInputNewton() ? 1 : 0;
	return ioBuffer.OutErrorStatus();
}

int CALL RadSetBInputHantila(int* n, int enabled)
{
	SetBInputHantila(enabled != 0);
	*n = 1;
	return ioBuffer.OutErrorStatus();
}

int CALL RadGetBInputHantila(int* enabled)
{
	*enabled = GetBInputHantila() ? 1 : 0;
	return ioBuffer.OutErrorStatus();
}

int CALL RadSetHantilaAlpha(int* n, double alpha)
{
	SetHantilaAlpha(alpha);
	*n = 1;
	return ioBuffer.OutErrorStatus();
}

int CALL RadGetHantilaAlpha(double* alpha)
{
	*alpha = GetHantilaAlpha();
	return ioBuffer.OutErrorStatus();
}

int CALL RadSetHantilaRelax(int* n, double relax)
{
	SetHantilaRelax(relax);
	*n = 1;
	return ioBuffer.OutErrorStatus();
}

int CALL RadGetHantilaRelax(double* relax)
{
	*relax = GetHantilaRelax();
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
                     double* points, int container_handle)
{
	ComputeFieldBatch(B_out, H_out, n_points, points, container_handle);
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


