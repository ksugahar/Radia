/*-------------------------------------------------------------------------
*
* File name:      radinter.cpp
*
* Project:        RADIA
*
* Description:    C-style interface functions
*
* Author(s):      Oleg Chubar
*
* First release:  1997
*
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include <string.h>
#include <array>
#include <Python.h>

#include "radappl.h"
#include "gmvect.h"
#include "radiobuf.h"

//-------------------------------------------------------------------------

extern "C" {

void RecMag( double,double,double, double,double,double, double,double,double );
void ExtrudedPolygon();
void ExtrudedPolygon2();
void ExtrudedPolygonOpt( double, double, double**, int, double* );
void ExtrudedPolygonDLL( double, double, double*, int, char, double* );
void PlanarPolygon();
void Polyhedron1();
void PolyhedronOpt( double**, int, int**, int*, int, double* );
void PolyhedronDLL( double*, int, int*, int*, int, double*, double*, double*, double* );
void Polyhedron2();
void MultGenExtrPolygon();
void MultGenExtrPolygonOpt( double**, double*, int*, int, double* );
void MultGenExtrPolygonDLL( double*, int*, double*, int, double* );
void MultGenExtrPolygonCur();
void MultGenExtrPolygonMag();
void MultGenExtrRectangle();
void MultGenExtrRectangleOpt( double**, int, double* );
void MultGenExtrRectangleDLL( double*, double*, int, double* );
void MultGenExtrTriangle();
//void MultGenExtrTriangleDLL( double, double, double*, double*, int, char, double*, const char*,const char*,const char* );
void MultGenExtrTriangleDLL( double, double, double*, double*, int, char, double*, const char*,const char*,const char*,const char* ); //OC30072018

//void ArcMag( double,double,double, double,double, double,double, double, int, double,double,double, char* );
void ArcMag( double,double,double, double,double, double,double, double, int, char*, double,double,double );
void ArcPolygon();
void ArcPolygonDLL( double,double, char, double*, int, double,double, int, char, double,double,double );
//void CylMag( double,double,double, double, double, int, double,double,double, char* );
void CylMag( double,double,double, double, double, int, char*, double,double,double );
void RecCur( double,double,double, double,double,double, double,double,double );
void ArcCur( double,double,double, double,double, double,double, double, int, double, char*, char* );
void RaceTrack( double,double,double, double,double, double,double, double, int, double, char*, char* );
void FlmCur();
void FlmCurOpt( double**, long, double );
void FlmCurDLL( double*, int, double );
void Rectngl( double,double,double, double,double );
void Group( int*, long );
void AddToGroup(int, int*, long);
void OutGroupSize( int );
void OutGroupSubObjectKeys( int );

void DuplicateElementG3D();
void DuplicateElementG3DOpt( int, const char* );
void CreateFromG3DObjectWithSymmetries( int );
void SubdivideElementG3D();
void SubdivideElementG3DOpt( int, double*, char, double*, int, const char*, const char*, const char* );
void CutElementG3D();
void CutElementG3DOpt1( int, double,double,double, double,double,double, const char* );
void CutElementG3DOpt0( int, double,double,double, double,double,double );
void CutElementG3DOpt( int, double,double,double, double,double,double, const char* );
void SubdivideElementG3DByParPlanes();
void GeometricalVolume( int );
void GeometricalLimits( int );
void FldCmpMetForSubdRecMag( int, int, int );
void SetLocMgnInSbdRecMag();
void RecMagsAsExtrPolygons( char* );
void RecMagsAsPolyhedrons( char* );
void RecognizeRecMags( char* );
void ExtPgnsAsPolyhedrons( char* );
void NumberOfDegOfFreedom( int );
void MagnOfObj( int );
void ObjField( int, char* );
void SetObjMagn( int, double,double,double );
void BackgroundFieldSource( double,double,double );
void CoefficientFunctionFieldSource( PyObject* );
void ScaleCurInObj( int,double );

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
void MaterialStd( char*, double );
void NonlinearIsotropMaterial( double,double,double, double,double,double );
void NonlinearIsotropMaterial2( double,double, double,double, double,double );
void NonlinearIsotropMaterial3();
void NonlinearIsotropMaterial3Opt( double**, long );
void NonlinearLaminatedMaterialML();
void NonlinearLaminatedMaterialFrm( double*,double*,double*, double, double* );
void NonlinearLaminatedMaterialTab( double*, int, double, double* );
void NonlinearAnisotropMaterial();
void NonlinearAnisotropMaterialOpt0( double*, int, double*, int );
void NonlinearAnisotropMaterialOpt1( double**, double** );
void NonlinearAnisotropMaterialOpt2( double**, double );
void NonlinearAnisotropMaterialOpt3( double, double** );
void ApplyMaterial( int, int );
void MvsH( int, char*, double,double,double );

void Field( int, char*, double,double,double, double,double,double, int, char*, double );
//void FieldArbitraryPointsStruct();
void FieldArbitraryPointsStruct( int, char* ); //OCTEST18042016
void FieldArbitraryPointsArray( long, const char*, double**, long );
void FieldInt( int, char*, char*, double,double,double, double,double,double );
void FieldForce( int, int );
void FieldEnergy( int, int, int,int,int );
void FieldForceThroughEnergy( int, int, char*, int,int,int );
void FieldTorqueThroughEnergy( int, int, char*, double,double,double, int,int,int );
void CompCriterium( double, double, double, double, double,double );
void CompPrecision();
void CompPrecisionOpt( const char*, const char*, const char*, const char*, const char*, const char*, const char*, const char* );
void MultipoleThresholds(double, double, double, double); // Maybe to be removed later
void PreRelax( int, int );
void ShowInteractMatrix(int);
void ShowInteractVector(int, char*);
void ManualRelax( int, int, int, double );
//void AutoRelax( int, double, int, int );
void AutoRelax();
void AutoRelaxOpt( int, double, int, int, const char* );
void UpdateSourcesForRelax( int );
void SolveGen( int, double, int, int );
void ParticleTrajectory( int, double, double,double,double,double, double,double, int );
void FocusingPotential( int, double,double,double, double,double,double, int );
//void FocusingKickPer( int, double,double,double, double,double,double, double,int, double,double,double, double,int,double,int, const char*, int,int,double,double, const char*, double );
void FocusingKickPer( int, double,double,double, double,double,double, double,double, double,double,double, double,int,double,int, const char*, int,int,double,double, const char*, double, const char* );
void FocusingKickPerFormStrRep( double*,double*,double*,double*,double*, int,int, double, int, const char* );
//void FocusingKick( int, double,double,double, double,double,double, double*,long,int, double,double,double, double,int,double,int, const char*, double,double );
void FocusingKickML();
void ShimSignature( int, char*, double,double,double, double,double,double, double,double,double, int, double,double,double );

void TolForConvergence( double, double, double );
void RandomizationOnOrOff( char* );
void PhysicalUnits();

void GraphicsForElemWithoutSymChilds(int);
//void GraphicsForElemWithSymChilds(int);
void GraphicsForElemWithSymChildsExt();
void GraphicsForAllWithoutSymChilds();
void GraphicsForAllWithSymChilds();
int GraphicsForElemVTK( int, const char*, const char*, const char* );

void ApplyDrawAttrToElem( int, double,double,double, double );
void ApplyColorToElem( int, double,double,double );
void RemoveDrawAttrFromElem(int);

void DeleteElement( int );
void DeleteAllElements1();
void InterruptTime( double );
void RadiaVersion();
//void DumpElem( int );
void DumpElem();
void DumpElemOpt( int*, int, const char* );
void DumpElemParse();
void DumpElemParseOpt( const unsigned char*, int );
void GenDump();
void DeleteAllElements2();

void StartProf(int, int, int);
void StopProf();
void OutCommandsInfo();
void ReturnInput(double, int);
void MemAllocMethForIntrctMatr(char*);

void ProcMPI( const char*, double*, long*, long*, long*);
//void ProcMPI( const char* );

//int AuxSetOptionNameAndValue(const char* OptTot, char* OptName, char** OptValue);
int AuxSetOptionNameAndValue(const char* OptTot, char* OptName, const char** OptValue);
}


//-------------------------------------------------------------------------

radTApplication rad;
radTIOBuffer ioBuffer;
radTYield radYield;
radTConvergRepair& radCR = rad.CnRep;

//-------------------------------------------------------------------------

//This is a "missing" function in VC++8 SDK;
//To allow linking DEBUG configuration with Libs compiled in release mode
//copied from Microsoft's "invarg.c"
//OC101015: Commented-out; no longer required for VC++14 (?)
//#ifdef WIN32
//#ifdef _DEBUG
//extern "C" _CRTIMP void __cdecl _invalid_parameter_noinfo(void)
//{
//	_invalid_parameter(nullptr, nullptr, nullptr, 0, 0);
//}
//#endif
//#endif

//-------------------------------------------------------------------------

int AuxSetOptionNameAndValue(const char* OptTot, char* OptName, const char** OptValue)
{
	if(*OptTot != '\0')
	{
		strncpy(OptName, OptTot, 199);
		OptName[199] = '\0';
		//char *pEndOptName = strrchr(OptName, '-');
		char *pEndOptName = strrchr(OptName, '>'); //OC19122019
		//if(pEndOptName == nullptr) { rad.Send.ErrorMessage("Radia::Error062"); return 0;}
		if((pEndOptName == nullptr) || (*(--pEndOptName) != '-')) { rad.Send.ErrorMessage("Radia::Error062"); return 0;} //OC19122019
		*pEndOptName = '\0';
		//*OptValue = strrchr(OptTot, '>') + 1;
		*OptValue = pEndOptName + 2; //OC19122019
	}
	return 1;
}

//-------------------------------------------------------------------------

//void AuxParseOptionNamesAndValues(int Nopt, const char** NonParsedOpts, const char** OptionNames, const char** OptionValues, int& OptionCount)
void AuxParseOptionNamesAndValues(const char** NonParsedOpts, const char** OptionNames, const char** OptionValues, int& OptionCount)
{
	int Nopt = OptionCount;
	OptionCount = 0;
	for(int i=0; i<Nopt; i++)
	{
		const char* Opt = NonParsedOpts[i];
		if(Opt != 0)
		{
			if(*Opt != '\0') 
			{
				//if(!AuxSetOptionNameAndValue(Opt, (char*)(OptionNames[OptionCount]), (char**)(OptionValues + OptionCount))) return;
				if(!AuxSetOptionNameAndValue(Opt, (char*)(OptionNames[OptionCount]), OptionValues + OptionCount)) return;
				OptionCount++;
			}
		}
	}
}

//-------------------------------------------------------------------------



//-------------------------------------------------------------------------

void RecMag(double xc, double yc, double zc, 
			double Lx, double Ly, double Lz, 
			double Mx, double My, double Mz)
{
	double J[] = {0.,0.,0.};
	double CPoi[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	//double CPoi[] = {xc, yc, zc};

	short OldActOnDoubles = radCR.ActOnDoubles;
	if(rad.TreatRecMagsAsExtrPolygons) radCR.ActOnDoubles = 0;

	double Dims[] = {radCR.Double(Lx), radCR.Double(Ly), radCR.Double(Lz)};
	//double Dims[] = {Lx, Ly, Lz};
	double Magn[] = {Mx, My, Mz};
	rad.SetRecMag(CPoi, 3, Dims, 3, Magn, 3, J, 3, 0);

	if(rad.TreatRecMagsAsExtrPolygons) radCR.ActOnDoubles = OldActOnDoubles;
}

//-------------------------------------------------------------------------

void RecCur(double xc, double yc, double zc, 
			double Lx, double Ly, double Lz, 
			double Jx, double Jy, double Jz)
{
	double Magn[] = {0.,0.,0.};
	double CPoi[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	double Dims[] = {radCR.Double(Lx), radCR.Double(Ly), radCR.Double(Lz)};
	double J[] = {Jx, Jy, Jz};
	rad.SetRecMag(CPoi, 3, Dims, 3, Magn, 3, J, 3, 1);
}

//-------------------------------------------------------------------------

void SetExtrPolygFirstPoint(double xc, double Lx, TVector2d& FirstPoint2d, char a, double* FirstPoi) //OC040306
{
	if((a == 'x') || (a == 'X'))
	{
		FirstPoi[0] = xc - 0.5*Lx;
		FirstPoi[1] = FirstPoint2d.x;
		FirstPoi[2] = FirstPoint2d.y;
	}
	else if((a == 'y') || (a == 'Y'))
	{
		FirstPoi[0] = FirstPoint2d.y;
		FirstPoi[1] = xc - 0.5*Lx;
		FirstPoi[2] = FirstPoint2d.x;
	}
	else
	{
		FirstPoi[0] = FirstPoint2d.x;
		FirstPoi[1] = FirstPoint2d.y;
		FirstPoi[2] = xc - 0.5*Lx;
	}
}

//-------------------------------------------------------------------------

void ExtrudedPolygon()
{

}

//-------------------------------------------------------------------------

void ExtrudedPolygon2()
{

}

//-------------------------------------------------------------------------

void MultGenExtrTriangle()
{

}

//-------------------------------------------------------------------------

void ExtrudedPolygonOpt(double xc, double Lx, double** Polygon, int AmOfVertices, double* M)
{
	std::vector<TVector2d> vArrayOfPoints2d(AmOfVertices);
	TVector2d* ArrayOfPoints2d = vArrayOfPoints2d.data();
	TVector2d* tArrayOfPoints2d = ArrayOfPoints2d;
	double** tPolygon = Polygon;

	for(int i=0; i<AmOfVertices; i++)
	{
		double* tPoint = *(tPolygon++);
		tArrayOfPoints2d->x = *(tPoint++);
		(tArrayOfPoints2d++)->y = *tPoint;
	}
	double FirstPoi[] = {xc - 0.5*Lx, ArrayOfPoints2d->x, ArrayOfPoints2d->y};

	rad.SetExtrudedPolygon(FirstPoi, 3, radCR.Double(Lx), ArrayOfPoints2d, AmOfVertices, M, 3, "x");
	// RAII: automatic cleanup
}

//-------------------------------------------------------------------------

void ExtrudedPolygonDLL(double xc, double Lx, double* Polygon, int AmOfVertices, char a, double* M)
{
	std::vector<TVector2d> vArrayOfPoints2d(AmOfVertices);
	TVector2d* ArrayOfPoints2d = vArrayOfPoints2d.data();
	TVector2d* tArrayOfPoints2d = ArrayOfPoints2d;
	double* tPolygon = Polygon;

	for(int i=0; i<AmOfVertices; i++)
	{
		tArrayOfPoints2d->x = *(tPolygon++);
		(tArrayOfPoints2d++)->y = *(tPolygon++);
	}

	//double FirstPoi[] = {xc - 0.5*Lx, ArrayOfPoints2d->x, ArrayOfPoints2d->y};
	double FirstPoi[3];
	SetExtrPolygFirstPoint(xc, Lx, ArrayOfPoints2d[0], a, FirstPoi); //OC040306

	for(int i1=0; i1<3; i1++) FirstPoi[i1] = radCR.Double(FirstPoi[i1]);

	rad.SetExtrudedPolygon(FirstPoi, 3, radCR.Double(Lx), ArrayOfPoints2d, AmOfVertices, M, 3, &a);
	// RAII: automatic cleanup
}

//-------------------------------------------------------------------------

void MultGenExtrTriangleDLL(double xc, double lx, double* pFlatVert, double* pFlatSubd, int nv, char a, double* pM, const char* sOpt1, const char* sOpt2, const char* sOpt3, const char* sOpt4) //OC30072018
//void MultGenExtrTriangleDLL(double xc, double lx, double* pFlatVert, double* pFlatSubd, int nv, char a, double* pM, const char* sOpt1, const char* sOpt2, const char* sOpt3)
{
	std::vector<TVector2d> vArrayOfPoints2d(nv);
	TVector2d* ArrayOfPoints2d = vArrayOfPoints2d.data();

	TVector2d* tArrayOfPoints2d = ArrayOfPoints2d;
	double* tPolygon = pFlatVert;

	for(int i=0; i<nv; i++)
	{
		tArrayOfPoints2d->x = *(tPolygon++);
		(tArrayOfPoints2d++)->y = *(tPolygon++);
	}

	double FirstPoi[3];
	SetExtrPolygFirstPoint(xc, lx, ArrayOfPoints2d[0], a, FirstPoi); //OC040306
	for(int i1=0; i1<3; i1++) FirstPoi[i1] = radCR.Double(FirstPoi[i1]);

	//std::array<char, 200> CharBuf1, CharBuf2, CharBuf3;
	std::array<char, 200> CharBuf1, CharBuf2, CharBuf3, CharBuf4; //OC30072018
	//const char* OptionNames[] = {CharBuf1.data(), CharBuf2.data(), CharBuf3.data()};
	const char* OptionNames[] = {CharBuf1.data(), CharBuf2.data(), CharBuf3.data(), CharBuf4.data()}; //OC30072018
	//const char* OptionValues[] = {0,0,0};
	const char* OptionValues[] = {0,0,0,0}; //OC30072018
	//const char* NonParsedOpts[] = {sOpt1, sOpt2, sOpt3};
	const char* NonParsedOpts[] = {sOpt1, sOpt2, sOpt3, sOpt4}; //OC30072018
	//int OptionCount = 3;
	int OptionCount = 4; //OC30072018
	//AuxParseOptionNamesAndValues(3, NonParsedOpts, OptionNames, OptionValues, OptionCount);
	AuxParseOptionNamesAndValues(NonParsedOpts, OptionNames, OptionValues, OptionCount);
	
	rad.SetMultGenExtrTriangle(FirstPoi, 3, radCR.Double(lx), ArrayOfPoints2d, nv, pFlatSubd, pM, 3, &a, OptionNames, OptionValues, OptionCount);
	// RAII: automatic cleanup
}

//-------------------------------------------------------------------------

void PlanarPolygon()
{

}

//-------------------------------------------------------------------------

void Polyhedron1()
{

}

//-------------------------------------------------------------------------

void PolyhedronOpt(double** Vertices, int AmOfVertices, int** Faces, int* AmOfPoInFaces, int AmOfFaces, double* M)
{
	std::vector<TVector3d> vArrayOfPoints(AmOfVertices);
	TVector3d* ArrayOfPoints = vArrayOfPoints.data();
	TVector3d* tArrayOfPoints = ArrayOfPoints;
	double** tVertices = Vertices;
	for(int i=0; i<AmOfVertices; i++)
	{
		double* tPoint = *(tVertices++);
		tArrayOfPoints->x = *(tPoint++);
		tArrayOfPoints->y = *(tPoint++);
		(tArrayOfPoints++)->z = *(tPoint++);
	}

	rad.SetPolyhedron1(ArrayOfPoints, AmOfVertices, Faces, AmOfPoInFaces, AmOfFaces, M, 0, 0, 0);
	//rad.SetPolyhedron1(ArrayOfPoints, AmOfVertices, Faces, AmOfPoInFaces, AmOfFaces, M, 3);
	// RAII: automatic cleanup
}

//-------------------------------------------------------------------------

void PolyhedronDLL(double* Vertices, int AmOfVertices, int* InFaces, int* AmOfPoInFaces, int AmOfFaces, double* M, double* M_LinCoef, double* J, double* J_LinCoef)
{
	bool MisDefined = false, M_LinCoefDefined = false; //declare at te top because of use of goto
	bool JisDefined = false, J_LinCoefDefined = false;
	bool LocMisDefined = false, LocM_LinCoefDefined = false;
	bool LocJisDefined = false, LocJ_LinCoefDefined = false;

	double *tVertices=0;
	TVector3d *tArrayOfPoints=0;
	int **Faces=0, **tFaces=0, *tInFaces=0;
	double *pJ=0, *pJ_LinCoef=0;

	std::vector<TVector3d> vArrayOfPoints(AmOfVertices);
	TVector3d* ArrayOfPoints = vArrayOfPoints.data();
	//TVector3d* tArrayOfPoints = ArrayOfPoints;
	tArrayOfPoints = ArrayOfPoints;
	//double* tVertices = Vertices;
	tVertices = Vertices;
	for(int i=0; i<AmOfVertices; i++)
	{
		tArrayOfPoints->x = *(tVertices++);
		tArrayOfPoints->y = *(tVertices++);
		(tArrayOfPoints++)->z = *(tVertices++);
	}

	//int** Faces = new int*[AmOfFaces];
	std::vector<int*> vFaces(AmOfFaces);
	Faces = vFaces.data();
	//int** tFaces = Faces;
	tFaces = Faces;
	//int* tInFaces = InFaces;
	tInFaces = InFaces;

	for(int k=0; k<AmOfFaces; k++)
	{
		*(tFaces++) = tInFaces;
		tInFaces += AmOfPoInFaces[k];
	}

	if(M != 0)
	{
		//bool LocMisDefined = false;
		LocMisDefined = false;
		for(int kk=0; kk<3; kk++) { if(M[kk] != 0.) { LocMisDefined = true; break;}}
		MisDefined = LocMisDefined;
	}
	if(M_LinCoef != 0)
	{
		//bool LocM_LinCoefDefined = false;
		LocM_LinCoefDefined = false;
		for(int kk=0; kk<9; kk++) { if(M_LinCoef[kk] != 0.) { LocM_LinCoefDefined = true; break;}}
		M_LinCoefDefined = LocM_LinCoefDefined;
	}

	if(J != 0)
	{
		//bool LocJisDefined = false;
		LocJisDefined = false;
		for(int kk=0; kk<3; kk++) { if(J[kk] != 0.) { LocJisDefined = true; break;}}
		JisDefined = LocJisDefined;
	}
	if(J_LinCoef != 0)
	{
		//bool LocJ_LinCoefDefined = false;
		LocJ_LinCoefDefined = false;
		for(int kk=0; kk<9; kk++) { if(J_LinCoef[kk] != 0.) { LocJ_LinCoefDefined = true; break;}}
		J_LinCoefDefined = LocJ_LinCoefDefined;
	}
	
	if(MisDefined && JisDefined) { rad.Send.ErrorMessage("Radia::Error120"); goto Finish;}
	if(M_LinCoefDefined) { rad.Send.ErrorMessage("Radia::Error121"); goto Finish;}

	//double *pJ = JisDefined? J : 0;
	//double *pJ_LinCoef = J_LinCoefDefined? J_LinCoef : 0;
	if(JisDefined) pJ = J;
	if(J_LinCoefDefined) pJ_LinCoef = J_LinCoef;

	rad.SetPolyhedron1(ArrayOfPoints, AmOfVertices, Faces, AmOfPoInFaces, AmOfFaces, M, 0, pJ, pJ_LinCoef);

Finish:
	// RAII: automatic cleanup via vArrayOfPoints and vFaces
	return;
}

//-------------------------------------------------------------------------

void Polyhedron2()
{

}

//-------------------------------------------------------------------------

void RecMagsAsExtrPolygons(char* OnOrOff)
{
	rad.RecMagsAsExtrPolygons(OnOrOff);
}

//-------------------------------------------------------------------------

void RecMagsAsPolyhedrons(char* OnOrOff)
{
	rad.RecMagsAsPolyhedrons(OnOrOff);
}

//-------------------------------------------------------------------------

void RecognizeRecMags(char* OnOrOff)
{
	rad.RecognizeRecMags(OnOrOff);
}

//-------------------------------------------------------------------------

void ExtPgnsAsPolyhedrons(char* OnOrOff)
{
	rad.ExtPgnsAsPolyhedrons(OnOrOff);
}

//-------------------------------------------------------------------------

void MultGenExtrPolygon()
{

}

//-------------------------------------------------------------------------

void MultGenExtrPolygonOpt(double** Layers, double* Heights, int* AmOfPointsInLayers, int AmOfLayers, double* M)
{
	// RAII: Use std::vector for all allocations
	std::vector<TVector2d*> vLayerPolygons(AmOfLayers);
	std::vector<std::vector<TVector2d>> vLayerStorage(AmOfLayers);
	TVector2d** LayerPolygons = vLayerPolygons.data();
	TVector2d** tLayerPolygons = LayerPolygons;
	double** tLayers = Layers;

	for(int i=0; i<AmOfLayers; i++)
	{
		int AmOfPointsInTheLayer = AmOfPointsInLayers[i];
		vLayerStorage[i].resize(AmOfPointsInTheLayer);
		*tLayerPolygons = vLayerStorage[i].data();

		TVector2d* tPoints = *(tLayerPolygons++);
		double* tLayer = *(tLayers++);
		for(int k=0; k<AmOfPointsInTheLayer; k++)
		{
			tPoints->x = *(tLayer++); (tPoints++)->y = *(tLayer++);
		}
	}
	rad.SetMultGenExtrPolygon(LayerPolygons, AmOfPointsInLayers, Heights, AmOfLayers, M, 3);

	// RAII: vLayerPolygons and vLayerStorage cleaned up automatically
}

//-------------------------------------------------------------------------

void MultGenExtrPolygonDLL(double* Layers, int* AmOfPointsInLayers, double* Heights, int AmOfLayers, double* M)
{
	// RAII: Use std::vector for all allocations
	std::vector<TVector2d*> vLayerPolygons(AmOfLayers);
	std::vector<std::vector<TVector2d>> vLayerStorage(AmOfLayers);
	TVector2d** LayerPolygons = vLayerPolygons.data();
	TVector2d** tLayerPolygons = LayerPolygons;

	double* tLayer = Layers;
	for(int i=0; i<AmOfLayers; i++)
	{
		int AmOfPointsInTheLayer = AmOfPointsInLayers[i];
		vLayerStorage[i].resize(AmOfPointsInTheLayer);
		*tLayerPolygons = vLayerStorage[i].data();

		TVector2d* tPoints = *(tLayerPolygons++);
		for(int k=0; k<AmOfPointsInTheLayer; k++)
		{
			tPoints->x = *(tLayer++); (tPoints++)->y = *(tLayer++);
		}
		Heights[i] = radCR.Double(Heights[i]);
	}

	rad.SetMultGenExtrPolygon(LayerPolygons, AmOfPointsInLayers, Heights, AmOfLayers, M, 3);

	// RAII: vLayerPolygons and vLayerStorage cleaned up automatically
}

//-------------------------------------------------------------------------



//-------------------------------------------------------------------------



//-------------------------------------------------------------------------

void MultGenExtrPolygonCurOrMag(char cur_or_mag)
//void MultGenExtrPolygonCur()
{//radObjMltExtPgnCur[z:0,a:"z",({{x1,y1},{x2,y2},...},{{R1,T1,H1},{R2,T2,H2},...}},I,Frame->Loc|Lab]
//or radObjMltExtPgnMag[z:0,a:"z",({{x1,y1},{x2,y2},...},{{R1,T1,H1},{R2,T2,H2},...}},{{mx1,my1,mz1},{mx2,my2,mz2},...},Frame->Loc|Lab]
//Ri: {{x,y,z},{vx,vy,vz},phi}
//Ti: {vx,vy,vz}
//Hi: k or {kx,ky} or {{kx,ky},phi}

}

//-------------------------------------------------------------------------

void MultGenExtrPolygonCur()
{
	MultGenExtrPolygonCurOrMag('c');
}
void MultGenExtrPolygonMag()
{
	MultGenExtrPolygonCurOrMag('m');
}

//-------------------------------------------------------------------------

void MultGenExtrRectangle()
{

}

//-------------------------------------------------------------------------

void MultGenExtrRectangleOpt(double** Layers, int AmOfLayers, double* M)
{
	std::vector<TVector3d> vRectCenPoints(AmOfLayers);
	std::vector<TVector2d> vRectDims(AmOfLayers);
	TVector3d* RectCenPoints = vRectCenPoints.data();
	TVector2d* RectDims = vRectDims.data();
	TVector3d* tRectCenPoints = RectCenPoints;
	TVector2d* tRectDims = RectDims;
	double** tLayers = Layers;

	for(int i=0; i<AmOfLayers; i++)
	{
		double* tCoord = *(tLayers++);
		tRectCenPoints->x = *(tCoord++); tRectCenPoints->y = *(tCoord++); (tRectCenPoints++)->z = *(tCoord++); 
		tRectDims->x = *(tCoord++); (tRectDims++)->y = *(tCoord++);
	}
	rad.SetMultGenExtrRectangle(RectCenPoints, RectDims, AmOfLayers, M, 3);
	// RAII: automatic cleanup
}

//-------------------------------------------------------------------------

void MultGenExtrRectangleDLL(double* pFlatCenPts, double* pFlatRtgSizes, int AmOfLayers, double* pM)
{
	std::vector<TVector3d> vRectCenPoints(AmOfLayers);
	std::vector<TVector2d> vRectDims(AmOfLayers);
	TVector3d* RectCenPoints = vRectCenPoints.data();
	TVector2d* RectDims = vRectDims.data();

	TVector3d* tRectCenPoints = RectCenPoints;
	TVector2d* tRectDims = RectDims;

	double* tCenPts = pFlatCenPts;
	double* tRtgSizes = pFlatRtgSizes;
	for(int i=0; i<AmOfLayers; i++)
	{
		tRectCenPoints->x = *(tCenPts++); tRectCenPoints->y = *(tCenPts++); (tRectCenPoints++)->z = *(tCenPts++); 
		tRectDims->x = *(tRtgSizes++); (tRectDims++)->y = *(tRtgSizes++);
	}
	rad.SetMultGenExtrRectangle(RectCenPoints, RectDims, AmOfLayers, pM, 3);
	// RAII: automatic cleanup
}

//-------------------------------------------------------------------------

void ArcMag(double xc, double yc, double zc, double rmin, double rmax, double phimin, double phimax, double Height, int NumberOfSegm, char* Orient, double mx, double my, double mz)
{
	double CPoi[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	double Radii[] = {fabs(radCR.Double(rmin)), fabs(radCR.Double(rmax))};
	double Angles[] = {phimin, phimax}; // Consider shrinking this
	double Magn[] = {mx, my, mz};
	rad.SetArcMag(CPoi, 3, Radii, 2, Angles, 2, Height, NumberOfSegm, Magn, 3, Orient);
}

//-------------------------------------------------------------------------

void ArcPolygon()
{

}

//-------------------------------------------------------------------------

void ArcPolygonDLL(double xc, double yc, char a, double* pFlatVert, int nv, double PhiMin, double PhiMax, int nseg, char sym_no, double mx, double my, double mz)
{
	double CenP[] = {radCR.Double(xc), radCR.Double(yc)};

	std::vector<TVector2d> vArrayOfPoints2d(nv);
	TVector2d* ArrayOfPoints2d = vArrayOfPoints2d.data();
	double *tFlatVert = pFlatVert;
	TVector2d *tArrayOfPoints2d = ArrayOfPoints2d;
	for(int i=0; i<nv; i++)
	{
		tArrayOfPoints2d->x = *(tFlatVert++);
		tArrayOfPoints2d->y = *(tFlatVert++);
		tArrayOfPoints2d++;
	}

	double Angles[] = {radCR.DoublePlus(PhiMin), radCR.DoubleMinus(PhiMax)};

	char SymOrNoSymStr[20];
	if((sym_no == 's') || (sym_no == 'S')) {
		strncpy(SymOrNoSymStr, "sym", 19);
		SymOrNoSymStr[19] = '\0';
	}
	else {
		strncpy(SymOrNoSymStr, "nosym", 19);
		SymOrNoSymStr[19] = '\0';
	}

	double Magn[] = {mx, my, mz};

	rad.SetArcPolygon(CenP, &a, ArrayOfPoints2d, nv, Angles, nseg, SymOrNoSymStr, Magn);
	// RAII: automatic cleanup
}

//-------------------------------------------------------------------------

//void CylMag(double xc, double yc, double zc, double r, double h, int NumberOfSegm, double mx, double my, double mz, char* Orient)
//{
//	double CPoi[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
//	double Magn[] = {mx, my, mz};
//	rad.SetCylMag(CPoi, 3, r, h, NumberOfSegm, Magn, 3, Orient);
//}

//-------------------------------------------------------------------------

void CylMag(double xc, double yc, double zc, double r, double h, int NumberOfSegm, char* Orient, double mx, double my, double mz)
{
	double CPoi[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	double Magn[] = {mx, my, mz};
	rad.SetCylMag(CPoi, 3, r, h, NumberOfSegm, Magn, 3, Orient);
}

//-------------------------------------------------------------------------

void ArcCur(double xc, double yc, double zc, double rmin, double rmax, double phimin, double phimax, double Height, int NumberOfSegm, double J_azim, char* ManOrAuto, char* Orient)
{
	double CPoi[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	double Radii[] = {fabs(radCR.Double(rmin)), fabs(radCR.Double(rmax))};
	double Angles[] = {phimin, phimax}; // Consider shrinking this
	rad.SetArcCur(CPoi, 3, Radii, 2, Angles, 2, Height, J_azim, NumberOfSegm, ManOrAuto, Orient);
}

//-------------------------------------------------------------------------

void RaceTrack(double xc, double yc, double zc, double rmin, double rmax, double Lx, double Ly, double Height, int NumberOfSegm, double J_azim, char* ManOrAuto, char* Orient)
{
	double CPoi[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	double Radii[] = {fabs(radCR.Double(rmin)), fabs(radCR.Double(rmax))};
	double StrPartDims[] = {(Lx==0.)? Lx : radCR.Double(Lx), (Ly==0.)? Ly : radCR.Double(Ly)};
	rad.SetRaceTrack(CPoi, 3, Radii, 2, StrPartDims, 2, Height, J_azim, NumberOfSegm, ManOrAuto, Orient);
}

//-------------------------------------------------------------------------

void FlmCur()
{

}

//-------------------------------------------------------------------------

void FlmCurOpt(double** Points, long LenPoints, double Cur)
{
	std::vector<TVector3d> vArrayOfPoints(LenPoints);
	TVector3d* ArrayOfPoints = vArrayOfPoints.data();
	TVector3d* tArrayOfPoints = ArrayOfPoints;
	double** tPoints = Points;

	for(long k=0; k<LenPoints; k++)
	{
		double *aPoint = *(tPoints++);
		*(tArrayOfPoints++) = TVector3d(aPoint[0], aPoint[1], aPoint[2]);
	}

	rad.SetFlmCur(Cur, ArrayOfPoints, (int)LenPoints);
	// RAII: automatic cleanup
}

//-------------------------------------------------------------------------

void FlmCurDLL(double* Points, int LenPoints, double Cur)
{
	std::vector<TVector3d> vArrayOfPoints(LenPoints);
	TVector3d* ArrayOfPoints = vArrayOfPoints.data();
	TVector3d* tArrayOfPoints = ArrayOfPoints;
	double* tPoints = Points;

	for(long k=0; k<LenPoints; k++)
	{
		double x = *(tPoints++);
		double y = *(tPoints++);
		double z = *(tPoints++);
		*(tArrayOfPoints++) = TVector3d(radCR.Double(x), radCR.Double(y), radCR.Double(z));
	}
	rad.SetFlmCur(Cur, ArrayOfPoints, (int)LenPoints);
	// RAII: automatic cleanup
}

//-------------------------------------------------------------------------

void Rectngl(double xc, double yc, double zc, double Lx, double Ly)
{
	double CPoi[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	double Dims[] = {radCR.Double(Lx), radCR.Double(Ly)};
	rad.SetRectangle(CPoi, 3, Dims, 2);
}

//-------------------------------------------------------------------------

void BackgroundFieldSource(double Bx, double By, double Bz)
{
	double B[] = {Bx, By, Bz};
	rad.SetBackgroundFieldSource(B, 3);
}

void CoefficientFunctionFieldSource(PyObject* callback)
{
	rad.SetCoefficientFunctionFieldSource(callback);
}

//-------------------------------------------------------------------------

void Group(int* ArrayOfKeys, long lenArrayOfKeys)
{
	rad.SetGroup(ArrayOfKeys, lenArrayOfKeys);
}

//-------------------------------------------------------------------------

void AddToGroup(int GroupKey, int* ArrayOfKeys, long lenArrayOfKeys)
{
	rad.AddToGroup(GroupKey, ArrayOfKeys, lenArrayOfKeys);
}

//-------------------------------------------------------------------------

void OutGroupSize(int ElemKey)
//void OutGroupSize(int ElemKey, char deep=0)
{
	rad.OutGroupSize(ElemKey);
	//rad.OutGroupSize(ElemKey, deep);
}

//-------------------------------------------------------------------------

void OutGroupSubObjectKeys(int ElemKey)
{
	rad.OutGroupSubObjectKeys(ElemKey);
}

//-------------------------------------------------------------------------

void DuplicateElementG3D()
{

}

//-------------------------------------------------------------------------

void DuplicateElementG3DOpt(int ElemKey, const char* Opt)
{
	const char* OptionNames[] = {0};
	const char* OptionValues[] = {0};
	int OptionCount = 0;

	std::array<char, 200> CharBuf;
	if((Opt != 0) && (*Opt != '\0'))
	{
		strncpy(CharBuf.data(), Opt, 199);
		CharBuf[199] = '\0';
		char *pEndOptName = strrchr(CharBuf.data(), '-');
		if(pEndOptName == nullptr) { rad.Send.ErrorMessage("Radia::Error062"); return;}
		*pEndOptName = '\0';
		OptionNames[0] = CharBuf.data();
		OptionValues[0] = strrchr(Opt, '>') + 1;
		OptionCount++;
	}

	rad.DuplicateElement_g3d(ElemKey, OptionNames, OptionValues, OptionCount);
}

//-------------------------------------------------------------------------

void CreateFromG3DObjectWithSymmetries(int ElemKey)
{
	rad.CreateFromObj_g3dWithSym(ElemKey);
}

//-------------------------------------------------------------------------

void NumberOfDegOfFreedom(int ElemKey)
{
	rad.ComputeNumberOfDegOfFreedom(ElemKey);
}

//-------------------------------------------------------------------------

void MagnOfObj(int ElemKey)
{
	//rad.ComputeMagnInCenter(ElemKey);
	rad.ComputeMagnOrJ_InCenter(ElemKey, 'M');
}

//-------------------------------------------------------------------------

void ObjField(int ElemKey, char* fldType)
{
	rad.ComputeMagnOrJ_InCenter(ElemKey, *fldType);
}

//-------------------------------------------------------------------------

void ScaleCurInObj(int ElemKey, double scaleCoef)
{
	rad.ScaleCurrent(ElemKey, scaleCoef);
}

//-------------------------------------------------------------------------

void SetObjMagn(int ElemKey, double Mx, double My, double Mz)
{
	rad.SetObjMagn(ElemKey, Mx, My, Mz);
}

//-------------------------------------------------------------------------

void SubdivideElementG3D()
{

}

//-------------------------------------------------------------------------

void SubdivideElementG3DOpt(int ElemKey, double* SubdivArray, char TypeExtraSpec, double* ExtraSpec, int LenExtraSpec, const char* Opt1, const char* Opt2, const char* Opt3)
{
	const char* OptionNames[] = {0,0,0};
	const char* OptionValues[] = {0,0,0};
	int OptionCount = 0;

	std::array<char, 200> CharBuf1, CharBuf2, CharBuf3;
	if(Opt1 != 0)
	{
		if(*Opt1 != '\0')
		{
			strncpy(CharBuf1.data(), Opt1, 199);
			CharBuf1[199] = '\0';
			char *pEndOptName = strrchr(CharBuf1.data(), '-');
			if(pEndOptName == nullptr) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			*pEndOptName = '\0';
			OptionNames[0] = CharBuf1.data();
			OptionValues[0] = strrchr(Opt1, '>') + 1;
			OptionCount++;
		}
	}
	if(Opt2 != 0)
	{
		if(*Opt2 != '\0')
		{
			strncpy(CharBuf2.data(), Opt2, 199);
			CharBuf2[199] = '\0';
			char *pEndOptName = strrchr(CharBuf2.data(), '-');
			if(pEndOptName == nullptr) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			*pEndOptName = '\0';
			OptionNames[OptionCount] = CharBuf2.data();
			OptionValues[OptionCount] = strrchr(Opt2, '>') + 1;
			OptionCount++;
		}
	}
	if(Opt3 != 0)
	{
		if(*Opt3 != '\0')
		{
			strncpy(CharBuf3.data(), Opt3, 199);
			CharBuf3[199] = '\0';
			char *pEndOptName = strrchr(CharBuf3.data(), '-');
			if(pEndOptName == nullptr) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			*pEndOptName = '\0';
			OptionNames[OptionCount] = CharBuf3.data();
			OptionValues[OptionCount] = strrchr(Opt3, '>') + 1;
			OptionCount++;
		}
	}

	rad.SubdivideElement_g3d(ElemKey, SubdivArray, 6, TypeExtraSpec, ExtraSpec, LenExtraSpec, OptionNames, OptionValues, OptionCount);
}

//-------------------------------------------------------------------------

void CutElementG3D()
{

}

//-------------------------------------------------------------------------

void CutElementG3DOpt1(int ElemKey, double x, double y, double z, double nx, double ny, double nz, const char* Opt)
{
	double PointOnPlane[] = {x,y,z};
	double PlaneNormal[] = {nx,ny,nz};

	char CharBuf[200];
	strncpy(CharBuf, Opt, 199);
	CharBuf[199] = '\0';
	char *pEndOptName = strrchr(CharBuf, '-');
	if(pEndOptName == nullptr) { rad.Send.ErrorMessage("Radia::Error062"); return;}
	*pEndOptName = '\0';

	const char* OptionNames[] = {CharBuf};
	const char* OptionValues[] = {0};
	OptionValues[0] = strrchr(Opt, '>') + 1;

	rad.CutElement_g3d(ElemKey, PointOnPlane, 3, PlaneNormal, 3, OptionNames, OptionValues, 1);
}

void CutElementG3DOpt0(int ElemKey, double x, double y, double z, double nx, double ny, double nz)
{
	double PointOnPlane[] = {x,y,z};
	double PlaneNormal[] = {nx,ny,nz};

	const char* OptionNames[] = {0};
	const char* OptionValues[] = {0};

	rad.CutElement_g3d(ElemKey, PointOnPlane, 3, PlaneNormal, 3, OptionNames, OptionValues, 0);
}

void CutElementG3DOpt(int ElemKey, double x, double y, double z, double nx, double ny, double nz, const char* Opt1)
{
	double PointOnPlane[] = {x,y,z};
	double PlaneNormal[] = {nx,ny,nz};

	const char* OptionNames[] = {0};
	const char* OptionValues[] = {0};
	int OptionCount = 0;

	std::array<char, 200> CharBuf1;
	if(*Opt1 != '\0')
	{
		strncpy(CharBuf1.data(), Opt1, 199);
		CharBuf1[199] = '\0';
		char *pEndOptName = strrchr(CharBuf1.data(), '-');
		if(pEndOptName == nullptr) { rad.Send.ErrorMessage("Radia::Error062"); return;}
		*pEndOptName = '\0';
		OptionNames[0] = CharBuf1.data();
		OptionValues[0] = strrchr(Opt1, '>') + 1;
		OptionCount++;
	}
	rad.CutElement_g3d(ElemKey, PointOnPlane, 3, PlaneNormal, 3, OptionNames, OptionValues, OptionCount);
}

//-------------------------------------------------------------------------

void SubdivideElementG3DByParPlanes()
{

}

//-------------------------------------------------------------------------

void GeometricalVolume(int ElemKey)
{
	rad.ComputeGeometricalVolume(ElemKey);
}

//-------------------------------------------------------------------------

void GeometricalLimits(int ElemKey)
{
	rad.ComputeGeometricalLimits(ElemKey);
}

//-------------------------------------------------------------------------

void FldCmpMetForSubdRecMag(int ElemKey, int Switch, int SubLevel)
{
	rad.FieldCompMethForSubdividedRecMag(ElemKey, Switch, SubLevel);
}

//-------------------------------------------------------------------------

void SetLocMgnInSbdRecMag() // May be removed
{

}

//-------------------------------------------------------------------------

void Translation(double vx, double vy, double vz)
{
	//double TranslArray[] = {radCR.Double(vx), radCR.Double(vy), radCR.Double(vz)};
	double TranslArray[] = {vx, vy, vz};
	rad.SetTranslation(TranslArray, 3);
}

//-------------------------------------------------------------------------

void Rotation(double xc, double yc, double zc, 
			  double vx, double vy, double vz, double Angle)
{
	//double PoiOnAxis[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	//double AxisVect[] = {radCR.Double(vx), radCR.Double(vy), radCR.Double(vz)};

	double PoiOnAxis[] = {xc, yc, zc};
	double AxisVect[] = {vx, vy, vz};

	rad.SetRotation(PoiOnAxis, 3, AxisVect, 3, Angle);
}

//-------------------------------------------------------------------------

void PlaneSym(double xc, double yc, double zc, 
			  double nx, double ny, double nz)
{
	//double PoiOnPlane[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	//double PlaneNormal[] = {radCR.Double(nx), radCR.Double(ny), radCR.Double(nz)};

	double PoiOnPlane[] = {xc, yc, zc};
	double PlaneNormal[] = {nx, ny, nz};

	rad.SetPlaneSym(PoiOnPlane, 3, PlaneNormal, 3, 1);
}

//-------------------------------------------------------------------------

void FieldInversion()
{
	rad.SetFieldInversion();
}

//-------------------------------------------------------------------------

void CombineTransformLeft(int ThisElemKey, int AnotherElemKey)
{
	rad.CombineTransformations(ThisElemKey, AnotherElemKey, 'L');
}

//-------------------------------------------------------------------------

void CombineTransformRight(int ThisElemKey, int AnotherElemKey)
{
	rad.CombineTransformations(ThisElemKey, AnotherElemKey, 'R');
}

//-------------------------------------------------------------------------

void ApplySymmetry(int g3dElemKey, int TransElemKey, int Multiplicity)
{
	rad.ApplySymmetry(g3dElemKey, TransElemKey, Multiplicity);
}

//-------------------------------------------------------------------------

void TransformObject(int g3dElemKey, int TransElemKey)
{
	rad.ApplySymmetry(g3dElemKey, TransElemKey, 1);	
}

//-------------------------------------------------------------------------

void LinearMaterial(double KsiPar, double KsiPer, double Mrx, double Mry, double Mrz)
{
	double KsiArray[] = {KsiPar, KsiPer};
	double RemMagnArray[] = {Mrx, Mry, Mrz};
	rad.SetLinearMaterial(KsiArray, 2, RemMagnArray, 3);
}

//-------------------------------------------------------------------------

void LinearMaterial2(double KsiPar, double KsiPer, double Mr)
{
	double KsiArray[] = {KsiPar, KsiPer};
	rad.SetLinearMaterial(KsiArray, 2, &Mr, 1);
}

//-------------------------------------------------------------------------

void MaterialStd(char* MatName, double Mr)
{
	rad.SetMaterialStd(MatName, Mr);
}

//-------------------------------------------------------------------------

void NonlinearIsotropMaterial(double Ms1, double Ms2, double Ms3, 
							  double ks1, double ks2, double ks3)
{
	double Ms[] = {Ms1, Ms2, Ms3};
	double ks[] = {ks1, ks2, ks3};
	rad.SetNonlinearIsotropMaterial(Ms, 3, ks, 3);
}

//-------------------------------------------------------------------------

void NonlinearIsotropMaterial2(double ks1, double Ms1, double ks2, double Ms2, double ks3, double Ms3)
{
	double Ms[] = {Ms1, Ms2, Ms3};
	double ks[] = {ks1, ks2, ks3};
	rad.SetNonlinearIsotropMaterial(Ms, 3, ks, 3);
}

//-------------------------------------------------------------------------

void NonlinearIsotropMaterial3()
{
	int lenArrayOfPoints2d;
	TVector2d* ArrayOfPoints2d = nullptr;
	if(!rad.Send.GetArrayOfVector2d(ArrayOfPoints2d, lenArrayOfPoints2d)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	rad.SetNonlinearIsotropMaterial(ArrayOfPoints2d, lenArrayOfPoints2d);

	if(ArrayOfPoints2d != nullptr) delete[] ArrayOfPoints2d;
}

//-------------------------------------------------------------------------

void NonlinearIsotropMaterial3Opt(double** HandM_Array, long LenHandM_Array)
{
	std::vector<TVector2d> vArrayOfPoints2d(LenHandM_Array);
	TVector2d* ArrayOfPoints2d = vArrayOfPoints2d.data();
	TVector2d* tArrayOfPoints2d = ArrayOfPoints2d;
	double** tHandM_Array = HandM_Array;

	for(long i=0; i<LenHandM_Array; i++)
	{
		tArrayOfPoints2d->x = **tHandM_Array;
		(tArrayOfPoints2d++)->y = (*tHandM_Array)[1];
		tHandM_Array++;
	}

	rad.SetNonlinearIsotropMaterial(ArrayOfPoints2d, (int)LenHandM_Array);
	// RAII: automatic cleanup
}

//-------------------------------------------------------------------------

void NonlinearLaminatedMaterialML()
{

}

//-------------------------------------------------------------------------

void NonlinearLaminatedMaterialFrm(double* pKsiMs1, double* pKsiMs2, double* pKsiMs3, double PackFactor, double* dN)
{
	int lenArrayOfPoints2d = 0;
	if(pKsiMs1 != 0) lenArrayOfPoints2d++;
	if(pKsiMs2 != 0) lenArrayOfPoints2d++;
	if(pKsiMs3 != 0) lenArrayOfPoints2d++;

	if(lenArrayOfPoints2d == 0) { rad.Send.ErrorMessage("Radia::Error000"); return;}

	std::vector<TVector2d> vArrayOfPoints2d(lenArrayOfPoints2d);
	TVector2d* ArrayOfPoints2d = vArrayOfPoints2d.data();

	if(pKsiMs1 != 0)
	{
		ArrayOfPoints2d->x = pKsiMs1[0];
		ArrayOfPoints2d->y = pKsiMs1[1];
	}
	if(pKsiMs2 != 0)
	{
		ArrayOfPoints2d[1].x = pKsiMs2[0];
		ArrayOfPoints2d[1].y = pKsiMs2[1];
	}
	if(pKsiMs3 != 0)
	{
		ArrayOfPoints2d[2].x = pKsiMs3[0];
		ArrayOfPoints2d[2].y = pKsiMs3[1];
	}
	rad.SetNonlinearLaminatedMaterial(ArrayOfPoints2d, lenArrayOfPoints2d, PackFactor, dN);
	// RAII: automatic cleanup
}

//-------------------------------------------------------------------------

void NonlinearLaminatedMaterialTab(double* pFlatMatDef, int AmOfMatPts, double PackFactor, double* dN)
{
	if((pFlatMatDef == 0) || (AmOfMatPts <= 3))
	{
		rad.Send.ErrorMessage("Radia::Error088");
	}
	std::vector<TVector2d> vArrayOfPoints2d(AmOfMatPts);
	TVector2d* ArrayOfPoints2d = vArrayOfPoints2d.data();

	TVector2d *tArrayOfPoints2d = ArrayOfPoints2d;
	double *tFlatMatDef = pFlatMatDef;
	for(int i=0; i<AmOfMatPts; i++)
	{
		tArrayOfPoints2d->x = *(tFlatMatDef++);
		(tArrayOfPoints2d++)->y = *(tFlatMatDef++);
	}
	rad.SetNonlinearLaminatedMaterial(ArrayOfPoints2d, AmOfMatPts, PackFactor, dN);
	// RAII: automatic cleanup
}

//-------------------------------------------------------------------------

void NonlinearAnisotropMaterial()
{

}

//-------------------------------------------------------------------------

void NonlinearAnisotropMaterialOpt0(double* pDataPar, int lenDataPar, double* pDataPer, int lenDataPer)
{
	rad.SetNonlinearAnisotropMaterial0(pDataPar, lenDataPar, pDataPer, lenDataPer);
}

//-------------------------------------------------------------------------

void NonlinearAnisotropMaterialOpt1(double** Par, double** Per)
{
	double KsiPar[4], KsiPer[4], MsPar[3], MsPer[3], Hc[2];
	double* Ksi[] = {KsiPar, KsiPer};
	double* Ms[] = {MsPar, MsPer};

	char DependenceIsNonlinear[] = {1,1};

	double **tPar = Par, **tPer = Per;
	for(int i=0; i<3; i++)
	{
		KsiPar[i] = (*tPar)[0]; MsPar[i] = (*tPar)[1];
		KsiPer[i] = (*tPer)[0]; MsPer[i] = (*tPer)[1];
		tPar++; tPer++;
	}
	KsiPar[3] = (*tPar)[0]; Hc[0] = (*tPar)[1];
	KsiPer[3] = (*tPer)[0]; Hc[1] = 0; //Hc[1] = (*tPer)[1];

	//rad.SetNonlinearAnisotropMaterial(Ksi, Ms, Hc, DependenceIsNonlinear);
	rad.SetNonlinearAnisotropMaterial(Ksi, Ms, Hc, 2, DependenceIsNonlinear);
}

//-------------------------------------------------------------------------

void NonlinearAnisotropMaterialOpt2(double** Par, double Per)
{
	double KsiPar[4], KsiPer[4], MsPar[3], MsPer[3], Hc[2];
	double* Ksi[] = {KsiPar, KsiPer};
	double* Ms[] = {MsPar, MsPer};

	char DependenceIsNonlinear[] = {1,0};

	double **tPar = Par;
	for(int i=0; i<3; i++)
	{
		KsiPar[i] = (*tPar)[0]; MsPar[i] = (*tPar)[1];
		tPar++;
	}
	KsiPar[3] = (*tPar)[0]; Hc[0] = (*tPar)[1];
	KsiPer[0] = Per; Hc[1] = 0;

	//rad.SetNonlinearAnisotropMaterial(Ksi, Ms, Hc, DependenceIsNonlinear);
	rad.SetNonlinearAnisotropMaterial(Ksi, Ms, Hc, 2, DependenceIsNonlinear);
}

//-------------------------------------------------------------------------

void NonlinearAnisotropMaterialOpt3(double Par, double** Per)
{
	double KsiPar[4], KsiPer[4], MsPar[3], MsPer[3], Hc[2];
	double* Ksi[] = {KsiPar, KsiPer};
	double* Ms[] = {MsPar, MsPer};

	char DependenceIsNonlinear[] = {0,1};

	double **tPer = Per;
	for(int i=0; i<3; i++)
	{
		KsiPer[i] = (*tPer)[0]; MsPer[i] = (*tPer)[1];
		tPer++;
	}
	KsiPar[0] = Par; Hc[0] = 0;
	KsiPer[3] = (*tPer)[0]; Hc[1] = 0; //Hc[1] = (*tPer)[1];

	//rad.SetNonlinearAnisotropMaterial(Ksi, Ms, Hc, DependenceIsNonlinear);
	rad.SetNonlinearAnisotropMaterial(Ksi, Ms, Hc, 2, DependenceIsNonlinear);
}

//-------------------------------------------------------------------------

void ApplyMaterial(int g3dRelaxElemKey, int MaterElemKey)
{
	rad.ApplyMaterial(g3dRelaxElemKey, MaterElemKey);
}

//-------------------------------------------------------------------------

void MvsH(int g3dRelaxOrMaterElemKey, char* MagnChar, double hx, double hy, double hz)
{
	double H[] = {hx, hy, hz};
	rad.ComputeMvsH(g3dRelaxOrMaterElemKey, MagnChar, H, 3);
}

//-------------------------------------------------------------------------

void PreRelax(int ElemKey, int SrcElemKey)
{
	rad.PreRelax(ElemKey, SrcElemKey);
}

//-------------------------------------------------------------------------

void ShowInteractMatrix(int InteractElemKey)
{
	rad.ShowInteractMatrix(InteractElemKey);
}

//-------------------------------------------------------------------------

void ShowInteractVector(int InteractElemKey, char* FieldVectID)
{
	rad.ShowInteractVector(InteractElemKey, FieldVectID);
}

//-------------------------------------------------------------------------

void ManualRelax(int InteractElemKey, int MethNo, int IterNumber, double RelaxParam)
{
	rad.MakeManualRelax(InteractElemKey, MethNo, IterNumber, RelaxParam);
}

//-------------------------------------------------------------------------

void AutoRelax()
{

}

//-------------------------------------------------------------------------

void AutoRelaxOpt(int InteractElemKey, double PrecOnMagnetiz, int MaxIterNumber, int MethNo, const char* Opt1)
{
	std::array<char, 200> CharBuf1;
	const char* OptionNames[] = {CharBuf1.data()};
	const char* OptionValues[] = {0};
	const char* NonParsedOpts[] = {Opt1};
	int OptionCount = 1;
	AuxParseOptionNamesAndValues(NonParsedOpts, OptionNames, OptionValues, OptionCount);

	rad.MakeAutoRelax(InteractElemKey, PrecOnMagnetiz, MaxIterNumber, MethNo, OptionNames, OptionValues, OptionCount);
}

//-------------------------------------------------------------------------

void UpdateSourcesForRelax(int InteractElemKey)
{
	rad.UpdateSourcesForRelax(InteractElemKey);
}

//-------------------------------------------------------------------------

void SolveGen(int ObjKey, double PrecOnMagnetiz, int MaxIterNumber, int MethNo)
{
	if(MethNo == 0) MethNo = 4; //Default method
	rad.SolveGen(ObjKey, PrecOnMagnetiz, MaxIterNumber, MethNo);
}

//-------------------------------------------------------------------------

void CompCriterium(double InAbsPrecB, double InAbsPrecA, double InAbsPrecB_int, 
				   double InAbsPrecFrc, double InAbsPrecTrjCoord, double InAbsPrecTrjAngle)
{
	rad.SetCompCriterium(InAbsPrecB, InAbsPrecA, InAbsPrecB_int, InAbsPrecFrc, InAbsPrecTrjCoord, InAbsPrecTrjAngle);
}

//-------------------------------------------------------------------------

void CompPrecision()
{

}

//-------------------------------------------------------------------------

void CompPrecisionOpt(const char* Opt1, const char* Opt2, const char* Opt3, const char* Opt4, const char* Opt5, const char* Opt6, const char* Opt7, const char* Opt8)
//void CompPrecisionOpt(const char* Opt1, const char* Opt2, const char* Opt3, const char* Opt4, const char* Opt5, const char* Opt6, const char* Opt7, const char* Opt8, const char* Opt9)
{
	const char* OptionNames[] = {0,0,0,0,0,0,0,0,0};
	double OptionValues[] = {0,0,0,0,0,0,0,0,0};
	int OptionCount = 0;

	int MaxAmOfOptions = 8;
	std::array<char, 200> CharBuf1, CharBuf2, CharBuf3, CharBuf4, CharBuf5, CharBuf6, CharBuf7, CharBuf8;
	char *TotCharBuf[] = {CharBuf1.data(), CharBuf2.data(), CharBuf3.data(), CharBuf4.data(), CharBuf5.data(), CharBuf6.data(), CharBuf7.data(), CharBuf8.data()};
	const char *InOpt[] = {Opt1, Opt2,  Opt3, Opt4, Opt5, Opt6, Opt7, Opt8};
	for(int i=0; i<MaxAmOfOptions; i++)
	{
		if(InOpt[i] != 0) //OC240108
		{
			if(*(InOpt[i]) != '\0')
			{
				strncpy(TotCharBuf[i], InOpt[i], 199);
				TotCharBuf[i][199] = '\0';
				//char *pEndOptName = strrchr(TotCharBuf[i], '-');
				char *pEndOptName = strrchr(TotCharBuf[i], '>'); //OC19122019
				//if(pEndOptName == nullptr) { rad.Send.ErrorMessage("Radia::Error062"); return;}
				if((pEndOptName == nullptr) || (*(--pEndOptName) != '-')) { rad.Send.ErrorMessage("Radia::Error062"); return;} //OC19122019
				*pEndOptName = '\0';
				OptionNames[i] = TotCharBuf[i];
				//OptionValues[i] = atof(strrchr(InOpt[i], '>') + 1);
				OptionValues[i] = atof(pEndOptName + 2); //OC19122019
				if(OptionValues[i] == 0) { rad.Send.ErrorMessage("Radia::Error057"); return;}
				OptionCount++;
			}
		}
	}
	rad.SetCompPrecisions(OptionNames, OptionValues, OptionCount);
}

//-------------------------------------------------------------------------

// Maybe to be removed later
void MultipoleThresholds(double a0, double a1, double a2, double a3)
{
	double MltplThresh[] = {a0, a1, a2, a3};
	rad.SetMltplThresh(MltplThresh); 
}

//-------------------------------------------------------------------------

void Field(int ElemKey, char* FieldChar, double x1, double y1, double z1, 
		   double x2, double y2, double z2, int Np, char* ShowArgFlag, double StrtArg)
{
	double StObsPoi[] = {radCR.Double(x1), radCR.Double(y1), radCR.Double(z1)};
	double FiObsPoi[] = {radCR.Double(x2), radCR.Double(y2), radCR.Double(z2)};
	rad.ComputeField(ElemKey, FieldChar, StObsPoi, 3, FiObsPoi, 3, Np, ShowArgFlag, radCR.Double(StrtArg));
}

//-------------------------------------------------------------------------
/**
void FieldArbitraryPointsStruct()
{
	int ElemKey;
	if(!rad.Send.GetInteger(ElemKey)) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	char* FieldChar;
	if(!rad.Send.GetString((const char*&)FieldChar)) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	//if(!rad.Send.GetString(FieldChar)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	radTVectorOfVector3d VectorOfVector3d;
	radTVectInputCell VectInputCell;
	//if(!rad.Send.GetArbitraryListOfVector3d(VectorOfVector3d, VectInputCell)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	int resListRead = rad.Send.GetArbitraryListOfVector3d(VectorOfVector3d, VectInputCell);
	if(!resListRead) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	else if(resListRead == 2) return; //OC29092010 calculation should not be performed because symbol was sent (?)

	for(auto& vec : VectorOfVector3d) //OC050504
	{
	TVector3d vSmallRand(radCR.Double(0), radCR.Double(0), radCR.Double(0));
	vec += vSmallRand;
	}

	rad.ComputeField(ElemKey, FieldChar, VectorOfVector3d, VectInputCell);
	rad.Send.DisownString(FieldChar);
}
**/

//This version is derived from the one of D. Hidas (18/04/2016).
//It works with Math. v.10 Plot[] (that eventually sends a symbolic variable to radFld).
//Attempts to place "getting" of a List of Points to a separate function in radTSend 
//lead to re-appearing the 
void FieldArbitraryPointsStruct(int ElemKey, char* FieldChar) 
{

}

//-------------------------------------------------------------------------

void FieldArbitraryPointsArray(long ElemKey, const char* FieldChar, double** Points, long LenPoints)
{
	char LocStr[50];
	strncpy(LocStr, FieldChar, 49);
	LocStr[49] = '\0';
	rad.ComputeField((int)ElemKey, LocStr, Points, LenPoints);
}

//-------------------------------------------------------------------------

void FieldInt(int ElemKey, char* IntID, char* FieldIntChar, double x1, double y1, double z1, double x2, double y2, double z2)
{
	double StPoi[] = {radCR.Double(x1), radCR.Double(y1), radCR.Double(z1)};
	double FiPoi[] = {radCR.Double(x2), radCR.Double(y2), radCR.Double(z2)};
	rad.ComputeFieldInt(ElemKey, IntID, FieldIntChar, StPoi, 3, FiPoi, 3);
}

//-------------------------------------------------------------------------

void FieldForce(int ElemKey, int ShapeElemKey)
{
	rad.ComputeFieldForce(ElemKey, ShapeElemKey);
}

//-------------------------------------------------------------------------

void FieldEnergy(int DestElemKey, int SourceElemKey, int kx, int ky, int kz)
{
	int SubdArray[] = {kx, ky, kz};
	rad.ComputeFieldEnergy(DestElemKey, SourceElemKey, SubdArray, 3);
}

//-------------------------------------------------------------------------

void FieldForceThroughEnergy(int DestElemKey, int SourceElemKey, char* ForceComponID, int kx, int ky, int kz)
{
	int SubdArray[] = {kx, ky, kz};
	rad.ComputeFieldForceThroughEnergy(DestElemKey, SourceElemKey, ForceComponID, SubdArray, 3);
}

//-------------------------------------------------------------------------

void FieldTorqueThroughEnergy(int DestElemKey, int SourceElemKey, char* TorqueComponID, double x0, double y0, double z0, int kx, int ky, int kz)
{
	int SubdArray[] = {kx, ky, kz};
	double TorqueCenPo[] = {x0, y0, z0};
	rad.ComputeFieldTorqueThroughEnergy(DestElemKey, SourceElemKey, TorqueComponID, SubdArray, 3, TorqueCenPo, 3);
}

//-------------------------------------------------------------------------

void ParticleTrajectory(int ElemKey, double E, double x0, double dxdy0, double z0, double dzdy0, double y0, double y1, int Np)
{
	rad.ComputeParticleTrajectory(ElemKey, E, radCR.Double(x0), dxdy0, radCR.Double(z0), dzdy0, radCR.Double(y0), radCR.Double(y1), Np);
}

//-------------------------------------------------------------------------

void FocusingPotential(int ElemKey, double x1, double y1, double z1, double x2, double y2, double z2, int Np)
{
	double StPoi[] = {radCR.Double(x1), radCR.Double(y1), radCR.Double(z1)};
	double FiPoi[] = {radCR.Double(x2), radCR.Double(y2), radCR.Double(z2)};
	rad.ComputeFocusPotent(ElemKey, StPoi, 3, FiPoi, 3, Np);
}

//-------------------------------------------------------------------------

void FocusingKickPer(int ElemKey, double x1, double y1, double z1, double nsx, double nsy, double nsz, double per, double nper, double n1x, double n1y, double n1z, double r1, int np1, double r2, int np2, const char* Comment, int nharm, int ns, double d1, double d2, const char* strKickUnit, double energyGeV, const char* strOutFormat)
{
	double P1[] = {radCR.Double(x1), radCR.Double(y1), radCR.Double(z1)};
	double Nlong[] = {radCR.Double(nsx), radCR.Double(nsy), radCR.Double(nsz)};
	double N1[] = {radCR.Double(n1x), radCR.Double(n1y), radCR.Double(n1z)};
	rad.ComputeFocusKickPer(ElemKey, P1, Nlong, per, nper, N1, r1, np1, r2, np2, Comment, nharm, ns, d1, d2, strKickUnit, energyGeV, strOutFormat);
}

//-------------------------------------------------------------------------

void FocusingKickPerFormStrRep(double* pKickData1, double* pKickData2, double* pBtE2Int, double* pCoordDir1, double* pCoordDir2, int np1, int np2, double per, int nper, const char* Comment)
{
	rad.ComposeFocusKickPerFormStrRep(pKickData1, pKickData2, pBtE2Int, pCoordDir1, pCoordDir2, np1, np2, per, nper, Comment);
}

//-------------------------------------------------------------------------

//void FocusingKick(int ElemKey, double x1,double y1,double z1, double nsx,double nsy,double nsz, double* dsArr,long Len_dsArr,int ns, double n1x,double n1y,double n1z, double r1,int np1,double r2,int np2, const char* Comment, double d1,double d2)
void FocusingKickML()
{

}

//-------------------------------------------------------------------------

void ShimSignature(int ElemKey, char* FldID, double vx, double vy, double vz, double x1, double y1, double z1, double x2, double y2, double z2, int Np, double vix, double viy, double viz)
{
	double V[] = {vx, vy, vz};
	double StPoi[] = {radCR.Double(x1), radCR.Double(y1), radCR.Double(z1)};
	double FiPoi[] = {radCR.Double(x2), radCR.Double(y2), radCR.Double(z2)};
	double Vi[] = {radCR.Double(vix), radCR.Double(viy), radCR.Double(viz)};
	rad.ComputeShimSignature(ElemKey, FldID, V, StPoi, FiPoi, Np, Vi);
}

//-------------------------------------------------------------------------

void TolForConvergence(double AbsRandMagnitude, double RelRandMagnitude, double ZeroRandMagnitude)
{
	rad.SetTolForConvergence(AbsRandMagnitude, RelRandMagnitude, ZeroRandMagnitude);
}

//-------------------------------------------------------------------------

void RandomizationOnOrOff(char* OnOrOff)
{
	rad.RandomizationOnOrOff(OnOrOff);
}

//-------------------------------------------------------------------------

void PhysicalUnits()
{
	rad.SetAndShowPhysUnits();
}

//-------------------------------------------------------------------------

//void DumpElem(int ElemKey, char* strFormat)
//{
//	rad.DumpElem(ElemKey, strFormat);
//}

void DumpElem()
{

}

//-------------------------------------------------------------------------

void DumpElemOpt(int* arKeys, int nKeys, const char* AscOrBin)
{
	rad.DumpElem(arKeys, nKeys, AscOrBin, true); //?
}

//-------------------------------------------------------------------------

void DumpElemParse()
{

}

//-------------------------------------------------------------------------

void DumpElemParseOpt(const unsigned char* sBytes, int nBytes)
{
	rad.DumpElemParse(sBytes, nBytes); //?
}

//-------------------------------------------------------------------------

void GenDump()
{
	rad.GenDump();
}

//-------------------------------------------------------------------------

void GraphicsForElemWithoutSymChilds(int ElemKey)
{
	rad.GraphicsForElem_g3d(ElemKey, 0);
}

//-------------------------------------------------------------------------

//void GraphicsForElemWithSymChilds(int ElemKey)
//{
//	rad.GraphicsForElem_g3d(ElemKey, 1);
//}

void GraphicsForElemWithSymChildsExt()
{

}

//-------------------------------------------------------------------------

int GraphicsForElemVTK(int ElemKey, const char* Opt1, const char* Opt2, const char* Opt3) //OC04112019 (requested by R.Nagler)
{
	std::array<char, 200> CharBuf1, CharBuf2, CharBuf3;
	const char* OptionNames[] ={ CharBuf1.data(), CharBuf2.data(), CharBuf3.data() };
	const char* OptionValues[] ={ 0, 0, 0 };
	const char* NonParsedOpts[] ={ Opt1, Opt2, Opt3 };
	int OptionCount = 3;
	AuxParseOptionNamesAndValues(NonParsedOpts, OptionNames, OptionValues, OptionCount);

	return rad.GraphicsForElem_g3d_VTK(ElemKey, OptionNames, OptionValues, OptionCount); //rad.GoObjGeometry in Rob's version
}

//-------------------------------------------------------------------------

void GraphicsForAllWithoutSymChilds()
{
	rad.GraphicsForAll_g3d(0);
}

//-------------------------------------------------------------------------

void GraphicsForAllWithSymChilds()
{
	rad.GraphicsForAll_g3d(1);
}

//-------------------------------------------------------------------------

void ApplyDrawAttrToElem(int ElemKey, double R, double G, double B, double Thickness)
{
	double RGB_col[] = {R, G, B};
	rad.ApplyDrawAttrToElem_g3d(ElemKey, RGB_col, 3, Thickness);
}

//-------------------------------------------------------------------------

void ApplyColorToElem(int ElemKey, double R, double G, double B)
{
	double RGB_col[] = {R, G, B};
	rad.ApplyDrawAttrToElem_g3d(ElemKey, RGB_col, 3);
}

//-------------------------------------------------------------------------

void RemoveDrawAttrFromElem(int ElemKey)
{
	rad.RemoveDrawAttrFromElem_g3d(ElemKey);
}

//-------------------------------------------------------------------------

void DeleteElement(int ElemKey)
{
	rad.DeleteElement(ElemKey);
}

//-------------------------------------------------------------------------

void DeleteAllElements1()
{
	rad.DeleteAllElements(1);
}

//-------------------------------------------------------------------------

void DeleteAllElements2()
{
	rad.DeleteAllElements(2);
}

//-------------------------------------------------------------------------

void InterruptTime(double t)
{
	radYield.YieldInit(t);
	rad.Send.Double(t);
}

//-------------------------------------------------------------------------

void RadiaVersion()
{
	rad.ReturnVersionID();
}

//-------------------------------------------------------------------------

void ReturnInput(double Input, int NumTimes)
{
	rad.ReturnInput(Input, NumTimes);
}

//-------------------------------------------------------------------------

void MemAllocMethForIntrctMatr(char* TotOrParts)
{
	rad.SetMemAllocMethForIntrctMatr(TotOrParts);
}

//-------------------------------------------------------------------------

void OutCommandsInfo()
{
	rad.Send.OrdinaryMessage("Radia::usage");
}

//-------------------------------------------------------------------------

void ProcMPI(const char* OnOrOff, double* arData=0, long* pnData=0, long* pRankFrom=0, long* pRankTo=0) //OC19032020
//void ProcMPI(const char* OnOrOff)
{
	rad.ProcMPI(OnOrOff, arData, pnData, pRankFrom, pRankTo); //OC19032020
	//rad.ProcMPI(OnOrOff);
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------