/*-------------------------------------------------------------------------
*
* File name:      rad_obj_subdivision.cpp
*
* Project:        RADIA
*
* Description:    Legacy Triangle-based subdivision entry points
*
* The Triangle library is no longer bundled with this repository.  These
* stubs keep the legacy C/Python API link-compatible while directing users to
* the Netgen/Cubit mesh workflows used by the maintained code path.
*
-------------------------------------------------------------------------*/

#include "rad_application.h"

namespace {
constexpr const char* kTriangleRemovedError = "Radia::Error126";
}

//-------------------------------------------------------------------------

int radTApplication::SetMultGenExtrTriangle(
	double* FirstPoi,
	long lenFirstPoi,
	double Lx,
	TVector2d* ArrayOfPoints2d,
	long lenArrayOfPoints2d,
	double* arSubdData,
	double* Magn,
	long lenMagn,
	const char* OrientStr,
	const char** OptionNames,
	const char** OptionValues,
	int OptionCount)
{
	(void)FirstPoi;
	(void)lenFirstPoi;
	(void)Lx;
	(void)ArrayOfPoints2d;
	(void)lenArrayOfPoints2d;
	(void)arSubdData;
	(void)Magn;
	(void)lenMagn;
	(void)OrientStr;
	(void)OptionNames;
	(void)OptionValues;
	(void)OptionCount;

	Send.ErrorMessage(kTriangleRemovedError);
	return 0;
}

//-------------------------------------------------------------------------

int radTApplication::TriangulatePolygon(
	TVector2d* ArrayOfPoints2d,
	long lenArrayOfPoints2d,
	double* arSubdData,
	char triSubdParamBorderCode,
	double triAngMin,
	double triAreaMax,
	const char* sTriExtOpt,
	TVector2d*& arTriVertPt,
	int& numTriVertPt,
	int*& arTriVertInd,
	int& numTri)
{
	(void)ArrayOfPoints2d;
	(void)lenArrayOfPoints2d;
	(void)arSubdData;
	(void)triSubdParamBorderCode;
	(void)triAngMin;
	(void)triAreaMax;
	(void)sTriExtOpt;

	arTriVertPt = nullptr;
	numTriVertPt = 0;
	arTriVertInd = nullptr;
	numTri = 0;

	Send.ErrorMessage(kTriangleRemovedError);
	return 0;
}

//-------------------------------------------------------------------------
