/*-------------------------------------------------------------------------
*
* File name:      radsend.cpp
*
* Project:        RADIA
*
* Description:    Interface functions
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radsend.h"
#include <stdio.h>
#include <string.h>
#include <vector>

//#ifdef __JAVA__
//#ifndef __SEND2JAVA_H
//#include "Send2Java.h"
//#endif
//extern CSendToJava gSendToJava;
//#endif

//#ifdef __DLLVBA__
//#ifndef __SEND2VBA_H
//#include "Send2VBA.h"
//#endif
//extern radTSendToVBA gSendToVBA;
//#endif

#include "rad_io_buffer.h"
extern radTIOBuffer ioBuffer;
#include "rad_graphics_3d.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTSend::ErrorMessage(const char* MessageString)
{
//#ifdef __JAVA__
//	gSendToJava.SendErrorMessage(MessageString);
//#endif
//#ifdef __DLLVBA__
//	gSendToVBA.SendErrorMessage(MessageString);
//#endif
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	ioBuffer.StoreErrorMessage(MessageString);
#endif
}

//-------------------------------------------------------------------------

void radTSend::WarningMessage(const char* MessageString)
{
//#ifdef __JAVA__
//	gSendToJava.SendWarningMessage(MessageString);
//#endif
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	ioBuffer.StoreWarningMessage(MessageString);
#endif
}

//-------------------------------------------------------------------------

void radTSend::OrdinaryMessage(const char* MessageString)
{
}

//-------------------------------------------------------------------------

void radTSend::String(const char* MessageString)
{
#ifdef __JAVA__
	gSendToJava.SendString(MessageString);
#endif
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	ioBuffer.StoreString(MessageString);
#endif
}

//-------------------------------------------------------------------------

void radTSend::ByteString(const unsigned char* MessageString, long len)
{
//#ifdef __JAVA__
//	gSendToJava.SendString(MessageString);
//#endif
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	ioBuffer.StoreByteString((const char*)MessageString, len);
#endif
}

//-------------------------------------------------------------------------

void radTSend::Double(double d)
{
#ifdef __JAVA__
	gSendToJava.SendDouble(d);
#endif
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	ioBuffer.StoreDouble(d);
#endif
}

//-------------------------------------------------------------------------

void radTSend::MyMLPutDouble(double d)
{
}

//-------------------------------------------------------------------------

void radTSend::DoubleList(double* ArrayOfDouble, int lenArrayOfDouble)
{
//#ifdef __JAVA__
#if defined __JAVA__ || defined ALPHA__DLL__ || defined ALPHA__LIB__
	int Dims[] = {lenArrayOfDouble};
	MultiDimArrayOfDouble(ArrayOfDouble, Dims, 1);
#endif
}

//-------------------------------------------------------------------------

void radTSend::Long(long LongIntValue)
{
#ifdef __JAVA__
	gSendToJava.SendLong(LongIntValue);
#endif
}

//-------------------------------------------------------------------------

void radTSend::Int(int IntValue)
{
#ifdef __JAVA__
	gSendToJava.SendInt(IntValue);
#endif
#ifdef __DLLVBA__
	gSendToVBA.SendInt(IntValue);
#endif
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	ioBuffer.StoreInt(IntValue);
#endif
}

//-------------------------------------------------------------------------

void radTSend::IntList(int* ArrayOfInt, int lenArrayOfInt)
{
#ifdef __JAVA__
	int Dims[] = { lenArrayOfInt};
	gSendToJava.SendMultiDimArrayOfInt(ArrayOfInt, Dims, 1);
#endif
#ifdef __DLLVBA__
	int Dims[] = { lenArrayOfInt};
	gSendToVBA.SendMultiDimArrayOfInt(ArrayOfInt, Dims, 1);
#endif
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	int Dims[] = { lenArrayOfInt};
	ioBuffer.StoreMultiDimArrayOfInt(ArrayOfInt, Dims, 1);
#endif
}

//-------------------------------------------------------------------------

void radTSend::InitOutList(int NumberOfElem, char DrawFacilityInd)
{
}

//-------------------------------------------------------------------------

void radTSend::Vector3d(const TVector3d* VectorPtr)
{
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	double TotOutArray[] = {VectorPtr->x, VectorPtr->y, VectorPtr->z};
	int Dims[] = {3};
	MultiDimArrayOfDouble(TotOutArray, Dims, 1);
#endif
}

//-------------------------------------------------------------------------

void radTSend::Vector3d(const TVector3df* VectorPtr)
{
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	double TotOutArray[] = {VectorPtr->x, VectorPtr->y, VectorPtr->z};
	int Dims[] = {3};
	MultiDimArrayOfDouble(TotOutArray, Dims, 1);
#endif
}

//-------------------------------------------------------------------------

void radTSend::ArrayOfVector3d(const TVector3d* ArrayOfVector3d, int lenArray)
{
}

//-------------------------------------------------------------------------

void radTSend::Matrix3d(const TMatrix3d* MatrixPtr)
{
}

//-------------------------------------------------------------------------

void radTSend::Matrix3d(const TMatrix3df* MatrixPtr)
{
}

//-------------------------------------------------------------------------

void radTSend::MatrixOfMatrix3d(TMatrix3d** MatrixOfMatrix3d, int AmOfStr, int AmOfCol)
{
}

//-------------------------------------------------------------------------

void radTSend::MatrixOfMatrix3d(TMatrix3df** MatrixOfMatrix3d, int AmOfStr, int AmOfCol)
{
}

//-------------------------------------------------------------------------

void radTSend::SubArbNestedArrays(double* Data, int* Dims, int Depth, int& CntData)
{
}

//-------------------------------------------------------------------------

void radTSend::ArbNestedArrays(double* Data, int* Dims, int Depth)
{
//#ifdef __JAVA__
#if defined __JAVA__ || defined ALPHA__DLL__ || defined ALPHA__LIB__
	MultiDimArrayOfDouble(Data, Dims, Depth);
#endif
}

//-------------------------------------------------------------------------

void radTSend::Color(const radRGB& RGB_color, char DrawFacilityInd)
{// This is not for Lines !!!
}

//-------------------------------------------------------------------------

void radTSend::GenInitDraw(char DrawFacilityInd)
{
}

//-------------------------------------------------------------------------

void radTSend::InitDrawSurfElem(int DrawAttrAreSet, radTDrawAttr& DrawAttr, int NumberOfSymChilds_pl_Orig, char DrawFacilityInd)
{
}

//-------------------------------------------------------------------------

void radTSend::InitDrawLinElem(int DrawAttrAreSet, radTDrawAttr& DrawAttr, int NumberOfSymChilds_pl_Orig, char DrawFacilityInd)
{
}

//-------------------------------------------------------------------------

void radTSend::InitDrawLineWithThickness(int DrawAttrAreSet, radTDrawAttr& DrawAttr, char DrawFacilityInd)
{
}

//-------------------------------------------------------------------------

void radTSend::DrawEdgesSuppression(char DrawFacilityInd)
{
}

//-------------------------------------------------------------------------

void radTSend::Polygon(const TVector3d* Side, int lenSide, char DrawFacilityInd)
{
	if(DrawFacilityInd == 2) // VTK export
	{
		if(!ShowFaces) return;
		AddGeomPolygon(Side, lenSide, GeomPolygons);
	}
}

//-------------------------------------------------------------------------

void radTSend::AddGeomPolygon(const TVector3d* Side, int lenSide, radTVectGeomPolygon& VectGeomPolygons)
{
	if((Side == 0) || (lenSide == 0)) return;

	TVector3d* tSide = (TVector3d*)Side;

	radTGeomPolygon aNewGeomPolygon;
	aNewGeomPolygon.vVertCoords.resize(lenSide*3);
	aNewGeomPolygon.VertCoords = aNewGeomPolygon.vVertCoords.data();

	double* tVertCoords = aNewGeomPolygon.VertCoords;

	for(int i = 0; i < lenSide; i++)
	{
		*(tVertCoords++) = tSide->x;
		*(tVertCoords++) = tSide->y;
		*(tVertCoords++) = tSide->z;

		if(tSide->x < Limits3D[0]) Limits3D[0] = tSide->x;
		else if(tSide->x > Limits3D[1]) Limits3D[1] = tSide->x;
		if(tSide->y < Limits3D[2]) Limits3D[2] = tSide->y;
		else if(tSide->y > Limits3D[3]) Limits3D[3] = tSide->y;
		if(tSide->z < Limits3D[4]) Limits3D[4] = tSide->z;
		else if(tSide->z > Limits3D[5]) Limits3D[5] = tSide->z;

		tSide++;
	}

	aNewGeomPolygon.Nv = lenSide;

	if(!(radTg3dGraphPresent::DrawAttrStack).empty())
	{
		radRGB& RGB = ((radTg3dGraphPresent::DrawAttrStack).begin())->RGB_col;
		aNewGeomPolygon.ColRGB[0] = float(RGB.Red); 
		aNewGeomPolygon.ColRGB[1] = float(RGB.Green); 
		aNewGeomPolygon.ColRGB[2] = float(RGB.Blue);
	}

	VectGeomPolygons.push_back(aNewGeomPolygon);
}

//-------------------------------------------------------------------------

void radTSend::Line(const TVector3d* EdgePoints, int lenEdgePoints, char DrawFacilityInd)
{
	if(DrawFacilityInd == 2)
	{
		if(!ShowLines)
		{
			return;
		}

		AddGeomPolygon(EdgePoints, lenEdgePoints, GeomLines);
	}
}

//-------------------------------------------------------------------------

int radTSend::GeomDataToBuffer() //OC04112019
//int radTSend::GeomDataToBuffer(int& nVP, int& nP, int& nVL, int& nL, int& key) //OC04112019
{
	int key = 0;
	int nP = (int)(GeomPolygons.size());
	int nL = (int)(GeomLines.size());
	if((nP > 0) || (nL > 0))
	{
		key = rand(); //to improve?
		if(nP > 0) ioBuffer.StoreGeomPolygData(GeomPolygons, key);
		if(nL > 0) ioBuffer.StoreGeomPolygData(GeomLines, key + 1);
	}
	return key;
}

//-------------------------------------------------------------------------
/**
void radTSend::FrameLines(char DrawFacilityInd)
{//Don't use this ! Use radTg3dGraphPresent::DrawFrameLines() instead !
	const double RelFrameLinesOffset = 0.07;

	const double RelArrowHeight = 0.08;
	const double ArrowBottomRadToHeightRatio = 0.15;

	const double RelCharHeight = 0.05;
	const double CharWidthToHeigthRatio = 0.6;

	radRGB FrameLinesColor(1,1,1); //Frame lines are white
	radTDrawAttr FrameLinesDrawAttr;
	FrameLinesDrawAttr.RGB_col = FrameLinesColor;
	radTg3dGraphPresent::DrawAttrStack.push_back(FrameLinesDrawAttr);

	double OrigLimits3D[6];
	double *tOrigLimits3D = OrigLimits3D, *tLimits3D = Limits3D;
	for(int k=0; k<6; k++) *(tOrigLimits3D++) = *(tLimits3D++);

	double Sx = OrigLimits3D[1] - OrigLimits3D[0];
	double Sy = OrigLimits3D[3] - OrigLimits3D[2];
	double Sz = OrigLimits3D[5] - OrigLimits3D[4];
	double Smax = (Sx > Sy)? ((Sx > Sz)? Sx : Sz) : ((Sy > Sz)? Sy : Sz);

	double Dx = Sx*RelFrameLinesOffset;
	double Dy = Sy*RelFrameLinesOffset;
	double Dz = Sz*RelFrameLinesOffset;

// Frame contour lines
	TVector3d Side[5], Segment[2];
	Side[0].x = OrigLimits3D[0] - Dx; Side[0].y = OrigLimits3D[2] - Dy; Side[0].z = OrigLimits3D[4] - Dz; 
	Side[1].x = Side[0].x; Side[1].y = Side[0].y; Side[1].z = OrigLimits3D[5] + Dz; 
	Side[2].x = Side[0].x; Side[2].y = OrigLimits3D[3] + Dy; Side[2].z = Side[1].z; 
	Side[3].x = Side[0].x; Side[3].y = Side[2].y; Side[3].z = Side[0].z;
	Side[4] = Side[0];
	Segment[0] = Side[0];
	Line(Side, 5, DrawFacilityInd);
	Side[0].x = OrigLimits3D[1] + Dx;
	Side[4].x = Side[3].x = Side[2].x = Side[1].x = Side[0].x;
	Line(Side, 5, DrawFacilityInd);

	Segment[1] = Side[0];
	Line(Segment, 2, DrawFacilityInd);
	Segment[0].y = Side[1].y; Segment[0].z = Side[1].z;
	Segment[1] = Side[1];
	Line(Segment, 2, DrawFacilityInd);
	Segment[0].y = Side[2].y; Segment[0].z = Side[2].z;
	Segment[1] = Side[2];
	Line(Segment, 2, DrawFacilityInd);
	Segment[0].y = Side[3].y; Segment[0].z = Side[3].z;
	Segment[1] = Side[3];
	Line(Segment, 2, DrawFacilityInd);

// Frame arrows
	double ArrowHeight = RelArrowHeight*Smax;
	double ArrowBottomRad = ArrowBottomRadToHeightRatio*ArrowHeight;

	TVector3d PyramidArrowInfo[4];
	TVector3d ArrowHeightVect(ArrowHeight,0.,0.);
	TVector3d ArrowBottomRad1(0., -ArrowBottomRad, 0.), ArrowBottomRad2(0., 0., -ArrowBottomRad);
	if(ArrowHeight < Sx)
	{
		PyramidArrowInfo[0].x = 0.5*(OrigLimits3D[0] + OrigLimits3D[1]) + 0.5*ArrowHeight;
		PyramidArrowInfo[0].y = OrigLimits3D[2] - Dy;
		PyramidArrowInfo[0].z = OrigLimits3D[4] - Dz;
		PyramidArrowInfo[1] = ArrowHeightVect;
		PyramidArrowInfo[2] = ArrowBottomRad1; PyramidArrowInfo[3] = ArrowBottomRad2;
		DrawPyramidArrow(PyramidArrowInfo, DrawFacilityInd);
	}
	if(ArrowHeight < Sy)
	{
		PyramidArrowInfo[0].x = OrigLimits3D[1] + Dx;
		PyramidArrowInfo[0].y = 0.5*(OrigLimits3D[2] + OrigLimits3D[3]) + 0.5*ArrowHeight;
		PyramidArrowInfo[0].z = OrigLimits3D[4] - Dz;
		PyramidArrowInfo[1].x = 0.; PyramidArrowInfo[1].y = ArrowHeight; PyramidArrowInfo[1].z = 0.;
		PyramidArrowInfo[2].x = ArrowBottomRad; PyramidArrowInfo[2].y = PyramidArrowInfo[2].z = 0.;
		PyramidArrowInfo[3].x = PyramidArrowInfo[3].y = 0.; PyramidArrowInfo[3].z = -ArrowBottomRad;
		DrawPyramidArrow(PyramidArrowInfo, DrawFacilityInd);
	}
	if(ArrowHeight < Sz)
	{
		PyramidArrowInfo[0].x = OrigLimits3D[1] + Dx;
		PyramidArrowInfo[0].y = OrigLimits3D[3] + Dy;
		PyramidArrowInfo[0].z = 0.5*(OrigLimits3D[4] + OrigLimits3D[5]) + 0.5*ArrowHeight;
		PyramidArrowInfo[1].x = PyramidArrowInfo[1].y = 0.; PyramidArrowInfo[1].z = ArrowHeight;
		PyramidArrowInfo[2].x = ArrowBottomRad; PyramidArrowInfo[2].y = PyramidArrowInfo[2].z = 0.;
		PyramidArrowInfo[3].x = 0.; PyramidArrowInfo[3].y = ArrowBottomRad; PyramidArrowInfo[3].z = 0.;
		DrawPyramidArrow(PyramidArrowInfo, DrawFacilityInd);
	}

// Frame characters
	double AbsCharHeigth = RelCharHeight*Smax;
	double AbsCharWidth = CharWidthToHeigthRatio*AbsCharHeigth;

	TVector3d CharInfo3D[3];
	CharInfo3D[0].x = 0.5*(OrigLimits3D[0] + OrigLimits3D[1]) - 0.5*AbsCharWidth;
	CharInfo3D[0].y = OrigLimits3D[2] - Dy;
	CharInfo3D[0].z = OrigLimits3D[4] - Dz - 1.5*AbsCharHeigth;
	CharInfo3D[1].x = 0.; CharInfo3D[1].y = -1.; CharInfo3D[1].z = 0.;
	CharInfo3D[2].x = 0.; CharInfo3D[2].y = 0.; CharInfo3D[2].z = AbsCharHeigth;
	DrawCharacter('X', CharWidthToHeigthRatio, CharInfo3D, DrawFacilityInd);

	CharInfo3D[0].x = OrigLimits3D[1] + Dx;
	CharInfo3D[0].y = 0.5*(OrigLimits3D[2] + OrigLimits3D[3]) - 0.5*AbsCharWidth;
	CharInfo3D[0].z = OrigLimits3D[4] - Dz - 1.5*AbsCharHeigth;
	CharInfo3D[1].x = 1.; CharInfo3D[1].y = 0.; CharInfo3D[1].z = 0.;
	CharInfo3D[2].x = 0.; CharInfo3D[2].y = 0.; CharInfo3D[2].z = AbsCharHeigth;
	DrawCharacter('Y', CharWidthToHeigthRatio, CharInfo3D, DrawFacilityInd);

	CharInfo3D[0].x = OrigLimits3D[1] + Dx;
	CharInfo3D[0].y = OrigLimits3D[3] + Dy + AbsCharWidth;
	CharInfo3D[0].z = 0.5*(OrigLimits3D[4] + OrigLimits3D[5]) - 0.5*AbsCharHeigth;
	CharInfo3D[1].x = 1.; CharInfo3D[1].y = 0.; CharInfo3D[1].z = 0.;
	CharInfo3D[2].x = 0.; CharInfo3D[2].y = 0.; CharInfo3D[2].z = AbsCharHeigth;
	DrawCharacter('Z', CharWidthToHeigthRatio, CharInfo3D, DrawFacilityInd);

	radTg3dGraphPresent::DrawAttrStack.pop_back();
}
**/
//-------------------------------------------------------------------------
/**
void radTSend::DrawPyramidArrow(TVector3d* PyramidArrowInfo, char DrawFacilityInd)
{//Don't use this ! Use radTg3dGraphPresent::DrawPyramidArrow() instead !
// PyramidArrowInfo[0] - main vertex of pyramid
// PyramidArrowInfo[1] - height vector of pyramid
// PyramidArrowInfo[2] - bottom radius_1 vect
// PyramidArrowInfo[3] - bottom radius_2 vect

	TVector3d ArrowVerices[5], Face[4];
	*ArrowVerices = *PyramidArrowInfo;
	TVector3d ArrowOrigin = PyramidArrowInfo[0] - PyramidArrowInfo[1];
	TVector3d& ArrowBottomRad1 = PyramidArrowInfo[2];
	TVector3d& ArrowBottomRad2 = PyramidArrowInfo[3];
	ArrowVerices[1] = ArrowOrigin + ArrowBottomRad1;
	ArrowVerices[2] = ArrowOrigin + ArrowBottomRad2;
	ArrowVerices[3] = ArrowOrigin - ArrowBottomRad1;
	ArrowVerices[4] = ArrowOrigin - ArrowBottomRad2;
	Face[0] = ArrowVerices[0]; Face[1] = ArrowVerices[1]; Face[2] = ArrowVerices[2];
	Polygon(Face, 3, DrawFacilityInd);
	Face[1] = ArrowVerices[2]; Face[2] = ArrowVerices[3];
	Polygon(Face, 3, DrawFacilityInd);
	Face[1] = ArrowVerices[3]; Face[2] = ArrowVerices[4];
	Polygon(Face, 3, DrawFacilityInd);
	Face[1] = ArrowVerices[4]; Face[2] = ArrowVerices[1];
	Polygon(Face, 3, DrawFacilityInd);
	Face[0] = ArrowVerices[1]; Face[1] = ArrowVerices[4]; Face[2] = ArrowVerices[3]; Face[3] = ArrowVerices[2];
	Polygon(Face, 4, DrawFacilityInd);
}
**/
//-------------------------------------------------------------------------
/**
void radTSend::DrawCharacter(char Ch, double Ratio, TVector3d* Info3D, char DrawFacilityInd)
{//Don't use this ! Use radTg3dGraphPresent::DrawCharacter() instead !
// Info3D[0], Info3D[1] - point and normal vector of a plane in which the character lies
// Info3D[2] - height vector of a character (the Height is the length of Info3D[2])
// Width = Ratio*Height

	TVector3d& Normal = Info3D[1];
	Normal = (1./sqrt(Normal.x*Normal.x + Normal.y*Normal.y + Normal.z*Normal.z))*Normal;

	TVector3d LowerLeft = *Info3D;
	TVector3d WidthVect = Ratio*(Info3D[2]^Normal);
	TVector3d LowerRight = LowerLeft + WidthVect;
	TVector3d UpperLeft = LowerLeft + Info3D[2];
	TVector3d UpperRight = UpperLeft + WidthVect;

	TVector3d Segment[4];

	if((Ch == 'X') || (Ch == 'x'))
	{
		Segment[0] = LowerLeft; Segment[1] = UpperRight;
		Line(Segment, 2, DrawFacilityInd);
		Segment[0] = UpperLeft; Segment[1] = LowerRight;
		Line(Segment, 2, DrawFacilityInd);
	}
	else if((Ch == 'Y') || (Ch == 'y'))
	{
		TVector3d Middle = 0.25*(LowerLeft + UpperLeft + LowerRight + UpperRight);
		Segment[0] = UpperLeft; Segment[1] = Middle; Segment[2] = UpperRight;
		Line(Segment, 3, DrawFacilityInd);
		Segment[0] = 0.5*(LowerLeft + LowerRight); Segment[1] = Middle;
		Line(Segment, 2, DrawFacilityInd);
	}
	else if((Ch == 'Z') || (Ch == 'z'))
	{
		Segment[0] = UpperLeft; Segment[1] = UpperRight; Segment[2] = LowerLeft; Segment[3] = LowerRight;
		Line(Segment, 4, DrawFacilityInd);
	}
}
**/
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

int radTSend::GetArrayOfDouble(double*& Data, long& lenData)
{
	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetArrayOfVector3d(TVector3d*& ArrayOfVector3d, int& lenArrayOfVector3d)
{
	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetVector3d(TVector3d& vect3d)
{
	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetVector2d(TVector2d& vect2d)
{
	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetArrayOfVector2d(TVector2d*& ArrayOfVector2d, int& lenArrayOfVector2d)
{
	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetArrayOfVector2dVersion2(TVector2d*& ArrayOfVector2d, int& lenArrayOfVector2d)
{
	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetArrayOfArrayOfVector3d(TVector3d**& ArrayOfArrayOfVector3d, int*& ArrayOfLengths, int& lenArrayOfArrayOfVector3d)
{
	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetArrayOfArrayOfInt(int**& ArrayOfArrayOfInt, int*& ArrayOfLengths, int& lenArrayOfArrayOfInt)
{
	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetInteger(int& Value)
{
	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetDouble(double& Value)
{
	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetString(const char*& Str)
//int radTSend::GetString(char*& Str)
{
	return 1;
}

//-------------------------------------------------------------------------

void radTSend::DisownString(char* Str)
{
}

//-------------------------------------------------------------------------

int radTSend::GetArbitraryListOfVector3d(radTVectorOfVector3d& VectorOfVector3d, radTVectInputCell& VectInputCell)
{
	return 1;
}

//-------------------------------------------------------------------------

void radTSend::MultiDimArrayOfDouble(double* Array, int* Dims, int NumDims)
{
#ifdef __JAVA__
	gSendToJava.SendMultiDimArrayOfDouble(Array, Dims, NumDims);
#endif
#ifdef __DLLVBA__
	gSendToVBA.SendMultiDimArrayOfDouble(Array, Dims, NumDims);
#endif
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	ioBuffer.StoreMultiDimArrayOfDouble(Array, Dims, NumDims);
#endif
}

//-------------------------------------------------------------------------

void radTSend::ArrayOfPairOfVect3d(radTVectPairOfVect3d* pVectPairOfVect3d)
{
//#ifdef __JAVA__
#if defined __JAVA__ || defined ALPHA__DLL__ || defined ALPHA__LIB__

	int AmOfPoints = (int)pVectPairOfVect3d->size();
	int NumDims = 3;
	int Dims[] = {3,2,AmOfPoints};

	long TotLen = Dims[0]*Dims[1]*Dims[2];
	std::vector<double> vTotArray(TotLen);
	double *TotArray = vTotArray.data();
	double *tTotArray = TotArray;
	for(int k=0; k<AmOfPoints; k++)
	{
		radTPairOfVect3d& aPair = (*pVectPairOfVect3d)[k];
		TVector3d &V1 = aPair.V1, &V2 = aPair.V2;
		*(tTotArray++) = V1.x; *(tTotArray++) = V1.y; *(tTotArray++) = V1.z;
		*(tTotArray++) = V2.x; *(tTotArray++) = V2.y; *(tTotArray++) = V2.z;
	}
	MultiDimArrayOfDouble(TotArray, Dims, NumDims);

#endif
}

//-------------------------------------------------------------------------

void radTSend::OutFieldForceOrTorqueThroughEnergyCompRes(char* ForceComponID, TVector3d& Vect, char ID)
{// This is only for Force and Torque!
	char* BufChar = ForceComponID;
	//char* EqEmptyStr = (ID=='f')? "FxFyFz" : "TxTyTz";
	//char EqEmptyStr[6];
	char EqEmptyStr[10]; //OC150505
	strncpy(EqEmptyStr, "TxTyTz", 9); EqEmptyStr[9] = '\0';
	if(ID=='f') { strncpy(EqEmptyStr, "FxFyFz", 9); EqEmptyStr[9] = '\0'; }

	char SmallID = ID;
	char CapitalID = (SmallID=='f')? 'F' : 'T';

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0')
		{
			char* BufChar_pl_1 = BufChar+1;
			if((((*BufChar==CapitalID) || (*BufChar==SmallID)) &&
			   (*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z')) ||
			   (*BufChar == 'X') || (*BufChar == 'x') ||
			   (*BufChar == 'Y') || (*BufChar == 'y') ||
			   (*BufChar == 'Z') || (*BufChar == 'z')) ItemCount++;
			BufChar++;
		}
		BufChar = ForceComponID;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 3;
	}

//#ifdef __JAVA__
#if defined __JAVA__ || defined ALPHA__DLL__ || defined ALPHA__LIB__

	double TotOutArray[10];
	double *t = TotOutArray;
	int nv = 0;

	while(*BufChar != '\0')
	{
		if((*(BufChar)==CapitalID) || (*(BufChar)==SmallID))
		{
			char* BufChar_pl_1 = BufChar+1;
			if((*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z'))
			{ *(t++) = Vect.x; *(t++) =Vect.y; *(t++) = Vect.z; nv += 3;}
		}
		else if((*(BufChar)=='X') || (*(BufChar)=='x')) { *(t++) = Vect.x; nv++;}
		else if((*(BufChar)=='Y') || (*(BufChar)=='y')) { *(t++) = Vect.y; nv++;}
		else if((*(BufChar)=='Z') || (*(BufChar)=='z')) { *(t++) = Vect.z; nv++;}
		BufChar++;
	}
	int Dims[] = { nv};
	MultiDimArrayOfDouble(TotOutArray, Dims, 1);
#endif
}

//-------------------------------------------------------------------------

void radTSend::OutFieldCompRes(char* FieldChar, radTField* FieldArray, double* ArgArray, int Np)
{
	char* BufChar = FieldChar;
	//char* EqEmptyStr = "BHAM";
	char EqEmptyStr[] = "BHAM"; //OC01052013

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while(*BufChar != '\0') 
		{
			if((*BufChar == 'B') || (*BufChar == 'b') || 
			   (*BufChar == 'H') || (*BufChar == 'h') ||
			   (*BufChar == 'A') || (*BufChar == 'a') ||
			   (*BufChar == 'M') || (*BufChar == 'm') ||
			   (*BufChar == 'J') || (*BufChar == 'j') ||
			   (*BufChar == 'P') || (*BufChar == 'p')) ItemCount++;
			BufChar++;
		}
		BufChar = FieldChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 4;
	}
	char* ActualInitCharPtr = BufChar;

//#ifdef __JAVA__
#if defined __JAVA__ || defined ALPHA__DLL__ || defined ALPHA__LIB__

	std::vector<double> vTotOutArray(14*Np);
	double *TotOutArray = vTotOutArray.data();
	double *t = TotOutArray;
	int nv = 0;

	radTField* FieldPtr = FieldArray;
	for(int i=0; i<Np; i++)
	{
		nv = 0;
		if(ArgArray != nullptr) // Argument Needed
		{
			*(t++) = ArgArray[i]; nv++;
		}

		while(*BufChar != '\0') 
		{
			char* BufChar_p_1 = BufChar+1;
			if(*(BufChar)=='B' || *(BufChar)=='b')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->B.x; nv++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->B.y; nv++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->B.z; nv++;}
				else { *(t++) = FieldPtr->B.x; *(t++) = FieldPtr->B.y; *(t++) = FieldPtr->B.z; nv += 3;}
			}
			else if(*(BufChar)=='H' || *(BufChar)=='h')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->H.x; nv++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->H.y; nv++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->H.z; nv++;}
				else { *(t++) = FieldPtr->H.x; *(t++) = FieldPtr->H.y; *(t++) = FieldPtr->H.z; nv += 3;}
			}
			else if(*(BufChar)=='A' || *(BufChar)=='a')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->A.x; nv++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->A.y; nv++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->A.z; nv++;}
				else { *(t++) = FieldPtr->A.x; *(t++) = FieldPtr->A.y; *(t++) = FieldPtr->A.z; nv += 3;}
			}
			else if(*(BufChar)=='M' || *(BufChar)=='m')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->M.x; nv++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->M.y; nv++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->M.z; nv++;}
				else { *(t++) = FieldPtr->M.x; *(t++) = FieldPtr->M.y; *(t++) = FieldPtr->M.z; nv += 3;}
			}
			else if(*(BufChar)=='J' || *(BufChar)=='j')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->J.x; nv++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->J.y; nv++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->J.z; nv++;}
				else { *(t++) = FieldPtr->J.x; *(t++) = FieldPtr->J.y; *(t++) = FieldPtr->J.z; nv += 3;}
			}
			else if(*(BufChar)=='P' || *(BufChar)=='p')	{ *(t++) = FieldPtr->Phi; nv++;}
			BufChar++;
		}
		FieldPtr++;
		BufChar = ActualInitCharPtr;
	}
	int Dims[] = { nv, Np};
	MultiDimArrayOfDouble(TotOutArray, Dims, 2);

	// RAII: automatic cleanup via vTotOutArray
#endif
}

//-------------------------------------------------------------------------

//void radTSend::OutFieldIntCompRes(char* FieldIntChar, radTField* FieldPtr, double* ArgArray, int Np)
void radTSend::OutFieldIntCompRes(char* FieldIntChar, radTField* FieldArray, double* ArgArray, int Np)
{
	char* BufChar = FieldIntChar;
	char* BufCharPrev = nullptr;
	//char* EqEmptyStr = "Ib";
	char EqEmptyStr[] = "Ib"; //OC01052013

	short I_used = 0;
	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			if(((*BufChar == 'B') || (*BufChar == 'b') || 
			    (*BufChar == 'H') || (*BufChar == 'h')) ||
			   (((*BufChar == 'X') || (*BufChar == 'x') ||
			     (*BufChar == 'Y') || (*BufChar == 'y') ||
				 (*BufChar == 'Z') || (*BufChar == 'z')) &&
				(*BufCharPrev != 'B') && (*BufCharPrev != 'b') &&
				(*BufCharPrev != 'H') && (*BufCharPrev != 'h'))) ItemCount++;

			if((*BufChar == 'I') || (*BufChar == 'i')) I_used = 1;
			BufCharPrev = BufChar;
			BufChar++;
		}
		BufChar = FieldIntChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 1;
	}
	if(I_used && (ItemCount == 0))
	{
		BufChar = EqEmptyStr;
		ItemCount = 1;
	}

	char* ActualInitCharPtr = BufChar;

#if defined __JAVA__ || defined ALPHA__DLL__ || defined ALPHA__LIB__

	//double TotOutArray[10];
	//double *t = TotOutArray;
	//int nv = 0;

	std::vector<double> vTotOutArray(10*Np);
	double *TotOutArray = vTotOutArray.data();
	double *t = TotOutArray;
	int nv = 0;

	radTField* FieldPtr = FieldArray;
	for(int i=0; i<Np; i++)
	{
		nv = 0;
		if(ArgArray != nullptr) // Argument Needed
		{
			*(t++) = ArgArray[i]; nv++;
		}

		while(*BufChar != '\0') 
		{
			char* BufChar_pl_1 = BufChar+1;
			char* BufChar_mi_1 = BufChar-1;

			if((*BufChar =='I') || (*BufChar == 'i'))
			{
				if((*BufChar_pl_1 == 'X') || (*BufChar_pl_1 == 'x')) { *(t++) = FieldPtr->Ib.x; nv++;}
				else if((*BufChar_pl_1 == 'Y') || (*BufChar_pl_1 == 'y')) { *(t++) = FieldPtr->Ib.y; nv++;}
				else if((*BufChar_pl_1 == 'Z') || (*BufChar_pl_1 == 'z')) { *(t++) = FieldPtr->Ib.z; nv++;}
				else if((*BufChar_pl_1 != 'B') && (*BufChar_pl_1 != 'b') &&
					(*BufChar_pl_1 != 'H') && (*BufChar_pl_1 != 'h') &&
					(*BufChar_pl_1 != 'X') && (*BufChar_pl_1 != 'x') &&
					(*BufChar_pl_1 != 'Y') && (*BufChar_pl_1 != 'y') &&
					(*BufChar_pl_1 != 'Z') && (*BufChar_pl_1 != 'z')) 
				{ *(t++) = FieldPtr->Ib.x; *(t++) = FieldPtr->Ib.y; *(t++) = FieldPtr->Ib.z; nv += 3; break;}
			}
			else if((*BufChar == 'B') || (*BufChar == 'b'))
			{
				if((*BufChar_pl_1 == 'X') || (*BufChar_pl_1 == 'x')) { *(t++) = FieldPtr->Ib.x; nv++;}
				else if((*BufChar_pl_1 == 'Y') || (*BufChar_pl_1 == 'y')) { *(t++) = FieldPtr->Ib.y; nv++;}
				else if((*BufChar_pl_1 == 'Z') || (*BufChar_pl_1 == 'z')) { *(t++) = FieldPtr->Ib.z; nv++;}
				else { *(t++) = FieldPtr->Ib.x; *(t++) = FieldPtr->Ib.y; *(t++) = FieldPtr->Ib.z; nv += 3;}
			}
			else if((*BufChar == 'H') || (*BufChar == 'h'))
			{
				if((*BufChar_pl_1 == 'X') || (*BufChar_pl_1 == 'x')) { *(t++) = FieldPtr->Ih.x; nv++;}
				else if((*BufChar_pl_1 == 'Y') || (*BufChar_pl_1 == 'y')) { *(t++) = FieldPtr->Ih.y; nv++;}
				else if((*BufChar_pl_1 == 'Z') || (*BufChar_pl_1 == 'z')) { *(t++) = FieldPtr->Ih.z; nv++;}
				else { *(t++) = FieldPtr->Ih.x; *(t++) = FieldPtr->Ih.y; *(t++) = FieldPtr->Ih.z; nv += 3;}
			}
			else if(((*BufChar == 'X') || (*BufChar == 'x')) &&
				(*BufChar_mi_1 != 'I') && (*BufChar_mi_1 != 'i') &&
				(*BufChar_mi_1 != 'B') && (*BufChar_mi_1 != 'b') &&
				(*BufChar_mi_1 != 'H') && (*BufChar_mi_1 != 'h')) { *(t++) = FieldPtr->Ib.x; nv++;}
			else if(((*BufChar == 'Y') || (*BufChar == 'y')) &&
				(*BufChar_mi_1 != 'I') && (*BufChar_mi_1 != 'i') &&
				(*BufChar_mi_1 != 'B') && (*BufChar_mi_1 != 'b') &&
				(*BufChar_mi_1 != 'H') && (*BufChar_mi_1 != 'h')) { *(t++) = FieldPtr->Ib.y; nv++;}
			else if(((*BufChar == 'Z') || (*BufChar == 'z')) &&
				(*BufChar_mi_1 != 'I') && (*BufChar_mi_1 != 'i') &&
				(*BufChar_mi_1 != 'B') && (*BufChar_mi_1 != 'b') &&
				(*BufChar_mi_1 != 'H') && (*BufChar_mi_1 != 'h')) { *(t++) = FieldPtr->Ib.z; nv++;}
				BufChar++;
		}
		FieldPtr++;
		BufChar = ActualInitCharPtr;
	}

	//int Dims[] = { nv};
	//MultiDimArrayOfDouble(TotOutArray, Dims, 1);

	int Dims[] = { nv, Np};
	MultiDimArrayOfDouble(TotOutArray, Dims, 2);

	// RAII: automatic cleanup via vTotOutArray
#endif
}

//-------------------------------------------------------------------------

void radTSend::OutRelaxResultsInfo(double* RelaxStatusParamArray, int lenRelaxStatusParamArray, int ActualIterNum)
{
//#ifdef __JAVA__
#if defined __JAVA__ || defined ALPHA__DLL__ || defined ALPHA__LIB__
	int TotOutElem = lenRelaxStatusParamArray + 1;
	std::vector<double> vTotOutArray(TotOutElem);
	double *TotOutArray = vTotOutArray.data();
	double *t = TotOutArray;
	double *tRelaxStatusParamArray = RelaxStatusParamArray;
	for(int i=0; i<lenRelaxStatusParamArray; i++) *(t++) = *(tRelaxStatusParamArray++);
	*t = ActualIterNum;

	int Dims[] = { TotOutElem};
	MultiDimArrayOfDouble(TotOutArray, Dims, 1);

	// RAII: automatic cleanup via vTotOutArray
#endif
}

//-------------------------------------------------------------------------

void radTSend::OutMagnetizCompRes(char* MagnChar, TVector3d& M_vect)
{
	char* BufChar = MagnChar;
	//char* EqEmptyStr = "MxMyMz";
	char EqEmptyStr[] = "MxMyMz";

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			char* BufChar_pl_1 = BufChar+1;
			if((((*BufChar == 'M') || (*BufChar == 'm')) && 
			   (*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z')) ||
			   (*BufChar == 'X') || (*BufChar == 'x') ||
			   (*BufChar == 'Y') || (*BufChar == 'y') ||
			   (*BufChar == 'Z') || (*BufChar == 'z')) ItemCount++;
			BufChar++;
		}
		BufChar = MagnChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 3;
	}

//#ifdef __JAVA__
#if defined __JAVA__ || defined ALPHA__DLL__ || defined ALPHA__LIB__

	double TotOutArray[10];
	double *t = TotOutArray;
	int nv = 0;

	while (*BufChar != '\0') 
	{
		if((*(BufChar)=='M') || (*(BufChar)=='m'))
		{
			char* BufChar_pl_1 = BufChar+1;
			if((*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z'))
			{ *(t++) = M_vect.x; *(t++) =M_vect.y; *(t++) = M_vect.z; nv += 3;}
		}
		else if((*(BufChar)=='X') || (*(BufChar)=='x')) { *(t++) = M_vect.x; nv++;}
		else if((*(BufChar)=='Y') || (*(BufChar)=='y')) { *(t++) = M_vect.y; nv++;}
		else if((*(BufChar)=='Z') || (*(BufChar)=='z')) { *(t++) = M_vect.z; nv++;}
		BufChar++;
	}
	int Dims[] = { nv};
	MultiDimArrayOfDouble(TotOutArray, Dims, 1);
#endif
}

//-------------------------------------------------------------------------

void radTSend::DeallocateGeomPolygonData()
{
	int AmOfGeomPolygons = (int)GeomPolygons.size();
	int AmOfGeomLines = (int)GeomLines.size();

	if((AmOfGeomPolygons == 0) && (AmOfGeomLines == 0)) return;

	int k;
	for(k=0; k<AmOfGeomPolygons; k++)
	{
		double* pCoords = GeomPolygons[k].VertCoords;
		if(pCoords != 0) delete[] pCoords;
	}
	GeomPolygons.clear();

	for(k=0; k<AmOfGeomLines; k++)
	{
		double* pCoords = GeomLines[k].VertCoords;
		if(pCoords != 0) delete[] pCoords;
	}
	GeomLines.clear();
}

//-------------------------------------------------------------------------
