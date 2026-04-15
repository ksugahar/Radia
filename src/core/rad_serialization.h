/*-------------------------------------------------------------------------
*
* File name:      radsend.h
*
* Project:        RADIA
*
* Description:    Interface functions (data input / output)
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#ifndef __RADSEND_H
#define __RADSEND_H

#include "rad_auxiliary_structures.h"

#include <map>
#include <vector>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// Drawing attribute types (kept for binary serialization compatibility)
struct radRGB {
	double Red, Green, Blue;
	radRGB(double InRed =0, double InGreen =0, double InBlue =0)
	{
		Red = InRed; Green = InGreen; Blue = InBlue;
	}
};

class radTrans;
class radTField;

//-------------------------------------------------------------------------

class radTSend {
public:

	radTSend() {}

	static void ErrorMessage(const char*);
	void OrdinaryMessage(const char*);
	static void WarningMessage(const char*);

	void String(const char*);
	void ByteString(const unsigned char* MessageString, long len);
	void Vector3d(const TVector3d*);
	void Vector3d(const TVector3df*);
	void ArrayOfVector3d(const TVector3d*, int);
	void Matrix3d(const TMatrix3d*);
	void Matrix3d(const TMatrix3df*);
	void MatrixOfMatrix3d(TMatrix3d**, int, int);
	void MatrixOfMatrix3d(TMatrix3df**, int, int);
	void Long(long);
	void Int(int);
	void IntList(int*, int);
	void Double(double);
	void DoubleList(double*, int);
	void ArbNestedArrays(double*, int*, int);
	void SubArbNestedArrays(double*, int*, int, int&);

	void MultiDimArrayOfDouble(double*, int*, int);

	void ArrayOfPairOfVect3d(radTVectPairOfVect3d* pVectPairOfVect3d);
	void OutFieldForceOrTorqueThroughEnergyCompRes(char* ForceComponID, TVector3d& Vect, char ID);
	void OutFieldIntCompRes(char* FieldIntChar, radTField* FieldPtr, double* ArgArray = 0, int Np = 1);
	void OutFieldCompRes(char* FieldChar, radTField* FieldArray, double* ArgArray, int Np);
	void OutRelaxResultsInfo(double* RelaxStatusParamArray, int lenRelaxStatusParamArray, int ActualIterNum);
	void OutMagnetizCompRes(char* MagnChar, TVector3d& M_vect);

	void MyMLPutDouble(double);

	void InitOutList(int);

	int GetInteger(int&);
	int GetDouble(double&);
	int GetString(const char*&);
	void DisownString(char* Str);
	int GetArbitraryListOfVector3d(radTVectorOfVector3d&, radTVectInputCell&);

	int GetVector3d(TVector3d& vect3d);
	int GetVector2d(TVector2d& vect2d);

	int GetArrayOfVector3d(TVector3d*&, int&);
	int GetArrayOfVector2d(TVector2d*&, int&);
	int GetArrayOfVector2dVersion2(TVector2d*&, int&);
	int GetArrayOfArrayOfVector3d(TVector3d**&, int*&, int&);
	int GetArrayOfArrayOfInt(int**&, int*&, int&);

	int GetArrayOfDouble(double*&, long&);
};

//-------------------------------------------------------------------------

#endif

