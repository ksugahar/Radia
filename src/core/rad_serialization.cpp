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
*
-------------------------------------------------------------------------*/

#include "rad_serialization.h"
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

void radTSend::InitOutList(int NumberOfElem)
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

// Graphics3D legacy code removed (FrameLines, DrawPyramidArrow, DrawCharacter)

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

// DeallocateGeomPolygonData removed (Graphics3D/VTK geometry code)

//-------------------------------------------------------------------------
