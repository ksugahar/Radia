/*-------------------------------------------------------------------------
*
* File name:      radopnam.h
*
* Project:        RADIA
*
* Description:    Option names for Radia / Mathematica functions
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#ifndef __RADOPNAM_H
#define __RADOPNAM_H

#include <string.h>
#include <vector>
#include <map>

//-------------------------------------------------------------------------

struct radTOptionNames {
	char B[25], A[25], BInt[25], Force[25], Torque[25], Energy[25], Coord[25], Angle[25]; // Precisions

	char Frame[25], FrameValues[3][25];
	char SubdParamCode[25], SubdParamBorderCode[25], SubdParamCodeValues[2][25];
	char SubdCoils[25], SubdCoilsValues[4][25];
	char TriAngMin[25], TriAreaMax[25], TriExtOpt[25]; //ki->Numb|Size,TriAngMin->...,TriAreaMax->...,TriExtOpt->\"...\"

	char ShowLines[25], ShowLinesValues[4][25];
	char ShowFaces[25], ShowFacesValues[4][25];
	char ShowFrameAxes[25], ShowFrameAxesValues[4][25];
	char ShowSymChilds[25], ShowSymChildsValues[4][25];

	char FreeSym[25], FreeSymValues[4][25];
	char ZeroM[25], ZeroM_Values[4][25];

	char LinTreat[25]; //, LinCoefValues[4][25];
	char Debug[25];

	//map<string, vector<pair<string, int> > > mOptData;
	map<string, map<string, int> > mOptData;

	radTOptionNames()
	{
		//partially obsolete definitions:
		strncpy(B, "PrcB", 24); B[24] = '\0';
		strncpy(A, "PrcA", 24); A[24] = '\0';
		strncpy(BInt, "PrcBInt", 24); BInt[24] = '\0';
		strncpy(Force, "PrcForce", 24); Force[24] = '\0';
		strncpy(Torque, "PrcTorque", 24); Torque[24] = '\0';
		strncpy(Energy, "PrcEnergy", 24); Energy[24] = '\0';
		strncpy(Coord, "PrcCoord", 24); Coord[24] = '\0';
		strncpy(Angle, "PrcAngle", 24); Angle[24] = '\0';

		strncpy(Frame, "Frame", 24); Frame[24] = '\0';
		strncpy(FrameValues[0], "Loc", 24); FrameValues[0][24] = '\0';
		strncpy(FrameValues[1], "LabTot", 24); FrameValues[1][24] = '\0';
		strncpy(FrameValues[2], "Lab", 24); FrameValues[2][24] = '\0';

		strncpy(SubdParamCode, "kxkykz", 24); SubdParamCode[24] = '\0';
		strncpy(SubdParamBorderCode, "ki", 24); SubdParamBorderCode[24] = '\0';
		strncpy(SubdParamCodeValues[0], "Numb", 24); SubdParamCodeValues[0][24] = '\0';
		strncpy(SubdParamCodeValues[1], "Size", 24); SubdParamCodeValues[1][24] = '\0';

		strncpy(SubdCoils, "DivCoils", 24); SubdCoils[24] = '\0';
		strncpy(SubdCoilsValues[0], "No", 24); SubdCoilsValues[0][24] = '\0';
		strncpy(SubdCoilsValues[1], "Yes", 24); SubdCoilsValues[1][24] = '\0';
		strncpy(SubdCoilsValues[2], "False", 24); SubdCoilsValues[2][24] = '\0';
		strncpy(SubdCoilsValues[3], "True", 24); SubdCoilsValues[3][24] = '\0';

		strncpy(ShowLines, "EdgeLines", 24); ShowLines[24] = '\0';
		strncpy(ShowLinesValues[0], "No", 24); ShowLinesValues[0][24] = '\0';
		strncpy(ShowLinesValues[1], "Yes", 24); ShowLinesValues[1][24] = '\0';
		strncpy(ShowLinesValues[2], "False", 24); ShowLinesValues[2][24] = '\0';
		strncpy(ShowLinesValues[3], "True", 24); ShowLinesValues[3][24] = '\0';

		strncpy(ShowFaces, "Faces", 24); ShowFaces[24] = '\0';
		strncpy(ShowFacesValues[0], "No", 24); ShowFacesValues[0][24] = '\0';
		strncpy(ShowFacesValues[1], "Yes", 24); ShowFacesValues[1][24] = '\0';
		strncpy(ShowFacesValues[2], "False", 24); ShowFacesValues[2][24] = '\0';
		strncpy(ShowFacesValues[3], "True", 24); ShowFacesValues[3][24] = '\0';

		strncpy(ShowFrameAxes, "Axes", 24); ShowFrameAxes[24] = '\0';
		strncpy(ShowFrameAxesValues[0], "No", 24); ShowFrameAxesValues[0][24] = '\0';
		strncpy(ShowFrameAxesValues[1], "Yes", 24); ShowFrameAxesValues[1][24] = '\0';
		strncpy(ShowFrameAxesValues[2], "False", 24); ShowFrameAxesValues[2][24] = '\0';
		strncpy(ShowFrameAxesValues[3], "True", 24); ShowFrameAxesValues[3][24] = '\0';

		strncpy(ShowSymChilds, "ShowSym", 24); ShowSymChilds[24] = '\0';
		strncpy(ShowSymChildsValues[0], "No", 24); ShowSymChildsValues[0][24] = '\0';
		strncpy(ShowSymChildsValues[1], "Yes", 24); ShowSymChildsValues[1][24] = '\0';
		strncpy(ShowSymChildsValues[2], "False", 24); ShowSymChildsValues[2][24] = '\0';
		strncpy(ShowSymChildsValues[3], "True", 24); ShowSymChildsValues[3][24] = '\0';

		strncpy(FreeSym, "FreeSym", 24); FreeSym[24] = '\0';
		strncpy(FreeSymValues[0], "No", 24); FreeSymValues[0][24] = '\0';
		strncpy(FreeSymValues[1], "Yes", 24); FreeSymValues[1][24] = '\0';
		strncpy(FreeSymValues[2], "False", 24); FreeSymValues[2][24] = '\0';
		strncpy(FreeSymValues[3], "True", 24); FreeSymValues[3][24] = '\0';

		strncpy(ZeroM, "ZeroM", 24); ZeroM[24] = '\0';
		strncpy(ZeroM_Values[0], "No", 24); ZeroM_Values[0][24] = '\0';
		strncpy(ZeroM_Values[1], "Yes", 24); ZeroM_Values[1][24] = '\0';
		strncpy(ZeroM_Values[2], "False", 24); ZeroM_Values[2][24] = '\0';
		strncpy(ZeroM_Values[3], "True", 24); ZeroM_Values[3][24] = '\0';

		strncpy(TriAngMin, "TriAngMin", 24); TriAngMin[24] = '\0';
		strncpy(TriAreaMax, "TriAreaMax", 24); TriAreaMax[24] = '\0';
		strncpy(TriExtOpt, "TriExtOpt", 24); TriExtOpt[24] = '\0';

		strncpy(LinTreat, "Lin", 24); LinTreat[24] = '\0';
		//strncpy(LinCoefValues[0], "Rel", 24); LinCoefValues[0][24] = '\0';
		//strncpy(LinCoefValues[1], "Loc", 24); LinCoefValues[1][24] = '\0';
		//strncpy(LinCoefValues[2], "Abs", 24); LinCoefValues[2][24] = '\0';
		//strncpy(LinCoefValues[3], "Lab", 24); LinCoefValues[3][24] = '\0';

		strncpy(Debug, "Debug", 24); Debug[24] = '\0';

		//newer definitions
		map<string, int> vRealVal;
		vRealVal["d"] = 0;
		mOptData[B] = vRealVal;
		mOptData[A] = vRealVal;
		mOptData[BInt] = vRealVal;
		mOptData[Force] = vRealVal;
		mOptData[Torque] = vRealVal;
		mOptData[Energy] = vRealVal;
		mOptData[Coord] = vRealVal;
		mOptData[Angle] = vRealVal;
		mOptData[TriAngMin] = vRealVal;
		mOptData[TriAreaMax] = vRealVal;

		map<string, int> vStringVal;
		vStringVal["s"] = 0;
		mOptData[TriExtOpt] = vStringVal;

		map<string, int> vFrameVal;
		vFrameVal["Loc"] = 0;
		vFrameVal["LabTot"] = 1;
		vFrameVal["Lab"] = 2;
		mOptData[Frame] = vFrameVal;

		map<string, int> vSubdPar;
		vSubdPar["Numb"] = 0;
		vSubdPar["Size"] = 1;
		mOptData[SubdParamCode] = vSubdPar;
		mOptData[SubdParamBorderCode] = vSubdPar;

		map<string, int> vNoYes;
		vNoYes["No"] = 0;
		vNoYes["Yes"] = 1;
		vNoYes["False"] = 0;
		vNoYes["True"] = 1;
		mOptData[SubdCoils] = vNoYes;
		mOptData[ShowLines] = vNoYes;
		mOptData[ShowFaces] = vNoYes;
		mOptData[ShowFrameAxes] = vNoYes;
		mOptData[ShowSymChilds] = vNoYes;
		mOptData[FreeSym] = vNoYes;
		mOptData[ZeroM] = vNoYes;
		mOptData[Debug] = vNoYes;

		map<string, int> vRelAbs;
		vRelAbs["Rel"] = 0;
		vRelAbs["Abs"] = 1;
		vRelAbs["Loc"] = 0;
		vRelAbs["Lab"] = 1;
		mOptData[LinTreat] = vRelAbs;
	}

	bool parseOption(const char* optName, const char* optValueIn, char& cOptValueOut, double& dOptValueOut, char* sOptValueOut, char& outCode)
	{
		map<string, map<string, int> >::const_iterator it = mOptData.find(optName);
		if(it == mOptData.end()) return false;

		const map<string, int> &mapCurValues = it->second;
		int numCurValues = (int)mapCurValues.size();
		if(numCurValues <= 0) return false;

		outCode = 0;

		if(numCurValues == 1)
		{
			const char* c_strKey = mapCurValues.begin()->first.c_str();
			if((*c_strKey == 'd') || (*c_strKey == 'r'))
			{
				if(optValueIn != 0) dOptValueOut = atof(optValueIn);

				outCode = 'd';
				return true;
			}
			else if(*c_strKey == 's')
			{//return without parsing optValueIn
				if(optValueIn != 0)
				{
					if(sOptValueOut != 0) strcpy(sOptValueOut, optValueIn);
				}

				outCode = 's';
				return true;
			}
		}

		if(optValueIn != 0)
		{
			map<string, int>::const_iterator itVal = mapCurValues.find(optValueIn);
			if(itVal == mapCurValues.end()) return false;
			cOptValueOut = (char)(itVal->second);
		}

		outCode = 'c';
		return true;
	}

	bool findParseOptionValues(const char** arAllOptionNamesIn, const char** arAllOptionValuesIn, int numAllOptIn, const char** arOptNamesToFind, int numOptToFind, char* cArOptValsFoundParsed, double* dArOptValsFoundParsed, char** sArOptValsFoundParsed)
	{
		if((numAllOptIn == 0) || (arAllOptionNamesIn == 0) || (arAllOptionValuesIn == 0)) return true;
		if((arOptNamesToFind == 0) || (numOptToFind == 0)) return true;
		if((cArOptValsFoundParsed == 0) && (dArOptValsFoundParsed == 0) && (sArOptValsFoundParsed)) return true;

		char cAux, sAux[1024], cOutCode;
		double dAux;
		sAux[0] = '\0';

		for(int i=0; i<numAllOptIn; i++)
		{
			//if(cArOptValsFoundParsed != 0) p_c = cArOptValsFoundParsed + countOptFound;
			//if(dArOptValsFoundParsed != 0) p_d = dArOptValsFoundParsed + countOptFound;
			//if(sArOptValsFoundParsed != 0) p_s = sArOptValsFoundParsed[countOptFound];

			const char *sCurOptionNameIn = arAllOptionNamesIn[i];
			cOutCode = 0;
			if(!parseOption(sCurOptionNameIn, arAllOptionValuesIn[i], cAux, dAux, sAux, cOutCode)) return false;

			int cOptCount=-1, dOptCount=-1, sOptCount=-1;

			for(int j=0; j<numOptToFind; j++)
			{
				char auxOutCode=0;
				if(!parseOption(arOptNamesToFind[j], 0, cAux, dAux, sAux, auxOutCode)) return false;
				if(auxOutCode == 'c') cOptCount++;
				else if(auxOutCode == 'd') dOptCount++;
				else if(auxOutCode == 's') sOptCount++;

				if(strcmp(sCurOptionNameIn, arOptNamesToFind[j]) == 0)
				{
					if(cArOptValsFoundParsed != 0) 
					{
						//if(cOutCode == 'c') cArOptValsFoundParsed[j] = cAux;
						if(cOutCode == 'c') cArOptValsFoundParsed[cOptCount] = cAux;
					}
					if(dArOptValsFoundParsed != 0) 
					{
						//if(cOutCode == 'd') dArOptValsFoundParsed[j] = dAux;
						if(cOutCode == 'd') dArOptValsFoundParsed[dOptCount] = dAux;
					}
					if(sArOptValsFoundParsed != 0) 
					{
						//if(cOutCode == 's') strcpy(sArOptValsFoundParsed[j], sAux);
						if(cOutCode == 's') strcpy(sArOptValsFoundParsed[sOptCount], sAux);
					}
				}
			}
		}
		return true;
	}
};

//-------------------------------------------------------------------------

#endif
