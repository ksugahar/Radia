/*-------------------------------------------------------------------------
*
* File name:      radauxst.h
*
* Project:        RADIA
*
* Description:    Auxiliary structures
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADAUXST_H
#define __RADAUXST_H

#ifdef WIN32
#pragma warning(disable : 4786)
#endif

#include "gmvect.h"
#include "gmvectf.h"

#include <vector>
#include <map>
#include <iostream>

#ifdef _MSC_VER
using namespace std; // Porting
#endif
#ifdef __MWERKS__
using namespace std; // Porting
#endif
#ifdef __GNUC__
//#define std
using namespace std; // Porting
#endif

//-------------------------------------------------------------------------

struct radTInputCell {

	char Type;
	long AuxNum;

	radTInputCell() { AuxNum = 0;}
	radTInputCell(char InType, long InAuxNum) { Type = InType; AuxNum = InAuxNum;}
};

//-------------------------------------------------------------------------

struct radTPairOfDouble {
	double First, Second;
	radTPairOfDouble(double InFirst =0., double InSecond =0.) { First = InFirst; Second = InSecond;}
};

//-------------------------------------------------------------------------

struct radTPairOfVect3d {
	TVector3d V1, V2;
	radTPairOfVect3d(TVector3d& In_V1, TVector3d& In_V2) { V1 = In_V1; V2 = In_V2;}
	radTPairOfVect3d(TVector3d& In_V1) { V1 = In_V1; V2.x = V2.y = V2.z = 0.;}
	radTPairOfVect3d() {}
};

//-------------------------------------------------------------------------

struct radTPairIntDouble {
	int mInt;
	double mDouble;

	radTPairIntDouble(int inInt, double inDouble)
	{
		mInt = inInt; mDouble = inDouble;
	}

	static bool less(const radTPairIntDouble& P1, const radTPairIntDouble& P2)
	{
		return (P1.mDouble < P2.mDouble);
	}

	//friend int operator <(const radTAuxIndNorm&, const radTAuxIndNorm&);
	//friend bool operator >(const radTAuxIndNorm&, const radTAuxIndNorm&);
	//friend int operator ==(const radTAuxIndNorm&, const radTAuxIndNorm&);
	//friend bool operator !=(const radTAuxIndNorm&, const radTAuxIndNorm&);
};

//-------------------------------------------------------------------------

using radTVectInputCell = vector<radTInputCell>;
using radTVectorOfVector3d = vector<TVector3d>;
using radTVectPairOfDouble = vector<radTPairOfDouble>;
using radTVectPairOfVect3d = vector<radTPairOfVect3d>;

#ifdef __MWERKS__
/*
null_template
struct iterator_traits <radTPairOfVect3d*> {
	 typedef ptrdiff_t difference_type;
	 typedef radTPairOfVect3d value_type;
	 typedef radTPairOfVect3d* pointer;
	 typedef radTPairOfVect3d& reference;
	 typedef random_access_iterator_tag iterator_category;
};
null_template
struct iterator_traits <radTPairOfDouble*> {
	 typedef ptrdiff_t difference_type;
	 typedef radTPairOfDouble value_type;
	 typedef radTPairOfDouble* pointer;
	 typedef radTPairOfDouble& reference;
	 typedef random_access_iterator_tag iterator_category;
};
*/
#endif

//-------------------------------------------------------------------------

struct radTGeomPolygon {
	std::vector<double> vVertCoords; // Vector for ownership
	double* VertCoords;
	int Nv;
	float ColRGB[3];

	radTGeomPolygon()
	{
		VertCoords = 0;
		Nv = 0;
		ColRGB[0] = ColRGB[1] = ColRGB[2] = -1;
	}

	// Destructor to ensure cleanup (though vector handles it automatically)
	~radTGeomPolygon()
	{
		// RAII: vector cleanup is automatic
	}
};

using radTVectGeomPolygon = vector<radTGeomPolygon>;

//-------------------------------------------------------------------------

#endif
