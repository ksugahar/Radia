/*-------------------------------------------------------------------------
*
* File name:      radgeom_types.h
*
* Project:        RADIA
*
* Description:    Geometric type definitions for polygons and vertices
*
* Author(s):      Oleg Chubar
*
* First release:  1997
*
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADGEOM_TYPES_H
#define __RADGEOM_TYPES_H

#ifndef __GMVECT_H
#include "gmvect.h"
#endif

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct radTIntPtrAndInt {
	std::vector<int> vInt;  // Vector for ownership
	int* pInt;
	int AnInt;
	radTIntPtrAndInt(int* In_pInt =0, int InAnInt =0) { pInt=In_pInt; AnInt=InAnInt;}
	// New constructor that takes ownership via vector
	radTIntPtrAndInt(std::vector<int>&& v, int InAnInt) : vInt(std::move(v)), AnInt(InAnInt) {
		pInt = vInt.data();
	}
};

//-------------------------------------------------------------------------

using radTVectIntPtrAndInt = vector<radTIntPtrAndInt>;

//-------------------------------------------------------------------------

struct radTPtrsToPgnAndVect2d {
	radTPolygon* pPgn;
	TVector2d* pVect2d;
	int AmOfPoints;
	radTPtrsToPgnAndVect2d() { pPgn = 0; pVect2d = 0; AmOfPoints = 0;}
};

//-------------------------------------------------------------------------

struct radTVertexPointLiaison {
#ifdef __GNUC__
	vector<int> FirstIndVect;
	vector<int> SecondIndVect;
#else
	vector<int, allocator<int> > FirstIndVect;
	vector<int, allocator<int> > SecondIndVect;
#endif
	char AdjSegmentUsed;

	radTVertexPointLiaison() { AdjSegmentUsed = 0;}
};

//-------------------------------------------------------------------------

#endif
