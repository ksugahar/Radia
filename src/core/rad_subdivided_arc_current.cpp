/*-------------------------------------------------------------------------
*
* File name:      radsbdac.h
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 subdivided rectangular cross-section arc with azimuthal current
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#include "rad_subdivided_arc_current.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------

void radTSubdividedArcCur::Dump(std::ostream& o, int ShortSign)
{
	((radTg3d*)((radTGroup*)this))->radTg3d::Dump(o, ShortSign);

	o << "Subdivided ArcCur";
	
	if(ShortSign==1) return;

	o << endl;
	radTArcCur::DumpPureObjInfo(o, ShortSign);

	o << endl;
	radTGroup::DumpPureObjInfo(o, ShortSign);

	o << endl;
	((radTg3d*)((radTGroup*)this))->radTg3d::DumpTransApplied(o);

	o << endl;
	o << "   Memory occupied (incl. the content): " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------
