/*-------------------------------------------------------------------------
*
* File name:      radcast.cpp
*
* Project:        RADIA
*
* Description:    "Dynamic cast" for RADIA classes
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#include "rad_type_cast.h"

#ifndef __RADINTRC_H
#include "rad_interaction.h"
#endif

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTInteraction* radTCast::InteractCast(radTg* gPtr)
{
	radTInteraction Interact;
	if(gPtr->Type_g()==Interact.Type_g()) return static_cast<radTInteraction*>(gPtr);
	else return nullptr;
	// Do not move it to the .h file: compilation problems may appear!
}
