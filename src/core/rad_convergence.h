/*-------------------------------------------------------------------------
*
* File name:      radcnvrg.h
*
* Project:        RADIA
*
* Description:    Small perturbation for numerical stability
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
*
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
* Modified:       2024 - Changed from random to deterministic perturbation
*                        for reproducible results
*
-------------------------------------------------------------------------*/

#ifndef __RADCNVRG_H
#define __RADCNVRG_H

#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "gmvect.h"

//-------------------------------------------------------------------------

class radTConvergRepair {
public:

	double AbsRand;
	double RelRand;
	double ZeroRand;

	short ActOnDoubles;

	radTConvergRepair()
	{
		ActOnDoubles = 1;

		AbsRand = 1.E-09;
		RelRand = 1.E-11;
		ZeroRand = AbsRand;
	}

	void SwitchActOnDoubles(short InActOnDoubles, double InAbsRand =0., double InRelRand =0., double InZeroRand =0.)
	{
		ActOnDoubles = InActOnDoubles;

		AbsRand = InAbsRand;
		RelRand = InRelRand;
		ZeroRand = InZeroRand;
	}

	// Deterministic perturbation functions
	// These add small fixed offsets to avoid numerical singularities
	// when observation points lie exactly on element boundaries.
	// The perturbation is ~1e-9 mm which is physically negligible.

	double Double(double A)
	{
		if(!ActOnDoubles) return A;
		// Use fixed offset instead of random for determinism
		// The 0.31830988618 factor (1/pi) provides an irrational offset
		// to avoid exact alignment with typical grid coordinates
		A += AbsRand * 0.31830988618;
		A *= (1. + RelRand * 0.27182818285);  // e/10
		if(A == 0.) A = ZeroRand * 0.31830988618;
		return A;
	}
	double DoublePlus(double A)
	{
		if(!ActOnDoubles) return A;
		A += AbsRand * 0.31830988618;
		A *= (1. + RelRand * 0.27182818285);
		if(A == 0.) A = ZeroRand * 0.31830988618;
		return A;
	}
	double DoubleMinus(double A)
	{
		if(!ActOnDoubles) return A;
		A -= AbsRand * 0.31830988618;
		A *= (1. - RelRand * 0.27182818285);
		if(A == 0.) A = -ZeroRand * 0.31830988618;
		return A;
	}

	double AbsRandMagnitude(double A)
	{
		if(!ActOnDoubles) return 0.;
		double AbsFromRel = RelRand*A;
		double AbsMax = (AbsFromRel<AbsRand)? AbsRand : AbsFromRel;
		return (A!=0.)? AbsMax : ((AbsMax<ZeroRand)? ZeroRand : AbsMax);
	}
};

//-------------------------------------------------------------------------

#endif

