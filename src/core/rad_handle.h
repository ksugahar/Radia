/*-------------------------------------------------------------------------
*
* File name:      radhandl.h
*
* Project:        RADIA
*
* Description:    Smart pointer for RADIA objects
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#ifndef __RADHANDLE_H
#define __RADHANDLE_H

#include <atomic>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

template<class T> class radTHandle {
public:
	T* rep;
	// ATOMIC reference count (2026-07-11): distinct radTHandle instances
	// sharing one count are copied/destroyed concurrently -- e.g. every
	// batch field evaluation copies the container handle out of
	// GlobalMapOfHandlers, and NGSolve assembly evaluates the RadiaField
	// CoefficientFunction from parallel TaskManager workers.  A plain int
	// count loses increments under that race, hits zero early, and deletes
	// live objects mid-evaluation (0xC0000374 heap corruption).
	std::atomic<int>* pcount;

	radTHandle () { rep=0; pcount=0;}
	radTHandle (T* pp) : rep(pp), pcount(new std::atomic<int>(1)) {}
	radTHandle (const radTHandle& r) : rep(r.rep), pcount(r.pcount)
	{
		if(pcount != 0) pcount->fetch_add(1, std::memory_order_relaxed);
	}

	void destroy()
	{
		if(pcount!=0)
			if(pcount->fetch_sub(1, std::memory_order_acq_rel)==1)
			{
				delete rep;
				delete pcount;
				rep=0; pcount=0;
			}
	}

	T* operator->() { return rep;}
	T* obj() { return rep;}
	void bind(const radTHandle& r)
	{
		if(rep!=r.rep)
		{
			if(r.rep!=0)
			{
				// Acquire the new reference BEFORE releasing the old one
				// (shared_ptr assignment order; rep!=r.rep excludes self-bind).
				r.pcount->fetch_add(1, std::memory_order_relaxed);
				destroy();
				rep = r.rep;
				pcount = r.pcount;
			}
			else
			{
				rep = 0;
				pcount = 0;
			}
		}
	}

	radTHandle& operator=(const radTHandle& r) 
	{ 
		bind(r); return *this;
	}

	// NOTE: operator< and operator== are provided as free function
	// templates below (lines ~91 / ~98). Member versions used to be
	// defined here as well, but C++20 rewriting candidates (P1185)
	// make MSVC report C2666 ambiguity between the member and the
	// free template when both exist. Removed 2026-04-12.

	~radTHandle()
	{
		destroy();
	}
};

//-------------------------------------------------------------------------

template<class T> inline int operator <(const radTHandle<T>& h1, const radTHandle<T>& h2)
{
	return (h1.rep < h2.rep); 
}

//-------------------------------------------------------------------------

template<class T> inline int operator ==(const radTHandle<T>& h1, const radTHandle<T>& h2)
{
	return (h1.rep == h2.rep); 
}

//-------------------------------------------------------------------------

#endif


