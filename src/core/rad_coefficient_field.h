/*-------------------------------------------------------------------------
*
* File name:      rad_coefficient_field.h
*
* Project:        RADIA
*
* Description:    CoefficientFunction-based background field source
*                 (Python callback support)
*
* First release:  2025
*
*-------------------------------------------------------------------------*/

#ifndef __RAD_COEFFICIENT_FIELD_H
#define __RAD_COEFFICIENT_FIELD_H

#include "rad_geometry_3d.h"

// Forward declaration to avoid Python.h dependency in header
struct _object;
using PyObject = _object;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTCoefficientFunctionFieldSource : public radTg3d {
public:
	PyObject* cf_callback;  // Python callable (CF or function)

	radTCoefficientFunctionFieldSource(PyObject* callback);

	radTCoefficientFunctionFieldSource(const radTCoefficientFunctionFieldSource& src);

	virtual ~radTCoefficientFunctionFieldSource();

	void B_comp(radTField* FieldPtr);

	void B_intComp(radTField* FieldPtr);


	// Dump / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)

	int DuplicateItself(radThg& hg, radTApplication*, char);

	int SizeOfThis() { return sizeof(radTCoefficientFunctionFieldSource); }
};

//-------------------------------------------------------------------------

#endif
