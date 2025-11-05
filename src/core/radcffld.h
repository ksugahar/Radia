/*-------------------------------------------------------------------------
*
* File name:      radcffld.h
*
* Project:        RADIA
*
* Description:    CoefficientFunction-based background field source
*                 (Python callback support)
*
* First release:  2025
*
*-------------------------------------------------------------------------*/

#ifndef __RADCFFLD_H
#define __RADCFFLD_H

#include "radg3d.h"

// Forward declaration to avoid Python.h dependency in header
struct _object;
using PyObject = _object;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTCoefficientFunctionFieldSource : public radTg3d {
public:
	PyObject* cf_callback;  // Python callable (CF or function)

	radTCoefficientFunctionFieldSource(PyObject* callback);

	radTCoefficientFunctionFieldSource(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers);

	radTCoefficientFunctionFieldSource(const radTCoefficientFunctionFieldSource& src);

	virtual ~radTCoefficientFunctionFieldSource();

	void B_comp(radTField* FieldPtr);

	void B_intComp(radTField* FieldPtr);

	radTg3dGraphPresent* CreateGraphPresent();

	void Dump(std::ostream& o, int ShortSign = 0);

	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut,
	             map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers,
	             int& gUniqueMapKey, int elemKey);

	int DuplicateItself(radThg& hg, radTApplication*, char);

	int SizeOfThis() { return sizeof(radTCoefficientFunctionFieldSource); }
};

//-------------------------------------------------------------------------

#endif
