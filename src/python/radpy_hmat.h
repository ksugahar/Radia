/*-------------------------------------------------------------------------
* H-Matrix field evaluation Python wrappers
-------------------------------------------------------------------------*/

#include "radentry_hmat.h"

//-------------------------------------------------------------------------
// FldBatch - Batch field evaluation with optional H-matrix
//-------------------------------------------------------------------------
static PyObject* radia_FldBatch(PyObject* self, PyObject* args)
{
	PyObject *oCmpnId=0, *oP=0, *oResB=0;
	double *arCoord=0, *arB=0;
	int use_hmatrix = 0;  // Default: direct calculation

	try
	{
		int ind=0;
		if(!PyArg_ParseTuple(args, "iOO|i:FldBatch", &ind, &oCmpnId, &oP, &use_hmatrix))
			throw CombErStr(strEr_BadFuncArg, ": FldBatch");
		if((ind == 0) || (oCmpnId == 0) || (oP == 0))
			throw CombErStr(strEr_BadFuncArg, ": FldBatch");

		char sCmpnId[256];
		CPyParse::CopyPyStringToC(oCmpnId, sCmpnId, 256);

		int nCoord = 0;
		char resP = CPyParse::CopyPyNestedListElemsToNumAr(oP, 'd', arCoord, nCoord);
		if(resP == 0) throw CombErStr(strEr_BadFuncArg, ": FldBatch: array / list of points");

		int nP = (int)round(nCoord/3.);
		if(nP*3 != nCoord) throw CombErStr(strEr_BadFuncArg, ": FldBatch: array / list of points");

		const int maxNumFldCmpn = 14;
		arB = new double[maxNumFldCmpn*nP];
		int nB = 0;

		// Call H-matrix batch evaluation
		g_pyParse.ProcRes(RadFldBatch(arB, &nB, ind, sCmpnId, arCoord, nP, use_hmatrix));

		if(nB <= 0) { nB = 1; *arB = 0;}

		if(nB == 1) oResB = Py_BuildValue("d", *arB);
		else if(nB > 1) oResB = CPyParse::SetDataListOfLists(arB, nB, nP);
		if(oResB) Py_XINCREF(oResB);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}
	if(arCoord != 0) delete[] arCoord;
	if(arB != 0) delete[] arB;
	return oResB;
}

//-------------------------------------------------------------------------
// SetHMatrixFieldEval - Enable/disable H-matrix acceleration
//-------------------------------------------------------------------------
static PyObject* radia_SetHMatrixFieldEval(PyObject* self, PyObject* args)
{
	try
	{
		int enabled = 0;
		double eps = 1e-6;

		if(!PyArg_ParseTuple(args, "i|d:SetHMatrixFieldEval", &enabled, &eps))
			throw CombErStr(strEr_BadFuncArg, ": SetHMatrixFieldEval");

		g_pyParse.ProcRes(RadSetHMatrixFieldEval(enabled, eps));

		Py_RETURN_NONE;
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		return NULL;
	}
}

//-------------------------------------------------------------------------
// ClearHMatrixCache - Clear H-matrix cache
//-------------------------------------------------------------------------
static PyObject* radia_ClearHMatrixCache(PyObject* self, PyObject* args)
{
	try
	{
		g_pyParse.ProcRes(RadClearHMatrixCache());
		Py_RETURN_NONE;
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		return NULL;
	}
}

//-------------------------------------------------------------------------
// GetHMatrixStats - Get H-matrix statistics
//-------------------------------------------------------------------------
static PyObject* radia_GetHMatrixStats(PyObject* self, PyObject* args)
{
	PyObject *oRes=0;
	double *arStats=0;

	try
	{
		arStats = new double[10];  // Reserve space for stats
		int nStats = 0;

		g_pyParse.ProcRes(RadGetHMatrixStats(arStats, &nStats));

		if(nStats > 0) {
			oRes = PyList_New(nStats);
			for(int i=0; i<nStats; i++) {
				PyList_SetItem(oRes, i, Py_BuildValue("d", arStats[i]));
			}
		} else {
			oRes = PyList_New(0);
		}

		if(oRes) Py_XINCREF(oRes);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}
	if(arStats != 0) delete[] arStats;
	return oRes;
}

//-------------------------------------------------------------------------
// UpdateHMatrixMagnetization - Fast magnetization update
//-------------------------------------------------------------------------
static PyObject* radia_UpdateHMatrixMagnetization(PyObject* self, PyObject* args)
{
	try
	{
		int obj = 0;

		if(!PyArg_ParseTuple(args, "i:UpdateHMatrixMagnetization", &obj))
			throw CombErStr(strEr_BadFuncArg, ": UpdateHMatrixMagnetization");

		int result = RadUpdateHMatrixMagnetization(obj);

		if(result == -2) {
			throw CombErStr("H-matrix not built yet. Call FldBatch with use_hmatrix=1 first.", "");
		}
		else if(result != 0) {
			throw CombErStr("Failed to update magnetization", "");
		}

		Py_RETURN_NONE;
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		return NULL;
	}
}

//-------------------------------------------------------------------------
