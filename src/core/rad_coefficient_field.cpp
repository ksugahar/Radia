/*-------------------------------------------------------------------------
*
* File name:      radcffld.cpp
*
* Project:        RADIA
*
* Description:    CoefficientFunction-based background field source
*                 (Python callback support)
*
* First release:  2025
*
*-------------------------------------------------------------------------*/

#include "rad_coefficient_field.h"
#include "rad_exception.h"
#include <Python.h>
#include <iostream>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTCoefficientFunctionFieldSource::radTCoefficientFunctionFieldSource(PyObject* callback)
	: cf_callback(callback)
{
	if(cf_callback) {
		Py_INCREF(cf_callback);
	}
}

//-------------------------------------------------------------------------

radTCoefficientFunctionFieldSource::radTCoefficientFunctionFieldSource(
	CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
{
	// Binary deserialization not supported for Python callbacks
	DumpBinParse_g3d(inStr, mKeysOldNew, gMapOfHandlers);
	cf_callback = 0;
}

//-------------------------------------------------------------------------

radTCoefficientFunctionFieldSource::radTCoefficientFunctionFieldSource(
	const radTCoefficientFunctionFieldSource& src)
	: radTg3d(src), cf_callback(src.cf_callback)
{
	if(cf_callback) {
		Py_INCREF(cf_callback);
	}
}

//-------------------------------------------------------------------------

radTCoefficientFunctionFieldSource::~radTCoefficientFunctionFieldSource()
{
	if(cf_callback) {
		// Check if Python is still initialized before trying to acquire GIL
		// During interpreter shutdown, Py_IsInitialized() returns 0
		if(Py_IsInitialized()) {
			PyGILState_STATE gstate = PyGILState_Ensure();
			Py_DECREF(cf_callback);
			PyGILState_Release(gstate);
		}
		// If Python is shutting down, let it handle cleanup
		// Attempting to DECREF during shutdown can cause crashes
	}
}

//-------------------------------------------------------------------------

void radTCoefficientFunctionFieldSource::B_comp(radTField* FieldPtr)
{
	if(!cf_callback) return;

	// Acquire Python GIL for callback
	PyGILState_STATE gstate = PyGILState_Ensure();

	try {
		// Prepare coordinates as Python list [x, y, z] in mm
		PyObject* coords = PyList_New(3);
		PyList_SetItem(coords, 0, PyFloat_FromDouble(FieldPtr->P.x));
		PyList_SetItem(coords, 1, PyFloat_FromDouble(FieldPtr->P.y));
		PyList_SetItem(coords, 2, PyFloat_FromDouble(FieldPtr->P.z));

		// Call Python callback: result = callback([x, y, z])
		// Result can be [Bx, By, Bz] or {'B': [Bx, By, Bz], 'A': [Ax, Ay, Az]}
		PyObject* result = PyObject_CallFunctionObjArgs(cf_callback, coords, nullptr);
		Py_DECREF(coords);

		if(!result) {
			// Python exception occurred
			PyErr_Print();
			PyGILState_Release(gstate);
			return;
		}

		// Check if result is a dictionary (for B and A)
		if(PyDict_Check(result)) {
			// Dictionary format: {'B': [Bx, By, Bz], 'A': [Ax, Ay, Az]}
			PyObject* B_obj = PyDict_GetItemString(result, "B");
			PyObject* A_obj = PyDict_GetItemString(result, "A");

			// Extract B field
			if(B_obj && (PyList_Check(B_obj) || PyTuple_Check(B_obj))) {
				if(PySequence_Size(B_obj) == 3) {
					PyObject* bx = PySequence_GetItem(B_obj, 0);
					PyObject* by = PySequence_GetItem(B_obj, 1);
					PyObject* bz = PySequence_GetItem(B_obj, 2);

					double Bx = PyFloat_AsDouble(bx);
					double By = PyFloat_AsDouble(by);
					double Bz = PyFloat_AsDouble(bz);

					Py_DECREF(bx); Py_DECREF(by); Py_DECREF(bz);

					if(!PyErr_Occurred()) {
						// B_from_cf is magnetic flux density in Tesla
						TVector3d B_from_cf(Bx, By, Bz);
						if(FieldPtr->FieldKey.B_) FieldPtr->B += B_from_cf;
						// Convert B to H using H = B/μ₀ where μ₀ = 1.25663706212e-6 T/(A/m)
						if(FieldPtr->FieldKey.H_) FieldPtr->H += B_from_cf * 795774.715459;  // 1/μ₀
					}
				}
			}

			// Extract A field (vector potential)
			if(A_obj && (PyList_Check(A_obj) || PyTuple_Check(A_obj))) {
				if(PySequence_Size(A_obj) == 3) {
					PyObject* ax = PySequence_GetItem(A_obj, 0);
					PyObject* ay = PySequence_GetItem(A_obj, 1);
					PyObject* az = PySequence_GetItem(A_obj, 2);

					double Ax = PyFloat_AsDouble(ax);
					double Ay = PyFloat_AsDouble(ay);
					double Az = PyFloat_AsDouble(az);

					Py_DECREF(ax); Py_DECREF(ay); Py_DECREF(az);

					if(!PyErr_Occurred()) {
						TVector3d A_from_cf(Ax, Ay, Az);
						if(FieldPtr->FieldKey.A_) FieldPtr->A += A_from_cf;
					}
				}
			}

			Py_DECREF(result);
		}
		// Backward compatibility: [Bx, By, Bz] format
		else if(PyList_Check(result) || PyTuple_Check(result)) {
			Py_ssize_t size = PySequence_Size(result);
			if(size != 3) {
				std::cerr << "[CFFieldSource] Callback must return 3 components [Bx, By, Bz], got "
				          << size << std::endl;
				Py_DECREF(result);
				PyGILState_Release(gstate);
				return;
			}

			PyObject* item0 = PySequence_GetItem(result, 0);
			PyObject* item1 = PySequence_GetItem(result, 1);
			PyObject* item2 = PySequence_GetItem(result, 2);

			double Bx = PyFloat_AsDouble(item0);
			double By = PyFloat_AsDouble(item1);
			double Bz = PyFloat_AsDouble(item2);

			Py_DECREF(item0);
			Py_DECREF(item1);
			Py_DECREF(item2);
			Py_DECREF(result);

			if(PyErr_Occurred()) {
				PyErr_Print();
				PyGILState_Release(gstate);
				return;
			}

			// B_from_cf is magnetic flux density in Tesla
			TVector3d B_from_cf(Bx, By, Bz);
			if(FieldPtr->FieldKey.B_) FieldPtr->B += B_from_cf;
			// Convert B to H using H = B/μ₀ where μ₀ = 1.25663706212e-6 T/(A/m)
			if(FieldPtr->FieldKey.H_) FieldPtr->H += B_from_cf * 795774.715459;  // 1/μ₀

			// For backward compatibility, A is not provided in list format
		}
		else {
			std::cerr << "[CFFieldSource] Callback must return [Bx, By, Bz] or {'B': [...], 'A': [...]}" << std::endl;
			Py_DECREF(result);
			PyGILState_Release(gstate);
			return;
		}

		if(PyErr_Occurred()) {
			PyErr_Print();
		}

	} catch (...) {
		std::cerr << "[CFFieldSource] Exception in B_comp" << std::endl;
	}

	PyGILState_Release(gstate);
}

//-------------------------------------------------------------------------

void radTCoefficientFunctionFieldSource::B_intComp(radTField* FieldPtr)
{
	if(!cf_callback) return;

	if(FieldPtr->FieldKey.FinInt_)
	{
		// For arbitrary field, use numerical integration
		// Simple trapezoidal rule: Integral ~= (B(P1) + B(P2))/2 * L

		TVector3d P1 = FieldPtr->P;
		TVector3d P2 = FieldPtr->NextP;
		TVector3d D = P2 - P1;
		double L = sqrt(D.x*D.x + D.y*D.y + D.z*D.z);

		// Evaluate field at both endpoints
		radTField F1 = *FieldPtr;
		F1.P = P1;
		F1.B.Zero();
		F1.H.Zero();
		F1.FieldKey.FinInt_ = 0;  // Disable infinite integral for evaluation
		B_comp(&F1);

		radTField F2 = *FieldPtr;
		F2.P = P2;
		F2.B.Zero();
		F2.H.Zero();
		F2.FieldKey.FinInt_ = 0;
		B_comp(&F2);

		// Average and multiply by length
		TVector3d B_avg = (F1.B + F2.B) * 0.5;
		TVector3d H_avg = (F1.H + F2.H) * 0.5;

		TVector3d BufIb = L * B_avg;
		TVector3d BufIh = L * H_avg;

		if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufIb;
		if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufIh;
	}

	// Infinite integral: set to zero (formally infinite for non-localized fields)
}

//-------------------------------------------------------------------------

radTg3dGraphPresent* radTCoefficientFunctionFieldSource::CreateGraphPresent()
{
	// No geometry to display for field source
	return 0;
}

//-------------------------------------------------------------------------

void radTCoefficientFunctionFieldSource::Dump(std::ostream& o, int ShortSign)
{
	radTg3d::Dump(o);
	o << "CoefficientFunction-based background field source";
	if(ShortSign==1) return;
	o << endl;
	o << "   Python callback: " << (cf_callback ? "registered" : "none");
	o << endl;
	o << "   Memory occupied: " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------

void radTCoefficientFunctionFieldSource::DumpBin(
	CAuxBinStrVect& oStr, vector<int>& vElemKeysOut,
	map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers,
	int& gUniqueMapKey, int elemKey)
{
	// Binary serialization not supported for Python callbacks
	// Silently do nothing - user won't be able to save/load CF field sources
}

//-------------------------------------------------------------------------

int radTCoefficientFunctionFieldSource::DuplicateItself(
	radThg& hg, radTApplication*, char)
{
	return FinishDuplication(new radTCoefficientFunctionFieldSource(*this), hg);
}

//-------------------------------------------------------------------------
