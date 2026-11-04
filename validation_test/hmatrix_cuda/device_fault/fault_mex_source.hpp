// Force-included into the unmodified radia_hdiv_device_mex.cpp of the fault build: the source keeps
// its text, its mexFunction becomes the inner entry of the shim and each mx allocation can fail.
#pragma once
#define mexFunction radia_hdiv_device_inner
#include "mex.h"
#include "fault_state.hpp"
inline mxArray* faultCreateNumericMatrix(mwSize m,mwSize n,mxClassID c,mxComplexity x){return radia_hdiv_fault::mxFails()?nullptr:mxCreateNumericMatrix(m,n,c,x);}
inline mxArray* faultCreateDoubleMatrix(mwSize m,mwSize n,mxComplexity x){return radia_hdiv_fault::mxFails()?nullptr:mxCreateDoubleMatrix(m,n,x);}
inline mxArray* faultCreateDoubleScalar(double v){return radia_hdiv_fault::mxFails()?nullptr:mxCreateDoubleScalar(v);}
inline mxArray* faultCreateString(const char* s){return radia_hdiv_fault::mxFails()?nullptr:mxCreateString(s);}
inline mxArray* faultCreateStructMatrix(mwSize m,mwSize n,int k,const char** names){return radia_hdiv_fault::mxFails()?nullptr:mxCreateStructMatrix(m,n,k,names);}
inline mxArray* faultCreateLogicalScalar(bool v){return radia_hdiv_fault::mxFails()?nullptr:mxCreateLogicalScalar(v);}
inline char* faultArrayToString(const mxArray* a){return radia_hdiv_fault::mxFails()?nullptr:mxArrayToString(a);}
inline int faultAddField(mxArray* a,const char* name){return radia_hdiv_fault::mxFails()?-1:mxAddField(a,name);}
#undef mxCreateNumericMatrix
#undef mxCreateDoubleMatrix
#undef mxCreateDoubleScalar
#undef mxCreateString
#undef mxCreateStructMatrix
#undef mxCreateLogicalScalar
#undef mxArrayToString
#undef mxAddField
#define mxCreateNumericMatrix faultCreateNumericMatrix
#define mxCreateDoubleMatrix faultCreateDoubleMatrix
#define mxCreateDoubleScalar faultCreateDoubleScalar
#define mxCreateString faultCreateString
#define mxCreateStructMatrix faultCreateStructMatrix
#define mxCreateLogicalScalar faultCreateLogicalScalar
#define mxArrayToString faultArrayToString
#define mxAddField faultAddField
