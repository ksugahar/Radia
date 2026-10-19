#pragma once
#include "mex.h"
#include <cstdio>
#include <exception>

namespace radia_mex {
using Dispatcher=void(*)(int,mxArray**,int,const mxArray**);
// This call frame must finish exception destruction before mexFunction asks
// MATLAB to raise an error. Do not inline it back into the error gateway.
#ifdef _MSC_VER
__declspec(noinline)
#else
__attribute__((noinline))
#endif
inline bool checked(Dispatcher dispatch,int nl,mxArray** out,int nr,const mxArray** in,
                    char* message,std::size_t capacity) noexcept {
    try{dispatch(nl,out,nr,in);return true;}
    catch(const std::exception& error){std::snprintf(message,capacity,"%s",error.what());}
    catch(...){std::snprintf(message,capacity,"Unknown native exception");}
    return false;
}
}
