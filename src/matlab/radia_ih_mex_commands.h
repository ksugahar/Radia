#pragma once

#include "mex.h"

#include <cstddef>
#include <string>

bool DispatchIHCommand(const std::string& command, int nlhs, mxArray* plhs[],
                       int nrhs, const mxArray* prhs[]);
void CleanupIHHandles();
std::size_t IHHandleCount();
