#pragma once

#include "mex.h"

#include <cstddef>
#include <string>

bool DispatchReactorCommand(const std::string& command, int nlhs,
                            mxArray* plhs[], int nrhs,
                            const mxArray* prhs[]);
void CleanupReactorHandles();
std::size_t ReactorHandleCount();
