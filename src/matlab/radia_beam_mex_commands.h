#pragma once

#include "mex.h"

#include <memory>
#include <string>

namespace ngcomp {
class GridFunction;
}

bool DispatchBeamCommand(const std::string& command, int nlhs,
                         mxArray* plhs[], int nrhs,
                         const mxArray* prhs[]);

void BeamTransferFromGridFunction(
    std::shared_ptr<ngcomp::GridFunction> field, int nlhs,
    mxArray* plhs[], int nrhs, const mxArray* prhs[]);
