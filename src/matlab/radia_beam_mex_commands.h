#pragma once

#include "mex.h"

#include <memory>
#include <string>

namespace ngcomp {
class GridFunction;
}

namespace rad_hdiv {
class HDivFieldEvaluator;
}

bool DispatchBeamCommand(const std::string& command, int nlhs,
                         mxArray* plhs[], int nrhs,
                         const mxArray* prhs[]);

void BeamTransferFromGridFunction(
    std::shared_ptr<ngcomp::GridFunction> field, int nlhs,
    mxArray* plhs[], int nrhs, const mxArray* prhs[]);

void BeamTrackGridFunction(
    std::shared_ptr<ngcomp::GridFunction> field, int nlhs,
    mxArray* plhs[], int nrhs, const mxArray* prhs[]);

void BeamTrackReferenceOrbitToPlane(
    std::shared_ptr<rad_hdiv::HDivFieldEvaluator> field, int nlhs,
    mxArray* plhs[], int nrhs, const mxArray* prhs[]);
