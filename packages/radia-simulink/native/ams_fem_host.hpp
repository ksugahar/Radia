/* SPDX-License-Identifier: MPL-2.0 */
#pragma once
#include <vector>
#include <string>
struct RadiaHCurlStage {
    int n=0,e=0,nq=0;
    std::vector<int> rp,ci,dofs,slots,incident_rp,incident,free,grp,gci;
    std::vector<double> shape,curls,weights,volume,coordinates,gvalues;
};
// Standard NGSolve geometry/basis APIs; no TaskManager regions or Python calls.
RadiaHCurlStage RadiaStageHCurl(const std::string&,const std::vector<std::string>&,
                              const std::string&);
