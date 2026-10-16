/* SPDX-License-Identifier: MPL-2.0 */
#pragma once
#include <vector>
#include <string>
struct RadiaHCurlStage {
    // Disjoint serial wall phases; element includes incident-row collection.
    double mesh_seconds=0,space_seconds=0,allocation_seconds=0;
    double element_seconds=0,csr_seconds=0,gradient_seconds=0,coordinate_seconds=0;
    int n=0,e=0,nq=0;
    std::vector<int> rp,ci,dofs,slots,incident_rp,incident,free,grp,gci;
    std::vector<double> shape,curls,weights,volume,coordinates,gvalues;
};
// Standard NGSolve geometry/basis APIs; no TaskManager regions or Python calls.
RadiaHCurlStage RadiaStageHCurl(const std::string&,const std::vector<std::string>&,
                              const std::string&);
