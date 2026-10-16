/* SPDX-License-Identifier: MPL-2.0 */
#include "ams_fem_host.hpp"
#include <hcurlhofespace.hpp>
#include <hcurlfe.hpp>
#include <meshaccess.hpp>
#include <algorithm>
#include <climits>
#include <cmath>
#include <chrono>
#include <set>
#include <stdexcept>

RadiaHCurlStage RadiaStageHCurl(const std::string& path,
    const std::vector<std::string>& labels,const std::string& dirichlet) {
    auto phaseStart=std::chrono::steady_clock::now();
    auto elapsed=[&](){auto now=std::chrono::steady_clock::now();
        double seconds=std::chrono::duration<double>(now-phaseStart).count();
        phaseStart=now;return seconds;};
    auto mesh=std::make_shared<ngcomp::MeshAccess>(path);
    auto ne=mesh->GetNE(ngfem::VOL);
    if(mesh->GetDimension()!=3||!ne||ne>size_t((INT_MAX-256)/36)||mesh->GetDeformation()||mesh->GetCurveOrder()>1)
        throw std::invalid_argument("HCurl GPU assembly requires undeformed affine 3D tetrahedra");
    std::set<std::string> expected(labels.begin(),labels.end()),actual;
    for(auto name:mesh->GetMaterials(ngfem::BND))actual.emplace(name);
    if(expected.empty()||expected.size()!=labels.size()||expected!=actual)
        throw std::invalid_argument("HCurl .vol boundary label contract mismatch");
    double meshSeconds=elapsed();
    ngcore::Flags flags;flags.SetFlag("order",1);flags.SetFlag("nograds");
    flags.SetFlag("dirichlet",dirichlet);
    auto space=std::make_shared<ngcomp::HCurlHighOrderFESpace>(mesh,flags);
    space->Update();space->FinalizeUpdate();
    auto nd=space->GetNDof();
    if(!nd||nd>INT_MAX-256||nd!=mesh->GetNEdges())
        throw std::overflow_error("HCurl assembly requires one int32 DOF per edge");
    RadiaHCurlStage s;s.mesh_seconds=meshSeconds;s.space_seconds=elapsed();
    s.n=int(nd);s.e=int(ne);
    ngfem::IntegrationRule rule(ngfem::ET_TET,2);s.nq=int(rule.Size());
    // Each numeric variable remains below 2 GiB; larger workflows require an
    // explicit HDF5/out-of-core design, not overflow or silent splitting.
    if(size_t(s.e)*s.nq*18>size_t(2ULL<<30)/sizeof(double))
        throw std::overflow_error("Mapped HCurl shape array exceeds 2 GiB capacity");
    s.dofs.resize(ne*6);s.slots.resize(ne*36);s.curls.resize(ne*18);
    s.shape.resize(ne*s.nq*18);s.weights.resize(ne*s.nq);s.volume.resize(ne);
    s.free.resize(nd);auto mask=space->GetFreeDofs();
    for(int i=0;i<s.n;++i)s.free[i]=mask->Test(i);
    std::vector<std::vector<int>> rows(nd),incident(nd);
    ngstd::LocalHeap heap(1<<22,"radia_cuda_hcurl_stage");
    s.allocation_seconds=elapsed();
    for(int e=0;e<s.e;++e){
        ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::VOL,e);
        if(mesh->GetElType(eid)!=ngfem::ET_TET||mesh->GetElement(eid).is_curved)
            throw std::invalid_argument("HCurl GPU assembly rejects curved or non-tetrahedral elements");
        auto& trafo=mesh->GetTrafo(eid,heap);
        auto& fe=dynamic_cast<const ngfem::HCurlFiniteElement<3>&>(space->GetFE(eid,heap));
        ngstd::Array<int> dofs;space->GetDofNrs(eid,dofs);
        if(dofs.Size()!=6||fe.GetNDof()!=6)throw std::invalid_argument("Unsupported HCurl element order");
        for(int a=0;a<6;++a){
            if(dofs[a]<0||dofs[a]>=s.n)throw std::invalid_argument("Unsupported HCurl constrained DOF");
            s.dofs[e*6+a]=dofs[a];incident[dofs[a]].push_back(e*6+a);
            for(int b=0;b<6;++b)rows[dofs[a]].push_back(dofs[b]);
        }
        ngbla::Matrix<> shape(6,3),curl(6,3);
        ngfem::IntegrationPoint center(.25,.25,.25);
        ngfem::MappedIntegrationPoint<3,3> cmip(center,trafo);fe.CalcMappedCurlShape(cmip,curl);
        for(int a=0;a<6;++a)for(int c=0;c<3;++c){double v=curl(a,c);
            if(!std::isfinite(v))throw std::runtime_error("Nonfinite NGSolve mapped curl");s.curls[e*18+a*3+c]=v;}
        double volume=0;
        for(int q=0;q<s.nq;++q){ngfem::MappedIntegrationPoint<3,3> mip(rule[q],trafo);
            fe.CalcMappedShape(mip,shape);double w=mip.GetWeight();
            if(!std::isfinite(w)||w<=0)throw std::runtime_error("Invalid NGSolve mapped measure");
            s.weights[e*s.nq+q]=w;volume+=w;
            for(int a=0;a<6;++a)for(int c=0;c<3;++c){double v=shape(a,c);
                if(!std::isfinite(v))throw std::runtime_error("Nonfinite NGSolve mapped shape");
                s.shape[(size_t(e)*s.nq+q)*18+a*3+c]=v;}
        }s.volume[e]=volume;
    }
    s.element_seconds=elapsed();
    s.rp.push_back(0);s.incident_rp.push_back(0);
    for(int i=0;i<s.n;++i){auto& row=rows[i];std::sort(row.begin(),row.end());row.erase(std::unique(row.begin(),row.end()),row.end());
        if(row.empty())throw std::runtime_error("Unused HCurl DOF");
        if(s.ci.size()+row.size()>INT_MAX-256)throw std::overflow_error("HCurl CSR exceeds int32 capacity");
        s.ci.insert(s.ci.end(),row.begin(),row.end());s.rp.push_back(int(s.ci.size()));
        s.incident.insert(s.incident.end(),incident[i].begin(),incident[i].end());s.incident_rp.push_back(int(s.incident.size()));}
    for(int e=0;e<s.e;++e)for(int a=0;a<6;++a)for(int b=0;b<6;++b){int i=s.dofs[e*6+a],j=s.dofs[e*6+b];
        auto begin=s.ci.begin()+s.rp[i],end=s.ci.begin()+s.rp[i+1];auto p=std::lower_bound(begin,end,j);
        if(p==end||*p!=j)throw std::runtime_error("HCurl CSR scatter contract failed");s.slots[e*36+a*6+b]=int(p-s.ci.begin());}
    s.csr_seconds=elapsed();
    auto grad=space->CreateGradient();
    if(grad->Height()!=s.n||grad->Width()!=mesh->GetNV())throw std::runtime_error("NGSolve gradient/vertex contract failed");
    s.grp.push_back(0);
    for(int i=0;i<s.n;++i){auto cols=grad->GetRowIndices(i);auto values=grad->GetRowValues(i);
        for(int k=0;k<cols.Size();++k){s.gci.push_back(cols[k]);s.gvalues.push_back(values[k]);}s.grp.push_back(int(s.gci.size()));}
    s.gradient_seconds=elapsed();
    const int nv=int(mesh->GetNV());s.coordinates.resize(size_t(nv)*3);
    for(int i=0;i<nv;++i){auto p=mesh->GetPoint<3>(i);for(int c=0;c<3;++c)s.coordinates[i+size_t(nv)*c]=p[c];}
    s.coordinate_seconds=elapsed();
    return s;
}
