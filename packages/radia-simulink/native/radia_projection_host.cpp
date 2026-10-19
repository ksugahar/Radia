#include "mex_error_boundary.hpp"
// Experimental mexcuda host bridge: original NGSolve SetValues semantics.
// MATLAB evaluates the captured mapped targets in one batch, then replays CF.
#include "mex.h"
#include "radia_projection_host.cuh"
#include <h1hofespace.hpp>
#include <hcurlhofespace.hpp>
#include <hdivhofespace.hpp>
#include <meshaccess.hpp>
#include <array>
#include <map>
#include <chrono>
#include <random>
#include <vector>
#include <limits>
#include <cmath>
#include <cstdint>
#include "hdiv_rt0_host.cuh"
#include "hdiv_element_setup.cuh"

namespace {
using Point=std::array<double,3>;
struct Targets final:ngfem::CoefficientFunction {
    mutable std::map<Point,std::size_t> indices;
    mutable std::vector<Point> points;
    std::vector<Point> values;
    bool collecting=true;
    Targets():CoefficientFunction(3){}
    double Evaluate(const ngfem::BaseMappedIntegrationPoint&) const override {
        throw std::runtime_error("Vector projection requires vector Evaluate");
    }
    void Evaluate(const ngfem::BaseMappedIntegrationPoint& mip,
                  ngbla::FlatVector<> result) const override {
        auto p=mip.GetPoint();
        Point key{p[0],p.Size()>1?p[1]:0,p.Size()>2?p[2]:0};
        for(double v:key)if(!std::isfinite(v))throw std::runtime_error("Nonfinite mapped target");
        auto found=indices.find(key);
        if(collecting){
            if(found==indices.end()){indices.emplace(key,points.size());points.push_back(key);}
            for(int j=0;j<3;++j)result(j)=0;
        }else{
            if(found==indices.end())throw std::runtime_error("Projection targets changed; rebuild context");
            for(int j=0;j<3;++j)result(j)=values.at(found->second)[j];
        }
    }
    void Evaluate(const ngfem::BaseMappedIntegrationRule& rule,
                  ngbla::BareSliceMatrix<> result) const override {
        for(std::size_t i=0;i<rule.Size();++i){
            ngbla::Vector<> v(3);Evaluate(rule[i],v);
            for(int j=0;j<3;++j)result(i,j)=v(j);
        }
    }
};
struct Context {
    std::shared_ptr<ngcomp::MeshAccess> mesh;
    std::shared_ptr<ngcomp::FESpace> space;
    std::shared_ptr<ngcomp::GridFunction> gf;
    std::shared_ptr<Targets> targets=std::make_shared<Targets>();
};
std::shared_ptr<ngcomp::GridFunction> createHeadlessGridFunction(
    const std::shared_ptr<ngcomp::FESpace>& space){
    // NGSolve's default visualization registers a SolutionData object that
    // strongly owns MeshAccess. A standalone numerical MEX has no Netgen UI.
    ngcore::Flags flags;flags.SetFlag("novisual",true);
    return ngcomp::CreateGridFunction(space,"radia_cuda_projection",flags);
}
std::map<std::uint64_t,std::unique_ptr<Context>> contexts;
std::map<std::uint64_t,std::unique_ptr<radia::hdiv_cuda::MaterialStagingContext>> materialContexts;
#ifdef RADIA_PROJECTION_LIFECYCLE_TEST
// Only the latest released generation is observed. Weak references do not
// keep FE payloads alive or accumulate a history of released objects.
struct ReleasedOwners {
    std::uint64_t handle=0;
    std::weak_ptr<ngcomp::MeshAccess> mesh;
    std::weak_ptr<ngcomp::FESpace> space;
    std::weak_ptr<ngcomp::GridFunction> gf;
    std::weak_ptr<Targets> targets;
} releasedOwners;

// Phase probes intentionally avoid the public handle registry. They use the
// same concrete FE/CF types as create(), with no Python or TaskManager region.
std::array<double,8> probeOwnerPhase(const std::string& filename,const std::string& phase){
    int level=-1;
    bool headless=false;
    if(phase=="mesh-only")level=0;
    else if(phase=="targets-cf")level=1;
    else if(phase=="fespace-construct")level=2;
    else if(phase=="fespace-update")level=3;
    else if(phase=="gridfunction-create")level=4;
    else if(phase=="gridfunction-update")level=5;
    else if(phase=="interpolate")level=6;
    else if(phase=="gridfunction-novisual-create"){level=4;headless=true;}
    else if(phase=="gridfunction-novisual-update"){level=5;headless=true;}
    else if(phase=="interpolate-novisual"){level=6;headless=true;}
    if(level<0)throw std::invalid_argument("Unknown projection ownership phase");
    if(!contexts.empty()||!materialContexts.empty())
        throw std::runtime_error("Ownership phase probe requires no live handles");
    auto mesh=std::make_shared<ngcomp::MeshAccess>(filename);
    if(mesh->GetDimension()!=3)throw std::runtime_error("Phase probe requires 3D mesh");
    std::shared_ptr<Targets> targets;
    std::shared_ptr<ngcomp::FESpace> space;
    std::shared_ptr<ngcomp::GridFunction> gf;
    if(level>=1)targets=std::make_shared<Targets>();
    if(level>=2){
        ngcore::Flags flags;flags.SetFlag("order",1);
        space=std::make_shared<ngcomp::HDivHighOrderFESpace>(mesh,flags);
    }
    if(level>=3){space->Update();space->FinalizeUpdate();}
    if(level>=4)gf=headless ? createHeadlessGridFunction(space) :
        ngcomp::CreateGridFunction(space,"radia_cuda_projection",ngcore::Flags());
    if(level>=5)gf->Update();
    if(level>=6)radia::cuda_projection::Interpolate(targets,*gf);
    std::weak_ptr<ngcomp::MeshAccess> weakMesh=mesh;
    std::weak_ptr<ngcomp::FESpace> weakSpace=space;
    std::weak_ptr<ngcomp::GridFunction> weakGf=gf;
    std::weak_ptr<Targets> weakTargets=targets;
    std::array<double,8> result{double(mesh.use_count()),double(space.use_count()),
        double(gf.use_count()),double(targets.use_count()),0,0,0,0};
    targets.reset();gf.reset();space.reset();mesh.reset();
    result[4]=double(weakMesh.use_count());result[5]=double(weakSpace.use_count());
    result[6]=double(weakGf.use_count());result[7]=double(weakTargets.use_count());
    return result;
}
#endif
std::uint64_t initialGeneration(){
    std::random_device entropy;
    const auto nonce=(std::uint64_t(entropy())<<32)^std::uint64_t(entropy());
    const auto ticks=std::uint64_t(std::chrono::steady_clock::now().time_since_epoch().count());
    return ((nonce^ticks)&0x00007fffffffffffULL)+1;
}
// Probabilistic reload separation; the checked registry still owns liveness.
std::uint64_t generation=initialGeneration();
std::string text(const mxArray* a){
    if(!mxIsChar(a))throw std::runtime_error("Expected char array");
    char* p=mxArrayToString(a);if(!p)throw std::runtime_error("Invalid text");
    std::string s(p);mxFree(p);return s;
}
std::uint64_t token(const mxArray* a){
    if(!mxIsUint64(a)||mxGetNumberOfElements(a)!=1)throw std::runtime_error("Expected scalar uint64 handle");
    auto t=*static_cast<const std::uint64_t*>(mxGetData(a));
    if(!contexts.count(t))throw std::runtime_error("Stale or wrong projection handle");return t;
}
std::uint64_t materialToken(const mxArray* a){
    if(!mxIsUint64(a)||mxGetNumberOfElements(a)!=1)throw std::runtime_error("Expected scalar uint64 material handle");
    const auto t=*static_cast<const std::uint64_t*>(mxGetData(a));
    if(!materialContexts.count(t))throw std::runtime_error("Stale or wrong material staging handle");return t;
}
void cleanup(){contexts.clear();materialContexts.clear();}
void arity(int nr,int expected,int nl,int outputs){
    if(nr!=expected||nl!=outputs)throw std::runtime_error("Invalid command arity");
}
}
void dispatchProjection(int nl,mxArray* out[],int nr,const mxArray* in[]){
    if(nr<1)throw std::runtime_error("Command required");
    auto command=text(in[0]);
    if(command=="assemble_hdiv_elements"||command=="assemble_hdiv_elements_sparse_mass"||command=="assemble_hdiv_elements_sparse"){
        arity(nr,4,nl,1);
        if(!mxIsDouble(in[2])||mxIsSparse(in[2])||mxIsComplex(in[2])||mxGetNumberOfElements(in[2])!=1)throw std::invalid_argument("HDiv order must be scalar double");
        const double order=mxGetScalar(in[2]);if(!std::isfinite(order)||order<1||order>3||std::floor(order)!=order)throw std::invalid_argument("HDiv order must be 1..3");
        if(!mxIsDouble(in[3])||mxIsSparse(in[3])||mxIsComplex(in[3])||mxGetNumberOfElements(in[3])!=3)throw std::invalid_argument("Applied field must contain three doubles");
        std::array<double,3> applied;for(int c=0;c<3;++c){applied[c]=mxGetDoubles(in[3])[c];if(!std::isfinite(applied[c]))throw std::invalid_argument("Nonfinite applied field");}
        out[0]=radia::hdiv_cuda::inspectElements(text(in[1]),int(order),applied,command!="assemble_hdiv_elements",command=="assemble_hdiv_elements_sparse");return;
    }
    if(command=="create_hdiv_material_context"){
        arity(nr,3,nl,1);
        if(!mxIsDouble(in[2])||mxIsComplex(in[2])||mxIsSparse(in[2])||mxGetNumberOfElements(in[2])!=1)throw std::invalid_argument("Material order must be scalar double");
        const double order=mxGetScalar(in[2]);if(!std::isfinite(order)||order<1||order>2||std::floor(order)!=order)throw std::invalid_argument("Material order must be 1..2");
        if(materialContexts.size()>=2)throw std::runtime_error("Material staging context count budget exceeded");
        auto c=std::make_unique<radia::hdiv_cuda::MaterialStagingContext>(text(in[1]),int(order));
        if(generation>=0x0000ffffffffffffULL)throw std::runtime_error("Projection generation exhausted");
        const auto t=0x4d41000000000000ULL|generation++;
        out[0]=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);*static_cast<std::uint64_t*>(mxGetData(out[0]))=t;
        materialContexts.emplace(t,std::move(c));mexLock();mexAtExit(cleanup);return;
    }
    if(command=="hdiv_material_context_chunk"){
        arity(nr,4,nl,1);const auto t=materialToken(in[1]);std::array<double,2> interval{};
        for(int k=0;k<2;++k){const auto* a=in[k+2];if(!mxIsDouble(a)||mxIsComplex(a)||mxIsSparse(a)||mxGetNumberOfElements(a)!=1)throw std::invalid_argument("Material interval must be scalar double");interval[k]=mxGetScalar(a);if(!std::isfinite(interval[k])||std::floor(interval[k])!=interval[k])throw std::invalid_argument("Material interval must be finite integer");}
        if(interval[0]<0||interval[0]>4096||interval[1]<1||interval[1]>4096)throw std::invalid_argument("Material interval range");
        out[0]=radia::hdiv_cuda::materialQuadrature(*materialContexts.at(t),std::size_t(interval[0]),std::size_t(interval[1]));return;
    }
    if(command=="destroy_hdiv_material_context"){
        arity(nr,2,nl,0);materialContexts.erase(materialToken(in[1]));mexUnlock();return;
    }
    if(command=="hdiv_material_quadrature_chunk"){
        arity(nr,5,nl,1);std::array<double,3> values{};
        for(int k=0;k<3;++k){const auto* a=in[k+2];if(!mxIsDouble(a)||mxIsComplex(a)||mxIsSparse(a)||mxGetNumberOfElements(a)!=1)throw std::invalid_argument("Material chunk arguments must be scalar doubles");values[k]=mxGetScalar(a);if(!std::isfinite(values[k])||std::floor(values[k])!=values[k])throw std::invalid_argument("Material chunk arguments must be finite integers");}
        if(values[0]<1||values[0]>2||values[1]<0||values[1]>4096||values[2]<1||values[2]>4096)throw std::invalid_argument("Material chunk argument range");
        out[0]=radia::hdiv_cuda::materialQuadrature(text(in[1]),int(values[0]),std::size_t(values[1]),std::size_t(values[2]));return;
    }
    if(command=="hdiv_material_quadrature"){
        arity(nr,3,nl,1);
        if(!mxIsDouble(in[2])||mxIsComplex(in[2])||mxIsSparse(in[2])||mxGetNumberOfElements(in[2])!=1)throw std::invalid_argument("Material order must be scalar double");
        const double order=mxGetScalar(in[2]);if(!std::isfinite(order)||order<1||order>2||std::floor(order)!=order)throw std::invalid_argument("Material order must be 1..2");
        out[0]=radia::hdiv_cuda::materialQuadrature(text(in[1]),int(order));return;
    }
    if(command=="inspect_hdiv_elements"){
        arity(nr,3,nl,1);
        if(!mxIsDouble(in[2])||mxIsComplex(in[2])||mxIsSparse(in[2])||mxGetNumberOfElements(in[2])!=1)throw std::invalid_argument("HDiv order must be scalar double");
        const double order=mxGetScalar(in[2]);if(!std::isfinite(order)||order<1||order>3||std::floor(order)!=order)throw std::invalid_argument("HDiv order must be 1..3");
        out[0]=radia::hdiv_cuda::inspectElements(text(in[1]),int(order));return;
    }
    if(command=="assemble_hdiv_rt0"){
        arity(nr,3,nl,1);
        if(!mxIsDouble(in[2])||mxIsComplex(in[2])||mxIsSparse(in[2])||mxGetNumberOfElements(in[2])!=3)throw std::invalid_argument("Applied field must contain three real doubles");
        std::array<double,3> applied;for(int i=0;i<3;++i){applied[i]=mxGetDoubles(in[2])[i];if(!std::isfinite(applied[i]))throw std::invalid_argument("Nonfinite applied field");}
        out[0]=radia::rt0_cuda::assemble(text(in[1]),applied);return;
    }
    if(command=="create"){
        arity(nr,4,nl,1);auto c=std::make_unique<Context>();
        if(!mxIsDouble(in[3])||mxIsComplex(in[3])||mxGetNumberOfElements(in[3])!=1)
            throw std::runtime_error("Expected scalar real order");
        double order=mxGetScalar(in[3]);
        if(!std::isfinite(order)||order<1||order>20||std::floor(order)!=order)
            throw std::runtime_error("Invalid order");
        c->mesh=std::make_shared<ngcomp::MeshAccess>(text(in[1]));
        if(c->mesh->GetDimension()!=3)throw std::runtime_error("Only 3D vector projection supported");
        ngcore::Flags flags;flags.SetFlag("order",int(order));auto kind=text(in[2]);
        if(kind=="hdiv")c->space=std::make_shared<ngcomp::HDivHighOrderFESpace>(c->mesh,flags);
        else if(kind=="hcurl")c->space=std::make_shared<ngcomp::HCurlHighOrderFESpace>(c->mesh,flags);
        else if(kind=="vectorh1")c->space=std::make_shared<ngcomp::CompoundFESpaceAllSame>(
            std::make_shared<ngcomp::H1HighOrderFESpace>(c->mesh,flags),3,flags);
        else throw std::runtime_error("Unsupported vector FE space");
        c->space->Update();c->space->FinalizeUpdate();
        c->gf=createHeadlessGridFunction(c->space);c->gf->Update();
        radia::cuda_projection::Interpolate(c->targets,*c->gf);
        c->targets->collecting=false;c->targets->values.resize(c->targets->points.size());
        if(generation>=0x0000ffffffffffffULL)throw std::runtime_error("Projection generation exhausted");
        auto t=0x5052000000000000ULL|generation++;
        out[0]=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
        *static_cast<std::uint64_t*>(mxGetData(out[0]))=t;
        contexts.emplace(t,std::move(c));mexLock();mexAtExit(cleanup);return;
    }
    if(command=="live_count"){
        arity(nr,1,nl,1);out[0]=mxCreateDoubleScalar(contexts.size()+materialContexts.size());return;
    }
#ifdef RADIA_PROJECTION_LIFECYCLE_TEST
    if(command=="diagnostic_owner_phase"){
        arity(nr,3,nl,1);
        const auto counts=probeOwnerPhase(text(in[1]),text(in[2]));
        out[0]=mxCreateDoubleMatrix(1,counts.size(),mxREAL);
        std::copy(counts.begin(),counts.end(),mxGetDoubles(out[0]));return;
    }
    if(command=="diagnostic_released_owners"){
        arity(nr,1,nl,1);
        const char* names[]={"last_released_handle","live_projection_contexts","live_material_contexts",
            "mesh_strong_count","space_strong_count","gridfunction_strong_count","targets_strong_count"};
        out[0]=mxCreateStructMatrix(1,1,7,names);
        auto* handle=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
        *static_cast<std::uint64_t*>(mxGetData(handle))=releasedOwners.handle;
        mxSetField(out[0],0,names[0],handle);
        const double counts[]={double(contexts.size()),double(materialContexts.size()),
            double(releasedOwners.mesh.use_count()),double(releasedOwners.space.use_count()),
            double(releasedOwners.gf.use_count()),double(releasedOwners.targets.use_count())};
        for(int i=0;i<6;++i)mxSetField(out[0],0,names[i+1],mxCreateDoubleScalar(counts[i]));
        return;
    }
#endif
    if(command=="destroy"){
        arity(nr,2,nl,0);const auto handle=token(in[1]);
#ifdef RADIA_PROJECTION_LIFECYCLE_TEST
        const auto& c=*contexts.at(handle);releasedOwners.handle=handle;
        releasedOwners.mesh=c.mesh;releasedOwners.space=c.space;
        releasedOwners.gf=c.gf;releasedOwners.targets=c.targets;
#endif
        contexts.erase(handle);mexUnlock();return;
    }
    if(command=="targets"){
        arity(nr,2,nl,1);const auto& p=contexts.at(token(in[1]))->targets->points;
        out[0]=mxCreateDoubleMatrix(p.size(),3,mxREAL);double* a=mxGetDoubles(out[0]);
        for(std::size_t i=0;i<p.size();++i)for(int j=0;j<3;++j)a[i+j*p.size()]=p[i][j];return;
    }
    if(command=="project"){
        arity(nr,3,nl,2);auto& c=*contexts.at(token(in[1]));auto& v=c.targets->values;
        if(!mxIsDouble(in[2])||mxIsSparse(in[2])||mxIsComplex(in[2])||mxGetNumberOfDimensions(in[2])!=2||mxGetM(in[2])!=v.size()||mxGetN(in[2])!=3)
            throw std::runtime_error("Values must be full real double target-count by 3");
        const double* a=mxGetDoubles(in[2]);
        for(std::size_t i=0;i<3*v.size();++i)if(!std::isfinite(a[i]))throw std::runtime_error("Nonfinite field");
        for(std::size_t i=0;i<v.size();++i)for(int j=0;j<3;++j)v[i][j]=a[i+j*v.size()];
        auto profile=radia::cuda_projection::Interpolate(c.targets,*c.gf);
        auto vector=c.gf->GetVector().FVDouble();out[0]=mxCreateDoubleMatrix(vector.Size(),1,mxREAL);
        for(std::size_t i=0;i<vector.Size();++i)mxGetDoubles(out[0])[i]=vector(i);
        out[1]=mxCreateDoubleScalar(profile.interpolation_seconds);return;
    }
    throw std::runtime_error("Unknown projection command");

}
extern "C" void RadiaProjectionGateway(int nl,mxArray** out,int nr,const mxArray** in){

    char errorMessage[2048]={0};
    if(!radia_mex::checked(dispatchProjection,nl,out,nr,in,errorMessage,sizeof(errorMessage))){
        mexErrMsgIdAndTxt("radia:cuda:projection","%s",errorMessage);
    }
}
