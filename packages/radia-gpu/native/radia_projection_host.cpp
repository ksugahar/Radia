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
#include <vector>
#include <limits>
#include <cmath>
#include <cstdint>
#include "hdiv_rt0_host.cuh"

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
std::map<std::uint64_t,std::unique_ptr<Context>> contexts;
std::uint64_t generation=1;
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
void cleanup(){contexts.clear();}
void arity(int nr,int expected,int nl,int outputs){
    if(nr!=expected||nl!=outputs)throw std::runtime_error("Invalid command arity");
}
}
extern "C" void RadiaProjectionGateway(int nl,mxArray* out[],int nr,const mxArray* in[]){
try{
    if(nr<1)throw std::runtime_error("Command required");
    auto command=text(in[0]);
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
        c->gf=ngcomp::CreateGridFunction(c->space,"radia_cuda_projection",ngcore::Flags());c->gf->Update();
        radia::cuda_projection::Interpolate(c->targets,*c->gf);
        c->targets->collecting=false;c->targets->values.resize(c->targets->points.size());
        if(generation>=0x0000ffffffffffffULL)throw std::runtime_error("Projection generation exhausted");
        auto t=0x5052000000000000ULL|generation++;
        out[0]=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
        *static_cast<std::uint64_t*>(mxGetData(out[0]))=t;
        contexts.emplace(t,std::move(c));mexLock();mexAtExit(cleanup);return;
    }
    if(command=="live_count"){
        arity(nr,1,nl,1);out[0]=mxCreateDoubleScalar(contexts.size());return;
    }
    if(command=="destroy"){
        arity(nr,2,nl,0);contexts.erase(token(in[1]));mexUnlock();return;
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
}catch(const std::exception& e){mexErrMsgIdAndTxt("radia:cuda:projection", "%s",e.what());}
}
