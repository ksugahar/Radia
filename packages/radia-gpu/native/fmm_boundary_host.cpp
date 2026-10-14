/* SPDX-License-Identifier: MPL-2.0 */
// Native checked flat-P1 boundary maps; shared SS integration owns near blocks.
#include <mex.h>
#include <h1hofespace.hpp>
#include <scalarfe.hpp>
#include <meshaccess.hpp>
#include "rad_bem_galerkin.h"
#include <array>
#include <map>
#include <set>
#include <vector>
#include <cmath>
#include <climits>
#include <stdexcept>
namespace {
using Own=std::unique_ptr<mxArray,decltype(&mxDestroyArray)>;
Own own(mxArray* p){if(!p)throw std::bad_alloc();return Own(p,mxDestroyArray);}
std::string text(const mxArray* p){if(!p||!mxIsChar(p))throw std::invalid_argument("Character boundary metadata required");char* c=mxArrayToString(p);if(!c)throw std::bad_alloc();std::string s(c);mxFree(c);return s;}
void put(mxArray* p,const char* name,Own v){mxSetField(p,0,name,v.release());}
Own vector(const std::vector<double>& v){auto p=own(mxCreateDoubleMatrix(v.size(),1,mxREAL));std::copy(v.begin(),v.end(),mxGetDoubles(p.get()));return p;}
Own matrix(const std::vector<double>& v,int cols){size_t n=v.size()/cols;auto p=own(mxCreateDoubleMatrix(n,cols,mxREAL));for(size_t i=0;i<n;++i)for(int c=0;c<cols;++c)mxGetDoubles(p.get())[i+n*c]=v[i*cols+c];return p;}
}
mxArray* RadiaPrepareBoundary(const mxArray* path,const mxArray* labels){
    auto mesh=std::make_shared<ngcomp::MeshAccess>(text(path));size_t nt=mesh->GetNE(ngfem::BND);
    if(mesh->GetDimension()!=3||!nt||nt>4096||mesh->GetCurveOrder()>1||mesh->GetDeformation())throw std::invalid_argument("Boundary FMM preparation supports 1..4096 flat exterior triangles only");
    if(!mxIsCell(labels))throw std::invalid_argument("Boundary label contract must be a cell array");
    std::set<std::string> expected,actual;for(size_t i=0;i<mxGetNumberOfElements(labels);++i)expected.insert(text(mxGetCell(labels,i)));
    for(auto name:mesh->GetMaterials(ngfem::BND))actual.emplace(name);
    if(expected.empty()||expected!=actual||expected.size()!=mxGetNumberOfElements(labels))throw std::invalid_argument("Boundary FMM label contract mismatch");
    ngcore::Flags flags;flags.SetFlag("order",1);auto space=std::make_shared<ngcomp::H1HighOrderFESpace>(mesh,flags);space->Update();space->FinalizeUpdate();
    ngstd::LocalHeap heap(1<<22,"radia_cuda_boundary_maps");ngfem::IntegrationRule rule(ngfem::ET_TRIG,16);int nq=int(rule.Size());
    std::vector<int> active;std::map<int,int> numbering;
    for(int t=0;t<int(nt);++t){ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::BND,t);ngstd::Array<int> dofs;space->GetDofNrs(eid,dofs);for(int d:dofs)if(!numbering.count(d)){if(d<0)throw std::invalid_argument("Constrained boundary DOF unsupported");numbering[d]=int(active.size());active.push_back(d);}}
    int nv=int(active.size());std::vector<double> vertices(size_t(nv)*3),nodes(nt*18),points,normals,mapRows,mapCols,mapValues,centers(nt*3),radius(nt),weights(nt*nq),shape(nt*nq*3);
    std::vector<int64_t> triangles(nt*3);std::vector<std::array<int,3>> localDofs(nt);
    for(int i=0;i<nv;++i){auto p=mesh->GetPoint<3>(active[i]);for(int c=0;c<3;++c)vertices[i*3+c]=p[c];}
    const double refs[3][2]={{0,0},{1,0},{0,1}};
    for(int t=0;t<int(nt);++t){ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::BND,t);int inside,outside;mesh->GetSElNeighbouringDomains(t,inside,outside);
        if(mesh->GetElType(eid)!=ngfem::ET_TRIG||mesh->GetElement(eid).is_curved||inside<=0||outside>0)throw std::invalid_argument("Boundary FMM rejects curved, internal or unoriented surface elements");
        auto& trafo=mesh->GetTrafo(eid,heap);auto& fe=dynamic_cast<const ngfem::BaseScalarFiniteElement&>(space->GetFE(eid,heap));ngstd::Array<int> dofs;space->GetDofNrs(eid,dofs);
        if(dofs.Size()!=3||fe.GetNDof()!=3)throw std::invalid_argument("Boundary FMM requires scalar P1 traces");ngbla::Vector<> values(3);
        for(int a=0;a<3;++a){ngfem::IntegrationPoint ip(refs[a][0],refs[a][1]);fe.CalcShape(ip,values);int match=-1;for(int j=0;j<3;++j)if(std::abs(values[j]-1)<1e-12)match=j;
            if(match<0)throw std::runtime_error("NGSolve P1 reference-corner contract failed");int d=numbering.at(dofs[match]);localDofs[t][a]=match;triangles[t*3+a]=d;
            ngfem::MappedIntegrationPoint<2,3> mip(ip,trafo);for(int c=0;c<3;++c)nodes[t*18+a*3+c]=mip.GetPoint()[c];}
        for(int a=0;a<3;++a)for(int c=0;c<3;++c)nodes[t*18+(a+3)*3+c]=.5*(nodes[t*18+a*3+c]+nodes[t*18+((a+1)%3)*3+c]);
        for(int c=0;c<3;++c)centers[t*3+c]=(nodes[t*18+c]+nodes[t*18+3+c]+nodes[t*18+6+c])/3;
        for(int a=0;a<3;++a){double d2=0;for(int c=0;c<3;++c){double d=nodes[t*18+a*3+c]-centers[t*3+c];d2+=d*d;}radius[t]=std::max(radius[t],std::sqrt(d2));}
        for(int q=0;q<nq;++q){ngfem::MappedIntegrationPoint<2,3> mip(rule[q],trafo);fe.CalcShape(rule[q],values);double w=mip.GetWeight();if(!std::isfinite(w)||w<=0)throw std::runtime_error("Invalid boundary mapped weight");weights[t*nq+q]=w;
            for(int c=0;c<3;++c){points.push_back(mip.GetPoint()[c]);normals.push_back(mip.GetNV()[c]);}
            for(int a=0;a<3;++a){double v=values[localDofs[t][a]];shape[(t*nq+q)*3+a]=v;mapRows.push_back(double(t*nq+q));mapCols.push_back(double(triangles[t*3+a]));mapValues.push_back(w*v);}}
    }
    radia::bem::P1EntryProvider provider(vertices.data(),nv,triangles.data(),int(nt),nodes.data(),17,10,0);
    std::map<std::pair<int,int>,std::array<double,2>> correction;size_t pairs=0;
    constexpr double scale=1./(4*3.141592653589793238462643383279502884);
    for(int a=0;a<int(nt);++a)for(int b=0;b<int(nt);++b){double d2=0;for(int c=0;c<3;++c){double d=centers[a*3+c]-centers[b*3+c];d2+=d*d;}double bound=2*(radius[a]+radius[b]);if(d2>bound*bound)continue;
        if(++pairs>2000000)throw std::overflow_error("Boundary exact-near pair budget exceeded");double sl[9],dl[9];provider.Pair(a,b,sl,dl);
        for(int i=0;i<nq;++i)for(int j=0;j<nq;++j){int ti=a*nq+i,sj=b*nq+j;double r[3],r2=0,dot=0;for(int c=0;c<3;++c){r[c]=points[ti*3+c]-points[sj*3+c];r2+=r[c]*r[c];dot+=r[c]*normals[sj*3+c];}if(r2<=1e-300)continue;
            double k=scale/std::sqrt(r2),kd=k*dot/r2,w=weights[ti]*weights[sj];for(int x=0;x<3;++x)for(int y=0;y<3;++y){double v=w*shape[ti*3+x]*shape[sj*3+y];sl[x*3+y]-=v*k;dl[x*3+y]-=v*kd;}}
        for(int x=0;x<3;++x)for(int y=0;y<3;++y){auto& entry=correction[{int(triangles[a*3+x]),int(triangles[b*3+y])}];entry[0]+=sl[x*3+y];entry[1]+=dl[x*3+y];}}
    std::vector<double> rows,cols,sl,dl,global;for(const auto& [key,v]:correction){rows.push_back(key.first);cols.push_back(key.second);sl.push_back(v[0]);dl.push_back(v[1]);}for(int d:active)global.push_back(d);
    const char* names[]={"schema","points","normals","map_rows","map_columns","map_values","correction_rows","correction_columns","correction_sl","correction_dl","global_dofs","dofs","near_pairs"};auto result=own(mxCreateStructMatrix(1,1,13,names));
    put(result.get(),"schema",own(mxCreateString("radia.flat-p1-boundary-fmm.v1")));put(result.get(),"points",matrix(points,3));put(result.get(),"normals",matrix(normals,3));
    put(result.get(),"map_rows",vector(mapRows));put(result.get(),"map_columns",vector(mapCols));put(result.get(),"map_values",vector(mapValues));put(result.get(),"correction_rows",vector(rows));put(result.get(),"correction_columns",vector(cols));put(result.get(),"correction_sl",vector(sl));put(result.get(),"correction_dl",vector(dl));put(result.get(),"global_dofs",vector(global));put(result.get(),"dofs",own(mxCreateDoubleScalar(nv)));put(result.get(),"near_pairs",own(mxCreateDoubleScalar(pairs)));return result.release();
}
