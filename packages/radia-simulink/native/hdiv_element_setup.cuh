// NGSolve-owned high-order basis staging. Not a charge map or material solver.
#pragma once
#include "hdiv_rt0_host.cuh"
#include "../../../src/core/gpu/hdiv/hdiv_material_staging.hpp"
#include <l2hofespace.hpp>
#include <scalarfe.hpp>
#include <calcinverse.hpp>
#include <map>
#include <algorithm>
#include <limits>
namespace radia::hdiv_cuda {
inline double power(double x,int n){double v=1.;while(n--)v*=x;return v;}
inline std::vector<std::array<int,3>> monomials(int p,bool volume){
    std::vector<std::array<int,3>> result;
    for(int i=0;i<=p;++i)for(int j=0;j<=p-i;++j)for(int k=0;k<=(volume?p-i-j:0);++k)result.push_back({i,j,k});return result;
}
inline std::vector<double> chargeMap(std::shared_ptr<ngcomp::MeshAccess> mesh,
    std::shared_ptr<ngcomp::HDivHighOrderFESpace> space,int order,bool curved,
    std::vector<double>& hosts,std::vector<double>& kinds,std::vector<double>& exponents,std::vector<std::map<std::size_t,double>>* sparseRows=nullptr){
    ngcore::Flags vf,sf;vf.SetFlag("order",order-1);sf.SetFlag("order",order);
    auto volumeSpace=std::make_shared<ngcomp::L2HighOrderFESpace>(mesh,vf);
    auto surfaceSpace=std::make_shared<ngcomp::L2SurfaceHighOrderFESpace>(mesh,sf);
    volumeSpace->Update();volumeSpace->FinalizeUpdate();surfaceSpace->Update();surfaceSpace->FinalizeUpdate();
    const auto nd=space->GetNDof();std::vector<double> result;
    ngstd::LocalHeap heap(1<<24,"radia_cuda_hdiv_charge_map");
    for(int kind=0;kind<2;++kind){
        const bool volume=kind==0;auto mons=monomials(volume?order-1:order,volume);
        const auto count=mesh->GetNE(volume?ngfem::VOL:ngfem::BND);
        if((hosts.size()+count*mons.size())>8192)throw std::runtime_error("HDiv charge budget exceeded");
        for(std::size_t e=0;e<count;++e){
            ngstd::HeapReset reset(heap);ngfem::ElementId eid(volume?ngfem::VOL:ngfem::BND,int(e));
            if(!volume){int in,out;mesh->GetSElNeighbouringDomains(int(e),in,out);if(in>0&&out>0)throw std::invalid_argument("Internal interface charge staging is not connected");}
            auto& trafo=mesh->GetTrafo(eid,heap);const auto& fe=space->GetFE(eid,heap);
            const auto& dg=dynamic_cast<const ngfem::BaseScalarFiniteElement&>(volume?volumeSpace->GetFE(eid,heap):surfaceSpace->GetFE(eid,heap));
            ngstd::Array<int> dofs;space->GetDofNrs(eid,dofs);const int nm=int(mons.size()),ns=dg.GetNDof(),nf=fe.GetNDof();
            if(nm!=ns||nf!=dofs.Size())throw std::runtime_error("Unexpected simplex DG basis dimensions");
            ngbla::Matrix<> mass(ns,ns),coupling(ns,nf),monoMass(nm,nm),cross(nm,ns);mass=0.;coupling=0.;monoMass=0.;cross=0.;
            ngbla::Vector<> scalar(ns),density(nf),mono(nm);ngbla::Matrix<> normalShape(nf,3);
            auto accumulate=[&](const ngfem::IntegrationPoint& ip,double w,int mode){
                dg.CalcShape(ip,scalar);
                if(mode==0){for(int i=0;i<ns;++i)for(int j=0;j<ns;++j)mass(i,j)+=w*scalar[i]*scalar[j];}
                if(mode==1){for(int i=0;i<ns;++i)for(int j=0;j<nf;++j)coupling(i,j)+=w*scalar[i]*density[j];}
                if(mode==2){
                    const double x=curved?ip(0):ip(1),y=curved?ip(1):(volume?ip(2):1.-ip(0)-ip(1));
                    const double z=curved?ip(2):1.-ip(0)-ip(1)-ip(2);
                    for(int i=0;i<nm;++i)mono[i]=power(x,mons[i][0])*power(y,mons[i][1])*power(z,mons[i][2]);
                    for(int i=0;i<nm;++i){for(int j=0;j<nm;++j)monoMass(i,j)+=w*mono[i]*mono[j];for(int j=0;j<ns;++j)cross(i,j)+=w*mono[i]*scalar[j];}
                }
            };
            for(int mode=0;mode<3;++mode){
                const int degree=mode==0?2*dg.Order():(mode==1?fe.Order()+dg.Order()-(volume?1:0):2*order+2);
                ngfem::IntegrationRule rule(volume?ngfem::ET_TET:ngfem::ET_TRIG,degree);
                for(const auto& ip:rule){
                    if(mode==2){accumulate(ip,ip.Weight(),mode);continue;}
                    if(volume){ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);
                        if(mode==1){dynamic_cast<const ngfem::HDivFiniteElement<3>&>(fe).CalcMappedDivShape(mip,density);density*=-1.;}
                        accumulate(ip,mip.GetWeight(),mode);
                    }else{ngfem::MappedIntegrationPoint<2,3> mip(ip,trafo);
                        if(mode==1){dynamic_cast<const ngfem::HDivNormalFiniteElement<2>&>(fe).CalcMappedShape(mip,normalShape);const auto normal=mip.GetNV();for(int i=0;i<nf;++i){density[i]=0.;for(int c=0;c<3;++c)density[i]+=normalShape(i,c)*normal[c];}}
                        accumulate(ip,mip.GetWeight(),mode);
                    }
                }
            }
            // The two original DG inversions and oriented change of basis,
            // retained separately to match _charge_basis[_curved].
            ngbla::CalcInverse(mass);ngbla::CalcInverse(monoMass);
            ngbla::Matrix<> densityMap(ns,nf),transform(nm,ns),local(nm,nf);
            densityMap=mass*coupling;transform=monoMass*cross;local=transform*densityMap;
            for(int i=0;i<nm;++i){hosts.push_back(double(e));kinds.push_back(double(kind));for(int c=0;c<3;++c)exponents.push_back(mons[i][c]);
                const auto offset=result.size();if(sparseRows)sparseRows->emplace_back();else result.resize(offset+nd,0.);
                for(int j=0;j<nf;++j){if(dofs[j]<0||std::size_t(dofs[j])>=nd)throw std::invalid_argument("Constrained HDiv charge DOF unsupported");if(sparseRows)sparseRows->back()[dofs[j]]=local(i,j);else result[offset+dofs[j]]=local(i,j);}
            }
        }
    }
    if(curved){double maximum=0.;
        if(sparseRows){for(const auto& row:*sparseRows)for(const auto& item:row)maximum=std::max(maximum,std::fabs(item.second));}
        else for(double v:result)maximum=std::max(maximum,std::fabs(v));
        const double floor=4096.*std::numeric_limits<double>::epsilon()*maximum;
        if(sparseRows){for(auto& row:*sparseRows)for(auto& item:row)if(std::fabs(item.second)<floor)item.second=0.;}
        else for(auto& v:result)if(std::fabs(v)<floor)v=0.;
    }
    return result;
}
// MATLAB owns this CSC layout. Domain assembly above remains per-element.
inline mxArray* sparseRowsMatrix(const std::vector<std::map<std::size_t,double>>& entries,std::size_t columnsCount){
    std::size_t nnz=0;for(const auto& row:entries)for(const auto& entry:row)if(entry.second!=0.)++nnz;
    auto* matrix=mxCreateSparse(entries.size(),columnsCount,nnz,mxREAL);auto* columns=mxGetJc(matrix);auto* rows=mxGetIr(matrix);auto* values=mxGetDoubles(matrix);
    for(const auto& row:entries)for(const auto& entry:row)if(entry.second!=0.)++columns[entry.first+1];
    for(std::size_t col=0;col<columnsCount;++col)columns[col+1]+=columns[col];
    std::vector<mwIndex> next(columns,columns+columnsCount);
    for(std::size_t row=0;row<entries.size();++row)for(const auto& entry:entries[row])if(entry.second!=0.){const auto offset=next[entry.first]++;rows[offset]=row;values[offset]=entry.second;}
    return matrix;
}
inline mxArray* inspectElements(const std::string& path,int order,const std::array<double,3>& applied={0.,0.,0.},bool sparseMass=false,bool sparseCharge=false){
    auto mesh=std::make_shared<ngcomp::MeshAccess>(path);const auto ne=mesh->GetNE(ngfem::VOL);
    if(mesh->GetCurveOrder()>2)throw std::invalid_argument("HDiv GPU geometry supports at most the existing P2 map");
    if(mesh->GetDimension()!=3||!ne||ne>4096||mesh->GetDeformation()||order<1||order>3)
        throw std::invalid_argument("HDiv staging requires orders 1..3 and 1..4096 undeformed tetrahedra");
    ngcore::Flags flags;flags.SetFlag("order",order);
    auto space=std::make_shared<ngcomp::HDivHighOrderFESpace>(mesh,flags);space->Update();space->FinalizeUpdate();
    const auto nd=space->GetNDof();if(!nd||nd>16384)throw std::runtime_error("HDiv staging DOF budget exceeded");
    // Count local contributions before allocating global storage. Duplicate mass
    // entries only reduce this bound; charge rows belong to individual elements.
    const auto checkedAdd=[](std::size_t a,std::size_t b){
        if(b>std::numeric_limits<std::size_t>::max()-a)throw std::overflow_error("HDiv staging size addition overflow");return a+b;};
    const auto checkedMul=[](std::size_t a,std::size_t b){
        if(a&&b>std::numeric_limits<std::size_t>::max()/a)throw std::overflow_error("HDiv staging size product overflow");return a*b;};
    const auto volumeMonomials=monomials(order-1,true).size(),surfaceMonomials=monomials(order,false).size();
    const auto chargeCount=checkedAdd(checkedMul(std::size_t(ne),volumeMonomials),checkedMul(std::size_t(mesh->GetNE(ngfem::BND)),surfaceMonomials));
    if(chargeCount>8192)throw std::runtime_error("HDiv charge budget exceeded before staging");
    std::size_t massEntries=0,chargeEntries=0;
    if(sparseMass||sparseCharge)for(bool volume:{true,false}){
        const auto region=volume?ngfem::VOL:ngfem::BND;
        for(std::size_t e=0;e<std::size_t(mesh->GetNE(region));++e){
            ngstd::Array<int> dofs;space->GetDofNrs(ngfem::ElementId(region,int(e)),dofs);
            for(auto d:dofs)if(d<0||std::size_t(d)>=nd)throw std::invalid_argument("Constrained HDiv staging unsupported");
            const auto count=std::size_t(dofs.Size());
            if(volume&&sparseMass)massEntries=checkedAdd(massEntries,checkedMul(count,count));
            if(sparseCharge)chargeEntries=checkedAdd(chargeEntries,checkedMul(count,volume?volumeMonomials:surfaceMonomials));
        }
    }
    // 128 bytes/contribution is a conservative admission allowance for the
    // MSVC map node plus concurrent CSC values/indices, not a measured peak.
    const auto sparseBytes=[&](std::size_t entries,std::size_t rows){
        return checkedAdd(checkedMul(entries,128),checkedAdd(checkedMul(rows,sizeof(std::map<std::size_t,double>)),checkedMul(checkedAdd(std::size_t(nd),1),2*sizeof(mwIndex))));};
    const auto massBytes=sparseMass?sparseBytes(massEntries,nd):checkedMul(checkedMul(std::size_t(nd),std::size_t(nd)),2*sizeof(double));
    const auto mapBytes=sparseCharge?sparseBytes(chargeEntries,chargeCount):checkedMul(checkedMul(chargeCount,std::size_t(nd)),2*sizeof(double));
    if(massBytes>(std::size_t(4)<<30)||mapBytes>(std::size_t(4)<<30)||checkedAdd(massBytes,mapBytes)>(std::size_t(6)<<30))
        throw std::runtime_error("HDiv staging allocation budget exceeded");
    // Reference-node order from _charge_basis_curved. Always use GetTrafo.
    const double refs[10][3]={{0,0,0},{1,0,0},{0,1,0},{0,0,1},{.5,0,0},{.5,.5,0},{0,.5,0},{0,0,.5},{.5,0,.5},{0,.5,.5}};
    std::vector<double> nodes(ne*30),mass(sparseMass?0:nd*nd,0.),rhs(nd,0.),mean(3*nd,0.),cellIds,faceIds,faceNodes;std::size_t curved=0;double volumeTotal=0.;
    std::vector<std::map<std::size_t,double>> massRows(sparseMass?nd:0);
    ngstd::LocalHeap heap(1<<24,"radia_cuda_hdiv_element_staging");
    for(std::size_t e=0;e<ne;++e){
        ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::VOL,int(e));
        if(mesh->GetElType(eid)!=ngfem::ET_TET)throw std::invalid_argument("HDiv staging rejects non-tetrahedra");
        const bool isCurved=mesh->GetElement(eid).is_curved;curved+=isCurved;
        if(isCurved&&order>2)throw std::invalid_argument("Existing curved charge semantics support orders 1..2");
        auto& trafo=mesh->GetTrafo(eid,heap);const auto& fe=dynamic_cast<const ngfem::HDivFiniteElement<3>&>(space->GetFE(eid,heap));
        for(auto id:mesh->GetElVertices(eid))cellIds.push_back(double(id));
        ngstd::Array<int> dofs;space->GetDofNrs(eid,dofs);
        if(dofs.Size()!=fe.GetNDof())throw std::runtime_error("HDiv local DOF mismatch");
        for(auto d:dofs)if(d<0||std::size_t(d)>=nd)throw std::invalid_argument("Constrained HDiv staging unsupported");
        for(int n=0;n<10;++n){ngfem::IntegrationPoint ip(refs[n][0],refs[n][1],refs[n][2]);ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);for(int c=0;c<3;++c)nodes[e*30+n*3+c]=mip.GetPoint()[c];}
        ngfem::IntegrationRule rule(ngfem::ET_TET,2*fe.Order());ngbla::Matrix<> shape(fe.GetNDof(),3);
        for(const auto& ip:rule){ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);fe.CalcMappedShape(mip,shape);const double w=mip.GetWeight();
            if(!std::isfinite(w)||w<=0)throw std::runtime_error("Invalid mapped measure");
            for(int i=0;i<dofs.Size();++i)for(int c=0;c<3;++c){rhs[dofs[i]]+=w*shape(i,c)*applied[c];mean[c*nd+dofs[i]]+=w*shape(i,c);}
            for(int i=0;i<dofs.Size();++i)for(int j=0;j<dofs.Size();++j){double dot=0.;for(int c=0;c<3;++c)dot+=shape(i,c)*shape(j,c);if(sparseMass)massRows[dofs[i]][dofs[j]]+=w*dot;else mass[dofs[i]*nd+dofs[j]]+=w*dot;}
        }
        // Existing CPU mean magnetization divides by ng.Integrate(1,mesh),
        // whose declared default order is 5, independently of the mass rule.
        ngfem::IntegrationRule volumeRule(ngfem::ET_TET,5);
        for(const auto& ip:volumeRule){ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);volumeTotal+=mip.GetWeight();}
    }
    const double triRefs[6][2]={{0,0},{1,0},{0,1},{.5,0},{.5,.5},{0,.5}};
    for(std::size_t f=0;f<mesh->GetNE(ngfem::BND);++f){ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::BND,int(f));
        if(mesh->GetElType(eid)!=ngfem::ET_TRIG)throw std::invalid_argument("HDiv curved skin requires triangles");
        for(auto id:mesh->GetElVertices(eid))faceIds.push_back(double(id));
        auto& trafo=mesh->GetTrafo(eid,heap);for(const auto& ref:triRefs){ngfem::IntegrationPoint ip(ref[0],ref[1]);ngfem::MappedIntegrationPoint<2,3> mip(ip,trafo);for(int c=0;c<3;++c)faceNodes.push_back(mip.GetPoint()[c]);}
    }
    if(volumeTotal<=0||!std::isfinite(volumeTotal))throw std::runtime_error("Invalid HDiv volume");for(auto& v:mean)v/=volumeTotal;
    std::vector<double> hosts,kinds,expo;std::vector<std::map<std::size_t,double>> chargeRows;auto B=chargeMap(mesh,space,order,curved!=0,hosts,kinds,expo,sparseCharge?&chargeRows:nullptr);
    const char* fields[]={"cell_nodes","mass","curved_cells","order","charge_map","host","kind","expo","face_nodes","cell_vertices","face_vertices","rhs","mean_operator"};auto* out=mxCreateStructMatrix(1,1,13,fields);
    mxSetField(out,0,"cell_nodes",radia::rt0_cuda::matrix(ne,30,nodes));if(sparseMass)mxSetField(out,0,"mass",sparseRowsMatrix(massRows,nd));else mxSetField(out,0,"mass",radia::rt0_cuda::matrix(nd,nd,mass));
    mxSetField(out,0,"curved_cells",mxCreateDoubleScalar(double(curved)));mxSetField(out,0,"order",mxCreateDoubleScalar(order));
    mxSetField(out,0,"face_nodes",radia::rt0_cuda::matrix(faceIds.size()/3,18,faceNodes));mxSetField(out,0,"cell_vertices",radia::rt0_cuda::matrix(ne,4,cellIds));mxSetField(out,0,"face_vertices",radia::rt0_cuda::matrix(faceIds.size()/3,3,faceIds));
    mxSetField(out,0,"rhs",radia::rt0_cuda::matrix(nd,1,rhs));mxSetField(out,0,"mean_operator",radia::rt0_cuda::matrix(3,nd,mean));
    if(sparseCharge)mxSetField(out,0,"charge_map",sparseRowsMatrix(chargeRows,nd));else mxSetField(out,0,"charge_map",radia::rt0_cuda::matrix(hosts.size(),nd,B));mxSetField(out,0,"host",radia::rt0_cuda::matrix(hosts.size(),1,hosts));mxSetField(out,0,"kind",radia::rt0_cuda::matrix(kinds.size(),1,kinds));mxSetField(out,0,"expo",radia::rt0_cuda::matrix(hosts.size(),3,expo));return out;
}
// The publisher owns MATLAB allocations; staging above it owns no MATLAB ABI.
inline mxArray* materialQuadrature(MaterialStagingContext& context,std::size_t first=0,std::size_t count=std::numeric_limits<std::size_t>::max()){
    const auto payload=stageMaterialQuadrature(context,first,count);
    const auto nq=payload.weights.size(),nd=payload.dofs;
    if(nd>std::numeric_limits<mwIndex>::max()||3*nq>std::numeric_limits<mwIndex>::max()||payload.values.size()>std::numeric_limits<mwIndex>::max())
        throw std::overflow_error("Material MATLAB index width exceeded");
    const char* names[]={"evaluation","weights","points","element","integration_order","dofs","scope","global_point_indices","total_points","total_elements","first_element","element_count","potential_shape_entries","tangent_nonzero_bound"};auto* out=mxCreateStructMatrix(1,1,14,names);
    auto* evaluation=mxCreateSparse(3*nq,nd,payload.values.size(),mxREAL);
    auto* jc=mxGetJc(evaluation);auto* ir=mxGetIr(evaluation);auto* values=mxGetDoubles(evaluation);
    for(std::size_t i=0;i<payload.columnOffsets.size();++i)jc[i]=mwIndex(payload.columnOffsets[i]);
    for(std::size_t i=0;i<payload.values.size();++i){ir[i]=mwIndex(payload.rowIndices[i]);values[i]=payload.values[i];}
    mxSetField(out,0,"evaluation",evaluation);mxSetField(out,0,"weights",radia::rt0_cuda::matrix(nq,1,payload.weights));
    mxSetField(out,0,"points",radia::rt0_cuda::matrix(nq,3,payload.points));mxSetField(out,0,"element",radia::rt0_cuda::matrix(nq,1,payload.elements));
    mxSetField(out,0,"integration_order",mxCreateDoubleScalar(payload.integrationOrder));mxSetField(out,0,"dofs",mxCreateDoubleScalar(double(nd)));
    std::vector<double> indices(payload.globalPointIndices.begin(),payload.globalPointIndices.end());
    mxSetField(out,0,"global_point_indices",radia::rt0_cuda::matrix(nq,1,indices));
    mxSetField(out,0,"total_points",mxCreateDoubleScalar(double(payload.totalPoints)));mxSetField(out,0,"total_elements",mxCreateDoubleScalar(double(payload.totalElements)));
    mxSetField(out,0,"first_element",mxCreateDoubleScalar(double(payload.firstElement)));mxSetField(out,0,"element_count",mxCreateDoubleScalar(double(payload.elementCount)));
    mxSetField(out,0,"potential_shape_entries",mxCreateDoubleScalar(double(payload.potentialShapeEntries)));mxSetField(out,0,"tangent_nonzero_bound",mxCreateDoubleScalar(double(payload.tangentNonzeroBound)));
    mxSetField(out,0,"scope",mxCreateString("NGSolve mapped material quadrature; component-major evaluation, zero-based elements; no nonlinear solve acceptance"));return out;
}

inline mxArray* materialQuadrature(const std::string& path,int order,std::size_t first=0,std::size_t count=std::numeric_limits<std::size_t>::max()){
    MaterialStagingContext context(path,order);return materialQuadrature(context,first,count);
}

}
