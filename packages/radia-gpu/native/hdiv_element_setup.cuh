// NGSolve-owned high-order basis staging. Not a charge map or material solver.
#pragma once
#include "hdiv_rt0_host.cuh"
namespace radia::hdiv_cuda {
inline mxArray* inspectElements(const std::string& path,int order){
    auto mesh=std::make_shared<ngcomp::MeshAccess>(path);const auto ne=mesh->GetNE(ngfem::VOL);
    if(mesh->GetDimension()!=3||!ne||ne>128||mesh->GetDeformation()||order<1||order>3)
        throw std::invalid_argument("HDiv staging requires orders 1..3 and 1..128 undeformed tetrahedra");
    ngcore::Flags flags;flags.SetFlag("order",order);
    auto space=std::make_shared<ngcomp::HDivHighOrderFESpace>(mesh,flags);space->Update();space->FinalizeUpdate();
    const auto nd=space->GetNDof();if(!nd||nd>1024)throw std::runtime_error("HDiv staging DOF budget exceeded");
    // Reference-node order from _charge_basis_curved. Always use GetTrafo.
    const double refs[10][3]={{0,0,0},{1,0,0},{0,1,0},{0,0,1},{.5,0,0},{.5,.5,0},{0,.5,0},{0,0,.5},{.5,0,.5},{0,.5,.5}};
    std::vector<double> nodes(ne*30),mass(nd*nd,0.);std::size_t curved=0;
    ngstd::LocalHeap heap(1<<24,"radia_cuda_hdiv_element_staging");
    for(std::size_t e=0;e<ne;++e){
        ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::VOL,int(e));
        if(mesh->GetElType(eid)!=ngfem::ET_TET)throw std::invalid_argument("HDiv staging rejects non-tetrahedra");
        const bool isCurved=mesh->GetElement(eid).is_curved;curved+=isCurved;
        if(isCurved&&order>2)throw std::invalid_argument("Existing curved charge semantics support orders 1..2");
        auto& trafo=mesh->GetTrafo(eid,heap);const auto& fe=dynamic_cast<const ngfem::HDivFiniteElement<3>&>(space->GetFE(eid,heap));
        ngstd::Array<int> dofs;space->GetDofNrs(eid,dofs);
        if(dofs.Size()!=fe.GetNDof())throw std::runtime_error("HDiv local DOF mismatch");
        for(auto d:dofs)if(d<0||std::size_t(d)>=nd)throw std::invalid_argument("Constrained HDiv staging unsupported");
        for(int n=0;n<10;++n){ngfem::IntegrationPoint ip(refs[n][0],refs[n][1],refs[n][2]);ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);for(int c=0;c<3;++c)nodes[e*30+n*3+c]=mip.GetPoint()[c];}
        ngfem::IntegrationRule rule(ngfem::ET_TET,2*fe.Order());ngbla::Matrix<> shape(fe.GetNDof(),3);
        for(const auto& ip:rule){ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);fe.CalcMappedShape(mip,shape);const double w=mip.GetWeight();
            if(!std::isfinite(w)||w<=0)throw std::runtime_error("Invalid mapped measure");
            for(int i=0;i<dofs.Size();++i)for(int j=0;j<dofs.Size();++j)for(int c=0;c<3;++c)mass[dofs[i]*nd+dofs[j]]+=w*shape(i,c)*shape(j,c);
        }
    }
    const char* fields[]={"cell_nodes","mass","curved_cells","order"};auto* out=mxCreateStructMatrix(1,1,4,fields);
    mxSetField(out,0,"cell_nodes",radia::rt0_cuda::matrix(ne,30,nodes));mxSetField(out,0,"mass",radia::rt0_cuda::matrix(nd,nd,mass));
    mxSetField(out,0,"curved_cells",mxCreateDoubleScalar(double(curved)));mxSetField(out,0,"order",mxCreateDoubleScalar(order));return out;
}
}
