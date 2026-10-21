// Owned NGSolve material staging, independent of Python and MATLAB ABIs.
#pragma once
#include <hdivfe.hpp>
#include <hdivhofespace.hpp>
#include <meshaccess.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
namespace radia::hdiv_cuda {
// CSC evaluation rows are component-major; coordinates are row-major triples.
// All indices are zero-based and owned. No NGSolve or publisher buffers escape.
struct MaterialQuadraturePayload {
    std::size_t dofs=0,totalPoints=0,totalElements=0,firstElement=0,elementCount=0;
    std::size_t potentialShapeEntries=0,tangentNonzeroBound=0;
    int integrationOrder=0;
    std::vector<std::size_t> columnOffsets,rowIndices,globalPointIndices;
    std::vector<double> values,weights,points,elements;
};
// NGSolve owns the material integration rule and mapped HDiv shapes.
// Numeric component-major evaluation map, not nodal interpolation.
struct MaterialStagingContext {
    std::shared_ptr<ngcomp::MeshAccess> mesh;
    std::shared_ptr<ngcomp::HDivHighOrderFESpace> space;
    std::shared_ptr<ngcomp::FESpace> materialSpace;
    std::size_t ne=0,nd=0,totalPoints=0;
    int degree=0;
    MaterialStagingContext(const std::string& path,int order){
        mesh=std::make_shared<ngcomp::MeshAccess>(path);ne=std::size_t(mesh->GetNE(ngfem::VOL));
        if(mesh->GetDimension()!=3||!ne||ne>4096||mesh->GetDeformation()||mesh->GetCurveOrder()>2||order<1||order>2)
            throw std::invalid_argument("Material quadrature requires bounded P1/P2 tetrahedral, hexahedral or pyramidal HDiv");
        ngcore::Flags flags;flags.SetFlag("order",order);space=std::make_shared<ngcomp::HDivHighOrderFESpace>(mesh,flags);space->Update();space->FinalizeUpdate();
        nd=std::size_t(space->GetNDof());if(!nd||nd>16384)throw std::invalid_argument("Material DOF budget exceeded");
        degree=std::max(6,2*order+2);ngcore::Flags irFlags;irFlags.SetFlag("order",degree);
        materialSpace=ngcomp::CreateFESpace("irspace",mesh,irFlags);if(!materialSpace)throw std::runtime_error("NGSolve integration space unavailable");materialSpace->Update();materialSpace->FinalizeUpdate();
        // NGSolve v6.2.2608 IntegrationRuleSpace uses IntegrationRule(type,2*order).
        // https://github.com/NGSolve/ngsolve/blob/v6.2.2608/comp/irspace.cpp
        // Retain the exported factory; its constructor is not a Windows DLL export.
        totalPoints=std::size_t(materialSpace->GetNDof());
        if(!totalPoints||totalPoints>2000000)throw std::invalid_argument("Material total-point budget exceeded");
        for(std::size_t e=0;e<ne;++e){
            const auto type=mesh->GetElType(ngfem::ElementId(ngfem::VOL,int(e)));
            if(type!=ngfem::ET_TET&&type!=ngfem::ET_HEX&&type!=ngfem::ET_PYRAMID)
                throw std::invalid_argument("Material staging supports tetrahedra, hexahedra and pyramids only; charge geometry acceptance is separate");
            ngfem::IntegrationRule rule(type,2*degree);
            if(rule.Size()<1||rule.Size()>4096)throw std::invalid_argument("Material element-point budget exceeded");
        }
    }
};
inline MaterialQuadraturePayload stageMaterialQuadrature(MaterialStagingContext& context,std::size_t first=0,std::size_t count=std::numeric_limits<std::size_t>::max()){
    const auto& mesh=context.mesh;const auto& space=context.space;const auto& materialSpace=context.materialSpace;
    const auto ne=context.ne,nd=context.nd,totalPoints=context.totalPoints;const int degree=context.degree;
    if(count==std::numeric_limits<std::size_t>::max())count=std::size_t(ne);
    if(first>=std::size_t(ne)||!count||count>std::size_t(ne)-first)throw std::invalid_argument("Invalid material element interval");
    std::size_t entries=0,tangentBound=0;std::vector<std::size_t> globalPoints;ngstd::LocalHeap heap(1<<24,"radia_hdiv_material_quadrature");
    for(std::size_t e=first;e<first+count;++e){ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::VOL,int(e));
        ngfem::IntegrationRule rule(mesh->GetElType(eid),2*degree);
        const auto& fe=dynamic_cast<const ngfem::HDivFiniteElement<3>&>(space->GetFE(eid,heap));
        if(fe.GetNDof()<1||fe.GetNDof()>256)throw std::invalid_argument("Material local shape budget exceeded");
        const auto local=std::size_t(fe.GetNDof()),points=std::size_t(rule.Size());
        if(local>(4000000-entries)/(points*3))throw std::invalid_argument("Material sparse staging budget exceeded");
        entries+=local*points*3;
        if(local*local>std::numeric_limits<std::size_t>::max()-tangentBound)throw std::overflow_error("Material tangent bound overflow");
        tangentBound+=local*local;
        ngstd::Array<int> sampleDofs;materialSpace->GetDofNrs(eid,sampleDofs);
        if(sampleDofs.Size()!=rule.Size())throw std::runtime_error("Material rule point count mismatch");
        if(std::size_t(sampleDofs.Size())>65536-globalPoints.size())throw std::invalid_argument("Material point budget exceeded");
        for(auto index:sampleDofs){if(index<0||std::size_t(index)>=totalPoints)throw std::runtime_error("Invalid material global sample index");globalPoints.push_back(std::size_t(index));}
    }
    std::sort(globalPoints.begin(),globalPoints.end());
    if(std::adjacent_find(globalPoints.begin(),globalPoints.end())!=globalPoints.end())throw std::runtime_error("Material samples are not element-disjoint");
    const auto nq=globalPoints.size();
    if(first==0&&count==std::size_t(ne)&&nq!=totalPoints)throw std::runtime_error("Material global sample coverage mismatch");
    // Each global material point belongs to one element. Its local HDiv DOFs
    // are unique, so no arithmetic reduction is needed for this evaluation
    // map. Append by column rather than allocating a tree node per nonzero.
    std::vector<std::vector<std::pair<std::size_t,double>>> columns(nd);
    std::vector<double> weights(nq),points(3*nq),elements(nq);
    for(std::size_t e=first;e<first+count;++e){ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::VOL,int(e));
        ngfem::IntegrationRule rule(mesh->GetElType(eid),2*degree);
        auto& trafo=mesh->GetTrafo(eid,heap);const auto& fe=dynamic_cast<const ngfem::HDivFiniteElement<3>&>(space->GetFE(eid,heap));
        ngstd::Array<int> dofs;space->GetDofNrs(eid,dofs);if(dofs.Size()!=fe.GetNDof())throw std::runtime_error("Material local DOF mismatch");
        for(auto d:dofs)if(d<0||std::size_t(d)>=nd)throw std::invalid_argument("Constrained material DOF unsupported");
        std::vector<int> uniqueDofs(dofs.begin(),dofs.end());std::sort(uniqueDofs.begin(),uniqueDofs.end());
        if(std::adjacent_find(uniqueDofs.begin(),uniqueDofs.end())!=uniqueDofs.end())throw std::runtime_error("Duplicate material local DOF unsupported");
        ngbla::Matrix<> shape(fe.GetNDof(),3);
        ngstd::Array<int> sampleDofs;materialSpace->GetDofNrs(eid,sampleDofs);if(sampleDofs.Size()!=rule.Size())throw std::runtime_error("Material rule point count mismatch");
        // Match NGSolve's pyramid mass-assembly mapping for parity, not greater
        // geometric accuracy: in 2607/2608 this rule path underestimates the
        // analytic straight-pyramid volume by ~1e-10; point mapping is closer.
        // Do not insert a Radia epsilon correction or change tet/hex mapping.
        const ngfem::BaseMappedIntegrationRule* pyramidRule=
            mesh->GetElType(eid)==ngfem::ET_PYRAMID?&trafo(rule,heap):nullptr;
        int localPoint=0;
        for(const auto& ip:rule){const auto index=sampleDofs[localPoint++];if(index<0)throw std::runtime_error("Invalid material sample index");
            const auto pos=std::lower_bound(globalPoints.begin(),globalPoints.end(),std::size_t(index));
            if(pos==globalPoints.end()||*pos!=std::size_t(index))throw std::runtime_error("Unmapped material sample index");const auto q=std::size_t(pos-globalPoints.begin());
            auto publishPoint=[&](const ngfem::MappedIntegrationPoint<3,3>& mip){fe.CalcMappedShape(mip,shape);
            const double w=mip.GetWeight();if(!std::isfinite(w)||w<=0)throw std::runtime_error("Invalid material measure");
            weights[q]=w;elements[q]=double(e);for(int c=0;c<3;++c){points[3*q+c]=mip.GetPoint()[c];for(int i=0;i<dofs.Size();++i){const double v=shape(i,c);if(!std::isfinite(v))throw std::runtime_error("Nonfinite material shape");if(v!=0.)columns[dofs[i]].emplace_back(std::size_t(c*nq+q),v);}}
            };
            if(pyramidRule)publishPoint(static_cast<const ngfem::MappedIntegrationPoint<3,3>&>((*pyramidRule)[localPoint-1]));
            else{ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);publishPoint(mip);}
        }
    }
    std::size_t nonzeros=0;
    for(auto& column:columns){
        std::sort(column.begin(),column.end(),[](const auto& a,const auto& b){return a.first<b.first;});
        if(std::adjacent_find(column.begin(),column.end(),[](const auto& a,const auto& b){return a.first==b.first;})!=column.end())throw std::runtime_error("Duplicate material evaluation entry");
        if(column.size()>entries-nonzeros)throw std::runtime_error("Material evaluation entry budget exceeded");nonzeros+=column.size();
    }
    MaterialQuadraturePayload out;
    out.dofs=nd;out.totalPoints=totalPoints;out.totalElements=ne;
    out.firstElement=first;out.elementCount=count;out.integrationOrder=degree;
    out.potentialShapeEntries=entries;out.tangentNonzeroBound=tangentBound;
    out.columnOffsets.resize(nd+1);out.rowIndices.reserve(nonzeros);out.values.reserve(nonzeros);
    for(std::size_t col=0;col<nd;++col){
        out.columnOffsets[col]=out.values.size();
        for(const auto& entry:columns[col]){out.rowIndices.push_back(entry.first);out.values.push_back(entry.second);}
    }
    out.columnOffsets[nd]=out.values.size();
    out.globalPointIndices=std::move(globalPoints);out.weights=std::move(weights);
    out.points=std::move(points);out.elements=std::move(elements);
    return out;
}
}
