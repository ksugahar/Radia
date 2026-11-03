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
    // curveOrder 2 applies NGSolve's isoparametric P2 curving to a linear checked
    // mesh, exactly as the Python FE side does with Mesh.Curve(2). A loaded P2 .vol
    // stays as stored: it has no CAD callback, so a second Curve would replace it.
    MaterialStagingContext(const std::string& path,int order,int curveOrder=1){
        mesh=std::make_shared<ngcomp::MeshAccess>(path);ne=std::size_t(mesh->GetNE(ngfem::VOL));
        // int32 element, DOF and point numbering is the structural bound; device
        // budgets are checked where the tables are consumed.
        constexpr std::size_t count=std::numeric_limits<int>::max();
        if(mesh->GetDimension()!=3||!ne||ne>count||mesh->GetDeformation()||mesh->GetCurveOrder()>curveOrder||order<0||order>2||(curveOrder!=1&&curveOrder!=2))
            throw std::invalid_argument("Material quadrature requires bounded RT0/P1/P2 tetrahedral, hexahedral or pyramidal HDiv on an affine or P2-curved mesh");
        if(curveOrder==2&&mesh->GetCurveOrder()!=2){mesh->Curve(2);if(mesh->GetCurveOrder()!=2)throw std::runtime_error("NGSolve did not curve the mesh to order 2");}
        ngcore::Flags flags;flags.SetFlag("order",order);space=std::make_shared<ngcomp::HDivHighOrderFESpace>(mesh,flags);space->Update();space->FinalizeUpdate();
        nd=std::size_t(space->GetNDof());if(!nd||nd>count)throw std::invalid_argument("Material DOF numbering exceeds int32");
        degree=std::max(6,2*order+2);ngcore::Flags irFlags;irFlags.SetFlag("order",degree);
        materialSpace=ngcomp::CreateFESpace("irspace",mesh,irFlags);if(!materialSpace)throw std::runtime_error("NGSolve integration space unavailable");materialSpace->Update();materialSpace->FinalizeUpdate();
        // NGSolve v6.2.2608 IntegrationRuleSpace uses IntegrationRule(type,2*order).
        // https://github.com/NGSolve/ngsolve/blob/v6.2.2608/comp/irspace.cpp
        // Retain the exported factory; its constructor is not a Windows DLL export.
        totalPoints=std::size_t(materialSpace->GetNDof());
        if(!totalPoints||totalPoints>count)throw std::invalid_argument("Material point numbering exceeds int32");
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
namespace p2_staging {
// Independent host form of the device P2 TET map (corners (0,0,0),(1,0,0),
// (0,1,0),(0,0,1), then edges 01,12,20,03,13,23): X(xi) and Jac[k][c]=dX_k/dxi_c.
constexpr double tetNodes[10][3]={{0,0,0},{1,0,0},{0,1,0},{0,0,1},{.5,0,0},{.5,.5,0},{0,.5,0},{0,0,.5},{.5,0,.5},{0,.5,.5}};
constexpr double triNodes[6][2]={{0,0},{1,0},{0,1},{.5,0},{.5,.5},{0,.5}};
inline double tetMap(const double* nodes,const double* xi,double X[3],double J[3][3]){
    const double L[4]={1-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]},dL[4][3]={{-1,-1,-1},{1,0,0},{0,1,0},{0,0,1}};
    constexpr int edge[6][2]={{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};double N[10],dN[10][3];
    for(int i=0;i<4;++i){N[i]=L[i]*(2*L[i]-1);for(int c=0;c<3;++c)dN[i][c]=(4*L[i]-1)*dL[i][c];}
    for(int k=0;k<6;++k){const int a=edge[k][0],b=edge[k][1];N[4+k]=4*L[a]*L[b];for(int c=0;c<3;++c)dN[4+k][c]=4*(dL[a][c]*L[b]+L[a]*dL[b][c]);}
    for(int k=0;k<3;++k){X[k]=0;for(int c=0;c<3;++c)J[k][c]=0;}
    for(int i=0;i<10;++i)for(int k=0;k<3;++k){X[k]+=N[i]*nodes[3*i+k];for(int c=0;c<3;++c)J[k][c]+=dN[i][c]*nodes[3*i+k];}
    return J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
}
// Physical P2 nodes of a volume element through NGSolve's element map.
inline void tetNodesOf(ngfem::ElementTransformation& trafo,std::vector<double>& out){
    for(const auto& n:tetNodes){ngfem::IntegrationPoint ip(n[0],n[1],n[2]);ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);
        for(int d=0;d<3;++d){const double v=mip.GetPoint()(d);if(!std::isfinite(v))throw std::runtime_error("Nonfinite P2 node");out.push_back(v);}}
}
// The P2 map must reproduce NGSolve's curved map (point and Jacobian) at
// points no rule symmetry aliases; it is exact for an order-2 curved mesh.
inline void checkTetMap(ngfem::ElementTransformation& trafo,const double* nodes){
    double size=0;for(int i=1;i<10;++i)for(int d=0;d<3;++d)size=std::max(size,std::abs(nodes[3*i+d]-nodes[d]));
    for(const auto& ip:{ngfem::IntegrationPoint(.137,.211,.303),ngfem::IntegrationPoint(.052,.618,.171)}){
        ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);const double xi[3]={ip(0),ip(1),ip(2)};double X[3],J[3][3];const double det=tetMap(nodes,xi,X,J);
        for(int c=0;c<3;++c)if(std::abs(X[c]-mip.GetPoint()(c))>1e-12*std::max(1.,size))throw std::runtime_error("NGSolve curved TET map differs from the device P2 convention");
        if(!(std::abs(det)>0)||std::abs(det-mip.GetJacobiDet())>1e-10*std::abs(det))throw std::runtime_error("NGSolve curved TET Jacobian differs from the device P2 convention");
    }
}
}
// TET evaluation tables for device E construction. NGSolve owns the
// reference HDiv shapes; they are tabulated once per vertex-rank class (the
// global vertex order that orients NGSolve's shapes) and every element is
// checked against its class before anything is published. The device maps
// shape(i,c) = J(c,:) shapeRef(i,:) / det J with J = [V0-V3, V1-V3, V2-V3],
// NGSolve's TET reference (1,0,0),(0,1,0),(0,0,1),(0,0,0); on a P2-curved mesh
// J is the P2 map Jacobian at each rule point and the element geometry is its
// 10 nodes. The convention is checked per element against the NGSolve map.
struct FlatEvaluationTables {
    std::size_t dofs=0,totalPoints=0,elements=0,local=0,rulePoints=0,nodes=4;
    int integrationOrder=0;
    std::vector<double> rule;          // rulePoints x (xi0,xi1,xi2,weight)
    std::vector<double> shapes;        // classes x rulePoints x local x 3 reference shapes
    std::vector<double> vertices;      // elements x nodes x 3: NGSolve vertex order, or curved P2 nodes
    std::vector<int> elementClass,elementDofs,elementPoints; // elements; elements x local; elements x rulePoints
};
inline FlatEvaluationTables stageFlatEvaluationTables(MaterialStagingContext& context){
    const auto& mesh=context.mesh;const auto& space=context.space;const auto& materialSpace=context.materialSpace;
    const bool curved=mesh->GetCurveOrder()==2;
    FlatEvaluationTables out;out.dofs=context.nd;out.totalPoints=context.totalPoints;out.elements=context.ne;out.integrationOrder=context.degree;out.nodes=curved?10:4;
    const ngfem::IntegrationRule rule(ngfem::ET_TET,2*context.degree);out.rulePoints=std::size_t(rule.Size());
    for(const auto& ip:rule){for(int d=0;d<3;++d)out.rule.push_back(ip(d));out.rule.push_back(ip.Weight());}
    // Asymmetric probe points: class membership and the affine convention are
    // checked at points no rule symmetry can alias.
    const ngfem::IntegrationPoint probes[2]={ngfem::IntegrationPoint(.137,.211,.303),ngfem::IntegrationPoint(.052,.618,.171)};
    std::vector<int> classKey;std::vector<double> classProbe;std::vector<int> covered(context.totalPoints,0);
    ngstd::LocalHeap heap(1<<24,"radia_hdiv_flat_evaluation");
    for(std::size_t e=0;e<context.ne;++e){ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::VOL,int(e));
        if(mesh->GetElType(eid)!=ngfem::ET_TET)throw std::invalid_argument("Device evaluation tables support tetrahedra only");
        const auto& fe=dynamic_cast<const ngfem::HDivFiniteElement<3>&>(space->GetFE(eid,heap));
        ngstd::Array<int> dofs;space->GetDofNrs(eid,dofs);
        if(!out.local)out.local=std::size_t(fe.GetNDof());
        if(std::size_t(fe.GetNDof())!=out.local||std::size_t(dofs.Size())!=out.local||out.local>64)
            throw std::invalid_argument("Device evaluation tables require one local HDiv size <= 64");
        for(auto d:dofs)if(d<0||std::size_t(d)>=context.nd)throw std::invalid_argument("Constrained material DOF unsupported");
        std::vector<int> sorted(dofs.begin(),dofs.end());std::sort(sorted.begin(),sorted.end());
        if(std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end())throw std::runtime_error("Duplicate material local DOF unsupported");
        out.elementDofs.insert(out.elementDofs.end(),dofs.begin(),dofs.end());
        ngstd::Array<int> samples;materialSpace->GetDofNrs(eid,samples);
        if(std::size_t(samples.Size())!=out.rulePoints)throw std::runtime_error("Material rule point count mismatch");
        for(auto q:samples){if(q<0||std::size_t(q)>=context.totalPoints||covered[q]++)throw std::runtime_error("Material samples are not element-disjoint");
            out.elementPoints.push_back(q);}
        const auto vertexNumbers=mesh->GetElement(eid).Vertices();
        int key=0;for(int v=0;v<4;++v){int rank=0;for(int w=0;w<4;++w)rank+=vertexNumbers[w]<vertexNumbers[v];key=4*key+rank;}
        auto& trafo=mesh->GetTrafo(eid,heap);
        if(curved){const auto at=out.vertices.size();p2_staging::tetNodesOf(trafo,out.vertices);p2_staging::checkTetMap(trafo,out.vertices.data()+at);}
        else{double V[4][3];
            for(int v=0;v<4;++v){const auto p=mesh->GetPoint<3>(vertexNumbers[v]);for(int d=0;d<3;++d){V[v][d]=p(d);out.vertices.push_back(p(d));}}
            double J[3][3];for(int c=0;c<3;++c)for(int d=0;d<3;++d)J[c][d]=V[d][c]-V[3][c];
            const double det=J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
            double size=0;for(int c=0;c<3;++c)for(int d=0;d<3;++d)size=std::max(size,std::abs(J[c][d]));
            for(const auto& ip:probes){ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);
                for(int c=0;c<3;++c){double x=V[3][c];for(int d=0;d<3;++d)x+=J[c][d]*ip(d);
                    if(std::abs(x-mip.GetPoint()(c))>1e-12*std::max(1.,size))throw std::runtime_error("NGSolve TET map differs from the device affine convention");}
                if(!(std::abs(det)>0)||std::abs(det-mip.GetJacobiDet())>1e-12*std::abs(det))throw std::runtime_error("NGSolve TET Jacobian differs from the device affine convention");}}
        ngbla::Matrix<> shape(out.local,3);std::vector<double> probe;
        for(const auto& ip:probes){fe.CalcShape(ip,shape);for(std::size_t i=0;i<out.local;++i)for(int c=0;c<3;++c)probe.push_back(shape(i,c));}
        auto found=std::find(classKey.begin(),classKey.end(),key);std::size_t cls=std::size_t(found-classKey.begin());
        if(found==classKey.end()){
            classKey.push_back(key);classProbe.insert(classProbe.end(),probe.begin(),probe.end());
            for(const auto& ip:rule){fe.CalcShape(ip,shape);for(std::size_t i=0;i<out.local;++i)for(int c=0;c<3;++c){
                if(!std::isfinite(shape(i,c)))throw std::runtime_error("Nonfinite reference HDiv shape");out.shapes.push_back(shape(i,c));}}
        }else{
            const double* reference=classProbe.data()+cls*probe.size();
            for(std::size_t i=0;i<probe.size();++i)if(std::abs(probe[i]-reference[i])>1e-12*(1+std::abs(reference[i])))
                throw std::runtime_error("NGSolve HDiv reference shapes are not determined by the element vertex order");
        }
        out.elementClass.push_back(int(cls));
    }
    for(auto c:covered)if(c!=1)throw std::runtime_error("Material global sample coverage mismatch");
    return out;
}
}
