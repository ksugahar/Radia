// Bounded affine-tetra RT0 setup. NGSolve owns mapped basis/orientation.
#pragma once
#include <hdivfe.hpp>
#include <hdivhofespace.hpp>
#include <meshaccess.hpp>
#include <map>
#include <array>
#include <algorithm>
namespace radia::rt0_cuda {
using Vec=ngbla::Vec<3>;
using Key=std::array<int,3>;
struct Facet {int element,opposite,count=1;};
inline mxArray* matrix(std::size_t rows,std::size_t cols,const std::vector<double>& values){
    auto* a=mxCreateDoubleMatrix(rows,cols,mxREAL);
    for(std::size_t i=0;i<rows;++i)for(std::size_t j=0;j<cols;++j)mxGetDoubles(a)[i+rows*j]=values[i*cols+j];return a;
}
inline mxArray* assemble(const std::string& path,const std::array<double,3>& field){
    auto mesh=std::make_shared<ngcomp::MeshAccess>(path);
    const auto ne=mesh->GetNE(ngfem::VOL);
    if(mesh->GetDimension()!=3||!ne||ne>16||mesh->GetDeformation())throw std::invalid_argument("RT0 setup requires 1..16 undeformed affine tetrahedra");
    ngcore::Flags flags;flags.SetFlag("order",0);
    auto space=std::make_shared<ngcomp::HDivHighOrderFESpace>(mesh,flags);space->Update();space->FinalizeUpdate();
    const auto nd=space->GetNDof();if(!nd||nd>256)throw std::runtime_error("RT0 setup DOF budget exceeded");
    const std::array<std::array<double,3>,4> refs={{{1,0,0},{0,1,0},{0,0,1},{0,0,0}}};
    std::vector<double> cells,faces,mass(nd*nd,0.),rhs(nd,0.),mean(3*nd,0.);double volume=0.;
    std::map<Key,Facet> topology;std::vector<std::array<Vec,4>> vertices;
    std::vector<std::array<int,4>> ids;
    ngstd::LocalHeap heap(1<<22,"radia_cuda_rt0_setup");
    for(std::size_t e=0;e<ne;++e){
        ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::VOL,int(e));
        if(mesh->GetElType(eid)!=ngfem::ET_TET||mesh->GetElement(eid).is_curved)throw std::invalid_argument("RT0 setup rejects non-affine or non-tetra elements");
        auto& trafo=mesh->GetTrafo(eid,heap);const auto& fe=dynamic_cast<const ngfem::HDivFiniteElement<3>&>(space->GetFE(eid,heap));
        ngstd::Array<int> dofs;space->GetDofNrs(eid,dofs);if(dofs.Size()!=4||fe.GetNDof()!=4)throw std::runtime_error("Unexpected RT0 local DOFs");
        std::array<Vec,4> points;std::array<int,4> numbers;auto raw=mesh->GetElVertices(eid);
        for(int i=0;i<4;++i){
            ngfem::IntegrationPoint ip(refs[i][0],refs[i][1],refs[i][2]);ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);points[i]=mip.GetPoint();
            double best=1.e100;int match=-1;
            for(auto id:raw){auto point=mesh->GetPoint<3>(id);double distance=ngbla::L2Norm(point-points[i]);if(distance<best){best=distance;match=id;}}
            if(match<0||best>1.e-10*std::max(1.,ngbla::L2Norm(points[i])))throw std::runtime_error("Mapped tetra vertex mismatch");numbers[i]=match;
            for(int c=0;c<3;++c)cells.push_back(points[i][c]);
        }
        vertices.push_back(points);ids.push_back(numbers);
        for(int opposite=0;opposite<4;++opposite){Key key;int k=0;for(int i=0;i<4;++i)if(i!=opposite)key[k++]=numbers[i];std::sort(key.begin(),key.end());auto found=topology.find(key);if(found==topology.end())topology.emplace(key,Facet{int(e),opposite});else if(++found->second.count>2)throw std::runtime_error("Nonmanifold RT0 facet");}
        ngfem::IntegrationRule rule(ngfem::ET_TET,2);
        for(const auto& ip:rule){
            ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);ngbla::Matrix<> shape(4,3);fe.CalcMappedShape(mip,shape);double w=mip.GetWeight();volume+=w;
            for(int i=0;i<4;++i){if(dofs[i]<0||std::size_t(dofs[i])>=nd)throw std::runtime_error("Constrained RT0 DOF unsupported");
                for(int c=0;c<3;++c){rhs[dofs[i]]+=w*field[c]*shape(i,c);mean[c*nd+dofs[i]]+=w*shape(i,c);}
                for(int j=0;j<4;++j)for(int c=0;c<3;++c)mass[dofs[i]*nd+dofs[j]]+=w*shape(i,c)*shape(j,c);
            }
        }
    }
    std::map<Key,int> skin;
    for(std::size_t f=0;f<mesh->GetNE(ngfem::BND);++f){
        auto raw=mesh->GetElVertices(ngfem::ElementId(ngfem::BND,int(f)));if(raw.Size()!=3)throw std::invalid_argument("RT0 skin must be triangles");Key key{int(raw[0]),int(raw[1]),int(raw[2])};std::sort(key.begin(),key.end());if(++skin[key]!=1)throw std::runtime_error("Duplicate mesh skin facet");
        if(!topology.count(key))throw std::runtime_error("Mesh skin does not match tetra topology");
    }
    std::size_t nf=0;for(const auto& [key,f]:topology)if(f.count==1){if(!skin.count(key))throw std::runtime_error("Missing exterior mesh skin");++nf;}
    if(ne+nf>128||volume<=0)throw std::runtime_error("RT0 geometry budget/volume invalid");
    std::vector<double> B((ne+nf)*nd,0.);
    for(std::size_t e=0;e<ne;++e){
        ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::VOL,int(e));auto& trafo=mesh->GetTrafo(eid,heap);const auto& fe=dynamic_cast<const ngfem::HDivFiniteElement<3>&>(space->GetFE(eid,heap));ngstd::Array<int> dofs;space->GetDofNrs(eid,dofs);
        ngfem::IntegrationPoint ip(.25,.25,.25);ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);ngbla::Vector<> div(4);fe.CalcMappedDivShape(mip,div);for(int i=0;i<4;++i)B[e*nd+dofs[i]]=-div[i];
    }
    std::size_t frow=ne;
    for(const auto& [key,f]:topology)if(f.count==1){
        ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::VOL,f.element);auto& trafo=mesh->GetTrafo(eid,heap);const auto& fe=dynamic_cast<const ngfem::HDivFiniteElement<3>&>(space->GetFE(eid,heap));ngstd::Array<int> dofs;space->GetDofNrs(eid,dofs);
        std::array<int,3> local;int k=0;Vec ref(0.);for(int i=0;i<4;++i)if(i!=f.opposite){local[k++]=i;for(int c=0;c<3;++c)ref[c]+=refs[i][c]/3.;}
        const auto& points=vertices[f.element];Vec u=points[local[1]]-points[local[0]],v=points[local[2]]-points[local[0]];Vec normal=ngbla::Cross(u,v);double norm=ngbla::L2Norm(normal);if(norm<=0)throw std::runtime_error("Degenerate RT0 skin");normal/=norm;
        if(ngbla::InnerProduct(normal,points[f.opposite]-points[local[0]])>0)normal*=-1.;
        ngfem::IntegrationPoint ip(ref[0],ref[1],ref[2]);ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);ngbla::Matrix<> shape(4,3);fe.CalcMappedShape(mip,shape);
        for(int i=0;i<4;++i)for(int c=0;c<3;++c)B[frow*nd+dofs[i]]+=shape(i,c)*normal[c];
        for(int i:local)for(int c=0;c<3;++c)faces.push_back(points[i][c]);++frow;
    }
    for(auto& v:mean)v/=volume;
    const char* names[]={"cells","faces","charge_map","mass","rhs","mean_operator"};auto* out=mxCreateStructMatrix(1,1,6,names);
    mxSetField(out,0,"cells",matrix(cells.size(),1,cells));mxSetField(out,0,"faces",matrix(faces.size(),1,faces));mxSetField(out,0,"charge_map",matrix(ne+nf,nd,B));mxSetField(out,0,"mass",matrix(nd,nd,mass));mxSetField(out,0,"rhs",matrix(nd,1,rhs));mxSetField(out,0,"mean_operator",matrix(3,nd,mean));return out;
}
}
