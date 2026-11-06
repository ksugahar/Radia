// NGSolve-owned reference tables for device construction of the TET charge map B.
#pragma once
#include "hdiv_material_staging.hpp"
#include <algorithm>
#include <array>
#include <map>
namespace radia::hdiv_cuda {
// Row order matches the CPU charge basis: cell monomials (cell-major), then
// monomials of the exterior faces in boundary-element order. No B value is
// produced here.
// Affine: device rows are T_cell[class]/det J on the cell DOFs and
// T_face[class]/|(F1-F0)x(F2-F0)| on the face DOFs; tables are L2 projections
// onto the Gram monomials in the cell_verts/face_verts frames.
// Curved (P2): monomials are NGSolve reference monomials (cells in xi, faces
// in the boundary element's own eta, as the CPU curved charge basis); tables
// are the reference moments R and the device solves G_e X = R with G_e the
// curved monomial Gram on NGSolve's default mass rules (cellRule/faceRule).
// The geometry is the 10/6 P2 nodes and the mesh vertex numbers.
struct FlatChargeTables {
    std::size_t dofs=0,cells=0,faces=0,cellMonomials=0,faceMonomials=0,cellWidth=0,faceWidth=0,cellClasses=0,faceClasses=0;
    bool curved=false,hex=false,mixed=false;
    // Mixed mesh: tables per class and hosts in order, cells then exterior faces; cellVertices and faceVertices hold the unit-cube lattices.
    std::vector<int> hostTypes,hostClass,classInfo,hostDofs,chargeHost,chargeKind,chargeExponents;std::vector<double> classTables;
    std::vector<int> cellClass,cellDofs,faceClass,faceDofs,cellExponents,faceExponents,cellNumbers,faceNumbers;
    std::vector<double> cellTables,faceTables,cellVertices,faceVertices,cellRule,faceRule;
};
namespace charge_staging {
inline int rankKey(const int* numbers,int count){
    int key=0;for(int v=0;v<count;++v){int rank=0;for(int w=0;w<count;++w)rank+=numbers[w]<numbers[v];key=count*key+rank;}return key;
}
// Solves the small symmetric positive monomial Gram M X = R in place (R is n x k, row-major).
inline void solveSmall(std::vector<double> M,std::vector<double>& R,std::size_t n,std::size_t k){
    for(std::size_t c=0;c<n;++c){std::size_t p=c;for(std::size_t r=c+1;r<n;++r)if(std::abs(M[r*n+c])>std::abs(M[p*n+c]))p=r;
        if(!(std::abs(M[p*n+c])>0))throw std::runtime_error("Singular charge monomial Gram");
        if(p!=c){for(std::size_t j=0;j<n;++j)std::swap(M[c*n+j],M[p*n+j]);for(std::size_t j=0;j<k;++j)std::swap(R[c*k+j],R[p*k+j]);}
        for(std::size_t r=0;r<n;++r){if(r==c)continue;const double f=M[r*n+c]/M[c*n+c];if(f==0.)continue;
            for(std::size_t j=c;j<n;++j)M[r*n+j]-=f*M[c*n+j];for(std::size_t j=0;j<k;++j)R[r*k+j]-=f*R[c*k+j];}}
    for(std::size_t c=0;c<n;++c)for(std::size_t j=0;j<k;++j)R[c*k+j]/=M[c*n+c];
}
// NGSolve TET reference point of cell vertex v (mesh vertex order).
constexpr double cellVertexXi[4][3]={{1,0,0},{0,1,0},{0,0,1},{0,0,0}};
// NGSolve's default orders of the CPU curved charge forms on Curve(2), pinned
// by tests/test_hdiv_device_curved_rule_contract.py: the mass rules
// (mv = 2 pv on cells, mb = 2 p on faces) are inexact on curved geometry and
// fix the projection; the moment rules (bv, bb) integrate polynomials exactly.
inline int curvedCellMassOrder(int order){return 2*std::max(order-1,0);}
inline int curvedFaceMassOrder(int order){return 2*order;}
inline int curvedCellMomentOrder(int order){return std::max(order+std::max(order-1,0)-1,0);}
inline int curvedFaceMomentOrder(int order){return 2*order+1;}
// Independent host form of the device P2 triangle map (corners (0,0),(1,0),(0,1), then edges 01,12,20).
inline void triMap(const double* nodes,double e0,double e1,double X[3],double& measure){
    const double L[3]={1-e0-e1,e0,e1},dL[3][2]={{-1,-1},{1,0},{0,1}};constexpr int edge[3][2]={{0,1},{1,2},{2,0}};
    double N[6],dN[6][2];
    for(int i=0;i<3;++i){N[i]=L[i]*(2*L[i]-1);for(int c=0;c<2;++c)dN[i][c]=(4*L[i]-1)*dL[i][c];}
    for(int k=0;k<3;++k){const int a=edge[k][0],b=edge[k][1];N[3+k]=4*L[a]*L[b];for(int c=0;c<2;++c)dN[3+k][c]=4*(dL[a][c]*L[b]+L[a]*dL[b][c]);}
    double T[2][3]={};for(int k=0;k<3;++k)X[k]=0;
    for(int i=0;i<6;++i)for(int k=0;k<3;++k){X[k]+=N[i]*nodes[3*i+k];for(int c=0;c<2;++c)T[c][k]+=dN[i][c]*nodes[3*i+k];}
    const double n[3]={T[0][1]*T[1][2]-T[0][2]*T[1][1],T[0][2]*T[1][0]-T[0][0]*T[1][2],T[0][0]*T[1][1]-T[0][1]*T[1][0]};
    measure=std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
}
inline void triNodesOf(ngfem::ElementTransformation& trafo,std::vector<double>& out){
    for(const auto& n:p2_staging::triNodes){ngfem::IntegrationPoint ip(n[0],n[1]);ngfem::MappedIntegrationPoint<2,3> mip(ip,trafo);
        for(int d=0;d<3;++d){const double v=mip.GetPoint()(d);if(!std::isfinite(v))throw std::runtime_error("Nonfinite P2 face node");out.push_back(v);}}
}
inline void checkTriMap(ngfem::ElementTransformation& trafo,const double* nodes){
    double size=0;for(int i=1;i<6;++i)for(int d=0;d<3;++d)size=std::max(size,std::abs(nodes[3*i+d]-nodes[d]));
    for(const auto& ip:{ngfem::IntegrationPoint(.137,.211),ngfem::IntegrationPoint(.618,.171)}){
        ngfem::MappedIntegrationPoint<2,3> mip(ip,trafo);double X[3],measure;triMap(nodes,ip(0),ip(1),X,measure);
        for(int c=0;c<3;++c)if(std::abs(X[c]-mip.GetPoint()(c))>1e-12*std::max(1.,size))throw std::runtime_error("NGSolve curved face map differs from the device P2 convention");
        if(!(measure>0)||std::abs(measure-mip.GetMeasure())>1e-10*measure)throw std::runtime_error("NGSolve curved face measure differs from the device P2 convention");
    }
}
inline void ruleRows(const ngfem::IntegrationRule& rule,int dimension,std::vector<double>& out){
    for(const auto& ip:rule){for(int d=0;d<3;++d)out.push_back(d<dimension?ip(d):0.);out.push_back(ip.Weight());}
}
}
inline FlatChargeTables stageFlatChargeTables(MaterialStagingContext& context){
    using namespace charge_staging;
    const auto& mesh=context.mesh;const auto& space=context.space;const int order=space->GetOrder();
    const bool curved=mesh->GetCurveOrder()==2;
    if(mesh->GetCurveOrder()>2||order<0||order>2||(curved&&order<1))
        throw std::invalid_argument("Device charge tables require affine TET RT0/P1/P2 HDiv or P2-curved TET P1/P2 HDiv");
    FlatChargeTables out;out.dofs=context.nd;out.cells=context.ne;out.curved=curved;
    // Cell charge degree max(order-1,0) and face degree order, as the CPU charge basis.
    const int pv=std::max(order-1,0);
    for(int i=0;i<=pv;++i)for(int j=0;j<=pv-i;++j)for(int k=0;k<=pv-i-j;++k){out.cellExponents.insert(out.cellExponents.end(),{i,j,k});}
    for(int i=0;i<=order;++i)for(int j=0;j<=order-i;++j){out.faceExponents.insert(out.faceExponents.end(),{i,j});}
    const std::size_t cm=out.cellExponents.size()/3,fm=out.faceExponents.size()/2;out.cellMonomials=cm;out.faceMonomials=fm;
    // Affine: exact for the monomial Gram (degree 2 order) and the projected div/trace.
    // Curved: NGSolve's moment orders; the curved Gram rules are published for the device.
    const ngfem::IntegrationRule tetRule(ngfem::ET_TET,curved?curvedCellMomentOrder(order):2*order+2),
        triRule(ngfem::ET_TRIG,curved?curvedFaceMomentOrder(order):2*order+2);
    if(curved){ruleRows(ngfem::IntegrationRule(ngfem::ET_TET,curvedCellMassOrder(order)),3,out.cellRule);
        ruleRows(ngfem::IntegrationRule(ngfem::ET_TRIG,curvedFaceMassOrder(order)),2,out.faceRule);}
    auto cellMonomial=[&](const double* xi,double* m){const double g[3]={xi[1],xi[2],1-xi[0]-xi[1]-xi[2]};const double* c=curved?xi:g;
        for(std::size_t a=0;a<cm;++a){const int* e=&out.cellExponents[3*a];m[a]=std::pow(c[0],e[0])*std::pow(c[1],e[1])*std::pow(c[2],e[2]);}};
    auto faceMonomial=[&](double g0,double g1,double* m){
        for(std::size_t a=0;a<fm;++a){const int* e=&out.faceExponents[2*a];m[a]=std::pow(g0,e[0])*std::pow(g1,e[1]);}};
    const ngfem::IntegrationPoint probes[2]={ngfem::IntegrationPoint(.137,.211,.303),ngfem::IntegrationPoint(.052,.618,.171)};
    ngstd::LocalHeap heap(1<<24,"radia_hdiv_flat_charge");
    std::vector<int> classKey;std::vector<double> classProbe;std::vector<int> cellVertexNumbers(4*context.ne);
    std::size_t local=0;
    for(std::size_t e=0;e<context.ne;++e){ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::VOL,int(e));
        if(mesh->GetElType(eid)!=ngfem::ET_TET)throw std::invalid_argument("Device charge tables support tetrahedra only");
        const auto& fe=dynamic_cast<const ngfem::HDivFiniteElement<3>&>(space->GetFE(eid,heap));
        ngstd::Array<int> dofs;space->GetDofNrs(eid,dofs);
        if(!local)local=std::size_t(fe.GetNDof());
        if(std::size_t(fe.GetNDof())!=local||std::size_t(dofs.Size())!=local||local>64)
            throw std::invalid_argument("Device charge tables require one local HDiv size <= 64");
        for(auto d:dofs)if(d<0||std::size_t(d)>=context.nd)throw std::invalid_argument("Constrained charge DOF unsupported");
        std::vector<int> sorted(dofs.begin(),dofs.end());std::sort(sorted.begin(),sorted.end());
        if(std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end())throw std::runtime_error("Duplicate charge local DOF unsupported");
        out.cellDofs.insert(out.cellDofs.end(),dofs.begin(),dofs.end());
        const auto vertexNumbers=mesh->GetElement(eid).Vertices();int numbers[4];
        for(int v=0;v<4;++v){numbers[v]=vertexNumbers[v];cellVertexNumbers[4*e+v]=numbers[v];}
        const int key=rankKey(numbers,4);
        auto& trafo=mesh->GetTrafo(eid,heap);
        if(curved){const auto at=out.cellVertices.size();p2_staging::tetNodesOf(trafo,out.cellVertices);p2_staging::checkTetMap(trafo,out.cellVertices.data()+at);
            out.cellNumbers.insert(out.cellNumbers.end(),numbers,numbers+4);}
        else{double V[4][3];
            for(int v=0;v<4;++v){const auto p=mesh->GetPoint<3>(vertexNumbers[v]);for(int d=0;d<3;++d){V[v][d]=p(d);out.cellVertices.push_back(p(d));}}
            double J[3][3];for(int c=0;c<3;++c)for(int d=0;d<3;++d)J[c][d]=V[d][c]-V[3][c];
            const double det=J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
            double size=0;for(int c=0;c<3;++c)for(int d=0;d<3;++d)size=std::max(size,std::abs(J[c][d]));
            for(const auto& ip:probes){ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);
                for(int c=0;c<3;++c){double x=V[3][c];for(int d=0;d<3;++d)x+=J[c][d]*ip(d);
                    if(std::abs(x-mip.GetPoint()(c))>1e-12*std::max(1.,size))throw std::runtime_error("NGSolve TET map differs from the device affine convention");}
                if(!(std::abs(det)>0)||std::abs(det-mip.GetJacobiDet())>1e-12*std::abs(det))throw std::runtime_error("NGSolve TET Jacobian differs from the device affine convention");}}
        ngbla::Vector<> div(local);ngbla::Matrix<> shape(local,3);std::vector<double> probe;
        for(const auto& ip:probes){fe.CalcDivShape(ip,div);fe.CalcShape(ip,shape);
            for(std::size_t i=0;i<local;++i){probe.push_back(div(i));for(int c=0;c<3;++c)probe.push_back(shape(i,c));}}
        auto found=std::find(classKey.begin(),classKey.end(),key);const std::size_t cls=std::size_t(found-classKey.begin());
        if(found==classKey.end()){
            classKey.push_back(key);classProbe.insert(classProbe.end(),probe.begin(),probe.end());
            // R = -int m(xi) div_ref phi(xi). Affine T = Mref^{-1} R is -div u projected in the
            // cell_verts frame; curved keeps R for the per-element device solve.
            std::vector<double> M(cm*cm,0.),R(cm*local,0.),m(cm);
            for(const auto& ip:tetRule){const double xi[3]={ip(0),ip(1),ip(2)};cellMonomial(xi,m.data());fe.CalcDivShape(ip,div);
                for(std::size_t a=0;a<cm;++a){for(std::size_t b=0;b<cm;++b)M[a*cm+b]+=ip.Weight()*m[a]*m[b];
                    for(std::size_t i=0;i<local;++i)R[a*local+i]-=ip.Weight()*m[a]*div(i);}}
            if(!curved)solveSmall(M,R,cm,local);
            for(double v:R)if(!std::isfinite(v))throw std::runtime_error("Nonfinite charge cell table");
            out.cellTables.insert(out.cellTables.end(),R.begin(),R.end());
        }else{
            const double* reference=classProbe.data()+cls*probe.size();
            for(std::size_t i=0;i<probe.size();++i)if(std::abs(probe[i]-reference[i])>1e-12*(1+std::abs(reference[i])))
                throw std::runtime_error("NGSolve HDiv reference shapes are not determined by the element vertex order");
        }
        out.cellClass.push_back(int(cls));
    }
    out.cellWidth=local;out.cellClasses=classKey.size();
    // Exterior faces: a boundary triangle owned by exactly one volume element
    // (sorted corner numbers, as the CPU charge basis); two owners is an
    // internal interface without single-sided charge; anything else is corrupt.
    std::map<std::array<int,3>,std::vector<std::pair<std::size_t,int>>> owners;
    for(std::size_t e=0;e<context.ne;++e)for(int skip=0;skip<4;++skip){std::array<int,3> key;int n=0;
        for(int v=0;v<4;++v)if(v!=skip)key[n++]=cellVertexNumbers[4*e+v];std::sort(key.begin(),key.end());owners[key].push_back({e,skip});}
    std::vector<int> faceKey;
    const std::size_t nse=std::size_t(mesh->GetNSE());
    for(std::size_t s=0;s<nse;++s){ngstd::HeapReset reset(heap);ngfem::ElementId sid(ngfem::BND,int(s));
        if(mesh->GetElType(sid)!=ngfem::ET_TRIG)throw std::invalid_argument("Device charge tables support triangular boundary elements only");
        const auto vertexNumbers=mesh->GetElement(sid).Vertices();int numbers[3];std::array<int,3> key;
        for(int v=0;v<3;++v)key[v]=numbers[v]=vertexNumbers[v];std::sort(key.begin(),key.end());
        const auto owner=owners.find(key);const std::size_t count=owner==owners.end()?0:owner->second.size();
        if(count==2)continue;
        if(count!=1)throw std::runtime_error("Charge boundary element matches "+std::to_string(count)+" volume facets; the mesh is non-conforming or corrupt");
        const std::size_t cell=owner->second[0].first;ngfem::ElementId eid(ngfem::VOL,int(cell));
        const auto& fe=dynamic_cast<const ngfem::HDivFiniteElement<3>&>(space->GetFE(eid,heap));
        int perm[3];for(int v=0;v<3;++v){perm[v]=-1;for(int w=0;w<4;++w)if(cellVertexNumbers[4*cell+w]==numbers[v])perm[v]=w;
            if(perm[v]<0)throw std::runtime_error("Charge face vertex is not a vertex of its cell");}
        ngstd::Array<int> faceDofs;space->GetDofNrs(sid,faceDofs);
        if(!out.faceWidth)out.faceWidth=std::size_t(faceDofs.Size());
        if(std::size_t(faceDofs.Size())!=out.faceWidth||!out.faceWidth||out.faceWidth>local)throw std::invalid_argument("Device charge tables require one face DOF count");
        std::vector<std::size_t> column(out.faceWidth);
        for(std::size_t k=0;k<out.faceWidth;++k){std::size_t i=0;while(i<local&&out.cellDofs[cell*local+i]!=faceDofs[k])++i;
            if(i==local)throw std::runtime_error("Charge face DOF is not a DOF of its cell");column[k]=i;}
        const double* a=cellVertexXi[perm[0]];const double* b=cellVertexXi[perm[1]];const double* c=cellVertexXi[perm[2]];
        // Face parameter xi(g) = o + g0 e0 + g1 e1. Affine: the face_verts frame
        // (o = F0). Curved: the boundary element's NGSolve reference, whose
        // vertices are (1,0),(0,1),(0,0) (o = F2).
        const double* o=curved?c:a;const double* p=curved?a:b;const double* q=curved?b:c;
        const double e0[3]={p[0]-o[0],p[1]-o[1],p[2]-o[2]},e1[3]={q[0]-o[0],q[1]-o[1],q[2]-o[2]};
        const double nu[3]={e0[1]*e1[2]-e0[2]*e1[1],e0[2]*e1[0]-e0[0]*e1[2],e0[0]*e1[1]-e0[1]*e1[0]};
        // R = int m(g) phi_ref(xi(g)).nu (the flux density u.n times the face
        // measure); affine T = Mtri^{-1} R, curved keeps R.
        std::vector<double> M(fm*fm,0.),R(fm*local,0.),m(fm);ngbla::Matrix<> shape(local,3);
        for(const auto& ip:triRule){const double g0=ip(0),g1=ip(1);faceMonomial(g0,g1,m.data());
            ngfem::IntegrationPoint xi(o[0]+g0*e0[0]+g1*e1[0],o[1]+g0*e0[1]+g1*e1[1],o[2]+g0*e0[2]+g1*e1[2]);fe.CalcShape(xi,shape);
            for(std::size_t i=0;i<local;++i){const double flux=shape(i,0)*nu[0]+shape(i,1)*nu[1]+shape(i,2)*nu[2];
                for(std::size_t r=0;r<fm;++r)R[r*local+i]+=ip.Weight()*m[r]*flux;}
            for(std::size_t r=0;r<fm;++r)for(std::size_t t=0;t<fm;++t)M[r*fm+t]+=ip.Weight()*m[r]*m[t];}
        if(!curved)solveSmall(M,R,fm,local);
        double scale=0;for(double v:R){if(!std::isfinite(v))throw std::runtime_error("Nonfinite charge face table");scale=std::max(scale,std::abs(v));}
        std::vector<unsigned char> onFace(local,0);for(auto i:column)onFace[i]=1;
        for(std::size_t r=0;r<fm;++r)for(std::size_t i=0;i<local;++i)if(!onFace[i]&&std::abs(R[r*local+i])>1e-12*scale)
            throw std::runtime_error("HDiv normal trace is not carried by the boundary-element DOFs");
        std::vector<double> table(fm*out.faceWidth);
        for(std::size_t r=0;r<fm;++r)for(std::size_t k=0;k<out.faceWidth;++k)table[r*out.faceWidth+k]=R[r*local+column[k]];
        const int fkey=rankKey(numbers,3);auto found=std::find(faceKey.begin(),faceKey.end(),fkey);const std::size_t cls=std::size_t(found-faceKey.begin());
        if(found==faceKey.end()){faceKey.push_back(fkey);out.faceTables.insert(out.faceTables.end(),table.begin(),table.end());}
        else{const double* reference=out.faceTables.data()+cls*table.size();
            for(std::size_t i=0;i<table.size();++i)if(std::abs(table[i]-reference[i])>1e-12*(1+std::abs(reference[i])))
                throw std::runtime_error("NGSolve HDiv face traces are not determined by the face vertex order");}
        out.faceClass.push_back(int(cls));out.faceDofs.insert(out.faceDofs.end(),faceDofs.begin(),faceDofs.end());
        if(curved){auto& strafo=mesh->GetTrafo(sid,heap);const auto at=out.faceVertices.size();triNodesOf(strafo,out.faceVertices);
            checkTriMap(strafo,out.faceVertices.data()+at);out.faceNumbers.insert(out.faceNumbers.end(),numbers,numbers+3);}
        else for(int v=0;v<3;++v){const auto pt=mesh->GetPoint<3>(vertexNumbers[v]);for(int d=0;d<3;++d)out.faceVertices.push_back(pt(d));}
        ++out.faces;
    }
    out.faceClasses=faceKey.size();
    // Charge rows are int32 indices downstream; the device budget, not a fixed count, bounds the size.
    constexpr std::size_t top=std::numeric_limits<int>::max();
    if(!out.faces||!cm||!fm||out.faces>top/fm||out.cells>(top-out.faces*fm)/cm)throw std::invalid_argument("Device charge map outside the bounded TET scope");
    return out;
}
}
#include "hdiv_hex_staging.hpp"
