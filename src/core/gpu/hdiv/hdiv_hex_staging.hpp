// NGSolve-owned reference tables and Q2 lattices of HEX meshes for device
// construction of the charge map B, the evaluation map E and the charge Gram.
// The CPU hex charge basis integrates the charges in the reference measure
// (the Piola Jacobians cancel up to sgn det J), so B rows depend on the
// geometry only through that sign: per vertex-rank class and orientation,
// cell rows are sgn det J times the reference-cube projection of
// -div_ref phi onto the tensor monomials, face rows the reference-square
// projection of phi_ref . nu over the boundary element's own (u,v), which
// carries no sign. No B, E or Gram value is produced here.
#pragma once
#include <algorithm>
#include <array>
#include <map>
namespace radia::hdiv_cuda::hex_staging {
// NGSolve ET_HEX reference vertices and faces; ET_QUAD vertices (0,0),(1,0),(1,1),(0,1).
constexpr double hexVertexXi[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
constexpr int hexFaces[6][4]={{0,1,2,3},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}};
inline void lag3(double t,double v[3],double d[3]){
    v[0]=2*(t-.5)*(t-1);v[1]=4*t*(1-t);v[2]=2*t*(t-.5);d[0]=4*t-3;d[1]=4-8*t;d[2]=4*t-1;
}
// Independent host form of the device Q2 lattice map: node ix + 3 iy + 9 iz at xi = (ix,iy,iz)/2.
inline double hexMap(const double* nodes,const double* xi,double X[3],double J[3][3]){
    double v[3][3],d[3][3];for(int c=0;c<3;++c)lag3(xi[c],v[c],d[c]);
    for(int k=0;k<3;++k){X[k]=0;for(int c=0;c<3;++c)J[k][c]=0;}
    for(int iz=0;iz<3;++iz)for(int iy=0;iy<3;++iy)for(int ix=0;ix<3;++ix){const double* n=nodes+3*(ix+3*iy+9*iz);
        const double s=v[0][ix]*v[1][iy]*v[2][iz],g[3]={d[0][ix]*v[1][iy]*v[2][iz],v[0][ix]*d[1][iy]*v[2][iz],v[0][ix]*v[1][iy]*d[2][iz]};
        for(int k=0;k<3;++k){X[k]+=s*n[k];for(int c=0;c<3;++c)J[k][c]+=g[c]*n[k];}}
    return J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
}
// Physical Q2 lattice of a volume (27 nodes) or boundary (9 nodes) element through NGSolve's map.
inline void latticeOf(ngfem::ElementTransformation& trafo,int dimension,std::vector<double>& out){
    for(int iz=0;iz<(dimension==3?3:1);++iz)for(int iy=0;iy<3;++iy)for(int ix=0;ix<3;++ix){
        double point[3];
        if(dimension==3){ngfem::IntegrationPoint ip(.5*ix,.5*iy,.5*iz);ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);for(int d=0;d<3;++d)point[d]=mip.GetPoint()(d);}
        else{ngfem::IntegrationPoint ip(.5*ix,.5*iy);ngfem::MappedIntegrationPoint<2,3> mip(ip,trafo);for(int d=0;d<3;++d)point[d]=mip.GetPoint()(d);}
        for(double v:point){if(!std::isfinite(v))throw std::runtime_error("Nonfinite HEX lattice node");out.push_back(v);}}
}
// The Q2 lattice reproduces NGSolve's affine, trilinear or Curve(2) hex map
// (all inside the tensor Q2 span). NGSolve orients hexes either way (its
// structured meshes are left-handed); det J must keep one sign per element,
// which is returned.
inline int checkHexMap(ngfem::ElementTransformation& trafo,const double* nodes){
    int sign=0;double size=0;for(int i=1;i<27;++i)for(int d=0;d<3;++d)size=std::max(size,std::abs(nodes[3*i+d]-nodes[d]));
    for(const auto& ip:{ngfem::IntegrationPoint(.137,.239,.361),ngfem::IntegrationPoint(.713,.211,.557),ngfem::IntegrationPoint(.319,.823,.173)}){
        ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);const double xi[3]={ip(0),ip(1),ip(2)};double X[3],J[3][3];const double det=hexMap(nodes,xi,X,J);
        for(int c=0;c<3;++c)if(std::abs(X[c]-mip.GetPoint()(c))>1e-12*std::max(1.,size))throw std::runtime_error("NGSolve HEX map differs from the device Q2 lattice convention");
        const int s=det>0?1:det<0?-1:0;
        if(!s||(sign&&s!=sign)||!(std::abs(det-mip.GetJacobiDet())<=1e-10*std::abs(det)))
            throw std::runtime_error("Device HEX route requires a Jacobian of one nonzero sign matching NGSolve's map");
        sign=s;
    }
    return sign;
}
inline int rankKey(const int* numbers,int count){
    int key=0;for(int v=0;v<count;++v){int rank=0;for(int w=0;w<count;++w)rank+=numbers[w]<numbers[v];key=count*key+rank;}return key;
}
inline void solveSmall(std::vector<double> M,std::vector<double>& R,std::size_t n,std::size_t k){
    for(std::size_t c=0;c<n;++c){std::size_t p=c;for(std::size_t r=c+1;r<n;++r)if(std::abs(M[r*n+c])>std::abs(M[p*n+c]))p=r;
        if(!(std::abs(M[p*n+c])>0))throw std::runtime_error("Singular HEX charge monomial Gram");
        if(p!=c){for(std::size_t j=0;j<n;++j)std::swap(M[c*n+j],M[p*n+j]);for(std::size_t j=0;j<k;++j)std::swap(R[c*k+j],R[p*k+j]);}
        for(std::size_t r=0;r<n;++r){if(r==c)continue;const double f=M[r*n+c]/M[c*n+c];if(f==0.)continue;
            for(std::size_t j=c;j<n;++j)M[r*n+j]-=f*M[c*n+j];for(std::size_t j=0;j<k;++j)R[r*k+j]-=f*R[c*k+j];}}
    for(std::size_t c=0;c<n;++c)for(std::size_t j=0;j<k;++j)R[c*k+j]/=M[c*n+c];
}
// The element of NGSolve's periodic wrapper: the plain element with its vertex order.
inline ngfem::FiniteElement& hexElement(MaterialStagingContext& context,ngfem::ElementId element,ngstd::LocalHeap& heap){
    auto& fe=context.space->GetFE(element,heap);
    if(!context.vertexMap.empty()){int numbers[8];const auto vertices=context.mesh->GetElement(element).Vertices();
        for(int v=0;v<8;++v)numbers[v]=context.orientedVertex(vertices[v]);
        fe.SetVertexNumbers(ngcore::FlatArray<int>(8,numbers));}
    return fe;
}
inline void requireHexMesh(const ngcomp::MeshAccess& mesh,std::size_t ne){
    for(std::size_t e=0;e<ne;++e)if(mesh.GetElType(ngfem::ElementId(ngfem::VOL,int(e)))!=ngfem::ET_HEX)
        throw std::invalid_argument("Device HEX tables require a pure hexahedral mesh");
}
}
namespace radia::hdiv_cuda {
inline bool isHexMesh(const MaterialStagingContext& context){
    return context.ne&&context.mesh->GetElType(ngfem::ElementId(ngfem::VOL,0))==ngfem::ET_HEX;
}
// Rows: cell monomials xi^i eta^j zeta^k (k outer, then j, i inner, each
// exponent <= order), then exterior boundary quads in boundary-element order
// with u^i v^j (j outer); the CPU _mons_hex/_mons_quad order.
inline FlatChargeTables stageHexChargeTables(MaterialStagingContext& context){
    using namespace hex_staging;
    const auto& mesh=context.mesh;const auto& space=context.space;const int order=space->GetOrder();
    if(order<1||order>2||mesh->GetCurveOrder()>2)throw std::invalid_argument("Device HEX charge tables require BDM1/BDM2 HDiv on an affine, trilinear or Curve(2) hex mesh");
    requireHexMesh(*mesh,context.ne);
    FlatChargeTables out;out.dofs=context.nd;out.cells=context.ne;out.hex=true;
    for(int k=0;k<=order;++k)for(int j=0;j<=order;++j)for(int i=0;i<=order;++i)out.cellExponents.insert(out.cellExponents.end(),{i,j,k});
    for(int j=0;j<=order;++j)for(int i=0;i<=order;++i)out.faceExponents.insert(out.faceExponents.end(),{i,j});
    const std::size_t cm=out.cellExponents.size()/3,fm=out.faceExponents.size()/2;out.cellMonomials=cm;out.faceMonomials=fm;
    // Exact for the monomial Gram and the projected div/trace of the HDiv shapes.
    const ngfem::IntegrationRule cubeRule(ngfem::ET_HEX,2*order+4),squareRule(ngfem::ET_QUAD,2*order+4);
    auto monomials=[](const int* exponents,std::size_t count,int dimension,const double* x,double* m){
        for(std::size_t a=0;a<count;++a){double v=1;for(int c=0;c<dimension;++c)for(int r=0;r<exponents[dimension*a+c];++r)v*=x[c];m[a]=v;}};
    const ngfem::IntegrationPoint probes[2]={ngfem::IntegrationPoint(.137,.239,.361),ngfem::IntegrationPoint(.713,.211,.557)};
    ngstd::LocalHeap heap(1<<24,"radia_hdiv_hex_charge");
    std::vector<int> classKey;std::vector<double> classProbe;std::vector<int> cellVertexNumbers(8*context.ne);std::size_t local=0;
    for(std::size_t e=0;e<context.ne;++e){ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::VOL,int(e));
        const auto& fe=dynamic_cast<const ngfem::HDivFiniteElement<3>&>(hexElement(context,eid,heap));
        ngstd::Array<int> dofs;space->GetDofNrs(eid,dofs);
        if(!local)local=std::size_t(fe.GetNDof());
        if(std::size_t(fe.GetNDof())!=local||std::size_t(dofs.Size())!=local||local>128)
            throw std::invalid_argument("Device HEX charge tables require one local HDiv size <= 128");
        for(auto d:dofs)if(d<0||std::size_t(d)>=context.nd)throw std::invalid_argument("Constrained charge DOF unsupported");
        std::vector<int> sorted(dofs.begin(),dofs.end());std::sort(sorted.begin(),sorted.end());
        if(std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end())throw std::runtime_error("Duplicate charge local DOF unsupported");
        out.cellDofs.insert(out.cellDofs.end(),dofs.begin(),dofs.end());
        const auto vertexNumbers=mesh->GetElement(eid).Vertices();int numbers[8];
        for(int v=0;v<8;++v){numbers[v]=context.orientedVertex(vertexNumbers[v]);cellVertexNumbers[8*e+v]=vertexNumbers[v];}
        auto& trafo=mesh->GetTrafo(eid,heap);
        const auto at=out.cellVertices.size();latticeOf(trafo,3,out.cellVertices);const int orientation=checkHexMap(trafo,out.cellVertices.data()+at);
        ngbla::Vector<> div(local);ngbla::Matrix<> shape(local,3);std::vector<double> probe;
        for(const auto& ip:probes){fe.CalcDivShape(ip,div);fe.CalcShape(ip,shape);
            for(std::size_t i=0;i<local;++i){probe.push_back(div(i));for(int c=0;c<3;++c)probe.push_back(shape(i,c));}}
        const int key=2*rankKey(numbers,8)+(orientation<0);auto found=std::find(classKey.begin(),classKey.end(),key);const std::size_t cls=std::size_t(found-classKey.begin());
        if(found==classKey.end()){
            classKey.push_back(key);classProbe.insert(classProbe.end(),probe.begin(),probe.end());
            // T = Mref^{-1} R with R = -sgn(det J) int m(xi) div_ref phi(xi) dxi on the unit cube,
            // as int (-div u) m dx = sgn(det J) int (-div_ref phi) m dxi under the Piola map.
            std::vector<double> M(cm*cm,0.),R(cm*local,0.),m(cm);
            for(const auto& ip:cubeRule){const double xi[3]={ip(0),ip(1),ip(2)};monomials(out.cellExponents.data(),cm,3,xi,m.data());fe.CalcDivShape(ip,div);
                for(std::size_t a=0;a<cm;++a){for(std::size_t b=0;b<cm;++b)M[a*cm+b]+=ip.Weight()*m[a]*m[b];
                    for(std::size_t i=0;i<local;++i)R[a*local+i]-=ip.Weight()*m[a]*div(i);}}
            solveSmall(M,R,cm,local);
            for(double& v:R){v*=orientation;if(!std::isfinite(v))throw std::runtime_error("Nonfinite HEX charge cell table");}
            out.cellTables.insert(out.cellTables.end(),R.begin(),R.end());
        }else{const double* reference=classProbe.data()+cls*probe.size();
            for(std::size_t i=0;i<probe.size();++i)if(std::abs(probe[i]-reference[i])>1e-12*(1+std::abs(reference[i])))
                throw std::runtime_error("NGSolve HDiv HEX reference shapes are not determined by the element vertex order");}
        out.cellClass.push_back(int(cls));
    }
    out.cellWidth=local;out.cellClasses=classKey.size();
    // Exterior quads (sorted corner numbers owned by one cell, as _exterior_bnd_elements); two owners is an internal interface.
    std::map<std::array<int,4>,std::vector<std::pair<std::size_t,int>>> owners;
    for(std::size_t e=0;e<context.ne;++e)for(int f=0;f<6;++f){std::array<int,4> key;
        for(int v=0;v<4;++v)key[v]=cellVertexNumbers[8*e+hexFaces[f][v]];std::sort(key.begin(),key.end());owners[key].push_back({e,f});}
    std::vector<int> faceKey;std::vector<double> faceMoments;const std::size_t nse=std::size_t(mesh->GetNSE());
    for(std::size_t s=0;s<nse;++s){ngstd::HeapReset reset(heap);ngfem::ElementId sid(ngfem::BND,int(s));
        if(mesh->GetElType(sid)!=ngfem::ET_QUAD)throw std::invalid_argument("Device HEX charge tables require quadrilateral boundary elements");
        if(context.excludedBoundary(sid))continue;
        const auto vertexNumbers=mesh->GetElement(sid).Vertices();int numbers[4];std::array<int,4> key;
        for(int v=0;v<4;++v)key[v]=numbers[v]=vertexNumbers[v];std::sort(key.begin(),key.end());
        const auto owner=owners.find(key);const std::size_t count=owner==owners.end()?0:owner->second.size();
        if(count==2)continue;
        if(count!=1)throw std::runtime_error("HEX boundary element matches "+std::to_string(count)+" volume facets; the mesh is non-conforming or corrupt");
        const std::size_t cell=owner->second[0].first;ngfem::ElementId eid(ngfem::VOL,int(cell));
        const auto& fe=dynamic_cast<const ngfem::HDivFiniteElement<3>&>(hexElement(context,eid,heap));
        const double* P[4];
        for(int v=0;v<4;++v){int w=0;while(w<8&&cellVertexNumbers[8*cell+w]!=numbers[v])++w;
            if(w==8)throw std::runtime_error("HEX face vertex is not a vertex of its cell");P[v]=hexVertexXi[w];}
        // Boundary (u,v) -> cell xi is affine on a cube face: P2 = P1 + P3 - P0.
        for(int c=0;c<3;++c)if(P[2][c]!=P[1][c]+P[3][c]-P[0][c])throw std::runtime_error("HEX boundary vertex order is not a reference-square traversal of its cell face");
        const double e0[3]={P[1][0]-P[0][0],P[1][1]-P[0][1],P[1][2]-P[0][2]},e1[3]={P[3][0]-P[0][0],P[3][1]-P[0][1],P[3][2]-P[0][2]};
        const double nu[3]={e0[1]*e1[2]-e0[2]*e1[1],e0[2]*e1[0]-e0[0]*e1[2],e0[0]*e1[1]-e0[1]*e1[0]};
        ngstd::Array<int> faceDofs;space->GetDofNrs(sid,faceDofs);
        if(!out.faceWidth)out.faceWidth=std::size_t(faceDofs.Size());
        if(std::size_t(faceDofs.Size())!=out.faceWidth||!out.faceWidth||out.faceWidth>local)throw std::invalid_argument("Device HEX charge tables require one face DOF count");
        std::vector<std::size_t> column(out.faceWidth);
        for(std::size_t k=0;k<out.faceWidth;++k){std::size_t i=0;while(i<local&&out.cellDofs[cell*local+i]!=faceDofs[k])++i;
            if(i==local)throw std::runtime_error("HEX face DOF is not a DOF of its cell");column[k]=i;}
        // R = int m(u,v) phi_ref(xi(u,v)).nu du dv: the flux u.n ds in the reference measure.
        std::vector<double> M(fm*fm,0.),R(fm*local,0.),m(fm);ngbla::Matrix<> shape(local,3);
        for(const auto& ip:squareRule){const double g[2]={ip(0),ip(1)};monomials(out.faceExponents.data(),fm,2,g,m.data());
            ngfem::IntegrationPoint xi(P[0][0]+g[0]*e0[0]+g[1]*e1[0],P[0][1]+g[0]*e0[1]+g[1]*e1[1],P[0][2]+g[0]*e0[2]+g[1]*e1[2]);fe.CalcShape(xi,shape);
            for(std::size_t i=0;i<local;++i){const double flux=shape(i,0)*nu[0]+shape(i,1)*nu[1]+shape(i,2)*nu[2];
                for(std::size_t r=0;r<fm;++r)R[r*local+i]+=ip.Weight()*m[r]*flux;}
            for(std::size_t r=0;r<fm;++r)for(std::size_t t=0;t<fm;++t)M[r*fm+t]+=ip.Weight()*m[r]*m[t];}
        // Faces of one class may lie on different reference faces, and the monomial Gram amplifies their roundoff:
        // the class check compares these moments, whose equality fixes the solved table of every face in the class.
        std::vector<double> moments(fm*out.faceWidth);
        for(std::size_t r=0;r<fm;++r)for(std::size_t k=0;k<out.faceWidth;++k)moments[r*out.faceWidth+k]=R[r*local+column[k]];
        solveSmall(M,R,fm,local);
        double scale=0;for(double v:R){if(!std::isfinite(v))throw std::runtime_error("Nonfinite HEX charge face table");scale=std::max(scale,std::abs(v));}
        std::vector<unsigned char> onFace(local,0);for(auto i:column)onFace[i]=1;
        for(std::size_t r=0;r<fm;++r)for(std::size_t i=0;i<local;++i)if(!onFace[i]&&std::abs(R[r*local+i])>1e-12*scale)
            throw std::runtime_error("HDiv HEX normal trace is not carried by the boundary-element DOFs");
        std::vector<double> table(fm*out.faceWidth);
        for(std::size_t r=0;r<fm;++r)for(std::size_t k=0;k<out.faceWidth;++k)table[r*out.faceWidth+k]=R[r*local+column[k]];
        int faceRanks[4];for(int v=0;v<4;++v)faceRanks[v]=context.orientedVertex(numbers[v]);
        const int fkey=rankKey(faceRanks,4);auto found=std::find(faceKey.begin(),faceKey.end(),fkey);const std::size_t cls=std::size_t(found-faceKey.begin());
        if(found==faceKey.end()){faceKey.push_back(fkey);out.faceTables.insert(out.faceTables.end(),table.begin(),table.end());faceMoments.insert(faceMoments.end(),moments.begin(),moments.end());}
        else{const double* reference=faceMoments.data()+cls*moments.size();
            for(std::size_t i=0;i<moments.size();++i)if(std::abs(moments[i]-reference[i])>1e-12*(1+std::abs(reference[i])))
                throw std::runtime_error("NGSolve HDiv HEX face traces are not determined by the face vertex order");}
        out.faceClass.push_back(int(cls));out.faceDofs.insert(out.faceDofs.end(),faceDofs.begin(),faceDofs.end());
        latticeOf(mesh->GetTrafo(sid,heap),2,out.faceVertices);
        ++out.faces;
    }
    out.faceClasses=faceKey.size();
    constexpr std::size_t top=std::numeric_limits<int>::max();
    if(!out.faces||out.faces>top/fm||out.cells>(top-out.faces*fm)/cm)throw std::invalid_argument("Device HEX charge map outside the bounded scope");
    return out;
}
// HEX evaluation tables: NGSolve reference shapes per vertex-rank class on
// IntegrationRule(ET_HEX, 2 degree) and the 27-node Q2 lattice per element;
// the device evaluates J at every rule point from the lattice.
inline FlatEvaluationTables stageHexEvaluationTables(MaterialStagingContext& context){
    using namespace hex_staging;
    const auto& mesh=context.mesh;const auto& space=context.space;const auto& materialSpace=context.materialSpace;
    requireHexMesh(*mesh,context.ne);
    FlatEvaluationTables out;out.dofs=context.nd;out.totalPoints=context.totalPoints;out.elements=context.ne;out.integrationOrder=context.degree;out.nodes=27;
    const ngfem::IntegrationRule rule(ngfem::ET_HEX,2*context.degree);out.rulePoints=std::size_t(rule.Size());
    for(const auto& ip:rule){for(int d=0;d<3;++d)out.rule.push_back(ip(d));out.rule.push_back(ip.Weight());}
    const ngfem::IntegrationPoint probes[2]={ngfem::IntegrationPoint(.137,.239,.361),ngfem::IntegrationPoint(.713,.211,.557)};
    std::vector<int> classKey;std::vector<double> classProbe;std::vector<int> covered(context.totalPoints,0);
    ngstd::LocalHeap heap(1<<24,"radia_hdiv_hex_evaluation");
    for(std::size_t e=0;e<context.ne;++e){ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::VOL,int(e));
        const auto& fe=dynamic_cast<const ngfem::HDivFiniteElement<3>&>(hexElement(context,eid,heap));
        ngstd::Array<int> dofs;space->GetDofNrs(eid,dofs);
        if(!out.local)out.local=std::size_t(fe.GetNDof());
        if(std::size_t(fe.GetNDof())!=out.local||std::size_t(dofs.Size())!=out.local||out.local>128)
            throw std::invalid_argument("Device HEX evaluation tables require one local HDiv size <= 128");
        for(auto d:dofs)if(d<0||std::size_t(d)>=context.nd)throw std::invalid_argument("Constrained material DOF unsupported");
        std::vector<int> sorted(dofs.begin(),dofs.end());std::sort(sorted.begin(),sorted.end());
        if(std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end())throw std::runtime_error("Duplicate material local DOF unsupported");
        out.elementDofs.insert(out.elementDofs.end(),dofs.begin(),dofs.end());
        ngstd::Array<int> samples;materialSpace->GetDofNrs(eid,samples);
        if(std::size_t(samples.Size())!=out.rulePoints)throw std::runtime_error("Material rule point count mismatch");
        for(auto q:samples){if(q<0||std::size_t(q)>=context.totalPoints||covered[q]++)throw std::runtime_error("Material samples are not element-disjoint");
            out.elementPoints.push_back(q);}
        const auto vertexNumbers=mesh->GetElement(eid).Vertices();int numbers[8];for(int v=0;v<8;++v)numbers[v]=context.orientedVertex(vertexNumbers[v]);
        auto& trafo=mesh->GetTrafo(eid,heap);
        const auto at=out.vertices.size();latticeOf(trafo,3,out.vertices);checkHexMap(trafo,out.vertices.data()+at);
        ngbla::Matrix<> shape(out.local,3);std::vector<double> probe;
        for(const auto& ip:probes){fe.CalcShape(ip,shape);for(std::size_t i=0;i<out.local;++i)for(int c=0;c<3;++c)probe.push_back(shape(i,c));}
        const int key=rankKey(numbers,8);auto found=std::find(classKey.begin(),classKey.end(),key);const std::size_t cls=std::size_t(found-classKey.begin());
        if(found==classKey.end()){classKey.push_back(key);classProbe.insert(classProbe.end(),probe.begin(),probe.end());
            for(const auto& ip:rule){fe.CalcShape(ip,shape);for(std::size_t i=0;i<out.local;++i)for(int c=0;c<3;++c){
                if(!std::isfinite(shape(i,c)))throw std::runtime_error("Nonfinite reference HDiv shape");out.shapes.push_back(shape(i,c));}}}
        else{const double* reference=classProbe.data()+cls*probe.size();
            for(std::size_t i=0;i<probe.size();++i)if(std::abs(probe[i]-reference[i])>1e-12*(1+std::abs(reference[i])))
                throw std::runtime_error("NGSolve HDiv HEX reference shapes are not determined by the element vertex order");}
        out.elementClass.push_back(int(cls));
    }
    for(auto c:covered)if(c!=1)throw std::runtime_error("Material global sample coverage mismatch");
    return out;
}
}
