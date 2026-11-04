// NGSolve-owned reference tables of straight tet, pyramid and hex meshes for device
// construction of the charge map B, the evaluation map E and the charge Gram.
// As on the pure hex route, charges are reference monomials in the reference measure,
// so B rows depend on the geometry only through sgn det J: per element type, vertex-rank
// class and orientation, cell rows are sgn det J times the reference projection of
// -div_ref phi and face rows the projection of phi_ref . nu over the boundary element's
// own coordinates. Every host is also stored as the Q2 lattice of its trilinear map
// from the unit cube (hdiv_typed_host_math.inc). No B, E or Gram value is produced here.
#pragma once
#include <algorithm>
#include <array>
#include <map>
#include <vector>
namespace radia::hdiv_cuda::mixed_staging {
enum HostType {hexHost=0,quadHost=1,tetHost=2,triHost=3,pyramidHost=4};
struct CellShape {ngfem::ELEMENT_TYPE type;int host,vertices,nodes;double xi[8][3];int faces,faceSize[6],face[6][4];};
// NGSolve reference vertices and faces.
inline const CellShape& cellShape(ngfem::ELEMENT_TYPE type){
    static const CellShape tet{ngfem::ET_TET,tetHost,4,4,{{1,0,0},{0,1,0},{0,0,1},{0,0,0}},4,{3,3,3,3},{{0,1,2},{0,1,3},{0,2,3},{1,2,3}}};
    static const CellShape pyramid{ngfem::ET_PYRAMID,pyramidHost,5,5,{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1}},5,{4,3,3,3,3},{{0,1,2,3},{0,1,4},{1,2,4},{2,3,4},{3,0,4}}};
    static const CellShape hex{ngfem::ET_HEX,hexHost,8,27,{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}},6,{4,4,4,4,4,4},
        {{0,1,2,3},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}}};
    if(type==ngfem::ET_TET)return tet;
    if(type==ngfem::ET_PYRAMID)return pyramid;
    if(type==ngfem::ET_HEX)return hex;
    throw std::invalid_argument("Mixed device tables support tetrahedra, pyramids and hexahedra only");
}
inline std::vector<std::array<int,3>> cellMonomials(int host,int order){
    std::vector<std::array<int,3>> out;
    if(host==hexHost){for(int k=0;k<=order;++k)for(int j=0;j<=order;++j)for(int i=0;i<=order;++i)out.push_back({i,j,k});}
    else if(host==tetHost){for(int total=0;total<order;++total)for(int k=0;k<=total;++k)for(int j=0;j<=total-k;++j)out.push_back({total-k-j,j,k});}
    else{out.push_back({0,0,0});if(order==2)for(const auto& e:{std::array<int,3>{1,0,0},std::array<int,3>{0,1,0},std::array<int,3>{0,0,1},std::array<int,3>{1,1,-1}})out.push_back(e);}
    return out;
}
inline std::vector<std::array<int,3>> faceMonomials(int host,int order){
    std::vector<std::array<int,3>> out;
    if(host==quadHost){for(int j=0;j<=order;++j)for(int i=0;i<=order;++i)out.push_back({i,j,0});}
    else{for(int total=0;total<=order;++total)for(int j=0;j<=total;++j)out.push_back({total-j,j,0});}
    return out;
}
inline double monomial(const std::array<int,3>& e,const double* x){
    double v=1;for(int c=0;c<2;++c)for(int r=0;r<e[c];++r)v*=x[c];
    if(e[2]<0)return v/(1-x[2]);
    for(int r=0;r<e[2];++r)v*=x[2];
    return v;
}
inline void vertexPoints(const ngcomp::MeshAccess& mesh,ngfem::ElementId eid,int count,std::vector<double>& out,int* numbers){
    const auto vertices=mesh.GetElement(eid).Vertices();
    if(int(vertices.Size())!=count)throw std::runtime_error("Mixed element vertex count disagrees with its type");
    for(int v=0;v<count;++v){numbers[v]=vertices[v];const auto point=mesh.GetPoint<3>(vertices[v]);for(int d=0;d<3;++d)out.push_back(point(d));}
}
// X and dX/dxi of the straight map at the reference point xi, from the element vertices.
inline double straightMap(int host,const double* V,const double* xi,double X[3],double J[3][3]){
    if(host==tetHost){for(int k=0;k<3;++k){X[k]=xi[0]*V[3*0+k]+xi[1]*V[3*1+k]+xi[2]*V[3*2+k]+(1-xi[0]-xi[1]-xi[2])*V[3*3+k];
            for(int c=0;c<3;++c)J[k][c]=V[3*c+k]-V[9+k];}}
    else{const double top=1-xi[2],u=xi[0]/top,v=xi[1]/top;
        for(int k=0;k<3;++k){const double e1=V[3+k]-V[k],e3=V[9+k]-V[k],b=V[6+k]-V[3+k]-V[9+k]+V[k],ez=V[12+k]-V[k];
            X[k]=V[k]+xi[0]*e1+xi[1]*e3+xi[2]*ez+xi[0]*v*b;J[k][0]=e1+v*b;J[k][1]=e3+u*b;J[k][2]=ez+u*v*b;}}
    return J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
}
// Netgen 6.2.2608 maps a pyramid point with z (1 - 1e-10), so the check on a pyramid allows that constant.
inline int checkStraightMap(ngfem::ElementTransformation& trafo,int host,const double* V){
    double size=0;for(int i=1;i<(host==tetHost?4:5);++i)for(int d=0;d<3;++d)size=std::max(size,std::abs(V[3*i+d]-V[d]));
    const double tolerance=(host==pyramidHost?1e-10:1e-12)*std::max(1.,size);int sign=0;
    for(const auto& ip:{ngfem::IntegrationPoint(.137,.239,.361),ngfem::IntegrationPoint(.211,.113,.557),ngfem::IntegrationPoint(.319,.173,.123)}){
        ngfem::MappedIntegrationPoint<3,3> mip(ip,trafo);const double xi[3]={ip(0),ip(1),ip(2)};double X[3],J[3][3];const double det=straightMap(host,V,xi,X,J);
        for(int c=0;c<3;++c)if(std::abs(X[c]-mip.GetPoint()(c))>tolerance)throw std::runtime_error("NGSolve tet or pyramid map differs from the straight device map");
        const int s=det>0?1:det<0?-1:0;
        if(!s||(sign&&s!=sign))throw std::runtime_error("Mixed device route requires a Jacobian of one nonzero sign");
        sign=s;
    }
    return sign;
}
// Q2 lattice of the trilinear map from the unit cube; lattices of tets, pyramids and triangles collapse corners.
inline void cubeLattice(int host,const double* V,std::vector<double>& out){
    if(host==triHost){
        for(int iv=0;iv<3;++iv)for(int iu=0;iu<3;++iu){const double u=.5*iu,v=.5*iv,x=u*(1-v);
            for(int k=0;k<3;++k)out.push_back(x*V[k]+v*V[3+k]+(1-x-v)*V[6+k]);}
        return;
    }
    for(int iz=0;iz<3;++iz)for(int iy=0;iy<3;++iy)for(int ix=0;ix<3;++ix){const double a=.5*ix,b=.5*iy,c=.5*iz;
        for(int k=0;k<3;++k){
            if(host==tetHost){const double x=a*(1-b)*(1-c),y=b*(1-c);out.push_back(x*V[k]+y*V[3+k]+c*V[6+k]+(1-x-y-c)*V[9+k]);}
            else out.push_back((1-c)*((1-a)*(1-b)*V[k]+a*(1-b)*V[3+k]+a*b*V[6+k]+(1-a)*b*V[9+k])+c*V[12+k]);}}
}
inline void uniqueDofs(const ngstd::Array<int>& dofs,std::size_t nd){
    for(auto d:dofs)if(d<0||std::size_t(d)>=nd)throw std::invalid_argument("Constrained mixed DOF unsupported");
    std::vector<int> sorted(dofs.begin(),dofs.end());std::sort(sorted.begin(),sorted.end());
    if(std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end())throw std::runtime_error("Duplicate mixed local DOF unsupported");
}
}
namespace radia::hdiv_cuda {
inline bool isMixedMesh(const MaterialStagingContext& context){
    bool tet=false,pyramid=false,hex=false;
    for(std::size_t e=0;e<context.ne;++e){const auto type=context.mesh->GetElType(ngfem::ElementId(ngfem::VOL,int(e)));
        tet=tet||type==ngfem::ET_TET;pyramid=pyramid||type==ngfem::ET_PYRAMID;hex=hex||type==ngfem::ET_HEX;}
    return pyramid||(tet&&hex);
}
inline FlatChargeTables stageMixedChargeTables(MaterialStagingContext& context){
    using namespace mixed_staging;
    const auto& mesh=context.mesh;const auto& space=context.space;const int order=space->GetOrder();
    if(order<1||order>2||mesh->GetCurveOrder()>1)throw std::invalid_argument("Mixed device charge tables require BDM1/BDM2 HDiv on a straight mesh");
    FlatChargeTables out;out.dofs=context.nd;out.cells=context.ne;out.mixed=true;
    ngstd::LocalHeap heap(1<<24,"radia_hdiv_mixed_charge");
    std::vector<const CellShape*> shapes(context.ne);
    struct ClassEntry {int key;std::size_t size;std::vector<double> probe;};std::vector<ClassEntry> classes;
    auto addCharges=[&](const std::vector<std::array<int,3>>& mons,int host,int kind){
        for(const auto& e:mons){out.chargeHost.push_back(host);out.chargeKind.push_back(kind);out.chargeExponents.insert(out.chargeExponents.end(),e.begin(),e.end());}};
    auto addClass=[&](int key,const std::vector<double>& table,std::size_t monomials,std::size_t width,const std::vector<double>& probe,const char* what){
        for(std::size_t c=0;c<classes.size();++c)if(classes[c].key==key){
            for(std::size_t i=0;i<probe.size();++i)if(std::abs(probe[i]-classes[c].probe[i])>1e-12*(1+std::abs(classes[c].probe[i])))
                throw std::runtime_error(std::string("NGSolve HDiv ")+what+" are not determined by the element vertex order");
            return int(c);}
        out.classInfo.insert(out.classInfo.end(),{int(out.classTables.size()),int(monomials),int(width)});
        out.classTables.insert(out.classTables.end(),table.begin(),table.end());classes.push_back({key,table.size(),probe});
        return int(classes.size()-1);};
    const ngfem::IntegrationPoint probes[2]={ngfem::IntegrationPoint(.137,.239,.361),ngfem::IntegrationPoint(.211,.113,.557)};
    for(std::size_t e=0;e<context.ne;++e){ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::VOL,int(e));
        const auto type=mesh->GetElType(eid);const CellShape& shape=cellShape(type);shapes[e]=&shape;
        const auto& fe=dynamic_cast<const ngfem::HDivFiniteElement<3>&>(space->GetFE(eid,heap));
        ngstd::Array<int> dofs;space->GetDofNrs(eid,dofs);const std::size_t local=std::size_t(fe.GetNDof());
        if(std::size_t(dofs.Size())!=local||local>128)throw std::invalid_argument("Mixed device charge tables require local HDiv size <= 128");
        uniqueDofs(dofs,context.nd);out.hostDofs.insert(out.hostDofs.end(),dofs.begin(),dofs.end());
        auto& trafo=mesh->GetTrafo(eid,heap);int numbers[8];int orientation;
        const auto vertexNumbers=mesh->GetElement(eid).Vertices();for(int v=0;v<shape.vertices;++v)numbers[v]=vertexNumbers[v];
        if(shape.host==hexHost){const auto at=out.cellVertices.size();hex_staging::latticeOf(trafo,3,out.cellVertices);orientation=hex_staging::checkHexMap(trafo,out.cellVertices.data()+at);}
        else{std::vector<double> V;int check[8];vertexPoints(*mesh,eid,shape.vertices,V,check);orientation=checkStraightMap(trafo,shape.host,V.data());cubeLattice(shape.host,V.data(),out.cellVertices);}
        out.hostTypes.push_back(shape.host);
        const auto mons=cellMonomials(shape.host,order);const std::size_t cm=mons.size();addCharges(mons,int(e),0);
        ngbla::Vector<> div(local);ngbla::Matrix<> values(local,3);std::vector<double> probe;
        for(const auto& ip:probes){fe.CalcDivShape(ip,div);fe.CalcShape(ip,values);
            for(std::size_t i=0;i<local;++i){probe.push_back(div(i));for(int c=0;c<3;++c)probe.push_back(values(i,c));}}
        const int key=2*((shape.host<<25)+hex_staging::rankKey(numbers,shape.vertices))+(orientation<0);
        const bool known=std::any_of(classes.begin(),classes.end(),[&](const ClassEntry& c){return c.key==key;});
        std::vector<double> table;
        if(!known){
            // T = Mref^{-1} R with R = -sgn(det J) int m div_ref phi on the reference element.
            const ngfem::IntegrationRule rule(type,2*order+6);std::vector<double> M(cm*cm,0.),R(cm*local,0.),m(cm);
            for(const auto& ip:rule){const double xi[3]={ip(0),ip(1),ip(2)};for(std::size_t a=0;a<cm;++a)m[a]=monomial(mons[a],xi);
                fe.CalcDivShape(ip,div);
                for(std::size_t a=0;a<cm;++a){for(std::size_t b=0;b<cm;++b)M[a*cm+b]+=ip.Weight()*m[a]*m[b];
                    for(std::size_t i=0;i<local;++i)R[a*local+i]-=ip.Weight()*m[a]*div(i);}}
            hex_staging::solveSmall(M,R,cm,local);
            for(double& v:R){v*=orientation;if(!std::isfinite(v))throw std::runtime_error("Nonfinite mixed charge cell table");}
            table=std::move(R);
        }
        out.hostClass.push_back(addClass(key,table,cm,local,probe,"reference shapes"));
    }
    // Exterior faces: sorted corner numbers owned by one cell, as _exterior_bnd_elements.
    std::map<std::vector<int>,std::vector<std::size_t>> owners;
    for(std::size_t e=0;e<context.ne;++e){const CellShape& shape=*shapes[e];ngfem::ElementId eid(ngfem::VOL,int(e));const auto numbers=mesh->GetElement(eid).Vertices();
        for(int f=0;f<shape.faces;++f){std::vector<int> key;for(int v=0;v<shape.faceSize[f];++v)key.push_back(numbers[shape.face[f][v]]);std::sort(key.begin(),key.end());owners[key].push_back(e);}}
    const std::size_t nse=std::size_t(mesh->GetNSE());
    for(std::size_t s=0;s<nse;++s){ngstd::HeapReset reset(heap);ngfem::ElementId sid(ngfem::BND,int(s));
        const auto type=mesh->GetElType(sid);
        if(type!=ngfem::ET_QUAD&&type!=ngfem::ET_TRIG)throw std::invalid_argument("Mixed device charge tables require triangular or quadrilateral boundary elements");
        const int host=type==ngfem::ET_QUAD?quadHost:triHost,corners=host==quadHost?4:3;
        const auto vertexNumbers=mesh->GetElement(sid).Vertices();int numbers[4];std::vector<int> key;
        for(int v=0;v<corners;++v){key.push_back(numbers[v]=vertexNumbers[v]);}std::sort(key.begin(),key.end());
        const auto owner=owners.find(key);const std::size_t count=owner==owners.end()?0:owner->second.size();
        if(count==2)continue;
        if(count!=1)throw std::runtime_error("Mixed boundary element matches "+std::to_string(count)+" volume facets; the mesh is non-conforming or corrupt");
        const std::size_t cell=owner->second[0];ngfem::ElementId eid(ngfem::VOL,int(cell));const CellShape& shape=*shapes[cell];
        const auto& fe=dynamic_cast<const ngfem::HDivFiniteElement<3>&>(space->GetFE(eid,heap));const std::size_t local=std::size_t(fe.GetNDof());
        const auto cellNumbers=mesh->GetElement(eid).Vertices();const double* P[4];
        for(int v=0;v<corners;++v){int w=0;while(w<shape.vertices&&cellNumbers[w]!=numbers[v])++w;
            if(w==shape.vertices)throw std::runtime_error("Mixed face vertex is not a vertex of its cell");P[v]=shape.xi[w];}
        // Boundary coordinates to cell xi are affine; the triangle has vertices (1,0),(0,1),(0,0), the quad (0,0),(1,0),(1,1),(0,1).
        double origin[3],e0[3],e1[3];
        for(int c=0;c<3;++c){if(host==quadHost){if(P[2][c]!=P[1][c]+P[3][c]-P[0][c])throw std::runtime_error("Mixed boundary vertex order is not a reference-square traversal of its cell face");
                origin[c]=P[0][c];e0[c]=P[1][c]-P[0][c];e1[c]=P[3][c]-P[0][c];}
            else{origin[c]=P[2][c];e0[c]=P[0][c]-P[2][c];e1[c]=P[1][c]-P[2][c];}}
        const double nu[3]={e0[1]*e1[2]-e0[2]*e1[1],e0[2]*e1[0]-e0[0]*e1[2],e0[0]*e1[1]-e0[1]*e1[0]};
        ngstd::Array<int> faceDofs,localDofs;space->GetDofNrs(sid,faceDofs);space->GetDofNrs(eid,localDofs);const std::size_t width=std::size_t(faceDofs.Size());
        if(!width||width>local)throw std::invalid_argument("Mixed device charge tables require boundary DOFs of the owning cell");
        std::vector<std::size_t> column(width);
        for(std::size_t k=0;k<width;++k){std::size_t i=0;while(i<local&&localDofs[i]!=faceDofs[k])++i;
            if(i==local)throw std::runtime_error("Mixed face DOF is not a DOF of its cell");column[k]=i;}
        const auto mons=faceMonomials(host,order);const std::size_t fm=mons.size();
        const ngfem::IntegrationRule rule(type,2*order+6);std::vector<double> M(fm*fm,0.),R(fm*local,0.),m(fm);ngbla::Matrix<> values(local,3);
        for(const auto& ip:rule){const double g[3]={ip(0),ip(1),0.};for(std::size_t a=0;a<fm;++a)m[a]=monomial(mons[a],g);
            ngfem::IntegrationPoint xi(origin[0]+g[0]*e0[0]+g[1]*e1[0],origin[1]+g[0]*e0[1]+g[1]*e1[1],origin[2]+g[0]*e0[2]+g[1]*e1[2]);fe.CalcShape(xi,values);
            for(std::size_t i=0;i<local;++i){const double flux=values(i,0)*nu[0]+values(i,1)*nu[1]+values(i,2)*nu[2];
                for(std::size_t r=0;r<fm;++r)R[r*local+i]+=ip.Weight()*m[r]*flux;}
            for(std::size_t r=0;r<fm;++r)for(std::size_t t=0;t<fm;++t)M[r*fm+t]+=ip.Weight()*m[r]*m[t];}
        std::vector<double> moments(fm*width);
        for(std::size_t r=0;r<fm;++r)for(std::size_t k=0;k<width;++k)moments[r*width+k]=R[r*local+column[k]];
        hex_staging::solveSmall(M,R,fm,local);
        double scale=0;for(double v:R){if(!std::isfinite(v))throw std::runtime_error("Nonfinite mixed charge face table");scale=std::max(scale,std::abs(v));}
        std::vector<unsigned char> onFace(local,0);for(auto i:column)onFace[i]=1;
        for(std::size_t r=0;r<fm;++r)for(std::size_t i=0;i<local;++i)if(!onFace[i]&&std::abs(R[r*local+i])>1e-10*scale)
            throw std::runtime_error("HDiv mixed normal trace is not carried by the boundary-element DOFs");
        std::vector<double> table(fm*width);
        for(std::size_t r=0;r<fm;++r)for(std::size_t k=0;k<width;++k)table[r*width+k]=R[r*local+column[k]];
        out.hostClass.push_back(addClass(-1-((host<<25)+hex_staging::rankKey(numbers,corners)),table,fm,width,moments,"face traces"));
        out.hostDofs.insert(out.hostDofs.end(),faceDofs.begin(),faceDofs.end());
        out.hostTypes.push_back(host);addCharges(mons,int(out.faces),1);
        if(host==quadHost)hex_staging::latticeOf(mesh->GetTrafo(sid,heap),2,out.faceVertices);
        else{std::vector<double> V;int check[8];vertexPoints(*mesh,sid,3,V,check);cubeLattice(triHost,V.data(),out.faceVertices);}
        ++out.faces;
    }
    if(!out.faces)throw std::invalid_argument("Mixed device charge tables require exterior boundary elements");
    return out;
}
// One evaluation table set per element type present, with global DOF and point numbers; the device
// evaluates J from the lattice (hex) or the vertices (tet, pyramid) at every rule point.
inline std::vector<FlatEvaluationTables> stageMixedEvaluationTables(MaterialStagingContext& context){
    using namespace mixed_staging;
    const auto& mesh=context.mesh;const auto& space=context.space;const auto& materialSpace=context.materialSpace;
    if(mesh->GetCurveOrder()>1)throw std::invalid_argument("Mixed device evaluation tables require a straight mesh");
    std::vector<FlatEvaluationTables> groups;std::vector<int> covered(context.totalPoints,0);
    ngstd::LocalHeap heap(1<<24,"radia_hdiv_mixed_evaluation");
    const ngfem::IntegrationPoint probes[2]={ngfem::IntegrationPoint(.137,.239,.361),ngfem::IntegrationPoint(.211,.113,.557)};
    for(const auto type:{ngfem::ET_TET,ngfem::ET_PYRAMID,ngfem::ET_HEX}){const CellShape& shape=cellShape(type);
        FlatEvaluationTables out;out.dofs=context.nd;out.totalPoints=context.totalPoints;out.integrationOrder=context.degree;out.nodes=std::size_t(shape.nodes);
        const ngfem::IntegrationRule rule(type,2*context.degree);out.rulePoints=std::size_t(rule.Size());
        for(const auto& ip:rule){for(int d=0;d<3;++d)out.rule.push_back(ip(d));out.rule.push_back(ip.Weight());}
        std::vector<int> classKey;std::vector<double> classProbe;
        for(std::size_t e=0;e<context.ne;++e){ngstd::HeapReset reset(heap);ngfem::ElementId eid(ngfem::VOL,int(e));
            if(mesh->GetElType(eid)!=type)continue;
            const auto& fe=dynamic_cast<const ngfem::HDivFiniteElement<3>&>(space->GetFE(eid,heap));
            ngstd::Array<int> dofs;space->GetDofNrs(eid,dofs);
            if(!out.local)out.local=std::size_t(fe.GetNDof());
            if(std::size_t(fe.GetNDof())!=out.local||std::size_t(dofs.Size())!=out.local||out.local>128)
                throw std::invalid_argument("Mixed device evaluation tables require one local HDiv size <= 128 per element type");
            uniqueDofs(dofs,context.nd);out.elementDofs.insert(out.elementDofs.end(),dofs.begin(),dofs.end());
            ngstd::Array<int> samples;materialSpace->GetDofNrs(eid,samples);
            if(std::size_t(samples.Size())!=out.rulePoints)throw std::runtime_error("Material rule point count mismatch");
            for(auto q:samples){if(q<0||std::size_t(q)>=context.totalPoints||covered[q]++)throw std::runtime_error("Material samples are not element-disjoint");
                out.elementPoints.push_back(q);}
            auto& trafo=mesh->GetTrafo(eid,heap);int numbers[8];
            if(shape.host==hexHost){const auto vertexNumbers=mesh->GetElement(eid).Vertices();for(int v=0;v<8;++v)numbers[v]=vertexNumbers[v];
                const auto at=out.vertices.size();hex_staging::latticeOf(trafo,3,out.vertices);hex_staging::checkHexMap(trafo,out.vertices.data()+at);}
            else{const auto at=out.vertices.size();vertexPoints(*mesh,eid,shape.vertices,out.vertices,numbers);checkStraightMap(trafo,shape.host,out.vertices.data()+at);}
            ngbla::Matrix<> values(out.local,3);std::vector<double> probe;
            for(const auto& ip:probes){fe.CalcShape(ip,values);for(std::size_t i=0;i<out.local;++i)for(int c=0;c<3;++c)probe.push_back(values(i,c));}
            const int key=hex_staging::rankKey(numbers,shape.vertices);auto found=std::find(classKey.begin(),classKey.end(),key);const std::size_t cls=std::size_t(found-classKey.begin());
            if(found==classKey.end()){classKey.push_back(key);classProbe.insert(classProbe.end(),probe.begin(),probe.end());
                for(const auto& ip:rule){fe.CalcShape(ip,values);for(std::size_t i=0;i<out.local;++i)for(int c=0;c<3;++c){
                    if(!std::isfinite(values(i,c)))throw std::runtime_error("Nonfinite reference HDiv shape");out.shapes.push_back(values(i,c));}}}
            else{const double* reference=classProbe.data()+cls*probe.size();
                for(std::size_t i=0;i<probe.size();++i)if(std::abs(probe[i]-reference[i])>1e-12*(1+std::abs(reference[i])))
                    throw std::runtime_error("NGSolve HDiv reference shapes are not determined by the element vertex order");}
            out.elementClass.push_back(int(cls));++out.elements;
        }
        if(out.elements)groups.push_back(std::move(out));
    }
    for(auto c:covered)if(c!=1)throw std::runtime_error("Material global sample coverage mismatch");
    return groups;
}
}
