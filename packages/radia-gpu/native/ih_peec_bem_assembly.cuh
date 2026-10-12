// Bounded flat-P1 staging reuses the existing Galerkin singular integrals and
// finite-segment field kernels. No Python or CPU linear solve is involved.
#include "rad_bem_galerkin.h"
#include "rad_biot_savart_filaments.h"
mxArray* realMatrixOut(const std::vector<double>& rowMajor,size_t rows,size_t columns){
 auto* out=mxCreateDoubleMatrix(rows,columns,mxREAL);for(size_t i=0;i<rows;++i)for(size_t j=0;j<columns;++j)mxGetDoubles(out)[i+rows*j]=rowMajor[i*columns+j];return out;
}
std::vector<double> readRows(const mxArray* a,size_t columns,size_t limit){
 if(mxGetNumberOfDimensions(a)!=2||mxGetN(a)!=columns||mxGetM(a)<1||mxGetM(a)>limit)throw std::invalid_argument("Bounded flat geometry shape mismatch");
 auto columnMajor=realVector(a,mxGetM(a)*columns);std::vector<double> rows(columnMajor.size());for(size_t i=0;i<mxGetM(a);++i)for(size_t j=0;j<columns;++j)rows[i*columns+j]=columnMajor[i+mxGetM(a)*j];return rows;
}
mxArray* assemblePEECBEM(const mxArray* cfg){
 auto vertices=readRows(field(cfg,"vertices_m"),3,128),triangleInput=readRows(field(cfg,"triangles"),3,256),segments=readRows(field(cfg,"segments_m"),6,1024);
 const size_t nv=vertices.size()/3,nt=triangleInput.size()/3,ns=segments.size()/6;
 for(size_t i=0;i<ns;++i){double squared=0.;for(int d=0;d<3;++d){double delta=segments[i*6+d+3]-segments[i*6+d];squared+=delta*delta;}if(!std::isfinite(squared)||squared<=1e-30)throw std::invalid_argument("Degenerate filament segment");}
 auto assignment=realVector(field(cfg,"segment_coil"),ns);size_t coils=0;for(double x:assignment){if(x<1||x>16||x!=int(x))throw std::invalid_argument("Coil IDs must be 1..16");coils=std::max(coils,size_t(x));}
 std::vector<int64_t> triangles(triangleInput.size());std::map<std::pair<int64_t,int64_t>,std::pair<int,int>> edges;
 std::vector<double> nodes(nt*18),area(nt),normal(nt*3),gradient(nt*9),centers(nt*3),mass(nv*nv,0.),stiffness(nv*nv,0.),field(nt*3*nv,0.);
 auto cross=[](const double* a,const double* b,double* c){c[0]=a[1]*b[2]-a[2]*b[1];c[1]=a[2]*b[0]-a[0]*b[2];c[2]=a[0]*b[1]-a[1]*b[0];};
 for(size_t e=0;e<nt;++e){
  for(size_t k=0;k<3;++k){double index=triangleInput[e*3+k];if(index<1||index>nv||index!=int(index))throw std::invalid_argument("Triangle index outside vertices");triangles[e*3+k]=int64_t(index)-1;}
  auto* tri=&triangles[e*3];if(tri[0]==tri[1]||tri[1]==tri[2]||tri[2]==tri[0])throw std::invalid_argument("Repeated triangle vertex");
  const double *p[3]={&vertices[tri[0]*3],&vertices[tri[1]*3],&vertices[tri[2]*3]};double d1[3],d2[3],raw[3];for(int d=0;d<3;++d){d1[d]=p[1][d]-p[0][d];d2[d]=p[2][d]-p[0][d];}cross(d1,d2,raw);double twoArea=std::sqrt(raw[0]*raw[0]+raw[1]*raw[1]+raw[2]*raw[2]);if(!std::isfinite(twoArea)||twoArea<1e-30)throw std::invalid_argument("Degenerate P1 triangle");area[e]=twoArea*.5;
  for(int d=0;d<3;++d){normal[e*3+d]=raw[d]/twoArea;centers[e*3+d]=(p[0][d]+p[1][d]+p[2][d])/3.;for(int k=0;k<3;++k)nodes[e*18+k*3+d]=p[k][d];nodes[e*18+9+d]=(p[0][d]+p[1][d])*.5;nodes[e*18+12+d]=(p[1][d]+p[2][d])*.5;nodes[e*18+15+d]=(p[2][d]+p[0][d])*.5;}
  for(int k=0;k<3;++k){int a=(k+1)%3,b=(k+2)%3;double opposite[3],g[3];for(int d=0;d<3;++d)opposite[d]=p[b][d]-p[a][d];cross(&normal[e*3],opposite,g);for(int d=0;d<3;++d){gradient[e*9+k*3+d]=g[d]/twoArea;field[(e*3+d)*nv+tri[k]]=-g[d]/twoArea;}
   auto edge=std::minmax(tri[k],tri[(k+1)%3]);auto& count=edges[{edge.first,edge.second}];++count.first;count.second+=tri[k]<tri[(k+1)%3]?1:-1;}
  for(int i=0;i<3;++i)for(int j=0;j<3;++j){mass[tri[i]*nv+tri[j]]+=area[e]*(i==j?2.:1.)/12.;double dot=0.;for(int d=0;d<3;++d)dot+=gradient[e*9+i*3+d]*gradient[e*9+j*3+d];stiffness[tri[i]*nv+tri[j]]+=area[e]*dot;}
 }
 for(const auto& edge:edges)if(edge.second.first!=2||edge.second.second!=0)throw std::invalid_argument("Closed consistently oriented genus-zero surface required");
 std::vector<bool> connected(nv,false);connected[0]=true;for(size_t pass=0;pass<nv;++pass)for(const auto& edge:edges){auto a=edge.first.first,b=edge.first.second;if(connected[a]||connected[b])connected[a]=connected[b]=true;}for(bool present:connected)if(!present)throw std::invalid_argument("Body surface must be connected and use every vertex");
 if(int64_t(nv)-int64_t(edges.size())+int64_t(nt)!=2)throw std::invalid_argument("Loop/cohomology body requires a separate unsupported frontend");
 double volume=0.;for(size_t e=0;e<nt;++e)for(int d=0;d<3;++d)volume+=area[e]*centers[e*3+d]*normal[e*3+d]/3.;if(!std::isfinite(volume)||volume<=0.)throw std::invalid_argument("Surface must have positive outward oriented volume");
 auto integerOption=[&](const char* name,int low,int high){double x=number(field(cfg,name));if(x<low||x>high||x!=int(x))throw std::invalid_argument("Invalid declared quadrature order");return int(x);};
 std::vector<double> sl(nv*nv),dl(nv*nv);radia::bem::AssembleSLDL(vertices.data(),int(nv),triangles.data(),int(nt),nodes.data(),integerOption("regular_quad_degree",5,31),integerOption("singular_order",2,12),1,sl.data(),dl.data());
 std::vector<double> h(nv*3*coils),a(nt*3*coils),weak(nv*coils,0.),currents(ns),zero(ns,0.),hr(nv*3),hi(nv*3),ar(nt*3),ai(nt*3);
 for(size_t coil=0;coil<coils;++coil){bool present=false;for(size_t s=0;s<ns;++s){currents[s]=assignment[s]==coil+1?1.:0.;present|=currents[s]!=0.;}if(!present)throw std::invalid_argument("Missing coil ID");
  radia::bs::HFromSegmentsComplex(segments.data(),int(ns),vertices.data(),int(nv),currents.data(),zero.data(),hr.data(),hi.data(),1);
  radia::bs::AFromSegmentsComplex(segments.data(),int(ns),centers.data(),int(nt),currents.data(),zero.data(),ar.data(),ai.data(),1);
  for(size_t i=0;i<nv*3;++i)h[i*coils+coil]=hr[i];for(size_t i=0;i<nt*3;++i)a[i*coils+coil]=ar[i];
  for(size_t e=0;e<nt;++e)for(int k=0;k<3;++k){double dot=0.;for(int d=0;d<3;++d){double mean=0.;for(int c=0;c<3;++c)mean+=hr[triangles[e*3+c]*3+d]/3.;dot+=gradient[e*9+k*3+d]*mean;}weak[triangles[e*3+k]*coils+coil]-=area[e]*dot;}
 }
 const char* names[]={"mass","surface_stiffness","single_layer","double_layer","field_gradient","face_area_m2","face_normals","unit_h","unit_a","incident_weak_map","body_volume_m3"};auto* out=mxCreateStructMatrix(1,1,11,names);
 auto set=[&](int index,const std::vector<double>& values,size_t rows,size_t columns){mxSetField(out,0,names[index],realMatrixOut(values,rows,columns));};
 set(0,mass,nv,nv);set(1,stiffness,nv,nv);set(2,sl,nv,nv);set(3,dl,nv,nv);set(4,field,nt*3,nv);set(5,area,nt,1);set(6,normal,nt,3);set(7,h,nv*3,coils);set(8,a,nt*3,coils);set(9,weak,nv,coils);mxSetField(out,0,names[10],mxCreateDoubleScalar(volume));return out;
}
