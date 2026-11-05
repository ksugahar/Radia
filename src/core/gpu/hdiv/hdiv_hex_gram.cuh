// HEX charge Gram generated on the device from Q2 host lattices (27 nodes
// per cell, 9 per boundary quad, node ix + 3 iy + 9 iz), reference-measure
// charge monomials and 1D Gauss rules. The directed host-pair dispatch
// mirrors the CPU QuadBlockHex without images: BDM1 touching conforming pairs
// take the pair-domain Duffy rule, non-touching near-band pairs the plain pair
// product, separated pairs the outer tensor far product; BDM2 affine near
// pairs take the analytic sub-simplex source potential on the outer tensor
// rule, other near pairs the composite radial rule. Anisotropic or
// non-conforming BDM1 near pairs take the CPU graded near tensor family:
// the smootherstep-graded outer tensor rule with the analytic inner on affine
// sources and the cone/fan inner on distorted ones. Every host pair of the pair table (all pairs for dense, the
// near candidates for ACA) is computed once per direction and stored as
// 0.5 (G_AB + G_BA^T), the CPU GetHexSymBlock; beyond the table ACA evaluates the
// far product once per host pair in a dense leaf and per entry on low-rank rows and columns. Only geometry and charges cross H2D,
// the 1D rules as kernel arguments; no Gram value leaves the device.
// A mixed tet, pyramid and hex mesh adds the host types and takes one family for
// every pair, as the CPU QuadBlockTyped: the graded near rule on self, touching and
// near pairs, the outer tensor rule beyond them, densities in the unit-cube measure.
// A pair of two different hosts integrates one direction, with the cheaper host as
// source, and the block is not symmetrized.
// Mirror and cyclic images (pure hex only) follow the CPU QuadBlockHex image dispatch: a directed
// image block moves its target points by T^-1, each image keeps its own pair table section and
// symmetrized blocks, and an entry folds the direct term with sign * 0.5 (G_T(a,b) + G_T(b,a)).
#pragma once
#include <chrono>
#include <climits>
#include <cmath>
#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_scan.cuh>
#include "hmatrix_resident.cuh"
#define RADIA_TYPED_KERNEL __device__
#define RADIA_TYPED_FIELD rad_hdiv_cuda
#include "hdiv_typed_host_math.inc"
namespace radia_hmatrix_cuda {
constexpr int hexRuleCap=16,hexThreads=128;
// Scratch of a graded block: source charges x 84 polynomial coefficients. Only a mixed mesh grades at order 2.
__device__ __host__ inline int hexGradedCoeff(int order){return (order==2?27:8)*84;}
constexpr double hexInvFourPi=0.07957747154594766788;
// CPU HEX_EXPONENTIAL_RULE_MAX_ANISOTROPY_DEFAULT, HEX_AFFINE_EXACT_NEAR_FACTOR and HexDistortedFarFactor().
constexpr double hexMaxAnisotropy=12.0,hexAffineExactNear=1.0,hexDistortedFar=1.0;
// CPU HOST_CONE_SINH_MIN_WIDTH and HOST_CONE_SINH_MAX_WIDTH.
constexpr double hexSinhMinWidth=1e-3,hexSinhMaxWidth=1.0;
struct HexRule {double x[hexRuleCap],w[hexRuleCap];int n;};
// One image T = M R(angle): mirror flips, cos/sin of the forward and of the inverse rotation, the
// folded sign and a nonzero flag when the image rotates. All float64: it crosses H2D as one ledger row.
struct HexImage {double flip[3],c,s,ci,si,sign,rotated;};
static_assert(sizeof(HexImage)==72,"hdiv_gpu._IMAGE_ROW_BYTES ledger contract");
constexpr int hexMaxImages=63;  // CPU GetHexBlock refuses a larger image index.
constexpr int hexFarLevels=3;
enum HexClass:int {hexFarClass=0,hexDuffyClass=1,hexProductClass=2,hexAffineClass=3,hexCompositeClass=4,hexGradedClass=5,hexSubPairClass=6};
// Global host g: cells 0..cells-1, then boundary quads. Charges of one host
// are contiguous (first[g], count[g]); exponents are reference xi powers.
struct HexDevice {
    const double* cellNodes;const double* faceNodes;int cells,hosts;
    const int* first;const int* count;const int* chargeHost;const int* exponents;
    const double* center;const double* size;const unsigned char* affine;
    // BDM2 affine hosts: physical polynomial coefficients of every reference
    // monomial (cells 27 x 84, quads 9 x 35) about affineCenter.
    const double* cellCoeff;const double* faceCoeff;const double* affineCenter;
    const double* farX;const double* farW;int farStride;
    const unsigned long long* keys;const std::size_t* offsets;const double* values;std::size_t pairs;
    // nearCell: plain outer rule of the cell targets of near pairs on a mixed mesh, empty when they take the graded nearOut.
    // innerSwitch: a near target point beyond (1 + innerSwitch) source radii reads the far rule; 0 keeps the near inner everywhere.
    // tangent: remainder rule of the tangent closed form of a warped source on a mixed mesh; empty keeps the radial and fan rules.
    HexRule out,in,pair,pairAffine,nearOut,inSelf,nearCell,tangent;double nearGrade,farInner,innerSwitch,candidate;int order;
    // Far rule by separation on a mixed mesh: a pair at or beyond farRatio[k] of its summed radii takes farLevel[k].
    HexRule farLevel[hexFarLevels];double farRatio[hexFarLevels];const double* farLevelX[hexFarLevels];const double* farLevelW[hexFarLevels];
    // Host types of a mixed mesh (radia_hdiv_typed::HostType); null on a pure hex mesh.
    const int* type;
    // Images of a pure hex mesh; image index img = 1..imageCount, 0 is the direct term.
    const HexImage* images;int imageCount;
};
using radia_hdiv_typed::hexLag3;using radia_hdiv_typed::hexMapX;using radia_hdiv_typed::hexMapJ;using radia_hdiv_typed::quadMapX;using radia_hdiv_typed::quadMapT;
__device__ inline bool hexCell(const HexDevice& d,int g){return g<d.cells;}
__device__ inline const double* hexNodes(const HexDevice& d,int g){
    return g<d.cells?d.cellNodes+81*std::size_t(g):d.faceNodes+27*std::size_t(g-d.cells);}
__device__ inline void hexPoint(const HexDevice& d,int g,const double* xi,double X[3]){
    if(hexCell(d,g))hexMapX(hexNodes(d,g),xi,X);else quadMapX(hexNodes(d,g),xi,X);}
__device__ inline double hexMono(const int* e,const double* xi){
    double v=1.0;for(int k=0;k<3;++k)for(int i=0;i<e[k];++i)v*=xi[k];return v;}
__device__ inline const int* hexExponent(const HexDevice& d,int g,int local){return d.exponents+3*std::size_t(d.first[g]+local);}
// Charge density per unit-cube measure of host g.
__device__ inline double hexDensity(const HexDevice& d,int g,const int* e,const double* xi){
    return d.type?radia_hdiv_typed::density(d.type[g],e,xi):hexMono(e,xi);}
__device__ inline void hexNeumaier(double& sum,double& comp,double x){
    const double t=sum+x;if(fabs(sum)>=fabs(x))comp+=(sum-t)+x;else comp+=(x-t)+sum;sum=t;}
// CPU ImageEvalPoint: T^-1 v, the mirror flips first, then the inverse rotation; img 0 is the identity.
__device__ inline void hexImagePoint(const HexImage& m,const double* v,double* o){
    double x=v[0]*m.flip[0],y=v[1]*m.flip[1];const double z=v[2]*m.flip[2];
    if(m.rotated!=0.0){const double a=x,b=y;x=m.ci*a-m.si*b;y=m.si*a+m.ci*b;}
    o[0]=x;o[1]=y;o[2]=z;
}
__device__ inline void hexImageVector(const HexImage& m,const double* v,double* o){
    double x=v[0],y=v[1];const double z=v[2]*m.flip[2];
    if(m.rotated!=0.0){const double a=x,b=y;x=m.c*a-m.s*b;y=m.s*a+m.c*b;}
    o[0]=x*m.flip[0];o[1]=y*m.flip[1];o[2]=z;
}
__device__ inline void hexImagePoint(const HexDevice& d,int img,const double* v,double* o){
    double x=v[0],y=v[1],z=v[2];
    if(img>0){const HexImage& m=d.images[img-1];x*=m.flip[0];y*=m.flip[1];z*=m.flip[2];
        if(m.rotated!=0.0){const double a=x,b=y;x=m.ci*a-m.si*b;y=m.si*a+m.ci*b;}}
    o[0]=x;o[1]=y;o[2]=z;
}
// CPU ImageApplyVector: T v, the forward rotation, then the mirror flips.
__device__ inline void hexImageVector(const HexDevice& d,int img,const double* v,double* o){
    double x=v[0],y=v[1],z=v[2];
    if(img>0){const HexImage& m=d.images[img-1];
        if(m.rotated!=0.0){const double a=x,b=y;x=m.c*a-m.s*b;y=m.s*a+m.c*b;}
        x*=m.flip[0];y*=m.flip[1];z*=m.flip[2];}
    o[0]=x;o[1]=y;o[2]=z;
}
// True separation of the target from the image of the source, |c_t - T c_s|: the integrand of a directed
// block is 1/|T^-1 x - y| = 1/|x - T y|. For a rotation this differs from |c_t - T^-1 c_s| whenever the centres
// are not symmetric about the sector mid-plane, and a touching pair would be classified far by that reading.
__device__ inline double hexSep(const HexDevice& d,int t,int s,int img=0){
    const double* a=d.center+3*std::size_t(t);double b[3];hexImageVector(d,img,d.center+3*std::size_t(s),b);
    const double dx=a[0]-b[0],dy=a[1]-b[1],dz=a[2]-b[2];return sqrt(dx*dx+dy*dy+dz*dz);}
// Level 0 is the outer rule; a level with no points is absent. The separation is that of the target from the image of the source.
__device__ inline int hexFarLevel(const HexDevice& d,int t,int s,int img=0){
    const double sep=hexSep(d,t,s,img),sz=d.size[t]+d.size[s];int level=0;
    for(int k=0;k<hexFarLevels;++k)if(d.farLevel[k].n&&sep>=d.farRatio[k]*sz)level=k+1;
    return level;
}
__device__ inline const HexRule& hexFarRule(const HexDevice& d,int level){return level?d.farLevel[level-1]:d.out;}
__device__ inline const double* hexFarPoints(const HexDevice& d,int level,int s){const int n=hexFarRule(d,level).n;
    return (level?d.farLevelX[level-1]:d.farX)+3*std::size_t(s)*(n*n*n);}
__device__ inline const double* hexFarWeights(const HexDevice& d,int level,int s){const int n=hexFarRule(d,level).n;
    return (level?d.farLevelW[level-1]:d.farW)+std::size_t(s)*(n*n*n);}
// CPU HexHostsTouch: a shared Q2 lattice node, the source nodes mapped by T, within 1e-9 of the radius sum.
__device__ inline bool hexTouch(const HexDevice& d,int t,int s,int img=0){
    if(img==0&&t==s)return true;
    const double rs=d.size[t]+d.size[s];if(hexSep(d,t,s,img)>1.05*rs+1e-12)return false;
    const double* a=hexNodes(d,t);const double* b=hexNodes(d,s);const int na=hexCell(d,t)?27:9,nb=hexCell(d,s)?27:9;
    const double tol=1e-9*rs+1e-14;
    for(int j=0;j<nb;++j){double q[3];hexImageVector(d,img,b+3*j,q);
        for(int i=0;i<na;++i){const double ex=q[0]-a[3*i],ey=q[1]-a[3*i+1],ez=q[2]-a[3*i+2];
            if(ex*ex+ey*ey+ez*ez<=tol*tol)return true;}}
    return false;
}
// CPU self_pair: the same host, invariant under the image (it lies on the mirror plane or the rotation axis).
__device__ inline bool hexSelfPair(const HexDevice& d,int t,int s,int img){
    if(t!=s)return false;if(img==0)return true;
    const double* c=d.center+3*std::size_t(t);double r[3];hexImagePoint(d,img,c,r);
    return fabs(r[0]-c[0])+fabs(r[1]-c[1])+fabs(r[2]-c[2])<1e-6*d.size[t]+1e-12;
}
struct HexAdjacency {int dim,permT[3],flipT[3],permS[3],flipS[3];};
__device__ inline void hexRefNode(int dim,int node,double* z){
    if(dim==3){z[0]=0.5*(node%3);z[1]=0.5*((node/3)%3);z[2]=0.5*(node/9);}else{z[0]=0.5*(node%3);z[1]=0.5*(node/3);z[2]=0.0;}}
// CPU HexPairAdjacencyOf: the shared entity (vertex 0, edge 1, face 2, cell 3)
// and reference frames that put it at zeta_k = 1 (target) / 0 (source) for
// k >= e with equal leading coordinates; -2 marks a non-conforming contact.
__device__ inline HexAdjacency hexAdjacency(const HexDevice& d,int t,int s,int img=0){
    HexAdjacency adj{-1,{0,1,2},{0,0,0},{0,1,2},{0,0,0}};
    const int dT=hexCell(d,t)?3:2,dS=hexCell(d,s)?3:2,na=dT==3?27:9,nb=dS==3?27:9;
    const double* a=hexNodes(d,t);const double* b=hexNodes(d,s);
    const double tol=1e-9*(d.size[t]+d.size[s])+1e-14;
    unsigned char pi[27],pj[27];int count=0;
    for(int j=0;j<nb;++j){double q[3];hexImageVector(d,img,b+3*j,q);
        for(int i=0;i<na;++i){const double ex=q[0]-a[3*i],ey=q[1]-a[3*i+1],ez=q[2]-a[3*i+2];
            if(ex*ex+ey*ey+ez*ez<=tol*tol){if(count<27){pi[count]=(unsigned char)i;pj[count]=(unsigned char)j;}++count;}}}
    if(!count)return adj;
    int e;
    if(count==1)e=0;else if(count==3)e=1;else if(count==9)e=2;else if(count==27&&dT==3&&dS==3)e=3;else{adj.dim=-2;return adj;}
    const int perms3[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}},perms2[2][3]={{0,1,2},{1,0,2}};
    const int nperT=dT==3?6:2,nperS=dS==3?6:2,nflipT=1<<dT,nflipS=1<<dS;const double eps=1e-12;
    for(int pT=0;pT<nperT;++pT){const int* permT=dT==3?perms3[pT]:perms2[pT];
        for(int fT=0;fT<nflipT;++fT){const int flipT[3]={fT&1,(fT>>1)&1,(fT>>2)&1};bool okT=true;
            for(int q=0;q<count&&okT;++q){double ref[3],zeta[3];hexRefNode(dT,pi[q],ref);
                for(int k=0;k<dT;++k)zeta[k]=flipT[k]?1.0-ref[permT[k]]:ref[permT[k]];
                for(int k=e;k<dT&&okT;++k)okT=fabs(zeta[k]-1.0)<eps;}
            if(!okT)continue;
            for(int pS=0;pS<nperS;++pS){const int* permS=dS==3?perms3[pS]:perms2[pS];
                for(int fS=0;fS<nflipS;++fS){const int flipS[3]={fS&1,(fS>>1)&1,(fS>>2)&1};bool ok=true;
                    for(int q=0;q<count&&ok;++q){double refT[3],zT[3],refS[3],zS[3];hexRefNode(dT,pi[q],refT);hexRefNode(dS,pj[q],refS);
                        for(int k=0;k<dT;++k)zT[k]=flipT[k]?1.0-refT[permT[k]]:refT[permT[k]];
                        for(int k=0;k<dS;++k)zS[k]=flipS[k]?1.0-refS[permS[k]]:refS[permS[k]];
                        for(int k=e;k<dS&&ok;++k)ok=fabs(zS[k])<eps;
                        for(int k=0;k<e&&ok;++k)ok=fabs(zT[k]-zS[k])<eps;}
                    if(!ok)continue;
                    adj.dim=e;for(int k=0;k<3;++k){adj.permT[k]=permT[k];adj.flipT[k]=flipT[k];adj.permS[k]=permS[k];adj.flipS[k]=flipS[k];}
                    return adj;}}}}
    adj.dim=-2;return adj;
}
__device__ inline void hexAxes(const HexDevice& d,int g,double& lo,double& hi){
    const double* n=hexNodes(d,g);const int tip[3]={2,6,18};lo=1e300;hi=0.0;
    for(int a=0;a<(hexCell(d,g)?3:2);++a){const double* p=n+3*tip[a];
        const double L=sqrt((p[0]-n[0])*(p[0]-n[0])+(p[1]-n[1])*(p[1]-n[1])+(p[2]-n[2])*(p[2]-n[2]));lo=fmin(lo,L);hi=fmax(hi,L);}
}
// CPU HexPairExponentialRulesResolve: axis aspect ratios and, off contact,
// the longest extent over the closest lattice gap stay within the bound.
__device__ inline bool hexExponential(const HexDevice& d,int t,int s,bool touching,int img=0){
    double loT,hiT,loS,hiS;hexAxes(d,t,loT,hiT);hexAxes(d,s,loS,hiS);
    double ratio=(loT>0.0&&loS>0.0)?fmax(hiT/loT,hiS/loS):1e300;
    if(!touching){const double* a=hexNodes(d,t);const double* b=hexNodes(d,s);const int na=hexCell(d,t)?27:9,nb=hexCell(d,s)?27:9;double gap2=1e300;
        for(int j=0;j<nb;++j){double q[3];hexImageVector(d,img,b+3*j,q);
            for(int i=0;i<na;++i){const double dx=a[3*i]-q[0],dy=a[3*i+1]-q[1],dz=a[3*i+2]-q[2];gap2=fmin(gap2,dx*dx+dy*dy+dz*dz);}}
        const double gap=sqrt(gap2),extent=fmax(hiT,hiS);ratio=fmax(ratio,gap>0.0?extent/gap:1e300);}
    return ratio<=hexMaxAnisotropy;
}
// Directed (target t, source s, image img) block family of the CPU QuadBlockHex.
__device__ inline int hexClassify(const HexDevice& d,int t,int s,HexAdjacency& adj,int img=0){
    adj.dim=-1;
    if(d.type)return hexTouch(d,t,s,img)||hexSep(d,t,s,img)<=d.nearGrade*(d.size[t]+d.size[s])?hexGradedClass:hexFarClass;
    const bool touching=hexTouch(d,t,s,img);const double sep=hexSep(d,t,s,img),sz=d.size[t]+d.size[s];
    const bool exactNear=touching||sep<=hexAffineExactNear*sz;
    if(d.order==1){
        bool ok=true;if(touching||sep<=d.nearGrade*sz)ok=hexExponential(d,t,s,touching,img);
        if(touching&&ok){adj=hexAdjacency(d,t,s,img);if(adj.dim>=0)return hexDuffyClass;}
        if(ok&&!touching&&sep<=d.nearGrade*sz)return hexProductClass;
    }
    if(d.affine[t]&&d.affine[s]&&(d.order==2||!exactNear))return exactNear?hexAffineClass:hexFarClass;
    if(!touching&&sep>hexDistortedFar*sz)return hexFarClass;
    if(d.order==2)return hexCompositeClass;
    // The CPU sub-simplex pair rule lies between the near band and the distorted far switch; the staged bands leave it empty.
    return (img==0&&t==s)||touching||sep<=d.nearGrade*sz?hexGradedClass:hexSubPairClass;
}
// Far tensor rule of source host s (CPU GetHexFarRule / PhiInnerHexTensorVec).
__device__ inline void hexFarInner(const HexDevice& d,int s,const double* p,int nS,double* inner){
    const int n1=d.out.n,np=hexCell(d,s)?n1*n1*n1:n1*n1;const double* X=d.farX+3*std::size_t(s)*d.farStride;const double* W=d.farW+std::size_t(s)*d.farStride;
    for(int q=0;q<np;++q){const double dx=p[0]-X[3*q],dy=p[1]-X[3*q+1],dz=p[2]-X[3*q+2];const double r=sqrt(dx*dx+dy*dy+dz*dz);
        if(r<1e-300)continue;const double w=W[q]/r;
        const double xi[3]={d.out.x[q%n1],d.out.x[(q/n1)%n1],hexCell(d,s)?d.out.x[q/(n1*n1)]:0.0};
        for(int l=0;l<nS;++l)inner[l]+=w*hexDensity(d,s,hexExponent(d,s,l),xi);}
}
// The far rule's source densities, once per block: every target point reads the same table.
__device__ inline void hexFarDensities(const HexDevice& d,int s,int nS,double* table,int level=0){
    const HexRule& R=hexFarRule(d,level);const int n1=R.n,np=hexCell(d,s)?n1*n1*n1:n1*n1;
    for(int q=threadIdx.x;q<np;q+=blockDim.x){const double xi[3]={R.x[q%n1],R.x[(q/n1)%n1],hexCell(d,s)?R.x[q/(n1*n1)]:0.0};
        for(int l=0;l<nS;++l)table[std::size_t(q)*nS+l]=hexDensity(d,s,hexExponent(d,s,l),xi);}
    __syncthreads();
}
__device__ inline void hexFarInnerTable(const HexDevice& d,int s,const double* p,int nS,const double* table,double* inner,int level=0){
    const int n1=hexFarRule(d,level).n,np=hexCell(d,s)?n1*n1*n1:n1*n1;const double* X=hexFarPoints(d,level,s);const double* W=hexFarWeights(d,level,s);
    for(int q=0;q<np;++q){const double dx=p[0]-X[3*q],dy=p[1]-X[3*q+1],dz=p[2]-X[3*q+2];const double r=sqrt(dx*dx+dy*dy+dz*dz);
        if(r<1e-300)continue;const double w=W[q]/r;const double* row=table+std::size_t(q)*nS;
        for(int l=0;l<nS;++l)inner[l]+=w*row[l];}
}
// Exact potentials of an affine source from its outer faces (CPU PhiInnerHexAffine{Cell,Face}Vec). BDM2 reads the
// host coefficient tables; a BDM1 graded pair passes the coefficients of its nS monomials in local.
__device__ inline void hexAffineInner(const HexDevice& d,int s,const double* p,int nS,double* inner,const double* local=nullptr){
    const double* nd=hexNodes(d,s);const double* c=d.affineCenter+3*std::size_t(s);const double lp[3]={p[0]-c[0],p[1]-c[1],p[2]-c[2]};
    if(hexCell(d,s)){
        // The corners as Q2 lattice nodes, one face in order and then the opposite one.
        const int corner[8]={0,2,8,6,18,20,26,24};
        // A BDM1 charge has degree at most 3 on a cell and 2 on a face; its other polynomial coefficients are zero.
        const int used=d.order==2?84:20;
        double C[8][3];for(int i=0;i<8;++i)for(int k=0;k<3;++k)C[i][k]=nd[3*corner[i]+k]-c[k];
        double moments[84];if(d.order==2)rad_hdiv_cuda::HexPotentialMomentsUpTo6(C,lp,moments);else rad_hdiv_cuda::HexPotentialMomentsUpTo3(C,lp,moments);
        for(int l=0;l<nS;++l){const int* e=hexExponent(d,s,l);const double* coeff=local?local+std::size_t(l)*84:d.cellCoeff+(std::size_t(s)*27+e[0]+3*e[1]+9*e[2])*84;
            double v=0.0;for(int i=0;i<used;++i)v+=coeff[i]*moments[i];inner[l]+=v;}
        return;
    }
    const int corner[4]={0,2,8,6};const int f=s-d.cells,used=d.order==2?35:10;
    double V[4][3];for(int i=0;i<4;++i)for(int k=0;k<3;++k)V[i][k]=nd[3*corner[i]+k]-c[k];
    double moments[35];if(d.order==2)rad_hdiv_cuda::QuadPotentialMomentsUpTo4(V,lp,moments);else rad_hdiv_cuda::QuadPotentialMomentsUpTo2(V,lp,moments);
    for(int l=0;l<nS;++l){const int* e=hexExponent(d,s,l);const double* coeff=local?local+std::size_t(l)*35:d.faceCoeff+(std::size_t(f)*9+e[0]+3*e[1])*35;
        double v=0.0;for(int i=0;i<used;++i)v+=coeff[i]*moments[i];inner[l]+=v;}
}
// CPU PhiInnerHexRadialHostVec: radial cones from the reference point over
// the complete host (six square faces of the cube, four edges of the quad).
__device__ inline void hexRadialInner(const HexDevice& d,int s,const double* p,const double* ref,int nS,double* inner){
    const double* nd=hexNodes(d,s);const HexRule& R=d.in;const int nq=R.n;
    if(hexCell(d,s)){
        for(int fixed=0;fixed<3;++fixed){const int aa=(fixed+1)%3,ab=(fixed+2)%3;
            for(int side=0;side<2;++side){const double distance=fabs(double(side)-ref[fixed]);if(distance<1e-300)continue;
                for(int ia=0;ia<nq;++ia)for(int ib=0;ib<nq;++ib){double face[3]={0.0,0.0,0.0};face[fixed]=double(side);face[aa]=R.x[ia];face[ab]=R.x[ib];
                    for(int ir=0;ir<nq;++ir){const double rp=R.x[ir];double sr[3];for(int k=0;k<3;++k)sr[k]=ref[k]+rp*(face[k]-ref[k]);
                        double X[3];hexMapX(nd,sr,X);const double dx=p[0]-X[0],dy=p[1]-X[1],dz=p[2]-X[2];const double r=sqrt(dx*dx+dy*dy+dz*dz);
                        if(r<1e-300)continue;const double w=R.w[ia]*R.w[ib]*R.w[ir]*rp*rp*distance/r;
                        for(int l=0;l<nS;++l)inner[l]+=w*hexMono(hexExponent(d,s,l),sr);}}}}
        return;
    }
    for(int fixed=0;fixed<2;++fixed){const int free=1-fixed;
        for(int side=0;side<2;++side)for(int it=0;it<nq;++it){double edge[2]={0.0,0.0};edge[fixed]=double(side);edge[free]=R.x[it];
            const double tangent[2]={free==0?1.0:0.0,free==1?1.0:0.0},ev[2]={edge[0]-ref[0],edge[1]-ref[1]};
            const double area=fabs(ev[0]*tangent[1]-ev[1]*tangent[0]);if(area<1e-300)continue;
            for(int ir=0;ir<nq;++ir){const double rp=R.x[ir];const double sr[3]={ref[0]+rp*ev[0],ref[1]+rp*ev[1],0.0};
                double X[3];quadMapX(nd,sr,X);const double dx=p[0]-X[0],dy=p[1]-X[1],dz=p[2]-X[2];const double r=sqrt(dx*dx+dy*dy+dz*dz);
                if(r<1e-300)continue;const double w=R.w[it]*R.w[ir]*rp*area/r;
                for(int l=0;l<nS;++l)inner[l]+=w*hexMono(hexExponent(d,s,l),sr);}}}
}
__device__ inline double hexClamp01(double v){return fmin(fmax(v,0.0),1.0);}
// CPU HexQ2Inverse: Newton on the Q2 map. A converged point stays unclamped; any other final point is clamped to the cube.
__device__ inline bool hexQ2Newton(const double* nodes,const double* target,double ref[3],bool& converged){
    converged=false;ref[0]=ref[1]=ref[2]=0.5;double last[3]={0.5,0.5,0.5};bool valid=false;double scale=0.0;
    for(int n=1;n<27;++n){const double dx=nodes[3*n]-nodes[0],dy=nodes[3*n+1]-nodes[1],dz=nodes[3*n+2]-nodes[2];scale=fmax(scale,sqrt(dx*dx+dy*dy+dz*dz));}
    const double tolerance=2e-11*fmax(scale,1.0)+1e-13;
    for(int it=0;it<20;++it){double X[3],J[3][3];hexMapJ(nodes,ref,X,J);
        const double res[3]={X[0]-target[0],X[1]-target[1],X[2]-target[2]};
        if(sqrt(res[0]*res[0]+res[1]*res[1]+res[2]*res[2])<=tolerance){converged=true;return true;}
        const double det=J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
        if(!isfinite(det)||fabs(det)<1e-18){if(!valid)return false;break;}
        valid=true;
        const double I[3][3]={{(J[1][1]*J[2][2]-J[1][2]*J[2][1])/det,(J[0][2]*J[2][1]-J[0][1]*J[2][2])/det,(J[0][1]*J[1][2]-J[0][2]*J[1][1])/det},
                              {(J[1][2]*J[2][0]-J[1][0]*J[2][2])/det,(J[0][0]*J[2][2]-J[0][2]*J[2][0])/det,(J[0][2]*J[1][0]-J[0][0]*J[1][2])/det},
                              {(J[1][0]*J[2][1]-J[1][1]*J[2][0])/det,(J[0][1]*J[2][0]-J[0][0]*J[2][1])/det,(J[0][0]*J[1][1]-J[0][1]*J[1][0])/det}};
        for(int row=0;row<3;++row)ref[row]-=I[row][0]*res[0]+I[row][1]*res[1]+I[row][2]*res[2];
        if(!isfinite(ref[0])||!isfinite(ref[1])||!isfinite(ref[2])){for(int k=0;k<3;++k)ref[k]=last[k];break;}
        for(int k=0;k<3;++k)last[k]=ref[k];
    }
    for(int k=0;k<3;++k)ref[k]=hexClamp01(ref[k]);
    return valid;
}
// The composite near rule clamps every inverse, converged or not.
__device__ inline bool hexQ2Inverse(const double* nodes,const double* target,double ref[3]){
    bool converged;const bool valid=hexQ2Newton(nodes,target,ref,converged);
    for(int k=0;k<3;++k)ref[k]=hexClamp01(ref[k]);
    return valid;
}
__device__ inline double quadDistance2(const double* nodes,const double* target,const double* uv){
    double X[3];quadMapX(nodes,uv,X);const double dx=X[0]-target[0],dy=X[1]-target[1],dz=X[2]-target[2];return dx*dx+dy*dy+dz*dz;}
// CPU QuadQ2ClosestReference: coarse grid, then monotone Gauss-Newton with a
// projected-gradient fallback on the unit square.
__device__ inline bool quadClosest(const double* nodes,const double* target,double ref[2]){
    double best=1e300;
    for(int iv=0;iv<=4;++iv)for(int iu=0;iu<=4;++iu){const double c[2]={0.25*iu,0.25*iv};const double dd=quadDistance2(nodes,target,c);
        if(isfinite(dd)&&dd<best){best=dd;ref[0]=c[0];ref[1]=c[1];}}
    if(!(best<1e300))return false;
    for(int it=0;it<40;++it){double X[3],T[3][2];quadMapT(nodes,ref,X,T);
        const double res[3]={X[0]-target[0],X[1]-target[1],X[2]-target[2]};
        const double current=res[0]*res[0]+res[1]*res[1]+res[2]*res[2];
        const double n00=T[0][0]*T[0][0]+T[1][0]*T[1][0]+T[2][0]*T[2][0],n01=T[0][0]*T[0][1]+T[1][0]*T[1][1]+T[2][0]*T[2][1],
                     n11=T[0][1]*T[0][1]+T[1][1]*T[1][1]+T[2][1]*T[2][1];
        const double r0=T[0][0]*res[0]+T[1][0]*res[1]+T[2][0]*res[2],r1=T[0][1]*res[0]+T[1][1]*res[1]+T[2][1]*res[2];
        const double det=n00*n11-n01*n01;if(!isfinite(det)||det<=1e-30)return false;
        const double d0=(n11*r0-n01*r1)/det,d1=(n00*r1-n01*r0)/det;bool accepted=false;double alpha=1.0;
        for(int ls=0;ls<16;++ls,alpha*=0.5){const double trial[2]={hexClamp01(ref[0]-alpha*d0),hexClamp01(ref[1]-alpha*d1)};
            const double s0=trial[0]-ref[0],s1=trial[1]-ref[1];if(sqrt(s0*s0+s1*s1)<=2e-12)return true;
            const double td=quadDistance2(nodes,target,trial);
            if(isfinite(td)&&td<=current+1e-15*fmax(current,1e-30)){ref[0]=trial[0];ref[1]=trial[1];accepted=true;
                if(current-td<=1e-13*fmax(current,1e-30))return true;break;}}
        if(!accepted){const double metric=fmax(n00+n11,1e-30);alpha=1.0;
            for(int ls=0;ls<24;++ls,alpha*=0.5){const double trial[2]={hexClamp01(ref[0]-alpha*r0/metric),hexClamp01(ref[1]-alpha*r1/metric)};
                const double s0=trial[0]-ref[0],s1=trial[1]-ref[1];if(sqrt(s0*s0+s1*s1)<=2e-12)return true;
                const double td=quadDistance2(nodes,target,trial);
                if(isfinite(td)&&td<=current+1e-15*fmax(current,1e-30)){ref[0]=trial[0];ref[1]=trial[1];accepted=true;
                    if(current-td<=1e-13*fmax(current,1e-30))return true;break;}}}
        if(!accepted)return true;
    }
    return true;
}
// CPU HexQ2ClosestReference: a converged interior inverse is the closest point; otherwise the nearest of the six Q2 faces.
__device__ inline bool hexQ2Closest(const double* nodes,const double* target,double ref[3]){
    bool converged;double inverse[3];const bool valid=hexQ2Newton(nodes,target,inverse,converged);
    const double slack=1e-9;bool inside=valid;for(int k=0;k<3;++k)inside=inside&&inverse[k]>=-slack&&inverse[k]<=1.0+slack;
    if(inside){for(int k=0;k<3;++k)ref[k]=hexClamp01(inverse[k]);return true;}
    double best=1e300;bool found=false;
    for(int fixed=0;fixed<3;++fixed){const int aa=(fixed+1)%3,ab=(fixed+2)%3;
        for(int side=0;side<2;++side){double face[27];
            for(int iv=0;iv<3;++iv)for(int iu=0;iu<3;++iu){int index[3];index[fixed]=2*side;index[aa]=iu;index[ab]=iv;
                const double* from=nodes+3*(index[0]+3*index[1]+9*index[2]);double* to=face+3*(iu+3*iv);to[0]=from[0];to[1]=from[1];to[2]=from[2];}
            double uv[2]={0.5,0.5};if(!quadClosest(face,target,uv))continue;
            const double dd=quadDistance2(face,target,uv);
            if(isfinite(dd)&&dd<best){best=dd;found=true;ref[fixed]=double(side);ref[aa]=uv[0];ref[ab]=uv[1];}}}
    return found;
}
// One-sided sinh nodes on [0,1] for a peak of width b at 0; a peak too narrow or too wide takes the plain fine rule.
__device__ inline void hexSinhOneSided(const HexRule& R,double b,double* x,double* w){
    if(b>=hexSinhMinWidth&&b<hexSinhMaxWidth){const double mu=0.5*asinh(1.0/b);
        for(int i=0;i<R.n;++i){const double arg=2.0*mu*R.x[i];x[i]=b*sinh(arg);w[i]=R.w[i]*2.0*b*mu*cosh(arg);}}
    else for(int i=0;i<R.n;++i){x[i]=R.x[i];w[i]=R.w[i];}
}
// Two-sided sinh nodes for a peak of width b at t0; only a peak wider than the interval takes the coarse plain rule.
__device__ inline int hexSinhTwoSided(const HexRule& coarse,const HexRule& R,double t0,double b,double* x,double* w){
    if(b>0.0&&b<hexSinhMaxWidth){const double s0=asinh(t0/b),s1=asinh((1.0-t0)/b),mu=0.5*(s0+s1),eta=0.5*(s0-s1);
        for(int i=0;i<R.n;++i){const double arg=mu*(2.0*R.x[i]-1.0)-eta;x[i]=t0+b*sinh(arg);w[i]=R.w[i]*2.0*b*mu*cosh(arg);}
        return R.n;}
    for(int i=0;i<coarse.n;++i){x[i]=coarse.x[i];w[i]=coarse.w[i];}
    return coarse.n;
}
__device__ inline double hexDot(const double* a,const double* b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
// CPU PhiInnerHexConeFanVec: cones from the reference apex over the complete host, a four-fan polar rule on
// every cone base, with sinh substitutions on the ray, the fan radius and the edge parameter.
__device__ inline void hexConeFanInner(const HexDevice& d,int s,const double* p,const double* ref,int nS,double* inner){
    const double* nd=hexNodes(d,s);const HexRule& fine=d.inSelf;
    double tn[hexRuleCap],tw[hexRuleCap],sn[hexRuleCap],sw[hexRuleCap],rn[hexRuleCap],rw[hexRuleCap],apex[3];
    if(hexCell(d,s)){double J[3][3];hexMapJ(nd,ref,apex,J);
        const double dpx=p[0]-apex[0],dpy=p[1]-apex[1],dpz=p[2]-apex[2],dist=sqrt(dpx*dpx+dpy*dpy+dpz*dpz);
        for(int fixed=0;fixed<3;++fixed){const int aa=(fixed+1)%3,ab=(fixed+2)%3;
            for(int side=0;side<2;++side){const double nsgn=double(side)-ref[fixed],fd=fabs(nsgn);if(fd<1e-300)continue;
                double N3[3],A[3],B[3];for(int k=0;k<3;++k){N3[k]=J[k][fixed]*nsgn;A[k]=J[k][aa];B[k]=J[k][ab];}
                const double a0=ref[aa],b0=ref[ab],AA=hexDot(A,A),AB=hexDot(A,B),BB=hexDot(B,B),AN=hexDot(A,N3),BN=hexDot(B,N3),det=AA*BB-AB*AB;
                double ca=a0,cb=b0;if(det>1e-24*AA*BB){ca=hexClamp01(a0-(BB*AN-AB*BN)/det);cb=hexClamp01(b0-(AA*BN-AB*AN)/det);}
                double G0[3];for(int k=0;k<3;++k)G0[k]=N3[k]+(ca-a0)*A[k]+(cb-b0)*B[k];
                const double peak=sqrt(hexDot(G0,G0)+dist*dist);
                for(int fixed2=0;fixed2<2;++fixed2)for(int side2=0;side2<2;++side2){
                    const double cFixed=fixed2==0?ca:cb,cFree=fixed2==0?cb:ca,s2=fabs(double(side2)-cFixed);if(s2<1e-300)continue;
                    const double* Efix=fixed2==0?A:B;const double* Efree=fixed2==0?B:A;
                    double Gedge[3];for(int k=0;k<3;++k)Gedge[k]=G0[k]+(double(side2)-cFixed)*Efix[k];
                    const double EE=hexDot(Efree,Efree);double t0=cFree;if(EE>0.0)t0=hexClamp01(cFree-hexDot(Gedge,Efree)/EE);
                    double R[3];for(int k=0;k<3;++k)R[k]=Gedge[k]+(t0-cFree)*Efree[k];
                    const double bt=EE>0.0?sqrt(hexDot(R,R)+dist*dist)/sqrt(EE):0.0;
                    const int nt=hexSinhTwoSided(d.in,fine,t0,bt,tn,tw);
                    for(int it=0;it<nt;++it){const double t=tn[it],wa=fixed2==0?double(side2)-ca:t-ca,wb=fixed2==0?t-cb:double(side2)-cb;
                        double JW[3];for(int k=0;k<3;++k)JW[k]=wa*A[k]+wb*B[k];const double jw=sqrt(hexDot(JW,JW));
                        hexSinhOneSided(fine,jw>0.0?peak/jw:0.0,sn,sw);
                        for(int is=0;is<fine.n;++is){const double rho=sn[is],wrho=tw[it]*sw[is]*rho*s2,a=ca+rho*wa,b=cb+rho*wb;
                            double ray[3];ray[fixed]=nsgn;ray[aa]=a-a0;ray[ab]=b-b0;
                            double v[3];for(int k=0;k<3;++k)v[k]=N3[k]+(a-a0)*A[k]+(b-b0)*B[k];const double vn=sqrt(hexDot(v,v));
                            hexSinhOneSided(fine,vn>0.0&&dist>0.0?dist/vn:0.0,rn,rw);
                            for(int ir=0;ir<fine.n;++ir){const double r=rn[ir];double y[3],X[3];for(int k=0;k<3;++k)y[k]=ref[k]+r*ray[k];
                                hexMapX(nd,y,X);const double dx=p[0]-X[0],dy=p[1]-X[1],dz=p[2]-X[2],pr=sqrt(dx*dx+dy*dy+dz*dz);if(pr<1e-300)continue;
                                const double weight=wrho*rw[ir]*r*r*fd/pr;for(int l=0;l<nS;++l)inner[l]+=weight*hexMono(hexExponent(d,s,l),y);}}}}}}
        return;}
    double T[3][2];quadMapT(nd,ref,apex,T);
    const double dpx=p[0]-apex[0],dpy=p[1]-apex[1],dpz=p[2]-apex[2],dist=sqrt(dpx*dpx+dpy*dpy+dpz*dpz);
    for(int fixed2=0;fixed2<2;++fixed2){const int free2=1-fixed2;
        for(int side2=0;side2<2;++side2){const double s2=fabs(double(side2)-ref[fixed2]);if(s2<1e-300)continue;
            double Efree[3],Gedge[3];for(int k=0;k<3;++k){Efree[k]=T[k][free2];Gedge[k]=(double(side2)-ref[fixed2])*T[k][fixed2];}
            const double EE=hexDot(Efree,Efree);double t0=ref[free2];if(EE>0.0)t0=hexClamp01(ref[free2]-hexDot(Gedge,Efree)/EE);
            double R[3];for(int k=0;k<3;++k)R[k]=Gedge[k]+(t0-ref[free2])*Efree[k];
            const double bt=EE>0.0?sqrt(hexDot(R,R)+dist*dist)/sqrt(EE):0.0;
            const int nt=hexSinhTwoSided(d.in,fine,t0,bt,tn,tw);
            for(int it=0;it<nt;++it){double w2[2];w2[fixed2]=double(side2)-ref[fixed2];w2[free2]=tn[it]-ref[free2];
                double Jw[3];for(int k=0;k<3;++k)Jw[k]=w2[0]*T[k][0]+w2[1]*T[k][1];const double jw=sqrt(hexDot(Jw,Jw));
                hexSinhOneSided(fine,jw>0.0&&dist>0.0?dist/jw:0.0,sn,sw);
                for(int is=0;is<fine.n;++is){const double rho=sn[is],uv[2]={ref[0]+rho*w2[0],ref[1]+rho*w2[1]};double X[3];
                    quadMapX(nd,uv,X);const double dx=p[0]-X[0],dy=p[1]-X[1],dz=p[2]-X[2],pr=sqrt(dx*dx+dy*dy+dz*dz);if(pr<1e-300)continue;
                    const double weight=tw[it]*sw[is]*rho*s2/pr,y3[3]={uv[0],uv[1],0.0};
                    for(int l=0;l<nS;++l)inner[l]+=weight*hexMono(hexExponent(d,s,l),y3);}}}}
}
// CPU QuadBlockHexNearTensor inner: the exact analytic potential on an affine source, the tangent rule on a warped
// source of a mixed mesh, otherwise the cone/fan rule, its apex the outer point itself (self) or the source's closest
// point. A self pair stays on the unmapped outer point p; every other pair evaluates at the image-mapped point pm.
__device__ inline bool hexGradedInner(const HexDevice& d,bool self,int s,const double* p,const double* pm,const double* xi,int nS,double* inner,const double* local){
    if(d.affine[s]){hexAffineInner(d,s,self?p:pm,nS,inner,local);return true;}
    if(d.tangent.n){const HexRule& Q=d.type[s]==radia_hdiv_typed::quadHost?d.inSelf:d.tangent;return radia_hdiv_typed::warpedPotentials(d.type[s],hexNodes(d,s),self?p:pm,self?xi:nullptr,hexExponent(d,s,0),nS,d.order,Q.x,Q.w,Q.n,inner);}
    if(self){hexConeFanInner(d,s,p,xi,nS,inner);return true;}
    double ref[3]={0.5,0.5,0.0};
    if(!(hexCell(d,s)?hexQ2Closest(hexNodes(d,s),pm,ref):quadClosest(hexNodes(d,s),pm,ref)))return false;
    hexConeFanInner(d,s,pm,ref,nS,inner);return true;
}
// CPU QuadBlockTyped inner of a tet, triangle or pyramid source at the target point p, mapped by T^-1 for an image: exact on the
// polynomial charges of an affine source, by its Taylor polynomial or on levels on its rational charge, by the tangent rule or radial
// from the anchor on a non-parallelogram pyramid. The target's own pre-image is used only when p is that point on the source itself.
__device__ inline bool hexTypedInner(const HexDevice& d,bool onSource,int s,const double* p,const double* xi,bool analytic,const double* center,const double* coeff,int nS,double* inner){
    using namespace radia_hdiv_typed;
    const int type=d.type[s];const double* nd=hexNodes(d,s);const int* e=hexExponent(d,s,0);const HexRule& R=d.inSelf;const HexRule& Q=d.tangent;
    if(analytic)affinePotentials(type,nd,center,p,nS,coeff,inner);
    if(type!=pyramidHost)return true;
    if(analytic){
        if(d.tangent.n)straightPyramidRational(nd,p,e,nS,Q.x,Q.w,Q.n,inner);
        else for(int l=0;l<nS;++l)if(e[3*l+2]<0)inner[l]+=pyramidRationalPotential(nd,center,p,R.x,R.w,R.n);
        return true;}
    if(d.tangent.n)return warpedPotentials(type,nd,p,onSource?xi:nullptr,e,nS,d.order,Q.x,Q.w,Q.n,inner);
    double anchor[3];if(onSource)collapse(type,xi,anchor);else pyramidAnchor(nd,p,anchor);
    pyramidRadial(nd,p,anchor,e,0,nS,R.x,R.w,R.n,inner);return true;
}
// CPU QuadBlockHexCompositeProduct inner: the exact self radial, the radial
// rule anchored at the (clamped) source pre-image near the source, the outer
// tensor rule beyond 1.5 far_inner host radii. A self pair stays on the unmapped point p, the others take pm.
__device__ inline bool hexCompositeInner(const HexDevice& d,bool self,int s,const double* p,const double* pm,const double* xi,int nS,double* inner){
    if(self){hexRadialInner(d,s,p,xi,nS,inner);return true;}
    const double* c=d.center+3*std::size_t(s);
    const double dx=pm[0]-c[0],dy=pm[1]-c[1],dz=pm[2]-c[2];
    if(sqrt(dx*dx+dy*dy+dz*dz)>1.5*d.farInner*d.size[s]){hexFarInner(d,s,pm,nS,inner);return true;}
    double ref[3]={0.5,0.5,0.5};
    if(!(hexCell(d,s)?hexQ2Inverse(hexNodes(d,s),pm,ref):quadClosest(hexNodes(d,s),pm,ref)))return false;
    hexRadialInner(d,s,pm,ref,nS,inner);return true;
}
// Host-level geometry (CPU constructor): corner mean and corner radius, the
// affine test of HexAffineInverseForms / QuadAffineInverseForms and, for
// BDM2, the physical monomial polynomials about the reference-centre image.
__device__ inline int hexPolyIdx(int ax,int ay,int az){
    const int deg=ax+ay+az;int idx=0;for(int k=0;k<deg;++k)idx+=(k+1)*(k+2)/2;for(int a=0;a<ax;++a)idx+=deg-a+1;return idx+ay;}
__device__ inline void hexPolyMulLinear(double* poly,int& deg,const double lin[4],int ncoeff){
    double tmp[84];for(int i=0;i<84;++i)tmp[i]=0.0;
    for(int total=0;total<=deg;++total)for(int ax=0;ax<=total;++ax)for(int ay=0;ay<=total-ax;++ay){const int az=total-ax-ay;
        const double c=poly[hexPolyIdx(ax,ay,az)];if(c==0.0)continue;
        tmp[hexPolyIdx(ax,ay,az)]+=c*lin[0];tmp[hexPolyIdx(ax+1,ay,az)]+=c*lin[1];tmp[hexPolyIdx(ax,ay+1,az)]+=c*lin[2];tmp[hexPolyIdx(ax,ay,az+1)]+=c*lin[3];}
    ++deg;for(int i=0;i<ncoeff;++i)poly[i]=tmp[i];
}
// Affine test and inverse linear forms of one host (CPU HexAffineInverseForms / QuadAffineInverseForms).
__device__ inline bool hexAffineForms(const double* nd,bool cell,double lin[3][4],double& scaleInv){
    bool ok=false;
    if(cell){const double* o=nd;const double* px=nd+6;const double* py=nd+18;const double* pz=nd+54;double A[3][3];
        for(int k=0;k<3;++k){A[k][0]=px[k]-o[k];A[k][1]=py[k]-o[k];A[k][2]=pz[k]-o[k];}
        const double det=A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
        double scale=0.0;for(int j=0;j<3;++j){double n2=0.0;for(int k=0;k<3;++k)n2+=A[k][j]*A[k][j];scale=fmax(scale,sqrt(n2));}
        ok=!(fabs(det)<1e-300||scale<1e-300);const double tol=1e-10*scale+1e-12;
        for(int iz=0;iz<3&&ok;++iz)for(int iy=0;iy<3&&ok;++iy)for(int ix=0;ix<3&&ok;++ix){const double* p=nd+3*(ix+3*iy+9*iz);double e2=0.0;
            for(int k=0;k<3;++k){const double v=p[k]-(o[k]+0.5*ix*A[k][0]+0.5*iy*A[k][1]+0.5*iz*A[k][2]);e2+=v*v;}ok=sqrt(e2)<=tol;}
        if(ok){const double id=1.0/det;double B[3][3];
            B[0][0]=(A[1][1]*A[2][2]-A[1][2]*A[2][1])*id;B[0][1]=-(A[0][1]*A[2][2]-A[0][2]*A[2][1])*id;B[0][2]=(A[0][1]*A[1][2]-A[0][2]*A[1][1])*id;
            B[1][0]=-(A[1][0]*A[2][2]-A[1][2]*A[2][0])*id;B[1][1]=(A[0][0]*A[2][2]-A[0][2]*A[2][0])*id;B[1][2]=-(A[0][0]*A[1][2]-A[0][2]*A[1][0])*id;
            B[2][0]=(A[1][0]*A[2][1]-A[1][1]*A[2][0])*id;B[2][1]=-(A[0][0]*A[2][1]-A[0][1]*A[2][0])*id;B[2][2]=(A[0][0]*A[1][1]-A[0][1]*A[1][0])*id;
            for(int q=0;q<3;++q){lin[q][0]=-(B[q][0]*o[0]+B[q][1]*o[1]+B[q][2]*o[2]);lin[q][1]=B[q][0];lin[q][2]=B[q][1];lin[q][3]=B[q][2];}
            scaleInv=1.0/fabs(det);}}
    else{const double* o=nd;const double* pu=nd+6;const double* pv=nd+18;double a[3],b[3];
        for(int k=0;k<3;++k){a[k]=pu[k]-o[k];b[k]=pv[k]-o[k];}
        const double aa=a[0]*a[0]+a[1]*a[1]+a[2]*a[2],ab=a[0]*b[0]+a[1]*b[1]+a[2]*b[2],bb=b[0]*b[0]+b[1]*b[1]+b[2]*b[2],det=aa*bb-ab*ab;
        ok=det>1e-300;const double tol=1e-10*fmax(sqrt(aa),sqrt(bb))+1e-12;
        for(int j=0;j<3&&ok;++j)for(int i=0;i<3&&ok;++i){const double* q=nd+3*(i+3*j);double e2=0.0;
            for(int k=0;k<3;++k){const double v=q[k]-(o[k]+0.5*i*a[k]+0.5*j*b[k]);e2+=v*v;}ok=e2<=tol*tol;}
        if(ok){const double i00=bb/det,i01=-ab/det,i11=aa/det;
            for(int k=0;k<3;++k){lin[0][k+1]=i00*a[k]+i01*b[k];lin[1][k+1]=i01*a[k]+i11*b[k];}
            lin[0][0]=-(lin[0][1]*o[0]+lin[0][2]*o[1]+lin[0][3]*o[2]);lin[1][0]=-(lin[1][1]*o[0]+lin[1][2]*o[1]+lin[1][3]*o[2]);
            scaleInv=1.0/sqrt(det);}}
    return ok;
}
// Physical polynomial of one reference monomial about the host centre.
__device__ inline bool hexMonomialCoeff(const double lin[3][4],double scaleInv,bool cell,const int* e,double* dst){
    const int axes=cell?3:2,ncoeff=cell?84:35;
    double poly[84];for(int i=0;i<84;++i)poly[i]=0.0;poly[0]=1.0;int deg=0;
    for(int a=0;a<axes;++a)for(int r=0;r<e[a];++r)hexPolyMulLinear(poly,deg,lin[a],ncoeff);
    bool finite=true;for(int i=0;i<ncoeff;++i){dst[i]=scaleInv*poly[i];finite=finite&&isfinite(dst[i]);}
    return finite;
}
__device__ inline void hexCenterForms(const double* nd,bool cell,double lin[3][4],double ac[3]){
    const double mid[3]={0.5,0.5,0.5};if(cell)hexMapX(nd,mid,ac);else quadMapX(nd,mid,ac);
    for(int a=0;a<(cell?3:2);++a)lin[a][0]+=lin[a][1]*ac[0]+lin[a][2]*ac[1]+lin[a][3]*ac[2];
}
// Beyond the corner radius a point lies outside a straight source, where the far tensor rule converges.
__device__ inline bool hexBeyondSource(const HexDevice& d,int s,const double* p){
    const double* c=d.center+3*std::size_t(s);const double dx=p[0]-c[0],dy=p[1]-c[1],dz=p[2]-c[2],limit=(1.0+d.innerSwitch)*d.size[s];
    return dx*dx+dy*dy+dz*dz>=limit*limit;
}
// Phase 1 tabulates, per target point, the inner source vector and the
// target weights; phase 2 sums every block entry over the points in the CPU
// accumulation order (Neumaier for the pair product, as the CPU).
__device__ inline void hexTensorBlock(const HexDevice& d,int t,int s,int img,int cls,int nT,int nS,double* scratch,double* out,int* invalid){
    const bool cellT=hexCell(d,t),cellS=hexCell(d,s),product=cls==hexProductClass,graded=cls==hexGradedClass,self=hexSelfPair(d,t,s,img);
    const int farLevel=cls==hexFarClass&&d.type?hexFarLevel(d,t,s,img):0;
    const HexRule& R=product?(d.affine[t]&&d.affine[s]?d.pairAffine:d.pair):graded?(cellT&&d.nearCell.n?d.nearCell:d.nearOut):hexFarRule(d,farLevel);
    const int N=R.n,nc=cellT?N:1,npT=N*N*nc;int npS=0;
    double* source=scratch;
    if(product){const int ncs=cellS?N:1;npS=N*N*ncs;
        for(int q=threadIdx.x;q<npS;q+=blockDim.x){const int a=q/(N*ncs),b=(q/ncs)%N,c=q%ncs;
            const double eta[3]={R.x[a],R.x[b],cellS?R.x[c]:0.0};double* row=source+std::size_t(q)*(4+nS);
            hexPoint(d,s,eta,row);row[3]=R.w[a]*R.w[b]*(cellS?R.w[c]:1.0);
            for(int l=0;l<nS;++l)row[4+l]=hexMono(hexExponent(d,s,l),eta);}
        __syncthreads();}
    // A BDM1 graded pair on an affine source keeps its monomial coefficients for this block only.
    const double* local=nullptr;
    const bool typedSource=graded&&d.type&&d.type[s]>=radia_hdiv_typed::tetHost;bool analytic=false;double typedCenter[3]={0.0,0.0,0.0};
    if(typedSource){using namespace radia_hdiv_typed;double lin[3][4],scaleInv=0.0;const double mid[3]={0.5,0.5,0.5};
        analytic=affineForms(d.type[s],hexNodes(d,s),lin,scaleInv);hexPoint(d,s,mid,typedCenter);
        if(!analytic&&d.type[s]!=pyramidHost)atomicExch(invalid,6);
        const int ncoeff=cellS?cellPolynomials:facePolynomials;
        if(analytic)for(int l=threadIdx.x;l<nS;l+=blockDim.x)if(!chargePolynomials(lin,scaleInv,typedCenter,cellS,hexExponent(d,s,l),1,scratch+std::size_t(l)*ncoeff))atomicExch(invalid,1);
        local=scratch;__syncthreads();}
    else if(graded&&d.affine[s]){double lin[3][4],scaleInv=0.0,ac[3];hexAffineForms(hexNodes(d,s),cellS,lin,scaleInv);hexCenterForms(hexNodes(d,s),cellS,lin,ac);
        const int ncoeff=cellS?84:35;
        for(int l=threadIdx.x;l<nS;l+=blockDim.x)if(!hexMonomialCoeff(lin,scaleInv,cellS,hexExponent(d,s,l),scratch+std::size_t(l)*ncoeff))atomicExch(invalid,1);
        local=scratch;__syncthreads();}
    const bool switched=graded&&d.innerSwitch>0.0;const int nFar=hexFarRule(d,farLevel).n;
    const int farPoints=cls==hexFarClass||switched?(cellS?nFar*nFar*nFar:nFar*nFar):0;
    const std::size_t gradedCoeff=graded?std::size_t(hexGradedCoeff(d.order)):std::size_t(0);
    double* farTable=scratch+gradedCoeff;if(farPoints)hexFarDensities(d,s,nS,farTable,farLevel);
    double* rows=scratch+std::size_t(npS)*(4+nS)+gradedCoeff+std::size_t(farPoints)*nS;const int width=nS+nT;
    for(int q=threadIdx.x;q<npT;q+=blockDim.x){double xi[3],wq;
        if(product){const int a=q/(N*nc),b=(q/nc)%N,c=q%nc;xi[0]=R.x[a];xi[1]=R.x[b];xi[2]=cellT?R.x[c]:0.0;wq=R.w[a]*R.w[b]*(cellT?R.w[c]:1.0);}
        else{const int ix=q%N,iy=(q/N)%N,iz=q/(N*N);xi[0]=R.x[ix];xi[1]=R.x[iy];xi[2]=cellT?R.x[iz]:0.0;wq=R.w[ix]*R.w[iy]*(cellT?R.w[iz]:1.0);}
        double X[3],Xm[3];hexPoint(d,t,xi,X);hexImagePoint(d,img,X,Xm);double* row=rows+std::size_t(q)*width;for(int l=0;l<nS;++l)row[l]=0.0;
        if(product){for(int p=0;p<npS;++p){const double* sp=source+std::size_t(p)*(4+nS);
                const double dx=Xm[0]-sp[0],dy=Xm[1]-sp[1],dz=Xm[2]-sp[2];const double r=sqrt(dx*dx+dy*dy+dz*dz);if(r<1e-300)continue;
                const double w=sp[3]/r;for(int l=0;l<nS;++l)row[l]+=w*sp[4+l];}}
        else if(cls==hexFarClass||(switched&&hexBeyondSource(d,s,Xm)))hexFarInnerTable(d,s,Xm,nS,farTable,row,farLevel);
        else if(cls==hexAffineClass)hexAffineInner(d,s,Xm,nS,row);
        else if(typedSource){const bool onSource=self&&fabs(X[0]-Xm[0])+fabs(X[1]-Xm[1])+fabs(X[2]-Xm[2])<=1e-12*d.size[s];
            if(!hexTypedInner(d,onSource,s,Xm,xi,analytic,typedCenter,local,nS,row))atomicExch(invalid,7);}
        else if(graded){if(!hexGradedInner(d,self,s,X,Xm,xi,nS,row,local))atomicExch(invalid,d.tangent.n?7:5);}
        else if(!hexCompositeInner(d,self,s,X,Xm,xi,nS,row))atomicExch(invalid,4);
        for(int l=0;l<nT;++l)row[nS+l]=wq*hexDensity(d,t,hexExponent(d,t,l),xi);
    }
    __syncthreads();
    for(int e=threadIdx.x;e<nT*nS;e+=blockDim.x){const int lt=e/nS,ls=e%nS;double sum=0.0,comp=0.0;
        for(int q=0;q<npT;++q){const double* row=rows+std::size_t(q)*width;const double v=row[nS+lt]*row[ls];
            if(product)hexNeumaier(sum,comp,v);else sum+=v;}
        out[e]=(sum+comp)*hexInvFourPi;}
}
// CPU QuadBlockHexPairDuffy for one BDM1 touching pair (at most 8 x 8
// charges): per-thread partial sums of each subdomain, a fixed-order block
// tree reduction, Neumaier over the subdomains.
__device__ inline void hexDuffyBlock(const HexDevice& d,int t,int s,int img,const HexAdjacency& adj,int nT,int nS,double* red,double* acc,double* cmp,double* out){
    for(int e=threadIdx.x;e<64;e+=blockDim.x){acc[e]=0.0;cmp[e]=0.0;}
    const int dT=hexCell(d,t)?3:2,dS=hexCell(d,s)?3:2,e=adj.dim,ntrT=dT-e,ntrS=dS-e,nrel=e,kcone=nrel+ntrT+ntrS,ndim=dT+dS;
    const HexRule& R=d.affine[t]&&d.affine[s]?d.pairAffine:d.pair;const int N=R.n;
    long long total=N;for(int k=1;k<ndim;++k)total*=N;
    const double* ndT=hexNodes(d,t);const double* ndS=hexNodes(d,s);
    const int ndomain=2*nrel+ntrT+ntrS,nsignsets=1<<(nrel>1?nrel-1:0);
    double part[64];__syncthreads();
    for(int sub=0;sub<ndomain*nsignsets;++sub){
        const int dsel=sub/nsignsets,signset=sub%nsignsets;int dom;double sgn;
        if(dsel<2*nrel){dom=dsel/2;sgn=dsel%2==0?1.0:-1.0;}else{dom=nrel+(dsel-2*nrel);sgn=1.0;}
        double relsign[3]={1.0,1.0,1.0};
        {int bit=0;for(int c=0;c<nrel;++c){if(c==dom)continue;if(bit<(nrel>1?nrel-1:0))relsign[c]=((signset>>bit)&1)?-1.0:1.0;++bit;}}
        const int extra=(dom>=nrel&&nrel>0)?2:1;
        for(int ex=0;ex<extra;++ex){
            if(extra==2)relsign[nrel-1]=ex==0?1.0:-1.0;
            for(int i=0;i<nT*nS;++i)part[i]=0.0;
            for(long long flat=threadIdx.x;flat<total;flat+=blockDim.x){
                int idx[6];long long rest=flat;for(int k=0;k<ndim;++k){idx[k]=int(rest%N);rest/=N;}
                const double wv=R.x[idx[0]];double weight=R.w[idx[0]],cone[6];int yi=1;
                for(int c=0;c<kcone;++c){if(c==dom){cone[c]=sgn*wv;continue;}
                    const double g=R.x[idx[yi]];weight*=R.w[idx[yi]];++yi;cone[c]=c<nrel?relsign[c]*wv*g:wv*g;}
                for(int c=1;c<kcone;++c)weight*=wv;
                double zT[3]={0.0,0.0,0.0},zS[3]={0.0,0.0,0.0};
                for(int i=0;i<nrel;++i){const double u=cone[i],lo=fmax(0.0,-u),hi=fmin(1.0,1.0-u),g=R.x[idx[yi]];
                    weight*=R.w[idx[yi]]*(hi-lo);++yi;const double f=lo+g*(hi-lo);zT[i]=f;zS[i]=f+u;}
                for(int j=0;j<ntrT;++j)zT[e+j]=1.0-cone[nrel+j];
                for(int j=0;j<ntrS;++j)zS[e+j]=cone[nrel+ntrT+j];
                double xi[3]={0.0,0.0,0.0},eta[3]={0.0,0.0,0.0};
                for(int k=0;k<dT;++k)xi[adj.permT[k]]=adj.flipT[k]?1.0-zT[k]:zT[k];
                for(int k=0;k<dS;++k)eta[adj.permS[k]]=adj.flipS[k]?1.0-zS[k]:zS[k];
                double XT[3],XTm[3],XS[3];
                if(dT==3)hexMapX(ndT,xi,XT);else quadMapX(ndT,xi,XT);
                hexImagePoint(d,img,XT,XTm);
                if(dS==3)hexMapX(ndS,eta,XS);else quadMapX(ndS,eta,XS);
                const double dx=XTm[0]-XS[0],dy=XTm[1]-XS[1],dz=XTm[2]-XS[2];const double r=sqrt(dx*dx+dy*dy+dz*dz);
                if(r>1e-300){weight/=r;double qS[8];for(int l=0;l<nS;++l)qS[l]=hexMono(hexExponent(d,s,l),eta);
                    for(int lt=0;lt<nT;++lt){const double wl=weight*hexMono(hexExponent(d,t,lt),xi);for(int ls=0;ls<nS;++ls)part[lt*nS+ls]+=wl*qS[ls];}}
            }
            for(int i=0;i<nT*nS;++i){red[threadIdx.x]=part[i];__syncthreads();
                for(int stride=blockDim.x/2;stride>0;stride>>=1){if(threadIdx.x<stride)red[threadIdx.x]+=red[threadIdx.x+stride];__syncthreads();}
                if(threadIdx.x==0)hexNeumaier(acc[i],cmp[i],red[0]);__syncthreads();}
        }
    }
    for(int i=threadIdx.x;i<nT*nS;i+=blockDim.x)out[i]=(acc[i]+cmp[i])*hexInvFourPi;
}
// One directed host-pair block per work item; t < 0 marks the absent second
// direction of a self pair. Codes: 3 sub-simplex pair rule, 4 failed composite anchor, 5 failed closest point.
struct HexWork {int t,s,img;std::size_t raw;};
__global__ void __launch_bounds__(hexThreads) hexBlocks(HexDevice d,const HexWork* work,std::size_t first,std::size_t count,
        double* scratch,std::size_t stride,double* raw,int* invalid){
    __shared__ double red[hexThreads],acc[64],cmp[64];__shared__ int cls;__shared__ HexAdjacency adj;
    for(std::size_t item=first+blockIdx.x;item<first+count;item+=gridDim.x){
        const HexWork w=work[item];if(w.t<0)continue;
        if(threadIdx.x==0){HexAdjacency a;cls=hexClassify(d,w.t,w.s,a,w.img);adj=a;}
        __syncthreads();
        const int nT=d.count[w.t],nS=d.count[w.s];
        if(cls==hexSubPairClass){if(threadIdx.x==0)atomicExch(invalid,3);}
        else if(cls==hexDuffyClass){if(nT>8||nS>8){if(threadIdx.x==0)atomicExch(invalid,2);}
            else hexDuffyBlock(d,w.t,w.s,w.img,adj,nT,nS,red,acc,cmp,raw+w.raw);}
        else hexTensorBlock(d,w.t,w.s,w.img,cls,nT,nS,scratch+std::size_t(blockIdx.x)*stride,raw+w.raw,invalid);
        __syncthreads();
    }
}
__global__ void hexPrepareHosts(HexDevice d,double* center,double* size,unsigned char* affine,double* affineCenter,double* cellCoeff,double* faceCoeff,int* invalid){
    const int g=int(blockIdx.x)*blockDim.x+threadIdx.x;if(g>=d.hosts)return;
    const bool cell=hexCell(d,g);const double* nd=hexNodes(d,g);
    const int corners[8]={0,2,6,8,18,20,24,26};const int nc=cell?8:4;
    double c[3]={0.0,0.0,0.0};for(int i=0;i<nc;++i)for(int k=0;k<3;++k)c[k]+=nd[3*corners[i]+k]/nc;
    double r=0.0;for(int i=0;i<nc;++i){const double* p=nd+3*corners[i];const double dx=p[0]-c[0],dy=p[1]-c[1],dz=p[2]-c[2];r=fmax(r,sqrt(dx*dx+dy*dy+dz*dz));}
    for(int k=0;k<3;++k)center[3*g+k]=c[k];size[g]=r;
    for(int k=0;k<3*(cell?27:9);++k)if(!isfinite(nd[k])){atomicExch(invalid,1);return;}
    if(!(r>0.0)){atomicExch(invalid,1);return;}
    double lin[3][4],scaleInv=0.0;const bool ok=hexAffineForms(nd,cell,lin,scaleInv);
    affine[g]=ok?1:0;
    if(!ok)return;
    double ac[3];hexCenterForms(nd,cell,lin,ac);
    for(int k=0;k<3;++k)affineCenter[3*g+k]=ac[k];
    if(d.order!=2||d.type)return;
    for(int mono=0;mono<(cell?27:9);++mono){const int e[3]={mono%3,(mono/3)%3,cell?mono/9:0};
        double* dst=cell?cellCoeff+(std::size_t(g)*27+mono)*84:faceCoeff+(std::size_t(g-d.cells)*9+mono)*35;
        if(!hexMonomialCoeff(lin,scaleInv,cell,e,dst))atomicExch(invalid,1);}
}
// Outer tensor rule points of every host (CPU GetHexFarRule): x fastest.
__global__ void hexFarRules(HexDevice d,HexRule R,double* X,double* W){
    const int g=int(blockIdx.x)*blockDim.x+threadIdx.x;if(g>=d.hosts)return;
    const int n1=R.n,np=hexCell(d,g)?n1*n1*n1:n1*n1;const std::size_t stride=std::size_t(n1)*n1*n1;
    for(int q=0;q<np;++q){const int ix=q%n1,iy=(q/n1)%n1,iz=q/(n1*n1);const double xi[3]={R.x[ix],R.x[iy],hexCell(d,g)?R.x[iz]:0.0};
        hexPoint(d,g,xi,X+3*(std::size_t(g)*stride+q));W[std::size_t(g)*stride+q]=R.w[ix]*R.w[iy]*(hexCell(d,g)?R.w[iz]:1.0);}
}
// A host pair enters the table of an image when either directed block can be non-far: both directions test their own separation.
__device__ inline bool hexCandidate(const HexDevice& d,int i,int j,int img=0){
    const double bound=d.candidate*(d.size[i]+d.size[j])+1e-12;
    return !(hexSep(d,i,j,img)>bound)||(img!=0&&!(hexSep(d,j,i,img)>bound));}
// Table rows: hosts per section, section 0 the direct term and section img each image.
__device__ inline unsigned long long hexKey(const HexDevice& d,int img,int lo,int hi){
    const auto hosts=(unsigned long long)d.hosts;return ((unsigned long long)img*hosts+(unsigned long long)lo)*hosts+(unsigned long long)hi;}
__device__ inline void hexKeyParts(const HexDevice& d,unsigned long long key,int& img,int& lo,int& hi){
    const auto hosts=(unsigned long long)d.hosts;img=int(key/(hosts*hosts));key%=hosts*hosts;lo=int(key/hosts);hi=int(key%hosts);}
__global__ void hexCountPairs(HexDevice d,std::size_t* counts){
    const auto row=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(row>=std::size_t(d.hosts)*(1+d.imageCount))return;
    const int img=int(row/d.hosts),i=int(row%d.hosts);
    std::size_t c=0;for(int j=i;j<d.hosts;++j)if(hexCandidate(d,i,j,img))++c;counts[row]=c;
}
// Keys (img, lo, hi) ascend because sections, hosts and their partners are visited in order.
__global__ void hexFillPairs(HexDevice d,const std::size_t* starts,unsigned long long* keys,std::size_t* sizes){
    const auto row=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(row>=std::size_t(d.hosts)*(1+d.imageCount))return;
    const int img=int(row/d.hosts),i=int(row%d.hosts);
    std::size_t k=starts[row];
    for(int j=i;j<d.hosts;++j)if(hexCandidate(d,i,j,img)){keys[k]=hexKey(d,img,i,j);
        sizes[k]=std::size_t(d.count[i])*std::size_t(d.count[j]);++k;}
}
__global__ void hexPairWork(HexDevice d,const unsigned long long* keys,const std::size_t* offsets,std::size_t pairs,HexWork* work){
    const auto p=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(p>=pairs)return;
    int img,lo,hi;hexKeyParts(d,keys[p],img,lo,hi);
    const std::size_t size=std::size_t(d.count[lo])*std::size_t(d.count[hi]);
    const bool oneSided=d.type&&lo!=hi;const bool lowerSource=oneSided&&radia_hdiv_typed::lowerHostIsSource(d.type[lo],d.type[hi]);
    work[2*p]={lowerSource?-1:lo,hi,img,2*offsets[p]};work[2*p+1]={lo==hi||(oneSided&&!lowerSource)?-1:hi,lo,img,2*offsets[p]+size};
}
// Cost class of one work item, cheapest first. The order only groups items: each block writes its own values.
__global__ void hexWorkCost(HexDevice d,const HexWork* work,std::size_t items,unsigned char* cost){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=items)return;
    const HexWork w=work[i];if(w.t<0){cost[i]=0;return;}
    HexAdjacency a;const int cls=hexClassify(d,w.t,w.s,a,w.img);
    if(cls==hexFarClass){cost[i]=1;return;}
    if(cls!=hexGradedClass){cost[i]=cls==hexAffineClass?2:3;return;}
    if(!d.type){cost[i]=4;return;}
    using namespace radia_hdiv_typed;
    if(d.type[w.s]==triHost||d.type[w.s]==tetHost){cost[i]=2;return;}
    double lin[3][4],scaleInv=0.0;
    cost[i]=d.type[w.s]==pyramidHost&&!affineForms(d.type[w.s],hexNodes(d,w.s),lin,scaleInv)?5:4;
}
__global__ void hexCostCount(const unsigned char* cost,std::size_t items,unsigned long long* counts){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<items)atomicAdd(counts+cost[i],1ULL);
}
// Host-side launch size. A launch runs one block per processor; a block takes more items while
// kernels stay short, so a slow item is shared out instead of holding every other processor.
// Scheduling only: each block writes its own values.
class HexLaunchPace {
    std::size_t blocks,items;std::chrono::steady_clock::time_point started;
public:
    explicit HexLaunchPace(std::size_t processors):blocks(processors),items(processors){}
    std::size_t begin(std::size_t remaining){started=std::chrono::steady_clock::now();return std::min(items,remaining);}
    void end(){const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        if(seconds<.02)items*=2;else if(seconds>.1)items=std::max(blocks,items/2);}
};
// 0.5 (G_AB + G_BA^T) into the lo-host row orientation, or the one integrated direction of a mixed mesh pair (CPU GetHexSymBlock).
// One direction of an image pair is G_T or G_T^-1 by the host types alone; the list closed under inversion with equal signs sums both.
__global__ void hexSymmetrize(HexDevice d,const unsigned long long* keys,const std::size_t* offsets,std::size_t pairs,const double* raw,double* values,int* invalid){
    const auto p=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(p>=pairs)return;
    int img,lo,hi;hexKeyParts(d,keys[p],img,lo,hi);
    const int nA=d.count[lo],nB=d.count[hi];const double* ab=raw+2*offsets[p];const double* ba=lo==hi?ab:ab+std::size_t(nA)*nB;
    const bool oneSided=d.type&&lo!=hi;const bool lowerSource=oneSided&&radia_hdiv_typed::lowerHostIsSource(d.type[lo],d.type[hi]);
    for(int a=0;a<nA;++a)for(int b=0;b<nB;++b){const double v=lowerSource?ba[b*nA+a]:oneSided?ab[a*nB+b]:0.5*(ab[a*nB+b]+ba[b*nA+a]);if(!isfinite(v))atomicExch(invalid,1);values[offsets[p]+std::size_t(a)*nB+b]=v;}
}
// Far tensor product of one charge pair beyond the pair table, target points mapped by T^-1.
__device__ inline double hexFarEntry(const HexDevice& d,std::size_t a,std::size_t b,int img=0){
    const int t=d.chargeHost[a],s=d.chargeHost[b];const int* ea=d.exponents+3*a;const int* eb=d.exponents+3*b;
    const int level=d.type?hexFarLevel(d,t,s,img):0;const HexRule& R=hexFarRule(d,level);
    const int n1=R.n,npT=hexCell(d,t)?n1*n1*n1:n1*n1,npS=hexCell(d,s)?n1*n1*n1:n1*n1;
    const double* XT=hexFarPoints(d,level,t);const double* WT=hexFarWeights(d,level,t);
    const double* XS=hexFarPoints(d,level,s);const double* WS=hexFarWeights(d,level,s);
    double sum=0.0;
    for(int qt=0;qt<npT;++qt){const double xt[3]={R.x[qt%n1],R.x[(qt/n1)%n1],hexCell(d,t)?R.x[qt/(n1*n1)]:0.0};double inner=0.0;
        double xm[3];hexImagePoint(d,img,XT+3*qt,xm);
        for(int qs=0;qs<npS;++qs){const double dx=xm[0]-XS[3*qs],dy=xm[1]-XS[3*qs+1],dz=xm[2]-XS[3*qs+2];
            const double r=sqrt(dx*dx+dy*dy+dz*dz);if(r<1e-300)continue;
            const double xs[3]={R.x[qs%n1],R.x[(qs/n1)%n1],hexCell(d,s)?R.x[qs/(n1*n1)]:0.0};inner+=(WS[qs]/r)*hexDensity(d,s,eb,xs);}
        sum+=(WT[qt]*hexDensity(d,t,ea,xt))*inner;}
    return sum*hexInvFourPi;
}
// The symmetrized entry of one term: 0 direct, img an image.
__device__ inline double hexTermEntry(const HexDevice& d,std::size_t a,std::size_t b,int img){
    const int ga=d.chargeHost[a],gb=d.chargeHost[b],lo=ga<gb?ga:gb,hi=ga<gb?gb:ga;
    const unsigned long long key=hexKey(d,img,lo,hi);
    std::size_t L=0,H=d.pairs;while(L<H){const std::size_t m=(L+H)/2;if(d.keys[m]<key)L=m+1;else H=m;}
    if(L<d.pairs&&d.keys[L]==key){const double* v=d.values+d.offsets[L];const int la=int(a)-d.first[ga],lb=int(b)-d.first[gb];
        return ga<=gb?v[la*d.count[gb]+lb]:v[lb*d.count[ga]+la];}
    // A pair missing from the table must be a far pair by the same test.
    if(hexCandidate(d,ga,gb,img))return NAN;
    // The direct far product is symmetric; an image pair folds both directed blocks like the CPU GetHexSymBlock.
    return img==0?hexFarEntry(d,a,b):0.5*(hexFarEntry(d,a,b,img)+hexFarEntry(d,b,a,img));
}
__device__ inline double hexEntry(const HexDevice& d,std::size_t a,std::size_t b){
    double value=hexTermEntry(d,a,b,0);
    for(int img=1;img<=d.imageCount;++img)value+=d.images[img-1].sign*hexTermEntry(d,a,b,img);
    return value;
}
// Device state of one HEX Gram: host lattices, charges, the pair table and its
// symmetrized blocks. Raw directed blocks, work items and scratch live only
// inside prepare(); the counts that size the table are the only D2H values.
struct HexPrepared {
    Buffer<double> cells,faces,center,size,affineCenter,cellCoeff,faceCoeff,farX,farW,values,levelX[hexFarLevels],levelW[hexFarLevels];
    Buffer<int> first,count,chargeHost,exponents,invalid,types;Buffer<unsigned char> affine;Buffer<HexImage> images;
    Buffer<unsigned long long> keys;Buffer<std::size_t> offsets;
    HexDevice d{};int order=1;
    // Work items and wall time per cost class, for the setup record.
    unsigned long long classItems[6]={};std::size_t classMicroseconds[6]={};
    static HexRule rule(const HostVector<double>& x,const HostVector<double>& w){
        HexRule r{};r.n=int(x.size());for(int k=0;k<r.n;++k){r.x[k]=x[k];r.w[k]=w[k];}return r;}
    // CPU HexGradedTensorCloud: the C2 smootherstep clusters the 1D near rule at both ends of every axis.
    static HexRule graded(HexRule r){for(int k=0;k<r.n;++k){const double s=r.x[k];
            r.x[k]=s*s*s*(10.0-15.0*s+6.0*s*s);r.w[k]*=30.0*s*s*(1.0-s)*(1.0-s);}return r;}
    template<class Record>int flag(cudaStream_t stream,Record& record,const char* payload){int value=0;
        check(cudaMemcpyAsync(&value,invalid.data,sizeof(int),cudaMemcpyDeviceToHost,stream));
        check(cudaStreamSynchronize(stream));record(payload,"D2H",sizeof(int));return value;}
    template<class Record>static std::size_t last(const Buffer<std::size_t>& a,const Buffer<std::size_t>& b,std::size_t m,cudaStream_t stream,
            Record& record,const char* payload){
        std::size_t x=0,y=0;check(cudaMemcpyAsync(&x,a.data+m-1,sizeof(x),cudaMemcpyDeviceToHost,stream));
        check(cudaMemcpyAsync(&y,b.data+m-1,sizeof(y),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));
        record(payload,"D2H",2*sizeof(std::size_t));return x+y;}
    // Host-side layout of the charge arrays: grouped by global host, cells before faces.
    template<class Block>static void layout(const Block& b,HostVector<int>& firstH,HostVector<int>& countH,HostVector<int>& globalH,int& order){
        const std::size_t cellsN=b.cellCount(),hosts=cellsN+b.faceCount(),n=b.charges();
        firstH.assign(hosts,0);countH.assign(hosts,0);globalH.assign(n,0);order=1;
        for(std::size_t a=0;a<n;++a){const int g=b.kind[a]==0?b.host[a]:int(cellsN)+b.host[a];
            if(a&&g<globalH[a-1])throw std::invalid_argument("HEX charges must be grouped by host, cells before faces");
            if(!countH[std::size_t(g)])firstH[std::size_t(g)]=int(a);
            ++countH[std::size_t(g)];globalH[a]=g;
            for(int c=0;c<3;++c)order=std::max(order,b.expo[3*a+c]);}
    }
    template<class Block,class Record>void prepare(const Block& b,bool dense,cudaStream_t stream,Record record,std::size_t& launches,int& invalidCode){
        prepareHosts(b,stream,record,launches,invalidCode);prepareTable(b,dense,stream,record,launches,invalidCode);}
    // Lattices, host data and far rules: everything a field evaluation reads, without the pair table.
    template<class Block,class Record>void prepareHosts(const Block& b,cudaStream_t stream,Record record,std::size_t& launches,int& invalidCode){
        const std::size_t cellsN=b.cellCount(),facesN=b.faceCount(),hosts=cellsN+facesN;
        HostVector<int> firstH{HostFaultAllocator<int>("hex.first")},countH{HostFaultAllocator<int>("hex.count")},globalH{HostFaultAllocator<int>("hex.global")};
        layout(b,firstH,countH,globalH,order);
        auto upload=[&](auto& buffer,const auto& v,const char* payload){if(v.empty())return;buffer.upload(v);record(payload,"H2D",checkedProduct(v.size(),sizeof(v[0])));};
        upload(cells,b.cells,"cell_q2_nodes_m");upload(faces,b.faces,"face_q2_nodes_m");upload(first,firstH,"hex_host_first");upload(count,countH,"hex_host_count");
        upload(chargeHost,globalH,"charge_host");upload(exponents,b.expo,"charge_exponents");upload(invalid,HostVector<int>{0},"hex_invalid_flag");
        d.cellNodes=cells.data;d.faceNodes=faces.data;d.cells=int(cellsN);d.hosts=int(hosts);d.first=first.data;d.count=count.data;
        d.chargeHost=chargeHost.data;d.exponents=exponents.data;d.order=order;d.nearGrade=b.nearGrade;d.farInner=b.farInner;d.innerSwitch=b.innerSwitch;
        d.out=rule(b.ruleOut,b.ruleOutW);d.in=rule(b.ruleIn,b.ruleInW);d.pair=rule(b.rulePair,b.rulePairW);d.pairAffine=rule(b.rulePairAffine,b.rulePairAffineW);
        d.nearOut=graded(rule(b.ruleNear,b.ruleNearW));d.inSelf=rule(b.ruleInSelf,b.ruleInSelfW);d.nearCell=rule(b.ruleNearCell,b.ruleNearCellW);d.tangent=rule(b.ruleTangent,b.ruleTangentW);
        for(int k=0;k<hexFarLevels;++k){d.farLevel[k]=rule(b.ruleFar[k],b.ruleFarW[k]);d.farRatio[k]=b.farRatio[k];}
        const bool typed=!b.hostTypes.empty();if(typed){upload(types,b.hostTypes,"host_types");d.type=types.data;}
        d.images=nullptr;d.imageCount=0;
        if(!b.images.empty()){upload(images,b.images,"image_transforms");d.images=images.data;d.imageCount=int(b.images.size());}
        center.allocate(3*hosts);size.allocate(hosts);affine.allocate(hosts);affineCenter.allocate(3*hosts);
        cellCoeff.allocate(order==2&&!typed?std::max<std::size_t>(1,checkedProduct(cellsN,27*84)):1);
        faceCoeff.allocate(order==2&&!typed?std::max<std::size_t>(1,checkedProduct(facesN,9*35)):1);
        d.center=center.data;d.size=size.data;d.affine=affine.data;d.affineCenter=affineCenter.data;d.cellCoeff=cellCoeff.data;d.faceCoeff=faceCoeff.data;
        const unsigned hostBlocks=unsigned((hosts+127)/128);
        hexPrepareHosts<<<hostBlocks,128,0,stream>>>(d,center.data,size.data,affine.data,affineCenter.data,cellCoeff.data,faceCoeff.data,invalid.data);
        check(cudaGetLastError());++launches;
        if((invalidCode=flag(stream,record,"prepare_flag"))!=0)throw std::invalid_argument("Nonfinite or degenerate HEX host geometry");
        d.farStride=d.out.n*d.out.n*d.out.n;
        farX.allocate(checkedProduct(hosts,3*std::size_t(d.farStride)));farW.allocate(checkedProduct(hosts,std::size_t(d.farStride)));
        d.farX=farX.data;d.farW=farW.data;
        hexFarRules<<<hostBlocks,128,0,stream>>>(d,d.out,farX.data,farW.data);check(cudaGetLastError());++launches;
        for(int k=0;k<hexFarLevels;++k)if(d.farLevel[k].n){const std::size_t points=checkedProduct(checkedProduct(std::size_t(d.farLevel[k].n),d.farLevel[k].n),d.farLevel[k].n);
            levelX[k].allocate(checkedProduct(hosts,3*points));levelW[k].allocate(checkedProduct(hosts,points));d.farLevelX[k]=levelX[k].data;d.farLevelW[k]=levelW[k].data;
            hexFarRules<<<hostBlocks,128,0,stream>>>(d,d.farLevel[k],levelX[k].data,levelW[k].data);check(cudaGetLastError());++launches;}
    }
    template<class Block,class Record>void prepareTable(const Block& b,bool dense,cudaStream_t stream,Record record,std::size_t& launches,int& invalidCode){
        const std::size_t hosts=b.cellCount()+b.faceCount();const bool typed=!b.hostTypes.empty();
        // Dense keeps every host pair; ACA every pair inside the HexHostsTouch radius, which holds every non-far family.
        d.candidate=dense?1e300:1.05;
        // One table section per term: the direct pairs, then those of each image.
        const std::size_t tableRows=checkedProduct(hosts,1+b.images.size());const unsigned rowBlocks=unsigned((tableRows+127)/128);
        if(tableRows>std::size_t(INT_MAX))throw std::invalid_argument("HEX image pair table exceeds the scan range");
        Buffer<std::size_t> hostPairs,hostStarts;hostPairs.allocate(tableRows);hostStarts.allocate(tableRows);
        hexCountPairs<<<rowBlocks,128,0,stream>>>(d,hostPairs.data);check(cudaGetLastError());++launches;
        {std::size_t bytes=0;check(cub::DeviceScan::ExclusiveSum(nullptr,bytes,hostPairs.data,hostStarts.data,int(tableRows),stream));
            Buffer<unsigned char> temporary;temporary.allocate(std::max<std::size_t>(1,bytes));
            check(cub::DeviceScan::ExclusiveSum(temporary.data,bytes,hostPairs.data,hostStarts.data,int(tableRows),stream));++launches;
            check(cudaStreamSynchronize(stream));}
        const std::size_t pairs=last(hostStarts,hostPairs,tableRows,stream,record,"hex_pair_count");
        if(!pairs)throw std::runtime_error("HEX pair table is empty");
        keys.allocate(pairs);Buffer<std::size_t> blockSizes;blockSizes.allocate(pairs);offsets.allocate(pairs);
        hexFillPairs<<<rowBlocks,128,0,stream>>>(d,hostStarts.data,keys.data,blockSizes.data);check(cudaGetLastError());++launches;
        {std::size_t bytes=0;check(cub::DeviceScan::ExclusiveSum(nullptr,bytes,blockSizes.data,offsets.data,int(pairs),stream));
            Buffer<unsigned char> temporary;temporary.allocate(std::max<std::size_t>(1,bytes));
            check(cub::DeviceScan::ExclusiveSum(temporary.data,bytes,blockSizes.data,offsets.data,int(pairs),stream));++launches;
            check(cudaStreamSynchronize(stream));}
        const std::size_t total=last(offsets,blockSizes,pairs,stream,record,"hex_block_values");
        // Phase-1 rows of the widest family: outer tensor (or pair product) points x (source + target) charges, plus the far density table.
        const std::size_t maxCount=order==2?27:8,nOut=std::size_t(d.out.n)*d.out.n*d.out.n;
        const std::size_t nPair=std::size_t(std::max(d.pair.n,d.pairAffine.n)),nPair3=nPair*nPair*nPair;
        const std::size_t nNear=std::size_t(d.nearOut.n),nNearCell=std::size_t(d.nearCell.n?d.nearCell.n:d.nearOut.n);
        const std::size_t nearPoints=d.nearCell.n?std::max(nNear*nNear,nNearCell*nNearCell*nNearCell):nNear*nNear*nNear;
        const std::size_t stride=std::max({nOut*3*maxCount,order==1&&!typed?nPair3*(4+maxCount)+nPair3*2*maxCount:std::size_t(0),
                                           order==1||typed?std::size_t(hexGradedCoeff(order))+(d.innerSwitch>0.0?nOut*maxCount:std::size_t(0))+nearPoints*2*maxCount:std::size_t(0)});
        // One block per processor: more would only share processors and lengthen each kernel.
        int gpu=0;check(cudaGetDevice(&gpu));cudaDeviceProp properties;check(cudaGetDeviceProperties(&properties,gpu));
        const std::size_t launchItems=std::size_t(std::max(properties.multiProcessorCount,1));
        const std::size_t items=checkedProduct(2,pairs);
        if(items>std::size_t(INT_MAX))throw std::runtime_error("HEX pair table exceeds the work sort");
        constexpr std::size_t costClasses=6;
        const std::size_t bytes=checkedProduct(3*total,sizeof(double))+checkedProduct(items,2*sizeof(HexWork)+2)+checkedProduct(launchItems*stride,sizeof(double))+costClasses*sizeof(unsigned long long);
        {std::size_t free,totalMemory;check(cudaMemGetInfo(&free,&totalMemory));
            if(!fitsAboveReserve(free,totalMemory,bytes))throw std::runtime_error("HEX Gram pair table would exceed the GPU reserve");}
        values.allocate(total);Buffer<double> raw,scratch;Buffer<HexWork> unordered,work;raw.allocate(2*total);unordered.allocate(items);work.allocate(items);scratch.allocate(launchItems*stride);
        const unsigned pairBlocks=unsigned((pairs+127)/128);
        hexPairWork<<<pairBlocks,128,0,stream>>>(d,keys.data,offsets.data,pairs,unordered.data);check(cudaGetLastError());++launches;
        // A launch waits for its slowest item, so costly near pairs run together, not one per launch among cheap ones.
        {Buffer<unsigned char> cost,sorted,temporary;cost.allocate(items);sorted.allocate(items);
            hexWorkCost<<<unsigned((items+127)/128),128,0,stream>>>(d,unordered.data,items,cost.data);check(cudaGetLastError());++launches;
            std::size_t sortBytes=0;check(cub::DeviceRadixSort::SortPairs(nullptr,sortBytes,cost.data,sorted.data,unordered.data,work.data,int(items),0,8,stream));
            temporary.allocate(std::max<std::size_t>(1,sortBytes));
            check(cub::DeviceRadixSort::SortPairs(temporary.data,sortBytes,cost.data,sorted.data,unordered.data,work.data,int(items),0,8,stream));++launches;
            Buffer<unsigned long long> classes;classes.allocate(costClasses);check(cudaMemsetAsync(classes.data,0,sizeof(classItems),stream));
            hexCostCount<<<unsigned((items+127)/128),128,0,stream>>>(sorted.data,items,classes.data);check(cudaGetLastError());++launches;
            check(cudaMemcpyAsync(classItems,classes.data,sizeof(classItems),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));
            record("hex_cost_classes","D2H",sizeof(classItems));}
        // Bounded launches keep every kernel short on display-attached GPUs. Class 0 holds the absent directions; the pace restarts where the cost jumps.
        std::size_t item=classItems[0];
        for(std::size_t c=1;c<costClasses;++c){HexLaunchPace pace(launchItems);const std::size_t end=item+classItems[c];const auto classStarted=std::chrono::steady_clock::now();
            while(item<end){const auto here=pace.begin(end-item);
                hexBlocks<<<unsigned(std::min(launchItems,here)),hexThreads,0,stream>>>(d,work.data,item,here,scratch.data,stride,raw.data,invalid.data);
                check(cudaGetLastError());check(cudaStreamSynchronize(stream));++launches;pace.end();item+=here;}
            classMicroseconds[c]=std::size_t(std::chrono::duration<double,std::micro>(std::chrono::steady_clock::now()-classStarted).count());}
        if((invalidCode=flag(stream,record,"hex_block_flag"))!=0)
            throw std::runtime_error(invalidCode==3?"HEX host pair needs the sub-simplex pair rule between a near band below the distorted far switch, which has no device route"
                                    :invalidCode==1?"HEX graded near rule produced a nonfinite source polynomial"
                                    :invalidCode==4?"HEX composite near rule could not invert the source host map"
                                    :invalidCode==5?"HEX graded near rule could not find the closest source point"
                                    :invalidCode==6?"Mixed charge Gram has a degenerate tet or triangle source"
                                    :invalidCode==7?"Mixed charge Gram has a degenerate tangent of a warped source"
                                    :"HEX pair-domain Duffy block exceeds 8 x 8 charges");
        hexSymmetrize<<<pairBlocks,128,0,stream>>>(d,keys.data,offsets.data,pairs,raw.data,values.data,invalid.data);check(cudaGetLastError());++launches;
        if((invalidCode=flag(stream,record,"hex_symmetrize_flag"))!=0)throw std::runtime_error("HEX device Gram produced a nonfinite block entry");
        d.keys=keys.data;d.offsets=offsets.data;d.values=values.data;d.pairs=pairs;
    }
};
} // namespace radia_hmatrix_cuda
