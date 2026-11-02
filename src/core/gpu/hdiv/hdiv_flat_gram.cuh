// TET BDM1/BDM2 charge Gram generated on the device from host geometry,
// charge monomials and reference rules. The entry dispatch mirrors the CPU
// HighOrderTetEntryStrategy without images: well-separated pairs take the low
// double rule. Affine: self and touching pairs of affine P1 hosts the shared
// graded outer rule (CPU P1TouchHosts/P1TouchEntryValue); other pairs the
// exact analytic moment potentials on the outer rule. Curved P2: touching
// hosts (shared edge, or a shared vertex with a face) the symmetrized curved
// Duffy potentials, other near pairs the one-sided high-rule tensor product
// (CPU direct curved block); P1TouchHosts is off for curved, as on the CPU. Only
// primitives cross to the device; the Gram is written into the resident
// deferred dense leaf (or the ACA leaves of hdiv_flat_aca.cuh) and never
// leaves the device.
#pragma once
#include <climits>
#include <cmath>
#include <cstdlib>
#include "hmatrix_resident.cuh"
#include "hdiv_p1_diagonal_input.hpp"
#define RADIA_P1_KERNEL __device__
#define RADIA_P1_FIELD rad_hdiv_cuda
#include "hdiv_p1_diagonal_math.inc"
#undef RADIA_P1_FIELD
#undef RADIA_P1_KERNEL
#include "hdiv_hex_gram.cuh"
namespace radia_hmatrix_cuda {
inline std::size_t flatSum(std::size_t a,std::size_t b){
    if(a>std::numeric_limits<std::size_t>::max()-b)throw std::overflow_error("Flat device byte count overflow");return a+b;}
struct FlatRulePoint {double xi[3],weight;};
static_assert(sizeof(FlatRulePoint)==32,"hdiv_gpu._RULE_ROW_BYTES ledger contract");
struct FlatCharge {
    double vertices[12],inverse[9],center[3],polyCenter[3],size,a,b[3],c[9];
    int kind,host,expo[3],degree,numbers[4];
};
struct FlatTransfer {const char* payload;const char* direction;std::size_t bytes;};
struct FlatBlock {
    HostVector<double> cells{HostFaultAllocator<double>("flat.cells")},faces{HostFaultAllocator<double>("flat.faces")};
    HostVector<int> host{HostFaultAllocator<int>("flat.host")},kind{HostFaultAllocator<int>("flat.kind")},expo{HostFaultAllocator<int>("flat.expo")};
    HostVector<FlatRulePoint> tet{HostFaultAllocator<FlatRulePoint>("flat.tet")},tri{HostFaultAllocator<FlatRulePoint>("flat.tri")},
        tetLow{HostFaultAllocator<FlatRulePoint>("flat.tet-low")},triLow{HostFaultAllocator<FlatRulePoint>("flat.tri-low")};
    double farFactor=1e30;
    // Curved P2 hosts: cells/faces hold 10/6 P2 nodes, mesh vertex numbers
    // decide touching and (gl,gw) is the [0,1] Duffy Gauss rule.
    bool curved=false;
    HostVector<double> gl{HostFaultAllocator<double>("flat.gl")},gw{HostFaultAllocator<double>("flat.gw")};
    HostVector<int> cellNumbers{HostFaultAllocator<int>("flat.cell-numbers")},faceNumbers{HostFaultAllocator<int>("flat.face-numbers")};
    // HEX hosts: cells/faces hold the 27/9-node Q2 lattices, charges are
    // reference monomials grouped by host, and the four [0,1] Gauss rules
    // (outer tensor, radial inner, pair product, affine pair product) with the
    // near band and the radial reach of the CPU hex Gram.
    bool hex=false;double nearGrade=1.0,farInner=1.0;
    HostVector<double> ruleOut{HostFaultAllocator<double>("hex.out")},ruleOutW{HostFaultAllocator<double>("hex.out-w")},
        ruleIn{HostFaultAllocator<double>("hex.in")},ruleInW{HostFaultAllocator<double>("hex.in-w")},
        rulePair{HostFaultAllocator<double>("hex.pair")},rulePairW{HostFaultAllocator<double>("hex.pair-w")},
        rulePairAffine{HostFaultAllocator<double>("hex.pair-affine")},rulePairAffineW{HostFaultAllocator<double>("hex.pair-affine-w")};
    std::size_t charges()const{return host.size();}
    std::size_t cellStride()const{return hex?81:curved?30:12;}
    std::size_t faceStride()const{return hex?27:curved?18:9;}
    std::size_t cellCount()const{return cells.size()/cellStride();}
    std::size_t faceCount()const{return faces.size()/faceStride();}
    void validate()const{
        const auto n=host.size();
        // Charge, cell and face indices are int32 on both sides of the binding.
        if(!n||n>std::size_t(INT_MAX)||kind.size()!=n||expo.size()!=3*n||cells.size()%cellStride()||faces.size()%faceStride()||
           cellCount()>std::size_t(INT_MAX)||faceCount()>std::size_t(INT_MAX))
            throw std::invalid_argument("Flat charge payload shape mismatch");
        if(hex){validateHex();return;}
        if(curved){
            if(gl.empty()||gl.size()>64||gw.size()!=gl.size()||cellNumbers.size()!=4*cellCount()||faceNumbers.size()!=3*faceCount())
                throw std::invalid_argument("Curved charge geometry requires 4/3 vertex numbers per host and one Duffy rule of 1..64 points");
            for(std::size_t k=0;k<gl.size();++k)if(!std::isfinite(gl[k])||!std::isfinite(gw[k])||!(gl[k]>0)||!(gl[k]<1)||!(gw[k]>0))
                throw std::invalid_argument("Curved Duffy rule requires interior [0,1] points and positive weights");
            for(const auto* numbers:{&cellNumbers,&faceNumbers})for(int v:*numbers)if(v<0)throw std::invalid_argument("Curved host vertex numbers must be nonnegative");
        }else if(!gl.empty()||!gw.empty()||!cellNumbers.empty()||!faceNumbers.empty())
            throw std::invalid_argument("Affine charge geometry takes no curved Duffy rule or vertex numbers");
        if(tet.empty()||tri.empty()||tet.size()>4096||tri.size()>4096||tetLow.size()>4096||triLow.size()>4096)
            throw std::invalid_argument("Flat reference rule size outside budget");
        if(!(farFactor>0)||std::isnan(farFactor))throw std::invalid_argument("Invalid far-pair factor");
        if(farFactor<1e29&&(tetLow.empty()||triLow.empty()))throw std::invalid_argument("Far split requires both low rules");
        for(const auto* rule:{&tet,&tri,&tetLow,&triLow})for(const auto& q:*rule)
            if(!std::isfinite(q.xi[0])||!std::isfinite(q.xi[1])||!std::isfinite(q.xi[2])||!std::isfinite(q.weight))
                throw std::invalid_argument("Nonfinite reference rule");
        for(const auto* coordinates:{&cells,&faces})for(double v:*coordinates)
            if(!std::isfinite(v))throw std::invalid_argument("Nonfinite flat charge geometry");
        // Per-component exponent bounds come before any sum, so no exponent can overflow a degree.
        const std::size_t cellTotal=cellCount(),faceTotal=faceCount();
        for(std::size_t a=0;a<n;++a){
            const int k=kind[a],h=host[a],limit=k==0?1:2;
            if((k!=0&&k!=1)||h<0||std::size_t(h)>=(k==0?cellTotal:faceTotal))
                throw std::invalid_argument("Flat charge kind or host index out of range");
            for(int c=0;c<3;++c)if(expo[3*a+c]<0||expo[3*a+c]>limit)
                throw std::invalid_argument("Flat charge exponent outside BDM1/BDM2 support");
            if((k==1&&expo[3*a+2]!=0)||expo[3*a]+expo[3*a+1]+expo[3*a+2]>limit)
                throw std::invalid_argument("Flat device Gram supports cell degree <= 1 and face degree <= 2 only");
        }
        // The device has no one-sided far rule; that diagnostic route stays CPU-only.
        const char* oneSided=std::getenv("RADIA_HDIV_HO_FAR_ONESIDED");
        if(oneSided&&oneSided[0]!='\0'&&oneSided[0]!='0')
            throw std::invalid_argument("RADIA_HDIV_HO_FAR_ONESIDED has no device Gram route");
    }
    void validateHex()const{
        if(curved||!tet.empty()||!tri.empty()||!tetLow.empty()||!triLow.empty()||farFactor<1e29||!gl.empty()||!gw.empty()||
           !cellNumbers.empty()||!faceNumbers.empty())
            throw std::invalid_argument("HEX charge geometry takes no TET rules, far split, curved Duffy rule or vertex numbers");
        for(const auto* rule:{&ruleOut,&ruleIn,&rulePair,&rulePairAffine}){const auto& w=rule==&ruleOut?ruleOutW:rule==&ruleIn?ruleInW:rule==&rulePair?rulePairW:rulePairAffineW;
            if(rule->empty()||rule->size()>std::size_t(hexRuleCap)||w.size()!=rule->size())throw std::invalid_argument("HEX rules require 1..16 Gauss points each");
            for(std::size_t k=0;k<rule->size();++k)if(!std::isfinite((*rule)[k])||!std::isfinite(w[k])||!((*rule)[k]>0)||!((*rule)[k]<1)||!(w[k]>0))
                throw std::invalid_argument("HEX rules require interior [0,1] points and positive weights");}
        if(rulePair.size()<2||rulePairAffine.size()<2)throw std::invalid_argument("HEX pair rules require at least 2 points");
        if(!std::isfinite(nearGrade)||!(nearGrade>0)||!std::isfinite(farInner)||!(farInner>0))throw std::invalid_argument("HEX near band and radial reach must be positive");
        for(const auto* coordinates:{&cells,&faces})for(double v:*coordinates)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite HEX host lattice");
        const std::size_t cellTotal=cellCount(),faceTotal=faceCount();int order=1;
        for(std::size_t a=0;a<charges();++a){const int k=kind[a],h=host[a];
            if((k!=0&&k!=1)||h<0||std::size_t(h)>=(k==0?cellTotal:faceTotal))throw std::invalid_argument("HEX charge kind or host index out of range");
            for(int c=0;c<3;++c){if(expo[3*a+c]<0||expo[3*a+c]>2)throw std::invalid_argument("HEX charge exponent outside BDM1/BDM2 support");order=std::max(order,expo[3*a+c]);}
            if(k==1&&expo[3*a+2]!=0)throw std::invalid_argument("HEX face charges carry two reference exponents");}
        // Per-host counts are bounded by the tensor monomials of the order: BDM1 Duffy blocks hold at most 8 x 8.
        HostVector<int> first{HostFaultAllocator<int>("hex.validate-first")},count{HostFaultAllocator<int>("hex.validate-count")},global{HostFaultAllocator<int>("hex.validate-global")};
        HexPrepared::layout(*this,first,count,global,order);
        const int monomials=(order+1)*(order+1);
        for(std::size_t g=0;g<count.size();++g)if(count[g]>(g<cellTotal?monomials*(order+1):monomials))
            throw std::invalid_argument("HEX host carries more charges than its tensor monomials");
        // CPU diagnostic knobs that change the HEX dispatch have no device route.
        for(const char* knob:{"RADIA_HDIV_HEX_PAIR_DUFFY","RADIA_HDIV_HEX_DISTORTED_FAR_FACTOR","RADIA_HDIV_HEX_MAX_ANISOTROPY",
                              "RADIA_HDIV_HEX_FAR_ONESIDED","RADIA_HDIV_HEX_GLPAIR_W_N"}){const char* value=std::getenv(knob);
            if(value&&value[0]!='\0')throw std::invalid_argument(std::string(knob)+" has no device HEX Gram route");}
    }
};
__device__ inline double flatPow(double base,int exponent){double r=1.0;for(int i=0;i<exponent;++i)r*=base;return r;}
// Physical-point monomial through the host reference coordinates (CPU EvalMono).
__device__ inline double flatMonomial(const FlatCharge& q,const double p[3]){
    const double* V=q.vertices;const double d[3]={p[0]-V[0],p[1]-V[1],p[2]-V[2]};
    if(q.kind==0){
        const double* I=q.inverse;
        const double l0=I[0]*d[0]+I[1]*d[1]+I[2]*d[2];
        const double l1=I[3]*d[0]+I[4]*d[1]+I[5]*d[2];
        const double l2=I[6]*d[0]+I[7]*d[1]+I[8]*d[2];
        return flatPow(l0,q.expo[0])*flatPow(l1,q.expo[1])*flatPow(l2,q.expo[2]);
    }
    const double a1d=(V[3]-V[0])*d[0]+(V[4]-V[1])*d[1]+(V[5]-V[2])*d[2];
    const double a2d=(V[6]-V[0])*d[0]+(V[7]-V[1])*d[1]+(V[8]-V[2])*d[2];
    const double* G=q.inverse;
    const double l0=G[0]*a1d+G[1]*a2d,l1=G[2]*a1d+G[3]*a2d;
    return flatPow(l0,q.expo[0])*flatPow(l1,q.expo[1]);
}
// Exact INT_host m(y)/|p-y| dy for cell degree <= 1 and face degree <= 2 (CPU PhiAtHO_Analytic).
__device__ inline double flatInnerPotential(const FlatCharge& q,const double p[3]){
    const double* c=q.polyCenter;const double lp[3]={p[0]-c[0],p[1]-c[1],p[2]-c[2]};
    if(q.kind==0){
        double V[4][3];for(int i=0;i<4;++i)for(int k=0;k<3;++k)V[i][k]=q.vertices[3*i+k]-c[k];
        const double I0=rad_hdiv_cuda::PhiTet(V,lp);
        if(q.degree==0)return q.a*I0;
        double M1[3];rad_hdiv_cuda::TetMoment1(V,lp,M1);
        return q.a*I0+q.b[0]*M1[0]+q.b[1]*M1[1]+q.b[2]*M1[2];
    }
    double V[3][3];for(int i=0;i<3;++i)for(int k=0;k<3;++k)V[i][k]=q.vertices[3*i+k]-c[k];
    const double I0=rad_hdiv_cuda::TriPotential(V,lp);
    if(q.degree==0)return q.a*I0;
    double M1[3];rad_hdiv_cuda::TriMoment1(V,lp,M1);
    double result=q.a*I0+q.b[0]*M1[0]+q.b[1]*M1[1]+q.b[2]*M1[2];
    if(q.degree>=2){double M2[3][3];rad_hdiv_cuda::TriMoment2(V,lp,M2);
        for(int k=0;k<3;++k)for(int l=0;l<3;++l)result+=q.c[3*k+l]*M2[k][l];}
    return result;
}
struct FlatDevice {
    const FlatCharge* charges;const double* points;const double* weights;const double* lowPoints;const double* lowWeights;
    std::size_t n,stride,lowStride;const int* outerCount;const int* lowCount;double farFactor;
    // Per-charge graded-rule source and per-host flags: nonzero when a host
    // carries any charge outside affine P1 (cell constant, face degree <= 1).
    const radia_hdiv_python::P1ChargeSource* p1;const int* cellNonP1;const int* faceNonP1;
    // Curved P2 hosts: node tables and the touching Duffy rule.
    bool curved;const double* cellNodes;const double* faceNodes;const double* gl;const double* gw;int nq;
    bool hex;HexDevice hexd;
};
__device__ inline double flatQuadDot(const FlatDevice& d,std::size_t target,std::size_t source){
    const auto* P=d.points+3*d.stride*target;const auto* W=d.weights+d.stride*target;double sum=0.0;
    for(int k=0;k<d.outerCount[target];++k)sum+=W[k]*flatInnerPotential(d.charges[source],P+3*k);
    return sum*0.07957747154594766788;
}
__device__ inline double flatQuadDotFar(const FlatDevice& d,std::size_t target,std::size_t source){
    const auto* Px=d.lowPoints+3*d.lowStride*target;const auto* Wx=d.lowWeights+d.lowStride*target;
    const auto* Py=d.lowPoints+3*d.lowStride*source;const auto* Wy=d.lowWeights+d.lowStride*source;double sum=0.0;
    for(int i=0;i<d.lowCount[target];++i){double inner=0.0;
        for(int j=0;j<d.lowCount[source];++j){
            const double dx=Px[3*i]-Py[3*j],dy=Px[3*i+1]-Py[3*j+1],dz=Px[3*i+2]-Py[3*j+2];
            const double r=sqrt(dx*dx+dy*dy+dz*dz);if(r<1e-300)continue;
            inner+=Wy[j]/r;
        }
        sum+=Wx[i]*inner;
    }
    return sum*0.07957747154594766788;
}
// CPU P1TouchHosts for direct pairs: both hosts carry only affine P1 charges
// and share a vertex with bitwise-equal coordinates (self pairs included).
__device__ inline bool flatP1Touch(const FlatDevice& d,const FlatCharge& a,const FlatCharge& b){
    if((a.kind==0?d.cellNonP1:d.faceNonP1)[a.host]||(b.kind==0?d.cellNonP1:d.faceNonP1)[b.host])return false;
    for(int i=0;i<(a.kind?3:4);++i)for(int j=0;j<(b.kind?3:4);++j)
        if(a.vertices[3*i]==b.vertices[3*j]&&a.vertices[3*i+1]==b.vertices[3*j+1]&&a.vertices[3*i+2]==b.vertices[3*j+2])return true;
    return false;
}
// CPU CurvedHostsTouch: two shared mesh vertices, or one when a face is involved.
__device__ inline bool curvedHostsTouch(const FlatCharge& a,const FlatCharge& b){
    int shared=0;for(int i=0;i<(a.kind?3:4);++i)for(int j=0;j<(b.kind?3:4);++j)if(a.numbers[i]==b.numbers[j]){++shared;break;}
    return shared>=2||(shared==1&&(a.kind||b.kind));
}
// Target outer rule against the singularity-resolving curved Duffy potential of the source host.
__device__ inline double curvedDuffyDot(const FlatDevice& d,std::size_t target,std::size_t source){
    const auto& s=d.charges[source];const auto* P=d.points+3*d.stride*target;const auto* W=d.weights+d.stride*target;double sum=0.0;
    const double* nodes=s.kind?d.faceNodes+18*std::size_t(s.host):d.cellNodes+30*std::size_t(s.host);
    for(int k=0;k<d.outerCount[target];++k){const double* p=P+3*k;
        const double inner=s.kind?rad_hdiv_cuda::CurvedTriPotential(reinterpret_cast<const double(*)[3]>(nodes),s.expo[0],s.expo[1],p,d.gl,d.gw,d.nq,true)
            :rad_hdiv_cuda::CurvedTetPotential(reinterpret_cast<const double(*)[3]>(nodes),s.expo[0],s.expo[1],s.expo[2],p,d.gl,d.gw,d.nq,true);
        sum+=W[k]*inner;}
    return sum*0.07957747154594766788;
}
// CPU QuadBlockHOCurvedDirect: high outer rule of both hosts, target-major sums.
__device__ inline double curvedDirectDot(const FlatDevice& d,std::size_t target,std::size_t source){
    const auto* Px=d.points+3*d.stride*target;const auto* Wx=d.weights+d.stride*target;
    const auto* Py=d.points+3*d.stride*source;const auto* Wy=d.weights+d.stride*source;double sum=0.0;
    for(int i=0;i<d.outerCount[target];++i){double inner=0.0;
        for(int j=0;j<d.outerCount[source];++j){
            const double dx=Px[3*i]-Py[3*j],dy=Px[3*i+1]-Py[3*j+1],dz=Px[3*i+2]-Py[3*j+2];
            const double r=sqrt(dx*dx+dy*dy+dz*dz);if(r<1e-300)continue;
            inner+=Wy[j]/r;
        }
        sum+=Wx[i]*inner;
    }
    return sum*0.07957747154594766788;
}
__device__ inline double flatEntry(const FlatDevice& d,std::size_t a,std::size_t b){
    if(d.hex)return hexEntry(d.hexd,a,b);
    const auto& qa=d.charges[a];const auto& qb=d.charges[b];
    if(a!=b&&d.farFactor<1e29){
        const double dx=qa.center[0]-qb.center[0],dy=qa.center[1]-qb.center[1],dz=qa.center[2]-qb.center[2];
        if(sqrt(dx*dx+dy*dy+dz*dz)>d.farFactor*(qa.size+qb.size))
            return 0.5*(flatQuadDotFar(d,a,b)+flatQuadDotFar(d,b,a));
    }
    if(d.curved){
        if(curvedHostsTouch(qa,qb))return a==b?curvedDuffyDot(d,a,a):0.5*(curvedDuffyDot(d,a,b)+curvedDuffyDot(d,b,a));
        return curvedDirectDot(d,a,b);
    }
    if(flatP1Touch(d,qa,qb))
        return a==b?radia_hdiv_p1_diagonal::pairEntry(d.p1[a],d.p1[a])
                   :0.5*(radia_hdiv_p1_diagonal::pairEntry(d.p1[a],d.p1[b])+radia_hdiv_p1_diagonal::pairEntry(d.p1[b],d.p1[a]));
    if(a==b)return flatQuadDot(d,a,a);
    return 0.5*(flatQuadDot(d,a,b)+flatQuadDot(d,b,a));
}
// One thread per charge: host maps, monomial-folded outer/low rules and the
// A + B.y + y^T C y physical expansion (CPU constructor + InitHOPolynomialCoefficients).
__global__ void prepareFlatCharges(std::size_t n,const double* cells,std::size_t cellCount,const double* faces,std::size_t faceCount,
        const int* hosts,const int* kinds,const int* exponents,
        const FlatRulePoint* tet,int nt,const FlatRulePoint* tri,int nr,const FlatRulePoint* tetLow,int ntl,const FlatRulePoint* triLow,int nrl,
        std::size_t stride,std::size_t lowStride,FlatCharge* charges,double* points,double* weights,int* outerCount,
        double* lowPoints,double* lowWeights,int* lowCount,radia_hdiv_python::P1ChargeSource* p1,
        int* cellNonP1,int* faceNonP1,int* invalid){
    const auto a=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(a>=n)return;
    FlatCharge q{};q.kind=kinds[a];q.host=hosts[a];
    for(int k=0;k<3;++k)q.expo[k]=exponents[3*a+k];
    if((q.kind!=0&&q.kind!=1)||q.host<0||q.expo[0]<0||q.expo[1]<0||q.expo[2]<0||
       (q.kind==0&&std::size_t(q.host)>=cellCount)||(q.kind==1&&(std::size_t(q.host)>=faceCount||q.expo[2]!=0))){atomicExch(invalid,1);return;}
    // BDM1/BDM2: cell degree <= 1, face degree <= 2. Higher modes need the Duffy/ladder kernels.
    // Bound each component before summing; this backstops FlatBlock::validate.
    const int limit=q.kind==0?1:2;
    if(q.expo[0]>limit||q.expo[1]>limit||q.expo[2]>limit){atomicExch(invalid,2);return;}
    q.degree=q.expo[0]+q.expo[1]+q.expo[2];
    if(q.degree>limit){atomicExch(invalid,2);return;}
    const int vertexCount=q.kind==0?4:3;const double* V=q.kind==0?cells+12*std::size_t(q.host):faces+9*std::size_t(q.host);
    for(int i=0;i<3*vertexCount;++i)q.vertices[i]=V[i];
    double measure,beta[3][3]={};int coordinates;
    if(q.kind==0){
        double E[9];for(int r=0;r<3;++r)for(int col=0;col<3;++col)E[r*3+col]=V[3*(col+1)+r]-V[r];
        const double det=E[0]*(E[4]*E[8]-E[5]*E[7])-E[1]*(E[3]*E[8]-E[5]*E[6])+E[2]*(E[3]*E[7]-E[4]*E[6]);
        const double iv=1.0/det;double* I=q.inverse;
        I[0]=(E[4]*E[8]-E[5]*E[7])*iv;I[1]=-(E[1]*E[8]-E[2]*E[7])*iv;I[2]=(E[1]*E[5]-E[2]*E[4])*iv;
        I[3]=-(E[3]*E[8]-E[5]*E[6])*iv;I[4]=(E[0]*E[8]-E[2]*E[6])*iv;I[5]=-(E[0]*E[5]-E[2]*E[3])*iv;
        I[6]=(E[3]*E[7]-E[4]*E[6])*iv;I[7]=-(E[0]*E[7]-E[1]*E[6])*iv;I[8]=(E[0]*E[4]-E[1]*E[3])*iv;
        measure=fabs(det);
        for(int i=0;i<3;++i)for(int k=0;k<3;++k)beta[i][k]=I[3*i+k];coordinates=3;
    }else{
        double a1[3],a2[3];for(int k=0;k<3;++k){a1[k]=V[3+k]-V[k];a2[k]=V[6+k]-V[k];}
        const double a1a2=a1[0]*a2[0]+a1[1]*a2[1]+a1[2]*a2[2];
        const double g[4]={a1[0]*a1[0]+a1[1]*a1[1]+a1[2]*a1[2],a1a2,a1a2,a2[0]*a2[0]+a2[1]*a2[1]+a2[2]*a2[2]};
        const double iv=1.0/(g[0]*g[3]-g[1]*g[2]);double* G=q.inverse;
        G[0]=g[3]*iv;G[1]=-g[1]*iv;G[2]=-g[2]*iv;G[3]=g[0]*iv;
        const double cr[3]={a1[1]*a2[2]-a1[2]*a2[1],a1[2]*a2[0]-a1[0]*a2[2],a1[0]*a2[1]-a1[1]*a2[0]};
        measure=2.0*(0.5*sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]));
        for(int k=0;k<3;++k){beta[0][k]=G[0]*a1[k]+G[1]*a2[k];beta[1][k]=G[2]*a1[k]+G[3]*a2[k];}coordinates=2;
    }
    if(!isfinite(measure)||!(measure>0)){atomicExch(invalid,1);return;}
    for(int i=0;i<9;++i)if(!isfinite(q.inverse[i])){atomicExch(invalid,1);return;}
    // Far/near centre and bounding radius; the polynomial centre sums before dividing.
    for(int i=0;i<vertexCount;++i)for(int k=0;k<3;++k)q.center[k]+=V[3*i+k]/double(vertexCount);
    for(int k=0;k<3;++k){for(int i=0;i<vertexCount;++i)q.polyCenter[k]+=V[3*i+k];q.polyCenter[k]/=vertexCount;}
    for(int i=0;i<vertexCount;++i){const double dx=V[3*i]-q.center[0],dy=V[3*i+1]-q.center[1],dz=V[3*i+2]-q.center[2];
        const double r=sqrt(dx*dx+dy*dy+dz*dz);if(r>q.size)q.size=r;}
    {
        double V0[3];for(int k=0;k<3;++k)V0[k]=V[k]-q.polyCenter[k];
        double factorA[2],factorB[2][3];int factors=0;
        for(int i=0;i<coordinates;++i)for(int c=0;c<q.expo[i];++c){
            for(int k=0;k<3;++k)factorB[factors][k]=beta[i][k];
            factorA[factors]=-(beta[i][0]*V0[0]+beta[i][1]*V0[1]+beta[i][2]*V0[2]);++factors;
        }
        double A=1.0,B[3]={0,0,0},C[3][3]={};
        for(int f=0;f<factors;++f){
            const double al=factorA[f];const double* be=factorB[f];
            const double nA=A*al;double nB[3],nC[3][3];
            for(int k=0;k<3;++k)nB[k]=A*be[k]+al*B[k];
            for(int k=0;k<3;++k)for(int l=0;l<3;++l)nC[k][l]=al*C[k][l]+0.5*(B[k]*be[l]+be[k]*B[l]);
            A=nA;for(int k=0;k<3;++k){B[k]=nB[k];for(int l=0;l<3;++l)C[k][l]=nC[k][l];}
        }
        q.a=A;for(int k=0;k<3;++k){q.b[k]=B[k];for(int l=0;l<3;++l)q.c[3*k+l]=C[k][l];}
    }
    auto map=[&](const FlatRulePoint& r,double* P){
        if(q.kind==0)for(int k=0;k<3;++k)P[k]=V[k]+r.xi[0]*(V[3+k]-V[k])+r.xi[1]*(V[6+k]-V[k])+r.xi[2]*(V[9+k]-V[k]);
        else for(int k=0;k<3;++k)P[k]=V[k]+r.xi[0]*(V[3+k]-V[k])+r.xi[1]*(V[6+k]-V[k]);
    };
    const FlatRulePoint* outer=q.kind==0?tet:tri;const int no=q.kind==0?nt:nr;
    for(int k=0;k<no;++k){double* P=points+3*(stride*a+k);map(outer[k],P);
        weights[stride*a+k]=outer[k].weight*measure*flatMonomial(q,P);}
    outerCount[a]=no;
    const FlatRulePoint* low=q.kind==0?tetLow:triLow;const int nl=q.kind==0?ntl:nrl;
    for(int k=0;k<nl;++k){double* P=lowPoints+3*(lowStride*a+k);map(low[k],P);
        lowWeights[lowStride*a+k]=low[k].weight*measure*flatMonomial(q,P);}
    lowCount[a]=nl;
    // CPU P1TouchHosts eligibility, folded per host; the graded source uses the shared CPU frame.
    if(q.expo[2]!=0||q.expo[0]+q.expo[1]>(q.kind==1?1:0))atomicExch((q.kind==0?cellNonP1:faceNonP1)+q.host,1);
    else if(!radia_hdiv_python::p1ChargeSourceInto(q.kind,q.vertices,q.expo,p1[a])||!isfinite(p1[a].constant)||
            !isfinite(p1[a].gradient[0])||!isfinite(p1[a].gradient[1])||!isfinite(p1[a].gradient[2])){atomicExch(invalid,1);return;}
    charges[a]=q;
}
// One thread per curved charge (CPU curved high-order constructor): centre is
// the corner mean, size the P2-node bounding radius, and every rule point is
// mapped with weight (rule weight x curved measure) x reference monomial.
__global__ void prepareCurvedCharges(std::size_t n,const double* cellNodes,std::size_t cellCount,const double* faceNodes,std::size_t faceCount,
        const int* cellNumbers,const int* faceNumbers,const int* hosts,const int* kinds,const int* exponents,
        const FlatRulePoint* tet,int nt,const FlatRulePoint* tri,int nr,const FlatRulePoint* tetLow,int ntl,const FlatRulePoint* triLow,int nrl,
        std::size_t stride,std::size_t lowStride,FlatCharge* charges,double* points,double* weights,int* outerCount,
        double* lowPoints,double* lowWeights,int* lowCount,int* invalid){
    const auto a=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(a>=n)return;
    FlatCharge q{};q.kind=kinds[a];q.host=hosts[a];
    for(int k=0;k<3;++k)q.expo[k]=exponents[3*a+k];
    if((q.kind!=0&&q.kind!=1)||q.host<0||q.expo[0]<0||q.expo[1]<0||q.expo[2]<0||
       (q.kind==0&&std::size_t(q.host)>=cellCount)||(q.kind==1&&(std::size_t(q.host)>=faceCount||q.expo[2]!=0))){atomicExch(invalid,1);return;}
    const int limit=q.kind==0?1:2;
    if(q.expo[0]>limit||q.expo[1]>limit||q.expo[2]>limit){atomicExch(invalid,2);return;}
    q.degree=q.expo[0]+q.expo[1]+q.expo[2];
    if(q.degree>limit){atomicExch(invalid,2);return;}
    const int corners=q.kind==0?4:3,nodes=q.kind==0?10:6;
    const double* X=q.kind==0?cellNodes+30*std::size_t(q.host):faceNodes+18*std::size_t(q.host);
    const int* numbers=q.kind==0?cellNumbers+4*std::size_t(q.host):faceNumbers+3*std::size_t(q.host);
    for(int i=0;i<corners;++i){q.numbers[i]=numbers[i];for(int k=0;k<3;++k)q.vertices[3*i+k]=X[3*i+k];}
    for(int i=0;i<corners;++i)for(int k=0;k<3;++k)q.center[k]+=X[3*i+k]/double(corners);
    for(int i=0;i<nodes;++i){const double dx=X[3*i]-q.center[0],dy=X[3*i+1]-q.center[1],dz=X[3*i+2]-q.center[2];
        const double r=sqrt(dx*dx+dy*dy+dz*dz);if(r>q.size)q.size=r;}
    auto mapped=[&](const FlatRulePoint& r,double* P,double& w){double measure;
        if(q.kind==0)rad_hdiv_cuda::CurvedTetMapMeasure(reinterpret_cast<const double(*)[3]>(X),r.xi[0],r.xi[1],r.xi[2],P,measure);
        else rad_hdiv_cuda::CurvedTriMapMeasure(reinterpret_cast<const double(*)[3]>(X),r.xi[0],r.xi[1],P,measure);
        double mono=flatPow(r.xi[0],q.expo[0])*flatPow(r.xi[1],q.expo[1]);if(q.kind==0)mono*=flatPow(r.xi[2],q.expo[2]);
        if(!isfinite(measure)||!(measure>0))atomicExch(invalid,1);
        w=r.weight*measure*mono;};
    const FlatRulePoint* outer=q.kind==0?tet:tri;const int no=q.kind==0?nt:nr;
    for(int k=0;k<no;++k)mapped(outer[k],points+3*(stride*a+k),weights[stride*a+k]);
    outerCount[a]=no;
    const FlatRulePoint* low=q.kind==0?tetLow:triLow;const int nl=q.kind==0?ntl:nrl;
    for(int k=0;k<nl;++k)mapped(low[k],lowPoints+3*(lowStride*a+k),lowWeights[lowStride*a+k]);
    lowCount[a]=nl;
    charges[a]=q;
}
__global__ void assembleFlatGram(FlatDevice d,std::size_t first,std::size_t count,double* matrix,int* invalid){
    const auto i=first+std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=first+count)return;
    const auto row=i/d.n,col=i%d.n;if(row>col)return;
    const double value=flatEntry(d,row,col);
    if(!isfinite(value))atomicExch(invalid,1);
    matrix[row*d.n+col]=value;if(row!=col)matrix[col*d.n+row]=value;
}
// Prepared charges shared by the dense and ACA routes: geometry, rules and
// the per-charge host maps live on the device; only primitives cross H2D.
struct FlatPrepared {
    Buffer<double> cells,faces;Buffer<int> hosts,kinds,exponents,invalid,outerCount,lowCount;
    Buffer<FlatRulePoint> tet,tri,tetLow,triLow;Buffer<FlatCharge> charges;Buffer<double> qp,qw,lp,lw;
    Buffer<radia_hdiv_python::P1ChargeSource> p1;Buffer<int> cellNonP1,faceNonP1;
    Buffer<double> gl,gw;Buffer<int> cellNumbers,faceNumbers;
    std::size_t n=0,stride=0,lowStride=0;double farFactor=1e30;bool curved=false,hex=false;HexPrepared hexPrepared;
    // Every buffer prepare() allocates, uploads included, so one reserve check
    // covers them; the HEX pair table is sized and checked inside prepare().
    static std::size_t workspaceBytes(const FlatBlock& block){
        const auto n=block.charges();
        if(block.hex){const std::size_t hosts=block.cellCount()+block.faceCount(),nOut=block.ruleOut.size();
            const auto lattices=checkedProduct(block.cells.size()+block.faces.size(),sizeof(double));
            const auto perHost=checkedProduct(hosts,checkedProduct(4*nOut*nOut*nOut+8,sizeof(double))+2*sizeof(int)+1+2*sizeof(std::size_t));
            const auto coefficients=checkedProduct(checkedProduct(block.cellCount(),27*84)+checkedProduct(block.faceCount(),9*35),sizeof(double));
            return flatSum(flatSum(lattices,perHost),flatSum(coefficients,checkedProduct(n,4*sizeof(int))+sizeof(int)));}
        const std::size_t stride=std::max(block.tet.size(),block.tri.size()),lowStride=std::max<std::size_t>(1,std::max(block.tetLow.size(),block.triLow.size()));
        const auto points=checkedProduct(checkedProduct(n,stride),4*sizeof(double)),low=checkedProduct(checkedProduct(n,lowStride),4*sizeof(double));
        const auto rules=checkedProduct(block.tet.size()+block.tri.size()+block.tetLow.size()+block.triLow.size(),sizeof(FlatRulePoint));
        const auto uploads=checkedProduct(block.cells.size()+block.faces.size(),sizeof(double))+checkedProduct(n,5*sizeof(int))+rules;
        const auto p1=flatSum(checkedProduct(n,sizeof(radia_hdiv_python::P1ChargeSource)),
            checkedProduct(block.cellCount()+block.faceCount(),sizeof(int)));
        const auto curved=flatSum(checkedProduct(block.gl.size()+block.gw.size(),sizeof(double)),
            checkedProduct(block.cellNumbers.size()+block.faceNumbers.size(),sizeof(int)));
        return flatSum(flatSum(flatSum(checkedProduct(n,sizeof(FlatCharge)+2*sizeof(int))+sizeof(int),flatSum(points,low)),uploads),flatSum(p1,curved));
    }
    int flag(cudaStream_t stream,HostVector<FlatTransfer>& transfers,const char* payload){int value=0;
        check(cudaMemcpyAsync(&value,invalid.data,sizeof(int),cudaMemcpyDeviceToHost,stream));
        check(cudaStreamSynchronize(stream));transfers.push_back({payload,"D2H",sizeof(int)});return value;}
    // allPairs: the dense route tabulates every HEX host pair; ACA and the
    // diagnostic only the near candidates and evaluate far entries on demand.
    void prepare(const FlatBlock& block,cudaStream_t stream,HostVector<FlatTransfer>& transfers,std::size_t& launches,int& invalidCode,bool allPairs){
        n=block.charges();farFactor=block.farFactor;curved=block.curved;hex=block.hex;
        if(hex){hexPrepared.prepare(block,allPairs,stream,[&](const char* payload,const char* direction,std::size_t bytes){
                transfers.push_back({payload,direction,bytes});},launches,invalidCode);
            invalid.upload(HostVector<int>{0});transfers.push_back({"invalid_flag","H2D",sizeof(int)});return;}
        stride=std::max(block.tet.size(),block.tri.size());lowStride=std::max<std::size_t>(1,std::max(block.tetLow.size(),block.triLow.size()));
        auto upload=[&](auto& buffer,const auto& values,const char* payload){
            if(values.empty())return;buffer.upload(values);
            transfers.push_back({payload,"H2D",checkedProduct(values.size(),sizeof(values[0]))});};
        upload(cells,block.cells,curved?"cell_p2_nodes_m":"cell_vertices_m");upload(faces,block.faces,curved?"face_p2_nodes_m":"face_vertices_m");
        upload(hosts,block.host,"charge_host");upload(kinds,block.kind,"charge_kind");upload(exponents,block.expo,"charge_exponents");
        upload(tet,block.tet,"tet_rule");upload(tri,block.tri,"tri_rule");upload(tetLow,block.tetLow,"tet_low_rule");upload(triLow,block.triLow,"tri_low_rule");
        upload(invalid,HostVector<int>{0},"invalid_flag");
        if(curved){upload(gl,block.gl,"curve_gauss_points");upload(gw,block.gw,"curve_gauss_weights");
            upload(cellNumbers,block.cellNumbers,"cell_vertex_numbers");upload(faceNumbers,block.faceNumbers,"face_vertex_numbers");}
        charges.allocate(n);outerCount.allocate(n);lowCount.allocate(n);
        qp.allocate(checkedProduct(n,3*stride));qw.allocate(checkedProduct(n,stride));
        lp.allocate(checkedProduct(n,3*lowStride));lw.allocate(checkedProduct(n,lowStride));
        // Every host flag is cleared on the device; a face-free or cell-free block keeps a one-int buffer.
        const std::size_t cellCount=block.cellCount(),faceCount=block.faceCount();
        p1.allocate(n);cellNonP1.allocate(std::max<std::size_t>(1,cellCount));faceNonP1.allocate(std::max<std::size_t>(1,faceCount));
        check(cudaMemsetAsync(cellNonP1.data,0,std::max<std::size_t>(1,cellCount)*sizeof(int),stream));
        check(cudaMemsetAsync(faceNonP1.data,0,std::max<std::size_t>(1,faceCount)*sizeof(int),stream));
        if(curved)prepareCurvedCharges<<<unsigned((n+127)/128),128,0,stream>>>(n,cells.data,cellCount,faces.data,faceCount,
            cellNumbers.data,faceNumbers.data,hosts.data,kinds.data,exponents.data,tet.data,int(block.tet.size()),tri.data,int(block.tri.size()),
            tetLow.data,int(block.tetLow.size()),triLow.data,int(block.triLow.size()),stride,lowStride,
            charges.data,qp.data,qw.data,outerCount.data,lp.data,lw.data,lowCount.data,invalid.data);
        else prepareFlatCharges<<<unsigned((n+127)/128),128,0,stream>>>(n,cells.data,cellCount,faces.data,faceCount,
            hosts.data,kinds.data,exponents.data,tet.data,int(block.tet.size()),tri.data,int(block.tri.size()),
            tetLow.data,int(block.tetLow.size()),triLow.data,int(block.triLow.size()),stride,lowStride,
            charges.data,qp.data,qw.data,outerCount.data,lp.data,lw.data,lowCount.data,p1.data,cellNonP1.data,faceNonP1.data,invalid.data);
        check(cudaGetLastError());++launches;
        if((invalidCode=flag(stream,transfers,"prepare_flag"))!=0)
            throw std::invalid_argument(invalidCode==2?"Flat device Gram supports cell degree <= 1 and face degree <= 2 only"
                                                       :"Invalid flat charge host, kind, exponent or degenerate geometry");
    }
    FlatDevice device()const{return {charges.data,qp.data,qw.data,lp.data,lw.data,n,stride,lowStride,outerCount.data,lowCount.data,farFactor,
        p1.data,cellNonP1.data,faceNonP1.data,curved,cells.data,faces.data,gl.data,gw.data,int(gl.size),hex,hexPrepared.d};}
};
inline void flatReserve(std::size_t bytes,const char* what){
    std::size_t free,total;check(cudaMemGetInfo(&free,&total));
    if(!fitsAboveReserve(free,total,bytes))throw std::runtime_error(what);
}
// Deferred-dense oracle for Resident::materializeDeferred: writes the
// normalized Gram and its sigma = sqrt(raw diagonal) on the resident stream.
class FlatGramOracle {
    const FlatBlock& block;
public:
    HostVector<FlatTransfer> transfers{HostFaultAllocator<FlatTransfer>("flat.transfers")};
    std::size_t ownedBytes=0,launches=0;int invalidCode=0;
    explicit FlatGramOracle(const FlatBlock& b):block(b){block.validate();}
    template<class Rule>static std::size_t ruleBytes(const Rule& r){return checkedProduct(r.size(),sizeof(FlatRulePoint));}
    void fillDense(const Plan& p,double* coefficients,const std::size_t*,double* scale,cudaStream_t stream){
        const auto n=block.charges();
        if(p.n!=n||p.directions.size()!=1||p.deferredDense.size()!=1||!p.deferredLowrank.empty()||p.coefficientCount!=checkedProduct(n,n))
            throw std::invalid_argument("Flat Gram requires one deferred dense charge block");
        const auto& span=p.deferredDense[0];
        if(span.row||span.column||span.rows!=n||span.columns!=n||span.offset)
            throw std::invalid_argument("Flat Gram dense span mismatch");
        for(std::size_t i=0;i<n;++i)if(p.permutation[i]!=i)throw std::invalid_argument("Flat Gram requires identity charge order");
        ownedBytes=FlatPrepared::workspaceBytes(block);flatReserve(ownedBytes,"Flat Gram GPU reserve would be exceeded");
        FlatPrepared prepared;prepared.prepare(block,stream,transfers,launches,invalidCode,true);
        const FlatDevice device=prepared.device();
        // Bounded launches keep each kernel short on display-attached GPUs.
        const auto entries=checkedProduct(n,n);constexpr std::size_t batch=std::size_t(1)<<16;
        if((entries+255)/256>std::size_t(INT_MAX))throw std::invalid_argument("Flat dense Gram exceeds the CUDA grid range");
        for(std::size_t first=0;first<entries;first+=batch){const auto here=std::min(batch,entries-first);
            assembleFlatGram<<<unsigned((here+127)/128),128,0,stream>>>(device,first,here,coefficients,prepared.invalid.data);
            check(cudaGetLastError());check(cudaStreamSynchronize(stream));++launches;}
        analyticPhysicalScale<<<unsigned((n+255)/256),256,0,stream>>>(n,coefficients,scale,prepared.invalid.data);check(cudaGetLastError());++launches;
        if((invalidCode=prepared.flag(stream,transfers,"energy_flag"))!=0)throw std::runtime_error("Flat device Gram produced a nonfinite entry or nonpositive self-energy");
        normalizeAnalytic<<<unsigned((entries+255)/256),256,0,stream>>>(n,coefficients,scale);check(cudaGetLastError());++launches;
        check(cudaStreamSynchronize(stream));
    }
};
} // namespace radia_hmatrix_cuda
