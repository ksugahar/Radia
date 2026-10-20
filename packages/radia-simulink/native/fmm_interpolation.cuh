/* SPDX-License-Identifier: MPL-2.0 */
// Independently authored tensor interpolation algebra for Laplace FMM.
// Mathematical reference: Fong & Darve, JCP 228 (2009), 8712-8725,
// doi:10.1016/j.jcp.2009.08.031. No upstream implementation is copied.
#pragma once
#include <cmath>
#ifdef __CUDACC__
#define FMM_HD __host__ __device__
#else
#define FMM_HD
#endif
namespace radia_fmm {
constexpr int MAX_ORDER=8;
struct Node {
    double center[3],half;
    int children[8],parent,octant,depth,sourceBegin,sourceCount,targetBegin,targetCount;
    FMM_HD bool leaf()const {for(int i=0;i<8;++i)if(children[i]>=0)return false;return true;}
};
FMM_HD inline double abscissa(int i,int p){return cos(3.141592653589793238462643383279502884*(i+.5)/p);}
FMM_HD inline void basis(double x,int p,double* v,double* d){
    double nodes[MAX_ORDER],w[MAX_ORDER];int hit=-1;
    for(int i=0;i<p;++i){double t=3.141592653589793238462643383279502884*(i+.5)/p;nodes[i]=cos(t);w[i]=(i%2?-1.:1.)*sin(t);
        if(fabs(x-nodes[i])<8*2.220446049250313e-16)hit=i;}
    if(hit>=0){double sum=0;for(int i=0;i<p;++i){v[i]=i==hit?1.:0.;
        d[i]=i==hit?0.:w[i]/(w[hit]*(nodes[hit]-nodes[i]));sum+=d[i];}d[hit]=-sum;return;}
    double sum=0,second=0;for(int i=0;i<p;++i){v[i]=w[i]/(x-nodes[i]);sum+=v[i];second+=v[i]/(x-nodes[i]);}
    for(int i=0;i<p;++i){v[i]/=sum;d[i]=v[i]*(second/sum-1./(x-nodes[i]));}
}
FMM_HD inline void point(const Node& node,int k,int p,double* x){
    x[0]=node.center[0]+node.half*abscissa(k%p,p);
    x[1]=node.center[1]+node.half*abscissa((k/p)%p,p);
    x[2]=node.center[2]+node.half*abscissa(k/(p*p),p);
}
FMM_HD inline void direct(const double* x,const double* y,double q,const double* dip,double& potential,double* field){
    double r[3]={x[0]-y[0],x[1]-y[1],x[2]-y[2]};double r2=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
    double inv=1./sqrt(r2),inv3=inv/r2,dr=dip[0]*r[0]+dip[1]*r[1]+dip[2]*r[2];
    potential+=q*inv+dr*inv3;
    for(int k=0;k<3;++k)field[k]+=q*r[k]*inv3+3*r[k]*dr*inv3/r2-dip[k]*inv3;
}
}
#undef FMM_HD
