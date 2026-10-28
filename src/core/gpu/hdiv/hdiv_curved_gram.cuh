// P2 reference-monomial Gram using the existing curved panel kernels.
// Touch classification and symmetric entry arithmetic mirror the CPU path.
#pragma once
#include <memory>
#include "hmatrix_host_fault.hpp"
namespace radia_hmatrix_cuda {
struct CurvedCharge {double nodes[30];int vertices[4],expo[3],kind,host;bool affine=false;};
struct ReferencePoint {double xi[3],weight;};
struct CurvedBlock {
    HostVector<CurvedCharge> charges{HostFaultAllocator<CurvedCharge>("geometry.charges")};
    HostVector<ReferencePoint> tet{HostFaultAllocator<ReferencePoint>("geometry.tet-rule")},tri{HostFaultAllocator<ReferencePoint>("geometry.tri-rule")};
    HostVector<double> gl{HostFaultAllocator<double>("geometry.gl")},gw{HostFaultAllocator<double>("geometry.gw")};
    bool hostVectorCache=false,hostMomentCache=false,profileAssembly=false;
};
struct HostTargetPoint {std::size_t representative;ReferencePoint quadrature;};
struct HostTargetSpan {std::size_t first,count;};
struct MappedHostPoint {double position[3];};
// Group only identical mapped hosts. Monomial weights stay charge-specific.
template<class PointAllocator,class SpanAllocator>
inline void curvedHostTargets(const CurvedBlock& block,std::vector<HostTargetPoint,PointAllocator>& points,std::vector<HostTargetSpan,SpanAllocator>& spans){
    using IndexAllocator=typename std::allocator_traits<PointAllocator>::template rebind_alloc<std::size_t>;
    std::vector<std::size_t,IndexAllocator> representatives{IndexAllocator(points.get_allocator())};spans.resize(block.charges.size());
    for(std::size_t i=0;i<block.charges.size();++i){const auto& a=block.charges[i];bool found=false;
        for(auto representative:representatives){const auto& b=block.charges[representative];if(a.kind!=b.kind||a.host!=b.host)continue;
            if(a.affine!=b.affine)throw std::invalid_argument("Inconsistent geometry for charge host");
            for(int k=0;k<(a.kind?18:30);++k)if(a.nodes[k]!=b.nodes[k])throw std::invalid_argument("Inconsistent mapped charge host");
            for(int k=0;k<(a.kind?3:4);++k)if(a.vertices[k]!=b.vertices[k])throw std::invalid_argument("Inconsistent host topology");
            spans[i]=spans[representative];found=true;break;
        }
        if(!found){const auto& rule=a.kind?block.tri:block.tet;spans[i]={points.size(),rule.size()};representatives.push_back(i);for(const auto& q:rule)points.push_back({i,q});}
    }
}
__device__ bool curvedTouch(const CurvedCharge& a,const CurvedCharge& b){
    int shared=0;for(int i=0;i<(a.kind?3:4);++i)for(int j=0;j<(b.kind?3:4);++j)if(a.vertices[i]==b.vertices[j]){++shared;break;}
    return shared>=2||(shared==1&&(a.kind||b.kind));
}
__device__ double curvedMonomial(const CurvedCharge& a,const double* xi){
    double g[3]={xi[0],xi[1],xi[2]};if(a.affine){g[0]=xi[1];g[1]=a.kind?1.-xi[0]-xi[1]:xi[2];g[2]=1.-xi[0]-xi[1]-xi[2];}
    double result=1.;for(int c=0;c<3;++c)for(int k=0;k<a.expo[c];++k)result*=g[c];return result;
}
__device__ double affineReferencePotential(const CurvedCharge& a,const double* p){
    const int degree=a.expo[0]+a.expo[1]+a.expo[2];double moments[20]={};
    if(!a.kind){double v[4][3];for(int i=0;i<4;++i)for(int c=0;c<3;++c)v[i][c]=a.nodes[3*((i+1)%4)+c];
        rad_hdiv_cuda::TetReferencePotentialMomentsUpTo(v,p,degree,moments);
    }else{
        double v[3][3],edge[2][3],inverse[3][3]={};for(int i=0;i<3;++i)for(int c=0;c<3;++c)v[i][c]=a.nodes[3*((i+1)%3)+c];
        double aa=0.,bb=0.,ab=0.;for(int c=0;c<3;++c){edge[0][c]=v[1][c]-v[0][c];edge[1][c]=v[2][c]-v[0][c];aa+=edge[0][c]*edge[0][c];bb+=edge[1][c]*edge[1][c];ab+=edge[0][c]*edge[1][c];}
        const double determinant=aa*bb-ab*ab;if(determinant<=0)return nan("");
        for(int c=0;c<3;++c){inverse[0][c]=(bb*edge[0][c]-ab*edge[1][c])/determinant;inverse[1][c]=(aa*edge[1][c]-ab*edge[0][c])/determinant;}
        rad_hdiv_cuda::SurfaceReferencePotentialMomentsUpTo(v,p,degree,v[0],inverse,moments);
    }
    return moments[rad_hdiv_cuda::PotentialMomentIndex(a.expo[0],a.expo[1],a.expo[2])];
}
__device__ void curvedMappedPoint(const CurvedCharge& a,const ReferencePoint& q,double* x,double& weight){
    double measure;
    if(a.kind)rad_hdiv_cuda::CurvedTriMapMeasure(reinterpret_cast<const double(*)[3]>(a.nodes),q.xi[0],q.xi[1],x,measure);
    else rad_hdiv_cuda::CurvedTetMapMeasure(reinterpret_cast<const double(*)[3]>(a.nodes),q.xi[0],q.xi[1],q.xi[2],x,measure);
    weight=q.weight*measure*curvedMonomial(a,q.xi);
}
__device__ double curvedSourcePotential(const CurvedCharge& target,const CurvedCharge& source,const double* p,
        const ReferencePoint* tet,std::size_t nt,const ReferencePoint* tri,std::size_t ns,
        const double* gl,const double* gw,int nq){
    const bool touching=curvedTouch(target,source);double inner=0.;
        if(source.affine)inner=affineReferencePotential(source,p);
        else if(touching){
            if(source.kind)inner=rad_hdiv_cuda::CurvedTriPotential(reinterpret_cast<const double(*)[3]>(source.nodes),source.expo[0],source.expo[1],p,gl,gw,nq,true);
            else inner=rad_hdiv_cuda::CurvedTetPotential(reinterpret_cast<const double(*)[3]>(source.nodes),source.expo[0],source.expo[1],source.expo[2],p,gl,gw,nq,true);
        }else{
            const auto* sq=source.kind?tri:tet;const auto m=source.kind?ns:nt;
            for(std::size_t s=0;s<m;++s){double x[3],w;curvedMappedPoint(source,sq[s],x,w);
                const double dx=p[0]-x[0],dy=p[1]-x[1],dz=p[2]-x[2],r=sqrt(dx*dx+dy*dy+dz*dz);if(r>=1e-300)inner+=w/r;
            }
        }
    return inner;
}
__device__ double curvedQuadDot(const CurvedCharge& target,const CurvedCharge& source,
        const ReferencePoint* tet,std::size_t nt,const ReferencePoint* tri,std::size_t ns,
        const double* gl,const double* gw,int nq){
    const auto* tq=target.kind?tri:tet;const auto n=target.kind?ns:nt;double sum=0.;
    for(std::size_t t=0;t<n;++t){double p[3],weight;curvedMappedPoint(target,tq[t],p,weight);
        sum+=weight*curvedSourcePotential(target,source,p,tet,nt,tri,ns,gl,gw,nq);
    }
    return sum*0.07957747154594766788;
}
#include "hdiv_curved_moments.cuh"
// Geometry and monomial weights depend only on a fixed mapped host/charge.
// Hoist them without changing quadrature, Duffy treatment or reduction order.
__global__ void mapCurvedHostPoints(const CurvedCharge* charges,const HostTargetPoint* points,std::size_t np,MappedHostPoint* mapped){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=np)return;
    double weight;curvedMappedPoint(charges[points[i].representative],points[i].quadrature,mapped[i].position,weight);
}
__global__ void mapCurvedChargeWeights(const CurvedCharge* charges,std::size_t n,const HostTargetPoint* points,const HostTargetSpan* spans,std::size_t stride,double* weights){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=n*stride)return;const auto row=i/stride,k=i%stride;
    if(k>=spans[row].count){weights[i]=0.;return;}double position[3],weight;
    curvedMappedPoint(charges[row],points[spans[row].first+k].quadrature,position,weight);weights[i]=weight;
}
__device__ double mappedCurvedSourcePotential(const CurvedCharge* charges,std::size_t target,std::size_t source,const double* p,
    const ReferencePoint* tet,std::size_t nt,const ReferencePoint* tri,std::size_t ns,const double* gl,const double* gw,int nq,
    const HostTargetSpan* spans,const MappedHostPoint* mapped,const double* weights,std::size_t stride){
    if(charges[source].affine||curvedTouch(charges[target],charges[source]))
        return curvedSourcePotential(charges[target],charges[source],p,tet,nt,tri,ns,gl,gw,nq);
    const auto span=spans[source];double inner=0.;
    for(std::size_t k=0;k<span.count;++k){const auto* x=mapped[span.first+k].position;const double w=weights[source*stride+k];
        const double dx=p[0]-x[0],dy=p[1]-x[1],dz=p[2]-x[2],r=sqrt(dx*dx+dy*dy+dz*dz);if(r>=1e-300)inner+=w/r;
    }
    return inner;
}
__global__ void cacheCurvedHostPotentials(const CurvedCharge* charges,std::size_t n,const HostTargetPoint* points,std::size_t np,
    const ReferencePoint* tet,std::size_t nt,const ReferencePoint* tri,std::size_t ns,const double* gl,const double* gw,int nq,
    std::size_t first,std::size_t count,std::size_t sourceFirst,const HostTargetSpan* spans,const MappedHostPoint* mapped,const double* weights,std::size_t stride,double* cache,int* invalid){
    const auto local=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(local>=count)return;const auto i=first+local;if(i>=n*np)return;
    const auto source=sourceFirst+i/np,point=i%np;if(source>=n)return;const auto& q=points[point];const auto* p=mapped[point].position;
    const double value=mappedCurvedSourcePotential(charges,q.representative,source,p,tet,nt,tri,ns,gl,gw,nq,spans,mapped,weights,stride);if(!isfinite(value))atomicExch(invalid,1);cache[i]=value;
}
template<class Allocator=std::allocator<MomentSourceGroup>>
inline std::vector<MomentSourceGroup,Allocator> curvedSourceMomentGroups(const CurvedBlock& block,const Allocator& allocator=Allocator{}){
    std::vector<MomentSourceGroup,Allocator> groups(allocator);
    for(std::size_t i=0;i<block.charges.size();++i){const auto& charge=block.charges[i];bool found=false;
        for(auto& group:groups){const auto& representative=block.charges[group.member[0]];if(charge.kind!=representative.kind||charge.host!=representative.host)continue;
            if(group.count==6)throw std::invalid_argument("Fused source moments support at most six monomials per host");
            group.member[group.count++]=i;found=true;break;
        }
        if(!found){MomentSourceGroup group{};group.member[0]=i;group.count=1;groups.push_back(group);}
    }
    return groups;
}
template<int SourceKind>
__global__ void cacheCurvedHostMomentPotentials(const CurvedCharge* charges,std::size_t n,const HostTargetPoint* points,std::size_t np,
    const HostTargetSpan* spans,const MappedHostPoint* mapped,const double* weights,std::size_t stride,
    const MomentSourceGroup* groups,const std::size_t* groupIds,std::size_t groupCount,const double* gl,const double* gw,int nq,
    std::size_t first,std::size_t count,std::size_t sourceFirst,std::size_t sourceCount,double* cache,int* invalid){
    const auto local=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(local>=count)return;const auto i=first+local;if(i>=groupCount*np)return;
    const auto& group=groups[groupIds[i/np]];const auto point=i%np;const auto target=points[point].representative;const auto* p=mapped[point].position;
    const auto firstSource=group.member[0];double values[6]={};
    if(charges[firstSource].affine){for(std::size_t m=0;m<group.count;++m)values[m]=affineReferencePotential(charges[group.member[m]],p);}
    else if(curvedTouch(charges[target],charges[firstSource])){
        if constexpr(SourceKind==1)curvedTriMomentPotentials(charges,group,p,gl,gw,nq,values);
        else curvedTetMomentPotentials(charges,group,p,gl,gw,nq,values);
    }else{
        const auto span=spans[firstSource];
        for(std::size_t k=0;k<span.count;++k){const auto* x=mapped[span.first+k].position;const double dx=p[0]-x[0],dy=p[1]-x[1],dz=p[2]-x[2],r=sqrt(dx*dx+dy*dy+dz*dz);if(r>=1e-300)
            for(std::size_t m=0;m<group.count;++m)values[m]+=weights[group.member[m]*stride+k]/r;
        }
    }
    for(std::size_t m=0;m<group.count;++m){const auto source=group.member[m];if(source<sourceFirst||source>=sourceFirst+sourceCount)continue;
        if(!isfinite(values[m]))atomicExch(invalid,1);cache[(source-sourceFirst)*np+point]=values[m];
    }
}
__device__ double curvedCachedDot(const CurvedCharge* charges,std::size_t target,std::size_t source,
    const HostTargetPoint* points,std::size_t np,const HostTargetSpan* spans,const double* cache,const double* weights,std::size_t stride){
    const auto span=spans[target];double sum=0.;for(std::size_t k=0;k<span.count;++k){const auto point=span.first+k;sum+=weights[target*stride+k]*cache[source*np+point];}
    return sum*0.07957747154594766788;
}
__global__ void assembleCachedCurvedGram(const CurvedCharge* charges,std::size_t n,const HostTargetPoint* points,std::size_t np,
    const HostTargetSpan* spans,const double* cache,const double* weights,std::size_t stride,double* matrix,int* invalid){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=n*n)return;const auto row=i/n,col=i%n;if(row>col)return;
    double value=curvedCachedDot(charges,row,col,points,np,spans,cache,weights,stride);if(row!=col)value=.5*(value+curvedCachedDot(charges,col,row,points,np,spans,cache,weights,stride));
    if(!isfinite(value))atomicExch(invalid,1);matrix[i]=value;if(row!=col)matrix[col*n+row]=value;
}
// Store each ordered entry before the original symmetric averaging. Source
// tiling changes cache residency only, not source integration or target sums.
__global__ void assembleCachedCurvedColumns(const CurvedCharge* charges,std::size_t n,const HostTargetPoint* points,std::size_t np,
    const HostTargetSpan* spans,const double* cache,const double* weights,std::size_t stride,std::size_t sourceFirst,std::size_t sourceCount,double* matrix,int* invalid){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=n*sourceCount)return;
    const auto row=i/sourceCount,localCol=i%sourceCount,col=sourceFirst+localCol;
    const auto value=curvedCachedDot(charges,row,localCol,points,np,spans,cache,weights,stride);
    if(!isfinite(value))atomicExch(invalid,1);matrix[row*n+col]=value;
}
__global__ void symmetrizeCachedCurved(std::size_t n,double* matrix,int* invalid){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=n*n)return;const auto row=i/n,col=i%n;
    if(row>=col)return;const double value=.5*(matrix[i]+matrix[col*n+row]);
    if(!isfinite(value))atomicExch(invalid,1);matrix[i]=value;matrix[col*n+row]=value;
}
__global__ void assembleCurvedGram(const CurvedCharge* charges,std::size_t n,
    const ReferencePoint* tet,std::size_t nt,const ReferencePoint* tri,std::size_t ns,
    const double* gl,const double* gw,int nq,std::size_t first,std::size_t count,double* matrix,int* invalid){
    const auto local=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(local>=count)return;const auto i=first+local;if(i>=n*n)return;
    const auto row=i/n,col=i%n;if(row>col)return;
    double value=curvedQuadDot(charges[row],charges[col],tet,nt,tri,ns,gl,gw,nq);
    if(row!=col)value=.5*(value+curvedQuadDot(charges[col],charges[row],tet,nt,tri,ns,gl,gw,nq));
    if(!isfinite(value))atomicExch(invalid,1);matrix[i]=value;if(row!=col)matrix[col*n+row]=value;
}
}
