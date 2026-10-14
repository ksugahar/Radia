// P2 reference-monomial Gram using the existing curved panel kernels.
// Touch classification and symmetric entry arithmetic mirror the CPU path.
#pragma once
namespace radia_hmatrix_cuda {
struct CurvedCharge {double nodes[30];int vertices[4],expo[3],kind,host;bool affine=false;};
struct ReferencePoint {double xi[3],weight;};
struct CurvedBlock {std::vector<CurvedCharge> charges;std::vector<ReferencePoint> tet,tri;std::vector<double> gl,gw;bool hostVectorCache=false;};
struct HostTargetPoint {std::size_t representative;ReferencePoint quadrature;};
struct HostTargetSpan {std::size_t first,count;};
// Group only identical mapped hosts. Monomial weights stay charge-specific.
inline void curvedHostTargets(const CurvedBlock& block,std::vector<HostTargetPoint>& points,std::vector<HostTargetSpan>& spans){
    std::vector<std::size_t> representatives;spans.resize(block.charges.size());
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
__global__ void cacheCurvedHostPotentials(const CurvedCharge* charges,std::size_t n,const HostTargetPoint* points,std::size_t np,
    const ReferencePoint* tet,std::size_t nt,const ReferencePoint* tri,std::size_t ns,const double* gl,const double* gw,int nq,
    std::size_t first,std::size_t count,double* cache,int* invalid){
    const auto local=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(local>=count)return;const auto i=first+local;if(i>=n*np)return;
    const auto source=i/np,point=i%np;const auto& q=points[point];const auto& target=charges[q.representative];double p[3],weight;curvedMappedPoint(target,q.quadrature,p,weight);
    const double value=curvedSourcePotential(target,charges[source],p,tet,nt,tri,ns,gl,gw,nq);if(!isfinite(value))atomicExch(invalid,1);cache[i]=value;
}
__device__ double curvedCachedDot(const CurvedCharge* charges,std::size_t target,std::size_t source,
    const HostTargetPoint* points,std::size_t np,const HostTargetSpan* spans,const double* cache){
    const auto span=spans[target];double sum=0.;for(std::size_t k=0;k<span.count;++k){const auto point=span.first+k;double p[3],weight;curvedMappedPoint(charges[target],points[point].quadrature,p,weight);sum+=weight*cache[source*np+point];}
    return sum*0.07957747154594766788;
}
__global__ void assembleCachedCurvedGram(const CurvedCharge* charges,std::size_t n,const HostTargetPoint* points,std::size_t np,
    const HostTargetSpan* spans,const double* cache,double* matrix,int* invalid){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=n*n)return;const auto row=i/n,col=i%n;if(row>col)return;
    double value=curvedCachedDot(charges,row,col,points,np,spans,cache);if(row!=col)value=.5*(value+curvedCachedDot(charges,col,row,points,np,spans,cache));
    if(!isfinite(value))atomicExch(invalid,1);matrix[i]=value;if(row!=col)matrix[col*n+row]=value;
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
