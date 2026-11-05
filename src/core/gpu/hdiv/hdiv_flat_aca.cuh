// Flat-TET charge Gram compressed on the device by partial-pivot adaptive
// cross approximation. The host builds only geometric metadata: a cluster
// tree over the charge host elements (centre plus support radius, longest-axis
// median bisection) and HACApK-style box-gap admissibility on support boxes.
// A cluster ends as a leaf at leafSize charges, the indices HACApK counts; a host is never split.
// It is an independent partition, not the CPU HACApK leaf set. Every Gram
// entry, sigma, dense leaf and low-rank factor is generated and kept on the
// device; only per-leaf integer ranks/status and step counts cross D2H.
#pragma once
#include <functional>
#include <numeric>
#include "hdiv_flat_gram.cuh"
#include "hdiv_flat_expansion.cuh"
namespace radia_hmatrix_cuda {
// rankCapDense: a leaf that reaches maxRank unconverged is stored as an exact
// dense leaf of the same operator (counted); otherwise construction fails.
// expansionOrder: Chebyshev nodes per axis of the kernel expansion of far blocks; 0 runs ACA on Gram entries.
// With it, admissibility is max(width) <= eta gap on tight boxes, because interpolation needs both boxes separated.
struct FlatAcaOptions {
    double eps=0,eta=0;std::size_t leafSize=0,maxRank=0;bool rankCapDense=false;int expansionOrder=0;
    void validate()const{
        if(!std::isfinite(eps)||!(eps>0)||!(eps<1))throw std::invalid_argument("ACA eps must be finite in (0,1)");
        if(!std::isfinite(eta)||!(eta>0))throw std::invalid_argument("ACA admissibility eta must be finite and positive");
        if(!leafSize||leafSize>4096)throw std::invalid_argument("ACA leaf size must be in [1,4096] charges");
        if(!maxRank||maxRank>1024)throw std::invalid_argument("ACA max rank must be in [1,1024]");
        if(expansionOrder&&(expansionOrder<3||expansionOrder>hcaOrderCap))throw std::invalid_argument("Kernel expansion order must be in [3,8]");
    }
};
struct FlatAcaBlock {std::size_t row,column,rows,columns;bool admissible;std::size_t rowBox,columnBox;};
struct FlatAcaPartition {
    HostVector<HcaBox> boxes{HostFaultAllocator<HcaBox>("aca.boxes")};
    HostVector<std::size_t> permutation{HostFaultAllocator<std::size_t>("aca.permutation")};
    HostVector<FlatAcaBlock> blocks{HostFaultAllocator<FlatAcaBlock>("aca.blocks")};
    std::size_t elements=0,clusters=0;
};
// Geometry-only partition. Charges of one host element stay contiguous and in
// ascending input order; equal centres or zero extent end the split as a leaf.
inline FlatAcaPartition flatAcaPartition(const FlatBlock& block,const FlatAcaOptions& options){
    block.validate();options.validate();
    const auto n=block.charges();const std::size_t cellCount=block.cellCount(),faceCount=block.faceCount();
    constexpr std::size_t none=std::numeric_limits<std::size_t>::max();
    HostVector<std::size_t> cellElement(cellCount,none,HostFaultAllocator<std::size_t>("aca.cell-element")),
        faceElement(faceCount,none,HostFaultAllocator<std::size_t>("aca.face-element"));
    for(std::size_t a=0;a<n;++a)(block.kind[a]==0?cellElement:faceElement)[std::size_t(block.host[a])]=0;
    HostVector<double> centre{HostFaultAllocator<double>("aca.centre")},radius{HostFaultAllocator<double>("aca.radius")},tight{HostFaultAllocator<double>("aca.tight")};
    std::size_t elements=0;
    // Centre: corner mean; support radius: every geometry node (curved P2 nodes bound the curved host).
    auto addElement=[&](const double* V,int vertices,int nodes){
        double c[3]={0,0,0};for(int i=0;i<vertices;++i)for(int k=0;k<3;++k)c[k]+=V[3*i+k];
        for(int k=0;k<3;++k)c[k]/=vertices;double r=0;
        for(int i=0;i<nodes;++i){const double dx=V[3*i]-c[0],dy=V[3*i+1]-c[1],dz=V[3*i+2]-c[2];r=std::max(r,std::sqrt(dx*dx+dy*dy+dz*dz));}
        centre.insert(centre.end(),c,c+3);radius.push_back(r);return elements++;};
    // HEX lattices: the corner mean of nodes {0,2,6,8,18,20,24,26} ({0,2,6,8} on a quad), radius over all nodes.
    auto addHex=[&](const double* L,int corners){static constexpr int node[8]={0,2,6,8,18,20,24,26};double C[24];
        for(int i=0;i<corners;++i)for(int k=0;k<3;++k)C[3*i+k]=L[3*node[i]+k];
        double c[3]={0,0,0};for(int i=0;i<corners;++i)for(int k=0;k<3;++k)c[k]+=C[3*i+k];for(int k=0;k<3;++k)c[k]/=corners;double r=0;
        for(int i=0;i<(corners==8?27:9);++i){const double dx=L[3*i]-c[0],dy=L[3*i+1]-c[1],dz=L[3*i+2]-c[2];r=std::max(r,std::sqrt(dx*dx+dy*dy+dz*dz));}
        if(options.expansionOrder)for(int k=0;k<3;++k){double lo=L[k],hi=L[k];
            for(int i=1;i<(corners==8?27:9);++i){lo=std::min(lo,L[3*i+k]);hi=std::max(hi,L[3*i+k]);}
            tight.push_back(lo);tight.push_back(hi);}
        centre.insert(centre.end(),c,c+3);radius.push_back(r);return elements++;};
    for(std::size_t h=0;h<cellCount;++h)if(cellElement[h]!=none)cellElement[h]=block.hex?addHex(block.cells.data()+81*h,8)
        :addElement(block.cells.data()+block.cellStride()*h,4,int(block.cellStride()/3));
    for(std::size_t h=0;h<faceCount;++h)if(faceElement[h]!=none)faceElement[h]=block.hex?addHex(block.faces.data()+27*h,4)
        :addElement(block.faces.data()+block.faceStride()*h,3,int(block.faceStride()/3));
    HostVector<std::size_t> offsets(elements+1,0,HostFaultAllocator<std::size_t>("aca.charge-offsets")),members(n,0,HostFaultAllocator<std::size_t>("aca.charge-members"));
    auto elementOf=[&](std::size_t a){return (block.kind[a]==0?cellElement:faceElement)[std::size_t(block.host[a])];};
    for(std::size_t a=0;a<n;++a)++offsets[elementOf(a)+1];
    for(std::size_t e=0;e<elements;++e)offsets[e+1]+=offsets[e];
    {auto next=offsets;for(std::size_t a=0;a<n;++a)members[next[elementOf(a)]++]=a;}
    struct Node {std::size_t first,count,charge=0,charges=0,left=0,right=0;double lo[3],hi[3],width=0,tlo[3],thi[3],twidth=0;bool leaf=true;};
    HostVector<Node> nodes{HostFaultAllocator<Node>("aca.nodes")};
    HostVector<std::size_t> order(elements,0,HostFaultAllocator<std::size_t>("aca.order"));std::iota(order.begin(),order.end(),std::size_t(0));
    std::function<std::size_t(std::size_t,std::size_t)> build=[&](std::size_t first,std::size_t count)->std::size_t{
        Node node{first,count};double clo[3],chi[3];std::size_t nodeCharges=0;
        for(int k=0;k<3;++k){node.lo[k]=clo[k]=std::numeric_limits<double>::infinity();node.hi[k]=chi[k]=-std::numeric_limits<double>::infinity();}
        for(std::size_t i=first;i<first+count;++i){const auto e=order[i];nodeCharges+=offsets[e+1]-offsets[e];
            for(int k=0;k<3;++k){const double c=centre[3*e+k];clo[k]=std::min(clo[k],c);chi[k]=std::max(chi[k],c);
                node.lo[k]=std::min(node.lo[k],c-radius[e]);node.hi[k]=std::max(node.hi[k],c+radius[e]);}}
        double w=0;for(int k=0;k<3;++k)w+=(node.hi[k]-node.lo[k])*(node.hi[k]-node.lo[k]);node.width=std::sqrt(w);
        if(options.expansionOrder){double tw=0;
            for(int k=0;k<3;++k){node.tlo[k]=std::numeric_limits<double>::infinity();node.thi[k]=-std::numeric_limits<double>::infinity();
                for(std::size_t i=first;i<first+count;++i){node.tlo[k]=std::min(node.tlo[k],tight[6*order[i]+2*k]);node.thi[k]=std::max(node.thi[k],tight[6*order[i]+2*k+1]);}
                tw+=(node.thi[k]-node.tlo[k])*(node.thi[k]-node.tlo[k]);}
            node.twidth=std::sqrt(tw);}
        int axis=0;for(int k=1;k<3;++k)if(chi[k]-clo[k]>chi[axis]-clo[axis])axis=k;
        const auto index=nodes.size();nodes.push_back(node);
        if(nodeCharges<=options.leafSize||!(chi[axis]>clo[axis]))return index;
        std::sort(order.begin()+std::ptrdiff_t(first),order.begin()+std::ptrdiff_t(first+count),[&](std::size_t a,std::size_t b){
            const double ca=centre[3*a+axis],cb=centre[3*b+axis];return ca<cb||(ca==cb&&a<b);});
        const auto half=count/2;const auto left=build(first,half),right=build(first+half,count-half);
        nodes[index].leaf=false;nodes[index].left=left;nodes[index].right=right;return index;};
    build(0,elements);
    FlatAcaPartition result;result.elements=elements;result.clusters=nodes.size();
    HostVector<std::size_t> start(elements+1,0,HostFaultAllocator<std::size_t>("aca.order-start"));
    for(std::size_t i=0;i<elements;++i){const auto e=order[i];start[i+1]=start[i]+offsets[e+1]-offsets[e];
        result.permutation.insert(result.permutation.end(),members.begin()+std::ptrdiff_t(offsets[e]),members.begin()+std::ptrdiff_t(offsets[e+1]));}
    if(result.permutation.size()!=n)throw std::runtime_error("ACA partition does not cover every charge exactly once");
    for(auto& node:nodes){node.charge=start[node.first];node.charges=start[node.first+node.count]-node.charge;}
    if(options.expansionOrder)for(const auto& node:nodes){HcaBox box{};for(int k=0;k<3;++k){box.lo[k]=node.tlo[k];box.hi[k]=node.thi[k];}
        box.charge=node.charge;box.charges=node.charges;result.boxes.push_back(box);}
    auto separated=[&](const Node& s,const double* lo,const double* hi,double width){
        double gap=0;for(int k=0;k<3;++k){const double d=std::max({0.0,lo[k]-s.hi[k],s.lo[k]-hi[k]});gap+=d*d;}
        gap=std::sqrt(gap);return gap>0&&std::min(s.width,width)<=options.eta*gap;};
    // Every term of a low-rank ACA block must be smooth: the column cluster and the box hull of each of its images are
    // admissible against the row cluster. The list is closed under inversion, so the image of the row cluster is covered.
    // The kernel expansion takes no images (FlatAcaOracle refuses them), so its branch reads the direct boxes alone.
    auto admissible=[&](const Node& s,const Node& t){
        if(options.expansionOrder){double tgap=0;for(int k=0;k<3;++k){const double d=std::max({0.0,t.tlo[k]-s.thi[k],s.tlo[k]-t.thi[k]});tgap+=d*d;}
            tgap=std::sqrt(tgap);return tgap>0&&std::max(s.twidth,t.twidth)<=options.eta*tgap;}
        if(!separated(s,t.lo,t.hi,t.width))return false;
        for(const auto& image:block.images){double lo[3],hi[3];
            for(int k=0;k<3;++k){lo[k]=std::numeric_limits<double>::infinity();hi[k]=-std::numeric_limits<double>::infinity();}
            for(int corner=0;corner<8;++corner){
                const double x=corner&1?t.hi[0]:t.lo[0],y=corner&2?t.hi[1]:t.lo[1],z=corner&4?t.hi[2]:t.lo[2];
                const double p[3]={image.flip[0]*(image.c*x-image.s*y),image.flip[1]*(image.s*x+image.c*y),image.flip[2]*z};
                for(int k=0;k<3;++k){lo[k]=std::min(lo[k],p[k]);hi[k]=std::max(hi[k],p[k]);}}
            double width=0;for(int k=0;k<3;++k)width+=(hi[k]-lo[k])*(hi[k]-lo[k]);
            if(!separated(s,lo,hi,std::sqrt(width)))return false;}
        return true;};
    // Upper/diagonal blocks only: s precedes t in charge order whenever s != t.
    std::function<void(std::size_t,std::size_t)> split=[&](std::size_t si,std::size_t ti){
        const Node s=nodes[si],t=nodes[ti];
        if(si==ti){
            if(s.leaf){result.blocks.push_back({s.charge,s.charge,s.charges,s.charges,false,si,ti});return;}
            split(s.left,s.left);split(s.left,s.right);split(s.right,s.right);return;}
        if(admissible(s,t)){result.blocks.push_back({s.charge,t.charge,s.charges,t.charges,true,si,ti});return;}
        if(s.leaf&&t.leaf){result.blocks.push_back({s.charge,t.charge,s.charges,t.charges,false,si,ti});return;}
        if(s.leaf){split(si,t.left);split(si,t.right);return;}
        if(t.leaf){split(s.left,ti);split(s.right,ti);return;}
        split(s.left,t.left);split(s.left,t.right);split(s.right,t.left);split(s.right,t.right);};
    split(0,0);
    // Coverage: the upper triangle including the diagonal, each entry once.
    std::size_t covered=0;for(const auto& b:result.blocks)covered+=b.row==b.column?b.rows*(b.rows+1)/2:b.rows*b.columns;
    if(covered!=n*(n+1)/2)throw std::runtime_error("ACA block partition does not cover the upper Gram exactly once");
    return result;
}
// Per-leaf ACA state; offsets index the batch scratch. Norms stay on the device.
struct AcaLeafState {
    std::size_t row,column,rows,columns,u,v,residual,used;
    double frobenius2,maxAbs;int rank,status,pivotRow,steps,cap,factorRank,probe;
};
enum AcaStatus:int {acaActive=0,acaConverged=1,acaRankCap=2,acaNonfinite=3,acaDenseSmaller=4};
constexpr int acaThreads=128;
// Rows that must agree with the factors before a leaf with a small last term stops.
constexpr int acaProbes=8;
constexpr std::size_t acaNone=~std::size_t(0);  // device-safe "no index"
__device__ inline double acaEntry(const FlatDevice& d,const std::size_t* perm,const double* sigma,std::size_t i,std::size_t j){
    const auto a=perm[i],b=perm[j];return flatEntry(d,a,b)*((1./sigma[a])*(1./sigma[b]));
}
__device__ inline double acaBlockSum(double x,double* shared){
    shared[threadIdx.x]=x;__syncthreads();
    for(int s=acaThreads/2;s>0;s>>=1){if(int(threadIdx.x)<s)shared[threadIdx.x]+=shared[threadIdx.x+s];__syncthreads();}
    const double r=shared[0];__syncthreads();return r;
}
// Deterministic argmax: larger value wins, equal values keep the lower index.
__device__ inline void acaBlockArgmax(double& value,std::size_t& index,double* sv,std::size_t* si){
    sv[threadIdx.x]=value;si[threadIdx.x]=index;__syncthreads();
    for(int s=acaThreads/2;s>0;s>>=1){
        if(int(threadIdx.x)<s){const double o=sv[threadIdx.x+s];const auto oi=si[threadIdx.x+s];
            if(o>sv[threadIdx.x]||(o==sv[threadIdx.x]&&oi<si[threadIdx.x])){sv[threadIdx.x]=o;si[threadIdx.x]=oi;}}
        __syncthreads();}
    value=sv[0];index=si[0];__syncthreads();
}
// Lowest unused row other than skip, or none.
__device__ inline std::size_t acaNextRow(const AcaLeafState& L,const double* used,std::size_t skip,double* sv,std::size_t* si){
    double value=-1;std::size_t index=acaNone;
    for(auto r=std::size_t(threadIdx.x);r<L.rows;r+=acaThreads)if(r!=skip&&used[r]==0.0){value=0;index=r;break;}
    acaBlockArgmax(value,index,sv,si);return index;
}
// The unused row nearest the evenly spaced target of one probe, or none; the same row in every thread.
__device__ inline std::size_t acaProbeRow(const AcaLeafState& L,const double* used,std::size_t skip,int probe,double* sv,std::size_t* si){
    const double target=double(2*probe-1)*double(L.rows)/double(2*acaProbes);
    double value=-double(L.rows)-1.0;std::size_t index=acaNone;
    for(auto r=std::size_t(threadIdx.x);r<L.rows;r+=acaThreads){const double distance=fabs(double(r)-target);
        if(r!=skip&&used[r]==0.0&&-distance>value){value=-distance;index=r;}}
    acaBlockArgmax(value,index,sv,si);return index;
}
// What one ACA step reads: Gram entries of the permuted charges.
struct AcaGramEntry {FlatDevice d;const std::size_t* perm;const double* sigma;
    __device__ double operator()(const AcaLeafState& L,std::size_t i,std::size_t j)const{return acaEntry(d,perm,sigma,L.row+i,L.column+j);}};
// The kernel between the Chebyshev nodes of two clusters; row and column name the clusters here.
struct AcaKernelEntry {HcaNodes e;const HcaBox* boxes;
    __device__ double operator()(const AcaLeafState& L,std::size_t i,std::size_t j)const{double a[3],b[3];
        hcaNode(e,boxes[L.row],int(i),a);hcaNode(e,boxes[L.column],int(j),b);
        const double dx=a[0]-b[0],dy=a[1]-b[1],dz=a[2]-b[2];return hexInvFourPi/sqrt(dx*dx+dy*dy+dz*dz);}};
// The product of low-rank factors, recompressed; row and column are the offsets of the factors here.
struct AcaFactorEntry {const double* factors;
    __device__ double operator()(const AcaLeafState& L,std::size_t i,std::size_t j)const{double sum=0.0;
        for(int r=0;r<L.factorRank;++r)sum+=factors[L.row+i+L.rows*r]*factors[L.column+j+L.columns*r];return sum;}};
// One ACA step for every active leaf of the batch (one block per leaf). A
// negligible pivot row is marked used and the next unused row is tried; the
// leaf converges by row exhaustion only after every row was tried. A small last
// term does not stop a leaf by itself: where clusters are symmetric, singular
// pairs come nearly equal and that term misses the next one. Probe rows must then
// agree with the factors, the mean row of the remainder at eps; the first that
// does not becomes the next pivot, or sends a leaf at its rank cap dense. These
// rows are samples: passing them is evidence, not a bound on the whole leaf.
template<class Entry>__global__ void acaStep(Entry entry,AcaLeafState* leaves,double* scratch,double eps,int* active){
    __shared__ double sv[acaThreads];__shared__ std::size_t si[acaThreads];
    AcaLeafState L=leaves[blockIdx.x];if(L.status!=acaActive)return;
    const std::size_t i=std::size_t(L.pivotRow),k=std::size_t(L.rank),rows=L.rows,columns=L.columns;
    double* U=scratch+L.u;double* V=scratch+L.v;double* R=scratch+L.residual;double* used=scratch+L.used;
    double best=-1,rowAbs=0;std::size_t jstar=acaNone;int bad=0;
    for(auto j=std::size_t(threadIdx.x);j<columns;j+=acaThreads){
        const double e=entry(L,i,j);if(!isfinite(e))bad=1;rowAbs=fmax(rowAbs,fabs(e));
        double r=e;for(std::size_t l=0;l<k;++l)r-=U[i+rows*l]*V[j+columns*l];
        R[j]=r;if(fabs(r)>best){best=fabs(r);jstar=j;}
    }
    acaBlockArgmax(best,jstar,sv,si);
    {double flag=bad;std::size_t ignored=0;acaBlockArgmax(flag,ignored,sv,si);bad=flag>0;}
    {std::size_t ignored=0;acaBlockArgmax(rowAbs,ignored,sv,si);}
    const double maxAbs=fmax(L.maxAbs,rowAbs);
    if(bad||!isfinite(best)){if(threadIdx.x==0){L.status=acaNonfinite;leaves[blockIdx.x]=L;}return;}
    if(L.probe>0){
        double row2=0;for(auto j=std::size_t(threadIdx.x);j<columns;j+=acaThreads)row2+=R[j]*R[j];
        row2=acaBlockSum(row2,sv);
        if(row2<=eps*eps*L.frobenius2/double(rows)){
            const auto next=acaProbeRow(L,used,i,L.probe-1,sv,si);
            if(threadIdx.x==0){used[i]=1.0;L.maxAbs=maxAbs;++L.steps;
                if(L.probe==1||next==acaNone)L.status=acaConverged;else{--L.probe;L.pivotRow=int(next);atomicAdd(active,1);}
                leaves[blockIdx.x]=L;}
            return;
        }
        if(L.rank>=L.cap){if(threadIdx.x==0){used[i]=1.0;L.maxAbs=maxAbs;++L.steps;L.status=acaRankCap;leaves[blockIdx.x]=L;}return;}
    }
    if(!(best>1e-15*maxAbs)){
        const auto next=acaNextRow(L,used,i,sv,si);
        if(threadIdx.x==0){used[i]=1.0;L.maxAbs=maxAbs;++L.steps;
            if(next==acaNone)L.status=acaConverged;else{L.pivotRow=int(next);atomicAdd(active,1);}
            leaves[blockIdx.x]=L;}
        return;
    }
    const double pivot=R[jstar];
    double nv2=0,nu2=0,columnAbs=0,nextBest=-1;std::size_t nextRow=acaNone;
    for(auto j=std::size_t(threadIdx.x);j<columns;j+=acaThreads){const double v=R[j]/pivot;V[j+columns*k]=v;nv2+=v*v;}
    for(auto r=std::size_t(threadIdx.x);r<rows;r+=acaThreads){
        const double e=entry(L,r,jstar);if(!isfinite(e))bad=1;columnAbs=fmax(columnAbs,fabs(e));
        double c=e;for(std::size_t l=0;l<k;++l)c-=U[r+rows*l]*V[jstar+columns*l];
        U[r+rows*k]=c;nu2+=c*c;
        if(r!=i&&used[r]==0.0&&fabs(c)>nextBest){nextBest=fabs(c);nextRow=r;}
    }
    nv2=acaBlockSum(nv2,sv);nu2=acaBlockSum(nu2,sv);
    {double flag=bad;std::size_t ignored=0;acaBlockArgmax(flag,ignored,sv,si);bad=flag>0;}
    {std::size_t ignored=0;acaBlockArgmax(columnAbs,ignored,sv,si);}
    acaBlockArgmax(nextBest,nextRow,sv,si);
    // ||S_k||_F^2 = ||S_{k-1}||_F^2 + 2 sum_l (u_l.u_k)(v_l.v_k) + ||u_k||^2 ||v_k||^2.
    double cross=0;
    for(std::size_t l=0;l<k;++l){double uu=0,vv=0;
        for(auto r=std::size_t(threadIdx.x);r<rows;r+=acaThreads)uu+=U[r+rows*l]*U[r+rows*k];
        for(auto j=std::size_t(threadIdx.x);j<columns;j+=acaThreads)vv+=V[j+columns*l]*V[j+columns*k];
        cross+=acaBlockSum(uu,sv)*acaBlockSum(vv,sv);}
    const double term=nu2*nv2,frobenius2=L.frobenius2+2*cross+term;
    const bool settled=sqrt(term)<=eps*sqrt(frobenius2);
    const auto probeRow=settled?acaProbeRow(L,used,i,acaProbes,sv,si):acaNone;
    if(threadIdx.x==0){
        used[i]=1.0;L.maxAbs=fmax(maxAbs,columnAbs);L.rank=int(k+1);++L.steps;L.frobenius2=frobenius2;
        if(bad||!isfinite(term)||!isfinite(frobenius2)||!(frobenius2>0))L.status=acaNonfinite;
        else if(std::size_t(L.rank)==(rows<columns?rows:columns))L.status=acaConverged;
        else if(settled&&probeRow==acaNone)L.status=acaConverged;
        else if(settled){L.probe=acaProbes;L.pivotRow=int(probeRow);atomicAdd(active,1);}
        else if(L.rank>=L.cap)L.status=acaRankCap;
        else if(nextRow==acaNone)L.status=acaConverged;
        // Factors of this rank already take as much as the dense leaf, and the rank only grows.
        else if(std::size_t(L.rank)*(rows+columns)>=rows*columns)L.status=acaDenseSmaller;
        else{L.probe=0;L.pivotRow=int(nextRow);atomicAdd(active,1);}
        leaves[blockIdx.x]=L;
    }
}
__global__ void acaSigma(FlatDevice d,double* sigma,int* invalid){
    const auto a=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(a>=d.n)return;
    double value=0.0;if(!flatChargeScale(d,a,flatEntry(d,a,a),value)){atomicExch(invalid,1);return;}
    sigma[a]=value;
}
__global__ void acaReport(const AcaLeafState* leaves,std::size_t count,int* report){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=count)return;
    report[2*i]=leaves[i].rank;report[2*i+1]=leaves[i].status;
}
struct AcaDenseSpan {std::size_t row,column,rows,columns,offset,first;};
__global__ void acaDenseFill(FlatDevice d,const std::size_t* perm,const double* sigma,const AcaDenseSpan* spans,std::size_t count,
                             std::size_t first,std::size_t entries,double* coefficients,int* invalid){
    const auto g=first+std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(g>=first+entries)return;
    std::size_t lo=0,hi=count;while(hi-lo>1){const auto mid=(lo+hi)/2;if(spans[mid].first<=g)lo=mid;else hi=mid;}
    const auto& s=spans[lo];const auto local=g-s.first,i=local/s.columns,j=local%s.columns;
    const double value=acaEntry(d,perm,sigma,s.row+i,s.column+j);if(!isfinite(value))atomicExch(invalid,1);
    coefficients[s.offset+i*s.columns+j]=value;
}
// One host of the permuted charge order, and a dense leaf as a grid of directed host pairs.
struct AcaHostStart {std::size_t start;int host;};
struct AcaPairSpan {std::size_t row,column,columns,offset,firstPair,rowHost,columnHost,columnHosts;};
// HEX and mixed hosts, one directed host pair per item. Each term of a pair, direct or image, reads the table when
// near; a far one evaluates its tensor block once, where entry by entry it would repeat the double sum for every
// charge pair. A far image term folds both directed blocks, like the entry route and the CPU GetHexSymBlock.
__global__ void __launch_bounds__(hexThreads) acaDensePairs(HexDevice d,const AcaHostStart* hosts,const AcaPairSpan* spans,std::size_t count,const double* sigma,
        std::size_t first,std::size_t items,double* scratch,std::size_t stride,std::size_t blockEntries,double* coefficients,int* invalid){
    for(std::size_t item=first+blockIdx.x;item<first+items;item+=gridDim.x){
        std::size_t lo=0,hi=count;while(hi-lo>1){const auto mid=(lo+hi)/2;if(spans[mid].firstPair<=item)lo=mid;else hi=mid;}
        const AcaPairSpan s=spans[lo];const auto local=item-s.firstPair;
        const AcaHostStart target=hosts[s.rowHost+local/s.columnHosts],source=hosts[s.columnHost+local%s.columnHosts];
        const int t=target.host,h=source.host,nT=d.count[t],nS=d.count[h];
        double* sum=scratch+std::size_t(blockIdx.x)*stride;double* block=sum+blockEntries;
        double* reverse=d.imageCount?block+blockEntries:block;double* work=reverse+blockEntries;
        for(int e=threadIdx.x;e<nT*nS;e+=blockDim.x)sum[e]=0.0;
        for(int img=0;img<=d.imageCount;++img){const double sign=img?d.images[img-1].sign:1.0;
            if(hexCandidate(d,t,h,img)){
                for(int e=threadIdx.x;e<nT*nS;e+=blockDim.x)sum[e]+=sign*hexTermEntry(d,std::size_t(d.first[t]+e/nS),std::size_t(d.first[h]+e%nS),img);
            }else{
                hexTensorBlock(d,t,h,img,hexFarClass,nT,nS,work,block,invalid);__syncthreads();
                if(img){hexTensorBlock(d,h,t,img,hexFarClass,nS,nT,work,reverse,invalid);__syncthreads();}
                for(int e=threadIdx.x;e<nT*nS;e+=blockDim.x)sum[e]+=sign*(img?0.5*(block[e]+reverse[(e%nS)*nT+e/nS]):block[e]);
            }
            __syncthreads();}
        const auto base=s.offset+(target.start-s.row)*s.columns+(source.start-s.column);
        for(int e=threadIdx.x;e<nT*nS;e+=blockDim.x){const int lt=e/nS,ls=e%nS;const auto a=std::size_t(d.first[t]+lt),b=std::size_t(d.first[h]+ls);
            const double value=sum[e]*((1./sigma[a])*(1./sigma[b]));if(!isfinite(value))atomicExch(invalid,1);
            coefficients[base+std::size_t(lt)*s.columns+ls]=value;}
        __syncthreads();
    }
}
struct FlatAcaStats {
    std::size_t elements=0,clusters=0,denseLeaves=0,lowrankLeaves=0,admissibleBlocks=0,rankCapPromotions=0,storagePromotions=0,zeroRankLeaves=0,
        rankSum=0,maxRank=0,steps=0,factorCoefficients=0,denseCoefficients=0,batches=0,farPairs=0,midPairs=0,gradedPairs=0,
        farMicroseconds=0,midMicroseconds=0,gradedMicroseconds=0,acaMicroseconds=0,denseMicroseconds=0,
        kernelRankSum=0,momentCoefficients=0,momentMicroseconds=0,expansionOrder=0;
};
// Compress before the resident plan exists (ranks size the plan), then
// materialize into the resident coefficients on its stream. Factors move D2D.
class FlatAcaOracle {
    const FlatBlock& block;FlatAcaOptions options;FlatAcaPartition partition;
    FlatPrepared prepared;Buffer<double> sigma;cudaStream_t stream=nullptr;
    struct Kept {std::size_t block,rank,chunk,u,v,state;};
    HostVector<Kept> kept{HostFaultAllocator<Kept>("aca.kept")};
    HostVector<HostUnique<Buffer<double>>> chunks{HostFaultAllocator<HostUnique<Buffer<double>>>("aca.chunks")};
    HostVector<std::size_t> rankOf{HostFaultAllocator<std::size_t>("aca.rank-of")};
    bool compressed=false;std::size_t activeBytes=0;
    Buffer<int> active;std::size_t statesBytes=0,reportBytes=0,activeCount=0;
    // Upload the leaves of one batch, step them to the end and read back their (rank, status) pairs.
    template<class Entry>HostVector<int> runLeaves(const Entry& entry,const HostVector<AcaLeafState>& states,double eps,Buffer<double>& work,std::size_t scratch){
        check(cudaMemsetAsync(work.data,0,scratch*sizeof(double),stream));
        Buffer<AcaLeafState> deviceStates;deviceStates.upload(states);statesBytes+=checkedProduct(states.size(),sizeof(AcaLeafState));
        int running=int(states.size());
        // A step adds rank or retires one row, so rows + cap steps bound every leaf.
        std::size_t bound=0;for(const auto& s:states)bound=std::max(bound,s.rows+std::size_t(s.cap));
        for(std::size_t step=0;running&&step<bound;++step){
            check(cudaMemsetAsync(active.data,0,sizeof(int),stream));
            acaStep<Entry><<<unsigned(states.size()),acaThreads,0,stream>>>(entry,deviceStates.data,work.data,eps,active.data);
            check(cudaGetLastError());++launches;
            check(cudaMemcpyAsync(&running,active.data,sizeof(int),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));
            activeCount+=sizeof(int);++stats.steps;
        }
        if(running)throw std::runtime_error("ACA did not terminate within its step bound");
        Buffer<int> report;report.allocate(2*states.size());
        acaReport<<<unsigned((states.size()+127)/128),128,0,stream>>>(deviceStates.data,states.size(),report.data);check(cudaGetLastError());++launches;
        HostVector<int> ranks(2*states.size(),0,HostFaultAllocator<int>("aca.report"));
        check(cudaMemcpyAsync(ranks.data(),report.data,ranks.size()*sizeof(int),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));
        reportBytes+=ranks.size()*sizeof(int);
        return ranks;
    }
    // Keep a converged leaf's factors, or count why it is stored dense. Returns whether it was kept.
    bool keepFactors(int rankValue,int status,std::size_t blockIndex,std::size_t state,std::size_t& chunkCount,HostVector<Kept>& batchKept){
        const auto rank=std::size_t(rankValue);const auto& b=partition.blocks[blockIndex];
        if(status==acaNonfinite)throw std::runtime_error("ACA produced a nonfinite pivot, entry or norm");
        if(status==acaRankCap){if(!options.rankCapDense)throw std::runtime_error("ACA leaf reached max_rank without converging");++stats.rankCapPromotions;return false;}
        if(status==acaDenseSmaller){++stats.storagePromotions;return false;}
        if(status!=acaConverged||rank>std::min(b.rows,b.columns))throw std::runtime_error("ACA leaf status or rank invalid");
        // Rank 0 means every row was negligible: the leaf is stored dense, never dropped.
        if(!rank){++stats.zeroRankLeaves;return false;}
        if(checkedProduct(rank,b.rows+b.columns)>=checkedProduct(b.rows,b.columns)){++stats.storagePromotions;return false;}
        batchKept.push_back({blockIndex,rank,chunks.size(),chunkCount,chunkCount+b.rows*rank,state});
        chunkCount+=rank*(b.rows+b.columns);rankOf[blockIndex]=rank;
        return true;
    }
    // U (rows x rank) and V (columns x rank) are leading slices of the scratch columns.
    void moveFactors(const HostVector<AcaLeafState>& states,const HostVector<Kept>& batchKept,const Buffer<double>& work,std::size_t chunkCount,std::size_t& live){
        if(!chunkCount)return;
        flatReserve(checkedProduct(chunkCount,sizeof(double)),"Flat ACA factor storage would exceed the GPU reserve");
        auto chunk=makeHostUnique<Buffer<double>>("aca.chunk");chunk->allocate(chunkCount);
        for(const auto& k:batchKept){const auto& st=states[k.state];const auto& b=partition.blocks[k.block];
            check(cudaMemcpyAsync(chunk->data+k.u,work.data+st.u,b.rows*k.rank*sizeof(double),cudaMemcpyDeviceToDevice,stream));
            check(cudaMemcpyAsync(chunk->data+k.v,work.data+st.v,b.columns*k.rank*sizeof(double),cudaMemcpyDeviceToDevice,stream));}
        check(cudaStreamSynchronize(stream));live+=checkedProduct(chunkCount,sizeof(double));chunks.push_back(std::move(chunk));
        kept.insert(kept.end(),batchKept.begin(),batchKept.end());
    }
    // The hosts of the permuted charge order, and the host index at each permuted charge.
    void permutedHosts(HostVector<AcaHostStart>& hosts,HostVector<std::size_t>& hostAt)const{
        const auto n=block.charges();const int cellHosts=int(block.cellCount());
        hostAt.assign(n+1,0);
        for(std::size_t i=0;i<n;++i){const auto a=partition.permutation[i];const int g=block.kind[a]==0?block.host[a]:cellHosts+block.host[a];
            if(hosts.empty()||hosts.back().host!=g)hosts.push_back({i,g});
            hostAt[i]=hosts.size()-1;}
        hostAt[n]=hosts.size();
    }
    std::size_t processors()const{int gpu=0;check(cudaGetDevice(&gpu));cudaDeviceProp properties;check(cudaGetDeviceProperties(&properties,gpu));
        return std::size_t(std::max(properties.multiProcessorCount,1));}
    // Far blocks from the kernel at the Chebyshev nodes of their clusters, batch by batch: ACA on K, the moments times
    // its factors, then ACA on that product, which also proves the block low rank at the requested tolerance.
    void compressExpansion(const HostVector<std::size_t>& admissible,std::size_t& live,std::size_t& peakScratch){
        const int n1=options.expansionOrder,q=n1*n1*n1;const HcaNodes nodes=hcaNodes(n1);const HexDevice hd=prepared.device().hexd;
        stats.expansionOrder=std::size_t(n1);
        HostVector<HcaBox> boxes=partition.boxes;HostVector<unsigned char> used(boxes.size(),0,HostFaultAllocator<unsigned char>("aca.box-used"));
        for(const auto b:admissible){used[partition.blocks[b].rowBox]=1;used[partition.blocks[b].columnBox]=1;}
        std::size_t tableSize=0;for(std::size_t b=0;b<boxes.size();++b)if(used[b]){boxes[b].table=tableSize;tableSize=flatSum(tableSize,checkedProduct(boxes[b].charges,std::size_t(q)));}
        stats.momentCoefficients=tableSize;
        HostVector<AcaHostStart> hosts{HostFaultAllocator<AcaHostStart>("aca.dense-hosts")};
        HostVector<std::size_t> hostAt(HostFaultAllocator<std::size_t>("aca.dense-host-index"));permutedHosts(hosts,hostAt);
        HostVector<HcaMomentItem> items{HostFaultAllocator<HcaMomentItem>("aca.moment-items")};
        for(std::size_t b=0;b<boxes.size();++b)if(used[b])for(auto h=hostAt[boxes[b].charge];h<hostAt[boxes[b].charge+boxes[b].charges];++h)items.push_back({b,hosts[h].start,hosts[h].host});
        const std::size_t blocks=processors(),maxCount=hd.order==2?27:8,nOut=std::size_t(hd.out.n)*hd.out.n*hd.out.n,stride=nOut*(maxCount+3*n1);
        const std::size_t expansionBytes=flatSum(flatSum(checkedProduct(tableSize,sizeof(double)),checkedProduct(boxes.size(),sizeof(HcaBox))),
            flatSum(checkedProduct(items.size(),sizeof(HcaMomentItem)),checkedProduct(checkedProduct(blocks,stride),sizeof(double))));
        flatReserve(expansionBytes,"Kernel expansion moments would exceed the GPU reserve");
        Buffer<double> moments;moments.allocate(tableSize);
        Buffer<HcaBox> deviceBoxes;deviceBoxes.upload(boxes);transfers.push_back({"hca_boxes","H2D",checkedProduct(boxes.size(),sizeof(HcaBox))});
        {const auto started=std::chrono::steady_clock::now();
            Buffer<HcaMomentItem> deviceItems;deviceItems.upload(items);transfers.push_back({"hca_moment_items","H2D",checkedProduct(items.size(),sizeof(HcaMomentItem))});
            Buffer<double> scratch;scratch.allocate(checkedProduct(blocks,stride));
            HexLaunchPace pace(blocks);
            for(std::size_t item=0;item<items.size();){const auto here=pace.begin(items.size()-item);
                hcaMoments<<<unsigned(std::min(blocks,here)),hexThreads,0,stream>>>(hd,nodes,deviceBoxes.data,deviceItems.data,sigma.data,item,here,scratch.data,stride,moments.data,prepared.invalid.data);
                check(cudaGetLastError());check(cudaStreamSynchronize(stream));++launches;pace.end();item+=here;}
            if((invalidCode=prepared.flag(stream,transfers,"hca_moment_flag"))!=0)throw std::runtime_error("Kernel expansion produced a nonfinite moment");
            stats.momentMicroseconds=std::size_t(std::chrono::duration<double,std::micro>(std::chrono::steady_clock::now()-started).count());}
        const std::size_t kernelRank=std::min(options.maxRank,std::size_t(q)/2+1),kernelNeed=2*std::size_t(q)*kernelRank+2*std::size_t(q),batchDoubles=batchScratchBytes/sizeof(double);
        std::size_t spanBytes=0;
        for(std::size_t first=0;first<admissible.size();){
            HostVector<AcaLeafState> states{HostFaultAllocator<AcaLeafState>("aca.states")};std::size_t scratch=0,factorBound=0;
            for(;first<admissible.size();++first){const auto& b=partition.blocks[admissible[first]];
                const auto reach=checkedProduct(b.rows+b.columns,std::min({kernelRank,b.rows,b.columns}));
                if(!states.empty()&&(scratch+kernelNeed>batchDoubles||factorBound+reach>batchDoubles))break;
                const auto u=scratch,v=u+std::size_t(q)*kernelRank,r=v+std::size_t(q)*kernelRank;
                states.push_back(AcaLeafState{b.rowBox,b.columnBox,std::size_t(q),std::size_t(q),u,v,r,r+q,0.,0.,0,acaActive,0,0,int(kernelRank),0});
                scratch+=kernelNeed;factorBound+=reach;}
            const auto firstLeaf=first-states.size();
            const auto kernelBytes=checkedProduct(scratch,sizeof(double))+checkedProduct(states.size(),sizeof(AcaLeafState)+2*sizeof(int));
            flatReserve(kernelBytes,"Kernel expansion batch scratch would exceed the GPU reserve");peakScratch=std::max(peakScratch,flatSum(kernelBytes,expansionBytes));
            Buffer<double> work;work.allocate(scratch);
            const auto kernelRanks=runLeaves(AcaKernelEntry{nodes,deviceBoxes.data},states,options.eps*hcaKernelEps,work,scratch);
            // Moments times the factors of K, for each leaf that converged; the others are stored dense and counted.
            HostVector<HcaProductSpan> spans{HostFaultAllocator<HcaProductSpan>("aca.product-spans")};
            HostVector<AcaLeafState> factorStates{HostFaultAllocator<AcaLeafState>("aca.factor-states")};HostVector<std::size_t> factorBlock{HostFaultAllocator<std::size_t>("aca.factor-block")};
            std::size_t productSize=0,factorScratch=0;
            for(std::size_t s=0;s<states.size();++s){const auto rank=std::size_t(kernelRanks[2*s]);const int status=kernelRanks[2*s+1];const auto index=admissible[firstLeaf+s];
                const auto& b=partition.blocks[index];
                // A leaf whose K did not compress is stored dense, for the reason keepFactors counts.
                if(status!=acaConverged||!rank){std::size_t none=0;HostVector<Kept> unused{HostFaultAllocator<Kept>("aca.batch-kept")};
                    keepFactors(int(rank),status,index,s,none,unused);continue;}
                stats.kernelRankSum+=rank;
                const auto outU=productSize,outV=outU+b.rows*rank;productSize=outV+b.columns*rank;
                spans.push_back({b.rows,rank,boxes[b.rowBox].table,states[s].u,outU});spans.push_back({b.columns,rank,boxes[b.columnBox].table,states[s].v,outV});
                const auto slots=std::min({rank,b.rows,b.columns}),u=factorScratch,v=u+b.rows*slots,r=v+b.columns*slots;
                factorStates.push_back(AcaLeafState{outU,outV,b.rows,b.columns,u,v,r,r+b.columns,0.,0.,0,acaActive,0,0,int(slots),int(rank)});
                factorBlock.push_back(index);factorScratch=r+b.columns+b.rows;
            }
            if(spans.empty()){++stats.batches;continue;}
            const auto factorBytes=flatSum(flatSum(checkedProduct(productSize,sizeof(double)),checkedProduct(factorScratch,sizeof(double))),checkedProduct(checkedProduct(spans.size(),sizeof(HcaProductSpan))+factorStates.size(),sizeof(AcaLeafState)+2*sizeof(int)));
            flatReserve(factorBytes,"Kernel expansion factors would exceed the GPU reserve");peakScratch=std::max(peakScratch,flatSum(flatSum(kernelBytes,factorBytes),expansionBytes));
            Buffer<double> product;product.allocate(productSize);
            {Buffer<HcaProductSpan> deviceSpans;deviceSpans.upload(spans);spanBytes+=checkedProduct(spans.size(),sizeof(HcaProductSpan));
                HexLaunchPace pace(blocks);
                for(std::size_t item=0;item<spans.size();){const auto here=pace.begin(spans.size()-item);
                    hcaProduct<<<unsigned(std::min(blocks,here)),hexThreads,0,stream>>>(deviceSpans.data,q,moments.data,work.data,product.data,item,here);
                    check(cudaGetLastError());check(cudaStreamSynchronize(stream));++launches;pace.end();item+=here;}}
            Buffer<double> factorWork;factorWork.allocate(factorScratch);
            const auto ranks=runLeaves(AcaFactorEntry{product.data},factorStates,options.eps,factorWork,factorScratch);
            std::size_t chunkCount=0;HostVector<Kept> batchKept{HostFaultAllocator<Kept>("aca.batch-kept")};
            // The cap here is the rank of the factors, which reproduce the product exactly when reached.
            for(std::size_t s=0;s<factorStates.size();++s)keepFactors(ranks[2*s],ranks[2*s+1]==acaRankCap?acaConverged:ranks[2*s+1],factorBlock[s],s,chunkCount,batchKept);
            moveFactors(factorStates,batchKept,factorWork,chunkCount,live);
            ++stats.batches;
        }
        if(spanBytes)transfers.push_back({"hca_product_spans","H2D",spanBytes});
    }
    void fillDenseEntries(const Plan& p,double* coefficients,const std::size_t* devicePermutation,cudaStream_t residentStream){
        HostVector<AcaDenseSpan> spans{HostFaultAllocator<AcaDenseSpan>("aca.dense-spans")};std::size_t entries=0;
        for(const auto& s:p.deferredDense){spans.push_back({s.row,s.column,s.rows,s.columns,s.offset,entries});entries+=checkedProduct(s.rows,s.columns);}
        flatReserve(checkedProduct(spans.size(),sizeof(AcaDenseSpan)),"Flat ACA dense-span upload would exceed the GPU reserve");
        Buffer<AcaDenseSpan> deviceSpans;deviceSpans.upload(spans);transfers.push_back({"aca_dense_spans","H2D",checkedProduct(spans.size(),sizeof(AcaDenseSpan))});
        const FlatDevice device=prepared.device();constexpr std::size_t batch=std::size_t(1)<<16;
        for(std::size_t first=0;first<entries;first+=batch){const auto here=std::min(batch,entries-first);
            acaDenseFill<<<unsigned((here+127)/128),128,0,residentStream>>>(device,devicePermutation,sigma.data,deviceSpans.data,spans.size(),first,here,coefficients,prepared.invalid.data);
            check(cudaGetLastError());check(cudaStreamSynchronize(residentStream));++launches;}
    }
    // Host metadata only: the permuted hosts and, per dense leaf, its host ranges. Their product enumerates the work.
    void fillDenseHosts(const Plan& p,double* coefficients,cudaStream_t residentStream){
        HostVector<AcaHostStart> hosts{HostFaultAllocator<AcaHostStart>("aca.dense-hosts")};
        HostVector<std::size_t> hostAt(HostFaultAllocator<std::size_t>("aca.dense-host-index"));permutedHosts(hosts,hostAt);
        HostVector<AcaPairSpan> spans{HostFaultAllocator<AcaPairSpan>("aca.dense-pair-spans")};std::size_t pairs=0;
        for(const auto& s:p.deferredDense){const auto rowHost=hostAt[s.row],columnHost=hostAt[s.column];
            if(hosts[rowHost].start!=s.row||hosts[columnHost].start!=s.column)throw std::invalid_argument("ACA dense leaf does not start at a host");
            const auto rowHosts=hostAt[s.row+s.rows]-rowHost,columnHosts=hostAt[s.column+s.columns]-columnHost;
            spans.push_back({s.row,s.column,s.columns,s.offset,pairs,rowHost,columnHost,columnHosts});pairs=flatSum(pairs,checkedProduct(rowHosts,columnHosts));}
        const std::size_t maxCount=prepared.hexPrepared.order==2?27:8,nOut=std::size_t(prepared.hexPrepared.d.out.n)*prepared.hexPrepared.d.out.n*prepared.hexPrepared.d.out.n;
        const std::size_t blockEntries=maxCount*maxCount,stride=blockEntries*(prepared.hexPrepared.d.imageCount?3:2)+nOut*3*maxCount;
        const std::size_t blocks=processors();
        flatReserve(flatSum(flatSum(checkedProduct(hosts.size(),sizeof(AcaHostStart)),checkedProduct(spans.size(),sizeof(AcaPairSpan))),checkedProduct(checkedProduct(blocks,stride),sizeof(double))),
            "Flat ACA dense host pairs would exceed the GPU reserve");
        Buffer<AcaHostStart> deviceHosts;deviceHosts.upload(hosts);transfers.push_back({"aca_dense_host_starts","H2D",checkedProduct(hosts.size(),sizeof(AcaHostStart))});
        Buffer<AcaPairSpan> deviceSpans;deviceSpans.upload(spans);transfers.push_back({"aca_dense_pair_spans","H2D",checkedProduct(spans.size(),sizeof(AcaPairSpan))});
        Buffer<double> scratch;scratch.allocate(checkedProduct(blocks,stride));
        HexLaunchPace pace(blocks);
        for(std::size_t item=0;item<pairs;){const auto here=pace.begin(pairs-item);
            acaDensePairs<<<unsigned(std::min(blocks,here)),hexThreads,0,residentStream>>>(prepared.hexPrepared.d,deviceHosts.data,deviceSpans.data,spans.size(),sigma.data,item,here,scratch.data,stride,blockEntries,coefficients,prepared.invalid.data);
            check(cudaGetLastError());check(cudaStreamSynchronize(residentStream));++launches;pace.end();item+=here;}
    }
public:
    HostVector<FlatTransfer> transfers{HostFaultAllocator<FlatTransfer>("aca.transfers")};
    std::size_t ownedBytes=0,launches=0;int invalidCode=0;FlatAcaStats stats;
    static constexpr std::size_t batchScratchBytes=std::size_t(256)<<20;
    FlatAcaOracle(const FlatBlock& b,const FlatAcaOptions& o):block(b),options(o){
        if(options.expansionOrder&&block.hostTypes.empty())throw std::invalid_argument("The kernel expansion of the charge Gram belongs to mixed meshes");
        if(options.expansionOrder&&!block.images.empty())throw std::invalid_argument("The kernel expansion of the charge Gram has no image route: it reads the direct term only");
        partition=flatAcaPartition(block,options);stats.elements=partition.elements;stats.clusters=partition.clusters;
        resourceFaultPoint("aca.stream.create");check(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));resourceAcquired();
    }
    FlatAcaOracle(const FlatAcaOracle&)=delete;FlatAcaOracle& operator=(const FlatAcaOracle&)=delete;
    ~FlatAcaOracle(){if(stream){cudaStreamSynchronize(stream);resourceReleased(cudaStreamDestroy(stream)==cudaSuccess);}}
    const HostVector<std::size_t>& permutation()const{return partition.permutation;}
    void compress(){
        if(compressed)throw std::runtime_error("ACA oracle already compressed");
        // Prepared charges, sigma, the permutation and the active counter stay live through compression.
        const auto n=block.charges();std::size_t live=flatSum(FlatPrepared::workspaceBytes(block),checkedProduct(n,2*sizeof(double))+sizeof(int));
        flatReserve(live,"Flat ACA GPU reserve would be exceeded");
        prepared.prepare(block,stream,transfers,launches,invalidCode,false);
        sigma.allocate(n);const FlatDevice device=prepared.device();
        acaSigma<<<unsigned((n+127)/128),128,0,stream>>>(device,sigma.data,prepared.invalid.data);check(cudaGetLastError());++launches;
        if((invalidCode=prepared.flag(stream,transfers,"energy_flag"))!=0)throw std::runtime_error("Flat device Gram produced a nonfinite or nonpositive self-energy");
        rankOf.assign(partition.blocks.size(),0);
        HostVector<std::size_t> admissible{HostFaultAllocator<std::size_t>("aca.admissible")};
        for(std::size_t b=0;b<partition.blocks.size();++b)if(partition.blocks[b].admissible)admissible.push_back(b);
        stats.admissibleBlocks=admissible.size();
        const auto maxRank=options.maxRank;std::size_t peakScratch=0;const auto acaStarted=std::chrono::steady_clock::now();
        active.allocate(1);
        if(options.expansionOrder){if(!admissible.empty())compressExpansion(admissible,live,peakScratch);}
        else{
        Buffer<std::size_t> perm;perm.upload(partition.permutation);
        transfers.push_back({"aca_permutation","H2D",checkedProduct(n,sizeof(std::size_t))});
        for(std::size_t first=0;first<admissible.size();){
            // Bounded batch: U, V, residual row and used flags per leaf.
            HostVector<AcaLeafState> states{HostFaultAllocator<AcaLeafState>("aca.states")};std::size_t scratch=0;
            for(;first<admissible.size();++first){const auto& b=partition.blocks[admissible[first]];
                // The rank of a leaf cannot pass its smaller side.
                const auto leafRank=std::min({maxRank,b.rows,b.columns});const auto need=checkedProduct(b.rows+b.columns,leafRank)+b.columns+b.rows;
                if(!states.empty()&&checkedProduct(scratch+need,sizeof(double))>batchScratchBytes)break;
                AcaLeafState s{b.row,b.column,b.rows,b.columns,scratch,scratch+b.rows*leafRank,scratch+(b.rows+b.columns)*leafRank,
                    scratch+(b.rows+b.columns)*leafRank+b.columns,0.,0.,0,acaActive,0,0,int(maxRank),0};
                states.push_back(s);scratch+=need;}
            const auto batchBytes=checkedProduct(scratch,sizeof(double))+checkedProduct(states.size(),sizeof(AcaLeafState)+2*sizeof(int));
            flatReserve(batchBytes,"Flat ACA batch scratch would exceed the GPU reserve");peakScratch=std::max(peakScratch,batchBytes);
            Buffer<double> work;work.allocate(scratch);
            const auto ranks=runLeaves(AcaGramEntry{device,perm.data,sigma.data},states,options.eps,work,scratch);
            const auto firstLeaf=first-states.size();std::size_t chunkCount=0;
            HostVector<Kept> batchKept{HostFaultAllocator<Kept>("aca.batch-kept")};
            for(std::size_t s=0;s<states.size();++s)keepFactors(ranks[2*s],ranks[2*s+1],admissible[firstLeaf+s],s,chunkCount,batchKept);
            moveFactors(states,batchKept,work,chunkCount,live);
            ++stats.batches;
        }
        }
        stats.acaMicroseconds=std::size_t(std::chrono::duration<double,std::micro>(std::chrono::steady_clock::now()-acaStarted).count())-stats.momentMicroseconds;
        if(prepared.hex){const auto& h=prepared.hexPrepared;
            stats.farPairs=h.classItems[1];stats.midPairs=h.classItems[2]+h.classItems[3];stats.gradedPairs=h.classItems[4]+h.classItems[5];
            stats.farMicroseconds=h.classMicroseconds[1];stats.midMicroseconds=h.classMicroseconds[2]+h.classMicroseconds[3];stats.gradedMicroseconds=h.classMicroseconds[4]+h.classMicroseconds[5];}
        if(statesBytes)transfers.push_back({"aca_leaf_spans","H2D",statesBytes});
        if(reportBytes)transfers.push_back({"aca_leaf_status","D2H",reportBytes});
        if(activeCount)transfers.push_back({"aca_active_count","D2H",activeCount});
        for(const auto& k:kept){stats.rankSum+=k.rank;stats.maxRank=std::max(stats.maxRank,k.rank);
            const auto& b=partition.blocks[k.block];stats.factorCoefficients+=k.rank*(b.rows+b.columns);}
        stats.lowrankLeaves=kept.size();stats.denseLeaves=partition.blocks.size()-kept.size();
        for(std::size_t b=0;b<partition.blocks.size();++b)if(!rankOf[b])stats.denseCoefficients+=partition.blocks[b].rows*partition.blocks[b].columns;
        ownedBytes=live+peakScratch;activeBytes=live;compressed=true;
    }
    template<class LeafType>HostVector<LeafType> leaves()const{
        if(!compressed)throw std::runtime_error("ACA oracle is not compressed");
        HostVector<LeafType> result;result.reserve(partition.blocks.size());
        for(std::size_t b=0;b<partition.blocks.size();++b){const auto& s=partition.blocks[b];
            LeafType leaf{};leaf.row=s.row;leaf.column=s.column;leaf.rows=s.rows;leaf.columns=s.columns;leaf.rank=rankOf[b];leaf.deferred=true;
            result.push_back(std::move(leaf));}
        return result;
    }
    void fillDense(const Plan& p,double* coefficients,const std::size_t* devicePermutation,double* scale,cudaStream_t residentStream){
        if(!compressed)throw std::runtime_error("ACA oracle is not compressed");
        const auto n=block.charges();
        if(p.n!=n||p.permutation.size()!=n||p.deferredLowrank.size()!=kept.size()||p.deferredDense.size()+kept.size()!=partition.blocks.size()||!p.uploads.empty())
            throw std::invalid_argument("ACA plan does not match its compressed partition");
        for(std::size_t i=0;i<n;++i)if(p.permutation[i]!=partition.permutation[i])throw std::invalid_argument("ACA plan permutation mismatch");
        // Plan spans follow leaf order; kept leaves are in admissible-block order.
        for(std::size_t k=0;k<kept.size();++k){const auto& span=p.deferredLowrank[k];const auto& b=partition.blocks[kept[k].block];
            if(span.row!=b.row||span.column!=b.column||span.rows!=b.rows||span.columns!=b.columns||span.rank!=kept[k].rank)
                throw std::invalid_argument("ACA low-rank span mismatch");
            const auto* chunk=chunks[kept[k].chunk]->data;
            // HACApK layout: a1 = V (columns x rank), a2 = U (rows x rank), column-major.
            check(cudaMemcpyAsync(coefficients+span.a,chunk+kept[k].v,span.columns*span.rank*sizeof(double),cudaMemcpyDeviceToDevice,residentStream));
            check(cudaMemcpyAsync(coefficients+span.b,chunk+kept[k].u,span.rows*span.rank*sizeof(double),cudaMemcpyDeviceToDevice,residentStream));
        }
        const auto denseStarted=std::chrono::steady_clock::now();
        if(!p.deferredDense.empty()){
            if(prepared.hex)fillDenseHosts(p,coefficients,residentStream);else fillDenseEntries(p,coefficients,devicePermutation,residentStream);
            if((invalidCode=prepared.flag(residentStream,transfers,"dense_flag"))!=0)throw std::runtime_error("Flat ACA dense leaf produced a nonfinite entry");
        }
        check(cudaMemcpyAsync(scale,sigma.data,n*sizeof(double),cudaMemcpyDeviceToDevice,residentStream));
        check(cudaStreamSynchronize(residentStream));
        stats.denseMicroseconds=std::size_t(std::chrono::duration<double,std::micro>(std::chrono::steady_clock::now()-denseStarted).count());
    }
};
} // namespace radia_hmatrix_cuda
