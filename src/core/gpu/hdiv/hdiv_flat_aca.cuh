// Flat-TET charge Gram compressed on the device by partial-pivot adaptive
// cross approximation. The host builds only geometric metadata: a cluster
// tree over the charge host elements (centre plus support radius, longest-axis
// median bisection) and HACApK-style box-gap admissibility on support boxes.
// It is an independent partition, not the CPU HACApK leaf set. Every Gram
// entry, sigma, dense leaf and low-rank factor is generated and kept on the
// device; only per-leaf integer ranks/status and step counts cross D2H.
#pragma once
#include <functional>
#include <numeric>
#include "hdiv_flat_gram.cuh"
namespace radia_hmatrix_cuda {
// rankCapDense: a leaf that reaches maxRank unconverged is stored as an exact
// dense leaf of the same operator (counted); otherwise construction fails.
struct FlatAcaOptions {
    double eps=0,eta=0;std::size_t leafSize=0,maxRank=0;bool rankCapDense=false;
    void validate()const{
        if(!std::isfinite(eps)||!(eps>0)||!(eps<1))throw std::invalid_argument("ACA eps must be finite in (0,1)");
        if(!std::isfinite(eta)||!(eta>0))throw std::invalid_argument("ACA admissibility eta must be finite and positive");
        if(!leafSize||leafSize>4096)throw std::invalid_argument("ACA leaf size must be in [1,4096] elements");
        if(!maxRank||maxRank>1024)throw std::invalid_argument("ACA max rank must be in [1,1024]");
    }
};
struct FlatAcaBlock {std::size_t row,column,rows,columns;bool admissible;};
struct FlatAcaPartition {
    HostVector<std::size_t> permutation{HostFaultAllocator<std::size_t>("aca.permutation")};
    HostVector<FlatAcaBlock> blocks{HostFaultAllocator<FlatAcaBlock>("aca.blocks")};
    std::size_t elements=0,clusters=0;
};
// Geometry-only partition. Charges of one host element stay contiguous and in
// ascending input order; equal centres or zero extent end the split as a leaf.
inline FlatAcaPartition flatAcaPartition(const FlatBlock& block,const FlatAcaOptions& options){
    block.validate();options.validate();
    const auto n=block.charges();const std::size_t cellCount=block.cells.size()/12,faceCount=block.faces.size()/9;
    constexpr std::size_t none=std::numeric_limits<std::size_t>::max();
    HostVector<std::size_t> cellElement(cellCount,none,HostFaultAllocator<std::size_t>("aca.cell-element")),
        faceElement(faceCount,none,HostFaultAllocator<std::size_t>("aca.face-element"));
    for(std::size_t a=0;a<n;++a)(block.kind[a]==0?cellElement:faceElement)[std::size_t(block.host[a])]=0;
    HostVector<double> centre{HostFaultAllocator<double>("aca.centre")},radius{HostFaultAllocator<double>("aca.radius")};
    std::size_t elements=0;
    auto addElement=[&](const double* V,int vertices){
        double c[3]={0,0,0};for(int i=0;i<vertices;++i)for(int k=0;k<3;++k)c[k]+=V[3*i+k];
        for(int k=0;k<3;++k)c[k]/=vertices;double r=0;
        for(int i=0;i<vertices;++i){const double dx=V[3*i]-c[0],dy=V[3*i+1]-c[1],dz=V[3*i+2]-c[2];r=std::max(r,std::sqrt(dx*dx+dy*dy+dz*dz));}
        centre.insert(centre.end(),c,c+3);radius.push_back(r);return elements++;};
    for(std::size_t h=0;h<cellCount;++h)if(cellElement[h]!=none)cellElement[h]=addElement(block.cells.data()+12*h,4);
    for(std::size_t h=0;h<faceCount;++h)if(faceElement[h]!=none)faceElement[h]=addElement(block.faces.data()+9*h,3);
    HostVector<std::size_t> offsets(elements+1,0,HostFaultAllocator<std::size_t>("aca.charge-offsets")),members(n,0,HostFaultAllocator<std::size_t>("aca.charge-members"));
    auto elementOf=[&](std::size_t a){return (block.kind[a]==0?cellElement:faceElement)[std::size_t(block.host[a])];};
    for(std::size_t a=0;a<n;++a)++offsets[elementOf(a)+1];
    for(std::size_t e=0;e<elements;++e)offsets[e+1]+=offsets[e];
    {auto next=offsets;for(std::size_t a=0;a<n;++a)members[next[elementOf(a)]++]=a;}
    struct Node {std::size_t first,count,charge=0,charges=0,left=0,right=0;double lo[3],hi[3],width=0;bool leaf=true;};
    HostVector<Node> nodes{HostFaultAllocator<Node>("aca.nodes")};
    HostVector<std::size_t> order(elements,0,HostFaultAllocator<std::size_t>("aca.order"));std::iota(order.begin(),order.end(),std::size_t(0));
    std::function<std::size_t(std::size_t,std::size_t)> build=[&](std::size_t first,std::size_t count)->std::size_t{
        Node node{first,count};double clo[3],chi[3];
        for(int k=0;k<3;++k){node.lo[k]=clo[k]=std::numeric_limits<double>::infinity();node.hi[k]=chi[k]=-std::numeric_limits<double>::infinity();}
        for(std::size_t i=first;i<first+count;++i){const auto e=order[i];
            for(int k=0;k<3;++k){const double c=centre[3*e+k];clo[k]=std::min(clo[k],c);chi[k]=std::max(chi[k],c);
                node.lo[k]=std::min(node.lo[k],c-radius[e]);node.hi[k]=std::max(node.hi[k],c+radius[e]);}}
        double w=0;for(int k=0;k<3;++k)w+=(node.hi[k]-node.lo[k])*(node.hi[k]-node.lo[k]);node.width=std::sqrt(w);
        int axis=0;for(int k=1;k<3;++k)if(chi[k]-clo[k]>chi[axis]-clo[axis])axis=k;
        const auto index=nodes.size();nodes.push_back(node);
        if(count<=options.leafSize||!(chi[axis]>clo[axis]))return index;
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
    auto admissible=[&](const Node& s,const Node& t){
        double gap=0;for(int k=0;k<3;++k){const double d=std::max({0.0,t.lo[k]-s.hi[k],s.lo[k]-t.hi[k]});gap+=d*d;}
        gap=std::sqrt(gap);return gap>0&&std::min(s.width,t.width)<=options.eta*gap;};
    // Upper/diagonal blocks only: s precedes t in charge order whenever s != t.
    std::function<void(std::size_t,std::size_t)> split=[&](std::size_t si,std::size_t ti){
        const Node s=nodes[si],t=nodes[ti];
        if(si==ti){
            if(s.leaf){result.blocks.push_back({s.charge,s.charge,s.charges,s.charges,false});return;}
            split(s.left,s.left);split(s.left,s.right);split(s.right,s.right);return;}
        if(admissible(s,t)){result.blocks.push_back({s.charge,t.charge,s.charges,t.charges,true});return;}
        if(s.leaf&&t.leaf){result.blocks.push_back({s.charge,t.charge,s.charges,t.charges,false});return;}
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
    double frobenius2,maxAbs;int rank,status,pivotRow,steps;
};
enum AcaStatus:int {acaActive=0,acaConverged=1,acaRankCap=2,acaNonfinite=3};
constexpr int acaThreads=128;
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
// One ACA step for every active leaf of the batch (one block per leaf). A
// negligible pivot row is marked used and the next unused row is tried; the
// leaf converges by row exhaustion only after every row was tried.
__global__ void acaStep(FlatDevice d,const std::size_t* perm,const double* sigma,AcaLeafState* leaves,double* scratch,
                        int maxRank,double eps,int* active){
    __shared__ double sv[acaThreads];__shared__ std::size_t si[acaThreads];
    AcaLeafState L=leaves[blockIdx.x];if(L.status!=acaActive)return;
    const std::size_t i=std::size_t(L.pivotRow),k=std::size_t(L.rank),rows=L.rows,columns=L.columns;
    double* U=scratch+L.u;double* V=scratch+L.v;double* R=scratch+L.residual;double* used=scratch+L.used;
    double best=-1,rowAbs=0;std::size_t jstar=acaNone;int bad=0;
    for(auto j=std::size_t(threadIdx.x);j<columns;j+=acaThreads){
        const double e=acaEntry(d,perm,sigma,L.row+i,L.column+j);if(!isfinite(e))bad=1;rowAbs=fmax(rowAbs,fabs(e));
        double r=e;for(std::size_t l=0;l<k;++l)r-=U[i+rows*l]*V[j+columns*l];
        R[j]=r;if(fabs(r)>best){best=fabs(r);jstar=j;}
    }
    acaBlockArgmax(best,jstar,sv,si);
    {double flag=bad;std::size_t ignored=0;acaBlockArgmax(flag,ignored,sv,si);bad=flag>0;}
    {std::size_t ignored=0;acaBlockArgmax(rowAbs,ignored,sv,si);}
    const double maxAbs=fmax(L.maxAbs,rowAbs);
    if(bad||!isfinite(best)){if(threadIdx.x==0){L.status=acaNonfinite;leaves[blockIdx.x]=L;}return;}
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
        const double e=acaEntry(d,perm,sigma,L.row+r,L.column+jstar);if(!isfinite(e))bad=1;columnAbs=fmax(columnAbs,fabs(e));
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
    if(threadIdx.x==0){
        const double term=nu2*nv2,frobenius2=L.frobenius2+2*cross+term;
        used[i]=1.0;L.maxAbs=fmax(maxAbs,columnAbs);L.rank=int(k+1);++L.steps;L.frobenius2=frobenius2;
        if(bad||!isfinite(term)||!isfinite(frobenius2)||!(frobenius2>0))L.status=acaNonfinite;
        else if(sqrt(term)<=eps*sqrt(frobenius2)||std::size_t(L.rank)==(rows<columns?rows:columns))L.status=acaConverged;
        else if(L.rank>=maxRank)L.status=acaRankCap;
        else if(nextRow==acaNone)L.status=acaConverged;
        else{L.pivotRow=int(nextRow);atomicAdd(active,1);}
        leaves[blockIdx.x]=L;
    }
}
__global__ void acaSigma(FlatDevice d,double* sigma,int* invalid){
    const auto a=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(a>=d.n)return;
    const double diagonal=flatEntry(d,a,a);if(!isfinite(diagonal)||!(diagonal>0)){atomicExch(invalid,1);return;}
    sigma[a]=sqrt(diagonal);
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
struct FlatAcaStats {
    std::size_t elements=0,clusters=0,denseLeaves=0,lowrankLeaves=0,admissibleBlocks=0,rankCapPromotions=0,storagePromotions=0,zeroRankLeaves=0,
        rankSum=0,maxRank=0,steps=0,factorCoefficients=0,denseCoefficients=0,batches=0;
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
public:
    HostVector<FlatTransfer> transfers{HostFaultAllocator<FlatTransfer>("aca.transfers")};
    std::size_t ownedBytes=0,launches=0;int invalidCode=0;FlatAcaStats stats;
    static constexpr std::size_t batchScratchBytes=std::size_t(256)<<20;
    FlatAcaOracle(const FlatBlock& b,const FlatAcaOptions& o):block(b),options(o){
        partition=flatAcaPartition(block,options);stats.elements=partition.elements;stats.clusters=partition.clusters;
        resourceFaultPoint("aca.stream.create");check(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));resourceAcquired();
    }
    FlatAcaOracle(const FlatAcaOracle&)=delete;FlatAcaOracle& operator=(const FlatAcaOracle&)=delete;
    ~FlatAcaOracle(){if(stream){cudaStreamSynchronize(stream);resourceReleased(cudaStreamDestroy(stream)==cudaSuccess);}}
    const HostVector<std::size_t>& permutation()const{return partition.permutation;}
    void compress(){
        if(compressed)throw std::runtime_error("ACA oracle already compressed");
        const auto n=block.charges();std::size_t live=FlatPrepared::workspaceBytes(block)+checkedProduct(n,sizeof(double));
        flatReserve(live,"Flat ACA GPU reserve would be exceeded");
        prepared.prepare(block,stream,transfers,launches,invalidCode);
        sigma.allocate(n);const FlatDevice device=prepared.device();
        acaSigma<<<unsigned((n+127)/128),128,0,stream>>>(device,sigma.data,prepared.invalid.data);check(cudaGetLastError());++launches;
        if((invalidCode=prepared.flag(stream,transfers,"energy_flag"))!=0)throw std::runtime_error("Flat device Gram produced a nonfinite or nonpositive self-energy");
        Buffer<std::size_t> perm;perm.upload(partition.permutation);
        transfers.push_back({"aca_permutation","H2D",checkedProduct(n,sizeof(std::size_t))});
        live+=checkedProduct(n,sizeof(std::size_t));
        rankOf.assign(partition.blocks.size(),0);
        HostVector<std::size_t> admissible{HostFaultAllocator<std::size_t>("aca.admissible")};
        for(std::size_t b=0;b<partition.blocks.size();++b)if(partition.blocks[b].admissible)admissible.push_back(b);
        stats.admissibleBlocks=admissible.size();
        const auto maxRank=options.maxRank;std::size_t statesBytes=0,reportBytes=0,activeCount=0,peakScratch=0;
        Buffer<int> active;active.allocate(1);
        for(std::size_t first=0;first<admissible.size();){
            // Bounded batch: U, V, residual row and used flags per leaf.
            HostVector<AcaLeafState> states{HostFaultAllocator<AcaLeafState>("aca.states")};std::size_t scratch=0;
            for(;first<admissible.size();++first){const auto& b=partition.blocks[admissible[first]];
                const auto need=checkedProduct(b.rows+b.columns,maxRank)+b.columns+b.rows;
                if(!states.empty()&&checkedProduct(scratch+need,sizeof(double))>batchScratchBytes)break;
                AcaLeafState s{b.row,b.column,b.rows,b.columns,scratch,scratch+b.rows*maxRank,scratch+(b.rows+b.columns)*maxRank,
                    scratch+(b.rows+b.columns)*maxRank+b.columns,0.,0.,0,acaActive,0,0};
                states.push_back(s);scratch+=need;}
            const auto batchBytes=checkedProduct(scratch,sizeof(double))+checkedProduct(states.size(),sizeof(AcaLeafState)+2*sizeof(int));
            flatReserve(batchBytes,"Flat ACA batch scratch would exceed the GPU reserve");peakScratch=std::max(peakScratch,batchBytes);
            Buffer<double> work;work.allocate(scratch);check(cudaMemsetAsync(work.data,0,scratch*sizeof(double),stream));
            Buffer<AcaLeafState> deviceStates;deviceStates.upload(states);statesBytes+=checkedProduct(states.size(),sizeof(AcaLeafState));
            int running=int(states.size());
            // A step adds rank or retires one row, so rows + maxRank steps bound every leaf.
            std::size_t bound=0;for(const auto& s:states)bound=std::max(bound,s.rows+maxRank);
            for(std::size_t step=0;running&&step<bound;++step){
                check(cudaMemsetAsync(active.data,0,sizeof(int),stream));
                acaStep<<<unsigned(states.size()),acaThreads,0,stream>>>(device,perm.data,sigma.data,deviceStates.data,work.data,int(maxRank),options.eps,active.data);
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
            const auto firstLeaf=first-states.size();std::size_t chunkCount=0;
            HostVector<Kept> batchKept{HostFaultAllocator<Kept>("aca.batch-kept")};
            for(std::size_t s=0;s<states.size();++s){const auto rank=std::size_t(ranks[2*s]);const int status=ranks[2*s+1];
                const auto& b=partition.blocks[admissible[firstLeaf+s]];
                if(status==acaNonfinite)throw std::runtime_error("ACA produced a nonfinite pivot, entry or norm");
                if(status==acaRankCap){if(!options.rankCapDense)throw std::runtime_error("ACA leaf reached max_rank without converging");++stats.rankCapPromotions;continue;}
                if(status!=acaConverged||rank>std::min(b.rows,b.columns))throw std::runtime_error("ACA leaf status or rank invalid");
                // Rank 0 means every row was negligible: the leaf is stored dense, never dropped.
                if(!rank){++stats.zeroRankLeaves;continue;}
                if(checkedProduct(rank,b.rows+b.columns)>=checkedProduct(b.rows,b.columns)){++stats.storagePromotions;continue;}
                batchKept.push_back({admissible[firstLeaf+s],rank,chunks.size(),chunkCount,chunkCount+b.rows*rank,s});
                chunkCount+=rank*(b.rows+b.columns);rankOf[admissible[firstLeaf+s]]=rank;
            }
            if(chunkCount){
                flatReserve(checkedProduct(chunkCount,sizeof(double)),"Flat ACA factor storage would exceed the GPU reserve");
                auto chunk=makeHostUnique<Buffer<double>>("aca.chunk");chunk->allocate(chunkCount);
                // U (rows x rank) and V (columns x rank) are leading slices of the scratch columns.
                for(const auto& k:batchKept){const auto& st=states[k.state];const auto& b=partition.blocks[k.block];
                    check(cudaMemcpyAsync(chunk->data+k.u,work.data+st.u,b.rows*k.rank*sizeof(double),cudaMemcpyDeviceToDevice,stream));
                    check(cudaMemcpyAsync(chunk->data+k.v,work.data+st.v,b.columns*k.rank*sizeof(double),cudaMemcpyDeviceToDevice,stream));}
                check(cudaStreamSynchronize(stream));live+=checkedProduct(chunkCount,sizeof(double));chunks.push_back(std::move(chunk));
                kept.insert(kept.end(),batchKept.begin(),batchKept.end());
            }
            ++stats.batches;
        }
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
        HostVector<AcaDenseSpan> spans{HostFaultAllocator<AcaDenseSpan>("aca.dense-spans")};std::size_t entries=0;
        for(const auto& s:p.deferredDense){spans.push_back({s.row,s.column,s.rows,s.columns,s.offset,entries});entries+=checkedProduct(s.rows,s.columns);}
        if(!spans.empty()){
            Buffer<AcaDenseSpan> deviceSpans;deviceSpans.upload(spans);transfers.push_back({"aca_dense_spans","H2D",checkedProduct(spans.size(),sizeof(AcaDenseSpan))});
            const FlatDevice device=prepared.device();constexpr std::size_t batch=std::size_t(1)<<16;
            for(std::size_t first=0;first<entries;first+=batch){const auto here=std::min(batch,entries-first);
                acaDenseFill<<<unsigned((here+127)/128),128,0,residentStream>>>(device,devicePermutation,sigma.data,deviceSpans.data,spans.size(),first,here,coefficients,prepared.invalid.data);
                check(cudaGetLastError());check(cudaStreamSynchronize(residentStream));++launches;}
            if((invalidCode=prepared.flag(residentStream,transfers,"dense_flag"))!=0)throw std::runtime_error("Flat ACA dense leaf produced a nonfinite entry");
        }
        check(cudaMemcpyAsync(scale,sigma.data,n*sizeof(double),cudaMemcpyDeviceToDevice,residentStream));
        check(cudaStreamSynchronize(residentStream));
    }
};
} // namespace radia_hmatrix_cuda
