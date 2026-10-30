// Resident symmetric HACApK leaf application; no FE or assembly substitution.
#pragma once
#include <cuda_runtime.h>
#include <vector>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <utility>
#include <memory>
#include "hmatrix_host_fault.hpp"
#include "hmatrix_leaf_coverage.hpp"

namespace radia_hmatrix_cuda {
inline void check(cudaError_t status) {
    if(status!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(status));
}
#include "hdiv_assembly_profile.cuh"
inline std::size_t checkedProduct(std::size_t a,std::size_t b) {
    if(b && a>std::numeric_limits<std::size_t>::max()/b)
        throw std::overflow_error("H-matrix allocation overflow");
    return a*b;
}
struct Direction {
    std::size_t input,output,inputs,outputs,rank,a,b,tmp,partial;
    bool transpose;
};
struct Leaf {
    std::size_t row,column,rows,columns,rank;
    // rank>0: a1[columns,rank], a2[rows,rank], column-major HACApK layout.
    // rank==0: a1[rows,columns], row-major dense layout.
    HostVector<double> a1{HostFaultAllocator<double>("leaf.a1")},a2{HostFaultAllocator<double>("leaf.a2")};
    bool deferred=false;
};
// Borrow only during synchronous host packing; Plan retains its own values.
struct BorrowedValues {
    const double* data=nullptr;std::size_t count=0;
    BorrowedValues()=default;
    BorrowedValues(const double* p,std::size_t n):data(p),count(n){if(n&&!p)throw std::invalid_argument("Missing borrowed factor values");}
    std::size_t size()const{return count;}
    bool empty()const{return !count;}
    const double* begin()const{return data;}
    const double* end()const{return count?data+count:data;}
};
struct BorrowedLeaf {
    std::size_t row=0,column=0,rows=0,columns=0,rank=0;
    BorrowedValues a1,a2;bool deferred=false;
};
struct DeferredDenseSpan {std::size_t row,column,rows,columns,offset;};
struct CoefficientUpload {std::size_t destination,source,count;};
struct Plan {
    std::size_t n=0,tmpCount=0,partialCount=0;
    std::size_t coefficientCount=0;
    HostVector<Direction> directions{HostFaultAllocator<Direction>("plan.directions")};
    HostVector<double> coefficients{HostFaultAllocator<double>("plan.coefficients")},scale{HostFaultAllocator<double>("plan.scale")};
    HostVector<std::size_t> permutation{HostFaultAllocator<std::size_t>("plan.permutation")},rowOffsets{HostFaultAllocator<std::size_t>("plan.row-offsets")},partialIndices{HostFaultAllocator<std::size_t>("plan.partial-indices")};
    HostVector<DeferredDenseSpan> deferredDense{HostFaultAllocator<DeferredDenseSpan>("plan.deferred-dense")};
    HostVector<CoefficientUpload> uploads{HostFaultAllocator<CoefficientUpload>("plan.uploads")};
};
template<class LeafList,class IndexList,class ScaleList>
inline Plan packLeaves(std::size_t n,const LeafList& leaves,
                 const IndexList& permutation,const ScaleList& scale) {
    if(!n || permutation.size()!=n || scale.size()!=n)
        throw std::invalid_argument("Invalid H-matrix dimension or scaling");
    Plan p;p.n=n;p.permutation.assign(permutation.begin(),permutation.end());p.scale.assign(scale.begin(),scale.end());
    HostVector<bool> seen(n,false,HostFaultAllocator<bool>("pack.seen"));
    for(auto i:permutation){if(i>=n||seen[i])throw std::invalid_argument("Invalid permutation");seen[i]=true;}
    for(double s:scale)if(!std::isfinite(s)||s<=0)throw std::invalid_argument("Invalid physical scale");
    HostVector<HostVector<std::size_t>> contributions(n,HostVector<std::size_t>{HostFaultAllocator<std::size_t>("pack.row-contribution")},HostFaultAllocator<HostVector<std::size_t>>("pack.contributions"));
    std::size_t copyCount=0;
    for(const auto& leaf:leaves) {
        if(!leaf.rows||!leaf.columns||leaf.row>=n||leaf.column>=n||
           leaf.rows>n-leaf.row||leaf.columns>n-leaf.column)
            throw std::invalid_argument("Leaf range outside matrix");
        if(leaf.row>leaf.column)throw std::invalid_argument("Pack upper/diagonal leaves only");
        if(leaf.row==leaf.column&&leaf.rows!=leaf.columns)
            throw std::invalid_argument("Diagonal leaf must be square");
        if(leaf.row<leaf.column&&leaf.rows>leaf.column-leaf.row)
            throw std::invalid_argument("Upper leaf overlaps diagonal");
        if(leaf.deferred&&(leaf.rank||!leaf.a1.empty()||!leaf.a2.empty()))throw std::invalid_argument("Invalid deferred dense storage");
        if(!leaf.deferred&&(leaf.a1.size()!=checkedProduct(leaf.columns,leaf.rank?leaf.rank:leaf.rows)||
           leaf.a2.size()!=(leaf.rank?checkedProduct(leaf.rows,leaf.rank):0))
           )
            throw std::invalid_argument("Invalid factor storage");
        for(double v:leaf.a1)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite leaf");
        for(double v:leaf.a2)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite leaf");
        if(leaf.a1.size()>std::numeric_limits<std::size_t>::max()-copyCount||
           leaf.a2.size()>std::numeric_limits<std::size_t>::max()-copyCount-leaf.a1.size())throw std::overflow_error("Host coefficient copy count overflow");
        copyCount+=leaf.a1.size()+leaf.a2.size();
    }
    // Only retained host factors/nondeferred dense values, never deferred G.
    // Reserving once avoids repeated coefficient-buffer growth during pack.
    p.coefficients.reserve(copyCount);
    for(const auto& leaf:leaves) {
        const auto firstCount=checkedProduct(leaf.columns,leaf.rank?leaf.rank:leaf.rows),secondCount=leaf.rank?checkedProduct(leaf.rows,leaf.rank):0;
        if(firstCount>std::numeric_limits<std::size_t>::max()-secondCount||firstCount+secondCount>std::numeric_limits<std::size_t>::max()-p.coefficientCount)throw std::overflow_error("Coefficient offset overflow");
        const auto a=p.coefficientCount,b=a+firstCount;p.coefficientCount+=firstCount+secondCount;
        if(leaf.deferred)p.deferredDense.push_back({leaf.row,leaf.column,leaf.rows,leaf.columns,a});
        else{
            p.uploads.push_back({a,p.coefficients.size(),firstCount+secondCount});
            if(!leaf.a1.empty())p.coefficients.insert(p.coefficients.end(),leaf.a1.begin(),leaf.a1.end());
            if(!leaf.a2.empty())p.coefficients.insert(p.coefficients.end(),leaf.a2.begin(),leaf.a2.end());
        }
        for(int transpose=0;transpose<(leaf.row==leaf.column?1:2);++transpose) {
            Direction d{transpose?leaf.row:leaf.column,transpose?leaf.column:leaf.row,
                transpose?leaf.rows:leaf.columns,transpose?leaf.columns:leaf.rows,
                leaf.rank,a,b,p.tmpCount,p.partialCount,transpose!=0};
            if(d.rank>std::numeric_limits<std::size_t>::max()-p.tmpCount||
               d.outputs>std::numeric_limits<std::size_t>::max()-p.partialCount)
                throw std::overflow_error("H-matrix workspace overflow");
            p.tmpCount+=d.rank;p.partialCount+=d.outputs;
            for(std::size_t i=0;i<d.outputs;++i)contributions[d.output+i].push_back(d.partial+i);
            p.directions.push_back(d);
        }
    }
    if(p.coefficients.size()!=copyCount)throw std::runtime_error("Host coefficient pack count mismatch");
    if(p.directions.empty())throw std::invalid_argument("Empty H-matrix");
    p.rowOffsets.push_back(0);
    for(const auto& row:contributions){p.partialIndices.insert(p.partialIndices.end(),row.begin(),row.end());p.rowOffsets.push_back(p.partialIndices.size());}
    return p;
}
template<class PermAllocator,class ScaleAllocator>
inline Plan pack(std::size_t n,const HostVector<Leaf>& leaves,
                 const std::vector<std::size_t,PermAllocator>& permutation,const std::vector<double,ScaleAllocator>& scale){
    return packLeaves(n,leaves,permutation,scale);
}
// Higher-order HDiv charge maps can amplify low-order dot-product roundoff.
// Keep the physical operator and stable leaf order; compensate scalar reductions.
// Stock builds use --fmad=false, so these correction terms are not contracted.
struct CompensatedDeviceSum {
    double value=0, correction=0;
    __device__ void add(double x) {
        const double t=value+x;
        correction += fabs(value)>=fabs(x) ? (value-t)+x : (x-t)+value;
        value=t;
    }
    __device__ void addProduct(double a,double b) {
        const double product=a*b;
        add(product);add(fma(a,b,-product));
    }
    __device__ double get() const {return value+correction;}
};
__global__ void gatherInput(std::size_t n,const std::size_t* permutation,const double* scale,
                            const double* state,double* input) {
    auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(i<n)input[i]=scale[permutation[i]]*state[permutation[i]];
}
__global__ void projectRanks(const Direction* directions,std::size_t count,
                            const double* coefficients,const double* input,double* tmp) {
    const auto task=std::size_t(blockIdx.x);if(task>=count)return;
    const auto d=directions[task];
    for(auto k=std::size_t(threadIdx.x);k<d.rank;k+=blockDim.x) {
        const auto base=d.transpose?d.b:d.a;CompensatedDeviceSum value;
        for(std::size_t j=0;j<d.inputs;++j)value.addProduct(coefficients[base+j+d.inputs*k],input[d.input+j]);
        tmp[d.tmp+k]=value.get();
    }
}
__global__ void applyLeaves(const Direction* directions,std::size_t count,
                           const double* coefficients,const double* input,const double* tmp,double* partial) {
    const auto task=std::size_t(blockIdx.x);if(task>=count)return;
    const auto d=directions[task];
    for(auto i=std::size_t(threadIdx.x);i<d.outputs;i+=blockDim.x) {
        CompensatedDeviceSum value;
        if(d.rank) {
            const auto base=d.transpose?d.a:d.b;
            for(std::size_t k=0;k<d.rank;++k)value.addProduct(coefficients[base+i+d.outputs*k],tmp[d.tmp+k]);
        }else {
            for(std::size_t j=0;j<d.inputs;++j) {
                const auto offset=d.transpose?j*d.outputs+i:i*d.inputs+j;
                value.addProduct(coefficients[d.a+offset],input[d.input+j]);
            }
        }
        partial[d.partial+i]=value.get();
    }
}
__global__ void reduceOutput(std::size_t n,const std::size_t* offsets,const std::size_t* indices,
        const std::size_t* permutation,const double* scale,const double* partial,double* state) {
    auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=n)return;
    CompensatedDeviceSum value;
    // Stable leaf order, no cross-block atomic updates.
    for(auto j=offsets[i];j<offsets[i+1];++j)value.add(partial[indices[j]]);
    state[permutation[i]]=scale[permutation[i]]*value.get();
}
template<class T> struct Buffer {
    T* data=nullptr;std::size_t size=0;int device=0;
#ifdef RADIA_HMATRIX_FAULT_INJECTION
    bool faultCounted=false;
#endif
    Buffer()=default;Buffer(const Buffer&)=delete;Buffer& operator=(const Buffer&)=delete;
    ~Buffer(){if(data){int previous=device;cudaGetDevice(&previous);if(previous!=device)cudaSetDevice(device);
        const auto status=cudaFree(data);
#ifdef RADIA_HMATRIX_FAULT_INJECTION
        if(faultCounted){if(status==cudaSuccess)deviceFaultState.liveBytes-=size*sizeof(T);else ++deviceFaultState.freeErrors;}
#else
        (void)status;
#endif
        if(previous!=device)cudaSetDevice(previous);}}
    void allocate(std::size_t n){
        if(data)throw std::runtime_error("Device buffer already allocated");
        const auto bytes=checkedProduct(n,sizeof(T));check(cudaGetDevice(&device));
#ifdef RADIA_HMATRIX_FAULT_INJECTION
        faultCounted=deviceFaultState.enabled;
        if(n && faultCounted){auto& s=deviceFaultState;++s.attempts;
            if(s.failAt && s.attempts==s.failAt)check(cudaErrorMemoryAllocation);
            if(bytes>std::numeric_limits<std::size_t>::max()-s.liveBytes)throw std::overflow_error("Device diagnostic byte overflow");}
#endif
        if(n)check(cudaMalloc(reinterpret_cast<void**>(&data),bytes));size=n;
#ifdef RADIA_HMATRIX_FAULT_INJECTION
        if(n && faultCounted){auto& s=deviceFaultState;s.liveBytes+=bytes;s.peakBytes=std::max(s.peakBytes,s.liveBytes);}
#endif
    }
    template<class Allocator>void upload(const std::vector<T,Allocator>& v){allocate(v.size());if(size){
        check(cudaMemcpy(data,v.data(),size*sizeof(T),cudaMemcpyHostToDevice));
        // Pageable H2D may return after staging but before DMA completion.
        // Consumers use nonblocking streams and do not inherit default-stream
        // ordering. Publish only after this upload's stream has completed.
        check(cudaStreamSynchronize(nullptr));
    }}
};
struct SparseMap {
    std::size_t rows=0,columns=0;
    HostVector<std::size_t> offsets,indices;
    HostVector<double> values;
    void validate()const {
        if(!rows||!columns||offsets.size()!=rows+1||offsets.front()!=0||offsets.back()!=values.size()||indices.size()!=values.size())
            throw std::invalid_argument("Invalid charge-map shape");
        for(std::size_t i=0;i<rows;++i)if(offsets[i]>offsets[i+1])throw std::invalid_argument("Invalid charge-map offsets");
        for(auto j:indices)if(j>=columns)throw std::invalid_argument("Invalid charge-map index");
        for(auto v:values)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite charge-map value");
    }
    SparseMap transpose()const {
        validate();SparseMap t;t.rows=columns;t.columns=rows;t.offsets.assign(columns+1,0);
        for(auto i:indices)++t.offsets[i+1];for(std::size_t i=0;i<columns;++i)t.offsets[i+1]+=t.offsets[i];
        t.indices.resize(indices.size());t.values.resize(values.size());auto next=t.offsets;
        for(std::size_t row=0;row<rows;++row)for(auto k=offsets[row];k<offsets[row+1];++k){auto q=next[indices[k]]++;t.indices[q]=row;t.values[q]=values[k];}return t;
    }
};
__global__ void sparseMapApply(std::size_t rows,const std::size_t* offsets,const std::size_t* indices,
                              const double* values,const double* x,double* y) {
    auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=rows)return;
    CompensatedDeviceSum result;for(auto k=offsets[i];k<offsets[i+1];++k)result.addProduct(values[k],x[indices[k]]);y[i]=result.get();
}
struct DeviceMap {
    Buffer<std::size_t> offsets,indices;Buffer<double> values;std::size_t rows=0;
    void upload(const SparseMap& map){rows=map.rows;offsets.upload(map.offsets);indices.upload(map.indices);values.upload(map.values);}
    void apply(const double* x,double* y,cudaStream_t stream){
        sparseMapApply<<<unsigned((rows+255)/256),256,0,stream>>>(rows,offsets.data,indices.data,values.data,x,y);check(cudaGetLastError());
    }
    void abandon(){offsets.data=nullptr;indices.data=nullptr;values.data=nullptr;}
};
} // namespace radia_hmatrix_cuda
#include "analytic_matrix_build.cuh"
#include "hdiv_curved_gram.cuh"
namespace radia_hmatrix_cuda {
class Resident {
    Buffer<Direction> directions;Buffer<double> coefficients,scale,input,tmp,partial,state;
    Buffer<std::size_t> permutation,rowOffsets,partialIndices;
    DeviceMap chargeMap,chargeTranspose;Buffer<double> dofState;std::size_t dofs=0;
    cudaStream_t stream=nullptr;std::size_t n,count;bool initialized=false,materialized=true;int device=0;
    void valid(bool incomplete=false){if(!materialized&&!incomplete)throw std::runtime_error("Deferred matrix is not materialized");int current;check(cudaGetDevice(&current));if(current!=device)throw std::runtime_error("H-matrix context belongs to another GPU");auto status=cudaStreamQuery(stream);if(status!=cudaSuccess&&status!=cudaErrorNotReady)check(status);}
public:
    std::size_t ownedBytes=0,setupPeakOwnedBytes=0,applications=0,kernelsPerApplication=4;
    bool assemblyProfileEnabled=false;double assemblyProfileSeconds[6]={};
    std::size_t assemblyFrameLeaves=0,assemblyFrameDenseCoefficients=0,assemblyFrameStorageCapCoefficients=0;
    std::size_t assemblyEntryBatches=0,assemblyEntryCount=0,assemblyLeafCount=0,assemblyLowrankCount=0;
    std::size_t assemblyCacheTileLoads=0,assemblyCacheTileHits=0,assemblyCacheFullTileBatches=0,assemblyCacheHostScratchBytes=0;
    std::size_t storedCoefficientCount=0,sameLeafDenseCoefficientCount=0,lowrankFactorCoefficientCount=0;
    std::size_t leafRankSum=0,maxLeafRank=0,storageSavingLowrankLeafCount=0;
    std::size_t deferredDenseLeaves=0,deferredDenseCoefficients=0,coefficientUploadBytes=0;
    explicit Resident(const Plan& p,const SparseMap* map=nullptr,std::size_t concurrentWorkspace=0):n(p.n),count(p.directions.size()) {
        check(cudaGetDevice(&device));
        materialized=p.deferredDense.empty();storedCoefficientCount=p.coefficientCount;coefficientUploadBytes=checkedProduct(p.coefficients.size(),sizeof(double));
        deferredDenseLeaves=p.deferredDense.size();for(const auto& span:p.deferredDense)deferredDenseCoefficients+=checkedProduct(span.rows,span.columns);
        for(const auto& direction:p.directions)if(!direction.transpose){
            const auto dense=checkedProduct(direction.inputs,direction.outputs);
            if(dense>std::numeric_limits<std::size_t>::max()-sameLeafDenseCoefficientCount)throw std::overflow_error("Leaf statistics overflow");
            sameLeafDenseCoefficientCount+=dense;
            if(direction.rank){
                const auto firstFactor=checkedProduct(direction.rank,direction.inputs),secondFactor=checkedProduct(direction.rank,direction.outputs);
                if(firstFactor>std::numeric_limits<std::size_t>::max()-secondFactor)throw std::overflow_error("Factor statistics overflow");
                const auto factors=firstFactor+secondFactor;
                if(factors>std::numeric_limits<std::size_t>::max()-lowrankFactorCoefficientCount||direction.rank>std::numeric_limits<std::size_t>::max()-leafRankSum)throw std::overflow_error("Rank statistics overflow");
                lowrankFactorCoefficientCount+=factors;leafRankSum+=direction.rank;maxLeafRank=std::max(maxLeafRank,direction.rank);
                if(factors<dense)++storageSavingLowrankLeafCount;
            }
        }
        if(count>2147483647u||n>std::size_t(2147483647u)*256)throw std::invalid_argument("CUDA grid range exceeded");
        auto bytes=[&](std::size_t size,std::size_t width){auto b=checkedProduct(size,width);if(b>std::numeric_limits<std::size_t>::max()-ownedBytes)throw std::overflow_error("Device allocation overflow");ownedBytes+=b;};
        bytes(count,sizeof(Direction));bytes(p.coefficientCount,sizeof(double));bytes(n,sizeof(double));
        bytes(n,2*sizeof(double));bytes(p.tmpCount,sizeof(double));bytes(p.partialCount,sizeof(double));
        bytes(n,sizeof(std::size_t));bytes(p.rowOffsets.size(),sizeof(std::size_t));bytes(p.partialIndices.size(),sizeof(std::size_t));
        std::size_t free,total;check(cudaMemGetInfo(&free,&total));const auto reserve=std::max(std::size_t(2)<<30,total/5);
        if(free<=reserve||ownedBytes>free-reserve)throw std::runtime_error("GPU memory reserve would be exceeded");
        SparseMap transpose;
        auto concurrentBudget=[&](){if(concurrentWorkspace&&
            (concurrentWorkspace>(std::size_t(2)<<30)||ownedBytes>(std::size_t(2)<<30)-concurrentWorkspace))throw std::runtime_error("Combined deferred setup workspace exceeds 2 GiB");};
        if(map){map->validate();if(map->rows!=n||map->columns>2147483647u)throw std::invalid_argument("Charge-map dimension mismatch");transpose=map->transpose();
            bytes(map->offsets.size(),sizeof(std::size_t));bytes(transpose.offsets.size(),sizeof(std::size_t));
            bytes(map->indices.size(),2*sizeof(std::size_t));bytes(map->values.size(),2*sizeof(double));bytes(map->columns,sizeof(double));
            if(ownedBytes>free-reserve)throw std::runtime_error("GPU charge-map memory reserve would be exceeded");
            concurrentBudget();
            dofs=map->columns;kernelsPerApplication=6;chargeMap.upload(*map);chargeTranspose.upload(transpose);dofState.allocate(dofs);
        }
        concurrentBudget();directions.upload(p.directions);
        if(materialized)coefficients.upload(p.coefficients);
        else{coefficients.allocate(p.coefficientCount);for(const auto& segment:p.uploads){
            if(segment.destination>p.coefficientCount||segment.count>p.coefficientCount-segment.destination||segment.source>p.coefficients.size()||segment.count>p.coefficients.size()-segment.source)throw std::invalid_argument("Invalid coefficient upload range");
            check(cudaMemcpy(coefficients.data+segment.destination,p.coefficients.data()+segment.source,segment.count*sizeof(double),cudaMemcpyHostToDevice));
        }}
        scale.upload(p.scale);
        permutation.upload(p.permutation);rowOffsets.upload(p.rowOffsets);partialIndices.upload(p.partialIndices);
        input.allocate(n);tmp.allocate(p.tmpCount);partial.allocate(p.partialCount);state.allocate(n);
        resourceFaultPoint("stream.create");
        check(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));resourceAcquired();
        try{resourceFaultPoint("stream.after-create");}
        catch(...){resourceReleased(cudaStreamDestroy(stream)==cudaSuccess);stream=nullptr;throw;}
    }
    template<class Oracle> void materializeDeferred(const Plan& p,Oracle& oracle){
        valid(true);if(materialized||p.deferredDense.empty()||p.n!=n||p.coefficientCount!=coefficients.size)throw std::invalid_argument("Invalid deferred materialization state");
        oracle.fillDense(p,coefficients.data,permutation.data,scale.data,stream);
        check(cudaStreamSynchronize(stream));materialized=true;
    }
    void assemble(const AnalyticBlock& block){
        valid();if(initialized||count!=1||coefficients.size!=checkedProduct(n,n)||block.sources.size()!=n)
            throw std::invalid_argument("Analytic assembly requires a fresh single dense block");
        std::size_t extra=checkedProduct(block.sources.size(),sizeof(AnalyticSource));
        const auto quadBytes=checkedProduct(block.quadrature.size(),sizeof(QuadraturePoint));
        if(quadBytes>std::numeric_limits<std::size_t>::max()-extra-sizeof(int))throw std::overflow_error("Assembly workspace overflow");extra+=quadBytes+sizeof(int);
        std::size_t free,total;check(cudaMemGetInfo(&free,&total));const auto reserve=std::max(std::size_t(2)<<30,total/5);
        if(free<=reserve||extra>free-reserve)throw std::runtime_error("Assembly GPU reserve would be exceeded");
        Buffer<AnalyticSource> sources;Buffer<QuadraturePoint> quadrature;Buffer<int> invalid;
        sources.upload(block.sources);quadrature.upload(block.quadrature);invalid.upload(HostVector<int>{0});setupPeakOwnedBytes=ownedBytes+extra;
        assembleAnalyticConstant<<<unsigned((n*n+63)/64),64,0,stream>>>(sources.data,quadrature.data,n,coefficients.data,invalid.data);check(cudaGetLastError());
        analyticPhysicalScale<<<unsigned((n+255)/256),256,0,stream>>>(n,coefficients.data,scale.data,invalid.data);check(cudaGetLastError());
        int flag=0;check(cudaMemcpyAsync(&flag,invalid.data,sizeof(int),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));
        if(flag)throw std::runtime_error("Analytic assembly produced invalid energy");
        normalizeAnalytic<<<unsigned((n*n+255)/256),256,0,stream>>>(n,coefficients.data,scale.data);check(cudaGetLastError());check(cudaStreamSynchronize(stream));
    }
    void assembleCurved(const CurvedBlock& block){
        valid();if(initialized||count!=1||coefficients.size!=checkedProduct(n,n)||block.charges.size()!=n)throw std::invalid_argument("Curved assembly requires a fresh dense block");
        assemblyProfileEnabled=block.profileAssembly;HDivAssemblyProfile profiler(block.profileAssembly,stream);
        HostVector<HostTargetPoint> hostPoints{HostFaultAllocator<HostTargetPoint>("resident.target-points")};HostVector<HostTargetSpan> hostSpans{HostFaultAllocator<HostTargetSpan>("resident.target-spans")};if(block.hostVectorCache)curvedHostTargets(block,hostPoints,hostSpans);
        const auto fullCacheCount=checkedProduct(n,hostPoints.size());
        const auto weightStride=std::max(block.tet.size(),block.tri.size());
        const auto mappedBytes=block.hostVectorCache?checkedProduct(hostPoints.size(),sizeof(MappedHostPoint))+checkedProduct(checkedProduct(n,weightStride),sizeof(double)):0;
        const auto cacheColumnBytes=checkedProduct(hostPoints.size(),sizeof(double));
        const auto cacheSources=hostPoints.empty()?std::size_t(0):std::min(n,std::max(std::size_t(1),(std::size_t(256)<<20)/cacheColumnBytes));
        const auto cacheCount=checkedProduct(cacheSources,hostPoints.size());
        const auto momentGroups=block.hostMomentCache?curvedSourceMomentGroups(block,HostFaultAllocator<MomentSourceGroup>("resident.source-groups")):HostVector<MomentSourceGroup>{HostFaultAllocator<MomentSourceGroup>("resident.source-groups")};
        HostVector<std::size_t> momentGroupIds{HostFaultAllocator<std::size_t>("resident.group-ids")},momentTileOffsets{HostFaultAllocator<std::size_t>("resident.tile-offsets")};
        if(block.hostMomentCache){for(std::size_t first=0;first<n;first+=cacheSources){const auto count=std::min(cacheSources,n-first);
            // Geometry-specific kernels keep the original per-component sums;
            // only independent source groups are scheduled separately.
            for(int kind=0;kind<2;++kind){momentTileOffsets.push_back(momentGroupIds.size());
                for(std::size_t g=0;g<momentGroups.size();++g){if(block.charges[momentGroups[g].member[0]].kind!=kind)continue;
                    bool active=false;for(std::size_t m=0;m<momentGroups[g].count;++m){const auto source=momentGroups[g].member[m];active|=source>=first&&source<first+count;}if(active)momentGroupIds.push_back(g);}
            }
        }momentTileOffsets.push_back(momentGroupIds.size());}

        const auto extra=checkedProduct(n,sizeof(CurvedCharge))+checkedProduct(block.tet.size()+block.tri.size(),sizeof(ReferencePoint))+checkedProduct(block.gl.size()+block.gw.size(),sizeof(double))+sizeof(int)
            +checkedProduct(hostPoints.size(),sizeof(HostTargetPoint))+checkedProduct(hostSpans.size(),sizeof(HostTargetSpan))+checkedProduct(cacheCount,sizeof(double))+mappedBytes+checkedProduct(momentGroups.size(),sizeof(MomentSourceGroup))+checkedProduct(momentGroupIds.size(),sizeof(std::size_t));
        if(extra>(std::size_t(2)<<30))throw std::invalid_argument("Curved host-vector workspace exceeds 2 GiB");
        std::size_t free,total;check(cudaMemGetInfo(&free,&total));const auto reserve=std::max(std::size_t(2)<<30,total/5);
        if(free<=reserve||extra>free-reserve)throw std::runtime_error("Curved assembly GPU reserve exceeded");
        Buffer<CurvedCharge> charges;Buffer<ReferencePoint> tet,tri;Buffer<double> gl,gw;Buffer<int> invalid;
        charges.upload(block.charges);tet.upload(block.tet);tri.upload(block.tri);gl.upload(block.gl);gw.upload(block.gw);invalid.upload(HostVector<int>{0});setupPeakOwnedBytes=ownedBytes+extra;
        Buffer<HostTargetPoint> targetPoints;Buffer<HostTargetSpan> targetSpans;Buffer<double> potentialCache,chargeWeights;Buffer<MappedHostPoint> mappedPoints;Buffer<MomentSourceGroup> deviceMomentGroups;Buffer<std::size_t> deviceMomentGroupIds;
        if(block.hostVectorCache){
            targetPoints.upload(hostPoints);targetSpans.upload(hostSpans);potentialCache.allocate(cacheCount);
            mappedPoints.allocate(hostPoints.size());chargeWeights.allocate(checkedProduct(n,weightStride));
            if(block.hostMomentCache){deviceMomentGroups.upload(momentGroups);deviceMomentGroupIds.upload(momentGroupIds);}
            profiler.begin(0);
            mapCurvedHostPoints<<<unsigned((hostPoints.size()+127)/128),128,0,stream>>>(charges.data,targetPoints.data,hostPoints.size(),mappedPoints.data);check(cudaGetLastError());
            mapCurvedChargeWeights<<<unsigned((n*weightStride+127)/128),128,0,stream>>>(charges.data,n,targetPoints.data,targetSpans.data,weightStride,chargeWeights.data);check(cudaGetLastError());
            profiler.end();
            // One stream orders cache production, consumption and reuse without host synchronization.
            for(std::size_t sourceFirst=0;sourceFirst<n;sourceFirst+=cacheSources){
                const auto sourcesHere=std::min(cacheSources,n-sourceFirst),entries=checkedProduct(sourcesHere,hostPoints.size());
                if(block.hostMomentCache){
                    const auto tile=sourceFirst/cacheSources;
                    for(int kind=0;kind<2;++kind){const auto offset=2*tile+std::size_t(kind),firstGroup=momentTileOffsets[offset],groupCount=momentTileOffsets[offset+1]-firstGroup,work=checkedProduct(groupCount,hostPoints.size());
                        profiler.begin(1+kind);
                        // Larger independent grids let near-field work from multiple
                        // source groups occupy the GPU without changing point sums.
                        constexpr std::size_t fusedBatch=std::size_t(1)<<22;
                        for(std::size_t first=0;first<work;first+=fusedBatch){const auto batch=std::min(fusedBatch,work-first);
                            if(kind==0)cacheCurvedHostMomentPotentials<0><<<unsigned((batch+127)/128),128,0,stream>>>(charges.data,n,targetPoints.data,targetPoints.size,targetSpans.data,mappedPoints.data,chargeWeights.data,weightStride,deviceMomentGroups.data,deviceMomentGroupIds.data+firstGroup,groupCount,gl.data,gw.data,int(gl.size),first,batch,sourceFirst,sourcesHere,potentialCache.data,invalid.data);
                            else cacheCurvedHostMomentPotentials<1><<<unsigned((batch+127)/128),128,0,stream>>>(charges.data,n,targetPoints.data,targetPoints.size,targetSpans.data,mappedPoints.data,chargeWeights.data,weightStride,deviceMomentGroups.data,deviceMomentGroupIds.data+firstGroup,groupCount,gl.data,gw.data,int(gl.size),first,batch,sourceFirst,sourcesHere,potentialCache.data,invalid.data);
                            check(cudaGetLastError());
                        }
                        profiler.end();
                    }
                }else{profiler.begin(3);for(std::size_t first=0;first<entries;first+=65536){const auto batch=std::min(std::size_t(65536),entries-first);
                    cacheCurvedHostPotentials<<<unsigned((batch+127)/128),128,0,stream>>>(charges.data,n,targetPoints.data,targetPoints.size,tet.data,tet.size,tri.data,tri.size,gl.data,gw.data,int(gl.size),first,batch,sourceFirst,targetSpans.data,mappedPoints.data,chargeWeights.data,weightStride,potentialCache.data,invalid.data);check(cudaGetLastError());
                }
                profiler.end();}
                profiler.begin(4);
                if(cacheCount==fullCacheCount){
                    assembleCachedCurvedGram<<<unsigned((n*n+127)/128),128,0,stream>>>(charges.data,n,targetPoints.data,targetPoints.size,targetSpans.data,potentialCache.data,chargeWeights.data,weightStride,coefficients.data,invalid.data);check(cudaGetLastError());
                }else{
                    assembleCachedCurvedColumns<<<unsigned((n*sourcesHere+127)/128),128,0,stream>>>(charges.data,n,targetPoints.data,targetPoints.size,targetSpans.data,potentialCache.data,chargeWeights.data,weightStride,sourceFirst,sourcesHere,coefficients.data,invalid.data);check(cudaGetLastError());
                }
                profiler.end();
            }
            profiler.begin(5);
            if(cacheCount!=fullCacheCount){symmetrizeCachedCurved<<<unsigned((n*n+127)/128),128,0,stream>>>(n,coefficients.data,invalid.data);check(cudaGetLastError());}
            check(cudaStreamSynchronize(stream));
        }else for(std::size_t first=0;first<n*n;first+=128){const auto batch=std::min(std::size_t(128),n*n-first);
            assembleCurvedGram<<<unsigned((batch+31)/32),32,0,stream>>>(charges.data,n,tet.data,tet.size,tri.data,tri.size,gl.data,gw.data,int(gl.size),first,batch,coefficients.data,invalid.data);check(cudaGetLastError());check(cudaStreamSynchronize(stream));
        }
        analyticPhysicalScale<<<unsigned((n+255)/256),256,0,stream>>>(n,coefficients.data,scale.data,invalid.data);check(cudaGetLastError());
        int flag=0;check(cudaMemcpyAsync(&flag,invalid.data,sizeof(int),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));if(flag)throw std::runtime_error("Curved assembly produced invalid energy");
        normalizeAnalytic<<<unsigned((n*n+255)/256),256,0,stream>>>(n,coefficients.data,scale.data);check(cudaGetLastError());profiler.end();check(cudaStreamSynchronize(stream));
        profiler.collect(assemblyProfileSeconds);
    }
    ~Resident(){if(stream){int previous=device;cudaGetDevice(&previous);if(previous!=device)cudaSetDevice(device);
        auto status=cudaStreamQuery(stream);
        if(status!=cudaSuccess&&status!=cudaErrorNotReady){directions.data=nullptr;coefficients.data=nullptr;scale.data=nullptr;input.data=nullptr;tmp.data=nullptr;partial.data=nullptr;state.data=nullptr;permutation.data=nullptr;rowOffsets.data=nullptr;partialIndices.data=nullptr;chargeMap.abandon();chargeTranspose.abandon();dofState.data=nullptr;}
        else{cudaStreamSynchronize(stream);resourceReleased(cudaStreamDestroy(stream)==cudaSuccess);}if(previous!=device)cudaSetDevice(previous);}}
    template<class Allocator>void upload(const std::vector<double,Allocator>& values){
        valid();const auto dimension=dofs?dofs:n;if(values.size()!=dimension)throw std::invalid_argument("State dimension mismatch");
        for(double v:values)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite state");
        check(cudaMemcpyAsync(dofs?dofState.data:state.data,values.data(),dimension*sizeof(double),cudaMemcpyHostToDevice,stream));
        check(cudaStreamSynchronize(stream));initialized=true;
    }
    std::size_t dimension()const{return dofs?dofs:n;}
    bool isDemag()const{return dofs!=0;}
    cudaStream_t executionStream(){valid();return stream;}
    // Internal device-to-device application for Krylov vectors. No host
    // vector transfer or synchronization; accepted uploaded state is untouched.
    void applyDevice(const double* x,double* y){
        valid();
        if(dofs)chargeMap.apply(x,state.data,stream);
        const double* source=dofs?state.data:x;
        gatherInput<<<unsigned((n+255)/256),256,0,stream>>>(n,permutation.data,scale.data,source,input.data);check(cudaGetLastError());
        projectRanks<<<unsigned(count),128,0,stream>>>(directions.data,count,coefficients.data,input.data,tmp.data);check(cudaGetLastError());
        applyLeaves<<<unsigned(count),128,0,stream>>>(directions.data,count,coefficients.data,input.data,tmp.data,partial.data);check(cudaGetLastError());
        double* target=dofs?state.data:y;
        reduceOutput<<<unsigned((n+255)/256),256,0,stream>>>(n,rowOffsets.data,partialIndices.data,permutation.data,scale.data,partial.data,target);check(cudaGetLastError());
        if(dofs)chargeTranspose.apply(state.data,y,stream);
    }
    void apply(std::size_t repeats){
        valid();if(!initialized)throw std::runtime_error("Upload initial state before apply");
        if(!repeats||repeats>10000)throw std::invalid_argument("Invalid repetition count");
        for(std::size_t k=0;k<repeats;++k){
            if(dofs)chargeMap.apply(dofState.data,state.data,stream);
            gatherInput<<<unsigned((n+255)/256),256,0,stream>>>(n,permutation.data,scale.data,state.data,input.data);check(cudaGetLastError());
            projectRanks<<<unsigned(count),128,0,stream>>>(directions.data,count,coefficients.data,input.data,tmp.data);check(cudaGetLastError());
            applyLeaves<<<unsigned(count),128,0,stream>>>(directions.data,count,coefficients.data,input.data,tmp.data,partial.data);check(cudaGetLastError());
            reduceOutput<<<unsigned((n+255)/256),256,0,stream>>>(n,rowOffsets.data,partialIndices.data,permutation.data,scale.data,partial.data,state.data);check(cudaGetLastError());
            if(dofs)chargeTranspose.apply(state.data,dofState.data,stream);
        }
        applications+=repeats;
    }
    HostVector<double> download(){
        valid();if(!initialized)throw std::runtime_error("No resident state");const auto dimension=dofs?dofs:n;HostVector<double> result(dimension);
        check(cudaMemcpyAsync(result.data(),dofs?dofState.data:state.data,dimension*sizeof(double),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));
        for(double v:result)if(!std::isfinite(v))throw std::runtime_error("Nonfinite resident result");return result;
    }
};
} // namespace radia_hmatrix_cuda
