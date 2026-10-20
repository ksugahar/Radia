// Bounded high-order/curved entry boundary. No dense Gram allocation.
#pragma once
#include "hmatrix_entry_oracle.cuh"
#include "hmatrix_curved_entry_cache.cuh"
namespace radia_hmatrix_cuda {
__global__ void evaluateCurvedEntries(const CurvedCharge* charges,std::size_t n,
    const ReferencePoint* tet,std::size_t nt,const ReferencePoint* tri,std::size_t ns,
    const double* gl,const double* gw,int nq,const AnalyticEntryIndex* indices,
    std::size_t count,double* values,int* invalid){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=count)return;
    auto row=indices[i].row,col=indices[i].column;
    if(row>=n||col>=n){atomicExch(invalid,1);values[i]=0.;return;}
    // Match the upper-triangle evaluation and ordered symmetric averaging of
    // assembleCurvedGram, including touching Duffy and curved reference frames.
    if(row>col){const auto tmp=row;row=col;col=tmp;}
    double value=curvedQuadDot(charges[row],charges[col],tet,nt,tri,ns,gl,gw,nq);
    if(row!=col)value=.5*(value+curvedQuadDot(charges[col],charges[row],tet,nt,tri,ns,gl,gw,nq));
    if(!isfinite(value))atomicExch(invalid,1);values[i]=value;
}
class CurvedEntryOracle {
    Buffer<CurvedCharge> charges;
    Buffer<ReferencePoint> tet,tri;
    Buffer<double> gl,gw,values;
    Buffer<AnalyticEntryIndex> indices;
    Buffer<int> invalid;
    std::size_t capacity,sourceCount;
    int device;
    std::mutex access;
    HostUnique<CurvedCachedEntryOracle> cached;
public:
    std::size_t ownedBytes=0,batches=0,entries=0,cacheTileLoads=0,cacheTileHits=0,cacheFullTileBatches=0,cacheHostScratchBytes=0;
    CurvedEntryOracle(const CurvedBlock& block,std::size_t batchCapacity,bool useCache=false,std::size_t cacheBudget=std::size_t(256)<<20,bool deviceDense=false)
        :capacity(batchCapacity),sourceCount(block.charges.size()){
        if(!sourceCount||sourceCount>8192||!capacity||capacity>EntryBatchCapacity||block.tet.empty()||block.tri.empty()||
           block.gl.empty()||block.gl.size()!=block.gw.size()||block.gl.size()>2147483647u)
            throw std::invalid_argument("Invalid curved entry oracle budget");
        check(cudaGetDevice(&device));
        if(deviceDense&&!useCache)throw std::invalid_argument("Device dense fill requires cached entries");
        if(useCache){cached=makeHostUnique<CurvedCachedEntryOracle>("oracle.cached-object",block,capacity,cacheBudget,deviceDense);ownedBytes=cached->ownedBytes;return;}
        auto add=[&](std::size_t count,std::size_t width){const auto bytes=checkedProduct(count,width);
            if(bytes>std::numeric_limits<std::size_t>::max()-ownedBytes)throw std::overflow_error("Curved oracle allocation overflow");
            ownedBytes+=bytes;};
        add(sourceCount,sizeof(CurvedCharge));add(block.tet.size()+block.tri.size(),sizeof(ReferencePoint));
        add(block.gl.size(),2*sizeof(double));add(capacity,sizeof(AnalyticEntryIndex)+sizeof(double));add(1,sizeof(int));
        std::size_t free,total;check(cudaMemGetInfo(&free,&total));const auto reserve=std::max(std::size_t(2)<<30,total/5);
        if(free<=reserve||ownedBytes>free-reserve)throw std::runtime_error("Curved entry oracle GPU reserve exceeded");
        charges.upload(block.charges);tet.upload(block.tet);tri.upload(block.tri);gl.upload(block.gl);gw.upload(block.gw);
        indices.allocate(capacity);values.allocate(capacity);invalid.allocate(1);
    }
    void fillDense(const Plan& p,double* coefficients,const std::size_t* permutation,const double* scale,cudaStream_t stream){
        if(!cached)throw std::invalid_argument("Missing device dense cache");cached->fillDense(p,coefficients,permutation,scale,stream);
        cacheTileLoads=cached->tileLoads;cacheTileHits=cached->tileHits;
    }
    template<class Allocator>HostVector<double> evaluate(const std::vector<AnalyticEntryIndex,Allocator>& request,cudaStream_t stream){
        std::lock_guard<std::mutex> lock(access);int current;check(cudaGetDevice(&current));
        if(current!=device)throw std::runtime_error("Curved entry oracle belongs to another GPU");
        if(request.size()>capacity)throw std::invalid_argument("Curved entry batch exceeds workspace budget");
        for(const auto& index:request)if(index.row>=sourceCount||index.column>=sourceCount)
            throw std::invalid_argument("Curved entry index outside charge range");
        HostVector<double> result(request.size(),0.,HostFaultAllocator<double>("oracle.result"));if(request.empty())return result;
        if(cached){result=cached->evaluate(request,stream);cacheTileLoads=cached->tileLoads;cacheTileHits=cached->tileHits;cacheFullTileBatches=cached->fullTileBatches;cacheHostScratchBytes=cached->hostScratchBytes;++batches;entries+=request.size();return result;}
        check(cudaMemsetAsync(invalid.data,0,sizeof(int),stream));
        check(cudaMemcpyAsync(indices.data,request.data(),request.size()*sizeof(AnalyticEntryIndex),cudaMemcpyHostToDevice,stream));
        evaluateCurvedEntries<<<unsigned((request.size()+31)/32),32,0,stream>>>(charges.data,sourceCount,
            tet.data,tet.size,tri.data,tri.size,gl.data,gw.data,int(gl.size),indices.data,request.size(),values.data,invalid.data);
        check(cudaGetLastError());int flag=0;
        check(cudaMemcpyAsync(result.data(),values.data,result.size()*sizeof(double),cudaMemcpyDeviceToHost,stream));
        check(cudaMemcpyAsync(&flag,invalid.data,sizeof(int),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));
        if(flag)throw std::runtime_error("Curved entry evaluation failed");
        ++batches;entries+=request.size();return result;
    }
};
}
