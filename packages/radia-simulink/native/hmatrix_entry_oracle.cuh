// Native assembly boundary: resident geometry, bounded ACA row/column batches.
// Include after hmatrix_resident.cuh; this is not an alternative ACA algorithm.
#pragma once
#include "../../../src/core/gpu/hdiv/hmatrix_resident.cuh"
#include "cHACApK_batch_size.h"
#include <mutex>
namespace radia_hmatrix_cuda {
constexpr std::size_t EntryBatchCapacity=RAD_HACAPK_BATCH_SIZE;
class AnalyticEntryOracle {
    Buffer<AnalyticSource> sources;
    Buffer<QuadraturePoint> quadrature;
    Buffer<AnalyticEntryIndex> indices;
    Buffer<double> values;
    Buffer<int> invalid;
    std::size_t capacity,sourceCount;
    int device;
    std::mutex access;
public:
    std::size_t ownedBytes=0,batches=0,entries=0,cacheTileLoads=0,cacheTileHits=0,cacheFullTileBatches=0,cacheHostScratchBytes=0;
    AnalyticEntryOracle(const AnalyticBlock& block,std::size_t batchCapacity)
        :capacity(batchCapacity),sourceCount(block.sources.size()){
        if(!sourceCount||!capacity||capacity>2147483647u)
            throw std::invalid_argument("Invalid analytic entry batch budget");
        check(cudaGetDevice(&device));
        auto add=[&](std::size_t count,std::size_t width){
            const auto bytes=checkedProduct(count,width);
            if(bytes>std::numeric_limits<std::size_t>::max()-ownedBytes)
                throw std::overflow_error("Entry oracle allocation overflow");
            ownedBytes+=bytes;
        };
        add(sourceCount,sizeof(AnalyticSource));add(block.quadrature.size(),sizeof(QuadraturePoint));
        add(capacity,sizeof(AnalyticEntryIndex)+sizeof(double));add(1,sizeof(int));
        // Validate metadata before a device kernel can dereference quadrature.
        for(const auto& source:block.sources){
            if(source.count!=(source.volume?64u:7u)||source.first>block.quadrature.size()||
               source.count>block.quadrature.size()-source.first)
                throw std::invalid_argument("Invalid analytic quadrature range");
        }
        std::size_t free,total;check(cudaMemGetInfo(&free,&total));
        const auto reserve=std::max(std::size_t(2)<<30,total/5);
        if(free<=reserve||ownedBytes>free-reserve)
            throw std::runtime_error("Entry oracle GPU memory reserve would be exceeded");
        sources.upload(block.sources);quadrature.upload(block.quadrature);
        indices.allocate(capacity);values.allocate(capacity);invalid.allocate(1);
    }
    template<class Allocator>HostVector<double> evaluate(const std::vector<AnalyticEntryIndex,Allocator>& request,cudaStream_t stream){
        // Host ACA leaf workers may share the geometry. Serialize this workspace;
        // no TaskManager, global callback state or per-entry CUDA launch is used.
        std::lock_guard<std::mutex> lock(access);
        int current;check(cudaGetDevice(&current));
        if(current!=device)throw std::runtime_error("Entry oracle belongs to another GPU");
        if(request.size()>capacity)throw std::invalid_argument("Entry batch exceeds workspace budget");
        for(const auto& index:request)if(index.row>=sourceCount||index.column>=sourceCount)
            throw std::invalid_argument("Entry index outside charge range");
        HostVector<double> result(request.size(),0.,HostFaultAllocator<double>("oracle.result"));if(request.empty())return result;
        check(cudaMemsetAsync(invalid.data,0,sizeof(int),stream));
        check(cudaMemcpyAsync(indices.data,request.data(),request.size()*sizeof(AnalyticEntryIndex),cudaMemcpyHostToDevice,stream));
        evaluateAnalyticEntries<<<unsigned((request.size()+63)/64),64,0,stream>>>(
            sources.data,quadrature.data,sourceCount,indices.data,request.size(),values.data,invalid.data);
        check(cudaGetLastError());int flag=0;
        check(cudaMemcpyAsync(result.data(),values.data,result.size()*sizeof(double),cudaMemcpyDeviceToHost,stream));
        check(cudaMemcpyAsync(&flag,invalid.data,sizeof(int),cudaMemcpyDeviceToHost,stream));
        check(cudaStreamSynchronize(stream));
        if(flag)throw std::runtime_error("Analytic entry evaluation failed");
        ++batches;entries+=request.size();return result;
    }
};
} // namespace radia_hmatrix_cuda
