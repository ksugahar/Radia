// Diagnostic-only charge-pair evaluation of the device flat Gram entry; never on
// the construction path. Explicit pairs cross H2D and their raw values D2H.
#pragma once
#include "hdiv_flat_gram.cuh"
namespace radia_hmatrix_cuda {
__global__ void flatPairEntries(FlatDevice d,const std::size_t* rows,const std::size_t* cols,std::size_t count,double* out){
    const std::size_t i=blockIdx.x*std::size_t(blockDim.x)+threadIdx.x;
    if(i<count)out[i]=flatEntry(d,rows[i],cols[i]);
}
// Diagnostic only: raw physical flatEntry values (the entry the dense and ACA
// routes write before normalization) for explicit charge pairs, from a fresh
// device preparation; the pairs cross H2D and the values D2H, no matrix does.
inline HostVector<double> flatDiagnosticEntries(const FlatBlock& block,const HostVector<std::size_t>& rows,const HostVector<std::size_t>& cols){
    block.validate();
    const auto n=block.charges(),count=rows.size();
    if(!count||count!=cols.size()||count>65536)throw std::invalid_argument("Flat diagnostic entries require 1..65536 matching charge pairs");
    for(std::size_t i=0;i<count;++i)if(rows[i]>=n||cols[i]>=n)throw std::invalid_argument("Flat diagnostic charge index out of range");
    flatReserve(flatSum(FlatPrepared::workspaceBytes(block),checkedProduct(count,2*sizeof(std::size_t)+sizeof(double))),
        "Flat diagnostic entries GPU reserve would be exceeded");
    cudaStream_t stream=nullptr;check(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
    struct Guard {cudaStream_t s;~Guard(){cudaStreamSynchronize(s);cudaStreamDestroy(s);}} guard{stream};
    HostVector<FlatTransfer> transfers{HostFaultAllocator<FlatTransfer>("flat.diagnostic-transfers")};std::size_t launches=0;int code=0;
    FlatPrepared prepared;prepared.prepare(block,stream,transfers,launches,code);
    Buffer<std::size_t> r,c;Buffer<double> values;r.upload(rows);c.upload(cols);values.allocate(count);
    flatPairEntries<<<unsigned((count+127)/128),128,0,stream>>>(prepared.device(),r.data,c.data,count,values.data);check(cudaGetLastError());
    HostVector<double> out(count);
    check(cudaMemcpyAsync(out.data(),values.data,count*sizeof(double),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));
    for(auto v:out)if(!std::isfinite(v))throw std::runtime_error("Flat diagnostic entry is nonfinite");
    return out;
}
} // namespace radia_hmatrix_cuda
