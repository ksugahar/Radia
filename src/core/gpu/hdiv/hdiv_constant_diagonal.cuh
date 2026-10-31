// Constant charge-space self entries, using caller-staged physical quadrature.
#pragma once
#include "hdiv_constant_diagonal_input.hpp"
#include "hmatrix_resident.cuh"
#include <mutex>

namespace radia_hmatrix_cuda {
__global__ void constantChargeDiagonal(const AnalyticSource* sources,
        const QuadraturePoint* quadrature,std::size_t n,double* diagonal,
        double* scale,int* invalid){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=n)return;
    const double value=analyticQuadDot(sources,quadrature,i,i);
    if(!isfinite(value)||value<=0.){atomicExch(invalid,1);return;}
    diagonal[i]=value;scale[i]=sqrt(value);
}
}
namespace radia_hdiv_python {
namespace constant_diagonal_detail {
using namespace radia_hmatrix_cuda;
struct Stream {
    cudaStream_t value=nullptr;int device=0;
    Stream(){
        check(cudaGetDevice(&device));resourceFaultPoint("constant-diagonal.stream.create");
        check(cudaStreamCreateWithFlags(&value,cudaStreamNonBlocking));resourceAcquired();
        try{resourceFaultPoint("constant-diagonal.stream.after-create");}
        catch(...){resourceReleased(cudaStreamDestroy(value)==cudaSuccess);value=nullptr;throw;}
    }
    ~Stream(){
        if(!value)return;
        int previous=device;cudaGetDevice(&previous);
        if(previous!=device)cudaSetDevice(device);
        resourceReleased(cudaStreamDestroy(value)==cudaSuccess);
        if(previous!=device)cudaSetDevice(previous);
    }
    Stream(const Stream&)=delete;Stream& operator=(const Stream&)=delete;
};
struct Workspace {
    Buffer<double> diagonal,scale;
    // Destroy/synchronize stream before releasing its buffers.
    Stream stream;
};
class Owner final:public ConstantDiagonal {
    std::mutex mutex;
    HostUnique<Workspace> workspace;
    ConstantDiagonalStats record;
    HostVector<double> download(bool scale){
        if(!workspace)throw std::runtime_error("Constant charge diagonal has been closed");
        int device;check(cudaGetDevice(&device));
        if(device!=workspace->stream.device)throw std::runtime_error("Constant charge diagonal belongs to another GPU");
        const auto& buffer=scale?workspace->scale:workspace->diagonal;
        HostVector<double> out(record.charges);
        if(record.diagnosticDownloadBytes>std::numeric_limits<std::size_t>::max()-out.size()*sizeof(double))
            throw std::overflow_error("Constant charge diagnostic byte counter overflow");
        check(cudaStreamSynchronize(workspace->stream.value));
        check(cudaMemcpy(out.data(),buffer.data,out.size()*sizeof(double),cudaMemcpyDeviceToHost));
        record.diagnosticDownloadBytes+=out.size()*sizeof(double);return out;
    }
public:
    explicit Owner(const ConstantDiagonalInput& input){
        // Validate before cudaMemGetInfo, allocation, stream creation or kernel launch.
        input.validate();record.charges=input.charges.size();
        const auto sourceBytes=checkedProduct(record.charges,sizeof(AnalyticSource));
        const auto quadratureBytes=checkedProduct(input.quadrature.size(),sizeof(QuadraturePoint));
        record.ownedDeviceBytes=checkedProduct(record.charges,2*sizeof(double));
        record.setupPeakOwnedBytes=record.ownedDeviceBytes+sourceBytes+quadratureBytes+sizeof(int);
        if(record.setupPeakOwnedBytes>(std::size_t(2)<<30))
            throw std::invalid_argument("Constant charge setup exceeds existing 2 GiB budget");
        std::size_t free,total;check(cudaMemGetInfo(&free,&total));
        const auto reserve=std::max(std::size_t(2)<<30,total/5);
        if(free<=reserve||record.setupPeakOwnedBytes>free-reserve)
            throw std::runtime_error("Constant charge diagonal exceeds GPU memory reserve");
        HostVector<AnalyticSource> sources;HostVector<QuadraturePoint> quadrature;
        sources.reserve(record.charges);quadrature.reserve(input.quadrature.size());
        for(const auto& charge:input.charges){
            AnalyticSource source{};std::copy(charge.vertices.begin(),charge.vertices.end(),source.vertices);
            source.first=charge.first;source.count=charge.count;source.volume=charge.kind==0;
            sources.push_back(source);
        }
        for(const auto& point:input.quadrature){
            QuadraturePoint q{};std::copy(point.point,point.point+3,q.point);q.weight=point.weight;
            quadrature.push_back(q);
        }
        Buffer<AnalyticSource> deviceSources;Buffer<QuadraturePoint> deviceQuadrature;Buffer<int> invalid;
        auto pending=makeHostUnique<Workspace>("constant-diagonal.workspace");
        deviceSources.upload(sources);deviceQuadrature.upload(quadrature);invalid.allocate(1);
        pending->diagonal.allocate(record.charges);pending->scale.allocate(record.charges);
        check(cudaMemsetAsync(invalid.data,0,sizeof(int),pending->stream.value));
        constantChargeDiagonal<<<unsigned((record.charges+127)/128),128,0,pending->stream.value>>>(
            deviceSources.data,deviceQuadrature.data,record.charges,pending->diagonal.data,pending->scale.data,invalid.data);
        check(cudaGetLastError());check(cudaStreamSynchronize(pending->stream.value));
        int failed=0;check(cudaMemcpy(&failed,invalid.data,sizeof(int),cudaMemcpyDeviceToHost));
        if(failed)throw std::runtime_error("Constant charge physical diagonal is nonpositive or nonfinite");
        record.setupUploadBytes=sourceBytes+quadratureBytes;record.setupDownloadBytes=sizeof(int);
        workspace=std::move(pending);
    }
    HostVector<double> diagnostic_diagonal()override{std::lock_guard<std::mutex> lock(mutex);return download(false);}
    HostVector<double> diagnostic_scale()override{std::lock_guard<std::mutex> lock(mutex);return download(true);}
    ConstantDiagonalStats stats()override{std::lock_guard<std::mutex> lock(mutex);return record;}
    void close()override{
        std::lock_guard<std::mutex> lock(mutex);workspace.reset();
        record.ownedDeviceBytes=0;record.closed=true;
    }
};
}
std::shared_ptr<ConstantDiagonal> create_constant_diagonal(const ConstantDiagonalInput& input){
    return std::allocate_shared<constant_diagonal_detail::Owner>(
        radia_hmatrix_cuda::HostFaultAllocator<constant_diagonal_detail::Owner>("constant-diagonal.owner"),input);
}
} // namespace radia_hdiv_python
