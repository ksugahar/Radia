// Diagnostic-only ACA leaf on a described matrix; never on the construction path.
#pragma once
#include "hdiv_flat_aca.cuh"
namespace radia_hmatrix_cuda {
// Diagnostic: diag(d) + w w^T from its two vectors, so a chosen stopping path can be driven without any kernel.
struct AcaDescribedEntry {const double* diagonal;const double* outer;
    __device__ double operator()(const AcaLeafState&,std::size_t i,std::size_t j)const{return (i==j?diagonal[i]:0.0)+outer[i]*outer[j];}};
// Diagnostic only, never on the construction path: rank, status and steps of one square leaf diag(d) + w w^T.
// The two vectors describe it and cross H2D; the entries are evaluated on the device and the factors stay there.
struct AcaLeafReport {int rank,status,steps;};
inline AcaLeafReport acaDiagnosticLeaf(const HostVector<double>& diagonal,const HostVector<double>& outer,double eps,int maxRank){
    const auto n=diagonal.size();
    if(!n||outer.size()!=n||maxRank<1||!(eps>0.0))throw std::invalid_argument("ACA diagnostic leaf requires matching vectors, a positive cap and a positive eps");
    const auto leafRank=std::min(std::size_t(maxRank),n),scratch=2*n*leafRank+2*n;
    flatReserve(checkedProduct(scratch+2*n,sizeof(double)),"ACA diagnostic leaf would exceed the GPU reserve");
    cudaStream_t stream=nullptr;check(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
    struct Guard {cudaStream_t s;~Guard(){cudaStreamSynchronize(s);cudaStreamDestroy(s);}} guard{stream};
    Buffer<double> deviceDiagonal,deviceOuter,work;deviceDiagonal.upload(diagonal);deviceOuter.upload(outer);work.allocate(scratch);
    check(cudaMemsetAsync(work.data,0,scratch*sizeof(double),stream));
    HostVector<AcaLeafState> states{HostFaultAllocator<AcaLeafState>("aca.diagnostic-state")};
    states.push_back(AcaLeafState{0,0,n,n,0,n*leafRank,2*n*leafRank,2*n*leafRank+n,0.,0.,0,acaActive,0,0,maxRank,0});
    Buffer<AcaLeafState> deviceStates;deviceStates.upload(states);Buffer<int> active;active.allocate(1);
    int running=1;
    for(std::size_t step=0;running&&step<n+leafRank+acaProbes;++step){
        check(cudaMemsetAsync(active.data,0,sizeof(int),stream));
        acaStep<AcaDescribedEntry><<<1u,acaThreads,0,stream>>>(AcaDescribedEntry{deviceDiagonal.data,deviceOuter.data},deviceStates.data,work.data,eps,active.data);
        check(cudaGetLastError());
        check(cudaMemcpyAsync(&running,active.data,sizeof(int),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));}
    if(running)throw std::runtime_error("ACA diagnostic leaf did not terminate within its step bound");
    AcaLeafState out;check(cudaMemcpyAsync(&out,deviceStates.data,sizeof(AcaLeafState),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));
    return {out.rank,out.status,out.steps};
}
} // namespace radia_hmatrix_cuda
