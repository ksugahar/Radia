#include "csr_pcg_cuda.cuh"
namespace radia_amg_kernels {
#include "ams_amg_kernels.cuh"
}
namespace radia_csr {
Result solve_gpu_host(int n,int nnz,const int* row,const int* col,const double* val,
                      const double* rhs,const double* dinv,double tol,int maxiter,double* x) {
    validate_dimensions(n,nnz,tol,maxiter);validate_host(n,nnz,row,col,val,rhs,dinv);
    DeviceBuffer<int> rp(n+1),ci(std::max(nnz,1));
    DeviceBuffer<double> av(std::max(nnz,1)),b(n),di(n),solution(n);
    check(cudaMemcpy(rp.data,row,(size_t(n)+1)*sizeof(int),cudaMemcpyHostToDevice));
    check(cudaMemcpy(ci.data,col,size_t(nnz)*sizeof(int),cudaMemcpyHostToDevice));
    check(cudaMemcpy(av.data,val,size_t(nnz)*sizeof(double),cudaMemcpyHostToDevice));
    check(cudaMemcpy(b.data,rhs,size_t(n)*sizeof(double),cudaMemcpyHostToDevice));
    check(cudaMemcpy(di.data,dinv,size_t(n)*sizeof(double),cudaMemcpyHostToDevice));
    auto result=solve_device(n,nnz,rp.data,ci.data,av.data,b.data,di.data,tol,maxiter,solution.data);
    check(cudaMemcpy(x,solution.data,size_t(n)*sizeof(double),cudaMemcpyDeviceToHost));
    return result;
}
}

namespace radia_csr {
void amg_smooth_gpu_host(int n,int nnz,const int* row,const int* col,const double* val,
                         const double* rhs,const double* l1,int sweeps,double* output) {
    validate_dimensions(n,nnz,1e-6,sweeps);
    // The common host guard checks finite positive weights; l1 is not inverted.
    validate_host(n,nnz,row,col,val,rhs,l1);
    DeviceBuffer<int> rp(size_t(n)+1),ci(std::max(nnz,1));
    DeviceBuffer<double> av(std::max(nnz,1)),b(n),weights(n),x(n),next(n);
    check(cudaMemcpy(rp.data,row,(size_t(n)+1)*sizeof(int),cudaMemcpyHostToDevice));
    check(cudaMemcpy(ci.data,col,size_t(nnz)*sizeof(int),cudaMemcpyHostToDevice));
    check(cudaMemcpy(av.data,val,size_t(nnz)*sizeof(double),cudaMemcpyHostToDevice));
    check(cudaMemcpy(b.data,rhs,size_t(n)*sizeof(double),cudaMemcpyHostToDevice));
    check(cudaMemcpy(weights.data,l1,size_t(n)*sizeof(double),cudaMemcpyHostToDevice));
    const size_t blocks=(size_t(n)+threads-1)/threads;
    radia_amg_kernels::initial_relax<<<blocks,threads>>>(n,b.data,weights.data,x.data);
    check(cudaGetLastError());
    for(int k=1;k<sweeps;++k) {
        radia_amg_kernels::multiply_relax<<<blocks,threads>>>(n,rp.data,ci.data,av.data,x.data,b.data,weights.data,next.data);
        check(cudaGetLastError());std::swap(x.data,next.data);
    }
    check(cudaMemcpy(output,x.data,size_t(n)*sizeof(double),cudaMemcpyDeviceToHost));
    for(int i=0;i<n;++i)if(!std::isfinite(output[i]))
        throw std::runtime_error("AMG smoothing produced nonfinite output");
}
}
