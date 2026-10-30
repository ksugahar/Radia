// Shared CUDA inverse-law evaluation; no FE, MEX, Python or Newton ownership.
#pragma once
#include "hdiv_inverse_bh.hpp"
#include "hmatrix_resident.cuh"
namespace radia_hdiv_bh {
__global__ void inverseLawKernel(const double* a,std::size_t n,const double* knots,
    std::size_t nk,const double* c,const double* base,double cap,double barrier,
    double last,double initial,double* result){
    const auto q=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(q>=n)return;
    const double amplitude=a[q],clipped=fmin(amplitude,cap);
    std::size_t lo=0,hi=nk-2;
    while(lo<hi){const auto mid=lo+(hi-lo+1)/2;if(clipped>=knots[mid])lo=mid;else hi=mid-1;}
    const double x=clipped-knots[lo];const auto* cc=c+4*lo;
    double h=((cc[0]*x+cc[1])*x+cc[2])*x+cc[3];
    double tangent=(3*cc[0]*x+2*cc[1])*x+cc[2];
    double energy=base[lo]+cc[0]*x*x*x*x/4+cc[1]*x*x*x/3+cc[2]*x*x/2+cc[3]*x;
    if(amplitude>cap){const double extra=amplitude-cap;h+=barrier*extra;tangent=barrier;energy+=last*extra+.5*barrier*extra*extra;}
    result[q]=amplitude>0?h/amplitude:initial;result[n+q]=tangent;result[2*n+q]=energy;
}
inline std::size_t validateInverseBH(const Law& law,const radia_hmatrix_cuda::HostVector<double>& amplitudes){
    using namespace radia_hmatrix_cuda;
    const auto nk=law.knots.size(),n=amplitudes.size();
    if(nk<2||nk>4096||!n||n>65536||law.coefficients.size()!=4*(nk-1)||law.integralBase.size()!=nk)
        throw std::invalid_argument("Inverse BH payload shape or budget exceeded");
    if(!std::isfinite(law.cap)||law.cap<=0||!std::isfinite(law.barrierSlope)||law.barrierSlope<=0||
        !std::isfinite(law.lastField)||law.lastField<=0||!std::isfinite(law.initialReluctivity)||law.initialReluctivity<=0||
        law.knots.front()!=0||law.cap>law.knots.back())throw std::invalid_argument("Invalid inverse BH scalar contract");
    for(std::size_t i=0;i<nk;++i)if(!std::isfinite(law.knots[i])||!std::isfinite(law.integralBase[i])||(i&&law.knots[i]<=law.knots[i-1]))throw std::invalid_argument("Invalid inverse BH knots or energy bases");
    for(double x:law.coefficients)if(!std::isfinite(x))throw std::invalid_argument("Nonfinite inverse BH coefficient");
    auto consistent=[](double a,double b){return std::abs(a-b)<=128*std::numeric_limits<double>::epsilon()*std::max(1.,std::max(std::abs(a),std::abs(b)));};
    if(!consistent(law.coefficients[3],0)||!consistent(law.integralBase[0],0)||!consistent(law.initialReluctivity,law.coefficients[2]))
        throw std::invalid_argument("Inverse BH origin contract mismatch");
    for(std::size_t i=0;i<nk-1;++i){const auto* c=law.coefficients.data()+4*i;const double x=law.knots[i+1]-law.knots[i];
        const double area=c[0]*x*x*x*x/4+c[1]*x*x*x/3+c[2]*x*x/2+c[3]*x;
        if(!consistent(law.integralBase[i+1],law.integralBase[i]+area)||(i+1<nk-1&&!consistent(((c[0]*x+c[1])*x+c[2])*x+c[3],law.coefficients[4*(i+1)+3])))
            throw std::invalid_argument("Inverse BH interval continuity mismatch");
    }
    const auto capInterval=std::min(std::size_t(std::upper_bound(law.knots.begin(),law.knots.end(),law.cap)-law.knots.begin()-1),nk-2);
    const auto* cc=law.coefficients.data()+4*capInterval;const double capX=law.cap-law.knots[capInterval];
    if(!consistent(law.lastField,((cc[0]*capX+cc[1])*capX+cc[2])*capX+cc[3]))throw std::invalid_argument("Inverse BH barrier continuity mismatch");
    for(double x:amplitudes)if(!std::isfinite(x)||x<0)throw std::invalid_argument("Inverse BH amplitude must be finite and nonnegative");
    const auto entries=n+nk+law.coefficients.size()+law.integralBase.size()+3*n;
    const auto bytes=checkedProduct(entries,sizeof(double));
    if(bytes>(std::size_t(8)<<20))throw std::invalid_argument("Inverse BH device payload budget exceeded");
    return bytes;
}
inline Samples evaluateInverseBH(const Law& law,const radia_hmatrix_cuda::HostVector<double>& amplitudes){
    using namespace radia_hmatrix_cuda;const auto n=amplitudes.size(),nk=law.knots.size();
    const auto bytes=validateInverseBH(law,amplitudes);
    Buffer<double> a,knots,c,base,result;a.upload(amplitudes);knots.upload(law.knots);c.upload(law.coefficients);base.upload(law.integralBase);result.allocate(3*n);
    inverseLawKernel<<<unsigned((n+127)/128),128>>>(a.data,n,knots.data,nk,c.data,base.data,law.cap,law.barrierSlope,law.lastField,law.initialReluctivity,result.data);
    check(cudaGetLastError());HostVector<double> values(3*n,HostFaultAllocator<double>("bh.scalar-download"));check(cudaMemcpy(values.data(),result.data,values.size()*sizeof(double),cudaMemcpyDeviceToHost));
    Samples out;out.ownedDeviceBytes=bytes;out.secant.assign(values.begin(),values.begin()+n);out.tangent.assign(values.begin()+n,values.begin()+2*n);out.energy.assign(values.begin()+2*n,values.end());
    for(std::size_t i=0;i<n;++i)if(!std::isfinite(out.secant[i])||!std::isfinite(out.tangent[i])||!std::isfinite(out.energy[i])||out.secant[i]<=0||out.tangent[i]<=0)throw std::runtime_error("Invalid inverse BH tangent or energy");
    return out;
}
}
