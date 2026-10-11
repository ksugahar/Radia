#include "optuna_cuda.h"
#include <cuda_runtime.h>
#include <cfloat>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace optuna_cuda {
namespace {
void check(cudaError_t status) {
    if (status != cudaSuccess)
        throw std::runtime_error(std::string("Optuna CUDA: ") + cudaGetErrorString(status));
}
struct Buffer {
    double* p = nullptr;
    explicit Buffer(std::size_t n) {
        if (n > std::numeric_limits<std::size_t>::max() / sizeof(double))
            throw std::overflow_error("Optuna CUDA buffer size overflow");
        if (n) check(cudaMalloc(reinterpret_cast<void**>(&p), n * sizeof(double)));
    }
    ~Buffer() { if (p) cudaFree(p); }
    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;
    void upload(const std::vector<double>& a) {
        if (!a.empty()) check(cudaMemcpy(p, a.data(), a.size()*sizeof(double), cudaMemcpyHostToDevice));
    }
};
__device__ double log_mass(double lo, double hi) {
    constexpr double inv = 0.70710678118654752440;
    double mass;
    if (hi <= 0) mass = .5*(erfc(-hi*inv)-erfc(-lo*inv));
    else if (lo >= 0) mass = .5*(erfc(lo*inv)-erfc(hi*inv));
    else mass = .5*(erf(hi*inv)-erf(lo*inv));
    return log(fmax(mass, DBL_MIN));
}
__device__ double component(double x, double w, double mu, double sigma,
                            double den, bool logarithmic, double step) {
    if (w == 0) return -INFINITY;
    double value;
    if (isfinite(step)) {
        double lo=x-.5*step, hi=x+.5*step;
        if (logarithmic) { lo=log(lo); hi=log(hi); }
        value=log_mass((lo-mu)/sigma,(hi-mu)/sigma)-den;
    } else {
        double z=((logarithmic?log(x):x)-mu)/sigma;
        value=-.5*z*z-0.91893853320467274178-log(sigma)-den;
    }
    return value+log(w);
}
__global__ void numerical_kernel(const double* x, std::size_t n,
    const double* w, const double* mu, const double* sigma, const double* den,
    std::size_t k, bool logarithmic, double step, double* result) {
    for (std::size_t i=blockIdx.x*std::size_t(blockDim.x)+threadIdx.x; i<n;
         i+=std::size_t(blockDim.x)*gridDim.x) {
        double maximum=-INFINITY;
        for (std::size_t j=0;j<k;++j)
            maximum=fmax(maximum,component(x[i],w[j],mu[j],sigma[j],den[j],logarithmic,step));
        double sum=0;
        for (std::size_t j=0;j<k;++j)
            sum+=exp(component(x[i],w[j],mu[j],sigma[j],den[j],logarithmic,step)-maximum);
        result[i]=maximum+log(sum);
    }
}
__global__ void categorical_kernel(const double* x,std::size_t n,
    const double* w,const double* probability,std::size_t k,std::size_t choices,double* result) {
    for (std::size_t i=blockIdx.x*std::size_t(blockDim.x)+threadIdx.x;i<n;
         i+=std::size_t(blockDim.x)*gridDim.x) {
        std::size_t choice=static_cast<std::size_t>(x[i])-1;
        double maximum=-INFINITY;
        for (std::size_t j=0;j<k;++j)
            maximum=fmax(maximum,log(fmax(probability[j*choices+choice],DBL_MIN))+log(w[j]));
        double sum=0;
        for (std::size_t j=0;j<k;++j)
            sum+=exp(log(fmax(probability[j*choices+choice],DBL_MIN))+log(w[j])-maximum);
        result[i]=maximum+log(sum);
    }
}
unsigned blocks(std::size_t n) { return static_cast<unsigned>(std::min<std::size_t>((n-1)/128+1,65535)); }
__device__ double joint_component(const double* x, std::size_t row, std::size_t kernel,
    double weight, const double* parameters, const Dimension* dimensions, std::size_t d) {
    double value=log(weight);
    for (std::size_t j=0;j<d;++j) {
        const Dimension dim=dimensions[j];
        const double sample=x[row*d+j];
        if (dim.categorical) {
            std::size_t choice=static_cast<std::size_t>(sample)-1;
            value+=log(fmax(parameters[dim.probabilities+kernel*dim.choices+choice],DBL_MIN));
        } else {
            value+=component(sample,1.0,parameters[dim.mu+kernel],parameters[dim.sigma+kernel],
                             parameters[dim.denominator+kernel],dim.log_scale,dim.step);
        }
    }
    return value;
}
__global__ void joint_kernel(const double* x,std::size_t n,const double* w,std::size_t k,
    const double* parameters,const Dimension* dimensions,std::size_t d,double* result) {
    for (std::size_t i=blockIdx.x*std::size_t(blockDim.x)+threadIdx.x;i<n;
         i+=std::size_t(blockDim.x)*gridDim.x) {
        double maximum=-INFINITY;
        for (std::size_t j=0;j<k;++j)
            maximum=fmax(maximum,joint_component(x,i,j,w[j],parameters,dimensions,d));
        double sum=0;
        for (std::size_t j=0;j<k;++j)
            sum+=exp(joint_component(x,i,j,w[j],parameters,dimensions,d)-maximum);
        result[i]=maximum+log(sum);
    }
}
std::vector<double> download(Buffer& b,std::size_t n) {
    std::vector<double> out(n);
    check(cudaMemcpy(out.data(),b.p,n*sizeof(double),cudaMemcpyDeviceToHost));
    for (double x:out) if (!std::isfinite(x))
        throw std::runtime_error("Optuna CUDA density returned a nonfinite value");
    return out;
}
}
std::vector<double> numerical(const std::vector<double>& x,const std::vector<double>& w,
    const std::vector<double>& mu,const std::vector<double>& sigma,
    const std::vector<double>& den,bool logarithmic,double step) {
    if (x.empty()) return {};
    Buffer dx(x.size()),dw(w.size()),dm(mu.size()),ds(sigma.size()),dd(den.size()),dy(x.size());
    dx.upload(x);dw.upload(w);dm.upload(mu);ds.upload(sigma);dd.upload(den);
    numerical_kernel<<<blocks(x.size()),128>>>(dx.p,x.size(),dw.p,dm.p,ds.p,dd.p,w.size(),logarithmic,step,dy.p);
    check(cudaGetLastError());
    return download(dy,x.size());
}
std::vector<double> categorical(const std::vector<double>& x,const std::vector<double>& w,
    const std::vector<double>& probabilities,std::size_t choices) {
    if (x.empty()) return {};
    Buffer dx(x.size()),dw(w.size()),dp(probabilities.size()),dy(x.size());
    dx.upload(x);dw.upload(w);dp.upload(probabilities);
    categorical_kernel<<<blocks(x.size()),128>>>(dx.p,x.size(),dw.p,dp.p,w.size(),choices,dy.p);
    check(cudaGetLastError());
    return download(dy,x.size());
}
std::vector<double> joint(const std::vector<double>& x,std::size_t n,
    const std::vector<double>& w,const std::vector<double>& parameters,
    const std::vector<Dimension>& dimensions) {
    if (!n) return {};
    // Allocate a byte-sized descriptor buffer with the same checked RAII owner.
    if (dimensions.size() > std::numeric_limits<std::size_t>::max()/sizeof(Dimension))
        throw std::overflow_error("Optuna CUDA descriptor size overflow");
    const std::size_t bytes=dimensions.size()*sizeof(Dimension);
    Buffer dx(x.size()),dw(w.size()),dp(parameters.size()),dy(n),dd((bytes-1)/sizeof(double)+1);
    dx.upload(x);dw.upload(w);dp.upload(parameters);
    check(cudaMemcpy(dd.p,dimensions.data(),bytes,cudaMemcpyHostToDevice));
    joint_kernel<<<blocks(n),128>>>(dx.p,n,dw.p,w.size(),dp.p,
        reinterpret_cast<const Dimension*>(dd.p),dimensions.size(),dy.p);
    check(cudaGetLastError());
    return download(dy,n);
}
}
