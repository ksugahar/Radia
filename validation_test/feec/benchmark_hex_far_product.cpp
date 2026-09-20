// Isolated numerical/parity benchmark for the production far contraction.
#include "../../src/core/rad_hex_far_product.h"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <random>
#include <stdexcept>

using Clock = std::chrono::steady_clock;

int benchmark()
{
    std::mt19937 generator(20260908);
    std::uniform_real_distribution<double> random(-1.0, 1.0);
    std::cout << std::setprecision(17) << "{\"cases\":[";
    bool first = true;
    for (const auto dimensions : {std::vector<int>{9,9,16,16}, {26,26,64,64},
                                 {26,9,64,16}, {9,26,16,64}, {26,26,125,125},
                                 {26,26,65,64}}) {
        const int nt=dimensions[0], ns=dimensions[1], qt=dimensions[2], qs=dimensions[3];
        std::vector<double> tx(3*qt), sx(3*qs), tw(qt), sw(qs), tv(qt*nt), sv(qs*ns);
        for (auto* data : {&tx,&sx,&tw,&sw,&tv,&sv})
            for (double& value : *data) value=random(generator);
        for (int i=0;i<qt;++i) tx[3*i] += 5.0;
        for (double& weight : tw) weight=std::abs(weight)/qt;
        for (double& weight : sw) weight=std::abs(weight)/qs;
        // Include a coincident point to preserve the scalar zero-radius guard.
        for (int axis=0;axis<3;++axis) tx[axis]=sx[axis];
        auto scalar = [&](double* output) {
            std::vector<double> inner(ns);
            for (int i=0;i<qt;++i) {
                std::fill(inner.begin(),inner.end(),0.0);
                for (int j=0;j<qs;++j) {
                    const double dx=tx[3*i]-sx[3*j], dy=tx[3*i+1]-sx[3*j+1], dz=tx[3*i+2]-sx[3*j+2];
                    const double radius=std::sqrt(dx*dx+dy*dy+dz*dz);
                    if (radius<1e-300) continue;
                    const double weight=sw[j]/radius;
                    for (int k=0;k<ns;++k) inner[k]+=weight*sv[j*ns+k];
                }
                for (int k=0;k<nt;++k) {
                    const double weight=tw[i]*tv[i*nt+k];
                    for (int j=0;j<ns;++j) output[k*ns+j]+=weight*inner[j];
                }
            }
        };
        auto candidate = [&](double* output) {
            radia_detail::HexFarProductBlas(nt,ns,qt,qs,tx.data(),tw.data(),tv.data(),
                                          sx.data(),sw.data(),sv.data(),output);
        };
        std::vector<double> reference(nt*ns), result(nt*ns);
        scalar(reference.data());
        const int previous=mkl_set_num_threads_local(2);
        candidate(result.data());
        const int restored=mkl_set_num_threads_local(previous);
        if(restored!=2) throw std::runtime_error("MKL local threads not restored");
        double scale=0.0,error=0.0;
        for (int i=0;i<nt*ns;++i) {
            scale=std::max(scale,std::abs(reference[i]));
            error=std::max(error,std::abs(reference[i]-result[i]));
        }
        if (!std::isfinite(error) || error>5e-13*scale) throw std::runtime_error("quadrature parity failed");
        auto measure = [&](auto&& kernel) {
            std::vector<double> samples;
            volatile double observed = 0.0;
            for (int repeat=0;repeat<7;++repeat) {
                const auto start=Clock::now();
                for (int k=0;k<200;++k) {
                    std::fill(result.begin(),result.end(),0.0);
                    kernel(result.data());
                    observed = result[k % result.size()];
                }
                samples.push_back(std::chrono::duration<double>(Clock::now()-start).count()/200);
            }
            std::sort(samples.begin(),samples.end());
            return samples[3];
        };
        const double before=measure(scalar),after=measure(candidate);
        if (!first) std::cout << ',';
        first=false;
        std::cout << "{\"nt\":" << nt << ",\"ns\":" << ns << ",\"qt\":" << qt
                  << ",\"qs\":" << qs << ",\"relative_max_error\":" << error/scale
                  << ",\"scalar_median_s\":" << before << ",\"blas_median_s\":" << after
                  << ",\"speedup\":" << before/after << '}';
    }
    std::cout << "]}\n";
    return 0;
}

#ifdef _WIN32
#define RADIA_BENCH_EXPORT __declspec(dllexport)
#else
#define RADIA_BENCH_EXPORT
#endif
extern "C" RADIA_BENCH_EXPORT int run_benchmark()
{
    try { return benchmark(); }
    catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
