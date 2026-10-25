// Analytic, provider-independent FFT validation. No stored Radia output oracle.
#include "fftw.h"
#include <algorithm>
#include <utility>
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>
#include <limits>
#include <stdexcept>
static double check(int nx,int ny,int sign,bool inplace,int batches=1) {
    const int n=nx*ny;
    std::vector<fftw_complex> data(n*batches),out(n*batches);
    for(int b=0;b<batches;++b) for(int i=0;i<n;++i) data[b*n+i]={float((i+2*b)%5-2),float((2*i+b)%7-3)};
    const auto original=data;
    if(nx==1) {
        auto plan=fftw_create_plan(ny,static_cast<fftw_direction>(sign),inplace?FFTW_IN_PLACE:FFTW_OUT_OF_PLACE);
        if(!plan) throw std::runtime_error("FFT plan allocation");
        fftw(plan,batches,data.data(),1,n,inplace?data.data():out.data(),1,n);
        fftw_destroy_plan(plan);
    } else {
        auto plan=fftw2d_create_plan(nx,ny,static_cast<fftw_direction>(sign),inplace?FFTW_IN_PLACE:FFTW_OUT_OF_PLACE);
        if(!plan) throw std::runtime_error("FFT plan allocation");
        fftwnd(plan,batches,data.data(),1,n,inplace?data.data():out.data(),1,n);
        fftwnd_destroy_plan(plan);
    }
    const auto& result=inplace?data:out;
    double worst=0,scale=0;
    for(int b=0;b<batches;++b) for(int kx=0;kx<nx;++kx) for(int ky=0;ky<ny;++ky) {
        std::complex<double> expected=0.;
        for(int x=0;x<nx;++x) for(int y=0;y<ny;++y) {
            auto v=original[b*n+x*ny+y];
            double angle=sign*2*std::acos(-1.)*(double(kx*x)/nx+double(ky*y)/ny);
            expected+=std::complex<double>(v.re,v.im)*std::polar(1.,angle);
        }
        auto got=result[b*n+kx*ny+ky];
        worst=std::max(worst,std::abs(std::complex<double>(got.re,got.im)-expected));
        scale=std::max(scale,std::abs(expected));
    }
    // A conservative analytical float roundoff envelope, not a changed solver gate.
    const double bound=32*n*std::numeric_limits<float>::epsilon()*std::max(1.,scale);
    if(!std::isfinite(worst) || worst>bound) throw std::runtime_error("FFT differs from independent analytic DFT");
    return worst;
}
int main() {
    double worst=0;
    for(auto dims: {std::pair<int,int>{1,8},{2,3},{3,4}})
        for(int sign:{-1,1}) for(bool inplace:{false,true})
            worst=std::max(worst,check(dims.first,dims.second,sign,inplace,2));
    std::printf("{\"analytic_dft\":\"PASS\",\"max_absolute_error\":%.17g}\n",worst);
}
