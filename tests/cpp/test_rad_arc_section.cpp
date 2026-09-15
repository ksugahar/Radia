// Standalone kernel test: only the C++ standard library is required.
#include "rad_arc_section.h"
#include <iostream>

using RadArcSection::Vec;
constexpr double pi = 3.14159265358979323846;

Vec midpoint(int n, double r, double z, double lo, double hi) {
    const double ri=.035, ro=.070, h=.105;
    Vec sum{};
    for(int k=0;k<n;++k) {
        const double phi=lo+(k+.5)*(hi-lo)/n;
        const double c=std::cos(phi), s=std::sin(phi);
        for(int i=0;i<n;++i) {
            const double radius=ri+(i+.5)*(ro-ri)/n;
            for(int j=0;j<n;++j) {
                const double dz=z-(-h/2+(j+.5)*h/n);
                const double dx=r-radius*c, dy=-radius*s;
                const double d2=dx*dx+dy*dy+dz*dz;
                const double w=radius/(d2*std::sqrt(d2));
                sum[0]+=c*dz*w; sum[1]+=s*dz*w;
                sum[2]+=(radius-r*c)*w;
            }
        }
    }
    for(double& v:sum) v*=(ro-ri)*h*(hi-lo)/(double(n)*n*n);
    return sum;
}

double relative(const Vec& a,const Vec& b) {
    double error=0, norm=0;
    for(int k=0;k<3;++k) { error+=(a[k]-b[k])*(a[k]-b[k]); norm+=b[k]*b[k]; }
    return std::sqrt(error/norm);
}

int main() {
    const double half=.105/2;
    const double exact=2*pi*.105*std::log((.070+std::hypot(.070,half))/(.035+std::hypot(.035,half)));
    const Vec axis=RadArcSection::Field(0,0,.035,.070,.105,0,2*pi);
    if(relative(axis,{0,0,exact})>1e-12) return 1;
    for(const auto& c: {std::array<double,4>{.012,.09,.2,4.8},
                       std::array<double,4>{.012,.09,0,2*pi},
                       std::array<double,4>{.012,.09,0,1e-5},
                       std::array<double,4>{40.,80.,.2,4.8}}) {
        const auto a=midpoint(64,c[0],c[1],c[2],c[3]);
        const auto b=midpoint(128,c[0],c[1],c[2],c[3]);
        const auto d=midpoint(256,c[0],c[1],c[2],c[3]);
        Vec coarse{}, fine{};
        for(int k=0;k<3;++k) { coarse[k]=(4*b[k]-a[k])/3; fine[k]=(4*d[k]-b[k])/3; }
        const auto actual=RadArcSection::Field(c[0],c[1],.035,.070,.105,c[2],c[3]);
        const double convergence=relative(coarse,fine), error=relative(actual,fine);
        std::cout << "reference_convergence=" << convergence << " kernel_error=" << error << '\n';
        if(!std::isfinite(error) || convergence>2e-7 || error>2e-7) return 2;
    }
    std::cout << "PASS: axis closed form and four independent volume references\n";
}
