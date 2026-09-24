// Standalone kernel test: only the C++ standard library is required.
#include "rad_arc_section.h"
#include <iostream>
#include <limits>

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

bool boundary_checks() {
    for(const auto& p: {std::array<double,2>{.035,0}, {.070,0},
                       {.035,.0525}, {.070,.0525}, {.050,.0525}}) {
        const auto full=RadArcSection::Field(p[0],p[1],.035,.070,.105,0,2*pi);
        const auto first=RadArcSection::Field(p[0],p[1],.035,.070,.105,0,1.7);
        const auto second=RadArcSection::Field(p[0],p[1],.035,.070,.105,1.7,2*pi);
        Vec sum{};
        for(int k=0;k<3;++k) sum[k]=first[k]+second[k];
        const double partition_error=relative(sum,full);
        if(!std::isfinite(partition_error) || partition_error>1e-8) return false;
        for(double angle: {1.7,2*pi}) {
            const auto actual=RadArcSection::Field(p[0],p[1],.035,.070,.105,0,angle);
            const auto left=RadArcSection::Field(p[0]-1e-8,p[1]-1e-8,.035,.070,.105,0,angle);
            const auto right=RadArcSection::Field(p[0]+1e-8,p[1]+1e-8,.035,.070,.105,0,angle);
            Vec average{};
            for(int k=0;k<3;++k) average[k]=(left[k]+right[k])/2;
            const double continuity_error=relative(actual,average);
            if(!std::isfinite(continuity_error) || continuity_error>3e-5) return false;
        }
    }
    return true;
}

int check_kernel() {
    const double half=.105/2;
    const double exact=2*pi*.105*std::log((.070+std::hypot(.070,half))/(.035+std::hypot(.035,half)));
    const Vec axis=RadArcSection::Field(0,0,.035,.070,.105,0,2*pi);
    const double axis_error=relative(axis,{0,0,exact});
    if(!std::isfinite(axis_error) || axis_error>1e-12) return 1;
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
    if(!boundary_checks()) return 3;
    std::cout << "PASS: five boundary partitions and ten two-sided limits\n";
    // Explicit default tolerance is bit-identical; out-of-range values throw.
    const auto implicit=RadArcSection::Field(.012,.09,.035,.070,.105,.2,4.8);
    const auto explicit_default=RadArcSection::Field(.012,.09,.035,.070,.105,.2,4.8,
                                                     RadArcSection::DefaultRelTol);
    if(implicit!=explicit_default) return 5;
    for(const double bad: {0., 1.e-13, 1.e-2, std::numeric_limits<double>::quiet_NaN()}) {
        bool threw=false;
        try { RadArcSection::Field(.012,.09,.035,.070,.105,.2,4.8,bad); }
        catch(const std::invalid_argument&) { threw=true; }
        if(!threw) return 6;
    }
    std::cout << "PASS: default tolerance identity and range checks\n";
    return 0;
}

int main() {
    try { return check_kernel(); }
    catch(const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 4;
    }
}
