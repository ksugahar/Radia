// Standalone kernel test: only the C++ standard library is required.
#include "rad_arc_section.h"
#include <iostream>
#include <limits>
#include <vector>

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

// Direct volume reference for one azimuth: composite 20-point Gauss
// in radius and height, independent of Section's corner or thin routes.
Vec section_reference(double phi, double r, double z, double ri, double ro, double h) {
    const auto& g=RadArcSection::Legendre(20);
    const double c=std::cos(phi), s=std::sin(phi);
    constexpr int panels=8;
    Vec sum{};
    for(int pu=0;pu<panels;++pu) for(int pv=0;pv<panels;++pv)
    for(int i=0;i<20;++i) for(int j=0;j<20;++j) {
        const double hu=(ro-ri)/(2*panels), hv=h/(2*panels);
        const double radius=ri+(2*pu+1)*hu+g.x[i]*hu;
        const double dz=z-(-h/2+(2*pv+1)*hv+g.x[j]*hv);
        const double dx=r-radius*c, dy=-radius*s, d2=dx*dx+dy*dy+dz*dz;
        const double k=g.w[i]*g.w[j]*hu*hv*radius/(d2*std::sqrt(d2));
        sum[0]+=c*dz*k; sum[1]+=s*dz*k; sum[2]+=(radius-r*c)*k;
    }
    return sum;
}

// Sections thin in one or both directions, from far away down to a few
// widths, and both sides of the 1/64 thin switch. Corner differences alone
// lose up to 4e-4 on these (anisotropic control: radial width 1e-12 m).
bool thin_section_checks() {
    struct Case { double ri, ro, h, r, z, phi; };
    std::vector<Case> cases={
        {.100,.101,.001,1.92,-2.60,.7}, {.100,.101,.001,2.89,1.45,0.},
        {.100,.101,.001,.95,.31,1.3}, {.100,.101,.001,.20,.40,.2},
        {.100,.100000000001,.001,.12,.02,.1}, {.100,.100000000001,.001,.101,.0004,.3},
        {.100,.101,1e-9,.13,1e-5,.05}, {.100,.101,1e-9,.1005,.003,0.}};
    // Radial thinness switch at distance 64*(ro-ri) straight outward.
    for(double side: {1-1e-6, 1+1e-6})
        cases.push_back({.100,.1001,.01,.1001+64e-4*side,0.,0.});
    // Height thinness switch at distance 64*h straight above.
    for(double side: {1-1e-6, 1+1e-6})
        cases.push_back({.100,.110,1e-4,.105,5e-5+64e-4*side,0.});
    for(const Case& k: cases) {
        const Vec actual=RadArcSection::Section(k.phi,k.r,k.z,k.ri,k.ro,k.h);
        const double error=relative(actual,section_reference(k.phi,k.r,k.z,k.ri,k.ro,k.h));
        std::cout << "thin_section_error=" << error << '\n';
        if(!std::isfinite(error) || error>1e-12) return false;
    }
    return true;
}

// Azimuthal integral of a thin section whose closest azimuth uses the
// corner route and whose far azimuths use the thin route.
bool thin_switch_integral_check() {
    const double ri=.100, ro=.1001, h=.01, r=.103, z=0., lo=0., hi=pi/2;
    const auto& g=RadArcSection::Legendre(20);
    constexpr int panels=200;
    Vec expected{};
    for(int p=0;p<panels;++p) {
        const double half=(hi-lo)/(2*panels), mid=lo+(2*p+1)*half;
        for(int i=0;i<20;++i) {
            const Vec v=section_reference(mid+half*g.x[i],r,z,ri,ro,h);
            for(int k=0;k<3;++k) expected[k]+=half*g.w[i]*v[k];
        }
    }
    const Vec actual=RadArcSection::Field(r,z,ri,ro,h,lo,hi);
    const double error=relative(actual,expected);
    std::cout << "thin_switch_integral_error=" << error << '\n';
    return std::isfinite(error) && error<=1e-9;
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
    if(!thin_section_checks()) return 5;
    std::cout << "PASS: twelve thin and anisotropic sections match their volume integrals\n";
    if(!thin_switch_integral_check()) return 6;
    std::cout << "PASS: azimuthal integral across the thin switch\n";
    return 0;
}

int main() {
    try { return check_kernel(); }
    catch(const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 4;
    }
}
