#ifndef RAD_ARC_SECTION_H
#define RAD_ARC_SECTION_H

#include <array>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace RadArcSection {
using Vec = std::array<double, 3>;

// Integrate J e_phi x (x-x') / |x-x'|^3 analytically in radius and z.
// With u=r'-r*cos(phi), v=z'-z, q=|r*sin(phi)|, the corner
// primitives follow by integrating 1/sqrt(u*u+v*v+q*q).
inline Vec Section(double phi, double r, double z, double ri, double ro, double h)
{
    const long double c = std::cos(phi), s = std::sin(phi);
    const long double a = r*c, q = std::abs(r*s);
    const long double us[] = {ri-a, ro-a}, vs[] = {-h/2-z, h/2-z};
    long double transverse = 0, axial = 0;
    for (int i=0; i<2; ++i) for (int j=0; j<2; ++j) {
        const long double u=us[i], v=vs[j];
        const long double d=std::sqrt(u*u+v*v+q*q);
        const long double uv=std::asinh(u/std::hypot(v,q));
        const long double vu=std::asinh(v/std::hypot(u,q));
        const int sign = (i==j) ? 1 : -1;
        transverse += sign*(d + (a==0 ? 0 : a*uv));
        axial += sign*((v==0 ? 0 : v*uv) - (a==0 ? 0 : a*vu)
                       - (q==0 ? 0 : q*std::atan2(u*v,q*d)));
    }
    return {double(c*transverse), double(s*transverse), double(axial)};
}

template<class F> Vec Gauss4(const F& f, double lo, double hi)
{
    static const double x[] = {-.8611363115940526,-.3399810435848563,
                              .3399810435848563,.8611363115940526};
    static const double w[] = {.3478548451374538,.6521451548625461,
                              .6521451548625461,.3478548451374538};
    Vec result{};
    const double mid=(lo+hi)/2, half=(hi-lo)/2;
    for (int i=0;i<4;++i) {
        const Vec value=f(mid+half*x[i]);
        for (int k=0;k<3;++k) result[k]+=half*w[i]*value[k];
    }
    return result;
}

template<class F> Vec Refine(const F& f, double lo, double hi,
                            const Vec& coarse, double atol, int depth)
{
    const double mid=(lo+hi)/2;
    const Vec left=Gauss4(f,lo,mid), right=Gauss4(f,mid,hi);
    Vec fine{};
    double error=0, magnitude=0;
    for(int k=0;k<3;++k) {
        fine[k]=left[k]+right[k];
        if (!std::isfinite(fine[k]))
            throw std::runtime_error("Arc section integral is non-finite");
        error=std::max(error,std::abs(fine[k]-coarse[k]));
        magnitude=std::max(magnitude,std::abs(fine[k]));
    }
    if(error<=atol+1.e-9*magnitude) return fine;
    if(depth==0) throw std::runtime_error("Arc section integral did not converge");
    const Vec a=Refine(f,lo,mid,left,atol/2,depth-1);
    const Vec b=Refine(f,mid,hi,right,atol/2,depth-1);
    for(int k=0;k<3;++k) fine[k]=a[k]+b[k];
    return fine;
}

inline Vec Field(double r, double z, double ri, double ro, double h,
                 double lo, double hi)
{
    const double pi=3.14159265358979323846;
    const double span=hi-lo;
    lo=std::remainder(lo,2*pi);
    hi=lo+span;
    auto f=[=](double phi){return Section(phi,r,z,ri,ro,h);};
    Vec result{};
    // Split at the closest azimuth and opposite azimuth; never sample endpoints.
    while(lo<hi) {
        const double next=(std::floor(lo/pi)+1)*pi;
        const double end=std::min(hi,std::min(lo+pi/2,next));
        if(end<=lo) throw std::runtime_error("Arc integration interval collapsed");
        const Vec value=Refine(f,lo,end,Gauss4(f,lo,end),
                              1.e-12*std::max(ro,h)*(end-lo),20);
        for(int k=0;k<3;++k) result[k]+=value[k];
        lo=end;
    }
    return result;
}
}
#endif
