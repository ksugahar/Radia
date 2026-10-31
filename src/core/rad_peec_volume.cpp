// Independent implementation of the Laplace-transform factorization of 1/r.
// No singular filament substitution or half-space surface impedance.
#include "rad_peec_volume.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>

namespace radia {
namespace {
constexpr double pi=3.141592653589793238462643383279502884;
// Six-point Gaussian integral is stable where four primitive values cancel.
constexpr double xp[6]={-.9324695142031520278,-.6612093864662645137,-.2386191860831969086,
                         .2386191860831969086,.6612093864662645137,.9324695142031520278};
constexpr double wp[6]={.1713244923791703450,.3607615730481386076,.4679139345726910474,
                       .4679139345726910474,.3607615730481386076,.1713244923791703450};
double gaussian_pair(double lo1,double hi1,double lo2,double hi2,double t) {
    const double w1=hi1-lo1,w2=hi2-lo2,h1=w1/2,h2=w2/2;
    const double shift=(lo1-lo2)+(w1-w2)/2;
    if(t*(h1+h2)<.25) {
        double value=0;
        for(int i=0;i<6;++i)for(int j=0;j<6;++j) {
            const double z=t*(shift+h1*xp[i]-h2*xp[j]);
            value+=wp[i]*wp[j]*std::exp(-z*z);
        }
        return h1*h2*value;
    }
    // F''(s)=exp(-t^2 s^2). Split off |s| analytically so disjoint
    // far boxes do not subtract four large, nearly equal linear terms.
    auto g=[](double z){return std::exp(-z*z)-std::sqrt(pi)*z*std::erfc(z);};
    const double overlap=std::max(0.,std::min(hi1,hi2)-std::max(lo1,lo2));
    const double tail=g(t*std::abs(hi1-lo2))+g(t*std::abs(lo1-hi2))
                     -g(t*std::abs(hi1-hi2))-g(t*std::abs(lo1-lo2));
    const double value=std::sqrt(pi)*overlap/t+tail/(2*t*t);
    if(!std::isfinite(value) || value<0)
        throw std::runtime_error("Nonfinite/negative Gaussian volume pairing; no clipping/fallback");
    return value;
}
double area(const VolumeBox& b) {
    if(b.axis<0 || b.axis>2 || (b.sign!=1 && b.sign!=-1))
        throw std::invalid_argument("Explicit oriented Cartesian current axis required");
    double result=1;
    for(int d=0;d<3;++d) {
        const double width=b.hi[d]-b.lo[d];
        if(!std::isfinite(b.lo[d]) || !std::isfinite(b.hi[d]) || !std::isfinite(width) || width<=0)
            throw std::invalid_argument("Finite positive volume-box dimensions required");
        if(d!=b.axis)result*=width;
    }
    if(!std::isfinite(result) || result<=0)throw std::invalid_argument("Finite positive cross-section required");
    return result;
}
}
VolumeGreenRule::VolumeGreenRule(int n) {
    if(n<8 || n>256)throw std::invalid_argument("Explicit volume quadrature order must be 8..256");
    u_.resize(n);weight_.resize(n);
    for(int i=0;i<(n+1)/2;++i) {
        double z=std::cos(pi*(i+.75)/(n+.5)),p=0,derivative=0;
        bool done=false;
        for(int it=0;it<100;++it) {
            double a=1,b=0;
            for(int k=1;k<=n;++k){double c=b;b=a;a=((2*k-1)*z*b-(k-1)*c)/k;}
            p=a;derivative=n*(z*a-b)/(z*z-1);
            const double next=z-p/derivative;
            if(std::abs(next-z)<=4*std::numeric_limits<double>::epsilon()) {z=next;done=true;break;}
            z=next;
        }
        if(!done)throw std::runtime_error("Legendre rule construction failed");
        // Recompute derivative at the converged root.
        double a=1,b=0;
        for(int k=1;k<=n;++k){double c=b;b=a;a=((2*k-1)*z*b-(k-1)*c)/k;}
        derivative=n*(z*a-b)/(z*z-1);
        const double w=1/((1-z*z)*derivative*derivative);
        u_[i]=(1-z)/2;u_[n-1-i]=(1+z)/2;
        weight_[i]=w;weight_[n-1-i]=w;
    }
}
double VolumeGreenRule::Inductance(const VolumeBox& a,const VolumeBox& b) const {
    const double aa=area(a),ab=area(b);
    if(a.axis!=b.axis)return 0.; // J_i dot J_j is exactly zero.
    double scale=0.;
    for(int d=0;d<3;++d) {
        scale=std::max(scale,std::max(a.hi[d],b.hi[d])-std::min(a.lo[d],b.lo[d]));
    }
    if(!std::isfinite(scale) || scale<=0)throw std::invalid_argument("Finite box-pair extent required");
    double sum=0.,correction=0.;
    for(size_t k=0;k<u_.size();++k) {
        const double v=1-u_[k],t=u_[k]/(scale*v);
        double product=1.;
        for(int d=0;d<3;++d)product*=gaussian_pair(a.lo[d],a.hi[d],b.lo[d],b.hi[d],t);
        const double term=weight_[k]*product/(scale*v*v);
        const double adjusted=term-correction,next=sum+adjusted;
        correction=(next-sum)-adjusted;sum=next;
    }
    const double result=1.e-7*(2/std::sqrt(pi))*sum/(aa*ab)*a.sign*b.sign;
    if(!std::isfinite(result))throw std::runtime_error("Nonfinite volume inductance");
    return result;
}
}
