#ifndef RAD_ARC_SECTION_H
#define RAD_ARC_SECTION_H

#include <array>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace RadArcSection {
using Vec = std::array<double, 3>;
template<class F> Vec Gauss4(const F& f, double lo, double hi);

inline Vec MomentSection(double phi, double r, double z, double ri, double ro, double h)
{
    const double c=std::cos(phi), s=std::sin(phi);
    const double rc=(ri+ro)/2, a=(ro-ri)/2, b=h/2;
    const double dx=r-rc*c, dy=-rc*s, dd=dx*dx+dy*dy+z*z;
    const double A=2*(rc-r*c)*a/dd, B=-2*z*b/dd, C=a*a/dd, D=b*b/dd;
    const double q=std::abs(A)+std::abs(B)+C+D;
    if(!(q<.125)) throw std::runtime_error("Arc moment expansion outside convergence region");
    // Expand (1+A*x+B*y+C*x*x+D*y*y)^(-3/2). All rectangular
    // moments are exact: mean(x^i*y^j) is zero for odd powers.
    constexpr int order=12, size=2*order+3;
    using Poly=std::array<std::array<double,size>,size>;
    Poly power{};
    power[0][0]=1.;
    auto moment=[](int i,int j) { return (i%2 || j%2) ? 0. : 1./((i+1.)*(j+1.)); };
    double coefficient=1., qpower=q, transverse=0., axial=0.;
    const double bound=(rc+a)*std::max(std::abs(z)+b,rc+a+std::abs(r*c));
    for(int n=0;n<=order;++n) {
        double t=0., v=0.;
        for(int i=0;i<=2*n;++i) for(int j=0;j<=2*n-i;++j) {
            const double p=power[i][j];
            if(p==0.) continue;
            t+=p*(rc*z*moment(i,j)+a*z*moment(i+1,j)
                  -rc*b*moment(i,j+1)-a*b*moment(i+1,j+1));
            v+=p*(rc*(rc-r*c)*moment(i,j)+a*(2*rc-r*c)*moment(i+1,j)
                  +a*a*moment(i+2,j));
        }
        transverse+=coefficient*t;
        axial+=coefficient*v;
        const double next_coefficient=-coefficient*(n+1.5)/(n+1.);
        const double tail=bound*std::abs(next_coefficient)*qpower/(1-q*(n+2.5)/(n+2.));
        if(tail<=1.e-12*std::max(std::abs(transverse),std::abs(axial))) {
            const double factor=4*a*b/(dd*std::sqrt(dd));
            return {factor*c*transverse,factor*s*transverse,factor*axial};
        }
        if(n==order) break;
        Poly next{};
        for(int i=0;i<=2*n;++i) for(int j=0;j<=2*n-i;++j) {
            const double p=power[i][j];
            next[i+1][j]+=A*p; next[i][j+1]+=B*p;
            next[i+2][j]+=C*p; next[i][j+2]+=D*p;
        }
        power=next;
        coefficient=next_coefficient;
        qpower*=q;
    }
    throw std::runtime_error("Arc moment expansion did not converge");
}

// Per unit azimuth, with u=r'-r*cos(phi), v=z'-z, q=|r*sin(phi)|,
// d^2=u^2+v^2+q^2 and a=r*cos(phi), the section integrands are
//   transverse T=-(u+a)*v/d^3 (times cos, sin)   axial X=(u+a)*u/d^3.
// Corner primitives difference nearly equal values when the section is
// thin compared with its distance, losing about eps*(distance/width)^2.
// Each direction is therefore judged separately: a direction thinner than
// 1/64 of the distance is integrated by 4-point Gauss (error about
// 128^-8), the other direction in closed form with rationalised
// differences. Everything is double on every platform.
constexpr double ThinSectionRatio = 1./64.;

namespace Detail {
// 1/d(p2)-1/d(p1), d(p)^2=p^2+rest, without subtracting nearly equal terms.
inline double DiffInv(double p1, double p2, double rest)
{
    const double d1=std::sqrt(p1*p1+rest), d2=std::sqrt(p2*p2+rest);
    return (p1-p2)*(p1+p2)/(d1*d2*(d1+d2));
}
// [p/d(p)] from p1 to p2, divided by rest.
inline double DiffRatioOver(double p1, double p2, double rest)
{
    const double d1=std::sqrt(p1*p1+rest), d2=std::sqrt(p2*p2+rest);
    if(p1*p2>0) return (p2-p1)*(p2+p1)/((p2*d1+p1*d2)*d1*d2);
    return (p2/d2-p1/d1)/rest;
}
// asinh(p2/rho)-asinh(p1/rho).
inline double DiffAsinh(double p1, double p2, double rho)
{
    const double d1=std::hypot(p1,rho), d2=std::hypot(p2,rho);
    if(p1>=0 && p2>=0) return std::log1p((p2-p1+d2-d1)/(p1+d1));
    if(p1<=0 && p2<=0) return std::log1p((p2-p1+d1-d2)/(d2-p2));
    return std::asinh(p2/rho)-std::asinh(p1/rho);
}
constexpr double GaussX[4]={-.8611363115940526,-.3399810435848563,
                            .3399810435848563,.8611363115940526};
constexpr double GaussW[4]={.3478548451374538,.6521451548625461,
                            .6521451548625461,.3478548451374538};
}

inline Vec Section(double phi, double r, double z, double ri, double ro, double h)
{
    using namespace Detail;
    const double c = std::cos(phi), s = std::sin(phi);
    const double a = r*c, q = std::abs(r*s);
    const double u1 = ri-a, u2 = ro-a, v1 = -h/2-z, v2 = h/2-z;
    const double du = u1>0 ? u1 : (u2<0 ? -u2 : 0.);
    const double dv = v1>0 ? v1 : (v2<0 ? -v2 : 0.);
    const double limit2 = ThinSectionRatio*ThinSectionRatio*(q*q+du*du+dv*dv);
    const bool thinU = (ro-ri)*(ro-ri) < limit2, thinV = h*h < limit2;
    // Thin widths come from the section data, not from endpoint differences.
    const double midU = (ri+ro)/2-a, halfU = (ro-ri)/2, midV = -z, halfV = h/2;
    double transverse = 0, axial = 0;
    if (thinU && thinV) {
        for (int i=0; i<4; ++i) for (int j=0; j<4; ++j) {
            const double u = midU+halfU*GaussX[i], v = midV+halfV*GaussX[j];
            const double dd = u*u+v*v+q*q, w = GaussW[i]*GaussW[j]/(dd*std::sqrt(dd));
            transverse -= w*(u+a)*v;
            axial += w*(u+a)*u;
        }
        transverse *= halfU*halfV; axial *= halfU*halfV;
    } else if (thinU) {
        for (int i=0; i<4; ++i) {
            const double u = midU+halfU*GaussX[i], rest = u*u+q*q;
            transverse += GaussW[i]*(u+a)*DiffInv(v1,v2,rest);
            axial += GaussW[i]*(u+a)*u*DiffRatioOver(v1,v2,rest);
        }
        transverse *= halfU; axial *= halfU;
    } else if (thinV) {
        for (int j=0; j<4; ++j) {
            const double v = midV+halfV*GaussX[j], rest = v*v+q*q;
            transverse += GaussW[j]*(-v)*(a*DiffRatioOver(u1,u2,rest)-DiffInv(u1,u2,rest));
            axial += GaussW[j]*(DiffAsinh(u1,u2,std::sqrt(rest))
                                -rest*DiffRatioOver(u1,u2,rest)-a*DiffInv(u1,u2,rest));
        }
        transverse *= halfV; axial *= halfV;
    } else {
        const double us[] = {u1, u2}, vs[] = {v1, v2};
        const double hu[] = {std::hypot(u1,q), std::hypot(u2,q)};
        const double hv[] = {std::hypot(v1,q), std::hypot(v2,q)};
        for (int i=0; i<2; ++i) for (int j=0; j<2; ++j) {
            const double u=us[i], v=vs[j];
            const double d=std::sqrt(u*u+v*v+q*q);
            const double uv=std::asinh(u/hv[j]);
            const double vu=std::asinh(v/hu[i]);
            const int sign = (i==j) ? 1 : -1;
            transverse += sign*(d + (a==0 ? 0 : a*uv));
            axial += sign*((v==0 ? 0 : v*uv) - (a==0 ? 0 : a*vu)
                           - (q==0 ? 0 : q*std::atan2(u*v,q*d)));
        }
    }
    return {c*transverse, s*transverse, axial};
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

// Fixed accuracy of the angular quadrature: relative 1e-9 (round-off level
// of the section kernel) plus an absolute part of 1e-12 x section size x
// interval.  It is not user-adjustable.
constexpr double RelTol = 1.e-9;
constexpr double AbsTolScale = 1.e-12;

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
    if(error<=atol+RelTol*magnitude) return fine;
    if(depth==0) throw std::runtime_error("Arc section integral did not converge");
    const Vec a=Refine(f,lo,mid,left,atol/2,depth-1);
    const Vec b=Refine(f,mid,hi,right,atol/2,depth-1);
    for(int k=0;k<3;++k) fine[k]=a[k]+b[k];
    return fine;
}

// Gauss-Legendre nodes and weights on [-1, 1] for 1..20 points, computed
// once by Newton iteration (thread-safe static initialisation).
struct LegendreRule { int n = 0; std::array<double, 20> x{}, w{}; };
inline const LegendreRule& Legendre(int n)
{
    static const std::array<LegendreRule, 21> rules = [] {
        std::array<LegendreRule, 21> table{};
        const double pi = 3.14159265358979323846;
        for (int m = 1; m <= 20; ++m) {
            table[m].n = m;
            for (int i = 0; i < m; ++i) {
                double t = std::cos(pi * (i + 0.75) / (m + 0.5)), derivative = 0.;
                for (int iteration = 0; iteration < 100; ++iteration) {
                    double p0 = 1., p1 = t;
                    for (int k = 2; k <= m; ++k) {
                        const double p2 = ((2 * k - 1) * t * p1 - (k - 1) * p0) / k;
                        p0 = p1; p1 = p2;
                    }
                    if (m == 1) { p1 = t; p0 = 1.; }
                    derivative = m * (t * p1 - p0) / (t * t - 1.);
                    const double step = p1 / derivative;
                    t -= step;
                    if (std::abs(step) < 1.e-16) break;
                }
                table[m].x[i] = t;
                table[m].w[i] = 2. / ((1. - t * t) * derivative * derivative);
            }
        }
        return table;
    }();
    return rules[n];
}

// Fixed-rule order for one azimuthal piece [lo, hi] of the section integral,
// or 0 when the adaptive rule must be kept.  The section is integrated in
// closed form, so only the azimuthal integrand's analyticity matters.  Its
// nearest complex singularity is at least the distance d from the observer
// to the arc solid; for a piece of half-length L (at the outer radius) the
// Bernstein-ellipse parameter is rho = t + sqrt(1 + t^2), t = d / L, and an
// n-point Gauss rule errs by O(rho^-2n).  n is chosen for 64 * rho^-2n below
// 1e-10, and points closer than L (n > 16) keep the adaptive rule.  The
// bound is relative to the integrand, not to an integral that cancels, so
// the caller also compares with n + 4 points and accepts only under the
// adaptive rule's own criterion (FarPiece).
inline int FarPieceOrder(double r, double z, double ri, double ro, double h,
                         double lo, double hi)
{
    const double pi = 3.14159265358979323846;
    const auto contains = [&](double angle) { return lo <= angle && angle <= hi; };
    const double c = (contains(0.) || contains(2 * pi) || contains(-2 * pi))
        ? 1. : std::max(std::cos(lo), std::cos(hi));
    // Closest point of the solid: nearest azimuth, then radius, then height.
    const double radius = c > 0. ? std::clamp(r * c, ri, ro) : ri;
    const double height = std::clamp(z, -h / 2, h / 2);
    const double d2 = r * r + radius * radius - 2 * r * radius * c
                      + (z - height) * (z - height);
    const double half = ro * (hi - lo) / 2;
    if (!(half > 0.) || !(d2 > 0.)) return 0;
    const double t = std::sqrt(d2) / half;
    const double rho = t + std::sqrt(1. + t * t);
    const int n = int(std::ceil(std::log(64. / 1.e-10) / (2. * std::log(rho))));
    return (n <= 16) ? std::max(n, 2) : 0;
}

template<class F> Vec GaussFixed(const F& f, double lo, double hi, int n)
{
    const LegendreRule& rule = Legendre(n);
    Vec result{};
    const double mid = (lo + hi) / 2, half = (hi - lo) / 2;
    for (int i = 0; i < n; ++i) {
        const Vec value = f(mid + half * rule.x[i]);
        for (int k = 0; k < 3; ++k) result[k] += half * rule.w[i] * value[k];
    }
    return result;
}

inline Vec AxisPrimitive(double v, double ri, double ro)
{
    // Radial antiderivative of the axial field and its first two z
    // derivatives. log1p and rationalized differences preserve thin widths.
    const double di=std::hypot(ri,v), d_o=std::hypot(ro,v);
    const double width=ro-ri, sum=ro+ri;
    const double logarithm=std::log1p(width*(1+sum/(di+d_o))/(ri+di));
    const double x=ri/di,y=ro/d_o;
    const double derivative_ratio=v*(ri-ro)*sum/((ri*d_o+ro*di)*di*d_o);
    return {v*logarithm,logarithm+v*derivative_ratio,
            derivative_ratio*(x*x+x*y+y*y)};
}

inline Vec FullCircleAxis(double r, double z, double ri, double ro, double h)
{
    if(std::abs(z)<=32*std::max(ro,h)) {
        const Vec lower=AxisPrimitive(z-h/2,ri,ro);
        const Vec upper=AxisPrimitive(z+h/2,ri,ro);
        const double two_pi=6.28318530717958647692;
        return {-two_pi*r*(upper[1]-lower[1])/2,0.,
                two_pi*(upper[0]-lower[0]-r*r*(upper[2]-lower[2])/4)};
    }
    // Axial source integration is exact. Radial integration gives B0 and
    // its z derivatives; the regular axis expansion avoids subtracting
    // nearly equal angular contributions to recover a tiny radial field.
    const double scale=std::max(ro,h);
    const double rho=r/scale, lower=(z-h/2)/scale, upper=(z+h/2)/scale;
    auto f=[=](double radius) {
        Vec value{};
        const double rr=radius*radius;
        if(std::abs(z)>32*scale) {
            // Rationalized endpoint differences, not axial quadrature.
            const double dl=std::hypot(radius,lower), du=std::hypot(radius,upper);
            const double gap=h/scale, sum=2*z/scale;
            const double delta=-gap*sum/(dl+du);
            const double l2=dl*dl,u2=du*du,l3=l2*dl,u3=u2*du;
            value[0]=rr*gap*sum/(du*dl*(upper*dl+lower*du));
            value[1]=rr*delta*(l2+dl*du+u2)/(l3*u3);
            const double diff5=delta*(l2*l2+l3*du+l2*u2+dl*u3+u2*u2)/(l3*l2*u3*u2);
            value[2]=-3*rr*(gap/(u3*u2)+lower*diff5);
            return value;
        }
        const double v[]={lower,upper};
        for(int i=0;i<2;++i) {
            const double dd=rr+v[i]*v[i], d=std::sqrt(dd);
            const double sign=i==0 ? -1. : 1.;
            value[0]+=sign*v[i]/d;
            value[1]+=sign*rr/(dd*d);
            value[2]-=sign*3*rr*v[i]/(dd*dd*d);
        }
        return value;
    };
    const double lo=ri/scale, hi=ro/scale;
    const Vec axis=Refine(f,lo,hi,Gauss4(f,lo,hi),1.e-13*(hi-lo),20);
    const double factor=6.28318530717958647692*scale;
    return {-factor*rho*axis[1]/2,0.,factor*(axis[0]-rho*rho*axis[2]/4)};
}

inline Vec Field(double r, double z, double ri, double ro, double h,
                 double lo, double hi)
{
    const double pi=3.14159265358979323846;
    const double span=hi-lo;
    lo=std::remainder(lo,2*pi);
    hi=lo+span;
    if(r==0.) {
        double transverse, axial;
        if(std::abs(z)>32*std::max(ro,h)) {
            const Vec section=MomentSection(0.,0.,z,ri,ro,h);
            transverse=section[0]; axial=section[2];
        } else {
            auto radial_distance=[=](double v) {
                return (ro-ri)*(ro+ri)/(std::hypot(ro,v)+std::hypot(ri,v));
            };
            transverse=radial_distance(z-h/2)-radial_distance(z+h/2);
            axial=AxisPrimitive(z+h/2,ri,ro)[0]-AxisPrimitive(z-h/2,ri,ro)[0];
        }
        return {transverse*(std::sin(hi)-std::sin(lo)),
                transverse*(std::cos(lo)-std::cos(hi)),axial*span};
    }
    Vec result{};
    // Split at the closest azimuth and opposite azimuth; never sample endpoints.
    while(lo<hi) {
        const double next=(std::floor(lo/pi)+1)*pi;
        const double end=std::min(hi,std::min(lo+pi/2,next));
        if(end<=lo) throw std::runtime_error("Arc integration interval collapsed");
        if(const int fixed = FarPieceOrder(r,z,ri,ro,h,lo,end)) {
            // Accept the fixed rule only under the adaptive rule's criterion,
            // estimated from the n and n + 4 point results; else refine below.
            auto smooth=[=](double phi){return Section(phi,r,z,ri,ro,h);};
            const Vec coarse=GaussFixed(smooth,lo,end,fixed);
            const Vec fine=GaussFixed(smooth,lo,end,fixed+4);
            double error=0, magnitude=0;
            bool finite=true;
            for(int k=0;k<3;++k) {
                finite=finite && std::isfinite(fine[k]) && std::isfinite(coarse[k]);
                error=std::max(error,std::abs(fine[k]-coarse[k]));
                magnitude=std::max(magnitude,std::abs(fine[k]));
            }
            if(finite && error<=AbsTolScale*std::max(ro,h)*(end-lo)+RelTol*magnitude) {
                for(int k=0;k<3;++k) result[k]+=fine[k];
                lo=end;
                continue;
            }
        }
        Vec value{};
        // Integrable logarithms and nearby exterior peaks occur at endpoints.
        // Map each half interval from its endpoint with phi = endpoint +/- L*t^4.
        const double half=(end-lo)/2;
        if(std::hypot(r,z)>32*std::max(ro,h)) {
            // The far moment kernel has no endpoint peak; retain its smooth rule.
            auto smooth=[=](double phi){return Section(phi,r,z,ri,ro,h);};
            value=Refine(smooth,lo,end,Gauss4(smooth,lo,end),
                         AbsTolScale*std::max(ro,h)*(end-lo),20);
        } else for(int side=0;side<2;++side) {
            const double endpoint=std::remainder(side==0 ? lo : end,2*pi);
            auto regular=[=](double t) {
                const double t3=t*t*t;
                Vec v=Section(endpoint+(side==0 ? 1 : -1)*half*t3*t,r,z,ri,ro,h);
                for(double& component:v) component*=4*half*t3;
                return v;
            };
            const Vec part=Refine(regular,0.,1.,Gauss4(regular,0.,1.),
                                  AbsTolScale*std::max(ro,h)*half,20);
            for(int k=0;k<3;++k) value[k]+=part[k];
        }
        for(int k=0;k<3;++k) result[k]+=value[k];
        lo=end;
    }
    return result;
}
}
#endif
