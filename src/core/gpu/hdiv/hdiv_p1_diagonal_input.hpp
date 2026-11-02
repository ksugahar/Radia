// Affine-TET P1 charge source shared by the CPU ChargeGram and the device flat
// Gram graded touching rule; not an HDiv operator or FE map.
#pragma once
#include <array>
#include <cfloat>
#include <cmath>
#include <stdexcept>
#ifdef __CUDACC__
#define RADIA_P1_HOST_DEVICE __host__ __device__
#else
#define RADIA_P1_HOST_DEVICE
#endif

namespace radia_hdiv_python {
struct P1Charge {
    std::array<double,12> vertices{}; // triangle uses the first nine entries
    std::array<int,3> exponent{};
    int kind=0; // 0: affine tetrahedron; 1: flat triangle
};
// POD consumed by the shared host/device pair-entry arithmetic. The polynomial
// is constant + gradient.(x-center).
struct P1ChargeSource {
    double vertices[12]{},center[3]{},gradient[3]{},constant=1.;
    int kind=0;
    bool linear=false;
};
// Radia's charge monomial frame, not an NGSolve FE/orientation map: the same
// centred linear polynomial as CPU InitHOPolynomialCoefficients. Returns false
// for a singular face frame; the caller checks finiteness. Every field is
// written: the device fills uninitialized buffer slots in place.
RADIA_P1_HOST_DEVICE inline bool p1ChargeSourceInto(int kind,const double* vertices,const int* exponent,
        P1ChargeSource& out){
    out=P1ChargeSource{};
    for(int i=0;i<12;++i)out.vertices[i]=vertices[i];
    out.kind=kind;
    out.linear=exponent[0]!=0||exponent[1]!=0;
    if(!out.linear)return true;
    double e1[3],e2[3],g00=0,g01=0,g11=0;
    for(int k=0;k<3;++k){
        out.center[k]=(vertices[k]+vertices[3+k]+vertices[6+k])/3.;
        e1[k]=vertices[3+k]-vertices[k];
        e2[k]=vertices[6+k]-vertices[k];
        g00+=e1[k]*e1[k];g01+=e1[k]*e2[k];g11+=e2[k]*e2[k];
    }
    const double determinant=g00*g11-g01*g01;
    if(!(determinant>0.&&determinant<=DBL_MAX))return false;
    const bool first=exponent[0]==1;
    const double a=(first?g11:-g01)/determinant;
    const double b=(first?-g01:g00)/determinant;
    out.constant=0.;
    for(int k=0;k<3;++k){
        out.gradient[k]=a*e1[k]+b*e2[k];
        out.constant-=out.gradient[k]*(vertices[k]-out.center[k]);
    }
    return true;
}
inline P1ChargeSource p1ChargeSource(const P1Charge& charge){
    P1ChargeSource out;
    if(!p1ChargeSourceInto(charge.kind,charge.vertices.data(),charge.exponent.data(),out))
        throw std::invalid_argument("P1 face monomial frame is singular or unrepresentable");
    for(int k=0;k<3;++k)
        if(!std::isfinite(out.center[k])||!std::isfinite(out.gradient[k]))
            throw std::invalid_argument("Nonfinite P1 face polynomial");
    if(!std::isfinite(out.constant))throw std::invalid_argument("Nonfinite P1 face polynomial");
    return out;
}
} // namespace radia_hdiv_python
