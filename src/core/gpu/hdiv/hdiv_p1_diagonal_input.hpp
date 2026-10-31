// Host-only affine-TET P1 charge contract, not an HDiv operator or FE map.
#pragma once
#include <array>
#include <algorithm>
#include <cmath>
#include "hmatrix_host_fault.hpp"

namespace radia_hdiv_python {
struct P1Charge {
    std::array<double,12> vertices{}; // triangle uses the first nine entries
    std::array<int,3> exponent{};
    std::size_t first=0,count=0;
    int kind=0; // 0: affine tetrahedron; 1: flat triangle
};
struct PhysicalChargePoint {double point[3];double weight;};
// POD consumed by the shared host/device self-entry arithmetic. The polynomial
// is A + gradient.(x-center); physical weights do NOT include this monomial.
struct P1ChargeSource {
    double vertices[12]{},center[3]{},gradient[3]{},constant=1.;
    std::size_t first=0,count=0;
    int kind=0;
    bool linear=false;
};
inline P1ChargeSource p1ChargeSource(const P1Charge& charge){
    P1ChargeSource out;
    std::copy(charge.vertices.begin(),charge.vertices.end(),out.vertices);
    out.kind=charge.kind;out.first=charge.first;out.count=charge.count;
    out.linear=charge.exponent[0]!=0||charge.exponent[1]!=0;
    if(!out.linear)return out;
    // Radia's charge monomial frame, not an NGSolve FE/orientation map.
    // Same centred linear polynomial as CPU InitHOPolynomialCoefficients.
    double e1[3],e2[3],g00=0,g01=0,g11=0;
    for(int k=0;k<3;++k){
        out.center[k]=(charge.vertices[k]+charge.vertices[3+k]+charge.vertices[6+k])/3.;
        e1[k]=charge.vertices[3+k]-charge.vertices[k];
        e2[k]=charge.vertices[6+k]-charge.vertices[k];
        g00+=e1[k]*e1[k];g01+=e1[k]*e2[k];g11+=e2[k]*e2[k];
    }
    const double determinant=g00*g11-g01*g01;
    if(!std::isfinite(determinant)||determinant<=0.)
        throw std::invalid_argument("P1 face monomial frame is singular or unrepresentable");
    const bool first=charge.exponent[0]==1;
    const double a=(first?g11:-g01)/determinant;
    const double b=(first?-g01:g00)/determinant;
    out.constant=0.;
    for(int k=0;k<3;++k){
        out.gradient[k]=a*e1[k]+b*e2[k];
        out.constant-=out.gradient[k]*(charge.vertices[k]-out.center[k]);
        if(!std::isfinite(out.center[k])||!std::isfinite(out.gradient[k]))
            throw std::invalid_argument("Nonfinite P1 face polynomial");
    }
    if(!std::isfinite(out.constant))throw std::invalid_argument("Nonfinite P1 face polynomial");
    return out;
}
struct P1DiagonalInput {
    int order=1;
    radia_hmatrix_cuda::HostVector<P1Charge> charges;
    radia_hmatrix_cuda::HostVector<PhysicalChargePoint> quadrature;

    void validate()const {
        if(order!=1||charges.empty()||charges.size()>4608)
            throw std::invalid_argument("P1 charge diagonal requires order 1 within bounded scope");
        std::size_t next=0;
        for(const auto& charge:charges){
            if(charge.kind!=0&&charge.kind!=1)
                throw std::invalid_argument("P1 charge diagonal requires affine TET or flat TRI supports");
            const auto& e=charge.exponent;
            if(e[0]<0||e[0]>1||e[1]<0||e[1]>1||e[2]!=0||e[0]+e[1]>(charge.kind?1:0))
                throw std::invalid_argument("P1 charge diagonal rejects unsupported monomials");
            if(charge.first!=next||!charge.count||charge.count>256||
               charge.first>quadrature.size()||charge.count>quadrature.size()-charge.first)
                throw std::invalid_argument("Invalid P1 charge quadrature partition");
            next+=charge.count;
            for(auto value:charge.vertices)if(!std::isfinite(value))
                throw std::invalid_argument("Nonfinite P1 charge geometry");
            if(charge.kind==1)for(int k=9;k<12;++k)if(charge.vertices[k]!=0.)
                throw std::invalid_argument("Triangle vertex padding must be zero");
            double edges[3][3]{};
            for(int i=0;i<(charge.kind?2:3);++i)for(int k=0;k<3;++k)
                edges[i][k]=charge.vertices[3*(i+1)+k]-charge.vertices[k];
            const double cross[3]={edges[0][1]*edges[1][2]-edges[0][2]*edges[1][1],
                                   edges[0][2]*edges[1][0]-edges[0][0]*edges[1][2],
                                   edges[0][0]*edges[1][1]-edges[0][1]*edges[1][0]};
            const double measure=charge.kind?
                std::hypot(cross[0],cross[1],cross[2]):
                std::fabs(cross[0]*edges[2][0]+cross[1]*edges[2][1]+cross[2]*edges[2][2]);
            if(!std::isfinite(measure)||measure<=0.)
                throw std::invalid_argument("Degenerate P1 charge support");
            (void)p1ChargeSource(charge); // Reject unrepresentable polynomial before CUDA.
        }
        if(next!=quadrature.size())throw std::invalid_argument("Unclaimed charge quadrature points");
        for(const auto& point:quadrature){
            for(auto value:point.point)if(!std::isfinite(value))
                throw std::invalid_argument("Nonfinite physical charge quadrature point");
            if(!std::isfinite(point.weight))
                throw std::invalid_argument("Nonfinite physical charge quadrature weight");
        }
    }
};
struct P1DiagonalStats {
    std::size_t charges=0,ownedDeviceBytes=0,setupPeakOwnedBytes=0;
    std::size_t setupUploadBytes=0,setupDownloadBytes=0,diagnosticDownloadBytes=0;
    bool closed=false;
};
class P1Diagonal {
public:
    virtual ~P1Diagonal()=default;
    // Explicit diagnostic downloads only; production construction stays resident.
    virtual radia_hmatrix_cuda::HostVector<double> diagnostic_diagonal()=0;
    virtual radia_hmatrix_cuda::HostVector<double> diagnostic_scale()=0;
    virtual P1DiagonalStats stats()=0;
    virtual void close()=0;
};
std::shared_ptr<P1Diagonal> create_p1_diagonal(const P1DiagonalInput&);
} // namespace radia_hdiv_python
