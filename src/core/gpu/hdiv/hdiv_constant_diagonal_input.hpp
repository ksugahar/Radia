// Host-only contract for the constant-charge component, not an HDiv operator.
#pragma once
#include <array>
#include <cmath>
#include "hmatrix_host_fault.hpp"

namespace radia_hdiv_python {
struct ConstantCharge {
    std::array<double,12> vertices{}; // triangle uses the first nine entries
    std::array<int,3> exponent{};
    std::size_t first=0,count=0;
    int kind=0; // 0: affine tetrahedron; 1: flat triangle
};
struct PhysicalChargePoint {double point[3];double weight;};
struct ConstantDiagonalInput {
    radia_hmatrix_cuda::HostVector<ConstantCharge> charges;
    radia_hmatrix_cuda::HostVector<PhysicalChargePoint> quadrature;

    void validate()const {
        if(charges.empty()||charges.size()>4608)
            throw std::invalid_argument("Constant charge diagonal dimensions exceed bounded scope");
        std::size_t next=0;
        for(const auto& charge:charges){
            if(charge.kind!=0&&charge.kind!=1)
                throw std::invalid_argument("Constant charge diagonal requires TET or TRI supports");
            for(auto power:charge.exponent)if(power!=0)
                throw std::invalid_argument("Constant charge diagonal rejects nonconstant monomials");
            if(charge.first!=next||!charge.count||charge.count>256||
               charge.first>quadrature.size()||charge.count>quadrature.size()-charge.first)
                throw std::invalid_argument("Invalid constant charge quadrature partition");
            next+=charge.count;
            for(auto value:charge.vertices)if(!std::isfinite(value))
                throw std::invalid_argument("Nonfinite constant charge geometry");
            if(charge.kind==1)for(int k=9;k<12;++k)if(charge.vertices[k]!=0.)
                throw std::invalid_argument("Triangle vertex padding must be zero");
            double e[3][3]{};
            for(int i=0;i<(charge.kind?2:3);++i)for(int k=0;k<3;++k)
                e[i][k]=charge.vertices[3*(i+1)+k]-charge.vertices[k];
            const double cross[3]={e[0][1]*e[1][2]-e[0][2]*e[1][1],
                                   e[0][2]*e[1][0]-e[0][0]*e[1][2],
                                   e[0][0]*e[1][1]-e[0][1]*e[1][0]};
            const double measure=charge.kind?
                std::hypot(cross[0],cross[1],cross[2]):
                std::fabs(cross[0]*e[2][0]+cross[1]*e[2][1]+cross[2]*e[2][2]);
            if(!std::isfinite(measure)||measure<=0.)
                throw std::invalid_argument("Degenerate constant charge support");
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
struct ConstantDiagonalStats {
    std::size_t charges=0,ownedDeviceBytes=0,setupPeakOwnedBytes=0;
    std::size_t setupUploadBytes=0,setupDownloadBytes=0,diagnosticDownloadBytes=0;
    bool closed=false;
};
class ConstantDiagonal {
public:
    virtual ~ConstantDiagonal()=default;
    // Explicit diagnostic downloads only; production construction stays resident.
    virtual radia_hmatrix_cuda::HostVector<double> diagnostic_diagonal()=0;
    virtual radia_hmatrix_cuda::HostVector<double> diagnostic_scale()=0;
    virtual ConstantDiagonalStats stats()=0;
    virtual void close()=0;
};
std::shared_ptr<ConstantDiagonal> create_constant_diagonal(const ConstantDiagonalInput&);
} // namespace radia_hdiv_python
