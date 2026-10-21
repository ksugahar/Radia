// MEX-free ownership boundary for the shared B^T G B resident CUDA operator.
#pragma once
#include <cstddef>
#include <memory>
#include <vector>
#include "hdiv_inverse_bh.hpp"
namespace radia_hdiv_python {
struct LeafInput {
    std::size_t row,column,rows,columns,rank;
    radia_hmatrix_cuda::HostVector<double> a1,a2;
    bool deferred=false;
};
struct MapInput {
    std::size_t rows,columns;
    radia_hmatrix_cuda::HostVector<std::size_t> offsets,indices;
    radia_hmatrix_cuda::HostVector<double> values;
};
struct MaterialResult {radia_hmatrix_cuda::HostVector<double> x;int iterations;double relative_residual;std::size_t owned_bytes,graph_launches;};
struct NewtonOptions {double tolerance=1e-11,linearTolerance=5e-14;int maximumIterations=40,maximumLinearIterations=10000,preconditioner=0,degree=32;double lowerFraction=.01;};
struct NewtonStep {double residual=0,energy=0,alpha=0,linearResidual=0;int linearIterations=0,backtracks=0;bool residualCertified=false;};
struct NewtonResult {radia_hdiv_bh::BoundaryVector<double> x;radia_hdiv_bh::BoundaryVector<NewtonStep> history;double residual=0,energy=0;std::size_t graphLaunches=0;};
struct Stats {std::size_t dimension,owned_bytes,applications;bool closed;};
class DeviceOperator {
public:
    virtual ~DeviceOperator()=default;
    virtual radia_hmatrix_cuda::HostVector<double> apply(const radia_hmatrix_cuda::HostVector<double>&)=0;
    virtual MaterialResult solve_linear(const MapInput&,const radia_hmatrix_cuda::HostVector<double>&,double,double,int,bool,bool,int,double)=0;
    virtual NewtonResult solve_nonlinear(radia_hdiv_bh::MaterialOperator&,const MapInput&,const radia_hmatrix_cuda::HostVector<double>&,const radia_hmatrix_cuda::HostVector<double>&,const NewtonOptions&)=0;
    virtual Stats stats()=0;
    virtual void close()=0;
};
std::shared_ptr<DeviceOperator> create(std::size_t,const radia_hmatrix_cuda::HostVector<LeafInput>&,
    const radia_hmatrix_cuda::HostVector<std::size_t>&,const radia_hmatrix_cuda::HostVector<double>&,const MapInput&);
} // namespace radia_hdiv_python
