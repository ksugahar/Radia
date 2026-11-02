// MEX-free ownership boundary for the shared B^T G B resident CUDA operator.
#pragma once
#include <cstddef>
#include <memory>
#include <string>
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
// Primitive flat-TET charge geometry; reference rules are (xi0,xi1,xi2,weight) rows.
struct FlatInput {
    radia_hmatrix_cuda::HostVector<double> cells,faces,tet,tri,tetLow,triLow;
    radia_hmatrix_cuda::HostVector<int> host,kind,expo;
    double farFactor=1e30;
};
// NGSolve charge-projection tables per vertex-rank class: B rows are generated
// on the device from these and the FlatInput vertices; no B value crosses.
// Exponents are host-side layout checks against FlatInput and never uploaded.
struct ChargeMapInput {
    std::size_t dofs=0,cellMonomials=0,faceMonomials=0,cellWidth=0,faceWidth=0;
    radia_hmatrix_cuda::HostVector<double> cellTables,faceTables;
    radia_hmatrix_cuda::HostVector<int> cellClass,cellDofs,faceClass,faceDofs,cellExponents,faceExponents;
};
struct Transfer {std::string payload,direction;std::size_t bytes;};
struct Assembly {std::string route;std::size_t charges,setupBytes,launches;radia_hmatrix_cuda::HostVector<Transfer> transfers;};
class DeviceOperator {
public:
    virtual ~DeviceOperator()=default;
    virtual radia_hmatrix_cuda::HostVector<double> apply(const radia_hmatrix_cuda::HostVector<double>&)=0;
    virtual MaterialResult solve_linear(const MapInput&,const radia_hmatrix_cuda::HostVector<double>&,double,double,int,bool,bool,int,double)=0;
    virtual NewtonResult solve_nonlinear(radia_hdiv_bh::MaterialOperator&,const radia_hmatrix_cuda::HostVector<double>&,const radia_hmatrix_cuda::HostVector<double>&,const NewtonOptions&)=0;
    // Host-evaluation contexts only (validation consumers until they migrate to device tables).
    virtual NewtonResult solve_nonlinear_host_reference(radia_hdiv_bh::MaterialOperator&,const MapInput&,const radia_hmatrix_cuda::HostVector<double>&,const radia_hmatrix_cuda::HostVector<double>&,const NewtonOptions&)=0;
    virtual Stats stats()=0;
    virtual Assembly assembly()=0;
    virtual void close()=0;
};
std::shared_ptr<DeviceOperator> create(std::size_t,const radia_hmatrix_cuda::HostVector<LeafInput>&,
    const radia_hmatrix_cuda::HostVector<std::size_t>&,const radia_hmatrix_cuda::HostVector<double>&,const MapInput&);
// The charge Gram, B and B^T are generated on the device; no Gram or B value crosses.
std::shared_ptr<DeviceOperator> createFlat(const FlatInput&,const ChargeMapInput&);
} // namespace radia_hdiv_python
