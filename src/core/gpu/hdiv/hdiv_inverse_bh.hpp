// Shared inverse-BH scalar contract, in SI A/m; common mu0 omitted from energy.
#pragma once
#include <cstddef>
#include <vector>
#include <memory>
#include <string>
#include "hmatrix_host_fault.hpp"
namespace radia_hdiv_bh {
// Shared host payload ownership participates in per-module fail-N diagnostics.
// Python/NumPy and MATLAB publication allocations remain runtime-owned.
template<class T>using BoundaryVector=radia_hmatrix_cuda::HostVector<T>;
struct Law {
    BoundaryVector<double> knots,coefficients,integralBase; // row-major cubic rows
    double cap=0,barrierSlope=0,lastField=0,initialReluctivity=0;
};
struct Samples {BoundaryVector<double> secant,tangent,energy;std::size_t ownedDeviceBytes=0;};
struct MaterialInput {
    std::size_t points=0,dofs=0;
    BoundaryVector<std::size_t> offsets,indices; // CSR component-major E, 3*points rows
    BoundaryVector<double> values,weights;
};
// Affine-TET evaluation tables: NGSolve reference shapes per vertex-rank class
// (classes x rulePoints x local x 3), rule rows (xi0,xi1,xi2,weight), element
// vertices in NGSolve order, DOF and material-point numbers. No E values.
struct FlatEvaluationInput {
    std::size_t elements=0,local=0,rulePoints=0,classes=0,points=0,dofs=0;
    BoundaryVector<double> rule,shapes,vertices;
    BoundaryVector<int> elementClass,elementDofs,elementPoints;
};
struct LedgerRow {std::string payload,direction;std::size_t bytes=0;};
struct MaterialEvaluation {
    BoundaryVector<double> gradient,tangentAction;
    double energy=0,maximumMagnetization=0;std::size_t ownedDeviceBytes=0,linearizationEpoch=0;
};
struct MaterialStats {std::size_t points=0,dofs=0,devicePayloadBound=0,setupUploadBytes=0,stateUploadBytes=0,evaluations=0,epoch=1,linearizationEpoch=0,tangentApplications=0;bool closed=false,linearizationValid=false;};
class MaterialOperator {
public:
    virtual ~MaterialOperator()=default;
    virtual MaterialEvaluation evaluate(const radia_hmatrix_cuda::HostVector<double>&,const radia_hmatrix_cuda::HostVector<double>&)=0;
    virtual MaterialEvaluation linearize(const radia_hmatrix_cuda::HostVector<double>&)=0;
    virtual radia_hmatrix_cuda::HostVector<double> tangent(const radia_hmatrix_cuda::HostVector<double>&,std::size_t)=0;
    virtual MaterialStats stats()=0;
    // Setup transfers of this context; the route name says how E was built.
    virtual std::string route()=0;
    virtual radia_hmatrix_cuda::HostVector<LedgerRow> transfers()=0;
    virtual void close()=0;
};
std::shared_ptr<MaterialOperator> createMaterial(const MaterialInput&,const Law&);
// E and E^T are generated on the device; no evaluation matrix crosses.
std::shared_ptr<MaterialOperator> createFlatMaterial(const FlatEvaluationInput&,const Law&);
MaterialEvaluation evaluateMaterial(const MaterialInput&,const Law&,const radia_hmatrix_cuda::HostVector<double>& state,const radia_hmatrix_cuda::HostVector<double>& direction);
Samples evaluate(const Law&,const radia_hmatrix_cuda::HostVector<double>& amplitudes);
}
