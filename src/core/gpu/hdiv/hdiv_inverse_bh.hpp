// Shared inverse-BH scalar contract, in SI A/m; common mu0 omitted from energy.
#pragma once
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <stdexcept>
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
// TET evaluation tables: NGSolve reference shapes per vertex-rank class
// (classes x rulePoints x local x 3), rule rows (xi0,xi1,xi2,weight), element
// geometry, DOF and material-point numbers. No E values. Geometry holds
// `nodes` points per element: 4 affine vertices in NGSolve order, the 10
// P2 nodes of a curved element (corners first, NGSolve reference frame), or
// the 27-node Q2 lattice of a HEX (node ix + 3 iy + 9 iz at xi = (ix,iy,iz)/2),
// or the 5 vertices of a straight pyramid in NGSolve order.
struct FlatEvaluationInput {
    std::size_t elements=0,local=0,rulePoints=0,classes=0,points=0,dofs=0,nodes=4;
    BoundaryVector<double> rule,shapes,vertices;
    BoundaryVector<int> elementClass,elementDofs,elementPoints;
};
// One table set for a pure mesh, one per element type of a mixed mesh. DOF and
// point numbers stay global: `points` of every group is the mesh total, and the
// groups together cover it once.
using FlatEvaluationGroups=BoundaryVector<FlatEvaluationInput>;
inline std::string flatEvaluationFamily(const FlatEvaluationGroups& in){
    return in.size()>1||in.at(0).nodes==5?"mixed":in[0].nodes==27?"hex":in[0].nodes==10?"curved-tet":"flat-tet";
}
// Device bytes of one matrix-free flat material build and its resident state,
// from table dimensions alone: no global E is stored. Builder transients are
// counted as live together, so the sum bounds the setup peak; CUB scratch is
// queried on the device and added there. Mass terms are upper bounds (one
// stored entry per element triple).
struct FlatEvaluationBudget {
    static constexpr std::size_t count=6;
    static constexpr const char* names[count]={"reference_tables","incidence","incidence_sort",
                                               "mass_triples","mass_csr","material_state"};
    std::size_t bytes[count]={};
    std::size_t total()const{
        std::size_t sum=0;for(auto b:bytes){if(b>std::numeric_limits<std::size_t>::max()-sum)throw std::overflow_error("Flat material budget overflow");sum+=b;}
        return sum;
    }
    std::size_t limiting()const{std::size_t at=0;for(std::size_t i=1;i<count;++i)if(bytes[i]>bytes[at])at=i;return at;}
};
inline std::size_t budgetProduct(std::size_t a,std::size_t b){
    if(b&&a>std::numeric_limits<std::size_t>::max()/b)throw std::overflow_error("Flat material budget overflow");return a*b;
}
inline std::size_t budgetSum(std::initializer_list<std::size_t> terms){
    std::size_t s=0;for(auto t:terms){if(t>std::numeric_limits<std::size_t>::max()-s)throw std::overflow_error("Flat material budget overflow");s+=t;}return s;
}
inline FlatEvaluationBudget flatEvaluationBudget(const FlatEvaluationGroups& groups,std::size_t lawDoubles){
    const auto mul=budgetProduct;const auto add=budgetSum;
    const auto nq=groups.at(0).points,nd=groups[0].dofs;
    constexpr std::size_t index=sizeof(std::size_t),real=sizeof(double),small=sizeof(int),key=sizeof(unsigned long long);
    FlatEvaluationBudget b;std::size_t triples=0;
    for(const auto& in:groups){const auto ne=in.elements,local=in.local,rule=in.rulePoints;
        const auto slots=mul(ne,local);triples=add({triples,mul(ne,mul(local,local))});
        // Resident element tables (geometry, shapes, rule, map table of a non-affine element, class, DOFs, points, slot order) plus the flag.
        b.bytes[0]=add({b.bytes[0],mul(mul(3*in.nodes,ne),real),in.nodes==4?std::size_t(0):mul(mul(10,mul(ne,rule)),real),mul(mul(mul(in.classes,rule),mul(local,3)),real),mul(mul(4,rule),real),
                        mul(add({ne,slots,mul(ne,rule),1}),small),slots});
        b.bytes[1]=add({b.bytes[1],mul(add({nd,1}),index),mul(slots,sizeof(unsigned))});
        b.bytes[2]=add({b.bytes[2],mul(slots,3*sizeof(unsigned))});
    }
    // The transpose of a mixed mesh sums one vector per group.
    if(groups.size()>1)b.bytes[1]=add({b.bytes[1],mul(nd,real)});
    b.bytes[3]=add({mul(triples,3*key+2*real+2*small),mul(nd,real),real,2*small});
    // The reference mass CSR and the frozen tangent values on its pattern.
    b.bytes[4]=add({mul(add({nd,1}),index),mul(triples,index+2*real),mul(nd,real)});
    b.bytes[5]=add({mul(add({mul(18,nq),mul(4,nd),2,lawDoubles}),real),small});
    return b;
}
struct LedgerRow {std::string payload,direction;std::size_t bytes=0;};
// Exact device allocation sizes: kept for the context lifetime vs held only during setup.
struct MaterialFootprint {std::size_t persistent=0,transient=0;};
// Diagnostic action y=Mx of the device-built reference mass with its sparsity count and Jacobi bound.
struct MassAction {radia_hmatrix_cuda::HostVector<double> y;std::size_t nonzeros=0;double bound=0;};
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
    virtual MaterialFootprint footprint()=0;
    virtual MassAction massAction(const radia_hmatrix_cuda::HostVector<double>&)=0;
    // Diagnostic E x or E^T y; vectors cross, no evaluation matrix does.
    virtual radia_hmatrix_cuda::HostVector<double> evaluationAction(const radia_hmatrix_cuda::HostVector<double>&,bool transpose)=0;
    virtual void close()=0;
};
// E and E^T act matrix-free on the device; no evaluation matrix crosses.
std::shared_ptr<MaterialOperator> createFlatMaterial(const FlatEvaluationGroups&,const Law&);
// Device-built reference mass alone, for the linear material solve: no law and
// no material point state. Vectors cross for diagnostics; the mass never does.
class MassOperator {
public:
    virtual ~MassOperator()=default;
    virtual MassAction massAction(const radia_hmatrix_cuda::HostVector<double>&)=0;
    virtual radia_hmatrix_cuda::HostVector<LedgerRow> transfers()=0;
    virtual MaterialFootprint footprint()=0;
    virtual std::size_t dofs()=0;
    virtual bool closed()=0;
    virtual void close()=0;
    // Geometry family of the device build: device-flat-tet-mass, device-curved-tet-mass or device-hex-mass.
    virtual std::string route()=0;
};
std::shared_ptr<MassOperator> createFlatMass(const FlatEvaluationGroups&);
// Host-only scope check of the flat tables (the same checks the context runs) and
// their device budget without CUB scratch or the law; no CUDA call.
FlatEvaluationBudget validateFlatEvaluationScope(const FlatEvaluationGroups&);
Samples evaluate(const Law&,const radia_hmatrix_cuda::HostVector<double>& amplitudes);
}
