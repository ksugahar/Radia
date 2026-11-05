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
struct NewtonOptions {double tolerance=1e-11,linearTolerance=5e-14;int maximumIterations=40,maximumLinearIterations=10000,preconditioner=0,degree=32;double lowerFraction=.01;};
struct NewtonStep {double residual=0,energy=0,alpha=0,linearResidual=0;int linearIterations=0,backtracks=0;bool residualCertified=false;};
struct NewtonResult {radia_hdiv_bh::BoundaryVector<double> x;radia_hdiv_bh::BoundaryVector<NewtonStep> history;double residual=0,energy=0;std::size_t graphLaunches=0;};
// Linear material solve A=B^T G B+inv_chi M over a device-built mass context.
// preconditioner: 0 mass-diagonal, 1 mass-cholesky, 2 mass-polynomial, 3 mass-polynomial-graph,
// 4 mass-chebyshev, 5 mass-chebyshev-graph.
struct LinearOptions {double invChi=0,tolerance=0;int maximumIterations=0,preconditioner=0,degree=8;double lowerFraction=.1;};
struct LinearResult {radia_hdiv_bh::BoundaryVector<double> x;int iterations=0;double relativeResidual=0;std::size_t ownedBytes=0,graphLaunches=0;};
struct Stats {std::size_t dimension,owned_bytes,applications;bool closed;};
// Primitive TET charge geometry; reference rules are (xi0,xi1,xi2,weight) rows.
// family "affine-tet": cells/faces hold 4/3 vertices per host. "curved-tet":
// cells/faces hold the 10/6 P2 nodes of the isoparametric map (corners
// first), the mesh vertex numbers decide touching hosts and (gl,gw) is the
// [0,1] Gauss rule of the touching Duffy potentials. "hex": cells/faces hold
// the 27/9-node Q2 lattices, charges are reference monomials and hexRules /
// hexBands carry the CPU hex Gram rules; TET rules stay empty. "mixed": as
// "hex" with the unit-cube lattices of tet, pyramid and hex hosts and hostTypes.
// compression is "dense" (one deferred dense block), "aca" (device ACA leaves on a geometric partition) or
// "hca" (mixed meshes: far blocks from the Laplace kernel expansion); the ACA parameters apply to the last two.
struct FlatInput {
    radia_hmatrix_cuda::HostVector<double> cells,faces,tet,tri,tetLow,triLow;
    radia_hmatrix_cuda::HostVector<int> host,kind,expo;
    double farFactor=1e30;
    std::string compression;
    double acaEps=0,acaEta=0;std::size_t acaLeafSize=0,acaMaxRank=0;bool acaRankCapDense=false;int acaExpansionOrder=0;
    std::string family="affine-tet";
    radia_hmatrix_cuda::HostVector<int> cellNumbers,faceNumbers;
    radia_hmatrix_cuda::HostVector<double> gl,gw;
    // HEX: rows (rule, point, weight) with rule 0 outer, 1 radial inner,
    // 2 pair product, 3 affine pair product, 4 near outer, 5 fine inner (BDM1 only), 6 plain near outer of cell targets
    // (mixed only, absent when they take rule 4), 7..9 far levels (mixed only) and 10 the tangent remainder of a warped
    // source (mixed only, absent when it keeps the radial and fan rules); bands = (near_grade, far_inner,
    // inner_switch, then the separation ratio of each far level, 0 when absent).
    radia_hmatrix_cuda::HostVector<double> hexRules,hexBands;
    radia_hmatrix_cuda::HostVector<int> hostTypes;
};
// NGSolve charge-projection tables per vertex-rank class: B rows are generated
// on the device from these and the FlatInput geometry; no B value crosses.
// Exponents check the layout against FlatInput on the host; only the curved map uploads them.
// Affine: tables are the reference projections, scaled per host on the device.
// Curved: tables are the reference moments R and the device solves the
// per-host monomial Gram G_e X = R on NGSolve's declared mass rules
// (cellRule/faceRule rows), so B = G_e^{-1} R carries the curved measure.
// HEX: tables are the reference-measure projections and B rows are the class
// tables unscaled (the Piola Jacobians cancel in the reference-frame charges).
// Mixed: as HEX with one table per class in classTables; classInfo holds
// (table start, monomials, width) per class, hostClass the class per host and
// hostDofs the DOFs of all hosts in host order. Other table fields stay empty.
struct ChargeMapInput {
    std::size_t dofs=0,cellMonomials=0,faceMonomials=0,cellWidth=0,faceWidth=0;
    radia_hmatrix_cuda::HostVector<double> cellTables,faceTables,cellRule,faceRule;
    radia_hmatrix_cuda::HostVector<int> cellClass,cellDofs,faceClass,faceDofs,cellExponents,faceExponents;
    bool curved=false,hex=false,mixed=false;
    radia_hmatrix_cuda::HostVector<double> classTables;
    radia_hmatrix_cuda::HostVector<int> hostClass,classInfo,hostDofs;
};
struct Transfer {std::string payload,direction;std::size_t bytes;};
struct Assembly {std::string route;std::size_t charges,setupBytes,launches;radia_hmatrix_cuda::HostVector<Transfer> transfers;
    std::string compression;radia_hmatrix_cuda::HostVector<std::pair<std::string,std::size_t>> counters;
    // Device-flat route only: the checked pre-allocation charge budget, in publication order.
    radia_hmatrix_cuda::HostVector<std::pair<std::string,std::size_t>> budget;};
class DeviceOperator {
public:
    virtual ~DeviceOperator()=default;
    virtual radia_hmatrix_cuda::HostVector<double> apply(const radia_hmatrix_cuda::HostVector<double>&)=0;
    virtual NewtonResult solve_nonlinear(radia_hdiv_bh::MaterialOperator&,const radia_hmatrix_cuda::HostVector<double>&,const radia_hmatrix_cuda::HostVector<double>&,const NewtonOptions&)=0;
    virtual LinearResult solve_linear(radia_hdiv_bh::MassOperator&,const radia_hmatrix_cuda::HostVector<double>&,const LinearOptions&)=0;
    virtual Stats stats()=0;
    virtual Assembly assembly()=0;
    virtual void close()=0;
};
// The charge Gram, B and B^T are generated on the device; no Gram or B value crosses.
std::shared_ptr<DeviceOperator> createFlat(const FlatInput&,const ChargeMapInput&);
// Gram-only constant-density block filled on the device from flat row-major linear
// tetrahedra (12 coordinates) and triangles (9); no charge map, no material solve.
std::shared_ptr<DeviceOperator> createAnalyticGram(const radia_hmatrix_cuda::HostVector<double>& cells,const radia_hmatrix_cuda::HostVector<double>& faces);
// Diagnostic: raw device flat-Gram entries for explicit charge pairs; no matrix crosses.
radia_hmatrix_cuda::HostVector<double> flatDiagnosticEntries(const FlatInput&,const radia_hmatrix_cuda::HostVector<std::size_t>& rows,
    const radia_hmatrix_cuda::HostVector<std::size_t>& cols);
} // namespace radia_hdiv_python
