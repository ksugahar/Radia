// Same packing, scaling, sparse maps and kernels as the standalone MEX gateway.
#include "hdiv_resident_python_api.hpp"
#include "hmatrix_material_solver.cuh"
#include "hmatrix_nonlinear_material.cuh"
#include "hdiv_flat_gram.cuh"
#include <mutex>
#include <atomic>
namespace radia_hdiv_python {
NewtonResult nonlinearBridge(radia_hmatrix_cuda::Resident&,radia_hdiv_bh::MaterialOperator&,const radia_hmatrix_cuda::HostVector<double>&,const radia_hmatrix_cuda::HostVector<double>&,const NewtonOptions&);
NewtonResult nonlinearHostReferenceBridge(radia_hmatrix_cuda::Resident&,radia_hdiv_bh::MaterialOperator&,const MapInput&,const radia_hmatrix_cuda::HostVector<double>&,const radia_hmatrix_cuda::HostVector<double>&,const NewtonOptions&);

namespace {
class Operator final:public DeviceOperator {
    std::mutex mutex;
    radia_hmatrix_cuda::HostUnique<radia_hmatrix_cuda::Resident> resident;
    const std::size_t dimension;
public:
    Assembly record;
    Operator(const radia_hmatrix_cuda::Plan& plan,const radia_hmatrix_cuda::SparseMap& map)
      :dimension(map.columns) {
        // Keep the existing 2 GiB combined native setup bound, not a new cap.
        resident=radia_hmatrix_cuda::makeHostUnique<radia_hmatrix_cuda::Resident>("python.resident-object",plan,&map,1);
        record={"host-leaves",plan.n,resident->ownedBytes,0,{}};
        record.transfers.push_back({"host_hmatrix_coefficients","H2D",plan.coefficients.size()*sizeof(double)});
        residentTransfers(plan,map);
    }
    Operator(const radia_hmatrix_cuda::Plan& plan,const radia_hmatrix_cuda::SparseMap& map,const radia_hmatrix_cuda::FlatBlock& block)
      :dimension(map.columns) {
        radia_hmatrix_cuda::FlatGramOracle oracle(block);
        resident=radia_hmatrix_cuda::makeHostUnique<radia_hmatrix_cuda::Resident>("python.resident-object",plan,&map,1);
        resident->materializeDeferred(plan,oracle);
        resident->setupPeakOwnedBytes=resident->ownedBytes+oracle.ownedBytes;
        record={"device-flat-tet",plan.n,resident->setupPeakOwnedBytes,oracle.launches,{}};
        residentTransfers(plan,map);
        for(const auto& t:oracle.transfers)record.transfers.push_back({t.payload,t.direction,t.bytes});
    }
    // Resident construction uploads: O(n) plan indices and the host-built charge map
    // (the remaining sparse-matrix upload of this slice, reported, not hidden).
    void residentTransfers(const radia_hmatrix_cuda::Plan& plan,const radia_hmatrix_cuda::SparseMap& map){
        const auto index=sizeof(std::size_t);
        record.transfers.push_back({"resident_plan","H2D",plan.directions.size()*sizeof(radia_hmatrix_cuda::Direction)+
            (plan.permutation.size()+plan.rowOffsets.size()+plan.partialIndices.size())*index+plan.scale.size()*sizeof(double)});
        record.transfers.push_back({"charge_map_csr_host_built","H2D",(map.offsets.size()+map.columns+1+2*map.indices.size())*index+2*map.values.size()*sizeof(double)});
    }
    radia_hmatrix_cuda::HostVector<double> apply(const radia_hmatrix_cuda::HostVector<double>& x)override {
        std::lock_guard<std::mutex> lock(mutex);
        if(!resident)throw std::runtime_error("HDiv resident has been closed");
        resident->upload(x);resident->apply(1);const auto y=resident->download();
        return {y.begin(),y.end()};
    }
    MaterialResult solve_linear(const MapInput& input,const radia_hmatrix_cuda::HostVector<double>& rhs,
        double invChi,double tolerance,int maxIterations,bool chebyshev,bool graph,int degree,double lower)override {
        std::lock_guard<std::mutex> lock(mutex);
        if(!resident)throw std::runtime_error("HDiv resident has been closed");
        if(input.rows!=dimension||input.columns!=dimension||rhs.size()!=dimension)
            throw std::invalid_argument("Material mass and RHS dimension mismatch");
        radia_hmatrix_cuda::SparseMap mass;mass.rows=input.rows;mass.columns=input.columns;
        mass.offsets.assign(input.offsets.begin(),input.offsets.end());
        mass.indices.assign(input.indices.begin(),input.indices.end());mass.values.assign(input.values.begin(),input.values.end());
        mass.validate();
        const auto result=radia_hmatrix_cuda::solveLinearMaterial(*resident,mass,rhs,invChi,tolerance,
            maxIterations,false,false,graph,chebyshev,degree,lower);
        return {{result.x.begin(),result.x.end()},result.iterations,result.relativeResidual,result.ownedBytes,result.preconditionerGraphLaunches};
    }
    NewtonResult solve_nonlinear(radia_hdiv_bh::MaterialOperator& material,const radia_hmatrix_cuda::HostVector<double>& rhs,const radia_hmatrix_cuda::HostVector<double>& initial,const NewtonOptions& options)override{
        std::lock_guard<std::mutex> lock(mutex);if(!resident)throw std::invalid_argument("HDiv resident has been closed");return nonlinearBridge(*resident,material,rhs,initial,options);
    }
    NewtonResult solve_nonlinear_host_reference(radia_hdiv_bh::MaterialOperator& material,const MapInput& reference,const radia_hmatrix_cuda::HostVector<double>& rhs,const radia_hmatrix_cuda::HostVector<double>& initial,const NewtonOptions& options)override{
        std::lock_guard<std::mutex> lock(mutex);if(!resident)throw std::invalid_argument("HDiv resident has been closed");return nonlinearHostReferenceBridge(*resident,material,reference,rhs,initial,options);
    }
    Assembly assembly()override {std::lock_guard<std::mutex> lock(mutex);return record;}
    Stats stats()override {
        std::lock_guard<std::mutex> lock(mutex);
        return {dimension,resident?resident->ownedBytes:0,resident?resident->applications:0,!resident};
    }
    void close()override {std::lock_guard<std::mutex> lock(mutex);resident.reset();}
};
}
std::shared_ptr<DeviceOperator> create(std::size_t n,const radia_hmatrix_cuda::HostVector<LeafInput>& leaves,
    const radia_hmatrix_cuda::HostVector<std::size_t>& permutation,const radia_hmatrix_cuda::HostVector<double>& scale,const MapInput& input) {
    if(!n||n>4608||!input.columns||input.columns>16384||input.rows!=n||leaves.empty()||leaves.size()>65536)
        throw std::invalid_argument("Unsupported bounded HDiv resident dimensions");
    std::size_t partials=0;
    for(const auto& leaf:leaves) {
        if(!leaf.rows||!leaf.columns||leaf.rows>n||leaf.columns>n)
            throw std::invalid_argument("Unsupported resident leaf dimensions");
        const auto count=leaf.rows+(leaf.row==leaf.column?0:leaf.columns);
        if(count>(std::size_t(1)<<24)-partials)
            throw std::invalid_argument("Resident host partial workspace budget exceeded");
        partials+=count;
    }
    auto plan=radia_hmatrix_cuda::packLeaves(n,leaves,permutation,scale);
    radia_hmatrix_cuda::SparseMap map;map.rows=input.rows;map.columns=input.columns;
    map.offsets.assign(input.offsets.begin(),input.offsets.end());
    map.indices.assign(input.indices.begin(),input.indices.end());map.values.assign(input.values.begin(),input.values.end());
    map.validate();return std::allocate_shared<Operator>(radia_hmatrix_cuda::HostFaultAllocator<Operator>("python.operator-owner"),plan,map);
}
std::shared_ptr<DeviceOperator> createFlat(const FlatInput& input,const MapInput& charges) {
    const auto n=input.host.size();
    if(!n||n>4608||charges.rows!=n||!charges.columns||charges.columns>16384)
        throw std::invalid_argument("Unsupported bounded flat HDiv resident dimensions");
    radia_hmatrix_cuda::FlatBlock block;
    block.cells.assign(input.cells.begin(),input.cells.end());block.faces.assign(input.faces.begin(),input.faces.end());
    block.host.assign(input.host.begin(),input.host.end());block.kind.assign(input.kind.begin(),input.kind.end());
    block.expo.assign(input.expo.begin(),input.expo.end());block.farFactor=input.farFactor;
    auto rule=[](const radia_hmatrix_cuda::HostVector<double>& rows,auto& out){
        if(rows.size()%4)throw std::invalid_argument("Reference rule rows require (xi0,xi1,xi2,weight)");
        for(std::size_t i=0;i<rows.size();i+=4)out.push_back({{rows[i],rows[i+1],rows[i+2]},rows[i+3]});};
    rule(input.tet,block.tet);rule(input.tri,block.tri);rule(input.tetLow,block.tetLow);rule(input.triLow,block.triLow);
    block.validate();
    // One deferred dense leaf: device storage is allocated, never uploaded.
    radia_hmatrix_cuda::HostVector<LeafInput> leaves;leaves.push_back(LeafInput{0,0,n,n,0,{},{},true});
    radia_hmatrix_cuda::HostVector<std::size_t> permutation(n);for(std::size_t i=0;i<n;++i)permutation[i]=i;
    auto plan=radia_hmatrix_cuda::packLeaves(n,leaves,permutation,radia_hmatrix_cuda::HostVector<double>(n,1.));
    radia_hmatrix_cuda::SparseMap map;map.rows=charges.rows;map.columns=charges.columns;
    map.offsets.assign(charges.offsets.begin(),charges.offsets.end());
    map.indices.assign(charges.indices.begin(),charges.indices.end());map.values.assign(charges.values.begin(),charges.values.end());
    map.validate();return std::allocate_shared<Operator>(radia_hmatrix_cuda::HostFaultAllocator<Operator>("python.operator-owner"),plan,map,block);
}
} // namespace radia_hdiv_python

radia_hdiv_bh::Samples radia_hdiv_bh::evaluate(const Law& law,const radia_hmatrix_cuda::HostVector<double>& a){return evaluateInverseBH(law,a);}

radia_hdiv_bh::MaterialEvaluation radia_hdiv_bh::evaluateMaterial(const MaterialInput& input,const Law& law,const radia_hmatrix_cuda::HostVector<double>& state,const radia_hmatrix_cuda::HostVector<double>& direction){return evaluateMappedMaterial(input,law,state,direction);}

namespace radia_hdiv_bh {
namespace {
std::atomic<unsigned> liveMaterialOwners{0};
class MaterialOwner final:public MaterialOperator {
    std::mutex mutex;radia_hmatrix_cuda::HostUnique<MappedMaterial> context;bool counted=false;
public:
    template<class Input>MaterialOwner(const Input& input,const Law& law){
        if(liveMaterialOwners.fetch_add(1)>=2){--liveMaterialOwners;throw std::invalid_argument("Material context count budget exceeded");}
        try{context=radia_hmatrix_cuda::makeHostUnique<MappedMaterial>("python.material-object",input,law);counted=true;}catch(...){--liveMaterialOwners;throw;}
    }
    ~MaterialOwner()override{close();}
    MaterialEvaluation evaluate(const radia_hmatrix_cuda::HostVector<double>& x,const radia_hmatrix_cuda::HostVector<double>& d)override{
        std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");return context->evaluate(x,d);
    }
    MaterialEvaluation linearize(const radia_hmatrix_cuda::HostVector<double>& x)override{std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");return context->linearize(x);}
    radia_hmatrix_cuda::HostVector<double> tangent(const radia_hmatrix_cuda::HostVector<double>& d,std::size_t epoch)override{std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");return context->tangentAction(d,epoch);}
    // The reference mass is the context's device-built M; no host mass crosses.
    radia_hdiv_python::NewtonResult solve(radia_hmatrix_cuda::Resident& op,const radia_hmatrix_cuda::HostVector<double>& rhs,const radia_hmatrix_cuda::HostVector<double>& initial,const radia_hdiv_python::NewtonOptions& options){
        std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");
        return radia_hmatrix_cuda::solveNonlinearMaterial(op,*context,rhs,initial,options);
    }
    // Host-evaluation contexts have no device mass; a device-table context never accepts a host reference.
    radia_hdiv_python::NewtonResult solveHostReference(radia_hmatrix_cuda::Resident& op,const radia_hdiv_python::MapInput& input,const radia_hmatrix_cuda::HostVector<double>& rhs,const radia_hmatrix_cuda::HostVector<double>& initial,const radia_hdiv_python::NewtonOptions& options){
        std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");
        if(context->deviceMass())throw std::invalid_argument("Device-table material owns its reference mass; host reference rejected");
        if(input.rows!=op.dimension()||input.columns!=op.dimension())throw std::invalid_argument("Reference mass dimension mismatch");
        radia_hmatrix_cuda::SparseMap reference;reference.rows=input.rows;reference.columns=input.columns;reference.offsets.assign(input.offsets.begin(),input.offsets.end());reference.indices.assign(input.indices.begin(),input.indices.end());reference.values.assign(input.values.begin(),input.values.end());reference.validate();
        return radia_hmatrix_cuda::solveNonlinearMaterial(op,*context,reference,rhs,initial,options);
    }
    MaterialStats stats()override{std::lock_guard<std::mutex> lock(mutex);if(context)return context->stats();MaterialStats out;out.closed=true;return out;}
    std::string route()override{std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");return context->route();}
    radia_hmatrix_cuda::HostVector<LedgerRow> transfers()override{std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");return context->transfers();}
    MaterialFootprint footprint()override{std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");return context->footprint();}
    void close()override{std::lock_guard<std::mutex> lock(mutex);context.reset();if(counted){--liveMaterialOwners;counted=false;}}
};
}
std::shared_ptr<MaterialOperator> createMaterial(const MaterialInput& input,const Law& law){return std::allocate_shared<MaterialOwner>(radia_hmatrix_cuda::HostFaultAllocator<MaterialOwner>("python.material-owner"),input,law);}
std::shared_ptr<MaterialOperator> createFlatMaterial(const FlatEvaluationInput& input,const Law& law){return std::allocate_shared<MaterialOwner>(radia_hmatrix_cuda::HostFaultAllocator<MaterialOwner>("python.material-owner"),input,law);}
}

radia_hdiv_python::NewtonResult radia_hdiv_python::nonlinearBridge(radia_hmatrix_cuda::Resident& op,radia_hdiv_bh::MaterialOperator& material,const radia_hmatrix_cuda::HostVector<double>& rhs,const radia_hmatrix_cuda::HostVector<double>& initial,const NewtonOptions& options){
    auto* owner=dynamic_cast<radia_hdiv_bh::MaterialOwner*>(&material);if(!owner)throw std::invalid_argument("Material owner ABI mismatch");return owner->solve(op,rhs,initial,options);
}
radia_hdiv_python::NewtonResult radia_hdiv_python::nonlinearHostReferenceBridge(radia_hmatrix_cuda::Resident& op,radia_hdiv_bh::MaterialOperator& material,const MapInput& reference,const radia_hmatrix_cuda::HostVector<double>& rhs,const radia_hmatrix_cuda::HostVector<double>& initial,const NewtonOptions& options){
    auto* owner=dynamic_cast<radia_hdiv_bh::MaterialOwner*>(&material);if(!owner)throw std::invalid_argument("Material owner ABI mismatch");return owner->solveHostReference(op,reference,rhs,initial,options);
}
