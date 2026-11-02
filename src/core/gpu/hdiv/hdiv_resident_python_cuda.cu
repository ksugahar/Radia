// Same packing, scaling, sparse maps and kernels as the standalone MEX gateway.
#include "hdiv_resident_python_api.hpp"
#include "hmatrix_material_solver.cuh"
#include "hmatrix_nonlinear_material.cuh"
#include "hdiv_flat_gram.cuh"
#include "hdiv_flat_aca.cuh"
#include "hdiv_flat_charge_map.cuh"
#include "hdiv_p1_diagonal.cuh"
#include <mutex>
#include <sstream>
#include <atomic>
namespace radia_hdiv_python {
NewtonResult nonlinearBridge(radia_hmatrix_cuda::Resident&,radia_hdiv_bh::MaterialOperator&,const radia_hmatrix_cuda::HostVector<double>&,const radia_hmatrix_cuda::HostVector<double>&,const NewtonOptions&);
NewtonResult nonlinearHostReferenceBridge(radia_hmatrix_cuda::Resident&,radia_hdiv_bh::MaterialOperator&,const MapInput&,const radia_hmatrix_cuda::HostVector<double>&,const radia_hmatrix_cuda::HostVector<double>&,const NewtonOptions&);

namespace {
// Named device bytes known before any allocation, summed as a conservative
// setup bound. ACA factor storage and the rank-sized plan buffers are known
// only after compression; their allocation sites keep live reserve checks.
struct ChargeBudget {
    static constexpr std::size_t count=6;
    static constexpr const char* names[count]={"charge_prepare","charge_map","charge_map_scratch","resident_vectors","gram_dense","aca_batch_scratch"};
    std::size_t bytes[count]{};
    std::size_t total()const{std::size_t t=0;for(auto b:bytes)t=radia_hmatrix_cuda::flatSum(t,b);return t;}
    std::size_t limiting()const{std::size_t k=0;for(std::size_t i=1;i<count;++i)if(bytes[i]>bytes[k])k=i;return k;}
};
ChargeBudget chargeBudget(const radia_hmatrix_cuda::FlatBlock& block,const FlatInput& input,const ChargeMapInput& charges,bool aca){
    using radia_hmatrix_cuda::checkedProduct;using radia_hmatrix_cuda::flatSum;
    const auto n=block.charges(),dofs=charges.dofs,nnz=radia_hdiv_bh::FlatChargeMapBuilder::nonzeros(charges);
    ChargeBudget b;
    b.bytes[0]=flatSum(radia_hmatrix_cuda::FlatPrepared::workspaceBytes(block),aca?checkedProduct(n,2*sizeof(double))+sizeof(int):0);
    b.bytes[1]=flatSum(checkedProduct(flatSum(checkedProduct(2,nnz),n+dofs+2),sizeof(std::size_t)),checkedProduct(flatSum(checkedProduct(2,nnz),dofs),sizeof(double)));
    b.bytes[2]=radia_hdiv_bh::FlatChargeMapBuilder::scratchBound(input,charges);
    b.bytes[3]=checkedProduct(flatSum(checkedProduct(5,n),1),sizeof(double));
    b.bytes[4]=aca?0:checkedProduct(checkedProduct(n,n),sizeof(double));
    b.bytes[5]=aca?radia_hmatrix_cuda::FlatAcaOracle::batchScratchBytes:0;
    std::size_t free,total;radia_hmatrix_cuda::check(cudaMemGetInfo(&free,&total));
    const auto reserve=radia_hmatrix_cuda::reserveBytes(total),available=free>reserve?free-reserve:0,need=b.total();
    if(need>available){const auto k=b.limiting();std::ostringstream m;
        m<<"Flat charge device budget exceeded: "<<need<<" bytes needed, "<<available<<" free above the reserve; limiting term "
         <<ChargeBudget::names[k]<<" = "<<b.bytes[k]<<" bytes";throw std::runtime_error(m.str());}
    return b;
}
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
        planTransfer(plan);
        // The host-leaf route still uploads its host-built charge map (reported, not hidden).
        const auto index=sizeof(std::size_t);
        record.transfers.push_back({"charge_map_csr_host_built","H2D",(map.offsets.size()+map.columns+1+2*map.indices.size())*index+2*map.values.size()*sizeof(double)});
    }
    // B and B^T are generated on the device from NGSolve tables; no CSR crosses.
    // The dense oracle allocates its workspace after B is built; ACA factors
    // and prepared charges are already live while the resident and B exist.
    template<class Oracle>
    // The checked charge budget and the live reserve checks replace the fixed
    // 2 GiB combined setup bound on this route.
    Operator(const radia_hmatrix_cuda::Plan& plan,Oracle& oracle,bool oracleLive,const FlatInput& flat,const ChargeMapInput& charges,const ChargeBudget& budget)
      :dimension(charges.dofs) {
        radia_hdiv_bh::FlatChargeMapBuilder builder;
        resident=radia_hmatrix_cuda::makeHostUnique<radia_hmatrix_cuda::Resident>("python.resident-object",plan,nullptr,0);
        resident->buildChargeMap(charges.dofs,builder.nonzeros(charges),radia_hdiv_bh::FlatChargeMapBuilder::scratchBound(flat,charges),
            [&](radia_hmatrix_cuda::DeviceMap& B,radia_hmatrix_cuda::DeviceMap& BT){builder.build(flat,charges,B,BT);});
        resident->materializeDeferred(plan,oracle);
        resident->setupPeakOwnedBytes=resident->ownedBytes+(oracleLive?oracle.ownedBytes+builder.scratchBytes:std::max(oracle.ownedBytes,builder.scratchBytes));
        record={"device-flat-tet",plan.n,resident->setupPeakOwnedBytes,oracle.launches+builder.launches,{}};
        record.compression=flat.compression;
        for(std::size_t i=0;i<ChargeBudget::count;++i)record.budget.push_back({ChargeBudget::names[i],budget.bytes[i]});
        planTransfer(plan);
        for(const auto& t:builder.transfers)record.transfers.push_back({t.payload,t.direction,t.bytes});
        for(const auto& t:oracle.transfers)record.transfers.push_back({t.payload,t.direction,t.bytes});
    }
    // Resident construction uploads O(n) plan indices.
    void planTransfer(const radia_hmatrix_cuda::Plan& plan){
        const auto index=sizeof(std::size_t);
        record.transfers.push_back({"resident_plan","H2D",plan.directions.size()*sizeof(radia_hmatrix_cuda::Direction)+
            (plan.permutation.size()+plan.rowOffsets.size()+plan.partialIndices.size())*index+plan.scale.size()*sizeof(double)});
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
std::shared_ptr<DeviceOperator> createFlat(const FlatInput& input,const ChargeMapInput& charges) {
    const auto n=input.host.size();
    if(!n||!charges.dofs)throw std::invalid_argument("Unsupported bounded flat HDiv resident dimensions");
    radia_hmatrix_cuda::FlatBlock block;
    block.cells.assign(input.cells.begin(),input.cells.end());block.faces.assign(input.faces.begin(),input.faces.end());
    block.host.assign(input.host.begin(),input.host.end());block.kind.assign(input.kind.begin(),input.kind.end());
    block.expo.assign(input.expo.begin(),input.expo.end());block.farFactor=input.farFactor;
    auto rule=[](const radia_hmatrix_cuda::HostVector<double>& rows,auto& out){
        if(rows.size()%4)throw std::invalid_argument("Reference rule rows require (xi0,xi1,xi2,weight)");
        for(std::size_t i=0;i<rows.size();i+=4)out.push_back({{rows[i],rows[i+1],rows[i+2]},rows[i+3]});};
    rule(input.tet,block.tet);rule(input.tri,block.tri);rule(input.tetLow,block.tetLow);rule(input.triLow,block.triLow);
    block.validate();radia_hdiv_bh::validateFlatChargeMap(input,charges);
    const auto budget=chargeBudget(block,input,charges,input.compression=="aca");
    const auto owner=radia_hmatrix_cuda::HostFaultAllocator<Operator>("python.operator-owner");
    if(input.compression=="dense"){
        if(input.acaEps!=0||input.acaEta!=0||input.acaLeafSize||input.acaMaxRank||input.acaRankCapDense)
            throw std::invalid_argument("Dense flat Gram takes no ACA parameters");
        // One deferred dense leaf: device storage is allocated, never uploaded.
        radia_hmatrix_cuda::HostVector<LeafInput> leaves;leaves.push_back(LeafInput{0,0,n,n,0,{},{},true});
        radia_hmatrix_cuda::HostVector<std::size_t> permutation(n);for(std::size_t i=0;i<n;++i)permutation[i]=i;
        auto plan=radia_hmatrix_cuda::packLeaves(n,leaves,permutation,radia_hmatrix_cuda::HostVector<double>(n,1.));
        radia_hmatrix_cuda::FlatGramOracle oracle(block);
        return std::allocate_shared<Operator>(owner,plan,oracle,false,input,charges,budget);
    }
    if(input.compression=="aca"){
        radia_hmatrix_cuda::FlatAcaOptions options;options.eps=input.acaEps;options.eta=input.acaEta;
        options.leafSize=input.acaLeafSize;options.maxRank=input.acaMaxRank;options.rankCapDense=input.acaRankCapDense;
        // Ranks are known only after device compression; they size the plan.
        radia_hmatrix_cuda::FlatAcaOracle oracle(block,options);oracle.compress();
        auto plan=radia_hmatrix_cuda::packLeaves(n,oracle.leaves<LeafInput>(),oracle.permutation(),radia_hmatrix_cuda::HostVector<double>(n,1.));
        auto result=std::allocate_shared<Operator>(owner,plan,oracle,true,input,charges,budget);
        const auto& s=oracle.stats;
        result->record.counters={{"elements",s.elements},{"clusters",s.clusters},{"dense_leaves",s.denseLeaves},
            {"lowrank_leaves",s.lowrankLeaves},{"admissible_blocks",s.admissibleBlocks},{"rank_cap_promotions",s.rankCapPromotions},
            {"storage_promotions",s.storagePromotions},{"zero_rank_leaves",s.zeroRankLeaves},{"rank_sum",s.rankSum},{"max_rank",s.maxRank},{"aca_steps",s.steps},
            {"aca_batches",s.batches},{"factor_coefficients",s.factorCoefficients},{"dense_coefficients",s.denseCoefficients}};
        return result;
    }
    throw std::invalid_argument("Flat charge Gram compression must be explicit 'dense' or 'aca'");
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
    MassAction massAction(const radia_hmatrix_cuda::HostVector<double>& x)override{std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");return context->massAction(x);}
    radia_hmatrix_cuda::HostVector<double> evaluationAction(const radia_hmatrix_cuda::HostVector<double>& x,bool transpose)override{std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");return context->evaluationAction(x,transpose);}
    void close()override{std::lock_guard<std::mutex> lock(mutex);context.reset();if(counted){--liveMaterialOwners;counted=false;}}
};
}
std::shared_ptr<MaterialOperator> createMaterial(const MaterialInput& input,const Law& law){return std::allocate_shared<MaterialOwner>(radia_hmatrix_cuda::HostFaultAllocator<MaterialOwner>("python.material-owner"),input,law);}
std::shared_ptr<MaterialOperator> createFlatMaterial(const FlatEvaluationInput& input,const Law& law){return std::allocate_shared<MaterialOwner>(radia_hmatrix_cuda::HostFaultAllocator<MaterialOwner>("python.material-owner"),input,law);}
FlatEvaluationBudget validateFlatEvaluationScope(const FlatEvaluationInput& input){validateFlatEvaluation(input);return flatEvaluationBudget(input,0);}
}

radia_hdiv_python::NewtonResult radia_hdiv_python::nonlinearBridge(radia_hmatrix_cuda::Resident& op,radia_hdiv_bh::MaterialOperator& material,const radia_hmatrix_cuda::HostVector<double>& rhs,const radia_hmatrix_cuda::HostVector<double>& initial,const NewtonOptions& options){
    auto* owner=dynamic_cast<radia_hdiv_bh::MaterialOwner*>(&material);if(!owner)throw std::invalid_argument("Material owner ABI mismatch");return owner->solve(op,rhs,initial,options);
}
radia_hdiv_python::NewtonResult radia_hdiv_python::nonlinearHostReferenceBridge(radia_hmatrix_cuda::Resident& op,radia_hdiv_bh::MaterialOperator& material,const MapInput& reference,const radia_hmatrix_cuda::HostVector<double>& rhs,const radia_hmatrix_cuda::HostVector<double>& initial,const NewtonOptions& options){
    auto* owner=dynamic_cast<radia_hdiv_bh::MaterialOwner*>(&material);if(!owner)throw std::invalid_argument("Material owner ABI mismatch");return owner->solveHostReference(op,reference,rhs,initial,options);
}
