// Device-generated charge Gram, B and B^T; resident kernels shared with the standalone MEX gateway.
#include "hdiv_resident_python_api.hpp"
#include "hmatrix_material_solver.cuh"
#include "hmatrix_nonlinear_material.cuh"
#include "hdiv_flat_gram.cuh"
#include "hdiv_flat_diagnostic.cuh"
#include "hdiv_flat_aca.cuh"
#include "hdiv_aca_diagnostic.cuh"
#include "hdiv_flat_charge_map.cuh"
#include "hdiv_flat_field.cuh"
#include <mutex>
#include <sstream>
#include <atomic>
namespace radia_hdiv_python {
NewtonResult nonlinearBridge(radia_hmatrix_cuda::Resident&,radia_hdiv_bh::MaterialOperator&,const radia_hmatrix_cuda::HostVector<double>&,const radia_hmatrix_cuda::HostVector<double>&,const NewtonOptions&);
LinearResult linearBridge(radia_hmatrix_cuda::Resident&,radia_hdiv_bh::MassOperator&,const radia_hmatrix_cuda::HostVector<double>&,const LinearOptions&);

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
radia_hmatrix_cuda::FlatBlock flatBlock(const FlatInput& input);
class Operator final:public DeviceOperator {
    std::mutex mutex;
    radia_hmatrix_cuda::HostUnique<radia_hmatrix_cuda::Resident> resident;
    const std::size_t dimension;
    const bool gramOnly=false;
    radia_hmatrix_cuda::HostUnique<radia_hmatrix_cuda::FlatBlock> fieldBlock;
public:
    Assembly record;
    Operator(const radia_hmatrix_cuda::Plan& plan,const radia_hmatrix_cuda::AnalyticBlock& block):dimension(plan.n),gramOnly(true) {
        resident=radia_hmatrix_cuda::makeHostUnique<radia_hmatrix_cuda::Resident>("python.resident-object",plan,nullptr,0);
        resident->assemble(block);
        record={"device-analytic-constant-block",plan.n,std::max(resident->ownedBytes,resident->setupPeakOwnedBytes),3,{}};
        record.compression="dense";
        planTransfer(plan);
        record.transfers.push_back({"analytic_geometry","H2D",block.sources.size()*sizeof(radia_hmatrix_cuda::AnalyticSource)+block.quadrature.size()*sizeof(radia_hmatrix_cuda::QuadraturePoint)});
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
        record={flat.family=="mixed"?"device-mixed":flat.family=="hex"?"device-hex":flat.family=="curved-tet"?"device-curved-tet":"device-flat-tet",plan.n,resident->setupPeakOwnedBytes,oracle.launches+builder.launches,{}};
        record.compression=flat.compression;
        if(flat.family!="mixed")fieldBlock=radia_hmatrix_cuda::makeHostUnique<radia_hmatrix_cuda::FlatBlock>("python.field-block",flatBlock(flat));
        for(std::size_t i=0;i<ChargeBudget::count;++i)record.budget.push_back({ChargeBudget::names[i],budget.bytes[i]});
        planTransfer(plan);
        for(const auto& t:builder.transfers)record.transfers.push_back({t.payload,t.direction,t.bytes});
        for(const auto& t:oracle.transfers)record.transfers.push_back({t.payload,t.direction,t.bytes});
    }
    // Resident construction uploads O(n) plan indices and any host-held coefficients; device-built plans hold none.
    void planTransfer(const radia_hmatrix_cuda::Plan& plan){
        const auto index=sizeof(std::size_t);
        record.transfers.push_back({"resident_plan","H2D",plan.directions.size()*sizeof(radia_hmatrix_cuda::Direction)+
            (plan.permutation.size()+plan.rowOffsets.size()+plan.partialIndices.size())*index+(plan.scale.size()+plan.coefficients.size())*sizeof(double)});
    }
    radia_hmatrix_cuda::HostVector<double> apply(const radia_hmatrix_cuda::HostVector<double>& x)override {
        std::lock_guard<std::mutex> lock(mutex);
        if(!resident)throw std::runtime_error("HDiv resident has been closed");
        resident->upload(x);resident->apply(1);const auto y=resident->download();
        return {y.begin(),y.end()};
    }
    NewtonResult solve_nonlinear(radia_hdiv_bh::MaterialOperator& material,const radia_hmatrix_cuda::HostVector<double>& rhs,const radia_hmatrix_cuda::HostVector<double>& initial,const NewtonOptions& options)override{
        std::lock_guard<std::mutex> lock(mutex);if(!resident)throw std::invalid_argument("HDiv resident has been closed");
        if(gramOnly)throw std::invalid_argument("The analytic Gram has no charge map or material solve");
        return nonlinearBridge(*resident,material,rhs,initial,options);
    }
    LinearResult solve_linear(radia_hdiv_bh::MassOperator& mass,const radia_hmatrix_cuda::HostVector<double>& rhs,const LinearOptions& options)override{
        std::lock_guard<std::mutex> lock(mutex);if(!resident)throw std::invalid_argument("HDiv resident has been closed");
        if(gramOnly)throw std::invalid_argument("The analytic Gram has no charge map or material solve");
        return linearBridge(*resident,mass,rhs,options);
    }
    FieldResult field(const radia_hmatrix_cuda::HostVector<double>& x,const radia_hmatrix_cuda::HostVector<double>& points)override{
        std::lock_guard<std::mutex> lock(mutex);if(!resident)throw std::invalid_argument("HDiv resident has been closed");
        if(!fieldBlock)throw std::invalid_argument("Device field evaluation has no route for this geometry family");
        auto computed=fieldBlock->hex?radia_hmatrix_cuda::hexDeviceField(*fieldBlock,*resident,x,points):radia_hmatrix_cuda::flatDeviceField(*fieldBlock,*resident,x,points);
        FieldResult out;out.field.assign(computed.field.begin(),computed.field.end());out.launches=computed.launches;out.ownedBytes=computed.ownedBytes;out.route=computed.route;
        for(const auto& t:computed.transfers)out.transfers.push_back({t.payload,t.direction,t.bytes});
        return out;
    }
    Assembly assembly()override {std::lock_guard<std::mutex> lock(mutex);return record;}
    Stats stats()override {
        std::lock_guard<std::mutex> lock(mutex);
        return {dimension,resident?resident->ownedBytes:0,resident?resident->applications:0,!resident};
    }
    void close()override {std::lock_guard<std::mutex> lock(mutex);resident.reset();}
};
}
namespace {
radia_hmatrix_cuda::FlatBlock flatBlock(const FlatInput& input){
    radia_hmatrix_cuda::FlatBlock block;
    block.cells.assign(input.cells.begin(),input.cells.end());block.faces.assign(input.faces.begin(),input.faces.end());
    block.host.assign(input.host.begin(),input.host.end());block.kind.assign(input.kind.begin(),input.kind.end());
    block.expo.assign(input.expo.begin(),input.expo.end());block.farFactor=input.farFactor;
    if(input.family!="affine-tet"&&input.family!="curved-tet"&&input.family!="hex"&&input.family!="mixed")
        throw std::invalid_argument("Charge geometry family must be 'affine-tet', 'curved-tet', 'hex' or 'mixed'");
    if(input.hostTypes.empty()==(input.family=="mixed"))throw std::invalid_argument("Host types are required for family 'mixed' and rejected otherwise");
    block.curved=input.family=="curved-tet";block.hex=input.family=="hex"||input.family=="mixed";
    block.hostTypes.assign(input.hostTypes.begin(),input.hostTypes.end());
    if(block.hex){
        if(input.hexRules.size()%3||input.hexBands.size()!=6)throw std::invalid_argument("HEX geometry requires (rule,point,weight) rows and six bands");
        radia_hmatrix_cuda::HostVector<double>* points[11]={&block.ruleOut,&block.ruleIn,&block.rulePair,&block.rulePairAffine,&block.ruleNear,&block.ruleInSelf,&block.ruleNearCell,&block.ruleFar[0],&block.ruleFar[1],&block.ruleFar[2],&block.ruleTangent};
        radia_hmatrix_cuda::HostVector<double>* weights[11]={&block.ruleOutW,&block.ruleInW,&block.rulePairW,&block.rulePairAffineW,&block.ruleNearW,&block.ruleInSelfW,&block.ruleNearCellW,&block.ruleFarW[0],&block.ruleFarW[1],&block.ruleFarW[2],&block.ruleTangentW};
        int previous=0;
        for(std::size_t i=0;i<input.hexRules.size();i+=3){const double r=input.hexRules[i];
            if(!(r>=0&&r<=10&&r==std::floor(r))||int(r)<previous)throw std::invalid_argument("HEX rule rows must be grouped by rule index 0..10");
            previous=int(r);points[previous]->push_back(input.hexRules[i+1]);weights[previous]->push_back(input.hexRules[i+2]);}
        block.nearGrade=input.hexBands[0];block.farInner=input.hexBands[1];block.innerSwitch=input.hexBands[2];
        for(int k=0;k<3;++k)block.farRatio[k]=input.hexBands[3+k];
    }else if(!input.hexRules.empty()||!input.hexBands.empty())throw std::invalid_argument("TET geometry takes no HEX rules");
    block.gl.assign(input.gl.begin(),input.gl.end());block.gw.assign(input.gw.begin(),input.gw.end());
    block.cellNumbers.assign(input.cellNumbers.begin(),input.cellNumbers.end());block.faceNumbers.assign(input.faceNumbers.begin(),input.faceNumbers.end());
    auto rule=[](const radia_hmatrix_cuda::HostVector<double>& rows,auto& out){
        if(rows.size()%4)throw std::invalid_argument("Reference rule rows require (xi0,xi1,xi2,weight)");
        for(std::size_t i=0;i<rows.size();i+=4)out.push_back({{rows[i],rows[i+1],rows[i+2]},rows[i+3]});};
    rule(input.tet,block.tet);rule(input.tri,block.tri);rule(input.tetLow,block.tetLow);rule(input.triLow,block.triLow);
    block.images=radia_hmatrix_cuda::flatImageRows(input.imageMasks,input.imageSigns,input.imageAngles);
    block.validate();return block;
}
}
void acaDiagnosticLeaf(const radia_hmatrix_cuda::HostVector<double>& diagonal,const radia_hmatrix_cuda::HostVector<double>& outer,double eps,int maxRank,int& rank,int& status,int& steps){
    const auto report=radia_hmatrix_cuda::acaDiagnosticLeaf(diagonal,outer,eps,maxRank);rank=report.rank;status=report.status;steps=report.steps;
}
std::size_t materialEnergyDiagnostic(const radia_hmatrix_cuda::HostVector<double>& energy,const radia_hmatrix_cuda::HostVector<double>& amplitude,double& sum,double& maximum){
    using namespace radia_hmatrix_cuda;
    if(energy.empty()||energy.size()!=amplitude.size())throw std::invalid_argument("Material energy diagnostic takes two matching nonempty vectors");
    const auto blocks=radia_hdiv_bh::materialEnergyReductionBlocks(energy.size());
    Buffer<double> e,a,partials,result;e.upload(energy);a.upload(amplitude);partials.allocate(3*blocks);result.allocate(2);
    radia_hdiv_bh::launchMaterialEnergyReduction(energy.size(),e.data,a.data,partials.data,result.data,nullptr);
    double out[2];check(cudaMemcpy(out,result.data,2*sizeof(double),cudaMemcpyDeviceToHost));sum=out[0];maximum=out[1];return blocks;
}
radia_hmatrix_cuda::HostVector<double> flatDiagnosticEntries(const FlatInput& input,const radia_hmatrix_cuda::HostVector<std::size_t>& rows,
    const radia_hmatrix_cuda::HostVector<std::size_t>& cols){
    if(input.host.empty())throw std::invalid_argument("Flat diagnostic entries require charges");
    return radia_hmatrix_cuda::flatDiagnosticEntries(flatBlock(input),rows,cols);
}
std::shared_ptr<DeviceOperator> createFlat(const FlatInput& input,const ChargeMapInput& charges) {
    const auto n=input.host.size();
    if(!n||!charges.dofs)throw std::invalid_argument("Unsupported bounded flat HDiv resident dimensions");
    const auto block=flatBlock(input);radia_hdiv_bh::validateFlatChargeMap(input,charges);
    const auto budget=chargeBudget(block,input,charges,input.compression=="aca"||input.compression=="hca");
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
    if(input.compression=="aca"||input.compression=="hca"){
        radia_hmatrix_cuda::FlatAcaOptions options;options.eps=input.acaEps;options.eta=input.acaEta;
        options.leafSize=input.acaLeafSize;options.maxRank=input.acaMaxRank;options.rankCapDense=input.acaRankCapDense;options.expansionOrder=input.acaExpansionOrder;
        // Ranks are known only after device compression; they size the plan.
        radia_hmatrix_cuda::FlatAcaOracle oracle(block,options);oracle.compress();
        auto plan=radia_hmatrix_cuda::packLeaves(n,oracle.leaves<LeafInput>(),oracle.permutation(),radia_hmatrix_cuda::HostVector<double>(n,1.));
        auto result=std::allocate_shared<Operator>(owner,plan,oracle,true,input,charges,budget);
        const auto& s=oracle.stats;
        result->record.counters={{"elements",s.elements},{"clusters",s.clusters},{"dense_leaves",s.denseLeaves},
            {"lowrank_leaves",s.lowrankLeaves},{"admissible_blocks",s.admissibleBlocks},{"rank_cap_promotions",s.rankCapPromotions},
            {"storage_promotions",s.storagePromotions},{"zero_rank_leaves",s.zeroRankLeaves},{"rank_sum",s.rankSum},{"max_rank",s.maxRank},{"aca_steps",s.steps},
            {"aca_batches",s.batches},{"factor_coefficients",s.factorCoefficients},{"dense_coefficients",s.denseCoefficients},
            {"far_pairs",s.farPairs},{"mid_pairs",s.midPairs},{"graded_pairs",s.gradedPairs},{"far_microseconds",s.farMicroseconds},
            {"mid_microseconds",s.midMicroseconds},{"graded_microseconds",s.gradedMicroseconds},{"aca_microseconds",s.acaMicroseconds},{"dense_microseconds",s.denseMicroseconds}};
        if(options.expansionOrder){result->record.counters.push_back({"expansion_order",s.expansionOrder});result->record.counters.push_back({"kernel_rank_sum",s.kernelRankSum});
            result->record.counters.push_back({"moment_coefficients",s.momentCoefficients});result->record.counters.push_back({"moment_microseconds",s.momentMicroseconds});}
        return result;
    }
    throw std::invalid_argument("Flat charge Gram compression must be explicit 'dense', 'aca' or 'hca'");
}
std::shared_ptr<DeviceOperator> createAnalyticGram(const radia_hmatrix_cuda::HostVector<double>& cells,const radia_hmatrix_cuda::HostVector<double>& faces){
    using namespace radia_hmatrix_cuda;
    const auto block=analyticSources(cells,faces);const auto n=block.sources.size();
    Leaf leaf{0,0,n,n,0};leaf.deferred=true;
    HostVector<std::size_t> permutation(n);for(std::size_t i=0;i<n;++i)permutation[i]=i;
    const auto plan=pack(n,{leaf},permutation,HostVector<double>(n,1.));
    return std::allocate_shared<Operator>(HostFaultAllocator<Operator>("python.operator-owner"),plan,block);
}
} // namespace radia_hdiv_python

radia_hdiv_bh::Samples radia_hdiv_bh::evaluate(const Law& law,const radia_hmatrix_cuda::HostVector<double>& a){return evaluateInverseBH(law,a);}

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
    MaterialStats stats()override{std::lock_guard<std::mutex> lock(mutex);if(context)return context->stats();MaterialStats out;out.closed=true;return out;}
    std::string route()override{std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");return context->route();}
    radia_hmatrix_cuda::HostVector<LedgerRow> transfers()override{std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");return context->transfers();}
    MaterialFootprint footprint()override{std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");return context->footprint();}
    MassAction massAction(const radia_hmatrix_cuda::HostVector<double>& x)override{std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");return context->massAction(x);}
    radia_hmatrix_cuda::HostVector<double> evaluationAction(const radia_hmatrix_cuda::HostVector<double>& x,bool transpose)override{std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");return context->evaluationAction(x,transpose);}
    void close()override{std::lock_guard<std::mutex> lock(mutex);context.reset();if(counted){--liveMaterialOwners;counted=false;}}
};
}
namespace {
std::atomic<unsigned> liveMassOwners{0};
// The same device builder as the material context; its evaluation tables and
// point measures are released after the build, only M and its diagonal stay.
class MassOwner final:public MassOperator {
    std::mutex mutex;radia_hmatrix_cuda::HostUnique<radia_hmatrix_cuda::DeviceMap> M;radia_hmatrix_cuda::HostUnique<radia_hmatrix_cuda::Buffer<double>> diagonal;
    std::size_t nd=0,nonzeros=0,massBytes=0;double bound=0;bool counted=false;std::string routeName;
    radia_hmatrix_cuda::HostVector<LedgerRow> ledger;MaterialFootprint bytes;int device=0;
    void live(){if(!M)throw std::invalid_argument("Mass context is closed");int current;radia_hmatrix_cuda::check(cudaGetDevice(&current));if(current!=device)throw std::runtime_error("Mass context belongs to another GPU");}
public:
    explicit MassOwner(const FlatEvaluationGroups& input){
        using namespace radia_hmatrix_cuda;
        if(liveMassOwners.fetch_add(1)>=2){--liveMassOwners;throw std::invalid_argument("Mass context count budget exceeded");}
        counted=true;
        try{
            validateFlatEvaluation(input);FlatEvaluationBudget budget;requireFlatEvaluationDevice(input,0,budget);
            check(cudaGetDevice(&device));nd=input[0].dofs;routeName="device-"+flatEvaluationFamily(input)+"-mass";
            auto map=makeHostUnique<DeviceMap>("python.mass-map");auto diag=makeHostUnique<Buffer<double>>("python.mass-diagonal");
            {FlatEvaluationOperator flat;Buffer<double> weights;FlatEvaluationBuilder builder;
                builder.build(input,flat,weights,*map,*diag);
                nonzeros=builder.massNonzeros;bound=builder.massBound;
                bytes.transient=builder.transientBytes+flat.bytes()+weights.size*sizeof(double);
                for(const auto& t:builder.transfers)ledger.push_back({t.payload,t.direction,t.bytes});}
            massBytes=(nd+1+nonzeros)*sizeof(std::size_t)+(nonzeros+nd)*sizeof(double);bytes.persistent=massBytes;
            M=std::move(map);diagonal=std::move(diag);
        }catch(...){--liveMassOwners;counted=false;throw;}
    }
    ~MassOwner()override{close();}
    radia_hmatrix_cuda::MassOperand operand(){live();return {M.get(),diagonal->data,bound,nd,massBytes,true};}
    std::mutex& guard(){return mutex;}
    MassAction massAction(const radia_hmatrix_cuda::HostVector<double>& x)override{
        using namespace radia_hmatrix_cuda;std::lock_guard<std::mutex> lock(mutex);live();
        if(x.size()!=nd)throw std::invalid_argument("Reference mass action dimension mismatch");
        for(double v:x)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite reference mass input");
        Buffer<double> in,out;in.upload(x);out.allocate(nd);M->apply(in.data,out.data,nullptr);check(cudaDeviceSynchronize());
        MassAction result;result.y.resize(nd);check(cudaMemcpy(result.y.data(),out.data,nd*sizeof(double),cudaMemcpyDeviceToHost));
        result.nonzeros=nonzeros;result.bound=bound;return result;
    }
    radia_hmatrix_cuda::HostVector<LedgerRow> transfers()override{std::lock_guard<std::mutex> lock(mutex);live();return ledger;}
    MaterialFootprint footprint()override{std::lock_guard<std::mutex> lock(mutex);live();return bytes;}
    std::size_t dofs()override{std::lock_guard<std::mutex> lock(mutex);return nd;}
    bool closed()override{std::lock_guard<std::mutex> lock(mutex);return !M;}
    std::string route()override{std::lock_guard<std::mutex> lock(mutex);return routeName;}
    void close()override{std::lock_guard<std::mutex> lock(mutex);diagonal.reset();M.reset();if(counted){--liveMassOwners;counted=false;}}
};
}
std::shared_ptr<MassOperator> createFlatMass(const FlatEvaluationGroups& input){return std::allocate_shared<MassOwner>(radia_hmatrix_cuda::HostFaultAllocator<MassOwner>("python.mass-owner"),input);}
std::shared_ptr<MaterialOperator> createFlatMaterial(const FlatEvaluationGroups& input,const Law& law){return std::allocate_shared<MaterialOwner>(radia_hmatrix_cuda::HostFaultAllocator<MaterialOwner>("python.material-owner"),input,law);}
FlatEvaluationBudget validateFlatEvaluationScope(const FlatEvaluationGroups& input){validateFlatEvaluation(input);return flatEvaluationBudget(input,0);}
}

radia_hdiv_python::NewtonResult radia_hdiv_python::nonlinearBridge(radia_hmatrix_cuda::Resident& op,radia_hdiv_bh::MaterialOperator& material,const radia_hmatrix_cuda::HostVector<double>& rhs,const radia_hmatrix_cuda::HostVector<double>& initial,const NewtonOptions& options){
    auto* owner=dynamic_cast<radia_hdiv_bh::MaterialOwner*>(&material);if(!owner)throw std::invalid_argument("Material owner ABI mismatch");return owner->solve(op,rhs,initial,options);
}
radia_hdiv_python::LinearResult radia_hdiv_python::linearBridge(radia_hmatrix_cuda::Resident& op,radia_hdiv_bh::MassOperator& mass,const radia_hmatrix_cuda::HostVector<double>& rhs,const LinearOptions& options){
    auto* owner=dynamic_cast<radia_hdiv_bh::MassOwner*>(&mass);if(!owner)throw std::invalid_argument("Mass owner ABI mismatch");
    if(options.preconditioner<0||options.preconditioner>5)throw std::invalid_argument("Unsupported material preconditioner");
    std::lock_guard<std::mutex> lock(owner->guard());
    const auto solution=radia_hmatrix_cuda::solveLinearMaterial(op,owner->operand(),rhs,options.invChi,options.tolerance,options.maximumIterations,
        radia_hmatrix_cuda::MassPreconditioner(options.preconditioner),options.degree,options.lowerFraction);
    LinearResult out;out.x.assign(solution.x.begin(),solution.x.end());out.iterations=solution.iterations;
    out.relativeResidual=solution.relativeResidual;out.ownedBytes=solution.ownedBytes;out.graphLaunches=solution.preconditionerGraphLaunches;return out;
}
