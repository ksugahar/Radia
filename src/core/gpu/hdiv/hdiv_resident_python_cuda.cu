// Same packing, scaling, sparse maps and kernels as the standalone MEX gateway.
#include "hdiv_resident_python_api.hpp"
#include "hmatrix_material_solver.cuh"
#include "hmatrix_nonlinear_material.cuh"
#include "hdiv_constant_diagonal.cuh"
#include <mutex>
#include <atomic>
namespace radia_hdiv_python {
NewtonResult nonlinearBridge(radia_hmatrix_cuda::Resident&,radia_hdiv_bh::MaterialOperator&,const MapInput&,const radia_hmatrix_cuda::HostVector<double>&,const radia_hmatrix_cuda::HostVector<double>&,const NewtonOptions&);

namespace {
class Operator final:public DeviceOperator {
    std::mutex mutex;
    radia_hmatrix_cuda::HostUnique<radia_hmatrix_cuda::Resident> resident;
    const std::size_t dimension;
public:
    Operator(const radia_hmatrix_cuda::Plan& plan,const radia_hmatrix_cuda::SparseMap& map)
      :dimension(map.columns) {
        // Keep the existing 2 GiB combined native setup bound, not a new cap.
        resident=radia_hmatrix_cuda::makeHostUnique<radia_hmatrix_cuda::Resident>("python.resident-object",plan,&map,1);
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
    NewtonResult solve_nonlinear(radia_hdiv_bh::MaterialOperator& material,const MapInput& reference,const radia_hmatrix_cuda::HostVector<double>& rhs,const radia_hmatrix_cuda::HostVector<double>& initial,const NewtonOptions& options)override{
        std::lock_guard<std::mutex> lock(mutex);if(!resident)throw std::invalid_argument("HDiv resident has been closed");return nonlinearBridge(*resident,material,reference,rhs,initial,options);
    }
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
} // namespace radia_hdiv_python

radia_hdiv_bh::Samples radia_hdiv_bh::evaluate(const Law& law,const radia_hmatrix_cuda::HostVector<double>& a){return evaluateInverseBH(law,a);}

radia_hdiv_bh::MaterialEvaluation radia_hdiv_bh::evaluateMaterial(const MaterialInput& input,const Law& law,const radia_hmatrix_cuda::HostVector<double>& state,const radia_hmatrix_cuda::HostVector<double>& direction){return evaluateMappedMaterial(input,law,state,direction);}

namespace radia_hdiv_bh {
namespace {
std::atomic<unsigned> liveMaterialOwners{0};
class MaterialOwner final:public MaterialOperator {
    std::mutex mutex;radia_hmatrix_cuda::HostUnique<MappedMaterial> context;bool counted=false;
public:
    MaterialOwner(const MaterialInput& input,const Law& law){
        if(liveMaterialOwners.fetch_add(1)>=2){--liveMaterialOwners;throw std::invalid_argument("Material context count budget exceeded");}
        try{context=radia_hmatrix_cuda::makeHostUnique<MappedMaterial>("python.material-object",input,law);counted=true;}catch(...){--liveMaterialOwners;throw;}
    }
    ~MaterialOwner()override{close();}
    MaterialEvaluation evaluate(const radia_hmatrix_cuda::HostVector<double>& x,const radia_hmatrix_cuda::HostVector<double>& d)override{
        std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");return context->evaluate(x,d);
    }
    MaterialEvaluation linearize(const radia_hmatrix_cuda::HostVector<double>& x)override{std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");return context->linearize(x);}
    radia_hmatrix_cuda::HostVector<double> tangent(const radia_hmatrix_cuda::HostVector<double>& d,std::size_t epoch)override{std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");return context->tangentAction(d,epoch);}
    radia_hdiv_python::NewtonResult solve(radia_hmatrix_cuda::Resident& op,const radia_hdiv_python::MapInput& input,const radia_hmatrix_cuda::HostVector<double>& rhs,const radia_hmatrix_cuda::HostVector<double>& initial,const radia_hdiv_python::NewtonOptions& options){
        std::lock_guard<std::mutex> lock(mutex);if(!context)throw std::invalid_argument("Material context is closed");
        if(input.rows!=op.dimension()||input.columns!=op.dimension())throw std::invalid_argument("Reference mass dimension mismatch");
        radia_hmatrix_cuda::SparseMap reference;reference.rows=input.rows;reference.columns=input.columns;reference.offsets.assign(input.offsets.begin(),input.offsets.end());reference.indices.assign(input.indices.begin(),input.indices.end());reference.values.assign(input.values.begin(),input.values.end());reference.validate();
        return radia_hmatrix_cuda::solveNonlinearMaterial(op,*context,reference,rhs,initial,options);
    }
    MaterialStats stats()override{std::lock_guard<std::mutex> lock(mutex);if(context)return context->stats();MaterialStats out;out.closed=true;return out;}
    void close()override{std::lock_guard<std::mutex> lock(mutex);context.reset();if(counted){--liveMaterialOwners;counted=false;}}
};
}
std::shared_ptr<MaterialOperator> createMaterial(const MaterialInput& input,const Law& law){return std::allocate_shared<MaterialOwner>(radia_hmatrix_cuda::HostFaultAllocator<MaterialOwner>("python.material-owner"),input,law);}
}

radia_hdiv_python::NewtonResult radia_hdiv_python::nonlinearBridge(radia_hmatrix_cuda::Resident& op,radia_hdiv_bh::MaterialOperator& material,const MapInput& reference,const radia_hmatrix_cuda::HostVector<double>& rhs,const radia_hmatrix_cuda::HostVector<double>& initial,const NewtonOptions& options){
    auto* owner=dynamic_cast<radia_hdiv_bh::MaterialOwner*>(&material);if(!owner)throw std::invalid_argument("Material owner ABI mismatch");return owner->solve(op,reference,rhs,initial,options);
}
