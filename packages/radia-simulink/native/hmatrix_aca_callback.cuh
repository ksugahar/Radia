// Module-local C ABI for the existing HACApK active-vector batch dispatch.
#pragma once
#include "hmatrix_entry_oracle.cuh"
#include <exception>

#include "hmatrix_leaf_bounds.h"
#include "cHACApK_base.h"
#include "cHACApK_frame_admission.h"
#define RAD_HACAPK_BATCH_ENTRIES
#include "cHACApK_calc_entry_ij.h"
namespace radia_hmatrix_cuda {
struct EntryCallbackBinding {
    void* oracle;
    HostVector<double> (*evaluate)(void*,const HostVector<AnalyticEntryIndex>&,cudaStream_t);
    cudaStream_t stream;
    std::exception_ptr failure;
    const HostVector<double>* physicalScale;
    bool deferDense;
    LeafCoefficientBounds frameBounds;
    bool frameAdmitted=false;
    std::size_t frameCoefficientBudget=std::size_t(2)<<30;
    std::size_t frameWorkingPayloadBudget=std::size_t(2)<<30;
};
inline thread_local EntryCallbackBinding* activeEntryBinding=nullptr;
class EntryCallbackScope {
    EntryCallbackBinding binding;
public:
    template<class Oracle>
    EntryCallbackScope(Oracle& oracle,cudaStream_t stream,const HostVector<double>* scale=nullptr,bool deferDense=false,std::size_t frameBudget=std::size_t(2)<<30,std::size_t workingBudget=std::size_t(2)<<30)
        :binding{&oracle,[](void* context,const HostVector<AnalyticEntryIndex>& request,cudaStream_t s){return static_cast<Oracle*>(context)->evaluate(request,s);},stream,{},scale,deferDense}{
        if(!frameBudget||frameBudget>(std::size_t(2)<<30))throw std::invalid_argument("Invalid coefficient frame allowance");
        binding.frameCoefficientBudget=frameBudget;
        if(!workingBudget||workingBudget>(std::size_t(2)<<30))throw std::invalid_argument("Invalid frame working payload allowance");
        binding.frameWorkingPayloadBudget=workingBudget;
        if(activeEntryBinding)throw std::runtime_error("Nested HACApK CUDA entry binding");
        activeEntryBinding=&binding;
    }
    ~EntryCallbackScope(){activeEntryBinding=nullptr;}
    EntryCallbackScope(const EntryCallbackScope&)=delete;
    EntryCallbackScope& operator=(const EntryCallbackScope&)=delete;
    void rethrow()const{if(binding.failure)std::rethrow_exception(binding.failure);}
    const LeafCoefficientBounds& admittedFrame()const{rethrow();if(!binding.frameAdmitted)throw std::runtime_error("HACApK frame was not admitted");return binding.frameBounds;}
};
}
#ifdef RAD_HACAPK_FRAME_ADMISSION
extern "C" int cHACApK_frame_admission(st_cHACApK_leafmtx const* leaves,int count,int dimension,int rankCap){
    auto* binding=radia_hmatrix_cuda::activeEntryBinding;
    try{
#ifdef RADIA_HMATRIX_FAULT_INJECTION
        if(radia_hmatrix_cuda::frameExceptionProbe)radia_hmatrix_cuda::frameExceptionProbe();
#endif
        if(!binding||binding->failure||binding->frameAdmitted||!binding->physicalScale||!leaves||count<1||dimension<1||rankCap<1)
            throw std::invalid_argument("Invalid HACApK frame admission state");
        const auto n=binding->physicalScale->size();
        if(std::size_t(dimension)!=n||std::size_t(count)>radia_hmatrix_cuda::checkedProduct(n,n))throw std::invalid_argument("HACApK frame count/dimension mismatch");
        radia_hmatrix_cuda::LeafCoefficientBounds bounds;
        for(int i=1;i<=count;++i){
            const auto* leaf=leaves[i];
            if(!leaf||leaf->a1||leaf->a2||leaf->kt||leaf->nstrtl<1||leaf->nstrtt<1||leaf->ndl<1||leaf->ndt<1||(leaf->ltmtx!=1&&leaf->ltmtx!=2))
                throw std::invalid_argument("HACApK frame must contain unfilled valid leaves");
            bounds.include(n,std::size_t(leaf->nstrtl-1),std::size_t(leaf->nstrtt-1),std::size_t(leaf->ndl),std::size_t(leaf->ndt),leaf->ltmtx==1,std::size_t(rankCap));
        }
        // Coefficients only. Existing pre-frame, oracle, Resident and reserve
        // checks still govern their separate overlapping allocation phases.
        if(bounds.storageCapBytes()>binding->frameCoefficientBudget)throw std::runtime_error("HACApK coefficient frame allowance exceeded");
        // Additional unchanged 2-GiB working-payload ceiling: completed native
        // coefficients, serial ACA scratch and native/owned-Plan coexistence.
        // Not total host admission: tree/control/packing metadata, allocator,
        // oracle/cache and later device/runtime phases have separate ownership.
        if(bounds.payloadLifecycleCapBytes()>binding->frameWorkingPayloadBudget)throw std::runtime_error("HACApK coefficient/ACA working payload allowance exceeded");
        binding->frameBounds=bounds;binding->frameAdmitted=true;return 1;
    }catch(...){if(binding&&!binding->failure)binding->failure=std::current_exception();return 0;}
}
#endif
#ifdef RAD_HACAPK_DENSE_LEAF_SINK
extern "C" int cHACApK_dense_leaf_sink(const st_cHACApK_leafmtx_t* leaf,const int* rows,const int* columns,int problem){
    auto* binding=radia_hmatrix_cuda::activeEntryBinding;if(!binding||!binding->deferDense)return 0;
    try{
        if(binding->failure)std::rethrow_exception(binding->failure);
        const auto n=binding->physicalScale?binding->physicalScale->size():0;
        if(!leaf||!rows||!columns||problem||!n||leaf->ltmtx!=2||leaf->kt||leaf->a1||leaf->a2||
            leaf->ndl<=0||leaf->ndt<=0||leaf->nstrtl<1||leaf->nstrtt<1||
            std::size_t(leaf->nstrtl)>n||std::size_t(leaf->nstrtt)>n||
            std::size_t(leaf->ndl)>n-std::size_t(leaf->nstrtl)+1||std::size_t(leaf->ndt)>n-std::size_t(leaf->nstrtt)+1)
            throw std::invalid_argument("Invalid deferred native leaf");
    }catch(...){if(!binding->failure)binding->failure=std::current_exception();}
    return 1;
}
#endif
extern "C" void cHACApK_entry_batch(int count,const int* rows,const int* columns,int problem,double* output){
    auto* binding=radia_hmatrix_cuda::activeEntryBinding;
    // No exception may cross the C compression frames. The owner must rethrow
    // before accepting leaves. A failed declared CUDA route never uses CPU.
    try{
        if(!binding)throw std::runtime_error("HACApK CUDA entry binding missing");
        if(binding->failure)std::rethrow_exception(binding->failure);
        if(count<1||count>RAD_HACAPK_BATCH_SIZE||problem!=0||!rows||!columns||!output)
            throw std::invalid_argument("Unsupported HACApK analytic batch request");
        radia_hmatrix_cuda::HostVector<radia_hmatrix_cuda::AnalyticEntryIndex> indices{radia_hmatrix_cuda::HostFaultAllocator<radia_hmatrix_cuda::AnalyticEntryIndex>("oracle.request")};
        for(int i=0;i<count;++i){
            if(rows[i]<1||columns[i]<1)throw std::invalid_argument("HACApK indices are one-based");
            indices.push_back({std::size_t(rows[i]-1),std::size_t(columns[i]-1)});
        }
        auto values=binding->evaluate(binding->oracle,indices,binding->stream);
        if(binding->physicalScale){
            const auto& scale=*binding->physicalScale;
            for(int i=0;i<count;++i){
                if(indices[i].row>=scale.size()||indices[i].column>=scale.size())throw std::invalid_argument("ACA physical scale shape mismatch");
                const auto rowScale=scale[indices[i].row],columnScale=scale[indices[i].column];
                if(!std::isfinite(rowScale)||!std::isfinite(columnScale)||rowScale<=0||columnScale<=0)
                    throw std::invalid_argument("Invalid ACA physical scale");
                values[i]*=(1./rowScale)*(1./columnScale);
            }
        }
        std::copy(values.begin(),values.end(),output);
    }catch(...){
        if(binding&&!binding->failure)binding->failure=std::current_exception();
        if(output&&count>0&&count<=RAD_HACAPK_BATCH_SIZE)
            std::fill(output,output+count,std::numeric_limits<double>::quiet_NaN());
    }
}
extern "C" double cHACApK_entry_ij(int row,int column,int problem){
    double value=0.;cHACApK_entry_batch(1,&row,&column,problem,&value);return value;
}
