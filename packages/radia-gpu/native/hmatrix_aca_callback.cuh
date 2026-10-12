// Module-local C ABI for the existing HACApK active-vector batch dispatch.
#pragma once
#include "hmatrix_entry_oracle.cuh"
#include <exception>
#define RAD_HACAPK_BATCH_ENTRIES
#include "cHACApK_calc_entry_ij.h"
namespace radia_hmatrix_cuda {
struct EntryCallbackBinding {
    AnalyticEntryOracle& oracle;
    cudaStream_t stream;
    std::exception_ptr failure;
};
inline thread_local EntryCallbackBinding* activeEntryBinding=nullptr;
class EntryCallbackScope {
    EntryCallbackBinding binding;
public:
    EntryCallbackScope(AnalyticEntryOracle& oracle,cudaStream_t stream):binding{oracle,stream,{}}{
        if(activeEntryBinding)throw std::runtime_error("Nested HACApK CUDA entry binding");
        activeEntryBinding=&binding;
    }
    ~EntryCallbackScope(){activeEntryBinding=nullptr;}
    EntryCallbackScope(const EntryCallbackScope&)=delete;
    EntryCallbackScope& operator=(const EntryCallbackScope&)=delete;
    void rethrow()const{if(binding.failure)std::rethrow_exception(binding.failure);}
};
}
extern "C" void cHACApK_entry_batch(int count,const int* rows,const int* columns,int problem,double* output){
    auto* binding=radia_hmatrix_cuda::activeEntryBinding;
    // No exception may cross the C compression frames. The owner must rethrow
    // before accepting leaves. A failed declared CUDA route never uses CPU.
    try{
        if(!binding)throw std::runtime_error("HACApK CUDA entry binding missing");
        if(binding->failure)std::rethrow_exception(binding->failure);
        if(count<1||count>256||problem!=0||!rows||!columns||!output)
            throw std::invalid_argument("Unsupported HACApK analytic batch request");
        std::vector<radia_hmatrix_cuda::AnalyticEntryIndex> indices;
        for(int i=0;i<count;++i){
            if(rows[i]<1||columns[i]<1)throw std::invalid_argument("HACApK indices are one-based");
            indices.push_back({std::size_t(rows[i]-1),std::size_t(columns[i]-1)});
        }
        const auto values=binding->oracle.evaluate(indices,binding->stream);
        std::copy(values.begin(),values.end(),output);
    }catch(...){
        if(binding&&!binding->failure)binding->failure=std::current_exception();
        if(output&&count>0&&count<=256)
            std::fill(output,output+count,std::numeric_limits<double>::quiet_NaN());
    }
}
