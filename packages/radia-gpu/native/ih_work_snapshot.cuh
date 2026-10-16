// Validation-only accepted-work banks. No CUDA or allocation in commit().
#pragma once
#include "ih_complex_solver.cuh"
#include "ih_esim_cell.cuh"
#include <vector>
#include <cstdint>
#include <string>
#include <chrono>
namespace radia_ih_cuda {
class WorkSnapshot {
public:
    enum class Kind {Real,Complex,Cells,Parameters};
    struct Slot {std::string name;std::size_t count,offset,bytes;Kind kind;};
    struct Meta {std::size_t step=0;double current=0,ambient=0,linear=0,certified=0,field=0,stageSeconds=0,gpuCopySeconds=0;};
private:
    std::vector<Slot> slots_;void* bank_[2]={nullptr,nullptr};std::size_t bankBytes_=0;
    unsigned accepted_=0;bool valid_=false;std::uint64_t epoch_=0;Meta metadata_[2];
    std::chrono::steady_clock::time_point stageStart_;
    cudaEvent_t start_=nullptr,stop_=nullptr;
public:
    WorkSnapshot()=default;
    WorkSnapshot(const WorkSnapshot&)=delete;WorkSnapshot& operator=(const WorkSnapshot&)=delete;
    ~WorkSnapshot(){for(auto p:bank_)if(p)cudaFree(p);if(start_)cudaEventDestroy(start_);if(stop_)cudaEventDestroy(stop_);}
    void add(const char* name,std::size_t count,Kind kind){
        const std::size_t width=kind==Kind::Real?sizeof(double):kind==Kind::Complex?sizeof(cuDoubleComplex):kind==Kind::Cells?sizeof(ESIMCellResult):sizeof(ESIMCellParameters);
        if(count>SIZE_MAX/width||bankBytes_>SIZE_MAX-15)throw std::overflow_error("Snapshot shape overflow");
        const auto offset=(bankBytes_+15)&~std::size_t(15);const auto bytes=count*width;
        if(bytes>SIZE_MAX-offset)throw std::overflow_error("Snapshot byte overflow");
        slots_.push_back({name,count,offset,bytes,kind});bankBytes_=offset+bytes;
    }
    void allocate(std::size_t otherRequired){
        if(!bankBytes_||bankBytes_>(2ULL<<30)||bankBytes_>SIZE_MAX/2)throw std::overflow_error("Snapshot direct payload budget exceeded");
        std::size_t available,total;cudaCheck(cudaMemGetInfo(&available,&total));const auto reserve=std::max<std::size_t>(2ULL<<30,total/5);
        if(available<=reserve||2*bankBytes_>available-reserve||otherRequired>available-reserve-2*bankBytes_)throw std::runtime_error("Atomicity and snapshot banks exceed GPU reserve");
        cudaCheck(cudaMalloc(&bank_[0],bankBytes_));cudaCheck(cudaMalloc(&bank_[1],bankBytes_));
        cudaCheck(cudaEventCreate(&start_));cudaCheck(cudaEventCreate(&stop_));
    }
    void begin(cudaStream_t stream){stageStart_=std::chrono::steady_clock::now();cudaCheck(cudaEventRecord(start_,stream));}
    void stage(const char* name,const void* source,cudaStream_t stream){
        for(const auto& slot:slots_)if(slot.name==name){if(slot.bytes&&!source)throw std::invalid_argument("Missing snapshot source");if(slot.bytes)cudaCheck(cudaMemcpyAsync(static_cast<char*>(bank_[1-accepted_])+slot.offset,source,slot.bytes,cudaMemcpyDeviceToDevice,stream));return;}
        throw std::invalid_argument("Unknown snapshot field");
    }
    void end(cudaStream_t stream,Meta meta){metadata_[1-accepted_]=meta;cudaCheck(cudaEventRecord(stop_,stream));}
    void finish(){float ms;cudaCheck(cudaEventElapsedTime(&ms,start_,stop_));metadata_[1-accepted_].gpuCopySeconds=double(ms)*1e-3;metadata_[1-accepted_].stageSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-stageStart_).count();}
    void commit(bool resetting)noexcept{if(resetting){valid_=false;++epoch_;}else{accepted_=1-accepted_;valid_=true;}}
    static_assert(noexcept(std::declval<WorkSnapshot&>().commit(false)),"Snapshot commit must not throw");
    void abandon()noexcept{bank_[0]=bank_[1]=nullptr;start_=stop_=nullptr;valid_=false;}
    bool valid()const noexcept{return valid_;}std::uint64_t epoch()const noexcept{return epoch_;}
    std::size_t bytes()const noexcept{return 2*bankBytes_;}
    const Meta& meta()const{return metadata_[accepted_];}
    const std::vector<Slot>& slots()const noexcept{return slots_;}
    template<class T>std::vector<T> read(const Slot& slot)const{
        if(!valid_||slot.bytes!=slot.count*sizeof(T))throw std::invalid_argument("Invalid snapshot read");
        std::vector<T> out(slot.count);if(slot.bytes)cudaCheck(cudaMemcpy(out.data(),static_cast<const char*>(bank_[accepted_])+slot.offset,slot.bytes,cudaMemcpyDeviceToHost));return out;
    }
};
}
