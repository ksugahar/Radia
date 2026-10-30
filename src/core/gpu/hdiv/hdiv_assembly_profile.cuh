// Opt-in stream-event diagnostics; event/runtime storage is not owned-buffer memory.
#pragma once
struct HDivAssemblyProfile {
    struct Sample {
        cudaEvent_t first=nullptr,last=nullptr;
        int stage;
        explicit Sample(int value):stage(value){
            try{
                resourceFaultPoint("profile.first.create");
                check(cudaEventCreate(&first));resourceAcquired();
                resourceFaultPoint("profile.first.after-create");
                resourceFaultPoint("profile.last.create");
                check(cudaEventCreate(&last));resourceAcquired();
                resourceFaultPoint("profile.last.after-create");
            }catch(...){release();throw;}
        }
        Sample(const Sample&)=delete;
        Sample& operator=(const Sample&)=delete;
        void release()noexcept{
            if(last){resourceReleased(cudaEventDestroy(last)==cudaSuccess);last=nullptr;}
            if(first){resourceReleased(cudaEventDestroy(first)==cudaSuccess);first=nullptr;}
        }
        ~Sample(){release();}
    };
    bool enabled;
    cudaStream_t stream;
    HostVector<HostUnique<Sample>> samples{HostFaultAllocator<HostUnique<Sample>>("profile.samples")};
    HDivAssemblyProfile(bool value,cudaStream_t owner):enabled(value),stream(owner){}
    void begin(int stage){if(enabled){samples.emplace_back(makeHostUnique<Sample>("profile.sample",stage));check(cudaEventRecord(samples.back()->first,stream));}}
    void end(){if(enabled)check(cudaEventRecord(samples.back()->last,stream));}
    // Caller synchronizes this stream after recording every sample's end.
    void collect(double* seconds){if(enabled)for(const auto& sample:samples){float ms=0;check(cudaEventElapsedTime(&ms,sample->first,sample->last));seconds[sample->stage]+=double(ms)*0.001;}}
};
