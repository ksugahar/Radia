// Test-only MEX entry over the reviewed kernel object, never shipped. These definitions take the place
// of the cudart import stubs when linking, so the kernel object itself is unchanged. CUDA library
// workspaces of cublas and cusolver do not pass through here and stay outside this count.
#include "mex.h"
#include "fault_state.hpp"
#include "hmatrix_host_fault.hpp"
#include <windows.h>
#include <cstring>
#include <map>
#include <string>
extern "C" void radia_hdiv_device_inner(int,mxArray*[],int,const mxArray*[]);
namespace {
constexpr int cudaSuccess=0,cudaErrorMemoryAllocation=2,cudaErrorUnknown=999;
enum Kind {None,Device,Sync,Mx,Host,Resource,Unknown};
struct State {
    Kind kind=None;bool mxPaused=false;std::size_t failAt=0,deviceAttempts=0,syncAttempts=0,mxAttempts=0,deviceLiveBytes=0,freeErrors=0;
    std::map<void*,std::size_t>* live=nullptr;
} state;
template<class F>F real(const char* name){
    static HMODULE module=GetModuleHandleA("cudart64_12.dll");
    if(!module)mexErrMsgIdAndTxt("radia:hdiv:fault","cudart64_12.dll is not loaded");
    const auto p=reinterpret_cast<F>(GetProcAddress(module,name));
    if(!p)mexErrMsgIdAndTxt("radia:hdiv:fault","Missing CUDA entry %s",name);
    return p;
}
bool hit(Kind kind,std::size_t attempts){return state.kind==kind&&state.failAt&&attempts==state.failAt;}
bool arm(const std::string& name,std::size_t at){
    using namespace radia_hmatrix_cuda;
    const Kind kind=name=="none"?None:name=="device"?Device:name=="sync"?Sync:name=="mx"?Mx:name=="host"?Host:name=="resource"?Resource:Unknown;
    if(kind==Unknown)return false;
    state.kind=kind;state.mxPaused=false;state.failAt=at;state.deviceAttempts=state.syncAttempts=state.mxAttempts=0;
    // Enabled once before any handle exists, so buffers released after a disarm stay counted.
    hostFaultState.enabled=resourceFaultState.enabled=true;
    hostFaultState.attempts=resourceFaultState.attempts=0;
    hostFaultState.failAt=kind==Host?at:0;resourceFaultState.failAt=kind==Resource?at:0;
    return true;
}
mxArray* field(double v){return mxCreateDoubleScalar(v);}
mxArray* report(){
    using namespace radia_hmatrix_cuda;
    const char* names[]={"device_attempts","sync_attempts","mx_attempts","host_attempts","resource_attempts","device_live_bytes",
        "device_live_count","free_errors","host_live_bytes","resource_live","resource_free_errors"};
    mxArray* r=mxCreateStructMatrix(1,1,11,names);
    const double values[]={double(state.deviceAttempts),double(state.syncAttempts),double(state.mxAttempts),double(hostFaultState.attempts),
        double(resourceFaultState.attempts),double(state.deviceLiveBytes),double(state.live?state.live->size():0),double(state.freeErrors),
        double(hostFaultState.liveBytes),double(resourceFaultState.live),double(resourceFaultState.freeErrors)};
    for(int i=0;i<11;++i)mxSetField(r,0,names[i],field(values[i]));
    return r;
}
}
namespace radia_hdiv_fault {
bool mxFails(){return !state.mxPaused&&hit(Mx,++state.mxAttempts);}
}
extern "C" int cudaMalloc(void** p,std::size_t bytes){
    static const auto forward=real<int(*)(void**,std::size_t)>("cudaMalloc");
    if(hit(Device,++state.deviceAttempts))return cudaErrorMemoryAllocation;
    const int status=forward(p,bytes);
    if(status==cudaSuccess&&*p){
        if(!state.live)state.live=new std::map<void*,std::size_t>();
        state.live->emplace(*p,bytes);state.deviceLiveBytes+=bytes;
    }
    return status;
}
extern "C" int cudaFree(void* p){
    static const auto forward=real<int(*)(void*)>("cudaFree");
    const int status=forward(p);
    if(p&&state.live){
        const auto it=state.live->find(p);
        if(it==state.live->end()||status!=cudaSuccess)++state.freeErrors;
        else{state.deviceLiveBytes-=it->second;state.live->erase(it);}
    }
    return status;
}
// The device work finishes before the error is reported, so cleanup never frees a buffer still in use.
extern "C" int cudaStreamSynchronize(void* stream){
    static const auto forward=real<int(*)(void*)>("cudaStreamSynchronize");
    const int status=forward(stream);
    if(hit(Sync,++state.syncAttempts))return cudaErrorUnknown;
    return status;
}
// No try block here: MATLAB raises its errors as C++ exceptions, and they must pass through this frame.
void mexFunction(int nl,mxArray* out[],int nr,const mxArray* in[]){
    if(nr>=1&&mxIsChar(in[0])){
        char command[32]={0};mxGetString(in[0],command,sizeof(command));
        if(!std::strcmp(command,"fault_arm")){
            if(nr!=3||!mxIsChar(in[1])||!mxIsDouble(in[2])||mxGetNumberOfElements(in[2])!=1)mexErrMsgIdAndTxt("radia:hdiv:fault","fault_arm(kind,at)");
            char kind[16]={0};mxGetString(in[1],kind,sizeof(kind));
            if(!arm(kind,std::size_t(mxGetScalar(in[2]))))mexErrMsgIdAndTxt("radia:hdiv:fault","Unknown fault kind %s",kind);
            return;
        }
        if(!std::strcmp(command,"fault_report")){out[0]=report();return;}
        if(!std::strcmp(command,"fault_mx_pause")){state.mxPaused=nr==2&&mxIsLogicalScalarTrue(in[1]);return;}
    }
    radia_hdiv_device_inner(nl,out,nr,in);
}
