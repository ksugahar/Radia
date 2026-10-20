#pragma once
#include <cstddef>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <vector>
#include <map>
#include <cstring>
#include <utility>
namespace radia_hmatrix_cuda {
struct HostFaultState {
    bool enabled=false;
    std::size_t failAt=0, attempts=0, liveBytes=0, peakBytes=0;
    const char* site="";
    char filter[96]={0};
};
inline thread_local HostFaultState hostFaultState;
inline thread_local void (*frameExceptionProbe)()=nullptr;
// Per-module diagnostic state; never hooks the process CUDA allocator.
struct DeviceFaultState {
    bool enabled=false;
    std::size_t failAt=0,attempts=0,liveBytes=0,peakBytes=0,freeErrors=0;
};
inline thread_local DeviceFaultState deviceFaultState;
inline void configureDeviceFault(std::size_t at){
    if(deviceFaultState.liveBytes)throw std::runtime_error("Cannot reset fault counters with live device buffers");
    deviceFaultState={true,at,0,0,0,0};
}

inline void configureHostFault(std::size_t at){
    if(hostFaultState.liveBytes)throw std::runtime_error("Cannot reset fault counters with live host buffers");
    hostFaultState={true,at,0,0,0,""};
}
inline void configureHostFaultSite(std::size_t at,const char* site){
    if(!site||std::strlen(site)>=sizeof(hostFaultState.filter))throw std::invalid_argument("Invalid host allocation site filter");
    configureHostFault(at);std::strcpy(hostFaultState.filter,site);
}
// Own resource handles only; opaque CUDA library workspaces remain sampled metrics.
struct ResourceFaultState {
    bool enabled=false;
    std::size_t failAt=0,attempts=0,live=0,freeErrors=0;
    const char* site="";
};
inline thread_local ResourceFaultState resourceFaultState;
inline void configureResourceFault(std::size_t at){
    if(resourceFaultState.live)throw std::runtime_error("Cannot reset resource diagnostics with live handles");
    resourceFaultState={true,at,0,0,0,""};
}
inline void resourceFaultPoint(const char* site){
    auto& s=resourceFaultState;if(!s.enabled)return;s.site=site;++s.attempts;
    if(s.failAt && s.attempts==s.failAt)throw std::bad_alloc();
}
inline void resourceAcquired(){if(resourceFaultState.enabled)++resourceFaultState.live;}
inline void resourceReleased(bool success)noexcept{
    if(resourceFaultState.enabled){if(success && resourceFaultState.live)--resourceFaultState.live;else ++resourceFaultState.freeErrors;}
}

template<class T> struct HostFaultAllocator {
    using value_type=T;
    using propagate_on_container_move_assignment=std::true_type;
    const char* site="host-vector";
    bool counted=hostFaultState.enabled;
    HostFaultAllocator() noexcept=default;
    explicit HostFaultAllocator(const char* name) noexcept:site(name){}
    template<class U>HostFaultAllocator(const HostFaultAllocator<U>& a) noexcept:site(a.site),counted(a.counted){}
    T* allocate(std::size_t n){
        if(n>std::numeric_limits<std::size_t>::max()/sizeof(T))throw std::bad_alloc();
        const auto bytes=n*sizeof(T);
        if(counted){
            auto& s=hostFaultState;
            const bool selected=!s.filter[0]||std::strcmp(site,s.filter)==0;
            if(selected){s.site=site;++s.attempts;if(s.failAt && s.attempts==s.failAt)throw std::bad_alloc();}
            if(bytes>std::numeric_limits<std::size_t>::max()-s.liveBytes)throw std::bad_alloc();
        }
        T* p=std::allocator<T>{}.allocate(n);
        if(counted){auto& s=hostFaultState;s.liveBytes+=bytes;if(s.liveBytes>s.peakBytes)s.peakBytes=s.liveBytes;}
        return p;
    }
    void deallocate(T* p,std::size_t n) noexcept{
        if(counted)hostFaultState.liveBytes-=n*sizeof(T);
        std::allocator<T>{}.deallocate(p,n);
    }
    template<class U>bool operator==(const HostFaultAllocator<U>& other)const noexcept{return counted==other.counted;}
    template<class U>bool operator!=(const HostFaultAllocator<U>& other)const noexcept{return !(*this==other);}
};
template<class T>struct HostObjectDeleter {
    HostFaultAllocator<T> allocator;
    void operator()(T* value)noexcept{
        if(value){std::allocator_traits<HostFaultAllocator<T>>::destroy(allocator,value);allocator.deallocate(value,1);}
    }
};
template<class T>using HostUnique=std::unique_ptr<T,HostObjectDeleter<T>>;
template<class T,class... Args>HostUnique<T> makeHostUnique(const char* site,Args&&... args){
    HostFaultAllocator<T> allocator(site);
    auto* value=allocator.allocate(1);
    try{std::allocator_traits<HostFaultAllocator<T>>::construct(allocator,value,std::forward<Args>(args)...);}
    catch(...){allocator.deallocate(value,1);throw;}
    return HostUnique<T>(value,HostObjectDeleter<T>{allocator});
}
template<class T>using HostVector=std::vector<T,HostFaultAllocator<T>>;
template<class Key,class Value>using HostMap=std::map<Key,Value,std::less<Key>,HostFaultAllocator<std::pair<const Key,Value>>>;
}
