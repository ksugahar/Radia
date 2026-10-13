// Host scheduling boundary. MATLAB builds never include NGSolve TaskManager.
#pragma once
#include <cstddef>
#ifdef RAD_MATLAB_HOST
#include <algorithm>
#include <exception>
#include <mutex>
#ifdef _OPENMP
#include <omp.h>
#endif
namespace radia { namespace host {
inline thread_local int requestedThreads=1;
struct IntRange {std::size_t size;explicit IntRange(std::size_t n):size(n){}};
inline int ThreadCount(){return requestedThreads;}
inline int ThreadId(){
#ifdef _OPENMP
 return omp_in_parallel()?omp_get_thread_num():0;
#else
 return 0;
#endif
}
struct ThreadScope {int previous;explicit ThreadScope(int n):previous(requestedThreads){
#ifdef _OPENMP
 requestedThreads=std::max(1,n);
#else
 (void)n;requestedThreads=1;
#endif
}~ThreadScope(){requestedThreads=previous;}};
template<class F>void ParallelFor(IntRange range,F&& function){
#ifdef _OPENMP
 std::exception_ptr failure;std::mutex mutex;const int count=requestedThreads;
#pragma omp parallel for num_threads(count) schedule(static)
 for(long long i=0;i<static_cast<long long>(range.size);++i){try{function(static_cast<std::size_t>(i));}catch(...){std::lock_guard<std::mutex> lock(mutex);if(!failure)failure=std::current_exception();}}
 if(failure)std::rethrow_exception(failure);
#else
 for(std::size_t i=0;i<range.size;++i)function(i);
#endif
}
}}
#else
#include <core/taskmanager.hpp>
namespace radia { namespace host {
using ngcore::IntRange;
inline int ThreadCount(){return ngcore::TaskManager::GetMaxThreads();}
inline int ThreadId(){return ngcore::TaskManager::GetThreadId();}
using ThreadScope=ngcore::RegionTaskManager;
template<class Range,class F>void ParallelFor(Range range,F&& function){ngcore::ParallelFor(range,function);}
}}
#endif
