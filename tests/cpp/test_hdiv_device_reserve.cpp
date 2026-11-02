// Guards: device buffers live at the same time are checked together against
// the GPU reserve (separate fit, combined fail is refused), with exact
// boundaries and overflow-checked sums.
// Assumes: host arithmetic only; the allocation sites are checked by the HDiv
// source contracts and the RTX smokes, not here.
#include "gpu/hdiv/hdiv_device_reserve.hpp"
#include <cstdio>
#include <limits>
#include <stdexcept>

namespace {
int failures=0;
void expect(bool ok,const char* what){
    if(!ok){std::fprintf(stderr,"FAILED: %s\n",what);++failures;}
}
}

int main(){
    using namespace radia_hmatrix_cuda;
    constexpr std::size_t MiB=std::size_t(1)<<20,GiB=std::size_t(1)<<30;
    const std::size_t total=16*GiB,reserve=reserveBytes(total);
    expect(reserve==total/5,"16 GiB device reserves total/5");
    expect(reserveBytes(4*GiB)==2*GiB,"small device reserves 2 GiB");
    const std::size_t free=reserve+150*MiB,map=100*MiB,scratch=100*MiB;
    expect(fitsAboveReserve(free,total,map),"map alone fits");
    expect(fitsAboveReserve(free,total,scratch),"scratch alone fits");
    expect(!fitsAboveReserve(free,total,checkedSum({map,scratch})),"map and scratch together must not fit");
    expect(fitsAboveReserve(free,total,checkedSum({map,50*MiB})),"exactly the free bytes above the reserve fit");
    expect(!fitsAboveReserve(free,total,checkedSum({map,50*MiB,1})),"one byte above the reserve is refused");
    expect(!fitsAboveReserve(reserve,total,0),"no free bytes above the reserve");
    bool overflow=false;
    try{checkedSum({std::numeric_limits<std::size_t>::max(),1});}catch(const std::overflow_error&){overflow=true;}
    expect(overflow,"overflowing sum throws");
    if(failures)return 1;
    std::puts("hdiv device reserve: ok");
    return 0;
}
