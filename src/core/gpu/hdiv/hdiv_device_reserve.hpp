#pragma once
#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <stdexcept>

namespace radia_hmatrix_cuda {
// Live GPU reserve rule of the device allocation sites. Buffers that are live
// at the same time are summed with overflow checks and must fit together above
// max(2 GiB, total/5): separate checks against one free figure can each pass
// while their sum crosses the reserve. Host-only, so it is tested without a GPU.
inline std::size_t reserveBytes(std::size_t total){return std::max(std::size_t(2)<<30,total/5);}
inline std::size_t checkedSum(std::initializer_list<std::size_t> parts){
    std::size_t sum=0;
    for(auto part:parts){
        if(part>std::numeric_limits<std::size_t>::max()-sum)throw std::overflow_error("Device allocation overflow");
        sum+=part;
    }
    return sum;
}
inline bool fitsAboveReserve(std::size_t free,std::size_t total,std::size_t bytes){
    const auto reserve=reserveBytes(total);
    return free>reserve&&bytes<=free-reserve;
}
}
