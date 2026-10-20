// Symbolic coefficient/serial ACA payload bounds; excludes metadata and runtime.
#pragma once
#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>
namespace radia_hmatrix_cuda {
struct LeafCoefficientBounds {
    std::size_t leaves=0,upperLeaves=0,lowerLeaves=0;
    std::size_t upperDenseCoefficients=0,upperFactorCapCoefficients=0;
    std::size_t upperStorageCapCoefficients=0,maximumSerialFactorCoefficients=0;
    std::size_t maximumSerialAcaAuxiliaryBytes=0;
    static std::size_t add(std::size_t a,std::size_t b){
        if(b>std::numeric_limits<std::size_t>::max()-a)throw std::overflow_error("Leaf bound sum overflow");return a+b;
    }
    static std::size_t multiply(std::size_t a,std::size_t b){
        if(a&&b>std::numeric_limits<std::size_t>::max()/a)throw std::overflow_error("Leaf bound product overflow");return a*b;
    }
    void include(std::size_t n,std::size_t row,std::size_t column,std::size_t rows,std::size_t columns,bool admissible,std::size_t rankCap){
        if(!n||!rows||!columns||!rankCap||row>=n||column>=n||rows>n-row||columns>n-column)
            throw std::invalid_argument("Invalid symbolic leaf extent/rank cap");
        if(row==column&&rows!=columns)throw std::invalid_argument("Non-square diagonal leaf");
        // Transactional: a rejected extent or arithmetic overflow leaves prior bounds intact.
        auto next=*this;next.leaves=add(next.leaves,1);
        if(row>column){next.lowerLeaves=add(next.lowerLeaves,1);*this=next;return;}
        const auto dense=multiply(rows,columns);
        const auto factors=admissible?multiply(rankCap,add(rows,columns)):0;
        next.upperLeaves=add(next.upperLeaves,1);
        next.upperDenseCoefficients=add(next.upperDenseCoefficients,dense);
        next.upperFactorCapCoefficients=add(next.upperFactorCapCoefficients,factors);
        // Cap exhaustion can promote an admissible block to exact dense.
        next.upperStorageCapCoefficients=add(next.upperStorageCapCoefficients,std::max(dense,factors));
        next.maximumSerialFactorCoefficients=std::max(next.maximumSerialFactorCoefficients,factors);
        if(admissible)next.maximumSerialAcaAuxiliaryBytes=std::max(next.maximumSerialAcaAuxiliaryBytes,serialAcaAuxiliaryBytes(rows,columns,rankCap));
        *this=next;
    }
    std::size_t storageCapBytes()const{return multiply(upperStorageCapCoefficients,sizeof(double));}
    // Serial fill owns completed output factors plus zaa/zab rank-cap scratch.
    // ACA masks/reference vectors coexist with that scratch. Dense promotion
    // frees scratch first; using maxima independently is conservative.
    std::size_t serialFillPayloadCapBytes()const{
        return add(add(storageCapBytes(),multiply(maximumSerialFactorCoefficients,sizeof(double))),maximumSerialAcaAuxiliaryBytes);
    }
    // Borrowed native values remain alive while the owned Plan is packed.
    // Deferred dense leaves make this conservative, not an exact live counter.
    std::size_t nativeAndPlanPayloadCapBytes()const{return multiply(storageCapBytes(),2);}
    std::size_t payloadLifecycleCapBytes()const{
        return std::max(serialFillPayloadCapBytes(),nativeAndPlanPayloadCapBytes());
    }
    // ACA masks + rank workspace + ACA row bounds OR ACA+ reference vectors.
    // Excludes output factors, callback buffers, allocator overhead and C tree/control.
    static std::size_t serialAcaAuxiliaryBytes(std::size_t rows,std::size_t columns,std::size_t rankCap){
        if(!rows||!columns||!rankCap)throw std::invalid_argument("Invalid ACA auxiliary extent");
        const auto extent=add(rows,columns);
        const auto masks=multiply(extent,sizeof(int));
        const auto values=multiply(add(rankCap,extent),sizeof(double));
        return add(masks,values);
    }
};
}
