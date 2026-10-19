// Analytical symbolic bounds and rejection atomicity; no matrix allocation.
#include "hmatrix_leaf_bounds.h"
#include <iostream>
using radia_hmatrix_cuda::LeafCoefficientBounds;
void require(bool value){if(!value)throw std::runtime_error("Leaf coefficient bound contract failed");}
int main(){try{
    LeafCoefficientBounds b;b.include(6,0,0,3,3,false,2);b.include(6,0,3,3,3,true,2);
    b.include(6,3,0,3,3,true,2);b.include(6,3,3,3,3,false,2);
    require(b.leaves==4&&b.upperLeaves==3&&b.lowerLeaves==1);
    require(b.upperDenseCoefficients==27&&b.upperFactorCapCoefficients==12&&b.upperStorageCapCoefficients==30&&b.maximumSerialFactorCoefficients==12&&b.storageCapBytes()==240);
    require(LeafCoefficientBounds::serialAcaAuxiliaryBytes(3,5,2)==112);
    bool auxOverflow=false;try{LeafCoefficientBounds::serialAcaAuxiliaryBytes(std::numeric_limits<std::size_t>::max(),1,1);}catch(const std::overflow_error&){auxOverflow=true;}require(auxOverflow);
    bool auxInvalid=false;try{LeafCoefficientBounds::serialAcaAuxiliaryBytes(0,5,2);}catch(const std::invalid_argument&){auxInvalid=true;}require(auxInvalid);
    require(b.maximumSerialAcaAuxiliaryBytes==88);
    require(b.serialFillPayloadCapBytes()==424&&b.nativeAndPlanPayloadCapBytes()==480&&b.payloadLifecycleCapBytes()==480);
    auto before=b;int rejected=0;
    try{b.include(6,5,0,2,1,false,1);}catch(const std::invalid_argument&){++rejected;}
    try{b.include(6,0,0,2,3,false,1);}catch(const std::invalid_argument&){++rejected;}
    try{b.include(6,0,0,1,1,false,0);}catch(const std::invalid_argument&){++rejected;}
    try{b.include(std::numeric_limits<std::size_t>::max(),0,0,std::numeric_limits<std::size_t>::max(),std::numeric_limits<std::size_t>::max(),false,1);}catch(const std::overflow_error&){++rejected;}
    require(rejected==4&&b.leaves==before.leaves&&b.upperStorageCapCoefficients==before.upperStorageCapCoefficients&&b.maximumSerialAcaAuxiliaryBytes==before.maximumSerialAcaAuxiliaryBytes);
    LeafCoefficientBounds promoted;promoted.include(3,0,0,3,3,true,1);
    require(promoted.upperStorageCapCoefficients==9&&promoted.upperFactorCapCoefficients==6);
    require(promoted.serialFillPayloadCapBytes()==200&&promoted.nativeAndPlanPayloadCapBytes()==144&&promoted.payloadLifecycleCapBytes()==200);
    LeafCoefficientBounds boundary;boundary.upperStorageCapCoefficients=(std::size_t(1)<<30)/sizeof(double);
    require(boundary.payloadLifecycleCapBytes()==(std::size_t(2)<<30));
    boundary.upperStorageCapCoefficients++;
    require(boundary.payloadLifecycleCapBytes()>(std::size_t(2)<<30));
    LeafCoefficientBounds overlap;overlap.upperStorageCapCoefficients=std::numeric_limits<std::size_t>::max()/sizeof(double);
    bool overlapOverflow=false;try{overlap.nativeAndPlanPayloadCapBytes();}catch(const std::overflow_error&){overlapOverflow=true;}require(overlapOverflow);
    LeafCoefficientBounds scratch;scratch.maximumSerialAcaAuxiliaryBytes=std::numeric_limits<std::size_t>::max();scratch.maximumSerialFactorCoefficients=1;
    bool scratchOverflow=false;try{scratch.serialFillPayloadCapBytes();}catch(const std::overflow_error&){scratchOverflow=true;}require(scratchOverflow);
    LeafCoefficientBounds bytes;bytes.upperStorageCapCoefficients=std::numeric_limits<std::size_t>::max();
    bool overflow=false;try{bytes.storageCapBytes();}catch(const std::overflow_error&){overflow=true;}require(overflow);
    std::cout<<"{\"accepted\":true,\"upper_dense_coefficients\":27,\"upper_storage_cap_coefficients\":30,\"lower_leaves_skipped\":1,\"invalid_or_overflow_rejected\":5,\"rejection_preserves_prior_bounds\":true,\"dense_promotion_budgeted\":true}"<<std::endl;
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
