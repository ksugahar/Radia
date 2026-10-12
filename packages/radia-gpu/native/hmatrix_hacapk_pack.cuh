// Direct in-process HACApK packing: never exposes pointers to MATLAB.
#pragma once
#include "hmatrix_resident.cuh"
extern "C" {
#include "cHACApK_base.h"
}
namespace radia_hmatrix_cuda {
inline Plan fromHacapk(const st_cHACApK_leafmtxp_t& matrix,
                      const st_cHACApK_lcontrol_t& control,
                      const std::vector<double>& physicalScale) {
    if(matrix.nd<=0||matrix.nlf<=0||!matrix.st_lf||!control.lod)
        throw std::invalid_argument("HACApK matrix is not assembled");
    std::vector<std::size_t> permutation(matrix.nd);
    for(int i=0;i<matrix.nd;++i){
        const auto original=control.lod[i+1];
        if(original<1||original>matrix.nd)throw std::invalid_argument("Invalid native HACApK permutation");
        permutation[i]=std::size_t(original-1);
    }
    std::vector<Leaf> leaves;
    for(int i=1;i<=matrix.nlf;++i){
        const auto* native=matrix.st_lf[i];if(!native)continue;
        if(native->nstrtl>native->nstrtt)continue; // CPU symmetric wrapper ignores these too.
        if(native->ndl<=0||native->ndt<=0||native->nstrtl<1||native->nstrtt<1||
           native->nstrtl>matrix.nd||native->nstrtt>matrix.nd||
           native->ndl>matrix.nd-native->nstrtl+1||native->ndt>matrix.nd-native->nstrtt+1)
            throw std::invalid_argument("Invalid native leaf range");
        if(native->ltmtx!=1&&native->ltmtx!=2)throw std::invalid_argument("Unsupported native leaf kind");
        Leaf leaf;leaf.row=native->nstrtl-1;leaf.column=native->nstrtt-1;leaf.rows=native->ndl;leaf.columns=native->ndt;
        if(native->ltmtx==1){
            if(native->kt<=0||!native->a1||!native->a2)throw std::invalid_argument("Empty native low-rank upper leaf");
            leaf.rank=native->kt;
            leaf.a1.assign(native->a1,native->a1+checkedProduct(leaf.columns,leaf.rank));
            leaf.a2.assign(native->a2,native->a2+checkedProduct(leaf.rows,leaf.rank));
        }else{
            if(!native->a1)throw std::invalid_argument("Empty native dense upper leaf");
            leaf.rank=0;leaf.a1.assign(native->a1,native->a1+checkedProduct(leaf.rows,leaf.columns));
        }
        leaves.push_back(std::move(leaf));
    }
    return pack(matrix.nd,leaves,permutation,physicalScale);
}
} // namespace radia_hmatrix_cuda
