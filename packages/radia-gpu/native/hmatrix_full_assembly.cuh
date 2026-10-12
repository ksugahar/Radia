// Existing HACApK hierarchy and compression; CUDA owns analytic entries.
#pragma once
#include "hmatrix_hacapk_pack.cuh"
#include "hmatrix_aca_callback.cuh"
#include <cstdlib>
extern "C" {
#include "cHACApK_pca_cluster.h"
int HACApK_build_hmatrix_wrapper(void*,void*,double*,int,int,int,double,int,double,int,int);
void HACApK_free_hmatrix_wrapper(void*,void*);
void HACApK_free_leafmtxp(void*);
}
namespace radia_hmatrix_cuda {
struct FullAssembly {
    Plan plan;
    std::size_t entryBatches=0,entryCount=0,oracleOwned=0,leafCount=0,lowrankCount=0;
};
inline FullAssembly assembleFull(const AnalyticBlock& block,double epsilon,int leafSize,double eta,int maxRank){
    const auto n=block.sources.size();
    // Bound worst-case host dense storage before existing C allocations.
    if(!n||n>4096||epsilon<=0||epsilon>=1||!std::isfinite(epsilon)||leafSize<1||
       leafSize>4096||eta<=0||!std::isfinite(eta)||maxRank<1||maxRank>4096)
        throw std::invalid_argument("Unsupported full analytic HACApK budget/options");
    static std::mutex buildMutex;std::lock_guard<std::mutex> lock(buildMutex);
    struct NativeOwner {
        st_cHACApK_leafmtxp matrix=static_cast<st_cHACApK_leafmtxp>(std::calloc(1,sizeof(st_cHACApK_leafmtxp_t)));
        st_cHACApK_lcontrol control=static_cast<st_cHACApK_lcontrol>(std::calloc(1,sizeof(st_cHACApK_lcontrol_t)));
        ~NativeOwner(){cHACApK_set_point_radius(nullptr);cHACApK_set_sym_fill(0);cHACApK_set_cluster_strategy(CHACAPK_CLUSTER_BBOX);
            HACApK_free_hmatrix_wrapper(matrix,control);HACApK_free_leafmtxp(matrix);std::free(control);}
    } native;
    if(!native.matrix||!native.control)throw std::bad_alloc();
    std::vector<double> centers(3*n),radii(n),scale(n);
    for(std::size_t i=0;i<n;++i){const auto& source=block.sources[i];const int vertices=source.volume?4:3;
        for(int v=0;v<vertices;++v)for(int k=0;k<3;++k)centers[3*i+k]+=source.vertices[3*v+k]/vertices;
        for(int v=0;v<vertices;++v){double d=0.;for(int k=0;k<3;++k){double x=source.vertices[3*v+k]-centers[3*i+k];d+=x*x;}radii[i]=std::max(radii[i],std::sqrt(d));}
    }
    AnalyticEntryOracle oracle(block,256);
    for(std::size_t start=0;start<n;start+=256){std::vector<AnalyticEntryIndex> indices;
        for(auto i=start;i<std::min(n,start+256);++i)indices.push_back({i,i});
        auto values=oracle.evaluate(indices,nullptr);
        for(std::size_t i=0;i<values.size();++i){if(values[i]<=0||!std::isfinite(values[i]))throw std::runtime_error("Invalid charge normalization diagonal");scale[start+i]=std::sqrt(values[i]);}
    }
    EntryCallbackScope callback(oracle,nullptr,&scale);
    cHACApK_set_point_radius(radii.data());cHACApK_set_sym_fill(1);cHACApK_set_cluster_strategy(CHACAPK_CLUSTER_BBOX);
    const int status=HACApK_build_hmatrix_wrapper(native.matrix,native.control,centers.data(),int(n),1,3,epsilon,leafSize,eta,maxRank,0);
    callback.rethrow();if(status)throw std::runtime_error("Native HACApK assembly failed");
    if(int(native.control->param[60])!=1)throw std::runtime_error("Unexpected native ACA dispatch");
    FullAssembly result;result.plan=fromHacapk(*native.matrix,*native.control,scale);
    result.entryBatches=oracle.batches;result.entryCount=oracle.entries;result.oracleOwned=oracle.ownedBytes;
    result.leafCount=native.matrix->nlf;result.lowrankCount=native.matrix->nlfkt;return result;
}
}
