// Existing HACApK hierarchy and compression; CUDA owns analytic entries.
#pragma once
#include "hmatrix_hacapk_pack.cuh"
#include "hmatrix_aca_callback.cuh"
#include "hmatrix_curved_entry_oracle.cuh"
#include <cstdlib>
#include "hmatrix_c_memory.h"
extern "C" {
#include "cHACApK_pca_cluster.h"
int HACApK_build_hmatrix_wrapper(void*,void*,double*,int,int,int,double,int,double,int,int);
void HACApK_free_hmatrix_wrapper(void*,void*);
void HACApK_free_leafmtxp(void*);
}
namespace radia_hmatrix_cuda {
struct FullAssembly {
    Plan plan;
    LeafCoefficientBounds frameBounds;
    std::shared_ptr<CurvedEntryOracle> retainedOracle;
    std::size_t entryBatches=0,entryCount=0,oracleOwned=0,leafCount=0,lowrankCount=0;
    std::size_t cacheTileLoads=0,cacheTileHits=0,cacheFullTileBatches=0,cacheHostScratchBytes=0;
};
inline std::mutex& fullAssemblyMutex(){static std::mutex mutex;return mutex;}
template<class Oracle>
inline FullAssembly assembleFullEntries(Oracle& oracle,HostVector<double> centers,const HostVector<double>& radii,double epsilon,int leafSize,double eta,int maxRank,bool deferDense=false,std::size_t frameBudget=std::size_t(2)<<30,std::size_t workingBudget=std::size_t(2)<<30){
    const auto n=radii.size();
    // C dense/factor products use signed int. The bounded extension has nffc=1;
    // n^2 <= 21,233,664 and 2*n^2 FP64 coefficients <= 339,738,624 bytes.
    // Restrict its tree/rank settings rather than exposing arbitrary leaf counts.
    if(!n||n>4608||centers.size()!=3*n||epsilon<=0||epsilon>=1||!std::isfinite(epsilon)||leafSize<1||
       leafSize>4096||eta<=0||!std::isfinite(eta)||maxRank<1||maxRank>4096)
        throw std::invalid_argument("Unsupported full analytic HACApK budget/options");
    const auto square=checkedProduct(n,n);
    if(square>std::size_t(std::numeric_limits<int>::max())/2||
       checkedProduct(square,2*sizeof(double))>(std::size_t(384)<<20)||
       (n>4096&&(leafSize<64||maxRank>256)))
        throw std::invalid_argument("Full HACApK C allocation/index budget exceeded");
    for(const auto x:centers)if(!std::isfinite(x))throw std::invalid_argument("Nonfinite cluster center");
    for(const auto r:radii)if(!std::isfinite(r)||r<0)throw std::invalid_argument("Invalid cluster support radius");
    std::lock_guard<std::mutex> lock(fullAssemblyMutex());
    struct NativeOwner {
        st_cHACApK_leafmtxp matrix=static_cast<st_cHACApK_leafmtxp>(std::calloc(1,sizeof(st_cHACApK_leafmtxp_t)));
        st_cHACApK_lcontrol control=static_cast<st_cHACApK_lcontrol>(std::calloc(1,sizeof(st_cHACApK_lcontrol_t)));
        ~NativeOwner(){cHACApK_set_point_radius(nullptr);cHACApK_set_sym_fill(0);cHACApK_set_cluster_strategy(CHACAPK_CLUSTER_BBOX);
            HACApK_free_hmatrix_wrapper(matrix,control);HACApK_free_leafmtxp(matrix);std::free(control);}
    } native;
    if(!native.matrix||!native.control)throw std::bad_alloc();
    HostVector<double> scale(n);
    for(std::size_t start=0;start<n;start+=EntryBatchCapacity){HostVector<AnalyticEntryIndex> indices;
        for(auto i=start;i<std::min(n,start+EntryBatchCapacity);++i)indices.push_back({i,i});
        auto values=oracle.evaluate(indices,nullptr);
        for(std::size_t i=0;i<values.size();++i){if(values[i]<=0||!std::isfinite(values[i]))throw std::runtime_error("Invalid charge normalization diagonal");scale[start+i]=std::sqrt(values[i]);}
    }
    EntryCallbackScope callback(oracle,nullptr,&scale,deferDense,frameBudget,workingBudget);
    cHACApK_set_point_radius(radii.data());cHACApK_set_sym_fill(1);cHACApK_set_cluster_strategy(CHACAPK_CLUSTER_BBOX);
    const int status=rad_hmatrix_guarded_build(native.matrix,native.control,centers.data(),int(n),1,3,epsilon,leafSize,eta,maxRank,0);
    callback.rethrow();if(status)throw std::runtime_error("Native HACApK assembly failed");
    if(int(native.control->param[60])!=1)throw std::runtime_error("Unexpected native ACA dispatch");
    FullAssembly result;
#ifdef RAD_HACAPK_FRAME_ADMISSION
    result.frameBounds=callback.admittedFrame();
#endif
    result.plan=fromHacapk(*native.matrix,*native.control,scale,deferDense);
    if(deferDense&&result.plan.deferredDense.empty())throw std::runtime_error("Declared device dense route did not defer any leaves");
    result.entryBatches=oracle.batches;result.entryCount=oracle.entries;result.oracleOwned=oracle.ownedBytes;
    result.cacheTileLoads=oracle.cacheTileLoads;result.cacheTileHits=oracle.cacheTileHits;result.cacheFullTileBatches=oracle.cacheFullTileBatches;result.cacheHostScratchBytes=oracle.cacheHostScratchBytes;
    for(int i=1;i<=native.matrix->nlf;++i){const auto* leaf=native.matrix->st_lf[i];
        if(leaf&&leaf->nstrtl<=leaf->nstrtt){++result.leafCount;if(leaf->ltmtx==1&&leaf->kt>0)++result.lowrankCount;}}
    return result;
}
inline FullAssembly assembleFull(const AnalyticBlock& block,double epsilon,int leafSize,double eta,int maxRank,std::size_t frameBudget=std::size_t(2)<<30,std::size_t workingBudget=std::size_t(2)<<30){
    const auto n=block.sources.size();if(!n||n>4096)throw std::invalid_argument("Unsupported analytic HACApK budget");
    HostVector<double> centers(3*n),radii(n);
    for(std::size_t i=0;i<n;++i){const auto& source=block.sources[i];const int vertices=source.volume?4:3;
        for(int v=0;v<vertices;++v)for(int k=0;k<3;++k)centers[3*i+k]+=source.vertices[3*v+k]/vertices;
        for(int v=0;v<vertices;++v){double d=0.;for(int k=0;k<3;++k){double x=source.vertices[3*v+k]-centers[3*i+k];d+=x*x;}radii[i]=std::max(radii[i],std::sqrt(d));}
    }
    AnalyticEntryOracle oracle(block,EntryBatchCapacity);return assembleFullEntries(oracle,std::move(centers),radii,epsilon,leafSize,eta,maxRank,false,frameBudget,workingBudget);
}
// P2 Lagrange geometry is inside the convex hull of its quadratic Bernstein
// controls. Edge control = 2*midpoint - (vertex_i+vertex_j)/2. Vertex-only boxes
// are unsafe for curved elements. Node order matches the shared P2 kernels.
inline void curvedSupportBox(const CurvedCharge& charge,double* center,double& radius){
    const int vertices=charge.kind?3:4,edges=charge.kind?3:6;
    const int pairs[6][2]={{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
    double controls[10][3]={},magnitude=1.;
    for(int i=0;i<vertices;++i)for(int c=0;c<3;++c)controls[i][c]=charge.nodes[3*i+c];
    for(int e=0;e<edges;++e)for(int c=0;c<3;++c)
        controls[vertices+e][c]=2.*charge.nodes[3*(vertices+e)+c]-.5*charge.nodes[3*pairs[e][0]+c]-.5*charge.nodes[3*pairs[e][1]+c];
    for(int i=0;i<vertices+edges;++i)for(int c=0;c<3;++c){
        if(!std::isfinite(controls[i][c]))throw std::invalid_argument("Curved support control overflow");
        magnitude=std::max(magnitude,std::abs(controls[i][c]));magnitude=std::max(magnitude,std::abs(charge.nodes[3*i+c]));}
    const double padding=64.*std::numeric_limits<double>::epsilon()*magnitude;
    radius=0.;for(int c=0;c<3;++c){double lo=controls[0][c],hi=lo;
        for(int i=1;i<vertices+edges;++i){lo=std::min(lo,controls[i][c]);hi=std::max(hi,controls[i][c]);}
        center[c]=.5*lo+.5*hi;
        radius=std::max(radius,std::max(std::abs(lo-center[c]),std::abs(hi-center[c]))+padding);}
    radius=std::nextafter(radius,std::numeric_limits<double>::infinity());
    if(!std::isfinite(radius))throw std::invalid_argument("Curved support radius overflow");
}
inline FullAssembly assembleFullCurved(const CurvedBlock& block,double epsilon,int leafSize,double eta,int maxRank,bool cache=false,std::size_t cacheBudget=std::size_t(256)<<20,bool deviceDense=false){
    const auto n=block.charges.size();if(!n||n>4608||
        (n>4096&&(leafSize<64||maxRank>256)))throw std::invalid_argument("Unsupported curved HACApK budget");
    HostVector<double> centers(3*n),radii(n);
    for(std::size_t i=0;i<n;++i)curvedSupportBox(block.charges[i],centers.data()+3*i,radii[i]);
    auto oracle=std::allocate_shared<CurvedEntryOracle>(HostFaultAllocator<CurvedEntryOracle>("oracle.shared"),block,EntryBatchCapacity,cache,cacheBudget,deviceDense);
    auto result=assembleFullEntries(*oracle,std::move(centers),radii,epsilon,leafSize,eta,maxRank,deviceDense);
    if(deviceDense)result.retainedOracle=std::move(oracle);return result;
}
}
