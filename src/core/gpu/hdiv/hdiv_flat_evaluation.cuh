// Device construction of the affine-TET HDiv material evaluation map E and E^T
// and of the reference mass M = E^T W E with its diagonal and Jacobi bound.
// NGSolve tabulates reference shapes per vertex-rank class on the host; only
// those tables, the element vertices, DOF numbers and point numbers cross.
#pragma once
#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_run_length_encode.cuh>
#include <cub/device/device_scan.cuh>
#include <cub/device/device_reduce.cuh>
#include <sstream>
#include "hmatrix_resident.cuh"
namespace radia_hdiv_bh {
// Joint (row,column) mass key row*nd+column: 64-bit, exact for any int32 DOF count.
using MassKey=unsigned long long;
__global__ void flatEvaluationEntries(std::size_t ne,std::size_t local,std::size_t rule,std::size_t nq,std::size_t classes,
    const double* vertices,const int* elementClass,const int* elementDofs,const int* elementPoints,const double* shapes,
    const double* ruleRows,std::size_t* indices,double* values,double* weights,int* invalid){
    const auto t=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(t>=ne*rule)return;
    const auto e=t/rule,k=t%rule;const double* V=vertices+12*e;
    double J[3][3];for(int c=0;c<3;++c)for(int d=0;d<3;++d)J[c][d]=V[3*d+c]-V[9+c];
    const double det=J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
    const int cls=elementClass[e],q=elementPoints[e*rule+k];
    if(!isfinite(det)||det==0.||cls<0||std::size_t(cls)>=classes||q<0||std::size_t(q)>=nq){atomicExch(invalid,1);return;}
    weights[q]=ruleRows[4*k+3]*fabs(det);
    // Rows keep ascending DOF columns; local sizes are bounded by 64 on the host.
    int order[64];const int* dofs=elementDofs+e*local;
    for(std::size_t i=0;i<local;++i){std::size_t j=i;while(j>0&&dofs[order[j-1]]>dofs[i]){order[j]=order[j-1];--j;}order[j]=int(i);}
    const double* S=shapes+(std::size_t(cls)*rule+k)*local*3;
    for(int c=0;c<3;++c){const auto base=(std::size_t(c)*nq+std::size_t(q))*local;
        for(std::size_t r=0;r<local;++r){const int i=order[r];
            const double v=(J[c][0]*S[3*i]+J[c][1]*S[3*i+1]+J[c][2]*S[3*i+2])/det;
            if(!isfinite(v))atomicExch(invalid,1);indices[base+r]=std::size_t(dofs[i]);values[base+r]=v;}}
}
__global__ void uniformOffsets(std::size_t rows,std::size_t width,std::size_t* offsets){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<=rows)offsets[i]=i*width;
}
__global__ void columnKeys(std::size_t nnz,const std::size_t* indices,unsigned* keys,unsigned* entries){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<nnz){keys[i]=unsigned(indices[i]);entries[i]=unsigned(i);}
}
// Stable radix order keeps rows ascending within each column (deterministic E^T).
__global__ void transposeEntries(std::size_t nnz,std::size_t width,const unsigned* entries,const double* values,
    std::size_t* indices,double* transposed){
    const auto j=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(j<nnz){indices[j]=entries[j]/width;transposed[j]=values[entries[j]];}
}
__global__ void columnOffsets(std::size_t columns,std::size_t nnz,const unsigned* keys,std::size_t* offsets){
    const auto col=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(col>columns)return;
    std::size_t lo=0,hi=nnz;while(lo<hi){const auto mid=lo+(hi-lo)/2;if(keys[mid]<col)lo=mid+1;else hi=mid;}offsets[col]=lo;
}
// One thread per element-local pair: M_e(r,s)=sum_k sum_c w_q (E_cqr E_cqs).
// The product E_cqr E_cqs commutes exactly, so M is symmetric by construction.
__global__ void flatMassEntries(std::size_t ne,std::size_t local,std::size_t rule,std::size_t nq,std::size_t nd,
    const int* elementPoints,const std::size_t* indices,const double* values,const double* weights,MassKey* keys,double* entries){
    const auto t=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(t>=ne*local*local)return;
    const auto e=t/(local*local),r=(t/local)%local,s=t%local;radia_hmatrix_cuda::CompensatedDeviceSum sum;
    for(std::size_t k=0;k<rule;++k){const auto q=std::size_t(elementPoints[e*rule+k]);const double w=weights[q];
        for(int c=0;c<3;++c){const auto base=(std::size_t(c)*nq+q)*local;sum.addProduct(w,values[base+r]*values[base+s]);}}
    const auto base=std::size_t(elementPoints[e*rule])*local;keys[t]=MassKey(indices[base+r])*nd+indices[base+s];entries[t]=sum.get();
}
// Sorted runs of equal (row,column) keys are summed in stable element order (deterministic).
__global__ void flatMassReduce(std::size_t nnz,std::size_t nd,const MassKey* keys,const int* starts,const int* counts,
    const double* entries,std::size_t* indices,double* values){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=nnz)return;radia_hmatrix_cuda::CompensatedDeviceSum sum;
    for(int k=0;k<counts[i];++k)sum.add(entries[starts[i]+k]);indices[i]=keys[i]%nd;values[i]=sum.get();
}
__global__ void flatMassRowOffsets(std::size_t nd,std::size_t nnz,const MassKey* keys,std::size_t* offsets){
    const auto row=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(row>nd)return;const auto key=MassKey(row)*nd;
    std::size_t lo=0,hi=nnz;while(lo<hi){const auto mid=lo+(hi-lo)/2;if(keys[mid]<key)lo=mid+1;else hi=mid;}offsets[row]=lo;
}
__device__ inline bool flatMassFind(const std::size_t* offsets,const std::size_t* indices,std::size_t row,std::size_t column,std::size_t& at){
    std::size_t lo=offsets[row],hi=offsets[row+1];while(lo<hi){const auto mid=lo+(hi-lo)/2;if(indices[mid]<column)lo=mid+1;else hi=mid;}
    at=lo;return lo<offsets[row+1]&&indices[lo]==column;
}
// Same symmetry tolerance and positive diagonal as the host mass check; row bound sum_k|M_rk|/M_rr.
__global__ void flatMassDiagnostics(std::size_t nd,const std::size_t* offsets,const std::size_t* indices,const double* values,
    double* diagonal,double* rowBound,int* invalid){
    const auto row=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(row>=nd)return;std::size_t at=0;
    if(!flatMassFind(offsets,indices,row,row,at)||!isfinite(values[at])||!(values[at]>0)){atomicExch(invalid,1);diagonal[row]=1.;rowBound[row]=0.;return;}
    const double d=values[at];double sum=0.;
    for(auto k=offsets[row];k<offsets[row+1];++k){const double v=values[k];std::size_t t=0;
        if(!isfinite(v)||!flatMassFind(offsets,indices,indices[k],row,t)||fabs(v-values[t])>1.e-13*fmax(1.,fmax(fabs(v),fabs(values[t]))))atomicCAS(invalid,0,2);
        sum+=fabs(v/d);}
    diagonal[row]=d;rowBound[row]=sum;
}
struct FlatEvaluationTransfer {const char* payload;const char* direction;std::size_t bytes;};
// Builds E (3*nq x nd, component-major rows), E^T, point measures and the reference
// mass CSR M with its diagonal into caller-owned device storage. Every transient
// buffer lives until build() returns, so transientBytes is their simultaneous total
// (driver-internal allocations excluded).
struct FlatEvaluationBuilder {
    template<class T>using Buffer=radia_hmatrix_cuda::Buffer<T>;
    radia_hmatrix_cuda::HostVector<FlatEvaluationTransfer> transfers;std::size_t launches=0,nonzeros=0,massNonzeros=0,transientBytes=0;double massBound=0.;
    void build(const FlatEvaluationInput& in,radia_hmatrix_cuda::DeviceMap& E,radia_hmatrix_cuda::DeviceMap& ET,Buffer<double>& weights,
               radia_hmatrix_cuda::DeviceMap& M,Buffer<double>& massDiagonal){
        using namespace radia_hmatrix_cuda;
        const auto ne=in.elements,local=in.local,rule=in.rulePoints,nq=in.points,nd=in.dofs;
        nonzeros=checkedProduct(checkedProduct(3*nq,local),1);
        Buffer<double> vertices,shapes,ruleRows;Buffer<int> classes,dofs,points,invalid;
        auto upload=[&](auto& buffer,const auto& values,const char* payload){buffer.upload(values);const auto bytes=checkedProduct(values.size(),sizeof(values[0]));
            transfers.push_back({payload,"H2D",bytes});transientBytes+=bytes;};
        auto scratch=[&](auto& buffer,std::size_t count){buffer.allocate(count);transientBytes+=checkedProduct(count,sizeof(*buffer.data));};
        upload(vertices,in.vertices,"element_vertices_m");upload(shapes,in.shapes,"reference_shapes");upload(ruleRows,in.rule,"material_rule");
        upload(classes,in.elementClass,"element_class");upload(dofs,in.elementDofs,"element_dofs");upload(points,in.elementPoints,"element_points");
        upload(invalid,HostVector<int>{0},"evaluation_invalid_flag");
        E.rows=3*nq;E.offsets.allocate(3*nq+1);E.indices.allocate(nonzeros);E.values.allocate(nonzeros);weights.allocate(nq);
        flatEvaluationEntries<<<unsigned((ne*rule+127)/128),128>>>(ne,local,rule,nq,in.classes,vertices.data,classes.data,dofs.data,points.data,
            shapes.data,ruleRows.data,E.indices.data,E.values.data,weights.data,invalid.data);check(cudaGetLastError());++launches;
        uniformOffsets<<<unsigned((3*nq+256)/256),256>>>(3*nq,local,E.offsets.data);check(cudaGetLastError());++launches;
        int failed=0;check(cudaMemcpy(&failed,invalid.data,sizeof(int),cudaMemcpyDeviceToHost));
        transfers.push_back({"evaluation_flag","D2H",sizeof(int)});
        if(failed)throw std::runtime_error("Device HDiv evaluation map has a degenerate element, invalid index or nonfinite shape");
        Buffer<unsigned> keys,sortedKeys,entries,sortedEntries;scratch(keys,nonzeros);scratch(sortedKeys,nonzeros);scratch(entries,nonzeros);scratch(sortedEntries,nonzeros);
        columnKeys<<<unsigned((nonzeros+255)/256),256>>>(nonzeros,E.indices.data,keys.data,entries.data);check(cudaGetLastError());++launches;
        int bits=1;while((std::size_t(1)<<bits)<nd)++bits;std::size_t bytes=0;
        check(cub::DeviceRadixSort::SortPairs(nullptr,bytes,keys.data,sortedKeys.data,entries.data,sortedEntries.data,int(nonzeros),0,bits));
        Buffer<unsigned char> temporary;scratch(temporary,bytes);
        check(cub::DeviceRadixSort::SortPairs(temporary.data,bytes,keys.data,sortedKeys.data,entries.data,sortedEntries.data,int(nonzeros),0,bits));++launches;
        ET.rows=nd;ET.offsets.allocate(nd+1);ET.indices.allocate(nonzeros);ET.values.allocate(nonzeros);
        transposeEntries<<<unsigned((nonzeros+255)/256),256>>>(nonzeros,local,sortedEntries.data,E.values.data,ET.indices.data,ET.values.data);check(cudaGetLastError());++launches;
        columnOffsets<<<unsigned((nd+256)/256),256>>>(nd,nonzeros,sortedKeys.data,ET.offsets.data);check(cudaGetLastError());++launches;
        buildMass(in,E,weights,points,M,massDiagonal,scratch);
        check(cudaDeviceSynchronize());
    }
private:
    template<class Scratch>void buildMass(const FlatEvaluationInput& in,radia_hmatrix_cuda::DeviceMap& E,Buffer<double>& weights,Buffer<int>& points,
                                          radia_hmatrix_cuda::DeviceMap& M,Buffer<double>& massDiagonal,Scratch& scratch){
        using namespace radia_hmatrix_cuda;
        const auto ne=in.elements,local=in.local,rule=in.rulePoints,nq=in.points,nd=in.dofs;
        const auto triples=checkedProduct(ne,checkedProduct(local,local));
        Buffer<MassKey> keys,sortedKeys,uniqueKeys;Buffer<double> entries,sortedEntries,rowBound,bound;Buffer<int> counts,starts,runs,invalid;
        scratch(keys,triples);scratch(sortedKeys,triples);scratch(uniqueKeys,triples);scratch(entries,triples);scratch(sortedEntries,triples);
        scratch(counts,triples);scratch(starts,triples);scratch(runs,1);scratch(rowBound,nd);scratch(bound,1);scratch(invalid,1);
        flatMassEntries<<<unsigned((triples+255)/256),256>>>(ne,local,rule,nq,nd,points.data,E.indices.data,E.values.data,weights.data,keys.data,entries.data);
        check(cudaGetLastError());++launches;
        int bits=1;while((std::size_t(1)<<bits)<nd*nd)++bits;std::size_t bytes=0;
        check(cub::DeviceRadixSort::SortPairs(nullptr,bytes,keys.data,sortedKeys.data,entries.data,sortedEntries.data,int(triples),0,bits));
        Buffer<unsigned char> sortScratch;scratch(sortScratch,bytes);
        check(cub::DeviceRadixSort::SortPairs(sortScratch.data,bytes,keys.data,sortedKeys.data,entries.data,sortedEntries.data,int(triples),0,bits));++launches;
        bytes=0;check(cub::DeviceRunLengthEncode::Encode(nullptr,bytes,sortedKeys.data,uniqueKeys.data,counts.data,runs.data,int(triples)));
        Buffer<unsigned char> encodeScratch;scratch(encodeScratch,bytes);
        check(cub::DeviceRunLengthEncode::Encode(encodeScratch.data,bytes,sortedKeys.data,uniqueKeys.data,counts.data,runs.data,int(triples)));++launches;
        int nnz=0;check(cudaMemcpy(&nnz,runs.data,sizeof(int),cudaMemcpyDeviceToHost));transfers.push_back({"mass_nonzeros","D2H",sizeof(int)});
        if(nnz<int(nd)||std::size_t(nnz)>triples)throw std::runtime_error("Device HDiv reference mass has an invalid sparsity pattern");
        massNonzeros=std::size_t(nnz);
        bytes=0;check(cub::DeviceScan::ExclusiveSum(nullptr,bytes,counts.data,starts.data,nnz));
        Buffer<unsigned char> scanScratch;scratch(scanScratch,bytes);
        check(cub::DeviceScan::ExclusiveSum(scanScratch.data,bytes,counts.data,starts.data,nnz));++launches;
        M.rows=nd;M.offsets.allocate(nd+1);M.indices.allocate(massNonzeros);M.values.allocate(massNonzeros);massDiagonal.allocate(nd);
        flatMassReduce<<<unsigned((massNonzeros+255)/256),256>>>(massNonzeros,nd,uniqueKeys.data,starts.data,counts.data,sortedEntries.data,M.indices.data,M.values.data);
        check(cudaGetLastError());++launches;
        flatMassRowOffsets<<<unsigned((nd+256)/256),256>>>(nd,massNonzeros,uniqueKeys.data,M.offsets.data);check(cudaGetLastError());++launches;
        check(cudaMemset(invalid.data,0,sizeof(int)));
        flatMassDiagnostics<<<unsigned((nd+255)/256),256>>>(nd,M.offsets.data,M.indices.data,M.values.data,massDiagonal.data,rowBound.data,invalid.data);
        check(cudaGetLastError());++launches;
        bytes=0;check(cub::DeviceReduce::Max(nullptr,bytes,rowBound.data,bound.data,int(nd)));
        Buffer<unsigned char> reduceScratch;scratch(reduceScratch,bytes);
        check(cub::DeviceReduce::Max(reduceScratch.data,bytes,rowBound.data,bound.data,int(nd)));++launches;
        int status=0;check(cudaMemcpy(&status,invalid.data,sizeof(int),cudaMemcpyDeviceToHost));transfers.push_back({"mass_flag","D2H",sizeof(int)});
        check(cudaMemcpy(&massBound,bound.data,sizeof(double),cudaMemcpyDeviceToHost));transfers.push_back({"mass_spectral_bound","D2H",sizeof(double)});
        if(status==1)throw std::runtime_error("Device HDiv reference mass has a missing or nonpositive diagonal");
        if(status)throw std::runtime_error("Device HDiv reference mass is not symmetric");
        if(!std::isfinite(massBound)||massBound<1.)throw std::runtime_error("Invalid device mass spectral bound");
    }
};
// Host bounds before any upload: every index is checked against its range.
// No element, DOF or point cap: int32 tables and the int item counts of the
// CUB sort/encode/scan calls are the structural bound, and the live device
// budget (requireFlatEvaluationDevice) decides what fits.
inline void validateFlatEvaluation(const FlatEvaluationInput& in){
    using radia_hmatrix_cuda::checkedProduct;
    const auto ne=in.elements,local=in.local,rule=in.rulePoints;constexpr std::size_t count=std::numeric_limits<int>::max();
    if(!ne||!local||local>64||!rule||rule>4096||!in.classes||in.classes>24||!in.dofs||in.points!=checkedProduct(ne,rule))
        throw std::invalid_argument("Flat evaluation tables outside the affine-TET scope");
    if(ne>count||in.dofs>count||in.points>count)throw std::invalid_argument("Flat evaluation tables exceed int32 element, DOF or point numbering");
    if(checkedProduct(checkedProduct(3,in.points),local)>count)throw std::invalid_argument("Flat evaluation entries exceed the checked int CUB item count");
    if(checkedProduct(ne,checkedProduct(local,local))>count)throw std::invalid_argument("Flat mass triples exceed the checked int CUB item count");
    if(in.vertices.size()!=12*ne||in.elementClass.size()!=ne||in.elementDofs.size()!=ne*local||in.elementPoints.size()!=ne*rule||
       in.rule.size()!=4*rule||in.shapes.size()!=in.classes*rule*local*3)throw std::invalid_argument("Flat evaluation table shape mismatch");
    for(double v:in.vertices)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite evaluation vertex");
    for(double v:in.shapes)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite reference shape");
    for(std::size_t k=0;k<rule;++k)if(!std::isfinite(in.rule[4*k])||!std::isfinite(in.rule[4*k+1])||!std::isfinite(in.rule[4*k+2])||!(in.rule[4*k+3]>0)||!std::isfinite(in.rule[4*k+3]))
        throw std::invalid_argument("Invalid material reference rule");
    for(int c:in.elementClass)if(c<0||std::size_t(c)>=in.classes)throw std::invalid_argument("Evaluation class index out of range");
    for(int d:in.elementDofs)if(d<0||std::size_t(d)>=in.dofs)throw std::invalid_argument("Evaluation DOF index out of range");
    std::vector<unsigned char> seen(in.points,0);
    for(int q:in.elementPoints){if(q<0||std::size_t(q)>=in.points||seen[q])throw std::invalid_argument("Evaluation points must cover each material point once");seen[q]=1;}
    for(std::size_t e=0;e<ne;++e)for(std::size_t i=0;i<local;++i)for(std::size_t j=0;j<i;++j)
        if(in.elementDofs[e*local+i]==in.elementDofs[e*local+j])throw std::invalid_argument("Duplicate evaluation local DOF");
}
// CUB scratch of FlatEvaluationBuilder, queried with null storage (no allocation).
// All of it is live together until build() returns; the scan is sized for its
// largest possible run count.
inline std::size_t flatEvaluationSortScratch(const FlatEvaluationInput& in){
    using namespace radia_hmatrix_cuda;
    const auto nnz=checkedProduct(checkedProduct(3,in.points),in.local),triples=checkedProduct(in.elements,checkedProduct(in.local,in.local)),nd=in.dofs;
    int bits=1;while((std::size_t(1)<<bits)<nd)++bits;int massBits=1;while((std::size_t(1)<<massBits)<nd*nd)++massBits;
    std::size_t total=0,bytes=0;auto add=[&]{if(bytes>std::numeric_limits<std::size_t>::max()-total)throw std::overflow_error("Flat material budget overflow");total+=bytes;bytes=0;};
    check(cub::DeviceRadixSort::SortPairs(nullptr,bytes,(const unsigned*)nullptr,(unsigned*)nullptr,(const unsigned*)nullptr,(unsigned*)nullptr,int(nnz),0,bits));add();
    check(cub::DeviceRadixSort::SortPairs(nullptr,bytes,(const MassKey*)nullptr,(MassKey*)nullptr,(const double*)nullptr,(double*)nullptr,int(triples),0,massBits));add();
    check(cub::DeviceRunLengthEncode::Encode(nullptr,bytes,(const MassKey*)nullptr,(MassKey*)nullptr,(int*)nullptr,(int*)nullptr,int(triples)));add();
    check(cub::DeviceScan::ExclusiveSum(nullptr,bytes,(const int*)nullptr,(int*)nullptr,int(triples)));add();
    check(cub::DeviceReduce::Max(nullptr,bytes,(const double*)nullptr,(double*)nullptr,int(nd)));add();
    return total;
}
// Live preflight before any material allocation: the dimension budget plus CUB
// scratch must fit the current free memory above the shared reserve. The
// failure names the limiting term; nothing spills to the host.
inline std::size_t requireFlatEvaluationDevice(const FlatEvaluationInput& in,std::size_t lawDoubles,FlatEvaluationBudget& budget){
    using namespace radia_hmatrix_cuda;
    budget=flatEvaluationBudget(in,lawDoubles);const auto scratch=flatEvaluationSortScratch(in),need=budget.total();
    if(scratch>std::numeric_limits<std::size_t>::max()-need)throw std::overflow_error("Flat material budget overflow");
    std::size_t free=0,total=0;check(cudaMemGetInfo(&free,&total));const auto reserve=std::max(std::size_t(2)<<30,total/5);
    const auto available=free>reserve?free-reserve:0;
    if(need+scratch>available){const auto at=budget.limiting();std::ostringstream message;
        message<<"Flat material device budget exceeded: "<<need+scratch<<" bytes needed, "<<available<<" free above the reserve; limiting term "
               <<FlatEvaluationBudget::names[at]<<" = "<<budget.bytes[at]<<" bytes";
        throw std::runtime_error(message.str());}
    return need+scratch;
}
}
