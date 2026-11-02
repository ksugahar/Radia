// Matrix-free device action of the affine-TET HDiv material evaluation map E
// and E^T, and device construction of the reference mass M = E^T W E with its
// diagonal and Jacobi bound. NGSolve tabulates reference shapes per vertex-rank
// class on the host; only those tables, the element vertices, DOF numbers and
// point numbers cross, and they stay resident: no global E is stored.
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
// Affine TET map in NGSolve vertex order: J columns are v_d - v_3 (d=0..2).
__device__ inline double affineTet(const double* V,double J[3][3]){
    for(int c=0;c<3;++c)for(int d=0;d<3;++d)J[c][d]=V[3*d+c]-V[9+c];
    return J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
}
// The one formula of E: contravariant Piola value of reference shape S, component c.
// The setup check, E x, E^T y and the mass all call it.
__device__ inline double piolaValue(const double J[3][3],double det,const double* S,int c){
    return (J[c][0]*S[0]+J[c][1]*S[1]+J[c][2]*S[2])/det;
}
// Local slots in ascending DOF order (stable insertion sort; local <= 64 on the host).
__global__ void flatElementOrder(std::size_t ne,std::size_t local,const int* elementDofs,unsigned char* order){
    const auto e=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(e>=ne)return;
    const int* dofs=elementDofs+e*local;unsigned char* o=order+e*local;
    for(std::size_t i=0;i<local;++i){std::size_t j=i;while(j>0&&dofs[o[j-1]]>dofs[i]){o[j]=o[j-1];--j;}o[j]=(unsigned char)i;}
}
// Point measures and the validity of every E value.
__global__ void flatEvaluationEntries(std::size_t ne,std::size_t local,std::size_t rule,std::size_t nq,std::size_t classes,
    const double* vertices,const int* elementClass,const int* elementPoints,const double* shapes,
    const double* ruleRows,double* weights,int* invalid){
    const auto t=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(t>=ne*rule)return;
    const auto e=t/rule,k=t%rule;double J[3][3];const double det=affineTet(vertices+12*e,J);
    const int cls=elementClass[e],q=elementPoints[e*rule+k];
    if(!isfinite(det)||det==0.||cls<0||std::size_t(cls)>=classes||q<0||std::size_t(q)>=nq){atomicExch(invalid,1);return;}
    weights[q]=ruleRows[4*k+3]*fabs(det);
    const double* S=shapes+(std::size_t(cls)*rule+k)*local*3;
    for(int c=0;c<3;++c)for(std::size_t i=0;i<local;++i)if(!isfinite(piolaValue(J,det,S+3*i,c)))atomicExch(invalid,1);
}
// v = E x, one thread per material point, local slots in ascending DOF order.
__global__ void flatEvaluationApply(std::size_t ne,std::size_t local,std::size_t rule,std::size_t nq,
    const double* vertices,const int* elementClass,const int* elementDofs,const int* elementPoints,const unsigned char* order,
    const double* shapes,const double* x,double* v){
    const auto t=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(t>=ne*rule)return;
    const auto e=t/rule,k=t%rule;double J[3][3];const double det=affineTet(vertices+12*e,J);
    const auto q=std::size_t(elementPoints[t]);const int* dofs=elementDofs+e*local;const unsigned char* o=order+e*local;
    const double* S=shapes+(std::size_t(elementClass[e])*rule+k)*local*3;
    for(int c=0;c<3;++c){radia_hmatrix_cuda::CompensatedDeviceSum sum;
        for(std::size_t r=0;r<local;++r){const int i=o[r];sum.addProduct(piolaValue(J,det,S+3*i,c),x[dofs[i]]);}
        v[std::size_t(c)*nq+q]=sum.get();}
}
// g = E^T y, one thread per DOF over its incident (element, slot) entries in
// ascending element order: deterministic, with no floating-point atomics. For
// each component the terms follow ascending (element, point).
__global__ void flatEvaluationTransposeApply(std::size_t nd,std::size_t local,std::size_t rule,std::size_t nq,
    const double* vertices,const int* elementClass,const int* elementPoints,const double* shapes,
    const std::size_t* offsets,const unsigned* incidence,const double* y,double* g){
    const auto d=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(d>=nd)return;radia_hmatrix_cuda::CompensatedDeviceSum sum;
    for(int c=0;c<3;++c)for(auto at=offsets[d];at<offsets[d+1];++at){
        const auto e=std::size_t(incidence[at])/local,i=std::size_t(incidence[at])%local;double J[3][3];const double det=affineTet(vertices+12*e,J);
        const double* S=shapes+std::size_t(elementClass[e])*rule*local*3+3*i;
        for(std::size_t k=0;k<rule;++k)sum.addProduct(piolaValue(J,det,S+k*local*3,c),y[std::size_t(c)*nq+std::size_t(elementPoints[e*rule+k])]);}
    g[d]=sum.get();
}
__global__ void slotKeys(std::size_t count,const int* elementDofs,unsigned* keys,unsigned* entries){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<count){keys[i]=unsigned(elementDofs[i]);entries[i]=unsigned(i);}
}
__global__ void columnKeys(std::size_t nnz,const std::size_t* indices,unsigned* keys,unsigned* entries){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<nnz){keys[i]=unsigned(indices[i]);entries[i]=unsigned(i);}
}
__global__ void columnOffsets(std::size_t columns,std::size_t nnz,const unsigned* keys,std::size_t* offsets){
    const auto col=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(col>columns)return;
    std::size_t lo=0,hi=nnz;while(lo<hi){const auto mid=lo+(hi-lo)/2;if(keys[mid]<col)lo=mid+1;else hi=mid;}offsets[col]=lo;
}
// One thread per element-local pair: M_e(r,s)=sum_k sum_c w_q (E_cqr E_cqs), the
// E values from piolaValue. The product E_cqr E_cqs commutes exactly, so M is
// symmetric by construction.
__global__ void flatMassEntries(std::size_t ne,std::size_t local,std::size_t rule,std::size_t nd,
    const double* vertices,const int* elementClass,const int* elementDofs,const int* elementPoints,const unsigned char* order,
    const double* shapes,const double* weights,MassKey* keys,double* entries){
    const auto t=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(t>=ne*local*local)return;
    const auto e=t/(local*local),r=(t/local)%local,s=t%local;double J[3][3];const double det=affineTet(vertices+12*e,J);
    const int i=order[e*local+r],j=order[e*local+s];const int* dofs=elementDofs+e*local;
    const double* S=shapes+std::size_t(elementClass[e])*rule*local*3;radia_hmatrix_cuda::CompensatedDeviceSum sum;
    for(std::size_t k=0;k<rule;++k){const double w=weights[elementPoints[e*rule+k]];const double* Sk=S+k*local*3;
        for(int c=0;c<3;++c)sum.addProduct(w,piolaValue(J,det,Sk+3*i,c)*piolaValue(J,det,Sk+3*j,c));}
    keys[t]=MassKey(dofs[i])*nd+MassKey(dofs[j]);entries[t]=sum.get();
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
// Resident tables of the matrix-free E: geometry-derived element state and the
// DOF incidence (entries e*local+slot ascending per DOF). No E value is stored.
struct FlatEvaluationOperator {
    template<class T>using Buffer=radia_hmatrix_cuda::Buffer<T>;
    std::size_t ne=0,local=0,rule=0,nq=0,nd=0;
    Buffer<double> vertices,shapes;Buffer<int> classes,dofs,points;Buffer<unsigned char> order;
    Buffer<std::size_t> incidenceOffsets;Buffer<unsigned> incidence;
    std::size_t bytes()const{
        return (vertices.size+shapes.size)*sizeof(double)+(classes.size+dofs.size+points.size)*sizeof(int)+order.size+
               incidenceOffsets.size*sizeof(std::size_t)+incidence.size*sizeof(unsigned);
    }
    void apply(const double* x,double* v,cudaStream_t stream){
        flatEvaluationApply<<<unsigned((ne*rule+127)/128),128,0,stream>>>(ne,local,rule,nq,vertices.data,classes.data,dofs.data,points.data,order.data,shapes.data,x,v);
        radia_hmatrix_cuda::check(cudaGetLastError());
    }
    void applyTranspose(const double* y,double* g,cudaStream_t stream){
        flatEvaluationTransposeApply<<<unsigned((nd+127)/128),128,0,stream>>>(nd,local,rule,nq,vertices.data,classes.data,points.data,shapes.data,
            incidenceOffsets.data,incidence.data,y,g);
        radia_hmatrix_cuda::check(cudaGetLastError());
    }
};
struct FlatEvaluationTransfer {const char* payload;const char* direction;std::size_t bytes;};
// Builds the matrix-free E tables, point measures and the reference mass CSR M
// with its diagonal into caller-owned device storage. Every transient buffer lives until build()
// returns, so transientBytes is their simultaneous total (driver-internal
// allocations excluded); the resident tables are the operator's bytes().
struct FlatEvaluationBuilder {
    template<class T>using Buffer=radia_hmatrix_cuda::Buffer<T>;
    radia_hmatrix_cuda::HostVector<FlatEvaluationTransfer> transfers;std::size_t launches=0,massNonzeros=0,transientBytes=0;double massBound=0.;
    void build(const FlatEvaluationInput& in,FlatEvaluationOperator& op,Buffer<double>& weights,radia_hmatrix_cuda::DeviceMap& M,Buffer<double>& massDiagonal){
        using namespace radia_hmatrix_cuda;
        const auto ne=in.elements,local=in.local,rule=in.rulePoints,nq=in.points,nd=in.dofs,slots=checkedProduct(ne,local);
        op.ne=ne;op.local=local;op.rule=rule;op.nq=nq;op.nd=nd;
        Buffer<double> ruleRows;Buffer<int> invalid;
        auto upload=[&](auto& buffer,const auto& values,const char* payload,bool resident){buffer.upload(values);const auto bytes=checkedProduct(values.size(),sizeof(values[0]));
            transfers.push_back({payload,"H2D",bytes});if(!resident)transientBytes+=bytes;};
        auto scratch=[&](auto& buffer,std::size_t count){buffer.allocate(count);transientBytes+=checkedProduct(count,sizeof(*buffer.data));};
        upload(op.vertices,in.vertices,"element_vertices_m",true);upload(op.shapes,in.shapes,"reference_shapes",true);upload(ruleRows,in.rule,"material_rule",false);
        upload(op.classes,in.elementClass,"element_class",true);upload(op.dofs,in.elementDofs,"element_dofs",true);upload(op.points,in.elementPoints,"element_points",true);
        upload(invalid,HostVector<int>{0},"evaluation_invalid_flag",false);
        op.order.allocate(slots);weights.allocate(nq);
        flatElementOrder<<<unsigned((ne+127)/128),128>>>(ne,local,op.dofs.data,op.order.data);check(cudaGetLastError());++launches;
        flatEvaluationEntries<<<unsigned((ne*rule+127)/128),128>>>(ne,local,rule,nq,in.classes,op.vertices.data,op.classes.data,op.points.data,
            op.shapes.data,ruleRows.data,weights.data,invalid.data);check(cudaGetLastError());++launches;
        int failed=0;check(cudaMemcpy(&failed,invalid.data,sizeof(int),cudaMemcpyDeviceToHost));
        transfers.push_back({"evaluation_flag","D2H",sizeof(int)});
        if(failed)throw std::runtime_error("Device HDiv evaluation map has a degenerate element, invalid index or nonfinite shape");
        int bits=1;while((std::size_t(1)<<bits)<nd)++bits;
        {   // DOF incidence: stable sort of the slot DOF numbers, so entries ascend per DOF.
            Buffer<unsigned> keys,sortedKeys,entries;scratch(keys,slots);scratch(sortedKeys,slots);scratch(entries,slots);
            op.incidence.allocate(slots);op.incidenceOffsets.allocate(nd+1);
            slotKeys<<<unsigned((slots+255)/256),256>>>(slots,op.dofs.data,keys.data,entries.data);check(cudaGetLastError());++launches;
            std::size_t bytes=0;check(cub::DeviceRadixSort::SortPairs(nullptr,bytes,keys.data,sortedKeys.data,entries.data,op.incidence.data,int(slots),0,bits));
            Buffer<unsigned char> temporary;scratch(temporary,bytes);
            check(cub::DeviceRadixSort::SortPairs(temporary.data,bytes,keys.data,sortedKeys.data,entries.data,op.incidence.data,int(slots),0,bits));++launches;
            columnOffsets<<<unsigned((nd+256)/256),256>>>(nd,slots,sortedKeys.data,op.incidenceOffsets.data);check(cudaGetLastError());++launches;
        }
        buildMass(in,op,weights,M,massDiagonal,scratch);
        check(cudaDeviceSynchronize());
    }
private:
    template<class Scratch>void buildMass(const FlatEvaluationInput& in,FlatEvaluationOperator& op,Buffer<double>& weights,
                                          radia_hmatrix_cuda::DeviceMap& M,Buffer<double>& massDiagonal,Scratch& scratch){
        using namespace radia_hmatrix_cuda;
        const auto ne=in.elements,local=in.local,rule=in.rulePoints,nd=in.dofs;
        const auto triples=checkedProduct(ne,checkedProduct(local,local));
        Buffer<MassKey> keys,sortedKeys,uniqueKeys;Buffer<double> entries,sortedEntries,rowBound,bound;Buffer<int> counts,starts,runs,invalid;
        scratch(keys,triples);scratch(sortedKeys,triples);scratch(uniqueKeys,triples);scratch(entries,triples);scratch(sortedEntries,triples);
        scratch(counts,triples);scratch(starts,triples);scratch(runs,1);scratch(rowBound,nd);scratch(bound,1);scratch(invalid,1);
        flatMassEntries<<<unsigned((triples+255)/256),256>>>(ne,local,rule,nd,op.vertices.data,op.classes.data,op.dofs.data,op.points.data,op.order.data,
                                                               op.shapes.data,weights.data,keys.data,entries.data);
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
    const auto slots=checkedProduct(in.elements,in.local),triples=checkedProduct(in.elements,checkedProduct(in.local,in.local)),nd=in.dofs;
    int bits=1;while((std::size_t(1)<<bits)<nd)++bits;int massBits=1;while((std::size_t(1)<<massBits)<nd*nd)++massBits;
    std::size_t total=0,bytes=0;auto add=[&]{if(bytes>std::numeric_limits<std::size_t>::max()-total)throw std::overflow_error("Flat material budget overflow");total+=bytes;bytes=0;};
    check(cub::DeviceRadixSort::SortPairs(nullptr,bytes,(const unsigned*)nullptr,(unsigned*)nullptr,(const unsigned*)nullptr,(unsigned*)nullptr,int(slots),0,bits));add();
    check(cub::DeviceRadixSort::SortPairs(nullptr,bytes,(const MassKey*)nullptr,(MassKey*)nullptr,(const double*)nullptr,(double*)nullptr,int(triples),0,massBits));add();
    check(cub::DeviceRunLengthEncode::Encode(nullptr,bytes,(const MassKey*)nullptr,(MassKey*)nullptr,(int*)nullptr,(int*)nullptr,int(triples)));add();
    check(cub::DeviceScan::ExclusiveSum(nullptr,bytes,(const int*)nullptr,(int*)nullptr,int(triples)));add();
    check(cub::DeviceReduce::Max(nullptr,bytes,(const double*)nullptr,(double*)nullptr,int(nd)));add();
    return total;
}
// Live preflight before any material allocation: the dimension budget plus CUB
// scratch must fit the current free memory
// above the shared reserve. The failure names the limiting term; nothing
// spills to the host.
inline std::size_t requireFlatEvaluationDevice(const FlatEvaluationInput& in,std::size_t lawDoubles,FlatEvaluationBudget& budget){
    using namespace radia_hmatrix_cuda;
    budget=flatEvaluationBudget(in,lawDoubles);
    const auto scratch=flatEvaluationSortScratch(in),need=budget.total();
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
