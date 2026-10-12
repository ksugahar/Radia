// Resident symmetric HACApK leaf application; no FE or assembly substitution.
#pragma once
#include <cuda_runtime.h>
#include <vector>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <utility>

namespace radia_hmatrix_cuda {
inline void check(cudaError_t status) {
    if(status!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(status));
}
inline std::size_t checkedProduct(std::size_t a,std::size_t b) {
    if(b && a>std::numeric_limits<std::size_t>::max()/b)
        throw std::overflow_error("H-matrix allocation overflow");
    return a*b;
}
struct Direction {
    std::size_t input,output,inputs,outputs,rank,a,b,tmp,partial;
    bool transpose;
};
struct Leaf {
    std::size_t row,column,rows,columns,rank;
    // rank>0: a1[columns,rank], a2[rows,rank], column-major HACApK layout.
    // rank==0: a1[rows,columns], row-major dense layout.
    std::vector<double> a1,a2;
};
struct Plan {
    std::size_t n=0,tmpCount=0,partialCount=0;
    std::vector<Direction> directions;
    std::vector<double> coefficients,scale;
    std::vector<std::size_t> permutation,rowOffsets,partialIndices;
};
inline Plan pack(std::size_t n,const std::vector<Leaf>& leaves,
                 const std::vector<std::size_t>& permutation,const std::vector<double>& scale) {
    if(!n || permutation.size()!=n || scale.size()!=n)
        throw std::invalid_argument("Invalid H-matrix dimension or scaling");
    Plan p;p.n=n;p.permutation=permutation;p.scale=scale;
    std::vector<bool> seen(n,false);
    for(auto i:permutation){if(i>=n||seen[i])throw std::invalid_argument("Invalid permutation");seen[i]=true;}
    for(double s:scale)if(!std::isfinite(s)||s<=0)throw std::invalid_argument("Invalid physical scale");
    std::vector<std::vector<std::size_t>> contributions(n);
    for(const auto& leaf:leaves) {
        if(!leaf.rows||!leaf.columns||leaf.row>=n||leaf.column>=n||
           leaf.rows>n-leaf.row||leaf.columns>n-leaf.column)
            throw std::invalid_argument("Leaf range outside matrix");
        if(leaf.row>leaf.column)throw std::invalid_argument("Pack upper/diagonal leaves only");
        if(leaf.row==leaf.column&&leaf.rows!=leaf.columns)
            throw std::invalid_argument("Diagonal leaf must be square");
        if(leaf.row<leaf.column&&leaf.rows>leaf.column-leaf.row)
            throw std::invalid_argument("Upper leaf overlaps diagonal");
        if(leaf.a1.size()!=checkedProduct(leaf.columns,leaf.rank?leaf.rank:leaf.rows)||
           leaf.a2.size()!=(leaf.rank?checkedProduct(leaf.rows,leaf.rank):0))
            throw std::invalid_argument("Invalid factor storage");
        for(double v:leaf.a1)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite leaf");
        for(double v:leaf.a2)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite leaf");
        const auto a=p.coefficients.size();p.coefficients.insert(p.coefficients.end(),leaf.a1.begin(),leaf.a1.end());
        const auto b=p.coefficients.size();p.coefficients.insert(p.coefficients.end(),leaf.a2.begin(),leaf.a2.end());
        for(int transpose=0;transpose<(leaf.row==leaf.column?1:2);++transpose) {
            Direction d{transpose?leaf.row:leaf.column,transpose?leaf.column:leaf.row,
                transpose?leaf.rows:leaf.columns,transpose?leaf.columns:leaf.rows,
                leaf.rank,a,b,p.tmpCount,p.partialCount,transpose!=0};
            if(d.rank>std::numeric_limits<std::size_t>::max()-p.tmpCount||
               d.outputs>std::numeric_limits<std::size_t>::max()-p.partialCount)
                throw std::overflow_error("H-matrix workspace overflow");
            p.tmpCount+=d.rank;p.partialCount+=d.outputs;
            for(std::size_t i=0;i<d.outputs;++i)contributions[d.output+i].push_back(d.partial+i);
            p.directions.push_back(d);
        }
    }
    if(p.directions.empty())throw std::invalid_argument("Empty H-matrix");
    p.rowOffsets.push_back(0);
    for(const auto& row:contributions){p.partialIndices.insert(p.partialIndices.end(),row.begin(),row.end());p.rowOffsets.push_back(p.partialIndices.size());}
    return p;
}
__global__ void gatherInput(std::size_t n,const std::size_t* permutation,const double* scale,
                            const double* state,double* input) {
    auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(i<n)input[i]=scale[permutation[i]]*state[permutation[i]];
}
__global__ void projectRanks(const Direction* directions,std::size_t count,
                            const double* coefficients,const double* input,double* tmp) {
    const auto task=std::size_t(blockIdx.x);if(task>=count)return;
    const auto d=directions[task];
    for(auto k=std::size_t(threadIdx.x);k<d.rank;k+=blockDim.x) {
        const auto base=d.transpose?d.b:d.a;double value=0;
        for(std::size_t j=0;j<d.inputs;++j)value+=coefficients[base+j+d.inputs*k]*input[d.input+j];
        tmp[d.tmp+k]=value;
    }
}
__global__ void applyLeaves(const Direction* directions,std::size_t count,
                           const double* coefficients,const double* input,const double* tmp,double* partial) {
    const auto task=std::size_t(blockIdx.x);if(task>=count)return;
    const auto d=directions[task];
    for(auto i=std::size_t(threadIdx.x);i<d.outputs;i+=blockDim.x) {
        double value=0;
        if(d.rank) {
            const auto base=d.transpose?d.a:d.b;
            for(std::size_t k=0;k<d.rank;++k)value+=coefficients[base+i+d.outputs*k]*tmp[d.tmp+k];
        }else {
            for(std::size_t j=0;j<d.inputs;++j) {
                const auto offset=d.transpose?j*d.outputs+i:i*d.inputs+j;
                value+=coefficients[d.a+offset]*input[d.input+j];
            }
        }
        partial[d.partial+i]=value;
    }
}
__global__ void reduceOutput(std::size_t n,const std::size_t* offsets,const std::size_t* indices,
        const std::size_t* permutation,const double* scale,const double* partial,double* state) {
    auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=n)return;
    double value=0;
    // Stable leaf order, no cross-block atomic updates.
    for(auto j=offsets[i];j<offsets[i+1];++j)value+=partial[indices[j]];
    state[permutation[i]]=scale[permutation[i]]*value;
}
template<class T> struct Buffer {
    T* data=nullptr;std::size_t size=0;int device=0;
    Buffer()=default;Buffer(const Buffer&)=delete;Buffer& operator=(const Buffer&)=delete;
    ~Buffer(){if(data){int previous=device;cudaGetDevice(&previous);if(previous!=device)cudaSetDevice(device);cudaFree(data);if(previous!=device)cudaSetDevice(previous);}}
    void allocate(std::size_t n){check(cudaGetDevice(&device));size=n;if(n)check(cudaMalloc(reinterpret_cast<void**>(&data),checkedProduct(n,sizeof(T))));}
    void upload(const std::vector<T>& v){allocate(v.size());if(size)check(cudaMemcpy(data,v.data(),size*sizeof(T),cudaMemcpyHostToDevice));}
};
struct SparseMap {
    std::size_t rows=0,columns=0;
    std::vector<std::size_t> offsets,indices;
    std::vector<double> values;
    void validate()const {
        if(!rows||!columns||offsets.size()!=rows+1||offsets.front()!=0||offsets.back()!=values.size()||indices.size()!=values.size())
            throw std::invalid_argument("Invalid charge-map shape");
        for(std::size_t i=0;i<rows;++i)if(offsets[i]>offsets[i+1])throw std::invalid_argument("Invalid charge-map offsets");
        for(auto j:indices)if(j>=columns)throw std::invalid_argument("Invalid charge-map index");
        for(auto v:values)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite charge-map value");
    }
    SparseMap transpose()const {
        validate();SparseMap t;t.rows=columns;t.columns=rows;t.offsets.assign(columns+1,0);
        for(auto i:indices)++t.offsets[i+1];for(std::size_t i=0;i<columns;++i)t.offsets[i+1]+=t.offsets[i];
        t.indices.resize(indices.size());t.values.resize(values.size());auto next=t.offsets;
        for(std::size_t row=0;row<rows;++row)for(auto k=offsets[row];k<offsets[row+1];++k){auto q=next[indices[k]]++;t.indices[q]=row;t.values[q]=values[k];}return t;
    }
};
__global__ void sparseMapApply(std::size_t rows,const std::size_t* offsets,const std::size_t* indices,
                              const double* values,const double* x,double* y) {
    auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=rows)return;
    double result=0;for(auto k=offsets[i];k<offsets[i+1];++k)result+=values[k]*x[indices[k]];y[i]=result;
}
struct DeviceMap {
    Buffer<std::size_t> offsets,indices;Buffer<double> values;std::size_t rows=0;
    void upload(const SparseMap& map){rows=map.rows;offsets.upload(map.offsets);indices.upload(map.indices);values.upload(map.values);}
    void apply(const double* x,double* y,cudaStream_t stream){
        sparseMapApply<<<unsigned((rows+255)/256),256,0,stream>>>(rows,offsets.data,indices.data,values.data,x,y);check(cudaGetLastError());
    }
    void abandon(){offsets.data=nullptr;indices.data=nullptr;values.data=nullptr;}
};
} // namespace radia_hmatrix_cuda
#include "analytic_matrix_build.cuh"
namespace radia_hmatrix_cuda {
class Resident {
    Buffer<Direction> directions;Buffer<double> coefficients,scale,input,tmp,partial,state;
    Buffer<std::size_t> permutation,rowOffsets,partialIndices;
    DeviceMap chargeMap,chargeTranspose;Buffer<double> dofState;std::size_t dofs=0;
    cudaStream_t stream=nullptr;std::size_t n,count;bool initialized=false;int device=0;
    void valid(){int current;check(cudaGetDevice(&current));if(current!=device)throw std::runtime_error("H-matrix context belongs to another GPU");auto status=cudaStreamQuery(stream);if(status!=cudaSuccess&&status!=cudaErrorNotReady)check(status);}
public:
    std::size_t ownedBytes=0,setupPeakOwnedBytes=0,applications=0,kernelsPerApplication=4;
    std::size_t assemblyEntryBatches=0,assemblyEntryCount=0,assemblyLeafCount=0,assemblyLowrankCount=0;
    explicit Resident(const Plan& p,const SparseMap* map=nullptr):n(p.n),count(p.directions.size()) {
        check(cudaGetDevice(&device));
        if(count>2147483647u||n>std::size_t(2147483647u)*256)throw std::invalid_argument("CUDA grid range exceeded");
        auto bytes=[&](std::size_t size,std::size_t width){auto b=checkedProduct(size,width);if(b>std::numeric_limits<std::size_t>::max()-ownedBytes)throw std::overflow_error("Device allocation overflow");ownedBytes+=b;};
        bytes(count,sizeof(Direction));bytes(p.coefficients.size(),sizeof(double));bytes(n,sizeof(double));
        bytes(n,2*sizeof(double));bytes(p.tmpCount,sizeof(double));bytes(p.partialCount,sizeof(double));
        bytes(n,sizeof(std::size_t));bytes(p.rowOffsets.size(),sizeof(std::size_t));bytes(p.partialIndices.size(),sizeof(std::size_t));
        std::size_t free,total;check(cudaMemGetInfo(&free,&total));const auto reserve=std::max(std::size_t(2)<<30,total/5);
        if(free<=reserve||ownedBytes>free-reserve)throw std::runtime_error("GPU memory reserve would be exceeded");
        SparseMap transpose;
        if(map){map->validate();if(map->rows!=n||map->columns>2147483647u)throw std::invalid_argument("Charge-map dimension mismatch");transpose=map->transpose();
            bytes(map->offsets.size(),sizeof(std::size_t));bytes(transpose.offsets.size(),sizeof(std::size_t));
            bytes(map->indices.size(),2*sizeof(std::size_t));bytes(map->values.size(),2*sizeof(double));bytes(map->columns,sizeof(double));
            if(ownedBytes>free-reserve)throw std::runtime_error("GPU charge-map memory reserve would be exceeded");
            dofs=map->columns;kernelsPerApplication=6;chargeMap.upload(*map);chargeTranspose.upload(transpose);dofState.allocate(dofs);
        }
        directions.upload(p.directions);coefficients.upload(p.coefficients);scale.upload(p.scale);
        permutation.upload(p.permutation);rowOffsets.upload(p.rowOffsets);partialIndices.upload(p.partialIndices);
        input.allocate(n);tmp.allocate(p.tmpCount);partial.allocate(p.partialCount);state.allocate(n);
        check(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
    }
    void assemble(const AnalyticBlock& block){
        valid();if(initialized||count!=1||coefficients.size!=checkedProduct(n,n)||block.sources.size()!=n)
            throw std::invalid_argument("Analytic assembly requires a fresh single dense block");
        std::size_t extra=checkedProduct(block.sources.size(),sizeof(AnalyticSource));
        const auto quadBytes=checkedProduct(block.quadrature.size(),sizeof(QuadraturePoint));
        if(quadBytes>std::numeric_limits<std::size_t>::max()-extra-sizeof(int))throw std::overflow_error("Assembly workspace overflow");extra+=quadBytes+sizeof(int);
        std::size_t free,total;check(cudaMemGetInfo(&free,&total));const auto reserve=std::max(std::size_t(2)<<30,total/5);
        if(free<=reserve||extra>free-reserve)throw std::runtime_error("Assembly GPU reserve would be exceeded");
        Buffer<AnalyticSource> sources;Buffer<QuadraturePoint> quadrature;Buffer<int> invalid;
        sources.upload(block.sources);quadrature.upload(block.quadrature);invalid.upload(std::vector<int>{0});setupPeakOwnedBytes=ownedBytes+extra;
        assembleAnalyticConstant<<<unsigned((n*n+63)/64),64,0,stream>>>(sources.data,quadrature.data,n,coefficients.data,invalid.data);check(cudaGetLastError());
        analyticPhysicalScale<<<unsigned((n+255)/256),256,0,stream>>>(n,coefficients.data,scale.data,invalid.data);check(cudaGetLastError());
        int flag=0;check(cudaMemcpyAsync(&flag,invalid.data,sizeof(int),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));
        if(flag)throw std::runtime_error("Analytic assembly produced invalid energy");
        normalizeAnalytic<<<unsigned((n*n+255)/256),256,0,stream>>>(n,coefficients.data,scale.data);check(cudaGetLastError());check(cudaStreamSynchronize(stream));
    }
    ~Resident(){if(stream){int previous=device;cudaGetDevice(&previous);if(previous!=device)cudaSetDevice(device);
        auto status=cudaStreamQuery(stream);
        if(status!=cudaSuccess&&status!=cudaErrorNotReady){directions.data=nullptr;coefficients.data=nullptr;scale.data=nullptr;input.data=nullptr;tmp.data=nullptr;partial.data=nullptr;state.data=nullptr;permutation.data=nullptr;rowOffsets.data=nullptr;partialIndices.data=nullptr;chargeMap.abandon();chargeTranspose.abandon();dofState.data=nullptr;}
        else{cudaStreamSynchronize(stream);cudaStreamDestroy(stream);}if(previous!=device)cudaSetDevice(previous);}}
    void upload(const std::vector<double>& values){
        valid();const auto dimension=dofs?dofs:n;if(values.size()!=dimension)throw std::invalid_argument("State dimension mismatch");
        for(double v:values)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite state");
        check(cudaMemcpyAsync(dofs?dofState.data:state.data,values.data(),dimension*sizeof(double),cudaMemcpyHostToDevice,stream));
        check(cudaStreamSynchronize(stream));initialized=true;
    }
    std::size_t dimension()const{return dofs?dofs:n;}
    bool isDemag()const{return dofs!=0;}
    cudaStream_t executionStream(){valid();return stream;}
    // Internal device-to-device application for Krylov vectors. No host
    // vector transfer or synchronization; accepted uploaded state is untouched.
    void applyDevice(const double* x,double* y){
        valid();
        if(dofs)chargeMap.apply(x,state.data,stream);
        const double* source=dofs?state.data:x;
        gatherInput<<<unsigned((n+255)/256),256,0,stream>>>(n,permutation.data,scale.data,source,input.data);check(cudaGetLastError());
        projectRanks<<<unsigned(count),128,0,stream>>>(directions.data,count,coefficients.data,input.data,tmp.data);check(cudaGetLastError());
        applyLeaves<<<unsigned(count),128,0,stream>>>(directions.data,count,coefficients.data,input.data,tmp.data,partial.data);check(cudaGetLastError());
        double* target=dofs?state.data:y;
        reduceOutput<<<unsigned((n+255)/256),256,0,stream>>>(n,rowOffsets.data,partialIndices.data,permutation.data,scale.data,partial.data,target);check(cudaGetLastError());
        if(dofs)chargeTranspose.apply(state.data,y,stream);
    }
    void apply(std::size_t repeats){
        valid();if(!initialized)throw std::runtime_error("Upload initial state before apply");
        if(!repeats||repeats>10000)throw std::invalid_argument("Invalid repetition count");
        for(std::size_t k=0;k<repeats;++k){
            if(dofs)chargeMap.apply(dofState.data,state.data,stream);
            gatherInput<<<unsigned((n+255)/256),256,0,stream>>>(n,permutation.data,scale.data,state.data,input.data);check(cudaGetLastError());
            projectRanks<<<unsigned(count),128,0,stream>>>(directions.data,count,coefficients.data,input.data,tmp.data);check(cudaGetLastError());
            applyLeaves<<<unsigned(count),128,0,stream>>>(directions.data,count,coefficients.data,input.data,tmp.data,partial.data);check(cudaGetLastError());
            reduceOutput<<<unsigned((n+255)/256),256,0,stream>>>(n,rowOffsets.data,partialIndices.data,permutation.data,scale.data,partial.data,state.data);check(cudaGetLastError());
            if(dofs)chargeTranspose.apply(state.data,dofState.data,stream);
        }
        applications+=repeats;
    }
    std::vector<double> download(){
        valid();if(!initialized)throw std::runtime_error("No resident state");const auto dimension=dofs?dofs:n;std::vector<double> result(dimension);
        check(cudaMemcpyAsync(result.data(),dofs?dofState.data:state.data,dimension*sizeof(double),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));
        for(double v:result)if(!std::isfinite(v))throw std::runtime_error("Nonfinite resident result");return result;
    }
};
} // namespace radia_hmatrix_cuda
