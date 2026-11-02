// Device construction of the affine-TET HDiv material evaluation map E and E^T.
// NGSolve tabulates reference shapes per vertex-rank class on the host; only
// those tables, the element vertices, DOF numbers and point numbers cross.
#pragma once
#include <cub/device/device_radix_sort.cuh>
#include "hmatrix_resident.cuh"
namespace radia_hdiv_bh {
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
struct FlatEvaluationTransfer {const char* payload;const char* direction;std::size_t bytes;};
// Builds E (3*nq x nd, component-major rows) and E^T into caller-owned device maps.
struct FlatEvaluationBuilder {
    template<class T>using Buffer=radia_hmatrix_cuda::Buffer<T>;
    radia_hmatrix_cuda::HostVector<FlatEvaluationTransfer> transfers;std::size_t launches=0,nonzeros=0;
    void build(const FlatEvaluationInput& in,radia_hmatrix_cuda::DeviceMap& E,radia_hmatrix_cuda::DeviceMap& ET,Buffer<double>& weights){
        using namespace radia_hmatrix_cuda;
        const auto ne=in.elements,local=in.local,rule=in.rulePoints,nq=in.points,nd=in.dofs;
        nonzeros=checkedProduct(checkedProduct(3*nq,local),1);
        Buffer<double> vertices,shapes,ruleRows;Buffer<int> classes,dofs,points,invalid;
        auto upload=[&](auto& buffer,const auto& values,const char* payload){buffer.upload(values);
            transfers.push_back({payload,"H2D",checkedProduct(values.size(),sizeof(values[0]))});};
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
        Buffer<unsigned> keys,sortedKeys,entries,sortedEntries;keys.allocate(nonzeros);sortedKeys.allocate(nonzeros);entries.allocate(nonzeros);sortedEntries.allocate(nonzeros);
        columnKeys<<<unsigned((nonzeros+255)/256),256>>>(nonzeros,E.indices.data,keys.data,entries.data);check(cudaGetLastError());++launches;
        int bits=1;while((std::size_t(1)<<bits)<nd)++bits;std::size_t scratch=0;
        check(cub::DeviceRadixSort::SortPairs(nullptr,scratch,keys.data,sortedKeys.data,entries.data,sortedEntries.data,int(nonzeros),0,bits));
        Buffer<unsigned char> temporary;temporary.allocate(scratch);
        check(cub::DeviceRadixSort::SortPairs(temporary.data,scratch,keys.data,sortedKeys.data,entries.data,sortedEntries.data,int(nonzeros),0,bits));++launches;
        ET.rows=nd;ET.offsets.allocate(nd+1);ET.indices.allocate(nonzeros);ET.values.allocate(nonzeros);
        transposeEntries<<<unsigned((nonzeros+255)/256),256>>>(nonzeros,local,sortedEntries.data,E.values.data,ET.indices.data,ET.values.data);check(cudaGetLastError());++launches;
        columnOffsets<<<unsigned((nd+256)/256),256>>>(nd,nonzeros,sortedKeys.data,ET.offsets.data);check(cudaGetLastError());++launches;
        check(cudaDeviceSynchronize());
    }
};
// Host bounds before any upload: every index is checked against its range.
inline void validateFlatEvaluation(const FlatEvaluationInput& in){
    const auto ne=in.elements,local=in.local,rule=in.rulePoints;
    if(!ne||ne>4096||!local||local>64||!rule||rule>4096||!in.classes||in.classes>24||!in.dofs||in.dofs>16384||in.points!=ne*rule||in.points>65536)
        throw std::invalid_argument("Flat evaluation tables outside the bounded affine-TET scope");
    if(radia_hmatrix_cuda::checkedProduct(3*in.points,local)>4000000)throw std::invalid_argument("Flat evaluation entry budget exceeded");
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
}
