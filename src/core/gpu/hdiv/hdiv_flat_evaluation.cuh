// Matrix-free device action of the TET HDiv material evaluation map E and
// E^T, and device construction of the reference mass M = E^T W E with its
// diagonal and Jacobi bound. NGSolve tabulates reference shapes per vertex-rank
// class on the host; only those tables, the element geometry (affine vertices
// or curved P2 nodes), the rule, DOF numbers and point numbers cross, and they
// stay resident: no global E is stored. J and det of a non-affine element are tabulated
// once per rule point at setup, since E^T needs them once per (DOF, slot, point, component).
#pragma once
#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_run_length_encode.cuh>
#include <cub/device/device_scan.cuh>
#include <cub/device/device_reduce.cuh>
#include <sstream>
#include "hmatrix_resident.cuh"
#include "hdiv_hex_gram.cuh"
namespace radia_hdiv_bh {
// Joint (row,column) mass key row*nd+column: 64-bit, exact for any int32 DOF count.
using MassKey=unsigned long long;
// Affine TET map in NGSolve vertex order: J columns are v_d - v_3 (d=0..2).
__device__ inline double affineTet(const double* V,double J[3][3]){
    for(int c=0;c<3;++c)for(int d=0;d<3;++d)J[c][d]=V[3*d+c]-V[9+c];
    return J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
}
// Element map Jacobian at rule row k: the affine TET (4 nodes, point
// independent), the isoparametric P2 TET (10 nodes) or the HEX Q2 lattice
// (27 nodes), Jac[k][c] = dX_k/dxi_c, the affine column convention. Returns det J.
__device__ inline double elementJacobian(const double* geometry,std::size_t nodes,std::size_t e,const double* ruleRow,double J[3][3]){
    if(nodes==4)return affineTet(geometry+12*e,J);
    if(nodes==27){double X[3];radia_hmatrix_cuda::hexMapJ(geometry+81*e,ruleRow,X,J);return rad_hdiv_cuda::det3(J);}
    if(nodes==5){const double* V=geometry+15*e;double X[3];radia_hdiv_typed::pyramidVertexMap(V,V+3,V+6,V+9,V+12,ruleRow,X,J);return rad_hdiv_cuda::det3(J);}
    double X[3];rad_hdiv_cuda::CurvedTetEval(reinterpret_cast<const double(*)[3]>(geometry+30*e),ruleRow[0],ruleRow[1],ruleRow[2],X,J);
    return rad_hdiv_cuda::det3(J);
}
// The one formula of E: contravariant Piola value of reference shape S, component c.
// The setup check, E x, E^T y and the mass all call it.
__device__ inline double piolaValue(const double J[3][3],double det,const double* S,int c){
    return (J[c][0]*S[0]+J[c][1]*S[1]+J[c][2]*S[2])/det;
}
// Element map table row of a non-affine element: J row by row, then det.
constexpr std::size_t elementMapDoubles=10;
__device__ inline const double (*tabulatedJacobian(const double* row))[3]{return reinterpret_cast<const double(*)[3]>(row);}
// Local slots in ascending DOF order (stable insertion sort; local <= 128 on the host).
__global__ void flatElementOrder(std::size_t ne,std::size_t local,const int* elementDofs,unsigned char* order){
    const auto e=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(e>=ne)return;
    const int* dofs=elementDofs+e*local;unsigned char* o=order+e*local;
    for(std::size_t i=0;i<local;++i){std::size_t j=i;while(j>0&&dofs[o[j-1]]>dofs[i]){o[j]=o[j-1];--j;}o[j]=(unsigned char)i;}
}
// Point measures and the validity of every E value.
__global__ void flatEvaluationEntries(std::size_t ne,std::size_t local,std::size_t rule,std::size_t nq,std::size_t classes,std::size_t nodes,
    const double* vertices,const int* elementClass,const int* elementPoints,const double* shapes,
    const double* ruleRows,double* weights,double* maps,int* invalid){
    const auto t=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(t>=ne*rule)return;
    const auto e=t/rule,k=t%rule;double J[3][3];const double det=elementJacobian(vertices,nodes,e,ruleRows+4*k,J);
    const int cls=elementClass[e],q=elementPoints[e*rule+k];
    if(!isfinite(det)||det==0.||cls<0||std::size_t(cls)>=classes||q<0||std::size_t(q)>=nq){atomicExch(invalid,1);return;}
    weights[q]=ruleRows[4*k+3]*fabs(det);
    if(nodes!=4){double* row=maps+elementMapDoubles*t;for(int c=0;c<3;++c)for(int d=0;d<3;++d)row[3*c+d]=J[c][d];row[9]=det;}
    const double* S=shapes+(std::size_t(cls)*rule+k)*local*3;
    for(int c=0;c<3;++c)for(std::size_t i=0;i<local;++i)if(!isfinite(piolaValue(J,det,S+3*i,c)))atomicExch(invalid,1);
}
// v = E x, one thread per material point, local slots in ascending DOF order.
__global__ void flatEvaluationApply(std::size_t ne,std::size_t local,std::size_t rule,std::size_t nq,std::size_t nodes,
    const double* vertices,const double* maps,const int* elementClass,const int* elementDofs,const int* elementPoints,const unsigned char* order,
    const double* shapes,const double* x,double* v){
    const auto t=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(t>=ne*rule)return;
    const auto e=t/rule,k=t%rule;double affine[3][3];double det=0.;const double (*J)[3]=affine;
    if(nodes==4)det=affineTet(vertices+12*e,affine);else{const double* row=maps+elementMapDoubles*t;J=tabulatedJacobian(row);det=row[9];}
    const auto q=std::size_t(elementPoints[t]);const int* dofs=elementDofs+e*local;const unsigned char* o=order+e*local;
    const double* S=shapes+(std::size_t(elementClass[e])*rule+k)*local*3;
    for(int c=0;c<3;++c){radia_hmatrix_cuda::CompensatedDeviceSum sum;
        for(std::size_t r=0;r<local;++r){const int i=o[r];sum.addProduct(piolaValue(J,det,S+3*i,c),x[dofs[i]]);}
        v[std::size_t(c)*nq+q]=sum.get();}
}
// g = E^T y, one thread per DOF over its incident (element, slot) entries in
// ascending element order: deterministic, with no floating-point atomics. For
// each component the terms follow ascending (element, point).
__global__ void flatEvaluationTransposeApply(std::size_t nd,std::size_t local,std::size_t rule,std::size_t nq,std::size_t nodes,
    const double* vertices,const double* maps,const int* elementClass,const int* elementPoints,const double* shapes,
    const std::size_t* offsets,const unsigned* incidence,const double* y,double* g){
    const auto d=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(d>=nd)return;radia_hmatrix_cuda::CompensatedDeviceSum sum;
    for(int c=0;c<3;++c)for(auto at=offsets[d];at<offsets[d+1];++at){
        const auto e=std::size_t(incidence[at])/local,i=std::size_t(incidence[at])%local;double affine[3][3];double det=0.;const double (*J)[3]=affine;
        if(nodes==4)det=affineTet(vertices+12*e,affine);
        const double* S=shapes+std::size_t(elementClass[e])*rule*local*3+3*i;
        for(std::size_t k=0;k<rule;++k){if(nodes!=4){const double* row=maps+elementMapDoubles*(e*rule+k);J=tabulatedJacobian(row);det=row[9];}
            sum.addProduct(piolaValue(J,det,S+k*local*3,c),y[std::size_t(c)*nq+std::size_t(elementPoints[e*rule+k])]);}}
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
__global__ void flatMassEntries(std::size_t ne,std::size_t local,std::size_t rule,std::size_t nd,std::size_t nodes,
    const double* vertices,const double* maps,const int* elementClass,const int* elementDofs,const int* elementPoints,const unsigned char* order,
    const double* shapes,const double* weights,MassKey* keys,double* entries){
    const auto t=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(t>=ne*local*local)return;
    const auto e=t/(local*local),r=(t/local)%local,s=t%local;double affine[3][3];double det=0.;const double (*J)[3]=affine;if(nodes==4)det=affineTet(vertices+12*e,affine);
    const int i=order[e*local+r],j=order[e*local+s];const int* dofs=elementDofs+e*local;
    const double* S=shapes+std::size_t(elementClass[e])*rule*local*3;radia_hmatrix_cuda::CompensatedDeviceSum sum;
    for(std::size_t k=0;k<rule;++k){if(nodes!=4){const double* row=maps+elementMapDoubles*(e*rule+k);J=tabulatedJacobian(row);det=row[9];}
        const double w=weights[elementPoints[e*rule+k]];const double* Sk=S+k*local*3;
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
// Rows [first, first+count) of the frozen material tangent on the reference mass pattern:
// T_rs = sum_q E_qr^T w_q (sec_q I + (dif_q - sec_q) n_q n_q^T) E_qs, the action of materialPointFields.
// One thread per row over its incident (element, slot) entries in ascending element order:
// deterministic, with no floating-point atomics. Groups add up, one launch after another.
__global__ void flatTangentRows(std::size_t first,std::size_t count,std::size_t local,std::size_t rule,std::size_t nq,std::size_t nodes,
    const double* vertices,const double* maps,const int* elementClass,const int* elementDofs,const int* elementPoints,const double* shapes,
    const std::size_t* incidenceOffsets,const unsigned* incidence,const std::size_t* massOffsets,const std::size_t* massIndices,
    const double* v,const double* amplitude,const double* law,const double* weights,double* values,int* invalid){
    const auto index=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(index>=count)return;const auto d=first+index;
    for(auto at=incidenceOffsets[d];at<incidenceOffsets[d+1];++at){
        const auto e=std::size_t(incidence[at])/local,i=std::size_t(incidence[at])%local;double affine[3][3];double det=0.;const double (*J)[3]=affine;
        if(nodes==4)det=affineTet(vertices+12*e,affine);
        const double* S=shapes+std::size_t(elementClass[e])*rule*local*3;const int* dofs=elementDofs+e*local;
        for(std::size_t j=0;j<local;++j){radia_hmatrix_cuda::CompensatedDeviceSum sum;
            for(std::size_t k=0;k<rule;++k){if(nodes!=4){const double* row=maps+elementMapDoubles*(e*rule+k);J=tabulatedJacobian(row);det=row[9];}
                const auto q=std::size_t(elementPoints[e*rule+k]);const double* Sk=S+k*local*3;const double w=weights[q],sec=law[q],dif=law[nq+q],a=amplitude[q];
                double direction[3],column[3];radia_hmatrix_cuda::CompensatedDeviceSum dot;
                for(int c=0;c<3;++c){direction[c]=a>0?v[std::size_t(c)*nq+q]/a:0.;column[c]=piolaValue(J,det,Sk+3*j,c);dot.addProduct(direction[c],column[c]);}
                for(int c=0;c<3;++c)sum.addProduct(piolaValue(J,det,Sk+3*i,c),w*(sec*column[c]+(dif-sec)*direction[c]*dot.get()));}
            std::size_t slot=0;if(!flatMassFind(massOffsets,massIndices,d,std::size_t(dofs[j]),slot)){atomicExch(invalid,1);return;}
            values[slot]+=sum.get();}
    }
}
// Sums one group's transpose into the mesh vector (mixed meshes only).
__global__ void flatEvaluationAdd(std::size_t nd,const double* partial,double* g){
    const auto d=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(d<nd)g[d]+=partial[d];
}
// Resident tables of the matrix-free E of one element type: geometry-derived element state and
// the DOF incidence (entries e*local+slot ascending per DOF). No E value is stored.
struct FlatEvaluationGroup {
    template<class T>using Buffer=radia_hmatrix_cuda::Buffer<T>;
    std::size_t ne=0,local=0,rule=0,nq=0,nd=0,nodes=4;
    Buffer<double> vertices,shapes,ruleRows,maps;Buffer<int> classes,dofs,points;Buffer<unsigned char> order;
    Buffer<std::size_t> incidenceOffsets;Buffer<unsigned> incidence;
    std::size_t bytes()const{
        return (vertices.size+shapes.size+ruleRows.size+maps.size)*sizeof(double)+(classes.size+dofs.size+points.size)*sizeof(int)+order.size+
               incidenceOffsets.size*sizeof(std::size_t)+incidence.size*sizeof(unsigned);
    }
    void apply(const double* x,double* v,cudaStream_t stream){
        flatEvaluationApply<<<unsigned((ne*rule+127)/128),128,0,stream>>>(ne,local,rule,nq,nodes,vertices.data,maps.data,classes.data,dofs.data,points.data,order.data,shapes.data,x,v);
        radia_hmatrix_cuda::check(cudaGetLastError());
    }
    void applyTranspose(const double* y,double* g,cudaStream_t stream){
        flatEvaluationTransposeApply<<<unsigned((nd+127)/128),128,0,stream>>>(nd,local,rule,nq,nodes,vertices.data,maps.data,classes.data,points.data,shapes.data,
            incidenceOffsets.data,incidence.data,y,g);
        radia_hmatrix_cuda::check(cudaGetLastError());
    }
    void addFrozenTangent(const radia_hmatrix_cuda::DeviceMap& pattern,const double* v,const double* amplitude,const double* law,const double* weights,double* values,int* invalid){
        // One launch per row batch: a row costs local*rule products per incident slot.
        constexpr std::size_t batch=8192;
        for(std::size_t first=0;first<nd;first+=batch){const auto rows=std::min(batch,nd-first);
            flatTangentRows<<<unsigned((rows+63)/64),64>>>(first,rows,local,rule,nq,nodes,vertices.data,maps.data,classes.data,dofs.data,points.data,shapes.data,
                incidenceOffsets.data,incidence.data,pattern.offsets.data,pattern.indices.data,v,amplitude,law,weights,values,invalid);
            radia_hmatrix_cuda::check(cudaGetLastError());}
    }
};
// Tet, pyramid and hex groups write disjoint points of E x; E^T y sums the groups in order.
struct FlatEvaluationOperator {
    static constexpr std::size_t maxGroups=3;
    FlatEvaluationGroup groups[maxGroups];std::size_t count=0;radia_hmatrix_cuda::Buffer<double> partial;
    std::size_t bytes()const{std::size_t sum=partial.size*sizeof(double);for(std::size_t i=0;i<count;++i)sum+=groups[i].bytes();return sum;}
    void apply(const double* x,double* v,cudaStream_t stream){for(std::size_t i=0;i<count;++i)groups[i].apply(x,v,stream);}
    void addFrozenTangent(const radia_hmatrix_cuda::DeviceMap& pattern,const double* v,const double* amplitude,const double* law,const double* weights,double* values,int* invalid){
        for(std::size_t i=0;i<count;++i)groups[i].addFrozenTangent(pattern,v,amplitude,law,weights,values,invalid);
    }
    void applyTranspose(const double* y,double* g,cudaStream_t stream){
        groups[0].applyTranspose(y,g,stream);
        for(std::size_t i=1;i<count;++i){groups[i].applyTranspose(y,partial.data,stream);
            const auto nd=groups[i].nd;flatEvaluationAdd<<<unsigned((nd+127)/128),128,0,stream>>>(nd,partial.data,g);radia_hmatrix_cuda::check(cudaGetLastError());}
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
    // One ledger row per payload: the groups of a mixed mesh add up.
    void transfer(const char* payload,const char* direction,std::size_t bytes){
        for(auto& row:transfers)if(std::string(row.payload)==payload&&std::string(row.direction)==direction){row.bytes+=bytes;return;}
        transfers.push_back({payload,direction,bytes});
    }
    void build(const FlatEvaluationGroups& in,FlatEvaluationOperator& op,Buffer<double>& weights,radia_hmatrix_cuda::DeviceMap& M,Buffer<double>& massDiagonal){
        using namespace radia_hmatrix_cuda;
        const auto nq=in.at(0).points,nd=in[0].dofs;
        Buffer<int> invalid;
        auto upload=[&](auto& buffer,const auto& values,const char* payload,bool resident){buffer.upload(values);const auto bytes=checkedProduct(values.size(),sizeof(values[0]));
            transfer(payload,"H2D",bytes);if(!resident)transientBytes+=bytes;};
        auto scratch=[&](auto& buffer,std::size_t count){buffer.allocate(count);transientBytes+=checkedProduct(count,sizeof(*buffer.data));};
        upload(invalid,HostVector<int>{0},"evaluation_invalid_flag",false);
        weights.allocate(nq);
        int bits=1;while((std::size_t(1)<<bits)<nd)++bits;
        op.count=in.size();
        if(in.size()>1)op.partial.allocate(nd);
        for(std::size_t g=0;g<in.size();++g){const auto& group=in[g];auto& part=op.groups[g];
            const auto ne=group.elements,local=group.local,rule=group.rulePoints,slots=checkedProduct(ne,local);
            part.ne=ne;part.local=local;part.rule=rule;part.nq=nq;part.nd=nd;part.nodes=group.nodes;
            upload(part.vertices,group.vertices,group.nodes==27?"element_q2_nodes_m":group.nodes==10?"element_p2_nodes_m":group.nodes==5?"element_pyramid_vertices_m":"element_vertices_m",true);
            upload(part.shapes,group.shapes,"reference_shapes",true);upload(part.ruleRows,group.rule,"material_rule",true);
            upload(part.classes,group.elementClass,"element_class",true);upload(part.dofs,group.elementDofs,"element_dofs",true);upload(part.points,group.elementPoints,"element_points",true);
            part.order.allocate(slots);if(part.nodes!=4)part.maps.allocate(checkedProduct(elementMapDoubles,checkedProduct(ne,rule)));
            flatElementOrder<<<unsigned((ne+127)/128),128>>>(ne,local,part.dofs.data,part.order.data);check(cudaGetLastError());++launches;
            flatEvaluationEntries<<<unsigned((ne*rule+127)/128),128>>>(ne,local,rule,nq,group.classes,part.nodes,part.vertices.data,part.classes.data,part.points.data,
                part.shapes.data,part.ruleRows.data,weights.data,part.maps.data,invalid.data);check(cudaGetLastError());++launches;
        }
        int failed=0;check(cudaMemcpy(&failed,invalid.data,sizeof(int),cudaMemcpyDeviceToHost));
        transfer("evaluation_flag","D2H",sizeof(int));
        if(failed)throw std::runtime_error("Device HDiv evaluation map has a degenerate element, invalid index or nonfinite shape");
        for(std::size_t g=0;g<in.size();++g){auto& part=op.groups[g];const auto slots=checkedProduct(part.ne,part.local);
            // DOF incidence: stable sort of the slot DOF numbers, so entries ascend per DOF.
            Buffer<unsigned> keys,sortedKeys,entries;scratch(keys,slots);scratch(sortedKeys,slots);scratch(entries,slots);
            part.incidence.allocate(slots);part.incidenceOffsets.allocate(nd+1);
            slotKeys<<<unsigned((slots+255)/256),256>>>(slots,part.dofs.data,keys.data,entries.data);check(cudaGetLastError());++launches;
            std::size_t bytes=0;check(cub::DeviceRadixSort::SortPairs(nullptr,bytes,keys.data,sortedKeys.data,entries.data,part.incidence.data,int(slots),0,bits));
            Buffer<unsigned char> temporary;scratch(temporary,bytes);
            check(cub::DeviceRadixSort::SortPairs(temporary.data,bytes,keys.data,sortedKeys.data,entries.data,part.incidence.data,int(slots),0,bits));++launches;
            columnOffsets<<<unsigned((nd+256)/256),256>>>(nd,slots,sortedKeys.data,part.incidenceOffsets.data);check(cudaGetLastError());++launches;
            check(cudaDeviceSynchronize());
        }
        buildMass(in,op,weights,M,massDiagonal,scratch);
        check(cudaDeviceSynchronize());
    }
private:
    template<class Scratch>void buildMass(const FlatEvaluationGroups& in,FlatEvaluationOperator& op,Buffer<double>& weights,
                                          radia_hmatrix_cuda::DeviceMap& M,Buffer<double>& massDiagonal,Scratch& scratch){
        using namespace radia_hmatrix_cuda;
        const auto nd=in[0].dofs;std::size_t triples=0;
        for(const auto& group:in)triples+=checkedProduct(group.elements,checkedProduct(group.local,group.local));
        Buffer<MassKey> keys,sortedKeys,uniqueKeys;Buffer<double> entries,sortedEntries,rowBound,bound;Buffer<int> counts,starts,runs,invalid;
        scratch(keys,triples);scratch(sortedKeys,triples);scratch(uniqueKeys,triples);scratch(entries,triples);scratch(sortedEntries,triples);
        scratch(counts,triples);scratch(starts,triples);scratch(runs,1);scratch(rowBound,nd);scratch(bound,1);scratch(invalid,1);
        std::size_t at=0;
        for(std::size_t g=0;g<op.count;++g){auto& part=op.groups[g];const auto count=part.ne*part.local*part.local;
            flatMassEntries<<<unsigned((count+255)/256),256>>>(part.ne,part.local,part.rule,nd,part.nodes,part.vertices.data,part.maps.data,part.classes.data,part.dofs.data,part.points.data,part.order.data,
                                                               part.shapes.data,weights.data,keys.data+at,entries.data+at);
            check(cudaGetLastError());++launches;at+=count;}
        int bits=1;while((std::size_t(1)<<bits)<nd*nd)++bits;std::size_t bytes=0;
        check(cub::DeviceRadixSort::SortPairs(nullptr,bytes,keys.data,sortedKeys.data,entries.data,sortedEntries.data,int(triples),0,bits));
        Buffer<unsigned char> sortScratch;scratch(sortScratch,bytes);
        check(cub::DeviceRadixSort::SortPairs(sortScratch.data,bytes,keys.data,sortedKeys.data,entries.data,sortedEntries.data,int(triples),0,bits));++launches;
        bytes=0;check(cub::DeviceRunLengthEncode::Encode(nullptr,bytes,sortedKeys.data,uniqueKeys.data,counts.data,runs.data,int(triples)));
        Buffer<unsigned char> encodeScratch;scratch(encodeScratch,bytes);
        check(cub::DeviceRunLengthEncode::Encode(encodeScratch.data,bytes,sortedKeys.data,uniqueKeys.data,counts.data,runs.data,int(triples)));++launches;
        int nnz=0;check(cudaMemcpy(&nnz,runs.data,sizeof(int),cudaMemcpyDeviceToHost));transfer("mass_nonzeros","D2H",sizeof(int));
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
        int status=0;check(cudaMemcpy(&status,invalid.data,sizeof(int),cudaMemcpyDeviceToHost));transfer("mass_flag","D2H",sizeof(int));
        check(cudaMemcpy(&massBound,bound.data,sizeof(double),cudaMemcpyDeviceToHost));transfer("mass_spectral_bound","D2H",sizeof(double));
        if(status==1)throw std::runtime_error("Device HDiv reference mass has a missing or nonpositive diagonal");
        if(status)throw std::runtime_error("Device HDiv reference mass is not symmetric");
        if(!std::isfinite(massBound)||massBound<1.)throw std::runtime_error("Invalid device mass spectral bound");
    }
};
// Host bounds before any upload: every index is checked against its range.
// No element, DOF or point cap: int32 tables and the int item counts of the
// CUB sort/encode/scan calls are the structural bound, and the live device
// budget (requireFlatEvaluationDevice) decides what fits.
inline void validateFlatEvaluation(const FlatEvaluationGroups& groups){
    using radia_hmatrix_cuda::checkedProduct;constexpr std::size_t count=std::numeric_limits<int>::max();
    if(groups.empty()||groups.size()>FlatEvaluationOperator::maxGroups)throw std::invalid_argument("Flat evaluation takes one to three element-type groups");
    const auto points=groups[0].points,dofs=groups[0].dofs;std::size_t covered=0,triples=0;
    for(const auto& in:groups){
        const auto ne=in.elements,local=in.local,rule=in.rulePoints;
        // TET: 24 vertex-rank classes, local <= 64. HEX: classes follow the mesh's vertex orders, local <= 128. Pyramid: 120 orders.
        const bool wide=in.nodes==27||in.nodes==5;
        if(!ne||!local||local>(wide?128u:64u)||!rule||rule>4096||!in.classes||in.classes>(wide?1024u:24u)||!in.dofs||
           (in.nodes!=4&&in.nodes!=10&&!wide)||in.points!=points||in.dofs!=dofs||(groups.size()>1&&in.nodes==10))
            throw std::invalid_argument("Flat evaluation tables outside the affine/curved P2 TET, straight pyramid and Q2 HEX scope");
        covered+=checkedProduct(ne,rule);triples+=checkedProduct(ne,checkedProduct(local,local));
        if(ne>count||in.dofs>count||in.points>count)throw std::invalid_argument("Flat evaluation tables exceed int32 element, DOF or point numbering");
        if(in.vertices.size()!=3*in.nodes*ne||in.elementClass.size()!=ne||in.elementDofs.size()!=ne*local||in.elementPoints.size()!=ne*rule||
           in.rule.size()!=4*rule||in.shapes.size()!=in.classes*rule*local*3)throw std::invalid_argument("Flat evaluation table shape mismatch");
        for(double v:in.vertices)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite evaluation vertex");
        for(double v:in.shapes)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite reference shape");
        for(std::size_t k=0;k<rule;++k)if(!std::isfinite(in.rule[4*k])||!std::isfinite(in.rule[4*k+1])||!std::isfinite(in.rule[4*k+2])||!(in.rule[4*k+3]>0)||!std::isfinite(in.rule[4*k+3]))
            throw std::invalid_argument("Invalid material reference rule");
        for(int c:in.elementClass)if(c<0||std::size_t(c)>=in.classes)throw std::invalid_argument("Evaluation class index out of range");
        for(int d:in.elementDofs)if(d<0||std::size_t(d)>=in.dofs)throw std::invalid_argument("Evaluation DOF index out of range");
        if(!in.periodicTrace)for(std::size_t e=0;e<ne;++e)for(std::size_t i=0;i<local;++i)for(std::size_t j=0;j<i;++j)
            if(in.elementDofs[e*local+i]==in.elementDofs[e*local+j])throw std::invalid_argument("Duplicate evaluation local DOF");
    }
    if(covered!=points)throw std::invalid_argument("Evaluation points must cover each material point once");
    if(triples>count)throw std::invalid_argument("Flat mass triples exceed the checked int CUB item count");
    std::vector<unsigned char> seen(points,0);
    for(const auto& in:groups)for(int q:in.elementPoints){if(q<0||std::size_t(q)>=points||seen[q])throw std::invalid_argument("Evaluation points must cover each material point once");seen[q]=1;}
}
// CUB scratch of FlatEvaluationBuilder, queried with null storage (no allocation).
// All of it is live together until build() returns; the scan is sized for its
// largest possible run count.
inline std::size_t flatEvaluationSortScratch(const FlatEvaluationGroups& groups){
    using namespace radia_hmatrix_cuda;
    const auto nd=groups[0].dofs;std::size_t triples=0;
    int bits=1;while((std::size_t(1)<<bits)<nd)++bits;int massBits=1;while((std::size_t(1)<<massBits)<nd*nd)++massBits;
    std::size_t total=0,bytes=0;auto add=[&]{if(bytes>std::numeric_limits<std::size_t>::max()-total)throw std::overflow_error("Flat material budget overflow");total+=bytes;bytes=0;};
    for(const auto& in:groups){const auto slots=checkedProduct(in.elements,in.local);triples+=checkedProduct(in.elements,checkedProduct(in.local,in.local));
        check(cub::DeviceRadixSort::SortPairs(nullptr,bytes,(const unsigned*)nullptr,(unsigned*)nullptr,(const unsigned*)nullptr,(unsigned*)nullptr,int(slots),0,bits));add();}
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
inline std::size_t requireFlatEvaluationDevice(const FlatEvaluationGroups& in,std::size_t lawDoubles,FlatEvaluationBudget& budget){
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
