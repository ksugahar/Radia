// Device magnetization field of a resident TET solution, with its mirror and cyclic images: the CPU direct route
// (HDivFieldEvaluator). Affine hosts take the closed-form cell and triangle fields of their physical polynomials about the
// host centroid; curved hosts integrate their P2 maps with the Gram's Duffy rule at every probe. E_T(x) = T E(T^-1 x),
// added with the image sign. Only probes, the state vector, geometry and host lists cross H2D and the requested field D2H;
// the charges are the resident charge map's product on the device.
#pragma once
#include "hdiv_hex_field.cuh"
namespace radia_hmatrix_cuda {
// One row per host: affine cell 20 cubic coefficients, affine face constant, slope and hessian, curved reference modes.
constexpr std::size_t flatFieldStride=20;
__device__ inline int flatCurvedCellMode(const int* e){return e[0]?1:e[1]?2:e[2]?3:0;}
__device__ inline int flatCurvedFaceMode(const int* e){
    if(e[0]==0)return e[1];
    return e[0]==1?3+e[1]:5;
}
__global__ void flatFieldPolynomials(FlatDevice d,std::size_t cellCount,std::size_t hosts,const std::size_t* first,const int* order,const double* charges,double* poly,int* invalid){
    const auto g=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(g>=hosts)return;
    const bool cell=g<cellCount;double* out=poly+g*flatFieldStride;
    for(std::size_t i=0;i<flatFieldStride;++i)out[i]=0.0;
    for(std::size_t k=first[g];k<first[g+1];++k){const int a=order[k];const double q=charges[a];const FlatCharge& c=d.charges[a];
        if(d.curved)out[cell?flatCurvedCellMode(c.expo):flatCurvedFaceMode(c.expo)]+=q;
        else if(cell){out[hexPolyIdx(0,0,0)]+=q*c.a;out[hexPolyIdx(1,0,0)]+=q*c.b[0];out[hexPolyIdx(0,1,0)]+=q*c.b[1];out[hexPolyIdx(0,0,1)]+=q*c.b[2];}
        else{out[0]+=q*c.a;for(int i=0;i<3;++i)out[1+i]+=q*c.b[i];for(int i=0;i<9;++i)out[4+i]+=q*c.c[i];}}
    for(std::size_t i=0;i<flatFieldStride;++i)if(!isfinite(out[i]))atomicExch(invalid,1);
}
__device__ inline void flatFieldHost(const FlatDevice& d,std::size_t cellCount,std::size_t g,const double* poly,const double* r,double* out){
    const double* p=poly+g*flatFieldStride;
    bool empty=true;for(std::size_t i=0;i<flatFieldStride&&empty;++i)empty=p[i]==0.0;if(empty)return;
    const bool cell=g<cellCount;
    if(d.curved){
        if(cell){rad_hdiv::FieldCurvedTetSource source;const double* X=d.cellNodes+30*g;
            for(int i=0;i<10;++i)for(int k=0;k<3;++k)source.nodes[i][k]=X[3*i+k];
            for(int i=0;i<4;++i)source.coefficient[i]=p[i];
            rad_hdiv_cuda::AddCurvedTet(source,r,out,d.gl,d.gw,d.nq);
        }else{rad_hdiv::FieldCurvedTriSource source;const double* X=d.faceNodes+18*(g-cellCount);
            for(int i=0;i<6;++i)for(int k=0;k<3;++k)source.nodes[i][k]=X[3*i+k];
            for(int i=0;i<6;++i)source.coefficient[i]=p[i];
            rad_hdiv_cuda::AddCurvedTriangle(source,r,out,d.gl,d.gw,d.nq);}
        return;
    }
    double value[3]={0.0,0.0,0.0};
    if(cell){const double* X=d.cellNodes+12*g;double c[3],V[4][3],lr[3];
        for(int k=0;k<3;++k){c[k]=0.25*(X[k]+X[3+k]+X[6+k]+X[9+k]);lr[k]=r[k]-c[k];for(int i=0;i<4;++i)V[i][k]=X[3*i+k]-c[k];}
        rad_hdiv_cuda::TetVolFieldCubic(V,lr,p,value);
    }else{const double* X=d.faceNodes+9*(g-cellCount);double c[3],V[3][3],lr[3],hessian[3][3];
        for(int k=0;k<3;++k){c[k]=(X[k]+X[3+k]+X[6+k])/3.0;lr[k]=r[k]-c[k];for(int i=0;i<3;++i)V[i][k]=X[3*i+k]-c[k];}
        for(int i=0;i<3;++i)for(int k=0;k<3;++k)hessian[i][k]=p[4+3*i+k];
        rad_hdiv_cuda::QuadTriField(V,lr,p[0],p+1,hessian,value);}
    for(int k=0;k<3;++k)out[k]+=value[k];
}
// One block per (probe, term); its threads share the hosts, then reduce.
__global__ void flatFieldRows(FlatDevice d,std::size_t cellCount,std::size_t hosts,const double* poly,const double* probes,std::size_t firstRow,double* rows){
    __shared__ double shared[hexFieldThreads];
    const auto row=firstRow+std::size_t(blockIdx.x);const std::size_t terms=1+std::size_t(d.imageCount),probe=row/terms;const int img=int(row%terms);
    double rm[3];if(img)hexImagePoint(d.images[img-1],probes+3*probe,rm);else for(int k=0;k<3;++k)rm[k]=probes[3*probe+k];
    double sum[3]={0.0,0.0,0.0};
    for(std::size_t g=threadIdx.x;g<hosts;g+=hexFieldThreads)flatFieldHost(d,cellCount,g,poly,rm,sum);
    for(int k=0;k<3;++k)sum[k]=hexFieldBlockSum(sum[k],shared);
    if(threadIdx.x==0){double image[3];
        if(img)hexImageVector(d.images[img-1],sum,image);else for(int k=0;k<3;++k)image[k]=sum[k];
        const double sign=img?d.images[img-1].sign:1.0;
        for(int k=0;k<3;++k)rows[3*row+k]=sign*image[k];}
}
// x is the prepared FE state; the Resident supplies q = B x on its own stream and device.
inline HexFieldResult flatDeviceField(const FlatBlock& block,Resident& resident,const HostVector<double>& x,const HostVector<double>& points){
    if(block.hex||!block.hostTypes.empty())throw std::invalid_argument("Device field evaluation has no route for this geometry family");
    if(points.empty()||points.size()%3||x.size()!=resident.fieldStateSize())throw std::invalid_argument("Device field evaluation requires (points,3) probes and the prepared state");
    for(double v:points)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite field probe");
    for(double v:x)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite field state");
    const std::size_t probes=points.size()/3,terms=1+block.images.size(),rowCount=checkedProduct(probes,terms),n=block.charges();
    if(rowCount>std::size_t(INT_MAX))throw std::invalid_argument("Device field probes exceed the launch range");
    const std::size_t cellCount=block.cellCount(),hosts=cellCount+block.faceCount();
    HexFieldResult result;result.route=block.curved?"device-curved-tet-direct":"device-flat-tet-polynomial";
    cudaStream_t stream=resident.executionStream();
    auto record=[&](const char* payload,const char* direction,std::size_t bytes){result.transfers.push_back({payload,direction,bytes});};
    // The charges of each host, in the order of the charge list: the coefficients sum deterministically.
    HostVector<std::size_t> first(hosts+1,0,HostFaultAllocator<std::size_t>("flat.field-first"));
    HostVector<int> order(n,0,HostFaultAllocator<int>("flat.field-order"));
    for(std::size_t a=0;a<n;++a){const std::size_t g=std::size_t(block.host[a])+(block.kind[a]?cellCount:0);
        if(block.host[a]<0||g>=hosts)throw std::invalid_argument("Invalid field charge host");++first[g+1];}
    for(std::size_t g=0;g<hosts;++g)first[g+1]+=first[g];
    {HostVector<std::size_t> next(first.begin(),first.end()-1,HostFaultAllocator<std::size_t>("flat.field-next"));
        for(std::size_t a=0;a<n;++a)order[next[std::size_t(block.host[a])+(block.kind[a]?cellCount:0)]++]=int(a);}
    result.ownedBytes=flatSum(FlatPrepared::workspaceBytes(block),flatSum(checkedProduct(rowCount,3*sizeof(double)),
        checkedProduct(checkedProduct(hosts,flatFieldStride)+n+x.size()+2*points.size(),sizeof(double))));
    flatReserve(result.ownedBytes,"Device field evaluation would exceed the GPU reserve");
    uploadHexFieldFarGauss();record("field_far_gauss","H2D",sizeof(rad_hdiv_cuda::GaussLegendreTable));
    FlatPrepared prepared;int code=0;prepared.prepare(block,stream,result.transfers,result.launches,code,false);
    const FlatDevice d=prepared.device();
    Buffer<std::size_t> hostFirst;Buffer<int> hostOrder;Buffer<double> state,charges,probe,poly,rows,field;
    hostFirst.upload(first);record("field_host_first","H2D",checkedProduct(first.size(),sizeof(std::size_t)));
    hostOrder.upload(order);record("field_host_order","H2D",checkedProduct(order.size(),sizeof(int)));
    state.upload(x);record("field_state","H2D",checkedProduct(x.size(),sizeof(double)));
    probe.upload(points);record("field_probes","H2D",checkedProduct(points.size(),sizeof(double)));
    charges.allocate(n);resident.chargeProducts(state.data,charges.data);
    poly.allocate(checkedProduct(hosts,flatFieldStride));
    flatFieldPolynomials<<<unsigned((hosts+127)/128),128,0,stream>>>(d,cellCount,hosts,hostFirst.data,hostOrder.data,charges.data,poly.data,prepared.invalid.data);check(cudaGetLastError());++result.launches;
    if(prepared.flag(stream,result.transfers,"field_polynomial_flag")!=0)throw std::runtime_error("Device field polynomial is nonfinite");
    rows.allocate(checkedProduct(rowCount,3));field.allocate(points.size());
    // A curved host integrates its map at every probe, so its launches stay one row long.
    const std::size_t launchRows=block.curved?1:64;
    for(std::size_t begin=0;begin<rowCount;begin+=launchRows){const auto here=std::min(launchRows,rowCount-begin);
        flatFieldRows<<<unsigned(here),hexFieldThreads,0,stream>>>(d,cellCount,hosts,poly.data,probe.data,begin,rows.data);
        check(cudaGetLastError());check(cudaStreamSynchronize(stream));++result.launches;}
    hexFieldReduce<<<unsigned((probes+63)/64),64,0,stream>>>(probes,terms,rows.data,field.data,prepared.invalid.data);check(cudaGetLastError());++result.launches;
    if(prepared.flag(stream,result.transfers,"field_flag")!=0)throw std::runtime_error("Device field evaluation produced a nonfinite value");
    result.field.resize(points.size());
    check(cudaMemcpyAsync(result.field.data(),field.data,points.size()*sizeof(double),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));
    record("field_values","D2H",points.size()*sizeof(double));
    return result;
}
} // namespace radia_hmatrix_cuda
