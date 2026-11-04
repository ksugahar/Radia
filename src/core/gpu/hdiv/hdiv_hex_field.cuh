// Device magnetization field of a resident HEX solution, with its mirror and cyclic images: the CPU direct route
// (HDivFieldEvaluator). All affine hosts take the exact 6-tet and 2-triangle polynomial fields of their physical
// charge densities; otherwise every host takes the shared tensor point cloud. E_T(x) = T E(T^-1 x), added with the
// image sign. Only probes, the state vector and geometry cross H2D and the requested field D2H; the charges are the
// resident charge map's product on the device.
#pragma once
#include "hdiv_flat_gram.cuh"
namespace radia_hmatrix_cuda {
constexpr int hexFieldThreads=64;
// CPU GaussLegendreTable of the far field rules; the unused row zero and columns stay zero.
inline void uploadHexFieldFarGauss(){
    rad_hdiv_cuda::GaussLegendreTable table{};
    for(int n=1;n<=rad_hdiv_cuda::FAR_MAX_NODES;++n)for(int i=0;i<n;++i){
        double t=std::cos(3.14159265358979323846*(i+0.75)/(n+0.5)),derivative=1.0;
        for(int iteration=0;iteration<100;++iteration){double p0=1.0,p1=t;
            for(int k=2;k<=n;++k){const double p2=((2.0*k-1.0)*t*p1-(k-1.0)*p0)/k;p0=p1;p1=p2;}
            derivative=n*(t*p1-p0)/(t*t-1.0);const double step=p1/derivative;t-=step;if(std::fabs(step)<1e-16)break;}
        table.x[n][i]=0.5*(1.0-t);table.w[n][i]=1.0/((1.0-t*t)*derivative*derivative);}
    check(cudaMemcpyToSymbol(rad_hdiv_cuda::DeviceFarGauss,&table,sizeof(table)));
}
// Physical density of one affine host about its centre, the CPU per-host polynomial.
__global__ void hexFieldPolynomials(HexDevice d,const double* charges,double* poly,int* invalid){
    const int g=int(blockIdx.x)*blockDim.x+threadIdx.x;if(g>=d.hosts)return;
    const bool cell=hexCell(d,g);const int ncoeff=cell?84:35;double lin[3][4],scaleInv=0.0,centre[3];double* out=poly+std::size_t(g)*84;
    if(!hexAffineForms(hexNodes(d,g),cell,lin,scaleInv)){atomicExch(invalid,1);return;}
    hexCenterForms(hexNodes(d,g),cell,lin,centre);
    for(int i=0;i<ncoeff;++i)out[i]=0.0;
    for(int l=0;l<d.count[g];++l){double mode[84];
        if(!hexMonomialCoeff(lin,scaleInv,cell,hexExponent(d,g,l),mode)){atomicExch(invalid,1);return;}
        const double q=charges[d.first[g]+l];for(int i=0;i<ncoeff;++i)out[i]+=q*mode[i];}
}
// CPU FromPolynomialTet (order 1: cubic and quadratic closed forms) and FromSexticQuarticPolynomialTet (order 2).
__device__ inline void hexFieldAffineHost(const HexDevice& d,int g,const double* poly,const double* r,double* out){
    const double* p=poly+std::size_t(g)*84;const bool cell=hexCell(d,g);const int ncoeff=cell?84:35;
    bool empty=true;for(int i=0;i<ncoeff&&empty;++i)empty=p[i]==0.0;if(empty)return;
    const double* nd=hexNodes(d,g);const double* c=d.affineCenter+3*std::size_t(g);const double lr[3]={r[0]-c[0],r[1]-c[1],r[2]-c[2]};
    if(cell){
        const int corner[8]={0,2,8,6,18,20,26,24},tets[6][4]={{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}};
        for(int sub=0;sub<6;++sub){double V[4][3];for(int i=0;i<4;++i)for(int k=0;k<3;++k)V[i][k]=nd[3*corner[tets[sub][i]]+k]-c[k];
            if(d.order==1){double value[3]={0.0,0.0,0.0};rad_hdiv_cuda::TetVolFieldCubic(V,lr,p,value);for(int k=0;k<3;++k)out[k]+=value[k];}
            else{double basis[84][3];rad_hdiv_cuda::TetVolFieldBasisUpTo6(V,lr,basis);
                for(int i=0;i<84;++i)for(int k=0;k<3;++k)out[k]+=p[i]*basis[i][k];}}
        return;
    }
    const int corner[4]={0,2,8,6},tris[2][3]={{0,1,2},{0,2,3}};
    for(int sub=0;sub<2;++sub){double V[3][3];for(int i=0;i<3;++i)for(int k=0;k<3;++k)V[i][k]=nd[3*corner[tris[sub][i]]+k]-c[k];
        if(d.order==1){double slope[3],hessian[3][3],value[3]={0.0,0.0,0.0};
            for(int a=0;a<3;++a){int e[3]={0,0,0};e[a]=1;slope[a]=p[hexPolyIdx(e[0],e[1],e[2])];}
            for(int a=0;a<3;++a)for(int b=0;b<3;++b){int e[3]={0,0,0};++e[a];++e[b];const double v=p[hexPolyIdx(e[0],e[1],e[2])];hessian[a][b]=a==b?v:0.5*v;}
            rad_hdiv_cuda::QuadTriField(V,lr,p[0],slope,hessian,value);for(int k=0;k<3;++k)out[k]+=value[k];}
        else{double basis[35][3];rad_hdiv_cuda::TriFieldBasisUpTo4(V,lr,basis);
            for(int i=0;i<35;++i)for(int k=0;k<3;++k)out[k]+=p[i]*basis[i][k];}}
}
// CPU shared tensor cloud: co-located charge modes combine at each point, the rule the host far table holds.
__device__ inline void hexFieldCloudHost(const HexDevice& d,int g,const double* charges,const double* r,double* out){
    const int n1=d.out.n,np=hexCell(d,g)?n1*n1*n1:n1*n1;const double* X=d.farX+3*std::size_t(g)*d.farStride;const double* W=d.farW+std::size_t(g)*d.farStride;
    for(int q=0;q<np;++q){const double xi[3]={d.out.x[q%n1],d.out.x[(q/n1)%n1],hexCell(d,g)?d.out.x[q/(n1*n1)]:0.0};
        double density=0.0;for(int l=0;l<d.count[g];++l)density+=charges[d.first[g]+l]*hexDensity(d,g,hexExponent(d,g,l),xi);
        rad_hdiv::AddCloudPoint(X+3*q,W[q]*density,r,out);}
}
__device__ inline double hexFieldBlockSum(double x,double* shared){
    shared[threadIdx.x]=x;__syncthreads();
    for(int s=hexFieldThreads/2;s>0;s>>=1){if(int(threadIdx.x)<s)shared[threadIdx.x]+=shared[threadIdx.x+s];__syncthreads();}
    const double r=shared[0];__syncthreads();return r;
}
// One block per (probe, term); its threads share the hosts, then reduce.
__global__ void hexFieldRows(HexDevice d,const double* charges,const double* poly,bool affine,const double* probes,std::size_t firstRow,double* rows){
    __shared__ double shared[hexFieldThreads];
    const auto row=firstRow+std::size_t(blockIdx.x);const std::size_t terms=1+std::size_t(d.imageCount),probe=row/terms;const int img=int(row%terms);
    double rm[3];hexImagePoint(d,img,probes+3*probe,rm);double sum[3]={0.0,0.0,0.0};
    for(int g=int(threadIdx.x);g<d.hosts;g+=hexFieldThreads){
        if(affine)hexFieldAffineHost(d,g,poly,rm,sum);else hexFieldCloudHost(d,g,charges,rm,sum);}
    for(int k=0;k<3;++k)sum[k]=hexFieldBlockSum(sum[k],shared);
    if(threadIdx.x==0){double image[3];hexImageVector(d,img,sum,image);const double sign=img?d.images[img-1].sign:1.0;
        for(int k=0;k<3;++k)rows[3*row+k]=sign*image[k];}
}
// Scales by 1/(4 pi) like the CPU field_from_solution.
__global__ void hexFieldReduce(std::size_t probes,std::size_t terms,const double* rows,double* field,int* invalid){
    const auto p=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(p>=probes)return;
    for(int k=0;k<3;++k){double sum=0.0;for(std::size_t t=0;t<terms;++t)sum+=rows[3*(p*terms+t)+k];
        const double value=sum/(4.0*3.14159265358979323846);if(!isfinite(value))atomicExch(invalid,1);field[3*p+k]=value;}
}
struct HexFieldResult {HostVector<double> field;HostVector<FlatTransfer> transfers{HostFaultAllocator<FlatTransfer>("hex.field-transfers")};std::size_t launches=0,ownedBytes=0;bool affine=false;};
// x is the prepared FE state; the Resident supplies q = B x on its own stream and device.
inline HexFieldResult hexDeviceField(const FlatBlock& block,Resident& resident,const HostVector<double>& x,const HostVector<double>& points){
    if(!block.hex||!block.hostTypes.empty())throw std::invalid_argument("Device field evaluation has no route for this geometry family");
    if(points.empty()||points.size()%3||x.size()!=resident.fieldStateSize())throw std::invalid_argument("Device field evaluation requires (points,3) probes and the prepared state");
    for(double v:points)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite field probe");
    for(double v:x)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite field state");
    const std::size_t probes=points.size()/3,terms=1+block.images.size(),rowCount=checkedProduct(probes,terms),n=block.charges();
    if(rowCount>std::size_t(INT_MAX))throw std::invalid_argument("Device field probes exceed the launch range");
    HexFieldResult result;cudaStream_t stream=resident.executionStream();
    auto record=[&](const char* payload,const char* direction,std::size_t bytes){result.transfers.push_back({payload,direction,bytes});};
    uploadHexFieldFarGauss();record("field_far_gauss","H2D",sizeof(rad_hdiv_cuda::GaussLegendreTable));
    HexPrepared prepared;int code=0;prepared.prepareHosts(block,stream,record,result.launches,code);
    const HexDevice d=prepared.d;const std::size_t hosts=block.cellCount()+block.faceCount();
    Buffer<double> state,charges,probe,poly,rows,field;
    state.upload(x);record("field_state","H2D",checkedProduct(x.size(),sizeof(double)));
    probe.upload(points);record("field_probes","H2D",checkedProduct(points.size(),sizeof(double)));
    charges.allocate(n);resident.chargeProducts(state.data,charges.data);
    // CPU route choice: one non-affine host sends every host to the cloud.
    {HostVector<unsigned char> flags(hosts,0,HostFaultAllocator<unsigned char>("hex.field-affine"));
        check(cudaMemcpyAsync(flags.data(),d.affine,hosts,cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));record("field_host_affine","D2H",hosts);
        result.affine=true;for(auto flag:flags)result.affine=result.affine&&flag!=0;}
    poly.allocate(result.affine?checkedProduct(hosts,84):1);
    if(result.affine){hexFieldPolynomials<<<unsigned((hosts+127)/128),128,0,stream>>>(d,charges.data,poly.data,prepared.invalid.data);check(cudaGetLastError());++result.launches;
        if(prepared.flag(stream,record,"field_polynomial_flag")!=0)throw std::runtime_error("Device field polynomial is nonfinite or its host stopped being affine");}
    rows.allocate(checkedProduct(rowCount,3));field.allocate(points.size());
    constexpr std::size_t launchRows=64;
    for(std::size_t first=0;first<rowCount;first+=launchRows){const auto here=std::min(launchRows,rowCount-first);
        hexFieldRows<<<unsigned(here),hexFieldThreads,0,stream>>>(d,charges.data,poly.data,result.affine,probe.data,first,rows.data);
        check(cudaGetLastError());check(cudaStreamSynchronize(stream));++result.launches;}
    hexFieldReduce<<<unsigned((probes+63)/64),64,0,stream>>>(probes,terms,rows.data,field.data,prepared.invalid.data);check(cudaGetLastError());++result.launches;
    if(prepared.flag(stream,record,"field_flag")!=0)throw std::runtime_error("Device field evaluation produced a nonfinite value");
    result.field.resize(points.size());
    check(cudaMemcpyAsync(result.field.data(),field.data,points.size()*sizeof(double),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));
    record("field_values","D2H",points.size()*sizeof(double));
    result.ownedBytes=checkedProduct(rowCount,3*sizeof(double))+checkedProduct(n+x.size()+2*points.size(),sizeof(double))+(result.affine?checkedProduct(hosts,84*sizeof(double)):0);
    return result;
}
} // namespace radia_hmatrix_cuda
