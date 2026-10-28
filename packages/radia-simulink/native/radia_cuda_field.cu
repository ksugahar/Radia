#include "mex_error_boundary.hpp"
// Experimental retained-cloud RadFld backend. No implicit element sampling,
// CPU evaluation fallback, Python runtime, or TaskManager execution.
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cuda_runtime.h>
#include "rad_hdiv_field_evaluator.h"
#include "../../../src/core/gpu/hdiv/radia_field_exact.cuh"
#include <algorithm>
#include <cstdint>
#include <chrono>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

using namespace rad_hdiv;
namespace {
void check(cudaError_t e) { if(e!=cudaSuccess) throw std::runtime_error(cudaGetErrorString(e)); }
std::size_t bytes(std::size_t n,std::size_t width) {
    if(n>std::numeric_limits<std::size_t>::max()/width) throw std::overflow_error("Allocation size overflow");
    return n*width;
}
struct DeviceBuffer {
    void* p=nullptr; std::size_t size=0;
    ~DeviceBuffer() { if(p) cudaFree(p); }
    DeviceBuffer()=default;
    DeviceBuffer(const DeviceBuffer&)=delete;
    void allocate(std::size_t n) { if(n) check(cudaMalloc(&p,n)); size=n; }
    template<class T> const T* as() const { return static_cast<const T*>(p); }
};
struct Context {
    DeviceBuffer points,order,nodes,images,sources,tets,triangles,cubic_triangles,curved_tets,curved_triangles,gauss_points,gauss_weights;
    FieldSourceView view{};
    cudaStream_t stream=nullptr;
    std::size_t npoints=0, nimages=0, owned=0;
    double theta=0;
    int device=0;
    ~Context() {
        if(!stream) return;
        int previous=device;
        cudaGetDevice(&previous);
        if(previous!=device) cudaSetDevice(device);
        const auto state=cudaStreamQuery(stream);
        if(state!=cudaSuccess && state!=cudaErrorNotReady) {
            // A device reset already released these allocations. Never free
            // stale addresses that another CUDA consumer may now own.
            for(auto* buffer:{&points,&order,&nodes,&images,&sources,&tets,&triangles,&cubic_triangles,&curved_tets,&curved_triangles,&gauss_points,&gauss_weights}) buffer->p=nullptr;
            if(previous!=device) cudaSetDevice(previous);
            return;
        }
        cudaStreamSynchronize(stream);
        for(auto* buffer:{&points,&order,&nodes,&images,&sources,&tets,&triangles,&cubic_triangles,&curved_tets,&curved_triangles,&gauss_points,&gauss_weights}) {
            if(buffer->p) cudaFree(buffer->p);
            buffer->p=nullptr;
        }
        cudaStreamDestroy(stream);
        if(previous!=device) cudaSetDevice(previous);
    }
};
std::map<std::uint64_t,std::unique_ptr<Context>> registry;
std::uint64_t last_generation=0;
void cleanup() { registry.clear(); }
void room(std::size_t need) {
    std::size_t available,total; check(cudaMemGetInfo(&available,&total));
    const std::size_t reserve=std::max<std::size_t>(2ULL<<30,total/5);
    if(available<=reserve || need>available-reserve)
        throw std::runtime_error("Field allocation would exceed GPU budget (reserve max(2 GiB,20%))");
}
const double* doubles(const mxArray* a,std::size_t cols,const char* name) {
    if(!mxIsDouble(a)||mxIsComplex(a)||mxIsSparse(a)||mxGetNumberOfDimensions(a)!=2||mxGetN(a)!=cols)
        throw std::invalid_argument(std::string(name)+" must be a full real double matrix with the declared column count");
    const double* p=mxGetPr(a);
    for(std::size_t i=0;i<mxGetNumberOfElements(a);++i)
        if(!std::isfinite(p[i])) throw std::invalid_argument(std::string(name)+" must be finite");
    return p;
}
std::vector<double> rows(const mxArray* a,std::size_t cols,const char* name) {
    const auto* p=doubles(a,cols,name); const std::size_t n=mxGetM(a);
    std::vector<double> r(bytes(n,cols));
    for(std::size_t i=0;i<n;++i) for(std::size_t k=0;k<cols;++k) r[cols*i+k]=p[i+n*k];
    return r;
}
std::vector<double> rows3(const mxArray* a,const char* name) { return rows(a,3,name); }
const mxArray* member(const mxArray* a,const char* name) {
    const auto* value=mxGetField(a,0,name);
    if(!value) throw std::invalid_argument(std::string("Source specification missing ")+name);
    return value;
}
std::vector<double> optionalRows(const mxArray* a,const char* name,std::size_t cols) {
    const auto* value=mxGetField(a,0,name);
    return value?rows(value,cols,name):std::vector<double>{};
}

double scalar(const mxArray* a) {
    if(!mxIsDouble(a)||mxIsComplex(a)||mxIsSparse(a)||mxGetNumberOfElements(a)!=1||!std::isfinite(mxGetScalar(a)))
        throw std::invalid_argument("Expected a finite real double scalar");
    return mxGetScalar(a);
}
std::string text(const mxArray* a) {
    if(!mxIsChar(a)||mxGetM(a)!=1) throw std::invalid_argument("Expected a character row");
    char* value=mxArrayToString(a); if(!value) throw std::bad_alloc();
    std::string s(value); mxFree(value); return s;
}
std::uint64_t token(const mxArray* a) {
    if(!mxIsUint64(a)||mxIsComplex(a)||mxGetNumberOfElements(a)!=1)
        throw std::invalid_argument("Handle must be a scalar uint64");
    return *static_cast<const std::uint64_t*>(mxGetData(a));
}
Context& context(const mxArray* a) {
    auto found=registry.find(token(a));
    if(found==registry.end()) throw std::invalid_argument("Invalid, stale, or foreign field handle");
    int device; check(cudaGetDevice(&device));
    if(device!=found->second->device) throw std::runtime_error("GPU device changed while field handle is live");
    check(cudaStreamQuery(found->second->stream));
    return *found->second;
}
template<class T> void upload(Context& c,DeviceBuffer& dst,const std::vector<T>& src) {
    dst.allocate(bytes(src.size(),sizeof(T)));
    if(dst.size) check(cudaMemcpyAsync(dst.p,src.data(),dst.size,cudaMemcpyHostToDevice,c.stream));
    c.owned+=dst.size;
}

// Explicit mixed policy: far-node arithmetic FP32, geometry/opening and
// exact near-element kernels FP64. Results are accumulated into FP64 output.
__device__ void multipoleMixed32(const FieldCloudNode& node,const double r[3],double out[3]) {
    float R[3],p[3],S[9];bool safe=true;
    for(int k=0;k<3;++k){R[k]=float(r[k]-node.center[k]);p[k]=float(node.dipole[k]);}
    for(int k=0;k<9;++k)S[k]=float(node.second[k]);
    float q=float(node.charge);
    for(int k=0;k<3;++k)safe=safe&&isfinite(R[k])&&isfinite(p[k])&&(node.dipole[k]==0||p[k]!=0);
    for(int k=0;k<9;++k)safe=safe&&isfinite(S[k])&&(node.second[k]==0||S[k]!=0);
    safe=safe&&isfinite(q)&&(node.charge==0||q!=0);
    const float r2=R[0]*R[0]+R[1]*R[1]+R[2]*R[2];
    safe=safe&&isfinite(r2)&&r2>1e-20f&&r2<1e20f;
    float value[3]{};
    if(safe){
        float inv=1.0f/sqrtf(r2),i3=inv/r2,i5=i3/r2,i7=i5/r2;
        float pr=p[0]*R[0]+p[1]*R[1]+p[2]*R[2],sr[3]{},trace=S[0]+S[4]+S[8];
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)sr[a]+=S[3*a+b]*R[b];
        float rsr=R[0]*sr[0]+R[1]*sr[1]+R[2]*sr[2];
        for(int a=0;a<3;++a){
            value[a]=q*R[a]*i3+(3.0f*R[a]*pr*i5-p[a]*i3)
                +(7.5f*R[a]*rsr*i7-1.5f*(trace*R[a]+2.0f*sr[a])*i5);
            safe=safe&&isfinite(value[a]);
        }
    }
    if(safe)for(int k=0;k<3;++k)out[k]+=double(value[k]);
    else AddFieldMultipole(node.center,node.charge,node.dipole,node.second,r,out);
}
__device__ void base(FieldSourceView view,const FieldSourceRef* sources,std::size_t npoints,
        const std::size_t* order,const FieldCloudNode* nodes,double theta,
        bool tree,bool mixed,const double r[3],double out[3]) {
    out[0]=out[1]=out[2]=0;
    if(!tree) {
        FieldCompensated3 sum;
        for(std::size_t i=0;i<npoints;++i) {
            double v[3]{}; rad_hdiv_cuda::AddFieldSource(view,sources[i],r,v); sum.Add(v);
        }
        sum.Store(out); return;
    }
    // The CPU tree is threaded on the host: no fixed device traversal stack.
    int index=0;
    while(index>=0) {
        const auto& node=nodes[index];
        const double x=r[0]-node.center[0], y=r[1]-node.center[1], z=r[2]-node.center[2];
        const double distance=sqrt(x*x+y*y+z*z);
        if(node.left>=0 && distance>0 && node.radius/distance<=theta) {
            if(mixed) multipoleMixed32(node,r,out);
            else AddFieldMultipole(node.center,node.charge,node.dipole,node.second,r,out);
            index=node.escape;
        } else if(node.left<0) {
            for(std::size_t p=node.begin;p<node.end;++p) {
                rad_hdiv_cuda::AddFieldSource(view,sources[order[p]],r,out);
            }
            index=node.escape;
        } else index=node.left;
    }
}
__global__ void evaluate(FieldSourceView view,const FieldSourceRef* sources,std::size_t npoints,
        const std::size_t* order,const FieldCloudNode* nodes,double theta,
        const FieldCloudImage* images,std::size_t nimages,const double* obs,
        std::size_t nobs,bool tree,bool mixed,double* result) {
    const std::size_t i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(i>=nobs) return;
    const double* r=obs+3*i;
    FieldCompensated3 sum;
    double v[3]; base(view,sources,npoints,order,nodes,theta,tree,mixed,r,v); sum.Add(v);
    for(std::size_t j=0;j<nimages;++j) {
        const auto& image=images[j]; double t[3]={r[0],r[1],r[2]};
        for(int k=0;k<3;++k) if(image.mask&(1<<k)) t[k]=-t[k];
        if(image.angle!=0) {
            const double c=cos(-image.angle),s=sin(-image.angle),x=t[0],y=t[1];
            t[0]=c*x-s*y; t[1]=s*x+c*y;
        }
        base(view,sources,npoints,order,nodes,theta,tree,mixed,t,v);
        if(image.angle!=0) {
            const double c=cos(image.angle),s=sin(image.angle),x=v[0],y=v[1];
            v[0]=c*x-s*y; v[1]=s*x+c*y;
        }
        for(int k=0;k<3;++k) { if(image.mask&(1<<k)) v[k]=-v[k]; sum.Add(k,image.sign*v[k]); }
    }
    sum.Store(result+3*i);
}
mxArray* metrics(Context& c,std::size_t temporary,std::size_t before,std::size_t after) {
    const char* fields[]={"owned_peak_bytes","persistent_owned_bytes","temporary_owned_bytes",
        "cuda_free_before_bytes","cuda_free_after_bytes","source_count","tree_nodes","image_count"};
    mxArray* a=mxCreateStructMatrix(1,1,8,fields);
    const double values[]={double(c.owned+temporary),double(c.owned),double(temporary),
        double(before),double(after),double(c.npoints),double(c.nodes.size/sizeof(FieldCloudNode)),double(c.nimages)};
    for(int i=0;i<8;++i) mxSetFieldByNumber(a,0,i,mxCreateDoubleScalar(values[i]));
    return a;
}
void dispatch(int nlhs,mxArray** out,int nrhs,const mxArray** in) {
    if(nrhs<1) throw std::invalid_argument("Command required");
    const auto command=text(in[0]);
    if(command=="live_count") {
        if(nrhs!=1||nlhs!=1) throw std::invalid_argument("live_count requires one output");
        out[0]=mxCreateDoubleScalar(double(registry.size())); return;
    }
    if(command=="destroy") {
        if(nrhs!=2||nlhs!=0) throw std::invalid_argument("destroy(handle) requires no output");
        auto found=registry.find(token(in[1]));
        if(found==registry.end()) throw std::invalid_argument("Invalid, stale, or foreign field handle");
        // Destruction remains possible after GPU reset or device selection.
        registry.erase(found); mexUnlock(); return;
    }
    if(command=="create_cloud" || command=="create_elements") {
        const bool elements=command=="create_elements";
        const int shift=elements?1:0;
        if(nrhs!=8-shift||nlhs!=2) throw std::invalid_argument("create_cloud requires xyz,q,masks,signs,angles,leaf,theta and two outputs");
        std::size_t n=0; std::vector<double> xyz,qv;
        if(!elements) {
            n=mxGetM(in[1]);
            if(!n || n>std::size_t(std::numeric_limits<int>::max()/2)) throw std::invalid_argument("Source count outside tree index range");
            room(bytes(n,sizeof(FieldCloudPoint)+sizeof(FieldSourceRef)+sizeof(std::size_t)+2*sizeof(FieldCloudNode)));
            xyz=rows3(in[1],"xyz");
            const auto* q=doubles(in[2],1,"strength");
            if(mxGetM(in[2])!=n) throw std::invalid_argument("strength length mismatch");
            if(n) qv.assign(q,q+n);
        } else if(!mxIsStruct(in[1]) || mxGetNumberOfElements(in[1])!=1) {
            throw std::invalid_argument("Element source must be a scalar versioned struct");
        }
        const auto* masks=doubles(in[3-shift],1,"masks"); const auto* signs=doubles(in[4-shift],1,"signs");
        const auto* angles=doubles(in[5-shift],1,"angles"); const std::size_t ni=mxGetM(in[3-shift]);
        if(mxGetM(in[4-shift])!=ni||mxGetM(in[5-shift])!=ni) throw std::invalid_argument("Image array length mismatch");
        std::vector<int> im(ni); std::vector<double> si(ni),an(ni);
        for(std::size_t i=0;i<ni;++i) {
            if(masks[i]!=std::floor(masks[i])||masks[i]<0||masks[i]>7) throw std::invalid_argument("Image masks must be integers in [0,7]");
            im[i]=int(masks[i]); si[i]=signs[i]; an[i]=angles[i];
        }
        const double leaf=scalar(in[6-shift]);
        if(leaf<1||leaf>65536||leaf!=std::floor(leaf)) throw std::invalid_argument("Invalid leaf size");
        FieldEvaluatorOptions options; options.leaf_size=int(leaf); options.theta=scalar(in[7-shift]);
        std::shared_ptr<HDivFieldEvaluator> host;
        if(!elements) {
            if(!n || n>std::size_t(std::numeric_limits<int>::max()/2)) throw std::invalid_argument("Source count outside tree index range");
            room(bytes(n,sizeof(FieldCloudPoint)+sizeof(FieldSourceRef)+sizeof(std::size_t)+2*sizeof(FieldCloudNode))+bytes(ni,sizeof(FieldCloudImage)));
            host=HDivFieldEvaluator::FromCloud(std::move(xyz),std::move(qv),im,si,options);
        } else {
            const auto* spec=in[1];
            if(scalar(member(spec,"version"))!=1) throw std::invalid_argument("Unsupported element source schema version");
            const auto kind=text(member(spec,"representation"));
            std::size_t vc=0,sc=0;
            if(kind=="tet") {vc=16;sc=22;}
            else if(kind=="polynomial_tet") {vc=32;sc=22;}
            else if(kind=="cubic_tet") {vc=32;sc=29;}
            else if(kind=="sextic_quartic_tet") {vc=96;sc=44;}
            else if(kind=="curved_tet") {vc=34;sc=24;}
            else throw std::invalid_argument("Unsupported element representation");
            const auto* volume=member(spec,"volume");const auto* surface=member(spec,"surface");
            n=mxGetM(volume)+mxGetM(surface);
            if(!n || n>std::size_t(std::numeric_limits<int>::max()/2)) throw std::invalid_argument("Source count outside tree index range");
            // Conservative admission precedes host tree/moment construction.
            room(bytes(n,sizeof(FieldTetSource)+sizeof(FieldSourceRef)+sizeof(std::size_t)+2*sizeof(FieldCloudNode))+bytes(ni,sizeof(FieldCloudImage)));
            auto v=rows(volume,vc,"volume"), t=rows(surface,sc,"surface");
            auto vo=optionalRows(spec,"volume_origins",3),so=optionalRows(spec,"surface_origins",3);
            if(kind=="curved_tet") {
                if(!vo.empty()||!so.empty()) throw std::invalid_argument("Curved sources use reference densities; origins are unsupported");
                auto gl=rows(member(spec,"gauss_points"),1,"gauss_points");
                auto gw=rows(member(spec,"gauss_weights"),1,"gauss_weights");
                if(gl.empty()||gl.size()>64||gl.size()!=gw.size()) throw std::invalid_argument("Gauss rule must contain 1..64 matched nodes and weights");
                for(std::size_t i=0;i<gl.size();++i)
                    if(gl[i]<0||gl[i]>1||gw[i]<=0) throw std::invalid_argument("Gauss nodes must lie in [0,1] and weights must be positive");
                host=HDivFieldEvaluator::FromCurvedTet(std::move(v),std::move(t),std::move(gl),std::move(gw),im,si,options);
            } else if(kind=="tet") host=HDivFieldEvaluator::FromTet(std::move(v),std::move(t),im,si,options,std::move(vo),std::move(so));
            else if(kind=="polynomial_tet") host=HDivFieldEvaluator::FromPolynomialTet(std::move(v),std::move(t),im,si,options,std::move(vo),std::move(so));
            else if(kind=="cubic_tet") host=HDivFieldEvaluator::FromCubicPolynomialTet(std::move(v),std::move(t),im,si,options,std::move(vo),std::move(so));
            else host=HDivFieldEvaluator::FromSexticQuarticPolynomialTet(std::move(v),std::move(t),im,si,options,std::move(vo),std::move(so));
        }
        host->SetImageRotations(an); auto packed=host->ExportSnapshot();
        for(const auto& node:packed.nodes) {
            if(!std::isfinite(node.radius)||!std::isfinite(node.charge))
                throw std::invalid_argument("Source extent/strength overflows tree moments");
            for(double x:node.center) if(!std::isfinite(x)) throw std::invalid_argument("Tree center overflow");
            for(double x:node.dipole) if(!std::isfinite(x)) throw std::invalid_argument("Tree dipole overflow");
            for(double x:node.second) if(!std::isfinite(x)) throw std::invalid_argument("Tree quadrupole overflow");
        }
        auto c=std::make_unique<Context>(); check(cudaGetDevice(&c->device));
        check(cudaStreamCreateWithFlags(&c->stream,cudaStreamNonBlocking));
        c->npoints=n; c->nimages=ni; c->theta=packed.theta;
        const std::size_t need=bytes(packed.points.size(),sizeof(FieldCloudPoint))+
            bytes(packed.order.size(),sizeof(std::size_t))+bytes(packed.nodes.size(),sizeof(FieldCloudNode))+
            bytes(packed.images.size(),sizeof(FieldCloudImage))+bytes(packed.sources.size(),sizeof(FieldSourceRef))+
            bytes(packed.tets.size(),sizeof(FieldTetSource))+bytes(packed.triangles.size(),sizeof(FieldTriSource))+
            bytes(packed.cubic_triangles.size(),sizeof(FieldCubicTriSource))+
            bytes(packed.curved_tets.size(),sizeof(FieldCurvedTetSource))+bytes(packed.curved_triangles.size(),sizeof(FieldCurvedTriSource))+
            bytes(packed.gauss_points.size()+packed.gauss_weights.size(),sizeof(double));
        room(need); std::size_t before,total,after; check(cudaMemGetInfo(&before,&total));
        check(UploadFieldFarGauss());
        upload(*c,c->sources,packed.sources);
        upload(*c,c->tets,packed.tets); upload(*c,c->triangles,packed.triangles);
        upload(*c,c->cubic_triangles,packed.cubic_triangles);
        upload(*c,c->curved_tets,packed.curved_tets); upload(*c,c->curved_triangles,packed.curved_triangles);
        upload(*c,c->gauss_points,packed.gauss_points); upload(*c,c->gauss_weights,packed.gauss_weights);
        upload(*c,c->points,packed.points); upload(*c,c->order,packed.order);
        upload(*c,c->nodes,packed.nodes); upload(*c,c->images,packed.images);
        c->view={c->tets.as<FieldTetSource>(),c->triangles.as<FieldTriSource>(),
            c->cubic_triangles.as<FieldCubicTriSource>(),c->curved_tets.as<FieldCurvedTetSource>(),
            c->curved_triangles.as<FieldCurvedTriSource>(),c->points.as<FieldCloudPoint>(),
            c->gauss_points.as<double>(),c->gauss_weights.as<double>(),static_cast<int>(packed.gauss_points.size())};
        check(cudaStreamSynchronize(c->stream)); check(cudaMemGetInfo(&after,&total));
        // Steady-clock generations survive module unload/reload without reusing
        // a counter starting at one. Creation is serialized by MATLAB.
        const auto generation=std::uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now().time_since_epoch()).count()) & ((1ULL<<56)-1);
        if(!generation || generation<=last_generation) throw std::overflow_error("Handle generation clock exhausted");
        const std::uint64_t id=(0xF1ULL<<56)|generation;
        last_generation=generation;
        out[0]=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
        *static_cast<std::uint64_t*>(mxGetData(out[0]))=id;
        out[1]=metrics(*c,0,before,after);
        registry.emplace(id,std::move(c)); mexLock(); return;
    }
    if(command=="field") {
        if(nrhs!=4||nlhs!=2) throw std::invalid_argument("field(handle,observations,algorithm) requires two outputs");
        auto& c=context(in[1]); const auto algorithm=text(in[3]);
        if(algorithm!="direct"&&algorithm!="tree"&&algorithm!="tree_mixed32") throw std::invalid_argument("Explicit direct, tree or tree_mixed32 algorithm required");
        const std::size_t n=mxGetM(in[2]);
        if(n>std::size_t(std::numeric_limits<int>::max())) throw std::invalid_argument("Observation count exceeds launch range");
        const std::size_t one=bytes(n,3*sizeof(double)), temporary=bytes(one,2);
        room(temporary);
        auto obs=rows3(in[2],"observations"); std::size_t before,total,after; check(cudaMemGetInfo(&before,&total));
        DeviceBuffer d_obs,d_out; d_obs.allocate(one); d_out.allocate(one);
        std::vector<double> values(obs.size());
        if(n) {
            check(cudaMemcpyAsync(d_obs.p,obs.data(),one,cudaMemcpyHostToDevice,c.stream));
            evaluate<<<unsigned((n+127)/128),128,0,c.stream>>>(c.view,c.sources.as<FieldSourceRef>(),c.npoints,
                c.order.as<std::size_t>(),c.nodes.as<FieldCloudNode>(),c.theta,c.images.as<FieldCloudImage>(),
                c.nimages,d_obs.as<double>(),n,algorithm!="direct",algorithm=="tree_mixed32",static_cast<double*>(d_out.p));
            check(cudaGetLastError());
            check(cudaMemcpyAsync(values.data(),d_out.p,one,cudaMemcpyDeviceToHost,c.stream));
            check(cudaStreamSynchronize(c.stream));
        }
        check(cudaMemGetInfo(&after,&total));
        for(double x:values) if(!std::isfinite(x)) throw std::runtime_error("Nonfinite CUDA field result");
        out[0]=mxCreateDoubleMatrix(n,3,mxREAL); auto* p=mxGetPr(out[0]);
        for(std::size_t i=0;i<n;++i) for(int k=0;k<3;++k) p[i+n*k]=values[3*i+k];
        out[1]=metrics(c,temporary,before,after); return;
    }
    throw std::invalid_argument("Unknown field command");
}
} // namespace
void mexFunction(int nl,mxArray** out,int nr,const mxArray** in){

    char errorMessage[2048]={0};
    if(!radia_mex::checked(dispatch,nl,out,nr,in,errorMessage,sizeof(errorMessage))){
        mexErrMsgIdAndTxt("radia:cuda:field","%s",errorMessage);
    }
}
