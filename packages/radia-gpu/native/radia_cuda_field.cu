// Experimental retained-cloud RadFld backend. No implicit element sampling,
// CPU evaluation fallback, Python runtime, or TaskManager execution.
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cuda_runtime.h>
#include "rad_hdiv_field_evaluator.h"
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
    DeviceBuffer points,order,nodes,images;
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
            points.p=order.p=nodes.p=images.p=nullptr;
            if(previous!=device) cudaSetDevice(previous);
            return;
        }
        cudaStreamSynchronize(stream);
        for(auto* buffer:{&points,&order,&nodes,&images}) {
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
std::vector<double> rows3(const mxArray* a,const char* name) {
    const auto* p=doubles(a,3,name); const std::size_t n=mxGetM(a);
    std::vector<double> r(bytes(n,3));
    for(std::size_t i=0;i<n;++i) for(int k=0;k<3;++k) r[3*i+k]=p[i+n*k];
    return r;
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

__device__ void base(const FieldCloudPoint* points,std::size_t npoints,
        const std::size_t* order,const FieldCloudNode* nodes,double theta,
        bool tree,const double r[3],double out[3]) {
    out[0]=out[1]=out[2]=0;
    if(!tree) {
        FieldCompensated3 sum;
        for(std::size_t i=0;i<npoints;++i) {
            double v[3]{}; AddCloudPoint(points[i].position,points[i].strength,r,v); sum.Add(v);
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
            AddFieldMultipole(node.center,node.charge,node.dipole,node.second,r,out);
            index=node.escape;
        } else if(node.left<0) {
            for(std::size_t p=node.begin;p<node.end;++p) {
                const auto& source=points[order[p]];
                AddCloudPoint(source.position,source.strength,r,out);
            }
            index=node.escape;
        } else index=node.left;
    }
}
__global__ void evaluate(const FieldCloudPoint* points,std::size_t npoints,
        const std::size_t* order,const FieldCloudNode* nodes,double theta,
        const FieldCloudImage* images,std::size_t nimages,const double* obs,
        std::size_t nobs,bool tree,double* result) {
    const std::size_t i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(i>=nobs) return;
    const double* r=obs+3*i;
    FieldCompensated3 sum;
    double v[3]; base(points,npoints,order,nodes,theta,tree,r,v); sum.Add(v);
    for(std::size_t j=0;j<nimages;++j) {
        const auto& image=images[j]; double t[3]={r[0],r[1],r[2]};
        for(int k=0;k<3;++k) if(image.mask&(1<<k)) t[k]=-t[k];
        if(image.angle!=0) {
            const double c=cos(-image.angle),s=sin(-image.angle),x=t[0],y=t[1];
            t[0]=c*x-s*y; t[1]=s*x+c*y;
        }
        base(points,npoints,order,nodes,theta,tree,t,v);
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
    if(command=="create_cloud") {
        if(nrhs!=8||nlhs!=2) throw std::invalid_argument("create_cloud requires xyz,q,masks,signs,angles,leaf,theta and two outputs");
        const std::size_t n=mxGetM(in[1]);
        if(!n || n>std::size_t(std::numeric_limits<int>::max()/2)) throw std::invalid_argument("Source count outside tree index range");
        auto xyz=rows3(in[1],"xyz"); const auto* q=doubles(in[2],1,"strength");
        if(mxGetM(in[2])!=n) throw std::invalid_argument("strength length mismatch");
        const auto* masks=doubles(in[3],1,"masks"); const auto* signs=doubles(in[4],1,"signs");
        const auto* angles=doubles(in[5],1,"angles"); const std::size_t ni=mxGetM(in[3]);
        if(mxGetM(in[4])!=ni||mxGetM(in[5])!=ni) throw std::invalid_argument("Image array length mismatch");
        std::vector<int> im(ni); std::vector<double> si(ni),an(ni);
        for(std::size_t i=0;i<ni;++i) {
            if(masks[i]!=std::floor(masks[i])||masks[i]<0||masks[i]>7) throw std::invalid_argument("Image masks must be integers in [0,7]");
            im[i]=int(masks[i]); si[i]=signs[i]; an[i]=angles[i];
        }
        const double leaf=scalar(in[6]);
        if(leaf<1||leaf>65536||leaf!=std::floor(leaf)) throw std::invalid_argument("Invalid leaf size");
        FieldEvaluatorOptions options; options.leaf_size=int(leaf); options.theta=scalar(in[7]);
        room(bytes(n,sizeof(FieldCloudPoint)+sizeof(std::size_t)+2*sizeof(FieldCloudNode))+
             bytes(ni,sizeof(FieldCloudImage)));
        auto host=HDivFieldEvaluator::FromCloud(std::move(xyz),std::vector<double>(q,q+n),im,si,options);
        host->SetImageRotations(an); auto packed=host->ExportCloudSnapshot();
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
            bytes(packed.images.size(),sizeof(FieldCloudImage));
        room(need); std::size_t before,total,after; check(cudaMemGetInfo(&before,&total));
        upload(*c,c->points,packed.points); upload(*c,c->order,packed.order);
        upload(*c,c->nodes,packed.nodes); upload(*c,c->images,packed.images);
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
        if(algorithm!="direct"&&algorithm!="tree") throw std::invalid_argument("Explicit direct or tree algorithm required");
        const std::size_t n=mxGetM(in[2]);
        if(n>std::size_t(std::numeric_limits<int>::max())) throw std::invalid_argument("Observation count exceeds launch range");
        auto obs=rows3(in[2],"observations"); const std::size_t one=bytes(n,3*sizeof(double)), temporary=bytes(one,2);
        room(temporary); std::size_t before,total,after; check(cudaMemGetInfo(&before,&total));
        DeviceBuffer d_obs,d_out; d_obs.allocate(one); d_out.allocate(one);
        std::vector<double> values(obs.size());
        if(n) {
            check(cudaMemcpyAsync(d_obs.p,obs.data(),one,cudaMemcpyHostToDevice,c.stream));
            evaluate<<<unsigned((n+127)/128),128,0,c.stream>>>(c.points.as<FieldCloudPoint>(),c.npoints,
                c.order.as<std::size_t>(),c.nodes.as<FieldCloudNode>(),c.theta,c.images.as<FieldCloudImage>(),
                c.nimages,d_obs.as<double>(),n,algorithm=="tree",static_cast<double*>(d_out.p));
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
void mexFunction(int nlhs,mxArray** out,int nrhs,const mxArray** in) {
    try { mxInitGPU(); static bool registered=false; if(!registered) { mexAtExit(cleanup); registered=true; }
        dispatch(nlhs,out,nrhs,in);
    } catch(const std::exception& e) { mexErrMsgIdAndTxt("radia:cuda:field", "%s",e.what()); }
}
