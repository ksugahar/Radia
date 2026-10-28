/* SPDX-License-Identifier: MPL-2.0 */
// Included inside the gateway namespace, after validated CSR/AMS definitions.
struct OwnedArray {
    mxArray* value=nullptr;
    explicit OwnedArray(const mxArray* source){value=mxDuplicateArray(source);
        if(!value)throw std::bad_alloc();mexMakeArrayPersistent(value);}
    ~OwnedArray(){if(value)mxDestroyArray(value);}
    OwnedArray(const OwnedArray&)=delete;
};
struct ResidentHierarchy {
    OwnedArray owner;Input input;AMS ams;
    explicit ResidentHierarchy(const mxArray* source):owner(source),ams(device_node(owner.value,input),input.validation){input.validation.check();}
};
struct ResidentAMS {
    std::unique_ptr<ResidentHierarchy> hierarchy;
    ResidentWork work;Flag validation;std::unique_ptr<AMSGraph> graph;
    int device;unsigned solves=0,rebuilds=1;
    explicit ResidentAMS(std::unique_ptr<ResidentHierarchy> h):hierarchy(std::move(h)),work(hierarchy->ams.a.n){ck(cudaGetDevice(&device));}
    void checkDevice(){int current;ck(cudaGetDevice(&current));if(current!=device)throw std::runtime_error("resident context belongs to another CUDA device");}
    void replace(const mxArray* source){checkDevice();auto next=std::make_unique<ResidentHierarchy>(source);
        if(next->ams.a.n!=hierarchy->ams.a.n)throw std::runtime_error("resident update cannot change dimensions");
        graph.reset();hierarchy.swap(next);++rebuilds;}
};
using Handle=uint64_t;
std::unordered_map<Handle,std::unique_ptr<ResidentAMS>> amsContexts;
constexpr Handle AMS_TAG=0xA65A;
Handle nextHandle(Handle tag){
    static Handle generation=Handle(std::chrono::high_resolution_clock::now().time_since_epoch().count())&0x0000ffffffffffffULL;
    if(generation==0x0000ffffffffffffULL)throw std::runtime_error("resident generation exhausted");
    return (tag<<48)|++generation;
}
Handle handleValue(const mxArray* a,Handle tag){
    if(mxGetClassID(a)!=mxUINT64_CLASS||mxIsComplex(a)||mxGetNumberOfElements(a)!=1)
        throw std::runtime_error("typed uint64 resident handle required");
    Handle h=*static_cast<const Handle*>(mxGetData(a));
    if((h>>48)!=tag)throw std::runtime_error("resident handle type mismatch");return h;
}
mxArray* handleArray(Handle h){auto a=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);*static_cast<Handle*>(mxGetData(a))=h;return a;}
ResidentAMS& amsContext(const mxArray* a){auto it=amsContexts.find(handleValue(a,AMS_TAG));
    if(it==amsContexts.end())throw std::runtime_error("stale or released AMS handle");it->second->checkDevice();return *it->second;}
void clearResidentContexts();
void registerResidentExit(){static bool registered=false;if(!registered){mexAtExit(clearResidentContexts);registered=true;}}
mxArray* gpuCopy(const double* p,int n){Output out(n);ck(cudaMemcpy(out.data(),p,size_t(n)*sizeof(double),cudaMemcpyDeviceToDevice));return mxGPUCreateMxArrayOnGPU(out.a);}
bool amsResidentCommand(const std::string& command,int nlhs,mxArray** out,int nrhs,const mxArray** args){
    if(command=="create"){
        if(nrhs!=2||nlhs!=1)throw std::runtime_error("handle=radia_cuda_ams('create',hierarchy)");
        auto h=std::make_unique<ResidentHierarchy>(args[1]);auto context=std::make_unique<ResidentAMS>(std::move(h));
        auto id=nextHandle(AMS_TAG);registerResidentExit();amsContexts.emplace(id,std::move(context));mexLock();out[0]=handleArray(id);return true;
    }
    if(command=="release"){
        if(nrhs!=2||nlhs!=0)throw std::runtime_error("radia_cuda_ams('release',handle)");
        amsContext(args[1]);amsContexts.erase(handleValue(args[1],AMS_TAG));mexUnlock();return true;
    }
    if(command=="update"){
        if(nrhs!=3||nlhs!=0)throw std::runtime_error("radia_cuda_ams('update',handle,hierarchy)");
        amsContext(args[1]).replace(args[2]);return true;
    }
    if(command!="solve")return false;
    if(nrhs!=7||nlhs!=2)throw std::runtime_error("[x,info]=radia_cuda_ams('solve',handle,operator,b,tol,maxiter,use_graph)");
    auto& c=amsContext(args[1]);auto allocationsBefore=nativeBufferAllocations.load(std::memory_order_relaxed);double tol=scalar(args[4]);int limit=integer(args[5]);
    if(!(tol>0&&tol<=1e-6)||limit<1||limit==INT_MAX||!mxIsLogicalScalar(args[6]))throw std::runtime_error("invalid resident solver options");
    ck(cudaMemset(c.validation.p,0,sizeof(int)));Input input(c.validation.p);CSR op(device_node(args[2],input),input.validation);int n=c.hierarchy->ams.a.n;
    if(op.n!=n||op.m!=n)throw std::runtime_error("resident operator dimension mismatch");
    const double* b;finite_vec(input,args[3],n,b);verify_mask<<<(n+T-1)/T,T>>>(n,c.hierarchy->ams.free,b,input.validation.p);input.validation.check();
    bool useGraph=mxIsLogicalScalarTrue(args[6]);
    if(useGraph&&!c.graph)c.graph=std::make_unique<AMSGraph>(c.hierarchy->ams,c.work.r.p,c.work.z.p);
    auto stats=residentPCG(c.hierarchy->ams,op,b,c.work,tol,limit,useGraph?c.graph.get():nullptr);++c.solves;
    ck(cudaGetLastError());ck(cudaDeviceSynchronize());out[0]=gpuCopy(c.work.x.p,n);
    const char* names[]={"iterations","true_relative_residual","workspace_bytes","workspace_allocations","solve_count","hierarchy_builds","graph_cached","native_buffer_allocations_during_solve","coarse_action_cache_bytes","iteration_host_syncs"};
    out[1]=mxCreateStructMatrix(1,1,10,names);
    double values[]={double(stats.iterations),stats.residual,double(c.work.bytes()),1.,double(c.solves),double(c.rebuilds)};
    for(int k=0;k<6;++k)mxSetField(out[1],0,names[k],mxCreateDoubleScalar(values[k]));
    mxSetField(out[1],0,names[6],mxCreateLogicalScalar(useGraph));
    mxSetField(out[1],0,names[7],mxCreateDoubleScalar(double(nativeBufferAllocations.load(std::memory_order_relaxed)-allocationsBefore)));
    mxSetField(out[1],0,names[8],mxCreateDoubleScalar(double(c.hierarchy->ams.coarseCacheBytes())));
    mxSetField(out[1],0,names[9],mxCreateDoubleScalar(1));return true;
}
