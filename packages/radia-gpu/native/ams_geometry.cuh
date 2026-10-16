/* SPDX-License-Identifier: MPL-2.0 */
// Included after the FEM array helpers. Own copied data, never NGSolve objects.
struct HCurlGeometry {
    RadiaHCurlStage stage;
    FEMArray mapped;
    int device;
    size_t deviceBytes=0,hostBytes=0;
    double creationSeconds=0;
    mutable bool invalidated=false;
    mutable double livenessCheckSeconds=0;
    mutable size_t livenessCheckCalls=0;
#ifdef RADIA_AMS_GEOMETRY_LIFECYCLE_TEST
    bool simulateDeviceMismatch=false;
#endif
    HCurlGeometry(const std::string& path,const std::vector<std::string>& labels,const std::string& dirichlet)
        :mapped(nullptr,mxDestroyArray){
        auto start=std::chrono::steady_clock::now();ck(cudaGetDevice(&device));
        stage=RadiaStageHCurl(path,labels,dirichlet);auto& s=stage;
        deviceBytes=(s.shape.size()+s.weights.size()+s.curls.size()+s.volume.size())*sizeof(double)+
            (s.dofs.size()+s.slots.size()+s.incident_rp.size()+s.incident.size()+s.free.size()+s.rp.size()+s.ci.size())*sizeof(int);
        size_t available,total;ck(cudaMemGetInfo(&available,&total));size_t reserve=std::max<size_t>(2ULL<<30,total/5);
        if(available<=reserve||deviceBytes>available-reserve)throw std::runtime_error("HCurl geometry exceeds CUDA memory reserve");
        mapped=femStruct({"shape","weights","curls","volume","dofs","slots","incident_rp","incident","free","rp","ci"});
        auto put=[&](const char* key,auto& values,mxClassID type){
            auto value=femVector(values,type,true);
            // Persist each child explicitly as well as its owning struct.
            // Destroying the struct owns recursive destruction of its children.
            mexMakeArrayPersistent(value.get());femPut(mapped.get(),key,std::move(value));
        };
        put("shape",s.shape,mxDOUBLE_CLASS);put("weights",s.weights,mxDOUBLE_CLASS);
        put("curls",s.curls,mxDOUBLE_CLASS);put("volume",s.volume,mxDOUBLE_CLASS);
        put("dofs",s.dofs,mxINT32_CLASS);put("slots",s.slots,mxINT32_CLASS);
        put("incident_rp",s.incident_rp,mxINT32_CLASS);put("incident",s.incident,mxINT32_CLASS);
        put("free",s.free,mxINT32_CLASS);put("rp",s.rp,mxINT32_CLASS);put("ci",s.ci,mxINT32_CLASS);
        ck(cudaDeviceSynchronize());mexMakeArrayPersistent(mapped.get());
        std::vector<double>().swap(s.shape);std::vector<double>().swap(s.weights);
        std::vector<double>().swap(s.curls);std::vector<double>().swap(s.volume);
        hostBytes=(s.coordinates.capacity()+s.gvalues.capacity())*sizeof(double)+
            (s.dofs.capacity()+s.slots.capacity()+s.incident_rp.capacity()+s.incident.capacity()+s.free.capacity()+s.rp.capacity()+s.ci.capacity()+s.grp.capacity()+s.gci.capacity())*sizeof(int);
        creationSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    }
    void checkDevice()const{
        int current;ck(cudaGetDevice(&current));
#ifdef RADIA_AMS_GEOMETRY_LIFECYCLE_TEST
        if(simulateDeviceMismatch)current=device+1;
#endif
        if(current!=device)throw std::runtime_error("HCurl geometry belongs to another CUDA device");
        if(invalidated)throw std::runtime_error("geometry invalidated by device reset");
    }
    FEMArray copy(const char* key)const{
        checkDevice();const mxArray* value=field(mapped.get(),key);
        mxArray* result=nullptr;mxArray* argument=const_cast<mxArray*>(value);
        const auto checkStart=std::chrono::steady_clock::now();
        auto failure=mexCallMATLABWithTrap(1,&result,1,&argument,"existsOnGPU");
        livenessCheckSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-checkStart).count();
        ++livenessCheckCalls;
        const bool live=!failure&&result&&mxIsLogicalScalarTrue(result)&&mxIsGPUArray(value);
        if(result)mxDestroyArray(result);if(failure)mxDestroyArray(failure);
        if(!live){invalidated=true;throw std::runtime_error("geometry invalidated by device reset");}
        // Read-only input clones only. All CUDA consumers use Input::get,
        // hence mxGPUGetDataReadOnly; kernels never write retained geometry.
        return femOwn(mxDuplicateArray(value));
    }
};
constexpr Handle GEOMETRY_TAG=0xA65C;
std::unordered_map<Handle,std::unique_ptr<HCurlGeometry>> geometryContexts;
size_t geometryContextCount(){return geometryContexts.size();}
HCurlGeometry& geometryContext(const mxArray* a){
    auto it=geometryContexts.find(handleValue(a,GEOMETRY_TAG));
    if(it==geometryContexts.end())throw std::runtime_error("stale or released HCurl geometry handle");
    it->second->checkDevice();return *it->second;
}
void clearGeometryContexts(){
    int previous=0;cudaGetDevice(&previous);
    for(auto& entry:geometryContexts){cudaSetDevice(entry.second->device);entry.second.reset();}
    geometryContexts.clear();cudaSetDevice(previous);
}
bool geometryCommand(const std::string& command,int nl,mxArray** out,int nr,const mxArray** in){
#ifdef RADIA_AMS_GEOMETRY_LIFECYCLE_TEST
    static bool failCreate=false;
    if(command=="geometry_test_device_mismatch"){
        if(nr!=2||nl!=0)throw std::invalid_argument("geometry mismatch test arity");
        auto it=geometryContexts.find(handleValue(in[1],GEOMETRY_TAG));
        if(it==geometryContexts.end())throw std::runtime_error("stale geometry test handle");
        it->second->simulateDeviceMismatch=true;return true;
    }
    if(command=="geometry_test_fail_create"){
        if(nr!=1||nl!=0)throw std::invalid_argument("geometry allocation test arity");
        failCreate=true;return true;
    }
#endif
    if(command=="geometry_release"){
        if(nr!=2||nl!=0)throw std::invalid_argument("radia_cuda_ams('geometry_release',handle)");
        auto it=geometryContexts.find(handleValue(in[1],GEOMETRY_TAG));
        if(it==geometryContexts.end())throw std::runtime_error("stale or released HCurl geometry handle");
        struct RestoreDevice {
            int previous;bool active=true;
            explicit RestoreDevice(int owner){ck(cudaGetDevice(&previous));ck(cudaSetDevice(owner));}
            cudaError_t restore(){auto status=cudaSetDevice(previous);if(status==cudaSuccess)active=false;return status;}
            ~RestoreDevice(){if(active)cudaSetDevice(previous);}
        } restore(it->second->device);
        geometryContexts.erase(it);mexUnlock(); // Balance lock even if restoration fails.
        ck(restore.restore());return true;
    }
    if(command!="geometry_create")return false;
    if(nr!=4||nl!=2)throw std::invalid_argument("[handle,info]=radia_cuda_ams('geometry_create',vol,labels,dirichlet)");
    if(!mxIsCell(in[2]))throw std::invalid_argument("Boundary label contract must be a cell array");
    std::vector<std::string> labels;
    for(size_t i=0;i<mxGetNumberOfElements(in[2]);++i){auto* label=mxGetCell(in[2],i);if(!label)throw std::invalid_argument("Empty boundary label cell");labels.push_back(word(label));}
    auto context=std::make_unique<HCurlGeometry>(word(in[1]),labels,word(in[3]));
    auto info=femStruct({"schema","creation_seconds","retained_device_bytes","retained_host_vector_capacity_bytes","dofs","elements","device"});
    femPut(info.get(),"schema",femOwn(mxCreateString("radia.hcurl-geometry.v1")));
    auto put=[&](const char* key,double value){femPut(info.get(),key,femOwn(mxCreateDoubleScalar(value)));};
    put("creation_seconds",context->creationSeconds);put("retained_device_bytes",double(context->deviceBytes));
    put("retained_host_vector_capacity_bytes",double(context->hostBytes));put("dofs",context->stage.n);
    put("elements",context->stage.e);put("device",context->device);
    auto id=nextHandle(GEOMETRY_TAG);auto handle=femOwn(handleArray(id));registerResidentExit();
#ifdef RADIA_AMS_GEOMETRY_LIFECYCLE_TEST
    if(failCreate){failCreate=false;throw std::bad_alloc();} // Simulate emplace allocation failure before lock.
#endif
    geometryContexts.emplace(id,std::move(context));mexLock();out[0]=handle.release();out[1]=info.release();return true;
}
