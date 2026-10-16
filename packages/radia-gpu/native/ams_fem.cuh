/* SPDX-License-Identifier: MPL-2.0 */
// Included in the AMS gateway: no NGSolve headers enter the CUDA translation unit.
using FEMArray=std::unique_ptr<mxArray,decltype(&mxDestroyArray)>;
FEMArray femOwn(mxArray* p){if(!p)throw std::bad_alloc();return FEMArray(p,mxDestroyArray);}
FEMArray femStruct(std::initializer_list<const char*> fields){std::vector<const char*> f(fields);return femOwn(mxCreateStructMatrix(1,1,int(f.size()),f.data()));}
void femPut(mxArray* p,const char* name,FEMArray value){mxSetField(p,0,name,value.release());}
template<class V> FEMArray femVector(const std::vector<V>& v,mxClassID type,bool gpu){
    if(!gpu){auto p=femOwn(mxCreateNumericMatrix(v.size(),1,type,mxREAL));std::copy(v.begin(),v.end(),static_cast<V*>(mxGetData(p.get())));return p;}
    mwSize n=v.size();auto* g=mxGPUCreateGPUArray(1,&n,type,mxREAL,MX_GPU_DO_NOT_INITIALIZE);
    try{ck(cudaMemcpy(mxGPUGetData(g),v.data(),v.size()*sizeof(V),cudaMemcpyHostToDevice));
        auto p=femOwn(mxGPUCreateMxArrayOnGPU(g));mxGPUDestroyGPUArray(g);return p;
    }catch(...){mxGPUDestroyGPUArray(g);throw;}
}
FEMArray femCSR(int n,int m,const std::vector<int>& rp,const std::vector<int>& ci,FEMArray values,bool gpu){
    auto p=femStruct({"rows","cols","rowptr","columns","values"});
    femPut(p.get(),"rows",femOwn(mxCreateDoubleScalar(n)));femPut(p.get(),"cols",femOwn(mxCreateDoubleScalar(m)));
    femPut(p.get(),"rowptr",femVector(rp,mxINT32_CLASS,gpu));femPut(p.get(),"columns",femVector(ci,mxINT32_CLASS,gpu));femPut(p.get(),"values",std::move(values));return p;
}
#include "ams_geometry.cuh"
__global__ void femLocal(int e,int nq,const double* shape,const double* weights,const double* curl,
    const double* volume,double massCoefficient,double jx,double jy,double jz,double* localMass,double* localCurl,double* localRhs){
    int k=blockIdx.x*blockDim.x+threadIdx.x;if(k>=e*36)return;int cell=k/36,a=(k%36)/6,b=k%6;
    double m=0,r=0;
    for(int q=0;q<nq;++q){const double* s=shape+(size_t(cell)*nq+q)*18;double dot=0;
        for(int c=0;c<3;++c)dot+=s[a*3+c]*s[b*3+c];m+=weights[cell*nq+q]*dot;
        if(b==0)r+=weights[cell*nq+q]*(s[a*3]*jx+s[a*3+1]*jy+s[a*3+2]*jz);}
    double dot=0;for(int c=0;c<3;++c)dot+=curl[cell*18+a*3+c]*curl[cell*18+b*3+c];
    localMass[k]=massCoefficient*m;localCurl[k]=volume[cell]*dot;if(b==0)localRhs[cell*6+a]=r;
}
// One row owner sums elements in mesh order; no floating-point atomics.
__global__ void femReduce(int n,const int* rp,const int* slots,const int* incidentRp,const int* incident,
    const int* free,const double* lm,const double* lc,const double* lr,double nu,double* mass,double* A,double* rhs){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    for(int p=rp[i];p<rp[i+1];++p){mass[p]=0;A[p]=0;}double r=0;
    for(int p=incidentRp[i];p<incidentRp[i+1];++p){int ea=incident[p],e=ea/6,a=ea%6;r+=lr[ea];
        for(int b=0;b<6;++b){int k=e*36+a*6+b,s=slots[k];mass[s]+=lm[k];A[s]+=lm[k]+nu*lc[k];}}
    rhs[i]=free[i]?r:0;
}
bool femCommand(const std::string& command,int nl,mxArray** out,int nr,const mxArray** in){
    if(geometryCommand(command,nl,out,nr,in))return true;
    const bool reused=command=="assemble_hcurl_geometry";
    if(command!="assemble_hcurl"&&!reused)return false;
    if(nr!=(reused?6:8)||nl!=3)throw std::invalid_argument("HCurl assembly requires model/fem/info outputs and explicit geometry or vol/labels/dirichlet inputs");
    auto totalStart=std::chrono::steady_clock::now();
    std::vector<std::string> labels;
    if(!reused){
        if(!mxIsCell(in[2]))throw std::invalid_argument("Boundary label contract must be a cell array of char rows");
        for(size_t i=0;i<mxGetNumberOfElements(in[2]);++i){auto* label=mxGetCell(in[2],i);if(!label)throw std::invalid_argument("Empty boundary label cell");labels.push_back(word(label));}
    }
    const int first=reused?2:4;
    double nu=scalar(in[first]),m=scalar(in[first+1]);
    if(nu<=0||m<=0)throw std::invalid_argument("Positive SI reluctivity and positive gauge mass required");
    if(!mxIsDouble(in[first+2])||mxIsComplex(in[first+2])||mxIsSparse(in[first+2])||mxGetNumberOfElements(in[first+2])!=3)
        throw std::invalid_argument("J must be a full three-component double vector in A/m^2");
    // Last input is the explicit polynomial nonlinearity coefficient, zero for linear runs.
    double gamma=scalar(in[first+3]);if(gamma<0)throw std::invalid_argument("Nonnegative constitutive gamma required");
    const double* j=mxGetPr(in[first+2]);for(int c=0;c<3;++c)if(!std::isfinite(j[c]))throw std::invalid_argument("Nonfinite source current");
    HCurlGeometry* geometry=reused?&geometryContext(in[1]):nullptr;
    const double checkSecondsBefore=geometry?geometry->livenessCheckSeconds:0;
    const size_t checkCallsBefore=geometry?geometry->livenessCheckCalls:0;
    RadiaHCurlStage temporary;
    if(!reused)temporary=RadiaStageHCurl(word(in[1]),labels,word(in[3]));
    const auto& s=reused?geometry->stage:temporary;
    double hostSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-totalStart).count();
    size_t nz=s.ci.size(),freeBytes=s.free.size()*sizeof(int);
    const size_t shapeCount=size_t(s.e)*s.nq*18,weightCount=size_t(s.e)*s.nq,curlCount=size_t(s.e)*18,volumeCount=size_t(s.e);
    size_t mappedBytes=(shapeCount+weightCount+curlCount+volumeCount)*sizeof(double);
    size_t topologyBytes=(s.dofs.size()+s.slots.size()+s.incident_rp.size()+s.incident.size()+s.rp.size()+s.ci.size())*sizeof(int)+freeBytes;
    size_t scratchBytes=size_t(s.e)*(36*2+6)*sizeof(double);
    size_t outputBytes=(2*nz+s.n)*sizeof(double);
    size_t ownedBytes=mappedBytes+topologyBytes+scratchBytes+outputBytes+
        (nz+s.n)*sizeof(double)+(s.rp.size()+s.ci.size())*sizeof(int)+sizeof(int);
    // Include the resident geometry even if mxDuplicateArray materializes
    // independent device storage for a returned/read-only array view.
    if(geometry)ownedBytes+=geometry->deviceBytes;
    size_t available,total;ck(cudaMemGetInfo(&available,&total));size_t reserve=std::max<size_t>(2ULL<<30,total/5);
    if(available<=reserve||ownedBytes>available-reserve)throw std::runtime_error("HCurl assembly exceeds CUDA memory reserve max(2 GiB,20%)");
    auto start=std::chrono::steady_clock::now();
    auto mapped=[&](const char* key,const auto& values,mxClassID type){return geometry?geometry->copy(key):femVector(values,type,true);};
    auto shape=mapped("shape",s.shape,mxDOUBLE_CLASS),weights=mapped("weights",s.weights,mxDOUBLE_CLASS);
    auto curls=mapped("curls",s.curls,mxDOUBLE_CLASS),volume=mapped("volume",s.volume,mxDOUBLE_CLASS);
    auto dofs=mapped("dofs",s.dofs,mxINT32_CLASS),slots=mapped("slots",s.slots,mxINT32_CLASS);
    auto irp=mapped("incident_rp",s.incident_rp,mxINT32_CLASS),incident=mapped("incident",s.incident,mxINT32_CLASS);
    auto free=mapped("free",s.free,mxINT32_CLASS),rp=mapped("rp",s.rp,mxINT32_CLASS);
    // ci is created only for the returned mass CSR below.
    Buffer lm(size_t(s.e)*36),lc(size_t(s.e)*36),lr(size_t(s.e)*6),mass(nz),A(nz),rhs(s.n);
    Input view;auto dp=[&](const FEMArray& a,size_t n){return view.get<double>(a.get(),mxDOUBLE_CLASS,n);};
    auto ip=[&](const FEMArray& a,size_t n){return view.get<int>(a.get(),mxINT32_CLASS,n);};
    femLocal<<<(s.e*36+T-1)/T,T>>>(s.e,s.nq,dp(shape,shapeCount),dp(weights,weightCount),dp(curls,curlCount),dp(volume,volumeCount),m,j[0],j[1],j[2],lm.p,lc.p,lr.p);
    femReduce<<<(s.n+T-1)/T,T>>>(s.n,ip(rp,s.rp.size()),ip(slots,s.slots.size()),ip(irp,s.incident_rp.size()),ip(incident,s.incident.size()),ip(free,s.free.size()),lm.p,lc.p,lr.p,nu,mass.p,A.p,rhs.p);
    verify_vec<<<(int(nz)+T-1)/T,T>>>(int(nz),mass.p,false,view.validation.p);
    verify_vec<<<(int(nz)+T-1)/T,T>>>(int(nz),A.p,false,view.validation.p);
    verify_vec<<<(s.n+T-1)/T,T>>>(s.n,rhs.p,false,view.validation.p);view.validation.check();ck(cudaGetLastError());ck(cudaDeviceSynchronize());
    double deviceSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    auto model=femStruct({"schema","A","G","coordinates","free","cycle"});
    femPut(model.get(),"schema",femOwn(mxCreateString("radia.ams-native-input.v1")));
    auto hostA=femOwn(mxCreateDoubleMatrix(nz,1,mxREAL));ck(cudaMemcpy(mxGetPr(hostA.get()),A.p,nz*sizeof(double),cudaMemcpyDeviceToHost));
    femPut(model.get(),"A",femCSR(s.n,s.n,s.rp,s.ci,std::move(hostA),false));
    int nv=int(s.coordinates.size()/3);
    femPut(model.get(),"G",femCSR(s.n,nv,s.grp,s.gci,femVector(s.gvalues,mxDOUBLE_CLASS,false),false));
    auto xyz=femOwn(mxCreateDoubleMatrix(nv,3,mxREAL));std::copy(s.coordinates.begin(),s.coordinates.end(),mxGetPr(xyz.get()));
    femPut(model.get(),"coordinates",std::move(xyz));femPut(model.get(),"free",femVector(s.free,mxINT32_CLASS,false));
    femPut(model.get(),"cycle",femOwn(mxCreateDoubleScalar(1)));
    auto fem=femStruct({"schema","mass","elements","dofs","slots","curls","volume","free","rhs","nu0","gamma"});
    femPut(fem.get(),"schema",femOwn(mxCreateString("radia.ams-nonlinear-tet6.v1")));
    auto massCSR=geometry?femStruct({"rows","cols","rowptr","columns","values"}):
        femCSR(s.n,s.n,s.rp,s.ci,femOwn(gpuCopy(mass.p,nz)),true);
    if(geometry){
        femPut(massCSR.get(),"rows",femOwn(mxCreateDoubleScalar(s.n)));
        femPut(massCSR.get(),"cols",femOwn(mxCreateDoubleScalar(s.n)));
        femPut(massCSR.get(),"rowptr",geometry->copy("rp"));femPut(massCSR.get(),"columns",geometry->copy("ci"));
        femPut(massCSR.get(),"values",femOwn(gpuCopy(mass.p,nz)));
    }
    femPut(fem.get(),"mass",std::move(massCSR));
    femPut(fem.get(),"elements",femOwn(mxCreateDoubleScalar(s.e)));femPut(fem.get(),"dofs",std::move(dofs));femPut(fem.get(),"slots",std::move(slots));
    femPut(fem.get(),"curls",std::move(curls));femPut(fem.get(),"volume",std::move(volume));femPut(fem.get(),"free",std::move(free));
    femPut(fem.get(),"rhs",femOwn(gpuCopy(rhs.p,s.n)));femPut(fem.get(),"nu0",femOwn(mxCreateDoubleScalar(nu)));femPut(fem.get(),"gamma",femOwn(mxCreateDoubleScalar(gamma)));
    auto info=femStruct({"schema","host_basis_seconds","host_staging_phases","upload_assembly_seconds","total_seconds","gpu_peak_owned_bytes_upper_bound","upload_bytes","download_bytes","python_adapter_used","assembly","geometry_reused","geometry_retained_device_bytes","geometry_retained_host_vector_capacity_bytes","geometry_liveness_check_seconds","geometry_liveness_check_calls"});
    femPut(info.get(),"schema",femOwn(mxCreateString("radia.hcurl-cuda-assembly.v1")));
    auto number=[&](const char* k,double x){femPut(info.get(),k,femOwn(mxCreateDoubleScalar(x)));};
    number("host_basis_seconds",hostSeconds);number("upload_assembly_seconds",deviceSeconds);
    auto phases=femStruct({"mesh_validation_seconds","space_setup_seconds","allocation_mask_seconds",
        "element_basis_incidence_seconds","csr_scatter_seconds","gradient_seconds","coordinates_seconds"});
    auto phase=[&](const char* k,double x){femPut(phases.get(),k,femOwn(mxCreateDoubleScalar(reused?0:x)));};
    phase("mesh_validation_seconds",s.mesh_seconds);phase("space_setup_seconds",s.space_seconds);
    phase("allocation_mask_seconds",s.allocation_seconds);phase("element_basis_incidence_seconds",s.element_seconds);
    phase("csr_scatter_seconds",s.csr_seconds);phase("gradient_seconds",s.gradient_seconds);
    phase("coordinates_seconds",s.coordinate_seconds);
    femPut(info.get(),"host_staging_phases",std::move(phases));
    number("total_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-totalStart).count());
    // Includes output copies concurrently with all setup buffers, and validation flag.
    number("gpu_peak_owned_bytes_upper_bound",ownedBytes);
    // ci is uploaded only with mass CSR; rp is uploaded for reduction and CSR.
    number("upload_bytes",reused?0:mappedBytes+topologyBytes+s.rp.size()*sizeof(int));number("download_bytes",nz*sizeof(double));
    femPut(info.get(),"geometry_reused",femOwn(mxCreateLogicalScalar(reused)));
    number("geometry_retained_device_bytes",geometry?double(geometry->deviceBytes):0);
    number("geometry_retained_host_vector_capacity_bytes",geometry?double(geometry->hostBytes):0);
    number("geometry_liveness_check_seconds",geometry?geometry->livenessCheckSeconds-checkSecondsBefore:0);
    number("geometry_liveness_check_calls",geometry?double(geometry->livenessCheckCalls-checkCallsBefore):0);
    femPut(info.get(),"python_adapter_used",femOwn(mxCreateLogicalScalar(false)));
    femPut(info.get(),"assembly",femOwn(mxCreateString("NGSolve affine HCurl mapped shapes; CUDA FP64 element quadrature and deterministic row reduction")));
    out[0]=model.release();out[1]=fem.release();out[2]=info.release();return true;
}
