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
    if(command!="assemble_hcurl")return false;
    if(nr!=8||nl!=3)throw std::invalid_argument("[model,fem,info]=radia_cuda_ams('assemble_hcurl',vol,labels,dirichlet,nu,mass,J,gamma)");
    auto totalStart=std::chrono::steady_clock::now();
    if(!mxIsCell(in[2]))throw std::invalid_argument("Boundary label contract must be a cell array of char rows");
    std::vector<std::string> labels;for(size_t i=0;i<mxGetNumberOfElements(in[2]);++i){auto* label=mxGetCell(in[2],i);if(!label)throw std::invalid_argument("Empty boundary label cell");labels.push_back(word(label));}
    double nu=scalar(in[4]),m=scalar(in[5]);
    if(nu<=0||m<=0)throw std::invalid_argument("Positive SI reluctivity and positive gauge mass required");
    if(!mxIsDouble(in[6])||mxIsComplex(in[6])||mxIsSparse(in[6])||mxGetNumberOfElements(in[6])!=3)
        throw std::invalid_argument("J must be a full three-component double vector in A/m^2");
    // in[7] is the explicit polynomial nonlinearity coefficient, zero for linear runs.
    double gamma=scalar(in[7]);if(gamma<0)throw std::invalid_argument("Nonnegative constitutive gamma required");
    const double* j=mxGetPr(in[6]);for(int c=0;c<3;++c)if(!std::isfinite(j[c]))throw std::invalid_argument("Nonfinite source current");
    auto s=RadiaStageHCurl(word(in[1]),labels,word(in[3]));
    double hostSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-totalStart).count();
    size_t nz=s.ci.size(),freeBytes=s.free.size()*sizeof(int);
    size_t mappedBytes=(s.shape.size()+s.weights.size()+s.curls.size()+s.volume.size())*sizeof(double);
    size_t topologyBytes=(s.dofs.size()+s.slots.size()+s.incident_rp.size()+s.incident.size()+s.rp.size()+s.ci.size())*sizeof(int)+freeBytes;
    size_t scratchBytes=size_t(s.e)*(36*2+6)*sizeof(double);
    size_t outputBytes=(2*nz+s.n)*sizeof(double);
    size_t ownedBytes=mappedBytes+topologyBytes+scratchBytes+outputBytes+
        (nz+s.n)*sizeof(double)+(s.rp.size()+s.ci.size())*sizeof(int)+sizeof(int);
    size_t available,total;ck(cudaMemGetInfo(&available,&total));size_t reserve=std::max<size_t>(2ULL<<30,total/5);
    if(available<=reserve||ownedBytes>available-reserve)throw std::runtime_error("HCurl assembly exceeds CUDA memory reserve max(2 GiB,20%)");
    auto start=std::chrono::steady_clock::now();
    auto shape=femVector(s.shape,mxDOUBLE_CLASS,true),weights=femVector(s.weights,mxDOUBLE_CLASS,true);
    auto curls=femVector(s.curls,mxDOUBLE_CLASS,true),volume=femVector(s.volume,mxDOUBLE_CLASS,true);
    auto dofs=femVector(s.dofs,mxINT32_CLASS,true),slots=femVector(s.slots,mxINT32_CLASS,true);
    auto irp=femVector(s.incident_rp,mxINT32_CLASS,true),incident=femVector(s.incident,mxINT32_CLASS,true);
    auto free=femVector(s.free,mxINT32_CLASS,true),rp=femVector(s.rp,mxINT32_CLASS,true);
    // ci is created only for the returned mass CSR below.
    Buffer lm(size_t(s.e)*36),lc(size_t(s.e)*36),lr(size_t(s.e)*6),mass(nz),A(nz),rhs(s.n);
    Input view;auto dp=[&](const FEMArray& a,size_t n){return view.get<double>(a.get(),mxDOUBLE_CLASS,n);};
    auto ip=[&](const FEMArray& a,size_t n){return view.get<int>(a.get(),mxINT32_CLASS,n);};
    femLocal<<<(s.e*36+T-1)/T,T>>>(s.e,s.nq,dp(shape,s.shape.size()),dp(weights,s.weights.size()),dp(curls,s.curls.size()),dp(volume,s.volume.size()),m,j[0],j[1],j[2],lm.p,lc.p,lr.p);
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
    femPut(fem.get(),"mass",femCSR(s.n,s.n,s.rp,s.ci,femOwn(gpuCopy(mass.p,nz)),true));
    femPut(fem.get(),"elements",femOwn(mxCreateDoubleScalar(s.e)));femPut(fem.get(),"dofs",std::move(dofs));femPut(fem.get(),"slots",std::move(slots));
    femPut(fem.get(),"curls",std::move(curls));femPut(fem.get(),"volume",std::move(volume));femPut(fem.get(),"free",std::move(free));
    femPut(fem.get(),"rhs",femOwn(gpuCopy(rhs.p,s.n)));femPut(fem.get(),"nu0",femOwn(mxCreateDoubleScalar(nu)));femPut(fem.get(),"gamma",femOwn(mxCreateDoubleScalar(gamma)));
    auto info=femStruct({"schema","host_basis_seconds","upload_assembly_seconds","total_seconds","gpu_peak_owned_bytes_upper_bound","upload_bytes","download_bytes","taskmanager_started","python_adapter_used","assembly"});
    femPut(info.get(),"schema",femOwn(mxCreateString("radia.hcurl-cuda-assembly.v1")));
    auto number=[&](const char* k,double x){femPut(info.get(),k,femOwn(mxCreateDoubleScalar(x)));};
    number("host_basis_seconds",hostSeconds);number("upload_assembly_seconds",deviceSeconds);
    number("total_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-totalStart).count());
    // Includes output copies concurrently with all setup buffers, and validation flag.
    number("gpu_peak_owned_bytes_upper_bound",ownedBytes);
    number("upload_bytes",mappedBytes+topologyBytes+(s.rp.size()+s.ci.size())*sizeof(int));number("download_bytes",nz*sizeof(double));
    femPut(info.get(),"taskmanager_started",femOwn(mxCreateLogicalScalar(false)));femPut(info.get(),"python_adapter_used",femOwn(mxCreateLogicalScalar(false)));
    femPut(info.get(),"assembly",femOwn(mxCreateString("NGSolve affine HCurl mapped shapes; CUDA FP64 element quadrature and deterministic row reduction")));
    out[0]=model.release();out[1]=fem.release();out[2]=info.release();return true;
}
