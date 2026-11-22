/* SPDX-License-Identifier: MPL-2.0 */
// Lowest-order affine-tetrahedron assembly from NGSolve-owned mapped curls.
__global__ void verify_elements(int elements,int n,int nz,const int* dofs,const int* slots,
    const double* curls,const double* volume,const int* rp,const int* ci,int* error){
    int e=blockIdx.x*blockDim.x+threadIdx.x;if(e>=elements)return;
    if(!isfinite(volume[e])||volume[e]<=0)atomicExch(error,1);
    for(int a=0;a<6;++a){int i=dofs[e*6+a];
        if(i<0||i>=n){atomicExch(error,1);continue;}
        for(int k=0;k<3;++k)if(!isfinite(curls[e*18+a*3+k]))atomicExch(error,1);
        for(int b=0;b<6;++b){int j=dofs[e*6+b],p=slots[e*36+a*6+b];
            if(p<0||p>=nz||p<rp[i]||p>=rp[i+1])atomicExch(error,1);
            else if(j<0||j>=n||ci[p]!=j)atomicExch(error,1);
        }
    }
}
__global__ void residual_seed(int n,const double* massx,const double* rhs,const int* free,double* r){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)r[i]=free[i]?massx[i]-rhs[i]:0;
}
template<bool Assemble> __global__ void nonlinear_tets(int elements,const int* dofs,const int* slots,
    const double* curls,const double* volume,const int* free,const double* x,double nu0,double gamma,
    double* residual,double* values,double* energy){
    int e=blockIdx.x*blockDim.x+threadIdx.x;if(e>=elements)return;
    double B[3]={0,0,0};
    for(int a=0;a<6;++a)for(int k=0;k<3;++k)B[k]+=curls[e*18+a*3+k]*x[dofs[e*6+a]];
    double b2=B[0]*B[0]+B[1]*B[1]+B[2]*B[2],nu=nu0+gamma*b2,w=volume[e];
    energy[e]=w*(0.5*nu0*b2+0.25*gamma*b2*b2);
    if(Assemble)for(int a=0;a<6;++a){
        const double* ca=curls+e*18+a*3;double ba=ca[0]*B[0]+ca[1]*B[1]+ca[2]*B[2];
        int i=dofs[e*6+a];if(free[i])atomicAdd(residual+i,w*nu*ba);
        for(int b=0;b<6;++b){const double* cb=curls+e*18+b*3;
            double bb=cb[0]*B[0]+cb[1]*B[1]+cb[2]*B[2];
            double dot=ca[0]*cb[0]+ca[1]*cb[1]+ca[2]*cb[2];
            atomicAdd(values+slots[e*36+a*6+b],w*(nu*dot+2*gamma*ba*bb));
        }
    }
}
__global__ void project_csr(int n,const int* rp,const int* ci,const int* free,const double* raw,double* values){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    for(int p=rp[i];p<rp[i+1];++p)values[p]=(free[i]&&free[ci[p]])?raw[p]:0;
}
__global__ void verify_edge_vertices(int edges,int vertices,const int* endpoints,int* invalid){
    const int edge=blockIdx.x*blockDim.x+threadIdx.x;
    if(edge>=edges)return;
    const int a=endpoints[2*edge],b=endpoints[2*edge+1];
    if(a<0||b<0||a>=vertices||b>=vertices||a==b)atomicExch(invalid,1);
}
struct NonlinearState {
    OwnedArray owner;Input input;CSR mass;int n,e,device;
    const int *dofs,*slots,*free;const double *curls,*volume,*rhs;double nu0,gamma;
    Buffer residual,raw,tangent,massx,energy;Blas blas;unsigned evaluations=0;
    radia_fem::Ledger amsLedger;
    std::unique_ptr<radia_fem::Ams> deviceAms;
    unsigned amsEvaluation=0;
    explicit NonlinearState(const mxArray* source):owner(source),mass(device_node(field(owner.value,"mass"),input),input.validation),
        n(mass.n),e(integer(field(owner.value,"elements"))),residual(n),raw(mass.nz),tangent(mass.nz),massx(n),energy(e){
        if(word(field(owner.value,"schema"))!="radia.ams-nonlinear-tet6.v1"||mass.m!=n||e<1||e>INT_MAX/36)
            throw std::runtime_error("unsupported nonlinear FE schema or dimensions");
        nu0=scalar(field(owner.value,"nu0"));gamma=scalar(field(owner.value,"gamma"));
        if(nu0<=0||gamma<0)throw std::runtime_error("constitutive tangent must be positive");
        dofs=input.get<int>(field(owner.value,"dofs"),mxINT32_CLASS,size_t(e)*6);
        slots=input.get<int>(field(owner.value,"slots"),mxINT32_CLASS,size_t(e)*36);
        curls=input.get<double>(field(owner.value,"curls"),mxDOUBLE_CLASS,size_t(e)*18);
        volume=input.get<double>(field(owner.value,"volume"),mxDOUBLE_CLASS,e);
        free=input.get<int>(field(owner.value,"free"),mxINT32_CLASS,n);
        finite_vec(input,field(owner.value,"rhs"),n,rhs);
        input.validation.check(); // CSR bounds checked before element slot dereferences.
        verify_elements<<<(e+T-1)/T,T>>>(e,n,mass.nz,dofs,slots,curls,volume,mass.rp,mass.ci,input.validation.p);
        verify_mask<<<(n+T-1)/T,T>>>(n,free,rhs,input.validation.p);input.validation.check();ck(cudaGetDevice(&device));
    }
    void checkDevice(){int d;ck(cudaGetDevice(&d));if(d!=device)throw std::runtime_error("nonlinear state belongs to another CUDA device");}
    const double* checkedX(const mxArray* a,Input& view){const double* x;finite_vec(view,a,n,x);
        verify_mask<<<(n+T-1)/T,T>>>(n,free,x,view.validation.p);view.validation.check();return x;}
    double evaluate(const double* x,bool assemble){
        mass.mult(x,massx.p);
        if(assemble){
            residual_seed<<<(n+T-1)/T,T>>>(n,massx.p,rhs,free,residual.p);
            ck(cudaMemcpy(raw.p,mass.a,size_t(mass.nz)*sizeof(double),cudaMemcpyDeviceToDevice));
            nonlinear_tets<true><<<(e+T-1)/T,T>>>(e,dofs,slots,curls,volume,free,x,nu0,gamma,residual.p,raw.p,energy.p);
            project_csr<<<(n+T-1)/T,T>>>(n,mass.rp,mass.ci,free,raw.p,tangent.p);++evaluations;
        }else nonlinear_tets<false><<<(e+T-1)/T,T>>>(e,dofs,slots,curls,volume,free,x,nu0,gamma,nullptr,nullptr,energy.p);
        cb(cublasSetPointerMode(blas.h,CUBLAS_POINTER_MODE_HOST));double integral,mx,bx;
        cb(cublasDasum(blas.h,e,energy.p,1,&integral));cb(cublasDdot(blas.h,n,x,1,massx.p,1,&mx));
        cb(cublasDdot(blas.h,n,x,1,rhs,1,&bx));double value=integral+0.5*mx-bx;
        ck(cudaGetLastError());if(!std::isfinite(value))throw std::runtime_error("nonlinear energy overflow");return value;
    }
};
constexpr Handle NONLINEAR_TAG=0xA65B;
std::unordered_map<Handle,std::unique_ptr<NonlinearState>> nonlinearContexts;
void clearResidentContexts(){amsContexts.clear();nonlinearContexts.clear();}
NonlinearState& nonlinearContext(const mxArray* a){auto it=nonlinearContexts.find(handleValue(a,NONLINEAR_TAG));
    if(it==nonlinearContexts.end())throw std::runtime_error("stale or released nonlinear handle");it->second->checkDevice();return *it->second;}
struct ResidentHierarchyTiming {
    bool built;
    bool refreshed;
    double seconds;
};

ResidentHierarchyTiming prepareNonlinearHierarchy(NonlinearState& c){
    if(!c.evaluations)throw std::runtime_error("evaluate the tangent before preparing AMS");
    Input view;
    const auto start=std::chrono::steady_clock::now();
    const bool built=!c.deviceAms;
    const bool refreshed=!built&&c.amsEvaluation!=c.evaluations;
    if(built){
        const int vertices=integer(field(c.owner.value,"vertices"));
        if(vertices<1)throw std::runtime_error("resident AMS needs mesh vertices");
        const int* edges=view.get<int>(field(c.owner.value,"edge_vertices"),mxINT32_CLASS,size_t(c.n)*2);
        if(vertices>INT_MAX/3)throw std::runtime_error("resident HCurl vertex count exceeds int32 capacity");
        const double* coords=view.get<double>(field(c.owner.value,"coordinates"),mxDOUBLE_CLASS,size_t(vertices)*3);
        verify_edge_vertices<<<(c.n+T-1)/T,T>>>(c.n,vertices,edges,view.validation.p);
        verify_vec<<<(vertices*3+T-1)/T,T>>>(vertices*3,coords,false,view.validation.p);
        view.validation.check();
        auto hierarchy=std::make_unique<radia_fem::Ams>();
        const radia_fem::CsrView op{c.n,c.n,c.mass.rp,c.mass.ci,c.raw.p};
        hierarchy->build(c.amsLedger,op,c.mass.nz,c.free,c.raw.p,c.n,vertices,edges,coords,64,0.25);
        c.deviceAms=std::move(hierarchy);
    }else if(refreshed){
        // If refresh fails, discard the partially updated hierarchy.
        auto hierarchy=std::move(c.deviceAms);
        hierarchy->refresh_values(c.amsLedger,c.raw.p);
        c.deviceAms=std::move(hierarchy);
    }
    c.amsEvaluation=c.evaluations;
    ck(cudaDeviceSynchronize());
    const double hierarchySeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return {built,refreshed,hierarchySeconds};
}

bool nonlinearCommand(const std::string& command,int nlhs,mxArray** out,int nrhs,const mxArray** args){
    if(command=="nonlinear_create"){
        if(nrhs!=2||nlhs!=1)throw std::runtime_error("handle=radia_cuda_ams('nonlinear_create',model)");
        auto context=std::make_unique<NonlinearState>(args[1]);Handle h=nextHandle(NONLINEAR_TAG);
        registerResidentExit();nonlinearContexts.emplace(h,std::move(context));mexLock();out[0]=handleArray(h);return true;
    }
    if(command=="nonlinear_release"){
        if(nrhs!=2||nlhs!=0)throw std::runtime_error("radia_cuda_ams('nonlinear_release',handle)");
        nonlinearContext(args[1]);nonlinearContexts.erase(handleValue(args[1],NONLINEAR_TAG));mexUnlock();return true;
    }
    if(command=="nonlinear_prepare_resident"){
        if(nrhs!=2||nlhs!=1)throw std::runtime_error("info=radia_cuda_ams('nonlinear_prepare_resident',handle)");
        const auto timing=prepareNonlinearHierarchy(nonlinearContext(args[1]));
        const char* fields[]={"hierarchy_built","hierarchy_refreshed","hierarchy_seconds"};
        out[0]=mxCreateStructMatrix(1,1,3,fields);
        mxSetField(out[0],0,fields[0],mxCreateLogicalScalar(timing.built));
        mxSetField(out[0],0,fields[1],mxCreateLogicalScalar(timing.refreshed));
        mxSetField(out[0],0,fields[2],mxCreateDoubleScalar(timing.seconds));
        return true;
    }
    if(command=="nonlinear_solve_resident"){
        if(nrhs!=5||nlhs!=2)
            throw std::runtime_error("[x,info]=radia_cuda_ams('nonlinear_solve_resident',handle,rhs,tol,maxiter)");
        auto& c=nonlinearContext(args[1]);
        const double tol=scalar(args[3]);
        const int maxiter=integer(args[4]);
        if(!c.evaluations||!(tol>0&&tol<=1e-6)||maxiter<1||maxiter==INT_MAX)
            throw std::runtime_error("evaluate the tangent before solving with valid CG controls");
        Input view;
        const double* rhs;
        finite_vec(view,args[2],c.n,rhs);
        verify_mask<<<(c.n+T-1)/T,T>>>(c.n,c.free,rhs,view.validation.p);
        view.validation.check();
        const auto timing=prepareNonlinearHierarchy(c);
        Output x(c.n);
        auto action=[&](const double* in,double* out){
            radia_fem::multiply<<<radia_fem::blocks_of(c.n),radia_fem::TH>>>(c.n,c.mass.rp,c.mass.ci,c.tangent.p,in,out);
        };
        auto precondition=[&](const double* in,double* out){c.deviceAms->apply(in,out);};
        radia_fem::Blas blas(c.amsLedger);
        const double rhsNormSquared=blas.dot(c.n,rhs,rhs);
        radia_fem::SolveReport report;
        double trueRelative=0;
        if(rhsNormSquared!=0){
            report=radia_fem::pcg(c.amsLedger,c.n,action,precondition,rhs,x.data(),tol,maxiter);
            radia_fem::Dev<double> residual(size_t(c.n),"HCurl true residual");
            action(x.data(),residual.p);
            radia_fem::k_axpy<<<radia_fem::blocks_of(c.n),radia_fem::TH>>>(c.n,-1.,rhs,residual.p);
            trueRelative=std::sqrt(blas.dot(c.n,residual.p,residual.p))/report.initial_norm;
            if(!(trueRelative<=tol))throw std::runtime_error("resident HCurl CG failed the true relative residual check");
        }
        out[0]=mxGPUCreateMxArrayOnGPU(x.a);
        const char* names[]={"iterations","true_relative_residual","hierarchy_built","hierarchy_refreshed","hierarchy_seconds","gpu_peak_owned_bytes"};
        out[1]=mxCreateStructMatrix(1,1,6,names);
        const double values[]={double(report.iterations),trueRelative,double(timing.built),double(timing.refreshed),
                               timing.seconds,double(radia_fem::device_memory().peak)};
        for(int k=0;k<6;++k)mxSetField(out[1],0,names[k],mxCreateDoubleScalar(values[k]));
        return true;
    }
#ifdef RADIA_AMS_VALIDATION
    if(command=="nonlinear_raw"){
        if(nrhs!=2||nlhs!=1)throw std::runtime_error("values=radia_cuda_ams('nonlinear_raw',handle)");
        auto& c=nonlinearContext(args[1]);if(!c.evaluations)throw std::runtime_error("evaluate tangent before downloading values");
        out[0]=mxCreateDoubleMatrix(c.mass.nz,1,mxREAL);
        ck(cudaMemcpy(mxGetPr(out[0]),c.raw.p,size_t(c.mass.nz)*sizeof(double),cudaMemcpyDeviceToHost));return true;
    }
#endif
    if(command=="resident_status"){
        if(nrhs!=1||nlhs!=1)throw std::runtime_error("status=radia_cuda_ams('resident_status')");
        const char* fields[]={"ams_contexts","nonlinear_contexts"};out[0]=mxCreateStructMatrix(1,1,2,fields);
        mxSetField(out[0],0,fields[0],mxCreateDoubleScalar(double(amsContexts.size())));
        mxSetField(out[0],0,fields[1],mxCreateDoubleScalar(double(nonlinearContexts.size())));return true;
    }
    if(command!="nonlinear_evaluate"&&command!="nonlinear_energy")return false;
    bool assemble=command=="nonlinear_evaluate";
    if(nrhs!=3||nlhs!=(assemble?3:1))throw std::runtime_error("nonlinear evaluation output count mismatch");
    auto& c=nonlinearContext(args[1]);ck(cudaMemset(c.input.validation.p,0,sizeof(int)));
    Input view(c.input.validation.p);const double* x=c.checkedX(args[2],view);double value=c.evaluate(x,assemble);
    if(assemble){out[0]=gpuCopy(c.residual.p,c.n);out[1]=gpuCopy(c.tangent.p,c.mass.nz);out[2]=mxCreateDoubleScalar(value);}
    else out[0]=mxCreateDoubleScalar(value);return true;
}
