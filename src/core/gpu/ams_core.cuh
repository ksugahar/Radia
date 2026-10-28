/* SPDX-License-Identifier: MPL-2.0 */
#pragma once
// Shared unchanged AMS/AMG cycle and PCG algebra, independent of gateways.
// Include inside a private kernel namespace after CUDA/standard headers and snapshot.hpp.
using radia_sparse::Node;
constexpr int T=256;
thread_local cudaStream_t launch_stream=nullptr;
void ck(cudaError_t s){if(s!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(s));}
void cb(cublasStatus_t s){if(s!=CUBLAS_STATUS_SUCCESS)throw std::runtime_error("cuBLAS failed: "+std::to_string(s));}
struct Flag {
    int* p=nullptr; bool owned=true;
    explicit Flag(int* borrowed):p(borrowed),owned(false){}
    explicit Flag(size_t count=1){size_t bytes=std::max<size_t>(1,count)*sizeof(int);ck(cudaMalloc(&p,bytes));ck(cudaMemset(p,0,bytes));}
    ~Flag(){if(p&&owned)cudaFree(p);}
    void check(){int value;ck(cudaMemcpy(&value,p,sizeof(int),cudaMemcpyDeviceToHost));if(value)throw std::runtime_error("invalid GPU numeric data or CSR indices");}
};
std::atomic<uint64_t> nativeBufferAllocations{0};
struct Buffer {
    double* p=nullptr; size_t count;
    explicit Buffer(size_t n):count(n){ ck(cudaMalloc(&p,std::max<size_t>(1,n)*sizeof(double)));++nativeBufferAllocations; }
    ~Buffer(){if(p)cudaFree(p);}
    Buffer(const Buffer&)=delete; Buffer& operator=(const Buffer&)=delete;
    void zero(){ck(cudaMemset(p,0,count*sizeof(double)));}
};
#include "ams_amg_kernels.cuh"
struct CSR {
    int n,m,nz;const int *rp,*ci;const double* a;
    CSR(const Node& v,Flag& validation):n(v.at("rows").integer()),m(v.at("cols").integer()) {
        if(n<1||m<1||n>INT_MAX-T||m>INT_MAX-T)throw std::runtime_error("operator dimensions exceed supported capacity");
        size_t count=v.at("values").count;
        if(count>INT_MAX-T)throw std::runtime_error("CSR index capacity exceeded");nz=int(count);
        rp=v.at("rowptr").data<int>(Node::integers,size_t(n)+1);
        ci=v.at("columns").data<int>(Node::integers,nz);
        a=v.at("values").data<double>(Node::reals,nz);
        verify_csr<<<(std::max(n,nz)+T-1)/T,T,0,launch_stream>>>(n,m,nz,rp,ci,a,validation.p);
    }
    void mult(const double* x,double* y)const{multiply<<<(n+T-1)/T,T,0,launch_stream>>>(n,rp,ci,a,x,y);}
    void add(const double* x,double alpha,double* y)const{multiply_add<<<(n+T-1)/T,T,0,launch_stream>>>(n,rp,ci,a,x,alpha,y);}
    void relax(Buffer& x,const double* b,const double* l1,Buffer& tmp)const{
        multiply_relax<<<(n+T-1)/T,T,0,launch_stream>>>(n,rp,ci,a,x.p,b,l1,tmp.p);
        std::swap(x.p,tmp.p);
    }
    void res(const double* x,const double* b,double* r)const{
        multiply_residual<<<(n+T-1)/T,T,0,launch_stream>>>(n,rp,ci,a,x,b,r);
    }
};
void finite_vec(Flag& validation,const Node& a,int n,const double*& result,bool positive=false){
    result=a.data<double>(Node::reals,n);
    verify_vec<<<(n+T-1)/T,T,0,launch_stream>>>(n,result,positive,validation.p);
}
struct Coarse {
    std::string kind;int n,used=0;const double *a=nullptr,*di=nullptr;const int *order=nullptr,*rp=nullptr,*ci=nullptr;
    std::unique_ptr<Buffer> work;
    std::unique_ptr<Flag> seen;
    std::unique_ptr<Buffer> inverse;
    std::unique_ptr<Flag> reverse;
    Coarse(const Node& v,Flag& validation,int rows):kind(v.at("kind").string()),n(rows){
        if(v.at("rows").integer()!=n)throw std::runtime_error("coarse size mismatch");
        if(kind=="dense_inverse") {
            a=v.at("values").data<double>(Node::reals,size_t(n)*n);
            if(size_t(n)*n>INT_MAX-T)throw std::runtime_error("dense coarse capacity exceeded");
            verify_vec<<<(n*n+T-1)/T,T,0,launch_stream>>>(n*n,a,false,validation.p);
        } else if(kind=="ngsolve_sparse_ldlt") {
            size_t count=v.at("diagonal_inverse").count;
            if(count>size_t(n))throw std::runtime_error("invalid used coarse rows");used=int(count);
            di=v.at("diagonal_inverse").data<double>(Node::reals,used);
            order=v.at("order").data<int>(Node::integers,n);
            count=v.at("values").count;
            if(count>INT_MAX-T)throw std::runtime_error("coarse capacity exceeded");
            rp=v.at("rowptr").data<int>(Node::integers,size_t(used)+1);
            ci=v.at("columns").data<int>(Node::integers,count);
            a=v.at("values").data<double>(Node::reals,count);
            seen=std::make_unique<Flag>(used);
            verify_csr<<<std::max(1,(std::max(used,int(count))+T-1)/T),T,0,launch_stream>>>(used,used,int(count),rp,ci,a,validation.p);
            verify_factor<<<(n+T-1)/T,T,0,launch_stream>>>(n,used,int(count),order,rp,ci,seen->p,validation.p);
            verify_permutation<<<std::max(1,(used+T-1)/T),T,0,launch_stream>>>(used,seen->p,validation.p);
            verify_vec<<<std::max(1,(used+T-1)/T),T,0,launch_stream>>>(used,di,false,validation.p);work=std::make_unique<Buffer>(used);
            // Validate every factor/index before construction dereferences it.
            validation.check();
            if(used>0&&used<=512){size_t available,total;ck(cudaMemGetInfo(&available,&total));size_t bytes=size_t(used)*used*sizeof(double)+size_t(used)*sizeof(int);
                size_t reserve=std::max<size_t>(2ULL<<30,total/5);if(available<=reserve||bytes>available-reserve)throw std::runtime_error("coarse action cache exceeds CUDA memory reserve");
                inverse=std::make_unique<Buffer>(size_t(used)*used);reverse=std::make_unique<Flag>(used);
                cache_inverse_order<<<(n+T-1)/T,T>>>(n,order,reverse->p);
                inverse_sparse_coarse<<<used,32>>>(used,rp,ci,a,di,inverse->p);
                verify_vec<<<(used*used+T-1)/T,T>>>(used*used,inverse->p,false,validation.p);validation.check();
            }
        } else if(kind!="ten_l1_jacobi_sweeps")throw std::runtime_error("unknown coarse solver");
    }
    void apply(const double* b,double* x)const{
        if(kind=="dense_inverse")dense_coarse<<<(n+T-1)/T,T,0,launch_stream>>>(n,a,b,x);
        else if(inverse)cached_sparse_coarse<<<(n+T-1)/T,T,0,launch_stream>>>(n,used,order,reverse->p,inverse->p,b,x);
        else sparse_coarse<<<1,32,0,launch_stream>>>(n,used,order,rp,ci,a,di,b,x,work->p);
    }
};
struct Level {
    CSR a;std::unique_ptr<CSR> p,pt;const double* l1;int smooth;
    Buffer x,b,r,tmp;std::unique_ptr<Coarse> coarse;
    Level(const Node& v,Flag& validation,bool last):a(v.at("A"),validation),smooth(v.at("num_smooth").integer()),x(a.n),b(a.n),r(a.n),tmp(a.n){
        if(a.n!=a.m||smooth<1)throw std::runtime_error("invalid AMG level");finite_vec(validation,v.at("l1"),a.n,l1,true);
        if(last)coarse=std::make_unique<Coarse>(v.at("coarse"),validation,a.n);
        else {p=std::make_unique<CSR>(v.at("P"),validation);pt=std::make_unique<CSR>(v.at("Pt"),validation);}
    }
    void relax(){a.relax(x,b.p,l1,tmp);}
    void initial(){initial_relax<<<(a.n+T-1)/T,T,0,launch_stream>>>(a.n,b.p,l1,x.p);}
};
struct Space {
    CSR p,pt;std::vector<std::unique_ptr<Level>> levels;Buffer rhs,out;
    Space(const Node& v,Flag& validation):p(v.at("P"),validation),pt(v.at("Pt"),validation),rhs(p.m),out(p.n){
        const auto& list=v.at("levels").items();
        if(list.empty()||list.size()>30)throw std::runtime_error("AMG levels require depth1..30");
        size_t count=list.size();
        for(size_t k=0;k<count;++k)levels.push_back(std::make_unique<Level>(list[k],validation,k+1==count));
        if(pt.n!=p.m||pt.m!=p.n||levels[0]->a.n!=p.m)throw std::runtime_error("auxiliary dimensions mismatch");
        for(size_t k=0;k+1<count;++k){auto& l=*levels[k];int nc=levels[k+1]->a.n;
            if(l.p->n!=l.a.n||l.p->m!=nc||l.pt->n!=nc||l.pt->m!=l.a.n)throw std::runtime_error("AMG transfer dimensions mismatch");}
    }
    void cycle(size_t k){
        auto& l=*levels[k];
        if(l.coarse){if(l.coarse->kind=="ten_l1_jacobi_sweeps"){l.initial();for(int s=1;s<10;++s)l.relax();}else l.coarse->apply(l.b.p,l.x.p);return;}
        l.initial();for(int s=1;s<l.smooth;++s)l.relax();
        l.a.res(l.x.p,l.b.p,l.r.p);
        l.pt->mult(l.r.p,levels[k+1]->b.p);cycle(k+1);l.p->add(levels[k+1]->x.p,1,l.x.p);
        for(int s=0;s<l.smooth;++s)l.relax();
    }
    void apply(const double* residual){pt.mult(residual,levels[0]->b.p);cycle(0);p.mult(levels[0]->x.p,out.p);}
    void add(const double* residual,double alpha,double* output){pt.mult(residual,levels[0]->b.p);cycle(0);p.add(levels[0]->x.p,alpha,output);}
};
struct AMS {
    CSR a;const double* l1;const int* free;int smooth,cycle,gradient;double omega;
    std::vector<std::unique_ptr<Space>> spaces;Buffer b,x,r,tmp;
    AMS(const Node& v,Flag& validation):a(v.at("A"),validation),smooth(v.at("num_smooth").integer()),cycle(v.at("cycle_type").integer()),
        gradient(0),omega(v.at("correction_weight").real()),b(a.n),x(a.n),r(a.n),tmp(a.n){
        if(v.at("schema").string()!="radia.ams-accelerator.v1"||a.n!=a.m||smooth<1||(cycle!=1&&cycle!=7))throw std::runtime_error("unsupported AMS snapshot");
        gradient=!v.at("beta_zero").flag();
        finite_vec(validation,v.at("l1"),a.n,l1,true);free=v.at("free").data<int>(Node::integers,a.n);
        const auto& list=v.at("spaces").items();if(list.size()!=size_t(3+gradient))throw std::runtime_error("invalid AMS spaces");
        for(size_t k=0;k<list.size();++k){spaces.push_back(std::make_unique<Space>(list[k],validation));
            if(spaces.back()->p.n!=a.n)throw std::runtime_error("AMS interpolation dimension mismatch");}
    }
    void res(){a.res(x.p,b.p,r.p);}
    size_t coarseCacheBytes()const{size_t bytes=0;for(const auto& s:spaces)for(const auto& l:s->levels)if(l->coarse&&l->coarse->inverse)bytes+=l->coarse->inverse->count*sizeof(double)+size_t(l->coarse->used)*sizeof(int);return bytes;}
    void relax(){for(int s=0;s<smooth;++s)a.relax(x,b.p,l1,tmp);}
    void initial(){initial_relax<<<(a.n+T-1)/T,T,0,launch_stream>>>(a.n,b.p,l1,x.p);for(int s=1;s<smooth;++s)a.relax(x,b.p,l1,tmp);}
    void grad(){if(!gradient)return;res();spaces[0]->add(r.p,omega,x.p);}
    void nodal(){res();for(int k=gradient;k<gradient+3;++k)spaces[k]->add(r.p,omega,x.p);}
    void apply(const double* rhs,double* output){ck(cudaMemcpyAsync(b.p,rhs,a.n*sizeof(double),cudaMemcpyDeviceToDevice,launch_stream));mask<<<(a.n+T-1)/T,T,0,launch_stream>>>(a.n,free,b.p);
        if(cycle==1){initial();grad();nodal();grad();relax();}
        else {initial();grad();relax();nodal();relax();grad();relax();}
        mask<<<(a.n+T-1)/T,T,0,launch_stream>>>(a.n,free,x.p);ck(cudaMemcpyAsync(output,x.p,a.n*sizeof(double),cudaMemcpyDeviceToDevice,launch_stream));}
};
// Record the fixed AMS cycle once per solve. All graph buffers remain owned by
// AMS and the PCG workspace until the graph is destroyed. Capture does not run
// kernels; each replay initializes its iterates before consuming them.
struct AMSGraph {
    cudaStream_t stream=nullptr;
    cudaGraph_t graph=nullptr;
    cudaGraphExec_t executable=nullptr;
    AMSGraph(AMS& ams,const double* rhs,double* output) {
        ck(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
        try {
            ck(cudaStreamBeginCapture(stream,cudaStreamCaptureModeThreadLocal));
            launch_stream=stream;
            try { ams.apply(rhs,output); }
            catch(...) { launch_stream=nullptr; cudaStreamEndCapture(stream,&graph); throw; }
            launch_stream=nullptr;
            ck(cudaStreamEndCapture(stream,&graph));
            ck(cudaGraphInstantiate(&executable,graph,nullptr,nullptr,0));
        } catch(...) { release(); throw; }
    }
    void release() noexcept {
        if(executable)cudaGraphExecDestroy(executable);
        if(graph)cudaGraphDestroy(graph);
        if(stream)cudaStreamDestroy(stream);
    }
    ~AMSGraph(){release();}
    AMSGraph(const AMSGraph&)=delete;
    AMSGraph& operator=(const AMSGraph&)=delete;
    void apply(){ck(cudaGraphLaunch(executable,nullptr));}
};
struct Blas {
    cublasHandle_t h;
    Blas(){cb(cublasCreate(&h));}
    ~Blas(){cublasDestroy(h);}
    void dot(int n,const double* a,const double* b,double* device_result){
        cb(cublasSetPointerMode(h,CUBLAS_POINTER_MODE_DEVICE));
        cb(cublasDdot(h,n,a,1,b,1,device_result));
    }
    double norm(int n,const double* a){
        cb(cublasSetPointerMode(h,CUBLAS_POINTER_MODE_HOST));
        double d;cb(cublasDnrm2(h,n,a,1,&d));return d;
    }
};
struct ResidentWork {
    Buffer r,p,z,ap,x,scalars;Blas blas;
    explicit ResidentWork(int n):r(n),p(n),z(n),ap(n),x(n),scalars(3){}
    size_t bytes()const{return (r.count+p.count+z.count+ap.count+x.count+3)*sizeof(double);}
};
struct SolveStats {int iterations=0;double residual=0;};
template<class Preconditioner> SolveStats residentPCG(Preconditioner& ams,const CSR& op,const double* b,ResidentWork& w,
    double tol,int maxiter,AMSGraph* graph) {
    int n=ams.a.n;w.x.zero();auto& r=w.r;auto& p=w.p;auto& z=w.z;auto& ap=w.ap;
    auto& scalars=w.scalars;auto& blas=w.blas;SolveStats stats;
    ck(cudaMemcpyAsync(r.p,b,n*sizeof(double),cudaMemcpyDeviceToDevice,launch_stream));
    double bnorm=blas.norm(n,b);
    if(!std::isfinite(bnorm)||bnorm<0)throw std::runtime_error("invalid RHS norm");
    if(bnorm==0)return stats;
    auto precondition=[&](){if(graph)graph->apply();else ams.apply(r.p,z.p);};
    precondition();ck(cudaMemcpyAsync(p.p,z.p,n*sizeof(double),cudaMemcpyDeviceToDevice,launch_stream));
    blas.dot(n,r.p,z.p,scalars.p);stats.residual=1;
    for(stats.iterations=1;stats.iterations<=maxiter;++stats.iterations){
        op.mult(p.p,ap.p);blas.dot(n,p.p,ap.p,scalars.p+1);
        pcg_update_device<<<(n+T-1)/T,T>>>(n,scalars.p,p.p,ap.p,w.x.p,r.p);
        stats.residual=blas.norm(n,r.p)/bnorm;
        if(!std::isfinite(stats.residual))throw std::runtime_error("nonfinite residual");
        if(stats.residual<=tol){op.res(w.x.p,b,r.p);stats.residual=blas.norm(n,r.p)/bnorm;
            if(stats.residual<=tol)break;
            precondition();ck(cudaMemcpyAsync(p.p,z.p,n*sizeof(double),cudaMemcpyDeviceToDevice,launch_stream));
            blas.dot(n,r.p,z.p,scalars.p);continue;
        }
        precondition();blas.dot(n,r.p,z.p,scalars.p+2);
        combine<<<(n+T-1)/T,T>>>(n,z.p,scalars.p,p.p);
        ck(cudaMemcpyAsync(scalars.p,scalars.p+2,sizeof(double),cudaMemcpyDeviceToDevice,launch_stream));
    }
    op.res(w.x.p,b,r.p);stats.residual=blas.norm(n,r.p)/bnorm;
    if(!std::isfinite(stats.residual)||stats.residual>tol||stats.iterations>maxiter)
        throw std::runtime_error("AMS PCG failed its FP64 true residual gate");
    return stats;
}
