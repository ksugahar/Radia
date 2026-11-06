/* SPDX-License-Identifier: MPL-2.0 */
#pragma once
// Device-resident AMS hierarchy refresh: the device counterpart of
// HypreBasedAMS::Update with reuse_hierarchy (CompactAMG::Refresh). The
// coarsening, the interpolation operators and every sparsity pattern stay as
// built; new operator values on the unchanged fine pattern are turned into the
// BC-modified fine matrix, the auxiliary Galerkin matrices, every coarse level,
// all l1 norms and a new coarsest solve, without any matrix leaving the device.
// Each Galerkin row is accumulated by one owner in the CPU loop order with
// explicitly rounded products and sums (no contraction, no atomics).
// The coarsest level must be an "ngsolve_sparse_ldlt" factor (the default
// hierarchy); it becomes a dense inverse on the rows the CPU sparse Cholesky
// factors (pinned gauge rows and empty rows solve to zero). A lean dense-inverse
// coarsest level, or any state the CPU would rebuild or treat differently,
// raises instead.
// Include inside a private kernel namespace after ams_core.cuh; the
// translation unit includes <cusolverDn.h> and links cuSOLVER.

void cs(cusolverStatus_t s){if(s!=CUSOLVER_STATUS_SUCCESS)throw std::runtime_error("cuSOLVER failed: "+std::to_string(int(s)));}

enum RefreshError : int { refresh_ok=0, refresh_pattern=1, refresh_unsorted=2, refresh_outside=3,
    refresh_nonfinite=4, refresh_coupled_zero=5 };

__global__ void refresh_same_pattern(int n,int nz,const int* rpa,const int* cia,const int* rpb,const int* cib,int* err){
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i<=n&&rpa[i]!=rpb[i])atomicExch(err,refresh_pattern);
    if(i<nz&&cia[i]!=cib[i])atomicExch(err,refresh_pattern);
}
__global__ void refresh_sorted_rows(int n,const int* rp,const int* ci,int* err){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    for(int p=rp[i]+1;p<rp[i+1];++p)if(ci[p]<=ci[p-1]){atomicExch(err,refresh_unsorted);return;}
}
// FillBCModifiedValues: identity rows for constrained dofs, zero couplings to them.
__global__ void refresh_bc_fill(int n,const int* rp,const int* ci,const int* free,const double* src,double* dst,int* err){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    for(int p=rp[i];p<rp[i+1];++p){
        if(!isfinite(src[p]))atomicExch(err,refresh_nonfinite);
        dst[p]=!free[i]?(ci[p]==i?1.0:0.0):(free[ci[p]]?src[p]:0.0);
    }
}
// HypreBasedAMS fine smoother norms.
__global__ void refresh_ams_l1(int n,const int* rp,const int* ci,const double* a,double* l1){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    double diag=0,off=0;
    for(int p=rp[i];p<rp[i+1];++p){if(ci[p]==i)diag=fabs(a[p]);else off=__dadd_rn(off,fabs(a[p]));}
    double temp=__dadd_rn(diag,__dmul_rn(0.5,off));
    l1[i]=temp<=__dmul_rn(4.0/3.0,diag)?(diag>0?diag:1.0):(temp>0?temp:1.0);
}
// GalerkinProductOnPattern: row u of C = sum_e Pt[u,e] sum_k A[e,k] P[k,:] on C's pattern.
__global__ void refresh_galerkin(int rows,const int* trp,const int* tci,const double* tv,
    const int* arp,const int* aci,const double* av,const int* prp,const int* pci,const double* pv,
    const int* crp,const int* cci,double* cv,int* err){
    int u=blockIdx.x*blockDim.x+threadIdx.x;if(u>=rows)return;
    const int begin=crp[u],end=crp[u+1];
    for(int p=begin;p<end;++p)cv[p]=0.0;
    for(int a=trp[u];a<trp[u+1];++a){
        const double pe=tv[a];if(pe==0.0)continue;
        const int e=tci[a];
        for(int b=arp[e];b<arp[e+1];++b){
            const double w=__dmul_rn(pe,av[b]);if(w==0.0)continue;
            const int k=aci[b];
            for(int c=prp[k];c<prp[k+1];++c){
                const int j=pci[c];int lo=begin,hi=end;
                while(lo<hi){int mid=(lo+hi)>>1;if(cci[mid]<j)lo=mid+1;else hi=mid;}
                if(lo==end||cci[lo]!=j){if(pv[c]!=0.0)atomicExch(err,refresh_outside);continue;}
                cv[lo]=__dadd_rn(cv[lo],__dmul_rn(w,pv[c]));
            }
        }
    }
    for(int p=begin;p<end;++p)if(!isfinite(cv[p]))atomicExch(err,refresh_nonfinite);
}
// HypreBasedAMS::FixZeroRows: unit diagonal on rows with zero l1 norm.
__global__ void refresh_fix_zero_rows(int n,const int* rp,const int* ci,double* a){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    double l1=0;int diag=-1;
    for(int p=rp[i];p<rp[i+1];++p){l1=__dadd_rn(l1,fabs(a[p]));if(ci[p]==i)diag=p;}
    if(l1==0.0&&diag>=0)a[diag]=1.0;
}
// CompactAMG::ComputeL1Norms.
__global__ void refresh_amg_l1(int n,const int* rp,const double* a,double* l1){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    double sum=0;for(int p=rp[i];p<rp[i+1];++p)sum=__dadd_rn(sum,fabs(a[p]));
    l1[i]=sum>0?sum:1.0;
}
// Coarsest-level row selection, one thread (coarse levels are small), following
// CompactAMG::SparseCoarsest: with pin_gauge one unknown of each connected
// zero-row-sum component is pinned (seed order); rows without stiffness solve
// to zero; a zero diagonal on a coupled row is rejected.
__global__ void refresh_coarse_rows(int n,const int* rp,const int* ci,const double* a,bool pin_gauge,
    int* keep,int* position,int* queue,int* seen,int* count,int* err){
    if(blockIdx.x||threadIdx.x)return;
    int m=0;
    for(int i=0;i<n;++i){
        double diag=0,off=0;
        for(int p=rp[i];p<rp[i+1];++p){if(!isfinite(a[p]))atomicExch(err,refresh_nonfinite);
            if(ci[p]==i)diag=__dadd_rn(diag,fabs(a[p]));else off=__dadd_rn(off,fabs(a[p]));}
        position[i]=-1;seen[i]=0;
        if(diag>0)keep[m++]=i;
        else if(off>0)atomicExch(err,refresh_coupled_zero);
    }
    if(pin_gauge){
        for(int seed=0;seed<n;++seed){
            if(seen[seed])continue;
            int head=0,tail=0;queue[tail++]=seed;seen[seed]=1;bool constant=true;
            while(head<tail){int i=queue[head++];double sum=0,scale=0;
                for(int p=rp[i];p<rp[i+1];++p){sum=__dadd_rn(sum,a[p]);scale=__dadd_rn(scale,fabs(a[p]));
                    int j=ci[p];if(j!=i&&a[p]!=0&&!seen[j]){seen[j]=1;queue[tail++]=j;}}
                if(fabs(sum)>__dmul_rn(1e-10,scale))constant=false;}
            if(constant)seen[seed]=2;
        }
        int kept=0;for(int k=0;k<m;++k)if(seen[keep[k]]!=2)keep[kept++]=keep[k];m=kept;
    }
    for(int k=0;k<m;++k)position[keep[k]]=k;
    count[0]=m;
}
__global__ void refresh_dense_gather(int m,const int* rp,const int* ci,const double* a,const int* keep,
    const int* position,double* dense){
    int k=blockIdx.x*blockDim.x+threadIdx.x;if(k>=m)return;
    double* column=dense+size_t(k)*m;
    for(int j=0;j<m;++j)column[j]=0;
    int i=keep[k];for(int p=rp[i];p<rp[i+1];++p){int q=position[ci[p]];if(q>=0)column[q]=a[p];}
}
__global__ void refresh_dense_place(int n,int m,const int* keep,const double* inverse,double* out,int* err){
    size_t t=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(t>=size_t(m)*m)return;
    int r=int(t/m),c=int(t%m);
    // potri leaves the lower triangle; the inverse is symmetric.
    double v=r>=c?inverse[size_t(c)*m+r]:inverse[size_t(r)*m+c];
    if(!isfinite(v))atomicExch(err,refresh_nonfinite);
    out[size_t(keep[r])*n+keep[c]]=v;
}

struct RefreshLevel { std::unique_ptr<Buffer> values,l1,inverse; };
struct RefreshSpace { std::unique_ptr<Buffer> auxiliary; std::vector<RefreshLevel> levels; };
struct AMSRefreshStats { double seconds=0; size_t owned_bytes=0; int coarse_rows_max=0; unsigned refreshes=0; };
// Owns every refreshed value buffer; the AMS hierarchy views point into it after
// the first refresh. Buffers are allocated once and rewritten in place.
struct AMSRefresh {
    std::unique_ptr<Buffer> fine,fineL1;std::vector<RefreshSpace> spaces;
    Flag error;cusolverDnHandle_t solver=nullptr;AMSRefreshStats stats;
    AMSRefresh(){cs(cusolverDnCreate(&solver));}
    ~AMSRefresh(){if(solver)cusolverDnDestroy(solver);}
    AMSRefresh(const AMSRefresh&)=delete;AMSRefresh& operator=(const AMSRefresh&)=delete;
    void raise(){
        int value;ck(cudaMemcpy(&value,error.p,sizeof(int),cudaMemcpyDeviceToHost));
        if(!value)return;
        ck(cudaMemset(error.p,0,sizeof(int)));
        const char* reason=value==refresh_pattern?"operator pattern differs from the hierarchy fine pattern":
            value==refresh_unsorted?"refresh requires ascending CSR columns":
            value==refresh_outside?"Galerkin contribution outside the frozen coarse pattern; a full AMS setup is required":
            value==refresh_coupled_zero?"coarsest level has a zero diagonal on a coupled row":
            "nonfinite refreshed AMS values";
        throw std::runtime_error(std::string("AMS device refresh: ")+reason);
    }
};
inline size_t refreshBudget(size_t bytes){
    size_t available,total;ck(cudaMemGetInfo(&available,&total));
    size_t reserve=std::max<size_t>(2ULL<<30,total/5);
    if(available<=reserve||bytes>available-reserve)throw std::runtime_error("AMS device refresh exceeds CUDA memory reserve max(2 GiB,20%)");
    return bytes;
}
inline int blocks(size_t n){return int(std::max<size_t>(1,(n+T-1)/T));}
// Dense inverse of the coarsest level into level.inverse (n x n, zero outside the factored rows).
inline void refreshCoarse(AMSRefresh& s,const CSR& a,Coarse& coarse,RefreshLevel& level,bool pinGauge){
    // After a first refresh the level holds this refresh's dense inverse.
    if(coarse.kind!="ngsolve_sparse_ldlt"&&!(level.inverse&&coarse.a==level.inverse->p))
        throw std::runtime_error("AMS device refresh supports an ngsolve_sparse_ldlt coarsest level only");
    const int n=a.n;
    if(size_t(n)*n>size_t(INT_MAX-T))throw std::runtime_error("AMS device refresh: coarsest level exceeds dense capacity");
    Flag keep(n),position(n),queue(n),seen(n),count(1);
    refresh_coarse_rows<<<1,1>>>(n,a.rp,a.ci,a.a,pinGauge,keep.p,position.p,queue.p,seen.p,count.p,s.error.p);
    ck(cudaGetLastError());s.raise();
    int m;ck(cudaMemcpy(&m,count.p,sizeof(int),cudaMemcpyDeviceToHost));
    if(!level.inverse){refreshBudget(size_t(n)*n*sizeof(double));level.inverse=std::make_unique<Buffer>(size_t(n)*n);s.stats.owned_bytes+=size_t(n)*n*sizeof(double);}
    level.inverse->zero();
    s.stats.coarse_rows_max=std::max(s.stats.coarse_rows_max,n);
    if(m>0){
        refreshBudget(size_t(m)*m*sizeof(double)*2);
        Buffer matrix(size_t(m)*m);
        refresh_dense_gather<<<blocks(m),T>>>(m,a.rp,a.ci,a.a,keep.p,position.p,matrix.p);
        ck(cudaGetLastError());
        int lwork=0;cs(cusolverDnDpotrf_bufferSize(s.solver,CUBLAS_FILL_MODE_LOWER,m,matrix.p,m,&lwork));
        int lworkInverse=0;cs(cusolverDnDpotri_bufferSize(s.solver,CUBLAS_FILL_MODE_LOWER,m,matrix.p,m,&lworkInverse));
        Buffer work(size_t(std::max(lwork,lworkInverse)));Flag info(1);
        cs(cusolverDnDpotrf(s.solver,CUBLAS_FILL_MODE_LOWER,m,matrix.p,m,work.p,lwork,info.p));
        int status;ck(cudaMemcpy(&status,info.p,sizeof(int),cudaMemcpyDeviceToHost));
        if(status!=0)throw std::runtime_error("AMS device refresh: coarsest level is not positive definite on its factored rows");
        cs(cusolverDnDpotri(s.solver,CUBLAS_FILL_MODE_LOWER,m,matrix.p,m,work.p,lworkInverse,info.p));
        ck(cudaMemcpy(&status,info.p,sizeof(int),cudaMemcpyDeviceToHost));
        if(status!=0)throw std::runtime_error("AMS device refresh: coarsest inverse failed");
        refresh_dense_place<<<blocks(size_t(m)*m),T>>>(n,m,keep.p,matrix.p,level.inverse->p,s.error.p);
        ck(cudaGetLastError());s.raise();
    }
    // The coarse action now reads the refreshed dense inverse.
    coarse.kind="dense_inverse";coarse.a=level.inverse->p;coarse.inverse.reset();coarse.reverse.reset();
}
// Refresh `ams` for the values of `op`, which must have the hierarchy fine pattern.
inline AMSRefreshStats refreshAMS(AMS& ams,AMSRefresh& s,const CSR& op){
    auto start=std::chrono::steady_clock::now();
    if(!ams.reuse)throw std::runtime_error("AMS device refresh requires a hierarchy built with reuse_hierarchy");
    const int n=ams.a.n;
    if(op.n!=n||op.m!=n||op.nz!=ams.a.nz)throw std::runtime_error("AMS device refresh: operator dimensions differ from the hierarchy");
    ck(cudaMemset(s.error.p,0,sizeof(int)));
    if(op.rp!=ams.a.rp||op.ci!=ams.a.ci)
        refresh_same_pattern<<<blocks(size_t(std::max(n+1,op.nz))),T>>>(n,op.nz,op.rp,op.ci,ams.a.rp,ams.a.ci,s.error.p);
    refresh_sorted_rows<<<blocks(n),T>>>(n,ams.a.rp,ams.a.ci,s.error.p);
    for(auto& space:ams.spaces)for(auto& level:space->levels)
        refresh_sorted_rows<<<blocks(level->a.n),T>>>(level->a.n,level->a.rp,level->a.ci,s.error.p);
    ck(cudaGetLastError());s.raise();
    const bool first=!s.fine;
    if(first){
        size_t bytes=size_t(ams.a.nz+n)*sizeof(double);
        for(auto& space:ams.spaces)for(auto& level:space->levels)bytes+=size_t(level->a.nz+level->a.n)*sizeof(double);
        refreshBudget(bytes);
        s.fine=std::make_unique<Buffer>(ams.a.nz);s.fineL1=std::make_unique<Buffer>(n);
        s.spaces.resize(ams.spaces.size());
        for(size_t k=0;k<ams.spaces.size();++k){
            auto& levels=ams.spaces[k]->levels;s.spaces[k].levels.resize(levels.size());
            for(size_t l=0;l<levels.size();++l){
                s.spaces[k].levels[l].values=std::make_unique<Buffer>(levels[l]->a.nz);
                s.spaces[k].levels[l].l1=std::make_unique<Buffer>(levels[l]->a.n);
            }
        }
        s.stats.owned_bytes+=bytes;
    }
    // Fine level: BC-modified values and the AMS smoother norms.
    refresh_bc_fill<<<blocks(n),T>>>(n,ams.a.rp,ams.a.ci,ams.free,op.a,s.fine->p,s.error.p);
    refresh_ams_l1<<<blocks(n),T>>>(n,ams.a.rp,ams.a.ci,s.fine->p,s.fineL1->p);
    for(size_t k=0;k<ams.spaces.size();++k){
        auto& space=*ams.spaces[k];auto& owned=s.spaces[k];
        // Auxiliary matrix Pt A_bc P on the frozen level-0 pattern, then FixZeroRows.
        auto& top=space.levels[0]->a;
        refresh_galerkin<<<blocks(top.n),T>>>(top.n,space.pt.rp,space.pt.ci,space.pt.a,ams.a.rp,ams.a.ci,s.fine->p,
            space.p.rp,space.p.ci,space.p.a,top.rp,top.ci,owned.levels[0].values->p,s.error.p);
        refresh_fix_zero_rows<<<blocks(top.n),T>>>(top.n,top.rp,top.ci,owned.levels[0].values->p);
        for(size_t l=0;l<space.levels.size();++l){
            auto& level=*space.levels[l];auto& values=owned.levels[l];
            if(l>0){auto& fine=*space.levels[l-1];
                refresh_galerkin<<<blocks(level.a.n),T>>>(level.a.n,fine.pt->rp,fine.pt->ci,fine.pt->a,fine.a.rp,fine.a.ci,owned.levels[l-1].values->p,
                    fine.p->rp,fine.p->ci,fine.p->a,level.a.rp,level.a.ci,values.values->p,s.error.p);}
            refresh_amg_l1<<<blocks(level.a.n),T>>>(level.a.n,level.a.rp,values.values->p,values.l1->p);
        }
        ck(cudaGetLastError());s.raise();
        auto& last=*space.levels.back();
        // Only the gradient space (first, when present) has the constant kernel.
        if(last.coarse&&last.coarse->kind!="ten_l1_jacobi_sweeps"){
            CSR view=last.a;view.a=owned.levels.back().values->p;
            refreshCoarse(s,view,*last.coarse,owned.levels.back(),ams.gradient&&k==0);
        }
    }
    // Point the hierarchy at the refreshed values (pattern pointers are unchanged).
    ams.a.a=s.fine->p;ams.l1=s.fineL1->p;
    for(size_t k=0;k<ams.spaces.size();++k)
        for(size_t l=0;l<ams.spaces[k]->levels.size();++l){
            auto& level=*ams.spaces[k]->levels[l];
            level.a.a=s.spaces[k].levels[l].values->p;level.l1=s.spaces[k].levels[l].l1->p;
        }
    ck(cudaGetLastError());ck(cudaDeviceSynchronize());s.raise();
    ++s.stats.refreshes;
    s.stats.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return s.stats;
}
