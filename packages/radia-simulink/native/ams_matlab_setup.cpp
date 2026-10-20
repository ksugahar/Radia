/* SPDX-License-Identifier: MPL-2.0 */
// MATLAB-owned AMS host setup, reused from the retired MEX constructor.
// SPARSESOLV_SERIAL_AMS_SETUP suppresses all internally owned TaskManager regions.
#ifndef SPARSESOLV_SERIAL_AMS_SETUP
#error MATLAB AMS must disable all internally owned TaskManager regions
#endif
#include "sparsesolv/preconditioners/compact_ams.hpp"
#include <mex.h>
#include <matrix.h>
#include "ams_setup_gpu.hpp"
#include <climits>
#include <filesystem>
#include <limits>
#include <memory>
#include <string>

namespace {
using Mx = std::unique_ptr<mxArray, decltype(&mxDestroyArray)>;
Mx own(mxArray* a) { if (!a) throw std::runtime_error("MAT array allocation failed"); return Mx(a, mxDestroyArray); }
Mx structure(std::initializer_list<const char*> names) {
    std::vector<const char*> fields(names);
    return own(mxCreateStructMatrix(1,1,int(fields.size()),fields.data()));
}
void put(mxArray* a,const char* key,Mx value) { mxSetField(a,0,key,value.release()); }
Mx number(double x) { return own(mxCreateDoubleScalar(x)); }
Mx text(const char* x) { return own(mxCreateString(x)); }
const mxArray* field(const mxArray* a,const char* key) {
    if(!a||!mxIsStruct(a)||mxGetNumberOfElements(a)!=1)throw std::runtime_error("Scalar input struct required");
    auto* v=mxGetField(a,0,key);if(!v)throw std::runtime_error(std::string("Missing input field: ")+key);return v;
}
double scalar(const mxArray* a) {
    if(!mxIsDouble(a)||mxIsComplex(a)||mxIsSparse(a)||mxGetNumberOfElements(a)!=1||!std::isfinite(mxGetScalar(a)))
        throw std::runtime_error("Finite real double scalar required");
    return mxGetScalar(a);
}
int dimension(const mxArray* a) {
    double x=scalar(a);if(x<1||x>INT_MAX-256||x!=std::floor(x))throw std::runtime_error("Invalid int32 dimension");return int(x);
}
std::string word(const mxArray* a) {
    if(!mxIsChar(a))throw std::runtime_error("Character schema required");
    char* p=mxArrayToString(a);if(!p)throw std::runtime_error("Invalid schema text");
    std::string s(p);mxFree(p);return s;
}
template<class T> const T* values(const mxArray* a,mxClassID type,size_t count) {
    if(!a||mxGetClassID(a)!=type||mxIsComplex(a)||mxIsSparse(a)||mxGetNumberOfElements(a)!=count)
        throw std::runtime_error("Input array class/shape mismatch");
    return static_cast<const T*>(mxGetData(a));
}
template<class T> Mx mxvector(const T* data,size_t n,mxClassID type) {
    auto a=own(mxCreateNumericMatrix(n,1,type,mxREAL));
    auto* out=static_cast<T*>(mxGetData(a.get()));
    ngla::CompactParallelFor(n,[&](size_t i){out[i]=data[i];});return a;
}
template<class T,class Range> Mx ngarray(const Range& source,mxClassID type) {
    auto a=own(mxCreateNumericMatrix(source.Size(),1,type,mxREAL));
    auto* out=static_cast<T*>(mxGetData(a.get()));
    ngla::CompactParallelFor(source.Size(),[&](size_t i){out[i]=source[i];});return a;
}
std::shared_ptr<ngla::SparseMatrix<double>> read_csr(const mxArray* a) {
    int n=dimension(field(a,"rows")),m=dimension(field(a,"cols"));
    auto* val=field(a,"values");size_t nz=mxGetNumberOfElements(val);
    if(nz>INT_MAX-256)throw std::overflow_error("Input CSR exceeds int32 capacity");
    auto* rp=values<int>(field(a,"rowptr"),mxINT32_CLASS,size_t(n)+1);
    auto* ci=values<int>(field(a,"columns"),mxINT32_CLASS,nz);
    auto* av=values<double>(val,mxDOUBLE_CLASS,nz);
    if(rp[0]!=0||rp[n]!=int(nz))throw std::runtime_error("Invalid CSR endpoints");
    ngcore::Array<int> counts(n);
    for(int i=0;i<n;++i) {
        if(rp[i]<0||rp[i+1]<rp[i]||rp[i+1]>int(nz))throw std::runtime_error("Invalid CSR row pointer");
        counts[i]=rp[i+1]-rp[i];
        for(int j=rp[i];j<rp[i+1];++j)
            if(ci[j]<0||ci[j]>=m||(j>rp[i]&&ci[j]<=ci[j-1])||!std::isfinite(av[j]))
                throw std::runtime_error("Invalid CSR column/value");
    }
    auto result=std::make_shared<ngla::SparseMatrix<double>>(counts,m);
    ngla::CompactParallelFor(n,[&](size_t i){
        auto cols=result->GetRowIndices(i);auto vals=result->GetRowValues(i);
        for(int j=0;j<counts[i];++j){cols[j]=ci[rp[i]+j];vals[j]=av[rp[i]+j];}
    });return result;
}
Mx csr(const ngla::SparseMatrix<double>& a) {
    if(a.NZE()>size_t(INT_MAX-256))throw std::overflow_error("Output CSR exceeds int32 capacity");
    auto result=structure({"rows","cols","rowptr","columns","values"});
    auto rows=own(mxCreateNumericMatrix(size_t(a.Height())+1,1,mxINT32_CLASS,mxREAL));
    auto cols=own(mxCreateNumericMatrix(a.NZE(),1,mxINT32_CLASS,mxREAL));
    auto vals=own(mxCreateDoubleMatrix(a.NZE(),1,mxREAL));
    auto* rp=static_cast<int*>(mxGetData(rows.get()));auto* ci=static_cast<int*>(mxGetData(cols.get()));
    auto* av=mxGetPr(vals.get());int k=0;
    for(int i=0;i<a.Height();++i){rp[i]=k;k+=int(a.GetRowIndices(i).Size());}rp[a.Height()]=k;
    ngla::CompactParallelFor(a.Height(),[&](size_t i){
        auto c=a.GetRowIndices(i);auto v=a.GetRowValues(i);
        for(int j=0;j<c.Size();++j){ci[rp[i]+j]=c[j];av[rp[i]+j]=v[j];}
    });
    put(result.get(),"rows",number(a.Height()));put(result.get(),"cols",number(a.Width()));
    put(result.get(),"rowptr",std::move(rows));put(result.get(),"columns",std::move(cols));put(result.get(),"values",std::move(vals));return result;
}
struct HostCSR {
    std::vector<int> rp,ci;std::vector<double> av;int n,m;
    explicit HostCSR(const ngla::SparseMatrix<double>& a):n(a.Height()),m(a.Width()) {
        if(n<1||m<1||n>INT_MAX-256||m>INT_MAX-256||a.NZE()>size_t(INT_MAX-256))
            throw std::overflow_error("CUDA setup CSR capacity");
        rp.resize(size_t(n)+1);ci.resize(a.NZE());av.resize(a.NZE());
        int offset=0;
        for(int i=0;i<n;++i){rp[i]=offset;offset+=int(a.GetRowIndices(i).Size());}rp[n]=offset;
        ngla::CompactParallelFor(n,[&](size_t i){auto cols=a.GetRowIndices(i);auto vals=a.GetRowValues(i);
            for(int k=0;k<cols.Size();++k){ci[rp[i]+k]=cols[k];av[rp[i]+k]=vals[k];}});
    }
    AMSSetupCSR view() const {return {n,m,int(av.size()),rp.data(),ci.data(),av.data()};}
};
struct GPUGeneration;
thread_local GPUGeneration* active_generation=nullptr;
struct GPUGeneration {
    void* device=nullptr;GPUGeneration* previous=nullptr;
    ngla::CompactGalerkinCallback previous_callback=nullptr;
    const ngla::SparseMatrix<double>* identity=nullptr;std::unique_ptr<HostCSR> a;
    explicit GPUGeneration(bool enabled) {
        if(!enabled)return;
        device=AMSSetupGPUCreate();previous=active_generation;
        previous_callback=ngla::compact_galerkin_callback;active_generation=this;
        ngla::compact_galerkin_callback=generate;
    }
    ~GPUGeneration(){if(device){ngla::compact_galerkin_callback=previous_callback;active_generation=previous;AMSSetupGPUDestroy(device);}}
    static std::shared_ptr<ngla::SparseMatrix<double>> generate(const ngla::SparseMatrix<double>& matrix,
        const ngla::SparseMatrix<double>& p,const ngla::SparseMatrix<double>& pt) {
        auto& scope=*active_generation;
        if(scope.identity!=&matrix){scope.a=std::make_unique<HostCSR>(matrix);scope.identity=&matrix;}
        HostCSR prolong(p),transpose(pt);
        auto product=AMSSetupGPURestrict(scope.device,&matrix,scope.a->view(),prolong.view(),transpose.view());
        ngcore::Array<int> counts(pt.Height());
        for(int i=0;i<pt.Height();++i)counts[i]=product.rowptr[i+1]-product.rowptr[i];
        auto out=std::make_shared<ngla::SparseMatrix<double>>(counts,p.Width());
        ngla::CompactParallelFor(pt.Height(),[&](size_t i){auto cols=out->GetRowIndices(i);auto vals=out->GetRowValues(i);
            for(int k=0;k<cols.Size();++k){cols[k]=product.columns[product.rowptr[i]+k];vals[k]=product.values[product.rowptr[i]+k];}});
        return out;
    }
    AMSSetupGPUStats stats() const {return device?AMSSetupGPUStatistics(device):AMSSetupGPUStats{};}
};
Mx coarse(const ngla::BaseMatrix* inverse,const std::vector<double>& dense,int n) {
    if(!dense.empty())throw std::runtime_error("Native host requires sparsecholesky coarse factors");
    auto* factor=dynamic_cast<const ngla::SparseCholesky<double>*>(inverse);
    if(!factor)throw std::runtime_error("Native host requires NGSolve sparsecholesky");
    auto first=factor->GetFirstInRow();auto ri=factor->GetFirstInRowRI();
    auto indices=factor->GetRowIndex2();auto lf=factor->GetLFact();
    if(lf.Size()>size_t(INT_MAX-256))throw std::overflow_error("Coarse factor exceeds int32 capacity");
    std::vector<int> rp(first.Size()),ci(lf.Size());
    for(int i=0;i<first.Size();++i){if(first[i]>size_t(INT_MAX-256))throw std::overflow_error("Coarse row capacity");rp[i]=int(first[i]);}
    for(int i=0;i<factor->GetNUsed();++i)
        for(size_t j=first[i];j<first[i+1];++j)ci[j]=indices[ri[i]+j-first[i]];
    auto result=structure({"kind","rows","rowptr","columns","values","order","diagonal_inverse"});
    put(result.get(),"kind",text("ngsolve_sparse_ldlt"));put(result.get(),"rows",number(n));
    put(result.get(),"rowptr",mxvector(rp.data(),rp.size(),mxINT32_CLASS));
    put(result.get(),"columns",mxvector(ci.data(),ci.size(),mxINT32_CLASS));
    put(result.get(),"values",ngarray<double>(lf,mxDOUBLE_CLASS));
    put(result.get(),"order",ngarray<int>(factor->GetOrder(),mxINT32_CLASS));
    put(result.get(),"diagonal_inverse",ngarray<double>(factor->GetDiag(),mxDOUBLE_CLASS));return result;
}
Mx amg(const ngla::CompactAMG& a) {
    auto result=own(mxCreateCellMatrix(a.NumLevels(),1));
    a.VisitAcceleratorLevels([&](size_t index,const auto& matrix,const auto* p,const auto* pt,
        const auto& l1,const auto* inverse,const auto& dense,int smooth){
        bool last=index+1==size_t(a.NumLevels());
        auto level=last?structure({"A","num_smooth","l1","coarse"}):structure({"A","num_smooth","l1","P","Pt"});
        put(level.get(),"A",csr(matrix));put(level.get(),"num_smooth",number(smooth));
        put(level.get(),"l1",mxvector(l1.data(),l1.size(),mxDOUBLE_CLASS));
        if(last)put(level.get(),"coarse",coarse(inverse,dense,matrix.Height()));
        else {if(!p||!pt)throw std::runtime_error("Missing prolongation");put(level.get(),"P",csr(*p));put(level.get(),"Pt",csr(*pt));}
        mxSetCell(result.get(),index,level.release());
    });return result;
}
Mx snapshot(const ngla::HypreBasedAMS& a) {
    auto out=structure({"schema","A","l1","cycle_type","num_smooth","correction_weight","beta_zero","free","spaces"});
    put(out.get(),"schema",text("radia.ams-accelerator.v1"));put(out.get(),"A",csr(a.GetAbc()));
    auto& l1=a.GetL1Norms();put(out.get(),"l1",mxvector(l1.data(),l1.size(),mxDOUBLE_CLASS));
    put(out.get(),"cycle_type",number(a.GetCycleType()));put(out.get(),"num_smooth",number(a.GetNumSmooth()));
    put(out.get(),"correction_weight",number(a.GetCorrectionWeight()));put(out.get(),"beta_zero",own(mxCreateLogicalScalar(a.GetBetaZero())));
    std::vector<int> free(a.GetNdofHC());auto mask=a.GetFreeDofs();
    for(int i=0;i<a.GetNdofHC();++i)free[i]=!mask||mask->Test(i);
    put(out.get(),"free",mxvector(free.data(),free.size(),mxINT32_CLASS));
    auto spaces=own(mxCreateCellMatrix(4,1));size_t index=0;
    auto append=[&](const auto& p,const auto& pt,const ngla::CompactAMG* auxiliary){
        if(!auxiliary)throw std::runtime_error("Missing auxiliary AMG");
        auto s=structure({"P","Pt","levels"});put(s.get(),"P",csr(p));put(s.get(),"Pt",csr(pt));put(s.get(),"levels",amg(*auxiliary));
        mxSetCell(spaces.get(),index++,s.release());
    };
    append(a.GetGrad(),a.GetGradT(),a.GetBGAsAMG());append(a.GetPix(),a.GetPixT(),a.GetBPixAsAMG());
    append(a.GetPiy(),a.GetPiyT(),a.GetBPiyAsAMG());append(a.GetPiz(),a.GetPizT(),a.GetBPizAsAMG());
    put(out.get(),"spaces",std::move(spaces));return out;
}
}
// This entry is linked into the CUDA gateway, never distributed as CPU MEX.
void RadiaAMSSetup(int nlhs,mxArray** output,int nrhs,const mxArray** input) {
    if((nrhs!=3&&nrhs!=4)||nlhs!=2)throw std::runtime_error("[snapshot,info]=radia_cuda_ams('setup',model,threads[,gpu_generation])");
    if(nrhs==4&&!mxIsLogicalScalar(input[3]))throw std::runtime_error("Logical gpu_generation required");
    const bool gpu_generation=nrhs==4&&mxIsLogicalScalarTrue(input[3]);
    const mxArray* model=input[1];
    if(word(field(model,"schema"))!="radia.ams-native-input.v1")throw std::runtime_error("Unsupported input schema");
    int threads=dimension(input[2]);
    ngla::CompactSetupThreadScope thread_scope(threads);
    ngla::HypreBasedAMS::RequireSerialSetup();
    // No TaskManager is activated. Match the serial CPU oracle's NGSolve
    // configuration and restore the caller's limit after this command.
    struct SerialNGSolve {
        int previous=ngcore::TaskManager::GetMaxThreads();
        SerialNGSolve(){ngcore::TaskManager::SetNumThreads(1);}
        ~SerialNGSolve(){ngcore::TaskManager::SetNumThreads(previous);}
    } serial_ngsolve;
        auto a=read_csr(field(model,"A"));auto g=read_csr(field(model,"G"));
        if(a->Height()!=a->Width()||g->Height()!=a->Height())throw std::runtime_error("Incompatible operator dimensions");
        auto* xyz=field(model,"coordinates");
        auto* points=values<double>(xyz,mxDOUBLE_CLASS,size_t(g->Width())*3);
        if(mxGetM(xyz)!=size_t(g->Width())||mxGetN(xyz)!=3)throw std::runtime_error("Coordinates must be vertices by three");
        std::vector<double> coords[3];for(int j=0;j<3;++j){coords[j].resize(g->Width());for(int i=0;i<g->Width();++i){
            double x=points[size_t(j)*g->Width()+i];if(!std::isfinite(x))throw std::runtime_error("Nonfinite coordinates");coords[j][i]=x;}}
        auto* mask=values<int>(field(model,"free"),mxINT32_CLASS,a->Height());
        auto free=std::make_shared<ngcore::BitArray>(a->Height());free->Clear();
        for(int i=0;i<a->Height();++i){if(mask[i]!=0&&mask[i]!=1)throw std::runtime_error("Invalid free mask");if(mask[i])free->SetBit(i);}
        double cycle_value=scalar(field(model,"cycle"));if(cycle_value!=1&&cycle_value!=7)throw std::runtime_error("Unsupported AMS cycle");
        GPUGeneration gpu(gpu_generation);
        auto start=std::chrono::steady_clock::now();
        ngla::HypreBasedAMS ams(a,g,free,coords[0],coords[1],coords[2],int(cycle_value));
        double setup_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        start=std::chrono::steady_clock::now();auto result=snapshot(ams);
        double export_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        auto diagnostic=structure({"schema","setup_seconds","export_seconds","setup_threads","python_adapter_used","coarse_backend","amg_phase_seconds","matrix_generation_backend","gpu_products","gpu_generation_seconds","gpu_upload_bytes","gpu_download_bytes","gpu_setup_peak_owned_bytes"});
        put(diagnostic.get(),"schema",text("radia.ams-matlab-setup.v1"));put(diagnostic.get(),"setup_seconds",number(setup_s));
        put(diagnostic.get(),"export_seconds",number(export_s));put(diagnostic.get(),"setup_threads",number(threads));
        put(diagnostic.get(),"python_adapter_used",own(mxCreateLogicalScalar(false)));put(diagnostic.get(),"coarse_backend",text("ngsolve_sparsecholesky"));
        std::array<double,7> phases{};
        for(const auto* auxiliary : {ams.GetBGAsAMG(),ams.GetBPixAsAMG(),ams.GetBPiyAsAMG(),ams.GetBPizAsAMG()})
            if(auxiliary)for(size_t i=0;i<phases.size();++i)phases[i]+=auxiliary->SetupPhases()[i];
        put(diagnostic.get(),"amg_phase_seconds",mxvector(phases.data(),phases.size(),mxDOUBLE_CLASS));
    auto stats=gpu.stats();
    put(diagnostic.get(),"matrix_generation_backend",text(gpu_generation?"CUDA deterministic Galerkin":"NGSolve CPU Galerkin"));
    put(diagnostic.get(),"gpu_products",number(stats.products));
    put(diagnostic.get(),"gpu_generation_seconds",number(stats.seconds));
    put(diagnostic.get(),"gpu_upload_bytes",number(stats.upload_bytes));
    put(diagnostic.get(),"gpu_download_bytes",number(stats.download_bytes));
    put(diagnostic.get(),"gpu_setup_peak_owned_bytes",number(stats.peak_owned_bytes));
    output[0]=result.release();output[1]=diagnostic.release();
}
