// Experimental numeric leaf boundary, not full HDiv assembly or solver acceptance.
#include "mex.h"
#include "hmatrix_hacapk_pack.cuh"
#include "hmatrix_entry_oracle.cuh"
#include "hmatrix_curved_entry_oracle.cuh"
#include "hmatrix_aca_callback.cuh"
#include "../../../src/core/gpu/hdiv/hmatrix_material_solver.cuh"
#include "../../../src/core/gpu/hdiv/hmatrix_nonlinear_material.cuh"
#ifdef RADIA_HMATRIX_FULL_HOST
#include "hmatrix_full_assembly.cuh"
#endif
#include <map>
#include <iterator>
#include <memory>
#include <chrono>
#include <string_view>
#include <cstdio>
#include <cstdlib>
using namespace radia_hmatrix_cuda;
namespace {
#ifdef RADIA_HMATRIX_FAULT_INJECTION
std::size_t mxFailureAt=0,mxAttempts=0;
bool mxProbeActive=false;
#endif
void mxAllocationProbe(){
#ifdef RADIA_HMATRIX_FAULT_INJECTION
    if(mxProbeActive){++mxAttempts;if(mxFailureAt&&mxAttempts==mxFailureAt)throw std::bad_alloc();}
#endif
}
template<class Factory,class... Args>mxArray* checkedArray(Factory factory,Args... args){
    mxAllocationProbe();auto* value=factory(args...);if(!value)throw std::bad_alloc();return value;
}
char* checkedText(const mxArray* value){mxAllocationProbe();auto* result=mxArrayToString(value);if(!result)throw std::bad_alloc();return result;}
// Opt-in fault-build trace only. Close every record before native work so an
// external watchdog can recover the last boundary even if CUDA blocks.
void diagnosticTrace(const char* phase) noexcept {
#ifdef RADIA_HMATRIX_FAULT_INJECTION
    const char* path=std::getenv("RADIA_HMATRIX_PROGRESS_PATH");
    if(path&&*path){if(auto* file=std::fopen(path,"a")){
        std::fprintf(file,"%s\n",phase);std::fclose(file);
    }}
#else
    (void)phase;
#endif
}
struct CommandTrace {
    ~CommandTrace(){diagnosticTrace("command.scope_exit");}
};
using ContextMap=HostMap<std::uint64_t,HostUnique<Resident>>;
HostUnique<ContextMap> contexts;
using MaterialContextMap=HostMap<std::uint64_t,HostUnique<radia_hdiv_bh::MappedMaterial>>;
HostUnique<MaterialContextMap> materialContexts;
std::uint64_t generation=std::uint64_t(std::chrono::steady_clock::now().time_since_epoch().count())&0x00007fffffffffffULL;
void cleanup(){materialContexts.reset();contexts.reset();}
std::size_t ownedPublicationArrays=0;
#ifdef RADIA_HMATRIX_FAULT_INJECTION
unsigned publicationFailureStage=0,publicationReachedStage=0;
#endif
void publicationProbe(unsigned stage){
#ifdef RADIA_HMATRIX_FAULT_INJECTION
    publicationReachedStage=stage;if(publicationFailureStage==stage)throw std::bad_alloc();
#else
    (void)stage;
#endif
}
struct PublicationArrayDeleter {
    void operator()(mxArray* value)const noexcept{mxDestroyArray(value);--ownedPublicationArrays;}
};
void publishContext(HostUnique<Resident> context,mxArray* output[]){
    diagnosticTrace("publication.enter");
    if(generation==0x0000ffffffffffffULL)throw std::overflow_error("Handle generation exhausted");
    const auto id=0x484d000000000000ULL|++generation;
    publicationProbe(1);
    auto* value=checkedArray(&mxCreateNumericMatrix,1,1,mxUINT64_CLASS,mxREAL);
    diagnosticTrace("publication.array_created");
    if(!value)throw std::bad_alloc();
    ++ownedPublicationArrays;
    std::unique_ptr<mxArray,PublicationArrayDeleter> array(value);
    *static_cast<std::uint64_t*>(mxGetData(value))=id;
    publicationProbe(2);
    // MSVC maps allocate a sentinel even when empty. Own the whole registry
    // only while handles exist, including its object, sentinel and nodes.
    try{
        if(!contexts)contexts=makeHostUnique<ContextMap>("registry.object",std::less<std::uint64_t>{},ContextMap::allocator_type("registry.node"));
        auto inserted=contexts->emplace(id,std::move(context));
        diagnosticTrace("publication.map_inserted");
        if(!inserted.second)throw std::runtime_error("Handle collision");
        try{publicationProbe(3);}
        catch(...){contexts->erase(inserted.first);throw;}
    }catch(...){if(contexts&&contexts->empty())contexts.reset();throw;}
    // Publish ownership and lock only after all throwing preparation is complete.
    mexLock();mexAtExit(cleanup);output[0]=array.release();--ownedPublicationArrays;
    diagnosticTrace("publication.return");
}
void publishMaterialContext(HostUnique<radia_hdiv_bh::MappedMaterial> context,mxArray* output[]){
    diagnosticTrace("publication.enter");
    if(generation==0x0000ffffffffffffULL)throw std::overflow_error("Handle generation exhausted");
    const auto id=0x4248000000000000ULL|++generation;
    publicationProbe(1);
    auto* value=checkedArray(&mxCreateNumericMatrix,1,1,mxUINT64_CLASS,mxREAL);
    diagnosticTrace("publication.array_created");
    if(!value)throw std::bad_alloc();
    ++ownedPublicationArrays;
    std::unique_ptr<mxArray,PublicationArrayDeleter> array(value);
    *static_cast<std::uint64_t*>(mxGetData(value))=id;
    publicationProbe(2);
    // MSVC maps allocate a sentinel even when empty. Own the whole registry
    // only while handles exist, including its object, sentinel and nodes.
    try{
        if(!materialContexts)materialContexts=makeHostUnique<MaterialContextMap>("bh.registry.object",std::less<std::uint64_t>{},MaterialContextMap::allocator_type("bh.registry.node"));
        auto inserted=materialContexts->emplace(id,std::move(context));
        diagnosticTrace("publication.map_inserted");
        if(!inserted.second)throw std::runtime_error("Handle collision");
        try{publicationProbe(3);}
        catch(...){materialContexts->erase(inserted.first);throw;}
    }catch(...){if(materialContexts&&materialContexts->empty())materialContexts.reset();throw;}
    // Publish ownership and lock only after all throwing preparation is complete.
    mexLock();mexAtExit(cleanup);output[0]=array.release();--ownedPublicationArrays;
    diagnosticTrace("publication.return");
}

#ifdef RADIA_HMATRIX_FAULT_INJECTION
void trappedFrameException(){
    mxArray* message=checkedArray(&mxCreateString,"Injected frame-admission exception");
    if(!message)throw std::bad_alloc();
    mxArray* failure=mexCallMATLABWithTrap(0,nullptr,1,&message,"error");
    mxDestroyArray(message);
    if(!failure)throw std::runtime_error("MATLAB exception injection did not raise");
    mxDestroyArray(failure);
    throw std::runtime_error("Injected trapped MATLAB exception at frame admission");
}
#endif
const mxArray* field(const mxArray* a,std::size_t i,const char* name){
    if(!a||!mxIsStruct(a)||i>=mxGetNumberOfElements(a))throw std::invalid_argument("Expected struct");
    auto* result=mxGetField(a,i,name);if(!result)throw std::invalid_argument("Missing required struct field");return result;
}
template<class Allocator>
std::vector<double,Allocator> numericValues(const mxArray* a,const Allocator& allocator){
    if(!a||!mxIsDouble(a)||mxIsSparse(a)||mxIsComplex(a)||mxGetNumberOfDimensions(a)>2)
        throw std::invalid_argument("Expected full real double array");
    const auto n=mxGetNumberOfElements(a);std::vector<double,Allocator> v(n,0.,allocator);auto* p=mxGetDoubles(a);
    for(std::size_t i=0;i<n;++i){if(!std::isfinite(p[i]))throw std::invalid_argument("Nonfinite input");v[i]=p[i];}return v;
}
HostVector<double> numbers(const mxArray* a){return numericValues(a,HostFaultAllocator<double>("normalization.numeric-input"));}
HostVector<double> hostNumbers(const mxArray* a,const char* site){return numericValues(a,HostFaultAllocator<double>(site));}
double scalarNumber(const mxArray* a){
    if(!a||!mxIsDouble(a)||mxIsSparse(a)||mxIsComplex(a)||mxGetNumberOfDimensions(a)>2||mxGetNumberOfElements(a)!=1)
        throw std::invalid_argument("Expected full real double scalar");
    const double value=*mxGetDoubles(a);
    if(!std::isfinite(value))throw std::invalid_argument("Nonfinite scalar input");
    return value;
}
HostVector<double> rowNumbers(const mxArray* a,std::size_t width){
    auto col=numericValues(a,HostFaultAllocator<double>("normalization.column-input"));if(mxGetN(a)!=width)throw std::invalid_argument("Curved row-major record width mismatch");
    HostVector<double> result(col.size(),0.,HostFaultAllocator<double>("normalization.row-input"));const auto rows=mxGetM(a);for(std::size_t r=0;r<rows;++r)for(std::size_t c=0;c<width;++c)result[r*width+c]=col[r+rows*c];return result;
}
class BoundedText {
    char value[128];
public:
    explicit BoundedText(const char* source){const auto n=std::strlen(source);if(n>=sizeof(value))throw std::invalid_argument("Text option exceeds 127-byte limit");std::memcpy(value,source,n+1);}
    bool operator==(const char* other)const noexcept{return std::strcmp(value,other)==0;}
    bool operator!=(const char* other)const noexcept{return !(*this==other);}
    const char* c_str()const noexcept{return value;}
    std::size_t rfind(const char* needle,std::size_t position)const noexcept{return std::string_view(value).rfind(needle,position);}
};
BoundedText text(const mxArray* a);
CurvedBlock curvedBlock(const mxArray* cfg){
    if(!mxIsStruct(cfg)||mxGetNumberOfElements(cfg)!=1)throw std::invalid_argument("Expected scalar curved configuration");
    const char* const allowed[]={"schema_version","geometry_mode","host_vector_cache","host_moment_cache","profile_assembly","cell_nodes","face_nodes","cell_vertices","face_vertices","host","kind","expo","tet_rule","tri_rule","gl","gw","charge_map"};
    for(int f=0;f<mxGetNumberOfFields(cfg);++f)if(std::none_of(std::begin(allowed),std::end(allowed),[&](const char* name){return std::strcmp(name,mxGetFieldNameByNumber(cfg,f))==0;}))throw std::invalid_argument("Unsupported curved option");
    const auto schema=scalarNumber(field(cfg,0,"schema_version"));if(schema!=1.)throw std::invalid_argument("Unsupported curved schema");
    const auto cells=rowNumbers(field(cfg,0,"cell_nodes"),30),faces=rowNumbers(field(cfg,0,"face_nodes"),18);
    const auto cv=rowNumbers(field(cfg,0,"cell_vertices"),4),fv=rowNumbers(field(cfg,0,"face_vertices"),3);
    const auto hosts=hostNumbers(field(cfg,0,"host"),"normalization.hosts"),kinds=hostNumbers(field(cfg,0,"kind"),"normalization.kinds");
    const auto expo=rowNumbers(field(cfg,0,"expo"),3);
    const auto n=hosts.size();if(!n||n>8192||kinds.size()!=n||expo.size()!=3*n||cv.size()/4!=cells.size()/30||fv.size()/3!=faces.size()/18)throw std::invalid_argument("Curved assembly shape/budget mismatch");
    const auto mode=text(field(cfg,0,"geometry_mode"));if(mode!="affine"&&mode!="curved")throw std::invalid_argument("Unsupported HDiv geometry mode");const bool affine=mode=="affine";
    CurvedBlock block;const auto gl=numbers(field(cfg,0,"gl")),gw=numbers(field(cfg,0,"gw"));block.gl.assign(gl.begin(),gl.end());block.gw.assign(gw.begin(),gw.end());
    if(const auto* cache=mxGetField(cfg,0,"host_vector_cache")){if(!mxIsLogicalScalar(cache))throw std::invalid_argument("Host-vector cache option must be scalar logical");block.hostVectorCache=mxIsLogicalScalarTrue(cache);}
    if(const auto* value=mxGetField(cfg,0,"host_moment_cache")){if(!mxIsLogicalScalar(value))throw std::invalid_argument("Host-moment option must be scalar logical");block.hostMomentCache=mxIsLogicalScalarTrue(value);}
    if(const auto* value=mxGetField(cfg,0,"profile_assembly")){if(!mxIsLogicalScalar(value))throw std::invalid_argument("Assembly profile option must be scalar logical");block.profileAssembly=mxIsLogicalScalarTrue(value);}
    if(block.profileAssembly&&!block.hostVectorCache)throw std::invalid_argument("Assembly profiling requires host-vector cache");
    if(block.hostMomentCache&&!block.hostVectorCache)throw std::invalid_argument("Fused moments require host-vector cache");
    if(block.gl.size()!=8||block.gw.size()!=8)throw std::invalid_argument("Curved assembly requires existing eight-point inner Gauss rule");
    for(int k=0;k<8;++k)if(block.gl[k]<=0||block.gl[k]>=1||block.gw[k]<=0||(k&&block.gl[k]<=block.gl[k-1]))throw std::invalid_argument("Invalid curved inner quadrature");
    auto rule=[&](const char* key,HostVector<ReferencePoint>& output){auto values=rowNumbers(field(cfg,0,key),4);if(values.empty()||values.size()>4*128)throw std::invalid_argument("Curved outer quadrature budget exceeded");
        for(std::size_t i=0;i<values.size();i+=4){ReferencePoint q{};double sum=0.;for(int c=0;c<3;++c){q.xi[c]=values[i+c];sum+=q.xi[c];if(q.xi[c]<0)throw std::invalid_argument("Invalid simplex quadrature");}q.weight=values[i+3];if(sum>1.+1e-14||q.weight<=0)throw std::invalid_argument("Invalid simplex quadrature");output.push_back(q);}};
    rule("tet_rule",block.tet);rule("tri_rule",block.tri);
    for(std::size_t i=0;i<n;++i){if(kinds[i]!=0&&kinds[i]!=1)throw std::invalid_argument("Invalid curved charge kind");CurvedCharge c{};c.kind=int(kinds[i]);const auto width=c.kind?18:30,vertices=c.kind?3:4;
        const auto& geometry=c.kind?faces:cells;const auto& ids=c.kind?fv:cv;
        if(hosts[i]<0||std::floor(hosts[i])!=hosts[i]||hosts[i]>=double(geometry.size()/width))throw std::invalid_argument("Invalid curved charge host");c.host=int(hosts[i]);
        for(int j=0;j<width;++j)c.nodes[j]=geometry[std::size_t(c.host)*width+j];
        for(int j=0;j<vertices;++j){double id=ids[std::size_t(c.host)*vertices+j];if(id<0||id>2147483647.||std::floor(id)!=id)throw std::invalid_argument("Invalid topological vertex");c.vertices[j]=int(id);}
        c.affine=affine;int degree=0;for(int j=0;j<3;++j){double e=expo[3*i+j];if(e<0||e>(affine?3:2)||std::floor(e)!=e)throw std::invalid_argument("Unsupported HDiv exponent");degree+=int(e);c.expo[j]=int(e);}
        if(degree>(c.kind?(affine?3:2):(affine?2:1))||(c.kind&&c.expo[2]))throw std::invalid_argument("Unsupported HDiv charge degree");block.charges.push_back(c);
    }
    return block;
}
std::size_t integer(const mxArray* a,bool positive=false){
    const auto value=scalarNumber(a);
    if(value<double(positive)||value>2147483647.||std::floor(value)!=value)
        throw std::invalid_argument("Expected bounded integer scalar");
    return std::size_t(value);
}
std::size_t cacheBytes(const mxArray* a){
    const auto value=scalarNumber(a);
    if(value<1.||value>double(std::size_t(2)<<30)||std::floor(value)!=value)
        throw std::invalid_argument("Expected bounded cache byte budget");
    return std::size_t(value);
}
struct MxTextDeleter {void operator()(char* p)const noexcept{mxFree(p);}};
BoundedText text(const mxArray* a){
    if(!mxIsChar(a)||mxGetNumberOfElements(a)>=128)throw std::invalid_argument("Expected char option of at most 127 characters");
    std::unique_ptr<char,MxTextDeleter> p(checkedText(a));
    if(!p)throw std::runtime_error("Text conversion failed");
    return BoundedText(p.get()); // The MATLAB conversion buffer releases on success or exception.
}
std::uint64_t token(const mxArray* a){if(!mxIsUint64(a)||mxIsComplex(a)||mxIsSparse(a)||mxGetNumberOfElements(a)!=1)throw std::invalid_argument("Expected checked uint64 handle");auto t=*static_cast<const std::uint64_t*>(mxGetData(a));if((!contexts||!contexts->count(t)))throw std::invalid_argument("Wrong type, stale generation or dead H-matrix handle");return t;}
void arity(int nr,int want,int nl,int outputs){if(nr!=want||nl!=outputs)throw std::invalid_argument("Invalid command arity");}
SparseMap readChargeMap(const mxArray* a,std::size_t rows){
    if(!mxIsDouble(a)||mxIsComplex(a)||mxGetNumberOfDimensions(a)!=2||mxGetM(a)!=rows||!mxGetN(a)||mxGetN(a)>2147483647u)
        throw std::invalid_argument("Charge map must be real double charge-count by DOF-count");
    SparseMap b;b.rows=rows;b.columns=mxGetN(a);b.offsets.assign(rows+1,0);
    using ChargeEntry=std::pair<std::size_t,double>;
    HostVector<HostVector<ChargeEntry>> entries{HostFaultAllocator<HostVector<ChargeEntry>>("normalization.rows")};
    entries.reserve(rows);
    for(std::size_t row=0;row<rows;++row)
        entries.emplace_back(HostFaultAllocator<ChargeEntry>("normalization.entries"));
    auto* values=mxGetDoubles(a);
    if(mxIsSparse(a)){
        auto* offsets=mxGetJc(a);auto* indices=mxGetIr(a);
        for(std::size_t col=0;col<b.columns;++col)for(auto k=offsets[col];k<offsets[col+1];++k){if(indices[k]>=rows)throw std::invalid_argument("Sparse map index out of bounds");entries[indices[k]].emplace_back(col,values[k]);}
    }else{
        for(std::size_t col=0;col<b.columns;++col)for(std::size_t row=0;row<rows;++row){double value=values[row+rows*col];if(value!=0)entries[row].emplace_back(col,value);}
    }
    for(std::size_t row=0;row<rows;++row){for(const auto& entry:entries[row]){b.indices.push_back(entry.first);b.values.push_back(entry.second);}b.offsets[row+1]=b.indices.size();}b.validate();return b;
}
void checkPayload(const mxArray* a,std::size_t& bytes,int depth=0){
    if(!a)return;if(depth>32)throw std::invalid_argument("Configuration nesting limit exceeded");
    if(mxIsStruct(a)){for(std::size_t i=0;i<mxGetNumberOfElements(a);++i)for(int f=0;f<mxGetNumberOfFields(a);++f)checkPayload(mxGetFieldByNumber(a,i,f),bytes,depth+1);return;}
    const auto limit=std::size_t(2)<<30;
    auto add=[&](std::size_t n,std::size_t width){auto count=checkedProduct(n,width);if(count>limit-bytes)throw std::invalid_argument("Payload exceeds 2 GiB: HDF5 frontend required (not connected to this experimental gateway)");bytes+=count;};
    if(mxIsSparse(a)){const auto columns=mxGetN(a);add(mxGetJc(a)[columns],mxGetElementSize(a)+sizeof(mwIndex));add(columns+1,sizeof(mwIndex));}
    else add(mxGetNumberOfElements(a),mxGetElementSize(a));
}
Plan read(const mxArray* cfg,bool normalized=false){
    if(mxGetNumberOfElements(cfg)!=1||integer(field(cfg,0,"schema_version"),true)!=(normalized?2:1))throw std::invalid_argument("Unsupported leaf schema");
    const auto n=integer(field(cfg,0,"dimension"),true);
    if(normalized&&n>4608)throw std::invalid_argument("Normalized dimension budget exceeded");
    auto values=numbers(field(cfg,0,"permutation"));
    if(values.size()!=n)throw std::invalid_argument("Permutation size mismatch");HostVector<std::size_t> permutation(n);
    for(std::size_t i=0;i<n;++i){if(values[i]<1||values[i]>n||std::floor(values[i])!=values[i])throw std::invalid_argument("Permutation is one-based");permutation[i]=std::size_t(values[i])-1;}
    auto* blocks=field(cfg,0,"leaves");if(!mxIsStruct(blocks))throw std::invalid_argument("Expected leaves struct array");
    HostVector<Leaf> leaves;
    for(std::size_t i=0;i<mxGetNumberOfElements(blocks);++i){
        Leaf leaf;leaf.row=integer(field(blocks,i,"row"),true)-1;leaf.column=integer(field(blocks,i,"column"),true)-1;
        leaf.rows=integer(field(blocks,i,"rows"),true);leaf.columns=integer(field(blocks,i,"columns"),true);leaf.rank=integer(field(blocks,i,"rank"));
        const auto a1=numbers(field(blocks,i,"a1")),a2=numbers(field(blocks,i,"a2"));leaf.a1.assign(a1.begin(),a1.end());leaf.a2.assign(a2.begin(),a2.end());leaves.push_back(std::move(leaf));
    }
    // Schema 2 already carries normalized factors. Preserve their values,
    // leaf order, permutation and physical scale through the Python packer.
    // Bypass the legacy cHACApK struct/int conversion, without changing math.
    // The existing fromHacapk adapter also preserves factor values.
    if(normalized){
        if(n>4608||leaves.empty()||leaves.size()>65536)throw std::invalid_argument("Normalized leaf budget exceeded");
        for(int f=0;f<mxGetNumberOfFields(cfg);++f){const std::string_view key=mxGetFieldNameByNumber(cfg,f);
            if(key!="schema_version"&&key!="dimension"&&key!="permutation"&&key!="scale"&&key!="leaves"&&key!="operator_kind"&&key!="charge_map"&&key!="assembly")throw std::invalid_argument("Unsupported normalized leaf option");
        }
        if(text(field(cfg,0,"assembly"))!="normalized_leaves")throw std::invalid_argument("Explicit normalized assembly required");
        std::size_t partials=0;
        for(const auto& leaf:leaves){
            if(!leaf.rows||!leaf.columns||leaf.rows>n||leaf.columns>n)throw std::invalid_argument("Invalid normalized leaf dimensions");
            const auto count=leaf.rows+(leaf.row==leaf.column?0:leaf.columns);
            if(count>(std::size_t(1)<<24)-partials)throw std::invalid_argument("Normalized partial workspace budget exceeded");
            partials+=count;
        }
        requireUpperLeafPartition(n,leaves,HostFaultAllocator<std::uint64_t>("normalized.coverage"));
        return pack(n,leaves,permutation,numbers(field(cfg,0,"scale")));
    }
    if(leaves.size()>2147483647u)throw std::overflow_error("Too many native leaves");
    HostVector<st_cHACApK_leafmtx_t> native(leaves.size());
    HostVector<st_cHACApK_leafmtx_t*> pointers(leaves.size()+1,nullptr);
    for(std::size_t i=0;i<leaves.size();++i){auto& leaf=leaves[i];if(leaf.row>leaf.column)throw std::invalid_argument("Numeric schema packs upper leaves only");
        auto& target=native[i];target.ltmtx=leaf.rank?1:2;target.kt=int(leaf.rank);target.nstrtl=int(leaf.row+1);target.nstrtt=int(leaf.column+1);target.ndl=int(leaf.rows);target.ndt=int(leaf.columns);target.a1=leaf.a1.data();target.a2=leaf.a2.data();pointers[i+1]=&target;
        // Validate lengths before the native adapter reads any pointers.
        if(leaf.a1.size()!=checkedProduct(leaf.columns,leaf.rank?leaf.rank:leaf.rows)||leaf.a2.size()!=(leaf.rank?checkedProduct(leaf.rows,leaf.rank):0))throw std::invalid_argument("Invalid native factor storage");
    }
    HostVector<int> lod(n+1);for(std::size_t i=0;i<n;++i)lod[i+1]=int(permutation[i]+1);
    st_cHACApK_leafmtxp_t matrix{};matrix.nd=int(n);matrix.nlf=int(leaves.size());matrix.st_lf=pointers.data();
    st_cHACApK_lcontrol_t control{};control.lod=lod.data();
    return fromHacapk(matrix,control,numbers(field(cfg,0,"scale")));
}
}

std::pair<radia_hdiv_bh::MaterialInput,radia_hdiv_bh::Law> readMaterialConfig(const mxArray* E,const mxArray* weightArray,const mxArray* cfg){
    const auto nq=mxGetNumberOfElements(weightArray),nd=mxGetN(E);
    if(!nq||nq>65536||!nd||nd>16384||!mxIsSparse(E)||mxGetM(E)!=3*nq||mxGetJc(E)[nd]>4000000)throw std::invalid_argument("Mapped material payload shape or budget exceeded");
        if(mxGetNumberOfElements(cfg)!=1)throw std::invalid_argument("Expected one inverse BH law");const auto* ka=field(cfg,0,"knots");const auto* ca=field(cfg,0,"coefficients");const auto* ba=field(cfg,0,"integral_base");const auto nk=mxGetNumberOfElements(ka);
        if(nk<2||nk>4096||mxGetM(ca)!=nk-1||mxGetN(ca)!=4||mxGetNumberOfElements(ba)!=nk)throw std::invalid_argument("Inverse BH law shape exceeded");
        const auto knots=numbers(ka),coefficients=rowNumbers(ca,4),base=numbers(ba);
        radia_hdiv_bh::Law law{{knots.begin(),knots.end()},{coefficients.begin(),coefficients.end()},{base.begin(),base.end()},scalarNumber(field(cfg,0,"cap")),scalarNumber(field(cfg,0,"barrier_slope")),scalarNumber(field(cfg,0,"last_field")),scalarNumber(field(cfg,0,"initial_reluctivity"))};
        const auto map=readChargeMap(E,3*nq);const auto w=numbers(weightArray);
        radia_hdiv_bh::MaterialInput input{nq,nd,{map.offsets.begin(),map.offsets.end()},{map.indices.begin(),map.indices.end()},{map.values.begin(),map.values.end()},{w.begin(),w.end()}};
    return {std::move(input),std::move(law)};
}
std::uint64_t materialToken(const mxArray* a){
    if(!mxIsUint64(a)||mxIsSparse(a)||mxIsComplex(a)||mxGetNumberOfElements(a)!=1)throw std::invalid_argument("Expected checked material uint64 handle");
    const auto id=*static_cast<const std::uint64_t*>(mxGetData(a));if(!materialContexts||!materialContexts->count(id))throw std::invalid_argument("Wrong type, stale generation or dead material handle");return id;
}
mxArray* publishMaterialEvaluation(const radia_hdiv_bh::MaterialEvaluation& result){
    const char* names[]={"gradient","tangent_action","energy","maximum_magnetization","device_payload_bound_bytes","linearization_epoch"};auto* out=checkedArray(&mxCreateStructMatrix,1,1,6,names);++ownedPublicationArrays;
    std::unique_ptr<mxArray,PublicationArrayDeleter> owner(out);
    const radia_hdiv_bh::BoundaryVector<double>* arrays[]={&result.gradient,&result.tangentAction};for(int i=0;i<2;++i){auto* a=checkedArray(&mxCreateDoubleMatrix,arrays[i]->size(),1,mxREAL);std::copy(arrays[i]->begin(),arrays[i]->end(),mxGetDoubles(a));mxSetField(out,0,names[i],a);}
    mxSetField(out,0,names[2],checkedArray(&mxCreateDoubleScalar,result.energy));mxSetField(out,0,names[3],checkedArray(&mxCreateDoubleScalar,result.maximumMagnetization));mxSetField(out,0,names[4],checkedArray(&mxCreateDoubleScalar,double(result.ownedDeviceBytes)));auto* stamp=checkedArray(&mxCreateNumericMatrix,1,1,mxUINT64_CLASS,mxREAL);*static_cast<std::uint64_t*>(mxGetData(stamp))=result.linearizationEpoch;mxSetField(out,0,names[5],stamp);--ownedPublicationArrays;return owner.release();
}
void dispatchCommand(int nl,mxArray* out[],int nr,const mxArray* in[]){
    diagnosticTrace("command.enter");CommandTrace trace;
    #ifdef RADIA_HMATRIX_FAULT_INJECTION
    mxProbeActive=false;
    #endif
    if(!nr)throw std::invalid_argument("Command required");auto command=text(in[0]);
    #ifdef RADIA_HMATRIX_FAULT_INJECTION
    // Enable before every numerical command, including the BH/Newton branches.
    mxProbeActive=command.rfind("diagnostic_",0)!=0;
    #endif
    if(command=="linearize_bh_material"){
        arity(nr,3,nl,1);auto& material=*materialContexts->at(materialToken(in[1]));if(mxGetNumberOfElements(in[2])!=material.dofs())throw std::invalid_argument("Material state dimension mismatch");
        const auto x=numbers(in[2]);out[0]=publishMaterialEvaluation(material.linearize(x));return;
    }
    if(command=="apply_bh_tangent"){
        arity(nr,4,nl,1);auto& material=*materialContexts->at(materialToken(in[1]));if(mxGetNumberOfElements(in[2])!=material.dofs()||!mxIsUint64(in[3])||mxIsSparse(in[3])||mxIsComplex(in[3])||mxGetNumberOfElements(in[3])!=1)throw std::invalid_argument("Expected direction and uint64 linearization epoch");
        const auto d=numbers(in[2]);const auto result=material.tangentAction({d.begin(),d.end()},*static_cast<const std::uint64_t*>(mxGetData(in[3])));
        out[0]=checkedArray(&mxCreateDoubleMatrix,result.size(),1,mxREAL);std::copy(result.begin(),result.end(),mxGetDoubles(out[0]));return;
    }
    if(command=="solve_nonlinear_material"){
        arity(nr,7,nl,1);auto& op=*contexts->at(token(in[1]));auto& material=*materialContexts->at(materialToken(in[2]));const auto n=op.dimension();
        std::size_t payload=0;for(int i=3;i<7;++i)checkPayload(in[i],payload);
        if(!mxIsSparse(in[3])||mxGetM(in[3])!=n||mxGetN(in[3])!=n||mxGetJc(in[3])[n]>4000000||mxGetNumberOfElements(in[4])!=n||mxGetNumberOfElements(in[5])!=n)throw std::invalid_argument("Newton input dimension or payload exceeded");
        const auto* cfg=in[6];if(mxGetNumberOfElements(cfg)!=1)throw std::invalid_argument("Expected one explicit Newton options struct");
        const auto pre=text(field(cfg,0,"preconditioner"));int mode;if(pre=="reference-diagonal")mode=0;else if(pre=="reference-chebyshev")mode=1;else if(pre=="reference-chebyshev-graph")mode=2;else throw std::invalid_argument("Unsupported explicit Newton reference preconditioner");
        radia_hdiv_python::NewtonOptions options{scalarNumber(field(cfg,0,"nonlinear_tolerance")),scalarNumber(field(cfg,0,"linear_tolerance")),int(integer(field(cfg,0,"maximum_iterations"),true)),int(integer(field(cfg,0,"maximum_linear_iterations"),true)),mode,int(integer(field(cfg,0,"degree"),true)),scalarNumber(field(cfg,0,"lower_fraction"))};
        const auto mass=readChargeMap(in[3],n);const auto rhs=numbers(in[4]),initial=numbers(in[5]);
        const auto result=solveNonlinearMaterial(op,material,mass,rhs,initial,options);
        const char* names[]={"x","relative_residual","energy","newton_iterations","preconditioner_graph_launches","history"};out[0]=checkedArray(&mxCreateStructMatrix,1,1,6,names);
        auto* x=checkedArray(&mxCreateDoubleMatrix,n,1,mxREAL);std::copy(result.x.begin(),result.x.end(),mxGetDoubles(x));mxSetField(out[0],0,names[0],x);
        const double scalars[]={result.residual,result.energy,double(result.history.size()),double(result.graphLaunches)};for(int i=1;i<=4;++i)mxSetField(out[0],0,names[i],checkedArray(&mxCreateDoubleScalar,scalars[i-1]));
        const char* fields[]={"relative_residual","energy","alpha","linear_relative_residual","linear_iterations","backtracks","original_residual_certified"};auto* history=checkedArray(&mxCreateStructMatrix,1,result.history.size(),7,fields);
        // Attach the root before filling it so dispatchChecked owns every partial row.
        mxSetField(out[0],0,names[5],history);
        for(std::size_t i=0;i<result.history.size();++i){const auto& row=result.history[i];const double values[]={row.residual,row.energy,row.alpha,row.linearResidual,double(row.linearIterations),double(row.backtracks),row.residualCertified?1.:0.};for(int f=0;f<7;++f)mxSetField(history,i,fields[f],checkedArray(&mxCreateDoubleScalar,values[f]));}
        return;
    }
    if(command=="create_bh_material"){
        arity(nr,4,nl,1);std::size_t payload=0;for(int i=1;i<4;++i)checkPayload(in[i],payload);
        if(materialContexts&&materialContexts->size()>=2)throw std::invalid_argument("Material context count budget exceeded");
        auto config=readMaterialConfig(in[1],in[2],in[3]);
        publishMaterialContext(makeHostUnique<radia_hdiv_bh::MappedMaterial>("bh.material.owner",config.first,config.second),out);return;
    }
    if(command=="evaluate_bh_material_context"){
        arity(nr,4,nl,1);auto& context=*materialContexts->at(materialToken(in[1]));
        if(mxGetNumberOfElements(in[2])!=context.dofs()||(mxGetNumberOfElements(in[3])&&mxGetNumberOfElements(in[3])!=context.dofs()))throw std::invalid_argument("Material state or direction mismatch");
        const auto x=numbers(in[2]),d=numbers(in[3]);out[0]=publishMaterialEvaluation(context.evaluate(x,d));return;
    }
    if(command=="destroy_bh_material"){
        arity(nr,2,nl,0);materialContexts->erase(materialToken(in[1]));if(materialContexts->empty())materialContexts.reset();mexUnlock();return;
    }
    if(command=="bh_material_stats"){
        arity(nr,2,nl,1);const auto stats=materialContexts->at(materialToken(in[1]))->stats();
        const char* names[]={"device_payload_bound_bytes","setup_upload_bytes","state_upload_bytes","evaluations","epoch"};out[0]=checkedArray(&mxCreateStructMatrix,1,1,5,names);
        const double values[]={double(stats.devicePayloadBound),double(stats.setupUploadBytes),double(stats.stateUploadBytes),double(stats.evaluations),double(stats.epoch)};
        for(int i=0;i<5;++i)mxSetField(out[0],0,names[i],checkedArray(&mxCreateDoubleScalar,values[i]));return;
    }
    if(command=="evaluate_bh_material"){
        arity(nr,6,nl,1);std::size_t payload=0;for(int i=1;i<6;++i)checkPayload(in[i],payload);
        auto config=readMaterialConfig(in[1],in[2],in[5]);if(mxGetNumberOfElements(in[3])!=config.first.dofs||(mxGetNumberOfElements(in[4])&&mxGetNumberOfElements(in[4])!=config.first.dofs))throw std::invalid_argument("Material state or direction mismatch");
        const auto state=numbers(in[3]),direction=numbers(in[4]);out[0]=publishMaterialEvaluation(radia_hdiv_bh::evaluateMappedMaterial(config.first,config.second,{state.begin(),state.end()},{direction.begin(),direction.end()}));return;
    }
    if(command=="inverse_bh_samples"){
        arity(nr,9,nl,1);std::size_t payload=0;for(int i=1;i<9;++i)checkPayload(in[i],payload);
        const auto nk=mxGetNumberOfElements(in[1]),ns=mxGetNumberOfElements(in[4]);
        if(nk<2||nk>4096||mxGetNumberOfElements(in[2])!=4*(nk-1)||mxGetNumberOfElements(in[3])!=nk||ns<1||ns>65536)
            throw std::invalid_argument("Inverse BH payload shape or budget exceeded");
        const auto k=numbers(in[1]),c=numbers(in[2]),b=numbers(in[3]),a=numbers(in[4]);
        radia_hdiv_bh::Law law{{k.begin(),k.end()},{c.begin(),c.end()},{b.begin(),b.end()},scalarNumber(in[5]),scalarNumber(in[6]),scalarNumber(in[7]),scalarNumber(in[8])};
        const auto result=radia_hdiv_bh::evaluateInverseBH(law,{a.begin(),a.end()});
        const char* names[]={"secant","tangent","energy","owned_device_payload_bytes"};
        out[0]=checkedArray(&mxCreateStructMatrix,1,1,4,names);
        const radia_hdiv_bh::BoundaryVector<double>* arrays[]={&result.secant,&result.tangent,&result.energy};
        for(int i=0;i<3;++i){auto* value=checkedArray(&mxCreateDoubleMatrix,a.size(),1,mxREAL);std::copy(arrays[i]->begin(),arrays[i]->end(),mxGetDoubles(value));mxSetField(out[0],0,names[i],value);}
        mxSetField(out[0],0,names[3],checkedArray(&mxCreateDoubleScalar,double(result.ownedDeviceBytes)));return;
    }
    #ifdef RADIA_HMATRIX_FAULT_INJECTION
    if(command=="diagnostic_mx_allocation_failure"){arity(nr,2,nl,0);mxFailureAt=integer(in[1]);mxAttempts=0;return;}
    if(command=="diagnostic_mx_allocation_report"){arity(nr,1,nl,1);const char* names[]={"attempts","fail_at"};out[0]=checkedArray(&mxCreateStructMatrix,1,1,2,names);mxSetField(out[0],0,names[0],checkedArray(&mxCreateDoubleScalar,double(mxAttempts)));mxSetField(out[0],0,names[1],checkedArray(&mxCreateDoubleScalar,double(mxFailureAt)));return;}
    #endif
#ifdef RADIA_HMATRIX_FULL_HOST
#ifdef RADIA_HMATRIX_FAULT_INJECTION
    if(command=="diagnostic_trapped_frame_exception"){
        arity(nr,2,nl,0);
        if(((contexts&&!contexts->empty())||(materialContexts&&!materialContexts->empty())))throw std::invalid_argument("Fault injection requires no live operators");
        const auto enabled=integer(in[1]);if(enabled>1)throw std::invalid_argument("Expected zero or one");
        frameExceptionProbe=enabled?trappedFrameException:nullptr;return;
    }
    if(command=="diagnostic_publication_failure"){
        arity(nr,2,nl,0);if(((contexts&&!contexts->empty())||(materialContexts&&!materialContexts->empty()))||ownedPublicationArrays)throw std::invalid_argument("Publication diagnostics require no live owners");
        const auto stage=integer(in[1]);if(stage>3)throw std::invalid_argument("Expected publication stage zero through three");
        publicationFailureStage=unsigned(stage);publicationReachedStage=0;return;
    }
    if(command=="diagnostic_publication_report"){
        arity(nr,1,nl,1);const char* names[]={"owned_output_arrays","live_handles","stage","reached_stage"};
        out[0]=checkedArray(&mxCreateStructMatrix,1,1,4,names);
        mxSetField(out[0],0,names[0],checkedArray(&mxCreateDoubleScalar,ownedPublicationArrays));
        mxSetField(out[0],0,names[1],checkedArray(&mxCreateDoubleScalar,((contexts?contexts->size():0)+(materialContexts?materialContexts->size():0))));
        mxSetField(out[0],0,names[2],checkedArray(&mxCreateDoubleScalar,publicationFailureStage));
        mxSetField(out[0],0,names[3],checkedArray(&mxCreateDoubleScalar,publicationReachedStage));return;
    }
    if(command=="diagnostic_resource_failure"){
        arity(nr,2,nl,0);if(((contexts&&!contexts->empty())||(materialContexts&&!materialContexts->empty())))throw std::invalid_argument("Resource diagnostics require no live operators");
        configureResourceFault(integer(in[1]));return;
    }
    if(command=="diagnostic_resource_report"){
        arity(nr,1,nl,1);const auto& r=resourceFaultState;
        const char* names[]={"attempts","live_resources","free_errors","fail_at","site"};
        out[0]=checkedArray(&mxCreateStructMatrix,1,1,5,names);
        const double values[]={double(r.attempts),double(r.live),double(r.freeErrors),double(r.failAt)};
        for(int i=0;i<4;++i)mxSetField(out[0],0,names[i],checkedArray(&mxCreateDoubleScalar,values[i]));
        mxSetField(out[0],0,"site",checkedArray(&mxCreateString,r.site));return;
    }
    if(command=="diagnostic_device_allocation_failure"){
        arity(nr,2,nl,0);
        if(((contexts&&!contexts->empty())||(materialContexts&&!materialContexts->empty())))throw std::invalid_argument("Fault injection requires no live operators");
        configureDeviceFault(integer(in[1]));return;
    }
    if(command=="diagnostic_device_allocation_report"){
        arity(nr,1,nl,1);const auto& r=deviceFaultState;
        const char* names[]={"attempts","live_bytes","peak_bytes","fail_at","free_errors"};
        out[0]=checkedArray(&mxCreateStructMatrix,1,1,5,names);
        const double values[]={double(r.attempts),double(r.liveBytes),double(r.peakBytes),double(r.failAt),double(r.freeErrors)};
        for(int i=0;i<5;++i)mxSetField(out[0],0,names[i],checkedArray(&mxCreateDoubleScalar,values[i]));return;
    }
    if(command=="diagnostic_host_allocation_failure"){
        if(nl||nr<2||nr>3)throw std::invalid_argument("Host allocation failure expects target and optional site");
        if(((contexts&&!contexts->empty())||(materialContexts&&!materialContexts->empty())))throw std::invalid_argument("Fault injection requires no live operators");
        const auto target=integer(in[1]);
        if(nr==3){const auto site=text(in[2]);configureHostFaultSite(target,site.c_str());}
        else configureHostFault(target);return;
    }
    if(command=="diagnostic_host_allocation_report"){
        arity(nr,1,nl,1);const auto& r=hostFaultState;
        const char* names[]={"attempts","live_bytes","peak_bytes","fail_at","site","filter"};
        out[0]=checkedArray(&mxCreateStructMatrix,1,1,6,names);
        const double values[]={double(r.attempts),double(r.liveBytes),double(r.peakBytes),double(r.failAt)};
        for(int i=0;i<4;++i)mxSetField(out[0],0,names[i],checkedArray(&mxCreateDoubleScalar,values[i]));
        mxSetField(out[0],0,"site",checkedArray(&mxCreateString,r.site));mxSetField(out[0],0,"filter",checkedArray(&mxCreateString,r.filter));return;
    }
    if(command=="diagnostic_c_allocation_failure"){
        arity(nr,2,nl,0);
        if(((contexts&&!contexts->empty())||(materialContexts&&!materialContexts->empty())))throw std::invalid_argument("Fault injection requires no live operators");
        rad_hmatrix_c_memory_configure(integer(in[1]));return;
    }
    if(command=="diagnostic_c_allocation_report"){
        arity(nr,1,nl,1);const auto r=rad_hmatrix_c_memory_last();
        const char* names[]={"allocations","frees","peak_bytes","live_bytes","fail_at","failure_code","last_file","last_line","last_kind","transferred_bytes"};
        out[0]=checkedArray(&mxCreateStructMatrix,1,1,10,names);
        const double values[]={double(r.allocations),double(r.frees),double(r.peak_bytes),double(r.live_bytes),double(r.fail_at),double(r.failure_code)};
        for(int i=0;i<6;++i)mxSetField(out[0],0,names[i],checkedArray(&mxCreateDoubleScalar,values[i]));
        mxSetField(out[0],0,"last_file",checkedArray(&mxCreateString,r.last_file?r.last_file:""));
        mxSetField(out[0],0,"last_line",checkedArray(&mxCreateDoubleScalar,r.last_line));
        mxSetField(out[0],0,"last_kind",checkedArray(&mxCreateString,r.last_kind?r.last_kind:""));
        mxSetField(out[0],0,"transferred_bytes",checkedArray(&mxCreateDoubleScalar,r.transferred_bytes));return;
    }
#endif
    if(command=="create_hdiv_hacapk"||command=="create_hdiv_hacapk_cached"||command=="create_hdiv_hacapk_device_dense"){
        const bool deviceDense=command=="create_hdiv_hacapk_device_dense";
        const bool cached=command!="create_hdiv_hacapk";
        if(cached&&nr==7)arity(nr,7,nl,1);else arity(nr,6,nl,1);
        const auto cacheBudget=nr==7?cacheBytes(in[6]):(std::size_t(256)<<20);
        if(!cacheBudget||cacheBudget>(std::size_t(2)<<30))throw std::invalid_argument("Invalid indexed cache budget");
        std::size_t bytes=0;for(int i=1;i<nr;++i)checkPayload(in[i],bytes);
        auto block=curvedBlock(in[1]);
        if(block.hostVectorCache||block.hostMomentCache||block.profileAssembly)
            throw std::invalid_argument("Curved HACApK entry oracle does not support cache/profile options");
        const auto epsilon=numbers(in[2]),eta=numbers(in[4]);
        if(epsilon.size()!=1||eta.size()!=1)throw std::invalid_argument("Expected scalar HACApK tolerance/admissibility");
        auto map=readChargeMap(field(in[1],0,"charge_map"),block.charges.size());
        auto assembled=assembleFullCurved(block,epsilon[0],int(integer(in[3],true)),eta[0],int(integer(in[5],true)),cached,cacheBudget,deviceDense);
        auto c=makeHostUnique<Resident>("resident.object",assembled.plan,&map,deviceDense?assembled.oracleOwned:0);
        if(deviceDense){
            if(!assembled.retainedOracle)throw std::runtime_error("Missing deferred assembly owner");
            c->materializeDeferred(assembled.plan,*assembled.retainedOracle);
            c->setupPeakOwnedBytes=c->ownedBytes+assembled.oracleOwned;
            assembled.cacheTileLoads=assembled.retainedOracle->cacheTileLoads;assembled.cacheTileHits=assembled.retainedOracle->cacheTileHits;
            assembled.retainedOracle.reset();
        }else c->setupPeakOwnedBytes=std::max(c->ownedBytes,assembled.oracleOwned);
        c->assemblyEntryBatches=assembled.entryBatches;c->assemblyEntryCount=assembled.entryCount;
        c->assemblyLeafCount=assembled.leafCount;c->assemblyLowrankCount=assembled.lowrankCount;
        c->assemblyFrameLeaves=assembled.frameBounds.leaves;c->assemblyFrameDenseCoefficients=assembled.frameBounds.upperDenseCoefficients;c->assemblyFrameStorageCapCoefficients=assembled.frameBounds.upperStorageCapCoefficients;
        c->assemblyCacheTileLoads=assembled.cacheTileLoads;c->assemblyCacheTileHits=assembled.cacheTileHits;
        c->assemblyCacheFullTileBatches=assembled.cacheFullTileBatches;c->assemblyCacheHostScratchBytes=assembled.cacheHostScratchBytes;
        publishContext(std::move(c),out);return;
    }
    if(command=="create_analytic_hacapk"){
        arity(nr,2,nl,1);if(!mxIsStruct(in[1])||mxGetNumberOfElements(in[1])!=1)throw std::invalid_argument("Expected scalar HACApK configuration");std::size_t bytes=0;checkPayload(in[1],bytes);
        if(integer(field(in[1],0,"schema_version"),true)!=1)throw std::invalid_argument("Unsupported HACApK schema");
        for(int f=0;f<mxGetNumberOfFields(in[1]);++f){std::string_view key=mxGetFieldNameByNumber(in[1],f);
            if(key!="schema_version"&&key!="cells"&&key!="faces"&&key!="operator_kind"&&key!="charge_map"&&key!="epsilon"&&key!="leaf_size"&&key!="eta"&&key!="max_rank"&&key!="frame_coefficient_budget_bytes"&&key!="frame_working_payload_budget_bytes")throw std::invalid_argument("Unsupported HACApK option");}
        auto epsilon=numbers(field(in[1],0,"epsilon")),eta=numbers(field(in[1],0,"eta"));
        if(epsilon.size()!=1||eta.size()!=1)throw std::invalid_argument("Expected scalar HACApK tolerance/admissibility");
        if(mxGetN(field(in[1],0,"cells"))>1||mxGetN(field(in[1],0,"faces"))>1)throw std::invalid_argument("Geometry uses flat coordinate column vectors");
        auto block=analyticSources(numbers(field(in[1],0,"cells")),numbers(field(in[1],0,"faces")),4096);
        const auto mode=text(field(in[1],0,"operator_kind"));SparseMap map;
        if(mode=="demag")map=readChargeMap(field(in[1],0,"charge_map"),block.sources.size());else if(mode!="gram")throw std::invalid_argument("Operator kind must be gram or demag");
        const auto* budgetField=mxGetField(in[1],0,"frame_coefficient_budget_bytes");
        const auto frameBudget=budgetField?cacheBytes(budgetField):(std::size_t(2)<<30);
        const auto* workingField=mxGetField(in[1],0,"frame_working_payload_budget_bytes");
        const auto workingBudget=workingField?cacheBytes(workingField):(std::size_t(2)<<30);
        auto assembled=assembleFull(block,epsilon[0],int(integer(field(in[1],0,"leaf_size"),true)),eta[0],int(integer(field(in[1],0,"max_rank"),true)),frameBudget,workingBudget);
        auto c=makeHostUnique<Resident>("resident.object",assembled.plan,mode=="demag"?&map:nullptr);
        c->setupPeakOwnedBytes=std::max(c->ownedBytes,assembled.oracleOwned);
        c->assemblyEntryBatches=assembled.entryBatches;c->assemblyEntryCount=assembled.entryCount;
        c->assemblyLeafCount=assembled.leafCount;c->assemblyLowrankCount=assembled.lowrankCount;
        c->assemblyFrameLeaves=assembled.frameBounds.leaves;c->assemblyFrameDenseCoefficients=assembled.frameBounds.upperDenseCoefficients;c->assemblyFrameStorageCapCoefficients=assembled.frameBounds.upperStorageCapCoefficients;
        publishContext(std::move(c),out);return;
    }
#endif
    if(command=="solve_linear_material"){
        if((nr!=7&&nr!=8&&nr!=10)||nl!=1)throw std::invalid_argument("Expected material solve and optional preconditioner");auto& op=*contexts->at(token(in[1]));std::size_t bytes=0;for(int i=2;i<7;++i)checkPayload(in[i],bytes);
        auto mass=readChargeMap(in[2],op.dimension());auto rhs=numbers(in[3]);auto invChi=numbers(in[4]);auto tolerance=numbers(in[5]);
        if(invChi.size()!=1||tolerance.size()!=1)throw std::invalid_argument("Expected scalar material/tolerance");
        const auto preconditioner=nr>=8?text(in[7]):BoundedText("mass-diagonal");
        if(preconditioner!="mass-diagonal"&&preconditioner!="mass-cholesky"&&preconditioner!="mass-polynomial"&&preconditioner!="mass-polynomial-graph"&&preconditioner!="mass-chebyshev"&&preconditioner!="mass-chebyshev-graph")throw std::invalid_argument("Unsupported material preconditioner");
        int chebyshevDegree=8;double chebyshevLowerFraction=.1;
        if(nr==10){
            checkPayload(in[8],bytes);checkPayload(in[9],bytes);
            if(preconditioner!="mass-chebyshev"&&preconditioner!="mass-chebyshev-graph")throw std::invalid_argument("Chebyshev settings require mass-chebyshev");
            const auto degree=integer(in[8],true);auto lower=numbers(in[9]);
            if(degree>64||lower.size()!=1||!std::isfinite(lower[0])||lower[0]<=0||lower[0]>=1)throw std::invalid_argument("Invalid fixed Chebyshev settings");
            chebyshevDegree=int(degree);chebyshevLowerFraction=lower[0];
        }
        auto result=solveLinearMaterial(op,mass,rhs,invChi[0],tolerance[0],int(integer(in[6],true)),preconditioner=="mass-cholesky",preconditioner=="mass-polynomial"||preconditioner=="mass-polynomial-graph",preconditioner=="mass-polynomial-graph"||preconditioner=="mass-chebyshev-graph",preconditioner=="mass-chebyshev"||preconditioner=="mass-chebyshev-graph",chebyshevDegree,chebyshevLowerFraction);
        const char* names[]={"x","iterations","true_relative_residual","solver_owned_device_bytes","preconditioner","preconditioner_graph_launches","chebyshev_degree","chebyshev_lower_fraction"};out[0]=checkedArray(&mxCreateStructMatrix,1,1,8,names);
        auto* x=checkedArray(&mxCreateDoubleMatrix,result.x.size(),1,mxREAL);std::copy(result.x.begin(),result.x.end(),mxGetDoubles(x));mxSetField(out[0],0,"x",x);
        mxSetField(out[0],0,"iterations",checkedArray(&mxCreateDoubleScalar,result.iterations));mxSetField(out[0],0,"true_relative_residual",checkedArray(&mxCreateDoubleScalar,result.relativeResidual));
        mxSetField(out[0],0,"solver_owned_device_bytes",checkedArray(&mxCreateDoubleScalar,result.ownedBytes));mxSetField(out[0],0,"preconditioner",checkedArray(&mxCreateString,preconditioner.c_str()));mxSetField(out[0],0,"preconditioner_graph_launches",checkedArray(&mxCreateDoubleScalar,result.preconditionerGraphLaunches));mxSetField(out[0],0,"chebyshev_degree",checkedArray(&mxCreateDoubleScalar,(preconditioner=="mass-chebyshev"||preconditioner=="mass-chebyshev-graph")?chebyshevDegree:0));mxSetField(out[0],0,"chebyshev_lower_fraction",checkedArray(&mxCreateDoubleScalar,(preconditioner=="mass-chebyshev"||preconditioner=="mass-chebyshev-graph")?chebyshevLowerFraction:0));return;
    }
    if(command=="analytic_aca_compress"||command=="analytic_aca_leaf"){
        arity(nr,6,nl,1);std::size_t payloadBytes=0;
        for(int i=1;i<6;++i)checkPayload(in[i],payloadBytes);
        auto block=analyticSources(numbers(in[1]),numbers(in[2]),2147483647u);
        const auto rows=integer(in[3],true),rankLimit=integer(in[5],true),n=block.sources.size();
        const auto epsilon=numbers(in[4]);
        if(rows>=n||rankLimit>std::min(rows,n-rows)||epsilon.size()!=1||epsilon[0]<=0||epsilon[0]>=1)
            throw std::invalid_argument("Invalid ACA block split/rank/tolerance");
        const auto columns=n-rows;
        // This bounded diagnostic returns factors directly to MATLAB. Reject
        // oversized generated payloads before host allocation (HDF5 pending).
        if(checkedProduct(n,rankLimit)>(std::size_t(2)<<30)/sizeof(double))
            throw std::runtime_error("ACA factor payload exceeds 2 GiB; HDF5 frontend required");
        HostVector<double> a2(checkedProduct(rows,rankLimit)),a1(checkedProduct(columns,rankLimit));
        HostVector<int> lod(n+1);for(std::size_t i=0;i<n;++i)lod[i+1]=int(i+1);
        double parameters[100]={};parameters[61]=1;parameters[62]=double(rankLimit);parameters[64]=1;
        AnalyticEntryOracle oracle(block,EntryBatchCapacity);HostVector<double> scale(n);
        for(std::size_t offset=0;offset<n;offset+=256){
            HostVector<AnalyticEntryIndex> diagonal;
            for(auto i=offset;i<std::min(n,offset+256);++i)diagonal.push_back({i,i});
            const auto values=oracle.evaluate(diagonal,nullptr);
            for(std::size_t i=0;i<values.size();++i){
                if(!std::isfinite(values[i])||values[i]<=0)throw std::runtime_error("Invalid charge energy diagonal");
                scale[offset+i]=std::sqrt(values[i]);
            }
        }
        EntryCallbackScope binding(oracle,nullptr,&scale);
        const int rank=cHACApK_acaplus(a2.data(),a1.data(),parameters,int(rows),int(columns),1,int(rows+1),lod.data(),0,int(rankLimit),epsilon[0],1.,epsilon[0]*1.e-3);
        binding.rethrow();if(rank<1||rank>int(rankLimit))throw std::runtime_error("ACA compression did not return valid factors");
        // Preserve the production CPU rank-budget rule. Exact entries remain
        // CUDA evaluated; this is not a switch to a CPU numerical backend.
        const bool promoted=rank>=int(rankLimit)&&rankLimit<rows&&rankLimit<columns;
        if(command=="analytic_aca_leaf"&&promoted){
            const auto count=checkedProduct(rows,columns);
            if(count>(std::size_t(2)<<30)/sizeof(double))
                throw std::runtime_error("Dense leaf exceeds 2 GiB; HDF5 frontend required");
            // Release discarded factors before allocating the exact dense leaf.
            HostVector<double>().swap(a1);HostVector<double>().swap(a2);
            std::unique_ptr<mxArray,decltype(&mxDestroyArray)> dense(checkedArray(&mxCreateDoubleMatrix,rows,columns,mxREAL),&mxDestroyArray);
            auto* values=mxGetDoubles(dense.get());
            for(std::size_t offset=0;offset<count;offset+=256){
                HostVector<AnalyticEntryIndex> indices;
                for(auto k=offset;k<std::min(count,offset+256);++k)
                    indices.push_back({k%rows,rows+k/rows});
                const auto batch=oracle.evaluate(indices,nullptr);
                std::copy(batch.begin(),batch.end(),values+offset);
            }
            const char* names[]={"rank","dense","promoted","aca_rank","entry_batches","entry_count","owned_device_bytes"};
            out[0]=checkedArray(&mxCreateStructMatrix,1,1,7,names);
            mxSetField(out[0],0,"rank",checkedArray(&mxCreateDoubleScalar,0));
            mxSetField(out[0],0,"dense",dense.release());
            mxSetField(out[0],0,"promoted",checkedArray(&mxCreateLogicalScalar,true));
            mxSetField(out[0],0,"aca_rank",checkedArray(&mxCreateDoubleScalar,rank));
            mxSetField(out[0],0,"entry_batches",checkedArray(&mxCreateDoubleScalar,oracle.batches));
            mxSetField(out[0],0,"entry_count",checkedArray(&mxCreateDoubleScalar,oracle.entries));
            mxSetField(out[0],0,"owned_device_bytes",checkedArray(&mxCreateDoubleScalar,oracle.ownedBytes));return;
        }
        // Return physical factors while keeping ACA's existing normalized entry
        // semantics and standard threshold multiplier during compression.
        for(int k=0;k<rank;++k){
            for(std::size_t i=0;i<rows;++i)a2[i+rows*k]*=scale[i];
            for(std::size_t i=0;i<columns;++i)a1[i+columns*k]*=scale[rows+i];
        }
        const char* names[]={"rank","a1","a2","entry_batches","entry_count","owned_device_bytes"};
        out[0]=checkedArray(&mxCreateStructMatrix,1,1,6,names);mxSetField(out[0],0,"rank",checkedArray(&mxCreateDoubleScalar,rank));
        std::unique_ptr<mxArray,decltype(&mxDestroyArray)> first(checkedArray(&mxCreateDoubleMatrix,columns,rank,mxREAL),&mxDestroyArray);
        std::unique_ptr<mxArray,decltype(&mxDestroyArray)> second(checkedArray(&mxCreateDoubleMatrix,rows,rank,mxREAL),&mxDestroyArray);
        if(!first||!second)throw std::bad_alloc();
        std::copy(a1.begin(),a1.begin()+columns*rank,mxGetDoubles(first.get()));std::copy(a2.begin(),a2.begin()+rows*rank,mxGetDoubles(second.get()));
        mxSetField(out[0],0,"a1",first.release());mxSetField(out[0],0,"a2",second.release());
        mxSetField(out[0],0,"entry_batches",checkedArray(&mxCreateDoubleScalar,oracle.batches));mxSetField(out[0],0,"entry_count",checkedArray(&mxCreateDoubleScalar,oracle.entries));
        mxSetField(out[0],0,"owned_device_bytes",checkedArray(&mxCreateDoubleScalar,oracle.ownedBytes));return;
    }
    if(command=="analytic_aca_vector"){
        arity(nr,7,nl,1);std::size_t payloadBytes=0;
        for(int i=1;i<7;++i)checkPayload(in[i],payloadBytes);
        auto block=analyticSources(numbers(in[1]),numbers(in[2]),2147483647u);
        const auto n=block.sources.size(),pivot=integer(in[4],true),direction=integer(in[6]);
        const auto permutation=numbers(in[3]),maskValues=numbers(in[5]);
        if(permutation.size()!=n||maskValues.size()!=n||pivot>n||direction>1)
            throw std::invalid_argument("Invalid ACA vector shape/pivot/direction");
        HostVector<int> lod(n+1),mask(n);HostVector<bool> seen(n,false);
        for(std::size_t i=0;i<n;++i){
            const auto value=permutation[i];
            if(value<1||value>n||std::floor(value)!=value||seen[std::size_t(value)-1])
                throw std::invalid_argument("Invalid ACA permutation");
            seen[std::size_t(value)-1]=true;lod[i+1]=int(value);
            if(maskValues[i]!=0&&maskValues[i]!=1)throw std::invalid_argument("Invalid ACA mask");
            mask[i]=int(maskValues[i]);
        }
        AnalyticEntryOracle oracle(block,EntryBatchCapacity);EntryCallbackScope binding(oracle,nullptr);
        HostVector<double> values(n,0.);
        // HACApK ip is zero-based relative to its one-based block start.
        cHACApK_active_entry_vector(int(n),int(pivot-1),1,1,lod.data(),0,mask.data(),int(direction),values.data());
        binding.rethrow();out[0]=checkedArray(&mxCreateDoubleMatrix,n,1,mxREAL);
        std::copy(values.begin(),values.end(),mxGetDoubles(out[0]));return;
    }
    if(command=="analytic_entry_batch"){
        // Native assembly diagnostic boundary; the full ACA adapter owns its
        // oracle across callbacks instead of uploading geometry on every row.
        arity(nr,5,nl,1);std::size_t payloadBytes=0;
        for(int i=1;i<5;++i)checkPayload(in[i],payloadBytes);
        const auto rows=numbers(in[3]),columns=numbers(in[4]);
        if(rows.size()!=columns.size())throw std::invalid_argument("Entry index shape mismatch");
        auto block=analyticSources(numbers(in[1]),numbers(in[2]),2147483647u);
        const std::size_t capacity=256;
        HostVector<AnalyticEntryIndex> request;
        for(std::size_t i=0;i<rows.size();++i){
            if(rows[i]<1||columns[i]<1||rows[i]>block.sources.size()||columns[i]>block.sources.size()||
               std::floor(rows[i])!=rows[i]||std::floor(columns[i])!=columns[i])
                throw std::invalid_argument("Entry indices must be one-based charge indices");
            request.push_back({std::size_t(rows[i])-1,std::size_t(columns[i])-1});
        }
        AnalyticEntryOracle oracle(block,capacity);HostVector<double> result;result.reserve(request.size());
        for(std::size_t offset=0;offset<request.size();offset+=capacity){
            const auto end=std::min(offset+capacity,request.size());
            const HostVector<AnalyticEntryIndex> batch(request.begin()+offset,request.begin()+end);
            const auto values=oracle.evaluate(batch,nullptr);result.insert(result.end(),values.begin(),values.end());
        }
        out[0]=checkedArray(&mxCreateDoubleMatrix,result.size(),1,mxREAL);std::copy(result.begin(),result.end(),mxGetDoubles(out[0]));return;
    }
    if(command=="curved_entry_batch"||command=="curved_entry_batch_cached"){
        const int outputs=command=="curved_entry_batch_cached"&&nl==2?2:1;
        if(command=="curved_entry_batch_cached"&&nr==5)arity(nr,5,nl,outputs);else arity(nr,4,nl,outputs);
        std::size_t payloadBytes=0;for(int i=1;i<nr;++i)checkPayload(in[i],payloadBytes);
        auto block=curvedBlock(in[1]);const auto rows=numbers(in[2]),columns=numbers(in[3]);
        if(rows.size()!=columns.size())throw std::invalid_argument("Curved entry index shape mismatch");
        HostVector<AnalyticEntryIndex> request;request.reserve(rows.size());
        for(std::size_t i=0;i<rows.size();++i){
            if(rows[i]<1||columns[i]<1||rows[i]>block.charges.size()||columns[i]>block.charges.size()||
               std::floor(rows[i])!=rows[i]||std::floor(columns[i])!=columns[i])
                throw std::invalid_argument("Curved entry indices must be one-based charge indices");
            request.push_back({std::size_t(rows[i])-1,std::size_t(columns[i])-1});
        }
        const auto cacheBudget=nr==5?cacheBytes(in[4]):(std::size_t(256)<<20);
        CurvedEntryOracle oracle(block,256,command=="curved_entry_batch_cached",cacheBudget);HostVector<double> result;result.reserve(request.size());
        for(std::size_t offset=0;offset<request.size();offset+=256){const auto end=std::min(offset+256,request.size());
            const HostVector<AnalyticEntryIndex> batch(request.begin()+offset,request.begin()+end);
            const auto values=oracle.evaluate(batch,nullptr);result.insert(result.end(),values.begin(),values.end());}
        out[0]=checkedArray(&mxCreateDoubleMatrix,result.size(),1,mxREAL);std::copy(result.begin(),result.end(),mxGetDoubles(out[0]));
        if(nl==2){const char* names[]={"cache_tile_loads","cache_tile_hits","owned_device_bytes","cache_budget_bytes"};
            out[1]=checkedArray(&mxCreateStructMatrix,1,1,4,names);const std::size_t data[]={oracle.cacheTileLoads,oracle.cacheTileHits,oracle.ownedBytes,cacheBudget};
            for(int i=0;i<4;++i)mxSetField(out[1],0,names[i],checkedArray(&mxCreateDoubleScalar,data[i]));}return;
    }
    if(command=="create_curved_hdiv"||command=="create_hdiv_elements"){
        arity(nr,2,nl,1);std::size_t payloadBytes=0;checkPayload(in[1],payloadBytes);auto block=curvedBlock(in[1]);const auto n=block.charges.size();
        Leaf leaf{0,0,n,n,0,HostVector<double>(checkedProduct(n,n),0.),{}};HostVector<std::size_t> permutation(n);for(std::size_t i=0;i<n;++i)permutation[i]=i;
        auto plan=pack(n,{leaf},permutation,HostVector<double>(n,1.));auto map=readChargeMap(field(in[1],0,"charge_map"),n);
        if(generation==0x0000ffffffffffffULL)throw std::overflow_error("Handle generation exhausted");
        auto c=makeHostUnique<Resident>("resident.object",plan,&map);c->assembleCurved(block);
        publishContext(std::move(c),out);return;
    }
    if(command=="create"||command=="create_normalized"||command=="create_analytic_constant_block"){
        arity(nr,2,nl,1);if(!mxIsStruct(in[1])||mxGetNumberOfElements(in[1])!=1)throw std::invalid_argument("Expected scalar H-matrix configuration");std::size_t payloadBytes=0;checkPayload(in[1],payloadBytes);Plan plan;AnalyticBlock block;
        if(command=="create"||command=="create_normalized")plan=read(in[1],command=="create_normalized");
        else{
            if(const auto* assembly=mxGetField(in[1],0,"assembly")){if(text(assembly)!="analytic_constant_block")throw std::invalid_argument("Unsupported analytic assembly mode");}
            if(integer(field(in[1],0,"schema_version"),true)!=1)throw std::invalid_argument("Unsupported analytic block schema");
            for(int f=0;f<mxGetNumberOfFields(in[1]);++f){std::string_view key=mxGetFieldNameByNumber(in[1],f);
                if(key!="schema_version"&&key!="cells"&&key!="faces"&&key!="operator_kind"&&key!="charge_map"&&key!="assembly")
                    throw std::invalid_argument("Unsupported analytic block option");
            }
            if(mxGetN(field(in[1],0,"cells"))>1||mxGetN(field(in[1],0,"faces"))>1)
                throw std::invalid_argument("Geometry uses flat row-major coordinate column vectors");
            // No images, curved/high-order interpretation or alternate near/far strategy.
            block=analyticSources(numbers(field(in[1],0,"cells")),numbers(field(in[1],0,"faces")));
            const auto n=block.sources.size();Leaf leaf{0,0,n,n,0,HostVector<double>(checkedProduct(n,n),0.),{}};
            HostVector<std::size_t> permutation(n);for(std::size_t i=0;i<n;++i)permutation[i]=i;
            plan=pack(n,{leaf},permutation,HostVector<double>(n,1.));
        }
        if(generation==0x0000ffffffffffffULL)throw std::overflow_error("Handle generation exhausted");
        const auto mode=text(field(in[1],0,"operator_kind"));SparseMap map;
        if(mode=="demag")map=readChargeMap(field(in[1],0,"charge_map"),plan.n);
        else if(mode!="gram")throw std::invalid_argument("Operator kind must be gram or demag");
        auto c=makeHostUnique<Resident>("resident.object",plan,mode=="demag"?&map:nullptr);if(command=="create_analytic_constant_block")c->assemble(block);
        publishContext(std::move(c),out);return;
    }
    if(command=="live_count"){arity(nr,1,nl,1);out[0]=checkedArray(&mxCreateDoubleScalar,((contexts?contexts->size():0)+(materialContexts?materialContexts->size():0)));return;}
    if(command=="destroy"){arity(nr,2,nl,0);diagnosticTrace("destroy.before_erase");contexts->erase(token(in[1]));if(contexts->empty())contexts.reset();diagnosticTrace("destroy.after_erase");mexUnlock();return;}
    if(command=="upload"){arity(nr,3,nl,0);contexts->at(token(in[1]))->upload(numbers(in[2]));return;}
    if(command=="apply"){arity(nr,3,nl,0);contexts->at(token(in[1]))->apply(integer(in[2],true));return;}
    if(command=="download"){
        arity(nr,2,nl,1);auto v=contexts->at(token(in[1]))->download();out[0]=checkedArray(&mxCreateDoubleMatrix,v.size(),1,mxREAL);std::copy(v.begin(),v.end(),mxGetDoubles(out[0]));return;
    }
    if(command=="stats"){
        arity(nr,2,nl,1);const auto& c=*contexts->at(token(in[1]));const char* names[]={"owned_device_bytes","applications","kernel_launches","host_synchronizations_per_apply","setup_peak_owned_device_bytes"};
        out[0]=checkedArray(&mxCreateStructMatrix,1,1,5,names);mxSetField(out[0],0,names[0],checkedArray(&mxCreateDoubleScalar,c.ownedBytes));mxSetField(out[0],0,names[1],checkedArray(&mxCreateDoubleScalar,c.applications));mxSetField(out[0],0,names[2],checkedArray(&mxCreateDoubleScalar,c.kernelsPerApplication*c.applications));mxSetField(out[0],0,names[3],checkedArray(&mxCreateDoubleScalar,0));mxSetField(out[0],0,names[4],checkedArray(&mxCreateDoubleScalar,std::max(c.ownedBytes,c.setupPeakOwnedBytes)));
        const char* extra[]={"assembly_entry_batches","assembly_entry_count","assembly_leaf_count","assembly_lowrank_count"};
        const std::size_t counts[]={c.assemblyEntryBatches,c.assemblyEntryCount,c.assemblyLeafCount,c.assemblyLowrankCount};
        for(int i=0;i<4;++i){mxAddField(out[0],extra[i]);mxSetField(out[0],0,extra[i],checkedArray(&mxCreateDoubleScalar,counts[i]));}
        const char* frameNames[]={"assembly_frame_leaf_count","assembly_frame_upper_dense_coefficients","assembly_frame_storage_cap_coefficients"};
        const std::size_t frameValues[]={c.assemblyFrameLeaves,c.assemblyFrameDenseCoefficients,c.assemblyFrameStorageCapCoefficients};
        for(int i=0;i<3;++i){mxAddField(out[0],frameNames[i]);mxSetField(out[0],0,frameNames[i],checkedArray(&mxCreateDoubleScalar,frameValues[i]));}
        mxAddField(out[0],"assembly_cache_tile_loads");mxSetField(out[0],0,"assembly_cache_tile_loads",checkedArray(&mxCreateDoubleScalar,c.assemblyCacheTileLoads));
        mxAddField(out[0],"assembly_cache_tile_hits");mxSetField(out[0],0,"assembly_cache_tile_hits",checkedArray(&mxCreateDoubleScalar,c.assemblyCacheTileHits));
        mxAddField(out[0],"assembly_cache_full_tile_batches");mxSetField(out[0],0,"assembly_cache_full_tile_batches",checkedArray(&mxCreateDoubleScalar,c.assemblyCacheFullTileBatches));
        mxAddField(out[0],"assembly_cache_host_scratch_bytes");mxSetField(out[0],0,"assembly_cache_host_scratch_bytes",checkedArray(&mxCreateDoubleScalar,c.assemblyCacheHostScratchBytes));
        mxAddField(out[0],"assembly_deferred_dense_leaves");mxSetField(out[0],0,"assembly_deferred_dense_leaves",checkedArray(&mxCreateDoubleScalar,c.deferredDenseLeaves));
        mxAddField(out[0],"assembly_deferred_dense_coefficients");mxSetField(out[0],0,"assembly_deferred_dense_coefficients",checkedArray(&mxCreateDoubleScalar,c.deferredDenseCoefficients));
        mxAddField(out[0],"assembly_coefficient_upload_bytes");mxSetField(out[0],0,"assembly_coefficient_upload_bytes",checkedArray(&mxCreateDoubleScalar,c.coefficientUploadBytes));
        const char* storageNames[]={"stored_coefficient_count","same_leaf_dense_coefficient_count","lowrank_factor_coefficient_count","leaf_rank_sum","max_leaf_rank","storage_saving_lowrank_leaf_count"};
        const std::size_t storageValues[]={c.storedCoefficientCount,c.sameLeafDenseCoefficientCount,c.lowrankFactorCoefficientCount,c.leafRankSum,c.maxLeafRank,c.storageSavingLowrankLeafCount};
        for(int i=0;i<6;++i){mxAddField(out[0],storageNames[i]);mxSetField(out[0],0,storageNames[i],checkedArray(&mxCreateDoubleScalar,storageValues[i]));}
        mxAddField(out[0],"assembly_profile_enabled");mxSetField(out[0],0,"assembly_profile_enabled",checkedArray(&mxCreateLogicalScalar,c.assemblyProfileEnabled));
        const char* stages[]={"profile_mapping_seconds","profile_tet_potential_seconds","profile_tri_potential_seconds","profile_scalar_potential_seconds","profile_contraction_seconds","profile_finalize_seconds"};
        for(int i=0;i<6;++i){mxAddField(out[0],stages[i]);mxSetField(out[0],0,stages[i],checkedArray(&mxCreateDoubleScalar,c.assemblyProfileSeconds[i]));}return;
    }
    throw std::invalid_argument("Unknown H-matrix command");
}
// Complete C++ exception handling in a separate call frame. The entry frame
// below has no native owners or catch scopes for MATLAB's error exception to
// unwind through (including nvcc-generated cleanup metadata).
#ifdef _MSC_VER
__declspec(noinline)
#else
__attribute__((noinline))
#endif
bool dispatchChecked(int nl,mxArray* out[],int nr,const mxArray* in[],char* message,std::size_t capacity) noexcept {
    for(int i=0;i<nl;++i)out[i]=nullptr;
    try{dispatchCommand(nl,out,nr,in);return true;}
    catch(const std::exception& e){diagnosticTrace("error.catch_enter");std::snprintf(message,capacity,"%s",e.what());}
    catch(...){std::snprintf(message,capacity,"Unknown native H-matrix exception");}
    // Assigned output roots own their attached fields. Withdraw all roots on
    // a caught C++ failure before crossing the POD-only MATLAB error gateway.
    for(int i=0;i<nl;++i){if(out[i]){mxDestroyArray(out[i]);out[i]=nullptr;}}
    diagnosticTrace("error.checked_return");return false;
}
void mexFunction(int nl,mxArray* out[],int nr,const mxArray* in[]){
    char errorMessage[2048]={0};
    if(!dispatchChecked(nl,out,nr,in,errorMessage,sizeof(errorMessage))){
        diagnosticTrace("error.before_matlab_exit");
        mexErrMsgIdAndTxt("radia:cuda:hmatrix","%s",errorMessage);
    }
}
