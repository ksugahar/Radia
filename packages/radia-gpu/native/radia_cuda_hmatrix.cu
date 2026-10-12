// Experimental numeric leaf boundary, not full HDiv assembly or solver acceptance.
#include "mex.h"
#include "hmatrix_hacapk_pack.cuh"
#include "hmatrix_entry_oracle.cuh"
#include "hmatrix_aca_callback.cuh"
#include "hmatrix_material_solver.cuh"
#ifdef RADIA_HMATRIX_FULL_HOST
#include "hmatrix_full_assembly.cuh"
#endif
#include <map>
#include <memory>
#include <chrono>
#include <string>
using namespace radia_hmatrix_cuda;
namespace {
std::map<std::uint64_t,std::unique_ptr<Resident>> contexts;
std::uint64_t generation=std::uint64_t(std::chrono::steady_clock::now().time_since_epoch().count())&0x00007fffffffffffULL;
void cleanup(){contexts.clear();}
const mxArray* field(const mxArray* a,std::size_t i,const char* name){
    if(!a||!mxIsStruct(a)||i>=mxGetNumberOfElements(a))throw std::invalid_argument("Expected struct");
    auto* result=mxGetField(a,i,name);if(!result)throw std::invalid_argument(std::string("Missing field: ")+name);return result;
}
std::vector<double> numbers(const mxArray* a){
    if(!a||!mxIsDouble(a)||mxIsSparse(a)||mxIsComplex(a)||mxGetNumberOfDimensions(a)>2)
        throw std::invalid_argument("Expected full real double array");
    const auto n=mxGetNumberOfElements(a);std::vector<double> v(n);auto* p=mxGetDoubles(a);
    for(std::size_t i=0;i<n;++i){if(!std::isfinite(p[i]))throw std::invalid_argument("Nonfinite input");v[i]=p[i];}return v;
}
std::vector<double> rowNumbers(const mxArray* a,std::size_t width){
    auto col=numbers(a);if(mxGetN(a)!=width)throw std::invalid_argument("Curved row-major record width mismatch");
    std::vector<double> result(col.size());const auto rows=mxGetM(a);for(std::size_t r=0;r<rows;++r)for(std::size_t c=0;c<width;++c)result[r*width+c]=col[r+rows*c];return result;
}
CurvedBlock curvedBlock(const mxArray* cfg){
    if(!mxIsStruct(cfg)||mxGetNumberOfElements(cfg)!=1)throw std::invalid_argument("Expected scalar curved configuration");
    const std::vector<std::string> allowed={"schema_version","geometry_mode","cell_nodes","face_nodes","cell_vertices","face_vertices","host","kind","expo","tet_rule","tri_rule","gl","gw","charge_map"};
    for(int f=0;f<mxGetNumberOfFields(cfg);++f)if(std::find(allowed.begin(),allowed.end(),mxGetFieldNameByNumber(cfg,f))==allowed.end())throw std::invalid_argument("Unsupported curved option");
    auto schema=numbers(field(cfg,0,"schema_version"));if(schema.size()!=1||schema[0]!=1.)throw std::invalid_argument("Unsupported curved schema");
    const auto cells=rowNumbers(field(cfg,0,"cell_nodes"),30),faces=rowNumbers(field(cfg,0,"face_nodes"),18);
    const auto cv=rowNumbers(field(cfg,0,"cell_vertices"),4),fv=rowNumbers(field(cfg,0,"face_vertices"),3);
    const auto hosts=numbers(field(cfg,0,"host")),kinds=numbers(field(cfg,0,"kind")),expo=rowNumbers(field(cfg,0,"expo"),3);
    const auto n=hosts.size();if(!n||n>512||kinds.size()!=n||expo.size()!=3*n||cv.size()/4!=cells.size()/30||fv.size()/3!=faces.size()/18)throw std::invalid_argument("Curved assembly shape/budget mismatch");
    const auto mode=text(field(cfg,0,"geometry_mode"));if(mode!="affine"&&mode!="curved")throw std::invalid_argument("Unsupported HDiv geometry mode");const bool affine=mode=="affine";
    CurvedBlock block;block.gl=numbers(field(cfg,0,"gl"));block.gw=numbers(field(cfg,0,"gw"));
    if(block.gl.size()!=8||block.gw.size()!=8)throw std::invalid_argument("Curved assembly requires existing eight-point inner Gauss rule");
    for(int k=0;k<8;++k)if(block.gl[k]<=0||block.gl[k]>=1||block.gw[k]<=0||(k&&block.gl[k]<=block.gl[k-1]))throw std::invalid_argument("Invalid curved inner quadrature");
    auto rule=[&](const char* key,std::vector<ReferencePoint>& output){auto values=rowNumbers(field(cfg,0,key),4);if(values.empty()||values.size()>4*128)throw std::invalid_argument("Curved outer quadrature budget exceeded");
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
    auto v=numbers(a);if(v.size()!=1||v[0]<double(positive)||v[0]>2147483647.||std::floor(v[0])!=v[0])
        throw std::invalid_argument("Expected bounded integer scalar");return std::size_t(v[0]);
}
std::string text(const mxArray* a){if(!mxIsChar(a))throw std::invalid_argument("Command must be char");char* p=mxArrayToString(a);if(!p)throw std::runtime_error("Text conversion failed");std::string s(p);mxFree(p);return s;}
std::uint64_t token(const mxArray* a){if(!mxIsUint64(a)||mxIsComplex(a)||mxIsSparse(a)||mxGetNumberOfElements(a)!=1)throw std::invalid_argument("Expected checked uint64 handle");auto t=*static_cast<const std::uint64_t*>(mxGetData(a));if(!contexts.count(t))throw std::invalid_argument("Wrong type, stale generation or dead H-matrix handle");return t;}
void arity(int nr,int want,int nl,int outputs){if(nr!=want||nl!=outputs)throw std::invalid_argument("Invalid command arity");}
SparseMap readChargeMap(const mxArray* a,std::size_t rows){
    if(!mxIsDouble(a)||mxIsComplex(a)||mxGetNumberOfDimensions(a)!=2||mxGetM(a)!=rows||!mxGetN(a)||mxGetN(a)>2147483647u)
        throw std::invalid_argument("Charge map must be real double charge-count by DOF-count");
    SparseMap b;b.rows=rows;b.columns=mxGetN(a);b.offsets.assign(rows+1,0);
    std::vector<std::vector<std::pair<std::size_t,double>>> entries(rows);
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
Plan read(const mxArray* cfg){
    if(mxGetNumberOfElements(cfg)!=1||integer(field(cfg,0,"schema_version"),true)!=1)throw std::invalid_argument("Unsupported leaf schema");
    const auto n=integer(field(cfg,0,"dimension"),true);auto values=numbers(field(cfg,0,"permutation"));
    if(values.size()!=n)throw std::invalid_argument("Permutation size mismatch");std::vector<std::size_t> permutation(n);
    for(std::size_t i=0;i<n;++i){if(values[i]<1||values[i]>n||std::floor(values[i])!=values[i])throw std::invalid_argument("Permutation is one-based");permutation[i]=std::size_t(values[i])-1;}
    auto* blocks=field(cfg,0,"leaves");if(!mxIsStruct(blocks))throw std::invalid_argument("Expected leaves struct array");
    std::vector<Leaf> leaves;
    for(std::size_t i=0;i<mxGetNumberOfElements(blocks);++i){
        Leaf leaf;leaf.row=integer(field(blocks,i,"row"),true)-1;leaf.column=integer(field(blocks,i,"column"),true)-1;
        leaf.rows=integer(field(blocks,i,"rows"),true);leaf.columns=integer(field(blocks,i,"columns"),true);leaf.rank=integer(field(blocks,i,"rank"));
        leaf.a1=numbers(field(blocks,i,"a1"));leaf.a2=numbers(field(blocks,i,"a2"));leaves.push_back(std::move(leaf));
    }
    if(leaves.size()>2147483647u)throw std::overflow_error("Too many native leaves");
    std::vector<st_cHACApK_leafmtx_t> native(leaves.size());
    std::vector<st_cHACApK_leafmtx_t*> pointers(leaves.size()+1,nullptr);
    for(std::size_t i=0;i<leaves.size();++i){auto& leaf=leaves[i];if(leaf.row>leaf.column)throw std::invalid_argument("Numeric schema packs upper leaves only");
        auto& target=native[i];target.ltmtx=leaf.rank?1:2;target.kt=int(leaf.rank);target.nstrtl=int(leaf.row+1);target.nstrtt=int(leaf.column+1);target.ndl=int(leaf.rows);target.ndt=int(leaf.columns);target.a1=leaf.a1.data();target.a2=leaf.a2.data();pointers[i+1]=&target;
        // Validate lengths before the native adapter reads any pointers.
        if(leaf.a1.size()!=checkedProduct(leaf.columns,leaf.rank?leaf.rank:leaf.rows)||leaf.a2.size()!=(leaf.rank?checkedProduct(leaf.rows,leaf.rank):0))throw std::invalid_argument("Invalid native factor storage");
    }
    std::vector<int> lod(n+1);for(std::size_t i=0;i<n;++i)lod[i+1]=int(permutation[i]+1);
    st_cHACApK_leafmtxp_t matrix{};matrix.nd=int(n);matrix.nlf=int(leaves.size());matrix.st_lf=pointers.data();
    st_cHACApK_lcontrol_t control{};control.lod=lod.data();
    return fromHacapk(matrix,control,numbers(field(cfg,0,"scale")));
}
}
void mexFunction(int nl,mxArray* out[],int nr,const mxArray* in[]){
try{
    if(!nr)throw std::invalid_argument("Command required");auto command=text(in[0]);
#ifdef RADIA_HMATRIX_FULL_HOST
    if(command=="create_analytic_hacapk"){
        arity(nr,2,nl,1);if(!mxIsStruct(in[1])||mxGetNumberOfElements(in[1])!=1)throw std::invalid_argument("Expected scalar HACApK configuration");std::size_t bytes=0;checkPayload(in[1],bytes);
        if(integer(field(in[1],0,"schema_version"),true)!=1)throw std::invalid_argument("Unsupported HACApK schema");
        for(int f=0;f<mxGetNumberOfFields(in[1]);++f){std::string key=mxGetFieldNameByNumber(in[1],f);
            if(key!="schema_version"&&key!="cells"&&key!="faces"&&key!="operator_kind"&&key!="charge_map"&&key!="epsilon"&&key!="leaf_size"&&key!="eta"&&key!="max_rank")throw std::invalid_argument("Unsupported HACApK option: "+key);}
        auto epsilon=numbers(field(in[1],0,"epsilon")),eta=numbers(field(in[1],0,"eta"));
        if(epsilon.size()!=1||eta.size()!=1)throw std::invalid_argument("Expected scalar HACApK tolerance/admissibility");
        if(mxGetN(field(in[1],0,"cells"))>1||mxGetN(field(in[1],0,"faces"))>1)throw std::invalid_argument("Geometry uses flat coordinate column vectors");
        auto block=analyticSources(numbers(field(in[1],0,"cells")),numbers(field(in[1],0,"faces")),4096);
        const auto mode=text(field(in[1],0,"operator_kind"));SparseMap map;
        if(mode=="demag")map=readChargeMap(field(in[1],0,"charge_map"),block.sources.size());else if(mode!="gram")throw std::invalid_argument("Operator kind must be gram or demag");
        auto assembled=assembleFull(block,epsilon[0],int(integer(field(in[1],0,"leaf_size"),true)),eta[0],int(integer(field(in[1],0,"max_rank"),true)));
        auto c=std::make_unique<Resident>(assembled.plan,mode=="demag"?&map:nullptr);
        c->setupPeakOwnedBytes=std::max(c->ownedBytes,assembled.oracleOwned);
        c->assemblyEntryBatches=assembled.entryBatches;c->assemblyEntryCount=assembled.entryCount;
        c->assemblyLeafCount=assembled.leafCount;c->assemblyLowrankCount=assembled.lowrankCount;
        if(generation==0x0000ffffffffffffULL)throw std::overflow_error("Handle generation exhausted");
        const auto t=0x484d000000000000ULL|++generation;
        auto* output=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);*static_cast<std::uint64_t*>(mxGetData(output))=t;
        contexts.emplace(t,std::move(c));mexLock();mexAtExit(cleanup);out[0]=output;return;
    }
#endif
    if(command=="solve_linear_material"){
        arity(nr,7,nl,1);auto& op=*contexts.at(token(in[1]));std::size_t bytes=0;for(int i=2;i<7;++i)checkPayload(in[i],bytes);
        auto mass=readChargeMap(in[2],op.dimension());auto rhs=numbers(in[3]);auto invChi=numbers(in[4]);auto tolerance=numbers(in[5]);
        if(invChi.size()!=1||tolerance.size()!=1)throw std::invalid_argument("Expected scalar material/tolerance");
        auto result=solveLinearMaterial(op,mass,rhs,invChi[0],tolerance[0],int(integer(in[6],true)));
        const char* names[]={"x","iterations","true_relative_residual","solver_owned_device_bytes","preconditioner"};out[0]=mxCreateStructMatrix(1,1,5,names);
        auto* x=mxCreateDoubleMatrix(result.x.size(),1,mxREAL);std::copy(result.x.begin(),result.x.end(),mxGetDoubles(x));mxSetField(out[0],0,"x",x);
        mxSetField(out[0],0,"iterations",mxCreateDoubleScalar(result.iterations));mxSetField(out[0],0,"true_relative_residual",mxCreateDoubleScalar(result.relativeResidual));
        mxSetField(out[0],0,"solver_owned_device_bytes",mxCreateDoubleScalar(result.ownedBytes));mxSetField(out[0],0,"preconditioner",mxCreateString("inv_chi mass diagonal"));return;
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
        std::vector<double> a2(checkedProduct(rows,rankLimit)),a1(checkedProduct(columns,rankLimit));
        std::vector<int> lod(n+1);for(std::size_t i=0;i<n;++i)lod[i+1]=int(i+1);
        double parameters[100]={};parameters[61]=1;parameters[62]=double(rankLimit);parameters[64]=1;
        AnalyticEntryOracle oracle(block,256);std::vector<double> scale(n);
        for(std::size_t offset=0;offset<n;offset+=256){
            std::vector<AnalyticEntryIndex> diagonal;
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
            std::vector<double>().swap(a1);std::vector<double>().swap(a2);
            std::unique_ptr<mxArray,decltype(&mxDestroyArray)> dense(mxCreateDoubleMatrix(rows,columns,mxREAL),&mxDestroyArray);
            auto* values=mxGetDoubles(dense.get());
            for(std::size_t offset=0;offset<count;offset+=256){
                std::vector<AnalyticEntryIndex> indices;
                for(auto k=offset;k<std::min(count,offset+256);++k)
                    indices.push_back({k%rows,rows+k/rows});
                const auto batch=oracle.evaluate(indices,nullptr);
                std::copy(batch.begin(),batch.end(),values+offset);
            }
            const char* names[]={"rank","dense","promoted","aca_rank","entry_batches","entry_count","owned_device_bytes"};
            out[0]=mxCreateStructMatrix(1,1,7,names);
            mxSetField(out[0],0,"rank",mxCreateDoubleScalar(0));
            mxSetField(out[0],0,"dense",dense.release());
            mxSetField(out[0],0,"promoted",mxCreateLogicalScalar(true));
            mxSetField(out[0],0,"aca_rank",mxCreateDoubleScalar(rank));
            mxSetField(out[0],0,"entry_batches",mxCreateDoubleScalar(oracle.batches));
            mxSetField(out[0],0,"entry_count",mxCreateDoubleScalar(oracle.entries));
            mxSetField(out[0],0,"owned_device_bytes",mxCreateDoubleScalar(oracle.ownedBytes));return;
        }
        // Return physical factors while keeping ACA's existing normalized entry
        // semantics and standard threshold multiplier during compression.
        for(int k=0;k<rank;++k){
            for(std::size_t i=0;i<rows;++i)a2[i+rows*k]*=scale[i];
            for(std::size_t i=0;i<columns;++i)a1[i+columns*k]*=scale[rows+i];
        }
        const char* names[]={"rank","a1","a2","entry_batches","entry_count","owned_device_bytes"};
        out[0]=mxCreateStructMatrix(1,1,6,names);mxSetField(out[0],0,"rank",mxCreateDoubleScalar(rank));
        auto* first=mxCreateDoubleMatrix(columns,rank,mxREAL);auto* second=mxCreateDoubleMatrix(rows,rank,mxREAL);
        std::copy(a1.begin(),a1.begin()+columns*rank,mxGetDoubles(first));std::copy(a2.begin(),a2.begin()+rows*rank,mxGetDoubles(second));
        mxSetField(out[0],0,"a1",first);mxSetField(out[0],0,"a2",second);
        mxSetField(out[0],0,"entry_batches",mxCreateDoubleScalar(oracle.batches));mxSetField(out[0],0,"entry_count",mxCreateDoubleScalar(oracle.entries));
        mxSetField(out[0],0,"owned_device_bytes",mxCreateDoubleScalar(oracle.ownedBytes));return;
    }
    if(command=="analytic_aca_vector"){
        arity(nr,7,nl,1);std::size_t payloadBytes=0;
        for(int i=1;i<7;++i)checkPayload(in[i],payloadBytes);
        auto block=analyticSources(numbers(in[1]),numbers(in[2]),2147483647u);
        const auto n=block.sources.size(),pivot=integer(in[4],true),direction=integer(in[6]);
        const auto permutation=numbers(in[3]),maskValues=numbers(in[5]);
        if(permutation.size()!=n||maskValues.size()!=n||pivot>n||direction>1)
            throw std::invalid_argument("Invalid ACA vector shape/pivot/direction");
        std::vector<int> lod(n+1),mask(n);std::vector<bool> seen(n,false);
        for(std::size_t i=0;i<n;++i){
            const auto value=permutation[i];
            if(value<1||value>n||std::floor(value)!=value||seen[std::size_t(value)-1])
                throw std::invalid_argument("Invalid ACA permutation");
            seen[std::size_t(value)-1]=true;lod[i+1]=int(value);
            if(maskValues[i]!=0&&maskValues[i]!=1)throw std::invalid_argument("Invalid ACA mask");
            mask[i]=int(maskValues[i]);
        }
        AnalyticEntryOracle oracle(block,256);EntryCallbackScope binding(oracle,nullptr);
        std::vector<double> values(n,0.);
        // HACApK ip is zero-based relative to its one-based block start.
        cHACApK_active_entry_vector(int(n),int(pivot-1),1,1,lod.data(),0,mask.data(),int(direction),values.data());
        binding.rethrow();out[0]=mxCreateDoubleMatrix(n,1,mxREAL);
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
        std::vector<AnalyticEntryIndex> request;
        for(std::size_t i=0;i<rows.size();++i){
            if(rows[i]<1||columns[i]<1||rows[i]>block.sources.size()||columns[i]>block.sources.size()||
               std::floor(rows[i])!=rows[i]||std::floor(columns[i])!=columns[i])
                throw std::invalid_argument("Entry indices must be one-based charge indices");
            request.push_back({std::size_t(rows[i])-1,std::size_t(columns[i])-1});
        }
        AnalyticEntryOracle oracle(block,capacity);std::vector<double> result;result.reserve(request.size());
        for(std::size_t offset=0;offset<request.size();offset+=capacity){
            const auto end=std::min(offset+capacity,request.size());
            const std::vector<AnalyticEntryIndex> batch(request.begin()+offset,request.begin()+end);
            const auto values=oracle.evaluate(batch,nullptr);result.insert(result.end(),values.begin(),values.end());
        }
        out[0]=mxCreateDoubleMatrix(result.size(),1,mxREAL);std::copy(result.begin(),result.end(),mxGetDoubles(out[0]));return;
    }
    if(command=="create_curved_hdiv"||command=="create_hdiv_elements"){
        arity(nr,2,nl,1);std::size_t payloadBytes=0;checkPayload(in[1],payloadBytes);auto block=curvedBlock(in[1]);const auto n=block.charges.size();
        Leaf leaf{0,0,n,n,0,std::vector<double>(checkedProduct(n,n),0.),{}};std::vector<std::size_t> permutation(n);for(std::size_t i=0;i<n;++i)permutation[i]=i;
        auto plan=pack(n,{leaf},permutation,std::vector<double>(n,1.));auto map=readChargeMap(field(in[1],0,"charge_map"),n);
        if(generation==0x0000ffffffffffffULL)throw std::overflow_error("Handle generation exhausted");
        auto c=std::make_unique<Resident>(plan,&map);c->assembleCurved(block);const auto t=0x484d000000000000ULL|++generation;
        auto* output=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);*static_cast<std::uint64_t*>(mxGetData(output))=t;
        contexts.emplace(t,std::move(c));mexLock();mexAtExit(cleanup);out[0]=output;return;
    }
    if(command=="create"||command=="create_analytic_constant_block"){
        arity(nr,2,nl,1);if(!mxIsStruct(in[1])||mxGetNumberOfElements(in[1])!=1)throw std::invalid_argument("Expected scalar H-matrix configuration");std::size_t payloadBytes=0;checkPayload(in[1],payloadBytes);Plan plan;AnalyticBlock block;
        if(command=="create")plan=read(in[1]);
        else{
            if(const auto* assembly=mxGetField(in[1],0,"assembly")){if(text(assembly)!="analytic_constant_block")throw std::invalid_argument("Unsupported analytic assembly mode");}
            if(integer(field(in[1],0,"schema_version"),true)!=1)throw std::invalid_argument("Unsupported analytic block schema");
            for(int f=0;f<mxGetNumberOfFields(in[1]);++f){std::string key=mxGetFieldNameByNumber(in[1],f);
                if(key!="schema_version"&&key!="cells"&&key!="faces"&&key!="operator_kind"&&key!="charge_map"&&key!="assembly")
                    throw std::invalid_argument("Unsupported analytic block option: "+key);
            }
            if(mxGetN(field(in[1],0,"cells"))>1||mxGetN(field(in[1],0,"faces"))>1)
                throw std::invalid_argument("Geometry uses flat row-major coordinate column vectors");
            // No images, curved/high-order interpretation or alternate near/far strategy.
            block=analyticSources(numbers(field(in[1],0,"cells")),numbers(field(in[1],0,"faces")));
            const auto n=block.sources.size();Leaf leaf{0,0,n,n,0,std::vector<double>(checkedProduct(n,n),0.),{}};
            std::vector<std::size_t> permutation(n);for(std::size_t i=0;i<n;++i)permutation[i]=i;
            plan=pack(n,{leaf},permutation,std::vector<double>(n,1.));
        }
        if(generation==0x0000ffffffffffffULL)throw std::overflow_error("Handle generation exhausted");
        const auto mode=text(field(in[1],0,"operator_kind"));SparseMap map;
        if(mode=="demag")map=readChargeMap(field(in[1],0,"charge_map"),plan.n);
        else if(mode!="gram")throw std::invalid_argument("Operator kind must be gram or demag");
        auto c=std::make_unique<Resident>(plan,mode=="demag"?&map:nullptr);if(command=="create_analytic_constant_block")c->assemble(block);const auto t=0x484d000000000000ULL|++generation;
        auto* output=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);*static_cast<std::uint64_t*>(mxGetData(output))=t;
        contexts.emplace(t,std::move(c));mexLock();mexAtExit(cleanup);out[0]=output;return;
    }
    if(command=="live_count"){arity(nr,1,nl,1);out[0]=mxCreateDoubleScalar(contexts.size());return;}
    if(command=="destroy"){arity(nr,2,nl,0);contexts.erase(token(in[1]));mexUnlock();return;}
    if(command=="upload"){arity(nr,3,nl,0);contexts.at(token(in[1]))->upload(numbers(in[2]));return;}
    if(command=="apply"){arity(nr,3,nl,0);contexts.at(token(in[1]))->apply(integer(in[2],true));return;}
    if(command=="download"){
        arity(nr,2,nl,1);auto v=contexts.at(token(in[1]))->download();out[0]=mxCreateDoubleMatrix(v.size(),1,mxREAL);std::copy(v.begin(),v.end(),mxGetDoubles(out[0]));return;
    }
    if(command=="stats"){
        arity(nr,2,nl,1);const auto& c=*contexts.at(token(in[1]));const char* names[]={"owned_device_bytes","applications","kernel_launches","host_synchronizations_per_apply","setup_peak_owned_device_bytes"};
        out[0]=mxCreateStructMatrix(1,1,5,names);mxSetField(out[0],0,names[0],mxCreateDoubleScalar(c.ownedBytes));mxSetField(out[0],0,names[1],mxCreateDoubleScalar(c.applications));mxSetField(out[0],0,names[2],mxCreateDoubleScalar(c.kernelsPerApplication*c.applications));mxSetField(out[0],0,names[3],mxCreateDoubleScalar(0));mxSetField(out[0],0,names[4],mxCreateDoubleScalar(std::max(c.ownedBytes,c.setupPeakOwnedBytes)));
        const char* extra[]={"assembly_entry_batches","assembly_entry_count","assembly_leaf_count","assembly_lowrank_count"};
        const std::size_t counts[]={c.assemblyEntryBatches,c.assemblyEntryCount,c.assemblyLeafCount,c.assemblyLowrankCount};
        for(int i=0;i<4;++i){mxAddField(out[0],extra[i]);mxSetField(out[0],0,extra[i],mxCreateDoubleScalar(counts[i]));}return;
    }
    throw std::invalid_argument("Unknown H-matrix command");
}catch(const std::exception& e){mexErrMsgIdAndTxt("radia:cuda:hmatrix","%s",e.what());}
}
