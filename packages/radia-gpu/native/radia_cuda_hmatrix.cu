// Experimental numeric leaf boundary, not full HDiv assembly or solver acceptance.
#include "mex.h"
#include "hmatrix_hacapk_pack.cuh"
#include "hmatrix_entry_oracle.cuh"
#include "hmatrix_aca_callback.cuh"
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
    if(command=="analytic_aca_compress"){
        arity(nr,6,nl,1);std::size_t payloadBytes=0;
        for(int i=1;i<6;++i)checkPayload(in[i],payloadBytes);
        auto block=analyticSources(numbers(in[1]),numbers(in[2]),2147483647u);
        const auto rows=integer(in[3],true),rankLimit=integer(in[5],true),n=block.sources.size();
        const auto epsilon=numbers(in[4]);
        if(rows>=n||rankLimit>std::min(rows,n-rows)||epsilon.size()!=1||epsilon[0]<=0||epsilon[0]>=1)
            throw std::invalid_argument("Invalid ACA block split/rank/tolerance");
        const auto columns=n-rows;
        std::vector<double> a2(checkedProduct(rows,rankLimit)),a1(checkedProduct(columns,rankLimit));
        std::vector<int> lod(n+1);for(std::size_t i=0;i<n;++i)lod[i+1]=int(i+1);
        double parameters[100]={};parameters[61]=1;parameters[62]=double(rankLimit);parameters[64]=1;
        AnalyticEntryOracle oracle(block,256);EntryCallbackScope binding(oracle,nullptr);
        const int rank=cHACApK_acaplus(a2.data(),a1.data(),parameters,int(rows),int(columns),1,int(rows+1),lod.data(),0,int(rankLimit),epsilon[0],1.,epsilon[0]);
        binding.rethrow();if(rank<1||rank>int(rankLimit))throw std::runtime_error("ACA compression did not return valid factors");
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
        out[0]=mxCreateStructMatrix(1,1,5,names);mxSetField(out[0],0,names[0],mxCreateDoubleScalar(c.ownedBytes));mxSetField(out[0],0,names[1],mxCreateDoubleScalar(c.applications));mxSetField(out[0],0,names[2],mxCreateDoubleScalar(c.kernelsPerApplication*c.applications));mxSetField(out[0],0,names[3],mxCreateDoubleScalar(0));mxSetField(out[0],0,names[4],mxCreateDoubleScalar(std::max(c.ownedBytes,c.setupPeakOwnedBytes)));return;
    }
    throw std::invalid_argument("Unknown H-matrix command");
}catch(const std::exception& e){mexErrMsgIdAndTxt("radia:cuda:hmatrix","%s",e.what());}
}
