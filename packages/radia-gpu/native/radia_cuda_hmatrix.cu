// Experimental numeric leaf boundary, not full HDiv assembly or solver acceptance.
#include "mex.h"
#include "hmatrix_hacapk_pack.cuh"
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
    return pack(n,leaves,permutation,numbers(field(cfg,0,"scale")));
}
}
void mexFunction(int nl,mxArray* out[],int nr,const mxArray* in[]){
try{
    if(!nr)throw std::invalid_argument("Command required");auto command=text(in[0]);
    if(command=="create"){
        arity(nr,2,nl,1);auto plan=read(in[1]);
        if(generation==0x0000ffffffffffffULL)throw std::overflow_error("Handle generation exhausted");
        const auto mode=text(field(in[1],0,"operator_kind"));SparseMap map;
        if(mode=="demag")map=readChargeMap(field(in[1],0,"charge_map"),plan.n);
        else if(mode!="gram")throw std::invalid_argument("Operator kind must be gram or demag");
        auto c=std::make_unique<Resident>(plan,mode=="demag"?&map:nullptr);const auto t=0x484d000000000000ULL|++generation;
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
        arity(nr,2,nl,1);const auto& c=*contexts.at(token(in[1]));const char* names[]={"owned_device_bytes","applications","kernel_launches","host_synchronizations_per_apply"};
        out[0]=mxCreateStructMatrix(1,1,4,names);mxSetField(out[0],0,names[0],mxCreateDoubleScalar(c.ownedBytes));mxSetField(out[0],0,names[1],mxCreateDoubleScalar(c.applications));mxSetField(out[0],0,names[2],mxCreateDoubleScalar(c.kernelsPerApplication*c.applications));mxSetField(out[0],0,names[3],mxCreateDoubleScalar(0));return;
    }
    throw std::invalid_argument("Unknown H-matrix command");
}catch(const std::exception& e){mexErrMsgIdAndTxt("radia:cuda:hmatrix","%s",e.what());}
}
