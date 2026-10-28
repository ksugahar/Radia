/* SPDX-License-Identifier: MPL-2.0 */
#include "snapshot.hpp"
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>
namespace {
#include "ams_core.cuh"

struct Allocation {
    void* p=nullptr;
    explicit Allocation(size_t bytes){ck(cudaMalloc(&p,std::max<size_t>(bytes,1)));}
    ~Allocation(){if(p)cudaFree(p);}
};
Node upload(const Node& input,std::vector<std::unique_ptr<Allocation>>& arrays){
    Node out;out.kind=input.kind;out.scalar=input.scalar;out.boolean=input.boolean;
    out.word=input.word;out.count=input.count;
    for(auto& item:input.fields)out.fields[item.first]=upload(item.second,arrays);
    for(auto& item:input.children)out.children.push_back(upload(item,arrays));
    if(input.kind==Node::integers||input.kind==Node::reals){
        if(input.count>size_t(INT_MAX-256))throw std::runtime_error("snapshot array capacity exceeded");
        const bool integers=input.kind==Node::integers;
        if(input.count!=(integers?input.int_values.size():input.real_values.size()))
            throw std::runtime_error("snapshot host array length mismatch");
        size_t bytes=input.count*(integers?sizeof(int):sizeof(double));
        auto storage=std::make_unique<Allocation>(bytes);
        if(bytes)ck(cudaMemcpy(storage->p,integers?static_cast<const void*>(input.int_values.data()):
                      static_cast<const void*>(input.real_values.data()),bytes,cudaMemcpyHostToDevice));
        out.device=storage->p;arrays.push_back(std::move(storage));
    }
    return out;
}
Node number(int value){Node n;n.kind=Node::number;n.scalar=value;return n;}
Node identity(int n){
    Node result;result.kind=Node::record;result.fields["rows"]=number(n);result.fields["cols"]=number(n);
    Node row,col,val;row.kind=col.kind=Node::integers;val.kind=Node::reals;
    row.count=size_t(n)+1;col.count=val.count=n;
    for(int i=0;i<n;++i){row.int_values.push_back(i);col.int_values.push_back(i);val.real_values.push_back(1.);}
    row.int_values.push_back(n);result.fields["rowptr"]=std::move(row);
    result.fields["columns"]=std::move(col);result.fields["values"]=std::move(val);return result;
}
int rows(const Node& node){
    int n=node.at("A").at("rows").integer();
    if(n<1||n>INT_MAX-256)throw std::runtime_error("invalid snapshot system dimension");
    return n;
}
struct Context {
    int device,n;
    std::vector<std::unique_ptr<Allocation>> arrays;
    Node data;
    Flag validation;
    Buffer rhs,output;
    std::unique_ptr<AMS> ams;
    std::unique_ptr<Space> amg;
    explicit Context(const Node& input):n(rows(input)),rhs(n),output(n){
        ck(cudaGetDevice(&device));
        const auto& schema=input.at("schema").string();
        if(schema=="radia.ams-accelerator.v1"){
            data=upload(input,arrays);ams=std::make_unique<AMS>(data,validation);
        }else if(schema=="radia.amg-accelerator.v1"){
            Node space;space.kind=Node::record;space.fields["P"]=identity(n);space.fields["Pt"]=identity(n);
            space.fields["levels"]=input.at("levels");
            if(space.at("levels").items().empty())throw std::runtime_error("empty AMG hierarchy");
            const auto& system=input.at("A");
            const auto& first=space.at("levels").items().front().at("A");
            if(system.at("rows").integer()!=first.at("rows").integer()
                ||system.at("cols").integer()!=first.at("cols").integer()
                ||system.at("rowptr").int_values!=first.at("rowptr").int_values
                ||system.at("columns").int_values!=first.at("columns").int_values
                ||system.at("values").real_values!=first.at("values").real_values)
                throw std::runtime_error("AMG operator and finest-level snapshot differ");
            data=upload(space,arrays);amg=std::make_unique<Space>(data,validation);
            if(amg->levels[0]->a.n!=n)throw std::runtime_error("AMG snapshot operator dimension mismatch");
        }else throw std::runtime_error("unsupported sparse accelerator schema");
        validation.check();ck(cudaGetLastError());
    }
    void selected(){
        int current;ck(cudaGetDevice(&current));
        if(current!=device)throw std::runtime_error("GPU hierarchy belongs to a different CUDA device");
        launch_stream=nullptr;
    }
    void input(const double* b,size_t count,bool solve){
        selected();if(count!=size_t(n))throw std::runtime_error("RHS dimension mismatch");
        ck(cudaMemset(validation.p,0,sizeof(int)));
        ck(cudaMemcpy(rhs.p,b,n*sizeof(double),cudaMemcpyHostToDevice));
        verify_vec<<<(n+T-1)/T,T>>>(n,rhs.p,false,validation.p);
        if(ams){
            output.zero();verify_mask<<<(n+T-1)/T,T>>>(n,ams->free,solve?rhs.p:output.p,validation.p);
        }
        validation.check();
    }
    const CSR& matrix()const{return ams?ams->a:amg->levels[0]->a;}
};
struct AMGAction {
    const CSR& a;Space& space;
    void apply(const double* b,double* x){
        space.apply(b);ck(cudaMemcpy(x,space.out.p,a.n*sizeof(double),cudaMemcpyDeviceToDevice));
    }
};
}
namespace radia_sparse {
void* create(const Node& node){return new Context(node);}
Stats solve(void* owner,const double* b,size_t n,double* out,double tol,int limit){
    if(!owner||!std::isfinite(tol)||tol<=0||tol>1e-6||limit<1||limit>100000)
        throw std::runtime_error("invalid sparse GPU solve controls");
    auto& c=*static_cast<Context*>(owner);c.input(b,n,true);
    ResidentWork work(c.n);SolveStats result;
    if(c.ams)result=residentPCG(*c.ams,c.matrix(),c.rhs.p,work,tol,limit,nullptr);
    else {AMGAction action{c.matrix(),*c.amg};result=residentPCG(action,c.matrix(),c.rhs.p,work,tol,limit,nullptr);}
    ck(cudaGetLastError());ck(cudaMemcpy(out,work.x.p,n*sizeof(double),cudaMemcpyDeviceToHost));
    return {result.iterations,result.residual};
}
void apply(void* owner,const double* b,size_t n,double* out){
    if(!owner)throw std::runtime_error("closed sparse GPU hierarchy");
    auto& c=*static_cast<Context*>(owner);c.input(b,n,false);
    if(c.ams)c.ams->apply(c.rhs.p,c.output.p);
    else {AMGAction action{c.matrix(),*c.amg};action.apply(c.rhs.p,c.output.p);}
    verify_vec<<<(c.n+T-1)/T,T>>>(c.n,c.output.p,false,c.validation.p);
    c.validation.check();ck(cudaGetLastError());
    ck(cudaMemcpy(out,c.output.p,n*sizeof(double),cudaMemcpyDeviceToHost));
}
void release(void* owner) noexcept{
    if(!owner)return;
    auto p=static_cast<Context*>(owner);int previous=-1;cudaGetDevice(&previous);
    cudaSetDevice(p->device);delete p;
    if(previous>=0)cudaSetDevice(previous);
}
}
