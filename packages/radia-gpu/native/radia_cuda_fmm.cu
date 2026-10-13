/* SPDX-License-Identifier: MPL-2.0 */
// Laplace scalar/dipole FMM: P2M/M2M/M2L/L2L/L2P plus exact point near field.
// Core HDiv compact interactions remain HACApK. No FE or BEM quadrature is
// inferred here: callers own density maps and analytic/singular corrections.
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cuda_runtime.h>
#include "fmm_interpolation.cuh"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cstdint>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using namespace radia_fmm;
mxArray* RadiaPrepareBoundary(const mxArray*,const mxArray*);
namespace {
void check(cudaError_t s){if(s!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(s));}
size_t product(size_t a,size_t b){if(b&&a>SIZE_MAX/b)throw std::overflow_error("FMM allocation overflow");return a*b;}
void room(size_t bytes){size_t free,total;check(cudaMemGetInfo(&free,&total));size_t reserve=std::max<size_t>(2ULL<<30,total/5);
    if(free<=reserve||bytes>free-reserve)throw std::runtime_error("FMM exceeds CUDA reserve max(2 GiB,20%)");}
struct Buffer {
    void* p=nullptr;size_t bytes=0;
    Buffer()=default;~Buffer(){if(p)cudaFree(p);}Buffer(const Buffer&)=delete;
    void allocate(size_t n){check(cudaMalloc(&p,std::max<size_t>(1,n)));bytes=n;}
    template<class V> V* as(){return static_cast<V*>(p);}
    template<class V> void upload(const std::vector<V>& v){allocate(product(v.size(),sizeof(V)));if(bytes)check(cudaMemcpy(p,v.data(),bytes,cudaMemcpyHostToDevice));}
};
std::string text(const mxArray* a){if(!a||!mxIsChar(a)||mxGetM(a)!=1)throw std::invalid_argument("Expected char row");char* p=mxArrayToString(a);if(!p)throw std::bad_alloc();std::string s(p);mxFree(p);return s;}
int number(const mxArray* a,int lo,int hi){if(!mxIsDouble(a)||mxIsSparse(a)||mxIsComplex(a)||mxGetNumberOfElements(a)!=1)throw std::invalid_argument("Expected integer double scalar");double v=mxGetScalar(a);
    if(!std::isfinite(v)||v!=floor(v)||v<lo||v>hi)throw std::invalid_argument("FMM scalar outside supported range");return int(v);}
double opening(const mxArray* a){if(!mxIsDouble(a)||mxIsSparse(a)||mxIsComplex(a)||mxGetNumberOfElements(a)!=1)throw std::invalid_argument("Expected separation scalar");double v=mxGetScalar(a);if(!std::isfinite(v)||v<1.1||v>8)throw std::invalid_argument("Separation must be in [1.1,8]");return v;}
std::vector<double> positions(const mxArray* a){if(!mxIsDouble(a)||mxIsComplex(a)||mxIsSparse(a)||mxGetNumberOfDimensions(a)!=2||mxGetN(a)!=3||!mxGetM(a)||mxGetM(a)>2000000)throw std::invalid_argument("Positions must be 1..2000000 by 3 full finite doubles");
    size_t n=mxGetM(a);auto* v=mxGetPr(a);std::vector<double> r(n*3);for(size_t i=0;i<n;++i)for(int c=0;c<3;++c){double x=v[i+n*c];if(!std::isfinite(x))throw std::invalid_argument("Nonfinite position");r[3*i+c]=x;}return r;}
struct Tree {
    std::vector<Node> nodes;std::vector<int> sources,targets,targetLeaf,farRp,far,nearRp,near;
    std::vector<std::vector<int>> levels,farRows,nearRows;int leafCapacity;double separation;
    const std::vector<double>& source;const std::vector<double>& target;
    Tree(const std::vector<double>& s,const std::vector<double>& t,int leaf,double sep):leafCapacity(leaf),separation(sep),source(s),target(t){
        std::vector<int> si(s.size()/3),ti(t.size()/3);std::iota(si.begin(),si.end(),0);std::iota(ti.begin(),ti.end(),0);targetLeaf.resize(ti.size());
        double lo[3]={s[0],s[1],s[2]},hi[3]={s[0],s[1],s[2]},center[3];
        for(auto* v:{&s,&t})for(size_t i=0;i<v->size()/3;++i)for(int c=0;c<3;++c){lo[c]=std::min(lo[c],(*v)[3*i+c]);hi[c]=std::max(hi[c],(*v)[3*i+c]);}
        double width=0;for(int c=0;c<3;++c){center[c]=lo[c]/2+hi[c]/2;width=std::max(width,hi[c]-lo[c]);}
        if(!std::isfinite(width)||width<=0)throw std::invalid_argument("Degenerate or overflowing FMM bounding box");
        build(si,ti,center,width*.500000000001,-1,0,0);farRows.resize(nodes.size());nearRows.resize(nodes.size());walk(0,0);
        flatten(farRows,farRp,far);flatten(nearRows,nearRp,near);
    }
    int build(const std::vector<int>& si,const std::vector<int>& ti,const double* center,double half,int parent,int octant,int depth){
        if(nodes.size()>=200000)throw std::overflow_error("FMM node capacity exceeded");int id=int(nodes.size());Node n{};
        std::copy(center,center+3,n.center);n.half=half;n.parent=parent;n.octant=octant;n.depth=depth;n.sourceBegin=int(sources.size());n.targetBegin=int(targets.size());n.sourceCount=int(si.size());n.targetCount=int(ti.size());std::fill(n.children,n.children+8,-1);nodes.push_back(n);
        if(levels.size()<=size_t(depth))levels.resize(depth+1);levels[depth].push_back(id);
        if(std::max(si.size(),ti.size())<=size_t(leafCapacity)||depth==12){sources.insert(sources.end(),si.begin(),si.end());targets.insert(targets.end(),ti.begin(),ti.end());for(int i:ti)targetLeaf[i]=id;return id;}
        std::vector<int> ss[8],tt[8];for(int i:si){int c=0;for(int k=0;k<3;++k)if(source[3*i+k]>=center[k])c|=1<<k;ss[c].push_back(i);}for(int i:ti){int c=0;for(int k=0;k<3;++k)if(target[3*i+k]>=center[k])c|=1<<k;tt[c].push_back(i);}
        for(int c=0;c<8;++c)if(!ss[c].empty()||!tt[c].empty()){double x[3];for(int k=0;k<3;++k)x[k]=center[k]+(c&(1<<k)?half/2:-half/2);int child=build(ss[c],tt[c],x,half/2,id,c,depth+1);nodes[id].children[c]=child;}return id;
    }
    void walk(int t,int s){const auto& a=nodes[t];const auto& b=nodes[s];if(!a.targetCount||!b.sourceCount)return;
        double d2=0;for(int k=0;k<3;++k){double d=a.center[k]-b.center[k];d2+=d*d;}double bound=separation*sqrt(3.)*(a.half+b.half);
        if(d2>bound*bound){farRows[t].push_back(s);return;}
        if(a.leaf()&&b.leaf()){nearRows[t].push_back(s);return;}
        if(!a.leaf()&&(b.leaf()||a.half>=b.half)){for(int child:a.children)if(child>=0)walk(child,s);}
        else for(int child:b.children)if(child>=0)walk(t,child);
    }
    static void flatten(const std::vector<std::vector<int>>& rows,std::vector<int>& rp,std::vector<int>& values){rp.push_back(0);for(const auto& row:rows){if(values.size()+row.size()>20000000)throw std::overflow_error("FMM interaction capacity exceeded");values.insert(values.end(),row.begin(),row.end());rp.push_back(int(values.size()));}}
};
__global__ void validate(int n,const double* v,int* error){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n&&!isfinite(v[i]))atomicExch(error,1);}
__global__ void p2m(int count,const Node* nodes,int p,int ns,int nc,const double* xyz,const int* order,
    const double* charge,const double* dip,double* moments){int q=p*p*p,index=blockIdx.x*blockDim.x+threadIdx.x;if(index>=count*q*nc)return;int k=index%q,c=(index/q)%nc,id=index/(q*nc);auto node=nodes[id];if(!node.leaf()||!node.sourceCount)return;
    double value=0;for(int i=node.sourceBegin;i<node.sourceBegin+node.sourceCount;++i){int s=order[i];double v[3][MAX_ORDER],d[3][MAX_ORDER];for(int a=0;a<3;++a)basis((xyz[s*3+a]-node.center[a])/node.half,p,v[a],d[a]);int a=k%p,b=(k/p)%p,e=k/(p*p);
        value+=charge[s+size_t(ns)*c]*v[0][a]*v[1][b]*v[2][e];
        value+=(dip[s+size_t(ns)*(3*c)]*d[0][a]*v[1][b]*v[2][e]+dip[s+size_t(ns)*(3*c+1)]*v[0][a]*d[1][b]*v[2][e]+dip[s+size_t(ns)*(3*c+2)]*v[0][a]*v[1][b]*d[2][e])/node.half;
    }moments[(size_t(id)*nc+c)*q+k]=value;
}
__global__ void m2m(int count,const int* level,const Node* nodes,int q,int nc,const double* transfer,double* m){int index=blockIdx.x*blockDim.x+threadIdx.x;if(index>=count*q*nc)return;int k=index%q,c=(index/q)%nc,id=level[index/(q*nc)];auto node=nodes[id];if(node.leaf())return;
    double value=0;for(int child:node.children)if(child>=0){auto ch=nodes[child];for(int j=0;j<q;++j)value+=transfer[(size_t(ch.octant)*q+k)*q+j]*m[(size_t(child)*nc+c)*q+j];}m[(size_t(id)*nc+c)*q+k]=value;
}
__global__ void m2l(int count,const Node* nodes,int p,int nc,const int* rp,const int* sources,const double* m,double* local){int q=p*p*p,index=blockIdx.x*blockDim.x+threadIdx.x;if(index>=count*q*nc)return;int k=index%q,c=(index/q)%nc,id=index/(q*nc);double x[3];point(nodes[id],k,p,x);double value=0;
    for(int i=rp[id];i<rp[id+1];++i){int s=sources[i];for(int j=0;j<q;++j){double y[3];point(nodes[s],j,p,y);double d2=0;for(int a=0;a<3;++a){double d=x[a]-y[a];d2+=d*d;}value+=m[(size_t(s)*nc+c)*q+j]/sqrt(d2);}}
    local[(size_t(id)*nc+c)*q+k]=value;
}
__global__ void l2l(int count,const int* level,const Node* nodes,int q,int nc,const double* transfer,double* local){int index=blockIdx.x*blockDim.x+threadIdx.x;if(index>=count*q*nc)return;int k=index%q,c=(index/q)%nc,id=level[index/(q*nc)];auto n=nodes[id];if(n.parent<0||!n.targetCount)return;double value=0;
    for(int j=0;j<q;++j)value+=transfer[(size_t(n.octant)*q+j)*q+k]*local[(size_t(n.parent)*nc+c)*q+j];local[(size_t(id)*nc+c)*q+k]+=value;
}
__global__ void evaluate(int nt,int ns,const double* target,const double* source,const int* targetLeaf,const Node* nodes,
    int p,int nc,const double* local,const int* nearRp,const int* near,const int* sourceOrder,const double* charge,const double* dip,bool omitSelf,
    double* potential,double* field,int* error){int t=blockIdx.x*blockDim.x+threadIdx.x;if(t>=nt)return;int id=targetLeaf[t],q=p*p*p;auto node=nodes[id];const double* x=target+3*t;
    double v[3][MAX_ORDER],d[3][MAX_ORDER];for(int a=0;a<3;++a)basis((x[a]-node.center[a])/node.half,p,v[a],d[a]);
    for(int c=0;c<nc;++c){double phi=0,E[3]={0,0,0};for(int k=0;k<q;++k){int a=k%p,b=(k/p)%p,e=k/(p*p);double l=local[(size_t(id)*nc+c)*q+k];phi+=l*v[0][a]*v[1][b]*v[2][e];E[0]-=l*d[0][a]*v[1][b]*v[2][e]/node.half;E[1]-=l*v[0][a]*d[1][b]*v[2][e]/node.half;E[2]-=l*v[0][a]*v[1][b]*d[2][e]/node.half;}
        for(int j=nearRp[id];j<nearRp[id+1];++j){auto n=nodes[near[j]];for(int k=n.sourceBegin;k<n.sourceBegin+n.sourceCount;++k){int s=sourceOrder[k];const double* y=source+3*s;double r2=0;for(int a=0;a<3;++a){double z=x[a]-y[a];r2+=z*z;}if(r2<=1e-300){if(!omitSelf)atomicExch(error,2);continue;}double dipole[3];for(int a=0;a<3;++a)dipole[a]=dip[s+size_t(ns)*(3*c+a)];direct(x,y,charge[s+size_t(ns)*c],dipole,phi,E);}}
        if(!isfinite(phi)||!isfinite(E[0])||!isfinite(E[1])||!isfinite(E[2]))atomicExch(error,3);
        potential[t+size_t(nt)*c]=phi;for(int a=0;a<3;++a)field[t+size_t(nt)*(3*c+a)]=E[a];}
}
struct Context {
    Buffer nodes,source,target,sourceOrder,targetLeaf,farRp,far,nearRp,near,transfer,moments,local,error;
    std::vector<std::unique_ptr<Buffer>> levels;std::vector<int> levelSizes;int ns,nt,nnode,p,q,device;bool omitSelf;size_t persistent=0,farPairs=0,nearPairs=0,nonleafFarPairs=0;cudaStream_t sentinel=nullptr;
    Context(const std::vector<double>& s,const std::vector<double>& t,int order,int leaf,double separation,bool omit):ns(int(s.size()/3)),nt(int(t.size()/3)),p(order),q(p*p*p),omitSelf(omit){
        Tree tree(s,t,leaf,separation);nnode=int(tree.nodes.size());check(cudaGetDevice(&device));
        farPairs=tree.far.size();nearPairs=tree.near.size();
        for(size_t i=0;i<tree.nodes.size();++i)if(!tree.nodes[i].leaf())nonleafFarPairs+=tree.farRows[i].size();
        std::vector<double> translation(size_t(8)*q*q);
        for(int oct=0;oct<8;++oct)for(int j=0;j<q;++j){double v[3][MAX_ORDER],d[MAX_ORDER];for(int a=0;a<3;++a){double x=.5*((oct&(1<<a)?1.:-1.)+abscissa(a==0?j%p:(a==1?(j/p)%p:j/(p*p)),p));basis(x,p,v[a],d);}for(int k=0;k<q;++k)translation[(size_t(oct)*q+k)*q+j]=v[0][k%p]*v[1][(k/p)%p]*v[2][k/(p*p)];}
        size_t bytes=s.size()*8+t.size()*8+tree.nodes.size()*sizeof(Node)+translation.size()*8+sizeof(int);
        for(auto* v:{&tree.sources,&tree.targetLeaf,&tree.farRp,&tree.far,&tree.nearRp,&tree.near})bytes+=v->size()*sizeof(int);
        for(const auto& v:tree.levels)bytes+=v.size()*sizeof(int);room(bytes);persistent=bytes;
        nodes.upload(tree.nodes);source.upload(s);target.upload(t);sourceOrder.upload(tree.sources);targetLeaf.upload(tree.targetLeaf);farRp.upload(tree.farRp);far.upload(tree.far);nearRp.upload(tree.nearRp);near.upload(tree.near);transfer.upload(translation);error.allocate(sizeof(int));
        for(auto& ids:tree.levels){auto b=std::make_unique<Buffer>();b->upload(ids);levelSizes.push_back(int(ids.size()));levels.push_back(std::move(b));}
        check(cudaStreamCreate(&sentinel));
    }
    ~Context(){int old=device;cudaGetDevice(&old);if(old!=device)cudaSetDevice(device);auto state=sentinel?cudaStreamQuery(sentinel):cudaErrorInvalidResourceHandle;
        if(state!=cudaSuccess&&state!=cudaErrorNotReady){for(auto* b:{&nodes,&source,&target,&sourceOrder,&targetLeaf,&farRp,&far,&nearRp,&near,&transfer,&moments,&local,&error})b->p=nullptr;for(auto& b:levels)b->p=nullptr;}
        else {cudaDeviceSynchronize();cudaStreamDestroy(sentinel);
            for(auto* b:{&nodes,&source,&target,&sourceOrder,&targetLeaf,&farRp,&far,&nearRp,&near,&transfer,&moments,&local,&error}){if(b->p)cudaFree(b->p);b->p=nullptr;}
            for(auto& b:levels){if(b->p)cudaFree(b->p);b->p=nullptr;}}
        if(old!=device)cudaSetDevice(old);
    }
    void valid(){int d;check(cudaGetDevice(&d));if(d!=device)throw std::runtime_error("FMM context belongs to another CUDA device");check(cudaStreamQuery(sentinel));}
};
std::map<uint64_t,std::unique_ptr<Context>> contexts;uint64_t generation=0;
void cleanup(){contexts.clear();}
uint64_t token(const mxArray* a){if(!mxIsUint64(a)||mxIsComplex(a)||mxGetNumberOfElements(a)!=1)throw std::invalid_argument("FMM handle must be scalar uint64");uint64_t h=*static_cast<const uint64_t*>(mxGetData(a));if(!contexts.count(h))throw std::invalid_argument("Stale or foreign FMM handle");return h;}
struct GPUView {const mxGPUArray* a;GPUView(const mxArray* x){if(!mxIsGPUArray(x))throw std::invalid_argument("FMM densities must be gpuArray doubles");a=mxGPUCreateFromMxArray(x);if(mxGPUGetClassID(a)!=mxDOUBLE_CLASS||mxGPUGetComplexity(a)!=mxREAL){mxGPUDestroyGPUArray(a);throw std::invalid_argument("FMM densities must be real double");}}
    ~GPUView(){mxGPUDestroyGPUArray(a);}const double* data(){return static_cast<const double*>(mxGPUGetDataReadOnly(a));}};
struct GPUOutput {mxGPUArray* a;GPUOutput(int n,int m){mwSize dims[2]={size_t(n),size_t(m)};a=mxGPUCreateGPUArray(2,dims,mxDOUBLE_CLASS,mxREAL,MX_GPU_DO_NOT_INITIALIZE);}~GPUOutput(){mxGPUDestroyGPUArray(a);}double* data(){return static_cast<double*>(mxGPUGetData(a));}};
void dispatch(int nl,mxArray** out,int nr,const mxArray** in){if(!nr)throw std::invalid_argument("FMM command required");std::string command=text(in[0]);
    if(command=="prepare_boundary"){if(nr!=3||nl!=1)throw std::invalid_argument("maps=radia_cuda_fmm('prepare_boundary',vol,boundaryLabels)");out[0]=RadiaPrepareBoundary(in[1],in[2]);return;}
    mxInitGPU();
    if(command=="live_count"){if(nr!=1||nl!=1)throw std::invalid_argument("live_count arity");out[0]=mxCreateDoubleScalar(contexts.size());return;}
    if(command=="destroy"){if(nr!=2||nl!=0)throw std::invalid_argument("destroy arity");contexts.erase(token(in[1]));mexUnlock();return;}
    if(command=="create"){if(nr!=7||nl!=1)throw std::invalid_argument("handle=radia_cuda_fmm('create',sources,targets,order,leafCapacity,separation,omitSelf)");if(!mxIsLogicalScalar(in[6]))throw std::invalid_argument("omitSelf must be explicit logical");
        auto source=positions(in[1]),target=positions(in[2]);auto context=std::make_unique<Context>(source,target,number(in[3],3,8),number(in[4],8,512),opening(in[5]),mxIsLogicalScalarTrue(in[6]));
        if(generation>=0x0000ffffffffffffULL)throw std::overflow_error("FMM handle generation exhausted");uint64_t h=0xF66A000000000000ULL|++generation;out[0]=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);*static_cast<uint64_t*>(mxGetData(out[0]))=h;contexts.emplace(h,std::move(context));mexAtExit(cleanup);mexLock();return;}
    if(command!="apply"||nr!=4||nl!=3)throw std::invalid_argument("[potential,field,info]=radia_cuda_fmm('apply',handle,charge,dipole)");
    auto& c=*contexts.at(token(in[1]));c.valid();GPUView charge(in[2]),dip(in[3]);size_t nq=mxGPUGetNumberOfElements(charge.a),nd=mxGPUGetNumberOfElements(dip.a);
    if(nq%size_t(c.ns)||nq/size_t(c.ns)<1||nq/size_t(c.ns)>3||nd!=nq*3)throw std::invalid_argument("FMM charge length ns*channels (1..3); dipole length ns*channels*3");
    int nc=int(nq/c.ns);size_t coefficients=product(product(size_t(c.nnode),c.q),nc),workBytes=product(coefficients,16),outputBytes=product(size_t(c.nt)*nc,32);
    if(coefficients>size_t(INT_MAX-256)||workBytes>size_t(2ULL<<30))throw std::overflow_error("FMM coefficient capacity exceeds bounded int32/2 GiB workflow");
    if(!c.moments.p||c.moments.bytes!=coefficients*8){room(workBytes+outputBytes);Buffer nextM,nextL;nextM.allocate(coefficients*8);nextL.allocate(coefficients*8);std::swap(c.moments.p,nextM.p);std::swap(c.moments.bytes,nextM.bytes);std::swap(c.local.p,nextL.p);std::swap(c.local.bytes,nextL.bytes);}
    else room(outputBytes);
    GPUOutput potential(c.nt,nc),field(c.nt,3*nc);check(cudaMemset(c.error.p,0,sizeof(int)));
    validate<<<(int(nq)+255)/256,256>>>(int(nq),charge.data(),c.error.as<int>());validate<<<(int(nd)+255)/256,256>>>(int(nd),dip.data(),c.error.as<int>());int error;check(cudaMemcpy(&error,c.error.p,sizeof(int),cudaMemcpyDeviceToHost));if(error)throw std::invalid_argument("Nonfinite FMM density");
    auto start=std::chrono::steady_clock::now();check(cudaMemset(c.moments.p,0,c.moments.bytes));check(cudaMemset(c.local.p,0,c.local.bytes));
    p2m<<<(int(coefficients)+127)/128,128>>>(c.nnode,c.nodes.as<Node>(),c.p,c.ns,nc,c.source.as<double>(),c.sourceOrder.as<int>(),charge.data(),dip.data(),c.moments.as<double>());
    for(int l=int(c.levels.size())-2;l>=0;--l){int n=c.levelSizes[l]*c.q*nc;m2m<<<(n+127)/128,128>>>(c.levelSizes[l],c.levels[l]->as<int>(),c.nodes.as<Node>(),c.q,nc,c.transfer.as<double>(),c.moments.as<double>());}
    m2l<<<(int(coefficients)+127)/128,128>>>(c.nnode,c.nodes.as<Node>(),c.p,nc,c.farRp.as<int>(),c.far.as<int>(),c.moments.as<double>(),c.local.as<double>());
    for(size_t l=1;l<c.levels.size();++l){int n=c.levelSizes[l]*c.q*nc;l2l<<<(n+127)/128,128>>>(c.levelSizes[l],c.levels[l]->as<int>(),c.nodes.as<Node>(),c.q,nc,c.transfer.as<double>(),c.local.as<double>());}
    evaluate<<<(c.nt+127)/128,128>>>(c.nt,c.ns,c.target.as<double>(),c.source.as<double>(),c.targetLeaf.as<int>(),c.nodes.as<Node>(),c.p,nc,c.local.as<double>(),c.nearRp.as<int>(),c.near.as<int>(),c.sourceOrder.as<int>(),charge.data(),dip.data(),c.omitSelf,potential.data(),field.data(),c.error.as<int>());
    check(cudaGetLastError());check(cudaMemcpy(&error,c.error.p,sizeof(int),cudaMemcpyDeviceToHost));if(error)throw std::runtime_error(error==2?"Coincident FMM source/target: select explicit self omission only with a declared correction":"FMM numerical overflow");check(cudaDeviceSynchronize());
    out[0]=mxGPUCreateMxArrayOnGPU(potential.a);out[1]=mxGPUCreateMxArrayOnGPU(field.a);
    const char* names[]={"apply_seconds","persistent_geometry_bytes","workspace_bytes","output_bytes","nodes","levels","order","python_adapter_used","far_pairs","near_pairs","nonleaf_far_pairs"};out[2]=mxCreateStructMatrix(1,1,11,names);
    double stats[]={std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count(),double(c.persistent),double(workBytes),double(outputBytes),double(c.nnode),double(c.levels.size()),double(c.p),0,double(c.farPairs),double(c.nearPairs),double(c.nonleafFarPairs)};for(int i=0;i<11;++i)mxSetField(out[2],0,names[i],mxCreateDoubleScalar(stats[i]));
}
}
void mexFunction(int nl,mxArray** out,int nr,const mxArray** in){try{dispatch(nl,out,nr,in);}catch(const std::exception& e){mexErrMsgIdAndTxt("radia:cuda:fmm","%s",e.what());}}
