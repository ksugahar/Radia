// Shared mapped material energy, gradient and exact Hessian-vector action.
// NGSolve supplies E and quadrature; no FE semantics or dense tangent fallback.
#pragma once
#include "hmatrix_inverse_bh.cuh"
#include "hdiv_flat_evaluation.cuh"
namespace radia_hdiv_bh {
__global__ void materialAmplitudes(std::size_t n,const double* v,double* a){
    const auto q=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(q<n)a[q]=hypot(hypot(v[q],v[n+q]),v[2*n+q]);
}
__global__ void materialPointFields(std::size_t n,const double* v,const double* dv,const double* a,
    const double* law,const double* weights,double* field,double* action,double* energy,int* invalid){
    const auto q=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(q>=n)return;
    const double sec=law[q],dif=law[n+q],w=weights[q];
    if(!isfinite(sec)||!isfinite(dif)||!isfinite(law[2*n+q])||sec<=0||dif<=0){atomicExch(invalid,1);return;}
    double direction[3];for(int c=0;c<3;++c)direction[c]=a[q]>0?v[c*n+q]/a[q]:0.;
    radia_hmatrix_cuda::CompensatedDeviceSum dot;if(dv)for(int c=0;c<3;++c)dot.addProduct(direction[c],dv[c*n+q]);
    for(int c=0;c<3;++c){field[c*n+q]=w*sec*v[c*n+q];if(dv)action[c*n+q]=w*(sec*dv[c*n+q]+(dif-sec)*direction[c]*dot.get());}
    energy[q]=w*law[2*n+q];
}
__global__ void materialEnergyReduction(std::size_t n,const double* energy,const double* amplitude,double* result){
    if(blockIdx.x||threadIdx.x)return;radia_hmatrix_cuda::CompensatedDeviceSum sum;double maximum=0;
    for(std::size_t q=0;q<n;++q){sum.add(energy[q]);maximum=fmax(maximum,amplitude[q]);}result[0]=sum.get();result[1]=maximum;
}
class MappedMaterial {
    using DeviceMap=radia_hmatrix_cuda::DeviceMap;
    template<class T>using Buffer=radia_hmatrix_cuda::Buffer<T>;
    std::size_t nq,nd,nk,payloadBound=0,setupBytes=0,stateBytes=0,calls=0;
    std::size_t stateEpoch=0,tangentCalls=0;bool frozen=false;
    int device=0;double cap,barrier,last,initial;
    // Explicit E/E^T: the host-evaluation route and the small-case oracle only.
    DeviceMap E,ET,M;FlatEvaluationOperator flat;bool matrixFree=false;
    void applyE(const double* x,double* y,cudaStream_t stream){if(matrixFree)flat.apply(x,y,stream);else E.apply(x,y,stream);}
    void applyET(const double* x,double* y,cudaStream_t stream){if(matrixFree)flat.applyTranspose(x,y,stream);else ET.apply(x,y,stream);}
    Buffer<double> m,delta,v,dv,amplitude,samples,weights,knots,c,base,field,action,pointEnergy,gradient,tangent,totals,massDiagonal;
    Buffer<int> invalid;
    std::string routeName;radia_hmatrix_cuda::HostVector<LedgerRow> ledger;
    // Reference mass exists only on the device-table route; persistent/transient bytes are exact allocation sizes.
    bool hasMass=false;double massBound=0.;std::size_t massBytes=0,massNonzeros=0,persistentBytes=0,transientBytes=0;FlatEvaluationBudget budget;
    void allocateState(const Law& law){
        m.allocate(nd);delta.allocate(nd);v.allocate(3*nq);dv.allocate(3*nq);amplitude.allocate(nq);samples.allocate(3*nq);
        knots.upload(law.knots);c.upload(law.coefficients);base.upload(law.integralBase);field.allocate(3*nq);action.allocate(3*nq);
        pointEnergy.allocate(nq);gradient.allocate(nd);tangent.allocate(nd);totals.allocate(2);invalid.allocate(1);
        ledger.push_back({"inverse_bh_law","H2D",(law.knots.size()+law.coefficients.size()+law.integralBase.size())*sizeof(double)});
    }
public:
    MappedMaterial(const MaterialInput& in,const Law& law):nq(in.points),nd(in.dofs),nk(law.knots.size()),cap(law.cap),barrier(law.barrierSlope),last(law.lastField),initial(law.initialReluctivity){
    using namespace radia_hmatrix_cuda;
    if(!nq||nq>65536||!nd||nd>16384||in.weights.size()!=nq||in.offsets.size()!=3*nq+1||in.indices.size()!=in.values.size()||in.values.size()>4000000)
        throw std::invalid_argument("Mapped material shape or payload budget exceeded");
    for(double w:in.weights)if(!std::isfinite(w)||w<=0)throw std::invalid_argument("Material measures must be positive finite");
    validateInverseBH(law,radia_hmatrix_cuda::HostVector<double>{0.});
    SparseMap e;e.rows=3*nq;e.columns=nd;e.offsets.assign(in.offsets.begin(),in.offsets.end());e.indices.assign(in.indices.begin(),in.indices.end());e.values.assign(in.values.begin(),in.values.end());e.validate();
    const auto transpose=e.transpose();
    std::size_t bytes=0;auto add=[&](std::size_t count,std::size_t width){const auto size=checkedProduct(count,width);if(size>(std::size_t(512)<<20)-bytes)throw std::invalid_argument("Mapped material device payload exceeds 512 MiB");bytes+=size;};
    add(e.offsets.size()+e.indices.size()+transpose.offsets.size()+transpose.indices.size(),sizeof(std::size_t));
    add(e.values.size()+transpose.values.size()+law.knots.size()+law.coefficients.size()+law.integralBase.size()+18*nq+4*nd+2,sizeof(double));add(1,sizeof(int));

        payloadBound=bytes;check(cudaGetDevice(&device));E.upload(e);ET.upload(transpose);
        setupBytes=(e.offsets.size()+e.indices.size()+transpose.offsets.size()+transpose.indices.size())*sizeof(std::size_t)+
            (e.values.size()+transpose.values.size()+in.weights.size()+law.knots.size()+law.coefficients.size()+law.integralBase.size())*sizeof(double);
        routeName="host-evaluation";
        ledger.push_back({"evaluation_csr_host_built","H2D",(e.offsets.size()+e.indices.size()+transpose.offsets.size()+transpose.indices.size())*sizeof(std::size_t)+
            (e.values.size()+transpose.values.size())*sizeof(double)});
        weights.upload(in.weights);ledger.push_back({"material_weights_host","H2D",in.weights.size()*sizeof(double)});allocateState(law);
    }
    // Matrix-free E/E^T from NGSolve reference tables resident on the device, with the
    // point measures and reference mass built there. explicitOracle additionally stores
    // the explicit E/E^T, applies them and builds M from them: small-case diagnostic only.
    MappedMaterial(const FlatEvaluationInput& in,const Law& law,bool explicitOracle=false):nq(in.points),nd(in.dofs),nk(law.knots.size()),cap(law.cap),barrier(law.barrierSlope),last(law.lastField),initial(law.initialReluctivity){
        using namespace radia_hmatrix_cuda;
        validateFlatEvaluation(in);validateInverseBH(law,radia_hmatrix_cuda::HostVector<double>{0.});
        // Dimension budget plus CUB scratch against live free memory; no fixed payload cap.
        payloadBound=requireFlatEvaluationDevice(in,law.knots.size()+law.coefficients.size()+law.integralBase.size(),explicitOracle,budget);
        check(cudaGetDevice(&device));
        FlatEvaluationBuilder builder;
        if(explicitOracle)builder.build(in,flat,weights,M,massDiagonal,&E,&ET);else builder.build(in,flat,weights,M,massDiagonal);
        matrixFree=!explicitOracle;routeName=explicitOracle?"device-flat-tet-explicit-oracle":"device-flat-tet-matrix-free";
        hasMass=true;massBound=builder.massBound;transientBytes=builder.transientBytes;
        massNonzeros=builder.massNonzeros;massBytes=(nd+1+massNonzeros)*sizeof(std::size_t)+(massNonzeros+nd)*sizeof(double);
        for(const auto& t:builder.transfers){ledger.push_back({t.payload,t.direction,t.bytes});if(std::string(t.direction)=="H2D")setupBytes+=t.bytes;}
        allocateState(law);setupBytes+=ledger.back().bytes;
        for(auto* map:{&E,&ET,&M})persistentBytes+=(map->offsets.size+map->indices.size)*sizeof(std::size_t)+map->values.size*sizeof(double);
        persistentBytes+=flat.bytes();
        for(auto* buffer:{&m,&delta,&v,&dv,&amplitude,&samples,&weights,&knots,&c,&base,&field,&action,&pointEnergy,&gradient,&tangent,&totals,&massDiagonal})
            persistentBytes+=buffer->size*sizeof(double);
        persistentBytes+=invalid.size*sizeof(int);
    }
    radia_hmatrix_cuda::MassOperand massOperand(){
        if(!hasMass)throw std::invalid_argument("Material context has no device-built reference mass");
        return {&M,massDiagonal.data,massBound,nd,massBytes,true};
    }
    MaterialFootprint footprint()const{return {persistentBytes,transientBytes};}
    bool deviceMass()const{return hasMass;}
    MassAction massAction(const radia_hmatrix_cuda::HostVector<double>& x){
        using namespace radia_hmatrix_cuda;
        if(!hasMass)throw std::invalid_argument("Material context has no device-built reference mass");
        if(x.size()!=nd)throw std::invalid_argument("Reference mass action dimension mismatch");
        for(double v:x)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite reference mass input");
        int current;check(cudaGetDevice(&current));if(current!=device)throw std::runtime_error("Mapped material belongs to another GPU");
        check(cudaMemcpy(delta.data,x.data(),nd*sizeof(double),cudaMemcpyHostToDevice));M.apply(delta.data,tangent.data,nullptr);
        MassAction out;out.y.resize(nd);check(cudaMemcpy(out.y.data(),tangent.data,nd*sizeof(double),cudaMemcpyDeviceToHost));
        out.nonzeros=massNonzeros;out.bound=massBound;return out;
    }
    // Diagnostic y=E x (transpose=false, x of size nd) or g=E^T y (transpose=true, y of size 3*nq).
    radia_hmatrix_cuda::HostVector<double> evaluationAction(const radia_hmatrix_cuda::HostVector<double>& x,bool transpose){
        using namespace radia_hmatrix_cuda;
        const auto in=transpose?3*nq:nd,out=transpose?nd:3*nq;
        if(x.size()!=in)throw std::invalid_argument("Evaluation action dimension mismatch");
        for(double value:x)if(!std::isfinite(value))throw std::invalid_argument("Nonfinite evaluation action input");
        int current;check(cudaGetDevice(&current));if(current!=device)throw std::runtime_error("Mapped material belongs to another GPU");
        double* input=transpose?field.data:delta.data;double* output=transpose?gradient.data:v.data;
        check(cudaMemcpy(input,x.data(),in*sizeof(double),cudaMemcpyHostToDevice));
        if(transpose)applyET(input,output,nullptr);else applyE(input,output,nullptr);
        frozen=false;HostVector<double> y(out);check(cudaMemcpy(y.data(),output,out*sizeof(double),cudaMemcpyDeviceToHost));return y;
    }
    const std::string& route()const{return routeName;}
    const radia_hmatrix_cuda::HostVector<LedgerRow>& transfers()const{return ledger;}
    std::size_t dofs()const{return nd;}
    MaterialStats stats()const{return {nq,nd,payloadBound,setupBytes,stateBytes,calls,1,stateEpoch,tangentCalls,false,frozen};}
    template<class StateAllocator,class DirectionAllocator=std::allocator<double>>
    MaterialEvaluation evaluate(const std::vector<double,StateAllocator>& state,const std::vector<double,DirectionAllocator>& direction){
        using namespace radia_hmatrix_cuda;
        if(calls==std::numeric_limits<std::size_t>::max())throw std::overflow_error("Material evaluation counter exhausted");
        frozen=false;if(stateEpoch==std::numeric_limits<std::size_t>::max())throw std::overflow_error("Material state epoch exhausted");++stateEpoch;
        int current;check(cudaGetDevice(&current));if(current!=device)throw std::runtime_error("Mapped material belongs to another GPU");
        if(state.size()!=nd||(!direction.empty()&&direction.size()!=nd))throw std::invalid_argument("Material state or direction dimension mismatch");
        for(double x:state)if(!std::isfinite(x))throw std::invalid_argument("Nonfinite material state");
        for(double x:direction)if(!std::isfinite(x))throw std::invalid_argument("Nonfinite material direction");
        const bool hasDirection=!direction.empty();check(cudaMemcpy(m.data,state.data(),nd*sizeof(double),cudaMemcpyHostToDevice));
        if(hasDirection)check(cudaMemcpy(delta.data,direction.data(),nd*sizeof(double),cudaMemcpyHostToDevice));check(cudaStreamSynchronize(nullptr));const auto uploaded=(state.size()+direction.size())*sizeof(double);
        if(uploaded>std::numeric_limits<std::size_t>::max()-stateBytes)throw std::overflow_error("Material upload counter exhausted");stateBytes+=uploaded;
        check(cudaMemset(invalid.data,0,sizeof(int)));applyE(m.data,v.data,nullptr);if(hasDirection)applyE(delta.data,dv.data,nullptr);
        const unsigned blocks=unsigned((nq+127)/128);materialAmplitudes<<<blocks,128>>>(nq,v.data,amplitude.data);check(cudaGetLastError());
        inverseLawKernel<<<blocks,128>>>(amplitude.data,nq,knots.data,nk,c.data,base.data,cap,barrier,last,initial,samples.data);check(cudaGetLastError());
        materialPointFields<<<blocks,128>>>(nq,v.data,hasDirection?dv.data:nullptr,amplitude.data,samples.data,weights.data,field.data,action.data,pointEnergy.data,invalid.data);check(cudaGetLastError());
        int failed=0;check(cudaMemcpy(&failed,invalid.data,sizeof(int),cudaMemcpyDeviceToHost));if(failed)throw std::runtime_error("Invalid mapped inverse BH tangent or energy");
        applyET(field.data,gradient.data,nullptr);if(hasDirection)applyET(action.data,tangent.data,nullptr);
        materialEnergyReduction<<<1,1>>>(nq,pointEnergy.data,amplitude.data,totals.data);check(cudaGetLastError());
        MaterialEvaluation out;out.gradient.resize(nd);check(cudaMemcpy(out.gradient.data(),gradient.data,nd*sizeof(double),cudaMemcpyDeviceToHost));
        if(hasDirection){out.tangentAction.resize(nd);check(cudaMemcpy(out.tangentAction.data(),tangent.data,nd*sizeof(double),cudaMemcpyDeviceToHost));}
        double scalar[2];check(cudaMemcpy(scalar,totals.data,2*sizeof(double),cudaMemcpyDeviceToHost));out.energy=scalar[0];out.maximumMagnetization=scalar[1];out.ownedDeviceBytes=payloadBound;
        for(double x:out.gradient)if(!std::isfinite(x))throw std::runtime_error("Nonfinite mapped material gradient");
        for(double x:out.tangentAction)if(!std::isfinite(x))throw std::runtime_error("Nonfinite mapped material tangent action");
        if(!std::isfinite(out.energy)||!std::isfinite(out.maximumMagnetization))throw std::runtime_error("Nonfinite mapped material energy");++calls;return out;
    }
    template<class Allocator>MaterialEvaluation linearize(const std::vector<double,Allocator>& state){auto out=evaluate(state,{});frozen=true;out.linearizationEpoch=stateEpoch;return out;}
    void applyFrozenDevice(const double* direction,double* output,cudaStream_t stream,std::size_t expectedEpoch){
        using namespace radia_hmatrix_cuda;
        if(!frozen||!expectedEpoch||expectedEpoch!=stateEpoch)throw std::invalid_argument("Stale material linearization epoch");
        int current;check(cudaGetDevice(&current));if(current!=device)throw std::runtime_error("Mapped material belongs to another GPU");
        if(tangentCalls==std::numeric_limits<std::size_t>::max())throw std::overflow_error("Material action counter exhausted");
        applyE(direction,dv.data,stream);
        materialPointFields<<<unsigned((nq+127)/128),128,0,stream>>>(nq,v.data,dv.data,amplitude.data,samples.data,weights.data,field.data,action.data,pointEnergy.data,invalid.data);check(cudaGetLastError());
        applyET(action.data,output,stream);++tangentCalls;
    }
    radia_hmatrix_cuda::HostVector<double> tangentAction(const radia_hmatrix_cuda::HostVector<double>& direction,std::size_t expectedEpoch){
        using namespace radia_hmatrix_cuda;
        if(direction.size()!=nd)throw std::invalid_argument("Material direction dimension mismatch");
        for(double x:direction)if(!std::isfinite(x))throw std::invalid_argument("Nonfinite material direction");
        if(!frozen||!expectedEpoch||expectedEpoch!=stateEpoch)throw std::invalid_argument("Stale material linearization epoch");
        check(cudaMemcpy(delta.data,direction.data(),nd*sizeof(double),cudaMemcpyHostToDevice));check(cudaStreamSynchronize(nullptr));
        const auto uploaded=nd*sizeof(double);if(uploaded>std::numeric_limits<std::size_t>::max()-stateBytes)throw std::overflow_error("Material upload counter exhausted");stateBytes+=uploaded;
        applyFrozenDevice(delta.data,tangent.data,nullptr,expectedEpoch);HostVector<double> out(nd,HostFaultAllocator<double>("bh.tangent-download"));check(cudaMemcpy(out.data(),tangent.data,nd*sizeof(double),cudaMemcpyDeviceToHost));
        for(double x:out)if(!std::isfinite(x))throw std::runtime_error("Nonfinite frozen material tangent action");return out;
    }
};
inline MaterialEvaluation evaluateMappedMaterial(const MaterialInput& input,const Law& law,const radia_hmatrix_cuda::HostVector<double>& state,const radia_hmatrix_cuda::HostVector<double>& direction){
    MappedMaterial context(input,law);return context.evaluate(state,direction);
}
}
