// Resident constitutive iteration for an explicitly assembled complete block.
// The assembler owns coil/body reaction and impedance-entry coefficients.
#pragma once
#include "ih_temperature_material.cuh"
#include "ih_complex_solver.cuh"
namespace radia_ih_cuda {
__global__ void esimOperatorValues(ComplexCSRView entryMap,const cuDoubleComplex* base,
        const ESIMCellResult* cells,cuDoubleComplex* values,int* invalid) {
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=entryMap.rows)return;
    auto value=base[i];
    for(int k=entryMap.rowptr[i];k<entryMap.rowptr[i+1];++k){const auto& cell=cells[entryMap.columns[k]];
        if(cell.status!=1||cell.trueResidual>1e-6)atomicExch(invalid,1);
        value=cuCadd(value,cuCmul(entryMap.values[k],cell.impedance));}
    if(!cellFinite(value))atomicExch(invalid,1);values[i]=value;
}
__global__ void esimFieldAmplitude(ComplexCSRView projection,const cuDoubleComplex* state,
        std::size_t samples,unsigned components,double floor,double* field,int* invalid) {
    const auto e=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(e>=samples)return;
    double squared=0.;for(unsigned c=0;c<components;++c){const auto row=e*components+c;auto h=make_cuDoubleComplex(0.,0.);
        for(int k=projection.rowptr[row];k<projection.rowptr[row+1];++k)h=cuCadd(h,cuCmul(projection.values[k],state[projection.columns[k]]));
        squared+=cuCreal(h)*cuCreal(h)+cuCimag(h)*cuCimag(h);}
    const double value=sqrt(squared);if(!isfinite(value))atomicExch(invalid,1);field[e]=fmax(value,floor);
}
__global__ void esimRelaxFields(double* field,const double* next,std::size_t n,double relaxation){
    const auto e=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(e<n)field[e]=(1.-relaxation)*field[e]+relaxation*next[e];
}
__global__ void esimLossWeights(const double* base,const ESIMCellResult* cells,std::size_t n,unsigned components,double* weights){
    const auto e=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(e<n)for(unsigned c=0;c<components;++c)weights[e*components+c]=base[e*components+c]*cuCreal(cells[e].impedance);
}
// Bounded research blocks use one ordered scalar reduction, not a vector download.
__global__ void esimCertify(const cuDoubleComplex* ax,const cuDoubleComplex* rhs,std::size_t n,
        const double* field,const double* next,const ESIMCellResult* cells,std::size_t samples,double* result,int* invalid){
    if(blockIdx.x||threadIdx.x)return;double rr=0.,bb=0.,change=0.;
    for(std::size_t i=0;i<n;++i){const auto r=cuCsub(ax[i],rhs[i]);rr+=cuCreal(r)*cuCreal(r)+cuCimag(r)*cuCimag(r);bb+=cuCreal(rhs[i])*cuCreal(rhs[i])+cuCimag(rhs[i])*cuCimag(rhs[i]);}
    for(std::size_t e=0;e<samples;++e){change=fmax(change,fabs(field[e]-next[e])/fmax(field[e],next[e]));if(cells[e].status!=1||!isfinite(cells[e].trueResidual)||cells[e].trueResidual>1e-6)atomicExch(invalid,1);}
    result[0]=bb==0.?sqrt(rr):sqrt(rr/bb);result[1]=change;if(!isfinite(result[0])||!isfinite(change))atomicExch(invalid,1);
}
struct ESIMCouplingView {
    RealCSRView temperatureEvaluation;TemperatureMaterialTable table;
    ComplexCSRView entryMap,fieldProjection,matrix;
    const cuDoubleComplex* baseValues;
    cuDoubleComplex* mutableValues;
    const ESIMCellParameters* geometry;
    const double *initialFields,*baseLossWeights;
    double* mutableLossWeights;
    std::size_t elements;unsigned components;int workspaceNodes,maxIterations;
    double tolerance,relaxation,fieldFloor;
};
class ESIMCoupling {
    ESIMCouplingView v_;
    ESIMCellParameters* parameters_=nullptr;ESIMCellResult* cells_=nullptr;ESIMCellResult* acceptedCells_=nullptr;
    double *fields_=nullptr,*next_=nullptr,*temperature_=nullptr,*coefficients_=nullptr,*scalar_=nullptr;
    cuDoubleComplex *workspace_=nullptr,*ax_=nullptr;int* invalid_=nullptr;
    std::size_t bytes_=0;unsigned candidateIterations_=0,acceptedIterations_=0;
    double candidateResidual_=0.,acceptedResidual_=0.,candidateFieldResidual_=0.,acceptedFieldResidual_=0.;
    template<class T>void allocate(T*& p,std::size_t n){cudaCheck(cudaMalloc(&p,n*sizeof(T)));bytes_+=n*sizeof(T);}
    void release()noexcept{for(auto* p:{fields_,next_,temperature_,coefficients_,scalar_})if(p)cudaFree(p);if(parameters_)cudaFree(parameters_);if(cells_)cudaFree(cells_);if(acceptedCells_)cudaFree(acceptedCells_);if(workspace_)cudaFree(workspace_);if(ax_)cudaFree(ax_);if(invalid_)cudaFree(invalid_);}
    void cells(cudaStream_t stream,const double* field){finiteElementESIMCells<<<unsigned((v_.elements+31)/32),32,0,stream>>>(parameters_,field,v_.elements,v_.workspaceNodes,workspace_,cells_);}
    void matrix(cudaStream_t stream){esimOperatorValues<<<unsigned((v_.entryMap.rows+127)/128),128,0,stream>>>(v_.entryMap,v_.baseValues,cells_,v_.mutableValues,invalid_);}
    void check(cudaStream_t stream){int bad;cudaCheck(cudaMemcpyAsync(&bad,invalid_,sizeof(int),cudaMemcpyDeviceToHost,stream));cudaCheck(cudaStreamSynchronize(stream));if(bad)throw std::runtime_error("Element ESIM material/operator candidate rejected");}
public:
    explicit ESIMCoupling(ESIMCouplingView v):v_(v){
        if(!v.elements||v.elements>4096||!v.matrix.rows||v.matrix.rows>4096||v.workspaceNodes<3||v.workspaceNodes>4096||
           v.table.temperatures<2||v.table.knotsCount<0||v.table.knotsCount==1||v.maxIterations<2||v.maxIterations>10000||
           !std::isfinite(v.tolerance)||v.tolerance<=0.||v.tolerance>1e-6||!std::isfinite(v.relaxation)||v.relaxation<=0.||v.relaxation>1.||
           !std::isfinite(v.fieldFloor)||v.fieldFloor<=0.||v.temperatureEvaluation.rows!=v.elements||v.fieldProjection.rows!=v.elements*v.components||!v.components)
            throw std::invalid_argument("Invalid bounded element ESIM coupling contract");
        const auto coefficientCount=v.table.knotsCount?std::size_t(v.table.knotsCount-1)*4*v.elements:0;
        if(coefficientCount>(2ULL<<30)/sizeof(double))throw std::overflow_error("ESIM coefficient workspace exceeds direct payload budget");
        const auto need=v.elements*(sizeof(ESIMCellParameters)+2*sizeof(ESIMCellResult)+3*sizeof(double)+std::size_t(v.workspaceNodes)*4*sizeof(cuDoubleComplex))+coefficientCount*sizeof(double)+v.matrix.rows*sizeof(cuDoubleComplex)+2*sizeof(double)+sizeof(int);
        std::size_t free,total;cudaCheck(cudaMemGetInfo(&free,&total));const auto reserve=std::max<std::size_t>(2ULL<<30,total/5);
        if(free<=reserve||need>free-reserve)throw std::runtime_error("Element ESIM context exceeds GPU reserve");
        try{allocate(parameters_,v.elements);allocate(cells_,v.elements);allocate(acceptedCells_,v.elements);cudaCheck(cudaMemset(acceptedCells_,0,v.elements*sizeof(ESIMCellResult)));allocate(fields_,v.elements);allocate(next_,v.elements);allocate(temperature_,v.elements);if(coefficientCount)allocate(coefficients_,coefficientCount);allocate(workspace_,v.elements*std::size_t(v.workspaceNodes)*4);allocate(ax_,v.matrix.rows);allocate(invalid_,1);allocate(scalar_,2);}catch(...){release();throw;}
    }
    ESIMCoupling(const ESIMCoupling&)=delete;ESIMCoupling& operator=(const ESIMCoupling&)=delete;
    ~ESIMCoupling(){release();}
    void abandonAfterDeviceReset()noexcept{parameters_=nullptr;cells_=acceptedCells_=nullptr;fields_=next_=temperature_=coefficients_=scalar_=nullptr;workspace_=ax_=nullptr;invalid_=nullptr;}
    std::size_t ownedBytes()const{return bytes_;}
    void validateTemperature(cudaStream_t stream,const double* temperature){
        cudaCheck(cudaMemsetAsync(invalid_,0,sizeof(int),stream));
        validateElementTemperatureRange<<<unsigned((v_.elements+127)/128),128,0,stream>>>(v_.temperatureEvaluation,temperature,v_.table.kelvin,v_.table.temperatures,invalid_);
        check(stream);
    }
    void prepare(cudaStream_t stream,const double* temperature,double){
        cudaCheck(cudaMemsetAsync(invalid_,0,sizeof(int),stream));
        updateElementESIMMaterials<<<unsigned((v_.elements+127)/128),128,0,stream>>>(v_.temperatureEvaluation,temperature,v_.table,v_.geometry,parameters_,temperature_,coefficients_,invalid_);
        cudaCheck(cudaMemcpyAsync(fields_,v_.initialFields,v_.elements*sizeof(double),cudaMemcpyDeviceToDevice,stream));check(stream);
    }
    ComplexSolver::Result solve(ComplexSolver& solver,double tolerance,int maxIterations){
        const auto stream=solver.stream();ComplexSolver::Result linear{};unsigned totalIterations=0;
        for(int iteration=0;iteration<v_.maxIterations;++iteration){
            cells(stream,fields_);matrix(stream);check(stream);linear=solver.solve(tolerance,maxIterations);totalIterations+=linear.iterations;
            esimFieldAmplitude<<<unsigned((v_.elements+127)/128),128,0,stream>>>(v_.fieldProjection,solver.solution(),v_.elements,v_.components,v_.fieldFloor,next_,invalid_);
            // Certify with direct cells at the returned total field, not lagged Z.
            cells(stream,next_);matrix(stream);check(stream);
            complexMatvec<<<unsigned((v_.matrix.rows+127)/128),128,0,stream>>>(v_.matrix,solver.solution(),ax_);
            esimCertify<<<1,1,0,stream>>>(ax_,solver.rhs(),v_.matrix.rows,fields_,next_,cells_,v_.elements,scalar_,invalid_);
            double residual[2];cudaCheck(cudaMemcpyAsync(residual,scalar_,sizeof(residual),cudaMemcpyDeviceToHost,stream));check(stream);
            if(iteration>0&&residual[0]<=v_.tolerance&&residual[1]<=v_.tolerance){
                esimLossWeights<<<unsigned((v_.elements+127)/128),128,0,stream>>>(v_.baseLossWeights,cells_,v_.elements,v_.components,v_.mutableLossWeights);
                candidateIterations_=unsigned(iteration+1);candidateResidual_=residual[0];candidateFieldResidual_=residual[1];return {totalIterations,linear.trueRelativeResidual};
            }
            esimRelaxFields<<<unsigned((v_.elements+127)/128),128,0,stream>>>(fields_,next_,v_.elements,v_.relaxation);
        }
        throw std::runtime_error("Element ESIM complete-block constitutive iteration did not converge; accepted state preserved");
    }
    void commit(){cudaCheck(cudaMemcpy(acceptedCells_,cells_,v_.elements*sizeof(ESIMCellResult),cudaMemcpyDeviceToDevice));acceptedIterations_=candidateIterations_;acceptedResidual_=candidateResidual_;acceptedFieldResidual_=candidateFieldResidual_;}
    void reset(){cudaCheck(cudaMemset(acceptedCells_,0,v_.elements*sizeof(ESIMCellResult)));acceptedIterations_=candidateIterations_=0;acceptedResidual_=candidateResidual_=acceptedFieldResidual_=candidateFieldResidual_=0.;}
    std::vector<ESIMCellResult> acceptedCells()const{std::vector<ESIMCellResult> out(v_.elements);cudaCheck(cudaMemcpy(out.data(),acceptedCells_,out.size()*sizeof(ESIMCellResult),cudaMemcpyDeviceToHost));return out;}
    unsigned iterations()const{return acceptedIterations_;}double residual()const{return acceptedResidual_;}double fieldResidual()const{return acceptedFieldResidual_;}
};
} // namespace radia_ih_cuda
