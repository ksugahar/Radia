// Persistent coupled IH state. Operator storage must outlive this context.
#pragma once
#include "ih_complex_solver.cuh"
#include "ih_thermal_kernels.cuh"
#include <memory>
#include <climits>
#include <vector>
namespace radia_ih_cuda {
struct CoupledOperators {
    ComplexCSRView eddy,lossProjection,thermalSystem;
    RealCSRView stiffness,convection,heatProjection;
    const cuDoubleComplex* rhsPerAmp;
    const double *lossWeights,*constantCoefficients,*initialTemperature;
    std::size_t heatSamples;
    unsigned lossComponents;
    double dt,conductivity,convectionScale;
    // Dynamic ESIM/material owners must supply their real update function.
    // It runs before each solve with the accepted temperature on device.
    bool dynamicMaterial;
    void (*updateMaterial)(void*,cudaStream_t,const double*,double);
    void* materialOwner;
};
class CoupledState {
    CoupledOperators op_;
    std::unique_ptr<ComplexSolver> eddy_,thermal_;
    double *heat_=nullptr,*candidateHeat_=nullptr,*load_=nullptr,*temperature_=nullptr,*candidate_=nullptr,*controls_=nullptr;
    int* invalid_=nullptr;
    cudaEvent_t heatReady_=nullptr;
    std::size_t steps_=0;
    void release()noexcept {
        // Solver destructors synchronize their own streams before state release.
        thermal_.reset();eddy_.reset();
        if(heatReady_)cudaEventDestroy(heatReady_);
        for(auto* p:{heat_,candidateHeat_,load_,temperature_,candidate_,controls_})if(p)cudaFree(p);
        if(invalid_)cudaFree(invalid_);
    }
public:
    CoupledState(const CoupledState&)=delete;
    CoupledState& operator=(const CoupledState&)=delete;
    explicit CoupledState(CoupledOperators op):op_(op) {
        if(!op.eddy.rows||op.eddy.rows>INT_MAX||!op.thermalSystem.rows||op.thermalSystem.rows>INT_MAX||
           !op.heatSamples||op.heatSamples>INT_MAX||!op.lossComponents||
           op.lossComponents>INT_MAX/op.heatSamples||
           (op.convection.values&&op.convection.rows!=op.thermalSystem.rows)||op.lossProjection.rows!=op.heatSamples*op.lossComponents||
           op.heatProjection.rows!=op.thermalSystem.rows||op.stiffness.rows!=op.thermalSystem.rows||
           !std::isfinite(op.dt)||op.dt<=0.||!std::isfinite(op.conductivity)||op.conductivity<0.||
           !std::isfinite(op.convectionScale)||op.convectionScale<0.||
           !op.rhsPerAmp||!op.lossWeights||!op.constantCoefficients||!op.initialTemperature)
            throw std::invalid_argument("Invalid coupled IH operator dimensions or parameters");
        if(op.dynamicMaterial&&!op.updateMaterial)throw std::invalid_argument("Dynamic ESIM/material update is missing; frozen substitution is forbidden");
        std::size_t available,total;cudaCheck(cudaMemGetInfo(&available,&total));
        const auto reserve=std::max<std::size_t>(2ULL<<30,total/5);
        // Includes solver workspaces, state and controls, not caller-owned operators.
        const auto need=(op.eddy.rows*9+op.thermalSystem.rows*9+2)*sizeof(cuDoubleComplex)+
            (2*op.heatSamples+3*op.thermalSystem.rows+2)*sizeof(double)+2*sizeof(BiCGState)+sizeof(int);
        if(available<=reserve||need>available-reserve)throw std::runtime_error("IH resident state exceeds GPU memory budget");
        try {
            eddy_=std::make_unique<ComplexSolver>(op.eddy);thermal_=std::make_unique<ComplexSolver>(op.thermalSystem);
            cudaCheck(cudaMalloc(&heat_,op.heatSamples*sizeof(double)));
            cudaCheck(cudaMalloc(&candidateHeat_,op.heatSamples*sizeof(double)));
            cudaCheck(cudaMalloc(&load_,op.thermalSystem.rows*sizeof(double)));
            cudaCheck(cudaMalloc(&temperature_,op.thermalSystem.rows*sizeof(double)));
            cudaCheck(cudaMalloc(&candidate_,op.thermalSystem.rows*sizeof(double)));
            cudaCheck(cudaMalloc(&controls_,2*sizeof(double)));cudaCheck(cudaMalloc(&invalid_,sizeof(int)));
            cudaCheck(cudaEventCreateWithFlags(&heatReady_,cudaEventDisableTiming));reset();
        }catch(...){release();throw;}
    }
    ~CoupledState(){release();}
    void reset() {
        cudaCheck(cudaMemcpyAsync(temperature_,op_.initialTemperature,op_.thermalSystem.rows*sizeof(double),cudaMemcpyDeviceToDevice,thermal_->stream()));
        cudaCheck(cudaMemsetAsync(heat_,0,op_.heatSamples*sizeof(double),thermal_->stream()));
        cudaCheck(cudaStreamSynchronize(thermal_->stream()));steps_=0;
    }
    // Called only at the Simulink discrete Update boundary. Outputs observes
    // accepted state; no thermal vector is returned merely to re-upload it.
    void update(double current,double ambient,double tolerance=1e-10,int maxIterations=2000) {
        if(!std::isfinite(current)||!std::isfinite(ambient))throw std::invalid_argument("Nonfinite IH controls");
        const double controls[2]={current,ambient};
        cudaCheck(cudaMemcpyAsync(controls_,controls,sizeof(controls),cudaMemcpyHostToDevice,eddy_->stream()));
        if(op_.updateMaterial)op_.updateMaterial(op_.materialOwner,eddy_->stream(),temperature_,current);
        scaleDrive<<<unsigned((op_.eddy.rows+127)/128),128,0,eddy_->stream()>>>(op_.rhsPerAmp,controls_,op_.eddy.rows,eddy_->rhs());
        eddy_->solve(tolerance,maxIterations);
        quadraticLoss<<<unsigned((op_.heatSamples+127)/128),128,0,eddy_->stream()>>>(op_.lossProjection,eddy_->solution(),op_.lossWeights,op_.heatSamples,op_.lossComponents,candidateHeat_);
        cudaCheck(cudaEventRecord(heatReady_,eddy_->stream()));
        cudaCheck(cudaStreamWaitEvent(thermal_->stream(),heatReady_,0));
        const auto blocks=unsigned((op_.thermalSystem.rows+127)/128);
        thermalLoad<<<blocks,128,0,thermal_->stream()>>>(op_.heatProjection,candidateHeat_,load_);
        thermalIncrementRHS<<<blocks,128,0,thermal_->stream()>>>(op_.stiffness,op_.convection,load_,temperature_,op_.constantCoefficients,controls_+1,op_.dt,op_.conductivity,op_.convectionScale,thermal_->rhs());
        thermal_->solve(tolerance,maxIterations);
        cudaCheck(cudaMemsetAsync(invalid_,0,sizeof(int),thermal_->stream()));
        thermalCandidate<<<blocks,128,0,thermal_->stream()>>>(temperature_,thermal_->solution(),candidate_,invalid_,op_.thermalSystem.rows);
        int invalid=0;cudaCheck(cudaMemcpyAsync(&invalid,invalid_,sizeof(int),cudaMemcpyDeviceToHost,thermal_->stream()));
        cudaCheck(cudaStreamSynchronize(thermal_->stream()));
        if(invalid)throw std::runtime_error("IH candidate temperature rejected; accepted state preserved");
        std::swap(temperature_,candidate_);std::swap(heat_,candidateHeat_);++steps_;
    }
    std::vector<double> temperature() {
        std::vector<double> out(op_.thermalSystem.rows);
        cudaCheck(cudaMemcpyAsync(out.data(),temperature_,out.size()*sizeof(double),cudaMemcpyDeviceToHost,thermal_->stream()));
        cudaCheck(cudaStreamSynchronize(thermal_->stream()));return out;
    }
    std::vector<double> heat() {
        std::vector<double> out(op_.heatSamples);
        cudaCheck(cudaMemcpyAsync(out.data(),heat_,out.size()*sizeof(double),cudaMemcpyDeviceToHost,thermal_->stream()));
        cudaCheck(cudaStreamSynchronize(thermal_->stream()));return out;
    }
    std::size_t steps()const{return steps_;}
};
} // namespace radia_ih_cuda
