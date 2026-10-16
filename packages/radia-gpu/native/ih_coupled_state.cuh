// Persistent coupled IH state. Operator storage must outlive this context.
#pragma once
#include "ih_complex_solver.cuh"
#include "ih_thermal_kernels.cuh"
#include "ih_work_snapshot.cuh"
#include <memory>
#include <climits>
#include <vector>
#include <type_traits>
#include <utility>
namespace radia_ih_cuda {
struct CoupledOperators {
    ComplexCSRView eddy,lossProjection,thermalSystem;
    ComplexPreconditioner eddyPreconditioner{};
    RealCSRView stiffness,convection,heatProjection;
    const cuDoubleComplex* rhsPerAmp;
    const double *lossWeights,*constantCoefficients,*initialTemperature;
    std::size_t heatSamples;
    unsigned lossComponents;
    double dt,conductivity,convectionScale;
    // Dynamic ESIM/material owners must supply their real update function.
    // It runs before each solve with the accepted temperature on device.
    WorkSnapshot* snapshot=nullptr;
    void (*stageWorkMaterial)(void*,WorkSnapshot&,cudaStream_t,WorkSnapshot::Meta&)=nullptr;
    bool dynamicMaterial;
    void (*updateMaterial)(void*,cudaStream_t,const double*,double);
    void* materialOwner;
    ComplexSolver::Result (*solveMaterial)(void*,ComplexSolver&,double,int)=nullptr;
    void (*stageMaterial)(void*,cudaStream_t,bool)=nullptr;
    void (*commitMaterial)(void*) noexcept=nullptr;

    void (*validateCandidateTemperature)(void*,cudaStream_t,const double*)=nullptr;
};
class CoupledState {
    CoupledOperators op_;
    std::unique_ptr<ComplexSolver> eddy_,thermal_;
    double *heat_=nullptr,*candidateHeat_=nullptr,*load_=nullptr,*temperature_=nullptr,*candidate_=nullptr,*controls_=nullptr;
    cuDoubleComplex* acceptedElectromagnetic_=nullptr;
    cuDoubleComplex* candidateElectromagnetic_=nullptr;
    int* invalid_=nullptr;
    cudaEvent_t heatReady_=nullptr;
    std::size_t steps_=0;
#ifdef RADIA_IH_TRANSACTION_TEST
    unsigned testFailureStage_=0;
#endif
    void testPoint(unsigned stage){
#ifdef RADIA_IH_TRANSACTION_TEST
        if(testFailureStage_==stage){testFailureStage_=0;throw std::runtime_error("Injected IH transaction rejection");}
#else
        (void)stage;
#endif
    }
    ComplexSolver::Result lastEddy_{},lastThermal_{};
    void release()noexcept {
        // Solver destructors synchronize their own streams before state release.
        thermal_.reset();eddy_.reset();
        if(heatReady_)cudaEventDestroy(heatReady_);
        for(auto* p:{heat_,candidateHeat_,load_,temperature_,candidate_,controls_})if(p)cudaFree(p);
        if(acceptedElectromagnetic_)cudaFree(acceptedElectromagnetic_);if(candidateElectromagnetic_)cudaFree(candidateElectromagnetic_);
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
        const auto need=((op.eddy.rows)*(op.eddyPreconditioner.precision==PreconditionerPrecision::None?10:12)+op.thermalSystem.rows*9+2)*sizeof(cuDoubleComplex)+
            (2*op.heatSamples+3*op.thermalSystem.rows+2)*sizeof(double)+2*sizeof(BiCGState)+sizeof(int)+op.eddy.rows*sizeof(cuDoubleComplex);
        if(available<=reserve||need>available-reserve)throw std::runtime_error("IH resident state exceeds GPU memory budget");
        try {
            eddy_=std::make_unique<ComplexSolver>(op.eddy,op.eddyPreconditioner);thermal_=std::make_unique<ComplexSolver>(op.thermalSystem);
            cudaCheck(cudaMalloc(&candidateElectromagnetic_,op.eddy.rows*sizeof(cuDoubleComplex)));
            cudaCheck(cudaMalloc(&acceptedElectromagnetic_,op.eddy.rows*sizeof(cuDoubleComplex)));
            cudaCheck(cudaMalloc(&heat_,op.heatSamples*sizeof(double)));
            cudaCheck(cudaMalloc(&candidateHeat_,op.heatSamples*sizeof(double)));
            cudaCheck(cudaMalloc(&load_,op.thermalSystem.rows*sizeof(double)));
            cudaCheck(cudaMalloc(&temperature_,op.thermalSystem.rows*sizeof(double)));
            cudaCheck(cudaMalloc(&candidate_,op.thermalSystem.rows*sizeof(double)));
            cudaCheck(cudaMalloc(&controls_,2*sizeof(double)));cudaCheck(cudaMalloc(&invalid_,sizeof(int)));
            cudaCheck(cudaEventCreateWithFlags(&heatReady_,cudaEventDisableTiming));reset();
        }catch(...){release();throw;}
    }
    void abandonAfterDeviceReset() noexcept {
        if(eddy_)eddy_->abandonAfterDeviceReset();
        if(thermal_)thermal_->abandonAfterDeviceReset();
        heat_=candidateHeat_=load_=temperature_=candidate_=controls_=nullptr;
        acceptedElectromagnetic_=candidateElectromagnetic_=nullptr;invalid_=nullptr;heatReady_=nullptr;
    }
    std::size_t ownedBytes()const {
        return eddy_->ownedBytes()+thermal_->ownedBytes()+2*op_.eddy.rows*sizeof(cuDoubleComplex)+
            (2*op_.heatSamples+3*op_.thermalSystem.rows+2)*sizeof(double)+sizeof(int);
    }
    cudaError_t drainCandidateWork()noexcept{const auto first=cudaStreamSynchronize(eddy_->stream());const auto second=cudaStreamSynchronize(thermal_->stream());return first!=cudaSuccess?first:second;}
#ifdef RADIA_IH_TRANSACTION_TEST
    void setTestFailure(unsigned stage){if(stage<1||stage>7)throw std::invalid_argument("Unknown transaction test stage");testFailureStage_=stage;}
    const cuDoubleComplex* testAcceptedEM()const noexcept{return acceptedElectromagnetic_;}
    const double* testAcceptedTemperature()const noexcept{return temperature_;}
    const double* testAcceptedHeat()const noexcept{return heat_;}
#endif
    std::size_t atomicityBytes()const noexcept{return op_.eddy.rows*sizeof(cuDoubleComplex);}
    ~CoupledState(){release();}
    void commitAccepted(const ComplexSolver::Result& em,const ComplexSolver::Result& thermal,bool resetting)noexcept{
        static_assert(noexcept(std::declval<CoupledOperators>().commitMaterial(nullptr)),"Material commit must be noexcept");
        if(op_.commitMaterial)op_.commitMaterial(op_.materialOwner);
        std::swap(acceptedElectromagnetic_,candidateElectromagnetic_);
        std::swap(temperature_,candidate_);std::swap(heat_,candidateHeat_);
        lastEddy_=em;lastThermal_=thermal;if(resetting)steps_=0;else ++steps_;
        if(op_.snapshot)op_.snapshot->commit(resetting);
    }
    static_assert(noexcept(std::declval<CoupledState&>().commitAccepted({}, {}, false)),"Commit must not throw");
    void reset() {
        cudaCheck(cudaMemcpyAsync(candidate_,op_.initialTemperature,op_.thermalSystem.rows*sizeof(double),cudaMemcpyDeviceToDevice,thermal_->stream()));
        cudaCheck(cudaMemsetAsync(candidateHeat_,0,op_.heatSamples*sizeof(double),thermal_->stream()));
        cudaCheck(cudaMemsetAsync(candidateElectromagnetic_,0,op_.eddy.rows*sizeof(cuDoubleComplex),thermal_->stream()));
        if(op_.stageMaterial)op_.stageMaterial(op_.materialOwner,thermal_->stream(),true);
        cudaCheck(cudaStreamSynchronize(thermal_->stream()));cudaCheck(cudaGetLastError());
        testPoint(7);commitAccepted({}, {}, true);
    }
    // Called only at the Simulink discrete Update boundary. Outputs observes
    // accepted state; no thermal vector is returned merely to re-upload it.
    void update(double current,double ambient,double tolerance=1e-10,int maxIterations=2000) {
        if(!std::isfinite(current)||!std::isfinite(ambient))throw std::invalid_argument("Nonfinite IH controls");
        const double controls[2]={current,ambient};
        cudaCheck(cudaMemcpyAsync(controls_,controls,sizeof(controls),cudaMemcpyHostToDevice,eddy_->stream()));
        if(op_.updateMaterial)op_.updateMaterial(op_.materialOwner,eddy_->stream(),temperature_,current);
        scaleDrive<<<unsigned((op_.eddy.rows+127)/128),128,0,eddy_->stream()>>>(op_.rhsPerAmp,controls_,op_.eddy.rows,eddy_->rhs());
        const auto eddyResult=op_.solveMaterial?op_.solveMaterial(op_.materialOwner,*eddy_,tolerance,maxIterations):eddy_->solve(tolerance,maxIterations);
        testPoint(1);
        quadraticLoss<<<unsigned((op_.heatSamples+127)/128),128,0,eddy_->stream()>>>(op_.lossProjection,eddy_->solution(),op_.lossWeights,op_.heatSamples,op_.lossComponents,candidateHeat_);
        cudaCheck(cudaEventRecord(heatReady_,eddy_->stream()));
        cudaCheck(cudaStreamWaitEvent(thermal_->stream(),heatReady_,0));
        const auto blocks=unsigned((op_.thermalSystem.rows+127)/128);
        thermalLoad<<<blocks,128,0,thermal_->stream()>>>(op_.heatProjection,candidateHeat_,load_);
        thermalIncrementRHS<<<blocks,128,0,thermal_->stream()>>>(op_.stiffness,op_.convection,load_,temperature_,op_.constantCoefficients,controls_+1,op_.dt,op_.conductivity,op_.convectionScale,thermal_->rhs());
        const auto thermalResult=thermal_->solve(tolerance,maxIterations);testPoint(2);
        cudaCheck(cudaMemsetAsync(invalid_,0,sizeof(int),thermal_->stream()));
        thermalCandidate<<<blocks,128,0,thermal_->stream()>>>(temperature_,thermal_->solution(),candidate_,invalid_,op_.thermalSystem.rows);
        int invalid=0;cudaCheck(cudaMemcpyAsync(&invalid,invalid_,sizeof(int),cudaMemcpyDeviceToHost,thermal_->stream()));
        cudaCheck(cudaStreamSynchronize(thermal_->stream()));
        if(invalid)throw std::runtime_error("IH candidate temperature rejected; accepted state preserved");
        if(op_.validateCandidateTemperature)op_.validateCandidateTemperature(op_.materialOwner,thermal_->stream(),candidate_);
        cudaCheck(cudaMemcpyAsync(candidateElectromagnetic_,eddy_->solution(),op_.eddy.rows*sizeof(cuDoubleComplex),cudaMemcpyDeviceToDevice,thermal_->stream()));
        testPoint(3);if(op_.stageMaterial)op_.stageMaterial(op_.materialOwner,thermal_->stream(),false);testPoint(4);
        if(op_.snapshot){
            auto& bank=*op_.snapshot;const auto stream=thermal_->stream();bank.begin(stream);
            bank.stage("x",candidateElectromagnetic_,stream);bank.stage("rhs",eddy_->rhs(),stream);
            bank.stage("T_material_eval",temperature_,stream);bank.stage("T_accepted_after",candidate_,stream);
            bank.stage("heat_Qe",candidateHeat_,stream);
            WorkSnapshot::Meta meta{};meta.step=steps_+1;meta.current=current;meta.ambient=ambient;meta.linear=eddyResult.trueRelativeResidual;meta.linearTolerance=tolerance;meta.linearIterations=eddyResult.iterations;meta.linearMaxIterations=unsigned(maxIterations);
            if(op_.stageWorkMaterial)op_.stageWorkMaterial(op_.materialOwner,bank,stream,meta);
            bank.norms(op_.eddy,eddy_->solution(),eddy_->rhs(),stream);bank.end(stream,meta);
        }
        testPoint(5);cudaCheck(cudaStreamSynchronize(thermal_->stream()));cudaCheck(cudaGetLastError());
        if(op_.snapshot)op_.snapshot->finish();testPoint(6);
        commitAccepted(eddyResult,thermalResult,false);
    }
    // Explicit observation only; never an intermediate Simulink transfer.
    std::vector<cuDoubleComplex> electromagneticState() {
        std::vector<cuDoubleComplex> out(op_.eddy.rows);
        cudaCheck(cudaMemcpyAsync(out.data(),acceptedElectromagnetic_,out.size()*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost,eddy_->stream()));
        cudaCheck(cudaStreamSynchronize(eddy_->stream()));return out;
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
    ComplexSolver::Result eddyResult()const{return lastEddy_;}
    ComplexSolver::Result thermalResult()const{return lastThermal_;}
    std::size_t steps()const{return steps_;}
};
} // namespace radia_ih_cuda
