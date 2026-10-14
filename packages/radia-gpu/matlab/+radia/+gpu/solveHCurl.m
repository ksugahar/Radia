function [solution,info] = solveHCurl(volPath,boundaryLabels,options)
%SOLVEHCURL Native .vol -> sparse magnetic FE -> CUDA AMS/Newton solve.
% H(B)=(Reluctivity+Gamma*|B|^2)*B with positive gauge mass. SI.
% CPU references are independent validation, never a runtime prerequisite.
arguments
    volPath (1,1) string
    boundaryLabels (1,:) string
    options.Dirichlet (1,1) string = ".*"
    options.Reluctivity (1,1) double {mustBeFinite,mustBePositive} = 1
    options.MassCoefficient (1,1) double {mustBeFinite,mustBePositive} = 1
    options.CurrentDensity (1,3) double {mustBeFinite} = [1 0 0]
    options.Gamma (1,1) double {mustBeFinite,mustBeNonnegative} = 0
    options.AMSReuse (1,1) string {mustBeMember(options.AMSReuse,["adaptive","frozen","rebuild"])} = "adaptive"
    options.SetupThreads (1,1) double {mustBeInteger,mustBeInRange(options.SetupThreads,1,64)} = 1
    options.LinearTolerance (1,1) double {mustBeFinite,mustBePositive,mustBeLessThanOrEqual(options.LinearTolerance,1e-6)} = 1e-10
    options.NonlinearTolerance (1,1) double {mustBeFinite,mustBePositive,mustBeLessThanOrEqual(options.NonlinearTolerance,1e-9)} = 1e-9
    options.ReturnOnGPU (1,1) logical = true
end
device=gpuDevice();wait(device);total=tic;
[model,fem,assembly]=radia.gpu.assembleHCurl(volPath,boundaryLabels, ...
    Dirichlet=options.Dirichlet,Reluctivity=options.Reluctivity,MassCoefficient=options.MassCoefficient, ...
    CurrentDensity=options.CurrentDensity,Gamma=options.Gamma);
nonlinear=radia_cuda_ams('nonlinear_create',fem);
nonlinearCleanup=onCleanup(@()radia_cuda_ams('nonlinear_release',nonlinear)); %#ok<NASGU>
context=[];contextCleanup=[]; %#ok<NASGU>
rhs=fem.rhs;x=zeros(size(rhs),'like',rhs);history=zeros(0,5);
normRhs=gather(norm(rhs));builds=0;firstIterations=0;rebuildNext=false;
setupSeconds=0;solveSeconds=0;assemblySeconds=0;peakOwned=assembly.gpu_peak_owned_bytes_upper_bound;
op=fem.mass;relative=inf;lastSolve=[];
for outer=0:30
    wait(device);timer=tic;[residual,op.values,oldEnergy]=radia_cuda_ams('nonlinear_evaluate',nonlinear,x);
    wait(device);assemblySeconds=assemblySeconds+toc(timer);
    if normRhs==0;relative=gather(norm(residual));else;relative=gather(norm(residual))/normRhs;end
    assert(isfinite(relative),'radia:cuda:nonlinear','Nonfinite nonlinear true residual');
    if relative<=options.NonlinearTolerance;break;end
    assert(outer<30,'radia:cuda:nonlinear','Newton did not converge');
    if isempty(context)||options.AMSReuse=="rebuild"||(options.AMSReuse=="adaptive"&&rebuildNext)
        wait(device);timer=tic;model.A.values=radia_cuda_ams('nonlinear_raw',nonlinear);
        [hierarchy,setup]=radia.gpu.prepareAMS(model,SetupThreads=options.SetupThreads);
        peakOwned=max(peakOwned,setup.gpu_setup_peak_owned_bytes);
        if isempty(context)
            context=radia.gpu.AMSContext(hierarchy);contextCleanup=onCleanup(@()delete(context));
        else;context.update(hierarchy);end
        clear hierarchy;wait(device);setupSeconds=setupSeconds+toc(timer);builds=builds+1;firstIterations=0;
    end
    wait(device);timer=tic;[delta,lastSolve]=context.solve(op,-residual,options.LinearTolerance,4000,true);
    wait(device);solveSeconds=solveSeconds+toc(timer);
    if firstIterations==0;firstIterations=lastSolve.iterations;end
    rebuildNext=lastSolve.iterations>max(30,2*firstIterations);
    slope=gather(residual'*delta);assert(isfinite(slope)&&slope<0,'radia:cuda:descent','Newton direction is not descent');alpha=1;
    while radia_cuda_ams('nonlinear_energy',nonlinear,x+alpha*delta)>oldEnergy+1e-4*alpha*slope
        if abs(alpha*slope)<=64*eps*max(1,abs(oldEnergy))
            [trialResidual,~,~]=radia_cuda_ams('nonlinear_evaluate',nonlinear,x+alpha*delta);
            trialRelative=gather(norm(trialResidual))/normRhs;
            if trialRelative<.9*relative||trialRelative<=options.NonlinearTolerance;break;end
        end
        alpha=alpha/2;assert(alpha>=2^-20,'radia:cuda:lineSearch','Newton line search failed');
    end
    x=x+alpha*delta;history(end+1,:)=[outer relative lastSolve.iterations lastSolve.true_relative_residual alpha]; %#ok<AGROW>
end
if options.ReturnOnGPU;solution=x;else;solution=gather(x);end
wait(device);
info=struct('schema','radia.native-hcurl-solve.v1','true_nonlinear_relative_residual',relative, ...
    'history',history,'hierarchy_builds',builds,'reuse_policy',options.AMSReuse,'assembly',assembly, ...
    'setup_seconds',setupSeconds,'assembly_seconds',assemblySeconds,'solve_seconds',solveSeconds, ...
    'total_seconds',toc(total),'max_phase_owned_bytes_upper_bound',peakOwned,'last_linear_solve',lastSolve, ...
    'memory_definition','max of assembly bound and AMS setup owned peak; not aggregate whole-process peak', ...
    'backend','MATLAB/mexcuda GPU; NGSolve host basis and serial sparsecholesky setup; no Python adapter');
end
