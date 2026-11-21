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
    options.LinearTolerance (1,1) double {mustBeFinite,mustBePositive,mustBeLessThanOrEqual(options.LinearTolerance,1e-6)} = 1e-10
    options.NonlinearTolerance (1,1) double {mustBeFinite,mustBePositive,mustBeLessThanOrEqual(options.NonlinearTolerance,1e-9)} = 1e-9
    options.ReturnOnGPU (1,1) logical = true
end
device=gpuDevice();wait(device);total=tic;
[fem,assembly]=radia.gpu.assembleHCurl(volPath,boundaryLabels, ...
    Dirichlet=options.Dirichlet,Reluctivity=options.Reluctivity,MassCoefficient=options.MassCoefficient, ...
    CurrentDensity=options.CurrentDensity,Gamma=options.Gamma);
nonlinear=radia_cuda_ams('nonlinear_create',fem);
nonlinearCleanup=onCleanup(@()radia_cuda_ams('nonlinear_release',nonlinear)); %#ok<NASGU>
rhs=fem.rhs;x=zeros(size(rhs),'like',rhs);history=zeros(0,5);
normRhs=gather(norm(rhs));builds=0;refreshes=0;refreshSeconds=0;
setupSeconds=0;solveSeconds=0;assemblySeconds=0;peakOwned=assembly.gpu_peak_owned_bytes_upper_bound;
relative=inf;lastSolve=[];
for outer=0:30
    wait(device);timer=tic;[residual,~,oldEnergy]=radia_cuda_ams('nonlinear_evaluate',nonlinear,x);
    wait(device);assemblySeconds=assemblySeconds+toc(timer);
    if normRhs==0;relative=gather(norm(residual));else;relative=gather(norm(residual))/normRhs;end
    assert(isfinite(relative),'radia:cuda:nonlinear','Nonfinite nonlinear true residual');
    if relative<=options.NonlinearTolerance;break;end
    assert(outer<30,'radia:cuda:nonlinear','Newton did not converge');
    wait(device);timer=tic;
    [delta,lastSolve]=radia_cuda_ams('nonlinear_solve_resident',nonlinear,-residual,options.LinearTolerance,4000);
    wait(device);phaseSeconds=toc(timer);
    builds=builds+lastSolve.hierarchy_built;refreshes=refreshes+lastSolve.hierarchy_refreshed;
    setupSeconds=setupSeconds+lastSolve.hierarchy_built*lastSolve.hierarchy_seconds;
    refreshSeconds=refreshSeconds+lastSolve.hierarchy_refreshed*lastSolve.hierarchy_seconds;
    solveSeconds=solveSeconds+max(0,phaseSeconds-lastSolve.hierarchy_seconds);
    peakOwned=max(peakOwned,lastSolve.gpu_peak_owned_bytes);
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
    'history',history,'hierarchy_builds',builds,'hierarchy_refreshes',refreshes,'refresh_seconds',refreshSeconds,'reuse_policy',"refresh",'assembly',assembly, ...
    'setup_seconds',setupSeconds,'assembly_seconds',assemblySeconds,'solve_seconds',solveSeconds, ...
    'total_seconds',toc(total),'max_phase_owned_bytes_upper_bound',peakOwned,'last_linear_solve',lastSolve, ...
    'memory_definition','max of assembly bound and AMS setup owned peak; not aggregate whole-process peak', ...
    'backend','MATLAB/mexcuda GPU; NGSolve host basis; no Python adapter');
end
