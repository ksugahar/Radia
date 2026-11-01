function result = validate_optuna_cuda(repoRoot, outputPath, executionMode)
%VALIDATE_OPTUNA_CUDA Explicit CPU/GPU execution and pinned upstream TPE acceptance.
arguments
    repoRoot (1,1) string
    outputPath (1,1) string = fullfile(repoRoot,'validation_test','optimization','results_optuna_cuda_100.json')
    executionMode (1,1) string {mustBeMember(executionMode,["cpu","gpu"])} = "gpu"
end
addpath(fullfile(repoRoot,'packages','radia-optuna','matlab'));
binary = string(java.io.File(which('optuna_mex')).getCanonicalPath());
expected = string(java.io.File(fullfile(repoRoot,'packages','radia-optuna','matlab',"optuna_mex."+mexext)).getCanonicalPath());
assert(strcmpi(binary,expected), 'Unexpected native gateway: %s',binary);
previousMode = radia.optuna.executionMode();
restoreMode = onCleanup(@()radia.optuna.executionMode(previousMode));
radia.optuna.executionMode(executionMode);
assert(radia.optuna.executionMode()==executionMode);
try
    rejected=optuna_mex('api.execution','auto');
    error('Unexpected execution mode accepted: %s',rejected);
catch exception
    assert(contains(string(exception.message),'execution mode must be cpu or gpu'));
end
assert(radia.optuna.executionMode()==executionMode);
if executionMode=="gpu"
    device = gpuDevice;
    gpuName=string(device.Name); capability=string(device.ComputeCapability); driver=string(device.DriverVersion);
else
    gpuName="not initialized"; capability=""; driver="";
end
fixturePath = fullfile(repoRoot,'packages','radia-optuna','oracles','optuna50_cuda_oracle.json');
oracle = jsondecode(fileread(fixturePath));
assert(string(oracle.optuna_version) == "5.0.0");
result = struct(schema="radia.optuna.cuda-acceptance.v1", ok=false, ...
    platform="windows-x64", source_sha=string(getenv('RADIA_OPTUNA_SOURCE_SHA')), ...
    execution_route="official-matlab-engine-or-batch", matlab=string(version), ...
    release=string(version('-release')), ...
    execution_mode=executionMode, gpu=gpuName, compute_capability=capability, ...
    driver=driver, binary="packages/radia-optuna/matlab/optuna_mex.mexw64", binary_sha256=sha256(binary), ...
    fixture_sha256=sha256(fixturePath), optuna_version=string(oracle.optuna_version), ...
    numpy_version=string(oracle.numpy_version), timing_class="diagnostic-only-not-performance-acceptance");
sources=["packages/radia-optuna/native/optuna_mex.cpp", ...
    "packages/radia-optuna/native/optuna_cuda.cu", "packages/radia-optuna/native/optuna_cuda.hpp", ...
    "packages/radia-optuna/tools/build_native.py", ...
    "packages/radia-optuna/matlab/+radia/+optuna/+internal/NativeKernels.m", ...
    "validation_test/optimization/validate_optuna_cuda.m"];
result.sources=struct(path=sources,sha256=arrayfun(@(p)sha256(fullfile(repoRoot,p)),sources));
before = optuna_mex('api.info');
assert(string(before.execution) == "matlab-mex-shared-cuda", 'A CUDA gateway is required.');
records = cell(0,1);
try
    for index = 1:numel(oracle.cases)
        entry = oracle.cases(index);
        lower = descriptors(entry.below);
        upper = descriptors(entry.above);
        flags = cellfun(@(e)e.categorical,lower);
        handle = optuna_mex('optuna.random_state.create',double(entry.seed));
        cleanup = onCleanup(@()optuna_mex('optuna.random_state.destroy',handle));
        actual = optuna_mex('optuna.tpe.best_joint',handle,double(entry.count),flags,lower,upper);
        errorBest = max(abs(actual(:)-entry.best(:)));
        assert(errorBest <= 5e-11, 'Upstream candidate differs in %s: %.17g',entry.name,errorBest);
        next = optuna_mex('optuna.random_state.rand',handle,8);
        assert(isequal(next(:),entry.next_random(:)), 'Random consumption differs in %s',entry.name);
        clear cleanup
        densityError = 0;
        if numel(lower) == 1
            samples = entry.samples(:);
            densityError = max(abs(densityCall(samples,lower{1})-entry.below_pdf(:)));
            densityError = max(densityError,max(abs(densityCall(samples,upper{1})-entry.above_pdf(:))));
            assert(densityError <= 1e-10,'Upstream density differs in %s: %.17g',entry.name,densityError);
            if ~flags
                handle = optuna_mex('optuna.random_state.create',double(entry.seed));
                cleanup = onCleanup(@()optuna_mex('optuna.random_state.destroy',handle));
                l=lower{1}; u=upper{1};
                value=optuna_mex('optuna.tpe.best_numerical',handle,double(entry.count), ...
                    l.weights,l.mu,l.sigma,u.weights,u.mu,u.sigma, ...
                    l.internal_low,l.internal_high,l.low,l.high,l.log,l.step);
                assert(abs(value-entry.best) <= 5e-11,'Fused numerical candidate differs in %s',entry.name);
                next=optuna_mex('optuna.random_state.rand',handle,8);
                assert(isequal(next(:),entry.next_random(:)));
                clear cleanup
            end
        end
        % Identical densities exercise first-candidate tie breaking on device.
        handle = optuna_mex('optuna.random_state.create',double(entry.seed));
        cleanup = onCleanup(@()optuna_mex('optuna.random_state.destroy',handle));
        tied=optuna_mex('optuna.tpe.best_joint',handle,double(entry.count),flags,lower,lower);
        first=entry.samples(1,:);
        assert(max(abs(tied(:)-first(:))) <= 5e-11,'Tie breaking differs');
        clear cleanup
        records{index}=struct(name=string(entry.name),density_max_abs_error=densityError, ...
            selected_candidate_max_abs_error=errorBest,random_state_exact=true,tie_breaking=true);
    end
    % History/grouping parity is covered once by validate_cuda.
    if executionMode=="gpu", result.density_diagnostics=profileDensity(); end
    result.ok=true;
catch exception
    result.error=string(getReport(exception,'extended','hyperlinks','off'));
end

result.cases=records;
result.native=optuna_mex('api.info');
result.density_launch_delta=result.native.density_launches-before.density_launches;
result.selection_launch_delta=result.native.selection_launches-before.selection_launches;
result.cpu_density_delta=result.native.cpu_density_calls-before.cpu_density_calls;
result.cpu_selection_delta=result.native.cpu_selection_calls-before.cpu_selection_calls;
if executionMode=="gpu"
    result.ok=result.ok && result.density_launch_delta>0 && result.selection_launch_delta>0 ...
        && result.cpu_density_delta==0 && result.cpu_selection_delta==0;
else
    result.ok=result.ok && result.density_launch_delta==0 && result.selection_launch_delta==0 ...
        && result.cpu_density_delta>0 && result.cpu_selection_delta>0 ...
        && result.native.host_to_device_bytes==before.host_to_device_bytes ...
        && result.native.device_to_host_bytes==before.device_to_host_bytes ...
        && result.native.cuda_device==before.cuda_device;
end
folder=fileparts(outputPath);
if ~isfolder(folder), mkdir(folder); end
fid=fopen(outputPath,'w');
assert(fid>=0,'Cannot write acceptance JSON');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
clear cleanup
assert(result.ok,'Optuna CUDA acceptance failed; see %s',outputPath);
end

function values=descriptors(raw)
if iscell(raw), values=raw; else, values=arrayfun(@(e)e,raw,UniformOutput=false); end
for index=1:numel(values)
    entry=values{index};
    entry.weights=entry.weights(:);
    if ~entry.categorical
        entry.mu=entry.mu(:); entry.sigma=entry.sigma(:);
        if isempty(entry.step), entry.step=NaN; end
    end
    values{index}=entry;
end
end

function result=densityCall(samples,e)
if e.categorical
    result=optuna_mex('optuna.parzen.log_pdf_categorical',samples,e.weights,e.probabilities);
else
    result=optuna_mex('optuna.parzen.log_pdf_numerical',samples,e.weights,e.mu,e.sigma, ...
        e.internal_low,e.internal_high,e.log,e.step);
end
end

function value=sha256(filename)
digest=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(filename,'rb'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
while ~feof(fid), digest.update(fread(fid,1048576,'*uint8')); end
value=lower(string(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[])));
end

function records=profileDensity()
% Diagnostic provenance on 100; formal timing acceptance remains separate.
shapes=[24 32; 4096 512];
records=cell(size(shapes,1),1);
for index=1:size(shapes,1)
    count=shapes(index,1); kernels=shapes(index,2);
    x=linspace(-3,3,count).'; mu=linspace(-2,2,kernels);
    sigma=linspace(0.4,1.3,kernels); weights=ones(1,kernels)/kernels;
    call=@()optuna_mex('optuna.parzen.log_pdf_numerical',x,weights,mu,sigma,-3,3,false,NaN);
    reference=@()continuousReference(x,weights,mu,sigma);
    start=tic; actual=call(); cold=toc(start);
    expected=reference();
    difference=max(abs(actual-expected));
    assert(difference<=1e-10,'Large density/reference mismatch');
    cudaTimes=zeros(1,5); hostTimes=zeros(1,5);
    for repetition=1:5
        start=tic; value=call(); cudaTimes(repetition)=toc(start);
        start=tic; hostValue=reference(); hostTimes(repetition)=toc(start);
    end
    records{index}=struct(candidates=count,kernels=kernels, ...
        first_call_seconds=cold,cuda_transfer_inclusive_median_seconds=median(cudaTimes), ...
        matlab_host_reference_median_seconds=median(hostTimes), ...
        host_over_cuda_ratio=median(hostTimes)/median(cudaTimes), ...
        max_abs_error=difference,classification="100-diagnostic-only");
end
end

function result=continuousReference(x,weights,mu,sigma)
mass=0.5*(erf((3-mu)./sigma/sqrt(2))-erf((-3-mu)./sigma/sqrt(2)));
terms=-0.5*((x-mu)./sigma).^2-log(sigma)-0.5*log(2*pi)-log(mass)+log(weights);
maximum=max(terms,[],2);
result=maximum+log(sum(exp(terms-maximum),2));
end
