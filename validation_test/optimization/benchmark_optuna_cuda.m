function result=benchmark_optuna_cuda(binaryDirectory,outputPath,variant,rounds)
%BENCHMARK_OPTUNA_CUDA Transfer-inclusive density and complete TPE calls.
% Run each variant in a fresh, owned MATLAB session on the same idle GPU host.
arguments
    binaryDirectory (1,1) string
    outputPath (1,1) string
    variant (1,1) string
    rounds (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(rounds,10)} = 31
end
addpath(binaryDirectory);
binary=string(which('optuna_mex'));
assert(binary==fullfile(binaryDirectory,"optuna_mex."+mexext));
previousThreads=maxNumCompThreads(1);
threadCleanup=onCleanup(@()maxNumCompThreads(previousThreads));
result=struct(schema="radia.optuna.cuda-performance.v1",ok=false,variant=variant, ...
    timing_class="formal-user-selected-idle-gpu-host",platform="windows-x64", ...
    matlab=string(version),matlab_threads=1,openmp_threads=string(getenv('OMP_NUM_THREADS')), ...
    rounds=rounds,binary_sha256=sha256(binary), ...
    benchmark_sha256=sha256(string(mfilename('fullpath'))+".m"), ...
    source_sha=string(getenv('RADIA_OPTUNA_SOURCE_SHA')));
try
    start=tic; device=gpuDevice; result.gpu_initialization_seconds=toc(start);
    result.gpu=string(device.Name); result.compute_capability=string(device.ComputeCapability);
    result.gpu_total_bytes=device.TotalMemory; result.gpu_available_before_bytes=device.AvailableMemory;
    shapes=[24 32 1;256 128 1;4096 512 1;24 32 5;256 128 5;4096 512 5];
    records=cell(size(shapes,1),1);
    for index=1:size(shapes,1)
        count=shapes(index,1); kernels=shapes(index,2); dimensions=shapes(index,3);
        x=repmat(linspace(-3,3,count).',1,dimensions);
        mu=linspace(-2,2,kernels); sigma=linspace(0.4,1.3,kernels);
        weights=ones(1,kernels)/kernels;
        e=struct(weights=weights,categorical=false,mu=mu,sigma=sigma, ...
            internal_low=-3,internal_high=3,low=-3,high=3,log=false,step=NaN);
        above=e; above.mu=mu+0.15;
        lower=repmat({e},1,dimensions); upper=repmat({above},1,dimensions);
        % Public density commands are univariate; dimension scaling is
        % measured through the joint candidate-evaluation command below.
        density=@()optuna_mex('optuna.parzen.log_pdf_numerical',x(:,1),weights,mu,sigma,-3,3,false,NaN);
        reference=@()hostDensity(x(:,1),weights,mu,sigma);
        start=tic; actual=density(); first=toc(start);
        difference=max(abs(actual-reference())); assert(difference<=1e-10);
        batch=max(1,floor(2400/count));
        gpuTimes=measure(density,rounds,batch); hostTimes=measure(reference,rounds,batch);
        expected=bestCall(count,dimensions,lower,upper);
        bestTimes=zeros(1,rounds);
        for repetition=1:rounds
            % RNG creation/destruction is outside the timed call. Each call
            % starts from the same seed; sampling and selection are timed.
            handle=optuna_mex('optuna.random_state.create',201);
            cleanup=onCleanup(@()optuna_mex('optuna.random_state.destroy',handle));
            start=tic;
            selected=optuna_mex('optuna.tpe.best_joint',handle,count,false(1,dimensions),lower,upper);
            bestTimes(repetition)=toc(start);
            assert(isequal(selected,expected),'Repeated candidate selection changed');
            clear cleanup
        end
        records{index}=struct(candidates=count,kernels=kernels,dimensions=dimensions, ...
            density_dimensions=1,density_first_call_seconds=first,density_batch=batch, ...
            density_seconds=gpuTimes,matlab_reference_seconds=hostTimes, ...
            density_median_seconds=median(gpuTimes),matlab_reference_median_seconds=median(hostTimes), ...
            host_over_cuda_ratio=median(hostTimes)/median(gpuTimes), ...
            tpe_seconds=bestTimes,tpe_median_seconds=median(bestTimes), ...
            selected=expected,max_abs_density_error=difference);
    end
    result.cases=vertcat(records{:}); result.native=optuna_mex('api.info');
    result.gpu_available_after_bytes=device.AvailableMemory;
    assert(result.native.density_launches>0 && result.native.selection_launches>0);
    result.ok=true;
catch exception
    result.error_identifier=string(exception.identifier);
    result.error_message=string(exception.message);
end
fid=fopen(outputPath,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true)); clear cleanup
assert(result.ok,'CUDA performance evaluation failed; see %s',outputPath);
end

function values=measure(call,rounds,batch)
value=call(); values=zeros(1,rounds); %#ok<NASGU>
for repetition=1:rounds
    start=tic;
    for iteration=1:batch, value=call(); end %#ok<NASGU>
    values(repetition)=toc(start)/batch;
end
end

function selected=bestCall(count,dimensions,lower,upper)
handle=optuna_mex('optuna.random_state.create',201);
cleanup=onCleanup(@()optuna_mex('optuna.random_state.destroy',handle));
selected=optuna_mex('optuna.tpe.best_joint',handle,count,false(1,dimensions),lower,upper);
end

function result=hostDensity(x,weights,mu,sigma)
mass=0.5*(erf((3-mu)./sigma/sqrt(2))-erf((-3-mu)./sigma/sqrt(2)));
terms=-0.5*((x-mu)./sigma).^2-log(sigma)-0.5*log(2*pi)-log(mass)+log(weights);
maximum=max(terms,[],2); result=maximum+log(sum(exp(terms-maximum),2));
end

function value=sha256(filename)
digest=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(filename,'rb'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
while ~feof(fid), digest.update(fread(fid,1048576,'*uint8')); end
value=lower(string(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[])));
end
