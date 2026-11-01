function tests = test_optuna_cuda
tests = functiontests(localfunctions);
end

function testDensityMatchesUpstream(testCase)
root=fileparts(mfilename('fullpath'));
oracle=jsondecode(fileread(fullfile(root,'..','..','packages','radia-optuna', ...
    'oracles','optuna_cuda_density_oracle.json')));
verifyEqual(testCase,string(oracle.optuna_version),"5.0.0");
previous=radia.optuna.executionMode();
cleanup=onCleanup(@()radia.optuna.executionMode(previous));
for mode=["cpu","gpu"]
radia.optuna.executionMode(mode);
before=optuna_mex('api.info');
for index=1:numel(oracle.cases)
    if iscell(oracle.cases), entry=oracle.cases{index}; else, entry=oracle.cases(index); end
    if string(entry.name)=="categorical"
        actual=optuna_mex('optuna.parzen.log_pdf_categorical',entry.samples, ...
            entry.weights,entry.probabilities);
    else
        step=entry.step; if isempty(step), step=NaN; end
        actual=optuna_mex('optuna.parzen.log_pdf_numerical',entry.samples, ...
            entry.weights,entry.mu,entry.sigma,entry.low,entry.high,entry.log_scale,step);
    end
    verifyEqual(testCase,actual,entry.expected(:),'AbsTol',5e-12,'RelTol',5e-12);
end
info=optuna_mex('api.info');
verifyEqual(testCase,string(info.execution_mode),mode);
if mode=="cpu"
    verifyEqual(testCase,info.density_launches,before.density_launches);
    verifyGreaterThan(testCase,info.cpu_density_calls,before.cpu_density_calls);
else
    verifyGreaterThan(testCase,info.density_launches,before.density_launches);
    verifyEqual(testCase,info.cpu_density_calls,before.cpu_density_calls);
end
end
end

function testCpuParallelDensityAndSwitchValidation(testCase)
previous=radia.optuna.executionMode();
cleanup=onCleanup(@()radia.optuna.executionMode(previous));
radia.optuna.executionMode("cpu");
before=optuna_mex('api.info');
x=linspace(-2,3,512)';
actual=optuna_mex('optuna.parzen.log_pdf_numerical',x,1,0,1,-2,3,false,NaN);
expected=-0.5*x.^2-0.5*log(2*pi)-log(0.5*(erf(3/sqrt(2))-erf(-2/sqrt(2))));
verifyEqual(testCase,actual,expected,'AbsTol',5e-12);
% Exercise the 24-candidate path with enough components to activate OpenMP.
% Identical components have the same analytical density, independently of
% candidate scheduling and the ordered mixture sum.
x=linspace(-2,3,24)';
weights=ones(1024,1)/1024;
beforeDefault=optuna_mex('api.info');
actual=optuna_mex('optuna.parzen.log_pdf_numerical',x,weights, ...
    zeros(1024,1),ones(1024,1),-2,3,false,NaN);
expected=-0.5*x.^2-0.5*log(2*pi)-log(0.5*(erf(3/sqrt(2))-erf(-2/sqrt(2))));
verifyEqual(testCase,actual,expected,'AbsTol',5e-12);
afterDefault=optuna_mex('api.info');
verifyGreaterThanOrEqual(testCase,afterDefault.cpu_max_density_threads,1);
if str2double(getenv('OMP_NUM_THREADS'))>1
    verifyEqual(testCase,afterDefault.cpu_parallel_density_calls, ...
        beforeDefault.cpu_parallel_density_calls+1);
end

verifyError(testCase,@()evaluateCommand('api.execution','auto'),'radia:optuna:mex:Exception');
verifyEqual(testCase,radia.optuna.executionMode(),"cpu");
after=optuna_mex('api.info');
verifyEqual(testCase,after.host_to_device_bytes,before.host_to_device_bytes);
verifyEqual(testCase,after.cuda_device,before.cuda_device);
end

function testGpuModelConstructionMatchesUpstream(testCase)
root=fileparts(mfilename('fullpath'));
oracle=jsondecode(fileread(fullfile(root,'..','..','packages','radia-optuna', ...
    'oracles','optuna50_cuda_oracle.json')));
verifyEqual(testCase,string(oracle.optuna_version),"5.0.0");
previous=radia.optuna.executionMode();
restoreMode=onCleanup(@()radia.optuna.executionMode(previous));
for mode=["cpu","gpu"]
    radia.optuna.executionMode(mode); before=optuna_mex('api.info');
    for entry=reshape(oracle.model_cases,1,[])
        good=entry.good; bad=entry.bad;
        if isempty(good), good=zeros(0,4); end
        if isempty(bad), bad=zeros(0,4); end
        for repeat=1:2
            state=optuna_mex('optuna.random_state.create',entry.seed);
            cleanup=onCleanup(@()optuna_mex('optuna.random_state.destroy',state));
            actual=optuna_mex('optuna.tpe.best_joint_observations',state,24, ...
                [false false false true],entry.lows,entry.highs, ...
                [false true false false],[NaN NaN 1 NaN],[0 0 0 3], ...
                good,bad,1,true,entry.endpoints,[],[]);
            verifyEqual(testCase,actual,reshape(entry.best,1,[]),'AbsTol',5e-11);
            next=optuna_mex('optuna.random_state.rand',state,8);
            verifyEqual(testCase,next(:),entry.next_random(:));
            clear cleanup
        end
    end
    after=optuna_mex('api.info');
    if mode=="gpu"
        verifyGreaterThan(testCase,after.model_build_launches,before.model_build_launches);
        verifyGreaterThan(testCase,after.model_reuses,before.model_reuses);
        verifyGreaterThan(testCase,after.model_download_bytes,before.model_download_bytes);
        verifyEqual(testCase,after.cpu_preparation_calls,before.cpu_preparation_calls);
        verifyGreaterThan(testCase,after.model_radix_sorts,before.model_radix_sorts);
    else
        verifyEqual(testCase,after.model_build_launches,before.model_build_launches);
        verifyEqual(testCase,after.host_to_device_bytes,before.host_to_device_bytes);
        verifyGreaterThan(testCase,after.cpu_cached_normalizer_values,before.cpu_cached_normalizer_values);
    end
end
% Differential guard for the existing stable large-array route, including
% signed zeros. This is separate from the pinned Python oracle cases above.
good=[zeros(129,1);-zeros(128,1)]; bad=good+0.2;
values=zeros(1,2); next=zeros(8,2);
for endpoints=[false true]
    for index=1:2
        modes=["cpu","gpu"];radia.optuna.executionMode(modes(index));
        state=optuna_mex('optuna.random_state.create',731);
        cleanup=onCleanup(@()optuna_mex('optuna.random_state.destroy',state));
        values(index)=optuna_mex('optuna.tpe.best_joint_observations',state,24, ...
            false,-2,2,false,NaN,0,good,bad,1,true,endpoints,[],[]);
        next(:,index)=optuna_mex('optuna.random_state.rand',state,8);
        clear cleanup
    end
    verifyEqual(testCase,values(1),values(2),'AbsTol',5e-12);
    verifyEqual(testCase,next(:,1),next(:,2));
end
end

function testCudaBuildBoundary(testCase)
info=optuna_mex('api.info');
verifyEqual(testCase,string(info.build_route),"nvcc-shared-library");
verifyFalse(testCase,info.taskmanager);
verifyNotEmpty(testCase,string(info.gpu_operations));
end

function testValidationAndRepeatedDeviceUse(testCase)
previous=radia.optuna.executionMode();
restoreMode=onCleanup(@()radia.optuna.executionMode(previous));
call=@()evaluateCommand('optuna.parzen.log_pdf_numerical',[0;1],[.3;.7], ...
    [0;1],[1;1],-2,3,false,NaN);
expected=call();
before=optuna_mex('api.info');
for index=1:20, verifyEqual(testCase,call(),expected,'AbsTol',0); end
after=optuna_mex('api.info');
if radia.optuna.executionMode()=="gpu"
    verifyEqual(testCase,after.device_allocation_calls,before.device_allocation_calls);
    verifyEqual(testCase,after.host_to_device_bytes,before.host_to_device_bytes);
    verifyGreaterThan(testCase,after.workspace_reuses,before.workspace_reuses);
    verifyGreaterThan(testCase,after.avoided_upload_bytes,before.avoided_upload_bytes);
end
verifyError(testCase,@()evaluateCommand('optuna.parzen.log_pdf_categorical', ...
    0,1,[.5,.5]),'radia:optuna:mex:Exception');
verifyError(testCase,@()evaluateCommand('optuna.parzen.log_pdf_numerical', ...
    1e308,1,0,1,-2,3,false,NaN),'radia:optuna:mex:Exception');
verifyEqual(testCase,call(),expected,'AbsTol',0);
radia.optuna.executionMode("cpu");
after=optuna_mex('api.info');
verifyEqual(testCase,after.cached_device_bytes,0);
radia.optuna.executionMode("gpu");
verifyEqual(testCase,call(),expected,'AbsTol',0);
end

function testExecutionFailureRestoresRandomState(testCase)
previous=radia.optuna.executionMode();
restoreMode=onCleanup(@()radia.optuna.executionMode(previous));
for mode=["cpu","gpu"]
radia.optuna.executionMode(mode);
lower=struct('categorical',false,'weights',1,'mu',0,'sigma',1, ...
    'low',-2,'high',3,'internal_low',-2,'internal_high',3,'log',false,'step',NaN);
upper=lower;
upper.sigma=1e-308; % Finite valid input whose evaluated tails underflow.
state=optuna_mex('optuna.random_state.create',42);
reference=optuna_mex('optuna.random_state.create',42);
cleanup=onCleanup(@()destroyStates(state,reference));
verifyError(testCase,@()evaluateCommand('optuna.tpe.best_joint',state,24,false,{lower},{upper}), ...
    'radia:optuna:mex:Exception');
actual=optuna_mex('optuna.random_state.rand',state,8);
expected=optuna_mex('optuna.random_state.rand',reference,8);
verifyEqual(testCase,actual,expected);
clear cleanup
state=optuna_mex('optuna.random_state.create',73);
reference=optuna_mex('optuna.random_state.create',73);
cleanup=onCleanup(@()destroyStates(state,reference));
for repeat=1:2
    before=optuna_mex('api.info');
    % Overflowing support produces an invalid model after GPU construction.
    verifyError(testCase,@()evaluateCommand('optuna.tpe.best_joint_observations', ...
        state,24,false,-1e308,1e308,false,NaN,0,0,1,1,true,false,[],[]), ...
        'radia:optuna:mex:Exception');
    actual=optuna_mex('optuna.random_state.rand',state,8);
    expected=optuna_mex('optuna.random_state.rand',reference,8);
    verifyEqual(testCase,actual,expected);
    after=optuna_mex('api.info');
    if mode=="gpu"
        verifyGreaterThan(testCase,after.model_build_launches,before.model_build_launches);
    end
end
clear cleanup
end
end

function destroyStates(state,reference)
optuna_mex('optuna.random_state.destroy',state);
optuna_mex('optuna.random_state.destroy',reference);
end

function value=evaluateCommand(varargin)
% Request the actual command output so errors cannot come from MEX nargout.
value=optuna_mex(varargin{:});
end
