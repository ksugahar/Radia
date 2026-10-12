function tests = test_optuna_cuda
tests = functiontests(localfunctions);
end

function testDensityMatchesUpstream(testCase)
root=fileparts(mfilename('fullpath'));
oracle=jsondecode(fileread(fullfile(root,'fixtures','optuna_cuda_density_oracle.json')));
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
verifyError(testCase,@()optuna_mex('api.execution','auto'),'radia:optuna:mex:Exception');
verifyEqual(testCase,radia.optuna.executionMode(),"cpu");
after=optuna_mex('api.info');
verifyEqual(testCase,after.host_to_device_bytes,before.host_to_device_bytes);
verifyEqual(testCase,after.cuda_device,before.cuda_device);
end

function testCudaBuildBoundary(testCase)
info=optuna_mex('api.info');
verifyEqual(testCase,string(info.build_route),"mexcuda");
verifyFalse(testCase,info.taskmanager);
verifyNotEmpty(testCase,string(info.gpu_operations));
end

function testValidationAndRepeatedDeviceUse(testCase)
call=@()optuna_mex('optuna.parzen.log_pdf_numerical',[0;1],[.3;.7], ...
    [0;1],[1;1],-2,3,false,NaN);
expected=call();
for index=1:20, verifyEqual(testCase,call(),expected,'AbsTol',0); end
verifyError(testCase,@()optuna_mex('optuna.parzen.log_pdf_categorical', ...
    0,1,[.5,.5]),'radia:optuna:mex:Exception');
verifyError(testCase,@()optuna_mex('optuna.parzen.log_pdf_numerical', ...
    1e308,1,0,1,-2,3,false,NaN),'radia:optuna:mex:Exception');
verifyEqual(testCase,call(),expected,'AbsTol',0);
end

function testExecutionFailureRestoresRandomState(testCase)
previous=radia.optuna.executionMode();
restoreMode=onCleanup(@()radia.optuna.executionMode(previous));
for mode=["cpu","gpu"]
radia.optuna.executionMode(mode);
lower=struct('categorical',true,'weights',1,'probabilities',[1,0]);
upper=struct('categorical',true,'weights',1,'probabilities',[0,1]);
state=optuna_mex('optuna.random_state.create',42);
reference=optuna_mex('optuna.random_state.create',42);
cleanup=onCleanup(@()destroyStates(state,reference));
verifyError(testCase,@()optuna_mex('optuna.tpe.best_joint',state,24,true,{lower},{upper}), ...
    'radia:optuna:mex:Exception');
actual=optuna_mex('optuna.random_state.rand',state,8);
expected=optuna_mex('optuna.random_state.rand',reference,8);
verifyEqual(testCase,actual,expected);
clear cleanup
end
end

function destroyStates(state,reference)
optuna_mex('optuna.random_state.destroy',state);
optuna_mex('optuna.random_state.destroy',reference);
end
