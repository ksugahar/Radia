function tests = test_optuna_cuda
tests = functiontests(localfunctions);
end

function testDensityMatchesUpstream(testCase)
root=fileparts(mfilename('fullpath'));
oracle=jsondecode(fileread(fullfile(root,'fixtures','optuna_cuda_density_oracle.json')));
verifyEqual(testCase,string(oracle.optuna_version),"5.0.0");
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
