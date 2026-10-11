function result = check_installed(payload, output, checkLTspice)
%CHECK_INSTALLED Exercise only the wheel's MATLAB payload in an owned session.
arguments
    payload (1,1) string
    output (1,1) string
    checkLTspice (1,1) logical = true
end
assert(isfolder(payload));
if ~isfolder(output), mkdir(output); end
addpath(payload);
resolved = string(which('radia.gpu.csrPCG'));
assert(startsWith(resolved,payload+filesep),'Checkout shadowed installed GPU code.');
for entryPoint = ["radia.simulink.buildCudaPCGBlock","radia_gpu_pcg_sfun","radia.ltspice.run"]
    assert(startsWith(string(which(entryPoint)),payload+filesep),'Installed adapter resolved outside the wheel.');
end
radia.gpu.requireRuntime();
system = struct('row_pointers',int32([0;2;5;7]), ...
    'columns',int32([0;1;0;1;2;1;2]), ...
    'values',[2;-1;-1;2;-1;-1;2],'inverse_diagonal',[.5;.5;.5]);
[x,info] = radia.gpu.csrPCG(gpuArray(system.row_pointers), ...
    gpuArray(system.columns),gpuArray(system.values),gpuArray([1;0;1]), ...
    gpuArray(system.inverse_diagonal),1e-12,100);
assert(norm(gather(x)-ones(3,1))<1e-12 && info.true_relative_residual<1e-12);
model = 'radia_gpu_installed_acceptance';
assert(~bdIsLoaded(model),'Acceptance model already belongs to another operation.');
cleanup = onCleanup(@() closeOwned(model));
blockOutput = fullfile(output,'simulink');
new_system(model);
workspace = get_param(model,'ModelWorkspace');
workspace.assignin('cudaSystem',system);
radia.simulink.buildCudaPCGBlock(model,OutputDirectory=blockOutput);
add_block('simulink/Sources/Constant',model+"/RHS",'Value','[1;0;1]');
add_block('simulink/Sinks/To Workspace',model+"/Solution", ...
    'VariableName','gpuSolution','SaveFormat','Array');
add_line(model,'RHS/1','CUDA PCG/1');
add_line(model,'CUDA PCG/1','Solution/1');
set_param(model,'SolverType','Fixed-step','Solver','FixedStepDiscrete', ...
    'FixedStep','0.001','StopTime','0.002','ReturnWorkspaceOutputs','on');
% Native edit/check/save/close/reopen on the exact path; never edit SLX ZIP data.
set_param(model,'SimulationCommand','update');
modelFile = fullfile(output,model+".slx");
save_system(model,modelFile); close_system(model,0); load_system(modelFile);
for index = 1:2
    simulation = sim(model);
    values = simulation.get('gpuSolution');
    assert(max(abs(values(:)-1))<1e-12);
    assert(radia.gpu.internal.systemStore('count','')==0,'GPU block leaked state.');
    evidence = jsondecode(fileread(fullfile(blockOutput,'result.json')));
    assert(evidence.passed && evidence.calls>=3 && evidence.max_true_residual<1e-12);
end
% A failed solve must leave failure evidence and release uploaded block state.
set_param(model+"/CUDA PCG",'Parameters', ...
    sprintf('cudaSystem,1e-12,1,0.001,''%s''',replace(blockOutput,"'","''")));
failed = false;
try
    sim(model);
catch
    failed = true;
end
assert(failed,'Nonconvergence unexpectedly succeeded.');
assert(radia.gpu.internal.systemStore('count','')==0,'Failed block leaked state.');
evidence = jsondecode(fileread(fullfile(blockOutput,'result.json')));
assert(~evidence.passed);
if checkLTspice
netlist = fullfile(output,'divider.cir');
file = fopen(netlist,'w'); assert(file>=0);
fprintf(file,'Radia independent resistor divider\nV1 in 0 2\nR1 in out 1000\nR2 out 0 1000\n.tran 0 0.001\n.end\n');
fclose(file);
circuit = radia.ltspice.run(netlist,OutputDirectory=fullfile(output,'ltspice'),Timeout_s=60);
assert(isfile(circuit.raw_file));
waveform = circuit.waveform;
voltage = waveform.values(:,waveform.names=="V(out)");
if isempty(voltage), voltage = waveform.values(:,strcmpi(waveform.names,"V(out)")); end
assert(~isempty(voltage) && max(abs(voltage(:)-1))<1e-8,'LTspice divider differs from analytical 1 V.');
end
result = struct('schema','radia.gpu.installed-acceptance.v1','passed',true, ...
    'matlab_release',version('-release'),'cuda_residual',info.true_relative_residual, ...
    'simulink_repeated_runs',2,'simulink_failure_cleanup',true, ...
    'ltspice_checked',checkLTspice);
file = fopen(fullfile(output,'result.json'),'w'); assert(file>=0);
fprintf(file,'%s\n',jsonencode(result)); fclose(file);
clear cleanup
end

function closeOwned(model)
if bdIsLoaded(model), close_system(model,0); end
end
