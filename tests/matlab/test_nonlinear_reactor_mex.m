function tests = test_nonlinear_reactor_mex
%TEST_NONLINEAR_REACTOR_MEX Verify the reduced HDiv-MMM native runtime.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
addpath(fullfile(repoRoot,"matlab"));
testCase.TestData.SetupInfo = radia.setup(Force=true);
end

function testCommandCatalogAndToroidalSaturation(testCase)
commands = string(radia.internal.callMex("api.commands"));
required = ["reactor.create","reactor.output","reactor.update", ...
    "reactor.snapshot","reactor.restore","reactor.reset", ...
    "reactor.info","reactor.destroy"];
verifyTrue(testCase,all(ismember(required,commands)));

before = radia.apiInfo();
config = radia.simulink.makeNonlinearReactorDemoConfig();
handle = radia.internal.callMex("reactor.create",config);
cleanup = onCleanup(@() destroyHandle(handle));
created = radia.apiInfo();
verifyEqual(testCase,created.reactor_handle_count, ...
    before.reactor_handle_count+1);
verifyEqual(testCase,created.handle_count,before.handle_count+1);

zero = radia.internal.callMex("reactor.output",handle,0.0);
moderate = radia.internal.callMex("reactor.output",handle,5.0);
saturated = radia.internal.callMex("reactor.output",handle,20.0);
verifyEqual(testCase,zero.peak_flux_density_T,0,"AbsTol",1e-15);
verifyGreaterThan(testCase,moderate.peak_flux_density_T,1.0);
verifyGreaterThan(testCase,saturated.peak_flux_density_T, ...
    moderate.peak_flux_density_T);
verifyLessThan(testCase,saturated.differential_inductance_H, ...
    moderate.differential_inductance_H);
verifyLessThan(testCase,moderate.differential_inductance_H, ...
    zero.differential_inductance_H);
verifyLessThan(testCase,saturated.residual_relative_norm,1e-10);
verifyLessThanOrEqual(testCase,saturated.nonlinear_iterations, ...
    config.max_iterations);
verifySize(testCase,saturated.flux_density_T,[config.n_samples,1]);
verifyTrue(testCase,all(isfinite(saturated.flux_density_T)));

radia.internal.callMex("reactor.destroy",handle);
clear cleanup
verifyError(testCase,@() radia.internal.callMex( ...
    "reactor.output",handle,1.0),"radia:mex:Exception");
after = radia.apiInfo();
verifyEqual(testCase,after.reactor_handle_count,before.reactor_handle_count);
verifyEqual(testCase,after.handle_count,before.handle_count);
end

function testLinearConstitutiveAnalyticResult(testCase)
mu0 = 4*pi*1e-7;
relativePermeability = 100;
airInductance = 2e-6;
resistance = 0.25;
sampleTime = 0.01;
fieldStrength = [0;1e3;1e6];
fluxDensity = mu0*relativePermeability*fieldStrength;
config = radia.simulink.makeNonlinearReactorConfig( ...
    0,reshape([1,0,0],[1,1,3]),1,1, ...
    [fieldStrength,fluxDensity], ...
    AirInductance_H=airInductance, ...
    WindingResistance_Ohm=resistance,SampleTime_s=sampleTime);
handle = radia.internal.callMex("reactor.create",config);
cleanup = onCleanup(@() destroyHandle(handle));

current = 3;
actual = radia.internal.callMex("reactor.output",handle,current);
expectedFlux = (airInductance+mu0*(relativePermeability-1))*current;
expectedInductance = airInductance+mu0*(relativePermeability-1);
expectedFluxDensity = mu0*relativePermeability*current;
expectedVoltage = resistance*current+expectedFlux/sampleTime;
verifyEqual(testCase,actual.flux_linkage_Wb_turn,expectedFlux, ...
    "RelTol",2e-13);
verifyEqual(testCase,actual.differential_inductance_H, ...
    expectedInductance,"RelTol",2e-13);
verifyEqual(testCase,actual.flux_density_T,expectedFluxDensity, ...
    "RelTol",2e-13);
verifyEqual(testCase,actual.voltage_V,expectedVoltage,"RelTol",2e-13);

radia.internal.callMex("reactor.update",handle,current);
steady = radia.internal.callMex("reactor.output",handle,current);
verifyEqual(testCase,steady.voltage_V,resistance*current, ...
    "AbsTol",2e-13);
end

function testOutputIsSideEffectFreeAndSnapshotRestores(testCase)
config = radia.simulink.makeNonlinearReactorDemoConfig();
handle = radia.internal.callMex("reactor.create",config);
cleanup = onCleanup(@() destroyHandle(handle));
initial = radia.internal.callMex("reactor.snapshot",handle);

first = radia.internal.callMex("reactor.output",handle,8.0);
second = radia.internal.callMex("reactor.output",handle,8.0);
unchanged = radia.internal.callMex("reactor.snapshot",handle);
verifyEqual(testCase,second,first);
verifyEqual(testCase,unchanged,initial);

radia.internal.callMex("reactor.update",handle,8.0);
advanced = radia.internal.callMex("reactor.snapshot",handle);
verifyEqual(testCase,advanced.accepted_steps,initial.accepted_steps+1);
verifyNotEqual(testCase,advanced.previous_flux_linkage_Wb_turn, ...
    initial.previous_flux_linkage_Wb_turn);
radia.internal.callMex("reactor.restore",handle,initial);
restored = radia.internal.callMex("reactor.snapshot",handle);
verifyEqual(testCase,restored,initial);
radia.internal.callMex("reactor.reset",handle);
verifyEqual(testCase,radia.internal.callMex("reactor.snapshot",handle),initial);
end

function testInvalidAbiFailsWithoutLeakingHandle(testCase)
before = radia.apiInfo();
config = radia.simulink.makeNonlinearReactorDemoConfig();

invalid = config;
invalid.schema = string(config.schema);
verifyError(testCase,@() radia.internal.callMex( ...
    "reactor.create",invalid),"radia:mex:Exception");

invalid = config;
invalid.demag_row_major = [invalid.demag_row_major,0];
verifyError(testCase,@() radia.internal.callMex( ...
    "reactor.create",invalid),"radia:mex:Exception");

invalid = config;
invalid.python_per_step = true;
verifyError(testCase,@() radia.internal.callMex( ...
    "reactor.create",invalid),"radia:mex:Exception");
verifyEqual(testCase,radia.apiInfo().reactor_handle_count, ...
    before.reactor_handle_count);
end

function testProductionSimulinkModel(testCase)
model = "radia_nonlinear_reactor";
if bdIsLoaded(model)
    close_system(model,0);
end
cleanup = onCleanup(@() closeModel(model));
before = radia.apiInfo();
in = Simulink.SimulationInput(model);
in = in.setModelParameter("StopTime","0.02","SimulationMode","normal");
out = sim(in);
current = out.logsout.get("current_A").Values;
fluxDensity = out.logsout.get("Bpeak_T").Values;
inductance = out.logsout.get("Ldiff_H").Values;
residual = out.logsout.get("residual").Values;
after = radia.apiInfo();

verifyEqual(testCase,max(abs(current.Data)),20,"AbsTol",1e-12);
verifyGreaterThan(testCase,max(fluxDensity.Data),2);
verifyLessThan(testCase,min(inductance.Data),1e-4);
verifyLessThan(testCase,max(residual.Data),1e-10);
verifyEqual(testCase,after.reactor_handle_count, ...
    before.reactor_handle_count);
verifyEqual(testCase,after.handle_count,before.handle_count);
end

function destroyHandle(handle)
try
    radia.internal.callMex("reactor.destroy",handle);
catch
end
end

function closeModel(model)
if bdIsLoaded(model)
    close_system(model,0);
end
end
