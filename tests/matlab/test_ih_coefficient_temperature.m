function tests=test_ih_coefficient_temperature
tests=functiontests(localfunctions);
end

function config=fixture()
config=radia.simulink.makeIHNativeSmokeConfig();
config.n_temperature=2;
config.rotation_mode='none';
config.temperature_representation='ngsolve-h1-coefficients';
config.temperature_constant_coefficients=[1;0];
config.temperature_cell_weights=[1;.5];
config.temperature_evaluation=struct('n_samples',3,'rows',[0;1;1;2;2], ...
    'cols',[0;0;1;0;1],'values',[1;1;.5;1;1]);
config.initial_temperature_K=[293.15;0];
config.mass_row_ptr=[0;2;4]; config.mass_col=[0;1;0;1];
config.mass_value=[1;.5;.5;1/3];
config.stiffness_row_ptr=config.mass_row_ptr; config.stiffness_col=config.mass_col;
config.stiffness_value=[0;0;0;1];
config.convection_row_ptr=config.mass_row_ptr; config.convection_col=config.mass_col;
config.convection_value=config.mass_value;
config.heat_to_temperature_projection=[1;.5];
config.convection_W_per_m2K=2;
config.thermal_tolerance=1e-13;
config.sample_time_s=.1;
end

function testConstantAmbientDoesNotHeatModalCoefficients(testCase)
config=radia.simulink.validateIHNativeConfig(fixture());
h=radia.internal.callMex('ih.thermal.create',config);
cleanup=onCleanup(@() radia.internal.callMex('ih.thermal.destroy',h));
for step=1:20
    radia.internal.callMex('ih.thermal.update',h,0,293.15,0);
end
actual=radia.internal.callMex('ih.thermal.output',h);
verifyEqual(testCase,actual,config.initial_temperature_K,AbsTol=1e-9);
verifyEqual(testCase,radia.simulink.ihTemperatureStatistics(config,actual), ...
    [293.15,293.15,293.15],AbsTol=1e-9);
end

function testCoefficientTransientMatchesIndependentSolve(testCase)
config=fixture();
h=radia.internal.callMex('ih.thermal.create',config);
cleanup=onCleanup(@() radia.internal.callMex('ih.thermal.destroy',h));
mass=[1,.5;.5,1/3]; stiffness=[0,0;0,1]; constant=[1;0];
reference=config.initial_temperature_K;
for step=1:100
    reference=(mass+.1*(stiffness+2*mass))\(mass*reference+.1*([10;5]+2*293.15*mass*constant));
    radia.internal.callMex('ih.thermal.update',h,10,293.15,0);
end
verifyEqual(testCase,radia.internal.callMex('ih.thermal.output',h),reference,AbsTol=1e-9);
end

function testModalRotationAndMissingEvaluationFail(testCase)
config=fixture();config.rotation_mode='periodic-uniform';
verifyError(testCase,@() radia.simulink.validateIHNativeConfig(config), ...
    'radia:simulink:IHConfigTemperatureRepresentation');
config=rmfield(fixture(),'temperature_evaluation');
verifyError(testCase,@() radia.simulink.validateIHNativeConfig(config), ...
    'radia:simulink:IHConfigTemperatureRepresentation');
end
