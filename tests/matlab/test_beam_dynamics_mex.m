function tests = test_beam_dynamics_mex
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabRoot = fullfile(repoRoot,"matlab");
if ~contains(path,matlabRoot)
    addpath(matlabRoot);
end
testCase.TestData.ElementaryChargeC = 1.602176634e-19;
end

function testCommandsArePublished(testCase)
commands = string(radia.internal.callMex("api.commands"));
expected = ["beam.species.proton";"beam.species.electron"; ...
    "beam.reference_particle.from_kinetic_energy_ev"; ...
    "beam.field.sample";"beam.equation.rhs";"beam.step";"beam.track"];
verifyTrue(testCase,all(ismember(expected,commands)));
end

function testReferenceParticleAndFieldSample(testCase)
species = radia.beam.proton();
reference = radia.beam.referenceParticle(species,220e6);
verifyEqual(testCase,reference.kinetic_energy_j, ...
    220e6*testCase.TestData.ElementaryChargeC,"RelTol",1e-15);
verifyEqual(testCase,reference.magnetic_rigidity_t_m, ...
    reference.momentum_kg_m_s/species.charge_c,"RelTol",1e-15);
electron = radia.beam.referenceParticle(radia.beam.electron(),1e6);
verifyLessThan(testCase,electron.magnetic_rigidity_t_m,0);

field = radia.beam.uniformField([0 0 0.7]);
sample = radia.beam.sampleField(field,[1 2 3]);
verifyEqual(testCase,sample.magnetic_t,[0 0 0.7]);
verifyEqual(testCase,sample.electric_v_m,[0 0 0]);
verifyEqual(testCase,string(sample.domain_status),"inside");
end

function testRHSAndStepAreIndependentlyInspectable(testCase)
species = radia.beam.proton();
field = radia.beam.uniformField([0 0 0.7]);
state = radia.beam.cartesianState([0 0 0],[2e-19 0 0]);
rhs = radia.beam.lorentzRHS(species,state,field);
verifyLessThan(testCase,rhs.dkinetic_momentum_kg_m_s(2),0);

result = radia.beam.step(species,state,field,1e-12, ...
    Stepper="classical-rk4");
verifyEqual(testCase,result.independent_after,1e-12,"AbsTol",1e-27);
verifyEqual(testCase,result.state_after.time_s,1e-12,"AbsTol",1e-27);
verifyEqual(testCase,string(result.backend),"native-cpp-mex");
end

function testPathLengthRHS(testCase)
species = radia.beam.proton();
state = radia.beam.cartesianState([0 0 0],[0 3e-19 4e-19]);
rhs = radia.beam.lorentzRHS(species,state,radia.beam.zeroField(), ...
    Independent="path_length");
verifyEqual(testCase,rhs.dposition_m,[0 0.6 0.8],"AbsTol",1e-15);
verifyEqual(testCase,rhs.dpath_length_m,1,"AbsTol",2e-15);
end

function testBorisUniformFieldTrajectory(testCase)
species = radia.beam.proton();
field = radia.beam.uniformField([0 0 0.7]);
momentum = 3e-19;
state = radia.beam.cartesianState([0 0 0],[momentum 0 0]);
stepProbe = radia.beam.step(species,state,field,1e-15);
gamma = stepProbe.invariants_before.relativistic_gamma;
omega = species.charge_c*0.7/(gamma*species.rest_mass_kg);
period = 2*pi/omega;
trajectory = radia.beam.track( ...
    species,state,field,0,period,period/2000,MaximumSteps=2100);

finalState = trajectory.samples{end};
verifyEqual(testCase,norm(finalState.kinetic_momentum_kg_m_s),momentum, ...
    "RelTol",2e-13);
verifyLessThan(testCase,norm(finalState.position_m(1:2)),2e-5);
verifyEqual(testCase,trajectory.summary.accepted_steps,2000);
verifyTrue(testCase,trajectory.summary.momentum_conservation_applicable);
verifyLessThan(testCase, ...
    trajectory.summary.maximum_relative_momentum_error,2e-13);
verifyEqual(testCase,string(trajectory.schema), ...
    "radia.beam.trajectory.result.v1");
end

function testElectricFieldDisablesMomentumConservationMetric(testCase)
species = radia.beam.proton();
field = radia.beam.uniformField([0 0 0],ElectricVM=[1e3 0 0]);
state = radia.beam.cartesianState([0 0 0],[1e-19 0 0]);
trajectory = radia.beam.track( ...
    species,state,field,0,1e-10,1e-11,Stepper="classical-rk4");

verifyFalse(testCase,trajectory.summary.momentum_conservation_applicable);
verifyTrue(testCase,isnan( ...
    trajectory.summary.maximum_relative_momentum_error));
end

function testInvalidStepperFailsLoudly(testCase)
species = radia.beam.proton();
state = radia.beam.cartesianState([0 0 0],[1e-19 0 0]);
verifyError(testCase,@() radia.beam.step( ...
    species,state,radia.beam.zeroField(),1e-3,Stepper="bad"), ...
    "MATLAB:validators:mustBeMember");
end
