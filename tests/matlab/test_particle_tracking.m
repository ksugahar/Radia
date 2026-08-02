function tests=test_particle_tracking
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabDirectory=fullfile(root,"matlab");
entries=string(strsplit(path,pathsep));
testCase.TestData.RemoveMatlabDirectory= ...
    ~any(strcmpi(entries,string(matlabDirectory)));
testCase.TestData.Path=matlabDirectory;
if testCase.TestData.RemoveMatlabDirectory
    addpath(matlabDirectory);
end
end

function teardownOnce(testCase)
if testCase.TestData.RemoveMatlabDirectory
    rmpath(testCase.TestData.Path);
end
end

function testRelativisticSpeedStaysBelowLightSpeed(testCase)
electronCharge=-1.602176634e-19;
electronMass=9.1093837139e-31;
speed=radia.tracking.speedFromKineticVoltage( ...
    electronCharge,electronMass,5e6,Relativistic=true);
verifyLessThan(testCase,speed,299792458.0);
end

function testStaticMagneticRhsUsesChargeOverGammaMass(testCase)
charge=1.602176634e-19;
mass=1.67262192595e-27;
velocity=[0.55*299792458.0;0.1*299792458.0;0.0];
gamma=1/sqrt(1-dot(velocity,velocity)/299792458.0^2);
normalizedMomentum=gamma*velocity/299792458.0;
magnetic=@(~,~,~)[0.0;0.0;1.2];
rhs=radia.tracking.relativisticLorentzRhs(0, ...
    [zeros(3,1);normalizedMomentum], ...
    charge,mass,@zeroField,magnetic,true);
expected=charge*cross(velocity,[0.0;0.0;1.2])/(mass*299792458.0);
verifyEqual(testCase,rhs(4:6),expected,"RelTol",5e-15);
verifyEqual(testCase,rhs(1:3),velocity,"RelTol",5e-15);
end

function testMagneticOnlyTrackingConservesRelativisticEnergy(testCase)
charge=-1.602176634e-19;
mass=9.1093837139e-31;
velocity=[0.4*299792458.0;0.0;0.0];
magnetic=@(~,~,~)[0.0;0.0;2e-3];
result=radia.tracking.trackLorentz(charge,mass,zeros(3,1),velocity, ...
    linspace(0,2e-8,101).',MagneticField=magnetic,Relativistic=true, ...
    RelativeTolerance=2e-11,AbsoluteTolerance=1e-14);
verifyLessThan(testCase,result.maximum_relative_kinetic_energy_drift,2e-9);
verifyTrue(testCase,result.success);
verifyEqual(testCase,result.units.velocity,"m/s");
verifySize(testCase,result.velocity_m_s,[101,3]);
verifySize(testCase,result.kinetic_energy_j,[101,1]);
verifyEmpty(testCase,result.stop_event);
end

function testTwoMomentaMeetSamePlanePointWithoutField(testCase)
plane=struct("point_m",[0.5;0;0],"normal",[1;0;0],"direction",1);
result=radia.tracking.twoMomentumExitDispersion( ...
    -1.602176634e-19,9.1093837139e-31,[0;0.1;0],[2e6;0;0], ...
    linspace(0,1e-6,101).',plane,RelativeMomentumOffset=1e-3, ...
    TransverseDirection=[0;1;0],Relativistic=true);
verifyEqual(testCase,result.eta_m,0.0,"AbsTol",1e-10);
verifyEqual(testCase,result.coincident_exit_error_m,0.0,"AbsTol",1e-12);
verifyEqual(testCase,result.minus_track.stop_event.face,"plane");
verifySize(testCase,result.minus_track.stop_event.velocity_m_s,[3,1]);
end

function testFiveMomentumFitRecoversQuadraticOptics(testCase)
offsets=[-1e-3;-5e-4;0;5e-4;1e-3];
positions=2e-4+3e-3*offsets+4*offsets.^2;
angles=-1e-4-5e-4*offsets+0.2*offsets.^2;
result=radia.tracking.fitFiveMomentumExitOptics( ...
    offsets,positions,angles);
verifyEqual(testCase,result.schema,"radia-five-momentum-exit-optics/v1");
verifyEqual(testCase,result.linear_regression_weights, ...
    [-400,-200,0,200,400],"AbsTol",2e-12);
verifyEqual(testCase,result.x0_m,2e-4,"AbsTol",1e-15);
verifyEqual(testCase,result.psi0_rad,-1e-4,"AbsTol",1e-15);
verifyEqual(testCase,result.eta_m,3e-3,"AbsTol",1e-12);
verifyEqual(testCase,result.eta_prime_rad,-5e-4,"AbsTol",1e-12);
verifyEqual(testCase,result.x_quadratic_m,4,"AbsTol",1e-9);
verifyEqual(testCase,result.psi_quadratic_rad,0.2,"AbsTol",1e-9);
verifyLessThan(testCase,result.max_x_residual_m,1e-14);
verifyLessThan(testCase,result.max_psi_residual_rad,1e-14);
verifyTrue(testCase,result.pass_all);
end

function testFiveMomentaHaveZeroExitOpticsWithoutField(testCase)
plane=struct("point_m",[0.5;0;0],"normal",[1;0;0],"direction",1);
result=radia.tracking.trackFiveMomentumExitOptics( ...
    -1.602176634e-19,9.1093837139e-31,[0;0.1;0],[2;0;0], ...
    linspace(0,1,21).',plane,ReferenceExitPointM=[0.5;0.1;0], ...
    TransverseDirection=[0;1;0],LongitudinalDirection=[1;0;0]);
verifyEqual(testCase,result.x0_m,0,"AbsTol",1e-12);
verifyEqual(testCase,result.psi0_rad,0,"AbsTol",1e-12);
verifyEqual(testCase,result.eta_m,0,"AbsTol",1e-12);
verifyEqual(testCase,result.eta_prime_rad,0,"AbsTol",1e-12);
verifyEqual(testCase,numel(result.tracks),5);
verifyTrue(testCase,result.pass_all);
for index=1:5
    verifyEqual(testCase,result.tracks(index).stop_event.face,"plane");
end
end

function testFiveMomentumFitRejectsInvalidContracts(testCase)
verifyError(testCase,@()radia.tracking.fitFiveMomentumExitOptics( ...
    [-1e-3;0;1e-3],zeros(3,1),zeros(3,1)), ...
    "radia:tracking:FiveMomentumSamples");
verifyError(testCase,@()radia.tracking.fitFiveMomentumExitOptics( ...
    [-1e-3;0;0;5e-4;1e-3],zeros(5,1),zeros(5,1)), ...
    "radia:tracking:MomentumOffsets");
verifyError(testCase,@()radia.tracking.fitFiveMomentumExitOptics( ...
    [-1e-3;-5e-4;0;5e-4;1e-3],zeros(5,1),zeros(5,1), ...
    EtaLimitM=NaN),"radia:tracking:AcceptanceLimits");
end

function testFiveMomentumTrackingHonorsStopBox(testCase)
plane=struct("point_m",[0.5;0;0],"normal",[1;0;0],"direction",1);
box=struct("minimum_m",[-1;-1;-1],"maximum_m",[0.25;1;1]);
operation=@()radia.tracking.trackFiveMomentumExitOptics( ...
    -1.602176634e-19,9.1093837139e-31,[0;0;0],[2;0;0], ...
    linspace(0,1,21).',plane,ReferenceExitPointM=[0.5;0;0], ...
    TransverseDirection=[0;1;0],StopBox=box);
verifyError(testCase,operation,"radia:tracking:ExitNotReached");
end

function testFiveMomentumTrackingRejectsParallelDirections(testCase)
plane=struct("point_m",[0.5;0;0],"normal",[1;0;0],"direction",1);
operation=@()radia.tracking.trackFiveMomentumExitOptics( ...
    -1.602176634e-19,9.1093837139e-31,[0;0;0],[2;0;0], ...
    linspace(0,1,21).',plane,ReferenceExitPointM=[0.5;0;0], ...
    TransverseDirection=[2;0;0]);
verifyError(testCase,operation,"radia:tracking:TransverseDirection");
end

function testComponentFieldAdapter(testCase)
field=radia.tracking.fieldFromComponents( ...
    @(x,y,z)x+y+z,@(x,y,z)x-y,@(x,y,z)z);
verifyEqual(testCase,field(1,2,3),[6;-1;3]);
end

function value=zeroField(~,~,~)
value=zeros(3,1);
end
