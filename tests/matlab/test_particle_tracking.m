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

function testComponentFieldAdapter(testCase)
field=radia.tracking.fieldFromComponents( ...
    @(x,y,z)x+y+z,@(x,y,z)x-y,@(x,y,z)z);
verifyEqual(testCase,field(1,2,3),[6;-1;3]);
end

function value=zeroField(~,~,~)
value=zeros(3,1);
end
