function tests=test_hcurl_topology_optimization
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(root,"matlab"));
testCase.TestData.Root=root;
testCase.TestData.MeshPath=writeUnitTetra();
end

function teardownOnce(testCase)
rmpath(fullfile(testCase.TestData.Root,"matlab"));
if isfile(testCase.TestData.MeshPath), delete(testCase.TestData.MeshPath); end
end

function testNativeHCurlOperatorAndShapeContractions(testCase)
[gram,operator,G,maps]=makeOperator();
cleanup=onCleanup(@()cleanupOperator(operator,gram));
L=zeros(2);
for component=1:3
    B=squeeze(maps(component,:,:));
    L=L+B'*G*B;
end
x=[1+0.2i;-0.4+0.3i];
verifyEqual(testCase,operator.matvec(x),L*x,"AbsTol",2e-11);
verifyEqual(testCase,operator.toDense(),L,"AbsTol",2e-11);

rho=0.7; power=1.5;
verifyEqual(testCase,operator.activationToDense(rho,Power=power), ...
    rho^(2*power)*L,"AbsTol",2e-11);
left=[0.2-0.1i;0.7+0.3i];
expectedActivation=2*power*rho^(2*power-1)*(left'*L*x);
verifyEqual(testCase,operator.activationContractions( ...
    rho,left,x,Power=power),expectedActivation,"AbsTol",3e-11);

cellVerts=unitTetraVertices();
velocity=zeros(1,1,4,3); velocity(1,1,:,:)=cellVerts;
verifyEqual(testCase,operator.directionalContractions(velocity,left,x), ...
    left'*L*x,"RelTol",2e-10,"AbsTol",2e-11);
clear cleanup
end

function testNativeMultifrequencyShapeAdjointMatchesFiniteDifference(testCase)
[gram,operator]=makeOperator();
cleanup=onCleanup(@()cleanupOperator(operator,gram));
L=operator.toDense();
R=[2.1,0.2;0.2,1.4];
resistance=struct("matrix",R,"jacobian",reshape(-R,[1,2,2]));
cellVerts=unitTetraVertices();
velocity=zeros(1,1,4,3); velocity(1,1,:,:)=cellVerts;
frequency=137; rhs=[1+0.2i;-0.3+0.1i];
result=radia.topopt.linearizeHCurlMultifrequencyJoule( ...
    operator,resistance,velocity,frequency,rhs);
step=1e-6;
plus=jouleObjective(R-step*R,L+step*L,frequency,rhs);
minus=jouleObjective(R+step*R,L-step*L,frequency,rhs);
finiteDifference=(plus-minus)/(2*step);
verifyEqual(testCase,result.gradient,finiteDifference, ...
    "RelTol",3e-6,"AbsTol",2e-9);
verifyEqual(testCase,result.objective,jouleObjective(R,L,frequency,rhs), ...
    "AbsTol",2e-13);
clear cleanup
end

function testNativeActivationAdjointMatchesFiniteDifference(testCase)
[gram,operator]=makeOperator();
cleanup=onCleanup(@()cleanupOperator(operator,gram));
cellGrams=zeros(1,2,2); cellGrams(1,:,:)=[1.8,0.15;0.15,1.1];
conductivity=struct("solid",5.0,"void",0.4,"power",3.0);
rho=0.63; frequency=91; rhs=[0.8+0.1i;-0.2+0.3i];
result=radia.topopt.linearizeHCurlActivationMultifrequencyJoule( ...
    operator,cellGrams,rho,frequency,rhs,conductivity,InductancePower=1.2);
step=1e-6;
plus=radia.topopt.linearizeHCurlActivationMultifrequencyJoule( ...
    operator,cellGrams,rho+step,frequency,rhs,conductivity,InductancePower=1.2);
minus=radia.topopt.linearizeHCurlActivationMultifrequencyJoule( ...
    operator,cellGrams,rho-step,frequency,rhs,conductivity,InductancePower=1.2);
finiteDifference=(plus.objective-minus.objective)/(2*step);
verifyEqual(testCase,result.gradient,finiteDifference, ...
    "RelTol",5e-6,"AbsTol",2e-9);
clear cleanup
end

function testActivationHCurlDriverUsesExistingTwoLevelLoop(testCase)
[gram,operator]=makeOperator();
operatorCleanup=onCleanup(@()cleanupOperator(operator,gram));
mesh=radia.ngsolve.Mesh.create(testCase.TestData.MeshPath);
meshCleanup=onCleanup(@()delete(mesh));
cellGrams=zeros(1,2,2); cellGrams(1,:,:)=[1.8,0.15;0.15,1.1];
conductivity=struct("solid",5.0,"void",0.4,"power",3.0);
frequency=91; rhs=[0.8+0.1i;-0.2+0.3i]; activation=0.63;
objective=activationObjective(operator,cellGrams,activation,frequency,rhs,conductivity);
initial=struct("mesh",mesh,"model",struct("generation",0), ...
    "normal_displacement",0,"thickness",1, ...
    "activation",activation,"objective",objective);
builder=@(state)activationInputs( ...
    state,operator,cellGrams,frequency,rhs,conductivity);
objectiveFcn=@(state)activationObjective( ...
    operator,cellGrams,state.activation,frequency,rhs,conductivity);
driverOptions=struct("InnerIterations",5,"MinimumInnerIterations",5, ...
    "MaxOuterIterations",1,"FinalizeTopology",false, ...
    "ActivationRemoveThreshold",0.1,"ActivationRestoreThreshold",0.9, ...
    "WorkDirectory","C:\temp\radia_hcurl_topopt_driver_test");
result=radia.topopt.optimizeHCurlEddyBubbleActivationHexSheet( ...
    initial,builder,@zeroDeformation,objectiveFcn, ...
    @(varargin)struct("generation",1), ...
    @(request)radia.ngsolve.Mesh.create(testCase.TestData.MeshPath),1,driverOptions);
verifyEqual(testCase,result.inner_iteration_count,5);
verifyEqual(testCase,result.cubit_rebuild_count,0);
verifyTrue(testCase,all(result.history.Route=="ngsolve_deform"));
clear meshCleanup operatorCleanup
end

function testNGSolveResistancePiolaTangentAndCellGrams(testCase)
mesh=radia.ngsolve.Mesh.create(testCase.TestData.MeshPath);
meshCleanup=onCleanup(@()delete(mesh));
space=radia.ngsolve.FESpace.create(mesh,"hcurl",1);
spaceCleanup=onCleanup(@()delete(space));
basis=eye(space.DofCount);
deformationSpace=radia.ngsolve.FESpace.create(mesh,"vectorh1",1,NoGrads=false);
deformationSpaceCleanup=onCleanup(@()delete(deformationSpace));
deformation=radia.ngsolve.GridFunction.fromFESpace( ...
    deformationSpace,Name="coordinate_deformation");
deformationCleanup=onCleanup(@()delete(deformation));
coordinates=radia.ngsolve.CoefficientFunction.coordinates(3);
coordinatesCleanup=onCleanup(@()delete(coordinates));
deformation.interpolate(coordinates);
referenceVertices=zeros(1,4,3);
referenceVertices(1,:,:)=[0,0,0;1,0,0;0,1,0;0,0,1];
sampled=radia.topopt.sampleHCurlSubtetVelocities( ...
    referenceVertices,int32(0),{deformation});
expectedMapped=[0,0,1;0,0,0;1,0,0;0,1,0];
verifyEqual(testCase,squeeze(sampled(1,1,:,:)),expectedMapped, ...
    "AbsTol",2e-14);
result=radia.topopt.assembleHCurlResistanceShapeTangents( ...
    space,basis,{deformation},Conductivity=2);
dR=squeeze(result.jacobian(1,:,:));
verifyEqual(testCase,dR,-result.matrix,"RelTol",2e-12,"AbsTol",2e-12);
cellGrams=radia.topopt.assembleHCurlCellCurlGrams(space,basis);
verifyEqual(testCase,squeeze(cellGrams(1,:,:)),2*result.matrix, ...
    "RelTol",2e-12,"AbsTol",2e-12);
clear coordinatesCleanup deformationCleanup deformationSpaceCleanup
clear spaceCleanup meshCleanup
end

function testHCurlCommandsArePublished(testCase)
commands=string(radia.internal.callMex('api.commands'));
expected=["hcurl.topopt.operator.create", ...
    "hcurl.topopt.operator.directional_contractions", ...
    "hcurl.topopt.resistance_shape_tangents", ...
    "hcurl.topopt.cell_curl_grams", ...
    "hcurl.topopt.sample_subtet_velocities", ...
    "hcurl.topopt.multifrequency_joule", ...
    "hcurl.topopt.activation_multifrequency_joule"];
verifyTrue(testCase,all(ismember(expected,commands)));
end

function [gram,operator,G,maps]=makeOperator()
cellVerts=unitTetraVertices();
a=0.5854101966249685; b=0.1381966011250105;
points=[a,b,b;b,a,b;b,b,a;b,b,b]; weights=ones(4,1)/24;
gram=radia.HACApKChargeGram.from_local_polynomials( ...
    cellVerts,1,int32([0;0]),[1,0,0,0;0,1,0,0], ...
    [0,0,0;1,0,0;0,1,0;0,0,1],points,weights, ...
    AcaEps=1e-12,LeafSize=4);
G=[gram.entry(1,1),gram.entry(1,2); ...
   gram.entry(2,1),gram.entry(2,2)];
maps=zeros(3,2,2);
maps(1,:,:)=[1.0,0.2;0.1,0.4];
maps(2,:,:)=[0.3,-0.1;0.2,0.7];
maps(3,:,:)=[-0.2,0.5;0.6,0.1];
cellTensor=zeros(1,4,3); cellTensor(1,:,:)=cellVerts;
operator=radia.topopt.HCurlTopologyOperator.create( ...
    gram,maps,cellTensor,int32([0;0]),int32(0),Mu=1);
end

function cleanupOperator(operator,gram)
delete(operator); delete(gram);
end

function value=jouleObjective(R,L,frequency,rhs)
state=(R+1i*2*pi*frequency*L)\rhs;
value=0.5*real(state'*R*state);
end

function inputs=activationInputs(state,operator,cellGrams,frequency,rhs,conductivity)
inputs=struct("operator",operator,"cellCurlGrams",cellGrams, ...
    "frequenciesHz",frequency,"rhs",rhs,"conductivity",conductivity, ...
    "area",1,"lpOptions",struct("VolumeMax",1, ...
    "ActivationMove",0.02,"InductancePower",1.2));
if numel(state.activation)~=1
    error("radia:test:ActivationShape","Expected one activation cell.");
end
end

function value=activationObjective(operator,cellGrams,activation,frequency,rhs,conductivity)
linearization=radia.topopt.linearizeHCurlActivationMultifrequencyJoule( ...
    operator,cellGrams,activation,frequency,rhs,conductivity,InductancePower=1.2);
value=linearization.objective;
end

function deformation=zeroDeformation(mesh,normalDisplacement) %#ok<INUSD>
space=radia.ngsolve.FESpace.create(mesh,"vectorh1",1,NoGrads=false);
deformation=radia.ngsolve.GridFunction.fromFESpace(space,Name="zero_deformation");
deformation.interpolate(radia.ngsolve.CoefficientFunction.constant([0,0,0]));
end

function vertices=unitTetraVertices()
vertices=[0,0,0;1,0,0;0,1,0;0,0,1];
end

function path=writeUnitTetra()
path=string(tempname("C:\temp"))+".vol";
file=fopen(path,"w");
if file<0, error("radia:test:MeshWrite","Could not create %s",path); end
cleanup=onCleanup(@()fclose(file));
fprintf(file,"mesh3d\ndimension\n3\ngeomtype\n0\n");
fprintf(file,"facedescriptors\n1\n1 1 0 1 1\n");
fprintf(file,"surfaceelements\n4\n");
fprintf(file,"1 1 1 0 3 1 2 3\n1 1 1 0 3 1 4 2\n");
fprintf(file,"1 1 1 0 3 2 4 3\n1 1 1 0 3 3 4 1\n");
fprintf(file,"volumeelements\n1\n1 4 1 2 3 4\n");
fprintf(file,"points\n4\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n");
fprintf(file,"pointelements\n0\nmaterials\n1\n1 conductor\n");
fprintf(file,"bcnames\n1\n1 outer\nendmesh\n");
end
