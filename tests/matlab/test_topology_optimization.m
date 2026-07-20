function tests=test_topology_optimization
tests=functiontests(localfunctions);
end
function setupOnce(testCase)
root=fileparts(fileparts(fileparts(mfilename("fullpath")))); path=fullfile(root,"matlab"); addpath(path); testCase.TestData.Path=path;
testCase.TestData.MeshPath=writeUnitTetra();
testCase.TestData.MultiMeshPath=writeTetraMesh(10);
end
function teardownOnce(testCase)
rmpath(testCase.TestData.Path);
if isfile(testCase.TestData.MeshPath), delete(testCase.TestData.MeshPath); end
if isfile(testCase.TestData.MultiMeshPath), delete(testCase.TestData.MultiMeshPath); end
end
function testVIMLinearization(testCase)
A=[3 -1;-1 2]; b=[1;0.5]; C=[1 2]; dA=zeros(2,2,2); dA(1,:,:)=[0.4 0;0 0]; dA(2,:,:)=[0 0;0 0.3];
r=radia.topopt.linearizeVIM(A,b,C,dA); epsilon=1e-7;
for k=1:2, shifted=(A+epsilon*squeeze(dA(k,:,:)))\b; observed=(C*shifted-r.response)/epsilon; verifyEqual(testCase,observed,r.response_jacobian(:,k),'RelTol',2e-6,'AbsTol',2e-8); end
end
function testLPAndCubitJournal(testCase)
r=radia.topopt.solveLPUpdate([0.5;0.5;0.5],[-3;-1;2],[1;1;1],1.5,MoveLimit=0.2);
verifyLessThanOrEqual(testCase,sum(r.density),1.5+1e-12); verifyLessThanOrEqual(testCase,max(abs(r.delta)),0.2+1e-12);
path="C:\temp\radia_topopt_density.jou"; info=radia.topopt.writeCubitJournal(path,[11;12;13;14],[0.9;0.1;0.7;0.2]);
verifyEqual(testCase,info.solid_count,2); verifyTrue(testCase,contains(string(fileread(path)),"add hex 11 13"));
end
function testSequentialVIMLP(testCase)
linearize=@localLinearize; r=radia.topopt.optimizeVIMLP(0.5*ones(3,1),ones(3,1),0.5,linearize, ...
 ObjectiveWeights=[1;0],MoveLimit=0.25,MaxIterations=5);
verifyLessThanOrEqual(testCase,sum(r.density),1.5+1e-12); verifyGreaterThanOrEqual(testCase,r.density(1),r.density(3));
end

function testSheetMetalLPAndMeshRouting(testCase)
n=3; u=zeros(n,1); t=ones(n,1); rho=ones(n,1); area=ones(n,1); L=[1,-2,1];
update=radia.topopt.solveSheetMetalLP(u,t,rho,[-ones(n,1);ones(n,1);ones(n,1)],area, ...
 VolumeMax=2.7,DisplacementMove=radia.topopt.localTrustRegion([1;2;1]),ThicknessMove=0.2,ActivationMove=0.1, ...
 ThicknessBounds=[0.8,1.2],Laplacian=L,CurvatureLimit=0.05);
verifyLessThanOrEqual(testCase,max(abs(L*update.normal_displacement)),0.05+1e-12);
verifyGreaterThanOrEqual(testCase,min(update.thickness),0.8-1e-12);
deform=radia.topopt.routeMeshUpdate([0.9;0.8],[2;3],[0.02;0.04]);
refine=radia.topopt.routeMeshUpdate([0.9;0.3],[2;10],[0.02;0.3]);
rebuild=radia.topopt.routeMeshUpdate([0.9;0.8],[2;3],[0.02;0.6]);
verifyEqual(testCase,deform.route,"ngsolve_deform");
verifyEqual(testCase,refine.route,"ngsolve_refine");
verifyEqual(testCase,refine.refine_elements,1);
verifyEqual(testCase,rebuild.route,"cubit_rebuild");
accepted=radia.topopt.acceptTrafoStep(@localQuality,0.08*ones(2,1));
verifyTrue(testCase,accepted.accepted); verifyEqual(testCase,accepted.scale,1);
end
function quality=localQuality(scale) %#ok<INUSD>
quality=struct("jacobian_determinants",[0.9;0.8],"jacobian_conditions",[2;3]);
end
function testMagneticShieldRMSGradient(testCase)
response=[3;4]; jacobian=[1,0,2;0,2,1];
gradient=radia.topopt.magneticShieldRMSGradient(response,jacobian);
expected=(response'*jacobian/(2*sqrt(mean(response.^2))))';
verifyEqual(testCase,gradient,expected,'AbsTol',1e-14);
end
function testVIMOperatorShapeDerivative(testCase)
M=[2,.2;.2,1.5]; B=[1,-1;.5,.25]; G=[.8,.1;.1,.6]; h=[2;-1];
dM=reshape([.1,.02;.02,-.03],[1,2,2]); dB=reshape([.03,0;-.01,.02],[1,2,2]); dG=reshape([.04,.01;.01,-.02],[1,2,2]);
r=radia.topopt.linearizeVIMOperator(M,B,G,h,.1,dM,dB,dG); epsilon=1e-7;
s=radia.topopt.linearizeVIMOperator(M+epsilon*squeeze(dM),B+epsilon*squeeze(dB),G+epsilon*squeeze(dG),h,.1,zeros(1,2,2),zeros(1,2,2),zeros(1,2,2));
verifyEqual(testCase,(s.matrix-r.matrix)/epsilon,squeeze(r.matrix_jacobian(1,:,:)),'RelTol',2e-6,'AbsTol',2e-9);
end

function testNativeVectorH1TrafoQuality(testCase)
mesh=radia.ngsolve.Mesh.create(testCase.TestData.MeshPath);
reference=mesh.trafoQuality(IntegrationOrder=2);
space=radia.ngsolve.FESpace.create(mesh,"vectorh1",1,NoGrads=false);
deformation=radia.ngsolve.GridFunction.fromFESpace(space,Name="deformation");
deformation.interpolate(radia.ngsolve.CoefficientFunction.constant([0,0,0]));
mesh.setDeformation(deformation);
quality=mesh.trafoQuality(IntegrationOrder=2, ...
 ReferenceDeterminants=reference.raw_jacobian_determinants);
verifyTrue(testCase,mesh.info().has_deformation);
verifyEqual(testCase,quality.jacobian_determinants,ones(1,1),'AbsTol',1e-13);
verifyEqual(testCase,quality.jacobian_conditions, ...
 reference.jacobian_conditions,'AbsTol',1e-13);
mesh.unsetDeformation();
verifyFalse(testCase,mesh.info().has_deformation);
end

function testTwoLevelTopologyCommitsHysteresisCrossing(testCase)
mesh=radia.ngsolve.Mesh.create(testCase.TestData.MeshPath);
initial=struct("mesh",mesh,"model",struct("generation",0), ...
 "normal_displacement",0,"thickness",1,"activation",0.8,"objective",0.04);
result=radia.topopt.optimizeHexSheetTopology(initial, ...
 @localRemovedActivationStep,@localZeroDeformation, ...
 @localActivationObjective,@localRebuildModel, ...
 @(request)radia.ngsolve.Mesh.create(testCase.TestData.MeshPath),1, ...
 InnerIterations=5,MinimumInnerIterations=5,MaxOuterIterations=1, ...
 WorkDirectory="C:\temp\radia_hex_topopt_matlab_test");
verifyEqual(testCase,result.inner_iteration_count,1);
verifyEqual(testCase,result.cubit_rebuild_count,1);
verifyEqual(testCase,result.hmatrix_rebuild_count,1);
verifyEqual(testCase,result.state.activation,0.2,'AbsTol',1e-14);
verifyEqual(testCase,result.history.Route(end),"cubit_rebuild");
verifyTrue(testCase,result.history.CubitRebuilt(end));
verifyTrue(testCase,result.history.HMatrixRebuilt(end));
verifyTrue(testCase,contains(result.history.CubitReason(end),"hysteresis"));
end

function testTwoLevelTopologyBatchesSparseCubitChanges(testCase)
n=10;
mesh=radia.ngsolve.Mesh.create(testCase.TestData.MultiMeshPath);
initial=struct("mesh",mesh,"model",struct("generation",0), ...
 "normal_displacement",zeros(n,1),"thickness",ones(n,1), ...
 "activation",ones(n,1),"objective",0);
result=radia.topopt.optimizeHexSheetTopology(initial, ...
 @localSparseRemovalStep,@localZeroDeformation, ...
 @localActivationObjective,@localRebuildModel, ...
 @(request)radia.ngsolve.Mesh.create(testCase.TestData.MultiMeshPath),ones(n,1), ...
 InnerIterations=5,MinimumInnerIterations=5,MaxOuterIterations=1, ...
 CubitBatchInterval=3,CubitBatchFraction=0.5, ...
 WorkDirectory="C:\temp\radia_hex_topopt_matlab_test");
verifyEqual(testCase,result.history.Route, ...
    ["ngsolve_deform";"ngsolve_deform";"cubit_rebuild"]);
verifyEqual(testCase,result.history.PendingTopologyChanges,[1;1;1]);
verifyEqual(testCase,result.cubit_rebuild_count,1);
verifyEqual(testCase,result.hmatrix_rebuild_count,1);
end

function testTwoLevelTopologySkipsUnneededCubit(testCase)
mesh=radia.ngsolve.Mesh.create(testCase.TestData.MeshPath);
initial=struct("mesh",mesh,"model",struct("generation",0), ...
 "normal_displacement",0,"thickness",1,"activation",0.8,"objective",0.04);
result=radia.topopt.optimizeHexSheetTopology(initial, ...
 @localStableActivationStep,@localZeroDeformation, ...
 @localActivationObjective,@localRebuildModel, ...
 @(request)radia.ngsolve.Mesh.create(testCase.TestData.MeshPath),1, ...
 InnerIterations=5,MinimumInnerIterations=5,MaxOuterIterations=1, ...
 WorkDirectory="C:\temp\radia_hex_topopt_matlab_test");
verifyEqual(testCase,result.inner_iteration_count,5);
verifyEqual(testCase,result.cubit_rebuild_count,0);
verifyEqual(testCase,result.hmatrix_rebuild_count,0);
verifyTrue(testCase,all(result.history.Route=="ngsolve_deform"));
end

function testTwoLevelTopologyRequiresFiveToTwentyInnerIterations(testCase)
mesh=radia.ngsolve.Mesh.create(testCase.TestData.MeshPath);
initial=struct("mesh",mesh,"model",struct,"normal_displacement",0, ...
 "thickness",1,"activation",0.8,"objective",0.04);
verifyError(testCase,@()radia.topopt.optimizeHexSheetTopology(initial, ...
 @localStableActivationStep,@localZeroDeformation,@localActivationObjective, ...
 @localRebuildModel,@(request)mesh,1,InnerIterations=4), ...
 "radia:topopt:InnerIterations");
end

function model=localLinearize(density)
model=struct("response",[density(1)-density(3);density(2)], ...
 "response_jacobian",[-1 0 1;0 1 0]);
end

function step=localRemovedActivationStep(state)
update=struct("normal_displacement",state.normal_displacement, ...
 "thickness",state.thickness,"activation",0.2*ones(size(state.activation)));
step=struct("update",update,"requires_cubit",false);
end

function step=localSparseRemovalStep(state)
activation=ones(size(state.activation)); activation(1)=0.2;
update=struct("normal_displacement",state.normal_displacement, ...
 "thickness",state.thickness,"activation",activation);
step=struct("update",update,"requires_cubit",false);
end

function step=localStableActivationStep(state)
update=struct("normal_displacement",state.normal_displacement, ...
 "thickness",state.thickness,"activation",min(0.85,state.activation+0.01));
step=struct("update",update,"requires_cubit",false);
end

function deformation=localZeroDeformation(mesh,normalDisplacement) %#ok<INUSD>
space=radia.ngsolve.FESpace.create(mesh,"vectorh1",1,NoGrads=false);
deformation=radia.ngsolve.GridFunction.fromFESpace(space,Name="topopt_deformation");
deformation.interpolate(radia.ngsolve.CoefficientFunction.constant([0,0,0]));
end

function objective=localActivationObjective(state)
objective=sum((1-state.activation).^2);
end

function model=localRebuildModel(mesh,normal,thickness,activation,route) %#ok<INUSD>
model=struct("generation",1,"route",string(route),"mesh_info",mesh.info());
end

function path=writeUnitTetra()
path=writeTetraMesh(1);
end

function path=writeTetraMesh(count)
path=string(tempname("C:\temp"))+".vol";
file=fopen(path,"w");
if file<0, error("radia:test:MeshWrite","Could not create %s",path); end
cleanup=onCleanup(@()fclose(file));
fprintf(file,"mesh3d\ndimension\n3\ngeomtype\n0\n");
fprintf(file,"facedescriptors\n1\n1 1 0 1 1\n");
fprintf(file,"surfaceelements\n%d\n",4*count);
faces=[1,2,3;1,4,2;2,4,3;3,4,1];
for cellIndex=0:count-1
    vertices=4*cellIndex+(1:4);
    for faceIndex=1:4
        face=vertices(faces(faceIndex,:));
        fprintf(file,"1 1 1 0 3 %d %d %d\n",face);
    end
end
fprintf(file,"volumeelements\n%d\n",count);
for cellIndex=0:count-1
    vertices=4*cellIndex+(1:4);
    fprintf(file,"1 4 %d %d %d %d\n",vertices);
end
fprintf(file,"points\n%d\n",4*count);
for cellIndex=0:count-1
    offset=2*cellIndex;
    fprintf(file,"%.17g 0 0\n",offset);
    fprintf(file,"%.17g 0 0\n",offset+1);
    fprintf(file,"%.17g 1 0\n",offset);
    fprintf(file,"%.17g 0 1\n",offset);
end
fprintf(file,"pointelements\n0\nmaterials\n1\n1 air\n");
fprintf(file,"bcnames\n1\n1 outer\nendmesh\n");
clear cleanup
end
