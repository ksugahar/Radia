function tests = test_radia_ngsolve_parity
% MATLAB equivalents of the important Python Radia-NGSolve numerical tests.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
addpath(fullfile(repoRoot, "matlab"));
radia.setup();
testCase.TestData.meshPath = writeUnitTetra();
end

function teardownOnce(testCase)
path = testCase.TestData.meshPath;
if isfile(path)
    delete(path);
end
end

function setup(~)
radia.UtiDelAll();
end

function teardown(~)
radia.UtiDelAll();
end

function testRadiaFieldApiAndAllFieldTypes(testCase)
% Python parity: test_ngsolve_integration and test_rad_ngsolve field APIs.
vertices = [0,0,0; 1,0,0; 0,1,0; 0,0,1];
object = radia.ObjTetrahedron(vertices, [0.2,-0.1,0.8]);
points = [0.15,0.20,0.10; 0.25,0.10,0.15];
fieldTypes = ["b", "h", "a", "m", "phi"];

for fieldType = fieldTypes
    field = radia.RadiaField(object, fieldType, Units="m");
    fieldMetadata = field.fieldInfo();
    coefficientMetadata = field.info();
    expectedDimension = 3;
    if fieldType == "phi"
        expectedDimension = 1;
    end

    verifyEqual(testCase, string(fieldMetadata.field_type), fieldType);
    verifyEqual(testCase, fieldMetadata.radia_obj, object);
    verifyEqual(testCase, coefficientMetadata.dimension, expectedDimension);
    actual = field.evaluate(testCase.TestData.meshPath, points);
    reference = radia.Fld(object, fieldType, points);
    verifyEqual(testCase, actual, reference, ...
        RelTol=3e-13, AbsTol=2e-14);
    verifyTrue(testCase, all(isfinite(actual), "all"));
    delete(field);
end
end

function testCurrentAndMagnetFieldRelations(testCase)
% Python parity: test_radia_ngsolve_fields coil and magnet checks.
mu0 = 4*pi*1e-7;
radius = 0.05;
width = 1e-3;
height = 1e-3;
currentDensity = 1000/(width*height);
coil = radia.ObjArcCur([0,0,0], ...
    [radius-width/2, radius+width/2], [-pi,pi], height, ...
    100, "m", "z", currentDensity);
axisPoint = [0,0,0.05];

B = radia.Fld(coil, "b", axisPoint);
H = radia.Fld(coil, "h", axisPoint);
A = radia.Fld(coil, "a", axisPoint);
M = radia.Fld(coil, "m", axisPoint);
phi = radia.Fld(coil, "phi", axisPoint);
verifyGreaterThan(testCase, abs(B(3)), 1e-6);
verifyLessThan(testCase, norm(B-mu0*H)/norm(B), 1e-2);
verifyLessThan(testCase, abs(A(3)), 1e-10);
verifyLessThan(testCase, norm(M), 1e-10);
verifyTrue(testCase, isfinite(phi));

radia.UtiDelAll();
vertices = [-0.005,-0.005,-0.005; 0.005,-0.005,-0.005; ...
    -0.005,0.005,-0.005; -0.005,-0.005,0.005];
magnet = radia.ObjTetrahedron(vertices, [0,0,795775]);
probe = [0,0,0.02];
magnetField = radia.Fld(magnet, "b", probe);
verifyGreaterThan(testCase, abs(magnetField(3)), 1e-8);
verifyLessThan(testCase, norm(radia.Fld(magnet, "m", probe)), 1e-6);
end

function testNGSolveSpaceGridFunctionCreation(testCase)
% Python parity: TestNGSolveFunctionSpaces for Radia-owned FE workflows.
mesh = radia.ngsolve.Mesh.create(testCase.TestData.meshPath);
spaces = ["h1", "hcurl", "hdiv"];
expectedDofs = [10, 14, 30];
for index = 1:numel(spaces)
    space = radia.ngsolve.FESpace.create(mesh, spaces(index), 2);
    gridFunction = radia.ngsolve.GridFunction.fromFESpace(space);
    verifyEqual(testCase, string(gridFunction.Space), spaces(index));
    verifyEqual(testCase, gridFunction.DofCount, expectedDofs(index));
    verifyEqual(testCase, gridFunction.vector(), ...
        zeros(expectedDofs(index),1));
    delete(gridFunction);
    delete(space);
end
delete(mesh);
end

function testHDivGridFunctionProjectionAndFieldAccuracy(testCase)
% Python parity: HDiv projection, direct evaluation, and far-field accuracy.
mesh = radia.ngsolve.Mesh.create(testCase.TestData.meshPath);
space = radia.ngsolve.FESpace.create(mesh, "hdiv", 2);
gridFunction = radia.ngsolve.GridFunction.fromFESpace(space);
constantSource = radia.ObjBckg([0.2,-0.1,0.8]);
constantField = radia.RadiaField(constantSource, "b");
gridFunction.interpolate(constantField);
projected = gridFunction.asCoefficient();
points = [0.10,0.10,0.10; 0.20,0.10,0.10; 0.10,0.20,0.10];
verifyEqual(testCase, projected.evaluate(testCase.TestData.meshPath, points), ...
    constantField.evaluate(testCase.TestData.meshPath, points), ...
    RelTol=2e-12, AbsTol=2e-13);

farVertices = [-2.05,-0.05,-0.05; -1.95,-0.05,-0.05; ...
    -2.05,0.05,-0.05; -2.05,-0.05,0.05];
farMagnet = radia.ObjTetrahedron(farVertices, [0,0,8e5]);
farField = radia.RadiaField(farMagnet, "b");
gridFunction.interpolate(farField);
farProjected = gridFunction.asCoefficient();
actual = farProjected.evaluate(testCase.TestData.meshPath, points);
reference = radia.Fld(farMagnet, "b", points);
relativeErrors = vecnorm(actual-reference,2,2) ./ vecnorm(reference,2,2);
verifyLessThan(testCase, max(relativeErrors), 0.1);

delete(farProjected);
delete(farField);
delete(projected);
delete(constantField);
delete(gridFunction);
delete(space);
delete(mesh);
end

function testRadiaFieldVoxelCoefficient(testCase)
% Python parity: as_voxel_cf returns an evaluable NGSolve coefficient.
vertices = [0,0,0; 1,0,0; 0,1,0; 0,0,1];
object = radia.ObjTetrahedron(vertices, [0.2,-0.1,0.8]);
field = radia.RadiaField(object, "b");
mesh = radia.ngsolve.Mesh.create(testCase.TestData.meshPath);
voxel = field.asVoxelCoefficient(mesh, 11);
point = [0.194,0.194,0.194];
verifyEqual(testCase, voxel.info().dimension, 3);
verifyEqual(testCase, voxel.evaluate(testCase.TestData.meshPath, point), ...
    radia.Fld(object, "b", point), RelTol=3e-12, AbsTol=3e-13);
delete(voxel);
delete(mesh);
delete(field);
end

function testTransformedContainerBatchField(testCase)
% Python parity: transformed-container batch evaluation crash regression.
rotationCenter = [0,-0.1,0];
rotationAngle = -0.349;
R = rotationMatrix([0,0,1], rotationAngle);
rng(0, "twister");
points = -0.2 + 0.4*rand(20000,3);

plain = radia.ObjCnt(buildCurrentContainer());
localPoints = (points-rotationCenter)*R + rotationCenter;
reference = radia.Fld(plain, "h", localPoints)*R.';

radia.UtiDelAll();
transformed = radia.ObjCnt(buildCurrentContainer());
radia.TrfOrnt(transformed, ...
    radia.TrfRot(rotationCenter, [0,0,1], rotationAngle));
actual = radia.Fld(transformed, "h", points);

for index = [1,1235,7778,20000]
    verifyRelativeVector(testCase, ...
        radia.Fld(transformed, "h", points(index,:)), actual(index,:), 1e-9);
end
scale = max(abs(reference), [], "all");
verifyGreaterThan(testCase, scale, 0);
verifyLessThan(testCase, max(abs(actual-reference), [], "all")/scale, 1e-9);
end

function testTransformedContainerLinearFormAssembly(testCase)
% Python parity: transformed RadiaField inside NGSolve assembly.
rotationCenter = [0,-0.1,0];
rotationAngle = -0.349;
container = radia.ObjCnt(buildCurrentContainer());
radia.TrfOrnt(container, ...
    radia.TrfRot(rotationCenter, [0,0,1], rotationAngle));
field = radia.RadiaField(container, "h");
mesh = radia.ngsolve.Mesh.create(testCase.TestData.meshPath);
space = radia.ngsolve.FESpace.create(mesh, "hdiv", 2);
form = radia.ngsolve.LinearForm.createFromCoefficient(space, field);
rhs = form.vector();
values = rhs.values();
verifyTrue(testCase, all(isfinite(values)));
verifyGreaterThan(testCase, norm(values), 0);

points = [0.05,0.02,0.01; 0.10,0.05,0.04; 0.02,0.15,0.02];
verifyEqual(testCase, field.evaluate(testCase.TestData.meshPath, points), ...
    radia.Fld(container, "h", points), RelTol=1e-9, AbsTol=1e-9);

delete(rhs);
delete(form);
delete(space);
delete(mesh);
delete(field);
end

function testTransformedPolyhedronAndContainerFields(testCase)
% Python parity: TrfOrnt must transform polyhedron and container fields.
vertices = [-0.005,-0.005,-0.005; 0.005,-0.005,-0.005; ...
    -0.005,0.005,-0.005; -0.005,-0.005,0.005];
magnet = radia.ObjTetrahedron(vertices, [0,0,954930]);
probe = [0.02,0.003,0.015];
reference = radia.Fld(magnet, "b", probe);
shift = [0.05,0,0];
radia.TrfOrnt(magnet, radia.TrfTrsl(shift));
verifyRelativeVector(testCase, ...
    radia.Fld(magnet, "b", probe+shift), reference, 1e-10);

radia.UtiDelAll();
offset = [0.02,0,0];
rotatedVertices = vertices + offset;
magnet = radia.ObjTetrahedron(rotatedVertices, [954930,0,0]);
probe = [0.05,0.004,0.006];
reference = radia.Fld(magnet, "b", probe);
R = rotationMatrix([0,0,1], pi/2);
radia.TrfOrnt(magnet, radia.TrfRot([0,0,0], [0,0,1], pi/2));
verifyRelativeVector(testCase, ...
    radia.Fld(magnet, "b", (R*probe.').'), (R*reference.').', 1e-10);

radia.UtiDelAll();
first = radia.ObjTetrahedron(vertices, [0,0,954930]);
second = radia.ObjTetrahedron(vertices+[0.02,0,0], [0,0,954930]);
container = radia.ObjCnt([first,second]);
probe = [0.01,0.005,0.02];
reference = radia.Fld(container, "b", probe);
shift = [0,0.06,0];
radia.TrfOrnt(container, radia.TrfTrsl(shift));
verifyRelativeVector(testCase, ...
    radia.Fld(container, "b", probe+shift), reference, 1e-10);
end

function container = buildCurrentContainer()
objects = zeros(1,4);
objects(1) = radia.TrfOrnt( ...
    radia.ObjRecCur([0,0,0], [0.01,0.2,0.01], [0,3e6,0]), ...
    radia.TrfTrsl([0.08,0,0.03]));
objects(2) = radia.TrfOrnt( ...
    radia.ObjRecCur([0,0,0], [0.01,0.2,0.01], [0,-3e6,0]), ...
    radia.TrfTrsl([-0.08,0,0.03]));
objects(3) = radia.TrfOrnt( ...
    radia.ObjArcCur([0,0,0], [0.075,0.085], [0,pi], ...
    0.01, 10, "a", "z", 3e6), radia.TrfTrsl([0,0.1,0.03]));
objects(4) = radia.TrfOrnt( ...
    radia.ObjArcCur([0,0,0], [0.075,0.085], [pi,2*pi], ...
    0.01, 10, "a", "z", 3e6), radia.TrfTrsl([0,-0.1,0.03]));
container = radia.ObjCnt(objects);
end

function R = rotationMatrix(axis, angle)
axis = axis(:)/norm(axis);
c = cos(angle);
s = sin(angle);
K = [0,-axis(3),axis(2); axis(3),0,-axis(1); ...
    -axis(2),axis(1),0];
R = eye(3)*c + s*K + (1-c)*(axis*axis.');
end

function verifyRelativeVector(testCase, actual, reference, tolerance)
scale = norm(reference);
verifyGreaterThan(testCase, scale, 0);
verifyLessThanOrEqual(testCase, norm(actual-reference)/scale, tolerance);
end

function path = writeUnitTetra()
if ispc && isfolder("C:\temp")
    scratch = "C:\temp";
else
    scratch = string(tempdir);
end
path = string(tempname(scratch)) + ".vol";
lines = [
    "mesh3d"
    "dimension"
    "3"
    "geomtype"
    "0"
    "facedescriptors"
    "1"
    "1 1 0 1 1"
    "surfaceelements"
    "4"
    "1 1 1 0 3 1 2 3"
    "1 1 1 0 3 1 4 2"
    "1 1 1 0 3 2 4 3"
    "1 1 1 0 3 3 4 1"
    "volumeelements"
    "1"
    "1 4 1 2 3 4"
    "points"
    "4"
    "0 0 0"
    "1 0 0"
    "0 1 0"
    "0 0 1"
    "pointelements"
    "0"
    "materials"
    "1"
    "1 air"
    "bcnames"
    "1"
    "1 outer"
    "endmesh"
    ];
file = fopen(path, "w");
if file < 0
    error("radia:test:MeshWrite", "Could not create %s", path);
end
cleanup = onCleanup(@() fclose(file));
fprintf(file, "%s\n", lines);
clear cleanup
end
