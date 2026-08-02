function tests = test_radia_mex
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
addpath(fullfile(repoRoot, "matlab"));
testCase.TestData.setupInfo = radia.setup(Force=true);
testCase.TestData.meshPath = writeUnitTetra();
end

function teardownOnce(testCase)
path = testCase.TestData.meshPath;
if isfile(path)
    delete(path);
end
end

function testApiAndTaskManager(testCase)
api = radia.apiInfo();
probe = radia.taskmanagerProbe(200000);
verifyEqual(testCase, api.api_version, 1);
verifyGreaterThanOrEqual(testCase, api.taskmanager_max_threads, 1);
verifyGreaterThanOrEqual(testCase, probe.used_threads, 1);
verifyTrue(testCase, isfinite(probe.checksum));
pathEntries = split(string(getenv("PATH")), pathsep);
excluded = string(testCase.TestData.setupInfo.excluded_openmp_runtime_dirs);
for index = 1:numel(excluded)
    verifyFalse(testCase, any(strcmpi(excluded(index), pathEntries)));
end
end

function testIHHandleCommandsAndLifecycle(testCase)
commands = string(radia.internal.callMex('api.commands'));
required = [ ...
    "ih.eddy.create", "ih.eddy.output", "ih.eddy.destroy", ...
    "ih.thermal.create", "ih.thermal.output", "ih.thermal.update", ...
    "ih.thermal.reset", "ih.thermal.destroy"];
verifyTrue(testCase, all(ismember(required, commands)));

before = radia.apiInfo();
config = radia.simulink.makeIHNativeSmokeConfig();
eddy = radia.internal.callMex('ih.eddy.create', config);
eddyCleanup = onCleanup(@() destroyIHHandle('ih.eddy.destroy', eddy));
created = radia.apiInfo();
verifyEqual(testCase, created.ih_handle_count, before.ih_handle_count + 1);
verifyEqual(testCase, created.handle_count, before.handle_count + 1);

heat = radia.internal.callMex( ...
    'ih.eddy.output', eddy, 2.0, 0.0, config.initial_temperature_K);
verifyEqual(testCase, heat, 4.0, "AbsTol", 1e-13);

thermal = radia.internal.callMex('ih.thermal.create', config);
thermalCleanup = onCleanup( ...
    @() destroyIHHandle('ih.thermal.destroy', thermal));
created = radia.apiInfo();
verifyEqual(testCase, created.ih_handle_count, before.ih_handle_count + 2);
verifyError(testCase, ...
    @() radia.internal.callMex('ih.thermal.output', eddy), ...
    "radia:mex:Exception");
verifyError(testCase, ...
    @() radia.internal.callMex( ...
        'ih.eddy.output', thermal, 1.0, 0.0, config.initial_temperature_K), ...
    "radia:mex:Exception");

initial = radia.internal.callMex('ih.thermal.output', thermal);
verifyEqual(testCase, initial, config.initial_temperature_K, ...
    "AbsTol", 1e-13);
radia.internal.callMex('ih.thermal.update', thermal, heat, 293.15, 0.0);
advanced = radia.internal.callMex('ih.thermal.output', thermal);
verifyEqual(testCase, advanced, 293.55, "AbsTol", 1e-12);
radia.internal.callMex('ih.thermal.reset', thermal);
verifyEqual(testCase, ...
    radia.internal.callMex('ih.thermal.output', thermal), ...
    config.initial_temperature_K, "AbsTol", 1e-13);

radia.internal.callMex('ih.eddy.destroy', eddy);
verifyError(testCase, ...
    @() radia.internal.callMex( ...
        'ih.eddy.output', eddy, 1.0, 0.0, config.initial_temperature_K), ...
    "radia:mex:Exception");
radia.internal.callMex('ih.thermal.destroy', thermal);
clear thermalCleanup eddyCleanup
after = radia.apiInfo();
verifyEqual(testCase, after.ih_handle_count, before.ih_handle_count);
verifyEqual(testCase, after.handle_count, before.handle_count);
end

function testIHHandleCreateRejectsInvalidConfig(testCase)
before = radia.apiInfo();
config = radia.simulink.makeIHNativeSmokeConfig();

invalid = config;
invalid.bh_mode = 'nonlinear';
verifyError(testCase, ...
    @() radia.internal.callMex('ih.eddy.create', invalid), ...
    "radia:mex:Exception");

invalid = config;
invalid.heat_cell_weights = 0;
verifyError(testCase, ...
    @() radia.internal.callMex('ih.eddy.create', invalid), ...
    "radia:mex:Exception");

invalid = config;
invalid.thermal_tolerance = 0;
verifyError(testCase, ...
    @() radia.internal.callMex('ih.thermal.create', invalid), ...
    "radia:mex:Exception");

invalid = config;
invalid.initial_temperature_K = 0;
verifyError(testCase, ...
    @() radia.internal.callMex('ih.thermal.create', invalid), ...
    "radia:mex:Exception");

invalid = config;
invalid.mass_row_ptr = [0; 2];
verifyError(testCase, ...
    @() radia.internal.callMex('ih.thermal.create', invalid), ...
    "radia:mex:Exception");

verifyEqual(testCase, radia.apiInfo().ih_handle_count, before.ih_handle_count);
end

function testNGSolveSpacesP1ToP6(testCase)
expectedHCurl = [6, 14, 29, 53, 88, 136];
expectedHDiv = [12, 30, 60, 105, 168, 252];
actualHCurl = zeros(1, 6);
actualHDiv = zeros(1, 6);
for order = 1:6
    info = radia.spaceInfo(testCase.TestData.meshPath, order);
    actualHCurl(order) = info.hcurl_ndof;
    actualHDiv(order) = info.hdiv_ndof;
end
verifyEqual(testCase, actualHCurl, expectedHCurl);
verifyEqual(testCase, actualHDiv, expectedHDiv);
canonical = radia.ngsolve.space_info(testCase.TestData.meshPath, 2);
verifyEqual(testCase, canonical.hcurl_ndof, expectedHCurl(2));
verifyEqual(testCase, canonical.hdiv_ndof, expectedHDiv(2));
end

function testNGSolveMatrixDump(testCase)
meshPath = testCase.TestData.meshPath;
[hcurlMass, hcurlInfo] = radia.ngsolveMatrix(meshPath, "hcurl", 2, "mass");
[canonicalMass, canonicalInfo] = radia.ngsolve.matrix_dump( ...
    meshPath, "hcurl", 2, "mass");
verifyEqual(testCase, canonicalMass, hcurlMass);
verifyEqual(testCase, canonicalInfo.dof_count, hcurlInfo.dof_count);
verifyEqual(testCase, size(hcurlMass), [14, 14]);
verifyEqual(testCase, hcurlInfo.dof_count, 14);
verifyEqual(testCase, hcurlInfo.shape, [14, 14]);
verifyEqual(testCase, hcurlInfo.nonzero_nze, nnz(hcurlMass));
verifyLessThan(testCase, norm(hcurlMass - hcurlMass', "fro"), 1e-12);
verifyTrue(testCase, all(isfinite(nonzeros(hcurlMass))));

hcurlCurlCurl = radia.ngsolveMatrix(meshPath, "hcurl", 1, "stiffness");
verifyEqual(testCase, size(hcurlCurlCurl), [6, 6]);
verifyTrue(testCase, issymmetric(hcurlCurlCurl));

hdivMass = radia.ngsolveMatrix(meshPath, "hdiv", 2, "mass");
verifyEqual(testCase, size(hdivMass), [30, 30]);
verifyLessThan(testCase, norm(hdivMass - hdivMass', "fro"), 1e-12);

hdivDivDiv = radia.ngsolveMatrix(meshPath, "hdiv", 1, "divdiv");
verifyEqual(testCase, size(hdivDivDiv), [12, 12]);
verifyLessThan(testCase, norm(hdivDivDiv - hdivDivDiv', "fro"), 1e-12);

h1Stiffness = radia.ngsolveMatrix(meshPath, "h1", 2, "stiffness");
verifyEqual(testCase, size(h1Stiffness), [10, 10]);
verifyLessThan(testCase, norm(h1Stiffness - h1Stiffness', "fro"), 1e-12);

h1Mass = radia.ngsolveMatrix(meshPath, "h1", 2, "mass");
verifyEqual(testCase, size(h1Mass), [10, 10]);
verifyLessThan(testCase, norm(h1Mass - h1Mass', "fro"), 1e-12);
end

function testPersistentNGSolveMeshSpaceFormMatrix(testCase)
before = radia.apiInfo();
mesh = radia.ngsolve.Mesh.create(testCase.TestData.meshPath);
space = radia.ngsolve.FESpace.create(mesh, "hcurl", 2);
form = radia.ngsolve.BilinearForm.create(space, "mass");
matrix = form.matrix();

meshInfo = mesh.info();
spaceInfo = space.info();
formInfo = form.info();
matrixInfo = matrix.info();
verifyEqual(testCase, meshInfo.dimension, 3);
verifyEqual(testCase, meshInfo.vertices, 4);
verifyEqual(testCase, meshInfo.elements, 1);
verifyEqual(testCase, spaceInfo.dof_count, 14);
verifyEqual(testCase, space.DofCount, 14);
verifyEqual(testCase, formInfo.rows, 14);
verifyEqual(testCase, matrixInfo.rows, 14);
verifyEqual(testCase, matrixInfo.cols, 14);
verifyTrue(testCase, matrixInfo.is_sparse);
verifyTrue(testCase, matrix.IsSparse);

triplets = matrix.values();
reference = radia.ngsolve.matrix_dump( ...
    testCase.TestData.meshPath, "hcurl", 2, "mass");
nativeSparse = matrix.sparse();
verifyEqual(testCase, nativeSparse, reference, "AbsTol", 1e-12);
verifyEqual(testCase, triplets.shape, [14, 14]);
verifyEqual(testCase, triplets.nonzero_nze, nnz(reference));
verifyEqual(testCase, sparse(triplets.row, triplets.col, triplets.values, ...
    14, 14), reference, "AbsTol", 1e-12);

input = matrix.vector();
inputValues = (1:14).';
input.setValues(inputValues);
output = matrix.matvec(input);
verifyEqual(testCase, output.values(), reference * inputValues, ...
    "AbsTol", 1e-11);

inverse = matrix.inverse();
solution = inverse.matvec(input);
verifyEqual(testCase, solution.values(), reference \ inputValues, ...
    "RelTol", 1e-9, "AbsTol", 1e-11);

delete(solution);
delete(inverse);
delete(output);
delete(input);
delete(matrix);
delete(form);
delete(space);
delete(mesh);
after = radia.apiInfo();
verifyEqual(testCase, after.handle_count, before.handle_count);
end

function testChargeGramDirectNGSolveMassRegistration(testCase)
mesh = radia.ngsolve.Mesh.create(testCase.TestData.meshPath);
meshCleanup = onCleanup(@() delete(mesh));
space = radia.ngsolve.FESpace.create(mesh, "hcurl", 2);
spaceCleanup = onCleanup(@() delete(space));
form = radia.ngsolve.BilinearForm.create(space, "mass");
formCleanup = onCleanup(@() delete(form));
matrix = form.matrix();
matrixCleanup = onCleanup(@() delete(matrix));

n = matrix.Rows;
manager = radia.HACApKChargeGram([(0:n-1).', zeros(n,2)], ones(n,1), zeros(n,1));
managerCleanup = onCleanup(@() delete(manager));
manager.configureGeometryMassMatrixNGSolve(matrix);
manager.configureMassMatrixNGSolve(matrix);

x = sin((1:n).');
nativeInput = matrix.vector();
inputCleanup = onCleanup(@() delete(nativeInput));
nativeInput.setValues(x);
nativeOutput = matrix.matvec(nativeInput);
outputCleanup = onCleanup(@() delete(nativeOutput));
expected = nativeOutput.values();
clear outputCleanup inputCleanup matrixCleanup formCleanup spaceCleanup meshCleanup

verifyEqual(testCase, manager.applyConfiguredGeometryMass(x), expected, ...
    "AbsTol", 2e-12);
verifyEqual(testCase, manager.applyConfiguredMassRiesz(expected), x, ...
    "RelTol", 2e-10, "AbsTol", 2e-11);
info = manager.operatorInfo();
verifyTrue(testCase, info.mass_matrix_configured);
verifyTrue(testCase, info.geometry_mass_configured);
verifyEqual(testCase, info.n_face, n);
clear managerCleanup
end

function testPersistentNGSolveCoefficientBilinearForms(testCase)
before = radia.apiInfo();
mesh = radia.ngsolve.Mesh.create(testCase.TestData.meshPath);
coefficient = radia.ngsolve.CoefficientFunction.constant(2.5);

h1 = radia.ngsolve.FESpace.create(mesh, "h1", 2);
hcurl = radia.ngsolve.FESpace.create(mesh, "hcurl", 1);
hdiv = radia.ngsolve.FESpace.create(mesh, "hdiv", 1);

h1Base = radia.ngsolve.BilinearForm.create(h1, "mass");
h1Weighted = radia.ngsolve.BilinearForm.createFromCoefficient( ...
    h1, "mass", coefficient);
hcurlBase = radia.ngsolve.BilinearForm.create(hcurl, "mass");
hcurlWeighted = radia.ngsolve.BilinearForm.createFromCoefficient( ...
    hcurl, "mass", coefficient);
hdivBase = radia.ngsolve.BilinearForm.create(hdiv, "mass");
hdivWeighted = radia.ngsolve.BilinearForm.createFromCoefficient( ...
    hdiv, "mass", coefficient);
h1Boundary = radia.ngsolve.BilinearForm.createBoundaryFromCoefficient( ...
    h1, coefficient, Label="h1_boundary_mass");
hcurlBoundary = radia.ngsolve.BilinearForm.createBoundaryFromCoefficient( ...
    hcurl, coefficient, Label="hcurl_boundary_mass");
hdivBoundary = radia.ngsolve.BilinearForm.createBoundaryFromCoefficient( ...
    hdiv, coefficient, Label="hdiv_boundary_mass");

h1BaseMatrix = h1Base.matrix();
h1WeightedMatrix = h1Weighted.matrix();
hcurlBaseMatrix = hcurlBase.matrix();
hcurlWeightedMatrix = hcurlWeighted.matrix();
hdivBaseMatrix = hdivBase.matrix();
hdivWeightedMatrix = hdivWeighted.matrix();
h1BoundaryMatrix = h1Boundary.matrix();
hcurlBoundaryMatrix = hcurlBoundary.matrix();
hdivBoundaryMatrix = hdivBoundary.matrix();
verifyEqual(testCase, string(h1Weighted.info().form), "mass");
verifyEqual(testCase, h1WeightedMatrix.sparse(), ...
    2.5 * h1BaseMatrix.sparse(), "AbsTol", 1e-12);
verifyEqual(testCase, hcurlWeightedMatrix.sparse(), ...
    2.5 * hcurlBaseMatrix.sparse(), "AbsTol", 1e-12);
verifyEqual(testCase, hdivWeightedMatrix.sparse(), ...
    2.5 * hdivBaseMatrix.sparse(), "AbsTol", 1e-12);
verifyEqual(testCase, string(h1Boundary.info().form), "boundary_coefficient");
verifyGreaterThan(testCase, norm(h1BoundaryMatrix.sparse(), "fro"), 0);
verifyGreaterThan(testCase, norm(hcurlBoundaryMatrix.sparse(), "fro"), 0);
verifyGreaterThan(testCase, norm(hdivBoundaryMatrix.sparse(), "fro"), 0);

delete(h1WeightedMatrix);
delete(h1BaseMatrix);
delete(hcurlWeightedMatrix);
delete(hcurlBaseMatrix);
delete(hdivWeightedMatrix);
delete(hdivBaseMatrix);
delete(h1BoundaryMatrix);
delete(hcurlBoundaryMatrix);
delete(hdivBoundaryMatrix);
delete(h1Weighted);
delete(h1Base);
delete(hcurlWeighted);
delete(hcurlBase);
delete(hdivWeighted);
delete(hdivBase);
delete(h1Boundary);
delete(hcurlBoundary);
delete(hdivBoundary);
delete(h1);
delete(hcurl);
delete(hdiv);
delete(coefficient);
delete(mesh);
after = radia.apiInfo();
verifyEqual(testCase, after.handle_count, before.handle_count);
end

function testPersistentNGSolveComplexBilinearForms(testCase)
before = radia.apiInfo();
mesh = radia.ngsolve.Mesh.create(testCase.TestData.meshPath);
coefficient = radia.ngsolve.CoefficientFunction.constant(2.5 + 0.5i);
space = radia.ngsolve.FESpace.create(mesh, "h1", 2, Complex=true);
base = radia.ngsolve.BilinearForm.create(space, "mass");
weighted = radia.ngsolve.BilinearForm.createFromCoefficient( ...
    space, "mass", coefficient);
boundary = radia.ngsolve.BilinearForm.createBoundaryFromCoefficient( ...
    space, coefficient, Label="complex_boundary_mass");
baseMatrix = base.matrix();
weightedMatrix = weighted.matrix();
boundaryMatrix = boundary.matrix();

verifyTrue(testCase, base.info().is_complex);
verifyTrue(testCase, baseMatrix.IsComplex);
verifyTrue(testCase, weightedMatrix.IsComplex);
verifyTrue(testCase, boundaryMatrix.IsComplex);
baseTriplets = baseMatrix.values();
weightedTriplets = weightedMatrix.values();
verifyTrue(testCase, ~isreal(baseTriplets.values));
verifyEqual(testCase, weightedTriplets.values, ...
    (2.5 + 0.5i) * baseTriplets.values, "AbsTol", 1e-12);
verifyGreaterThan(testCase, norm(boundaryMatrix.sparse(), "fro"), 0);

A = baseMatrix.sparse();
x = baseMatrix.vector();
xValues = (1:baseMatrix.Cols).' + 1i * (baseMatrix.Cols:-1:1).';
x.setValues(xValues);
y = baseMatrix.matvec(x);
verifyEqual(testCase, y.values(), A * xValues, "AbsTol", 1e-11);

complexSolver = radia.ngsolve.Solver.create( ...
    baseMatrix, "gmres", Tolerance=1e-11, MaxSteps=100);
complexRhs = baseMatrix.vector();
complexRhs.setValues(xValues);
complexSolution = complexSolver.solve(complexRhs);
complexCheck = baseMatrix.matvec(complexSolution);
verifyEqual(testCase, complexCheck.values(), xValues, ...
    "RelTol", 1e-8, "AbsTol", 1e-10);

constantForm = radia.ngsolve.LinearForm.create( ...
    space, "constant", Value=1.25 + 0.75i);
coefficientForm = radia.ngsolve.LinearForm.createFromCoefficient( ...
    space, coefficient);
constantRhs = constantForm.vector();
coefficientRhs = coefficientForm.vector();
verifyTrue(testCase, constantForm.info().is_complex);
verifyTrue(testCase, coefficientForm.info().is_complex);
verifyTrue(testCase, ~isreal(constantRhs.values()));
verifyTrue(testCase, ~isreal(coefficientRhs.values()));
verifyGreaterThan(testCase, norm(constantRhs.values()), 0);
verifyGreaterThan(testCase, norm(coefficientRhs.values()), 0);

delete(coefficientRhs);
delete(constantRhs);
delete(coefficientForm);
delete(constantForm);
delete(y);
delete(x);
delete(complexCheck);
delete(complexSolution);
delete(complexRhs);
delete(complexSolver);
delete(weightedMatrix);
delete(baseMatrix);
delete(boundaryMatrix);
delete(weighted);
delete(base);
delete(boundary);
delete(space);
delete(coefficient);
delete(mesh);
after = radia.apiInfo();
verifyEqual(testCase, after.handle_count, before.handle_count);
end

function testPersistentNGSolveNativeKrylovSolver(testCase)
before = radia.apiInfo();
mesh = radia.ngsolve.Mesh.create(testCase.TestData.meshPath);
space = radia.ngsolve.FESpace.create(mesh, "h1", 2);
form = radia.ngsolve.BilinearForm.create(space, "mass");
matrix = form.matrix();
rhs = matrix.vector();
rhsValues = (1:matrix.Rows).';
rhs.setValues(rhsValues);
solver = radia.ngsolve.Solver.create( ...
    matrix, "cg", Tolerance=1e-11, MaxSteps=100);
solverInfo = solver.info();
verifyEqual(testCase, string(solverInfo.method), "cg");
verifyFalse(testCase, solverInfo.is_complex);
solution = solver.solve(rhs);
check = matrix.matvec(solution);
verifyEqual(testCase, check.values(), rhsValues, ...
    "RelTol", 1e-8, "AbsTol", 1e-10);
verifyGreaterThan(testCase, solver.info().steps, 0);

preconditioner = matrix.inverse();
preconditionedSolver = radia.ngsolve.Solver.create( ...
    matrix, "cg", Tolerance=1e-11, MaxSteps=20, ...
    Preconditioner=preconditioner);
preconditionedInfo = preconditionedSolver.info();
verifyTrue(testCase, preconditionedInfo.has_preconditioner);
preconditionedSolution = preconditionedSolver.solve(rhs);
preconditionedCheck = matrix.matvec(preconditionedSolution);
verifyEqual(testCase, preconditionedCheck.values(), rhsValues, ...
    "RelTol", 1e-8, "AbsTol", 1e-10);

delete(preconditionedCheck);
delete(preconditionedSolution);
delete(preconditionedSolver);
delete(preconditioner);
delete(check);
delete(solution);
delete(solver);
delete(rhs);
delete(matrix);
delete(form);
delete(space);
delete(mesh);
after = radia.apiInfo();
verifyEqual(testCase, after.handle_count, before.handle_count);
end

function testPersistentNGSolveGridFunctionAndLinearForm(testCase)
before = radia.apiInfo();
mesh = radia.ngsolve.Mesh.create(testCase.TestData.meshPath);
space = radia.ngsolve.FESpace.create(mesh, "h1", 2);
gridFunction = radia.ngsolve.GridFunction.fromFESpace( ...
    space, Name="temperature");
gridInfo = gridFunction.info();
verifyEqual(testCase, string(gridInfo.name), "temperature");
verifyEqual(testCase, gridInfo.dof_count, space.DofCount);

linearForm = radia.ngsolve.LinearForm.create( ...
    space, "constant", Value=2.0, Label="rhs");
linearInfo = linearForm.info();
verifyEqual(testCase, string(linearInfo.source), "constant");
verifyEqual(testCase, string(linearInfo.label), "rhs");
verifyEqual(testCase, linearInfo.size, space.DofCount);
rhs = linearForm.vector();
verifyTrue(testCase, rhs.IsView);
rhsValues = rhs.values();
verifySize(testCase, rhsValues, [space.DofCount, 1]);
verifyTrue(testCase, all(isfinite(rhsValues)));
verifyGreaterThan(testCase, norm(rhsValues), 0);

gridValues = (1:space.DofCount).';
gridFunction.setVector(gridValues);
verifyEqual(testCase, gridFunction.vector(), gridValues);
delete(linearForm);
verifyEqual(testCase, rhs.values(), rhsValues, "AbsTol", 1e-12);

delete(rhs);
delete(gridFunction);
delete(space);
delete(mesh);
after = radia.apiInfo();
verifyEqual(testCase, after.handle_count, before.handle_count);
end

function testNGSolveCoefficientAndGridFunctionObjects(testCase)
meshPath = testCase.TestData.meshPath;
before = radia.apiInfo();

scalar = radia.ngsolve.CoefficientFunction.constant(2.5);
scalarInfo = scalar.info();
verifyEqual(testCase, scalarInfo.dimension, 1);
verifyEmpty(testCase, scalarInfo.dimensions);
verifyFalse(testCase, scalarInfo.is_complex);
samplePoints = [0.1, 0.2, 0.3; 0.2, 0.1, 0.1];
verifyEqual(testCase, scalar.evaluate(meshPath, samplePoints), ...
    2.5 * ones(2, 1), "AbsTol", 1e-12);

vector = radia.ngsolve.CoefficientFunction.constant([1; 2; 3]);
vectorInfo = vector.info();
verifyEqual(testCase, vectorInfo.dimension, 3);
verifyEqual(testCase, vectorInfo.dimensions, 3);
verifyFalse(testCase, vectorInfo.is_complex);

sumCoefficient = scalar + scalar;
scaledCoefficient = 4.0 * scalar;
verifyEqual(testCase, sumCoefficient.info().dimension, 1);
verifyEqual(testCase, scaledCoefficient.info().dimension, 1);

gridFunction = radia.ngsolve.GridFunction.create(meshPath, "h1", 2);
gridInfo = gridFunction.info();
verifyEqual(testCase, gridInfo.dof_count, 10);
verifyEqual(testCase, gridFunction.DofCount, 10);
verifyEqual(testCase, string(gridFunction.Space), "h1");
verifySize(testCase, gridFunction.vector(), [10, 1]);

values = (1:10).';
gridFunction.setVector(values);
verifyEqual(testCase, gridFunction.vector(), values);

vectorView = gridFunction.vectorHandle();
vectorInfo = vectorView.info();
verifyEqual(testCase, vectorInfo.size, 10);
verifyTrue(testCase, vectorInfo.is_view);
verifyEqual(testCase, vectorInfo.component, 1);
verifyEqual(testCase, vectorView.values(), values);

workVector = vectorView.copy();
verifyFalse(testCase, workVector.IsView);
workVector.setZero();
verifyEqual(testCase, workVector.norm(), 0, "AbsTol", 1e-14);
workVector.axpy(2.0, vectorView);
verifyEqual(testCase, workVector.values(), 2.0 * values, "AbsTol", 1e-12);
verifyEqual(testCase, workVector.norm(), 2.0 * norm(values), "AbsTol", 1e-12);
verifyEqual(testCase, workVector.dot(vectorView), 2.0 * sum(values.^2), ...
    "AbsTol", 1e-12);
vectorView.scale(0.5);
verifyEqual(testCase, gridFunction.vector(), 0.5 * values, "AbsTol", 1e-12);

gridFunction.interpolate(scalar);
% NGSolve's H1 interpolation assigns the nodal DoFs; the remaining
% high-order DoFs are intentionally left at their previous zero state.
interpolated = gridFunction.vector();
verifyEqual(testCase, interpolated(1:4), 2.5 * ones(4, 1), ...
    "AbsTol", 1e-12);
verifyEqual(testCase, interpolated(5:10), zeros(6, 1), ...
    "AbsTol", 1e-12);

complexGridFunction = radia.ngsolve.GridFunction.create( ...
    meshPath, "h1", 1, Complex=true);
complexInfo = complexGridFunction.info();
verifyTrue(testCase, complexInfo.is_complex);
complexValues = (1:complexInfo.dof_count).' + ...
    1i * (complexInfo.dof_count:-1:1).';
complexGridFunction.setVector(complexValues);
verifyEqual(testCase, complexGridFunction.vector(), complexValues);
complexVector = complexGridFunction.vectorHandle();
complexWork = complexVector.copy();
verifyEqual(testCase, complexWork.values(), complexValues);
verifyEqual(testCase, complexWork.norm(), norm(complexValues), "AbsTol", 1e-12);
verifyEqual(testCase, complexWork.dot(complexWork), ...
    sum(complexValues .* complexValues), "AbsTol", 1e-12);
verifyEqual(testCase, complexWork.dot(complexWork, Conjugate=true), ...
    sum(abs(complexValues).^2), "AbsTol", 1e-12);

gridCoefficient = gridFunction.asCoefficient();
gridCoefficientInfo = gridCoefficient.info();
verifyEqual(testCase, gridCoefficientInfo.dimension, 1);
verifyFalse(testCase, gridCoefficientInfo.is_complex);

delete(gridCoefficient);
delete(workVector);
delete(vectorView);
delete(complexWork);
delete(complexVector);
delete(complexGridFunction);
retainedView = gridFunction.vectorHandle();
delete(gridFunction);
verifyEqual(testCase, retainedView.values(), interpolated, "AbsTol", 1e-12);
delete(retainedView);
delete(scaledCoefficient);
delete(sumCoefficient);
delete(vector);
delete(scalar);
after = radia.apiInfo();
verifyEqual(testCase, after.handle_count, before.handle_count);
end

function testNativeHCurlEddyCLNBasis(testCase)
meshPath = testCase.TestData.meshPath;
[~, parent] = radia.ngsolve.matrix_dump(meshPath, "hcurl", 2, "mass");
ports = zeros(parent.dof_count, 3);
ports(1, 1) = 1;
ports(2, 2) = 1;

result = radia.ngsolve.hcurl_eddy_cln_native_basis( ...
    meshPath, 2, ports, 2);
verifyEqual(testCase, result.dof_count, parent.dof_count);
verifySize(testCase, result.free_dofs, [1, parent.dof_count]);
verifyEqual(testCase, size(result.vectors, 1), parent.dof_count);
verifyEqual(testCase, result.rank, size(result.vectors, 2));
verifyGreaterThan(testCase, result.rank, 0);
verifyLessThan(testCase, result.orthogonality_error, 1e-10);
verifyEqual(testCase, string(result.operator), "mass+curlcurl");
verifyTrue(testCase, all(isfinite(result.vectors), "all"));
end

function testNativeHCurlEddyCLNModelProjection(testCase)
meshPath = testCase.TestData.meshPath;
order = 6;
[parentMass, parentMassInfo] = radia.ngsolveMatrix( ...
    meshPath, "hcurl", order, "mass");
parentCurlCurl = radia.ngsolveMatrix(meshPath, "hcurl", order, "stiffness");
ports = zeros(parentMassInfo.dof_count, 2);
ports(1, 1) = 1;
ports(2, 2) = 1;

model = radia.ngsolve.hcurl_eddy_cln_model( ...
    meshPath, order, ports, 2, Conductivity=2.5, Reluctivity=3.0);
basis = model.native_basis;
V = basis.vectors;

verifyEqual(testCase, string(model.assembly_schema), ...
    "radia.hcurl.eddy_cln.native_diffusion.v1");
verifyEqual(testCase, model.state_order, basis.rank);
verifyEqual(testCase, model.port_count, size(ports, 2));
expectedMassGram = V' * parentMass * V;
expectedCurlCurlGram = V' * parentCurlCurl * V;
expectedPortRHS = V' * ports;
verifyLessThan(testCase, norm(basis.mass_gram - expectedMassGram, "fro"), 1e-10);
verifyLessThan(testCase, norm(basis.curlcurl_gram - expectedCurlCurlGram, "fro"), 1e-10);
verifyLessThan(testCase, norm(basis.port_rhs - expectedPortRHS, "fro"), 1e-12);
verifyLessThan(testCase, norm(model.resistance - 3.0 * basis.curlcurl_gram, "fro"), 1e-12);
verifyLessThan(testCase, norm(model.inductance - 2.5 * basis.mass_gram, "fro"), 1e-12);
verifyTrue(testCase, model.passive);
verifyTrue(testCase, all(isfinite(model.Ad), "all"));
verifyTrue(testCase, all(isfinite(model.Bd), "all"));
end

function testNativeSimulinkStateSpaceHandle(testCase)
before = radia.apiInfo();
A = [1, 0.1; 0, 1];
B = [0; 0.1];
C = [1, 0];
D = 0;
x0 = [2; 0];
handle = radia.internal.callMex( ...
    'simulink.state_space.create', A, B, C, D, x0);
cleanup = onCleanup(@() destroyNativeStateSpace(handle));
info = radia.internal.callMex('simulink.state_space.info', handle);
verifyEqual(testCase, info.state_size, 2);
verifyEqual(testCase, info.input_size, 1);
verifyEqual(testCase, info.output_size, 1);
verifyEqual(testCase, radia.internal.callMex( ...
    'simulink.state_space.output', handle, 3), 2, "AbsTol", 1e-12);
verifyEqual(testCase, radia.internal.callMex( ...
    'simulink.state_space.output', handle, 3), 2, "AbsTol", 1e-12);
info = radia.internal.callMex('simulink.state_space.info', handle);
verifyEqual(testCase, info.step_count, 0);
radia.internal.callMex('simulink.state_space.update', handle, 3);
snapshot = radia.internal.callMex('simulink.state_space.snapshot', handle);
verifyEqual(testCase, snapshot.state, [2; 0.3], "AbsTol", 1e-12);
verifyEqual(testCase, snapshot.step_count, 1);
radia.internal.callMex('simulink.state_space.update', handle, 0);
radia.internal.callMex('simulink.state_space.restore', handle, snapshot);
restored = radia.internal.callMex('simulink.state_space.snapshot', handle);
verifyEqual(testCase, restored.state, snapshot.state, "AbsTol", 0);
verifyEqual(testCase, restored.step_count, snapshot.step_count, "AbsTol", 0);
radia.internal.callMex('simulink.state_space.reset', handle);
verifyEqual(testCase, radia.internal.callMex( ...
    'simulink.state_space.step', handle, 0), 2, "AbsTol", 1e-12);
clear cleanup
destroyNativeStateSpace(handle);
after = radia.apiInfo();
verifyEqual(testCase, after.handle_count, before.handle_count);
end

function destroyNativeStateSpace(handle)
try
    radia.internal.callMex('simulink.state_space.destroy', handle);
catch
end
end

function testNativePeriodicMotorAngleFamilyHandle(testCase)
before = radia.apiInfo();
grid = [0; pi];
period = 2*pi;
A = reshape([1.0, 0.5], 1, 1, 2);
B = reshape([1.0, 3.0], 1, 1, 2);
C = reshape([1.0, 3.0], 1, 1, 2);
D = zeros(1, 1, 2);
Q = reshape([2.0, 4.0], 1, 1, 2);
R = reshape([1.0, 2.0], 1, 1, 2);
S = reshape([0.0, 2.0], 1, 1, 2);
handle = radia.internal.callMex( ...
    'simulink.state_space.create', grid, period, A, B, C, D, Q, R, S, 2.0);
cleanup = onCleanup(@() destroyNativeStateSpace(handle));

info = radia.internal.callMex('simulink.state_space.info', handle);
verifyEqual(testCase, string(info.model_kind), "periodic_angle_family");
verifyEqual(testCase, info.snapshot_count, 2);
verifyEqual(testCase, info.state_size, 1);
verifyEqual(testCase, info.input_size, 1);
verifyEqual(testCase, info.linear_output_size, 1);
verifyEqual(testCase, info.output_size, 2);
verifyEqual(testCase, info.period, period, "AbsTol", 0);

first = radia.internal.callMex( ...
    'simulink.state_space.output', handle, pi/2, 3.0);
verifyEqual(testCase, first, [4.0; 19.5], "AbsTol", 1e-12);
repeated = radia.internal.callMex( ...
    'simulink.state_space.output', handle, pi/2, 3.0);
verifyEqual(testCase, repeated, first, "AbsTol", 0);
info = radia.internal.callMex('simulink.state_space.info', handle);
verifyEqual(testCase, info.step_count, 0);

radia.internal.callMex( ...
    'simulink.state_space.update', handle, pi/2, 3.0);
snapshot = radia.internal.callMex('simulink.state_space.snapshot', handle);
verifyEqual(testCase, snapshot.state, 7.5, "AbsTol", 1e-12);
verifyEqual(testCase, snapshot.step_count, 1);
second = radia.internal.callMex( ...
    'simulink.state_space.output', handle, pi/2, 3.0);
verifyEqual(testCase, second, [15.0; 122.625], "AbsTol", 1e-12);
info = radia.internal.callMex('simulink.state_space.info', handle);
verifyEqual(testCase, info.step_count, 1);
verifyEqual(testCase, info.last_coordinate, pi/2, "AbsTol", 1e-12);

radia.internal.callMex( ...
    'simulink.state_space.update', handle, pi/2, 3.0);
radia.internal.callMex('simulink.state_space.restore', handle, snapshot);
verifyEqual(testCase, radia.internal.callMex( ...
    'simulink.state_space.output', handle, pi/2, 3.0), ...
    second, "AbsTol", 0);

radia.internal.callMex('simulink.state_space.reset', handle);
wrapped = radia.internal.callMex( ...
    'simulink.state_space.step', handle, 2*pi + pi/2, 3.0);
verifyEqual(testCase, wrapped, first, "AbsTol", 1e-12);
clear cleanup
destroyNativeStateSpace(handle);
after = radia.apiInfo();
verifyEqual(testCase, after.handle_count, before.handle_count);
end

function testComplexMixedGalerkinKernels(testCase)
Kkk = [4 + 1i, 0.3 - 0.2i; 0.3 + 0.2i, 3 - 0.5i];
Kke = [0.2 + 0.1i; 0.4 - 0.3i];
Kek = [0.5 - 0.2i, 0.1 + 0.4i];
Kee = 2.5 + 0.7i;
actual = radia.schurComplement(Kkk, Kke, Kek, Kee);
expected = Kkk - Kke * (Kee \ Kek);
verifyEqual(testCase, actual, expected, "AbsTol", 1e-14);

A = [3 + 1i, 0.2; -0.1i, 2 - 0.4i];
b = [1 + 0.5i; -0.2 + 1i];
verifyEqual(testCase, radia.denseSolve(A, b), A \ b, "AbsTol", 1e-14);
end

function testCLNReductionKernels(testCase)
K = [2, 0.1; 0.1, 1];
N = [3, 0.2; 0.2, 1.5];
reduced = radia.clnLanczos(K, N, 2, 1e-30);
verifyEqual(testCase, reduced.n_input, 2);
verifyEqual(testCase, reduced.n_output, 2);
verifySize(testCase, reduced.Q, [2, 2]);
verifySize(testCase, reduced.R_diag, [2, 2]);
verifySize(testCase, reduced.L_tridiag, [2, 2]);
verifyTrue(testCase, all(isfinite([reduced.Q, reduced.R_diag, reduced.L_tridiag]), "all"));

diagValues = [4, 2, 1];
T = radia.clnBuildTridiagonal(diagValues);
verifyEqual(testCase, diag(T), [6; 3; 1], "AbsTol", 1e-14);
verifyEqual(testCase, T(1, 2), -diagValues(2), "AbsTol", 1e-14);
verifyEqual(testCase, T(2, 3), -diagValues(3), "AbsTol", 1e-14);

frequency = 1000;
s = 1i * 2 * pi * frequency;
rhs = [1; 0];
solution = (reduced.R_diag + s * reduced.L_tridiag) \ rhs;
expectedZ = 1 / solution(1);
verifyEqual(testCase, radia.clnImpedance( ...
    reduced.R_diag, reduced.L_tridiag, frequency), expectedZ, "RelTol", 1e-12);
verifySize(testCase, radia.clnImpedanceSweep( ...
    reduced.R_diag, reduced.L_tridiag, [0, frequency]), [2, 1]);

Q = [1, 0; 0, 1; 1, 1];
M = [1, 2; 3, 4; 5, 6];
verifyEqual(testCase, radia.clnTransformCoupling(Q, M), Q' * M, "AbsTol", 1e-14);
verifyEqual(testCase, radia.clnTransformPort(Q, [2; 3; 4]), Q' * [2; 3; 4], ...
    "AbsTol", 1e-14);

P = [1, 2, 3; 2, 4, 6; 3, 6, 9];
aca = radia.clnAcaCompress(P, 1e-10, 3);
verifyEqual(testCase, aca.n, 3);
verifyGreaterThanOrEqual(testCase, aca.k, 1);
verifyLessThanOrEqual(testCase, aca.k, 3);
verifyTrue(testCase, aca.converged);
verifyEqual(testCase, aca.U * aca.V', P, "RelTol", 1e-8);
end

function testEVRSTMethodAlgebra(testCase)
C = [1, 0; 0, 1; 1, -1];
D = zeros(1, 3);
G = zeros(2, 1);
Q = [1; 0];
MR = [4, 0.2, 0.1; 0.2, 3, 0.4; 0.1, 0.4, 2];
ML = [2, 0.1, 0; 0.1, 5, 0.3; 0, 0.3, 4];
P = [1; 2; -1];

result = radia.evrsTMethod(C, D, G, Q, MR, ML, P);
CQ = C * Q;
verifyEqual(testCase, result.current_evrs, CQ, "AbsTol", 1e-14);
verifyEqual(testCase, result.resistance_t, C' * MR * C, "AbsTol", 1e-14);
verifyEqual(testCase, result.inductance_t, C' * ML * C, "AbsTol", 1e-14);
verifyEqual(testCase, result.resistance_evrs, CQ' * MR * CQ, "AbsTol", 1e-14);
verifyEqual(testCase, result.port_evrs, Q' * C' * P, "AbsTol", 1e-14);
verifyEqual(testCase, result.diagnostics.div_curl_norm, 0, "AbsTol", 1e-14);
verifyEqual(testCase, result.diagnostics.evrs_resistance_galerkin_residual, ...
    0, "AbsTol", 1e-14);
end

function testSIBCKernels(testCase)
s = 1i * 2 * pi * 750;
sigma = 5.8e7;
mu = 4 * pi * 1e-7;
surfaceMeasure = 0.03;
kSibc = 2.4;
d = 0.2;
verifyEqual(testCase, radia.skinImpedance(s, sigma, mu), ...
    sqrt(mu * s / sigma), "AbsTol", 1e-18);
verifyEqual(testCase, radia.sibcAdmittanceTail(s, surfaceMeasure, sigma, mu), ...
    surfaceMeasure * sqrt(sigma / (mu * s)), "RelTol", 1e-14);
z = (s + d) / (kSibc * sqrt(s));
verifyEqual(testCase, radia.sibcTerminationImpedance(s, kSibc, d), ...
    z, "RelTol", 1e-14);
verifyEqual(testCase, radia.sibcTerminationAdmittance(s, kSibc, d), ...
    1 / z, "RelTol", 1e-14);
end

function testTetHCurlReducedGram(testCase)
cellVerts = [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1];
exponents = [0, 0, 0];
coefficients = [1, 0, 0];
refPoints = [0.25, 0.25, 0.25];
refWeights = 1;

gram = radia.tetHCurlReducedGram(cellVerts, exponents, coefficients, ...
    1, refPoints, refWeights);
canonicalGram = radia.hcurl.tet_reduced_gram(cellVerts, exponents, ...
    coefficients, 1, refPoints, refWeights);
verifySize(testCase, gram, [1, 1]);
verifyEqual(testCase, canonicalGram, gram, "AbsTol", 1e-14);
verifyTrue(testCase, isfinite(gram));
verifyGreaterThan(testCase, gram, 0);
end

function testComplexBiotSavartKernels(testCase)
segments = zeros(1, 2, 3);
segments(1, 2, 3) = 1;
obs = [0.3, 0.4, 0.5; -0.2, 0.1, 0.7];
currentRe = 1;
currentIm = 2;

[hRe, hIm] = radia.hFromSegmentsComplex(segments, obs, currentRe, currentIm);
[aRe, aIm] = radia.aFromSegmentsComplex(segments, obs, currentRe, currentIm);
verifySize(testCase, hRe, [2, 3]);
verifySize(testCase, aRe, [2, 3]);
verifyEqual(testCase, hIm, 2 * hRe, "RelTol", 1e-13);
verifyEqual(testCase, aIm, 2 * aRe, "RelTol", 1e-13);
verifyTrue(testCase, all(isfinite([hRe, hIm, aRe, aIm]), "all"));

vertices = reshape([0, 0, 0; 1, 0, 0; 0, 1, 0], 1, 3, 3);
surfaceCurrentRe = [0, 0, 1];
surfaceCurrentIm = [0, 1, 0];
[bRe, bIm] = radia.bFromTrianglesComplex( ...
    vertices, surfaceCurrentRe, surfaceCurrentIm, [0.2, 0.3, 1.0]);
[surfaceARe, surfaceAIm] = radia.aFromTrianglesComplex( ...
    vertices, surfaceCurrentRe, surfaceCurrentIm, [0.2, 0.3, 1.0]);
verifySize(testCase, bRe, [1, 3]);
verifySize(testCase, surfaceARe, [1, 3]);
verifyTrue(testCase, all(isfinite([bRe, bIm, surfaceARe, surfaceAIm]), "all"));
end

function testBEMGalerkinKernels(testCase)
vertices = [0, 0, 0; 1, 0, 0; 0, 1, 0];
triangles = int64([0, 1, 2]);
p2Nodes = [vertices(1, :), vertices(2, :), vertices(3, :), ...
    0.5 * (vertices(1, :) + vertices(2, :)), ...
    0.5 * (vertices(2, :) + vertices(3, :)), ...
    0.5 * (vertices(3, :) + vertices(1, :))];

[SL, DL] = radia.assembleSldlGalerkin(vertices, triangles, p2Nodes, 3, 2, 1);
verifySize(testCase, SL, [3, 3]);
verifySize(testCase, DL, [3, 3]);
verifyTrue(testCase, all(isfinite([SL, DL]), "all"));
verifyGreaterThan(testCase, trace(SL), 0);

dofs = int64(0:5);
[SLp2, DLp2] = radia.assembleSldlGalerkinP2( ...
    vertices, triangles, p2Nodes, dofs, 6, 3, 2, 1);
verifySize(testCase, SLp2, [6, 6]);
verifySize(testCase, DLp2, [6, 6]);
verifyTrue(testCase, all(isfinite([SLp2, DLp2]), "all"));
verifyGreaterThan(testCase, trace(SLp2), 0);
end

function testHACApKBEMLifecycle(testCase)
coordinates = [0, 0, 0; 1, 0, 0; 0, 1, 0];
entries = [4, 0.2, 0.1; 0.2, 3, 0.4; 0.1, 0.4, 2];
before = radia.apiInfo();
manager = radia.HACApKBEMManager(coordinates, entries);
created = radia.apiInfo();
verifyEqual(testCase, created.handle_count, before.handle_count + 1);
verifyTrue(testCase, manager.build(AcaEps=1e-10, LeafSize=2, Eta=2, ...
    MaxRank=20, PrintLevel=0));
info = manager.info();
verifyTrue(testCase, info.valid);
verifyEqual(testCase, info.n_dof, 3);
x = [1; -2; 0.5];
verifyEqual(testCase, manager.matvec(x), entries * x, "AbsTol", 1e-10);
delete(manager);
after = radia.apiInfo();
verifyEqual(testCase, after.handle_count, before.handle_count);
end

function testHACApKPEECLifecycle(testCase)
centers = [0, 0, 0; 0.01, 0, 0; 0.02, 0, 0];
directions = repmat([0, 0, 1], 3, 1);
lengths = 0.1 * ones(3, 1);
widths = 1e-3 * ones(3, 1);
heights = 1e-3 * ones(3, 1);
sigmas = 5.8e7 * ones(3, 1);
before = radia.apiInfo();
manager = radia.HACApKPEECManager( ...
    centers, directions, lengths, widths, heights, sigmas);
verifyTrue(testCase, manager.build(AcaEps=1e-10, LeafSize=2, Eta=3, ...
    MaxRank=20, PrintLevel=0));
created = radia.apiInfo();
verifyEqual(testCase, created.handle_count, before.handle_count + 1);
info = manager.info();
verifyTrue(testCase, info.valid);
verifyEqual(testCase, info.n_dof, 3);
y = manager.matvec([1; -2; 0.5]);
verifySize(testCase, y, [3, 1]);
verifyTrue(testCase, all(isfinite(y)));
verifyGreaterThan(testCase, norm(y), 0);
delete(manager);
after = radia.apiInfo();
verifyEqual(testCase, after.handle_count, before.handle_count);
end

function testHACApKChargeGramLifecycle(testCase)
centroids = [0, 0, 0; 1, 0, 0; 0, 1, 0];
measures = [0.5; 0.5; 0.5];
selfEnergy = [0.2; 0.2; 0.2];
before = radia.apiInfo();
manager = radia.HACApKChargeGram(centroids, measures, selfEnergy);
verifyTrue(testCase, manager.build(AcaEps=1e-10, LeafSize=2, Eta=2, ...
    MaxRank=20, PrintLevel=0));
created = radia.apiInfo();
verifyEqual(testCase, created.handle_count, before.handle_count + 1);
info = manager.info();
verifyTrue(testCase, info.valid);
verifyEqual(testCase, info.n_dof, 3);
verifyEqual(testCase, manager.entry(1, 1), selfEnergy(1), "AbsTol", 1e-14);
x = [1; -2; 0.5];
y = manager.matvec(x);
verifySize(testCase, y, [3, 1]);
verifyTrue(testCase, all(isfinite(y)));
ySym = manager.matvecSym(x);
yTranspose = manager.matvecTranspose(x);
verifySize(testCase, ySym, [3, 1]);
verifySize(testCase, yTranspose, [3, 1]);
verifyTrue(testCase, all(isfinite([ySym; yTranspose])));

manager.configureChargeMap(int32([0; 1; 2; 3]), int32([0; 1; 0]), [1; 1; 1], 2);
manager.configureGeometryMassMatrix(int32([0; 1]), int32([0; 1]), [1; 1], 2);
manager.configureMassMatrix(int32([0; 1]), int32([0; 1]), [2; 3], 2);
operatorInfo = manager.operatorInfo();
verifyTrue(testCase, operatorInfo.operator_configured);
verifyEqual(testCase, manager.applyConfiguredGeometryMass([1; 2]), [1; 2], ...
    "AbsTol", 1e-14);
verifyEqual(testCase, manager.applyConfiguredMassRiesz([2; 3]), [2; 3], ...
    "AbsTol", 1e-14);
expectedDemag = manager.applyConfiguredDemag([1; 2]);
verifyTrue(testCase, all(isfinite(expectedDemag)));
invChi = 0.4;
expectedOperator = expectedDemag + invChi * [2; 6];
verifyEqual(testCase, ...
    manager.applyConfiguredLinearMaterialOperator(invChi, [1; 2]), ...
    expectedOperator, "AbsTol", 2e-12);
batchInput = [1, 2; -0.5, 0.25];
batchOutput = manager.applyConfiguredLinearMaterialOperatorMany( ...
    invChi, batchInput);
verifyEqual(testCase, batchOutput(1,:).', expectedOperator, "AbsTol", 2e-12);
verifyEqual(testCase, batchOutput(2,:).', ...
    manager.applyConfiguredLinearMaterialOperator(invChi, batchInput(2,:).'), ...
    "AbsTol", 2e-12);

manager.setConfiguredConstraints(int32(0));
verifyEqual(testCase, manager.operatorInfo().constraint_count, 1);
constrained = manager.applyConfiguredLinearMaterialOperator( ...
    invChi, [1; 2]);
verifyEqual(testCase, constrained(1), 0, "AbsTol", 1e-14);
verifyEqual(testCase, ...
    manager.applyConfiguredLinearMaterialOperatorMany(invChi, batchInput), ...
    [constrained.'; manager.applyConfiguredLinearMaterialOperator( ...
        invChi, batchInput(2,:).').'], "AbsTol", 2e-12);
manager.setConfiguredConstraints(int32.empty);

rhsBatch = [1, -0.25; 0.2, 0.75];
solvedBatch = manager.solveConfiguredLinearMaterialAutoPrecMany( ...
    invChi, rhsBatch, Tol=1e-11, MaxIt=5000);
verifySize(testCase, solvedBatch.m, size(rhsBatch));
verifySize(testCase, solvedBatch.iters, [size(rhsBatch,1), 1]);
for k = 1:size(rhsBatch,1)
    residual = manager.applyConfiguredLinearMaterialOperator( ...
        invChi, solvedBatch.m(k,:).') - rhsBatch(k,:).';
    verifyLessThan(testCase, norm(residual), 2e-9);
end

manager.setConfiguredConstraints(int32(0));
activeRhs = [0; 0.75];
activeState = manager.solveConfiguredLinearMaterialAutoPrecMany( ...
    invChi, activeRhs.', Tol=1e-11, MaxIt=5000).m.';
responseMatrix = [0, 1];
adjoints = manager.solveConfiguredLinearMaterialAutoPrecMany( ...
    invChi, responseMatrix, Tol=1e-11, MaxIt=5000).m;
reduced = manager.reduceConfiguredCandidateSchur( ...
    invChi, int32(0), activeRhs, activeState, responseMatrix, adjoints, ...
    Tol=1e-11, MaxIt=5000, SolveBatchSize=1, MassRiesz=false);
verifySize(testCase, reduced.schur, [1, 1]);
verifySize(testCase, reduced.rhs, [1, 1]);
verifySize(testCase, reduced.response, [1, 1]);
verifyTrue(testCase, all(isfinite([reduced.schur(:); reduced.rhs(:); ...
    reduced.response(:)])));
manager.setConfiguredConstraints(int32.empty);
demag = manager.demagMatrix();
nativeInput = demag.vector();
nativeInput.setValues([1; 2]);
nativeOutput = demag.matvec(nativeInput);
verifyEqual(testCase, nativeOutput.values(), expectedDemag, "AbsTol", 2e-13);
verifyTrue(testCase, isfield(manager.stats(), "n_dof"));
delete(manager);
retainedOutput = demag.matvec(nativeInput);
verifyEqual(testCase, retainedOutput.values(), expectedDemag, "AbsTol", 2e-13);
delete(retainedOutput); delete(nativeOutput); delete(nativeInput); delete(demag);
after = radia.apiInfo();
verifyEqual(testCase, after.handle_count, before.handle_count);
end

function testHACApKSampledChargeGramFactories(testCase)
points = [0, 0, 0; 1, 0, 0; 0, 2, 0];
weights = [0.5; 1.25; 2.0];
kernelEpsilon = 0.2;

laplace = radia.HACApKChargeGram.from_sampled_laplace( ...
    points, weights, kernelEpsilon, AcaEps=1e-12, LeafSize=4);
laplaceCleanup = onCleanup(@() delete(laplace));
expectedLaplace = zeros(3);
for i = 1:3
    for j = 1:3
        distance = sqrt(sum((points(i,:) - points(j,:)).^2) + kernelEpsilon^2);
        expectedLaplace(i,j) = weights(i) * weights(j) / (4 * pi * distance);
        verifyEqual(testCase, laplace.entry(i, j), expectedLaplace(i,j), ...
            "RelTol", 2e-14);
    end
end
x = [1; -0.5; 0.25];
verifyEqual(testCase, laplace.matvec(x), expectedLaplace * x, ...
    "AbsTol", 2e-12);

referenceLength = 3.5;
planar = radia.HACApKChargeGram.from_sampled_planar_log( ...
    points, weights, kernelEpsilon, referenceLength, Build=false);
planarCleanup = onCleanup(@() delete(planar));
verifyEqual(testCase, planar.NDOF, 3);
expectedPlanar = zeros(3);
for i = 1:3
    for j = 1:3
        distance = sqrt(sum((points(i,1:2) - points(j,1:2)).^2) + ...
            kernelEpsilon^2);
        expectedPlanar(i,j) = -weights(i) * weights(j) * ...
            log(distance / referenceLength) / (2 * pi);
        verifyEqual(testCase, planar.entry(i, j), expectedPlanar(i,j), ...
            "AbsTol", 2e-14);
    end
end
verifyTrue(testCase, planar.build(AcaEps=1e-12, LeafSize=4));
verifyEqual(testCase, planar.matvec(x), expectedPlanar * x, ...
    "AbsTol", 2e-12);
clear planarCleanup laplaceCleanup
end

function testHACApKLocalPolynomialFactory(testCase)
cellVerts = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];
manager = radia.HACApKChargeGram.from_local_polynomials( ...
    cellVerts, 1, int32(0), 1, [0, 0, 0], ...
    [0.25, 0.25, 0.25], 1/6, AcaEps=1e-12, LeafSize=2);
cleanup = onCleanup(@() delete(manager));
verifyEqual(testCase, manager.NDOF, 1);
verifyGreaterThan(testCase, manager.entry(1, 1), 0);
verifyTrue(testCase, isfinite(manager.matvec(1)));
clear cleanup
end

function testHACApKHighOrderTetMatchesLocalPolynomial(testCase)
cellVerts = [0, 0, 0; 0.7, 0.1, 0; 0.1, 0.8, 0.05; 0.05, 0.2, 0.9];
exponents = int32([0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1]);
a = 0.5854101966249685;
b = 0.1381966011250105;
tetPoints = [a, b, b; b, a, b; b, b, a; b, b, b];
tetWeights = ones(4, 1) / 24;
trianglePoints = [1/6, 1/6; 2/3, 1/6; 1/6, 2/3];
triangleWeights = ones(3, 1) / 6;

highOrder = radia.HACApKChargeGram.from_high_order_tet( ...
    cellVerts, zeros(0,3), 1, int32(zeros(4,1)), int32(zeros(4,1)), ...
    exponents, tetPoints, tetWeights, trianglePoints, triangleWeights, ...
    InnerTetPoints=tetPoints, InnerTetWeights=tetWeights, ...
    InnerTrianglePoints=trianglePoints, InnerTriangleWeights=triangleWeights, ...
    Build=false);
highOrderCleanup = onCleanup(@() delete(highOrder));
local = radia.HACApKChargeGram.from_local_polynomials( ...
    cellVerts, 1, int32(zeros(4,1)), eye(4), double(exponents), ...
    tetPoints, tetWeights, Build=false);
localCleanup = onCleanup(@() delete(local));
verifyEqual(testCase, highOrder.NDOF, 4);
highOrderGram = zeros(4);
localGram = zeros(4);
for i = 1:4
    for j = 1:4
        highOrderGram(i,j) = highOrder.entry(i, j);
        localGram(i,j) = local.entry(i, j);
    end
end
verifyEqual(testCase, highOrderGram, localGram, "RelTol", 3e-14, "AbsTol", 1e-18);
verifyEqual(testCase, highOrderGram, highOrderGram.', "RelTol", 2e-14);
verifyTrue(testCase, highOrder.build(AcaEps=1e-12, LeafSize=8));
x = [1; -0.5; 0.25; -0.125];
verifyEqual(testCase, highOrder.matvec(x), highOrderGram * x, "AbsTol", 2e-12);
clear localCleanup highOrderCleanup
end

function testHACApKCurvedHighOrderTetReferenceDensity(testCase)
vertices = [0, 0, 0; 2, 0, 0; 0, 3, 0; 0, 0, 4];
cellNodes = [vertices; ...
    (vertices(1,:) + vertices(2,:)) / 2; ...
    (vertices(2,:) + vertices(3,:)) / 2; ...
    (vertices(3,:) + vertices(1,:)) / 2; ...
    (vertices(1,:) + vertices(4,:)) / 2; ...
    (vertices(2,:) + vertices(4,:)) / 2; ...
    (vertices(3,:) + vertices(4,:)) / 2];
a = 0.5854101966249685;
b = 0.1381966011250105;
tetPoints = [a, b, b; b, a, b; b, b, a; b, b, b];
tetWeights = ones(4, 1) / 24;
trianglePoints = [1/6, 1/6; 2/3, 1/6; 1/6, 2/3];
triangleWeights = ones(3, 1) / 6;
curvePoints = 0.5 * ([-0.906179845938664; -0.538469310105683; 0; ...
    0.538469310105683; 0.906179845938664] + 1);
curveWeights = 0.5 * [0.236926885056189; 0.478628670499367; ...
    0.568888888888889; 0.478628670499367; 0.236926885056189];
common = {cellNodes, zeros(0,3), int32([0, 1, 2, 3]), int32(zeros(0,3)), ...
    1, int32(0), int32(0), int32([0, 0, 0]), tetPoints, tetWeights, ...
    trianglePoints, triangleWeights, curvePoints, curveWeights};
physical = radia.HACApKChargeGram.from_curved_high_order_tet( ...
    common{:}, ReferenceDensity=false, Build=false);
physicalCleanup = onCleanup(@() delete(physical));
reference = radia.HACApKChargeGram.from_curved_high_order_tet( ...
    common{:}, ReferenceDensity=true, Build=false);
referenceCleanup = onCleanup(@() delete(reference));
physicalEntry = physical.entry(1, 1);
referenceEntry = reference.entry(1, 1);
verifyGreaterThan(testCase, physicalEntry, 0);
verifyGreaterThan(testCase, referenceEntry, 0);
verifyEqual(testCase, physicalEntry / referenceEntry, 24^2, "RelTol", 3e-13);
verifyTrue(testCase, physical.build(AcaEps=1e-12, LeafSize=2));
verifyEqual(testCase, physical.matvec(1), physicalEntry, "RelTol", 2e-12);
clear referenceCleanup physicalCleanup
end

function testHACApKHexFarBlock(testCase)
hexCellNodes = zeros(54, 3);
row = 1;
for offset = [0, 4]
    for iz = 0:2
        for iy = 0:2
            for ix = 0:2
                hexCellNodes(row,:) = [offset + ix/2, iy/2, iz/2];
                row = row + 1;
            end
        end
    end
end
gaussPoints = 0.5 * ([-0.906179845938664; -0.538469310105683; 0; ...
    0.538469310105683; 0.906179845938664] + 1);
gaussWeights = 0.5 * [0.236926885056189; 0.478628670499367; ...
    0.568888888888889; 0.478628670499367; 0.236926885056189];
manager = radia.HACApKChargeGram.from_hex( ...
    hexCellNodes, zeros(0,3), 2, 0, int32([0; 1]), int32([0; 0]), ...
    int32(zeros(2,3)), [0.25, 0.25, 0.25], 1/6, [1/3, 1/3], 1/2, ...
    gaussPoints, gaussWeights, gaussPoints, gaussWeights, ...
    [0.25, 0.25, 0.25], 1/6, [1/3, 1/3], 1/2, ...
    NearGrade=0.5, FarInnerFactor=1, Build=false);
cleanup = onCleanup(@() delete(manager));
verifyEqual(testCase, manager.NDOF, 2);
farEntry = manager.entry(1, 2);
verifyEqual(testCase, farEntry, manager.entry(2, 1), "RelTol", 2e-14);

[qx, qy, qz] = ndgrid(gaussPoints, gaussPoints, gaussPoints);
[wx, wy, wz] = ndgrid(gaussWeights, gaussWeights, gaussWeights);
points = [qx(:), qy(:), qz(:)];
weights = wx(:) .* wy(:) .* wz(:);
reference = 0;
source = points + [4, 0, 0];
for target = 1:size(points,1)
    distance = vecnorm(source - points(target,:), 2, 2);
    reference = reference + weights(target) * sum(weights ./ distance);
end
reference = reference / (4*pi);
verifyEqual(testCase, farEntry, reference, "RelTol", 3e-12);

gram = [manager.entry(1,1), farEntry; farEntry, manager.entry(2,2)];
verifyTrue(testCase, manager.build(AcaEps=1e-12, LeafSize=4));
x = [1; -0.25];
verifyEqual(testCase, manager.matvecSym(x), gram*x, "AbsTol", 2e-11);
clear cleanup
end

function testHACApKHexSelfBlockDirectionalDerivative(testCase)
hexCellNodes = zeros(27, 3);
row = 1;
for iz = 0:2
    for iy = 0:2
        for ix = 0:2
            hexCellNodes(row,:) = [ix/2, iy/2, iz/2];
            row = row + 1;
        end
    end
end
gaussPoints = 0.5 * ([-0.906179845938664; -0.538469310105683; 0; ...
    0.538469310105683; 0.906179845938664] + 1);
gaussWeights = 0.5 * [0.236926885056189; 0.478628670499367; ...
    0.568888888888889; 0.478628670499367; 0.236926885056189];
factory = @(nodes) radia.HACApKChargeGram.from_hex( ...
    nodes, zeros(0,3), 1, 0, int32(zeros(4,1)), int32(zeros(4,1)), ...
    int32([0,0,0; 1,0,0; 0,1,0; 0,0,1]), ...
    [0.25,0.25,0.25], 1/6, [1/3,1/3], 1/2, ...
    gaussPoints, gaussWeights, gaussPoints, gaussWeights, ...
    [0.25,0.25,0.25], 1/6, [1/3,1/3], 1/2, Build=false);
manager = factory(hexCellNodes);
cleanup = onCleanup(@() delete(manager));
selfBlock = zeros(4);
for row = 1:4
    for col = 1:4
        selfBlock(row,col) = manager.entry(row,col);
    end
end
derivative = manager.hexVolumeSelfBlockDirectionalDerivative(1, hexCellNodes);
verifySize(testCase, derivative, [4,4]);
verifyEqual(testCase, derivative, -selfBlock, "RelTol", 3e-13, "AbsTol", 3e-14);
verifyEqual(testCase, derivative, derivative.', "AbsTol", 2e-14);

translation = repmat([0.25, -0.5, 0.75], 27, 1);
translationDerivative = manager.hexVolumeSelfBlockDirectionalDerivative( ...
    1, translation);
verifyEqual(testCase, translationDerivative, zeros(4), "AbsTol", 2e-13);

step = 1e-4;
plus = factory((1 + step) * hexCellNodes);
plusCleanup = onCleanup(@() delete(plus));
minus = factory((1 - step) * hexCellNodes);
minusCleanup = onCleanup(@() delete(minus));
plusBlock = zeros(4);
minusBlock = zeros(4);
for row = 1:4
    for col = 1:4
        plusBlock(row,col) = plus.entry(row,col);
        minusBlock(row,col) = minus.entry(row,col);
    end
end
finiteDifference = (plusBlock - minusBlock) / (2*step);
verifyLessThan(testCase, norm(derivative - finiteDifference, "fro") / ...
    norm(finiteDifference, "fro"), 2e-8);

affineVelocity = zeros(27, 3);
affineVelocity(:,1) = 0.15*hexCellNodes(:,1) + 0.07*hexCellNodes(:,2);
affineVelocity(:,2) = -0.04*hexCellNodes(:,1) + 0.11*hexCellNodes(:,2);
affineVelocity(:,3) = 0.05*hexCellNodes(:,1) - 0.09*hexCellNodes(:,3);
affineDerivative = manager.hexVolumeSelfBlockDirectionalDerivative( ...
    1, affineVelocity);
plusAffine = factory(hexCellNodes + step*affineVelocity);
plusAffineCleanup = onCleanup(@() delete(plusAffine));
minusAffine = factory(hexCellNodes - step*affineVelocity);
minusAffineCleanup = onCleanup(@() delete(minusAffine));
plusAffineBlock = zeros(4);
minusAffineBlock = zeros(4);
for row = 1:4
    for col = 1:4
        plusAffineBlock(row,col) = plusAffine.entry(row,col);
        minusAffineBlock(row,col) = minusAffine.entry(row,col);
    end
end
affineFiniteDifference = (plusAffineBlock - minusAffineBlock) / (2*step);
verifyLessThan(testCase, norm(affineDerivative - affineFiniteDifference, "fro") / ...
    norm(affineFiniteDifference, "fro"), 6e-6);
clear minusAffineCleanup plusAffineCleanup
clear minusCleanup plusCleanup cleanup
end

function testHACApKNonAffineHexSelfBlockDirectionalDerivative(testCase)
hexCellNodes = zeros(27, 3);
row = 1;
for iz = 0:2
    for iy = 0:2
        for ix = 0:2
            hexCellNodes(row,:) = [ix/2, iy/2, iz/2];
            row = row + 1;
        end
    end
end
hexCellNodes(14,:) = hexCellNodes(14,:) + [0.03, -0.02, 0.04];
gaussPoints = 0.5 * ([-0.906179845938664; -0.538469310105683; 0; ...
    0.538469310105683; 0.906179845938664] + 1);
gaussWeights = 0.5 * [0.236926885056189; 0.478628670499367; ...
    0.568888888888889; 0.478628670499367; 0.236926885056189];
factory = @(nodes) radia.HACApKChargeGram.from_hex( ...
    nodes, zeros(0,3), 1, 0, int32(0), int32(0), int32([0,0,0]), ...
    [0.25,0.25,0.25], 1/6, [1/3,1/3], 1/2, ...
    gaussPoints, gaussWeights, gaussPoints, gaussWeights, ...
    [0.25,0.25,0.25], 1/6, [1/3,1/3], 1/2, Build=false);
manager = factory(hexCellNodes);
cleanup = onCleanup(@() delete(manager));
selfEntry = manager.entry(1,1);
derivative = manager.hexVolumeSelfBlockDirectionalDerivative(1, hexCellNodes);
verifyEqual(testCase, derivative, -selfEntry, "RelTol", 3e-11);

step = 1e-3;
plus = factory((1 + step) * hexCellNodes);
plusCleanup = onCleanup(@() delete(plus));
minus = factory((1 - step) * hexCellNodes);
minusCleanup = onCleanup(@() delete(minus));
finiteDifference = (plus.entry(1,1) - minus.entry(1,1)) / (2*step);
verifyEqual(testCase, derivative, finiteDifference, "RelTol", 2e-6);
clear minusCleanup plusCleanup cleanup
end

function testHACApKHexFaceSelfBlockDirectionalDerivative(testCase)
hexCellNodes = zeros(27,3);
row = 1;
for ix = 0:2
    for iy = 0:2
        for iz = 0:2
            hexCellNodes(row,:) = [ix/2,iy/2,iz/2];
            row = row + 1;
        end
    end
end
quadFaceNodes = zeros(54,3);
row = 1;
for face = 1:6
    for iu = 0:2
        for iv = 0:2
            u = iu/2;
            v = iv/2;
            if face == 1
                point = [0,u,v];
            elseif face == 2
                point = [u,0,v];
            elseif face == 3
                point = [1,1-u,1-v];
            elseif face == 4
                point = [1-u,1,1-v];
            elseif face == 5
                point = [u,v,0];
            else
                point = [1-u,1-v,1];
            end
            quadFaceNodes(row,:) = point;
            row = row + 1;
        end
    end
end
gaussPoints = 0.5 * ([-0.906179845938664; -0.538469310105683; 0; ...
    0.538469310105683; 0.906179845938664] + 1);
gaussWeights = 0.5 * [0.236926885056189; 0.478628670499367; ...
    0.568888888888889; 0.478628670499367; 0.236926885056189];
volumeExponents = zeros(27,3);
row = 1;
for iz = 0:2
    for iy = 0:2
        for ix = 0:2
            volumeExponents(row,:) = [ix,iy,iz];
            row = row + 1;
        end
    end
end
oneFaceExponents = zeros(9,3);
row = 1;
for iy = 0:2
    for ix = 0:2
        oneFaceExponents(row,:) = [ix,iy,0];
        row = row + 1;
    end
end
faceExponents = repmat(oneFaceExponents,6,1);
manager = radia.HACApKChargeGram.from_hex( ...
    hexCellNodes, quadFaceNodes, 1, 6, ...
    int32([zeros(27,1);repelem((0:5).',9)]), ...
    int32([zeros(27,1);ones(54,1)]), ...
    int32([volumeExponents;faceExponents]), ...
    [0.25,0.25,0.25], 1/6, [1/3,1/3], 1/2, ...
    gaussPoints, gaussWeights, gaussPoints, gaussWeights, ...
    [0.25,0.25,0.25], 1/6, [1/3,1/3], 1/2, Build=false);
cleanup = onCleanup(@() delete(manager));
translation = repmat([0.25,-0.5,0.75],9,1);
derivative = manager.hexFaceSelfBlockDirectionalDerivative(5,translation);
verifySize(testCase, derivative, [9,9]);
verifyEqual(testCase, derivative, zeros(9), "AbsTol", 3e-12);
translation3 = reshape([0.25,-0.5,0.75],1,1,3);
    completeDerivative = manager.hexChargeGramDirectionalDerivative( ...
        repmat(translation3,1,27,1), repmat(translation3,6,9,1));
    verifySize(testCase, completeDerivative, [81,81]);
    verifyEqual(testCase, completeDerivative, zeros(81), "AbsTol", 2e-11);
    cellScaling = reshape(hexCellNodes,[1,27,3]);
    faceScaling = permute(reshape(quadFaceNodes,[9,6,3]),[2,1,3]);
    scalingDerivative = manager.hexChargeGramDirectionalDerivative( ...
        cellScaling,faceScaling);
    verifyTrue(testCase,manager.build( ...
        AcaEps=1e-10,LeafSize=8,Eta=1,MaxRank=81));
    cellBatch = cat(1,reshape(cellScaling,[1,1,27,3]), ...
        reshape(repmat(translation3,1,27,1),[1,1,27,3]));
    faceBatch = cat(1,reshape(faceScaling,[1,6,9,3]), ...
        reshape(repmat(translation3,6,9,1),[1,6,9,3]));
    left = linspace(-0.3,0.8,81).';
    right = linspace(0.7,-0.2,81).';
    contractions = manager.directionalDerivativeContractions( ...
        "hex",cellBatch,faceBatch,left,right);
    verifyEqual(testCase,contractions, ...
        [left.'*scalingDerivative*right;left.'*completeDerivative*right], ...
        "RelTol",2e-6,"AbsTol",2e-9);
    derivativeOperator = manager.directionalDerivativeOperator("hex", ...
    repmat(translation3,1,27,1), repmat(translation3,6,9,1), ...
    AcaEps=1e-12, LeafSize=128, Eta=2);
operatorCleanup = onCleanup(@() delete(derivativeOperator));
verifyEqual(testCase, derivativeOperator.NDOF, 81);
operatorInfo = derivativeOperator.info();
verifyTrue(testCase, operatorInfo.valid);
verifyTrue(testCase, isfinite(derivativeOperator.entry(1,1)));
delete(manager);
verifyEqual(testCase, derivativeOperator.matvecSym(ones(81,1)), ...
    completeDerivative*ones(81,1), "AbsTol", 2e-10);
clear operatorCleanup cleanup
end

function testHACApKWedgeFarBlock(testCase)
triangleNodes = [1, 0; 0, 1; 0, 0; 0.5, 0.5; 0, 0.5; 0.5, 0];
wedgeCellNodes = zeros(36, 3);
row = 1;
for offset = [0, 4]
    for z = [0, 0.5, 1]
        for node = 1:6
            wedgeCellNodes(row,:) = [offset + triangleNodes(node,1), ...
                triangleNodes(node,2), z];
            row = row + 1;
        end
    end
end
gaussPoints = 0.5 * ([-0.906179845938664; -0.538469310105683; 0; ...
    0.538469310105683; 0.906179845938664] + 1);
gaussWeights = 0.5 * [0.236926885056189; 0.478628670499367; ...
    0.568888888888889; 0.478628670499367; 0.236926885056189];
trianglePoints = [1/6, 1/6; 2/3, 1/6; 1/6, 2/3];
triangleWeights = ones(3,1) / 6;
manager = radia.HACApKChargeGram.from_wedge( ...
    wedgeCellNodes, zeros(0,3), int32.empty, 2, 0, ...
    int32([0; 1]), int32([0; 0]), int32(zeros(2,3)), ...
    [0.25, 0.25, 0.25], 1/6, [1/3, 1/3], 1/2, ...
    trianglePoints, triangleWeights, gaussPoints, gaussWeights, ...
    gaussPoints, gaussWeights, [0.25, 0.25, 0.25], 1/6, ...
    [1/3, 1/3], 1/2, NearGrade=0.5, FarInnerFactor=1, Build=false);
cleanup = onCleanup(@() delete(manager));
verifyEqual(testCase, manager.NDOF, 2);
farEntry = manager.entry(1, 2);
verifyEqual(testCase, farEntry, manager.entry(2, 1), "RelTol", 2e-14);

points = [0.25, 0.5, 0.25; 0.5, 0.25, 0.5; 0.25, 0.25, 0.75];
weights = ones(3,1) / 6;
source = points + [4, 0, 0];
reference = 0;
for target = 1:size(points,1)
    distance = vecnorm(source - points(target,:), 2, 2);
    reference = reference + weights(target) * sum(weights ./ distance);
end
reference = reference / (4*pi);
verifyEqual(testCase, farEntry, reference, "RelTol", 3e-12);

gram = [manager.entry(1,1), farEntry; farEntry, manager.entry(2,2)];
verifyTrue(testCase, manager.build(AcaEps=1e-12, LeafSize=4));
x = [1; -0.25];
verifyEqual(testCase, manager.matvecSym(x), gram*x, "AbsTol", 2e-11);
clear cleanup
end

function testHACApKTetAndWedgeSelfBlockDirectionalDerivatives(testCase)
a = 0.5854101966249685;
b = 0.1381966011250105;
tetPoints = [a,b,b; b,a,b; b,b,a; b,b,b];
tetWeights = ones(4,1)/24;
trianglePoints = [1/6,1/6; 2/3,1/6; 1/6,2/3];
triangleWeights = ones(3,1)/6;
gaussPoints = 0.5 * ([-0.906179845938664; -0.538469310105683; 0; ...
    0.538469310105683; 0.906179845938664] + 1);
gaussWeights = 0.5 * [0.236926885056189; 0.478628670499367; ...
    0.568888888888889; 0.478628670499367; 0.236926885056189];

cellVerts = [0,0,0; 1,0,0; 0,1,0; 0,0,1];
faceVerts = cellVerts([1,3,2],:);
tet = radia.HACApKChargeGram.from_high_order_tet( ...
    cellVerts, faceVerts, 1, int32([0;0]), int32([0;1]), ...
    int32(zeros(2,3)), tetPoints, tetWeights, trianglePoints, ...
    triangleWeights, InnerTetPoints=tetPoints, ...
    InnerTetWeights=tetWeights, InnerTrianglePoints=trianglePoints, ...
    InnerTriangleWeights=triangleWeights, Build=false);
tetCleanup = onCleanup(@() delete(tet));
tetVolume = tet.tetVolumeSelfBlockDirectionalDerivative( ...
    1, repmat([0.25,-0.5,0.75],4,1));
tetFace = tet.tetFaceSelfBlockDirectionalDerivative( ...
    1, repmat([0.25,-0.5,0.75],3,1));
verifySize(testCase, tetVolume, [1,1]);
verifySize(testCase, tetFace, [1,1]);
verifyEqual(testCase, tetVolume, 0, "AbsTol", 2e-13);
verifyEqual(testCase, tetFace, 0, "AbsTol", 2e-13);
tetCellVelocity = reshape(repmat([0.25,-0.5,0.75],4,1),1,4,3);
tetFaceVelocity = reshape(repmat([0.25,-0.5,0.75],3,1),1,3,3);
tetComplete = tet.tetChargeGramDirectionalDerivative( ...
    tetCellVelocity,tetFaceVelocity);
tetRates = tet.tetChargeMapRowDirectionalRates( ...
    tetCellVelocity,tetFaceVelocity);
verifySize(testCase, tetComplete, [2,2]);
verifyEqual(testCase, tetComplete, zeros(2), "AbsTol", 3e-13);
verifyEqual(testCase, tetRates, zeros(2,1), "AbsTol", 3e-14);

triangleNodes = [1,0; 0,1; 0,0; 0.5,0.5; 0,0.5; 0.5,0];
wedgeCellNodes = zeros(18,3);
row = 1;
for z = [0,0.5,1]
    for node = 1:6
        wedgeCellNodes(row,:) = [triangleNodes(node,:),z];
        row = row + 1;
    end
end
wedgeFaceNodes = zeros(9,3);
wedgeFaceNodes(1:6,:) = [triangleNodes,zeros(6,1)];
wedge = radia.HACApKChargeGram.from_wedge( ...
    wedgeCellNodes, wedgeFaceNodes, int32(0), 1, 1, ...
    int32([0;0]), int32([0;1]), int32(zeros(2,3)), ...
    tetPoints, tetWeights, trianglePoints, triangleWeights, ...
    trianglePoints, triangleWeights, gaussPoints, gaussWeights, ...
    gaussPoints, gaussWeights, tetPoints, tetWeights, ...
    trianglePoints, triangleWeights, Build=false);
wedgeCleanup = onCleanup(@() delete(wedge));
wedgeVolume = wedge.wedgeVolumeSelfBlockDirectionalDerivative( ...
    1, repmat([0.25,-0.5,0.75],18,1));
wedgeFace = wedge.wedgeFaceSelfBlockDirectionalDerivative( ...
    1, repmat([0.25,-0.5,0.75],6,1));
verifySize(testCase, wedgeVolume, [1,1]);
verifySize(testCase, wedgeFace, [1,1]);
verifyEqual(testCase, wedgeVolume, 0, "AbsTol", 3e-12);
verifyEqual(testCase, wedgeFace, 0, "AbsTol", 3e-12);
wedgeCellVelocity = reshape(repmat([0.25,-0.5,0.75],18,1),1,18,3);
wedgeFaceVelocity = reshape(repmat([0.25,-0.5,0.75],9,1),1,9,3);
wedgeComplete = wedge.wedgeChargeGramDirectionalDerivative( ...
    wedgeCellVelocity,wedgeFaceVelocity);
verifySize(testCase, wedgeComplete, [2,2]);
verifyEqual(testCase, wedgeComplete, zeros(2), "AbsTol", 4e-12);
clear wedgeCleanup tetCleanup
end

function testHACApKConfiguredFieldRowsAndBatchedDerivatives(testCase)
a = 0.5854101966249685;
b = 0.1381966011250105;
tetPoints = [a,b,b; b,a,b; b,b,a; b,b,b];
tetWeights = ones(4,1)/24;
trianglePoints = [1/6,1/6; 2/3,1/6; 1/6,2/3];
triangleWeights = ones(3,1)/6;
cellVerts = [0,0,0; 1,0,0; 0,1,0; 0,0,1];
faceVerts = cellVerts([1,3,2],:);
manager = radia.HACApKChargeGram.from_high_order_tet( ...
    cellVerts, faceVerts, 1, int32([0;0]), int32([0;1]), ...
    int32(zeros(2,3)), tetPoints, tetWeights, trianglePoints, ...
    triangleWeights, InnerTetPoints=tetPoints, ...
    InnerTetWeights=tetWeights, InnerTrianglePoints=trianglePoints, ...
    InnerTriangleWeights=triangleWeights, Build=false);
cleanup = onCleanup(@() delete(manager));
verifyTrue(testCase, manager.build(AcaEps=1e-12, LeafSize=4));
manager.configureChargeMap(int32([0;1;2]), int32([0;1]), [1;1], 2);

observations = [2,2,2; -0.5,0.25,0.3];
weights = zeros(2,2,3);
weights(1,:,1) = [1,-0.5];
weights(2,:,2) = [0.25,1.5];
rows = manager.configuredFieldFunctionalRows(observations, weights);
verifySize(testCase, rows, [2,2]);
verifyTrue(testCase, all(isfinite(rows), "all"));

cellVelocity = zeros(1,1,4,3);
faceVelocity = zeros(1,1,3,3);
derivative = manager.configuredFieldFunctionalRowsDirectionalDerivative( ...
    observations, weights, cellVelocity, faceVelocity);
verifySize(testCase, derivative, [1,2,2]);
verifyEqual(testCase, derivative, zeros(1,2,2), "AbsTol", 1e-14);

left = [1,-0.25; 0.5,2];
right = [0.75;-1];
contractions = manager.directionalDerivativeContractionsMany( ...
    "tet", cellVelocity, faceVelocity, left, right);
verifySize(testCase, contractions, [2,1]);
verifyEqual(testCase, contractions, zeros(2,1), "AbsTol", 1e-14);
clear cleanup
end

function destroyIHHandle(command, handle)
try
    if handle ~= 0
        radia.internal.callMex(command, handle);
    end
catch
end
end

function testHACApKPlanar2DCellEdgeFarBlock(testCase)
cellMap = [0, 0; 1, 0; 0, 1; 0, 0];
edgeMap = [4, 0; 0, 1];
trianglePoints = [1/6, 1/6; 2/3, 1/6; 1/6, 2/3];
triangleWeights = ones(3,1) / 6;
gaussPoints = 0.5 * ([-0.906179845938664; -0.538469310105683; 0; ...
    0.538469310105683; 0.906179845938664] + 1);
gaussWeights = 0.5 * [0.236926885056189; 0.478628670499367; ...
    0.568888888888889; 0.478628670499367; 0.236926885056189];
manager = radia.HACApKChargeGram.from_planar_2d( ...
    2, 1, cellMap, int32(0), edgeMap, 1, 1, int32([0; 0]), ...
    int32([0; 1]), int32(zeros(2,3)), trianglePoints, triangleWeights, ...
    gaussPoints, gaussWeights, gaussPoints, gaussWeights, ...
    gaussPoints, gaussWeights, trianglePoints, triangleWeights, ...
    NearGrade=0.5, FarInnerFactor=1, Build=false);
cleanup = onCleanup(@() delete(manager));
verifyEqual(testCase, manager.NDOF, 2);
farEntry = manager.entry(1, 2);
verifyEqual(testCase, farEntry, manager.entry(2, 1), "RelTol", 2e-14);

referenceTrianglePoints = [2/3, 1/6; 1/6, 2/3; 1/6, 1/6];
edgePoints = [4 * ones(size(gaussPoints)), gaussPoints];
reference = 0;
for target = 1:size(referenceTrianglePoints,1)
    distance = vecnorm(edgePoints - referenceTrianglePoints(target,:), 2, 2);
    reference = reference - triangleWeights(target) * ...
        sum(gaussWeights .* log(distance));
end
reference = reference / (2*pi);
verifyEqual(testCase, farEntry, reference, "RelTol", 3e-12);

gram = [manager.entry(1,1), farEntry; farEntry, manager.entry(2,2)];
verifyTrue(testCase, manager.build(AcaEps=1e-12, LeafSize=4));
x = [0.75; -0.4];
verifyEqual(testCase, manager.matvecSym(x), gram*x, "AbsTol", 2e-11);
clear cleanup
end

function testHACApKCurvedPolytopeMatchesCurvedTet(testCase)
vertices = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];
cellVerticesGeometry = cat(3, vertices, ...
    [0, 0, 0; 1, 0, 0; 0, -1, 0; 0, 0, -1]);
cellNodes = zeros(20, 3);
for cell = 1:2
    v = cellVerticesGeometry(:,:,cell);
    rows = 10*(cell-1) + (1:10);
    cellNodes(rows,:) = [v; ...
        (v(1,:) + v(2,:))/2; (v(2,:) + v(3,:))/2; ...
        (v(3,:) + v(1,:))/2; (v(1,:) + v(4,:))/2; ...
        (v(2,:) + v(4,:))/2; (v(3,:) + v(4,:))/2];
end
a = 0.5854101966249685;
b = 0.1381966011250105;
tetPoints = [a, b, b; b, a, b; b, b, a; b, b, b];
tetWeights = ones(4, 1) / 24;
trianglePoints = [1/6, 1/6; 2/3, 1/6; 1/6, 2/3];
triangleWeights = ones(3, 1) / 6;
curvePoints = 0.5 * ([-0.906179845938664; -0.538469310105683; 0; ...
    0.538469310105683; 0.906179845938664] + 1);
curveWeights = 0.5 * [0.236926885056189; 0.478628670499367; ...
    0.568888888888889; 0.478628670499367; 0.236926885056189];

polytope = radia.HACApKChargeGram.from_curved_polytope( ...
    cellNodes, int32([0; 1; 2]), [0.25, 0.25, 0.25; 0.25, -0.25, -0.25], ...
    [1/6; 1/6], zeros(0,3), int32(0), zeros(0,3), zeros(0,1), ...
    tetPoints, tetWeights, trianglePoints, triangleWeights, ...
    curvePoints, curveWeights, 2, Build=false);
polytopeCleanup = onCleanup(@() delete(polytope));
highOrder = radia.HACApKChargeGram.from_curved_high_order_tet( ...
    cellNodes, zeros(0,3), int32([0,1,2,3; 0,1,4,5]), ...
    int32(zeros(0,3)), 2, int32([0;1]), int32([0;0]), ...
    int32(zeros(2,3)), tetPoints, tetWeights, trianglePoints, ...
    triangleWeights, curvePoints, curveWeights, Build=false);
highOrderCleanup = onCleanup(@() delete(highOrder));
verifyEqual(testCase, polytope.NDOF, 2);
polytopeGram = zeros(2);
highOrderGram = zeros(2);
for i = 1:2
    for j = 1:2
        polytopeGram(i,j) = polytope.entry(i,j);
        highOrderGram(i,j) = highOrder.entry(i,j);
    end
end
verifyEqual(testCase, polytopeGram, highOrderGram, "RelTol", 3e-13, ...
    "AbsTol", 2e-15);
verifyTrue(testCase, polytope.build(AcaEps=1e-12, LeafSize=4));
x = [1; -0.3];
verifyEqual(testCase, polytope.matvecSym(x), polytopeGram*x, "AbsTol", 2e-11);
clear highOrderCleanup polytopeCleanup
end

function testHACApKAnalyticTetFactory(testCase)
cellVerts = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];
faceVerts = [ ...
    0, 0, 0; 0, 1, 0; 1, 0, 0; ...
    0, 0, 0; 1, 0, 0; 0, 0, 1; ...
    0, 0, 0; 0, 0, 1; 0, 1, 0; ...
    1, 0, 0; 0, 1, 0; 0, 0, 1];
manager = radia.HACApKChargeGram.from_analytic_tet( ...
    cellVerts, faceVerts, 1, Build=false);
cleanup = onCleanup(@() delete(manager));
verifyEqual(testCase, manager.NDOF, 5);
gram = zeros(5);
for i = 1:5
    for j = 1:5
        gram(i,j) = manager.entry(i, j);
    end
end
verifyTrue(testCase, all(isfinite(gram), "all"));
verifyGreaterThan(testCase, min(diag(gram)), 0);
verifyEqual(testCase, gram, gram.', "RelTol", 2e-12);
verifyTrue(testCase, manager.build(AcaEps=1e-10, LeafSize=8));
x = [1; -0.5; 0.25; -0.125; 0.0625];
verifyEqual(testCase, manager.matvec(x), gram * x, "AbsTol", 2e-10);
clear cleanup
end

function testHACApKAnalyticPolytopeMatchesTet(testCase)
cellVerts = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];
faceVerts = [ ...
    0, 0, 0; 0, 1, 0; 1, 0, 0; ...
    0, 0, 0; 1, 0, 0; 0, 0, 1; ...
    0, 0, 0; 0, 0, 1; 0, 1, 0; ...
    1, 0, 0; 0, 1, 0; 0, 0, 1];
faceCentroids = zeros(4, 3);
faceMeasures = zeros(4, 1);
for face = 1:4
    triangle = faceVerts((3*face-2):(3*face), :);
    faceCentroids(face,:) = mean(triangle, 1);
    faceMeasures(face) = norm(cross( ...
        triangle(2,:) - triangle(1,:), triangle(3,:) - triangle(1,:))) / 2;
end
tet = radia.HACApKChargeGram.from_analytic_tet( ...
    cellVerts, faceVerts, 1, Build=false);
tetCleanup = onCleanup(@() delete(tet));
polytope = radia.HACApKChargeGram.from_analytic_polytope( ...
    faceVerts, int32([0; 4]), [0.25, 0.25, 0.25], 1/6, ...
    faceVerts, int32((0:4).'), faceCentroids, faceMeasures, 1, Build=false);
polytopeCleanup = onCleanup(@() delete(polytope));
tetGram = zeros(5);
polytopeGram = zeros(5);
for i = 1:5
    for j = 1:5
        tetGram(i,j) = tet.entry(i, j);
        polytopeGram(i,j) = polytope.entry(i, j);
    end
end
verifyEqual(testCase, polytopeGram(2:end,2:end), tetGram(2:end,2:end), ...
    "RelTol", 5e-12, "AbsTol", 5e-14);
relativeDifference = abs(polytopeGram - tetGram) ./ max(abs(tetGram), 1e-14);
verifyLessThan(testCase, max(relativeDifference, [], "all"), 5e-4);
verifyTrue(testCase, polytope.build(AcaEps=1e-10, LeafSize=8));
x = [1; -0.5; 0.25; -0.125; 0.0625];
verifyEqual(testCase, polytope.matvec(x), polytopeGram * x, "AbsTol", 2e-10);
clear polytopeCleanup tetCleanup
end

function testPersistentHDivFieldEvaluators(testCase)
cloud = radia.HDivFieldEvaluator.fromCloud([0, 0, 0; 1, 0, 0], [1; -0.25]);
cleanup = onCleanup(@() delete(cloud));
field = cloud.field([0, 1, 0; 0, 0, 2], Algorithm="direct");
verifySize(testCase, field, [2, 3]);
verifyTrue(testCase, all(isfinite(field), "all"));
stats = cloud.stats();
verifyEqual(testCase, stats.source_count, 2);
verifyEqual(testCase, cloud.candidateAlgorithm(2), "direct");
verifyEqual(testCase, cloud.lastAlgorithm(), "direct");
volumePoints = [0.15,0.20,0.10; 0.25,0.10,0.15];
fieldCoefficient = cloud.coefficientFunction(Algorithm="direct");
coefficientCleanup = onCleanup(@() delete(fieldCoefficient));
fieldCoefficientInfo = fieldCoefficient.info();
verifyEqual(testCase, string(fieldCoefficientInfo.kind), "hdiv_field");
verifyEqual(testCase, string(fieldCoefficientInfo.algorithm), "direct");
verifyEqual(testCase, ...
    fieldCoefficient.evaluate(testCase.TestData.meshPath, volumePoints), ...
    cloud.field(volumePoints, Algorithm="direct")/(4*pi), ...
    "RelTol", 3e-13, "AbsTol", 3e-15);

planar = radia.PlanarFieldEvaluator.create([0, 0; 1, 0], [1; -0.5]);
planarCleanup = onCleanup(@() delete(planar));
planarField = planar.field([0.25, 0.5; 2, 0.2]);
verifySize(testCase, planarField, [2, 2]);
verifySize(testCase, planar.az([0.25, 0.5; 2, 0.2]), [2, 1]);
verifyTrue(testCase, all(isfinite(planarField), "all"));
verifyEqual(testCase, planar.stats().base_source_count, 2);
planarCoefficient = planar.coefficientFunction();
planarCoefficientCleanup = onCleanup(@() delete(planarCoefficient));
planarCoefficientInfo = planarCoefficient.info();
verifyEqual(testCase, string(planarCoefficientInfo.kind), "planar_hdiv_field");
verifyEqual(testCase, planarCoefficientInfo.source_angle, 0);
verifyEqual(testCase, planarCoefficientInfo.target_angle, 0);
verifyEqual(testCase, ...
    planarCoefficient.evaluate(testCase.TestData.meshPath, volumePoints), ...
    planar.field(volumePoints(:,1:2)), "RelTol", 3e-14, "AbsTol", 2e-15);
clear planarCoefficientCleanup planarCleanup coefficientCleanup cleanup
end

function testChargeGramConfiguredFieldEvaluator(testCase)
cellVerts = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];
a = 0.5854101966249685;
b = 0.1381966011250105;
tetPoints = [a, b, b; b, a, b; b, b, a; b, b, b];
tetWeights = ones(4, 1) / 24;
trianglePoints = [1/6, 1/6; 2/3, 1/6; 1/6, 2/3];
triangleWeights = ones(3, 1) / 6;
manager = radia.HACApKChargeGram.from_high_order_tet( ...
    cellVerts, zeros(0,3), 1, int32(0), int32(0), int32([0,0,0]), ...
    tetPoints, tetWeights, trianglePoints, triangleWeights, Build=false);
managerCleanup = onCleanup(@() delete(manager));
verifyTrue(testCase, manager.build(AcaEps=1e-12, LeafSize=2));
manager.configureChargeMap(int32([0;1]), int32(0), 1, 1);
configured = manager.createFieldEvaluator(1);
configuredCleanup = onCleanup(@() delete(configured));

volume = [reshape(cellVerts.', 1, []), 1, 0, 0, 0];
reference = radia.HDivFieldEvaluator.fromTet(volume, zeros(0,22));
referenceCleanup = onCleanup(@() delete(reference));
observations = [2, 2, 2; -0.5, 0.25, 0.3; 0.3, 0.2, 1.5];
verifyEqual(testCase, configured.field(observations, Algorithm="direct"), ...
    reference.field(observations, Algorithm="direct"), ...
    "RelTol", 3e-13, "AbsTol", 2e-15);
clear referenceCleanup configuredCleanup managerCleanup
end

function testChargeGramConfiguredPlanarFieldEvaluator(testCase)
cellMap = [0, 0; 1, 0; 0, 1; 0, 0];
edgeMap = [4, 0; 0, 1];
trianglePoints = [1/6, 1/6; 2/3, 1/6; 1/6, 2/3];
triangleWeights = ones(3,1) / 6;
gaussPoints = 0.5 * ([-0.906179845938664; -0.538469310105683; 0; ...
    0.538469310105683; 0.906179845938664] + 1);
gaussWeights = 0.5 * [0.236926885056189; 0.478628670499367; ...
    0.568888888888889; 0.478628670499367; 0.236926885056189];
manager = radia.HACApKChargeGram.from_planar_2d( ...
    2, 1, cellMap, int32(0), edgeMap, 1, 1, int32([0;0]), ...
    int32([0;1]), int32(zeros(2,3)), trianglePoints, triangleWeights, ...
    gaussPoints, gaussWeights, gaussPoints, gaussWeights, ...
    gaussPoints, gaussWeights, trianglePoints, triangleWeights, Build=false);
managerCleanup = onCleanup(@() delete(manager));
verifyTrue(testCase, manager.build(AcaEps=1e-12, LeafSize=4));
manager.configureChargeMap(int32([0;1;2]), int32([0;1]), [1;1], 2);
configured = manager.createPlanarFieldEvaluator([1; -0.5]);
configuredCleanup = onCleanup(@() delete(configured));

vertices = [1, 0; 0, 1; 0, 0];
mid01 = (vertices(1,:) + vertices(2,:)) / 2;
mid12 = (vertices(2,:) + vertices(3,:)) / 2;
mid20 = (vertices(3,:) + vertices(1,:)) / 2;
subtriangles = cat(3, [vertices(1,:); mid01; mid20], ...
    [mid01; vertices(2,:); mid12], [mid20; mid12; vertices(3,:)], ...
    [mid01; mid12; mid20]);
cellPoints = zeros(12,2);
cellWeights = zeros(12,1);
row = 1;
for subtriangle = 1:4
    tri = subtriangles(:,:,subtriangle);
    edge1 = tri(2,:) - tri(1,:);
    edge2 = tri(3,:) - tri(1,:);
    scale = abs(det([edge1; edge2]));
    for q = 1:3
        cellPoints(row,:) = tri(1,:) + trianglePoints(q,1)*edge1 + ...
            trianglePoints(q,2)*edge2;
        cellWeights(row) = scale * triangleWeights(q);
        row = row + 1;
    end
end
edgePoints = [4*ones(size(gaussPoints)), gaussPoints];
reference = radia.PlanarFieldEvaluator.create( ...
    [cellPoints; edgePoints], [cellWeights; -0.5*gaussWeights]);
referenceCleanup = onCleanup(@() delete(reference));
observations = [2, 1.5; -0.5, 0.25; 5, 2];
verifyEqual(testCase, configured.field(observations), reference.field(observations), ...
    "RelTol", 3e-14, "AbsTol", 2e-15);
verifyEqual(testCase, configured.az(observations), reference.az(observations), ...
    "RelTol", 3e-14, "AbsTol", 2e-15);
clear referenceCleanup configuredCleanup managerCleanup
end

function testHACApKGeometryMassWithoutChargeMap(testCase)
manager = radia.HACApKChargeGram( ...
    [0, 0, 0; 1, 0, 0; 0, 1, 0], [1; 1; 1], [0; 0; 0]);
cleanup = onCleanup(@() delete(manager));
manager.configureGeometryMassMatrix( ...
    int32([0; 1; 2]), int32([0; 1; 2]), [2; 3; 4], 3);
verifyEqual(testCase, manager.applyConfiguredGeometryMass([1; 2; 3]), ...
    [2; 6; 12], "AbsTol", 1e-14);
verifyEqual(testCase, manager.operatorInfo().constraint_count, 0);
clear cleanup
end

function testEnergyStopLifecycle(testCase)
before = radia.apiInfo();
material = radia.EnergyStopMaterial( ...
    0.2, {[0, 0; 0.2, 1000]}, Alpha=5, Gamma=0, BMax=1);
created = radia.apiInfo();
verifyEqual(testCase, created.handle_count, before.handle_count + 1);
verifyEqual(testCase, material.StateSize, 6);

B = [0.05, 0, 0; 0.15, 0.02, 0; 0.25, 0.03, 0.01];
states = repmat(material.state0(), size(B, 1), 1);
H = material.forward(B, states);
newStates = material.commit(B, states);
energy = material.storedEnergy(B, states);
verifySize(testCase, H, [3, 3]);
verifySize(testCase, newStates, [3, 6]);
verifySize(testCase, energy, [3, 1]);
verifyTrue(testCase, all(isfinite(H), "all"));
verifyTrue(testCase, all(isfinite(newStates), "all"));
verifyGreaterThanOrEqual(testCase, min(energy), 0);

delete(material);
after = radia.apiInfo();
verifyEqual(testCase, after.handle_count, before.handle_count);
end

function testLegacyTetraFieldRoute(testCase)
radia.UtiDelAll();
cleanup = onCleanup(@() radia.UtiDelAll());
vertices = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];
object = radia.ObjTetrahedron(vertices, [0, 0, 1]);
container = radia.ObjCnt(object);

verifyEqual(testCase, radia.ObjGeoVol(object), 1 / 6, "AbsTol", 1e-14);
verifyGreaterThanOrEqual(testCase, radia.ObjDegFre(container), 0);

points = [2, 2, 2; 3, 2, 1];
B = radia.Fld(container, "b", points);
H = radia.Fld(container, "h", points);
A = radia.Fld(container, "a", points);
phi = radia.Fld(container, "phi", points);
verifySize(testCase, B, [2, 3]);
verifySize(testCase, H, [2, 3]);
verifySize(testCase, A, [2, 3]);
verifySize(testCase, phi, [2, 1]);
verifyTrue(testCase, all(isfinite([B, H, A, phi]), "all"));
verifyGreaterThan(testCase, norm(B, "fro"), 0);
clear cleanup
end

function testLegacyObjectsAndTransforms(testCase)
radia.UtiDelAll();
cleanup = onCleanup(@() radia.UtiDelAll());
vertices = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];
object = radia.ObjTetrahedron(vertices, [0, 0, 1]);
state = radia.ObjM(object);
verifySize(testCase, state.center, [1, 3]);
verifyEqual(testCase, state.magnetization, [0, 0, 1], "AbsTol", 1e-14);

radia.ObjSetM(object, [0.1, 0.2, 0.3]);
state = radia.ObjM(object);
verifyEqual(testCase, state.magnetization, [0.1, 0.2, 0.3], ...
    "AbsTol", 1e-14);

copy = radia.ObjDpl(object);
container = radia.ObjCnt(object);
radia.ObjAddToCnt(container, copy);
verifyEqual(testCase, radia.ObjCntSize(container), 2);
verifyEqual(testCase, sort(radia.ObjCntStuf(container)), sort([object, copy]));

translation = radia.TrfTrsl([1, 2, 3]);
rotation = radia.TrfRot([0, 0, 0], [0, 0, 1], pi / 3);
inversion = radia.TrfInv();
verifyGreaterThan(testCase, radia.TrfCmbL(translation, rotation), 0);
verifyGreaterThan(testCase, radia.TrfCmbR(rotation, inversion), 0);
moved = radia.TrfOrnt(copy, translation);
verifyEqual(testCase, radia.ObjGeoVol(moved), 1 / 6, "AbsTol", 1e-14);

verifyGreaterThan(testCase, radia.MatPM(1.2, 900000, [0, 0, 1]), 0);
verifyGreaterThan(testCase, radia.UtiVer(), 0);
clear cleanup
end

function testExtendedLegacyObjectsAndFields(testCase)
radia.UtiDelAll();
cleanup = onCleanup(@() radia.UtiDelAll());
% A repeated "off" call can leave the saved global perturbation at zero.
% Restore the documented defaults before integrating through a symmetry line.
radia.FldLenTol(1e-9, 1e-11, 1e-9);
radia.FldLenRndSw("on");

polygonObject = radia.ObjThckPgn(0, 0.2, [0, 0; 1, 0; 0, 1], "x", [0, 0, 1]);
cylinder = radia.ObjCylMag([0, 0, 0], 0.2, 0.4, 8, "z", [0, 0, 1]);
rectangularCurrent = radia.ObjRecCur([0, 0, 0], [0.4, 0.4, 0.1], [0, 0, 1]);
arcCurrent = radia.ObjArcCur([0, 0, 0], [0.2, 0.3], [0, pi / 2], ...
    0.1, 8, "m", "z", 1e5);
raceTrack = radia.ObjRaceTrk([0, 0, 0], [0.2, 0.3], [0.5, 0.5], ...
    0.1, 8, "m", "z", 1e5);
filament = radia.ObjFlmCur([0, 0, 0; 1, 0, 0; 1, 1, 0], 1);
arcPolygonMagnet = radia.ObjArcPgnMag([0, 0], "z", ...
    [0, 0; 0.2, 0; 0.2, 0.2], [0, pi / 4], 4, "nosym", [0, 0, 1]);
background = radia.ObjBckg([0, 0, 1]);

handles = [polygonObject, cylinder, rectangularCurrent, arcCurrent, raceTrack, ...
    filament, arcPolygonMagnet, background];
verifyTrue(testCase, all(handles > 0));
verifyGreaterThan(testCase, radia.ObjGeoVol(polygonObject), 0);
verifyGreaterThan(testCase, radia.ObjGeoVol(cylinder), 0);
verifyGreaterThan(testCase, radia.ObjGeoVol(arcPolygonMagnet), 0);
verifyEqual(testCase, exist("radia.ObjMltExtPgn", "file"), 0);
verifyEqual(testCase, exist("radia.ObjMltExtRtg", "file"), 0);
verifyEqual(testCase, exist("radia.ObjMltExtTri", "file"), 0);

shape = radia.FldFrcShpRtg([0, 0, 1], [1, 1]);
force = radia.FldFrc(cylinder, shape);
verifySize(testCase, force, [1, 3]);
verifyTrue(testCase, all(isfinite(force)));

points = radia.FldLst(cylinder, "b", [-1, 0, 0], [1, 0, 0], 5);
pointsWithArg = radia.FldLst(cylinder, "bz", [-1, 0, 0], [1, 0, 0], 5, "arg", 10);
verifySize(testCase, points, [5, 3]);
verifySize(testCase, pointsWithArg, [5, 2]);
verifyTrue(testCase, all(isfinite([points, pointsWithArg]), "all"));

integral = radia.FldInt(cylinder, "fin", "ibz", [-1, 0, 0], [1, 0, 0]);
centerField = radia.ObjCenFld(cylinder, "B");
verifyTrue(testCase, isscalar(integral) && isfinite(integral));
verifySize(testCase, centerField.center, [1, 3]);
verifySize(testCase, centerField.field, [1, 3]);
verifyTrue(testCase, all(isfinite([centerField.center, centerField.field])));

radia.FldCmpCrt(1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8);
radia.FldCmpPrc("PrcB->1e-8");
radia.FldLenTol(1e-12, 1e-11, 1e-12);
radia.FldLenRndSw("off");
radia.FldLenRndSw("on");
clear cleanup
end

function testExtendedMaterialsAndNonlinearSolve(testCase)
radia.UtiDelAll();
cleanup = onCleanup(@() radia.UtiDelAll());

isoFormula = radia.MatSatIsoFrm([0.01, 1e5; 0.1, 2e5; 0, 0]);
anisotropic = radia.MatSatAniso(1e-3, [0.01, 1e5, 0.01, 5e4, 0.01, 2e4, 1e-3]);
laminatedTable = radia.MatSatLamTab( ...
    [0, 0; 0.1, 1e4; 0.5, 5e4; 1, 1e5], 0.95, [0, 0, 1]);
laminatedFormula = radia.MatSatLamFrm([0.01, 1e5], [0, 0], [0, 0], ...
    0.95, [0, 0, 1]);
verifyTrue(testCase, all([isoFormula, anisotropic, laminatedTable, laminatedFormula] > 0));

mz = radia.MatMvsH(isoFormula, "mz", [0, 0, 0.1]);
verifyTrue(testCase, isscalar(mz) && isfinite(mz));

vertices = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];
object = radia.ObjTetrahedron(vertices, [0, 0, 1]);
radia.MatApl(object, isoFormula);
verifyError(testCase, @() radia.SolveNonl(object, 1e-4, 10, 0, 0), ...
    "radia:mex:Exception");
stats = radia.solveStats();
verifyTrue(testCase, isfield(stats, "t_matrix_build"));
verifyTrue(testCase, isempty(stats.t_matrix_build));
clear cleanup
end

function testHysteresisMaterialsAndState(testCase)
radia.UtiDelAll();
cleanup = onCleanup(@() radia.UtiDelAll());
tables = {{[0, 1, 2], [0, 1, 2]}};

play = radia.MatPlayHysteresis(1, 0.1, tables);
energy = radia.MatEnergyHysteresis(1, 0.1, tables);
verifyTrue(testCase, all([play, energy] > 0));

magnetization = radia.MatMvsH(play, "m", [100, 0, 0]);
verifySize(testCase, magnetization, [1, 3]);
verifyTrue(testCase, all(isfinite(magnetization)));
state = radia.MatHysSaveState(play);
verifySize(testCase, state, [16, 1]);
nuRev = radia.MatHysGetNuRev(play);
hIrr = radia.MatHysIrreversible(play, [0.1, 0, 0]);
verifyTrue(testCase, isfinite(nuRev));
verifySize(testCase, hIrr, [1, 3]);
verifyTrue(testCase, all(isfinite(hIrr)));
batchB = [0.1, 0, 0; 0.2, 0.01, 0];
batchStates = repmat(state.', size(batchB, 1), 1);
batchH = radia.MatHysForwardBatch(play, batchB, batchStates);
batchStatesOut = radia.MatHysCommitBatch(play, batchB, batchStates);
verifySize(testCase, batchH, [2, 3]);
verifySize(testCase, batchStatesOut, size(batchStates));
verifyTrue(testCase, all(isfinite([batchH, batchStatesOut]), "all"));
radia.MatHysCommitState(play);
radia.MatHysRestoreState(play, state);
clear cleanup
end

function testSolverConfigurationParity(testCase)
original = radia.getSolverConfig();
cleanup = onCleanup(@() radia.solverConfig(original));
radia.solverConfig(struct( ...
    "relax_param", 0.25, ...
    "newton_method", true, ...
    "newton_damping", true, ...
    "newton_damping_max_iter", 4, ...
    "newton_damping_min_omega", 0.02, ...
    "b_input_newton", false, ...
    "b_input_hantila", true, ...
    "hantila_alpha", 0.3, ...
    "hantila_relax", 0.4, ...
    "keep_magnetization", true));
changed = radia.getSolverConfig();
verifyEqual(testCase, changed.relax_param, 0.25, "AbsTol", 1e-14);
verifyTrue(testCase, changed.newton_method);
verifyTrue(testCase, changed.newton_damping);
verifyEqual(testCase, changed.newton_damping_max_iter, 4);
verifyEqual(testCase, changed.newton_damping_min_omega, 0.02, "AbsTol", 1e-14);
verifyFalse(testCase, changed.b_input_newton);
verifyTrue(testCase, changed.b_input_hantila);
verifyEqual(testCase, changed.hantila_alpha, 0.3, "AbsTol", 1e-14);
verifyEqual(testCase, changed.hantila_relax, 0.4, "AbsTol", 1e-14);
verifyTrue(testCase, changed.keep_magnetization);
verifyEqual(testCase, radia.getClusterStrategy(), radia.hluClusterStrategy());
clear cleanup
end

function testPlanarChargeKernels(testCase)
Xq = [0, 0; 1, 0];
Q = [1; -0.5];
P = [0.25, 0.5; 2, 0.2];
H = radia.planarChargeField(Xq, Q, P);
Az = radia.planarChargeAz(Xq, Q, P);
torque = radia.planarMaxwellTorqueCircle(Xq, Q, 3, [0.1, 0.2], [0, 0], 64);
force = radia.planarMaxwellForceCircle(Xq, Q, 3, [0.1, 0.2], [0, 0], 64);
verifySize(testCase, H, [2, 2]);
verifySize(testCase, Az, [2, 1]);
verifySize(testCase, force, [1, 2]);
verifyTrue(testCase, all(isfinite([H(:); Az(:); torque; force(:)])));
end

function testClosedFormCuboidAverages(testCase)
M = [0, 0, 1e5];
sourceMin = [0, 0, 0];
sourceMax = [1, 1, 1];
targetMin = [2, 0, 0];
targetMax = [3, 1, 1];
A = radia.averageDemagTensor(sourceMin, sourceMax, targetMin, targetMax);
B = radia.averageBInBox(M, sourceMin, sourceMax, targetMin, targetMax);
verifySize(testCase, A, [3, 3]);
verifySize(testCase, B, [1, 3]);
verifyTrue(testCase, all(isfinite([A(:); B(:)])));
end

function testEquivalenceSourceKernels(testCase)
centroids = [0, 0, 0];
normals = [0, 0, 1];
areas = 1;
Hsurf = [0, 0, 1];
obs = [0, 0, 2; 0.5, 0, 1.5];
H = radia.equivalenceSourceStaticH(centroids, normals, areas, Hsurf, obs);
verifySize(testCase, H, [2, 3]);
verifyTrue(testCase, all(isfinite(H), "all"));

EreSurf = [1, 0, 0];
EimSurf = [0.1, 0, 0];
HreSurf = [0, 1, 0];
HimSurf = [0, 0.1, 0];
[Ere, Eim, Hre, Him] = radia.equivalenceSourceHarmonic( ...
    centroids, normals, areas, EreSurf, EimSurf, HreSurf, HimSurf, obs, 1e6);
verifySize(testCase, Ere, [2, 3]);
verifySize(testCase, Eim, [2, 3]);
verifySize(testCase, Hre, [2, 3]);
verifySize(testCase, Him, [2, 3]);
verifyTrue(testCase, all(isfinite([Ere, Eim, Hre, Him]), "all"));
end

function testHACApKHLUParity(testCase)
oldTol = radia.hluGetTruncTol();
oldParallel = radia.hluGetParallel();
oldCap = radia.hluGetAccumCap();
cleanup = onCleanup(@() restoreHLU(oldTol, oldParallel, oldCap));
radia.hluSetTruncTol(oldTol);
radia.hluSetParallel(oldParallel);
radia.hluSetAccumCap(oldCap);
radia.hluSetParCutoff(0);

verifyGreaterThanOrEqual(testCase, radia.hluMaxThreads(), 1);
verifyTrue(testCase, isscalar(radia.hluClusterStrategy()));
timings = radia.hluLastTimings();
materialize = radia.hluMaterializeStats();
breakdown = radia.hluMixedBreakdown();
verifyTrue(testCase, isfield(timings, "t_decomp_sec"));
verifyTrue(testCase, isfield(materialize, "n_calls"));
verifyTrue(testCase, isfield(breakdown, "addmul"));
verifyTrue(testCase, isfield(breakdown.addmul, "rkxrk"));

errors = [radia.hluSelfTest(1, 10), ...
    radia.hluSelfTestRk(10, 2), ...
    radia.hluSelfTestAddmulRkRk(8, 8, 8, 2, 2, 2), ...
    radia.hluSelfTestRadiaExactWithMatrix(eye(162), ones(162, 1)), ...
    radia.hluSelfTestRadiaExactDiag(2), ...
    radia.hluSelfTestDepth3Asymmetric(2), ...
    radia.hluSelfTestMixedSiblingViaConversion(3), ...
    radia.hluSelfTestMixedSiblingNonUniform(3, 4, 1, 2), ...
    radia.hluSelfTestMixedSibling(3), ...
    radia.hluSelfTestRkDeep(10, 2)];
verifyTrue(testCase, all(isfinite(errors)));
verifyLessThan(testCase, errors, 1e-8);
clear cleanup
end

function testPersistentNGSolveCoefficientLinearForms(testCase)
before = radia.apiInfo();
mesh = radia.ngsolve.Mesh.create(testCase.TestData.meshPath);
scalar = radia.ngsolve.CoefficientFunction.constant(2.5);
vector = radia.ngsolve.CoefficientFunction.constant([1; 2; 3]);
h1 = radia.ngsolve.FESpace.create(mesh, "h1", 2);
hcurl = radia.ngsolve.FESpace.create(mesh, "hcurl", 1);
hdiv = radia.ngsolve.FESpace.create(mesh, "hdiv", 1);

h1Form = radia.ngsolve.LinearForm.createFromCoefficient( ...
    h1, scalar, Label="h1_coefficient_rhs");
hcurlForm = radia.ngsolve.LinearForm.createFromCoefficient( ...
    hcurl, vector, Label="hcurl_coefficient_rhs");
hdivForm = radia.ngsolve.LinearForm.createFromCoefficient( ...
    hdiv, vector, Label="hdiv_coefficient_rhs");
h1BoundaryForm = radia.ngsolve.LinearForm.createBoundaryFromCoefficient( ...
    h1, scalar, Label="h1_boundary_rhs");
hcurlBoundaryForm = radia.ngsolve.LinearForm.createBoundaryFromCoefficient( ...
    hcurl, vector, Label="hcurl_boundary_rhs");
hdivBoundaryForm = radia.ngsolve.LinearForm.createBoundaryFromCoefficient( ...
    hdiv, vector, Label="hdiv_boundary_rhs");

verifyEqual(testCase, string(h1Form.info().source), "coefficient");
verifyEqual(testCase, h1Form.info().size, h1.DofCount);
verifyEqual(testCase, hcurlForm.info().size, hcurl.DofCount);
verifyEqual(testCase, hdivForm.info().size, hdiv.DofCount);
verifyEqual(testCase, string(h1BoundaryForm.info().source), "boundary_coefficient");
verifyEqual(testCase, h1BoundaryForm.info().size, h1.DofCount);
verifyEqual(testCase, hcurlBoundaryForm.info().size, hcurl.DofCount);
verifyEqual(testCase, hdivBoundaryForm.info().size, hdiv.DofCount);
verifyTrue(testCase, all(isfinite(h1Form.vector().values())));
verifyTrue(testCase, all(isfinite(hcurlForm.vector().values())));
verifyTrue(testCase, all(isfinite(hdivForm.vector().values())));
verifyTrue(testCase, all(isfinite(h1BoundaryForm.vector().values())));
verifyTrue(testCase, all(isfinite(hcurlBoundaryForm.vector().values())));
verifyTrue(testCase, all(isfinite(hdivBoundaryForm.vector().values())));
verifyGreaterThan(testCase, norm(h1Form.vector().values()), 0);
verifyGreaterThan(testCase, norm(hcurlForm.vector().values()), 0);
verifyGreaterThan(testCase, norm(hdivForm.vector().values()), 0);
verifyGreaterThan(testCase, norm(h1BoundaryForm.vector().values()), 0);
verifyGreaterThan(testCase, norm(hcurlBoundaryForm.vector().values()), 0);
verifyGreaterThan(testCase, norm(hdivBoundaryForm.vector().values()), 0);

delete(h1Form);
delete(hcurlForm);
delete(hdivForm);
delete(h1BoundaryForm);
delete(hcurlBoundaryForm);
delete(hdivBoundaryForm);
delete(h1);
delete(hcurl);
delete(hdiv);
delete(vector);
delete(scalar);
delete(mesh);
after = radia.apiInfo();
verifyEqual(testCase, after.handle_count, before.handle_count);
end

function testAffineCellSelfEnergyShapeDerivative(testCase)
nodes = [0,0,0; 1,0,0; 0,1,0; 0,0,1];
velocities = zeros(2,4,3);
velocities(1,:,1) = nodes(:,1).';
velocities(1,:,2) = nodes(:,2).';
velocities(1,:,3) = nodes(:,3).';
velocities(2,:,1) = 0.25;
velocities(2,:,2) = -0.5;
velocities(2,:,3) = 0.75;
result = radia.topopt.affineCellSelfEnergyShapeDerivative( ...
    "tet", nodes, velocities);
verifyGreaterThan(testCase, result.value, 0);
verifySize(testCase, result.derivative, [2,1]);
verifyTrue(testCase, all(isfinite(result.derivative)));
verifyEqual(testCase, result.derivative(2), 0, "AbsTol", 2e-13);
end

function testStreamACATSVD(testCase)
leftCoordinate = linspace(-0.8, 0.9, 7).';
rightCoordinate = linspace(-0.7, 0.85, 5).';
left = [ones(7,1), leftCoordinate, leftCoordinate.^2];
right = [ones(5,1), rightCoordinate, rightCoordinate.^2];
matrix = left * right.';
[U, singularValues, V, acaRank] = radia.stream.acaTsvd( ...
    matrix, 3, MaxRank=5, AcaTolerance=1e-12);
reconstructed = U * diag(singularValues) * V.';
verifyGreaterThanOrEqual(testCase, acaRank, 3);
verifySize(testCase, U, [7,3]);
verifySize(testCase, V, [5,3]);
verifyLessThan(testCase, norm(reconstructed-matrix, "fro") / ...
    norm(matrix, "fro"), 2e-11);
end

function testProjectedAndReducedNativeMatrices(testCase)
mesh = radia.ngsolve.Mesh.create(testCase.TestData.meshPath);
space = radia.ngsolve.FESpace.create(mesh, "h1", 1);
form = radia.ngsolve.BilinearForm.create(space, "mass");
parent = form.matrix();
parentValues = full(parent.sparse());
projection = [1, 0.1i; 0.2, 1; -0.3i, 0.4; 0.5, -0.2i];
projected = radia.ngsolve.Matrix.projected(parent, projection);
projectedValues = projection' * parentValues * projection;

x = projected.vector();
xValues = [1+0.2i; -0.4+0.7i];
x.setValues(xValues);
y = projected.matvec(x);
verifyEqual(testCase, y.values(), projectedValues*xValues, ...
    "AbsTol", 2e-13);

dense = diag([2+0.5i, 3-0.25i, 4]);
reduced = radia.ngsolve.Matrix.reducedBlock( ...
    dense, {projected}, [1,2], 0.35-0.1i);
verifyEqual(testCase, reduced.termCount(), 1);
reducedX = reduced.vector();
reducedValues = [0.5+0.25i; -0.75+0.1i; 1.25-0.2i];
reducedX.setValues(reducedValues);
reducedY = reduced.matvec(reducedX);
expected = dense*reducedValues;
expected(1:2) = expected(1:2) + ...
    (0.35-0.1i)*projectedValues*reducedValues(1:2);
verifyEqual(testCase, reducedY.values(), expected, "AbsTol", 3e-13);

preconditioner = reduced.diagonalPreconditioner();
preconditioned = preconditioner.matvec(reducedX);
verifyEqual(testCase, preconditioned.values(), ...
    diag(dense).^(-1).*reducedValues, "AbsTol", 2e-14);

delete(preconditioned); delete(preconditioner); delete(reducedY);
delete(reducedX); delete(reduced); delete(y); delete(x); delete(projected);
delete(parent); delete(form); delete(space); delete(mesh);
end

function testRadiaFieldCoefficientCacheAndVoxel(testCase)
radia.UtiDelAll();
radiaCleanup = onCleanup(@() radia.UtiDelAll());
vertices = [0,0,0; 1,0,0; 0,1,0; 0,0,1];
object = radia.ObjTetrahedron(vertices, [0.2,-0.1,0.8]);
points = [0.15,0.20,0.10; 0.25,0.10,0.15];

field = radia.RadiaField(object, "b");
fieldCleanup = onCleanup(@() delete(field));
verifyEqual(testCase, field.evaluate(testCase.TestData.meshPath, points), ...
    radia.Fld(object, "b", points), "RelTol", 3e-13, "AbsTol", 2e-14);
field.prepareCache(points);
verifyEqual(testCase, field.cacheStats().size, 2);
field.evaluate(testCase.TestData.meshPath, points);
verifyGreaterThanOrEqual(testCase, field.cacheStats().hits, 2);
field.clearCache();
verifyFalse(testCase, field.cacheStats().enabled);

potential = radia.RadiaField(object, "phi");
potentialCleanup = onCleanup(@() delete(potential));
verifyEqual(testCase, potential.evaluate(testCase.TestData.meshPath, points), ...
    radia.Fld(object, "phi", points), "RelTol", 3e-13, "AbsTol", 2e-14);

mesh = radia.ngsolve.Mesh.create(testCase.TestData.meshPath);
meshCleanup = onCleanup(@() delete(mesh));
voxel = field.asVoxelCoefficient(mesh, 11);
voxelCleanup = onCleanup(@() delete(voxel));
gridPoint = [0.194,0.194,0.194];
verifyEqual(testCase, voxel.evaluate(testCase.TestData.meshPath, gridPoint), ...
    radia.Fld(object, "b", gridPoint), "RelTol", 3e-12, "AbsTol", 3e-13);
clear voxelCleanup meshCleanup potentialCleanup fieldCleanup radiaCleanup
end

function restoreHLU(tol, parallel, cap)
radia.hluSetTruncTol(tol);
radia.hluSetParallel(parallel);
radia.hluSetAccumCap(cap);
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
