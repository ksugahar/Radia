function tests = test_ngsolve_mex_pybind_parity
%TEST_NGSOLVE_MEX_PYBIND_PARITY Compare MATLAB MEX NGSolve with Python.
%   The expected values are created by the public Python NGSolve boundary on
%   the same .vol mesh.  They are never handwritten MATLAB goldens.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
addpath(fullfile(repoRoot, "matlab"));
% Exercise a just-linked gateway before MATLAB releases the staged MEX. The
% production lane leaves this unset and resolves matlab/radia_mex.mexw64.
mexTestDir = string(getenv("RADIA_MEX_TEST_DIR"));
if strlength(mexTestDir) > 0
    assert(isfolder(mexTestDir), ...
        "RADIA_MEX_TEST_DIR must name a directory containing radia_mex.");
    addpath(mexTestDir, "-begin");
    clear radia_mex
end
testCase.TestData.SetupInfo = radia.setup(Force=true);
if strlength(mexTestDir) > 0
    addpath(mexTestDir, "-begin");
    clear radia_mex
    assert(contains(which("radia_mex"), mexTestDir), ...
        "RADIA_MEX_TEST_DIR did not take precedence for radia_mex.");
    testCase.TestData.MexTestDir = mexTestDir;
end
testCase.TestData.MeshPath = writeUnitTetra();
[path, reference] = loadPythonReference(repoRoot, testCase.TestData.MeshPath);
testCase.TestData.PythonReferencePath = path;
testCase.TestData.PythonReference = reference;
end

function teardownOnce(testCase)
for path = [testCase.TestData.PythonReferencePath; testCase.TestData.MeshPath].'
    if isfile(path)
        delete(path);
    end
end
end

function testCoordinateCoefficientComponents(testCase)
coordinates = radia.ngsolve.CoefficientFunction.coordinates(3);
one = radia.ngsolve.CoefficientFunction.constant(1.0);
x = coordinates.component(1);
y = coordinates.component(2);
realScaled = 0.5 * x;
complexScaled = complex(0.5, 0.1) * x;
verifyFalse(testCase, realScaled.info().is_complex);
verifyTrue(testCase, complexScaled.info().is_complex);
field = x * (one - x) * y;
points = [0.10, 0.20, 0.30; 0.20, 0.25, 0.15];
verifyEqual(testCase, ...
    field.evaluate(testCase.TestData.MeshPath, points), ...
    points(:,1) .* (1.0 - points(:,1)) .* points(:,2), ...
    AbsTol=2e-14);
delete(field);
delete(complexScaled);
delete(realScaled);
delete(y);
delete(x);
delete(one);
delete(coordinates);
end

function testMeshAndFiniteElementSpaceMatchPython(testCase)
reference = testCase.TestData.PythonReference;
mesh = radia.ngsolve.Mesh.create(testCase.TestData.MeshPath);
cleanup = onCleanup(@() delete(mesh));
info = mesh.info();
verifyEqual(testCase, info.dimension, double(reference.mesh_dimension));
verifyEqual(testCase, info.vertices, double(reference.mesh_vertices));
verifyEqual(testCase, info.elements, double(reference.mesh_elements));

specifications = { ...
    "h1", 2, "h1_mass_ndof"; ...
    "hcurl", 2, "hcurl_mass_ndof"; ...
    "hdiv", 2, "hdiv_mass_ndof" ...
    };
for row = 1:size(specifications, 1)
    space = radia.ngsolve.FESpace.create(mesh, specifications{row, 1}, ...
        specifications{row, 2});
    spaceCleanup = onCleanup(@() delete(space));
    verifyEqual(testCase, space.DofCount, reference.(specifications{row, 3}));
    clear spaceCleanup
end
clear cleanup
end

function testCoefficientAndGridFunctionEvaluationMatchPython(testCase)
reference = testCase.TestData.PythonReference;
points = reference.points;
coordinates = radia.ngsolve.CoefficientFunction.coordinates(3);
coordinateCleanup = onCleanup(@() delete(coordinates));
constant = radia.ngsolve.CoefficientFunction.constant([1.25; -0.5; 0.75]);
constantCleanup = onCleanup(@() delete(constant));
verifyEqual(testCase, coordinates.evaluate(testCase.TestData.MeshPath, points), ...
    reference.coordinates, AbsTol=2e-14);
verifyEqual(testCase, constant.evaluate(testCase.TestData.MeshPath, points), ...
    reference.constant_vector, AbsTol=2e-14);

mesh = radia.ngsolve.Mesh.create(testCase.TestData.MeshPath);
meshCleanup = onCleanup(@() delete(mesh));
space = radia.ngsolve.FESpace.create(mesh, "h1", 2);
spaceCleanup = onCleanup(@() delete(space));
grid = radia.ngsolve.GridFunction.fromFESpace(space, Name="pybind_oracle");
gridCleanup = onCleanup(@() delete(grid));
grid.setVector(reference.h1_grid_values);
gridCoefficient = grid.asCoefficient();
coefficientCleanup = onCleanup(@() delete(gridCoefficient));
verifyEqual(testCase, grid.vector(), reference.h1_grid_values, AbsTol=0);
verifyEqual(testCase, gridCoefficient.evaluate(testCase.TestData.MeshPath, points), ...
    reference.h1_grid_evaluation, RelTol=3e-13, AbsTol=3e-14);
clear coefficientCleanup gridCleanup spaceCleanup meshCleanup
clear constantCleanup coordinateCleanup
end

function testAssembledSparseMatricesMatchPython(testCase)
reference = testCase.TestData.PythonReference;
cases = { ...
    "h1_mass", "h1", 2, "mass"; ...
    "h1_stiffness", "h1", 2, "stiffness"; ...
    "hcurl_mass", "hcurl", 2, "mass"; ...
    "hcurl_curlcurl", "hcurl", 1, "curlcurl"; ...
    "hdiv_mass", "hdiv", 2, "mass"; ...
    "hdiv_divdiv", "hdiv", 1, "divdiv" ...
    };
for row = 1:size(cases, 1)
    prefix = cases{row, 1};
    actual = radia.ngsolve.matrix_dump(testCase.TestData.MeshPath, ...
        cases{row, 2}, cases{row, 3}, cases{row, 4});
    expected = reference.(prefix + "_matrix");
    verifyEqual(testCase, size(actual), size(expected));
    verifyEqual(testCase, nnz(abs(actual) > 4e-14), ...
        nnz(abs(expected) > 4e-14));
    verifyEqual(testCase, actual, expected, RelTol=4e-13, AbsTol=4e-14);
end
end

function testPersistentFormsNativeMatvecAndSolverMatchPython(testCase)
reference = testCase.TestData.PythonReference;
mesh = radia.ngsolve.Mesh.create(testCase.TestData.MeshPath);
meshCleanup = onCleanup(@() delete(mesh));
space = radia.ngsolve.FESpace.create(mesh, "h1", 2);
spaceCleanup = onCleanup(@() delete(space));
form = radia.ngsolve.BilinearForm.create(space, "mass");
formCleanup = onCleanup(@() delete(form));
matrix = form.matrix();
matrixCleanup = onCleanup(@() delete(matrix));
verifyEqual(testCase, matrix.sparse(), reference.h1_mass_matrix, ...
    RelTol=4e-13, AbsTol=4e-14);

rhs = matrix.vector();
rhsCleanup = onCleanup(@() delete(rhs));
rhs.setValues(reference.h1_solver_rhs);
matvec = matrix.matvec(rhs);
matvecCleanup = onCleanup(@() delete(matvec));
verifyEqual(testCase, matvec.values(), ...
    reference.h1_mass_matrix * reference.h1_solver_rhs, ...
    RelTol=4e-13, AbsTol=4e-14);

solver = radia.ngsolve.Solver.create(matrix, "cg", ...
    Tolerance=1e-12, MaxSteps=100);
solverCleanup = onCleanup(@() delete(solver));
solution = solver.solve(rhs);
solutionCleanup = onCleanup(@() delete(solution));
verifyEqual(testCase, solution.values(), reference.h1_solver_solution, ...
    RelTol=2e-10, AbsTol=2e-11);
residual = matrix.matvec(solution);
residualCleanup = onCleanup(@() delete(residual));
verifyEqual(testCase, residual.values(), reference.h1_solver_rhs, ...
    RelTol=2e-10, AbsTol=2e-11);
verifyLessThan(testCase, norm(reference.h1_solver_residual), 2e-12);

constant = radia.ngsolve.CoefficientFunction.constant(2.5);
constantCleanup = onCleanup(@() delete(constant));
linear = radia.ngsolve.LinearForm.createFromCoefficient(space, constant);
linearCleanup = onCleanup(@() delete(linear));
linearVector = linear.vector();
linearCleanup2 = onCleanup(@() delete(linearVector));
verifyEqual(testCase, linearVector.values(), reference.h1_constant_rhs, ...
    RelTol=4e-13, AbsTol=4e-14);
clear linearCleanup2 linearCleanup constantCleanup residualCleanup solutionCleanup
clear solverCleanup matvecCleanup rhsCleanup matrixCleanup formCleanup
clear spaceCleanup meshCleanup
end

function testDirichletPoissonSolveUsesNGSolveFreeDofs(testCase)
mesh = radia.ngsolve.Mesh.create(testCase.TestData.MeshPath);
meshCleanup = onCleanup(@() delete(mesh));
space = radia.ngsolve.FESpace.create(mesh, "h1", 4, Dirichlet="outer");
spaceCleanup = onCleanup(@() delete(space));
free = space.freeDofs();
verifyEqual(testCase, space.Dirichlet, "outer");
verifyEqual(testCase, nnz(free), space.FreeDofCount);
verifyGreaterThan(testCase, nnz(free), 0);
verifyGreaterThan(testCase, nnz(~free), 0);

form = radia.ngsolve.BilinearForm.create(space, "stiffness");
formCleanup = onCleanup(@() delete(form));
matrix = form.matrix();
matrixCleanup = onCleanup(@() delete(matrix));
loadForm = radia.ngsolve.LinearForm.create(space, "constant", Value=1.0);
loadCleanup = onCleanup(@() delete(loadForm));
rhs = loadForm.vector();
rhsCleanup = onCleanup(@() delete(rhs));
inverse = matrix.inverse();
inverseCleanup = onCleanup(@() delete(inverse));
grid = radia.ngsolve.GridFunction.fromFESpace(space, Name="poisson_solution");
gridCleanup = onCleanup(@() delete(grid));
solution = grid.vectorHandle();
solutionCleanup = onCleanup(@() delete(solution));
inverse.matvecInto(rhs, solution);
applied = matrix.matvec(solution);
appliedCleanup = onCleanup(@() delete(applied));

values = solution.values();
residual = applied.values() - rhs.values();
verifyEqual(testCase, values(~free), zeros(nnz(~free), 1), AbsTol=1e-14);
verifyLessThan(testCase, norm(residual(free)), 1e-11);
verifyGreaterThan(testCase, norm(values(free)), 0);

clear appliedCleanup solutionCleanup gridCleanup inverseCleanup rhsCleanup loadCleanup
clear matrixCleanup formCleanup spaceCleanup meshCleanup
end

function [path, reference] = loadPythonReference(repoRoot, meshPath)
path = string(tempname("C:\temp")) + ".mat";
pythonExecutable = string(getenv("RADIA_PYTHON_EXECUTABLE"));
if strlength(pythonExecutable) == 0
    pythonExecutable = "python";
end
helper = fullfile(repoRoot, "tests", "matlab", ...
    "ngsolve_mex_python_reference.py");
command = quoteCommandArgument(pythonExecutable) + " " + ...
    quoteCommandArgument(helper) + " " + quoteCommandArgument(meshPath) + ...
    " " + quoteCommandArgument(path);
[status, output] = radia.internal.runPythonProcess(command);
if status ~= 0
    error("radia:test:PythonReference", ...
        "Python NGSolve reference failed (%d): %s", status, output);
end
reference = load(path);
end

function value = quoteCommandArgument(value)
value = '"' + replace(string(value), '"', '""') + '"';
end

function path = writeUnitTetra()
path = string(tempname("C:\temp")) + ".vol";
lines = [ ...
    "mesh3d"; "dimension"; "3"; "geomtype"; "0"; ...
    "facedescriptors"; "1"; "1 1 0 1 1"; ...
    "surfaceelements"; "4"; ...
    "1 1 1 0 3 1 2 3"; "1 1 1 0 3 1 4 2"; ...
    "1 1 1 0 3 2 4 3"; "1 1 1 0 3 3 4 1"; ...
    "volumeelements"; "1"; "1 4 1 2 3 4"; ...
    "points"; "4"; "0 0 0"; "1 0 0"; "0 1 0"; "0 0 1"; ...
    "pointelements"; "0"; "materials"; "1"; "1 air"; ...
    "bcnames"; "1"; "1 outer"; "endmesh" ...
    ];
file = fopen(path, "w");
if file < 0
    error("radia:test:MeshWrite", "Could not create %s", path);
end
cleanup = onCleanup(@() fclose(file));
fprintf(file, "%s\n", lines);
clear cleanup
end
