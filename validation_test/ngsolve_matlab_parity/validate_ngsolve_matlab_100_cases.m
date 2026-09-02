function report = validate_ngsolve_matlab_100_cases(options)
%VALIDATE_NGSOLVE_MATLAB_100_CASES Compare native MATLAB MEX with Python.
%   Python's public NGSolve API generates the oracle and .vol meshes. MATLAB
%   then independently assembles and solves every case through checked
%   radia.ngsolve MEX handles. A JSON validation report is always written.

arguments
    options.OutputJson (1,1) string = ""
    options.WorkDirectory (1,1) string = ""
    options.KeepArtifacts (1,1) logical = false
    options.PythonExecutable (1,1) string = "python"
end

started = tic;
testDirectory = string(fileparts(mfilename("fullpath")));
validationDirectory = string(fileparts(testDirectory));
repoRoot = string(fileparts(validationDirectory));
addpath(fullfile(repoRoot, "matlab"), "-begin");
radia.setup(Force=true);

if strlength(options.OutputJson) == 0
    options.OutputJson = fullfile(testDirectory, ...
        "results_ngsolve_matlab_100_case_parity.json");
end
if strlength(options.WorkDirectory) == 0
    options.WorkDirectory = string(tempname("C:\temp"));
end
if isfolder(options.WorkDirectory)
    rmdir(options.WorkDirectory, "s");
end
mkdir(options.WorkDirectory);
workCleanup = onCleanup(@() cleanupWorkDirectory( ...
    options.WorkDirectory, options.KeepArtifacts));

generator = fullfile(testDirectory, "generate_python_oracle.py");
command = quoteCommandArgument(options.PythonExecutable) + " " + ...
    quoteCommandArgument(generator) + " " + ...
    quoteCommandArgument(options.WorkDirectory);
[status, generatorOutput] = radia.internal.runPythonProcess(command);
if status ~= 0
    error("radia:validation:NGSolvePythonOracle", ...
        "Python NGSolve oracle failed (%d): %s", status, generatorOutput);
end

manifestPath = fullfile(options.WorkDirectory, "python_oracle_manifest.json");
oraclePath = fullfile(options.WorkDirectory, "python_oracle.mat");
manifest = jsondecode(fileread(manifestPath));
oracle = load(oraclePath);
assert(manifest.case_count == 100, "The oracle must contain exactly 100 cases.");
assert(manifest.dimension_counts.three_d > 0, "The oracle must include 3D cases.");
assert(manifest.maximum_dofs < manifest.dof_limit, ...
    "The Python oracle exceeded its DoF limit.");

caseResults = repmat(emptyCaseResult(), manifest.case_count, 1);
for index = 1:manifest.case_count
    specification = manifest.cases(index);
    try
        caseResults(index) = runOneCase(specification, oracle, manifest.dof_limit);
    catch exception
        failed = emptyCaseResult();
        failed.case_id = string(specification.case_id);
        failed.mesh_id = string(specification.mesh_id);
        failed.dimension = specification.dimension;
        failed.space = string(specification.space);
        failed.order = specification.order;
        failed.form = string(specification.form);
        failed.solve = logical(specification.solve);
        failed.failure = string(getReport(exception, "extended", ...
            "hyperlinks", "off"));
        caseResults(index) = failed;
    end
    if mod(index, 10) == 0
        fprintf("NGSolve MATLAB parity: %d/%d cases complete\n", ...
            index, manifest.case_count);
    end
end

apiInfo = radia.apiInfo();
passed = [caseResults.passed];
report = struct( ...
    "schema", "radia.ngsolve-python-matlab-100-case-validation.v1", ...
    "generated_utc", string(datetime("now", "TimeZone", "UTC", ...
        "Format", "yyyy-MM-dd'T'HH:mm:ss'Z'")), ...
    "all_passed", all(passed), ...
    "case_count", numel(caseResults), ...
    "passed_count", nnz(passed), ...
    "failed_count", nnz(~passed), ...
    "solve_case_count", nnz([caseResults.solve]), ...
    "two_d_case_count", nnz([caseResults.dimension] == 2), ...
    "three_d_case_count", nnz([caseResults.dimension] == 3), ...
    "dof_limit", manifest.dof_limit, ...
    "maximum_dofs", max([caseResults.dofs]), ...
    "maximum_matrix_relative_error", ...
        max([caseResults.matrix_relative_error]), ...
    "maximum_matvec_relative_error", ...
        max([caseResults.matvec_relative_error]), ...
    "maximum_solution_relative_error", ...
        max([caseResults.solution_relative_error]), ...
    "maximum_free_residual_relative", ...
        max([caseResults.free_residual_relative]), ...
    "native_handle_count_after", apiInfo.handle_count, ...
    "matlab_version", string(version), ...
    "matlab_release", string(version("-release")), ...
    "python_version", string(manifest.python_version), ...
    "numpy_version", string(manifest.numpy_version), ...
    "scipy_version", string(manifest.scipy_version), ...
    "ngsolve_version", string(manifest.ngsolve_version), ...
    "netgen_mesher_version", string(manifest.netgen_mesher_version), ...
    "python_total_duration_s", manifest.python_total_duration_s, ...
    "matlab_total_duration_s", toc(started), ...
    "case_results", caseResults);

writeJson(options.OutputJson, report);
fprintf("NGSolve MATLAB parity: %d/%d passed; report: %s\n", ...
    report.passed_count, report.case_count, options.OutputJson);
if ~report.all_passed
    first = caseResults(find(~passed, 1));
    error("radia:validation:NGSolveMatlabParity", ...
        "NGSolve MATLAB parity failed at %s: %s", ...
        first.case_id, first.failure);
end
assert(report.maximum_dofs < report.dof_limit, ...
    "A realized MATLAB case exceeded the DoF limit.");
assert(report.native_handle_count_after == 0, ...
    "Native NGSolve handles leaked after validation.");
clear workCleanup
end


function result = runOneCase(specification, oracle, dofLimit)
started = tic;
result = emptyCaseResult();
result.case_id = string(specification.case_id);
result.mesh_id = string(specification.mesh_id);
result.dimension = specification.dimension;
result.space = string(specification.space);
result.order = specification.order;
result.form = string(specification.form);
result.weight = specification.weight;
result.solve = logical(specification.solve);

mesh = radia.ngsolve.Mesh.create(string(specification.mesh_path));
meshCleanup = onCleanup(@() delete(mesh));
space = radia.ngsolve.FESpace.create(mesh, result.space, result.order, ...
    Dirichlet=string(specification.dirichlet));
spaceCleanup = onCleanup(@() delete(space));
coefficient = radia.ngsolve.CoefficientFunction.constant(result.weight);
coefficientCleanup = onCleanup(@() delete(coefficient));
form = radia.ngsolve.BilinearForm.createFromCoefficient( ...
    space, result.form, coefficient, Label=result.case_id);
formCleanup = onCleanup(@() delete(form));
matrix = form.matrix();
matrixCleanup = onCleanup(@() delete(matrix));

key = string(specification.oracle_key);
expectedDofs = double(oracle.(char(key + "_ndof")));
expectedNnz = double(oracle.(char(key + "_nnz")));
expectedFree = logical(oracle.(char(key + "_free_dofs")));
expectedFree = expectedFree(:);
actualFree = space.freeDofs();
actualMatrix = matrix.sparse();
expectedMatrix = oracle.(char(key + "_matrix"));
matrixDelta = actualMatrix - expectedMatrix;

result.dofs = space.DofCount;
result.free_dofs = nnz(actualFree);
result.nonzeros = nnz(abs(actualMatrix) > 5e-13);
result.matrix_absolute_error = norm(matrixDelta, "fro");
result.matrix_relative_error = relativeError(actualMatrix, expectedMatrix);
result.dof_match = result.dofs == expectedDofs && result.dofs < dofLimit;
result.free_dof_match = isequal(actualFree, expectedFree);
result.nonzero_match = result.nonzeros == nnz(abs(expectedMatrix) > 5e-13) ...
    && expectedNnz >= result.nonzeros;

input = matrix.vector();
inputCleanup = onCleanup(@() delete(input));
input.setValues(oracle.(char(key + "_input")));
matvec = matrix.matvec(input);
matvecCleanup = onCleanup(@() delete(matvec));
expectedMatvec = oracle.(char(key + "_matvec"));
result.matvec_relative_error = relativeError( ...
    matvec.values(), expectedMatvec);

if result.solve
    rhs = matrix.vector();
    rhsCleanup = onCleanup(@() delete(rhs));
    rhsValues = oracle.(char(key + "_rhs"));
    rhs.setValues(rhsValues);
    inverse = matrix.inverse();
    inverseCleanup = onCleanup(@() delete(inverse));
    solution = matrix.vector();
    solutionCleanup = onCleanup(@() delete(solution));
    inverse.matvecInto(rhs, solution);
    solutionValues = solution.values();
    expectedSolution = oracle.(char(key + "_solution"));
    result.solution_relative_error = relativeError( ...
        solutionValues, expectedSolution);
    applied = matrix.matvec(solution);
    appliedCleanup = onCleanup(@() delete(applied));
    residual = applied.values() - rhsValues;
    result.free_residual_relative = norm(residual(actualFree)) / ...
        max(norm(rhsValues(actualFree)), 1.0);
    clear appliedCleanup solutionCleanup inverseCleanup rhsCleanup
end

result.duration_s = toc(started);
result.passed = result.dof_match && result.free_dof_match && ...
    result.nonzero_match && result.matrix_relative_error <= 5e-12 && ...
    result.matvec_relative_error <= 2e-11 && ...
    result.solution_relative_error <= 5e-10 && ...
    result.free_residual_relative <= 2e-10;
if ~result.passed
    result.failure = "numerical tolerance or structural parity check failed";
end
clear matvecCleanup inputCleanup matrixCleanup formCleanup
clear coefficientCleanup spaceCleanup meshCleanup
end


function result = emptyCaseResult()
result = struct( ...
    "case_id", "", "mesh_id", "", "dimension", 0, "space", "", ...
    "order", 0, "form", "", "weight", 0, "solve", false, ...
    "dofs", 0, "free_dofs", 0, "nonzeros", 0, ...
    "dof_match", false, "free_dof_match", false, ...
    "nonzero_match", false, "matrix_absolute_error", 0, ...
    "matrix_relative_error", 0, "matvec_relative_error", 0, ...
    "solution_relative_error", 0, "free_residual_relative", 0, ...
    "duration_s", 0, "passed", false, "failure", "");
end


function value = relativeError(actual, expected)
value = norm(actual - expected, "fro") / max(norm(expected, "fro"), 1.0);
end


function value = quoteCommandArgument(value)
value = '"' + replace(string(value), '"', '""') + '"';
end


function writeJson(path, value)
directory = string(fileparts(path));
if strlength(directory) > 0 && ~isfolder(directory)
    mkdir(directory);
end
file = fopen(path, "w", "n", "UTF-8");
if file < 0
    error("radia:validation:JsonWrite", "Could not write %s", path);
end
cleanup = onCleanup(@() fclose(file));
fprintf(file, "%s\n", jsonencode(value, PrettyPrint=true));
clear cleanup
end


function cleanupWorkDirectory(path, keepArtifacts)
if ~keepArtifacts && isfolder(path)
    rmdir(path, "s");
end
end

