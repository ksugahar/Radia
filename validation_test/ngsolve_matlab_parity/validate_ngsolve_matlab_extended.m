function report = validate_ngsolve_matlab_extended(options)
%VALIDATE_NGSOLVE_MATLAB_EXTENDED Run breadth, scale, and analytic gates.
%   Breadth compares 500 small Python/MATLAB NGSolve assemblies. Scale
%   compares 20 matrix-free applications above 10,000 DoFs. Manufactured
%   solutions compare MATLAB/MEX finite-element solves with analytic fields.

arguments
    options.OutputDirectory (1,1) string = ""
    options.WorkDirectory (1,1) string = ""
    options.KeepArtifacts (1,1) logical = false
    options.ReuseOracle (1,1) logical = false
    options.PythonExecutable (1,1) string = "python"
    options.RunBreadth (1,1) logical = true
    options.RunScale (1,1) logical = true
    options.RunManufactured (1,1) logical = true
end

started = tic;
testDirectory = string(fileparts(mfilename("fullpath")));
validationDirectory = string(fileparts(testDirectory));
repoRoot = string(fileparts(validationDirectory));
addpath(fullfile(repoRoot, "matlab"), "-begin");
radia.setup(Force=true);
if strlength(options.OutputDirectory) == 0
    options.OutputDirectory = testDirectory;
end
if strlength(options.WorkDirectory) == 0
    options.WorkDirectory = string(tempname("C:\temp"));
end
if isfolder(options.WorkDirectory) && ~options.ReuseOracle
    rmdir(options.WorkDirectory, "s");
end
if ~isfolder(options.WorkDirectory)
    mkdir(options.WorkDirectory);
end
cleanup = onCleanup(@() cleanupWorkDirectory( ...
    options.WorkDirectory, options.KeepArtifacts));

generator = fullfile(testDirectory, "generate_extended_oracle.py");
command = quoteCommandArgument(options.PythonExecutable) + " " + ...
    quoteCommandArgument(generator) + " " + ...
    quoteCommandArgument(options.WorkDirectory);
[status, generatorOutput] = radia.internal.runPythonProcess(command);
if status ~= 0
    error("radia:validation:NGSolveExtendedOracle", ...
        "Extended Python oracle failed (%d): %s", status, generatorOutput);
end
manifest = jsondecode(fileread(fullfile( ...
    options.WorkDirectory, "python_extended_manifest.json")));
oracle = load(fullfile(options.WorkDirectory, "python_extended_oracle.mat"));
assert(manifest.breadth_case_count == 500);
assert(manifest.scale_case_count == 20);

breadthPath = fullfile(options.OutputDirectory, ...
    "results_ngsolve_matlab_breadth_500.json");
if options.RunBreadth
    breadth = runBreadth(manifest, oracle);
    writeJson(breadthPath, breadth);
else
    breadth = readJson(breadthPath);
end
scalePath = fullfile(options.OutputDirectory, ...
    "results_ngsolve_matlab_scale_20.json");
if options.RunScale
    scale = runScale(manifest, oracle);
    writeJson(scalePath, scale);
else
    scale = readJson(scalePath);
end
manufacturedPath = fullfile(options.OutputDirectory, ...
    "results_ngsolve_matlab_manufactured.json");
if options.RunManufactured
    manufactured = runManufactured(manifest);
    writeJson(manufacturedPath, manufactured);
else
    manufactured = readJson(manufacturedPath);
end

apiInfo = radia.apiInfo();
mexPath = string(which("radia_mex"));
tierDuration = breadth.duration_s + scale.duration_s + manufactured.duration_s;
report = struct( ...
    "schema", "radia.ngsolve-python-matlab-extended-validation.v1", ...
    "generated_utc", utcNow(), ...
    "all_passed", breadth.all_passed && scale.all_passed && ...
        manufactured.all_passed && apiInfo.handle_count == 0, ...
    "breadth", summaryWithoutCases(breadth), ...
    "scale", summaryWithoutCases(scale), ...
    "manufactured", summaryWithoutCases(manufactured), ...
    "executed_tiers", struct("breadth", options.RunBreadth, ...
        "scale", options.RunScale, "manufactured", options.RunManufactured), ...
    "native_handle_count_after", apiInfo.handle_count, ...
    "host_name", hostName(), ...
    "computer_architecture", string(computer("arch")), ...
    "mex_path", mexPath, ...
    "mex_sha256", fileSha256(mexPath), ...
    "matlab_version", string(version), ...
    "matlab_release", string(version("-release")), ...
    "ngsolve_version", string(manifest.ngsolve_version), ...
    "python_version", string(manifest.python_version), ...
    "total_duration_s", tierDuration, ...
    "invocation_duration_s", toc(started));
writeJson(fullfile(options.OutputDirectory, ...
    "results_ngsolve_matlab_extended.json"), report);
if ~report.all_passed
    error("radia:validation:NGSolveExtendedParity", ...
        "An extended NGSolve MATLAB validation tier failed.");
end
clear cleanup
end


function report = runBreadth(manifest, oracle)
started = tic;
count = manifest.breadth_case_count;
results = repmat(emptyBreadthResult(), count, 1);
mesh = [];
meshId = "";
for index = 1:count
    specification = manifest.breadth_cases(index);
    if string(specification.mesh_id) ~= meshId
        if ~isempty(mesh)
            delete(mesh);
        end
        mesh = radia.ngsolve.Mesh.create(string(specification.mesh_path));
        meshId = string(specification.mesh_id);
    end
    try
        results(index) = runBreadthCase(mesh, specification, oracle, ...
            manifest.dof_limit);
    catch exception
        results(index) = failedBreadthResult(specification, exception);
    end
    if mod(index, 25) == 0
        fprintf("NGSolve breadth: %d/%d cases complete\n", index, count);
    end
end
if ~isempty(mesh)
    delete(mesh);
end
passed = [results.passed];
report = struct( ...
    "schema", "radia.ngsolve-python-matlab-breadth-500.v1", ...
    "generated_utc", utcNow(), ...
    "all_passed", all(passed), ...
    "case_count", count, ...
    "passed_count", nnz(passed), ...
    "failed_count", nnz(~passed), ...
    "mesh_count", numel(manifest.small_meshes), ...
    "two_d_case_count", nnz([results.dimension] == 2), ...
    "three_d_case_count", nnz([results.dimension] == 3), ...
    "complex_case_count", nnz([results.is_complex]), ...
    "spatial_coefficient_case_count", ...
        nnz([results.coefficient_kind] == "affine_x"), ...
    "boundary_case_count", nnz([results.boundary]), ...
    "maximum_dofs", max([results.dofs]), ...
    "dof_limit", manifest.dof_limit, ...
    "maximum_matrix_relative_error", max([results.matrix_relative_error]), ...
    "maximum_matvec_relative_error", max([results.matvec_relative_error]), ...
    "maximum_solution_relative_error", max([results.solution_relative_error]), ...
    "duration_s", toc(started), ...
    "case_results", results);
end


function result = runBreadthCase(mesh, specification, oracle, dofLimit)
started = tic;
result = emptyBreadthResult();
result.case_id = string(specification.case_id);
result.mesh_id = string(specification.mesh_id);
result.dimension = specification.dimension;
result.space = string(specification.space);
result.order = specification.order;
result.form = string(specification.form);
result.coefficient_kind = string(specification.coefficient_kind);
result.boundary = logical(specification.boundary);
result.is_complex = specification.weight_imag ~= 0;
result.solve = logical(specification.solve);

space = makeSpace(mesh, specification);
spaceCleanup = onCleanup(@() delete(space));
coefficient = makeCoefficient(specification, result.dimension);
coefficientCleanup = onCleanup(@() delete(coefficient));
form = makeForm(space, coefficient, specification);
formCleanup = onCleanup(@() delete(form));
matrix = form.matrix();
matrixCleanup = onCleanup(@() delete(matrix));
key = string(specification.oracle_key);
actualMatrix = matrix.sparse();
expectedMatrix = oracle.(char(key + "_matrix"));
actualFree = space.freeDofs();
expectedFree = logical(oracle.(char(key + "_free_dofs")));
expectedFree = expectedFree(:);

result.dofs = space.DofCount;
result.nonzeros = nnz(abs(actualMatrix) > 5e-13);
result.dof_match = result.dofs == double(oracle.(char(key + "_ndof"))) ...
    && result.dofs < dofLimit;
result.free_dof_match = isequal(actualFree, expectedFree);
result.nonzero_match = result.nonzeros == ...
    nnz(abs(expectedMatrix) > 5e-13);
result.matrix_relative_error = relativeError(actualMatrix, expectedMatrix);

input = matrix.vector();
inputCleanup = onCleanup(@() delete(input));
input.setValues(oracle.(char(key + "_input")));
output = matrix.matvec(input);
outputCleanup = onCleanup(@() delete(output));
result.matvec_relative_error = relativeError( ...
    output.values(), oracle.(char(key + "_matvec")));
if result.solve
    rhs = matrix.vector();
    rhsCleanup = onCleanup(@() delete(rhs));
    rhs.setValues(oracle.(char(key + "_rhs")));
    inverse = matrix.inverse();
    inverseCleanup = onCleanup(@() delete(inverse));
    solution = inverse.matvec(rhs);
    solutionCleanup = onCleanup(@() delete(solution));
    result.solution_relative_error = relativeError( ...
        solution.values(), oracle.(char(key + "_solution")));
end
result.duration_s = toc(started);
result.passed = result.dof_match && result.free_dof_match && ...
    result.nonzero_match && result.matrix_relative_error < 2e-11 && ...
    result.matvec_relative_error < 2e-10 && ...
    (~result.solve || result.solution_relative_error < 2e-8);
end


function report = runScale(manifest, oracle)
started = tic;
count = manifest.scale_case_count;
results = repmat(emptyScaleResult(), count, 1);
mesh = [];
meshId = "";
for index = 1:count
    specification = manifest.scale_cases(index);
    if string(specification.mesh_id) ~= meshId
        if ~isempty(mesh)
            delete(mesh);
        end
        mesh = radia.ngsolve.Mesh.create(string(specification.mesh_path));
        meshId = string(specification.mesh_id);
    end
    try
        results(index) = runScaleCase(mesh, specification, oracle, ...
            manifest.scale_dof_floor, manifest.dof_limit);
    catch exception
        results(index) = failedScaleResult(specification, exception);
    end
    fprintf("NGSolve scale: %d/%d cases complete\n", index, count);
end
if ~isempty(mesh)
    delete(mesh);
end
passed = [results.passed];
report = struct( ...
    "schema", "radia.ngsolve-python-matlab-scale-20.v1", ...
    "generated_utc", utcNow(), ...
    "all_passed", all(passed), ...
    "case_count", count, ...
    "passed_count", nnz(passed), ...
    "failed_count", nnz(~passed), ...
    "minimum_dofs", min([results.dofs]), ...
    "maximum_dofs", max([results.dofs]), ...
    "dof_floor", manifest.scale_dof_floor, ...
    "dof_limit", manifest.dof_limit, ...
    "maximum_matvec_relative_error", max([results.matvec_relative_error]), ...
    "maximum_energy_relative_error", max([results.energy_relative_error]), ...
    "duration_s", toc(started), ...
    "case_results", results);
end


function result = runScaleCase(mesh, specification, oracle, dofFloor, dofLimit)
started = tic;
result = emptyScaleResult();
result.case_id = string(specification.case_id);
result.mesh_id = string(specification.mesh_id);
result.dimension = specification.dimension;
result.space = string(specification.space);
result.order = specification.order;
result.form = string(specification.form);
space = makeSpace(mesh, specification);
spaceCleanup = onCleanup(@() delete(space));
coefficient = makeCoefficient(specification, result.dimension);
coefficientCleanup = onCleanup(@() delete(coefficient));
form = makeForm(space, coefficient, specification);
formCleanup = onCleanup(@() delete(form));
matrix = form.matrix();
matrixCleanup = onCleanup(@() delete(matrix));
key = string(specification.oracle_key);
input = matrix.vector();
inputCleanup = onCleanup(@() delete(input));
input.setValues(oracle.(char(key + "_input")));
output = matrix.vector();
outputCleanup = onCleanup(@() delete(output));
matrix.matvecInto(input, output);
actual = output.values();
expected = oracle.(char(key + "_matvec"));
actualEnergy = real(input.dot(output));
expectedEnergy = double(oracle.(char(key + "_energy")));
result.dofs = space.DofCount;
result.dof_match = result.dofs == double(oracle.(char(key + "_ndof"))) ...
    && result.dofs >= dofFloor && result.dofs < dofLimit;
result.matvec_relative_error = relativeError(actual, expected);
result.energy_relative_error = abs(actualEnergy - expectedEnergy) / ...
    max(abs(expectedEnergy), 1.0);
result.input_norm_relative_error = abs(input.norm() - ...
    double(oracle.(char(key + "_input_norm")))) / ...
    max(double(oracle.(char(key + "_input_norm"))), 1.0);
result.matvec_norm_relative_error = abs(output.norm() - ...
    double(oracle.(char(key + "_matvec_norm")))) / ...
    max(double(oracle.(char(key + "_matvec_norm"))), 1.0);
result.duration_s = toc(started);
result.passed = result.dof_match && result.matvec_relative_error < 5e-10 ...
    && result.energy_relative_error < 5e-10 ...
    && result.input_norm_relative_error < 5e-12 ...
    && result.matvec_norm_relative_error < 5e-10;
end


function report = runManufactured(manifest)
started = tic;
families = ["tri_meshed", "quad_single", "tet_meshed", ...
    "hex_single", "wedge_single"];
orders = {[1, 2, 4], [1, 2, 4], [2, 4, 6], [2, 4, 6], [2, 4, 6]};
results = repmat(emptyManufacturedResult(), 15, 1);
row = 0;
for familyIndex = 1:numel(families)
    meshRow = manifest.small_meshes( ...
        find(string({manifest.small_meshes.id}) == families(familyIndex), 1));
    mesh = radia.ngsolve.Mesh.create(string(meshRow.path));
    meshCleanup = onCleanup(@() delete(mesh));
    for order = orders{familyIndex}
        row = row + 1;
        try
            results(row) = runManufacturedCase( ...
                mesh, families(familyIndex), meshRow.dimension, order, ...
                string(meshRow.path));
        catch exception
            failed = emptyManufacturedResult();
            failed.case_id = families(familyIndex) + "_p" + order;
            failed.mesh_id = families(familyIndex);
            failed.dimension = meshRow.dimension;
            failed.order = order;
            failed.failure = string(getReport(exception, "extended", ...
                "hyperlinks", "off"));
            results(row) = failed;
        end
        fprintf("NGSolve manufactured: %d/15 cases complete\n", row);
    end
    clear meshCleanup
end
convergence = repmat(struct("mesh_id", "", "first_error", 0, ...
    "final_error", 0, "monotone", false, "final_accurate", false), ...
    numel(families), 1);
for index = 1:numel(families)
    selected = results([results.mesh_id] == families(index));
    errors = [selected.sample_max_absolute_error];
    convergence(index).mesh_id = families(index);
    convergence(index).first_error = errors(1);
    convergence(index).final_error = errors(end);
    convergence(index).monotone = all(diff(errors) <= 1e-10);
    convergence(index).final_accurate = errors(end) < 5e-9;
end
passed = [results.passed];
report = struct( ...
    "schema", "radia.ngsolve-matlab-manufactured-solution.v1", ...
    "generated_utc", utcNow(), ...
    "all_passed", all(passed) && all([convergence.monotone]) && ...
        all([convergence.final_accurate]), ...
    "case_count", numel(results), ...
    "passed_count", nnz(passed), ...
    "failed_count", nnz(~passed), ...
    "mesh_family_count", numel(families), ...
    "maximum_free_residual_relative", ...
        max([results.free_residual_relative]), ...
    "maximum_final_sample_error", max([convergence.final_error]), ...
    "duration_s", toc(started), ...
    "convergence", convergence, ...
    "case_results", results);
end


function result = runManufacturedCase(mesh, meshId, dimension, order, meshPath)
started = tic;
result = emptyManufacturedResult();
result.case_id = meshId + "_p" + order;
result.mesh_id = meshId;
result.dimension = dimension;
result.order = order;
specification = struct("space", "h1", "order", order, ...
    "dirichlet", ".*", "weight_imag", 0);
space = makeSpace(mesh, specification);
spaceCleanup = onCleanup(@() delete(space));
coordinates = radia.ngsolve.CoefficientFunction.coordinates(dimension);
coordinatesCleanup = onCleanup(@() delete(coordinates));
one = radia.ngsolve.CoefficientFunction.constant(1.0);
oneCleanup = onCleanup(@() delete(one));
x = coordinates.component(1);
xCleanup = onCleanup(@() delete(x));
y = coordinates.component(2);
yCleanup = onCleanup(@() delete(y));
xx = x * (one - x);
yy = y * (one - y);
if dimension == 2
    exact = xx * yy;
    source = 2.0 * (xx + yy);
    points = [0.15, 0.20; 0.30, 0.25; 0.45, 0.20];
    expected = points(:,1) .* (1.0 - points(:,1)) .* ...
        points(:,2) .* (1.0 - points(:,2));
else
    z = coordinates.component(3);
    zCleanup = onCleanup(@() delete(z));
    zz = z * (one - z);
    points = [0.12, 0.18, 0.20; 0.22, 0.16, 0.45; 0.35, 0.20, 0.30];
    if meshId == "wedge_single"
        triangleBubble = x * y * (one - x - y);
        exact = triangleBubble * zz;
        source = 2.0 * (x + y) * zz + 2.0 * triangleBubble;
        expected = points(:,1) .* points(:,2) .* ...
            (1.0 - points(:,1) - points(:,2)) .* ...
            points(:,3) .* (1.0 - points(:,3));
    else
        exact = xx * yy * zz;
        source = 2.0 * (yy * zz + xx * zz + xx * yy);
        expected = points(:,1) .* (1.0 - points(:,1)) .* ...
            points(:,2) .* (1.0 - points(:,2)) .* ...
            points(:,3) .* (1.0 - points(:,3));
    end
end
exactCleanup = onCleanup(@() delete(exact));
sourceCleanup = onCleanup(@() delete(source));
unit = radia.ngsolve.CoefficientFunction.constant(1.0);
unitCleanup = onCleanup(@() delete(unit));
form = radia.ngsolve.BilinearForm.createFromCoefficient( ...
    space, "stiffness", unit, Label=result.case_id);
formCleanup = onCleanup(@() delete(form));
linear = radia.ngsolve.LinearForm.createFromCoefficient( ...
    space, source, Label=result.case_id + "_rhs");
linearCleanup = onCleanup(@() delete(linear));
matrix = form.matrix();
matrixCleanup = onCleanup(@() delete(matrix));
rhs = linear.vector();
rhsCleanup = onCleanup(@() delete(rhs));
inverse = matrix.inverse();
inverseCleanup = onCleanup(@() delete(inverse));
solution = inverse.matvec(rhs);
solutionCleanup = onCleanup(@() delete(solution));
residual = matrix.matvec(solution);
residualCleanup = onCleanup(@() delete(residual));
residual.axpy(-1.0, rhs);
free = space.freeDofs();
residualValues = residual.values();
rhsValues = rhs.values();
result.free_residual_relative = norm(residualValues(free)) / ...
    max(norm(rhsValues(free)), 1.0);
field = radia.ngsolve.GridFunction.fromFESpace(space, Name="solution");
fieldCleanup = onCleanup(@() delete(field));
field.setVector(solution.values());
fieldCoefficient = field.asCoefficient();
fieldCoefficientCleanup = onCleanup(@() delete(fieldCoefficient));
actual = fieldCoefficient.evaluate(meshPath, points);
result.dofs = space.DofCount;
result.sample_max_absolute_error = max(abs(actual - expected));
result.duration_s = toc(started);
result.passed = isfinite(result.sample_max_absolute_error) && ...
    result.sample_max_absolute_error < 0.1 && ...
    result.free_residual_relative < 2e-8;
end


function space = makeSpace(mesh, specification)
space = radia.ngsolve.FESpace.create(mesh, string(specification.space), ...
    specification.order, Dirichlet=string(specification.dirichlet), ...
    Complex=specification.weight_imag ~= 0);
end


function coefficient = makeCoefficient(specification, dimension)
if specification.weight_imag == 0
    weight = specification.weight_real;
else
    weight = complex(specification.weight_real, specification.weight_imag);
end
if string(specification.coefficient_kind) == "constant"
    coefficient = radia.ngsolve.CoefficientFunction.constant(weight);
    return
end
coordinates = radia.ngsolve.CoefficientFunction.coordinates(dimension);
x = coordinates.component(1);
base = radia.ngsolve.CoefficientFunction.constant(weight);
coefficient = base + 0.2 * x;
delete(base);
delete(x);
delete(coordinates);
end


function form = makeForm(space, coefficient, specification)
if logical(specification.boundary)
    form = radia.ngsolve.BilinearForm.createBoundaryFromCoefficient( ...
        space, coefficient, Label=string(specification.case_id));
else
    form = radia.ngsolve.BilinearForm.createFromCoefficient( ...
        space, string(specification.form), coefficient, ...
        Label=string(specification.case_id));
end
end


function result = emptyBreadthResult()
result = struct("case_id", "", "mesh_id", "", "dimension", 0, ...
    "space", "", "order", 0, "form", "", "coefficient_kind", "", ...
    "boundary", false, "is_complex", false, "solve", false, ...
    "dofs", 0, "nonzeros", 0, "dof_match", false, ...
    "free_dof_match", false, "nonzero_match", false, ...
    "matrix_relative_error", Inf, "matvec_relative_error", Inf, ...
    "solution_relative_error", 0, "duration_s", 0, ...
    "passed", false, "failure", "");
end


function result = failedBreadthResult(specification, exception)
result = emptyBreadthResult();
result.case_id = string(specification.case_id);
result.mesh_id = string(specification.mesh_id);
result.dimension = specification.dimension;
result.space = string(specification.space);
result.order = specification.order;
result.form = string(specification.form);
result.coefficient_kind = string(specification.coefficient_kind);
result.boundary = logical(specification.boundary);
result.is_complex = specification.weight_imag ~= 0;
result.solve = logical(specification.solve);
result.failure = string(getReport(exception, "extended", "hyperlinks", "off"));
end


function result = emptyScaleResult()
result = struct("case_id", "", "mesh_id", "", "dimension", 0, ...
    "space", "", "order", 0, "form", "", "dofs", 0, ...
    "dof_match", false, "matvec_relative_error", Inf, ...
    "energy_relative_error", Inf, "input_norm_relative_error", Inf, ...
    "matvec_norm_relative_error", Inf, "duration_s", 0, ...
    "passed", false, "failure", "");
end


function result = failedScaleResult(specification, exception)
result = emptyScaleResult();
result.case_id = string(specification.case_id);
result.mesh_id = string(specification.mesh_id);
result.dimension = specification.dimension;
result.space = string(specification.space);
result.order = specification.order;
result.form = string(specification.form);
result.failure = string(getReport(exception, "extended", "hyperlinks", "off"));
end


function result = emptyManufacturedResult()
result = struct("case_id", "", "mesh_id", "", "dimension", 0, ...
    "order", 0, "dofs", 0, "sample_max_absolute_error", Inf, ...
    "free_residual_relative", Inf, "duration_s", 0, ...
    "passed", false, "failure", "");
end


function result = summaryWithoutCases(report)
result = rmfield(report, "case_results");
end


function value = relativeError(actual, expected)
value = norm(actual - expected, "fro") / max(norm(expected, "fro"), 1.0);
end


function value = utcNow()
value = string(datetime("now", "TimeZone", "UTC", ...
    "Format", "yyyy-MM-dd'T'HH:mm:ss'Z'"));
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


function value = readJson(path)
if ~isfile(path)
    error("radia:validation:MissingTierResult", ...
        "A skipped validation tier requires an existing result: %s", path);
end
value = jsondecode(fileread(path));
end


function value = hostName()
value = string(getenv("COMPUTERNAME"));
if strlength(value) == 0
    value = string(java.net.InetAddress.getLocalHost().getHostName());
end
end


function value = fileSha256(path)
file = fopen(path, "rb");
assert(file >= 0, "Unable to open MEX file for hashing: %s", path);
cleanup = onCleanup(@() fclose(file));
bytes = fread(file, Inf, "*uint8");
digest = java.security.MessageDigest.getInstance("SHA-256");
digest.update(typecast(bytes, "int8"));
hashBytes = typecast(digest.digest(), "uint8");
value = lower(join(string(dec2hex(hashBytes, 2)), ""));
clear cleanup
end


function cleanupWorkDirectory(path, keepArtifacts)
if ~keepArtifacts && isfolder(path)
    rmdir(path, "s");
end
end
