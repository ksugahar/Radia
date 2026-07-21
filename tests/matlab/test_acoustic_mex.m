function tests = test_acoustic_mex
% Cross-language parity and error-contract tests for native acoustic MEX APIs.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
addpath(fullfile(repoRoot, "matlab"));
radia.setup(Force=true);
[path, reference] = loadPythonReference(repoRoot);
testCase.TestData.referencePath = path;
testCase.TestData.reference = reference;
end

function teardownOnce(testCase)
if isfile(testCase.TestData.referencePath)
    delete(testCase.TestData.referencePath);
end
end

function testRealSphereWrappersMatchPython(testCase)
reference = testCase.TestData.reference;
soft = radia.acoustic.softSphereScattering(2.3, 0.9, reference.points, Terms=18);
rigid = radia.acoustic.rigidSphereScattering(2.3, 0.9, reference.points, Terms=18);
fluid = radia.acoustic.fluidSphereScattering(2.3, 0.9, reference.mixed_points, ...
    InteriorWavenumber=1.4, DensityRatio=1.75, Terms=18);
elastic = radia.acoustic.elasticSphereScattering(2.3, 0.9, reference.points, ...
    LongitudinalSpeed=2.4, ShearSpeed=1.1, DensityRatio=1.8, Terms=15);

verifyEqual(testCase, soft.scattered, reference.soft_scattered, ...
    RelTol=2e-13, AbsTol=2e-14);
verifyEqual(testCase, soft.total, reference.soft_total, ...
    RelTol=2e-13, AbsTol=2e-14);
verifyEqual(testCase, rigid.scattered, reference.rigid_scattered, ...
    RelTol=2e-13, AbsTol=2e-14);
verifyEqual(testCase, fluid.total, reference.fluid_total, ...
    RelTol=2e-13, AbsTol=2e-14);
verifyEqual(testCase, fluid.inside_mask, logical(reference.fluid_inside_mask));
verifyEqual(testCase, elastic.scattered, reference.elastic_scattered, ...
    RelTol=3e-13, AbsTol=3e-14);
verifyEqual(testCase, string(soft.backend), "native-mex");
verifyEqual(testCase, string(elastic.backend), "native-mex");
end

function testComplexSphereBdfAndCqGridMatchPython(testCase)
reference = testCase.TestData.reference;
actual = radia.acoustic.softSphereScatteringComplexK( ...
    reference.complex_wavenumber, 0.9, reference.points, Terms=18);
verifyEqual(testCase, actual, reference.complex_scattered, ...
    RelTol=3e-13, AbsTol=3e-14);
verifyEqual(testCase, radia.acoustic.bdfDelta(reference.zeta, "BDF1"), ...
    reference.bdf1, RelTol=0, AbsTol=2e-15);
verifyEqual(testCase, radia.acoustic.bdfDelta(reference.zeta, "BDF2"), ...
    reference.bdf2, RelTol=0, AbsTol=2e-15);

verifyGridMatchesPython(testCase, reference, "bdf1", "BDF1", 15);
verifyGridMatchesPython(testCase, reference, "bdf2", "BDF2", 16);
end

function verifyGridMatchesPython(testCase, reference, prefix, method, numTime)
grid = radia.acoustic.cqGrid(numTime, 0.08, SoundSpeed=1.2, Method=method);
verifyEqual(testCase, grid.cq_radius, reference.(prefix + "_cq_radius"), ...
    RelTol=0, AbsTol=2e-15);
verifyEqual(testCase, grid.zeta, reference.(prefix + "_cq_zeta"), ...
    RelTol=0, AbsTol=2e-15);
verifyEqual(testCase, grid.cq_nodes, reference.(prefix + "_cq_nodes"), ...
    RelTol=2e-15, AbsTol=2e-14);
verifyEqual(testCase, grid.cq_wavenumbers, ...
    reference.(prefix + "_cq_wavenumbers"), ...
    RelTol=2e-15, AbsTol=2e-14);
verifyEqual(testCase, string(grid.backend), "native-mex");
end

function testExteriorContractFailsLoudly(testCase)
call = @() radia.acoustic.softSphereScattering(2.0, 1.0, [0, 0, 0.5]);
verifyError(testCase, call, "radia:mex:Exception");
invalidTerms = @() radia_mex( ...
    'acoustic.soft_sphere', 2.0, 1.0, [0, 0, 1.5], -2);
verifyError(testCase, invalidTerms, "radia:mex:Exception");
end

function [path, reference] = loadPythonReference(repoRoot)
path = string(tempname("C:\temp")) + ".mat";
pythonExecutable = string(getenv("RADIA_PYTHON_EXECUTABLE"));
if strlength(pythonExecutable) == 0
    pythonExecutable = "python";
end
helper = fullfile(repoRoot, "tests", "matlab", "acoustic_python_reference.py");
command = quoteCommandArgument(pythonExecutable) + " " + ...
    quoteCommandArgument(helper) + " " + quoteCommandArgument(path);
[status, output] = radia.internal.runPythonProcess(command);
if status ~= 0
    error("radia:test:PythonReference", ...
        "Python acoustic reference failed (%d): %s", status, output);
end
reference = load(path);
end

function value = quoteCommandArgument(value)
value = '"' + replace(string(value), '"', '""') + '"';
end
