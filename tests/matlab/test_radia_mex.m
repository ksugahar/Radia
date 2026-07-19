function tests = test_radia_mex
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

function testApiAndTaskManager(testCase)
api = radia.apiInfo();
probe = radia.taskmanagerProbe(200000);
verifyEqual(testCase, api.api_version, 1);
verifyGreaterThanOrEqual(testCase, api.taskmanager_max_threads, 1);
verifyGreaterThanOrEqual(testCase, probe.used_threads, 1);
verifyTrue(testCase, isfinite(probe.checksum));
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
verifySize(testCase, gram, [1, 1]);
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
delete(manager);
after = radia.apiInfo();
verifyEqual(testCase, after.handle_count, before.handle_count);
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
