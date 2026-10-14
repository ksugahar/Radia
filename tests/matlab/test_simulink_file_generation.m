function tests = test_simulink_file_generation
%TEST_SIMULINK_FILE_GENERATION Verify repository-external generated folders.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
addpath(fullfile(repoRoot, "matlab"), "-begin");
testCase.TestData.OriginalConfig = Simulink.fileGenControl("getConfig");
testCase.TestData.RootA = string(tempname("C:\temp"));
testCase.TestData.RootB = string(tempname("C:\temp"));
testCase.TestData.RepoRoot = repoRoot;
end

function teardownOnce(testCase)
Simulink.fileGenControl("setConfig", ...
    config=testCase.TestData.OriginalConfig);
roots = [testCase.TestData.RootA, testCase.TestData.RootB];
for root = roots
    if isfolder(root)
        rmdir(root, "s");
    end
end
end

function testConfigurePreserveAndForce(testCase)
Simulink.fileGenControl("set", ...
    CacheFolder=testCase.TestData.RepoRoot, ...
    CodeGenFolder=testCase.TestData.RepoRoot);
first = radia.simulink.configureFileGeneration( ...
    RootDirectory=testCase.TestData.RootA);
verifyTrue(testCase, first.available);
verifyTrue(testCase, first.changed);
verifyEqual(testCase, first.cache_folder, ...
    fullfile(testCase.TestData.RootA, "cache"));
verifyEqual(testCase, first.codegen_folder, ...
    fullfile(testCase.TestData.RootA, "codegen"));
verifyTrue(testCase, isfolder(first.cache_folder));
verifyTrue(testCase, isfolder(first.codegen_folder));

preserved = radia.simulink.configureFileGeneration( ...
    RootDirectory=testCase.TestData.RootB);
verifyFalse(testCase, preserved.changed);
verifyEqual(testCase, preserved.cache_folder, first.cache_folder);
verifyEqual(testCase, preserved.codegen_folder, first.codegen_folder);

forced = radia.simulink.configureFileGeneration( ...
    RootDirectory=testCase.TestData.RootB, Force=true);
verifyTrue(testCase, forced.changed);
verifyEqual(testCase, forced.cache_folder, ...
    fullfile(testCase.TestData.RootB, "cache"));
verifyEqual(testCase, forced.codegen_folder, ...
    fullfile(testCase.TestData.RootB, "codegen"));
end
