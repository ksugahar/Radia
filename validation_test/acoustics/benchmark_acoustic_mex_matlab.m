function result = benchmark_acoustic_mex_matlab(outputPath, matlabRoot, options)
%BENCHMARK_ACOUSTIC_MEX_MATLAB Benchmark the shared acoustic kernel via MEX.
arguments
    outputPath (1,1) string
    matlabRoot (1,1) string
    options.PointCount (1,1) double {mustBeInteger,mustBePositive} = 20000
    options.Terms (1,1) double {mustBeInteger,mustBeNonnegative} = 28
    options.Warmup (1,1) double {mustBeInteger,mustBeNonnegative} = 5
    options.Repeats (1,1) double {mustBeInteger,mustBePositive} = 31
end

setupTimer = tic;
addpath(matlabRoot);
radia.setup(Force=true);
setupSeconds = toc(setupTimer);

index = (0:options.PointCount-1).';
points = [ ...
    1.2 + mod(index,997)/997, ...
    (mod(index*17,991)-495)/4000, ...
    (mod(index*31,983)-491)/3500];
zeta = linspace(-0.8,0.8,options.PointCount) + ...
    1i*linspace(0.4,-0.4,options.PointCount);

scattering = @() radia_mex('acoustic.soft_sphere', ...
    3.1, 1.0, points, options.Terms);
transfer = @() radia_mex('acoustic.bdf_delta', zeta, 'BDF2');
firstTimer = tic;
first = scattering();
firstSeconds = toc(firstTimer);
for index = 1:options.Warmup
    warmScattering = scattering(); %#ok<NASGU>
    warmTransfer = transfer(); %#ok<NASGU>
end
scatteringSamples = zeros(options.Repeats,1);
transferSamples = zeros(options.Repeats,1);
for index = 1:options.Repeats
    timer = tic;
    sample = scattering(); %#ok<NASGU>
    scatteringSamples(index) = toc(timer);
end
for index = 1:options.Repeats
    timer = tic;
    sample = transfer(); %#ok<NASGU>
    transferSamples(index) = toc(timer);
end

result = struct( ...
    "schema", "radia.acoustic-backend-benchmark/v1", ...
    "backend", "matlab-mex", ...
    "host", hostName(), ...
    "matlab_version", string(version), ...
    "platform", string(computer), ...
    "point_count", options.PointCount, ...
    "terms", options.Terms, ...
    "warmup", options.Warmup, ...
    "repeats", options.Repeats, ...
    "setup_s", setupSeconds, ...
    "first_scattering_s", firstSeconds, ...
    "scattering_median_s", median(scatteringSamples), ...
    "scattering_min_s", min(scatteringSamples), ...
    "bdf_transfer_median_s", median(transferSamples), ...
    "bdf_transfer_min_s", min(transferSamples), ...
    "checksum_real", sum(real(first.scattered)), ...
    "checksum_imag", sum(imag(first.scattered)));

file = fopen(outputPath, "w");
if file < 0
    error("radia:benchmark:Output", "Could not open %s", outputPath);
end
cleanup = onCleanup(@() fclose(file));
fwrite(file, jsonencode(result, PrettyPrint=true), "char");
fprintf("%s\n", jsonencode(result, PrettyPrint=true));
end

function value = hostName()
value = string(getenv("COMPUTERNAME"));
if strlength(value) == 0
    [status, output] = system("hostname");
    if status == 0
        value = strtrim(string(output));
    end
end
end
