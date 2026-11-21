function report = rotating_magnet_eddy_gmsh_replay(sourceH5, resampledH5, steps, outJson)
%ROTATING_MAGNET_EDDY_GMSH_REPLAY MATLAB replay of the Gmsh-resampled record.
%   report = rotating_magnet_eddy_gmsh_replay(sourceH5, resampledH5, steps, outJson)
%   reads both HDF5 files with h5read only and checks, for the 0-based step
%   indices in STEPS (default [5 18 36], nonadjacent):
%     - time_s, solver_step and magnet pose equal the source entries;
%     - sample ids and cell-centre coordinates of the query points;
%     - signed J (each component) against an independent MATLAB evaluation of
%       the Quad4 interpolant at the cell centres (mean of the four corners of
%       the source grid), not the norm alone.
%   h5read returns arrays in reversed dimension order: J is 3 x m x n.
%   Writes a JSON report and errors on any mismatch.
if nargin < 3 || isempty(steps), steps = [5 18 36]; end
if nargin < 4, outJson = 'matlab_replay.json'; end
tol = 1e-12;

tS = h5read(sourceH5, '/time_s');            % n x 1
JS = h5read(sourceH5, '/field/J_A_per_m2');  % 3 x m x n
idS = double(h5read(sourceH5, '/samples/id'));
pS = h5read(sourceH5, '/samples/point_m');   % 3 x m
stepS = double(h5read(sourceH5, '/solver_step'));
cS = h5read(sourceH5, '/magnet/center_m');   % 3 x n
rS = h5read(sourceH5, '/magnet/rotation');   % 3 x 3 x n (transposed per step)

tR = h5read(resampledH5, '/time_s');
JR = h5read(resampledH5, '/field/J_A_per_m2');
idR = double(h5read(resampledH5, '/samples/id'));
pR = h5read(resampledH5, '/samples/point_m');
srcIdx = double(h5read(resampledH5, '/source_step_index'));
stepR = double(h5read(resampledH5, '/solver_step'));
cR = h5read(resampledH5, '/magnet/center_m');
rR = h5read(resampledH5, '/magnet/rotation');

% Source grid: id = i*n2 + k on the plate mid-plane (x = axis 1, z = axis 3).
[~, order] = sort(idS);
P = pS(:, order);
xs = unique(P(1, :)); zs = unique(P(3, :));
n1 = numel(xs); n2 = numel(zs);
assert(n1 * n2 == numel(idS), 'source samples are not a %d x %d grid', n1, n2);
J = JS(:, order, :);
G = reshape(J, 3, n2, n1, []);               % 3 x k x i x step
xc = (xs(1:end-1) + xs(2:end)) / 2; zc = (zs(1:end-1) + zs(2:end)) / 2;
[Zc, Xc] = ndgrid(zc, xc);                   % k fastest, as id = i*(n2-1) + k
expectId = (0:(n1-1)*(n2-1)-1)';
assert(isequal(idR(:), expectId), 'query ids differ from i*(n2-1)+k');
dp = max(abs([pR(1, :) - Xc(:)'; pR(3, :) - Zc(:)'; pR(2, :) - P(2, 1)]), [], 'all');
assert(dp < 1e-15, 'query coordinates differ from cell centres by %g m', dp);

rows = struct([]);
for s = steps
    j = s + 1; src = srcIdx(j) + 1;          % 1-based
    o = 0.25 * (G(:, 1:end-1, 1:end-1, src) + G(:, 2:end, 1:end-1, src) ...
              + G(:, 2:end, 2:end, src) + G(:, 1:end-1, 2:end, src));
    o = reshape(o, 3, []);
    v = JR(:, :, j);
    err = max(abs(v - o), [], 'all') / max(abs(o), [], 'all');
    perComp = max(abs(v - o), [], 2)' ./ max(abs(o), [], 2)';
    dt = abs(tR(j) - tS(src));
    dpose = max([abs(cR(:, j) - cS(:, src)); reshape(abs(rR(:, :, j) - rS(:, :, src)), [], 1)]);
    rows(end+1).step = s; %#ok<AGROW>
    rows(end).time_s = tR(j);
    rows(end).solver_step = stepR(j);
    rows(end).time_err_s = dt;
    rows(end).pose_err = dpose;
    rows(end).J_max_rel_err = err;
    rows(end).J_component_rel_err = perComp;
    rows(end).J_sample = v(:, 1)';
    fprintf('step %d t=%.6f s solver_step=%d: |dt|=%.1e pose=%.1e J rel=%.2e [%s]\n', ...
        s, tR(j), stepR(j), dt, dpose, err, num2str(perComp, '%.1e '));
    assert(dt == 0 && stepR(j) == stepS(src) && dpose == 0, 'time/step/pose mismatch at step %d', s);
    assert(err < tol, 'signed J mismatch at step %d: %g', s, err);
end
report = struct('shape_J_matlab', size(JR), 'unit', 'A/m^2', 'grid', [n1 n2], ...
                'steps', rows, 'tolerance', tol, 'matlab', version);
fid = fopen(outJson, 'w'); fwrite(fid, jsonencode(report, 'PrettyPrint', true)); fclose(fid);
fprintf('MATLAB replay passed: J %s (3 x m x n), report %s\n', mat2str(size(JR)), outJson);
end
