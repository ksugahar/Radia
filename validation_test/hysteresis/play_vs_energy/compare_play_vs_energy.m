%COMPARE_PLAY_VS_ENERGY  MATLAB companion to compare_play_vs_energy.py.
%
%   Loads the Radia C++ comparison results saved by the Python script and
%   overlays the MATLAB forensic reference implementations:
%     - bie.BInputPlayModel        (Type 6 mirror, scalar 1-D)
%     - bie.BInputEnergyEggerModel (Type 5 mirror, scalar 1-D)
%
%   Outputs an overlay figure that confirms (or exposes) any C++ vs MATLAB
%   discrepancy on the same H trajectories.

clear; close all;

this_dir = fileparts(mfilename('fullpath'));
matlab_pkg_dir = ['w:\02_学会資料\2026年度\2026_09_IGTE_Symposium\菅原\matlab'];
if ~exist(matlab_pkg_dir, 'dir')
    error('compare_play_vs_energy:matlabPkg', ...
          'MATLAB +bie/ package not found at: %s', matlab_pkg_dir);
end
addpath(matlab_pkg_dir);

%% Load Python comparison results
npz_path = fullfile(this_dir, 'comparison_results.npz');
if ~isfile(npz_path)
    error('compare_play_vs_energy:noNpz', ...
          'Run compare_play_vs_energy.py first to generate %s', npz_path);
end
% MATLAB doesn't read .npz natively; use Python via system call to convert,
% or load from an equivalent .mat. The Python script could be extended to
% also save .mat — for now, load from .npz via a tiny Python helper.
mat_path = fullfile(this_dir, 'comparison_results.mat');
if ~isfile(mat_path)
    fprintf('Converting .npz to .mat via Python helper...\n');
    cmd = sprintf(['python -c "import numpy as np, scipy.io as sio; ' ...
                   'd=dict(np.load(r''%s'')); sio.savemat(r''%s'', d)"'], ...
                  npz_path, mat_path);
    [status, out] = system(cmd);
    if status ~= 0
        error('compare_play_vs_energy:convertFail', 'npz->mat conversion failed: %s', out);
    end
end

S = load(mat_path);
K = double(S.K);
eta = double(S.eta(:));
% Reconstruct f_tables from flat arrays
f_table_sizes = double(S.f_table_sizes(:));
f_r_flat = double(S.f_r_flat(:));
f_f_flat = double(S.f_f_flat(:));
f_tables = cell(K, 1);
cursor = 0;
for k = 1:K
    n = f_table_sizes(k);
    r = f_r_flat(cursor+1:cursor+n);
    f = f_f_flat(cursor+1:cursor+n);
    f_tables{k} = [r, f];
    cursor = cursor + n;
end

%% Test 1: sinusoidal H trajectory
fprintf('\n=== Test 1: MATLAB Type 6 Play on sinusoidal H ===\n');
H_sin = double(S.H_sin(:));
n = numel(H_sin);

modelP = bie.BInputPlayModel(eta, f_tables);
B_play_ml = zeros(n, 1);
B_prev = 0;
tic;
for i = 1:n
    [B_play_ml(i), ~] = modelP.inverse(H_sin(i), 'BInit', B_prev);
    modelP.commit();
    B_prev = B_play_ml(i);
end
t_play_ml = toc;

err_play = abs(B_play_ml - double(S.B_play_sin(:)));
fprintf('  MATLAB Type 6 Play vs Radia C++ Type 6:\n');
fprintf('    max|B_ML - B_Radia| = %.3e T  (%.3e %% of |B|max)\n', ...
    max(err_play), 100*max(err_play)/max(abs(B_play_ml)));
fprintf('    MATLAB wallclock: %.2f ms/step  (Radia: %.3f ms/step)\n', ...
    1e3*t_play_ml/n, 1e3*double(S.t_play_sin)/n);

%% Test 1b: MATLAB Egger Type 5 on same H trajectory
fprintf('\n=== Test 1b: MATLAB Egger Type 5 on sinusoidal H ===\n');
modelE = bie.BInputEnergyEggerModel(eta, f_tables);
B_energy_ml = zeros(n, 1);
tic;
for i = 1:n
    try
        [B_energy_ml(i), ~] = modelE.inverse(H_sin(i));
        modelE.commit();
    catch
        B_energy_ml(i) = NaN;
    end
end
t_energy_ml = toc;
fprintf('  MATLAB Type 5 Egger range: B in [%.3e, %.3e]\n', ...
    min(B_energy_ml), max(B_energy_ml));

%% Plots — overlay C++ and MATLAB
fig = figure('Position', [60 60 1200 800], 'Color', 'w');

subplot(2,2,1); hold on; grid on;
plot(H_sin, double(S.B_play_sin(:)),  'b-',  'LineWidth', 1.4, 'DisplayName', 'Radia Type 6');
plot(H_sin, B_play_ml,                 'b:',  'LineWidth', 1.0, 'DisplayName', 'MATLAB Type 6');
plot(H_sin, double(S.B_energy_sin(:)), 'r-',  'LineWidth', 1.4, 'DisplayName', 'Radia Type 5');
xlabel('H [A/m]'); ylabel('B [T]');
title('BH loop: Play vs Energy (sin H)'); legend('Location','best');

subplot(2,2,2);
semilogy(max(err_play, 1e-16), 'b-'); grid on;
xlabel('step n'); ylabel('|B_{MATLAB} - B_{Radia}| [T]');
title('Type 6 cross-validation error');

subplot(2,2,3); hold on; grid on;
H_pwm = double(S.H_pwm(:));
plot(H_pwm, double(S.B_play_pwm(:)),  'b-', 'DisplayName', 'Radia Type 6');
plot(H_pwm, double(S.B_energy_pwm(:)), 'r-', 'DisplayName', 'Radia Type 5');
xlabel('H [A/m]'); ylabel('B [T]');
title('BH loop: PWM-like H'); legend('Location','best');

subplot(2,2,4);
bar(categorical({'C++ Play','MATLAB Play','C++ Energy'}), ...
    [1e3*double(S.t_play_sin)/n, 1e3*t_play_ml/n, 1e3*double(S.t_energy_sin)/n]);
ylabel('wallclock per step [ms]'); set(gca, 'YScale', 'log');
title('Per-step cost (sin)'); grid on;

fig_dir = fullfile(this_dir, 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
out_png = fullfile(fig_dir, 'compare_play_vs_energy_matlab.png');
exportgraphics(fig, out_png, 'Resolution', 130);
fprintf('\nSaved overlay figure: %s\n', out_png);
