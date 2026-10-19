function radia_bh_sfun(block)
% Temperature-dependent BH block. Supports a formula or a bilinear LUT.
setup(block);
end

function setup(block)
block.NumDialogPrms = 1;
block.NumInputPorts = 2;
block.NumOutputPorts = 2;
for k = 1:2
    block.InputPort(k).Dimensions = 1;
    block.InputPort(k).DatatypeID = 0;
    block.InputPort(k).DirectFeedthrough = true;
end
block.OutputPort(1).Dimensions = 1;
block.OutputPort(2).Dimensions = 1;
block.SampleTimes = [0 0];
block.SimStateCompliance = 'DefaultSimState';
block.RegBlockMethod('Outputs', @outputs);
end

function outputs(block)
cfg = block.DialogPrm(1).Data;
T = block.InputPort(1).Data;
H = block.InputPort(2).Data;
if ~isstruct(cfg) || ~isfield(cfg, 'mode')
    error('radia:bh:Config', 'BH configuration must contain mode.');
end
mode = lower(string(cfg.mode));
if mode == "formula"
    mu0 = 4*pi*1e-7;
    mu_r = cfg.mu_r_ref .* (1 + cfg.mu_r_temperature_slope .* (T - cfg.reference_temperature_K));
    if ~(isfinite(mu_r) && mu_r > 0)
        error('radia:bh:Material', 'Formula produced an invalid permeability.');
    end
    B = mu0 .* mu_r .* H;
    dBdH = mu0 .* mu_r;
elseif mode == "lut"
    [B, dBdH] = lut_eval(cfg, T, H);
else
    error('radia:bh:Mode', 'BH mode must be formula or lut.');
end
block.OutputPort(1).Data = B;
block.OutputPort(2).Data = dBdH;
end

function [B, dBdH] = lut_eval(cfg, T, H)
temps = cfg.temperature_K(:); fields = cfg.H_A_per_m(:);
table = cfg.B_T;
if ~isequal(size(table), [numel(temps), numel(fields)])
    error('radia:bh:LUT', 'B_T must be temperature-by-field.');
end
T = min(max(T, temps(1)), temps(end));
H = min(max(H, fields(1)), fields(end));
rows = zeros(numel(fields), 1);
for j = 1:numel(fields)
    rows(j) = interp1(temps, table(:, j), T, 'linear');
end
B = interp1(fields, rows, H, 'linear');
dBdH = (interp1(fields, rows, H + eps(max(1, abs(H))), 'linear', 'extrap') - ...
    interp1(fields, rows, H - eps(max(1, abs(H))), 'linear', 'extrap')) / (2*eps(max(1, abs(H))));
end
