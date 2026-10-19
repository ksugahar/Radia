function [rFlat, fFlat, tableSizes] = flattenHysteresisTables(K, tables)
%FLATTENHYSTERESISTABLES Convert MATLAB tables to the native flat contract.

if ~iscell(tables) || numel(tables) ~= K
    error("radia:hysteresis:Tables", "tables must be a cell array of length K.");
end
rFlat = zeros(0, 1);
fFlat = zeros(0, 1);
tableSizes = zeros(K, 1);
for k = 1:K
    table = tables{k};
    if iscell(table) && numel(table) == 2
        r = double(table{1}(:));
        f = double(table{2}(:));
    elseif isnumeric(table) && size(table, 2) == 2
        r = double(table(:, 1));
        f = double(table(:, 2));
    else
        error("radia:hysteresis:TableShape", ...
            "Each table must be an N-by-2 matrix or a {r,f} cell.");
    end
    if isempty(r) || numel(r) ~= numel(f)
        error("radia:hysteresis:TableLength", ...
            "Each r/f table must be non-empty and have matching lengths.");
    end
    tableSizes(k) = numel(r);
    rFlat = [rFlat; r]; %#ok<AGROW>
    fFlat = [fFlat; f]; %#ok<AGROW>
end
end
