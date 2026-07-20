function material = MatEnergyHysteresis(K, chi, tables, eps)
%MATENERGYHYSTERESIS Create an energy-based vector hysteresis material.

if nargin < 4
    eps = 1e-8;
end
[rFlat, fFlat, tableSizes] = radia.internal.flattenHysteresisTables(K, tables);
material = radia.internal.callMex('radia.MatEnergyHysteresis', double(K), ...
    double(chi), rFlat, fFlat, tableSizes, double(eps));
end
