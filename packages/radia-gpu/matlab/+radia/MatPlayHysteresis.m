function material = MatPlayHysteresis(K, eta, tables)
%MATPLAYHYSTERESIS Create a direct B-input play hysteresis material.

[rFlat, fFlat, tableSizes] = radia.internal.flattenHysteresisTables(K, tables);
material = radia.internal.callMex('radia.MatPlayHysteresis', double(K), ...
    double(eta), rFlat, fFlat, tableSizes);
end
