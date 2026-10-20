function material = MatSatLamTab(mhData, packingFactor, normal)
%MATSATLAMTAB Create a laminated nonlinear material from an M-H table.

material = radia.internal.callMex('radia.MatSatLamTab', double(mhData), ...
    double(packingFactor), double(normal));
end
