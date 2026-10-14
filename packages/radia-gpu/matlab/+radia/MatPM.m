function material = MatPM(Br, Hc, easyAxis)
%MATPM Create a permanent-magnet material.

material = radia.internal.callMex( ...
    'radia.MatPM', double(Br), double(Hc), double(easyAxis));
end
