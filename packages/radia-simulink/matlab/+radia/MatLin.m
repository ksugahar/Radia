function material = MatLin(muR, easyAxis)
%MATLIN Create an isotropic or uniaxial linear magnetic material.

if nargin < 2
    material = radia.internal.callMex('radia.MatLin', double(muR));
else
    material = radia.internal.callMex( ...
        'radia.MatLin', double(muR), double(easyAxis));
end
end
