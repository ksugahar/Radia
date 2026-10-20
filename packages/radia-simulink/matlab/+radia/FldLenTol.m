function FldLenTol(absValue, relValue, zeroValue)
%FLDLENTOL Set length randomization tolerances.

if nargin < 3
    zeroValue = 0;
end
radia.internal.callMex('radia.FldLenTol', double(absValue), double(relValue), ...
    double(zeroValue));
end
