function values = FldLst(object, fieldId, p1, p2, np, argOpt, start)
%FLDLST Evaluate a field along an equidistant line.

arguments
    object (1,1) double
    fieldId (1,1) string
    p1 (1,3) double
    p2 (1,3) double
    np (1,1) double {mustBeInteger, mustBePositive}
    argOpt (1,1) string = "noarg"
    start (1,1) double = 0
end

values = radia.internal.callMex('radia.FldLst', double(object), char(fieldId), ...
    double(p1), double(p2), double(np), char(argOpt), double(start));
end
