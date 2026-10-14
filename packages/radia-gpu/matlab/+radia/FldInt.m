function integral = FldInt(object, infFin, fieldId, p1, p2)
%FLDINT Integrate a field component along a straight line.

integral = radia.internal.callMex('radia.FldInt', double(object), ...
    char(string(infFin)), char(string(fieldId)), double(p1), double(p2));
end
