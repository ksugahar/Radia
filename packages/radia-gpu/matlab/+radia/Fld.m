function field = Fld(object, fieldType, points)
%FLD Evaluate a Radia field at one or more points.

field = radia.internal.callMex( ...
    'radia.Fld', double(object), char(string(fieldType)), double(points));
end
