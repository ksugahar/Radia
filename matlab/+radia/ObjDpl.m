function object = ObjDpl(source, option)
%OBJDPL Duplicate a Radia object.

if nargin < 2
    object = radia.internal.callMex('radia.ObjDpl', double(source));
else
    object = radia.internal.callMex( ...
        'radia.ObjDpl', double(source), char(string(option)));
end
end
