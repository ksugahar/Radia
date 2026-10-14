function ObjScaleCur(object, scale)
%OBJSCALECUR Scale the current of a Radia current object.

radia.internal.callMex( ...
    'radia.ObjScaleCur', double(object), double(scale));
end
