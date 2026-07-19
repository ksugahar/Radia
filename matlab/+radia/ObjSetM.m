function ObjSetM(object, magnetization)
%OBJSETM Set the magnetization of a Radia object.

radia.internal.callMex( ...
    'radia.ObjSetM', double(object), double(magnetization));
end
