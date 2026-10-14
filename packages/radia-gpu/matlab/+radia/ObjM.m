function state = ObjM(object)
%OBJM Return center and magnetization arrays for a Radia object.

state = radia.internal.callMex('radia.ObjM', double(object));
end
