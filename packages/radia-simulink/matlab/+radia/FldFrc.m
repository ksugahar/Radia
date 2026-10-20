function force = FldFrc(object, shape)
%FLDFRC Compute force using a Maxwell-tensor surface or object shape.

force = radia.internal.callMex('radia.FldFrc', double(object), double(shape));
end
