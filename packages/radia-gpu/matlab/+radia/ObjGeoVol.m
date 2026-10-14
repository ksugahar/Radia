function volume = ObjGeoVol(object)
%OBJGEOVOL Return the geometric volume of a Radia object.

volume = radia.internal.callMex('radia.ObjGeoVol', double(object));
end
