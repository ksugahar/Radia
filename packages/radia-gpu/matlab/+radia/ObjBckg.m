function object = ObjBckg(fieldVector)
%OBJBCKG Create a uniform background magnetic-field source.

object = radia.internal.callMex('radia.ObjBckg', double(fieldVector));
end
