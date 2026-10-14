function object = ObjRecCur(center, dimensions, currentDensity)
%OBJRECCUR Create a rectangular current-carrying block.

object = radia.internal.callMex('radia.ObjRecCur', double(center), ...
    double(dimensions), double(currentDensity));
end
