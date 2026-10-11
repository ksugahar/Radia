function object = ObjFlmCur(points, current)
%OBJFLMCUR Create a polygonal filament current.

object = radia.internal.callMex('radia.ObjFlmCur', double(points), double(current));
end
