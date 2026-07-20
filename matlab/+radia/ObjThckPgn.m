function object = ObjThckPgn(xc, lx, polygon, axis, magnetization)
%OBJTHCKPGN Create a uniformly magnetized extruded polygon.

object = radia.internal.callMex('radia.ObjThckPgn', double(xc), double(lx), ...
    double(polygon), char(string(axis)), double(magnetization));
end
