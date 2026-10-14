function object = ObjTetrahedron(vertices, magnetization)
%OBJTETRAHEDRON Create a magnetized tetrahedral Radia object.

object = radia.internal.callMex( ...
    'radia.ObjTetrahedron', double(vertices), double(magnetization));
end
