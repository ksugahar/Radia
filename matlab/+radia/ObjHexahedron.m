function object = ObjHexahedron(vertices, magnetization)
%OBJHEXAHEDRON Create a magnetized hexahedral Radia object.

object = radia.internal.callMex( ...
    'radia.ObjHexahedron', double(vertices), double(magnetization));
end
