function object = ObjWedge(vertices, magnetization)
%OBJWEDGE Create a magnetized wedge Radia object.

object = radia.internal.callMex( ...
    'radia.ObjWedge', double(vertices), double(magnetization));
end
