function object = ObjPyramid(vertices, magnetization)
%OBJPYRAMID Create a magnetized pyramid Radia object.

object = radia.internal.callMex( ...
    'radia.ObjPyramid', double(vertices), double(magnetization));
end
