function object = ObjTetrahedronCurrent(vertices, currentDensity)
%OBJTETRAHEDRONCURRENT Create a tetrahedral constant-current-density source.

object = radia.internal.callMex( ...
    'radia.ObjTetrahedronCurrent', double(vertices), double(currentDensity));
end
