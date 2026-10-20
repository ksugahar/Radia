function ObjAddToCnt(container, objects)
%OBJADDTOCNT Add object handles to an existing Radia container.

radia.internal.callMex( ...
    'radia.ObjAddToCnt', double(container), double(objects));
end
