function count = ObjCntSize(container)
%OBJCNTSize Return the number of objects in a Radia container.

count = radia.internal.callMex('radia.ObjCntSize', double(container));
end
