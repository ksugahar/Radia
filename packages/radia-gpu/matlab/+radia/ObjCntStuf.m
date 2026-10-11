function objects = ObjCntStuf(container)
%OBJCNTSTUF Return the object handles in a Radia container.

objects = radia.internal.callMex('radia.ObjCntStuf', double(container));
end
