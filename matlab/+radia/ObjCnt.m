function container = ObjCnt(objects)
%OBJCNT Create a Radia container from object handles.

container = radia.internal.callMex('radia.ObjCnt', double(objects));
end
