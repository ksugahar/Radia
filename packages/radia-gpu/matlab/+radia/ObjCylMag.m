function object = ObjCylMag(center, radius, height, nseg, axis, magnetization)
%OBJCYLMAG Create a uniformly magnetized cylinder.

object = radia.internal.callMex('radia.ObjCylMag', double(center), double(radius), ...
    double(height), double(nseg), char(string(axis)), double(magnetization));
end
