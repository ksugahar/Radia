function object = ObjArcCur(center, radii, phi, height, nseg, manAuto, axis, currentDensity)
%OBJARCCUR Create a finite-length arc current coil.

object = radia.internal.callMex('radia.ObjArcCur', double(center), double(radii), ...
    double(phi), double(height), double(nseg), char(string(manAuto)), ...
    char(string(axis)), double(currentDensity));
end
