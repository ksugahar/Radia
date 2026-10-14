function object = ObjRaceTrk(center, radii, lengths, height, nseg, manAuto, axis, currentDensity)
%OBJRACETRK Create a racetrack current coil.

object = radia.internal.callMex('radia.ObjRaceTrk', double(center), double(radii), ...
    double(lengths), double(height), double(nseg), char(string(manAuto)), ...
    char(string(axis)), double(currentDensity));
end
