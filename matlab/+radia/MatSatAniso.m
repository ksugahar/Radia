function material = MatSatAniso(dataPar, dataPer)
%MATSATANISO Create a nonlinear anisotropic material from two M-H curves.

material = radia.internal.callMex('radia.MatSatAniso', double(dataPar), double(dataPer));
end
