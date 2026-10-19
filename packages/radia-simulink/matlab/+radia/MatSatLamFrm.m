function material = MatSatLamFrm(ksiMs1, ksiMs2, ksiMs3, packing, normal)
%MATSATLAMFRM Create a laminated nonlinear formula material.

material = radia.internal.callMex('radia.MatSatLamFrm', double(ksiMs1), ...
    double(ksiMs2), double(ksiMs3), double(packing), double(normal));
end
