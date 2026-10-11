function material = MatSatIsoTab(bhData)
%MATSATISOTAB Create a nonlinear isotropic material from [H,B] rows.

material = radia.internal.callMex('radia.MatSatIsoTab', double(bhData));
end
