function checkCoreVolume(coreVolume_m3,coreArea_m2,ironPath_m)
%CHECKCOREVOLUME Reject an iron volume that contradicts the lumped circuit.
%   The hysteretic magnetic circuit uses one cross-section A and one iron
%   path l, so its iron volume is A*l by construction. A separately supplied
%   volume is redundant; when it disagrees, the reported iron work would be
%   rescaled by V/(A*l) and, in a gapped core where the air gap stores almost
%   all of the energy, nothing downstream would notice. NaN means "not
%   supplied" and derives the volume.
if isnan(coreVolume_m3),return,end
ironVolume=coreArea_m2*ironPath_m;
if ~(coreVolume_m3>0)||abs(coreVolume_m3-ironVolume)>1e-9*ironVolume
    error("radia:simulink:HystereticCoreVolumeMismatch", ...
        "CoreVolume_m3=%.17g contradicts CoreArea_m2*MagneticPath_m=%.17g. The lumped circuit fixes the iron volume to A*l; omit CoreVolume_m3 or make it consistent.", ...
        coreVolume_m3,ironVolume);
end
end
