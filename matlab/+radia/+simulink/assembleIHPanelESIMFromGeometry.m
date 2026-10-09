function config = assembleIHPanelESIMFromGeometry(wpVol, coilFile, configFile, options)
%ASSEMBLEIHPANELESIMFROMGEOMETRY Freeze certified panel ESIM at assembly current.
% Python outer Picard runs once at geometry update/initialization. Native time
% stepping uses this fixed material state; changing reference current rebuilds.
% options.python_executable selects the same explicit process-Python boundary
% as assembleIHOperatorsFromGeometry; the ESIM material solver is not duplicated.
arguments
    wpVol (1,1) string
    coilFile (1,1) string
    configFile (1,1) string
    options (1,1) struct = struct()
end
options.zs_mode = "per-panel-esim";
config = radia.simulink.assembleIHOperatorsFromGeometry(wpVol, coilFile, configFile, options);
end
