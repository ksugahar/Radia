function modelPath = buildIHSampleModel(options)
%BUILDIHSAMPLEMODEL Build the source-driven Radia IH sample object.
%   The sample contains separate masked IH Parameters, Eddy Current, and
%   Thermal subsystems. A Step block supplies the RMS current envelope, a
%   Ramp block supplies mechanical rotation angle, and a Constant supplies
%   ambient temperature. The object uses MATLAB and standard Simulink blocks
%   only; it does not launch Python or require a MEX backend.

arguments
    options.ModelName (1,1) string = "radia_ih_sample"
    options.OutputDirectory (1,1) string = ""
    options.Open (1,1) logical = false
end

plant = radia.simulink.makeIHPlant( ...
    HeatCapacity_J_per_K=150, ...
    ThermalConductance_W_per_K=3, ...
    SampleTime_s=0.05, ...
    InitialTemperature_K=293.15);

rotationAngle_rad = [0; pi / 2; pi; 3 * pi / 2];
coilCurrentRms_A = [0; 50; 100];
angleScale = [1.0; 1.2; 0.8; 1.1];
currentScale = (coilCurrentRms_A / 100).^2;
heatDensity_W_per_m3 = 4.0e6 * angleScale * currentScale.';
eddyLut = radia.simulink.makeIHEddyHeatDensityLUT( ...
    rotationAngle_rad, coilCurrentRms_A, heatDensity_W_per_m3, ...
    RegionVolumes_m3=2.5e-5, RegionNames="workpiece", ...
    CarrierFrequency_Hz=50e3, Source="Radia IH sample LUT");

radia.simulink.buildIHControlModel( ...
    options.ModelName, plant, eddyLut, StopTime_s=2.0, ...
    PlantBlock="standard", Save=false, Open=false);
root = options.ModelName;

delete_line(root, "coil_current_rms_A/1", "Eddy Current/1");
delete_block(root + "/coil_current_rms_A");
add_block("simulink/Sources/Step", root + "/coil_current_rms_A", ...
    "Time", "0.5", "Before", "40", "After", "80", ...
    "SampleTime", "0.05", "Position", [340 75 445 115]);
add_line(root, "coil_current_rms_A/1", "Eddy Current/1", ...
    "autorouting", "smart");

delete_line(root, "rotation_angle_rad/1", "Eddy Current/2");
delete_block(root + "/rotation_angle_rad");
add_block("simulink/Sources/Ramp", root + "/rotation_angle_rad", ...
    "slope", "2*pi", "start", "0", "InitialOutput", "0", ...
    "Position", [340 135 445 175]);
add_line(root, "rotation_angle_rad/1", "Eddy Current/2", ...
    "autorouting", "smart");

delete_line(root, "ambient_temperature_K/1", "Thermal/2");
delete_block(root + "/ambient_temperature_K");
add_block("simulink/Sources/Constant", root + "/ambient_temperature_K", ...
    "Value", "293.15", "Position", [555 205 640 245]);
add_line(root, "ambient_temperature_K/1", "Thermal/2", ...
    "autorouting", "smart");

set_param(root, "SaveOutput", "on", "OutputSaveName", "yout", ...
    "SaveTime", "on", "TimeSaveName", "tout");

outputDirectory = options.OutputDirectory;
if strlength(outputDirectory) == 0
    if ispc
        outputDirectory = "C:\temp\radia_ih_samples";
    else
        outputDirectory = fullfile(tempdir, "radia_ih_samples");
    end
end
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
modelPath = string(fullfile(outputDirectory, options.ModelName + ".slx"));
set_param(root, "SimulationCommand", "update");
save_system(root, modelPath);
if options.Open
    open_system(root);
end
end
