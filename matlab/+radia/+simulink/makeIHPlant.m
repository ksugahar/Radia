function plant = makeIHPlant(options)
%MAKEIHPLANT Create a discrete lumped thermal plant for Simulink.
%   PLANT = radia.simulink.makeIHPlant(...) returns the matrices and port
%   metadata used by a Simulink Discrete State-Space block. The input is
%   [power_W; ambient_temperature_K]. The outputs are
%   [temperature_K; heat_loss_W; energy_input_J; temperature_rate_K_per_s].
%
%   The electromagnetic calculation is intentionally outside this plant.
%   Radia VIM/FEM/SIBC/ESIM supplies power_W (or a position/temperature LUT
%   supplies it) and this model handles the slower thermal/control dynamics.

arguments
    options.HeatCapacity_J_per_K (1,1) double {mustBePositive} = 1.0
    options.ThermalConductance_W_per_K (1,1) double {mustBeNonnegative} = 0.0
    options.SampleTime_s (1,1) double {mustBePositive} = 1.0e-3
    options.InitialTemperature_K (1,1) double {mustBeFinite} = 293.15
    options.InitialEnergy_J (1,1) double {mustBeFinite} = 0.0
    options.AmbientTemperature_K (1,1) double {mustBeFinite} = 293.15
    options.Emissivity (1,1) double {mustBeNonnegative, mustBeLessThanOrEqual(options.Emissivity, 1)} = 0.0
    options.RadiatingArea_m2 (1,1) double {mustBeNonnegative} = 0.0
    options.RadiationLinearizationTemperature_K (1,1) double {mustBePositive} = 293.15
    options.UseRadiationLinearization (1,1) logical = false
end

C = options.HeatCapacity_J_per_K;
G_conv = options.ThermalConductance_W_per_K;
sigmaSB = 5.670374419e-8;
if options.UseRadiationLinearization && options.Emissivity > 0 && options.RadiatingArea_m2 > 0
    Tref = options.RadiationLinearizationTemperature_K;
    G_rad = 4 * options.Emissivity * sigmaSB * options.RadiatingArea_m2 * Tref^3;
else
    G_rad = 0.0;
end
G = G_conv + G_rad;
Ts = options.SampleTime_s;

% Exact zero-order-hold discretization of
% C*dT/dt = P - G*(T - T_ambient).
if G > 0
    a = exp(-G * Ts / C);
    bPower = (1 - a) / G;
    bAmbient = 1 - a;
else
    a = 1.0;
    bPower = Ts / C;
    bAmbient = 0.0;
end

% The second state accumulates the supplied heat input. This is useful for
% heating-time objectives and gives the controller an auditable energy port.
A = [a, 0; 0, 1];
B = [bPower, bAmbient; Ts, 0];
Cout = [1, 0; G, 0; 0, 1; -G / C, 0];
D = [0, 0; 0, -G; 0, 0; 1 / C, G / C];
x0 = [options.InitialTemperature_K; options.InitialEnergy_J];

plant = struct( ...
    "schema", "radia.ih.simulink.plant.v1", ...
    "A", A, "B", B, "C", Cout, "D", D, "x0", x0, ...
    "sample_time_s", Ts, ...
    "heat_capacity_J_per_K", C, ...
    "thermal_conductance_W_per_K", G, ...
    "convective_conductance_W_per_K", G_conv, ...
    "radiative_conductance_W_per_K", G_rad, ...
    "ambient_temperature_K", options.AmbientTemperature_K, ...
    "input_names", ["power_W"; "ambient_temperature_K"], ...
    "output_names", ["temperature_K"; "heat_loss_W"; ...
                      "energy_input_J"; "temperature_rate_K_per_s"], ...
    "state_names", ["temperature_K"; "energy_input_J"], ...
    "radiation_linearized", logical(G_rad > 0), ...
    "radiation_reference_temperature_K", options.RadiationLinearizationTemperature_K, ...
    "notes", "Power is supplied by a Radia EM/VIM/FEM/ESIM reduced model.");

% The matrices are sufficient for Simulink. Add an ss object when the
% Control System Toolbox is available, but keep the plant usable without it.
if exist("ss", "file") == 2
    plant.sys = ss(A, B, Cout, D, Ts);
    plant.sys.InputName = cellstr(plant.input_names);
    plant.sys.OutputName = cellstr(plant.output_names);
    plant.sys.StateName = cellstr(plant.state_names);
else
    plant.sys = [];
end
end
