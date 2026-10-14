function config = makeNonlinearReactorDemoConfig(options)
%MAKENONLINEARREACTORDEMOCONFIG Build the deterministic toroidal demo model.
%   The single retained mode is the divergence-free toroidal HDiv-MMM mode.
%   Quadrature samples preserve the radial 1/r field variation and provide
%   the distributed B output used by the public Simulink demonstration.

arguments
    options.MajorRadius_m (1,1) double {mustBeFinite,mustBePositive} = 0.05
    options.MinorRadius_m (1,1) double {mustBeFinite,mustBePositive} = 0.01
    options.Turns (1,1) double {mustBeInteger,mustBePositive} = 80
    options.WindingResistance_Ohm (1,1) double {mustBeFinite,mustBeNonnegative} = 0.30
    options.SampleTime_s (1,1) double {mustBeFinite,mustBePositive} = 1.0e-4
end
R = options.MajorRadius_m;
a = options.MinorRadius_m;
if a >= 0.5*R
    error("radia:simulink:NonlinearReactorGeometry", ...
        "MinorRadius_m must be less than half MajorRadius_m.");
end

nRadial = 4;
nSectionAngle = 12;
nAzimuth = 6;
nSamples = nRadial*nSectionAngle*nAzimuth;
points = zeros(nSamples,3);
modes = zeros(1,nSamples,3);
weights = zeros(nSamples,1);
areaWeight = pi*a^2/(nRadial*nSectionAngle);
excitation = 0;
airInductance = 0;
sample = 0;
for radialIndex = 1:nRadial
    rho = a*sqrt((radialIndex-0.5)/nRadial);
    for sectionIndex = 1:nSectionAngle
        alpha = 2*pi*(sectionIndex-0.5)/nSectionAngle;
        cylindricalRadius = R+rho*cos(alpha);
        z = rho*sin(alpha);
        for azimuthIndex = 1:nAzimuth
            phi = 2*pi*(azimuthIndex-1)/nAzimuth;
            sample = sample+1;
            tangent = [-sin(phi),cos(phi),0];
            points(sample,:) = [cylindricalRadius*cos(phi), ...
                cylindricalRadius*sin(phi),z];
            modes(1,sample,:) = (R/cylindricalRadius)*tangent;
            weights(sample) = (2*pi/nAzimuth)*cylindricalRadius*areaWeight;
            hPerAmp = options.Turns/(2*pi*cylindricalRadius);
            excitation = excitation+weights(sample)* ...
                dot(reshape(modes(1,sample,:),1,3),hPerAmp*tangent);
            airInductance = airInductance+ ...
                4*pi*1.0e-7*options.Turns^2/(2*pi)* ...
                areaWeight/cylindricalRadius/nAzimuth;
        end
    end
end

bhCurve = [ ...
    0,0; 50,0.55; 100,0.95; 200,1.30; 500,1.55; ...
    1000,1.68; 2000,1.78; 5000,1.90; 10000,1.98; ...
    50000,2.08; 200000,2.28];
config = radia.simulink.makeNonlinearReactorConfig( ...
    0,modes,weights,excitation,bhCurve, ...
    SamplePoints_m=points,AirInductance_H=airInductance, ...
    WindingResistance_Ohm=options.WindingResistance_Ohm, ...
    SampleTime_s=options.SampleTime_s,Source="toroidal-hdiv-mmm-mode");
config.geometry = struct("kind","toroid","major_radius_m",R, ...
    "minor_radius_m",a,"turns",options.Turns);
end
