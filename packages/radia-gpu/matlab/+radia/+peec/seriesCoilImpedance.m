function result = seriesCoilImpedance(coil, frequency_Hz, options)
%SERIESCOILIMPEDANCE Evaluate a PEEC coil at one AC frequency.
%
% The HACApK matrix stores the uniform-current (DC) partial inductance.
% This function preserves that low-frequency limit and replaces its internal
% term with a shape-specific AC model: Dowell for a strict rectangle and the
% Bessel solution for a physical round wire. The resulting R_eff(f) and
% L_eff(f) are suitable for a self-consistent narrow-band resonance solve.
arguments
    coil (1,1) struct
    frequency_Hz (1,1) double {mustBePositive,mustBeFinite}
    options.DCProperties (1,1) struct = struct()
    options.RelativePermeability (1,1) double {mustBePositive,mustBeFinite} = 1
    options.InternalImpedanceModel (1,1) string = "auto"
end

crossSectionKind = "unknown";
if isfield(coil,"cross_section_kind")
    crossSectionKind = lower(string(coil.cross_section_kind));
end
resolvedModel = lower(options.InternalImpedanceModel);
if resolvedModel == "auto"
    if crossSectionKind == "circle"
        resolvedModel = "round-bessel";
    elseif crossSectionKind == "equivalent-circle"
        resolvedModel = "equivalent-round-bessel";
    else
        % The MATLAB PEEC schema's widths/heights describe rectangular bars.
        % Older structs without the shape field therefore remain compatible.
        resolvedModel = "rectangular-dowell";
    end
end
supportedModels = ["rectangular-dowell","round-bessel", ...
    "equivalent-round-bessel"];
if ~any(resolvedModel == supportedModels)
    error("radia:peec:InternalImpedanceModel", ...
        "Unsupported internal impedance model: %s", ...
        options.InternalImpedanceModel);
end
required = ["lengths_m","widths_m","heights_m", ...
    "conductivities_S_per_m"];
missing = required(~isfield(coil,required));
if ~isempty(missing)
    error("radia:peec:CoilSchema", ...
        "Coil geometry is missing field(s): %s",strjoin(missing,", "));
end
if isempty(fieldnames(options.DCProperties))
    dc = radia.peec.seriesCoilProperties(coil);
else
    dc = options.DCProperties;
end
if ~all(isfield(dc,["inductance_H","resistance_Ohm"]))
    error("radia:peec:DCProperties", ...
        "DCProperties must contain inductance_H and resistance_Ohm.");
end

lengths = double(coil.lengths_m(:));
widths = double(coil.widths_m(:));
heights = double(coil.heights_m(:));
conductivities = double(coil.conductivities_S_per_m(:));
if any([numel(widths),numel(heights),numel(conductivities)] ~= numel(lengths))
    error("radia:peec:CoilShape", ...
        "Conductor geometry arrays must have the same number of entries.");
end

mu0 = 4*pi*1e-7;
omega = 2*pi*frequency_Hz;
areas = widths.*heights;
equivalentRadii = sqrt(areas/pi);
dcPerLength = 1./(conductivities.*areas);

acPerLength = zeros(size(lengths));
dcInternalInductancePerLength = zeros(size(lengths));
for index = 1:numel(lengths)
    if resolvedModel == "rectangular-dowell"
        broadWidth = max(widths(index),heights(index));
        diffusionThickness = min(widths(index),heights(index));
        acPerLength(index) = rectangularDowellImpedancePerLength( ...
            broadWidth,diffusionThickness,conductivities(index),omega, ...
            options.RelativePermeability);
        dcInternalInductancePerLength(index) = ...
            mu0*options.RelativePermeability*diffusionThickness/ ...
            (12*broadWidth);
    else
        acPerLength(index) = cylinderImpedancePerLength( ...
            equivalentRadii(index),conductivities(index),omega, ...
            options.RelativePermeability);
        dcInternalInductancePerLength(index) = ...
            mu0*options.RelativePermeability/(8*pi);
    end
end
uniformInternal = dcPerLength + ...
    1i*omega*dcInternalInductancePerLength;
internalCorrection = sum(lengths.*(acPerLength-uniformInternal));
impedance = double(dc.resistance_Ohm) + ...
    1i*omega*double(dc.inductance_H) + internalCorrection;
effectiveInductance = imag(impedance)/omega;
effectiveResistance = real(impedance);
if ~(isfinite(effectiveInductance) && effectiveInductance > 0 && ...
        isfinite(effectiveResistance) && effectiveResistance > 0)
    error("radia:peec:ACImpedance", ...
        "The frequency-dependent coil impedance is not passive and finite.");
end

skinDepths = sqrt(2./(omega*mu0*options.RelativePermeability.*conductivities));
result = struct( ...
    "frequency_Hz",frequency_Hz, ...
    "omega_rad_per_s",omega, ...
    "impedance_Ohm",impedance, ...
    "resistance_Ohm",effectiveResistance, ...
    "inductance_H",effectiveInductance, ...
    "reactance_Ohm",imag(impedance), ...
    "dc_resistance_Ohm",double(dc.resistance_Ohm), ...
    "dc_inductance_H",double(dc.inductance_H), ...
    "internal_impedance_correction_Ohm",internalCorrection, ...
    "equivalent_radius_m",mean(equivalentRadii), ...
    "skin_depth_m",mean(skinDepths), ...
    "internal_impedance_model",resolvedModel, ...
    "cross_section_kind",crossSectionKind, ...
    "proximity_effect_included",false, ...
    "dc_properties",dc);
end

function impedance = rectangularDowellImpedancePerLength( ...
    width,thickness,sigma,omega,muR)
mu0 = 4*pi*1e-7;
dcResistance = 1/(sigma*width*thickness);
skinDepth = sqrt(2/(omega*mu0*muR*sigma));
q = (1+1i)*thickness/(2*skinDepth);
if abs(q) < 1e-3
    q2 = q*q;
    factor = 1 + q2/3 - q2*q2/45 + 2*q2*q2*q2/945;
elseif real(q) > 30
    factor = q;
else
    factor = q/tanh(q);
end
impedance = dcResistance*factor;
end

function impedance = cylinderImpedancePerLength(radius,sigma,omega,muR)
mu0 = 4*pi*1e-7;
dcResistance = 1/(pi*radius^2*sigma);
dcInternalInductance = mu0*muR/(8*pi);
ka = radius*sqrt(-1i*omega*mu0*muR*sigma);
if abs(ka) < 1e-3
    impedance = dcResistance + 1i*omega*dcInternalInductance;
    return
end
% Scaled Bessel functions share the same scale factor, so their ratio is
% stable for both the low- and high-frequency regimes.
j0 = besselj(0,ka,1);
j1 = besselj(1,ka,1);
if abs(j1) <= realmin
    error("radia:peec:BesselImpedance", ...
        "The cylindrical conductor Bessel denominator vanished.");
end
impedance = ka*j0/(2*pi*radius^2*sigma*j1);
end
