function tests = test_peec_coil_geometry_resonance
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot,"matlab"));
testCase.TestData.netlist = fullfile(repoRoot,"matlab","samples", ...
    "ltspice_rlc_resonance_ac.cir");
testCase.TestData.outputRoot = "C:\temp\radia_peec_coil_geometry_test";
end

function testCircularTurnStackHasClosedOrientedTurns(testCase)
coil = radia.peec.circularTurnStack(0.03,12,SegmentsPerTurn=24);
verifyEqual(testCase,coil.segment_count,288);
verifySize(testCase,coil.centers_m,[288 3]);
verifySize(testCase,coil.polylines_m,[25 3 12]);
verifyEqual(testCase,vecnorm(coil.directions,2,2), ...
    ones(288,1),"AbsTol",1e-13);
for turn = 1:coil.turn_count
    points = coil.polylines_m(:,:,turn);
    verifyEqual(testCase,points(1,:),points(end,:),"AbsTol",1e-14);
end
end

function testPEECInductanceConvergesAndMatchesWheelerScale(testCase)
coarse = radia.peec.circularTurnStack(0.03,12,SegmentsPerTurn=24);
refined = radia.peec.circularTurnStack(0.03,12,SegmentsPerTurn=48);
coarseProperties = radia.peec.seriesCoilProperties(coarse);
refinedProperties = radia.peec.seriesCoilProperties(refined);
relativeRefinement = abs(coarseProperties.inductance_H- ...
    refinedProperties.inductance_H)/refinedProperties.inductance_H;
verifyLessThan(testCase,relativeRefinement,0.01);
verifyLessThan(testCase,coarseProperties.wheeler_relative_difference,0.12);
expectedResistance = coarse.wire_length_m/( ...
    coarse.conductivity_S_per_m*coarse.wire_width_m*coarse.wire_height_m);
verifyEqual(testCase,coarseProperties.resistance_Ohm, ...
    expectedResistance,"RelTol",1e-13);
end

function testFrequencyDependentImpedancePreservesDCAndAddsSkinEffect(testCase)
coil = radia.peec.circularTurnStack(0.03,12,SegmentsPerTurn=24);
dc = radia.peec.seriesCoilProperties(coil);
low = radia.peec.seriesCoilImpedance(coil,1,DCProperties=dc);
high = radia.peec.seriesCoilImpedance(coil,1e6,DCProperties=dc);
verifyEqual(testCase,low.inductance_H,dc.inductance_H,"RelTol",1e-9);
verifyEqual(testCase,low.resistance_Ohm,dc.resistance_Ohm,"RelTol",1e-9);
verifyGreaterThan(testCase,high.resistance_Ohm,4*dc.resistance_Ohm);
verifyLessThan(testCase,high.inductance_H,dc.inductance_H);
verifyFalse(testCase,high.proximity_effect_included);
verifyEqual(testCase,high.internal_impedance_model, ...
    "equivalent-round-bessel");
end

function testBesselCorrectionApproachesSurfaceImpedance(testCase)
frequency = 1e9;
coil = radia.peec.circularTurnStack(0.03,12,SegmentsPerTurn=24);
dc = radia.peec.seriesCoilProperties(coil);
ac = radia.peec.seriesCoilImpedance(coil,frequency,DCProperties=dc);
mu0 = 4*pi*1e-7;
omega = 2*pi*frequency;
area = coil.wire_width_m*coil.wire_height_m;
radius = sqrt(area/pi);
sigma = coil.conductivity_S_per_m;
skinDepth = sqrt(2/(omega*mu0*sigma));
surfacePerLength = (1+1i)/(sigma*skinDepth*2*pi*radius);
uniformPerLength = 1/(sigma*area) + 1i*omega*mu0/(8*pi);
expectedCorrection = coil.wire_length_m* ...
    (surfacePerLength-uniformPerLength);
verifyEqual(testCase,ac.internal_impedance_correction_Ohm, ...
    expectedCorrection,"RelTol",5e-3);
end

function testSelfConsistentResonanceUsesInductanceAtTheRoot(testCase)
capacitance = 3.3e-3;
coil = radia.peec.circularTurnStack(0.03,12,SegmentsPerTurn=24);
dc = radia.peec.seriesCoilProperties(coil);
root = radia.peec.solveSeriesResonance( ...
    coil,capacitance,DCProperties=dc);
verifyTrue(testCase,root.converged);
verifyLessThan(testCase,root.relative_reactance_residual,1e-8);
verifyEqual(testCase,root.coil_reactance_Ohm, ...
    -root.capacitor_reactance_Ohm,"RelTol",1e-9);
verifyEqual(testCase,root.inductance_H, ...
    imag(root.impedance_Ohm)/(2*pi*root.frequency_Hz),"RelTol",1e-13);

runResult = radia.ltspice.run(testCase.TestData.netlist, ...
    Parameters=struct( ...
        "Rval",root.resistance_Ohm, ...
        "Lval",root.inductance_H, ...
        "Cval",capacitance), ...
    OutputDirectory=fullfile(testCase.TestData.outputRoot,"ac"), ...
    RawFormat="binary");
[measuredFrequency,details] = ...
    radia.ltspice.measureSeriesResonance(runResult);
verifyEqual(testCase,measuredFrequency,root.frequency_Hz,"RelTol",5e-5);
verifyEqual(testCase,details.input_resistance_at_resonance_Ohm, ...
    root.resistance_Ohm,"RelTol",1e-10);
verifyLessThan(testCase,details.zero_vs_peak_relative_difference,5e-4);
end
