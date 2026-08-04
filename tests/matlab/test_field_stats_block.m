function tests = test_field_stats_block
%TEST_FIELD_STATS_BLOCK Numeric lock of the .m-based field reduction.
%   The Field Stats block computes [min mean max] in readable .m code
%   (radia.simulink.fieldStatsSFunction behind radia_field_stats_sfun)
%   instead of a primitive-block web; this locks the numbers and the
%   dynamic input sizing end-to-end through a real simulation.
tests = functiontests(localfunctions);
end

function testReducesVectorToMinMeanMax(testCase)
model = "field_stats_" + erase(string(java.util.UUID.randomUUID), "-");
new_system(model);
closer = onCleanup(@() close_system(model, 0)); %#ok<NASGU>
add_block("simulink/Sources/Constant", model + "/Field", ...
    Value="[4 1 2 5]", Position=[40 40 100 80]);
radia.simulink.addFieldStatsBlock(model, Position=[160 30 260 90]);
add_block("simulink/Ports & Subsystems/Out1", model + "/stats", ...
    Position=[320 50 350 70]);
add_line(model, "Field/1", "Field Stats/1", "autorouting", "smart");
add_line(model, "Field Stats/1", "stats/1", "autorouting", "smart");
set_param(model, "StopTime", "0", "SaveOutput", "on", ...
    "OutputSaveName", "yout");
out = sim(model);
y = out.yout;
if isa(y, "Simulink.SimulationData.Dataset")
    data = y.getElement(1).Values.Data;
else
    data = y;
end
verifyEqual(testCase, data(end, :), [1 3 5], "AbsTol", 1e-12);
end

function testScalarInputDegenerates(testCase)
model = "field_stats_" + erase(string(java.util.UUID.randomUUID), "-");
new_system(model);
closer = onCleanup(@() close_system(model, 0)); %#ok<NASGU>
add_block("simulink/Sources/Constant", model + "/Field", ...
    Value="7", Position=[40 40 100 80]);
radia.simulink.addFieldStatsBlock(model, Position=[160 30 260 90]);
add_block("simulink/Ports & Subsystems/Out1", model + "/stats", ...
    Position=[320 50 350 70]);
add_line(model, "Field/1", "Field Stats/1", "autorouting", "smart");
add_line(model, "Field Stats/1", "stats/1", "autorouting", "smart");
set_param(model, "StopTime", "0", "SaveOutput", "on", ...
    "OutputSaveName", "yout");
out = sim(model);
y = out.yout;
if isa(y, "Simulink.SimulationData.Dataset")
    data = y.getElement(1).Values.Data;
else
    data = y;
end
verifyEqual(testCase, data(end, :), [7 7 7], "AbsTol", 1e-12);
end
