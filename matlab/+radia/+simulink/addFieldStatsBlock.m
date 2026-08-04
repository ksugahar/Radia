function blockPath = addFieldStatsBlock(parentPath, options)
%ADDFIELDSTATSBLOCK Place the [min mean max] field-vector reduction block.
%   blockPath = addFieldStatsBlock(parentPath) adds a one-in / one-out
%   subsystem that reduces an N-wide field vector to the 3-wide
%   [min mean max] signal with named channels.  The computation lives
%   in readable .m code (radia.simulink.fieldStatsSFunction behind the
%   radia_field_stats_sfun Level-2 wrapper) rather than a
%   primitive-block web.  Wire it between a distributed-field port (temperature, heat
%   density) and a scope: a real IH configuration carries thousands of
%   DOFs (measured 2026-08-04: 3122 temperature DOFs) and a scope fed
%   the raw vector draws thousands of overlapping lines.  The mean is
%   the arithmetic mean over DOFs -- a display aid; volume-weighted
%   means stay owned by the result artifacts.
%
%   Used by buildIHNativeModel for the radia_ih scopes and published in
%   the Radia library under Utilities/Field Stats for user models that
%   wire the library IH block's q / T ports to their own scopes.

arguments
    parentPath (1,1) string
    options.BlockName (1,1) string = "Field Stats"
    options.Position (1,4) double = [45 35 285 105]
end

blockPath = parentPath + "/" + options.BlockName;
add_block("simulink/Ports & Subsystems/Subsystem", blockPath, ...
    Position=options.Position);
delete_line(blockPath, "In1/1", "Out1/1");
set_param(blockPath + "/In1", "Position", [40 128 70 142]);
set_param(blockPath + "/Out1", "Position", [420 138 450 152]);
% The reduction itself is readable .m code (fieldStatsSFunction via the
% radia_field_stats_sfun wrapper); the Mux only routes and names the
% three channels for scope legends.
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
    blockPath + "/Stats", FunctionName="radia_field_stats_sfun", ...
    Position=[170 80 300 210]);
add_block("simulink/Signal Routing/Mux", blockPath + "/Mux", ...
    Inputs="3", Position=[340 40 350 250]);
add_line(blockPath, "In1/1", "Stats/1", "autorouting", "smart");
set_param(add_line(blockPath, "Stats/1", "Mux/1", "autorouting", "smart"), ...
    "Name", "min");
set_param(add_line(blockPath, "Stats/2", "Mux/2", "autorouting", "smart"), ...
    "Name", "mean");
set_param(add_line(blockPath, "Stats/3", "Mux/3", "autorouting", "smart"), ...
    "Name", "max");
add_line(blockPath, "Mux/1", "Out1/1", "autorouting", "smart");

set_param(blockPath, "Tag", "RadiaFieldStats");
mask = Simulink.Mask.create(blockPath);
mask.Description = "Reduces an N-wide field vector (temperature, " + ...
    "heat density) to the 3-wide [min mean max] signal for scopes " + ...
    "and logging.  Base Simulink blocks only; the mean is the " + ...
    "arithmetic mean over DOFs (display aid -- volume-weighted " + ...
    "means stay owned by the result artifacts).";
mask.Display = "disp('[min mean max]');";
end
