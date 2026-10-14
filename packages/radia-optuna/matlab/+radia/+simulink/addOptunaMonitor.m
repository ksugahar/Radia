function monitorPath = addOptunaMonitor(parent, name, position)
%ADDOPTUNAMONITOR Add browser-free Scope and Pareto visualization sinks.
arguments
    parent (1,1) string
    name (1,1) string = "Optuna Monitor"
    position (1,4) double = [420 70 650 210]
end
monitorPath = parent + "/" + name;
add_block("simulink/Ports & Subsystems/Subsystem", monitorPath, ...
    Position=position);
delete_line(monitorPath, "In1/1", "Out1/1");
delete_block(monitorPath + "/In1");
delete_block(monitorPath + "/Out1");

labels = ["best","last","trials","status","pareto_x","pareto_y", ...
    "pareto_count","best_updated","pareto_revision"];
for index = 1:numel(labels)
    y = 25 + (index - 1) * 35;
    add_block("simulink/Ports & Subsystems/In1", ...
        monitorPath + "/" + labels(index), Port=string(index), ...
        Position=[25 y 55 y+20]);
end

add_block("simulink/Sinks/Scope", monitorPath + "/Optimization History", ...
    NumInputPorts="4", Position=[145 30 275 125]);
add_block("simulink/Sinks/XY Graph", monitorPath + "/Pareto Front", ...
    Position=[145 165 275 235]);
add_block("simulink/Sinks/Scope", monitorPath + "/Update Events", ...
    NumInputPorts="3", Position=[145 265 275 335]);
add_line(monitorPath, "best/1", "Optimization History/1");
add_line(monitorPath, "last/1", "Optimization History/2");
add_line(monitorPath, "trials/1", "Optimization History/3");
add_line(monitorPath, "status/1", "Optimization History/4");
add_line(monitorPath, "pareto_x/1", "Pareto Front/1");
add_line(monitorPath, "pareto_y/1", "Pareto Front/2");
add_line(monitorPath, "pareto_count/1", "Update Events/1");
add_line(monitorPath, "best_updated/1", "Update Events/2");
add_line(monitorPath, "pareto_revision/1", "Update Events/3");

mask = Simulink.Mask.create(monitorPath);
mask.Description = "Simulink-native Optuna history, best-update, and Pareto-front monitor. No browser is used.";
display = "disp('Optuna Monitor');";
for index = 1:numel(labels)
    display = display + "port_label('input'," + index + ",'" + labels(index) + "');";
end
mask.Display = display;
set_param(monitorPath, "UserData", struct( ...
    "visualization", "simulink-scope-xy", ...
    "browser_required", false, ...
    "best_update_policy", "pulse-on-primary-incumbent-change"), ...
    "UserDataPersistent", "on");
end
