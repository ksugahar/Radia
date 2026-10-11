function configureElectromagnetOptimizationSubsystem(blockPath)
%CONFIGUREELECTROMAGNETOPTIMIZATIONSUBSYSTEM Apply mask values internally.
arguments
    blockPath (1,1) string
end
if getSimulinkBlockHandle(blockPath) < 0
    error("radia:simulink:ElectromagnetTopologyBlock", ...
        "Block does not exist: %s",blockPath);
end
configFile = string(get_param(blockPath,"config_file"));
runRoot = string(get_param(blockPath,"run_root"));
timeout = string(get_param(blockPath,"timeout_s"));
python = string(get_param(blockPath,"python_executable"));
runnerExpression = string(get_param( ...
    blockPath,"radia_electromagnet_adjoint_runner"));
sampleTime = string(get_param(blockPath,"sample_time_s"));
historyCapacity = string(get_param(blockPath,"history_capacity"));
set_param(blockPath + "/Electromagnet Analysis","Parameters", ...
    "'em'," + configFile + "," + runRoot + "," + timeout + "," + python);
set_param(blockPath + "/Density Topology Optimization","Parameters", ...
    runnerExpression + "," + sampleTime + "," + historyCapacity);

runner = radia.simulink.resolveElectromagnetTopologyRunner( ...
    bdroot(blockPath),runnerExpression);
runner.setSolver(string(get_param(blockPath,"solver")));
end
