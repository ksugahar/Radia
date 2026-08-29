function value = ltspicePIDObjective(trial,modelFile)
%LTSPICEPIDOBJECTIVE Evaluate PID gains with the LTspice Simulink plant.
arguments
    trial (1,1) radia.optuna.Trial
    modelFile (1,1) string = ""
end

if strlength(modelFile) == 0
    modelFile = siblingPlantFile();
end

gains = trial.suggestVector(["Kp","Ki","Kd"], ...
    [0.05, 0.5, 1e-6], [8, 800, 5e-3], Log=[true,true,true]);
if ~isfile(modelFile)
    error("radia:simulink:PIDExampleModel", ...
        "The LTspice PID plant model does not exist: %s",modelFile);
end
runner = radia.optuna.SimulinkRunner(modelFile, ...
    ConfigureFcn=@(input,~) configureTrial(input,gains), ...
    ScoreFcn=@scoreTrial, StopTime="0.025");
value = runner.evaluate(trial);
end

function modelFile = siblingPlantFile()
harness = "radia_ltspice_pid_optuna";
if bdIsLoaded(harness)
    harnessFile = string(get_param(harness,"FileName"));
    if strlength(harnessFile) > 0
        modelFile = fullfile(fileparts(harnessFile), ...
            "radia_ltspice_pid_plant.slx");
        return
    end
end
modelFile = fullfile(radia.simulink.exampleDirectory(), ...
    "radia_ltspice_pid_plant.slx");
end

function input = configureTrial(input,gains)
input = input.setVariable("Kp",gains(1));
input = input.setVariable("Ki",gains(2));
input = input.setVariable("Kd",gains(3));
end

function value = scoreTrial(output,~)
errorSignal = output.get("pid_error");
controlSignal = output.get("pid_control");
if isempty(errorSignal) || isempty(controlSignal)
    error("radia:simulink:PIDExampleLog", ...
        "The PID example did not produce its error and control logs.");
end
tracking = mean(abs(errorSignal.Data(:)));
effort = mean(controlSignal.Data(:).^2);
value = tracking + 1e-3*effort;
end
