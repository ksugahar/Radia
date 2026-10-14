function setupSimulinkWorker(modelFolder)
%SETUPSIMULINKWORKER Make a file-backed trial model visible on a worker.
if strlength(modelFolder) > 0
    addpath(modelFolder);
end
end
