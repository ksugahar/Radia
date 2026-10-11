function info = taskmanagerProbe(n)
%TASKMANAGERPROBE Exercise NGSolve TaskManager inside MATLAB.

if nargin == 0
    info = radia.internal.callMex('taskmanager.probe');
else
    info = radia.internal.callMex('taskmanager.probe', double(n));
end
end
