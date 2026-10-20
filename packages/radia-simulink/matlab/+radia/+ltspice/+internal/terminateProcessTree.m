function terminated=terminateProcessTree(process)
%TERMINATEPROCESSTREE Stop an owned Windows process and its descendants.
terminated=true;
try
    if process.HasExited,return,end
    if ispc
        system(sprintf('taskkill /T /F /PID %d >NUL 2>&1',double(process.Id)));
    end
    if ~process.HasExited,process.Kill();end
    process.WaitForExit(5000);
    terminated=process.HasExited;
catch
    terminated=false;
end
end
