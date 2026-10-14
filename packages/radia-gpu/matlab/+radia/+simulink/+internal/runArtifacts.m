function out=runArtifacts(action,varargin)
%RUNARTIFACTS Bounded run-artifact retention for the LTspice coupling blocks.
%   The blocks write one heavy artifact folder per simulation step. Keeping
%   every folder grows without bound during a long run; deleting all of them
%   at Terminate destroys the evidence a failed or suspicious run needs.
%   This helper splits the two concerns:
%
%     * an always-on single-file run log, written to a deterministic path and
%       overwritten once per run, that survives regardless of exit code;
%     * a bounded ring of the most recent step folders, pruned as the run
%       advances so peak disk usage stays flat.
%
%   ACTIONS
%     path = runArtifacts("logpath",root,blockHandle)
%     runArtifacts("begin",logPath,headerText)
%     runArtifacts("append",logPath,lineText)
%     kept = runArtifacts("prune",kept,keepCount)
%     runArtifacts("finish",folder,failed)   % keeps folder when failed
arguments
    action (1,1) string
end
arguments (Repeating)
    varargin
end
out=[];
switch action
    case "logpath"
        root=string(varargin{1});name=string(getfullname(varargin{2}));
        out=fullfile(root,regexprep(name,'[^\w]','_')+".log");
    case "begin"
        writeLine(string(varargin{1}),string(varargin{2}),'w');
    case "append"
        writeLine(string(varargin{1}),string(varargin{2}),'a');
    case "prune"
        kept=string(varargin{1});keepCount=varargin{2};
        while numel(kept)>keepCount
            removeFolder(kept(1));kept(1)=[];
        end
        out=kept;
    case "finish"
        if ~varargin{2},removeFolder(string(varargin{1}));end
    otherwise
        error("radia:simulink:RunArtifactAction","Unsupported run-artifact action: %s",action);
end
end

function writeLine(path,text,mode)
folder=fileparts(path);
if strlength(folder)>0&&~isfolder(folder),mkdir(folder);end
file=fopen(path,mode);
if file<0,error("radia:simulink:RunArtifactLog","Could not write the run log %s.",path);end
cleanup=onCleanup(@()fclose(file));
fprintf(file,'%s\n',text);
clear cleanup
end

function removeFolder(folder)
if strlength(folder)>0&&isfolder(folder)
    try
        rmdir(folder,'s');
    catch cause
        warning("radia:simulink:RunArtifactCleanup","Could not remove run folder %s: %s",folder,cause.message);
    end
end
end
