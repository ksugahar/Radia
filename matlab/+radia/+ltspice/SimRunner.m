classdef SimRunner < handle
    %SIMRUNNER MATLAB-native LTspice run manager.
    properties, OutputFolder (1,1) string="C:\temp\radia_ltspice_runs"; Executable (1,1) string=""; CommandLineSwitches (1,:) string=strings(1,0); end
    properties (SetAccess=private), CompletedRuns (1,:) cell={}; Tasks (1,:) cell={}; end
    methods
        function obj=SimRunner(options)
            arguments, options.OutputFolder (1,1) string="C:\temp\radia_ltspice_runs"; options.Executable (1,1) string=""; end
            obj.OutputFolder=options.OutputFolder; obj.Executable=options.Executable;
        end
        function result=runNow(obj,netlist,options)
            arguments, obj; netlist (1,1) string {mustBeFile}; options.Parameters (1,1) struct=struct(); options.RunName (1,1) string="run"; options.Timeout_s (1,1) double {mustBePositive}=600; end
            obj.requireSupportedSwitches(strings(1,0));
            result=radia.ltspice.run(netlist,Parameters=options.Parameters, ...
                Executable=obj.Executable, ...
                OutputDirectory=fullfile(obj.OutputFolder,options.RunName), ...
                Timeout_s=options.Timeout_s);
            result.raw=radia.ltspice.RawRead(result.raw_file); result.log_reader=radia.ltspice.LogReader(result.log_file);
            obj.CompletedRuns{end+1}=result;
        end
        function results=runMany(obj,netlist,parameterSets,options)
            arguments, obj; netlist (1,1) string {mustBeFile}; parameterSets (1,:) cell; options.UseParallel (1,1) logical=false, end
            results=cell(size(parameterSets));
            if options.UseParallel
                if isempty(ver('parallel')), error("radia:ltspice:ParallelToolbox","Parallel Computing Toolbox is required."); end
                root=obj.OutputFolder; exe=obj.Executable;
                parfor k=1:numel(parameterSets)
                    results{k}=radia.ltspice.run(netlist,Parameters=parameterSets{k},Executable=exe,OutputDirectory=fullfile(root,sprintf("run_%06d",k)));
                end
            else
                for k=1:numel(parameterSets), results{k}=obj.runNow(netlist,Parameters=parameterSets{k},RunName=sprintf("run_%06d",k)); end
            end
        end
        function task=run(obj,netlist,options)
            arguments,obj;netlist (1,1) string {mustBeFile};options.WaitResource (1,1) logical=true;options.Callback=[];options.Timeout (1,1) double {mustBePositive}=600;options.RunFilename (1,1) string="";options.Switches (1,:) string=strings(1,0);end
            if ~options.WaitResource
                error("radia:ltspice:UnsupportedCompatibility", ...
                    "WaitResource=false is not implemented by the MATLAB task scheduler.");
            end
            obj.requireSupportedSwitches(options.Switches);
            if strlength(options.RunFilename)==0,name="run_"+string(numel(obj.CompletedRuns)+1);else,[~,name]=fileparts(options.RunFilename);end
            output=fullfile(obj.OutputFolder,name);
            if isempty(ver('parallel')),result=localRunTask(string(netlist),obj.Executable,output,options.Timeout);task=radia.ltspice.RunTask(result,options.Callback);
            else,pool=gcp('nocreate');if isempty(pool),pool=parpool('Processes');end,future=parfeval(pool,@localRunTask,1,string(netlist),obj.Executable,output,options.Timeout);task=radia.ltspice.RunTask(future,options.Callback);end
            obj.Tasks{end+1}=task;
        end
        function [raw,log]=run_now(obj,netlist,options)
            arguments
                obj
                netlist (1,1) string {mustBeFile}
                options.switches (1,:) string = strings(1,0)
                options.run_filename (1,1) string = ""
                options.timeout (1,1) double {mustBePositive} = 600
                options.exe_log (1,1) logical = false
            end
            if options.exe_log
                error("radia:ltspice:UnsupportedCompatibility", ...
                    "exe_log is not implemented by the MATLAB runner.");
            end
            obj.requireSupportedSwitches(options.switches);
            if strlength(options.run_filename)==0
                name="run";
            else
                [~,name]=fileparts(options.run_filename);
            end
            result=obj.runNow(netlist,RunName=name,Timeout_s=options.timeout);
            raw=string(result.raw_file);log=string(result.log_file);
        end
        function addCommandLineSwitch(obj,switchName,path),if nargin<3,path="";end,obj.CommandLineSwitches(end+1)=string(switchName)+string(path);end
        function clearCommandLineSwitches(obj),obj.CommandLineSwitches=strings(1,0);end
        function n=activeThreads(obj),n=sum(cellfun(@(x)x.isAlive(),obj.Tasks));end
        function answer=waitCompletion(obj,timeout,varargin)
            if nargin<2||isempty(timeout),timeout=Inf;end,t=tic;answer=true;for k=1:numel(obj.Tasks),remaining=max(0,timeout-toc(t));if ~obj.Tasks{k}.wait(remaining),answer=false;if ~isempty(varargin)&&varargin{1},obj.killAllSpice();end,break,end,end
        end
        function info=simInfo(obj),info=struct('completed',sum(cellfun(@(x)~x.isAlive(),obj.Tasks)),'active',obj.activeThreads(),'output_folder',obj.OutputFolder);end
        function runs=tasks(obj),runs=obj.Tasks;end
        function killAllSpice(obj),for k=1:numel(obj.Tasks),obj.Tasks{k}.cancel();end,end
        function cleanupFiles(obj),if isfolder(obj.OutputFolder),files=dir(fullfile(obj.OutputFolder,'run_*'));for k=1:numel(files),if files(k).isdir,rmdir(fullfile(files(k).folder,files(k).name),'s');end,end,end,end
        function path=createNetlist(~,ascFile,varargin)
            if ~isempty(varargin)
                error("radia:ltspice:UnsupportedCompatibility", ...
                    "create_netlist compatibility options are not implemented yet.");
            end
            result=radia.ltspice.schematicToNetlist(string(ascFile));
            path=string(result.netlist);
        end
        function answer=createRawFileWith(~,rawFilename,saveNames,conditions)
            answer=false; %#ok<NASGU>
            context=sprintf("destination=%s, traces=%d, conditions=%d", ...
                string(rawFilename),numel(saveNames),~isempty(conditions));
            error("radia:ltspice:UnsupportedCompatibility", ...
                "create_raw_file_with requires multi-run trace tagging and " + ...
                "condition filtering, which are not implemented yet (%s).", ...
                 context);
        end
        function exportSimLog(obj,path),rows=table();for k=1:numel(obj.Tasks),task=obj.Tasks{k};if ~task.isAlive(),rows=[rows;table(k,task.Status,task.RawFile,task.LogFile,'VariableNames',{'Task','Status','RawFile','LogFile'})];end,end,writetable(rows,path);end %#ok<AGROW>
        function setSimulator(obj,simulator),if isstruct(simulator)&&isfield(simulator,'executable'),obj.Executable=string(simulator.executable);elseif isstring(simulator)||ischar(simulator),obj.Executable=string(simulator);end,end
        function updateCompleted(~),end
        function args=validateCallbackArgs(~,callback,args),if isempty(callback),args=[];elseif nargin<3||isempty(args),args={};elseif ~iscell(args)&&~isstruct(args),error("radia:ltspice:CallbackArgs","Callback arguments must be a cell or struct.");end,end
        function add_command_line_switch(obj,varargin),obj.addCommandLineSwitch(varargin{:});end
        function clear_command_line_switches(obj),obj.clearCommandLineSwitches();end
        function x=active_threads(obj),x=obj.activeThreads();end
        function x=wait_completion(obj,varargin),x=obj.waitCompletion(varargin{:});end
        function x=sim_info(obj),x=obj.simInfo();end
        function kill_all_spice(obj),obj.killAllSpice();end
        function kill_all_ltspice(obj),obj.killAllSpice();end
        function cleanup_files(obj),obj.cleanupFiles();end
        function file_cleanup(obj),obj.cleanupFiles();end
        function x=create_netlist(obj,varargin),x=obj.createNetlist(varargin{:});end
        function x=create_raw_file_with(obj,varargin),x=obj.createRawFileWith(varargin{:});end
        function export_sim_log(obj,varargin),obj.exportSimLog(varargin{:});end
        function set_simulator(obj,varargin),obj.setSimulator(varargin{:});end
        function update_completed(obj),obj.updateCompleted();end
        function x=validate_callback_args(obj,varargin),x=obj.validateCallbackArgs(varargin{:});end
    end
    methods (Access=private)
        function requireSupportedSwitches(obj,overrides)
            if ~isempty(overrides)||~isempty(obj.CommandLineSwitches)
                error("radia:ltspice:UnsupportedCompatibility", ...
                    "LTspice command-line switches are not wired into the " + ...
                    "MATLAB process launcher yet.");
            end
        end
    end
end
function result=localRunTask(netlist,executable,output,timeout)
result=radia.ltspice.run(netlist,Executable=executable, ...
    OutputDirectory=output,Timeout_s=timeout);
result.raw=radia.ltspice.RawRead(result.raw_file);
result.log_reader=radia.ltspice.LogReader(result.log_file);
end
