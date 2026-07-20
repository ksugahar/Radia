function libraryPath=buildLibrary(options)
%BUILDLIBRARY Build the single distributable Radia Simulink block library.
arguments, options.OutputDirectory (1,1) string=""; end
matlabRoot=fileparts(fileparts(fileparts(mfilename("fullpath"))));
if strlength(options.OutputDirectory)==0,options.OutputDirectory=matlabRoot;end
if ~isfolder(options.OutputDirectory),mkdir(options.OutputDirectory);end
name="radia_simulink_library"; libraryPath=fullfile(options.OutputDirectory,name+".slx");
if bdIsLoaded(name),close_system(name,0);end
new_system(name,"Library");
applications=addEmptySubsystem(name,"Applications",[70 40 310 265]);
addApplicationBlock(applications,"Electromagnet","em",[45 30 230 75]);
addApplicationBlock(applications,"PCB PEEC","pcb",[45 90 230 135]);
addApplicationBlock(applications,"Motor","motor",[45 150 230 195]);
addApplicationBlock(applications,"Stream Function","streamfunction",[45 210 230 255]);
addApplicationBlock(applications,"Induction Heating","ih",[45 270 230 315]);

addEmptySubsystem(name,"LTspice",[70 305 310 405]);
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function",name+"/LTspice/LTspice Circuit", ...
 "FunctionName","radia_ltspice_sfun","Parameters","'', {'control'}, {'V(out)'}, 1e-3, 'C:\temp\radia_ltspice_block', inf, 300, ''",Position=[45 35 200 85]);
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function",name+"/LTspice/Hysteretic LTspice Plant", ...
 "FunctionName","radia_hysteretic_ltspice_sfun","Parameters","'', 1e-3, 'C:\temp\radia_hysteretic_ltspice_block'",Position=[45 105 200 155]);
addEmptySubsystem(name,"Optimization",[70 445 310 535]);
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function",name+"/Optimization/Optuna Optimization", ...
 "FunctionName","radia_optuna_sfun","Parameters","'', 20, 'minimize', '', 1, true",Position=[45 35 200 85]);
set_param(name,"Lock","on"); save_system(name,libraryPath); close_system(name,0);
end

function path=addEmptySubsystem(parent,name,position)
path=parent+"/"+name;
add_block("simulink/Ports & Subsystems/Subsystem",path,Position=position);
delete_line(path,"In1/1","Out1/1");
delete_block(path+"/In1");
delete_block(path+"/Out1");
end

function addApplicationBlock(parent,label,application,position)
path=parent+"/"+label;
parameters="'"+application+"', config_file, run_root, timeout_s, python_executable";
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function",path, ...
    FunctionName="radia_application_sfun", ...
    Parameters="'"+application+"', '', 'C:\temp\radia_simulink', 3600, 'python'", ...
    Position=position);
mask=Simulink.Mask.create(path);
mask.Description=label+" batch analysis. Create settings with radia.simulink.writeApplicationConfig; a rising trigger executes the validated Radia headless CLI once.";
mask.addParameter(Type="edit",Name="config_file", ...
    Prompt="Configuration JSON",Value="''",Evaluate="on");
mask.addParameter(Type="edit",Name="run_root", ...
    Prompt="Run artifact root", ...
    Value="'C:\temp\radia_simulink'",Evaluate="on");
mask.addParameter(Type="edit",Name="timeout_s", ...
    Prompt="Timeout (s)",Value="3600",Evaluate="on");
mask.addParameter(Type="edit",Name="python_executable", ...
    Prompt="Python executable",Value="'python'",Evaluate="on");
mask.Display="disp('"+label+"');" + ...
    "port_label('input',1,'run');" + ...
    "port_label('output',1,'status');" + ...
    "port_label('output',2,'primary');" + ...
    "port_label('output',3,'elapsed_s');";
set_param(path,"Parameters",parameters);
end
