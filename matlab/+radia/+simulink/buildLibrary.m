function libraryPath=buildLibrary(options)
%BUILDLIBRARY Build the single distributable Radia Simulink block library.
arguments, options.OutputDirectory (1,1) string=""; end
matlabRoot=fileparts(fileparts(fileparts(mfilename("fullpath"))));
if strlength(options.OutputDirectory)==0,options.OutputDirectory=matlabRoot;end
if ~isfolder(options.OutputDirectory),mkdir(options.OutputDirectory);end
name="radia_simulink_library"; libraryPath=fullfile(options.OutputDirectory,name+".slx");
if bdIsLoaded(name),close_system(name,0);end
new_system(name,"Library");
add_block("simulink/Ports & Subsystems/Subsystem",name+"/LTspice",Position=[70 45 260 135]);
delete_line(name+"/LTspice","In1/1","Out1/1"); delete_block(name+"/LTspice/In1"); delete_block(name+"/LTspice/Out1");
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function",name+"/LTspice/LTspice Circuit", ...
 "FunctionName","radia_ltspice_sfun","Parameters","'', {'control'}, {'V(out)'}, 1e-3, 'C:\temp\radia_ltspice_block', inf, 300, ''",Position=[45 35 200 85]);
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function",name+"/LTspice/Hysteretic LTspice Plant", ...
 "FunctionName","radia_hysteretic_ltspice_sfun","Parameters","'', 1e-3, 'C:\temp\radia_hysteretic_ltspice_block'",Position=[45 105 200 155]);
add_block("simulink/Ports & Subsystems/Subsystem",name+"/Optimization",Position=[70 175 260 265]);
delete_line(name+"/Optimization","In1/1","Out1/1"); delete_block(name+"/Optimization/In1"); delete_block(name+"/Optimization/Out1");
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function",name+"/Optimization/Optuna Optimization", ...
 "FunctionName","radia_optuna_sfun","Parameters","'', 20, 'minimize', '', 1, true",Position=[45 35 200 85]);
set_param(name,"Lock","on"); save_system(name,libraryPath); close_system(name,0);
end
