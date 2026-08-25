function copy_study(options)
%COPY_STUDY Copy one MAT-backed study to another storage file.
arguments
    options.from_study_name (1,1) string
    options.from_storage (1,1) string
    options.to_storage (1,1) string
    options.to_study_name (1,1) string = ""
end
source=radia.optuna.load_study(study_name=options.from_study_name, ...
    storage=options.from_storage);
targetName=options.to_study_name;
if strlength(targetName)==0, targetName=source.Name; end
if isfile(options.to_storage) || isfile(options.to_storage+".bak")
    error("radia:optuna:DuplicatedStudy", ...
        "Destination storage '%s' already exists.",options.to_storage);
end
loaded=builtin("load",source.StoragePath,"StudyData");
loaded.StudyData.Name=targetName;
folder=fileparts(options.to_storage);
if strlength(folder)==0
    folder=string(pwd);
elseif ~isfolder(folder)
    mkdir(folder);
end
temporary=string(tempname(folder))+".mat";
cleanup=onCleanup(@()removeTemporary(temporary));
StudyData=loaded.StudyData; %#ok<NASGU>
builtin("save",temporary,"StudyData","-mat");
radia.optuna.Study(StoragePath=temporary,AutoSave=false);
movefile(temporary,options.to_storage,"f");
clear cleanup
end

function removeTemporary(path)
if isfile(path), delete(path); end
end
