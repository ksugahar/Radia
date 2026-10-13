function receipt=smoke_published_modes(installedMatlabPath,outputPath)
arguments
    installedMatlabPath (1,1) string
    outputPath (1,1) string
end
addpath(installedMatlabPath,'-begin');
assert(startsWith(string(which('optuna_mex')),installedMatlabPath+filesep));
proposals=cell(1,2);values=cell(1,2);statistics=cell(1,2);
for index=1:2
    modes=["cpu","gpu"];radia.optuna.executionMode(modes(index));
    before=optuna_mex('api.info');
    study=radia.optuna.Study(Sampler=radia.optuna.TPESampler(Seed=37,Multivariate=true),AutoSave=false);
    x=zeros(40,2);y=zeros(40,1);
    for k=1:40
        trial=study.ask();x(k,1)=trial.suggest_float('x1',-5,10);x(k,2)=trial.suggest_float('x2',0,15);
        y(k)=(x(k,2)-5.1/(4*pi^2)*x(k,1)^2+5/pi*x(k,1)-6)^2+10*(1-1/(8*pi))*cos(x(k,1))+10;
        study.tell(trial,y(k));
    end
    proposals{index}=x;values{index}=y;after=optuna_mex('api.info');
    if index==1
        assert(after.density_launches==before.density_launches && after.host_to_device_bytes==before.host_to_device_bytes);
        assert(after.cpu_density_calls>before.cpu_density_calls);
    else
        assert(after.model_build_launches>before.model_build_launches && after.density_launches>before.density_launches);
        assert(after.cpu_preparation_calls==before.cpu_preparation_calls);
    end
    statistics{index}=struct('before',before,'after',after);clear trial study
end
assert(isequal(proposals{1},proposals{2}) && isequal(values{1},values{2}));
receipt=struct('ok',true,'trials_per_mode',40,'proposals_exact',true,'values_exact',true, ...
    'proposals',proposals{1},'values',values{1},'cpu',statistics{1},'gpu',statistics{2});
fid=fopen(outputPath,'w');assert(fid~=-1);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(receipt,PrettyPrint=true));
radia.optuna.executionMode('cpu');clear optuna_mex
end
