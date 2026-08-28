function result=validate_native_sobol_max_dimension()
%VALIDATE_NATIVE_SOBOL_MAX_DIMENSION Exercise the full Joe--Kuo table.
root=fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(root,"matlab"));
sampler=radia.optuna.QMCSampler(QMCType="sobol", ...
    Scramble=false,Seed=157);
tic
point=sampler.unitPoints(21201,1);
elapsed=toc;
assert(isequal(size(point),[1,21201]));
assert(all(point==0));
result=struct("schema","radia.validation.native-sobol-max.v1", ...
    "dimension",21201,"first_point_zero",true, ...
    "elapsed_seconds",elapsed);
disp(jsonencode(result,PrettyPrint=true));
end
