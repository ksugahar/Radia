function result=validate_binding_parity(pythonReport,outputPath)
%VALIDATE_BINDING_PARITY Same-input bitwise density through the two adapters.
% Model/optimizer RNG contracts remain the full upstream validation gates.
python=jsondecode(fileread(pythonReport));assert(python.ok);
package=fileparts(fileparts(mfilename('fullpath')));
oracle=jsondecode(fileread(fullfile(package,'oracles','optuna50_cuda_oracle.json')));
previous=radia.optuna.executionMode();
restore=onCleanup(@()radia.optuna.executionMode(previous));
radia.optuna.executionMode(string(python.backend));
count=0;
for i=1:numel(oracle.cases)
    entry=oracle.cases(i);
    if numel(entry.below)~=1,continue;end
    reference=python.cases(i);assert(strcmp(reference.name,entry.name));
    for side=["below","above"]
        d=entry.(side);x=entry.samples(:);
        if d.categorical
            actual=optuna_mex('optuna.parzen.log_pdf_categorical',x,d.weights,d.probabilities);
        else
            step=d.step;if isempty(step),step=NaN;end
            actual=optuna_mex('optuna.parzen.log_pdf_numerical',x,d.weights,d.mu,d.sigma, ...
                d.internal_low,d.internal_high,d.log,step);
        end
        expected=string(reference.(side+"_bits"));
        observed=string(num2hex(actual(:)));
        assert(isequal(observed(:),expected(:)), ...
            'radia:optuna:BitwiseParity','Same-input binding bits differ: %s/%s',entry.name,side);
        count=count+1;
    end
end
assert(count>0,'No same-input parity cases executed');
result=struct(ok=true,backend=string(python.backend),compared_density_arrays=count, ...
    scope="stored single-dimension density cases, bitwise; optimizer/model/RNG separate", ...
    archive_sha256=string(python.archive_sha256));
fid=fopen(outputPath,'w','n','UTF-8');assert(fid>=0);clean=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end
