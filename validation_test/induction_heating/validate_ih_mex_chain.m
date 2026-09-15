function report = validate_ih_mex_chain(inputDirectory, outputFile, caseNames)
% Consume real BEM/thermal assembly artifacts and compare native time stepping.
if nargin < 3, caseNames = ["solid","bored"]; end
before = radia.apiInfo(); cases = {};
for name = string(caseNames)
    data = jsondecode(fileread(fullfile(inputDirectory,name,"chain.json")));
    cfg = radia.simulink.validateIHNativeConfig(data.config);
    e = radia.internal.callMex('ih.eddy.create',cfg);
    ec = onCleanup(@() radia.internal.callMex('ih.eddy.destroy',e));
    t = radia.internal.callMex('ih.thermal.create',cfg);
    tc = onCleanup(@() radia.internal.callMex('ih.thermal.destroy',t));
    errorK = 0; relativeRise = 0; heatError = 0;
    for step = 1:size(data.reference_K,1)
        temp = radia.internal.callMex('ih.thermal.output',t);
        heat = radia.internal.callMex('ih.eddy.output',e,data.current_A,0,temp);
        heatError = max(heatError,norm(heat-data.current_A^2*cfg.heat_projection)/norm(heat));
        power = dot(cfg.heat_cell_weights,heat);
        assert(abs(power/(data.current_A^2*data.unit_power_W)-1)<1e-8);
        radia.internal.callMex('ih.thermal.update',t,heat,293.15,0);
        actual = radia.internal.callMex('ih.thermal.output',t);
        reference = data.reference_K(step,:).';
        if isfield(cfg,"temperature_constant_coefficients")
            sample = cfg.temperature_evaluation;
            evaluate = sparse(double(sample.rows)+1,double(sample.cols)+1,double(sample.values), ...
                sample.n_samples,cfg.n_temperature);
            delta = evaluate*(actual-reference);
            rise = evaluate*reference-293.15;
        else
            delta = actual-reference;
            rise = reference-293.15;
        end
        errorK = max(errorK,norm(delta,inf));
        relativeRise = max(relativeRise,norm(delta)/norm(rise));
    end
    assert(relativeRise < .02 && errorK < 1e-7 && heatError < 1e-12);
    clear tc ec
    cases{end+1} = struct('name',name,'nodes',data.nodes,'steps',size(data.reference_K,1), ...
        'unit_power_W',data.unit_power_W,'assembled_power_W',data.assembled_power_W, ...
        'max_temperature_difference_K',errorK,'max_relative_temperature_rise_error',relativeRise, ...
        'heat_vector_relative_error',heatError,'source_sha256',data.source_sha256,'scope',data.scope);
end
after = radia.apiInfo(); assert(after.ih_handle_count == before.ih_handle_count);
report = struct('passed',true,'scope','IH preassembled operator MEX parity; see per-case scope for source and geometry', ...
    'matlab_version',version,'cases',{cases},'native_build',jsondecode(fileread(which('radia_mex')+".build.json")), ...
    'handles_before',before.ih_handle_count,'handles_after',after.ih_handle_count);
fid=fopen(outputFile,'w');assert(fid>=0);fc=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
end
