function report = validate_ih_mex_chain(inputDirectory, outputFile)
% Consume real BEM/thermal assembly artifacts and compare native time stepping.
before = radia.apiInfo(); cases = {};
for name = ["solid","bored"]
    data = jsondecode(fileread(fullfile(inputDirectory,name,"chain.json")));
    cfg = data.config;
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
        errorK = max(errorK,norm(actual-reference,inf));
        relativeRise = max(relativeRise,norm(actual-reference)/norm(reference-293.15));
    end
    assert(relativeRise < .02 && errorK < 1e-7 && heatError < 1e-12);
    clear tc ec
    cases{end+1} = struct('name',name,'nodes',data.nodes,'steps',size(data.reference_K,1), ...
        'unit_power_W',data.unit_power_W,'assembled_power_W',data.assembled_power_W, ...
        'max_temperature_difference_K',errorK,'max_relative_temperature_rise_error',relativeRise, ...
        'heat_vector_relative_error',heatError,'source_sha256',data.source_sha256);
end
after = radia.apiInfo(); assert(after.ih_handle_count == before.ih_handle_count);
report = struct('passed',true,'scope','real weak BEM -> production thermal assembly -> MEX; in-memory coarse axisymmetric geometry represented by 3D volume', ...
    'matlab_version',version,'cases',{cases},'native_build',jsondecode(fileread(which('radia_mex')+".build.json")), ...
    'handles_before',before.ih_handle_count,'handles_after',after.ih_handle_count);
fid=fopen(outputFile,'w');assert(fid>=0);fc=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
end
