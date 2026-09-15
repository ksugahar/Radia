function report = validate_ih_mex_radial(outputFile)
% Axisymmetric radial FEM operators versus native Thermal MEX and exact steady heat.
% Independent MATLAB assembly; no CAD/export or EM-source accuracy claim.
arguments
    outputFile (1,1) string
end
before = radia.apiInfo();
cases = {};
for bore = [0,0.0075]
    r = linspace(bore,0.025,33).'; n = numel(r);
    k = 10; capacity = 1000; h = 100; ambient = 293.15; q = 1e5;
    M = zeros(n); K = zeros(n); C = zeros(n);
    for cell = 1:n-1
        ids = cell:cell+1; dr = r(cell+1)-r(cell);
        for xi = [-sqrt(3/5),0,sqrt(3/5)]
            w = 5/9; if xi == 0, w = 8/9; end
            N = [(1-xi)/2;(1+xi)/2]; radius = N.'*r(ids);
            dv = 2*pi*radius*dr*w/2;
            M(ids,ids) = M(ids,ids)+capacity*dv*(N*N.');
            D = [-1;1]/dr;
            K(ids,ids) = K(ids,ids)+k*dv*(D*D.');
        end
    end
    C(end,end) = 2*pi*r(end); % Only outer wall exchanges heat. Inner/axis insulated.
    volumes = sum(M,2)/capacity;
    cfg = radia.simulink.makeIHNativeSmokeConfig();
    cfg.n_temperature = n; cfg.n_heat = n; cfg.rotation_mode = 'none';
    cfg.initial_temperature_K = ambient*ones(n,1);
    cfg.temperature_cell_weights = sum(M,2);
    cfg.heat_cell_weights = volumes;
    cfg.heat_to_temperature_projection = reshape(diag(volumes).',[],1);
    cfg.heat_projection = q*ones(n,1);
    cfg.sample_time_s = 0.01; cfg.convection_W_per_m2K = h;
    cfg.thermal_tolerance = 1e-13; cfg.thermal_max_iterations = 500;
    for item = {'mass','stiffness','convection'}
        key = item{1};
        if strcmp(key,'mass'), A=M; elseif strcmp(key,'stiffness'), A=K; else, A=C; end
        % Shared full CSR pattern, including zeros, zero-based indices.
        cfg.([key '_row_ptr']) = (0:n:n*n).';
        cfg.([key '_col']) = repmat((0:n-1).',n,1);
        cfg.([key '_value']) = reshape(A.',[],1);
    end
    eddy = radia.internal.callMex('ih.eddy.create',cfg);
    eddyCleanup = onCleanup(@() radia.internal.callMex('ih.eddy.destroy',eddy));
    thermal = radia.internal.callMex('ih.thermal.create',cfg);
    cleanup = onCleanup(@() radia.internal.callMex('ih.thermal.destroy',thermal));
    reference = cfg.initial_temperature_K; maxRiseError = 0; maxBalance = 0;
    system = M + cfg.sample_time_s*(K+h*C);
    for step = 1:1000
        heat = radia.internal.callMex('ih.eddy.output',eddy,1,0,reference);
        assert(norm(heat-q,inf) < 1e-9);
        previous = radia.internal.callMex('ih.thermal.output',thermal);
        source = volumes.*heat;
        reference = system\(M*reference+cfg.sample_time_s*(source+h*C*ones(n,1)*ambient));
        radia.internal.callMex('ih.thermal.update',thermal,heat,ambient,0);
        actual = radia.internal.callMex('ih.thermal.output',thermal);
        maxRiseError = max(maxRiseError,norm(actual-reference,inf));
        balance = sum(M*(actual-previous))/cfg.sample_time_s + h*sum(C*(actual-ambient))-sum(source);
        maxBalance = max(maxBalance,abs(balance)/sum(source));
    end
    R = r(end); a = r(1);
    exact = ambient + q*(R*R-a*a)/(2*R*h) + q*(R*R-r.^2)/(4*k);
    if a > 0, exact = exact + q*a*a/(2*k)*log(r/R); end
    steadyRelative = norm(actual-exact)/norm(exact-ambient);
    assert(maxRiseError < 1e-7 && maxBalance < 1e-7 && steadyRelative < 0.02);
    clear cleanup eddyCleanup
    cases{end+1} = struct('bore_m',bore,'nodes',n,'steps',1000, ...
        'mex_vs_direct_max_K',maxRiseError,'power_balance_relative',maxBalance, ...
        'steady_temperature_rise_relative_L2',steadyRelative);
end
after = radia.apiInfo(); assert(after.ih_handle_count == before.ih_handle_count);
report = struct('scope','IH radial P1 thermal MEX; prescribed uniform source; not EM or production VOL validation', ...
    'matlab_version',version,'mex_path',which('radia_mex'),'cases',{cases}, ...
    'passed',true,'application','IH','steady_tolerance',0.02, ...
    'handle_count_before',before.ih_handle_count,'handle_count_after',after.ih_handle_count, ...
    'native_build',jsondecode(fileread(which('radia_mex') + ".build.json")));
fid = fopen(outputFile,'w'); assert(fid >= 0); fileCleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
end
