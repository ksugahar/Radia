function paths=recordReactionGPU(baseConfig,assembly,wheel,sourceCommit,outdir,options)
%RECORDREACTIONGPU Acquire actual accepted port work for two dedicated cases.
% Call only after compute admission and formula review. No physical pass claim.
% baseConfig/assembly must come from the installed public temperature-BH sample.
arguments
    baseConfig (1,1) struct
    assembly (1,1) struct
    wheel (1,:) char
    sourceCommit (1,:) char
    outdir (1,:) char
    options.CurrentA (1,1) double {mustBeFinite}=10000
    options.AmbientK (1,1) double {mustBeFinite,mustBePositive}=300
    options.Updates (1,1) double {mustBeInteger,mustBePositive}=5
end
assert(options.CurrentA~=0 && options.Updates<=1000,'Invalid bounded acquisition');
assert(~isempty(regexp(sourceCommit,'^[0-9a-f]{40}$','once')),'Full source SHA required');
assert(~isfolder(outdir) && ~isfile(outdir),'Output directory must be new');
assert(isfile(wheel),'Wheel required');
assert(strcmp(baseConfig.material_mode,'element-esim-operators'),'Element material required');
assert(isfield(baseConfig,'bh_knots') && isfield(baseConfig,'bh_values_T'),'Public BH sample required');
indices=assembly.unknown_indices;coil=indices.coil;port=indices.port;
n=size(baseConfig.eddy,1);
assert(isscalar(coil) && isscalar(port) && coil~=port && ...
    all([coil,port]>=1) && all([coil,port]<=n) && ...
    all(mod([coil,port],1)==0),'Single-coil checked indices required');
loaded=which('radia_cuda_ih');assert(~isempty(loaded),'Installed IH MEX required');
loaded=char(java.io.File(loaded).getCanonicalPath());
assert(endsWith(loaded,'.mexw64'),'Native Windows MEX required');
[~,resident]=inmem('-completenames');
for j=1:numel(resident)
    [~,name]=fileparts(resident{j});
    if strcmp(name,'radia_cuda_ih')
        actual=char(java.io.File(resident{j}).getCanonicalPath());
        assert(strcmpi(actual,loaded),'Another IH MEX is already loaded; reconnect the owning session');
    end
end
planPath=fullfile(fileparts(mfilename('fullpath')),'reaction_gate_plan.json');
plan=jsondecode(fileread(planPath));
identity=struct('source_commit',sourceCommit,'source_dirty',false,...
    'wheel_sha256',fileHash(wheel),'loaded_mex_sha256',fileHash(loaded),...
    'loaded_mex_path',loaded,'plan_sha256',fileHash(planPath));
% The caller must verify the clean source/build receipt before calling.
% source_dirty=false is an acquisition declaration, not a Git inspection here.
mkdir(outdir);paths=cell(2,1);
for caseIndex=1:2
    cfg=baseConfig;
    % Dedicated zero-coil-impedance experiment, not a post-hoc correction to V.
    cfg.eddy(coil,coil)=0;
    kind='temperature-bh-zero-resistance';
    if caseIndex==2
        kind='nonmagnetic-reactive';
        cfg=rmfield(cfg,{'bh_knots','bh_values_T'});
        cfg.initial_permeability_H_per_m=...
            (4*pi*1e-7)*ones(size(cfg.initial_permeability_H_per_m));
        cfg.saturation_magnetization_A_per_m=...
            zeros(size(cfg.saturation_magnetization_A_per_m));
        assert(all(isfinite(cfg.resistivity_ohm_m),'all') && ...
            all(cfg.resistivity_ohm_m>0,'all'),'Positive finite resistivity required');
    end
    assert(cfg.eddy(coil,coil)==0,'Zero supplied coil impedance required');
    caseDir=fullfile(outdir,kind);mkdir(caseDir);
    save(fullfile(caseDir,'input.mat'),'cfg','assembly','identity');
    h=radia_cuda_ih('create',cfg);closer=onCleanup(@()radia_cuda_ih('destroy',h));
    rows=cell(options.Updates,1);
    for k=1:options.Updates
        radia_cuda_ih('update',h,options.CurrentA,options.AmbientK);
        [temperature,heat,info]=radia_cuda_ih('outputs',h);
        em=radia_cuda_ih('electromagnetic_state',h);
        z=radia_cuda_ih('element_impedance',h);
        assert(info.steps==k && all(isfinite(em)) && all(isfinite(heat)),...
            'Incomplete accepted state');
        assert(abs(em(coil)-options.CurrentA)/abs(options.CurrentA)<=...
            plan.linear_true_relative_residual_limit,'Prescribed current constraint failed');
        row=struct('update_index',k,'voltage_update_index',k,...
            'current_update_index',k,'heat_update_index',k,'state','accepted',...
            'current_A',[options.CurrentA,0],...
            'port_voltage_V',[real(em(port)),imag(em(port))],...
            'face_heat_W',reshape(heat,1,[]),...
            'linear_true_relative_residual',info.eddy_true_residual,...
            'constitutive_true_relative_residual',info.constitutive_true_residual,...
            'actual_coil_current_A',[real(em(coil)),imag(em(coil))],...
            'temperature_K',reshape(temperature,1,[]),...
            'element_impedance_ohm',[real(z(:)),imag(z(:))],...
            'eddy_iterations',info.eddy_iterations,...
            'constitutive_iterations',info.constitutive_iterations);
        rows{k}=row;
        % Preserve each completed update before attempting the next one.
        writeJSON(fullfile(caseDir,sprintf('accepted-%04d.json',k)),row);
    end
    record=struct('phasor_convention',plan.phasor_convention,...
        'case_kind',kind,'supplied_coil_impedance_ohm',[0,0],...
        'artifact_identity',identity,'accepted_updates',{rows},...
        'acquisition_driver_sha256',fileHash([mfilename('fullpath'),'.m']),...
        'b1_closed',false,'release_accepted',false);
    if caseIndex==2
        record.material_contract='constant-mu0-positive-resistivity-no-remanence-no-bh';
    end
    paths{caseIndex}=fullfile(caseDir,'records.json');writeJSON(paths{caseIndex},record);
    clear closer;
end
end

function writeJSON(path,value)
f=fopen(path,'w');assert(f>=0,'Cannot write acquisition evidence');
closer=onCleanup(@()fclose(f));fprintf(f,'%s\n',jsonencode(value));
end

function value=fileHash(file)
f=fopen(file,'rb');assert(f>=0,'Cannot read evidence input');
closer=onCleanup(@()fclose(f));
md=java.security.MessageDigest.getInstance('SHA-256');md.update(fread(f,Inf,'*uint8'));
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
