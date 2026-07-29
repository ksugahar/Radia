function contract = compileFieldStudy(study, materialContract, options)
%COMPILEFIELDSTUDY Compile setup dictionaries into solver and fixed-Bus forms.

arguments
    study (1,1) struct
    materialContract (1,1) struct
    options.WindingContract (1,1) struct = struct()
    options.MaxBoundaries (1,1) double {mustBeInteger,mustBePositive} = 64
end
if ~isfield(study,"schema") || string(study.schema)~="radia.simulink.field-study-spec.v1"
    error("radia:simulink:FieldStudySpec","study must come from makeFieldStudySpec.");
end
if ~isfield(materialContract,"schema") || ...
        string(materialContract.schema)~="radia.simulink.material-dictionary.v1" || ...
        ~isfield(materialContract,"mesh") || ...
        strlength(string(materialContract.mesh.mesh_file))==0
    error("radia:simulink:FieldStudyMaterials", ...
        "materialContract must be compiled against a MeshFile.");
end
meshFile=string(materialContract.mesh.mesh_file);
assertPortableNames(string(materialContract.region_names(:)),".vol region");
boundaries=radia.simulink.inspectVolBoundaries(meshFile);
if double(boundaries.raw_boundary_id_count)>options.MaxBoundaries
    error("radia:simulink:FieldStudyBoundaryCapacity", ...
        "Raw boundary id count exceeds MaxBoundaries=%d.",options.MaxBoundaries);
end
dirichlet=study.dirichlet_values; robin=study.robin_boundaries;
dirichletNames=portableFields(dirichlet,"DirichletValues");
robinNames=portableFields(robin,"RobinBoundaries");
conductorNames=string(study.conductors(:));
forceNames=string(study.force_boundaries(:));
unknown=setdiff([dirichletNames;robinNames;conductorNames;forceNames], ...
    boundaries.boundary_names,"stable");
if ~isempty(unknown)
    error("radia:simulink:FieldStudyBoundary", ...
        "Study refers to unknown .vol boundaries: %s.",join(unknown,","));
end
if ~isempty(intersect(dirichletNames,robinNames))
    error("radia:simulink:FieldStudyBoundary", ...
        "A boundary cannot be both Dirichlet and Robin.");
end
if isempty(dirichletNames) && isempty(robinNames) && ...
        ~ismember(string(study.physics),["transient_heat","electrostatic_system"])
    error("radia:simulink:FieldStudyNullspace", ...
        "Specify a Dirichlet boundary or a positive thermal Robin boundary.");
end
physics=string(study.physics); scalar=physics~="harmonic_eddy";
if ~scalar && ~isempty(robinNames)
    error("radia:simulink:FieldStudyRobin", ...
        "Robin boundaries are available only for steady_heat.");
end
if ~ismember(physics,["steady_heat","transient_heat"]) && ~isempty(robinNames)
    error("radia:simulink:FieldStudyRobin", ...
        "Robin boundaries are available only for heat studies.");
end

runtime=emptyRuntime(options.MaxBoundaries);
runtime.physics_code=uint16(find(["electrostatic","current_flow", ...
    "steady_heat","harmonic_eddy","electrostatic_system","transient_heat"]==physics,1)-1);
runtime.formulation_code=uint16(study.formulation=="axisymmetric");
runtime.element_order=uint16(1+contains(string(study.element_family),"2"));
runtime.frequency_Hz=double(study.frequency_Hz);
runtime.model_depth_m=double(study.model_depth_m);
runtime.boundary_count=boundaries.raw_boundary_id_count;
runtime.time_sample_count=uint16(numel(study.time_s));
runtime.theta=double(study.theta);
if isnumeric(study.initial_temperature_K) && isscalar(study.initial_temperature_K)
    runtime.initial_temperature_K=double(study.initial_temperature_K);
else
    runtime.initial_temperature_K=NaN;
end
runtime.conductor_count=uint16(numel(conductorNames));
runtime.force_boundary_count=uint16(numel(forceNames));
for runtimeIndex=1:double(boundaries.raw_boundary_id_count)
        name=boundaries.raw_boundary_names(runtimeIndex);
        runtime.boundary_id(runtimeIndex)=boundaries.raw_boundary_ids(runtimeIndex);
        conductorIndex=find(conductorNames==name,1);
        if ~isempty(conductorIndex)
            runtime.boundary_kind(runtimeIndex)=uint16(3);
            runtime.boundary_value_real(runtimeIndex)=study.applied_voltages_V(conductorIndex);
        elseif isfield(dirichlet,char(name))
            value=dirichlet.(char(name));
            if ~isscalar(value) || ~isfinite(real(value)) || ~isfinite(imag(value))
                error("radia:simulink:FieldStudyBoundaryValue", ...
                    "Dirichlet value for '%s' must be a finite scalar.",name);
            end
            if physics~="current_flow" && imag(value)~=0
                error("radia:simulink:FieldStudyBoundaryValue", ...
                    "Only current_flow accepts complex boundary values.");
            end
            if physics=="harmonic_eddy" && real(value)~=0
                error("radia:simulink:FieldStudyHarmonicBoundary", ...
                    "harmonic_eddy accepts only explicit zero Dirichlet values.");
            end
            runtime.boundary_kind(runtimeIndex)=uint16(1);
            runtime.boundary_value_real(runtimeIndex)=real(value);
            runtime.boundary_value_imag(runtimeIndex)=imag(value);
        elseif isfield(robin,char(name))
            row=robin.(char(name));
            if ~isstruct(row) || ~isscalar(row) || ...
                    ~all(isfield(row,["transfer_w_per_m2_k","ambient_k"])) || ...
                    ~isfinite(row.transfer_w_per_m2_k) || row.transfer_w_per_m2_k<=0 || ...
                    ~isfinite(row.ambient_k)
                error("radia:simulink:FieldStudyRobin", ...
                    "Robin '%s' needs positive transfer_w_per_m2_k and finite ambient_k.",name);
            end
            runtime.boundary_kind(runtimeIndex)=uint16(2);
            runtime.robin_transfer_W_per_m2K(runtimeIndex)=row.transfer_w_per_m2_k;
            runtime.robin_ambient_K(runtimeIndex)=row.ambient_k;
        end
end

request=struct("physics",physics,"vol_text",fileread(meshFile), ...
    "source_name",string(meshFile),"element_family",string(study.element_family), ...
    "frequency_hz",double(study.frequency_Hz),"materials",struct());
if scalar
    request.formulation=string(study.formulation);
    if study.formulation=="planar",request.model_depth_m=double(study.model_depth_m);end
    request.dirichlet_values=pairStruct(dirichlet);
    request.robin_boundaries=robin;
    if ~isempty(fieldnames(study.terminal_pair)),request.terminal_pair=study.terminal_pair;end
    if physics=="transient_heat"
        request.time_s=study.time_s(:).';
        request.theta=study.theta;
        request.initial_temperature_k=study.initial_temperature_K;
        request.dirichlet_history_k=study.dirichlet_history_K;
        request.volumetric_source_history_w_per_m3= ...
            study.volumetric_source_history_W_per_m3;
    elseif physics=="electrostatic_system"
        request.conductors=cellstr(conductorNames);
        request.applied_voltages_v=study.applied_voltages_V(:).';
        request.force_boundaries=cellstr(forceNames);
    end
else
    if study.formulation=="axisymmetric"
        request.formulation="axisymmetric_henrotte";
    else
        request.formulation="planar";
    end
    request.dirichlet_boundaries=cellstr(dirichletNames);
end
request.materials=compileMaterials(physics,study,materialContract);
if ~scalar
    [request.branches,request.branch_current_a]=compileBranches( ...
        study,materialContract,options.WindingContract);
end
contract=struct("schema","radia.simulink.field-study.v1", ...
    "physics",physics,"mesh_sha256",materialContract.mesh.mesh_sha256, ...
    "boundary_names",boundaries.boundary_names,"request",request, ...
    "runtime",runtime,"runtime_policy",struct("fixed_width",true, ...
        "strings_per_step",false,"dictionary_lookup_per_step",false, ...
        "python_per_step",false,"batch_python_per_trigger",true));
end

function runtime=emptyRuntime(maxBoundaries)
runtime=struct("schema_version",uint16(1),"physics_code",uint16(0), ...
    "formulation_code",uint16(0),"element_order",uint16(1), ...
    "frequency_Hz",0,"model_depth_m",1,"boundary_count",uint16(0), ...
    "boundary_id",zeros(maxBoundaries,1,"uint32"), ...
    "boundary_kind",zeros(maxBoundaries,1,"uint16"), ...
    "boundary_value_real",zeros(maxBoundaries,1), ...
    "boundary_value_imag",zeros(maxBoundaries,1), ...
    "robin_transfer_W_per_m2K",zeros(maxBoundaries,1), ...
    "robin_ambient_K",zeros(maxBoundaries,1), ...
    "time_sample_count",uint16(0),"theta",1,"initial_temperature_K",293.15, ...
    "conductor_count",uint16(0),"force_boundary_count",uint16(0));
end

function names=portableFields(value,label)
names=string(fieldnames(value)); names=names(:);
assertPortableNames(names,label);
end

function assertPortableNames(names,label)
portable=cellfun(@(name)~isempty(regexp(name, ...
    "^[A-Za-z][A-Za-z0-9_]*$","once")),cellstr(names));
if any(~portable)
    error("radia:simulink:FieldStudyPortableName", ...
        "%s names must be portable MATLAB/JSON identifiers.",label);
end
end

function output=pairStruct(input)
output=struct(); names=string(fieldnames(input));
for k=1:numel(names)
    value=input.(char(names(k)));
    output.(char(names(k)))=[real(value),imag(value)];
end
end

function output=compileMaterials(physics,study,contract)
runtime=contract.runtime; output=struct(); sources=study.volumetric_sources;
sourceNames=portableFields(sources,"VolumetricSources");
unknown=setdiff(sourceNames,contract.region_names,"stable");
if ~isempty(unknown)
    error("radia:simulink:FieldStudySource", ...
        "VolumetricSources contains unknown regions: %s.",join(unknown,","));
end
mu0=4*pi*1e-7; eps0=8.8541878128e-12;
for k=1:numel(contract.region_names)
    region=contract.region_names(k); index=double(runtime.region_material_index(k));
    source=0;if isfield(sources,char(region)),source=sources.(char(region));end
    if ~isscalar(source) || ~isfinite(source)
        error("radia:simulink:FieldStudySource","Source for '%s' must be finite.",region);
    end
    if physics=="harmonic_eddy"
        if runtime.hysteresis_model_code(index)>0
            error("radia:simulink:FieldStudyHarmonicHysteresis", ...
                "harmonic_eddy does not treat hysteresis as a single-valued B-H curve.");
        end
        if runtime.bh_count(index)>0
            count=double(runtime.bh_count(index));
            B=runtime.bh_B_T(1:count,index);
            H=runtime.bh_H_A_per_m(1:count,index);
            if count<3 || B(1)~=0 || H(1)~=0
                error("radia:simulink:FieldStudyHarmonicBH", ...
                    "Nonlinear harmonic B-H data needs at least three rows starting at [0,0].");
            end
            rows=repmat(struct("b_t",0,"h_a_per_m",0),count,1);
            for j=1:count
                rows(j)=struct("b_t",B(j),"h_a_per_m",H(j));
            end
            row=struct("bh_curve",rows, ...
                "conductivity_s_per_m",runtime.conductivity_S_per_m(index));
        else
            row=struct("permeability_h_per_m",mu0*runtime.mu_r(index), ...
                "conductivity_s_per_m",runtime.conductivity_S_per_m(index));
        end
    elseif ismember(physics,["electrostatic","electrostatic_system"])
        row=struct("coefficient_si",eps0*runtime.relative_permittivity(index), ...
            "volumetric_source_si",source);
    elseif physics=="current_flow"
        row=struct("coefficient_si",runtime.conductivity_S_per_m(index), ...
            "volumetric_source_si",source);
        if study.frequency_Hz>0
            row.relative_permittivity=runtime.relative_permittivity(index);
        end
    else
        conductivity=runtime.thermal_conductivity_W_per_mK(index);
        if conductivity<=0
            error("radia:simulink:FieldStudyThermalMaterial", ...
                "Heat studies require positive thermal conductivity in region '%s'.",region);
        end
        row=struct("coefficient_si",conductivity,"volumetric_source_si",source);
        if physics=="transient_heat"
            density=runtime.density_kg_per_m3(index);
            specificHeat=runtime.specific_heat_J_per_kgK(index);
            if density<=0 || specificHeat<=0
                error("radia:simulink:FieldStudyThermalCapacity", ...
                    "transient_heat requires positive density and specific heat in region '%s'.",region);
            end
            row.volumetric_heat_capacity_j_per_m3_k=density*specificHeat;
        end
    end
    output.(char(region))=row;
end
end

function [branches,currents]=compileBranches(study,materials,windings)
if ~isfield(windings,"schema") || ...
        string(windings.schema)~="radia.simulink.winding-dictionary.v1"
    error("radia:simulink:FieldStudyWinding", ...
        "harmonic_eddy requires a compiled WindingContract.");
end
count=double(windings.runtime.winding_count); input=study.branch_current_A(:);
if numel(input)~=count || ~any(abs(input)>0)
    error("radia:simulink:FieldStudyCurrent", ...
        "BranchCurrent_A must contain one value per winding and one must be nonzero.");
end
branchCount=sum(double(windings.runtime.region_count(1:count)));
branches=repmat(struct("name","","material","","turns",0),branchCount,1);
currents=zeros(branchCount,2);
branchIndex=0;
for k=1:count
    n=double(windings.runtime.region_count(k));
    for j=1:n
        id=windings.runtime.region_id(j,k);
        regionIndex=find(materials.mesh.region_ids==id,1);
        if isempty(regionIndex)
            error("radia:simulink:FieldStudyWinding","Winding region id is absent from mesh.");
        end
        region=materials.mesh.region_names(regionIndex);
        turns=windings.runtime.effective_turns(k)* ...
            double(windings.runtime.region_polarity(j,k));
        branchIndex=branchIndex+1;
        branches(branchIndex)=struct("name",char(windings.winding_names(k)+"_"+region), ...
            "material",char(region),"turns",turns);
        currents(branchIndex,:)=[real(input(k)),imag(input(k))];
    end
end
end
