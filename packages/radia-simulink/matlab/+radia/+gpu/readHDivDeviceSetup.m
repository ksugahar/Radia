function setup=readHDivDeviceSetup(value)
% Checked radia.hdiv-device-setup.v1 payload for radia_hdiv_device.
% Integer tables are restored to int32 (exact integers only) and text to
% char; the MEX recomputes every dimension from the flat counts.
% Repeated solves pass the same payload, so only a changed payload is checked
% again; the last checked payload stays referenced until the next one replaces it.
persistent checkedInput checkedReal checkedSetup
if isstruct(value) && ~isempty(checkedInput) && isequal(value,checkedInput) && isequal(realFlags(value),checkedReal)
    setup=checkedSetup;return
end
assert(isstruct(value) && isscalar(value),'radia:gpu:DeviceSetup','Device setup must be a scalar struct');
required={'schema';'backend';'family';'order';'dofs';'mesh_sha256';'label_contract_sha256'; ...
    'charges';'evaluation';'rhs_operator';'mean_operator'};
assert(isequal(sort(fieldnames(value)),sort(required)),'radia:gpu:DeviceSetup','Invalid device setup fields');
setup=struct();
for name={'schema','backend','family','mesh_sha256','label_contract_sha256'}
    setup.(name{1})=char(value.(name{1}));
end
assert(strcmp(setup.schema,'radia.hdiv-device-setup.v1') && strcmp(setup.backend,'gpu') && ...
    any(strcmp(setup.family,{'affine-tet','curved-tet','hex'})),'radia:gpu:DeviceSetup','Unsupported device setup schema, backend or family');
curved=strcmp(setup.family,'curved-tet');hex=strcmp(setup.family,'hex');
assert(~isempty(regexp(setup.mesh_sha256,'^[0-9a-f]{64}$','once')) && ...
    ~isempty(regexp(setup.label_contract_sha256,'^[0-9a-f]{64}$','once')), ...
    'radia:gpu:DeviceSetup','Missing checked mesh/label identity');
setup.order=plain(value.order);setup.dofs=plain(value.dofs);nd=setup.dofs;
assert(isscalar(setup.order) && ismember(setup.order,[0 1 2]) && isscalar(nd) && nd>=1 && fix(nd)==nd, ...
    'radia:gpu:DeviceSetup','Invalid device setup order or dimension');
% Curved P2 TET and HEX geometry carry orders 1 and 2 only, matching the CPU routes.
assert(~(curved || hex) || ismember(setup.order,[1 2]),'radia:gpu:DeviceSetup','Curved or HEX device setup requires order 1 or 2');
setup.rhs_operator=double(value.rhs_operator);setup.mean_operator=double(value.mean_operator);
assert(isequal(size(setup.rhs_operator),[nd 3]) && isequal(size(setup.mean_operator),[3 nd]) && ...
    isreal(setup.rhs_operator) && isreal(setup.mean_operator) && ...
    all(isfinite(setup.rhs_operator(:))) && all(isfinite(setup.mean_operator(:))), ...
    'radia:gpu:DeviceSetup','Invalid NGSolve excitation/mean numbering');
setup.rhs_operator=setup.rhs_operator+0;setup.mean_operator=setup.mean_operator+0;
c=value.charges;
setup.charges=tables(c,{'cell_vertices','face_vertices','tet_rule','tri_rule','tet_low_rule','tri_low_rule', ...
    'cell_tables','face_tables'},{'charge_host','charge_kind','charge_exponents','cell_class','cell_dofs', ...
    'cell_exponents','face_class','face_dofs','face_exponents'});
curvedFields={'cell_numbers','face_numbers','curve_gl','curve_gw','cell_mass_rule','face_mass_rule'};
assert(all(isfield(c,curvedFields)==curved),'radia:gpu:DeviceSetup', ...
    'Curved geometry tables are required for curved-tet and rejected for affine-tet');
hexFields={'hex_rules','hex_bands'};
assert(all(isfield(c,hexFields)==hex),'radia:gpu:DeviceSetup', ...
    'HEX Gram rules and bands are required for hex and rejected otherwise');
extra=struct();
if curved
    extra=tables(c,{'curve_gl','curve_gw','cell_mass_rule','face_mass_rule'},{'cell_numbers','face_numbers'});
elseif hex
    extra=tables(c,hexFields,{});
end
for name=fieldnames(extra)'
    setup.charges.(name{1})=extra.(name{1});
end
setup.charges.family=setup.family;
setup.charges.far_factor=plain(c.far_factor);setup.charges.compression=char(c.compression);
setup.charges.dofs=plain(c.dofs);
assert(setup.charges.dofs==nd,'radia:gpu:DeviceSetup','Charge numbering differs from the FE dimension');
if strcmp(setup.charges.compression,'dense')
    assert(isempty(c.aca),'radia:gpu:DeviceSetup','Dense compression carries no ACA parameters');
    setup.charges.aca=zeros(0,0);
else
    a=c.aca;assert(isstruct(a) && isscalar(a),'radia:gpu:DeviceSetup','ACA compression requires its parameters');
    setup.charges.aca=struct('eps',plain(a.eps),'eta',plain(a.eta),'leaf_size',plain(a.leaf_size), ...
        'max_rank',plain(a.max_rank),'rank_cap',char(a.rank_cap));
end
setup.evaluation=tables(value.evaluation,{'rule','reference_shapes','vertices'}, ...
    {'element_class','element_dofs','element_points'});
checkedInput=value;checkedReal=realFlags(value);checkedSetup=setup;
end

function flags=realFlags(value)
% isequal ignores complexity. It also ignores class, sparsity and the sign of zero, which the read
% below normalizes, so equal payloads with equal flags read identically. Builtins only: this runs
% on every solve, and MATLAB code run this often shows up in the memory slope.
flags=[cellfun('isreal',struct2cell(value));cellfun('isreal',struct2cell(value.charges)); ...
    cellfun('isreal',struct2cell(value.evaluation))];
if isstruct(value.charges.aca)
    flags=[flags;cellfun('isreal',struct2cell(value.charges.aca))];
end
end

function v=plain(v)
v=double(v)+0;
end

function out=tables(s,reals,integers)
assert(isstruct(s) && isscalar(s) && all(isfield(s,[reals integers])), ...
    'radia:gpu:DeviceSetup','Missing device table fields');
out=struct();
for name=reals
    v=double(s.(name{1}));
    assert(isreal(v) && all(isfinite(v(:))),'radia:gpu:DeviceSetup','Non-finite device table: %s',name{1});
    out.(name{1})=v(:)+0;
end
for name=integers
    v=s.(name{1});
    if ~isa(v,'int32')
        v=full(double(v));
        assert(isreal(v) && all(v(:)==fix(v(:))) && all(abs(v(:))<=2^31-1), ...
            'radia:gpu:DeviceSetup','Integer device table is not exact int32: %s',name{1});
    end
    out.(name{1})=int32(v(:));
end
end
