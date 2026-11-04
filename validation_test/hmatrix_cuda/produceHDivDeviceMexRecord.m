function produceHDivDeviceMexRecord(inputPath,outPath,commonReceiptPath)
% Plain-MEX producer of the device HDiv frontend parity (PYBIND_MEX_BITWISE.md).
% Reads only the canonical input HDF5 written by hdiv_device_parity.py
% make-input, calls radia_hdiv_device on the shared nvcc kernel object and
% writes the MEX record. Frontend parity only; no tolerance and no CPU oracle.
arguments
    inputPath (1,:) char {mustBeFile}
    outPath (1,:) char {mustBeNonempty}
    commonReceiptPath (1,:) char {mustBeFile}
end
assert(~isfile(outPath),'radia:parity:Exists','Refusing to overwrite %s',outPath);
common=jsondecode(fileread(commonReceiptPath));
assert(isequal(sort(fieldnames(common)),sort({'source_sha';'kernel_object_sha256';'cuda_compiler'; ...
    'python_binary_sha256';'mex_binary_sha256'})),'radia:parity:Receipt','Unexpected common receipt fields');
assert(strcmp(textAttribute(inputPath,'/','schema'),'radia.hdiv-device-parity-input.v1') && ...
    any(strcmp(textAttribute(inputPath,'/','case'),{'device-flat-tet','device-mixed'})),'radia:parity:Input','Not a device parity input');
input=readHDivParityInput(inputPath);
nd=input.dofs;charges=input.charges;evaluation=input.evaluation;law=input.law;applyInput=input.apply;
rhs=input.rhs;initial=input.initial;newton=input.newton;
op=radia_hdiv_device('create_operator',charges);opOwner=onCleanup(@()radia_hdiv_device('destroy',op));
mass=radia_hdiv_device('create_mass',evaluation,nd);massOwner=onCleanup(@()radia_hdiv_device('destroy',mass));
material=radia_hdiv_device('create_material',evaluation,nd,law);materialOwner=onCleanup(@()radia_hdiv_device('destroy',material));
counters=radia_hdiv_device('assembly',op).hmatrix;
assert(counters.lowrank_leaves>=1,'radia:parity:ACA','Parity input did not exercise device ACA low-rank leaves');
applied=radia_hdiv_device('apply',op,applyInput);
lin=radia_hdiv_device('solve_linear',op,mass,rhs,input.inv_chi,input.tolerance,input.maximum_iterations,input.preconditioner);
new=radia_hdiv_device('solve_nonlinear',op,material,rhs,initial,newton);
clear materialOwner massOwner opOwner
assert(radia_hdiv_device('live_count')==0,'radia:parity:Leak','Device handles remain after the parity run');
[~,gpu]=system('nvidia-smi --query-gpu=uuid,compute_cap,name --format=csv,noheader');
rows=strtrim(splitlines(strtrim(gpu)));assert(isscalar(rows),'radia:parity:Device','Parity requires exactly one visible GPU');
parts=strtrim(split(rows{1},','));
receipt=common;receipt.frontend='mex';receipt.case=input.name;receipt.input_sha256=sha256(inputPath);
receipt.output_names={'apply';'linear_solution';'nonlinear_solution';'nonlinear_residual_history'};
receipt.artifact_sha256=sha256(which('radia_hdiv_device'));receipt.hmatrix=counters;
receipt.linear_relative_residual=lin.relative_residual;receipt.nonlinear_relative_residual=new.relative_residual;
receipt.device_uuid=parts{1};receipt.architecture=['sm_' strrep(parts{2},'.','')];receipt.device_name=parts{3};
receipt.matlab_release=version('-release');
outputs={applied,lin.x,new.x,[new.history.relative_residual]'};
for k=1:numel(outputs)
    name=['/outputs/' receipt.output_names{k}];value=double(outputs{k}(:));
    h5create(outPath,name,numel(value),'Datatype','double');h5write(outPath,name,value);
end
h5writeatt(outPath,'/','receipt_json',jsonencode(receipt));
end

function value=textAttribute(path,location,name)
value=char(string(h5readatt(path,location,name)));
end

function digest=sha256(path)
fid=fopen(path,'r');assert(fid>0,'radia:parity:Read','Cannot read %s',path);
bytes=fread(fid,inf,'*uint8');fclose(fid);
engine=java.security.MessageDigest.getInstance('SHA-256');engine.update(typecast(bytes,'int8'));
digest=lower(reshape(dec2hex(typecast(engine.digest(),'uint8'),2)',1,[]));
end
