function setup=loadHDivDeviceSetup(path)
% Read a radia.hdiv-device-setup.v1 HDF5 artifact into the checked struct that
% radia_hdiv_device consumes. The artifact comes from prepareHDivDeviceSetup or
% from `python -m radia_simulink.hdiv request.json setup.h5` run beforehand, so
% a MATLAB session consuming it needs no Python.
arguments
    path (1,:) char {mustBeNonempty}
end
setup=radia.gpu.readHDivDeviceSetup(readSetupFile(path));
end

function value=readSetupFile(path)
% Row-major HDF5 tables read as reversed MATLAB arrays; the flat tables are vectors.
attribute=@(location,name)h5readatt(path,location,name);
text=@(location,name)char(string(attribute(location,name)));
assert(strcmp(text('/','units'),'SI') && double(attribute('/','index_base'))==0,'radia:gpu:DeviceSetup', ...
    'Device setup file requires SI units and 0-based indices');
value=struct();
for name={'schema','backend','family','mesh_sha256','label_contract_sha256'}
    value.(name{1})=text('/',name{1});
end
value.order=double(attribute('/','order'));value.dofs=double(attribute('/','dofs'));
value.rhs_operator=h5read(path,'/rhs_operator')';value.mean_operator=h5read(path,'/mean_operator')';
value.charges=readTables(path,'/charges');value.evaluation=readTables(path,'/evaluation');
value.charges.far_factor=double(attribute('/charges','far_factor'));
value.charges.compression=text('/charges','compression');value.charges.dofs=double(attribute('/charges','dofs'));
if strcmp(value.charges.compression,'dense')
    value.charges.aca=zeros(0,0);
else
    value.charges.aca=struct('eps',double(attribute('/charges/aca','eps')),'eta',double(attribute('/charges/aca','eta')), ...
        'leaf_size',double(attribute('/charges/aca','leaf_size')),'max_rank',double(attribute('/charges/aca','max_rank')), ...
        'rank_cap',text('/charges/aca','rank_cap'));
end
end

function tables=readTables(path,group)
tables=struct();info=h5info(path,group);
for k=1:numel(info.Datasets)
    name=info.Datasets(k).Name;tables.(name)=h5read(path,[group '/' name]);
end
end
