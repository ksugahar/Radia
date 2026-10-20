function writeHDivH5(path, fields, metadata)
%WRITEHDIVH5 Typed dense HDiv output for independent CPU validation.
% This minimal writer rejects sparse input; CPU fixtures use canonical CSR.
arguments
    path (1,1) string
    fields (1,1) struct
    metadata (1,1) struct
end
[~,~,ext] = fileparts(path);
assert(any(strcmpi(ext,{'.h5','.hdf5'})), 'radia:HDF5', 'HDF5 only');
assert(~isfile(path), 'radia:HDF5', 'Refusing existing output');
assert(isfield(metadata,'units') && strcmp(metadata.units,'SI') && ...
    isfield(metadata,'labels') && iscellstr(metadata.labels), ...
    'radia:HDF5','Explicit SI units and text labels required');
[~,~,endian] = computer;
assert(endian=='L','radia:HDF5','Little-endian host required');
names = sort(fieldnames(fields));
assert(~isempty(names),'radia:HDF5','Empty fields');
% Validate before creating any file.
for k=1:numel(names)
    x=fields.(names{k});
    assert(~isempty(regexp(names{k},'^[A-Za-z][A-Za-z0-9_]*$','once')) && ...
        (isnumeric(x)||islogical(x)) && ~issparse(x) && ~isempty(x) && ismatrix(x), ...
        'radia:HDF5','Only nonempty dense rank-one/two numerical fields');
    assert(any(strcmp(class(x),{'double','single','int64','uint64','int32','uint32','uint8','logical'})), ...
        'radia:HDF5','Unsupported dtype');
end
for k=1:numel(names)
    x=fields.(names{k}); base='/fields/'+string(names{k});
    if isvector(x), shape=numel(x); x=x(:); else, shape=size(x); end
    cls=class(x); dtype=cls;
    if strcmp(cls,'double'), dtype='float64'; elseif strcmp(cls,'single'), dtype='float32';
    elseif strcmp(cls,'logical'), dtype='bool'; end
    parts={'real'}; representation='real';
    if ~isreal(x), parts={'real','imag'}; representation='complex-split'; end
    for j=1:numel(parts)
        if j==1, v=real(x); else, v=imag(x); end
        if islogical(v), v=uint8(v); end
        % MATLAB reverses HDF5 axes: nonconjugate transpose gives file [row,col].
        if numel(shape)==2, v=v.'; dims=fliplr(shape); else, dims=shape; end
        ds=base+'/'+parts{j};
        h5create(path,ds,dims,'Datatype',class(v)); h5write(path,ds,v);
        bytes=typecast(v(:),'uint8');
        md=java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
        digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
        h5writeatt(path,ds,'sha256',digest,'TextEncoding','UTF-8');
    end
    h5writeatt(path,base,'kind','dense');
    if numel(shape)==1, shapeText=sprintf('[%d]',shape); else, shapeText=jsonencode(shape); end
    h5writeatt(path,base,'shape_json',shapeText);
    h5writeatt(path,base,'dtype',dtype);
    h5writeatt(path,base,'representation',representation);
end
h5writeatt(path,'/','schema','radia.hdiv-hdf5.v1');
h5writeatt(path,'/','units','SI'); h5writeatt(path,'/','index_base',uint8(0));
h5writeatt(path,'/','layout','file_axes_row_column_C_order');
h5writeatt(path,'/','text_encoding','UTF-8');
h5writeatt(path,'/','metadata_json',jsonencode(metadata),'TextEncoding','UTF-8');
h5writeatt(path,'/','field_names_json',jsonencode(names));
end
