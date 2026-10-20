function [fields, metadata] = readHDivH5(path, maxBytes)
%READHDIVH5 Strict typed independent CPU fixture reader; no MAT fallback.
arguments
    path (1,1) string
    maxBytes (1,1) double {mustBePositive,mustBeFinite} = 512*2^20
end
[~,~,ext]=fileparts(path);
assert(any(strcmpi(ext,{'.h5','.hdf5'})),'radia:HDF5','HDF5 only');
assert(strcmp(h5readatt(path,'/','schema'),'radia.hdiv-hdf5.v1') && ...
    strcmp(h5readatt(path,'/','units'),'SI') && h5readatt(path,'/','index_base')==0 && ...
    strcmp(h5readatt(path,'/','layout'),'file_axes_row_column_C_order') && ...
    strcmp(h5readatt(path,'/','text_encoding'),'UTF-8'),'radia:HDF5','Contract mismatch');
metadata=jsondecode(h5readatt(path,'/','metadata_json'));
assert(strcmp(metadata.units,'SI') && (iscellstr(metadata.labels)||isempty(metadata.labels)), ...
    'radia:HDF5','Invalid physical metadata');
if isempty(metadata.labels), metadata.labels={}; else, metadata.labels=reshape(metadata.labels,1,[]); end
names=jsondecode(h5readatt(path,'/','field_names_json'));
assert(iscellstr(names) && numel(unique(names))==numel(names),'radia:HDF5','Invalid field inventory');
info=h5info(path,'/fields'); actual=cellfun(@(x) extractAfter(x,'/fields/'),{info.Groups.Name},'UniformOutput',false);
assert(isequal(sort(names(:)),sort(actual(:))),'radia:HDF5','Field inventory mismatch');
fields=struct; remaining=maxBytes;
for k=1:numel(names)
    name=names{k}; assert(~isempty(regexp(name,'^[A-Za-z][A-Za-z0-9_]*$','once')),'radia:HDF5','Invalid field');
    base='/fields/'+string(name); kind=h5readatt(path,base,'kind');
    if strcmp(kind,'dense')
        fields.(name)=values(base);
    elseif strcmp(kind,'csr')
        matrixShape=jsondecode(h5readatt(path,base,'matrix_shape_json'));
        assert(numel(matrixShape)==2 && all(matrixShape>=1 & matrixShape==fix(matrixShape) & matrixShape<flintmax),'radia:HDF5','Invalid CSR matrixShape');
        v=values(base+'/values'); assert(isvector(v),'radia:HDF5','CSR values must be vector');
        ptr=component(base+'/indptr',[matrixShape(1)+1], 'uint64');
        idx=component(base+'/indices',[numel(v)], 'uint64');
        assert(ptr(1)==0 && ptr(end)==numel(v) && all(ptr(2:end)>=ptr(1:end-1)) && ...
            all(idx<uint64(matrixShape(2))),'radia:HDF5','Invalid CSR indices');
        for row=1:matrixShape(1)
            a=double(ptr(row))+1; b=double(ptr(row+1));
            assert(all(idx(a+1:b)>idx(a:b-1)),'radia:HDF5','Noncanonical CSR');
        end
        rows=repelem((1:matrixShape(1))',double(diff(ptr)));
        fields.(name)=sparse(rows,double(idx)+1,double(v),matrixShape(1),matrixShape(2));
    else
        error('radia:HDF5','Unknown field kind');
    end
end
    function v=values(base)
        valueShape=jsondecode(h5readatt(path,base,'shape_json'));
        assert(numel(valueShape)<=2 && ~isempty(valueShape) && all(valueShape>=1 & valueShape==fix(valueShape)),'radia:HDF5','Invalid valueShape');
        dtype=h5readatt(path,base,'dtype');
        switch dtype
            case 'float64', cls='double';
            case 'float32', cls='single';
            case 'bool', cls='uint8';
            otherwise, cls=dtype;
        end
        v=component(base+'/real',valueShape,cls);
        representation=h5readatt(path,base,'representation');
        if strcmp(representation,'complex-split')
            assert(any(strcmp(cls,{'double','single'})),'radia:HDF5','Invalid complex dtype');
            v=complex(v,component(base+'/imag',valueShape,cls));
        else
            assert(strcmp(representation,'real'),'radia:HDF5','Invalid representation');
            if strcmp(dtype,'bool')
                assert(all(v(:)==0|v(:)==1),'radia:HDF5','Invalid logical data'); v=logical(v);
            end
        end
    end
    function v=component(ds,shape,cls)
        types={'double','single','int64','uint64','int32','uint32','uint8'};
        htypes={'H5T_IEEE_F64LE','H5T_IEEE_F32LE','H5T_STD_I64LE','H5T_STD_U64LE','H5T_STD_I32LE','H5T_STD_U32LE','H5T_STD_U8LE'};
        widths=[8 4 8 8 4 4 1]; at=find(strcmp(cls,types));
        assert(~isempty(at),'radia:HDF5','Unsupported dtype');
        info=h5info(path,ds); dims=double(info.Dataspace.Size);
        assert(isequal(dims(:),flip(shape(:))) && strcmp(info.Datatype.Type,htypes{at}), ...
            'radia:HDF5','Shape/dtype mismatch');
        bytes=prod(double(shape))*widths(at);
        assert(isfinite(bytes) && bytes<=remaining,'radia:HDF5','Read budget exceeded'); remaining=remaining-bytes;
        raw=h5read(path,ds); assert(strcmp(class(raw),cls),'radia:HDF5','Loaded dtype mismatch');
        md=java.security.MessageDigest.getInstance('SHA-256'); md.update(typecast(raw(:),'uint8'));
        digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
        assert(strcmp(digest,h5readatt(path,ds,'sha256')),'radia:HDF5','SHA256 mismatch');
        if numel(shape)==2, v=raw.'; else, v=raw(:); end
    end
end
