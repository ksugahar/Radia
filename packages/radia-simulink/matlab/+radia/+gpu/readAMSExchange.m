function data=readAMSExchange(path)
%READAMSEXCHANGE Read radia.ams-hdf5.v1 numeric data without MAT size limits.
manifest=jsondecode(native2unicode(uint8(h5read(path,'/manifest'))','UTF-8'));
assert(strcmp(manifest.schema,'radia.ams-hdf5.v1') && manifest.index_base==0 && ...
    strcmp(manifest.precision,'FP64') && strcmp(manifest.array_order,'numpy'), ...
    'radia:gpu:exchangeSchema','Unsupported AMS HDF5 exchange contract');
data=read_node(path,manifest.tree);
end

function value=read_node(path,node)
switch node.kind
    case 'struct'
        value=struct(); names=fieldnames(node.fields);
        for k=1:numel(names)
            value.(names{k})=read_node(path,node.fields.(names{k}));
        end
    case 'cell'
        value=cell(1,numel(node.items));
        for k=1:numel(value)
            if iscell(node.items); item=node.items{k}; else; item=node.items(k); end
            value{k}=read_node(path,item);
        end
    case 'array'
        assert(~isempty(regexp(node.path,'^/arrays/a[0-9]+$','once')), ...
            'radia:gpu:exchangePath','Invalid AMS dataset path');
        assert(any(strcmp(node.dtype,{'double','float64','int32','uint8'})), ...
            'radia:gpu:exchangeType','Unsupported AMS numeric type');
        value=h5read(path,node.path);
        expected=node.dtype;if strcmp(expected,'float64');expected='double';end
        assert(strcmp(class(value),expected),'radia:gpu:exchangeType','AMS dataset type mismatch');
        shape=double(node.shape(:)');
        assert(numel(shape)>=1 && numel(shape)<=2 && all(shape>=0) && ...
            all(shape==fix(shape)) && numel(value)==prod(shape), ...
            'radia:gpu:exchangeShape','AMS dataset shape mismatch');
        if numel(shape)==1
            value=reshape(value,[],1);
        else
            assert(isequal(size(value),fliplr(shape)), ...
                'radia:gpu:exchangeShape','AMS dataset dimensions mismatch');
            value=permute(value,[2 1]);
        end
    case 'logical'
        value=logical(node.value);
    case {'number','string'}
        value=node.value;
    otherwise
        error('radia:gpu:exchangeNode','Unsupported AMS HDF5 node');
end
end
