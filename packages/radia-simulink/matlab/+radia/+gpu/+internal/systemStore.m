function result = systemStore(action, key, value)
%SYSTEMSTORE Own only the GPU arrays uploaded by an active CUDA PCG block.
persistent systems
if isempty(systems), systems = containers.Map('KeyType','char','ValueType','any'); end
result = [];
switch action
    case 'create'
        if isKey(systems,key), error('radia:gpu:LiveBlock','CUDA block is already active: %s',key); end
        systems(key) = value;
    case 'update'
        if ~isKey(systems,key), error('radia:gpu:DeadBlock','CUDA block is not active: %s',key); end
        systems(key) = value;
    case 'get'
        if ~isKey(systems,key), error('radia:gpu:DeadBlock','CUDA block is not active: %s',key); end
        result = systems(key);
    case 'remove'
        if isKey(systems,key), result = systems(key); remove(systems,key); end
    case 'count'
        result = systems.Count;
    otherwise
        error('radia:gpu:StoreAction','Unknown CUDA block lifecycle action.');
end
end
