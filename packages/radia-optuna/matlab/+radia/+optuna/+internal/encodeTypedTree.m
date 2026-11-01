function tree = encodeTypedTree(value, label)
%ENCODETYPEDTREE Flatten one MATLAB value into typed node/pool arrays.
%   Nodes are emitted in preorder. Each node records its class code, parent
%   (zero for the root, nodes one based), MATLAB shape (shape pool) and the
%   offset/count of its elements in the pool owned by its class. Children:
%     struct  1 string node of field names, then numel*F values ordered
%             element-major (linear element index, then field order)
%     cell    numel values in linear order
%     table   1 string node of column names, then one node per column
%     datetime 2 scalar string nodes: TimeZone, Format
%   char/string text is UTF-8 in text_bytes; text_lengths holds the byte
%   length of each element (-1 for a missing string). datetime elements are
%   the exact double-double milliseconds since 1970-01-01 (UTC for zoned
%   values) in datetime_ms_hi/datetime_ms_lo. Unsupported classes,
%   function handles, objects, complex and sparse values raise
%   radia:optuna:StorageUnsupported before any file is written.
arguments
    value
    label (1,1) string = "value"
end
classes = radia.optuna.internal.TypedHdf5.ClassNames;
pools = radia.optuna.internal.TypedHdf5.PoolNames;
capacity = 64;
nodeClass = zeros(capacity, 1);
nodeParent = zeros(capacity, 1);
nodeShapeOffset = zeros(capacity, 1);
nodeShapeCount = zeros(capacity, 1);
nodeDataOffset = zeros(capacity, 1);
nodeDataCount = zeros(capacity, 1);
nodeCount = 0;
poolData = cell(numel(pools), 1);
for q = 1:numel(pools)
    poolData{q} = {};
end
poolSize = zeros(numel(pools), 1);
shapeChunks = {};
shapeSize = 0;
textLengths = {};
textLengthSize = 0;

emit(value, 0);

tree = struct();
tree.node_class = uint8(nodeClass(1:nodeCount));
tree.node_parent = int64(nodeParent(1:nodeCount));
tree.node_shape_offset = int64(nodeShapeOffset(1:nodeCount));
tree.node_ndims = int64(nodeShapeCount(1:nodeCount));
tree.node_data_offset = int64(nodeDataOffset(1:nodeCount));
tree.node_data_count = int64(nodeDataCount(1:nodeCount));
tree.shape = int64(vertcat(zeros(0,1), shapeChunks{:}));
tree.text_lengths = int64(vertcat(zeros(0,1), textLengths{:}));
for q = 1:numel(pools)
    poolClass = radia.optuna.internal.TypedHdf5.PoolClasses(q);
    tree.(pools(q)) = cast(vertcat(zeros(0,1,poolClass), poolData{q}{:}), poolClass);
end

    function emit(item, parent)
        index = newNode(parent, size(item));
        if isa(item, "double") && isreal(item) && ~issparse(item)
            setData(index, "double", "double", item(:));
        elseif isa(item, "function_handle")
            unsupported(label, index, "function handle");
        elseif islogical(item) && ~issparse(item)
            setData(index, "logical", "logical", uint8(item(:)));
        elseif isnumeric(item) && ~isobject(item)
            className = string(class(item));
            if ~any(classes(1:10) == className) || ~isreal(item) || issparse(item)
                unsupported(label, index, className + " (complex/sparse/unknown)");
            end
            setData(index, className, className, item(:));
        elseif ischar(item)
            if ~(ismatrix(item) && (size(item,1) == 1 || isempty(item)))
                unsupported(label, index, "char matrix " + mat2str(size(item)));
            end
            bytes = utf8(item);
            setText(index, "char", {bytes}, numel(bytes));
        elseif isstring(item)
            bytes = cell(numel(item), 1);
            lengths = zeros(numel(item), 1);
            for k = 1:numel(item)
                if ismissing(item(k))
                    bytes{k} = zeros(0, 1, "uint8");
                    lengths(k) = -1;
                else
                    bytes{k} = utf8(item(k));
                    lengths(k) = numel(bytes{k});
                end
            end
            setText(index, "string", bytes, lengths);
        elseif isdatetime(item)
            [hi, lo] = datetimeParts(item);
            restored = radia.optuna.internal.TypedHdf5.datetimeFromParts( ...
                hi, lo, string(item.TimeZone), size(item));
            if ~isequaln(restored, item)
                unsupported(label, index, "datetime without an exact representation");
            end
            setData(index, "datetime", "datetime_ms_hi", hi);
            appendPool("datetime_ms_lo", lo);
            emit(string(item.TimeZone), index);
            emit(string(item.Format), index);
        elseif isstruct(item)
            nodeClass(index) = classCode("struct");
            fields = string(fieldnames(item));
            emit(reshape(fields, 1, []), index);
            for k = 1:numel(item)
                for field = reshape(fields, 1, [])
                    emit(item(k).(field), index);
                end
            end
        elseif iscell(item)
            nodeClass(index) = classCode("cell");
            for k = 1:numel(item)
                emit(item{k}, index);
            end
        elseif istable(item)
            properties = item.Properties;
            if ~isempty(properties.RowNames) || ...
                    any(strlength(string(properties.VariableUnits)) > 0) || ...
                    any(strlength(string(properties.VariableDescriptions)) > 0)
                unsupported(label, index, "table row names/units/descriptions");
            end
            nodeClass(index) = classCode("table");
            columns = string(properties.VariableNames);
            emit(reshape(columns, 1, []), index);
            for k = 1:numel(columns)
                emit(item.(columns(k)), index);
            end
        else
            unsupported(label, index, string(class(item)));
        end
    end

    function index = newNode(parent, shape)
        nodeCount = nodeCount + 1;
        index = nodeCount;
        if index > numel(nodeClass)
            grow = numel(nodeClass);
            nodeClass(end+grow) = 0;
            nodeParent(end+grow) = 0;
            nodeShapeOffset(end+grow) = 0;
            nodeShapeCount(end+grow) = 0;
            nodeDataOffset(end+grow) = 0;
            nodeDataCount(end+grow) = 0;
        end
        nodeParent(index) = parent;
        nodeShapeOffset(index) = shapeSize;
        nodeShapeCount(index) = numel(shape);
        shapeChunks{end+1, 1} = shape(:);
        shapeSize = shapeSize + numel(shape);
        nodeDataOffset(index) = -1;
        nodeDataCount(index) = 0;
    end

    function setData(index, className, poolName, data)
        nodeClass(index) = classCode(className);
        nodeDataOffset(index) = appendPool(poolName, data);
        nodeDataCount(index) = numel(data);
    end

    function setText(index, className, bytes, lengths)
        % data offset/count address text_lengths; byte offsets follow from
        % the cumulative lengths of the preceding text elements.
        nodeClass(index) = classCode(className);
        nodeDataOffset(index) = textLengthSize;
        nodeDataCount(index) = numel(lengths);
        textLengths{end+1, 1} = lengths(:);
        textLengthSize = textLengthSize + numel(lengths);
        appendPool("text_bytes", vertcat(zeros(0,1,"uint8"), bytes{:}));
    end

    function offset = appendPool(poolName, data)
        p = find(pools == poolName, 1);
        offset = poolSize(p);
        poolData{p}{end+1, 1} = data(:);
        poolSize(p) = poolSize(p) + numel(data);
    end

    function code = classCode(className)
        code = find(classes == className, 1);
    end
end

function unsupported(label, index, what)
error("radia:optuna:StorageUnsupported", ...
    "Value '%s' node %d (%s) has no HDF5 storage mapping.", ...
    label, index, what);
end

function bytes = utf8(text)
bytes = reshape(uint8(unicode2native(char(text), 'UTF-8')), [], 1);
end

function [hi, lo] = datetimeParts(value)
% datetime keeps a double-double millisecond count; the public accessors
% round it, so read both words to store the exact value.
state = warning('off', 'MATLAB:structOnObject');
restore = onCleanup(@() warning(state));
raw = struct(value).data;
hi = reshape(real(raw), [], 1);
lo = reshape(imag(raw), [], 1);
end
