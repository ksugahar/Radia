function value = decodeTypedTree(tree, label)
%DECODETYPEDTREE Rebuild a MATLAB value from encodeTypedTree arrays.
%   Every structural inconsistency raises radia:optuna:StorageFormat.
arguments
    tree (1,1) struct
    label (1,1) string = "value"
end
classes = radia.optuna.internal.TypedHdf5.ClassNames;
count = numel(tree.node_class);
lengthsOk = numel(tree.node_parent) == count && ...
    numel(tree.node_shape_offset) == count && ...
    numel(tree.node_ndims) == count && ...
    numel(tree.node_data_offset) == count && ...
    numel(tree.node_data_count) == count;
if count == 0 || ~lengthsOk
    invalid("node arrays are empty or have inconsistent lengths");
end
textLengths = double(tree.text_lengths);
if any(textLengths < -1)
    invalid("text_lengths contains values below -1");
end
textByteOffsets = [0; cumsum(max(textLengths, 0))];
if textByteOffsets(end) ~= numel(tree.text_bytes)
    invalid("text_bytes length does not match text_lengths");
end
if numel(tree.datetime_ms_hi) ~= numel(tree.datetime_ms_lo)
    invalid("datetime pools differ in length");
end
cursor = 0;
value = take(0);
if cursor ~= count
    invalid(sprintf("%d trailing nodes", count - cursor));
end

    function item = take(parent)
        cursor = cursor + 1;
        index = cursor;
        if index > count
            invalid("node list ends early");
        end
        if double(tree.node_parent(index)) ~= parent
            invalid(sprintf("node %d has parent %d; expected %d", index, ...
                tree.node_parent(index), parent));
        end
        code = double(tree.node_class(index));
        if code < 1 || code > numel(classes)
            invalid(sprintf("node %d has class code %d", index, code));
        end
        className = classes(code);
        shapeOffset = double(tree.node_shape_offset(index));
        ndims = double(tree.node_ndims(index));
        if ndims < 2 || shapeOffset < 0 || shapeOffset + ndims > numel(tree.shape)
            invalid(sprintf("node %d has an invalid shape reference", index));
        end
        shape = double(reshape(tree.shape(shapeOffset+1:shapeOffset+ndims), 1, []));
        if any(shape < 0)
            invalid(sprintf("node %d has a negative dimension", index));
        end
        elements = prod(shape);
        offset = double(tree.node_data_offset(index));
        dataCount = double(tree.node_data_count(index));
        switch className
            case {'double','single','int8','uint8','int16','uint16', ...
                    'int32','uint32','int64','uint64','logical'}
                data = slice(className, offset, dataCount, elements, index);
                if className == "logical"
                    if any(data > 1)
                        invalid(sprintf("logical node %d is not 0/1", index));
                    end
                    data = logical(data);
                end
                item = reshape(data, shape);
            case {'char','string'}
                if offset < 0 || offset + dataCount > numel(textLengths) || ...
                        (className == "char" && dataCount ~= 1) || ...
                        (className == "string" && dataCount ~= elements)
                    invalid(sprintf("text node %d has an invalid range", index));
                end
                if className == "char"
                    text = textAt(offset + 1, index);
                    if numel(text) ~= elements
                        invalid(sprintf("char node %d decodes to %d characters", ...
                            index, numel(text)));
                    end
                    item = reshape(text, shape);
                else
                    item = strings(shape);
                    for k = 1:elements
                        position = offset + k;
                        if textLengths(position) < 0
                            item(k) = missing;
                        else
                            item(k) = string(textAt(position, index));
                        end
                    end
                end
            case 'datetime'
                hi = slice("datetime_ms_hi", offset, dataCount, elements, index);
                lo = slice("datetime_ms_lo", offset, dataCount, elements, index);
                timezone = take(index);
                format = take(index);
                if ~isstring(timezone) || ~isscalar(timezone) || ...
                        ~isstring(format) || ~isscalar(format)
                    invalid(sprintf("datetime node %d lacks TimeZone/Format", index));
                end
                item = radia.optuna.internal.TypedHdf5.datetimeFromParts( ...
                    hi, lo, timezone, shape);
                item.Format = format;
            case 'struct'
                requireNoData(index, offset, dataCount);
                fields = take(index);
                if ~isstring(fields) || any(ismissing(fields)) || ...
                        numel(unique(fields)) ~= numel(fields) || ...
                        ~all(arrayfun(@isvarname, fields))
                    invalid(sprintf("struct node %d has invalid field names", index));
                end
                if isempty(fields)
                    item = repmat(struct(), shape);
                else
                    item = cell2struct(cell([numel(fields), shape]), ...
                        cellstr(fields), 1);
                end
                for k = 1:elements
                    for field = reshape(fields, 1, [])
                        item(k).(field) = take(index);
                    end
                end
            case 'cell'
                requireNoData(index, offset, dataCount);
                item = cell(shape);
                for k = 1:elements
                    item{k} = take(index);
                end
            case 'table'
                requireNoData(index, offset, dataCount);
                columns = take(index);
                if numel(shape) ~= 2 || ~isstring(columns) || ...
                        numel(columns) ~= shape(2)
                    invalid(sprintf("table node %d has an invalid column list", index));
                end
                data = cell(1, numel(columns));
                for k = 1:numel(columns)
                    data{k} = take(index);
                    if size(data{k}, 1) ~= shape(1)
                        invalid(sprintf("table node %d column %s has %d rows", ...
                            index, columns(k), size(data{k}, 1)));
                    end
                end
                item = table(data{:}, 'VariableNames', cellstr(columns));
        end
    end

    function data = slice(poolName, offset, dataCount, elements, index)
        if dataCount ~= elements
            invalid(sprintf("node %d has %d elements; shape needs %d", ...
                index, dataCount, elements));
        end
        pool = tree.(poolName);
        if elements == 0
            data = pool([]);
            return
        end
        if offset < 0 || offset + dataCount > numel(pool)
            invalid(sprintf("node %d exceeds pool %s", index, poolName));
        end
        data = pool(offset+1:offset+dataCount);
    end

    function text = textAt(position, index)
        if textLengths(position) < 0
            invalid(sprintf("char node %d is marked missing", index));
        end
        bytes = tree.text_bytes(textByteOffsets(position)+1:textByteOffsets(position+1));
        text = reshape(native2unicode(reshape(bytes, 1, []), 'UTF-8'), 1, []);
    end

    function requireNoData(index, offset, dataCount)
        if offset ~= -1 || dataCount ~= 0
            invalid(sprintf("container node %d carries element data", index));
        end
    end

    function invalid(message)
        error("radia:optuna:StorageFormat", ...
            "Typed HDF5 tree '%s' is invalid: %s.", label, message);
    end
end
