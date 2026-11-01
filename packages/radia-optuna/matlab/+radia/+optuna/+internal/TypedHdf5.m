classdef TypedHdf5
    %TYPEDHDF5 Explicit typed mapping between MATLAB data and HDF5 files.
    %
    %   A file carries root attributes schema (ASCII), schema_version
    %   (int64) and encoding. Each stored value is one HDF5 group holding
    %   flat one-dimensional datasets produced by encodeTypedTree:
    %
    %     node_class         uint8  1-based code into attribute class_names
    %     node_parent        int64  parent node (0 = root), preorder, 1 based
    %     node_shape_offset  int64  zero-based offset into shape
    %     node_ndims         int64  number of MATLAB dimensions
    %     node_data_offset   int64  zero-based offset into the class pool
    %                               (text_lengths for char/string, -1 for
    %                               struct/cell/table containers)
    %     node_data_count    int64  element count in that pool
    %     shape              int64  MATLAB dimensions (column-major data)
    %     text_lengths       int64  UTF-8 byte length per text element,
    %                               -1 marks a missing string
    %     text_bytes         uint8  concatenated UTF-8 bytes
    %     double ... uint64  native-dtype element pools, exact bits
    %     logical            uint8  0/1
    %     datetime_ms_hi/lo  float64 double-double milliseconds since
    %                               1970-01-01T00:00:00 (UTC when zoned)
    %
    %   No MATLAB serialized object, MAT payload or function handle is ever
    %   written; unsupported values fail with radia:optuna:StorageUnsupported.

    properties (Constant)
        Encoding = "radia.optuna.typed-tree.v1"
        ClassNames = ["double","single","int8","uint8","int16","uint16", ...
            "int32","uint32","int64","uint64","logical","char","string", ...
            "datetime","struct","cell","table"]
        PoolNames = ["double","single","int8","uint8","int16","uint16", ...
            "int32","uint32","int64","uint64","logical","text_bytes", ...
            "datetime_ms_hi","datetime_ms_lo"]
        PoolClasses = ["double","single","int8","uint8","int16","uint16", ...
            "int32","uint32","int64","uint64","uint8","uint8", ...
            "double","double"]
        NodeFields = ["node_class","node_parent","node_shape_offset", ...
            "node_ndims","node_data_offset","node_data_count", ...
            "shape","text_lengths"]
        NodeClasses = ["uint8","int64","int64","int64","int64","int64", ...
            "int64","int64"]
    end

    methods (Static)
        function fid = createFile(path, schema, version)
            fid = H5F.create(char(path), 'H5F_ACC_TRUNC', ...
                'H5P_DEFAULT', 'H5P_DEFAULT');
            try
                radia.optuna.internal.TypedHdf5.writeAttribute( ...
                    fid, "/", "schema", string(schema));
                radia.optuna.internal.TypedHdf5.writeAttribute( ...
                    fid, "/", "schema_version", int64(version));
                radia.optuna.internal.TypedHdf5.writeAttribute(fid, "/", ...
                    "encoding", radia.optuna.internal.TypedHdf5.Encoding);
            catch exception
                H5F.close(fid);
                rethrow(exception)
            end
        end

        function fid = openFile(path, schema, version)
            % Open read-only and fail loudly on any header mismatch.
            import radia.optuna.internal.TypedHdf5
            path = string(path);
            if ~isfile(path)
                error("radia:optuna:StorageFormat", ...
                    "Storage '%s' does not exist.", path);
            end
            if ~H5F.is_hdf5(char(path))
                error("radia:optuna:StorageFormat", ...
                    "Storage '%s' is not an HDF5 file.", path);
            end
            fid = H5F.open(char(path), 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
            try
                foundSchema = TypedHdf5.readAttribute(fid, "/", "schema");
                foundVersion = TypedHdf5.readAttribute(fid, "/", "schema_version");
                foundEncoding = TypedHdf5.readAttribute(fid, "/", "encoding");
            catch exception
                H5F.close(fid);
                failure = MException("radia:optuna:StorageVersion", ...
                    "Storage '%s' has no readable schema header: %s", ...
                    path, exception.message);
                throw(addCause(failure, exception))
            end
            if ~isstring(foundSchema) || ~isa(foundVersion, "int64") || ...
                    ~isscalar(foundVersion) || ~isstring(foundEncoding) || ...
                    foundSchema ~= schema || foundVersion ~= int64(version) || ...
                    foundEncoding ~= TypedHdf5.Encoding
                H5F.close(fid);
                error("radia:optuna:StorageVersion", ...
                    "Storage '%s' has unsupported schema/version/encoding " + ...
                    "'%s'/%s/'%s'; expected '%s'/%d/'%s'.", path, ...
                    join(string(foundSchema)), join(string(foundVersion)), ...
                    join(string(foundEncoding)), schema, version, ...
                    TypedHdf5.Encoding);
            end
        end

        function writeTree(fid, name, value)
            %WRITETREE Encode VALUE fully, then write it as group NAME.
            import radia.optuna.internal.TypedHdf5
            tree = radia.optuna.internal.encodeTypedTree(value, name);
            group = H5G.create(fid, char(name), 'H5P_DEFAULT', ...
                'H5P_DEFAULT', 'H5P_DEFAULT');
            closer = onCleanup(@() H5G.close(group));
            TypedHdf5.writeAttribute(fid, name, "class_names", ...
                strjoin(TypedHdf5.ClassNames, ","));
            for field = [TypedHdf5.NodeFields, TypedHdf5.PoolNames]
                TypedHdf5.writeDataset(group, field, tree.(field));
            end
        end

        function value = readTree(fid, name)
            %READTREE Read group NAME with exact dtype checks and decode it.
            import radia.optuna.internal.TypedHdf5
            if ~H5L.exists(fid, char(name), 'H5P_DEFAULT')
                error("radia:optuna:StorageFormat", ...
                    "HDF5 group '%s' is missing.", name);
            end
            if TypedHdf5.readAttribute(fid, name, "class_names") ~= ...
                    strjoin(TypedHdf5.ClassNames, ",")
                error("radia:optuna:StorageFormat", ...
                    "HDF5 group '%s' uses an unknown class table.", name);
            end
            group = H5G.open(fid, char(name));
            closer = onCleanup(@() H5G.close(group));
            tree = struct();
            names = [TypedHdf5.NodeFields, TypedHdf5.PoolNames];
            classes = [TypedHdf5.NodeClasses, TypedHdf5.PoolClasses];
            for index = 1:numel(names)
                tree.(names(index)) = TypedHdf5.readDataset( ...
                    group, names(index), classes(index));
            end
            value = radia.optuna.internal.decodeTypedTree(tree, name);
        end

        function writeAttribute(fid, objectPath, name, value)
            % ASCII text or int64 attribute on an existing group/dataset.
            object = H5O.open(fid, char(objectPath), 'H5P_DEFAULT');
            closer = onCleanup(@() H5O.close(object));
            if isstring(value) || ischar(value)
                payload = char(value);
                if any(double(payload) > 127)
                    error("radia:optuna:StorageUnsupported", ...
                        "HDF5 attribute '%s' must be ASCII.", name);
                end
                type = H5T.copy('H5T_C_S1');
                H5T.set_size(type, max(1, numel(payload)));
                H5T.set_strpad(type, 'H5T_STR_NULLTERM');
                space = H5S.create('H5S_SCALAR');
                if isempty(payload)
                    payload = char(0);
                end
            else
                type = H5T.copy('H5T_NATIVE_INT64');
                payload = int64(value(:));
                space = H5S.create_simple(1, numel(payload), []);
            end
            typeCloser = onCleanup(@() H5T.close(type));
            spaceCloser = onCleanup(@() H5S.close(space));
            attribute = H5A.create(object, char(name), type, space, ...
                'H5P_DEFAULT');
            attributeCloser = onCleanup(@() H5A.close(attribute));
            H5A.write(attribute, 'H5ML_DEFAULT', payload);
        end

        function value = readAttribute(fid, objectPath, name)
            object = H5O.open(fid, char(objectPath), 'H5P_DEFAULT');
            closer = onCleanup(@() H5O.close(object));
            try
                attribute = H5A.open(object, char(name));
            catch
                error("radia:optuna:StorageFormat", ...
                    "HDF5 object '%s' lacks the '%s' attribute.", ...
                    objectPath, name);
            end
            attributeCloser = onCleanup(@() H5A.close(attribute));
            type = H5A.get_type(attribute);
            typeCloser = onCleanup(@() H5T.close(type));
            raw = H5A.read(attribute);
            if H5T.get_class(type) == H5ML.get_constant_value('H5T_STRING')
                raw = reshape(char(raw), 1, []);
                terminator = find(raw == char(0), 1);
                if ~isempty(terminator)
                    raw = raw(1:terminator-1);
                end
                value = string(raw);
            else
                value = raw;
            end
        end

        function value = datetimeFromParts(hi, lo, timezone, shape)
            if strlength(timezone) == 0
                epoch = datetime(1970, 1, 1);
            else
                epoch = datetime(1970, 1, 1, 'TimeZone', 'UTC');
            end
            value = datetime(reshape(hi, [], 1), 'ConvertFrom', 'epochtime', ...
                'Epoch', epoch, 'TicksPerSecond', 1000, ...
                'TimeZone', char(timezone)) + milliseconds(reshape(lo, [], 1));
            value = reshape(value, shape);
        end
    end

    methods (Static, Access=private)
        function writeDataset(group, name, data)
            data = data(:);
            typeName = radia.optuna.internal.TypedHdf5.h5Type( ...
                string(class(data)));
            space = H5S.create_simple(1, numel(data), []);
            spaceCloser = onCleanup(@() H5S.close(space));
            dataset = H5D.create(group, char(name), typeName, space, ...
                'H5P_DEFAULT');
            datasetCloser = onCleanup(@() H5D.close(dataset));
            if ~isempty(data)
                H5D.write(dataset, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', ...
                    'H5P_DEFAULT', data);
            end
        end

        function data = readDataset(group, name, expectedClass)
            if ~H5L.exists(group, char(name), 'H5P_DEFAULT')
                error("radia:optuna:StorageFormat", ...
                    "HDF5 dataset '%s' is missing.", name);
            end
            dataset = H5D.open(group, char(name));
            closer = onCleanup(@() H5D.close(dataset));
            space = H5D.get_space(dataset);
            spaceCloser = onCleanup(@() H5S.close(space));
            [rank, dims] = H5S.get_simple_extent_dims(space);
            if rank ~= 1
                error("radia:optuna:StorageFormat", ...
                    "HDF5 dataset '%s' must be one dimensional.", name);
            end
            type = H5D.get_type(dataset);
            typeCloser = onCleanup(@() H5T.close(type));
            expectedType = H5T.copy(radia.optuna.internal.TypedHdf5.h5Type( ...
                string(expectedClass)));
            expectedCloser = onCleanup(@() H5T.close(expectedType));
            if ~H5T.equal(type, expectedType)
                error("radia:optuna:StorageFormat", ...
                    "HDF5 dataset '%s' does not have dtype %s.", ...
                    name, expectedClass);
            end
            if dims == 0
                data = zeros(0, 1, char(expectedClass));
            else
                data = reshape(H5D.read(dataset, 'H5ML_DEFAULT', 'H5S_ALL', ...
                    'H5S_ALL', 'H5P_DEFAULT'), [], 1);
            end
            if ~isa(data, expectedClass)
                error("radia:optuna:StorageFormat", ...
                    "HDF5 dataset '%s' decoded as %s; expected %s.", ...
                    name, class(data), expectedClass);
            end
        end

        function typeName = h5Type(className)
            switch className
                case "double", typeName = 'H5T_NATIVE_DOUBLE';
                case "single", typeName = 'H5T_NATIVE_FLOAT';
                case "int8", typeName = 'H5T_NATIVE_INT8';
                case "uint8", typeName = 'H5T_NATIVE_UINT8';
                case "int16", typeName = 'H5T_NATIVE_INT16';
                case "uint16", typeName = 'H5T_NATIVE_UINT16';
                case "int32", typeName = 'H5T_NATIVE_INT32';
                case "uint32", typeName = 'H5T_NATIVE_UINT32';
                case "int64", typeName = 'H5T_NATIVE_INT64';
                case "uint64", typeName = 'H5T_NATIVE_UINT64';
                otherwise
                    error("radia:optuna:StorageUnsupported", ...
                        "No HDF5 dtype for MATLAB class %s.", className);
            end
        end
    end
end
