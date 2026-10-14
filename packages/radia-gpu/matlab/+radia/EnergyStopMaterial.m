classdef EnergyStopMaterial < handle
    %ENERGYSTOPMATERIAL Stateful B-input vector Stop material.

    properties (SetAccess = private)
        Alpha
        BMax
        NuBound
        Eta
        Gamma
        StateSize
    end

    properties (Access = private)
        NativeHandle = uint64(0)
    end

    methods
        function obj = EnergyStopMaterial(eta, gTables, options)
            arguments
                eta double
                gTables cell
                options.Alpha (1,1) double = 5.0
                options.Gamma double = 0.0
                options.BMax (1,1) double = 0.0
            end

            eta = double(eta(:).');
            if isempty(eta) || numel(gTables) ~= numel(eta)
                error("radia:EnergyStopMaterial:Tables", ...
                    "gTables must contain one table per eta value.");
            end
            gamma = double(options.Gamma(:).');
            if isscalar(gamma)
                gamma = repmat(gamma, size(eta));
            end
            if ~isequal(size(gamma), size(eta))
                error("radia:EnergyStopMaterial:Gamma", ...
                    "Gamma must be scalar or have the same size as eta.");
            end

            tableR = zeros(1, 0);
            tableG = zeros(1, 0);
            offsets = int32(0);
            for i = 1:numel(gTables)
                table = gTables{i};
                if iscell(table) && numel(table) == 2
                    r = double(table{1}(:).');
                    g = double(table{2}(:).');
                elseif isnumeric(table) && ismatrix(table) && size(table, 2) == 2
                    r = double(table(:, 1).');
                    g = double(table(:, 2).');
                else
                    error("radia:EnergyStopMaterial:TableShape", ...
                        "Each table must be an N-by-2 array or a {r,g} cell pair.");
                end
                if isempty(r) || ~isequal(size(r), size(g))
                    error("radia:EnergyStopMaterial:TableShape", ...
                        "Each r and g table must be nonempty and have matching sizes.");
                end
                tableR = [tableR, r]; %#ok<AGROW>
                tableG = [tableG, g]; %#ok<AGROW>
                offsets(end + 1) = int32(numel(tableR)); %#ok<AGROW>
            end

            obj.NativeHandle = radia.internal.callMex( ...
                'energy_stop.create', eta, tableR, tableG, offsets, gamma, ...
                options.Alpha, options.BMax);
            info = radia.internal.callMex('energy_stop.info', obj.NativeHandle);
            obj.Alpha = info.alpha;
            obj.BMax = info.b_max;
            obj.NuBound = info.nu_bound;
            obj.Eta = info.eta;
            obj.Gamma = info.gamma;
            obj.StateSize = info.state_size;
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex('energy_stop.destroy', obj.NativeHandle);
                catch
                end
                obj.NativeHandle = uint64(0);
            end
        end

        function state = state0(obj)
            obj.assertAlive();
            state = radia.internal.callMex('energy_stop.state0', obj.NativeHandle);
        end

        function H = forward(obj, B, states)
            obj.assertAlive();
            H = radia.internal.callMex( ...
                'energy_stop.forward', obj.NativeHandle, double(B), double(states));
        end

        function states = commit(obj, B, oldStates)
            obj.assertAlive();
            states = radia.internal.callMex( ...
                'energy_stop.commit', obj.NativeHandle, double(B), double(oldStates));
        end

        function energy = storedEnergy(obj, B, states)
            obj.assertAlive();
            energy = radia.internal.callMex( ...
                'energy_stop.stored_energy', obj.NativeHandle, double(B), double(states));
        end
    end

    methods (Access = private)
        function assertAlive(obj)
            if obj.NativeHandle == 0
                error("radia:EnergyStopMaterial:Deleted", ...
                    "The native EnergyStopMaterial has been deleted.");
            end
        end
    end
end
