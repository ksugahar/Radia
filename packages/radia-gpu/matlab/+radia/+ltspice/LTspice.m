classdef LTspice
    %LTSPICE Simulator discovery and batch-execution compatibility API.
    methods (Static)
        function answer = isAvailable()
            try
                radia.ltspice.findExecutable();
                answer = true;
            catch
                answer = false;
            end
        end

        function answer = is_available()
            answer = radia.ltspice.LTspice.isAvailable();
        end

        function path = createNetlist(circuitFile, varargin)
            result = radia.ltspice.schematicToNetlist(string(circuitFile));
            path = string(result.netlist);
        end

        function path = create_netlist(varargin)
            path = radia.ltspice.LTspice.createNetlist(varargin{:});
        end

        function code = run(netlistFile, varargin)
            radia.ltspice.run(string(netlistFile));
            code = 0;
        end

        function paths = getDefaultLibraryPaths()
            paths = strings(0, 1);
        end

        function paths = get_default_library_paths()
            paths = radia.ltspice.LTspice.getDefaultLibraryPaths();
        end

        function answer = using_macos_native_sim()
            answer = false;
        end

        function simulator = create_from(path, varargin)
            simulator = struct('executable', string(path), 'process_name', "");
            if ~isempty(varargin)
                simulator.process_name = string(varargin{1});
            end
        end

        function path = expand_and_check_local_dir(path, varargin)
            path = string(path);
            if ~isfolder(path)
                path = "";
            end
        end

        function name = guess_process_name(executable)
            [~, name, extension] = fileparts(executable);
            name = name + extension;
        end

        function switches = valid_switch(switchName, path)
            if nargin < 2
                path = "";
            end
            switches = string(switchName);
            if strlength(path) > 0
                switches(2) = string(path);
            end
        end
    end
end
