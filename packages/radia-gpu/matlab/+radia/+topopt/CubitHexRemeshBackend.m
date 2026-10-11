classdef CubitHexRemeshBackend
    %CUBITHEXREMESHBACKEND Batch Cubit adapter for HEX topology commits.
    %   Journal generation and mesh loading remain application-owned callbacks.

    properties (SetAccess=private)
        Executable string
        WorkDirectory string
        WriteJournalFcn function_handle
        LoadMeshFcn function_handle
        Arguments string
    end

    methods
        function obj=CubitHexRemeshBackend(executable,workDirectory,writeJournalFcn,loadMeshFcn,options)
            arguments
                executable (1,1) string
                workDirectory (1,1) string
                writeJournalFcn (1,1) function_handle
                loadMeshFcn (1,1) function_handle
                options.Arguments (1,:) string=["-batch"]
            end
            obj.Executable=executable;
            obj.WorkDirectory=workDirectory;
            obj.WriteJournalFcn=writeJournalFcn;
            obj.LoadMeshFcn=loadMeshFcn;
            obj.Arguments=options.Arguments;
        end

        function mesh=rebuild(obj,request)
            if ~isfolder(obj.WorkDirectory), mkdir(obj.WorkDirectory); end
            obj.WriteJournalFcn(request);
            journal=string(request.journal_path);
            meshPath=string(request.mesh_path);
            if ~isfile(journal)
                error("radia:topopt:CubitJournal", ...
                    "Cubit journal was not created: %s",journal);
            end
            previous=pwd;
            cleanup=onCleanup(@()cd(previous));
            cd(obj.WorkDirectory);
            command=strjoin([quoteArgument(obj.Executable), ...
                arrayfun(@quoteArgument,obj.Arguments),quoteArgument(journal)]," ");
            [status,output]=system(command);
            if status~=0
                lines=splitlines(string(output));
                firstLine=max(1,numel(lines)-29);
                lines=lines(firstLine:end);
                error("radia:topopt:CubitFailed", ...
                    "Cubit remesh failed (%d):\n%s",status,strjoin(lines,newline));
            end
            if ~isfile(meshPath)
                error("radia:topopt:CubitMesh", ...
                    "Cubit did not create mesh: %s",meshPath);
            end
            mesh=obj.LoadMeshFcn(meshPath);
            if ~isa(mesh,"radia.ngsolve.Mesh")
                error("radia:topopt:CubitMesh", ...
                    "LoadMeshFcn must return a radia.ngsolve.Mesh.");
            end
            clear cleanup
        end
    end
end

function value=quoteArgument(value)
quote=string(char(34));
value=quote+replace(string(value),quote,quote+quote)+quote;
end
