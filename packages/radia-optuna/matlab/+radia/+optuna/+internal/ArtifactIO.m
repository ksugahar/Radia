classdef ArtifactIO
    %ARTIFACTIO Binary-array boundary shared by MATLAB artifact stores.

    methods (Static)
        function bytes=readFile(path)
            [file,message]=fopen(path,"rb");
            if file<0
                error("radia:optuna:ArtifactRead","%s",message);
            end
            cleanup=onCleanup(@()fclose(file));
            bytes=fread(file,Inf,"*uint8");
            clear cleanup
        end

        function writeFile(path,bytes,exclusive)
            if nargin<3, exclusive=false; end
            if exclusive && isfile(path)
                error("radia:optuna:ArtifactFileExists", ...
                    "File already exists: %s",path);
            end
            [file,message]=fopen(path,"wb");
            if file<0
                error("radia:optuna:ArtifactWrite","%s",message);
            end
            cleanup=onCleanup(@()fclose(file));
            count=fwrite(file,reshape(uint8(bytes),[],1),"uint8");
            if count~=numel(bytes)
                error("radia:optuna:ArtifactWrite", ...
                    "The complete artifact could not be written.");
            end
            clear cleanup
        end

        function bytes=contentBytes(content)
            if isa(content,"uint8")
                bytes=reshape(content,[],1);
            elseif ischar(content) || (isstring(content) && isscalar(content))
                bytes=radia.optuna.internal.ArtifactIO.readFile(string(content));
            else
                error("radia:optuna:ArtifactContent", ...
                    "Artifact content must be uint8 data or a scalar file path.");
            end
        end
    end
end
