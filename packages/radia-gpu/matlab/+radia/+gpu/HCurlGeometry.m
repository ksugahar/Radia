classdef HCurlGeometry < handle
    %HCURLGEOMETRY Explicit immutable affine p=1 FE snapshot and GPU mapped data.
    % Existing objects retain their snapshot if the input file later changes.
    % Construct another object for new geometry, labels or constraints.
    % Snapshot SHA256 uses MATLAB's JVM; -nojvm sessions fail explicitly.
    properties (SetAccess=private)
        BoundaryLabels string
        Dirichlet string
        MeshSHA256 string
        Info struct
    end
    properties (Access=private)
        NativeHandle uint64
        Released logical = true
    end
    methods
        function self=HCurlGeometry(volPath,boundaryLabels,options)
            arguments
                volPath (1,1) string
                boundaryLabels (1,:) string
                options.Dirichlet (1,1) string = ".*"
            end
            timer=tic;
            before=meshHash(volPath);
            [h,info]=radia_cuda_ams('geometry_create',char(volPath),cellstr(boundaryLabels),char(options.Dirichlet));
            try
                after=meshHash(volPath);
                assert(before==after,'radia:cuda:geometryChanged','Mesh changed during geometry construction');
                self.BoundaryLabels=boundaryLabels;self.Dirichlet=options.Dirichlet;
                self.MeshSHA256=after;info.mesh_sha256=after;
                info.total_creation_seconds=toc(timer);self.Info=info;
                self.NativeHandle=h;self.Released=false;
            catch failure
                radia_cuda_ams('geometry_release',h);rethrow(failure);
            end
        end
        function [model,fem,info]=assemble(self,options)
            arguments
                self
                options.Reluctivity (1,1) double {mustBeFinite,mustBePositive} = 1
                options.MassCoefficient (1,1) double {mustBeFinite,mustBePositive} = 1
                options.CurrentDensity (1,3) double {mustBeFinite} = [1 0 0]
                options.Gamma (1,1) double {mustBeFinite,mustBeNonnegative} = 0
            end
            assert(~self.Released,'radia:cuda:released','HCurl geometry is released');
            [model,fem,info]=radia_cuda_ams('assemble_hcurl_geometry',self.NativeHandle, ...
                options.Reluctivity,options.MassCoefficient,options.CurrentDensity,options.Gamma);
            info.mesh_sha256=self.MeshSHA256;
        end
        function delete(self)
            if ~self.Released
                radia_cuda_ams('geometry_release',self.NativeHandle);
                self.Released=true;
            end
        end
    end
end
function digest=meshHash(path)
assert(isfile(path),'radia:cuda:vol','Checked .vol path does not exist');
assert(usejava('jvm'),'radia:cuda:hash','Geometry snapshot hashing requires MATLAB JVM');
file=fopen(path,'rb');assert(file>=0,'radia:cuda:hash','Cannot read mesh snapshot');
cleanup=onCleanup(@()fclose(file));
md=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
assert(fseek(file,0,'eof')==0,'radia:cuda:hash','Cannot seek mesh snapshot');
remaining=ftell(file);
assert(remaining>=0 && fseek(file,0,'bof')==0,'radia:cuda:hash','Cannot size mesh snapshot');
while remaining>0
    requested=min(remaining,2^20);
    [bytes,count]=fread(file,requested,'*uint8');
    [~,readError]=ferror(file);
    assert(count==requested && readError==0,'radia:cuda:hash','Mesh snapshot read failed');
    md.update(typecast(bytes,'int8'));
    remaining=remaining-count;
end
digest=string(lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[])));
end
