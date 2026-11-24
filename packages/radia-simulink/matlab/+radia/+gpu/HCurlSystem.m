classdef HCurlSystem < handle
    %HCURLSYSTEM Reuse an affine linear FE operator and its AMS hierarchy.
    % Geometry, labels, reluctivity and gauge mass are fixed for this object.
    % CurrentScale changes the declared source current without rebuilding.
    % Nonlinear material changes use solveHCurl and its resident hierarchy.
    properties (SetAccess=private)
        SetupInfo struct
    end
    properties (Access=private)
        Context uint64 = uint64(0)
        RHS
        Released logical = false
        Solves double = 0
    end
    methods
        function self=HCurlSystem(volPath,boundaryLabels,options)
            arguments
                volPath (1,1) string
                boundaryLabels (1,:) string
                options.Dirichlet (1,1) string = ".*"
                options.Reluctivity (1,1) double {mustBeFinite,mustBePositive} = 1
                options.MassCoefficient (1,1) double {mustBeFinite,mustBePositive} = 1
                options.CurrentDensity (1,3) double {mustBeFinite} = [1 0 0]
            end
            g=gpuDevice;wait(g);timer=tic;
            [fem,assembly]=radia.gpu.assembleHCurl(volPath,boundaryLabels, ...
                Dirichlet=options.Dirichlet,Reluctivity=options.Reluctivity, ...
                MassCoefficient=options.MassCoefficient,CurrentDensity=options.CurrentDensity,Gamma=0);
            h=radia_cuda_ams('nonlinear_create',fem);
            try
                [~,~,~]=radia_cuda_ams('nonlinear_evaluate',h,zeros(size(fem.rhs),'like',fem.rhs));
                setup=radia_cuda_ams('nonlinear_prepare_resident',h);
            catch e
                radia_cuda_ams('nonlinear_release',h);
                rethrow(e);
            end
            self.Context=h;self.RHS=fem.rhs;wait(g);
            self.SetupInfo=struct('assembly',assembly,'ams_setup',setup,'total_seconds',toc(timer), ...
                'backend','MATLAB/mexcuda GPU','geometry_and_coefficients_fixed',true);
        end
        function [x,info]=solve(self,options)
            arguments
                self
                options.CurrentScale (1,1) double {mustBeFinite} = 1
                options.Tolerance (1,1) double {mustBeFinite,mustBePositive,mustBeLessThanOrEqual(options.Tolerance,1e-6)} = 1e-10
                options.ReturnOnGPU (1,1) logical = true
            end
            assert(~self.Released,'radia:cuda:released','HCurl system is released');
            g=gpuDevice;wait(g);timer=tic;
            [x,linear]=radia_cuda_ams('nonlinear_solve_resident',self.Context, ...
                options.CurrentScale*self.RHS,options.Tolerance,4000);
            if ~options.ReturnOnGPU;x=gather(x);end
            wait(g);self.Solves=self.Solves+1;
            info=struct('linear',linear,'total_seconds',toc(timer),'solve_count',self.Solves, ...
                'matrix_reused',true,'hierarchy_reused',true,'current_scale',options.CurrentScale, ...
                'backend','MATLAB/mexcuda GPU');
        end
        function delete(self)
            if ~self.Released
                self.Released=true;
                if self.Context~=0;radia_cuda_ams('nonlinear_release',self.Context);end
                self.Context=uint64(0);self.RHS=[];
            end
        end
    end
end
