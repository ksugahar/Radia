classdef HCurlICSystem < handle
    %HCURLICSYSTEM Fixed affine HCurl FE matrix with retained GPU IC(0).
    properties (SetAccess=private)
        SetupInfo struct
    end
    properties (Access=private)
        Context
        RHS
        Released logical = false
    end
    methods
        function self=HCurlICSystem(volPath,boundaryLabels,options)
            arguments
                volPath (1,1) string
                boundaryLabels (1,:) string
                options.Dirichlet (1,1) string = ".*"
                options.Reluctivity (1,1) double {mustBeFinite,mustBePositive} = 1
                options.MassCoefficient (1,1) double {mustBeFinite,mustBePositive} = 1
                options.CurrentDensity (1,3) double {mustBeFinite} = [1 0 0]
                options.Shift (1,1) double {mustBeFinite,mustBeInRange(options.Shift,1,100)} = 1.05
            end
            g=gpuDevice;wait(g);timer=tic;
            [~,fem,assembly]=radia.gpu.assembleHCurl(volPath,boundaryLabels, ...
                Dirichlet=options.Dirichlet,Reluctivity=options.Reluctivity,MassCoefficient=options.MassCoefficient, ...
                CurrentDensity=options.CurrentDensity,Gamma=0);
            h=radia_cuda_ams('nonlinear_create',fem);cleanup=onCleanup(@()radia_cuda_ams('nonlinear_release',h));
            [~,values,~]=radia_cuda_ams('nonlinear_evaluate',h,zeros(size(fem.rhs),'like',fem.rhs));
            operator=fem.mass;operator.values=values;operator.free=fem.free;self.RHS=fem.rhs;clear cleanup;
            wait(g);factorTimer=tic;self.Context=radia.gpu.ICContext(operator,Shift=options.Shift);wait(g);
            self.SetupInfo=struct('assembly',assembly,'ic_setup',self.Context.SetupInfo, ...
                'ic_setup_seconds',toc(factorTimer),'total_seconds',toc(timer), ...
                'backend','MATLAB/mexcuda GPU IC(0)','geometry_and_coefficients_fixed',true);
        end
        function [x,info]=solve(self,options)
            arguments
                self
                options.CurrentScale (1,1) double {mustBeFinite} = 1
                options.Tolerance (1,1) double {mustBeFinite,mustBePositive,mustBeLessThanOrEqual(options.Tolerance,1e-6)} = 1e-10
                options.MaxIterations (1,1) double {mustBeInteger,mustBeInRange(options.MaxIterations,1,100000)} = 20000
                options.ReturnOnGPU (1,1) logical = true
            end
            assert(~self.Released,'radia:cuda:released','HCurl IC system is released');
            g=gpuDevice;wait(g);timer=tic;
            [x,linear]=self.Context.solve(options.CurrentScale*self.RHS,Tolerance=options.Tolerance, ...
                MaxIterations=options.MaxIterations,ReturnOnGPU=options.ReturnOnGPU);wait(g);
            info=struct('linear',linear,'total_seconds',toc(timer),'matrix_reused',true, ...
                'factor_reused',true,'current_scale',options.CurrentScale,'backend','MATLAB/mexcuda GPU IC(0)');
        end
        function delete(self)
            if ~self.Released
                self.Released=true;
                if ~isempty(self.Context);delete(self.Context);end
                self.Context=[];self.RHS=[];
            end
        end
    end
end
