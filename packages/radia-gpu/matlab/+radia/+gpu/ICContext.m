classdef ICContext < handle
    %ICCONTEXT Retain a CUDA IC(0) factor and FP64 PCG workspace.
    properties (SetAccess=private)
        SetupInfo struct
    end
    properties (Access=private)
        NativeHandle uint64 = uint64(0)
    end
    methods
        function self=ICContext(operator,options)
            arguments
                operator (1,1) struct
                options.Shift (1,1) double {mustBeFinite,mustBeInRange(options.Shift,1,100)} = 1.05
            end
            [self.NativeHandle,self.SetupInfo]=radia_cuda_ic('create',operator,options.Shift);
        end
        function [x,info]=solve(self,rhs,options)
            arguments
                self
                rhs
                options.Tolerance (1,1) double {mustBeFinite,mustBePositive,mustBeLessThanOrEqual(options.Tolerance,1e-6)} = 1e-10
                options.MaxIterations (1,1) double {mustBeInteger,mustBeInRange(options.MaxIterations,1,100000)} = 20000
                options.ReturnOnGPU (1,1) logical = true
            end
            assert(self.NativeHandle~=0,'radia:cuda:released','IC context is released');
            [x,info]=radia_cuda_ic('solve',self.NativeHandle,rhs,options.Tolerance,options.MaxIterations);
            if ~options.ReturnOnGPU;x=gather(x);end
        end
        function x=apply(self,rhs)
            assert(self.NativeHandle~=0,'radia:cuda:released','IC context is released');
            x=radia_cuda_ic('apply',self.NativeHandle,rhs);
        end
        function delete(self)
            if self.NativeHandle~=0
                radia_cuda_ic('release',self.NativeHandle);self.NativeHandle=uint64(0);
            end
        end
    end
end
