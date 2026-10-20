classdef AMSContext < handle
    %AMSCONTEXT Own a checked native AMS handle and persistent CUDA PCG workspace.
    % Explicitly update the hierarchy; changing the operator does not rebuild it.
    properties (Access=private)
        NativeHandle uint64 = uint64(0)
    end
    methods
        function self=AMSContext(hierarchy)
            self.NativeHandle=radia_cuda_ams('create',hierarchy);
        end
        function update(self,hierarchy)
            assert(self.NativeHandle~=0,'radia:cuda:released','AMS context is released');
            radia_cuda_ams('update',self.NativeHandle,hierarchy);
        end
        function [x,info]=solve(self,operator,rhs,tolerance,maxIterations,useGraph)
            if nargin<4;tolerance=1e-9;end
            if nargin<5;maxIterations=2000;end
            if nargin<6;useGraph=true;end
            assert(self.NativeHandle~=0,'radia:cuda:released','AMS context is released');
            [x,info]=radia_cuda_ams('solve',self.NativeHandle,operator,rhs,tolerance,maxIterations,useGraph);
        end
        function delete(self)
            if self.NativeHandle~=0
                radia_cuda_ams('release',self.NativeHandle);
                self.NativeHandle=uint64(0);
            end
        end
    end
end
