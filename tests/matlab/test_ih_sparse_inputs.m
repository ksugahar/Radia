classdef test_ih_sparse_inputs < matlab.unittest.TestCase
    % Sparse storage must never be read as numel contiguous doubles.
    methods (Test)
        function validatorRejectsSparseOperator(testCase)
            cfg = radia.simulink.makeIHNativeSmokeConfig();
            cfg.eddy_matrix_real = sparse(cfg.eddy_matrix_real);
            testCase.verifyError(@() radia.simulink.validateIHNativeConfig(cfg), ...
                'radia:simulink:IHConfigMatrix');
        end

        function validatorRejectsSparseCSRValues(testCase)
            cfg = radia.simulink.makeIHNativeSmokeConfig();
            cfg.stiffness_value = sparse(cfg.stiffness_value);
            testCase.verifyError(@() radia.simulink.validateIHNativeConfig(cfg), ...
                'radia:simulink:IHConfigVector');
        end

        function mexRejectsSparseOperatorWithoutValidator(testCase)
            cfg = radia.simulink.makeIHNativeSmokeConfig();
            cfg.eddy_matrix_real = sparse(cfg.eddy_matrix_real);
            testCase.verifyError(@() createAndDestroy('eddy',cfg), ...
                'radia:mex:Exception');
        end

        function mexRejectsEmptySparseStorageWithoutValidator(testCase)
            cfg = radia.simulink.makeIHNativeSmokeConfig();
            cfg.stiffness_value = sparse(1,1);
            testCase.verifyError(@() createAndDestroy('thermal',cfg), ...
                'radia:mex:Exception');
        end
    end
end

function createAndDestroy(kind,cfg)
% Request the required output so an arity error cannot satisfy the test.
handle = radia_mex(['ih.' kind '.create'],cfg);
radia_mex(['ih.' kind '.destroy'],handle);
end
