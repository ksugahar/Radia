classdef test_ih_sparse_inputs < matlab.unittest.TestCase
    % Sparse storage must never be read as numel contiguous doubles.
    methods (Test)
        function allocatorSurvivesOrdinaryClear(testCase)
            first = radia.internal.nextIHNativeHandle();
            clear('radia.internal.nextIHNativeHandle');
            second = radia.internal.nextIHNativeHandle();
            testCase.verifyClass(second,'uint64');
            testCase.verifyGreaterThan(second,first);
            testCase.verifyTrue(mislocked('radia.internal.nextIHNativeHandle'));
        end

        function destroyedIdentityStaysStaleAfterReload(testCase)
            % Dedicated CI only: never unload a borrowed session's MEX.
            testCase.assumeEqual(getenv('RADIA_MEX_RELOAD_TEST'),'1');
            cfg = radia.simulink.makeIHNativeSmokeConfig();
            old = radia_mex('ih.thermal.create',cfg);
            radia_mex('ih.thermal.destroy',old);
            testCase.assumeFalse(mislocked('radia_mex'), ...
                'Other live native handles prevent an isolated reload test.');
            clear radia_mex
            fresh = radia_mex('ih.thermal.create',cfg);
            cleanup = onCleanup(@() radia_mex('ih.thermal.destroy',fresh)); %#ok<NASGU>
            testCase.verifyNotEqual(fresh,old);
            testCase.verifyError(@() radia_mex('ih.thermal.reset',old), ...
                'radia:mex:Exception');
            testCase.verifyError(@() radia_mex('ih.thermal.destroy',old), ...
                'radia:mex:Exception');
            radia_mex('ih.thermal.reset',fresh);
        end

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
