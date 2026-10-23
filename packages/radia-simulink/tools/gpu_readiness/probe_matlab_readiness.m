function encoded = probe_matlab_readiness(expectedRoot, expectedGPU)
% Read-only MATLAB/license/GPU readiness in the actual runner account.
% Starts no pool, kernel, model, build or benchmark; never resets the GPU.
record = struct('schema','radia.runner-matlab-readiness.v1','passed',false, ...
    'matlab_version',version,'release',version('-release'),'matlabroot',matlabroot, ...
    'pid',feature('getpid'),'username',getenv('USERNAME'), ...
    'products',ver,'licenses',struct,'error','');
try
    actual = char(java.io.File(matlabroot).getCanonicalPath());
    expected = char(java.io.File(char(expectedRoot)).getCanonicalPath());
    assert(strcmpi(actual,expected),'Unexpected MATLAB installation.');
    assert(strcmp(version('-release'),'2026a'),'MATLAB R2026a required.');
    features = {'MATLAB','Simulink','Distrib_Computing_Toolbox'};
    for index = 1:numel(features)
        featureName = features{index};
        exists = license('test',featureName);
        [checkedOut, message] = license('checkout',featureName);
        record.licenses.(featureName) = struct('exists',logical(exists), ...
            'checkout',logical(checkedOut),'message',message);
        assert(exists && checkedOut,'License checkout failed: %s',featureName);
    end
    device = gpuDevice; % NO index argument: gpuDevice(1) resets the device.
    record.gpu = struct('name',device.Name,'uuid',device.UUID, ...
        'compute_capability',device.ComputeCapability, ...
        'graphics_driver',device.GraphicsDriverVersion, ...
        'toolkit',device.ToolkitVersion,'total_memory_bytes',device.TotalMemory);
    assert(contains(device.Name,char(expectedGPU),'IgnoreCase',true),'Unexpected GPU.');
    record.passed = true;
catch failure
    record.error = struct('identifier',failure.identifier,'message',failure.message);
end
encoded = jsonencode(record);
fprintf('RADIA_READINESS:%s\n',encoded);
end
