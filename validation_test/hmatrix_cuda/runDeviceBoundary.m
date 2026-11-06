function runDeviceBoundary(mode,checkout,mexDir,work,python,faultDir)
% One step of run_device_boundary.py in its own batch MATLAB: the same-kernel MEX
% record, the setup read controls, the MATLAB device workflow test, or the fault sweep.
addpath(mexDir,fullfile(checkout,'packages','radia-simulink','matlab'), ...
    fullfile(checkout,'validation_test','hmatrix_cuda'),fullfile(checkout,'tests','matlab'));
loaded=strrep(which('radia_hdiv_device'),'\','/');
assert(startsWith(loaded,strrep(mexDir,'\','/'),IgnoreCase=true),'Wrong MEX: %s',loaded);
switch mode
    case 'parity'
        produceHDivDeviceMexRecord(fullfile(work,'input.h5'),fullfile(work,'mex.h5'),fullfile(work,'common.json'));
    case 'parity-mixed'
        produceHDivDeviceMexRecord(fullfile(work,'input-mixed.h5'),fullfile(work,'mex-mixed.h5'),fullfile(work,'common.json'));
    case 'setup-read'
        results=runtests('test_hdiv_gpu_workflow_config');
        assert(~isempty(results) && all([results.Passed]),'Setup read controls did not all pass');
        names=string({results.Name});
        assert(any(endsWith(names,'testWarmAndColdReadsAgree')),'The warm and cold read control did not run');
        writeReport(fullfile(work,'setup-read.json'),struct('accepted',true,'cases',numel(results)));
    case 'device-test'
        writeReport(fullfile(work,'device-test.json'),test_prepared_hdiv_device_setup(checkedSetup(checkout,python,'hex_pair'),softIron(),fullfile(work,'device-test.h5')));
    case 'device-test-mixed'
        writeReport(fullfile(work,'device-test-mixed.json'),test_prepared_hdiv_device_setup(checkedSetup(checkout,python,'tet_pyramid_hex'),softIron(),fullfile(work,'device-test-mixed.h5')));
    case 'device-fault'
        addpath(faultDir,fullfile(checkout,'validation_test','hmatrix_cuda','device_fault'));
        flat=readHDivParityInput(fullfile(work,'fault-input.h5'));
        mixed=readHDivParityInput(fullfile(work,'fault-input-mixed.h5'));
        setup=checkedSetup(checkout,python,'hex_pair');
        hex=struct('name','hex-pair','dofs',setup.dofs,'charges',setup.charges,'evaluation',setup.evaluation, ...
            'law',radia.gpu.makeHDivBHMaterial(softIron()),'apply',sin((1:setup.dofs)'),'rhs',setup.rhs_operator*[1000;-2000;5000], ...
            'initial',zeros(setup.dofs,1),'inv_chi',.1,'tolerance',1e-11,'maximum_iterations',2000, ...
            'preconditioner','mass-cholesky','newton',flat.newton);
        writeReport(fullfile(work,'device-fault.json'),checkDeviceFaults([orderfields(flat) orderfields(hex) orderfields(mixed)],fullfile(work,'device-fault-trace.txt')));
    otherwise
        error('radia:gpu:DeviceBoundary','Unknown mode %s',mode);
end
end

function setup=checkedSetup(checkout,python,name)
fixtures=fullfile(checkout,'tests','fixtures','hdiv');
receipt=fullfile(fixtures,[name '.receipt.json']);
setup=radia.gpu.prepareHDivDeviceSetup(fullfile(fixtures,[name '.vol']),2, ...
    LabelContract=fullfile(fixtures,[name '.labels.json']),MeshCheckReceipt=receipt, ...
    MeshCheckReceiptSha256=sha256File(receipt), ...
    MeshCheckProvenance=jsondecode(fileread([receipt '.provenance.json'])), ...
    PythonExecutable=string(python));
end

function bh=softIron()
% Same soft-iron rows as the parity case.
bh=[0 0;100 .5;300 1;1000 1.4;5000 1.7;20000 1.9;100000 2.1];
end

function writeReport(path,report)
stream=fopen(path,'w');assert(stream>0,'radia:gpu:DeviceBoundary','Cannot write %s',path);
closer=onCleanup(@()fclose(stream));fwrite(stream,jsonencode(report),'char');
end

function digest=sha256File(path)
stream=fopen(path,'r');assert(stream>0,'radia:gpu:DeviceBoundary','Cannot read %s',path);
closer=onCleanup(@()fclose(stream));bytes=fread(stream,inf,'*uint8');
engine=java.security.MessageDigest.getInstance('SHA-256');engine.update(typecast(bytes,'int8'));
digest=lower(reshape(dec2hex(typecast(engine.digest(),'uint8'),2)',1,[]));
end
