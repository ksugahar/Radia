function runDeviceBoundary(mode,checkout,mexDir,work,python)
% One step of run_device_boundary.py in its own batch MATLAB: the same-kernel MEX
% record, the setup read controls, or the MATLAB device workflow test.
addpath(mexDir,fullfile(checkout,'packages','radia-simulink','matlab'), ...
    fullfile(checkout,'validation_test','hmatrix_cuda'),fullfile(checkout,'tests','matlab'));
loaded=strrep(which('radia_hdiv_device'),'\','/');
assert(startsWith(loaded,strrep(mexDir,'\','/'),IgnoreCase=true),'Wrong MEX: %s',loaded);
switch mode
    case 'parity'
        produceHDivDeviceMexRecord(fullfile(work,'input.h5'),fullfile(work,'mex.h5'),fullfile(work,'common.json'));
    case 'setup-read'
        results=runtests('test_hdiv_gpu_workflow_config');
        assert(~isempty(results) && all([results.Passed]),'Setup read controls did not all pass');
        names=string({results.Name});
        assert(any(endsWith(names,'testWarmAndColdReadsAgree')),'The warm and cold read control did not run');
        writeReport(fullfile(work,'setup-read.json'),struct('accepted',true,'cases',numel(results)));
    case 'device-test'
        fixtures=fullfile(checkout,'tests','fixtures','hdiv');
        receipt=fullfile(fixtures,'hex_pair.receipt.json');
        setup=radia.gpu.prepareHDivDeviceSetup(fullfile(fixtures,'hex_pair.vol'),2, ...
            LabelContract=fullfile(fixtures,'hex_pair.labels.json'),MeshCheckReceipt=receipt, ...
            MeshCheckReceiptSha256=sha256File(receipt), ...
            MeshCheckProvenance=jsondecode(fileread([receipt '.provenance.json'])), ...
            PythonExecutable=string(python));
        % Same soft-iron rows as the parity case.
        bh=[0 0;100 .5;300 1;1000 1.4;5000 1.7;20000 1.9;100000 2.1];
        writeReport(fullfile(work,'device-test.json'),test_prepared_hdiv_device_setup(setup,bh,fullfile(work,'device-test.h5')));
    otherwise
        error('radia:gpu:DeviceBoundary','Unknown mode %s',mode);
end
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
