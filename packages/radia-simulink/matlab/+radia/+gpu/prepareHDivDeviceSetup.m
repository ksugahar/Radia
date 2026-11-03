function setup=prepareHDivDeviceSetup(meshPath,order,options)
% Initialization boundary: a separate Radia Python process stages the checked
% mesh's NGSolve tables (python -m radia_simulink.hdiv) and returns them in HDF5.
% Staging loads the pip CUDA runtime, which must not share a process with
% MATLAB's own. No matrix is built on the host and nothing runs per Simulink
% step. Returns the checked radia.hdiv-device-setup.v1 struct consumed by
% radia_hdiv_device.
arguments
    meshPath (1,:) char {mustBeNonempty}
    order (1,1) double {mustBeMember(order,[0 1 2])}
    options.LabelContract (1,:) char {mustBeNonempty}
    options.MeshCheckReceipt (1,:) char {mustBeNonempty}
    options.MeshCheckReceiptSha256 (1,:) char {mustBeNonempty}
    options.MeshCheckProvenance (1,1) struct
    options.GramACA (1,1) struct=struct()
    options.GramCompression (1,:) char {mustBeMember(options.GramCompression,{'aca','dense'})}='aca'
    % 1: linear geometry; 2: NGSolve Curve(2) (orders 1 and 2). A pure-TET mesh
    % is affine-tet or curved-tet by curve order; a pure-HEX mesh is hex at either.
    options.CurveOrder (1,1) double {mustBeMember(options.CurveOrder,[1 2])}=1
    options.PythonExecutable (1,1) string=string(pyenv().Executable)
end
assert(strlength(options.PythonExecutable)>0 && ~contains(options.PythonExecutable,'"'),'radia:gpu:DeviceSetup', ...
    'Select the Python that has radia[gpu] through PythonExecutable or pyenv');
request=struct('mesh_path',meshPath,'backend','gpu','label_contract',options.LabelContract, ...
    'mesh_check_receipt',options.MeshCheckReceipt,'mesh_check_receipt_sha256',options.MeshCheckReceiptSha256, ...
    'mesh_check_provenance',options.MeshCheckProvenance,'order',order, ...
    'gram_compression',options.GramCompression,'curve_order',options.CurveOrder);
if ~isempty(fieldnames(options.GramACA))
    request.gram_aca=options.GramACA;
end
work=tempname;mkdir(work);workOwner=onCleanup(@()rmdir(work,'s'));
requestPath=fullfile(work,'request.json');outputPath=fullfile(work,'setup.h5');
stream=fopen(requestPath,'w','n','UTF-8');assert(stream>0,'radia:gpu:DeviceSetup','Cannot write %s',requestPath);
fwrite(stream,unicode2native(jsonencode(request),'UTF-8'));fclose(stream);
command=sprintf('"%s" -m radia_simulink.hdiv "%s" "%s"',options.PythonExecutable,requestPath,outputPath);
[status,output]=radia.internal.runPythonProcess(string(command));
if status~=0
    error('radia:gpu:DeviceSetup','HDiv device setup process failed with status %d:\n%s',status,output);
end
setup=radia.gpu.loadHDivDeviceSetup(outputPath);
families={'affine-tet','curved-tet'};
assert(setup.order==order && any(strcmp(setup.family,{families{options.CurveOrder},'hex'})),'radia:gpu:DeviceSetup', ...
    'Device setup order or geometry family differs from the request');
end
