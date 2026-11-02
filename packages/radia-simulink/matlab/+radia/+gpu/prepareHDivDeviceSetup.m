function setup=prepareHDivDeviceSetup(meshPath,order,options)
% Initialization boundary: Radia Python stages the checked mesh's NGSolve
% tables (radia_simulink.hdiv.prepare_hdiv_device_setup); no matrix is built
% on the host and nothing runs per Simulink step. Returns the checked
% radia.hdiv-device-setup.v1 struct consumed by radia_hdiv_device.
arguments
    meshPath (1,:) char {mustBeNonempty}
    order (1,1) double {mustBeMember(order,[0 1 2])}
    options.LabelContract (1,:) char {mustBeNonempty}
    options.MeshCheckReceipt (1,:) char {mustBeNonempty}
    options.MeshCheckReceiptSha256 (1,:) char {mustBeNonempty}
    options.MeshCheckProvenance (1,1) struct
    options.GramACA (1,1) struct=struct()
    options.GramCompression (1,:) char {mustBeMember(options.GramCompression,{'aca','dense'})}='aca'
    % 1: affine TET; 2: NGSolve Curve(2) isoparametric P2 TET (orders 1 and 2).
    options.CurveOrder (1,1) double {mustBeMember(options.CurveOrder,[1 2])}=1
end
keywords=struct('backend','gpu','label_contract',options.LabelContract, ...
    'mesh_check_receipt',options.MeshCheckReceipt,'mesh_check_receipt_sha256',options.MeshCheckReceiptSha256, ...
    'mesh_check_provenance',options.MeshCheckProvenance,'order',int64(order), ...
    'gram_compression',options.GramCompression,'curve_order',int64(options.CurveOrder));
if ~isempty(fieldnames(options.GramACA))
    aca=options.GramACA;aca.leaf_size=int64(aca.leaf_size);aca.max_rank=int64(aca.max_rank);
    keywords.gram_aca=aca;
end
result=radia.internal.callPython("radia_simulink.hdiv","prepare_hdiv_device_setup",{meshPath},Keywords=keywords);
setup=radia.gpu.readHDivDeviceSetup(result.value);
families={'affine-tet','curved-tet'};
assert(setup.order==order && strcmp(setup.family,families{options.CurveOrder}),'radia:gpu:DeviceSetup', ...
    'Device setup order or geometry family differs from the request');
end
