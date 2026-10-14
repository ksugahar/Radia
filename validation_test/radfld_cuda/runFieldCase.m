function runFieldCase(inputPath,algorithm,outputPath)
% Input read is outside the API timing. Geometry setup, all GPU transfers,
% evaluation, result recovery and native destruction are inside it.
kind=h5readatt(inputPath,'/','representation');
if strcmp(kind,'cloud')
    xyz=h5read(inputPath,'/xyz');q=h5read(inputPath,'/strength');
else
    source=struct('version',1,'representation',kind, ...
        'volume',h5read(inputPath,'/volume'),'surface',h5read(inputPath,'/surface'));
    if strcmp(kind,'curved_tet')
        source.gauss_points=h5read(inputPath,'/gauss_points');
        source.gauss_weights=h5read(inputPath,'/gauss_weights');
    else
        source.volume_origins=h5read(inputPath,'/volume_origins');
        source.surface_origins=h5read(inputPath,'/surface_origins');
    end
end
obs=h5read(inputPath,'/observations');
masks=h5read(inputPath,'/masks');signs=h5read(inputPath,'/signs');
angles=h5read(inputPath,'/angles');
leaf=h5read(inputPath,'/leaf');theta=h5read(inputPath,'/theta');
device=gpuDevice();wait(device);
memoryBefore=device.AvailableMemory;
hostBefore=memory;
timer=tic;
if strcmp(kind,'cloud')
    object=radia.gpu.FieldCloud(xyz,q,LeafSize=leaf,Theta=theta, ...
        ImageMasks=masks,ImageSigns=signs,ImageAngles=angles);
else
    object=radia.gpu.FieldElements(source,LeafSize=leaf,Theta=theta, ...
        ImageMasks=masks,ImageSigns=signs,ImageAngles=angles);
end
owner=onCleanup(@()delete(object));
[value,info]=object.field(obs,algorithm);
setupInfo=object.SetupInfo;
delete(object);clear owner;
wait(device);totalSeconds=toc(timer);
memoryAfter=device.AvailableMemory;
hostAfter=memory;
hostBytesBefore=hostBefore.MemUsedMATLAB;hostBytesAfter=hostAfter.MemUsedMATLAB;
save(outputPath,'value','info','setupInfo','totalSeconds','memoryBefore','memoryAfter', ...
    'hostBytesBefore','hostBytesAfter','-v7');
end
