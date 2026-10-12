function runFieldCase(inputPath,algorithm,outputPath)
% Input read is outside the API timing. Geometry setup, all GPU transfers,
% evaluation, result recovery and native destruction are inside it.
xyz=h5read(inputPath,'/xyz');q=h5read(inputPath,'/strength');
obs=h5read(inputPath,'/observations');
masks=h5read(inputPath,'/masks');signs=h5read(inputPath,'/signs');
angles=h5read(inputPath,'/angles');
leaf=h5read(inputPath,'/leaf');theta=h5read(inputPath,'/theta');
device=gpuDevice();wait(device);
memoryBefore=device.AvailableMemory;
timer=tic;
object=radia.gpu.FieldCloud(xyz,q,LeafSize=leaf,Theta=theta, ...
    ImageMasks=masks,ImageSigns=signs,ImageAngles=angles);
owner=onCleanup(@()delete(object));
[value,info]=object.field(obs,algorithm);
setupInfo=object.SetupInfo;
delete(object);clear owner;
wait(device);totalSeconds=toc(timer);
memoryAfter=device.AvailableMemory;
save(outputPath,'value','info','setupInfo','totalSeconds','memoryBefore','memoryAfter','-v7');
end
