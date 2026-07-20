function blkStruct=slblocks
%SLBLOCKS Register the packaged Radia block library in Library Browser.
Browser.Library="radia_simulink_library";
Browser.Name="Radia";
Browser.IsFlat=0;
blkStruct.Browser=Browser;
end
