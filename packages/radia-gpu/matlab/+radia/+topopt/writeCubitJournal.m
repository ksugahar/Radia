function info=writeCubitJournal(path,elementIds,density,options)
%WRITECUBITJOURNAL Write Cubit material block assignment for optimized cells.
arguments
 path (1,1) string; elementIds (:,1) double {mustBeInteger,mustBePositive}; density (:,1) double {mustBeBetween(density,0,1)}
 options.Threshold (1,1) double {mustBeBetween(options.Threshold,0,1)}=0.5
 options.SolidBlock (1,1) double {mustBeInteger,mustBePositive}=1001
 options.VoidBlock (1,1) double {mustBeInteger,mustBePositive}=1002
end
if numel(elementIds)~=numel(density)||numel(unique(elementIds))~=numel(elementIds), error("radia:topopt:Elements","Element IDs must be unique and match density."); end
solid=elementIds(density>=options.Threshold); void=elementIds(density<options.Threshold);
folder=fileparts(path); if strlength(folder)>0&&~isfolder(folder), mkdir(folder); end
file=fopen(path,'w'); if file<0, error("radia:topopt:Write","Cannot write %s",path); end
cleanup=onCleanup(@()fclose(file)); fprintf(file,"# Radia VIM linearized topology material assignment\nset echo off\n");
writeGroup(file,"radia_topopt_solid",solid,options.SolidBlock); writeGroup(file,"radia_topopt_void",void,options.VoidBlock); clear cleanup
info=struct("path",path,"solid_count",numel(solid),"void_count",numel(void),"threshold",options.Threshold);
end
function writeGroup(file,name,ids,block)
fprintf(file,"group '%s'",name); if ~isempty(ids), fprintf(file," add hex%s",sprintf(" %d",ids)); end
fprintf(file,"\nblock %d hex in group '%s'\n",block,name);
end
