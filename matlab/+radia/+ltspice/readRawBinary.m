function data = readRawBinary(rawFile)
%READRAWBINARY Read LTspice real transient binary RAW files.
arguments, rawFile (1,1) string {mustBeFile}, end
bytes=readBytes(rawFile); marker=unicode2native("Binary:"+newline,'UTF-16LE'); at=findPattern(bytes,marker);
if isempty(at), error("radia:ltspice:RawFormat","Binary marker not found."); end
headerEnd=at+numel(marker)-1; header=native2unicode(bytes(1:headerEnd),'UTF-16LE'); lines=splitlines(string(header));
isComplex=any(contains(lines,"Flags:") & contains(lines,"complex"));
nvar=headerInt(lines,"No. Variables:"); npoint=headerInt(lines,"No. Points:");
vline=find(strtrim(lines)=="Variables:",1); names=strings(1,nvar); types=strings(1,nvar);
for k=1:nvar
 tok=regexp(char(lines(vline+k)),'^\s*\d+\s+(\S+)\s+(\S+)\s*$','tokens','once'); names(k)=string(tok{1}); types(k)=string(tok{2});
end
payload=bytes(headerEnd+1:end);
if isComplex, stride=16*nvar; else, stride=8+4*(nvar-1); end
if numel(payload)<stride*npoint, error("radia:ltspice:RawTruncated","Binary RAW payload is truncated."); end
if isComplex, values=complex(zeros(npoint,nvar)); else, values=zeros(npoint,nvar); end
for p=1:npoint
 base=(p-1)*stride;
 if isComplex
  for k=1:nvar
   q=base+16*(k-1); pair=typecast(uint8(payload(q+(1:16))),'double'); values(p,k)=complex(pair(1),pair(2));
  end
 else
  values(p,1)=typecast(uint8(payload(base+(1:8))),'double');
  for k=2:nvar, q=base+8+4*(k-2); values(p,k)=double(typecast(uint8(payload(q+(1:4))),'single')); end
 end
end
data=struct("schema","radia.ltspice.raw.binary.v1","path",rawFile,"names",names,"types",types,"values",values,"is_complex",isComplex);
end
function bytes=readBytes(path), f=fopen(path,'rb'); c=onCleanup(@()fclose(f)); bytes=fread(f,inf,'*uint8').'; clear c, end
function at=findPattern(bytes,pattern), at=strfind(bytes,uint8(pattern)); if ~isempty(at),at=at(1);end, end
function value=headerInt(lines,label), row=find(startsWith(strtrim(lines),label),1); value=sscanf(char(extractAfter(strtrim(lines(row)),label)),'%d'); end
