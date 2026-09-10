function data = readRawBinary(rawFile)
%READRAWBINARY Read LTspice binary RAW files without guessing storage width.
arguments, rawFile (1,1) string {mustBeFile}, end
bytes=readBytes(rawFile); marker=unicode2native("Binary:"+newline,'UTF-16LE'); at=findPattern(bytes,marker);
if isempty(at), error("radia:ltspice:RawFormat","Binary marker not found."); end
headerEnd=at+numel(marker)-1; header=native2unicode(bytes(1:headerEnd),'UTF-16LE'); lines=splitlines(string(header));
flagRow=find(startsWith(lower(strtrim(lines)),"flags:"),1);
if isempty(flagRow),error("radia:ltspice:RawFlags","Binary RAW header has no Flags entry.");end
flags=lower(split(strtrim(extractAfter(strtrim(lines(flagRow)),":"))));flags(flags=="")=[];
known=["real","complex","double","forward","stepped","log","fastaccess"];
unknown=setdiff(flags,known);
if ~isempty(unknown),error("radia:ltspice:RawFlags","Unsupported binary RAW flag(s): %s.",join(unknown,", "));end
isComplex=any(flags=="complex");isReal=any(flags=="real");allDouble=any(flags=="double");
if isComplex==isReal,error("radia:ltspice:RawFlags","Binary RAW must declare exactly one of real or complex.");end
nvar=headerInt(lines,"No. Variables:"); npoint=headerInt(lines,"No. Points:");
vline=find(strtrim(lines)=="Variables:",1); names=strings(1,nvar); types=strings(1,nvar);
for k=1:nvar
 tok=regexp(char(lines(vline+k)),'^\s*\d+\s+(\S+)\s+(\S+)\s*$','tokens','once'); names(k)=string(tok{1}); types(k)=string(tok{2});
end
payload=bytes(headerEnd+1:end);
if isComplex,stride=16*nvar;elseif allDouble,stride=8*nvar;else,stride=8+4*(nvar-1);end
if numel(payload)<stride*npoint, error("radia:ltspice:RawTruncated","Binary RAW payload is truncated."); end
payload=payload(1:stride*npoint);
if isComplex
 doubles=reshape(typecast(uint8(payload),'double'),2*nvar,npoint).';
 values=complex(doubles(:,1:2:end),doubles(:,2:2:end));
elseif allDouble
 values=reshape(typecast(uint8(payload),'double'),nvar,npoint).';
else
 records=reshape(payload,stride,npoint);values=zeros(npoint,nvar);
 values(:,1)=typecast(reshape(records(1:8,:),1,[]),'double').';
 if nvar>1,values(:,2:end)=double(reshape(typecast(reshape(records(9:end,:),1,[]),'single'),nvar-1,npoint).');end
end
rawProperties=parseHeaderProperties(lines);
data=struct("schema","radia.ltspice.raw.binary.v1","path",rawFile,"names",names,"types",types,"values",values,"is_complex",isComplex,"raw_properties",rawProperties,"flags",flags);
end
function bytes=readBytes(path), f=fopen(path,'rb'); c=onCleanup(@()fclose(f)); bytes=fread(f,inf,'*uint8').'; clear c, end
function at=findPattern(bytes,pattern), at=strfind(bytes,uint8(pattern)); if ~isempty(at),at=at(1);end, end
function value=headerInt(lines,label), row=find(startsWith(strtrim(lines),label),1); value=sscanf(char(extractAfter(strtrim(lines(row)),label)),'%d'); end
function properties=parseHeaderProperties(lines)
properties=struct();
for k=1:numel(lines)
 line=strtrim(lines(k));colon=strfind(line,":");if isempty(colon),continue,end
 label=strtrim(extractBefore(line,colon(1)));if label=="Binary",continue,end
 key=matlab.lang.makeValidName(char(label));properties.(key)=strtrim(extractAfter(line,colon(1)));
end
end
