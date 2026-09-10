function data = readRawBinary(rawFile)
%READRAWBINARY Read explicitly supported little-endian LTspice RAW layouts.
arguments, rawFile (1,1) string {mustBeFile}, end
f=fopen(rawFile,'rb');
if f<0,error("radia:ltspice:Read","Cannot open %s",rawFile);end
c=onCleanup(@()fclose(f)); bytes=fread(f,inf,'*uint8').'; clear c
marker=unicode2native("Binary:"+newline,'UTF-16LE'); at=strfind(bytes,marker);
if isempty(at)
    marker=unicode2native("Binary:"+char(13)+newline,'UTF-16LE'); at=strfind(bytes,marker);
end
if isempty(at),error("radia:ltspice:RawFormat","Binary marker not found.");end
headerEnd=at(1)+numel(marker)-1;
lines=splitlines(string(native2unicode(bytes(1:headerEnd),'UTF-16LE')));
flags=split(lower(headerValue(lines,"Flags:")));
supported=["real","complex","double","forward","stepped","log"];
if any(~ismember(flags,supported)) || sum(ismember(flags,["real","complex"]))~=1
    error("radia:ltspice:RawFlags","Unsupported or conflicting RAW flags: %s",join(flags," "));
end
isComplex=any(flags=="complex"); allDouble=any(flags=="double");
nvar=headerInt(lines,"No. Variables:"); npoint=headerInt(lines,"No. Points:");
vline=find(strtrim(lines)=="Variables:");
if numel(vline)~=1 || vline+nvar>numel(lines)
    error("radia:ltspice:RawFormat","Invalid variable table.");
end
names=strings(1,nvar); types=strings(1,nvar);
for k=1:nvar
    tok=regexp(char(lines(vline+k)),'^\s*(\d+)\s+(\S+)\s+(\S+)\s*$','tokens','once');
    if isempty(tok) || str2double(tok{1})~=k-1
        error("radia:ltspice:RawFormat","Invalid variable row %d.",k);
    end
    names(k)=string(tok{2}); types(k)=string(tok{3});
end
properties=struct();
for k=1:vline-1
    tok=regexp(char(lines(k)),'^\s*([^:]+):\s*(.*)$','tokens','once');
    if ~isempty(tok),properties.(matlab.lang.makeValidName(tok{1}))=string(strtrim(tok{2}));end
end
if isComplex,stride=16*nvar;elseif allDouble,stride=8*nvar;else,stride=8+4*(nvar-1);end
payload=bytes(headerEnd+1:end);
if numel(payload)<stride*npoint,error("radia:ltspice:RawTruncated","Binary RAW payload is truncated.");end
if numel(payload)~=stride*npoint
    error("radia:ltspice:RawFormat","Payload size does not match the declared layout.");
end
if isComplex
    parts=reshape(decode(payload,'double'),2*nvar,npoint);
    values=complex(parts(1:2:end,:),parts(2:2:end,:)).';
elseif allDouble
    values=reshape(decode(payload,'double'),nvar,npoint).';
else
    records=reshape(payload,stride,npoint);
    timeBytes=records(1:8,:); signalBytes=records(9:end,:);
    values=[decode(timeBytes(:),'double').', ...
        reshape(double(decode(signalBytes(:),'single')),nvar-1,npoint).'];
end
data=struct("schema","radia.ltspice.raw.binary.v1","path",rawFile, ...
    "names",names,"types",types,"values",values,"is_complex",isComplex, ...
    "raw_properties",properties);
end

function value=headerValue(lines,label)
rows=find(startsWith(strtrim(lines),label));
if numel(rows)~=1,error("radia:ltspice:RawFormat","Expected one %s header.",label);end
value=strtrim(extractAfter(strtrim(lines(rows)),label));
end

function value=headerInt(lines,label)
value=str2double(headerValue(lines,label));
if ~isscalar(value)||~isfinite(value)||value<1||fix(value)~=value
    error("radia:ltspice:RawFormat","Invalid %s count.",label);
end
end

function values=decode(bytes,precision)
values=typecast(uint8(bytes(:).'),precision);
[~,~,endian]=computer;
if endian=='B',values=swapbytes(values);end
end
