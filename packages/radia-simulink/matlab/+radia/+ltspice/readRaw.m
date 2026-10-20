function data = readRaw(rawFile)
%READRAW Auto-detect ASCII or LTspice binary RAW data.
arguments, rawFile (1,1) string {mustBeFile}, end
f=fopen(rawFile,'r'); c=onCleanup(@()fclose(f)); prefix=fread(f,256,'*uint8'); clear c
if any(prefix==0), data=radia.ltspice.readRawBinary(rawFile); else, data=radia.ltspice.readRawAscii(rawFile); end
data=addSteps(data);
data=normalizeContract(data);
end
function data=addSteps(data)
t=real(data.values(:,1)); starts=[1;find(diff(t)<0)+1]; stops=[starts(2:end)-1;numel(t)];
data.step_ranges=[starts,stops]; data.step_count=numel(starts);
data.step_index=zeros(numel(t),1);
for k=1:numel(starts), data.step_index(starts(k):stops(k))=k; end
end
function data=normalizeContract(data)
data.contract_schema="radia.ltspice.raw.v2";
data.is_complex=~isreal(data.values);
data.axis_name=data.names(1);
if data.axis_name=="frequency", data.analysis="ac";
elseif data.axis_name=="time", data.analysis="transient";
else, data.analysis="other";
end
valid=matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(cellstr(data.names),"ReplacementStyle","hex"));
signals=struct();
for k=1:numel(valid), signals.(valid{k})=data.values(:,k); end
data.signals=signals;
if ~isfield(data,'raw_properties'),data.raw_properties=struct();end
if ~isfield(data,'aliases'),data.aliases=struct();end
if ~isfield(data,'alias_names'),data.alias_names=string(fieldnames(data.aliases));end
if ~isfield(data,'alias_formulas')
    keys=fieldnames(data.aliases);
    data.alias_formulas=strings(numel(keys),1);
    for k=1:numel(keys)
        data.alias_formulas(k)=string(data.aliases.(keys{k}));
    end
end
if ~isfield(data,'backannotations'),data.backannotations=strings(0,1);end
if isfield(data.raw_properties,'Plotname')
    data.plot_name=string(data.raw_properties.Plotname);
elseif data.analysis=="ac"
    data.plot_name="AC Analysis";
elseif data.analysis=="transient"
    data.plot_name="Transient Analysis";
else
    data.plot_name="";
end
end
