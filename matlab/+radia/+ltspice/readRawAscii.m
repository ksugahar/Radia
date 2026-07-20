function data = readRawAscii(rawFile)
%READRAWASCII Read an LTspice ASCII .raw waveform file.
arguments
    rawFile (1,1) string {mustBeFile}
end

text = fileread(rawFile);
lines = splitlines(string(text));
nVariables = readHeaderInteger(lines, "No. Variables:");
nPoints = readHeaderInteger(lines, "No. Points:");
variablesLine = find(strtrim(lines) == "Variables:", 1);
valuesLine = find(strtrim(lines) == "Values:", 1);
if isempty(variablesLine) || isempty(valuesLine) || valuesLine <= variablesLine
    error("radia:ltspice:RawFormat", ...
        "The file is not an LTspice ASCII RAW waveform: %s", rawFile);
end

names = strings(1, nVariables);
types = strings(1, nVariables);
for k = 1:nVariables
    tokens = regexp(char(lines(variablesLine + k)), ...
        '^\s*\d+\s+(\S+)\s+(\S+)\s*$', 'tokens', 'once');
    if isempty(tokens)
        error("radia:ltspice:RawFormat", "Invalid variable row %d.", k);
    end
    names(k) = string(tokens{1});
    types(k) = string(tokens{2});
end

isComplex=any(contains(lines,"Flags:") & contains(lines,"complex"));
if isComplex, values=complex(zeros(nPoints,nVariables)); else, values=zeros(nPoints,nVariables); end
cursor = valuesLine + 1;
for point = 1:nPoints
    while cursor <= numel(lines) && strlength(strtrim(lines(cursor))) == 0
        cursor = cursor + 1;
    end
    firstText=strtrim(lines(cursor)); firstParts=regexp(char(firstText),'^\s*\d+\s+(.+)$','tokens','once');
    if isempty(firstParts)
        error("radia:ltspice:RawFormat", "Invalid point row %d.", point);
    end
    values(point, 1) = parseValue(firstParts{1});
    cursor = cursor + 1;
    for variable = 2:nVariables
        entry = parseValue(strtrim(lines(cursor)));
        if ~isfinite(entry)
            error("radia:ltspice:RawFormat", ...
                "Missing variable %d at point %d.", variable, point);
        end
        values(point, variable) = entry;
        cursor = cursor + 1;
    end
end

function value=parseValue(text)
text=strtrim(string(text)); pair=regexp(char(text),'^\(?\s*([+\-0-9.eE]+)\s*,\s*([+\-0-9.eE]+)\s*\)?$','tokens','once');
if ~isempty(pair), value=complex(str2double(pair{1}),str2double(pair{2})); return, end
tokens=regexp(char(text),'[+\-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+\-]?\d+)?','match');
if isempty(tokens), value=NaN; else, value=str2double(tokens{end}); end
end

validNames = matlab.lang.makeUniqueStrings( ...
    matlab.lang.makeValidName(cellstr(names), "ReplacementStyle", "hex"));
signals = struct();
for k = 1:nVariables
    signals.(validNames{k}) = values(:, k);
end
data = struct( ...
    "schema", "radia.ltspice.raw.v1", ...
    "path", rawFile, ...
    "names", names, ...
    "types", types, ...
    "values", values, ...
    "signals", signals);
end

function value = readHeaderInteger(lines, label)
row = find(startsWith(strtrim(lines), label), 1);
if isempty(row)
    error("radia:ltspice:RawFormat", "Missing RAW header: %s", label);
end
value = sscanf(char(extractAfter(strtrim(lines(row)), label)), '%d');
if isempty(value) || value < 1
    error("radia:ltspice:RawFormat", "Invalid RAW header: %s", label);
end
end
