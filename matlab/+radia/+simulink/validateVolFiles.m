function reports = validateVolFiles(volFiles, options)
%VALIDATEVOLFILES Run check-vol before native solver initialization.
arguments
    volFiles (:,1) string
    options.ReportDirectory (1,1) string = ""
    options.Contract (1,1) string = ""
end
if isempty(volFiles), reports = strings(0,1); return; end
if strlength(options.ReportDirectory)==0
    options.ReportDirectory = fullfile("C:\temp","radia_vol_checks");
end
if ~isfolder(options.ReportDirectory), mkdir(options.ReportDirectory); end
reports = strings(size(volFiles));
for k=1:numel(volFiles)
    if ~isfile(volFiles(k)), error("radia:simulink:VolMissing", "Missing .vol file: %s",volFiles(k)); end
    [~,stem] = fileparts(volFiles(k));
    report = fullfile(options.ReportDirectory,stem+".vol-check.json");
    command = "check-vol " + quote(volFiles(k)) + " --format json --report-json " + quote(report);
    if strlength(options.Contract)>0
        command = command + " --contract " + quote(options.Contract) + " --strict-labels";
    end
    [status,output] = system(command);
    if status ~= 0
        error("radia:simulink:VolCheckFailed", "check-vol failed for %s:\n%s",volFiles(k),output);
    end
    reports(k) = report;
end
end

function value = quote(path)
value = '"' + replace(string(path),'"','""') + '"';
end
