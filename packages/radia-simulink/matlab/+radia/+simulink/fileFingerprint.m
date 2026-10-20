function fingerprint = fileFingerprint(paths)
%FILEFINGERPRINT Content fingerprints for a list of files.
%   fingerprint = fileFingerprint(paths) returns a struct array with one
%   row per input path: path (as given), bytes, and sha256 (lower-case
%   hex of the file content).  Every file must exist -- a missing
%   geometry input is a user error to surface immediately, not a state
%   to fingerprint around.

arguments
    paths (:,1) string
end

fingerprint = struct("path", {}, "bytes", {}, "sha256", {});
for index = 1:numel(paths)
    path = paths(index);
    if ~isfile(path)
        error("radia:simulink:FileFingerprintMissing", ...
            "File does not exist: %s", path);
    end
    identifier = fopen(path, "r");
    if identifier < 0
        error("radia:simulink:FileFingerprintUnreadable", ...
            "File cannot be read: %s", path);
    end
    closer = onCleanup(@() fclose(identifier));
    data = fread(identifier, Inf, "*uint8");
    clear closer
    digest = java.security.MessageDigest.getInstance("SHA-256");
    digest.update(data);
    hashBytes = typecast(digest.digest(), "uint8");
    fingerprint(index, 1) = struct( ...
        "path", char(path), ...
        "bytes", numel(data), ...
        "sha256", lower(reshape(dec2hex(hashBytes, 2).', 1, [])));
end
end
