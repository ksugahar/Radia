function mexPath = radia_motor_rom_sfun_build(repoRoot, outputDir)
% Build the Radia motor-ROM C-MEX S-function for Simulink.
%
% The numerical C ABI source is compiled into the MEX file directly.  This
% avoids a platform-specific import-library dependency and keeps the MEX
% adapter tied to the same source as radia_motor_rom.dll.

if nargin < 1 || isempty(repoRoot)
    here = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(fileparts(fileparts(here)));
end
if nargin < 2 || isempty(outputDir)
    outputDir = fullfile(repoRoot, 'build', 'simulink');
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

sfunSource = fullfile(repoRoot, 'src', 'radia', 'simulink', ...
                      'radia_motor_rom_sfun.cpp');
coreSource = fullfile(repoRoot, 'src', 'core', 'rad_motor_rom_c.cpp');
includeDir = fullfile(repoRoot, 'src', 'core');
assert(isfile(sfunSource), 'Missing S-function source: %s', sfunSource);
assert(isfile(coreSource), 'Missing Radia C ABI source: %s', coreSource);

mex('-R2018a', '-O', ['-I', includeDir], '-outdir', outputDir, ...
    sfunSource, coreSource);
mexPath = fullfile(outputDir, ['radia_motor_rom_sfun.', mexext]);
assert(isfile(mexPath), 'MEX build did not produce %s', mexPath);
fprintf('Built %s\n', mexPath);
end
