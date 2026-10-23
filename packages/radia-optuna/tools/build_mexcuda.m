function binary = build_mexcuda
%BUILD_MEXCUDA Invoke the canonical recipe from this package's location.
packageRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(packageRoot);
binary = buildOptunaMex;
end
