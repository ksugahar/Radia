function value = version()
%VERSION Return the upstream Optuna compatibility version.
%   MATLAB identifiers cannot begin with two underscores, so this function
%   is the language-equivalent surface for both optuna.__version__ and
%   optuna.version.__version__.

value = "5.0.0";
end
