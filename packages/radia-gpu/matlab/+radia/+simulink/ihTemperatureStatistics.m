function statistics = ihTemperatureStatistics(config, temperature)
%IHTEMPERATURESTATISTICS Sampled extrema and integrated mean, never modal min/max.
temperature = double(temperature(:));
weights = double(config.temperature_cell_weights(:));
if numel(temperature) ~= config.n_temperature || any(~isfinite(temperature))
    error("radia:simulink:IHConfigTemperatureRepresentation", "Invalid temperature state.");
end
if isfield(config,"temperature_constant_coefficients")
    if ~isfield(config,"temperature_evaluation")
        error("radia:simulink:IHConfigTemperatureRepresentation", "FE temperature evaluation operator is required.");
    end
    data = config.temperature_evaluation;
    rows = double(data.rows(:)); columns = double(data.cols(:)); values = double(data.values(:));
    n = double(data.n_samples);
    if ~isscalar(n) || ~isfinite(n) || n < 1 || n ~= fix(n) || ...
            numel(rows) ~= numel(columns) || numel(rows) ~= numel(values) || ...
            any(~isfinite([rows;columns;values])) || ...
            any(rows < 0 | rows >= n | rows ~= fix(rows)) || ...
            any(columns < 0 | columns >= config.n_temperature | columns ~= fix(columns))
        error("radia:simulink:IHConfigTemperatureRepresentation", "Invalid sparse temperature evaluation operator.");
    end
    evaluation = sparse(rows+1,columns+1,values,n,config.n_temperature);
    constant = double(config.temperature_constant_coefficients(:));
    if numel(constant) ~= config.n_temperature || any(~isfinite(constant)) || ...
            max(abs(evaluation*constant-1)) > 1e-9
        error("radia:simulink:IHConfigTemperatureRepresentation", "Temperature evaluation must preserve constants.");
    end
    sampled = evaluation*temperature;
    denominator = dot(weights,constant);
else
    sampled = temperature;
    denominator = sum(weights);
end
if ~isfinite(denominator) || denominator <= 0 || any(~isfinite(sampled))
    error("radia:simulink:IHConfigTemperatureRepresentation", "Invalid heat capacity or evaluated temperature.");
end
statistics = [min(sampled), dot(weights,temperature)/denominator, max(sampled)];
end
