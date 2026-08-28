classdef UpstreamNumerics
    %UPSTREAMNUMERICS Numeric primitives that must match Optuna 4.9 exactly.
    %   MATLAB and NumPy/Python disagree on three primitives that appear in
    %   every Optuna untransform path, so every sampler routes through this
    %   class instead of the MATLAB built-ins:
    %
    %     round        MATLAB rounds halves away from zero; NumPy rounds
    %                  halves to even (np.round).
    %     upper bound  Optuna returns np.nextafter(high, -inf) for a
    %                  non-single FloatDistribution, i.e. the interval is
    %                  half open, not clamped to high.
    %     step high    Optuna adjusts a non-divisible high with
    %                  decimal.Decimal(str(x)) arithmetic, not with binary
    %                  floating point.

    properties (Constant, Access=private)
        MaxSignificantDigits = 18
        MaxDecimalExponent = 22
    end

    methods (Static)
        function value = roundTiesToEven(value)
            %ROUNDTIESTOEVEN Elementwise numpy.round (banker's rounding).
            value = double(value);
            lower = floor(value);
            ties = (value - lower) == 0.5;
            value = round(value);
            if any(ties(:))
                even = lower(ties);
                odd = mod(even, 2) ~= 0;
                even(odd) = even(odd) + 1;
                value(ties) = even;
            end
        end

        function value = nextDown(value)
            %NEXTDOWN Elementwise np.nextafter(value, -inf) for finite input.
            value = double(value);
            finite = isfinite(value);
            if ~any(finite(:))
                return
            end
            selected = value(finite);
            bits = typecast(selected(:), "uint64");
            positive = selected(:) > 0;
            negative = selected(:) < 0;
            bits(positive) = bits(positive) - uint64(1);
            bits(negative) = bits(negative) + uint64(1);
            % +0 and -0 both step down to the largest negative subnormal.
            bits(~positive & ~negative) = ...
                bitor(bitshift(uint64(1), 63), uint64(1));
            selected(:) = typecast(bits, "double");
            value(finite) = selected;
        end

        function [high, adjusted] = adjustDiscreteUniformHigh(low, high, step)
            %ADJUSTDISCRETEUNIFORMHIGH Port of Optuna's decimal high fix-up.
            %   Optuna's _adjust_discrete_uniform_high compares
            %   Decimal(str(high)) - Decimal(str(low)) against
            %   Decimal(str(step)), so [0, 0.7] with Step=0.1 is divisible
            %   and high is preserved bit for bit. A binary floating-point
            %   test reports the same range as non-divisible and replaces
            %   high with 0.7000000000000001, i.e. above the requested
            %   bound, and emits a warning upstream never emits.
            low = double(low);
            high = double(high);
            step = double(step);
            adjusted = false;
            if ~isfinite(step) || step <= 0 || ~isfinite(low) || ~isfinite(high)
                return
            end
            helper = radia.optuna.internal.UpstreamNumerics;
            [lowMantissa, lowExponent] = helper.decimalParts(low);
            [highMantissa, highExponent] = helper.decimalParts(high);
            [stepMantissa, stepExponent] = helper.decimalParts(step);
            exponent = min([lowExponent, highExponent, stepExponent]);
            scaledLow = helper.rescale(lowMantissa, lowExponent - exponent);
            scaledHigh = helper.rescale(highMantissa, highExponent - exponent);
            scaledStep = helper.rescale(stepMantissa, stepExponent - exponent);
            range = scaledHigh - scaledLow;
            if mod(range, scaledStep) == 0
                return
            end
            scaledHigh = idivide(range, scaledStep, "floor") * scaledStep + ...
                scaledLow;
            high = helper.unscale(scaledHigh, exponent);
            adjusted = true;
        end

        function text = shortestDecimal(value)
            %SHORTESTDECIMAL Python repr(float): shortest round-trip decimal.
            value = double(value);
            text = string(sprintf("%.17g", value));
            for precision = 1:16
                candidate = string(sprintf("%.*g", precision, value));
                if str2double(candidate) == value
                    text = candidate;
                    return
                end
            end
        end
    end

    methods (Static, Access=private)
        function [mantissa, exponent] = decimalParts(value)
            % Exact (mantissa, exponent) of the shortest decimal that
            % round-trips to VALUE, mirroring Decimal(str(value)). Only the
            % numeric value matters downstream, so an exponent-form
            % representation such as "1e+02" is as good as "100.0".
            helper = radia.optuna.internal.UpstreamNumerics;
            text = char(helper.shortestDecimal(value));
            exponent = 0;
            marker = find(text == 'e' | text == 'E', 1);
            if ~isempty(marker)
                exponent = str2double(text(marker+1:end));
                text = text(1:marker-1);
            end
            negative = ~isempty(text) && text(1) == '-';
            if negative || (~isempty(text) && text(1) == '+')
                text = text(2:end);
            end
            point = find(text == '.', 1);
            if ~isempty(point)
                exponent = exponent - (numel(text) - point);
                text(point) = [];
            end
            digits = double(text) - 48;
            if isempty(digits) || any(digits < 0 | digits > 9)
                error("radia:optuna:StepPrecision", ...
                    "Cannot read a decimal representation of %.17g.", value);
            end
            trailing = find(digits ~= 0, 1, "last");
            if isempty(trailing)
                mantissa = int64(0);
                exponent = 0;
                return
            end
            % Normalize so the digit-count guard measures real precision.
            exponent = exponent + (numel(digits) - trailing);
            digits = digits(find(digits ~= 0, 1):trailing);
            if numel(digits) > helper.MaxSignificantDigits
                error("radia:optuna:StepPrecision", ...
                    "%.17g needs more than %d significant decimal digits.", ...
                    value, helper.MaxSignificantDigits);
            end
            mantissa = int64(0);
            for index = 1:numel(digits)
                mantissa = mantissa * int64(10) + int64(digits(index));
            end
            if negative
                mantissa = -mantissa;
            end
        end

        function value = rescale(mantissa, shift)
            % mantissa * 10^shift, refusing to saturate silently.
            helper = radia.optuna.internal.UpstreamNumerics;
            if mantissa == 0
                value = int64(0);
                return
            end
            if shift < 0 || shift > helper.MaxSignificantDigits
                error("radia:optuna:StepPrecision", ...
                    "The distribution decimal grid exceeds int64 precision.");
            end
            multiplier = int64(10)^int64(shift);
            if mantissa ~= 0 && ...
                    abs(double(mantissa)) * double(multiplier) > ...
                    double(intmax("int64")) / 4
                error("radia:optuna:StepPrecision", ...
                    "The distribution decimal grid exceeds int64 precision.");
            end
            value = mantissa * multiplier;
        end

        function value = unscale(scaled, exponent)
            % double(Decimal(scaled) * 10^exponent), correctly rounded.
            helper = radia.optuna.internal.UpstreamNumerics;
            if abs(scaled) > int64(flintmax) || ...
                    abs(exponent) > helper.MaxDecimalExponent
                error("radia:optuna:StepPrecision", ...
                    "Cannot represent the aligned high bound exactly; " + ...
                    "rescale the distribution so its decimal grid needs " + ...
                    "fewer digits.");
            end
            value = double(scaled);
            if exponent < 0
                value = value / 10^(-exponent);
            elseif exponent > 0
                value = value * 10^exponent;
            end
        end
    end
end
