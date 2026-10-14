classdef LiveMonitor < handle
    %LIVEMONITOR Live objective, best-value, duration, and parameter display.

    properties (SetAccess=private)
        Figure
        ObjectiveAxes
        DurationAxes
        ParameterAxes
        StatusText
        UpdateCount (1,1) double = 0
    end

    methods
        function obj = LiveMonitor(options)
            arguments
                options.Name (1,1) string = "Radia MATLAB Optuna"
                options.Visible (1,1) logical = true
            end
            visibility = "off";
            if options.Visible, visibility = "on"; end
            obj.Figure = figure("Name", options.Name, "NumberTitle", "off", ...
                "Visible", visibility, "Color", "white");
            layout = tiledlayout(obj.Figure, 2, 2, "TileSpacing", "compact");
            obj.ObjectiveAxes = nexttile(layout, 1, [1 2]);
            obj.DurationAxes = nexttile(layout, 3);
            obj.ParameterAxes = nexttile(layout, 4);
            obj.StatusText = annotation(obj.Figure, "textbox", ...
                [0.01 0.955 0.98 0.04], "EdgeColor", "none", ...
                "HorizontalAlignment", "right", "String", "waiting");
        end

        function update(obj, snapshot)
            obj.UpdateCount = obj.UpdateCount + 1;
            trials = snapshot.trial_table;
            complete = trials.State == "COMPLETE" & isfinite(trials.Value);
            cla(obj.ObjectiveAxes);
            pareto = snapshot.pareto_front;
            isMultiObjective = height(pareto) > 0 && numel(pareto.Values{1}) > 1;
            if isMultiObjective
                matrix = vertcat(pareto.Values{:});
                scatter(obj.ObjectiveAxes, matrix(:,1), matrix(:,2), 48, ...
                    pareto.TrialNumber, "filled");
                xlabel(obj.ObjectiveAxes, "Objective 1");
                ylabel(obj.ObjectiveAxes, "Objective 2");
                title(obj.ObjectiveAxes, "Live Pareto front");
                colorbar(obj.ObjectiveAxes);
            elseif any(complete)
                numbers = trials.TrialNumber(complete);
                values = trials.Value(complete);
                plot(obj.ObjectiveAxes, numbers, values, "o-", ...
                    "DisplayName", "objective");
                hold(obj.ObjectiveAxes, "on");
                if height(snapshot.best_trial) > 0
                    yline(obj.ObjectiveAxes, snapshot.best_trial.Value(1), "--", ...
                        "DisplayName", "best");
                end
                hold(obj.ObjectiveAxes, "off");
                legend(obj.ObjectiveAxes, "Location", "best");
                xlabel(obj.ObjectiveAxes, "Trial");
                ylabel(obj.ObjectiveAxes, "Objective");
            else
                xlabel(obj.ObjectiveAxes, "Trial");
                ylabel(obj.ObjectiveAxes, "Objective");
            end
            grid(obj.ObjectiveAxes, "on");

            cla(obj.DurationAxes);
            finished = isfinite(trials.Duration_s);
            if any(finished)
                bar(obj.DurationAxes, trials.TrialNumber(finished), ...
                    trials.Duration_s(finished));
            end
            xlabel(obj.DurationAxes, "Trial");
            ylabel(obj.DurationAxes, "Duration (s)");
            grid(obj.DurationAxes, "on");

            cla(obj.ParameterAxes);
            params = snapshot.param_table;
            numeric = isfinite(params.ValueNumeric);
            if any(numeric)
                names = unique(params.Name(numeric), "stable");
                hold(obj.ParameterAxes, "on");
                for k = 1:numel(names)
                    rows = numeric & params.Name == names(k);
                    plot(obj.ParameterAxes, params.TrialNumber(rows), ...
                        params.ValueNumeric(rows), ".-", "DisplayName", names(k));
                end
                hold(obj.ParameterAxes, "off");
                legend(obj.ParameterAxes, "Location", "best");
            end
            xlabel(obj.ParameterAxes, "Trial");
            ylabel(obj.ParameterAxes, "Parameter value");
            grid(obj.ParameterAxes, "on");
            obj.StatusText.String = sprintf("%s | trial %d | %s", ...
                snapshot.event, snapshot.trial_number, snapshot.trial_state);
            drawnow limitrate;
        end

        function delete(obj)
            if ~isempty(obj.Figure) && isgraphics(obj.Figure)
                delete(obj.Figure);
            end
        end
    end
end
