classdef (Abstract) ExceptionValue < handle
    %EXCEPTIONVALUE MATLAB value object for Python-style Optuna exceptions.

    properties (SetAccess=protected)
        args (1,:) cell = cell(1,0)
        identifier (1,1) string = "radia:optuna:OptunaError"
    end

    properties (SetAccess=private)
        notes (1,:) string = strings(1,0)
    end

    properties (Dependent, SetAccess=private)
        message
    end

    methods
        function obj=ExceptionValue(identifier,varargin)
            obj.identifier=string(identifier);
            obj.args=reshape(varargin,1,[]);
        end

        function value=get.message(obj)
            value=radia.optuna.ExceptionValue.formatArgs(obj.args);
        end

        function add_note(obj,note)
            if ~((ischar(note) && (isrow(note) || isempty(note))) || ...
                    (isstring(note) && isscalar(note)))
                error("radia:optuna:ExceptionNoteType", ...
                    "Exception notes must be scalar text.");
            end
            obj.notes(end+1)=string(note);
        end

        function result=with_traceback(obj,traceback) %#ok<INUSD>
            % MATLAB captures the stack when throw is called. Returning the
            % same handle is the language-equivalent identity contract.
            result=obj;
        end

        function throw(obj)
            exception=MException(char(obj.identifier),"%s",obj.message);
            for index=1:numel(obj.notes)
                exception=addCause(exception,MException( ...
                    "radia:optuna:ExceptionNote","%s",obj.notes(index)));
            end
            throw(exception);
        end
    end

    methods (Static, Access=private)
        function value=formatArgs(values)
            if isempty(values)
                value="";
                return
            end
            if isscalar(values)
                value=radia.optuna.ExceptionValue.pythonString(values{1});
                return
            end
            parts=strings(1,numel(values));
            for index=1:numel(values)
                parts(index)=radia.optuna.ExceptionValue.pythonRepr( ...
                    values{index});
            end
            value="("+strjoin(parts,", ")+")";
        end

        function value=pythonString(input)
            if ischar(input) || (isstring(input) && isscalar(input))
                value=string(input);
            elseif islogical(input) && isscalar(input)
                value=string(input);
                value=upper(extractBefore(value,2))+extractAfter(value,1);
            elseif isnumeric(input) && isscalar(input)
                value=string(input);
            else
                value=string(jsonencode(input));
            end
        end

        function value=pythonRepr(input)
            if ischar(input) || (isstring(input) && isscalar(input))
                value="'"+replace(string(input),"'","\'")+"'";
            else
                value=radia.optuna.ExceptionValue.pythonString(input);
            end
        end
    end
end
