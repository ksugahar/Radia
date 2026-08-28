classdef TrialRowIndex < handle
    %TRIALROWINDEX Trial-number buckets over one append-mostly history store.
    %   Optimization history is written once per suggest / report / attribute
    %   and read back per trial by freezeTrial, the samplers, and the
    %   pruners. Answering each read with a full scan of the store's trial
    %   column costs O(rows) per read, so a study that freezes or samples
    %   once per trial spends O(rows^2) overall. This keeps one bucket of row
    %   indices per trial number, turning a read into O(rows for that trial).
    %
    %   The index is a cache, never an authority: every lookup is handed the
    %   store's own trial column and rebuilds whenever the two disagree, so a
    %   caller that forgets to announce a write loses speed, never
    %   correctness. Rows come back in ascending order, exactly as find()
    %   would return them.

    properties (Access=private)
        Buckets cell = cell(0,1)
        IndexedRows (1,1) double = 0
        Usable (1,1) logical = false
    end

    methods
        function append(obj, trialNumber, row)
            %APPEND Record that ROW was appended for TRIALNUMBER.
            if ~obj.Usable
                return
            end
            slot = double(trialNumber) + 1;
            if ~(isscalar(slot) && isfinite(slot) && slot >= 1 && ...
                    slot == floor(slot))
                obj.invalidate();
                return
            end
            if slot > numel(obj.Buckets)
                obj.Buckets{slot, 1} = zeros(0, 1);
            end
            obj.Buckets{slot} = [obj.Buckets{slot}; double(row)];
            obj.IndexedRows = obj.IndexedRows + 1;
        end

        function invalidate(obj)
            %INVALIDATE Drop the cache after a deletion or a bulk rewrite.
            obj.Buckets = cell(0, 1);
            obj.IndexedRows = 0;
            obj.Usable = false;
        end

        function rows = lookup(obj, keys, trialNumber)
            %LOOKUP Ascending row indices of TRIALNUMBER within KEYS.
            if ~obj.refresh(keys)
                rows = find(reshape(double(keys), [], 1) == double(trialNumber));
                return
            end
            slot = double(trialNumber) + 1;
            if ~(isscalar(slot) && isfinite(slot) && slot >= 1 && ...
                    slot == floor(slot)) || slot > numel(obj.Buckets)
                rows = zeros(0, 1);
                return
            end
            rows = obj.Buckets{slot};
            if isempty(rows)
                rows = zeros(0, 1);
            end
        end

        function present = has(obj, keys, trialNumber)
            %HAS True when TRIALNUMBER owns at least one row.
            present = ~isempty(obj.lookup(keys, trialNumber));
        end
    end

    methods (Access=private)
        function usable = refresh(obj, keys)
            keys = reshape(double(keys), [], 1);
            if obj.Usable && obj.IndexedRows == numel(keys)
                usable = true;
                return
            end
            obj.Buckets = cell(0, 1);
            obj.IndexedRows = 0;
            obj.Usable = false;
            if isempty(keys)
                obj.Usable = true;
                obj.IndexedRows = 0;
                usable = true;
                return
            end
            slots = keys + 1;
            if any(~isfinite(slots)) || any(slots < 1) || ...
                    any(slots ~= floor(slots))
                % Trial numbers outside the nonnegative integers cannot be
                % bucketed; the caller falls back to scanning the column.
                usable = false;
                return
            end
            % accumarray does not promise the order it hands rows to the
            % reducer, so sort each bucket back into find() order.
            obj.Buckets = accumarray(slots, (1:numel(keys))', ...
                [max(slots), 1], @(rows) {sort(rows)}, {zeros(0, 1)});
            obj.IndexedRows = numel(keys);
            obj.Usable = true;
            usable = true;
        end
    end
end
