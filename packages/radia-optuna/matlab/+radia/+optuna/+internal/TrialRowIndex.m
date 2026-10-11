classdef TrialRowIndex < handle
    %TRIALROWINDEX Trial-number buckets over one append-mostly history store.
    %   Per-trial scans make freezing a whole study quadratic in its history
    %   size. This cache keeps ascending row buckets for nonnegative trial
    %   numbers. The source key column remains authoritative: LOOKUP rebuilds
    %   the cache whenever its row count no longer matches the source.

    properties (Access=private)
        Buckets cell = cell(0,1)
        IndexedRows (1,1) double = 0
        Usable (1,1) logical = false
    end

    methods
        function append(obj,trialNumber,row)
            if ~obj.Usable
                return
            end
            slot=double(trialNumber)+1;
            if ~(isscalar(slot) && isfinite(slot) && slot>=1 && ...
                    slot==floor(slot) && row==obj.IndexedRows+1)
                obj.invalidate();
                return
            end
            if slot>numel(obj.Buckets)
                obj.Buckets{slot,1}=zeros(0,1);
            end
            obj.Buckets{slot}=[obj.Buckets{slot};double(row)];
            obj.IndexedRows=obj.IndexedRows+1;
        end

        function invalidate(obj)
            obj.Buckets=cell(0,1);
            obj.IndexedRows=0;
            obj.Usable=false;
        end

        function rows=lookup(obj,keys,trialNumber)
            if ~obj.refresh(keys)
                rows=find(reshape(double(keys),[],1)==double(trialNumber));
                return
            end
            slot=double(trialNumber)+1;
            if ~(isscalar(slot) && isfinite(slot) && slot>=1 && ...
                    slot==floor(slot)) || slot>numel(obj.Buckets)
                rows=zeros(0,1);
                return
            end
            rows=obj.Buckets{slot};
            if isempty(rows)
                rows=zeros(0,1);
            end
        end
    end

    methods (Access=private)
        function usable=refresh(obj,keys)
            keys=reshape(double(keys),[],1);
            if obj.Usable && obj.IndexedRows==numel(keys)
                usable=true;
                return
            end
            obj.invalidate();
            if isempty(keys)
                obj.Usable=true;
                usable=true;
                return
            end
            slots=keys+1;
            if any(~isfinite(slots)) || any(slots<1) || ...
                    any(slots~=floor(slots))
                usable=false;
                return
            end
            obj.Buckets=accumarray(slots,(1:numel(keys))', ...
                [max(slots),1],@(rows) {sort(rows)},{zeros(0,1)});
            obj.IndexedRows=numel(keys);
            obj.Usable=true;
            usable=true;
        end
    end
end
