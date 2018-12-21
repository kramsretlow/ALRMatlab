function newY = CheckY(obj,Y)
%Checks any specified vector/matrix Y to see if it's reasonable to use it as a
%response.  Fixes problems if possible and throws warning/error if necessary. The
%supplied object Y can be either a numeric with at most 2 distinct values, a
%categorical with at most 2 categories, or a logical vector. Need to handle the case
%where all Y are the same separately to make sure values are properly matched to
%(low, high).
%Goal is to convert supplied Y into a vector/matrix taking values
%matching the coding. 

%- If Y is empty, just leave it empty and stop.
if isempty(Y)
    newY = Y;
    return
end

tmp = unique(Y);
newY = zeros(size(Y));
L = obj.Coding(1);
H = obj.Coding(2);

%- If Y is logical, just assign the low/high coding values to 0 and 1.
if islogical(Y)
    newY(~Y) = L;
    newY(Y) = H;
    return
end

%- If Y is numeric with one value:
%   * assign the right value if it matches one of L or H
%   * otherwise assign low.
%- If Y is numeric with two values, assign L to the lower ones, H to the
%  higher ones.
if isnumeric(Y)
    switch length(tmp)
        case 2
            newY(Y==tmp(1)) = L;
            newY(Y==tmp(2)) = H;
            if tmp(1)~=L || tmp(2)~=H
                warning(['Assigned Y using values different than the '...
                    'coding. The configuration has been updated '...
                    'but the coding has not been changed.'])
            end
        case 1
            if tmp==L || tmp==H
                newY(:) = tmp;
            else
                newY(:) = L;
                warning(['The supplied Y has only one value, not '...
                    'matching the coding. The low value was '...
                    'used instead.'])
            end
        otherwise
            error('Setting Y: Y must have exactly 1 or 2 distinct values')
    end
    return
end

%- If Y is categorical with one category:
%   * assign the right value if that category matches a state
%   * otherwise assign low
%- If Y is categorical with two categories:
%   * do the assignment if they match States
%   * otherwise assign the first category to L, the second to H.
if iscategorical(Y)
    switch length(tmp)
        case 2
            if tmp(1)==obj.States{1} && tmp(2)==obj.States{2}
                newY(Y==obj.States{1}) = L;
                newY(Y==obj.States{2}) = H;
            else
                newY(Y==tmp(1)) = L;
                newY(Y==tmp(2)) = H;
                warning(['Assigned Y using states that differ from '...
                    'the model''s states. The configuration has '...
                    'been changed but the states have not.'])
            end
        case 1
            if tmp==obj.States{1}
                newY(:) = L;
            elseif tmp==obj.States{2}
                newY(:) = H;
            else
                newY(:) = L;
                warning(['The supplied Y has only one state, not '...
                    'matching the existing states. The low value '...
                    'was used instead.'])
            end
        otherwise
            error('Setting Y: Y must have exactly 1 or 2 distinct values')
    end
    return
end
end

