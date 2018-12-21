function out = Mu(obj,nodes,reps)
%MU A function to compute the vector of centering adjustments.
%   TODO: Detailed explanation goes here
%   TODO: input checking (e.g. invalid nodes, reps... worried about speed though)
%   obj   An autologistic model object
%   nodes A vector of indices of the nodes to use (if not supplied, return all).
%   reps  a vector of indices of the replicates to use (if not supplied, return all).
%
% Notes:
% - For centered models, Mu is sized to match the size of Alpha (which depends on X
%   and Beta). For standard models, Mu is sized to be obj.N by obj.M (all zeros).
%   These should be the same, but might not be if certain parts of object aren't
%   initialized yet.

if obj.Centered
    L = obj.Coding(1);
    H = obj.Coding(2);
    out = (L*exp(obj.Alpha*L) + H*exp(obj.Alpha*H)) ./ ...
          (exp(obj.Alpha*L) + exp(obj.Alpha*H));
else
    out = zeros(obj.N,obj.M);
end

if nargin==2
    out = out(nodes,:);
end
if nargin==3
    out = out(nodes,reps);
end

end

