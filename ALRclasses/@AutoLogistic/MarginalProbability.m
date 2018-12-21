function out = MarginalProbability(obj, nodes)
%MARGINALPROBABILITY  Obtains the marginal probs of specified nodes being 'high'.
%   out = MarginalProbability(obj)  Gets the marginal probabilities of all nodes.
%   out = MarginalProbability(obj,nodes)  Gets marginal probabilities of nodes in the
%       supplied vector.
%
%   If obj.N is greater than 20, this is done by simulation.
%
%   obj     An AutoLogistic object.
%   nodes   A vector of node indices.
%   out     A vector of probabilities.

%TODO: input checking

if nargin==1
    nodes = 1:obj.N;
end
if obj.N <= 20
    T = obj.ProbabilityTable();
    probs = table2array(T(:,end));
    highs = table2array(T(:,1:obj.N)) == obj.Coding(2);
    weights = bsxfun(@times,highs,probs);
    marg = sum(weights);
    out = marg(nodes)';
else
    %TODO: implement the sampling way
end

end

