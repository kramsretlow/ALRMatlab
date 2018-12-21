function out = ConditionalProbability(obj, nodes, reps)
%CONDITIONALPROBABILITY  Computes the conditional probability that variables take
%their observed values, for the given autologistic model.  Optional arguments can 
%control which nodes and replicates have probabilities returned. 
%   out = ConditionalProbability(obj)  Computes the marginal probabilities for all 
%         nodes. out is an nn-by-mm matrix, where nn is either obj.N or 
%         length(nodes), and mm is either obj.M or length(reps).
%   out = ConditionalProbability(obj, nodes) Only returns probabilities for the nodes
%         with specified indices.
%   out = ConditionalProbability(obj, nodes, reps) Only returns probabilities for the
%         nodes and replicates specified.

% TODO: input checking

% Initialize objects
hi = obj.Coding(2);
out = zeros(obj.N, obj.M);

% Loop through replicates and get the conditional probabilities.
for m = 1:obj.M
    Ym = obj.Y(:,m);                          %-The current replicate's observations.
    alpha = obj.Alpha(:,m);                   %-Current replicate's unary parameters.
    mu = obj.Mu(1:obj.N, m);                   %-Current replicate's centering terms.
    
    % Get the (lambda-weighted) neighbour sums and add to the unary parameters.
    s = obj.AssociationMatrix * (Ym - mu);
    a_plus_s = alpha + s;
    
    % Get this replicate's conditional probabilities
    Yterm = exp(Ym.*a_plus_s);
    loterm = exp(obj.Coding(1)*a_plus_s);
    hiterm = exp(obj.Coding(2)*a_plus_s);
    out(:,m) = Yterm ./ (loterm + hiterm);
end

if nargin==2
    out = out(nodes,:);
end
if nargin==3
    out = out(nodes,reps);
end


end

