function out = PseudoLikelihood(obj)
%PSEUDOLIKELIHOOD  Computes the negative log pseudolikelihood for the given
%autologistic model.
%   out = PseudoLikelihood(obj)  Computes the neg log pseudolikelihood.
%
% TODO: check: need to worry about logPL = -Inf if loterm + hiterm almost zero?


% Get needed quantities and initialize things.
if ~obj.DimensionsOK
    error('Can''t compute pseudolikelihood: dimensions inconsistent.')
end
out = 0;

% Loop through replicates
for m = 1:obj.M

    Ym = obj.Y(:,m);                          %-The current replicate's observations.
    alpha = obj.Alpha(:,m);                   %-Current replicate's unary parameters.
    mu = obj.Mu(1:obj.N, m);                  %-Current replicate's centering terms.

    % Get the (lambda-weighted) neighbour sums and add to the unary parameters.
    s = obj.AssociationMatrix * (Ym - mu);
    a_plus_s = alpha + s;
    
    % Get this replicate's log pseudolikelihood
    loterm = exp(obj.Coding(1)*a_plus_s);
    hiterm = exp(obj.Coding(2)*a_plus_s);
    logPL = sum(Ym.*a_plus_s - log(loterm + hiterm));
    
    % Subtract this replicate's log PL from the total.
    out = out - logPL;

end

end

