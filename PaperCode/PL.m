function [f, g] = PL(theta, M)
%PLOBJ Compute the negative log pseudolikelihood function (and optionally its
%gradient) for ALR model M at parameter values theta.  Only works
%for cases with no replicates.

P = length(theta) - 1;
M.Beta = theta(1:P);
lam = theta(end);
M.Lambda = lam;
Y = M.Y;                              %-The observations.
X = M.X;                              %-The covariates.
alpha = M.Alpha;                      %-The unary parameters.
mu = M.Mu;                            %-The centering terms.
lo = M.Coding(1);
hi = M.Coding(2);

% Get the (lambda-weighted) neighbour sums and add to the unary parameters.
A = M.AssociationMatrix;  %NB: this includes the lambda, it's not the adj matrix.
Adjmat = M.Graph.adjacency;
s = A * (Y - mu);
a_plus_s = alpha + s;

% Get the log pseudolikelihood
loterm = exp(lo*a_plus_s);
hiterm = exp(hi*a_plus_s);
f = -sum(Y.*a_plus_s - log(loterm + hiterm));


if nargout > 1 %gradient required
    
    g = zeros(P + 1, 1);
    
    m = (lo*exp(lo*a_plus_s) + hi*exp(hi*a_plus_s)) ...
        ./ (exp(lo*a_plus_s) + exp(hi*a_plus_s));
    r = exp((lo+hi)*alpha) ./ (exp(lo*alpha) + exp(hi*alpha)).^2;
    
%     for i = 1:P
%         g(i) = X(:,i)'*(m-Y) + (hi-lo)^2 * (Y-q)'*A*(X(:,i).*r);
%     end
    
    g(1:P) = X'*(m-Y) + M.Centered*(hi-lo)^2 * (X.*repmat(r,[1 P]))'*A*(Y-m);
    
    g(end) = (m-Y)'*Adjmat*(Y-mu);
    
end


end

