function out = PLobj(param, mod)
% Objective function for maximum pseudolikelihood estimation of ALR models. 
%   param = parameter vector [beta', lambda]' where beta is a p-vector and lambda 
%           is a scalar. 
%   mod = the ALRsimple object.  

    mod.Beta = param(1:end-1);
    mod.Lambda = param(end);
    out = mod.PseudoLikelihood;

end