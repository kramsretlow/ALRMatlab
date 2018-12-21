function out = howclose(param,mod, probvec)
% Objective function for optimizing an ALRsimple model to match a probability table.
% Only works for 2-variable case.
%   param = parameter vector [alpha_1, alpha_2, lambda]'
%   mod = the ALRsimple object.  Needs to have Beta = 1 and the right graph
%   probvec = probability vector [phh phl plh pll]'

    mod.X = param(1:2);
    mod.Lambda = param(3);
    tbl = mod.ProbabilityTable;
    out = norm(probvec-tbl.Probability);


end