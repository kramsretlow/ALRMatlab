function [ TVdist ] = GetTV(theta, M, prob)
%GETTV Compute total variation distance between the model M with parameters theta, 
%         and probability vec prob. (for the small-scale network-data example) 
% Inputs:
%   theta   A column vector of parameters [beta0 beta1 lambda] to use in the CZO model.
%   CZO     A centered, zero/one type ALRsimple model with the appropriate graph and 
%           covariates.
%   prob    A probability vector such as generated by .ProbabiltyTable method of 
%           the ALR class.  Should be same length as rows of CZO.ProbabilityTable.


M.Beta = theta(1:2);
M.Lambda = theta(3);
tbl = M.ProbabilityTable;

TVdist = sum(abs(tbl.Probability-prob))/2;


end


