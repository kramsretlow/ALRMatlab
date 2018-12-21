function tbl = ProbabilityTable(obj,force)
%PROBABILITYTABLE Tabulates probabilities of all configurations.
%   ProbabilityTable(obj)  Tabulates the PMF of AutoLogistic model obj. If obj.N>20, 
%       an error is thrown to avoid run time or memory issues.
%   ProbabilityTable(obj,force)  If force evaluates to true, carries out the
%       tabulation regardless of the value of obj.N.

%Check for N too large
n = obj.N;
nr = 2^n;
if (nargin==1 || ~force) && n>20
    error(['Attempting to tabulate the PMF with more than 2^20 configurations. ' ...
           'You probably do not want to do this. Set argument ''force'' to true ' ...
           'to do it anyway.'])
elseif nargin==2 && force && n>20
    warning(['Tabulating the PMF with more than 2^20 configurations. This ' ...
             'could exceed time or memory limitations.'])
end

%Create a matrix holding all configurations. Include an extra column for
%probabilties. .
Y = zeros(nr,n+1);
L = obj.Coding(1);
H = obj.Coding(2);
for i=1:n
    Y(:,i) = repmat(repelem([H L]',nr/2^i), [2^(i-1) 1]);
end

%Get exp(negpotential) values for each configuration
Y(:,end) = exp(obj.Negpotential(Y(:,1:n)'));

%Normalize and put everything into a table
Y(:,end) = Y(:,end)/sum(Y(:,end));
tbl = array2table(Y);
tbl.Properties.VariableNames{end} = 'Probability';

end