function Z = standardize(X)
% Function to standardize a matrix of data.  Data is assumed to have variables
% in columns, cases in rows.  Each variable is standardized by subtracting the
% mean from all observations, and dividing by the standard deviation.
n = size(X,1);
Xbar = mean(X)';                      
StDev = std(X)';                      
Z = (X-repmat(Xbar',[n 1]))./repmat(StDev',[n,1]);
end