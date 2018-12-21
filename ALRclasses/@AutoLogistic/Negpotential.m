function out = Negpotential(obj, Y)
%NEGPOTENTIAL Computes the negpotential function Q() of an AL or ALR model at a
%particular configuration.
%   Negpotential(obj)  computes Q() at the value of Y stored in the object.
%   Negpotential(obj,Y)  computes Q() at the given Y.
%
%   obj is an ALmodel or ALRmodel object.
%   Y is an N-by-M matrix.  M is the number of replicate observations (commonly 1).
%   out is an M-by-1 vector of negpotential values.

if nargin==1
    Y = obj.Y;
elseif nargin==2
    %TODO: check that supplied Y matches the coding!
else
    error('Too many arguments.')
end

nvals = size(Y,2);
out = zeros(nvals,1);
mu = Mu(obj);
a = obj.Alpha;
AM = obj.AssociationMatrix;

for i=1:nvals
    %TODO: it might be possible to replace this loop using functions employing 
    %sum(bsxfun(@times,blah1,blah2), 2).  Look into possible speed improvements if
    %neccesary later.
    y = Y(:,i);
    out(i) = y'*a - y'*AM*mu + y'*AM*y/2;
end

end
