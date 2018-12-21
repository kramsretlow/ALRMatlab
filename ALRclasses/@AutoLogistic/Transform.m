function out = Transform(obj,coding,centered)
%TRANSFORM Perform parameter transformation between AL (not ALR) models. 
%   out = Transform(obj,coding,centered) 
%
%   Objects of classes descended from AutoLogistic have properties Alpha
%   and AssociationMatrix.  This function computes the values of Alpha and
%   AssociationMatrix for an equivalent AL model with a different coding and
%   centering. For the case of transformation from a centered model, the solution is
%   found using a fixed-point iteration scheme; otherwise it is computed explicitly.
%
%   inputs:
%     obj       An object with class descended from AutoLogistic.  Its members Alpha,
%               AssociationMatrix, Coding, and Centered are used to determine the
%               characteristics of the model being transformed from.
%     coding    A 2-vector giving the [low high] variable coding to transform into.
%     centered  A logical value determining the type of model to transform into.  
%               True for centered model and false for a standard one.
%   output out is a struct with elements:
%     alpha     A vector the same lenghth as obj.Alpha, giving the new unary
%               parameters.
%     assoc     If obj is of class ALRsimple, a scalar giving the new pairwise
%               parameter. Otherwise, a new association matrix.
%     niter     Number of iterations used to find the solution
%     a, b      The values of a = (H-L)/(hi-lo) and b = L - a*lo.

%Input checking
validateattributes(coding,{'numeric'},{'vector','numel',2,'increasing'},'Transform')
validateattributes(centered,{'logical'},{'size',[1 1]},'Transform')

% Here we're imagining transforming from an AL model with parameters phi and Alpha,
% having coding (L,H), to one with parameters alpha, assoc, and coding (lo,hi).
phi = obj.Alpha;
Omega = obj.AssociationMatrix;
L = obj.Coding(1);
H = obj.Coding(2);
lo = coding(1);
hi = coding(2);

% Compute required intermediate values
a = (H-L)/(hi-lo);
b = L - a*lo;
mu = obj.Mu;
one = ones(obj.N,1);

% Create output struct
out = struct('alpha',[],'assoc',[],'niter',0,'a',a,'b',b);

% Get the unary parameter vector. When transforming to standard model, just get the
% answer. Otherwise do the fixed point iteration.
if ~centered
    out.alpha = a * (phi + Omega*(b*one - mu));
else
    const = a*phi + a*Omega*(b*one - mu);
    newmu = @(u) (lo*exp(u*lo) + hi*exp(u*hi)) ./ (exp(u*lo) + exp(u*hi));
    g = @(u) a^2*Omega*newmu(u) + const;
    u_old = phi;
    tol = 1e-6;
    maxiter = 250;
    for i=1:maxiter
        u_new = g(u_old);
        err = max(abs(u_new-u_old));  %-OK if 1 is typical magnitude (Lange10, p. 70)
        u_old = u_new;
        if err < tol
            out.niter = i;
            break;
        end
        if i==maxiter
            out.niter = i;
            warning('Search for unary parameter vector did not converge.')
        end
    end
    out.alpha = u_old;
end

% Put in the pairwise parameter or association matrix
if isa(obj,'ALRsimple')
    out.assoc = a^2*obj.Lambda;
else
    out.assoc = a^2*Omega;
end    


end


