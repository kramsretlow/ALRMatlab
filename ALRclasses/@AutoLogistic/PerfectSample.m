function out = PerfectSample(obj)
%PERFECTSAMPLE draws a perfect sample from an autologistic model using CFTP.
%   out = PerfectSample(obj)  Draws one sample.
%   out = GibbsSample(obj, 'param', 'value')  Does Gibbs sampling with options
%         specified by parameter-value pairs (see below)
%
% Required arguments: 
%   obj   An object with class AutoLogistic.
%
% Output:
%   out   An obj.N-by-1 vector giving the draw.
%
%TODO:


%---Get needed objects---------------------------

% Parts of the AutoLogistic object
A = obj.AssociationMatrix;
mu = obj.Mu;
alpha = obj.Alpha;
lo = obj.Coding(1);
hi = obj.Coding(2);
N = obj.N;

% Initialize perfect sampler objects
T = 1;  %-Starting number of time steps back to go.
U = rand(N,2);  %-Initial pool of uniform random variates.
coalesce = false;

%---Do the sampler-------------------------------
%Note: we use matrix U to hold all the uniform random numbers needed to compute the
%lower and upper chains as stochastic recursive sequences. In this matrix each column
%holds the N random variates needed to do a full Gibbs sampling update of the
%variables in the graph. We think of the columns as going backwards in time to the
%left: column one is time 0, column 2 is time -1, ... column T+1 is time -T.  So to
%run the chains in forward time we go from right to left.
while ~coalesce
    
    L = lo*ones(N,1);                             %-Initialize the lower bound chain.
    H = hi*ones(N,1);                             %-Initialize the upper bound chain.   
    U = [U rand(N,T)];           %-Add T more columns of random numbers to the right.
    T = 2*T;                                                              %-Double T.
    for t = T+1:-1:1                           %-column t corresponds to time -(t-1).
        L = GibbsStep(L,t);
        H = GibbsStep(H,t);
    end
    coalesce = all(L==H);
    disp(['T = ', num2str(T), ': ', num2str(sum(H~=L)), ' elements different.'])
end
%disp(['Finished at T = ', num2str(T)])

out = L;

    %===Subfunction to do GS updates=================================================
    function Vout = GibbsStep(V,j)
    %Performs a single update step of all variables in obj, starting from
    %configuration V, and using random variables U(:,j)
    Vout = zeros(N,1);
    for i = 1:N;
        nbrsum = A(i,:)*(V-mu); 
        etalo = exp(lo*(alpha(i) + nbrsum));
        etahi = exp(hi*(alpha(i) + nbrsum));
        p_i = etahi/(etalo + etahi);
%         if U(i,j)<1-p_i
%             Vout(i) = lo;
%         else
%             Vout(i) = hi;
%         end
        if U(i,j)<p_i
            Vout(i) = hi;
        else
            Vout(i) = lo;
        end
    end        
    end
    %================================================================================

end




