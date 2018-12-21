function out = PerfectSample2(obj)
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
T = 2;                                    %-Starting number of time steps back to go.
seeds = randi(1e6,101,1);   %-Hold the seeds needed to reproduce the rand num stream.
coalesce = false;

%---Do the sampler-------------------------------
%Note: In this version we will not save the random numbers in memory.  Instead
%keep track of the seeds used.  seeds(j) will hold the seed used to generate samples
%from time -2^(j-1)T to -2^(j-2)+1.  We'll cap j at 100 as T*2^100 is a huge number
%of time steps. We'll call j the "epoch index." It tells us how far back in time we
%need to start (specifically, we start at time -T*2^(j-1)).
j = 0;     
while ~coalesce
    j = j + 1;
    L = lo*ones(N,1);
    H = hi*ones(N,1);
    [H, L] = DoSample(j,L,H);
    coalesce = all(L==H);
    disp(['Start from ', num2str(-T*2^(j-1)), ': ', num2str(sum(H~=L)), ' elements different.'])
end
if coalesce
    out = L;
else
    out = NaN(N,1);
    warning('Sampler did not coalesce. Returning NaNs.')
end


    %===Subfunction to run sampler forward for j epochs==============================
    function [L, H] = DoSample(j,L,H)
    % j is the number of epochs to go back.
    % L, H are the low and high chains.
    rng(seeds(j))
    if j==1
        for t = -T:0
            [L, H] = GibbsStep(L,H);
        end
    else
        for t = -T*2^(j-1) : -T*2^(j-2)-1
            [L, H] = GibbsStep(L,H);
        end
        [L,H] = DoSample(j-1,L,H);
    end
    end
    %================================================================================


    %===Subfunction to do GS updates=================================================
    function [L, H] = GibbsStep(L,H)
    %Performs a single update step of both the lower and upper chains, using common 
    %random variables.
    U = rand(N,1);
    for i = 1:N;
        %Do the lower chain update.
        nbrsum = A(i,:)*(L-mu); 
        etalo = exp(lo*(alpha(i) + nbrsum));
        etahi = exp(hi*(alpha(i) + nbrsum));
        p_i = etahi/(etalo + etahi);
        if U(i)<p_i
            L(i) = hi;
        else
            L(i) = lo;
        end
        %Do the lower chain update.
        nbrsum = A(i,:)*(H-mu); 
        etalo = exp(lo*(alpha(i) + nbrsum));
        etahi = exp(hi*(alpha(i) + nbrsum));
        p_i = etahi/(etalo + etahi);
        if U(i)<p_i
            H(i) = hi;
        else
            H(i) = lo;
        end
    end        
    end
    %================================================================================

end




