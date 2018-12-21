function out = GibbsSample(obj,varargin)
%GIBBSSAMPLE Does k iterations of a Gibbs sampler for autologistic model obj.
%   out = GibbsSample(obj)  Performs one GS step.
%   out = GibbsSample(obj, k)  Carries out k iterations of the sampler and returns 
%         the average of the k samples.
%   out = GibbsSample(obj, 'param', 'value')  Does Gibbs sampling with options
%         specified by parameter-value pairs (see below)
%
% Required arguments: 
%   obj   An object with class AutoLogistic.
% Optional arguments:
%   k     A positive integer number of iterations (default: 1)
% Param/value pairs:
%   'avg'     Logical. Should the average be returned, or all samples? Default: true.
%             NOTE: if avg==true, the proportion of 'high' draws is reported, not the 
%             numerical average (which is coding dependent).
%   'start'   Logical, categorical, or numeric vector. The starting configuration to 
%             use for sampling. Default: random.
%   'burnin'  Scalar.  The number of draws to discard at the beginning of the run.
%             Default: 0.
%   'verbose' Logical. Should progress be output to the console as computation
%             proceeds?
%
% Output:
%   out   If k==1 or avg==true, an obj.N-by-1 vector giving the draw or the average 
%         of k draws. Otherise, an obj.N-by-k array of samples, one per column.
%
%TODO:
%- Look for speedups, possibly even for certain cases like the 4-connected grid,
%  where we can use the checkerboard trick to speed up.

%---Parse the inputs-----------------------------
IP = inputParser;
IP.FunctionName = 'GibbSample';
addRequired(IP,'obj',@(x) isa(x,'AutoLogistic'))
addOptional(IP,'k',1,@(x) isscalar(x) && mod(x,1)==0)
addParameter(IP,'avg',true,@(x) validateattributes(x,{'logical'},{'scalar'}))
addParameter(IP,'start',[], ...
             @(x) validateattributes(x,{'numeric','logical','categorical'},...
             {'vector','numel',obj.N}))
addParameter(IP,'burnin',0,@(x) isscalar(x) && mod(x,1)==0)
addParameter(IP,'verbose',false,@(x) validateattributes(x,{'logical'},{'scalar'}))
parse(IP,obj,varargin{:})            

%---Get needed objects---------------------------
% Parts of the AutoLogistic object
A = IP.Results.obj.AssociationMatrix;
mu = IP.Results.obj.Mu;
alpha = IP.Results.obj.Alpha;
lo = IP.Results.obj.Coding(1);
hi = IP.Results.obj.Coding(2);
N = IP.Results.obj.N;
k = IP.Results.k;
% The starting value (call it Y)
if isempty(IP.Results.start)
    pick = binornd(1,0.5,[obj.N,1]);
    Y = pick;
    Y(pick==0) = lo;
    Y(pick==1) = hi;
else
    Y = IP.Results.start;
end
if isrow(Y)
    Y = Y';
end
Y = CheckY(IP.Results.obj,Y);  %Make Y match the object's coding.
% Other inputs
avg = IP.Results.avg;
burnin = IP.Results.burnin;
verbose = IP.Results.verbose;

%---Do the sampler-------------------------------
% Loop through k samples. In each sample, loop through all nodes in the graph.
temp = zeros(N,k);
for j = 1:k
    for i = 1:N;
        nbrsum = A(i,:)*(Y-mu); 
        etalo = exp(lo*(alpha(i) + nbrsum));
        etahi = exp(hi*(alpha(i) + nbrsum));
        if etahi==Inf
            p_i = 1;   %In rare cases, hi*(alpha(i) + nbrsum) could be too large.
        else
            p_i = etahi/(etalo + etahi);
        end
        u = rand();
        if u<p_i
            Y(i) = hi;
        else
            Y(i) = lo;
        end
    end
    temp(:,j) = Y;
    if verbose
        disp(['finished draw ', num2str(j), ' of ', num2str(k)])
    end
end

%---Finish up------------------------------------
temp(:,1:burnin) = [];   %-Does nothing if burnin==0.
if avg
    out = sum(temp==hi,2)/(k-burnin);
else
    out = temp;
end


end

