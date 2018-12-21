function ax = mysubplot(n,m,p,h,v)
%=============== SUBPLOT FUNCTION WITH USER-CONTROLLED SPACING =================
%
%Example Calls:
%   mysubplot(n,m,p)        Creates a subplot in the current figure in the same
%                           way as built-in subplot function, but with 1mm  
%                           spacing between figures and at borders.
%   mysubplot(n,m,p,h,v)    Creates a subplot with horizontal spacings given by
%                           vector h (specified left to right) and vertical
%                           spacings given by vector v (specified top to
%                           bottom).  All spacings are in cm.
%   ax = mysubplot(...)     Also returns the handle to the new axes.
%
%Notes:
%   * mysubplot operates on the figure currently active in the calling 
%     workspace.  If there is no such figure, one is created.
%   *TOFIX/TODO:  -Make it handle specification of unequal-sized figures correctly.
%   *****USE subplot('Position',[l b w h]) TO DO THIS*****
%
%KEYWORDS: MarkWolters Keyword1 Keyword2 ...
%
%===============================================================================

%-Get or create the parent figure-----------------------------------------------
fig = gcf;
set(fig,'Units','centimeters');
P = get(fig,'Position');
W = P(3);
H = P(4);
%-Input checking----------------------------------------------------------------
ok1 = max(p)<=n*m;
if nargin==3
    h = ones(m+1,1);
    v = ones(n+1,1);
    ok2 = true;
elseif nargin==5
    ok2 = length(h)==m+1 && length(v)==n+1;
else 
    ok2 = false;
end
ok3 = sum(h)<W && sum(v)<H;
if ~(ok1&&ok2&&ok3)
    error('Inputs are incorrectly specified')
end
%--Create the subplot-----------------------------------------------------------
% Note 1: we are using cm as units, but subplot('Position',[l b w h]) requires 
%   normalized units--otherwise get subplots overwriting each other.
% Note 2: the whole plot is envisioned as an n-by-m grid of plot areas.  We need
%   to handle the case where a subplot is to cover more than one of the areas.
h = h/W;                                %-Convert h to normalized units.
v = v/H;                                %-Convert v to normalized units.
plw = (1-sum(h))/m;                     %-Individual plot width.
plh = (1-sum(v))/n;                     %-Individual plot height.

r = ceil(p/m);                          %-Row each specified plot is in.
c = p-m*(r-1);                          %-Col each specified plot is in.
minr = min(r);                          %-Minimum row.
maxr = max(r);                          %-Max row.
minc = min(c);                          %-Min col.
maxc = max(c);                          %-Max col.
left = (minc-1)*plw + sum(h(1:minc));
bot = (n-maxr)*plh + sum(v(maxr+1:end));
if maxc==minc
    wid = plw;
else
    wid = (maxc-minc+1)*plw + sum(h(minc:maxc-1));
end
if maxr==minr;
    ht = plh;
else
    ht = (maxr-minr+1)*plh + sum(v(minr+1:maxr));
end
ax = subplot('Position',[left bot wid ht]);

end

