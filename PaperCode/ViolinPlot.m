function p = ViolinPlot( data, k, x0, w, color )
%VIOLINPLOT Adds a type of violin plot to the current axes.
%   This violin plot is based on an unequal-bin histogram.  There are k bins, with
%   bin edges placed at k+1 percentiles of data from 0 to 100%.  The plot is drawn 
%   vertically, placed on the current axes as a patch object at horizontal position 
%   x0, and it has width 2*w.
%
% Inputs:
%   data    The data vector
%   k       The number of bins to use for the histogram
%   x0      The horizontal location of the violin plot
%   w       The half-width of the plot
%   color   A [r g b] vector used in drawing the plot
%
% Output:
%   p       A handle to the patch object.

% Create the edges for the histogram
g = prctile(data, linspace(0,100,k+1));

% Get the histogram heights (scaled to max height w)
N = histcounts(data,g);
f = N./sum(N)./diff(g);  %-pdf normalization
f1 = x0 + w*f/max(f);
f2 = x0 - w*f/max(f);

% Make vectors needed to draw the patch
ixf1 = repelem(1:k,2);
ixf2 = fliplr(ixf1);
ixg1 = [1 repelem(2:k,2) k+1];
ixg2 = fliplr(ixg1);
xdata = [f1(ixf1) f2(ixf2)];
ydata = [g(ixg1) g(ixg2)];

% Plot the patch
p = patch(xdata,ydata,color);
hold on;
plot(x0,median(data),'k.')

end

