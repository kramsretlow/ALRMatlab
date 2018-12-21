function fig = myfigure(num,Dimensions)
% Utility function to create a figure that is well prepared for printing.  
% 
% fig = myfigure(num,Dimensions) creates a new figure using figure(num), then sets it
% to have size given by Dimensions, in cm, with the paper size set to be 3mm
% larger on each side.  It is then a simple matter to print the figure using a
% command like print(fig,'-r600','-depsc2','FigName') to create an eps file of with
% the correct dimensions.

fig = figure(num);
Position = [1 1 Dimensions(1) Dimensions(2)];
PaperSize = [Dimensions(1)+0.6 Dimensions(2)+0.6];
set(fig,'Units','centimeters',...
        'Position',Position,...
        'PaperUnits','centimeters',...
        'PaperSize',PaperSize,...
        'PaperPosition',[0.3 0.3 Dimensions])



end

