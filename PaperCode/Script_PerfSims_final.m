%% FINAL numeric comparisons for the Autologistic Plus/Minus Paper
% Building on what was learned in Script_Perfsims and Script_Perfsims2, include just
% the code needed for creating the publication figures.
% 
% It was decided that for the publication I'll only consider the straightforward
% "case 3" where the predictors are equal to the easting and northing coordinates.
% This is simple, agrees with previous papers' examples, and still exhibits all of
% the behaviours I want to demonstrate.

%% 1) Create objects to use later

% Misc variables
d = 30;            %-Square lattice dimension
beta = [-2 2 2]';  %-Regression coefficients
G = AutoLogistic.MakeGraph(d,d,1);  %-The graph.

% The covariates
X = ones(d^2,3);  %-Create the X matrix.
ivals = repmat(1:d,1,d)';  %-Easting index of nodes 1:d^2  
jvals = repelem(1:d,d)';   %-Northing index of nodes 1:d^2
Ei = 1/(d-1) * (ivals-1);  %-Easting coordinates of nodes on unit square.
Ni = 1/(d-1) * (jvals-1);  %-Northing coordinates of nodes on unit square.
X(:,2:3) = [Ei Ni];       %-Formerly called Case 3: pattern, no noise

% Create independence model with responses = high; then get endogenous probabilities.
% For plotting "large-scale structure," i.e., the indep-model probabilities.
ILM = ALRsimple;  %-"Independence Logistic Model"
ILM.Coding = [0 1];
ILM.Graph = G;
ILM.Beta = beta;
ILM.Lambda = 0;
ILM.X = X;
ILM.Y = ones(d^2,1);
Endo = ILM.ConditionalProbability;  %-P(Y==high| lambda = 0)

% Create SZO model
SZO = ALRsimple;
SZO.Coding = [0 1];
SZO.Graph = G;
SZO.Beta = beta;
SZO.X = X;

% Create SPM model 
% (recall we're using +/- 0.5 coding to make lambdas commensurate)
SPM = SZO;
SPM.Coding = [-1/2 1/2];

% Create CZO model
CZO = SZO;
CZO.Centered = true;


%% 2) Plot the graph and the endogenous structure

fig1 = myfigure(1,[8.5 7]);
fig1.Color = 'white';
clf
ax = axes;
imagesc([0 1], [0 1], reshape(Endo,d,d))
axis([-0.1 1.1 -0.1 1.1])
ax.Position = [0.1 0.186 0.775 0.815];
colormap(gray)
axis square
caxis([0 1])
set(gca,'YDir','normal')
cb = colorbar;
cb.Position = [0.8125 0.243 0.05 0.7];
cb.Label.String = 'Endogenous probability';
hold on
pl = plot(G,'-ko','XData',Ei,'YData',Ni);
pl.MarkerSize = 1;
pl.LineWidth = 0.25;
xlabel('x_1')
ylabel('x_2')
ax.Box = 'off';

%print(fig1,'-r600','-dpdf','-opengl','Sec4-EndoProbMap')
%print(fig1,'-dpdf','-painters','Sec4-EndoProbMap2')
%print(fig1,'-r600','-depsc2','Sec4-EndoProbMap3')
%print(fig1, 'Sec4-EndoProbMap', '-dtiff', '-opengl', '-r600')

%% 3) Plot average neighbour affinity and distance to endogenous vs. lambda
% Use the previously calculated values from 500 samples to do this. 

load('PefSimResults2.mat')  %-Retrieves object called expmt with the results.
expmt(expmt.scenario==4,:) = [];  %-Discard the "Case 4" runs we won't use here.

hspc = [1.2 3 3 1.2];     %-Horizontal spacings for mysubplot()
vspc = [0.7 1.2];         %-Vertical spacings for mysubplot()
fig2 = myfigure(2,[30 8]);
clf
models = {'SZO','SPM','CZO'};

for mm = 1:3  %-Loop through models
    
    %-Get the right rows of the table & right lambda, rho, helli values.
    ix = strcmp(expmt.model,models{mm});
    lams = expmt.lambda(ix);
    rhos = expmt.rhohat(ix);
    hellis = expmt.hellinger(ix);
    absdists = expmt.absdist(ix);
    
    %-Make the plot
    ax = mysubplot(1, 3, mm, hspc, vspc);
    plot(lams,rhos,'.k-')
    axis([0 lams(end) 0.5 1])
    xlabel('\lambda')
    ylabel('\rho')
    title(models{mm})
    ax2 = axes('Position',ax.Position);
    axis([0 lams(end) 0 1])
    ax2.Color = 'none';
    ax2.XAxis.Visible = 'off';
    ax2.YAxisLocation = 'right';
    ax2.FontWeight = 'bold';
    %*** plot only one of the two distance measures ***
    %line(lams,1-hellis,'Color','k','Linewidth',2,'Marker','.');
    line(lams,absdists,'Color','k','Linewidth',2,'Marker','.');
    %*********
    ylabel('\delta')
end


%% 3b) Same plot as abouve, but grouped by response

load('PefSimResults2.mat')  %-Retrieves object called expmt with the results.
expmt(expmt.scenario==4,:) = [];  %-Discard the "Case 4" runs we won't use here.

hspc = [1.5 3 3.5];     %-Horizontal spacings for mysubplot()
vspc = [0.5 1.25];         %-Vertical spacings for mysubplot()
fig3 = myfigure(3,[21 8]);
clf
models = {'SZO','SPM','CZO'};
% lnmrk = {'o','.','+'};
% lnwid = [0.5, 1 , 1.5];
% lnsty = {':','-','--'};
%lnmrk = {'0','s','x'};
% lnwid = [1, 1.5 , 1]*0.75;
% lnsty = {'-','-','-'};
lnmrk = {'.','.','.'};
lnwid = [0.5, 1.5 , 2];
lnsty = {'-','-.','-'};
mksz = [8 8 8];

lams = unique(expmt.lambda);

% Create the axes
ax1 = mysubplot(1, 2, 1, hspc, vspc);
    axis([0 lams(end) 0.5 1])
    xlabel('\lambda')
    ylabel('F_{\fontname{cmsy10}E}')
    hold on
ax2 = mysubplot(1, 2, 2, hspc, vspc);
    axis([0 lams(end) 0 1])
    xlabel('\lambda')
    ylabel('D_0')
    hold on

% Plot the lines into the axes
for mm = 1:3
    ix = strcmp(expmt.model,models{mm});
    rhos = expmt.rhohat(ix);
    absdists = expmt.absdist(ix);
    plot(ax1,lams,rhos,'Color','black','Marker',lnmrk{mm},'LineWidth',lnwid(mm), ...
         'LineStyle',lnsty{mm},'MarkerSize',mksz(mm))
    plot(ax2,lams,absdists,'Color','black','Marker',lnmrk{mm},'LineWidth',lnwid(mm), ...
         'LineStyle',lnsty{mm},'MarkerSize',mksz(mm));
end

% Add legend
leg = legend('S_{0,1}','S_{-1/2,1/2}','C_{0,1}');
leg.Position = [0.8621 0.7086 0.1211 0.2307];

%print(fig3,'-r600','-dpdf','-opengl','Sec4-AvgEffOfLambda')
%print(fig3,'-dpdf','-painters','Sec4-AvgEffOfLambda2')




%% 4) Plot samples drawn from the three models with different lambdas
% See my notebook 2016-05-03 for layout diagram

fig4 = myfigure(4,[20 11]);
clf
vspc = [0.05 0.01 0.05 0.05];
hspc = [0.1 0.005 0.05];
reps = 2;  %-number of replicates to draw for each case.
nlam = 11;
lam_max = 2.5;

lams = linspace(0,lam_max,nlam);
vgap = (1 - vspc(1) - 2*vspc(3) - vspc(4))/3;  
hgap = (1 - hspc(1) - (nlam-1)*hspc(2) - hspc(3))/nlam;
ht = (vgap - (reps-1)*vspc(2))/reps;
wd = hgap;

% %--- Generate all the needed perfect samples ---
% TMP = cell(reps,nlam,3);   %-Rows = reps, Cols = samples, pages = models
% for i = 1:nlam
%     SZO.Lambda = lams(i);
%     SPM.Lambda = lams(i);
%     CZO.Lambda = lams(i);
%     for j = 1:reps
%         TMP{j,i,1} = reshape(SZO.PerfectSample2,d,d);
%         TMP{j,i,2} = reshape(SPM.PerfectSample2,d,d) + 1/2;
%         TMP{j,i,3} = reshape(CZO.PerfectSample2,d,d);
%     end
% end
% %------------------------------------------------

%--- Make the plot ---
for k = 1:3
    %-LB is the (left,bottom) coords of the array of plots for model k.
    LB = [hspc(1), vspc(4) + (3-k)*(vgap+vspc(3))];
    
    for i = 1:nlam
        
        for j = 1:reps
            L = LB(1) + (i-1)*(wd+hspc(2));
            B = LB(2) + (j-1)*(ht+vspc(2));
            pos = [L, B, wd, ht];
            ax = axes('Position',pos);
            axis([0 1 0 1])
            imagesc([0 1],[0 1],TMP{j,i,k})
            caxis([0 1])
            colormap(gray)
            ax.YDir = 'normal';
            ax.Box = 'on';
            ax.TickLength = [0 0];
            ax.XTickLabel = '';
            ax.YTickLabel = '';
        end
        
        if k==1  %-Add lambda values annotations
            annotation('textbox',[L, B+ht+0.01, wd, 0.05], ...
                       'VerticalAlignment','bottom', ...
                       'HorizontalAlignment','center',...
                       'Margin', 0, ...
                       'LineStyle','none', ...
                       'String', ['\lambda = ' num2str(lams(i))], ...
                       'FontSize', 9)
        end
        
    end
end

%Add model names as annotations
strs = cell(3,1);
posns = cell(3,1);
strs{1,1} = {'S_{0,1}', 'model'};
strs{2,1} = {'S_{-1/2,1/2}', 'model'};
strs{3,1} = {'C_{0,1}', 'model'};
posns{1,1} = [0, 1-vspc(1)-vgap/2, hspc(1), 0];
posns{2,1} = [0, vspc(4)+vspc(3)+1.5*vgap, hspc(1), 0];
posns{3,1} = [0, vspc(4)+0.5*vgap, hspc(1), 0];
for i = 1:3
    annotation('textbox',posns{i,1}, ...
               'VerticalAlignment','middle', ...
               'HorizontalAlignment','center',...
               'Margin', 0, ...
               'LineStyle','none', ...
               'String', strs{i,1}, ...
               'FontSize', 12)
end 
%---------------------

%Save the figure
%print(fig4,'-r600','-dpdf','-opengl','Sec4-Realizations')


%% 4b) Estimate the marginal probability maps for the three models 
% Get the probability maps at different lambdas for each model with beta = [-2 2 2]';

nlam = 11;
lam_max = 2.5;
lams = linspace(0,lam_max,nlam);
nsamp = 500;   %-number of samples to draw for each case.
stor = zeros(SPM.N,length(lams),3);  %-rows = variables, cols = lambdas, pages = models
tmp = zeros(SPM.N,nsamp);  %-for holding sampling results temporarily.

tic
for i = 1:length(lams)
    
    % Do the SPM model
    SPM.Lambda = lams(i);
    SPM.Beta = [-2 2 2]';
    for s = 1:nsamp
        clc
        disp(['Starting SPM draw ' num2str(s) ' with lam = ' num2str(lams(i))])
        tmp(:,s) = SPM.PerfectSample2==SPM.Coding(2);
    end
    stor(:,i,1) = mean(tmp,2);
    
    % Do the SZO model
    SZO.Lambda = lams(i);
    SZO.Beta = [-2 2 2]';
    for s = 1:nsamp
        clc
        disp(['Starting SZO draw ' num2str(s) ' with lam = ' num2str(lams(i))])
        tmp(:,s) = SZO.PerfectSample2==SZO.Coding(2);
    end
    stor(:,i,2) = mean(tmp,2);
    
    % Do the CZO model
    CZO.Lambda = lams(i);
    CZO.Beta = [-2 2 2]';
    for s = 1:nsamp
        clc
        disp(['Starting CZO draw ' num2str(s) ' with lam = ' num2str(lams(i))])
        tmp(:,s) = CZO.PerfectSample2==CZO.Coding(2);
    end
    stor(:,i,3) = mean(tmp,2);
        
end
toc

% Save images of the probability maps
for i = 1:length(lams)
    str1 = ['SPMprobmap-222-lam' num2str(i) '.png'];
    str2 = ['SZOprobmap-222-lam' num2str(i) '.png'];
    str3 = ['CZOprobmap-222-lam' num2str(i) '.png'];
    imwrite(flipud(reshape(stor(:,i,1),d,d)), str1)
    imwrite(flipud(reshape(stor(:,i,2),d,d)), str2)
    imwrite(flipud(reshape(stor(:,i,3),d,d)), str3)
end



%% 4c) Plot samples AND marginals for the three models with different lambdas
% Follow the code of 4a) but put the marginals in the 1st row.

%needed for plotting from file...
% tst = imread('CZOprobmap-222-lam5.png',1)
% imagesc(double(flipud(tst))/255)
% set(gca,'YDir','normal')
% caxis([0 1])
%***************

%fig4 = myfigure(4,[20 14]);
fig4 = myfigure(4,[18 14]);
clf
vspc = [0.05 0.005 0.05 0.05];
hspc = [0.15 0.0035 0.05];
reps = 3;  %-number of replicates to draw for each case.
nlam = 11;
lam_max = 2.5;

lams = linspace(0,lam_max,nlam);
vgap = (1 - vspc(1) - 2*vspc(3) - vspc(4))/3;  
hgap = (1 - hspc(1) - (nlam-1)*hspc(2) - hspc(3))/nlam;
ht = (vgap - (reps-1)*vspc(2))/reps;
wd = hgap;

% %--- Generate all the needed perfect samples ---
% TMP = cell(reps,nlam,3);   %-Rows = reps, Cols = samples, pages = models
% for i = 1:nlam
%     SZO.Lambda = lams(i);
%     SPM.Lambda = lams(i);
%     CZO.Lambda = lams(i);
%     
%     nm1 = ['SZOprobmap-222-lam', num2str(i), '.png'];
%     nm2 = ['SPMprobmap-222-lam', num2str(i), '.png'];
%     nm3 = ['CZOprobmap-222-lam', num2str(i), '.png'];
%     tmp = imread(nm1);
%     TMP{1,i,1} = double(flipud(tmp))/255;
%     tmp = imread(nm2);
%     TMP{1,i,2} = double(flipud(tmp))/255;
%     tmp = imread(nm3);
%     TMP{1,i,3} = double(flipud(tmp))/255;
%     
%     for j = 2:reps
%         TMP{j,i,1} = reshape(SZO.PerfectSample2,d,d);
%         TMP{j,i,2} = reshape(SPM.PerfectSample2,d,d) + 1/2;
%         TMP{j,i,3} = reshape(CZO.PerfectSample2,d,d);
%     end
% end
% %2017-03-13: rearrange model order to suit desired plot order for figure (changed the
% %plotting code to reflect this change below)
% TMPP = TMP;
% TMPP(:,:,1) = TMP(:,:,1);
% TMPP(:,:,2) = TMP(:,:,3);
% TMPP(:,:,3) = TMP(:,:,2);
% TMP = TMPP;
% %------------------------------------------------

%--- Make the plot ---
for k = 1:3
    %-LB is the (left,bottom) coords of the array of plots for model k.
    LB = [hspc(1), vspc(4) + (3-k)*(vgap+vspc(3))];
    
    for i = 1:nlam
        
        for j = 1:reps
            L = LB(1) + (i-1)*(wd+hspc(2));
            B = LB(2) + (j-1)*(ht+vspc(2));
            pos = [L, B, wd, ht];
            ax = axes('Position',pos);
            axis([0 1 0 1])
            imagesc([0 1],[0 1],TMP{j,i,k})
            caxis([0 1])
            colormap(gray)
            ax.YDir = 'normal';
            ax.Box = 'on';
            ax.TickLength = [0 0];
            ax.XTickLabel = '';
            ax.YTickLabel = '';
        end
        
        if k==1  %-Add lambda values annotations
            annotation('textbox',[L, B+ht+0.01, wd, 0.05], ...
                       'VerticalAlignment','bottom', ...
                       'HorizontalAlignment','center',...
                       'Margin', 0, ...
                       'LineStyle','none', ...
                       'String', ['\lambda = ' num2str(lams(i))], ...
                       'FontSize', 9)
        end
        
    end
end

%Add model names as annotations
strs = cell(3,1);
posns = cell(3,1);
strs{1,1} = {'S_{0,1}', 'model'};
strs{2,1} = {'C_{0,1}', 'model'};
strs{3,1} = {'S_{-1/2,1/2}', 'model'};
posns{1,1} = [0, 1-vspc(1)-vgap/2, hspc(1), 0];
posns{2,1} = [0, vspc(4)+vspc(3)+1.5*vgap, hspc(1), 0];
posns{3,1} = [0, vspc(4)+0.5*vgap, hspc(1), 0];
for i = 1:3
    pos = posns{i,1};
    annotation('textbox',pos, ...
               'VerticalAlignment','middle', ...
               'HorizontalAlignment','left',...
               'Margin', 2, ...
               'LineStyle','none', ...
               'String', strs{i,1}, ...
               'FontSize', 12)
    pos1 = pos + [0 ht+vspc(2) 0 0];
    annotation('textbox',pos1, ...
               'VerticalAlignment','middle', ...
               'HorizontalAlignment','right',...
               'Margin', 2, ...
               'LineStyle','none', ...
               'String', 'draw 1 ', ...
               'FontSize', 9)
    pos2 = pos;
    annotation('textbox',pos2, ...
               'VerticalAlignment','middle', ...
               'HorizontalAlignment','right',...
               'Margin', 2, ...
               'LineStyle','none', ...
               'String', 'draw 2 ', ...
               'FontSize', 9)
    pos3 = pos - [0 ht+vspc(2) 0 0];
    annotation('textbox',pos3, ...
               'VerticalAlignment','middle', ...
               'HorizontalAlignment','right',...
               'Margin', 2, ...
               'LineStyle','none', ...
               'String', 'map ', ...
               'FontSize', 9)
end 
%---------------------

%Save the figure
%print(fig4,'-r600','-dpdf','-opengl','Sec4-Realizations')
%print(fig4, 'Sec4-RealizationsB', '-dtiff', '-opengl', '-r600')








%% 5) Get results and samples for a table of "nearest models"
% Use the standard (-1/2,1/2) model as reference and find the SPM and CZO models
% closest to them in "least squares unary parameter distance" sense.  Then draw a few
% perfect samples and save them as image files.

lams = 0:0.1:2.5; %-Lambda values to try.
params = zeros(length(lams),3,2);  %-Page 1 is for SZO, page 2 is for CZO
dists = zeros(length(lams),2);  %-To hold norms

for i = 1:length(lams)
    
    % Convert to nearest SZO model
    disp(['lambda = ' num2str(lams(i)) ' SZO'])
    SPM.Lambda = lams(i);
    tf = SPM.Transform([0 1],false);  %-transform to SZO (alphas)
    ls_beta = SPM.X\tf.alpha;   %-The regression coefficients closest to SPM's in 
                                % least-squares sense.
    params(i,:,1) = ls_beta;
    dists(i,1) = norm(ls_beta - SPM.Beta,1);
    
%     % Get draws from nearest SZO model
%     SZO.Lambda = lams(i);
%     SZO.Beta = ls_beta;
%     for j = 1:3
%         tmp = reshape(SZO.PerfectSample2,d,d);
%         str = ['SZOdraw-lam' num2str(lams(i)) '-' num2str(j) '.png'];
%         imwrite(flipud(tmp),str)
%     end
    
    % Convert to nearest CZO model
    disp(['lambda = ' num2str(lams(i)) ' CZO'])
    tf = SPM.Transform([0 1],true);  %-transform to CZO (alphas)
    ls_beta = SPM.X\tf.alpha;
    params(i,:,2) = ls_beta;
    dists(i,2) = norm(ls_beta - SPM.Beta,1);

%     % Get draws from nearest CZO model
%     CZO.Lambda = lams(i);
%     CZO.Beta = ls_beta;
%     for j = 1:3
%         tmp = reshape(CZO.PerfectSample2,d,d);
%         str = ['CZOdraw-lam' num2str(lams(i)) '-' num2str(j) '.png'];
%         imwrite(flipud(tmp),str)
%     end

end

%display results
disp(params)
disp(dists)

%---saved results---
% params(:,:,1) =
%    -2.0000    2.0000    2.0000
%    -2.9667    2.0000    2.0000
%    -3.9333    2.0000    2.0000
%    -4.9000    2.0000    2.0000
%    -5.8667    2.0000    2.0000
%    -6.8333    2.0000    2.0000
% params(:,:,2) =
%    -2.0000    2.0000    2.0000
%    -3.1707    3.1707    3.1707
%    -5.0882    5.0882    5.0882
%    -7.1973    7.1973    7.1973
%    -9.1838    9.1838    9.1838
%   -11.0929   11.0929   11.0929
% dists = 
%     0.0000    0.0000
%     0.9667    3.5120
%     1.9333    9.2646
%     2.9000   15.5918
%     3.8667   21.5513
%     4.8333   27.2786
%-------------------


%% 5b) For CZO case, get nearest models properly

% Put previous params into 3rd plane of params
params(:,:,3) = params(:,:,2);

% Make needed variables
X = SPM.X;
A = G.adjacency;
one = ones(SPM.N,1);

% mymu calculates the mu vector for given regression coeffs
mymu = @(theta) exp(X*theta)./(1+exp(X*theta));

%nlam = 20;
nlam = length(lams);
lams = linspace(0,2.5,nlam);
allbest = zeros(nlam,3);
allOF = zeros(nlam,1);
allflags = zeros(nlam,1);
for i = 1:length(lams)
    % objective function
    lambda = lams(i);
    OF = @(theta) norm(X*(theta-beta) + lambda/2*A*(one-2*mymu(theta)) ,2);

    % optimize and save
    %*** starting from previous solution gives worse results!!!***
    %*** starting from [-2, 2, 2] gives better results, but not smooth as f(lambda)
    %*** So try with random restarts****
    keepOF = 1e6;
    keepbest = zeros(3,1);
    keepflag = 0;
    for j = 1:25
        beta0 = beta + normrnd(0,3,3,1);
        [best, OFval, flag] = fminunc(OF,beta0);
        if OFval < keepOF
            keepOF = OFval;
            keepbest = best;
            keepflag = flag;
        end
    end
    params(i,:,2) = keepbest';

    allbest(i,:) = keepbest';
    allOF(i) = keepOF;
    allflags(i) = keepflag;
end
    
disp(params)
[lams' allbest allOF allflags]

% RESULTS
% [lams' allbest allOF allflags]
% ans =
%          0   -2.0000    2.0000    2.0000    0.0000    2.0000
%     0.1000   -2.1761    2.1761    2.1761    0.3747    2.0000
%     0.2000   -2.3743    2.3743    2.3743    0.8356    1.0000
%     0.3000   -2.5949    2.5949    2.5949    1.4042    1.0000
%     0.4000   -2.8368    2.8368    2.8368    2.1027    1.0000
%     0.5000   -3.0977    3.0977    3.0977    2.9517    1.0000
%     0.6000   -3.3740    3.3740    3.3740    3.9667    1.0000
%     0.7000   -3.6619    3.6619    3.6619    5.1563    1.0000
%     0.8000   -3.9582    3.9582    3.9582    6.5217    1.0000
%     0.9000   -4.2602    4.2602    4.2602    8.0580    1.0000
%     1.0000   -4.5660    4.5660    4.5660    9.7560    1.0000
%     1.1000   -4.8743    4.8743    4.8743   11.6041    1.0000
%     1.2000   -5.1846    5.1846    5.1846   13.5894    1.0000
%     1.3000   -5.4965    5.4965    5.4965   15.6994    1.0000
%     1.4000    0.0805    2.2350    2.2350   15.7726    1.0000
%     1.5000    0.2765    2.2567    2.2567   14.7039    1.0000
%     1.6000    0.4557    2.2905    2.2905   13.7124    1.0000
%     1.7000    0.6334    2.3247    2.3247   12.8765    1.0000
%     1.8000    0.8205    2.3511    2.3511   12.2563    1.0000
%     1.9000    1.0250    2.3635    2.3635   11.8852    1.0000
%     2.0000    2.0746   -2.0746   -2.0746   11.0352    1.0000
%     2.1000    2.0082   -2.0082   -2.0082    9.7284    1.0000
%     2.2000    1.9057   -1.9057   -1.9057    8.6766    2.0000
%     2.3000    1.7837   -1.7837   -1.7837    7.8624    2.0000
%     2.4000    1.6590   -1.6590   -1.6590    7.2410    1.0000
%     2.5000    1.5417   -1.5417   -1.5417    6.7634    1.0000

%% 5c) For CZO case, get nearest models properly and PLOT (more lambda vals)

% Make needed variables
X = SPM.X;
A = G.adjacency;
one = ones(SPM.N,1);

% mymu calculates the mu vector for given regression coeffs
mymu = @(theta) exp(X*theta)./(1+exp(X*theta));

nlam = 25;
lam_max = 20;
lams = linspace(0,lam_max,nlam);
allbest = zeros(nlam,3);
allOF = zeros(nlam,1);
allflags = zeros(nlam,1);
for i = 1:length(lams)
    % objective function
    lambda = lams(i);
    OF = @(theta) norm(X*(theta-beta) + lambda/2*A*(one-2*mymu(theta)) ,2);

    % optimize and save
    %*** starting from previous solution gives worse results!!!***
    %*** starting from [-2, 2, 2] gives better results, but not smooth as f(lambda)
    %*** So try with random restarts****
    keepOF = 1e6;
    keepbest = zeros(3,1);
    keepflag = 0;
    for j = 1:25
        beta0 = beta + normrnd(0,3,3,1);
        [best, OFval, flag] = fminunc(OF,beta0);
        if OFval < keepOF
            keepOF = OFval;
            keepbest = best;
            keepflag = flag;
        end
    end
    params(i,:,2) = keepbest';

    allbest(i,:) = keepbest';
    allOF(i) = keepOF;
    allflags(i) = keepflag;
end
    
disp(params)
[lams' allbest allOF allflags]

% quick and dirty figure for now, not sure if I need/want it for the paper...
fig5 = figure(5);
plot(lams,allOF)
line(lams,allbest(:,1),'Color','red')
line(lams,allbest(:,2),'Color','blue')


%% 5d) Plot least-squared distances for the CZO and SZO models (for publication)

% Make needed variables
X = SPM.X;
A = G.adjacency;
one = ones(SPM.N,1);
beta = [-2 2 2]';

% mymu calculates the mu vector for given regression coeffs
mymu = @(theta) exp(X*theta)./(1+exp(X*theta));

nlam = 150;
lam_max = 2.5;
lams = linspace(0,lam_max,nlam);
allbest = zeros(nlam,3,2);  %-Plane 1: CZO; plane 2: SZO
allOF = zeros(nlam,2);   %-col 1: CZO; col 2: SZO
allflags = zeros(nlam,1);
for i = 1:length(lams)
    % objective function (CZO)
    lambda = lams(i);
    OF = @(theta) norm(X*(theta-beta) + lambda/2*A*(one-2*mymu(theta)) ,2);

    % optimize and save
    %*** starting from previous solution gives worse results!!!***
    %*** starting from [-2, 2, 2] gives better results, but not smooth as f(lambda)
    %*** So try with random restarts****
    keepOF = 1e6;
    keepbest = zeros(3,1);
    keepflag = 0;
    for j = 1:25
        beta0 = beta + normrnd(0,3,3,1);
        [best, OFval, flag] = fminunc(OF,beta0);
        if OFval < keepOF
            keepOF = OFval;
            keepbest = best;
            keepflag = flag;
        end
    end

    allbest(i,:,1) = keepbest';
    allOF(i,1) = keepOF;
    allflags(i) = keepflag;
    
    % Get SZO nearest model
    T = X'*X;
    H = T\X'; 
    ls_beta = beta - lambda/2*H*A*one;
    allbest(i,:,2) = (ls_beta)';
    allOF(i,2) = norm(X*beta - lambda/2*A*one - X*ls_beta ,2);
    
end

%% 5d) continued--make the figure

fig6 = myfigure(6,[18,8]);
clf
plot(lams,allOF(:,1),'k-','LineWidth',2)
line(lams,allOF(:,2),'Color','black')
legend('C_{0,1}','S_{0,1}','Location','northwest')
xlabel('\lambda')
ylabel('Least-squares distance')

%print(fig6,'-r600','-dpdf','-opengl','Sec4-NearestOFs')
%print(fig6,'-dpdf','-painters','Sec4-NearestOFs')

%% 6) For the "nearest models" cases, estimate the marginal probability maps.
% requires the previous section(s) to be run so we have the params array, lams, etc.
% need to run the SPM model here too, get the probability maps.
% *** RUN CHUNKS 5 and 5b BEFORE RUNNING THIS ***

nsamp = 500;   %-number of samples to draw for each case.
stor = zeros(SPM.N,length(lams),3);  %-rows = variables, cols = lambdas, pages = models
tmp = zeros(SPM.N,nsamp);  %-for holding sampling results temporarily.

tic
for i = 1:length(lams)
    
    % Do the SPM model
    SPM.Lambda = lams(i);
    SPM.Beta = [-2 2 2]';
    for s = 1:nsamp
        clc
        disp(['Starting SPM draw ' num2str(s) ' with lam = ' num2str(lams(i))])
        tmp(:,s) = SPM.PerfectSample2==SPM.Coding(2);
    end
    stor(:,i,1) = mean(tmp,2);
    
    % Do the SZO model
    SZO.Lambda = lams(i);
    SZO.Beta = params(i,:,1)';
    for s = 1:nsamp
        clc
        disp(['Starting SZO draw ' num2str(s) ' with lam = ' num2str(lams(i))])
        tmp(:,s) = SZO.PerfectSample2==SZO.Coding(2);
    end
    stor(:,i,2) = mean(tmp,2);
    
    % Do the CZO model
    CZO.Lambda = lams(i);
    CZO.Beta = params(i,:,2)';
    for s = 1:nsamp
        clc
        disp(['Starting CZO draw ' num2str(s) ' with lam = ' num2str(lams(i))])
        tmp(:,s) = CZO.PerfectSample2==CZO.Coding(2);
    end
    stor(:,i,3) = mean(tmp,2);
        
end
toc

% Save images of the probability maps
for i = 1:length(lams)
    str1 = ['SPMprobmap-lam' num2str(i) '.png'];
    str2 = ['SZOprobmap-lam' num2str(i) '.png'];
    str3 = ['CZOprobmap-lam' num2str(i) '.png'];
    imwrite(flipud(reshape(stor(:,i,1),d,d)), str1)
    imwrite(flipud(reshape(stor(:,i,2),d,d)), str2)
    imwrite(flipud(reshape(stor(:,i,3),d,d)), str3)
end

















