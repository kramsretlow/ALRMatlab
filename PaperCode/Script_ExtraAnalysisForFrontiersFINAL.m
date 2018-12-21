%% n=2 example: showing param effects on P(Z1=h) and P(Z1=Z2)

%% Input values

% User inputs
G = AutoLogistic.MakeGraph(1,2,1);  %Graph to use.
ix1 = 1;  %-Index of Z1 in the graph.
ix2 = 2;  %-Index of Z2 in the graph.
ngrid1d = 250;  %-Number of gridpoints to use for 1D plots
ngrid2d = 50;  %-Number of gripdpoints to use per axis for 2D plots

% Derived objects
n = G.numnodes;
                                     %-assign_alpha is used for setting unary params.
assign_alpha = @(a1, a2) [a2*ones(ix1-1,1); a1; a2*ones(n-ix1,1)];

%% Set up AL model objects

% Here we'll let beta=1 and control the unary parameters by setting the n-vector X to
% the values we want (using assign_alpha). Initialize Lambda to zero and change it's
% value as desired later on.  Create an ALRsimple object for each of four models:
% CZO, CPM, SZO, SPM.

CZO = ALRsimple();
CZO.Coding = [0 1];
CZO.Centered = true;
CZO.Beta = 1;
CZO.Lambda = 0;
CZO.Graph = G;

CPM = CZO; 
CPM.Coding = [-1 1];

SPM = CPM;
SPM.Centered = false;

SZO = SPM;
SZO.Coding = [0 1];


%% Try plotting P(Z1=h) as fcn of p1, p2 for fixed lambda
%p1 and p2 are the endogenous probabilities (antilogit of the alphas)

%***inputs***
pmin = 0.01;      %-minimum endo prob to consider
pmax = 0.99;      %-maximum endo prob to consider
variant = CZO;    %-The model variant to use.
hspc = [1.3 0.4 0.4 0.4 0.4];  %-Horizontal spacings for mysubplot()
vspc = [0.7 1.3];                %-Vertical spacings for mysubplot()
checkp1 = 0.5;   %-For choosing max lambda: look at P(Z1=h) at p1=checkp1 & p2=checkp2
checkp2 = 0.95;
checkprob = 0.9;   %-For picking max lambda.
%proposed way to pick max lambda is to find the minimum lambda s.t. when p1=0.5 and
%p2=0.95, P(Z1=hi) = 0.9.  That is, set Z1's endo prob to 1/2, make Z2 very likely to
% take value hi, and increase lambda until the neighbor effect pulls Z1's prob up to
% 0.9
%************

lo = variant.Coding(1);
hi = variant.Coding(2);
if variant.Centered
    mufun = @(a) (lo*exp(lo*a)+hi*exp(hi*a))/(exp(lo*a)+exp(hi*a));
else
    mufun = @(a) 0;
end
prng = linspace(pmin,pmax,ngrid2d);
arng = log(prng./(1-prng))/(hi-lo);

%find lambda_max.  Will use four lambdas equally spaced in [0,lambda_max]
gfun = @(z,a1,a2,lam) exp(a1*z(1)+a2*z(2)+lam*(z(1)*z(2)-mufun(a2)*z(1)-mufun(a1)*z(2)));
partition = @(a1,a2,lam) gfun([lo lo],a1,a2,lam) + gfun([lo hi],a1,a2,lam) + gfun([hi lo],a1,a2,lam) + gfun([hi hi],a1,a2,lam);
pmf = @(z,a1,a2,lam) gfun(z,a1,a2,lam) / partition(a1,a2,lam);
checka1 = log(checkp1/(1-checkp1))/(hi-lo);
checka2 = log(checkp2/(1-checkp2))/(hi-lo);
PZ1eq1 = @(a1,a2,lam) pmf([hi lo],a1,a2,lam)+pmf([hi hi],a1,a2,lam);
tstfun = @(lam) (abs(PZ1eq1(checka1,checka2,lam)-0.5) - (checkprob-0.5)).^2;
lambda_max = fminbnd(tstfun,0,20);
disp(lambda_max)

% Make the graph
lams = linspace(0,lambda_max,4);
[A1, A2] = meshgrid(arng,arng);
[P1, P2] = meshgrid(prng,prng);
P1h = zeros(ngrid2d);

fig = myfigure(5,[18 5.75]);
clf

for k = 1:4
    for i = 1:ngrid2d
        for j = 1:ngrid2d
            P1h(i,j) = PZ1eq1(A1(i,j),A2(i,j),lams(k));
        end
    end
    disp(k)
    
    ax = mysubplot(1,4,k,hspc,vspc);
    contourf(P1,P2,P1h,[0 0.05 0.1:0.1:0.9 0.95 1])
    colormap(gray)
    caxis([0 1])
    line([0 1],[0.5 0.5],'Color','black','LineStyle','--')
    line([0.5 0.5],[0 1],'Color','black','LineStyle','--')
    axis square
    xlim([0 1])
    ylim([0 1])
    if k==1
        ylabel('p_2')
    else
        ax.YAxis.TickLabels = '';
    end
    xlabel('p_1')
    titstr = ['\lambda = ', num2str(lams(k),3)];
    title(titstr)
    %try overlaying thicker contour at probability 0.5
    hold on
    contour(P1,P2,P1h,[0.5 0.5],'k-','LineWidth',3)

end

%print(fig, 'PZ1eq1vsP1P2-CZO', '-dtiff', '-opengl', '-r600')


%% Plot P(Z1=Z2) as a function of lambda, for different p1, p2 combos

%***inputs***
p2rng = [0.05 0.25 0.5 0.75 0.95];  %-Will plot lines on each plot for each value in p2rng
p1rng = [0.05 0.25 0.5 0.75 0.95];    %-Make one plot for each p1 in p1rng
variant = CZO;    %-The model variant to use.
figno = 1;   %-figure number to use.
hgaps = [1.5 0.4 0.4];  %-Left, interior, and right spacings for hspace in mysubplot()
vspc = [0.7 1.3];                %-Vertical spacings for mysubplot()
logitrng = [-2.5 8];
%************

lo = variant.Coding(1);
hi = variant.Coding(2);
if variant.Centered
    mufun = @(a) (lo*exp(lo*a)+hi*exp(hi*a))/(exp(lo*a)+exp(hi*a));
else
    mufun = @(a) 0;
end
a2rng = log(p2rng./(1-p2rng))/(hi-lo);
a1rng = log(p1rng./(1-p1rng))/(hi-lo);

% Set the maximum lambda value.  Choose it such that logit(P(Z1==Z2)) = logitrng(2)
% when p1 = p2 = 0.5
gfun = @(z,a1,a2,lam) exp(a1*z(1)+a2*z(2)+lam*(z(1)*z(2)-mufun(a2)*z(1)-mufun(a1)*z(2)));
partition = @(a1,a2,lam) gfun([lo lo],a1,a2,lam) + gfun([lo hi],a1,a2,lam) + gfun([hi lo],a1,a2,lam) + gfun([hi hi],a1,a2,lam);
pmf = @(z,a1,a2,lam) gfun(z,a1,a2,lam) / partition(a1,a2,lam);
PZ1eqZ2 = @(a1,a2,lam) pmf([lo lo],a1,a2,lam)+pmf([hi hi],a1,a2,lam);
tstfun = @(lam) (log(PZ1eqZ2(0,0,lam)/(1-PZ1eqZ2(0,0,lam))) - logitrng(2)).^2;
lambda_max = fminbnd(tstfun,0,20);
disp(lambda_max)

%logit function
logit = @(p) log(p./(1-p));

%Make the graph
lams = linspace(0,lambda_max,ngrid1d);
nplot = length(p1rng);
hspc = [hgaps(1) hgaps(2)*ones(1,nplot-1) hgaps(3)];  %-Horizontal spacings for mysubplot()
fcnval = zeros(ngrid1d,1);

fig = myfigure(figno,[18 5]);
clf

for k = 1:nplot
    a1 = a1rng(k);

    ax = mysubplot(1,nplot,k,hspc,vspc);
    axis([0 lams(end) logitrng]);
    
    titstr = ['p_1 = ', num2str(p1rng(k),3)];
    title(titstr)
    if k==1
        ylabel('log odds of Z_1 = Z_2')
    else
        ax.YAxis.TickLabels = '';
    end
    xlabel('\lambda')
    grid on
    %axis square

    for i = 1:length(a2rng)
        for j = 1:ngrid1d
            fcnval(j) = logit(PZ1eqZ2(a1rng(k),a2rng(i),lams(j)));
        end
        
        ln = line(lams,fcnval,'Color','black');
        ln.LineWidth = 0.5 + 2*(i-1)/(length(a2rng)-1);
    end

end

%print(fig, 'PZ1eqZ2vsP1P2-SPM', '-dtiff', '-opengl', '-r600')



%% 30x30 spatial example: tabulate edge match probs vs. edge endo probs
% Generate perfect samples from SPM, SZO, and CZO models and count the proportion of
% edges "matched" in each of two regions of the graph:
% Region 1: nodes with endo probability <= 1/3 or >= 2/3
% Region 2: nodes with intermediate endo prob, in (1/3, 2/3)

% Get the X data and the endogenous probabilities
clear all
load('PerfSimXs.mat')  %X3, endo3 contain the case used in the paper.
X = X3;
endo = Endo3;
d = 30;
beta = [-2 2 2]';  %-Regression coefficients

%Set things up
ngen = 500;  %-Number of samples to generate per model
lams = linspace(0,1,9);
nlams = length(lams);
mods = {'CZO','SPM'};
nmods = length(mods);

%-Use mod as a temporary model object for doing runs.
mod = ALRsimple;   
G = AutoLogistic.MakeGraph(d,d,1);
mod.Graph = G;
mod.Beta = beta;
mod.X = X;
nedge = G.numedges;

%Create groups of nodes and also group edges (edge must have both nodes in the same
%group).
ix1 = find(endo<=1/3 | endo>=2/3);
ix2 = find(endo>1/3 & endo<2/3);
edg1 = zeros(nedge,1);
edg2 = zeros(nedge,1);
for i = 1:nedge
    thisedge = G.Edges.EndNodes(i,:);
    edg1(i) = all(ismember(thisedge, ix1));
    edg2(i) = all(ismember(thisedge, ix2));
end
edg1 = find(edg1);
edg2 = find(edg2);

% A function for counting matched edges among nodes in a group
nmatch = @(samp,edg) sum(samp(G.Edges.EndNodes(edg,1))==samp(G.Edges.EndNodes(edg,2)));

%Create matrices for holding output.  Match1 and Match2 are nlams-by-nmods, and hold
%the average of (proportion matched edges)
Match1 = zeros(nlams,nmods);
Match2 = zeros(nlams,nmods);

%Loop through models and lambdas and generate samples to get matched frac.
tic
for j = 1:nmods
    switch mods{j}
        case 'SZO'
            mod.Coding = [0 1];
            mod.Centered = false;
        case 'SPM'
            mod.Coding = [-1/2 1/2];
            mod.Centered = false;
        case 'CZO'
            mod.Coding = [0 1];
            mod.Centered = true;
    end
    for i = 1:nlams
        mod.Lambda = lams(i);
        for g = 1:ngen
            tmp = mod.PerfectSample2;
            Match1(i,j) = Match1(i,j) + nmatch(tmp,edg1);
            Match2(i,j) = Match2(i,j) + nmatch(tmp,edg2);
            disp([i j g])
        end
    end
end
Match1 = Match1/length(edg1)/ngen;
Match2 = Match2/length(edg2)/ngen;
toc

% Save the MF array to file
%save('EdgeMatchDataFinal.mat','Match1','Match2','lams','mods','endo')


%% Plot the results of the above

load('EdgeMatchDataFinal.mat')

fig11 = figure(11);
clf;
fig11.Units = 'centimeters';
fig11.Position = [1 1 8.5 8];
fig11.PaperPositionMode = 'auto';  %for printing at screen size
fig11.Color = [1 1 1];
ax = axes();
ax.Position = [0.1667    0.1442    0.7014    0.7808];
ax.YLim = [0.5 0.9];
ax.XLim = [0 1];
ax.YTick = 0.5:0.1:1;
%ax.XTick = vploc;
%ax.XTickLabel = repmat(lams,[4 1]);
%ax.XTickLabelRotation = 90;
%ax.YGrid = 'on';
ylabel('Estimated concordance probability')
xlabel('\lambda')
hold on

% plot(lams,Match1(:,2),'ok-','MarkerSize',4,'MarkerFaceColor',[0 0 0])
% plot(lams,Match2(:,2),'--ok','MarkerSize',4,'MarkerFaceColor',[0 0 0])
% plot(lams,Match1(:,1),'dk-','MarkerSize',4)
% plot(lams,Match2(:,1),'dk--','MarkerSize',4)

% Match1 = logit(Match1);
% Match2 = logit(Match2);
% ax.YLim = [0 2.1];
% ax.YTick = 0:0.5:2;

plot(lams,Match1(:,2),'ok-','MarkerSize',4,'MarkerFaceColor',[0 0 0],'LineWidth',1.5)
plot(lams,Match2(:,2),'-ok','MarkerSize',4,'MarkerFaceColor',[0 0 0],'LineWidth',1.5)
plot(lams,Match1(:,1),'dk-','MarkerSize',4)
plot(lams,Match2(:,1),'dk-','MarkerSize',4)

% legend('S_{-1/2,1/2} region 1','S_{-1/2,1/2} region 2', ...
%        'C_{0,1} region 1','C_{0,1} region 2','Location','NorthWest')

gap = 0.035;
text(lams(end)+gap,Match1(end,2),'Sy1')
text(lams(end)+gap,Match2(end,2),'Sy2')
text(lams(end)+gap,Match1(end,1),'C1')
text(lams(end)+gap,Match2(end,1),'C2')
   
%print(fig11, 'SpatialSimMatchFrac', '-dtiff', '-opengl', '-r600')
   
   
%% (Not used) 30x30 spatial example: tabulate edge match probs vs. edge endo probs
% Generate perfect samples from SPM, SZO, and CZO models and for each edge in the 
% graph, estimate its probability of being "matched" (i.e., having the same state at
% both nodes of the edge).  Tabulate how this matching probability depends on the
% unary parameters of the two nodes for each model.

% Get the X data and the endogenous probabilities
clear all
load('PerfSimXs.mat')  %X3, endo3 contain the case used in the paper.
X = X3;
endo = Endo3;
d = 30;
beta = [-2 2 2]';  %-Regression coefficients

%Set things up
ngen = 500;  %-Number of samples to generate per model
lams = [0 0.25 0.5 0.75 1.0];
nlams = length(lams);
mods = {'SZO','CZO','SPM'};
nmods = length(mods);

%-Use mod as a temporary model object for doing runs.
mod = ALRsimple;   
G = AutoLogistic.MakeGraph(d,d,1);
mod.Graph = G;
mod.Beta = beta;
mod.X = X;

%ismatch is a function for checking matched edges
ismatch = @(samp) samp(G.Edges.EndNodes(:,1))==samp(G.Edges.EndNodes(:,2));

%Create array for holding edge match fractions. Rows are edges; columns are lambdas;
%pages are models.
nedge = G.numedges;
MF = zeros(nedge,nlams,nmods);

%Loop through models and lambdas and generate samples to get matched frac.
tic
for k = 1:nmods
    switch mods{k}
        case 'SZO'
            mod.Coding = [0 1];
            mod.Centered = false;
        case 'SPM'
            mod.Coding = [-1/2 1/2];
            mod.Centered = false;
        case 'CZO'
            mod.Coding = [0 1];
            mod.Centered = true;
    end
    for j = 1:nlams
        mod.Lambda = lams(j);
        for g = 1:ngen
            tmp = mod.PerfectSample2;
            MF(:,j,k) = MF(:,j,k) + ismatch(tmp);
            disp([k j g])
        end
    end
end
MF = MF/ngen;
toc

% Save the MF array to file
%save('EdgeMatchData.mat','MF')


%% Tabulate output
% Group edges by their combined endogenous probability

load('EdgeMatchData.mat')

EdgeEndo = zeros(nedge,1);
for i = 1:nedge
    ix = G.Edges.EndNodes(i,:);
    EdgeEndo(i) = (endo(ix(1)) + endo(ix(2)))/2;
end

grp1 = EdgeEndo <= 1/3;
grp2 = (EdgeEndo > 1/3) & (EdgeEndo <= 2/3);
grp3 = ~(grp1|grp2);

%Create array to hold output.  one row per group, one page per lambda.  Columns are
%like:
%    SZO          CZO        SPM
%-----------  ----------  --------
% Q1 Q2 Q3     Q1 Q2 Q3   Q1 Q2 Q3

out = zeros(3,3*nmods,nlams);

for l = 1:nlams
    for m = 1:nmods
        cols = (3*(m-1)+1):(3*m);
        out(1,cols,l) = prctile(MF(grp1,l,m),[25 50 75]);
        out(2,cols,l) = prctile(MF(grp2,l,m),[25 50 75]);
        out(3,cols,l) = prctile(MF(grp3,l,m),[25 50 75]);
    end
end

disp(out)

%OR, use only two groups, and get quartiles for CZO and SPM.
% group 1 is edges with avg endo prob in the upper or lower thirds.
% group 2 is edges with avg endo prob between 1/3 and 2/3.
%Q1, Q2, Q3 are nlams-by-4.  First 2 columns are for CZO model, last 2 for SPM
%they hold the quartiles of the MF data for each group.
Q1 = zeros(nlams,4);
Q2 = zeros(nlams,4);
Q3 = zeros(nlams,4);

% Create the 2 groups
g1 = grp1|grp3;
g2 = grp2;

% Fill in Q1, Q2, Q3 using the array out
for l = 1:nlams
    %get the CZO values
    qqq1 = prctile(MF(g1,l,2),[25,50,75]);
    qqq2 = prctile(MF(g2,l,2),[25,50,75]);
    Q1(l,1) = qqq1(1);
    Q2(l,1) = qqq1(2);
    Q3(l,1) = qqq1(3);
    Q1(l,2) = qqq2(1);
    Q2(l,2) = qqq2(2);
    Q3(l,2) = qqq2(3);
    %get the SPM values
    qqq1 = prctile(MF(g1,l,3),[25,50,75]);
    qqq2 = prctile(MF(g2,l,3),[25,50,75]);
    Q1(l,3) = qqq1(1);
    Q2(l,3) = qqq1(2);
    Q3(l,3) = qqq1(3);
    Q1(l,4) = qqq2(1);
    Q2(l,4) = qqq2(2);
    Q3(l,4) = qqq2(3);
end

disp(' ')
disp(Q1)
disp(' ')
disp(Q2)
disp(' ')
disp(Q3)




%% Network Data Simulation


%% Create a table to hold the experimental run settings and the output.

%--- Simulation control variables ---
nreps = 200;  %-replications to do for each level combination.
lams = (0.25:0.25:1.5)';  %-lambda values to use.
models = {'SPM'; 'SZO'; 'CZO'};  %-names of models we're comparing.
cases = {'case 1'; 'case 2'};  %-graph types
objective = 'Hellinger';  %'Hellinger' or 'TV'

%--- Set up the expmt table ---
nlams = length(lams);
nmodels = length(models);
ncases = length(cases);
nruns = ncases*nlams;
expmt = table();
expmt.run = (1:(nruns))';
expmt.scenario = repelem(cases,nlams);
expmt.lambda = repmat(lams,ncases,1);
% SPM, SZO, CZO, Distances, endo1, endo2, endo3, marg1, marg2, marg3 to be filled
% in when we conduct the experiment.  Each cell of these table elements will
% contain a cell array or matrix with the nreps results for that experimental setting.
for i = 1:nmodels
    expmt.(models{i}) = cell(nruns,1);
end
expmt.Distance12 = cell(nruns,1);
expmt.Distance13 = cell(nruns,1);
expmt.endo1 = cell(nruns,1);
expmt.endo2 = cell(nruns,1);
expmt.endo3 = cell(nruns,1);
expmt.marg1 = cell(nruns,1);
expmt.marg2 = cell(nruns,1);
expmt.marg3 = cell(nruns,1);



%% Populate the SPM column
% This involves setting up the graphs, parameters, and covariates for the reference
% SPM model.

% Values used for all runs
beta = [0 1];
n = 16;
mu = 0;  %mean for X
sd = 1;  %sd for X

for r = 1:nruns
    
    %Initialize cell arrays to hold models etc.
    expmt.SPM{r} = cell(nreps,1);
    expmt.SZO{r} = cell(nreps,1);
    expmt.CZO{r} = cell(nreps,1);
    expmt.endo1{r} = zeros(n,nreps);
    expmt.endo2{r} = zeros(n,nreps);
    expmt.endo3{r} = zeros(n,nreps);
    expmt.marg1{r} = zeros(n,nreps);
    expmt.marg2{r} = zeros(n,nreps);
    expmt.marg3{r} = zeros(n,nreps);
    expmt.Distance12{r} = zeros(1,nreps);
    expmt.Distance13{r} = zeros(1,nreps);
    
    %Create nreps replicates of the model with different graphs and covariates.
    if strcmp(expmt.scenario{r},'case 1')
        m0 = 4;
        m = 1;
    else
        m0 = 2;
        m = 2;
    end
    A = triu(ones(m0),1);
    G0 = graph(A,'upper');
    G0.Edges.Weight = [];
    for i = 1:nreps
        G = G0;
        for j = (m0+1):n
            p = G.degree/sum(G.degree);
            ok = false;
            while ~ok  %Do multinomial draws until get m different nodes chosen.
                pick = find(mnrnd(m,p));
                ok = length(pick)==m;
            end
            G = addedge(G,j*ones(1,m),pick);
        end
        M = ALRsimple;
        M.Coding = [-1, 1];
        M.Graph = G;
        M.X = [ones(n,1) normrnd(mu,sd,n,1)];  
        M.Beta = beta;
        M.Lambda = 0;  %-only for getting endo
        expmt.endo1{r}(:,i) = M.MarginalProbability;
        M.Lambda = expmt.lambda(r);  %-the real lambda
        expmt.marg1{r}(:,i) = M.MarginalProbability;
        expmt.SPM{r}{i} = M;
        
        if mod(i,10)==0
            disp(['run ' num2str(r) ' rep ' num2str(i) ' complete.'])
        end
        
    end
            
end


%% Run the Simulation
tic

for r = 1:nruns
    
    for i = 1:nreps
        
        % Get probability table of reference model; create SZO model
        M = expmt.SPM{r}{i};
        pt = M.ProbabilityTable;
        M.Coding = [0 1];

        
        % Find the optimal SZO model
        if strcmp(objective, 'Hellinger')
            fcn = @(theta) GetHelli(theta, M, pt.Probability);
        else
            fcn = @(theta) GetTV(theta, M, pt.Probability);
        end
        opts = optimoptions(@fminunc, 'Display','none', ...
               'Algorithm','quasi-newton', 'MaxFunEvals',500);
        startval = [M.Beta; M.Lambda];
        [par, fval] = fminunc(fcn,startval,opts);

        % Add results to the table
        expmt.Distance12{r}(i) = fval;
        M.Beta = par(1:2);
        M.Lambda = 0;  %-only for getting endo
        expmt.endo2{r}(:,i) = M.MarginalProbability;
        M.Lambda = par(3);  %-the real lambda
        expmt.marg2{r}(:,i) = M.MarginalProbability;
        expmt.SZO{r}{i} = M;
        
        % Create the CZO model (add centering to M)
        M.Centered = true;
        
        % Find the optimal CZO model
        if strcmp(objective, 'Hellinger')
            fcn = @(theta) GetHelli(theta, M, pt.Probability);
        else
            fcn = @(theta) GetTV(theta, M, pt.Probability);
        end
        [par, fval] = fminunc(fcn,startval,opts);
        
        % Add results to the table 
        expmt.Distance13{r}(i) = fval;
        M.Beta = par(1:2);
        M.Lambda = 0;  %-only for getting endo
        expmt.endo3{r}(:,i) = M.MarginalProbability;
        M.Lambda = par(3);  %-the real lambda
        expmt.marg3{r}(:,i) = M.MarginalProbability;
        expmt.CZO{r}{i} = M;
        
        % Update progress to console
        if mod(i,10)==0
            disp(['run ' num2str(r) ' rep ' num2str(i) ' complete.'])
            disp(['total run time ' num2str(round(toc/60)) ' minutes so far.'])
        end
        
    end
    
end

% Save the experimental data to file
%save('DistanceExpmt200.mat','expmt')

%% Output analysis

%TODO: make a figure with the 4 boxplots for Helli dist and 4 for endo diff in a 2x4 
% arrangement. Then pick example(s) to demonstrate the case(s) we're considering.
% Then write it up!
% -to decide: include the SZO model, or just focus on CZO vs SPM?

%TODO: consider using "violin plots" instead of box plots? (some of our dists are
%bimodal)

load('DistanceExpmt200.mat')

nruns = length(expmt.run);
nreps = length(expmt.Distance12{1});
lams = unique(expmt.lambda);
nlams = length(lams);

%% Plot violin plots of hellinger distance for all cases, 
% grouped appropriately, on one axis. 

% Violin plot parameters.  we'll space the violin vertical centers 1 unit apart.
vpbins = 5;
vpcolor = 0.65*[1 1 1];
vpwidth = 0.4;

%Overall plot spacing control.  Horizontal axis will start at zero.
gap1 = 1;  %-gap between axis ends and first/last violin plot centers.
gap2 = 3;  %-gap between violin plot groups.
figwid = 7;  %-Figure width (inches)
fight = 2.5;  %-Figure height (inches)
labelpos = 1.1;  %-y location of group labels

%Locations for violin plots
vploc = [gap1+(0:5) gap1+5+gap2+(0:5) gap1+10+2*gap2+(0:5) gap1+15+3*gap2+(0:5)];


%Create the figure and axes
fig3 = figure(3);
clf
fig3.Units = 'inches';
fig3.Position = [1 1 figwid fight];
fig3.PaperPositionMode = 'auto';  %for printing at screen size
ax = axes();
ax.Position = [0.075 0.23 0.9 0.65];
ax.YLim = [0 1];
ax.XLim = [0 2*gap1 + 3*gap2 + 4*(nlams-1)];
ax.XTick = vploc;
ax.XTickLabel = repmat(lams,[4 1]);
ax.XTickLabelRotation = 90;
ax.YGrid = 'on';
ylabel('Hellinger Distance')
xlabel('Pairwise parameter of the baseline model')

%Put all Distance12's and Distanc13's into matrices
D12 = zeros(nreps,nruns);
D13 = zeros(nreps,nruns);
for r = 1:nruns
    D12(:,r) = expmt.Distance12{r};
    D13(:,r) = expmt.Distance13{r};
end

%Plot the violin plots.
%D12, case 1
for r = 1:nlams
    ViolinPlot(D12(:,r),vpbins,vploc(r),vpwidth,vpcolor);
end
text(gap1+2.5, labelpos, 'standard model, case 1', 'HorizontalAlignment', 'center')
%D12, case 2
for r = 1:nlams
    ViolinPlot(D12(:,r+nlams),vpbins,vploc(r+nlams),vpwidth,vpcolor);
end
text(gap1+gap2+7.5, labelpos, 'standard model, case 2', 'HorizontalAlignment', 'center')
%D13, case 1
for r = 1:nlams
    ViolinPlot(D13(:,r),vpbins,vploc(r+2*nlams),vpwidth,vpcolor);
end
text(gap1+2*gap2+12.5, labelpos, 'centered model, case 1', 'HorizontalAlignment', 'center')
%D13, case 2
for r = 1:nlams
    ViolinPlot(D13(:,r+nlams),vpbins,vploc(r+3*nlams),vpwidth,vpcolor);
end
text(gap1+3*gap2+17.5, labelpos, 'centered model, case 2', 'HorizontalAlignment', 'center')

%print(fig3, 'SimHelliDist5bin', '-dtiff', '-opengl', '-r600')



%% Make a similar plot for maximum endo difference
% Use the plot parameters etc. specified above.

% Note, there were two CZO case 2 runs (lambda=1.25, rep 11; lambda=1.5, rep 69) that
% produced NaNs when computing endogenous probabilities.  The cause of this was very
% large regression coefficients (absolute value approximately 100), that caused
% numerical problems.  So in the last group of violin plots be sure to remove NaNs.

%Put all max endo2 and endo3 values into matrices.
E12 = zeros(nreps,nruns);
E13 = zeros(nreps,nruns);
for r = 1:nruns
    d1 = abs(expmt.endo1{r} - expmt.endo2{r});
    E12(:,r) = max(d1,[],1,'omitnan');
    d2 = abs(expmt.endo1{r} - expmt.endo3{r});
    E13(:,r) = max(d2,[],1,'omitnan');    
end

%Create the figure and axes
fig4 = figure(4);
clf
fig4.Units = 'inches';
fig4.Position = [1 1 figwid fight];
fig4.PaperPositionMode = 'auto';  %for printing at screen size
ax = axes();
ax.Position = [0.075 0.23 0.9 0.65];
ax.YLim = [0 1];
ax.XLim = [0 2*gap1 + 3*gap2 + 4*(nlams-1)];
ax.XTick = vploc;
ax.XTickLabel = repmat(lams,[4 1]);
ax.XTickLabelRotation = 90;
ax.YGrid = 'on';
ylabel('Max. Endogenous Difference')
xlabel('Pairwise parameter, \lambda')

%Plot the violin plots.
%E12, case 1
for r = 1:nlams
    ViolinPlot(E12(:,r),vpbins,vploc(r),vpwidth,vpcolor);
end
text(gap1+2.5, labelpos, 'standard model, case 1', 'HorizontalAlignment', 'center')
%D12, case 2
for r = 1:nlams
    ViolinPlot(E12(:,r+nlams),vpbins,vploc(r+nlams),vpwidth,vpcolor);
end
text(gap1+gap2+7.5, labelpos, 'standard model, case 2', 'HorizontalAlignment', 'center')
%D13, case 1
for r = 1:nlams
    ViolinPlot(E13(:,r),vpbins,vploc(r+2*nlams),vpwidth,vpcolor);
end
text(gap1+2*gap2+12.5, labelpos, 'centered model, case 1', 'HorizontalAlignment', 'center')
%D13, case 2
for r = 1:nlams
    data = E13(:,r+nlams);
    data(isnan(data)) = [];
    ViolinPlot(data,vpbins,vploc(r+3*nlams),vpwidth,vpcolor);
end
text(gap1+3*gap2+17.5, labelpos, 'centered model, case 2', 'HorizontalAlignment', 'center')

%print(fig4, 'SimEndoDist5bin', '-dtiff', '-opengl', '-r600')


%% Make a similar plot for average endo difference
% small tweak of the above

%Put all max endo2 and endo3 values into matrices.
E12 = zeros(nreps,nruns);
E13 = zeros(nreps,nruns);
for r = 1:nruns
    d1 = abs(expmt.endo1{r} - expmt.endo2{r});
    E12(:,r) = mean(d1,1,'omitnan');
    d2 = abs(expmt.endo1{r} - expmt.endo3{r});
    E13(:,r) = mean(d2,1,'omitnan');    
end

%Create the figure and axes
fig4 = figure(4);
clf
fig4.Units = 'inches';
fig4.Position = [1 1 figwid fight];
fig4.PaperPositionMode = 'auto';  %for printing at screen size
ax = axes();
ax.Position = [0.075 0.23 0.9 0.65];
ax.YLim = [0 1];
ax.XLim = [0 2*gap1 + 3*gap2 + 4*(nlams-1)];
ax.XTick = vploc;
ax.XTickLabel = repmat(lams,[4 1]);
ax.XTickLabelRotation = 90;
ax.YGrid = 'on';
ylabel('Avg. |Endogenous Difference|')
xlabel('Pairwise parameter of the baseline model')

%Plot the violin plots.
%E12, case 1
for r = 1:nlams
    ViolinPlot(E12(:,r),vpbins,vploc(r),vpwidth,vpcolor);
end
text(gap1+2.5, labelpos, 'standard model, case 1', 'HorizontalAlignment', 'center')
%D12, case 2
for r = 1:nlams
    ViolinPlot(E12(:,r+nlams),vpbins,vploc(r+nlams),vpwidth,vpcolor);
end
text(gap1+gap2+7.5, labelpos, 'standard model, case 2', 'HorizontalAlignment', 'center')
%D13, case 1
for r = 1:nlams
    ViolinPlot(E13(:,r),vpbins,vploc(r+2*nlams),vpwidth,vpcolor);
end
text(gap1+2*gap2+12.5, labelpos, 'centered model, case 1', 'HorizontalAlignment', 'center')
%D13, case 2
for r = 1:nlams
    data = E13(:,r+nlams);
    data(isnan(data)) = [];
    ViolinPlot(data,vpbins,vploc(r+3*nlams),vpwidth,vpcolor);
end
text(gap1+3*gap2+17.5, labelpos, 'centered model, case 2', 'HorizontalAlignment', 'center')

%print(fig4, 'SimEndoDist5bin-avg', '-dtiff', '-opengl', '-r600')


%% Plot an example
% To illustrate the case, and to show the lambda effect for SPM model.
% We'll plot the graph with colored circles for the nodes.  Node color (grayscale)
% gives marginal probability of being high.  Since graph plot doesn't seem to allow
% outlines on node markers, we need to plot points with circular markers at the node 
% locations instead.
%
% Graphs are plotted with there centers at (0,0).  So we can scale and translate the
% XData properties of the graph plots to position several graphs on an invisible
% axes.  Here we'll plot 3 graphs, each with width 1.0.  Leave space of 0.1 between
% each graph.


%Choose a run an replicate to graph.
%For the paper: case 1 example was (4, 21). case 2 example was (7,5)
pickrun = 4;
pickrep = 21;

%Lambdas to use for the plot (3 of them)
trylams = [0 0.75 1.5];

%Figure and axes setup.
figwid = 7;  %-Figure width (inches)
fight = 2.5;  %-Figure height (inches)
lrgap = 0.15;  %-gap at either end of plot
midgap = 0.25;  %-gap between graphs
markersz = 8;   %-marker size
titlegap = 0.5;  %-gap between highest vertex and title.
tbgap = 0.5;  %-gap at top and bottom of the plot

%Create the figure and (invisible) axes
fig5 = figure(5);
clf
fig5.Units = 'inches';
fig5.Position = [1 1 figwid fight];
fig5.PaperPositionMode = 'auto';  %for printing at screen size
fig5.Color = [1 1 1];  %white background.
ax = axes();
ax.Position = [0 0 1 1];

%Get the model, the graph, and the marginal probs.
%***change SPM to SZO or CZO in the next line to see different models***
M = expmt.CZO{pickrun}{pickrep};
G = M.Graph;
M.Lambda = trylams(1);
marg1 = M.MarginalProbability;
M.Lambda = trylams(2);
marg2 = M.MarginalProbability;
M.Lambda = trylams(3);
marg3 = M.MarginalProbability;

%Plot the graph three times, in the desired locations
pl1 = plot(G,'EdgeColor','k','Nodecolor','none');
xmin = min(pl1.XData);
xmax = max(pl1.XData);
pl1.XData = (pl1.XData - xmin)/(xmax-xmin);
hold on
pl2 = plot(G,'EdgeColor','k','Nodecolor','none');
pl2.XData = pl1.XData + 1 + midgap;
pl3 = plot(G,'EdgeColor','k','Nodecolor','none');
pl3.XData = pl1.XData + 2 + 2*midgap;

%Set axis limits and make axes invisible
ax.XLim = [-lrgap 3+2*midgap+lrgap];
ax.Visible = 'off';

%Plot circles colored to show marginal probabilities
for i = 1:16
    plot(pl1.XData(i),pl1.YData(i),'ko','MarkerSize',markersz,...
         'MarkerEdgeColor',0.5*[1 1 1],'MarkerFaceColor',marg1(i)*[1 1 1]);
end
for i = 1:16
    plot(pl2.XData(i),pl2.YData(i),'ko','MarkerSize',markersz,...
         'MarkerEdgeColor',0.5*[1 1 1],'MarkerFaceColor',marg2(i)*[1 1 1]);
end
for i = 1:16
    plot(pl3.XData(i),pl3.YData(i),'ko','MarkerSize',markersz,...
         'MarkerEdgeColor',0.5*[1 1 1],'MarkerFaceColor',marg3(i)*[1 1 1]);
end

%Add titles and fix vertical ax limits
yloc = max(pl1.YData) + titlegap;
text(0.5,yloc,['\lambda = ' num2str(trylams(1))], 'HorizontalAlignment', 'center')
text(1.5 + midgap, yloc,['\lambda = ' num2str(trylams(2))], 'HorizontalAlignment', 'center')
text(2.5 + 2*midgap ,yloc,['\lambda = ' num2str(trylams(3))], 'HorizontalAlignment', 'center')
ax.YLim = [min(pl1.YData)-tbgap yloc+tbgap];

%print(fig5, 'SimGraphEx-7-5b', '-dtiff', '-opengl', '-r600')





%% Compare coefficients 

r = 12;  %-run to consider

%storage for coeffs (1s col: SZO; 2nd col: CZO)
beta0 = zeros(nreps,2);
beta1 = zeros(nreps,2);
lambda = zeros(nreps,2);
sumx = zeros(nreps,2);

for i = 1:nreps
    beta0(i,1) = expmt.SZO{r}{i}.Beta(1);
    beta0(i,2) = expmt.CZO{r}{i}.Beta(1);
    beta1(i,1) = expmt.SZO{r}{i}.Beta(2);
    beta1(i,2) = expmt.CZO{r}{i}.Beta(2);
    lambda(i,1) = expmt.SZO{r}{i}.Lambda;
    lambda(i,2) = expmt.CZO{r}{i}.Lambda;
    sumx(i,1) = sum(expmt.SZO{r}{i}.X(:,2));
    sumx(i,2) = sum(expmt.CZO{r}{i}.X(:,2));
end

figure(1)
plot(beta0(:,1),beta1(:,1),'k.')
xlim(prctile(beta0(:,1),[3 97]))
ylim(prctile(beta1(:,1),[3 97]))
figure(2)
plot(beta0(:,2),beta1(:,2),'k.')
xlim(prctile(beta0(:,2),[3 97]))
ylim(prctile(beta1(:,2),[3 97]))

figure(5)
clf
scatter3(beta0(:,2),beta1(:,2),lambda(:,2))
hold on
scatter3(0,1,1.5,10)
grid on
xlim(prctile(beta0(:,2),[3 97]))
ylim(prctile(beta1(:,2),[3 97]))
zlim(prctile(lambda(:,2),[3 97]))
xlabel('\beta_0')
ylabel('\beta_1')
zlabel('\lambda')

fig99 = figure(99);
clf
fig99.Units = 'centimeters';
fig5.Position = [1 1 8.5 7];
fig5.PaperPositionMode = 'auto';  %for printing at screen size
fig5.Color = [1 1 1];  %white background.
subplot(2,2,1)
scatter(lambda(:,2),beta0(:,2),4,'black','filled')
xlim(prctile(lambda(:,2),[2 96]))
ylim(prctile(beta0(:,2),[3 97]))
xlabel('\lambda')
ylabel('\beta_0')
grid on
subplot(2,2,2)
scatter(beta1(:,2),beta0(:,2),4,'black','filled')
xlim(prctile(beta1(:,2),[4 97]))
ylim(prctile(beta0(:,2),[3 97]))
xlabel('\beta_1')
ylabel('\beta_0')
grid on
subplot(2,2,4)
scatter(beta1(:,2),lambda(:,2),4,'black','filled')
xlim(prctile(beta1(:,2),[4 97]))
ylim(prctile(lambda(:,2),[2 96]))
xlabel('\beta_1')
ylabel('\lambda')
grid on

%print(fig99, 'NearestCZOParamsCase2Lam15', '-dtiff', '-opengl', '-r600')


%% compare the r most probable outcomes of SPM and CZO models.

runno = 12;  %-run to consider
i = 100;   %-replicate to consider
cumprob = 0.99;  %-we'll take the most likely models that add up to cumprob.

M1 = expmt.SPM{runno}{i};
M3 = expmt.CZO{runno}{i};

tbl1s = sortrows(M1.ProbabilityTable,17,'descend');
tbl2s = sortrows(M3.ProbabilityTable,17,'descend');
pick1 = cumsum(tbl1s.Probability)<=cumprob;
pick2 = cumsum(tbl2s.Probability)<=cumprob;
top1 = tbl1s(pick1,:);
top2 = tbl2s(pick2,:);
confs1 = top1{:,1:16}==1;
confs2 = top2{:,1:16}==1;

U = union(confs1,confs2,'rows');

disp([size(confs1,1) size(confs2,1) size(U,1)])


























