%% Making plots for "parameter interpretation" section of ALR +/- paper.

% Make the code general enough to handle a variety of simple graphs:
%   - Designate one vertex as "Z1" the different one with its own alpha_1
%   - Let all other vertices take the same alpha_2
%   - Simple smoothing, parameter lambda
%   - Designate another vertex in the common group as "Z2" to compare probabilities
%     with Z1.

%% Input values

% User inputs
G = AutoLogistic.MakeGraph(1,2,1);  %Graph to use.
ix1 = 1;  %-Index of Z1 in the graph.
ix2 = 2;  %-Index of Z2 in the graph.
ngrid1d = 250;  %-Number of gridpoints to use for 1D plots
ngrid2d = 25;  %-Number of gripdpoints to use per axis for 2D plots

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


%% Try plotting P(Z1=h) as fcn of a1, a2 for fixed lambda

%***inputs***
lambda_max = 12;   %-Will use four lambdas equally spaced in [0,lambda_max]
a1_max = 8;       %-Will vary alpha1 between [-a1_max, a1_max]
a2_max = 8;       %-Will vary alpha2 between [-a2_max, a2_max]
variant = SPM;    %-The model variant to use.
hspc = [1.2 0.4 0.4 0.4 0.4];  %-Horizontal spacings for mysubplot()
vspc = [0.7 1.3];                %-Vertical spacings for mysubplot()
%************

lams = linspace(0,lambda_max,4);
a1rng = linspace(-a1_max,a1_max,ngrid2d);
a2rng = linspace(-a2_max,a2_max,ngrid2d);
[A1, A2] = meshgrid(a1rng,a2rng);
P1h = zeros(ngrid2d);

fig = myfigure(1,[31 9.5]);
clf

for k = 1:4
    variant.Lambda = lams(k);
    for i = 1:ngrid2d
        for j = 1:ngrid2d
            variant.X = assign_alpha(A1(i,j),A2(i,j));
            P1h(i,j) = variant.MarginalProbability(ix1);
            disp([i j])
        end
    end

    ax = mysubplot(1,4,k,hspc,vspc);
    %contourf(A1,A2,P1h,linspace(0,1,15))
    contourf(A1,A2,P1h,[0 0.05 0.1:0.1:0.9 0.95 1])
    colormap(gray)
    caxis([0 1])
    line([a1rng(1) a1rng(end)],[0 0],'Color','black','LineStyle','--')
    line([0 0],[a2rng(1) a2rng(end)],'Color','black','LineStyle','--')
    axis square
    if k==1
        ylabel('\alpha_2')
    else
        ax.YAxis.TickLabels = '';
    end
    xlabel('\alpha_1')
    titstr = ['\lambda = ', num2str(lams(k),3)];
    title(titstr)

end

% print(fig,'-r600','-dpdf','-opengl','CPM-P1h-vs-a1a2-neq2-final2')


%% Try plotting P(Z1=Z2) as fcn of a1, a2 for fixed lambda

%***inputs***
lambda_max = 1;   %-Will use four lambdas equally spaced in [0,lambda_max]
a1_max = 4;       %-Will vary alpha1 between [-a1_max, a1_max]
a2_max = 4;       %-Will vary alpha2 between [-a2_max, a2_max]
variant = SPM;    %-The model variant to use.
hspc = [1.2 0.4 0.4 0.4 0.4];  %-Horizontal spacings for mysubplot()
vspc = [0.7 1.3];                %-Vertical spacings for mysubplot()
%************

lams = linspace(0,lambda_max,4);
a1rng = linspace(-a1_max,a1_max,ngrid2d);
a2rng = linspace(-a2_max,a2_max,ngrid2d);
[A1, A2] = meshgrid(a1rng,a2rng);
P1eq2 = zeros(ngrid2d);

fig = myfigure(1,[31 9.5]);
clf

for k = 1:4
    variant.Lambda = lams(k);
    for i = 1:ngrid2d
        for j = 1:ngrid2d
            variant.X = assign_alpha(A1(i,j),A2(i,j));
            tbl = variant.ProbabilityTable;
            rowidx = tbl{:,ix1}==tbl{:,ix2};
            P1eq2(i,j) = sum(tbl.Probability(rowidx));
            disp([i j])
        end
    end

    ax = mysubplot(1,4,k,hspc,vspc);
    contourf(A1,A2,P1eq2,linspace(0,1,15))
    colormap(gray)
    caxis([0 1])
    line([a1rng(1) a1rng(end)],[0 0],'Color','black','LineStyle','--')
    line([0 0],[a2rng(1) a2rng(end)],'Color','black','LineStyle','--')
    axis square
    if k==1
        ylabel('\alpha_2')
    else
        ax.YAxis.TickLabels = '';
    end
    xlabel('\alpha_1')
    titstr = ['\lambda = ', num2str(lams(k),3)];
    title(titstr)

end

% print(fig,'-r600','-dpdf','-opengl','CZO-P1eq2-vs-a1a2-neq2-final')


%% Plot P(Z1=Z2) as a function of lambda, for different a1, a2 combos

%***inputs***
lambda_max = 30;   %-Will vary lambda  between [0 lambda_max]
a1rng = [-4 -2 0 2 4];       %-Will plot lines on each plot for each value in a1rng
a2_max = 4;       %-Will make four plots, with a2 varying from 0 to a2_max
variant = CZO;    %-The model variant to use.
hspc = [2 0.4 0.4 0.4 0.4];  %-Horizontal spacings for mysubplot()
vspc = [0.7 1.3];                %-Vertical spacings for mysubplot()
%************

lams = linspace(0,lambda_max,ngrid1d);
a2rng = linspace(0,a2_max,4);
fcnval = zeros(ngrid1d,1);

fig = myfigure(1,[31 9.5]);
clf

for k = 1:4
    
    a2 = a2rng(k);

    ax = mysubplot(1,4,k,hspc,vspc);
    axis([0 lams(end) 0 1])
    titstr = ['\alpha_2 = ', num2str(a2rng(k),3)];
    title(titstr)
    if k==1
        ylabel('P(Z_1 = Z_2)')
    else
        ax.YAxis.TickLabels = '';
    end
    xlabel('\lambda')
    grid on
    axis square

    for i = 1:length(a1rng)
        a1 = a1rng(i);
        variant.X = assign_alpha(a1,a2);
        for j = 1:ngrid1d
            variant.Lambda = lams(j);
            tbl = variant.ProbabilityTable;
            rowidx = tbl{:,ix1}==tbl{:,ix2};
            fcnval(j) = sum(tbl.Probability(rowidx));
            disp([i j])
        end
        line(lams,fcnval,'Color','black')
        text(lams(ngrid1d/2),fcnval(ngrid1d/2),['\alpha_1 = ' num2str(a1rng(i))])
    end

end

% print(fig,'-r600','-dpdf','-opengl','SPM-P1eq2-vs-lam-neq2')


%% Try plotting P(Z1=Z2) as fcn of lambda, a1 for fixed a2

%***inputs***
lambda_max = 5;   %-Will vary lambda over [0,lambda_max]
a1_min = -6;     %-Will vary alpha1 between [-a1_min, a1_max]
a1_max = 6;       %-Will vary alpha1 between [-a1_min, a1_max]
a2_max = 4;       %-Will make four plots, with a2 varying from 0 to a2_max
variant = SPM;    %-The model variant to use.
hspc = [1.3 0.5 0.5 0.5 0.4];  %-Horizontal spacings for mysubplot()
vspc = [0.7 1.3];                %-Vertical spacings for mysubplot()
%************

lams = linspace(0,lambda_max,ngrid2d);
a1rng = linspace(a1_min,a1_max,ngrid2d);
a2rng = linspace(0,a2_max,4);
[LAM, A1] = meshgrid(lams,a1rng);
P1eq2 = zeros(ngrid2d);

fig = myfigure(1,[31 9.5]);
clf
fig.Color = 'white';

for k = 1:4
    for i = 1:ngrid2d
        for j = 1:ngrid2d
            variant.X = assign_alpha(A1(i,j),a2rng(k));
            variant.Lambda = LAM(i,j);
            tbl = variant.ProbabilityTable;
            rowidx = tbl{:,ix1}==tbl{:,ix2};
            P1eq2(i,j) = sum(tbl.Probability(rowidx));
            disp([i j])
        end
    end

    ax = mysubplot(1,4,k,hspc,vspc);
    %contourf(LAM,A1,P1eq2,linspace(0,1,15))
    %contourf(LAM,A1,P1eq2,0:0.05:1)
    contourf(LAM,A1,P1eq2,[0 0.05 0.1:0.1:0.9 0.95 1])
    colormap(gray)
    caxis([0 1])
    line([0 lambda_max],[0 0],'Color','black','LineStyle','--')
    %line([0 0],[a1rng(1) a1rng(end)],'Color','black','LineStyle','--')
    axis square
    if k==1
        ylabel('\alpha_1')
        ax.YAxis.FontSize = 12;
    else
        ax.YAxis.TickLabels = '';
    end
    xlabel('\lambda')
    ax.XAxis.FontSize = 12;
    titstr = ['\alpha_2 = ', num2str(a2rng(k),3)];
    title(titstr)
    ax.Title.FontSize = 13;

end

% print(fig,'-r600','-dpdf','-opengl','SPM-P1eq2-vs-LamA1-neq2')


%% Compute equivalent models for a given n=2 case.

%Set up a the joint pmf of Z1, Z2
phh = 0.88;
phl = 0.009;
plh = 0.001;
pll = 0.11;
probvec = [phh phl plh pll]';

%Get corr for this case
p1 = phh+phl;
p2 = phh+plh;
correl = (phh - p1*p2)/sqrt(p1*p2*(1-p1)*(1-p2));

%Create an SPM model
mod = ALRsimple();
mod.Graph = AutoLogistic.MakeGraph(1,2,1);  %Graph to use.
mod.Coding = [-1 1];
mod.Centered = false;
mod.Beta = 1;
mod.X = [0 0]';  %initial value

%Optimize
OF = @(par) howclose(par,mod,probvec);
options = optimset('TolFun',1e-5,'TolX',1e-5);
[out, fval, exitflag, output] = fminsearch(OF,[0 0 0]',options);

%Put best params into SPM and check probability table
mod.X = out(1:2);
mod.Lambda = out(3);
mod.ProbabilityTable

%Transform this SPM model into the others; show results
correl
out
czo = mod.Transform([0 1],true);
czo = [czo.alpha; czo.assoc]
cpm = mod.Transform([-1 1],true);
cpm = [cpm.alpha; cpm.assoc]
szo = mod.Transform([0 1],false);
szo = [szo.alpha; szo.assoc]

















