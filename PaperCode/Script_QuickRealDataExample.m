%% Looking for a simple demo of limitations of centered zero/one ALR model.


%% *** Hydrocotyle data ***
clear all 
clc

% Load the data in
load('hydrocotyle')  %-loads table hydro.
[n, ~] = size(hydro);

% horz and vert are mixed up in this data set (it's Germany), so fix it...
vert = max(hydro.horz)-hydro.horz+1;
horz = hydro.vert;

% Make adjacency matrices.
% A4: 4-nearest neighbor.  A8: 8-nearest neighbor.
% Just do this by dumb looping for simplicity.
A4 = zeros(n);
A8 = zeros(n);
for i = 1:n
    E = hydro.horz(i);
    N = hydro.vert(i);
    leftright = abs(hydro.horz-E)==1 & hydro.vert==N;
    updown = abs(hydro.vert-N)==1 & hydro.horz==E;
    corner = abs(hydro.horz-E)==1 & abs(hydro.vert-N)==1;
    A4(i,:) = leftright | updown;
    A8(i,:) = leftright | updown | corner;
end

% Make graphs
G4 = graph(A4);
G8 = graph(A8);


%% Plot the presensce/absence data and the altitude covariate
figwd = 5*4.5;
fight = 5*4.5;
hmax = max(horz);
vmax = max(vert);

% First find all vertices that aren't connected to four neighbors (for drawing the
% outline of the map under the patches later)
edgeix = find(G4.degree~=4);  

%==== Draw the presence/absence map ====
fig1 = myfigure(1,[figwd fight]);
clf
for i = 1:length(edgeix)
    ptch = patch(horz(edgeix(i))+[-1 1 1 -1]/2, vert(edgeix(i))+[-1 -1 1 1]/2, 'white');
    ptch.FaceColor = 'none';
    ptch.LineWidth = 2;
    ptch.EdgeColor = 0.5*[1 1 1];    
end
hold on
for i = 1:n
    centx = horz(i);
    centy = vert(i);  
    col = hydro.obs(i)*[1 1 1];
    ptch = patch(centx+[-1 1 1 -1]/2, centy+[-1 -1 1 1]/2, col);
    ptch.EdgeColor = 'none';
end
xlim([0 hmax+1])
ylim([0 vmax+1])
fig1.CurrentAxes.Visible = 'off';
fig1.CurrentAxes.Position = [0.01 0.01 .98 .98];

%==== Draw the altitude map ====
fig2 = myfigure(2,[figwd fight]);
clf
for i = 1:length(edgeix)
    ptch = patch(horz(edgeix(i))+[-1 1 1 -1]/2, vert(edgeix(i))+[-1 -1 1 1]/2, col);
    ptch.FaceColor = 'none';
    ptch.LineWidth = 2;
    ptch.EdgeColor = 0.5*[1 1 1];    
end
hold on
for i = 1:n
    centx = horz(i);
    centy = vert(i);  
    col = hydro.alt(i)/max(hydro.alt)*[1 1 1];
    ptch = patch(centx+[-1 1 1 -1]/2, centy+[-1 -1 1 1]/2, col);
    ptch.EdgeColor = 'none';
end
xlim([0 hmax+1])
ylim([0 vmax+1])
fig2.CurrentAxes.Visible = 'off';
fig2.CurrentAxes.Position = [0.01 0.01 .98 .98];



%% Make ALR models

% Turn the hydro table into an X matrix.   Just use altitude as predictor. 
% Use G4 as the graph.
X = [ones(n,1) hydro.alt];
p = 2;  % only using one covariate + intercept.

% Make SPM model (model 2). ***Choose 4- or 8-neighbor graph***
SPM = ALRsimple();
SPM.Y = hydro.obs;
SPM.X = X;
SPM.Graph = G4;
% SPM.Graph = G8;
SPM.Beta = [1 0]';
SPM.Coding = [-1 1]/2;  %-plus/minus one half to make magnitudes more comparable.

% Make the corresponding CZO model (model 1).
CZO = SPM;
CZO.Coding = [0 1];
CZO.Centered = 1;

% Make the corresponding SZO model (model 0).
SZO = CZO;
SZO.Centered = 0;


%% Try max PSlik for the ALR model.
%Note: By inspection I found three different local optima for the CZO model, reached
%from different starting points.  They are:
%
%             (1)        (2)       (3)
%icept     -1.7415   -0.8151    2.7787      
%altitude  -0.1693    0.1481   -0.1955      
%lambda     1.5062    1.5554    1.4279      
%                                
%PSlik:   846.8727  854.6500  890.3448
%
% Solution (1) is the best, presumably the global optimum.  It is used here.
% Solution (2) fits okay but has opposite signs of the best solution (!)
% Solution (3) is the one you arrive at if you use the logistic regression param 
%              estimates (with lambda=0) as the starting point.
%

%------- First get logistic regression fit and predictions
fit = fitglm(hydro(:,[4 6]),'Distribution','binomial');
disp(fit.Coefficients)
%--------

% Note, if re-running the code, no need to do optimization, just load the data object
% that has the fitted models inside.  If you do run the optimization, check the CZO
% solution to make sure it's the right one...

OF0 = @(theta) PL(theta,SZO);
OF1 = @(theta) PL(theta,CZO);
OF2 = @(theta) PL(theta,SPM);

options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on', ...
          'TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',10^5);
%startval = [1 zeros(1,p)]';
%startval = [fit.Coefficients.Estimate; 0];
startval = 2*rand(p+1,1) - 1;
[out0, fval0, exitflag0, output0] = fminunc(OF0,startval,options);
[out1, fval1, exitflag1, output1] = fminunc(OF1,startval,options);
[out2, fval2, exitflag2, output2] = fminunc(OF2,startval,options);

disp(' ')
disp([out0 out1 out2])
disp(' ')
disp([fval0 fval1 fval2])

SZO.Beta = out0(1:p);
SZO.Lambda = out0(end);
CZO.Beta = out1(1:p);
CZO.Lambda = out1(end);
SPM.Beta = out2(1:p);
SPM.Lambda = out2(end);


%% Draw a bunch of perfect samples from the fitted models for future use

%load the data to avoid long run times...
load('hydroSamples2')

nsamp = 250;
samp0 = zeros(n,nsamp);
samp1 = zeros(n,nsamp);
samp2 = zeros(n,nsamp);

tic
parfor i = 1:nsamp
    disp(i)
    samp0(:,i) = SZO.PerfectSample2;
    samp1(:,i) = CZO.PerfectSample2;
    samp2(:,i) = SPM.PerfectSample2;
end
toc

% hydroSamples2 is for the best MPL parameter settings found.
% original file hydroSamples was for a parameter estimate that was a local optimum.
%save('hydroSamples2','hydro','SZO','CZO','SPM','samp0','samp1','samp2')

%% Compare marginal prediction maps (including logistic)
%load('hydroSamples2')

pred00 = fit.Fitted.Probability;
pred0 = sum(samp0==SZO.Coding(2), 2)/nsamp;
pred1 = sum(samp1==CZO.Coding(2), 2)/nsamp;
pred2 = sum(samp2==SPM.Coding(2), 2)/nsamp;

%==== Logistic regression map ====
fig3 = myfigure(3,[figwd fight]);
clf
for i = 1:length(edgeix)
    ptch = patch(horz(edgeix(i))+[-1 1 1 -1]/2, vert(edgeix(i))+[-1 -1 1 1]/2, 'white');
    ptch.FaceColor = 'none';
    ptch.LineWidth = 2;
    ptch.EdgeColor = 0.5*[1 1 1];    
end
hold on
for i = 1:n
    col = pred00(i)*[1 1 1];
    ptch = patch(horz(i)+[-1 1 1 -1]/2, vert(i)+[-1 -1 1 1]/2, col);
    ptch.EdgeColor = 'none';
end
xlim([0 hmax+1])
ylim([0 vmax+1])
fig3.CurrentAxes.Visible = 'off';
fig3.CurrentAxes.Position = [0.01 0.01 .98 .98];

%==== SZO map ====
fig4 = myfigure(4,[figwd fight]);
clf
for i = 1:length(edgeix)
    ptch = patch(horz(edgeix(i))+[-1 1 1 -1]/2, vert(edgeix(i))+[-1 -1 1 1]/2, 'white');
    ptch.FaceColor = 'none';
    ptch.LineWidth = 2;
    ptch.EdgeColor = 0.5*[1 1 1];    
end
hold on
for i = 1:n
    col = pred0(i)*[1 1 1];
    ptch = patch(horz(i)+[-1 1 1 -1]/2, vert(i)+[-1 -1 1 1]/2, col);
    ptch.EdgeColor = 'none';
end
xlim([0 hmax+1])
ylim([0 vmax+1])
fig4.CurrentAxes.Visible = 'off';
fig4.CurrentAxes.Position = [0.01 0.01 .98 .98];

%==== CZO map ====
fig5 = myfigure(5,[figwd fight]);
clf
for i = 1:length(edgeix)
    ptch = patch(horz(edgeix(i))+[-1 1 1 -1]/2, vert(edgeix(i))+[-1 -1 1 1]/2, 'white');
    ptch.FaceColor = 'none';
    ptch.LineWidth = 2;
    ptch.EdgeColor = 0.5*[1 1 1];    
end
hold on
for i = 1:n
    col = pred1(i)*[1 1 1];
    ptch = patch(horz(i)+[-1 1 1 -1]/2, vert(i)+[-1 -1 1 1]/2, col);
    ptch.EdgeColor = 'none';
end
xlim([0 hmax+1])
ylim([0 vmax+1])
fig5.CurrentAxes.Visible = 'off';
fig5.CurrentAxes.Position = [0.01 0.01 .98 .98];

%==== SPM map ====
fig6 = myfigure(6,[figwd fight]);
clf
for i = 1:length(edgeix)
    ptch = patch(horz(edgeix(i))+[-1 1 1 -1]/2, vert(edgeix(i))+[-1 -1 1 1]/2, 'white');
    ptch.FaceColor = 'none';
    ptch.LineWidth = 2;
    ptch.EdgeColor = 0.5*[1 1 1];    
end
hold on
for i = 1:n
    col = pred2(i)*[1 1 1];
    ptch = patch(horz(i)+[-1 1 1 -1]/2, vert(i)+[-1 -1 1 1]/2, col);
    ptch.EdgeColor = 'none';
end
xlim([0 hmax+1])
ylim([0 vmax+1])
fig6.CurrentAxes.Visible = 'off';
fig6.CurrentAxes.Position = [0.01 0.01 .98 .98];


%% DELETE THIS CELL LATER
%==== Logistic regression map ====
fig99 = myfigure(88,[figwd fight]);
clf
for i = 1:length(edgeix)
    ptch = patch(horz(edgeix(i))+[-1 1 1 -1]/2, vert(edgeix(i))+[-1 -1 1 1]/2, 'white');
    ptch.FaceColor = 'none';
    ptch.LineWidth = 2;
    ptch.EdgeColor = 0.5*[1 1 1];    
end
hold on
for i = 1:n
    col = pred2_droplam(i)*[1 1 1];
    ptch = patch(horz(i)+[-1 1 1 -1]/2, vert(i)+[-1 -1 1 1]/2, col);
    ptch.EdgeColor = 'none';
end
xlim([0 hmax+1])
ylim([0 vmax+1])
fig99.CurrentAxes.Visible = 'off';
fig99.CurrentAxes.Position = [0.01 0.01 .98 .98];





%% Plot the ROC curves as a summary of GOF. Also display AUC.

[xx00, yy00, tt00, auc00] = perfcurve(hydro.obs,fit.Fitted.Probability,1);
[xx0, yy0, tt0, auc0] = perfcurve(hydro.obs,pred0,1);
[xx1, yy1, tt1, auc1] = perfcurve(hydro.obs,pred1,1);
[xx2, yy2, tt2, auc2] = perfcurve(hydro.obs,pred2,1);

disp([auc00 auc0 auc1 auc2])

fig7 = myfigure(7,[9 7]);
clf
plot(xx00,yy00,'k-.')
hold on
plot(xx0,yy0,'k-')
plot(xx1,yy1,'k-','color',0.6*[1 1 1],'linewidth',2)
plot(xx2,yy2,'k-','linewidth',2)

xlabel('False positive rate')
ylabel('True positive rate')

legend('logistic regression','traditional model','centered model','symmetric model')

%print(fig7, 'ROCcurve', '-dtiff', '-opengl', '-r600')


%% Draw samples from the fitted models to use for parametric bootstrap
% (which was done in HHC11E; they used 2000 samples)

nboot = 2000;
bootsamp0 = zeros(n,nboot);
bootsamp1 = zeros(n,nboot);
bootsamp2 = zeros(n,nboot);

tic
parfor i = 1:nboot
    disp(i)
    bootsamp0(:,i) = SZO.PerfectSample2;
    bootsamp1(:,i) = CZO.PerfectSample2;
    bootsamp2(:,i) = SPM.PerfectSample2;
end
toc

%save('hydroSamples2boot','bootsamp0','bootsamp1','bootsamp2')


%% Get MPL parameter estimates for each bootstrap sample.

%load('hydroSamples2boot')

THETA0 = zeros(nboot,3);  %-to hold bootstrap param estimates
THETA1 = zeros(nboot,3);
THETA2 = zeros(nboot,3);

psLik = zeros(nboot,3);  %-to hold pslikelihood minimum values (one column per model)
flag = zeros(nboot,3);  %-to hold optimization exit flags.

mod0 = SZO;
mod1 = CZO;
mod2 = SPM;

%start optimization from the original estimate
start0 = [SZO.Beta; SZO.Lambda];
start1 = [CZO.Beta; CZO.Lambda];
start2 = [SPM.Beta; SPM.Lambda];

options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on', ...
          'TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',10^5,'Display','none');

for i = 1:nboot
    
    disp(i)
    
    mod0.Y = bootsamp0(:,i);
    mod1.Y = bootsamp1(:,i);
    mod2.Y = bootsamp2(:,i);
    
    OF0 = @(theta) PL(theta,mod0);
    OF1 = @(theta) PL(theta,mod1);
    OF2 = @(theta) PL(theta,mod2);

    [out0, fval0, exitflag0 ] = fminunc(OF0,start0,options);
    [out1, fval1, exitflag1 ] = fminunc(OF1,start1,options);
    [out2, fval2, exitflag2] = fminunc(OF2,start2,options);

    THETA0(i,:) = out0;
    THETA1(i,:) = out1;
    THETA2(i,:) = out2;
    
    psLik(i,:) = [fval0 fval1 fval2];
    flag(i,:) = [exitflag0 exitflag1 exitflag2];
    
end

%save('hydroBootResults','THETA0','THETA1','THETA2','psLik','flag')


%% Use the bootstrap results to get precision estimates for params.
%Notes:
% - All optimization runs gave exit flag 3.
% - All param estimates had fairly smooth-looking unimodal distributions.

%load('hydroBootResults')
disp(' ')
disp('Original estimates')
disp([start0'; start1'; start2';])
disp(' ')
disp('Means of bootstrap estimates')
disp([mean(THETA0,1); mean(THETA1,1); mean(THETA2,1)])
disp(' ')
disp('Standard deviations of bootstrap estimates')
disp([std(THETA0,1); std(THETA1,1); std(THETA2,1)])
disp(' ')
disp('[2.5th; 97.5th] percentile of bootstrap estimates')
disp([prctile(THETA0,[2.5 97.5],1); prctile(THETA1,[2.5 97.5],1); prctile(THETA2,[2.5 97.5],1)])
disp(' ')
disp('width of 95% confidence interval')
disp([prctile(THETA0,97.5,1); prctile(THETA1,97.5,1); prctile(THETA2,97.5,1)] ...
     - [prctile(THETA0,2.5,1); prctile(THETA1,2.5,1); prctile(THETA2,2.5,1)])
 



%% Get samples for "covariate impact" as in BGW15ArXiVb
% To get the impact measures, need to drop out (set to zero) individual terms in the
% model and then re-evaluate the model predictions/fitted values.  Getting model 
% predictions requires sampling.
%
% For completeness' sake, drop out all three parameters in turn (intercept, altitude,
% AND association parameter).  Then draw nsamp samples from the model.  Do this for
% all three models.

nsamp = 250;

samp0_dropint = zeros(n,nsamp);
samp1_dropint = zeros(n,nsamp);
samp2_dropint = zeros(n,nsamp);
samp0_dropalt = zeros(n,nsamp);
samp1_dropalt = zeros(n,nsamp);
samp2_dropalt = zeros(n,nsamp);
samp0_droplam = zeros(n,nsamp);
samp1_droplam = zeros(n,nsamp);
samp2_droplam = zeros(n,nsamp);

%drop-intercept models
mod0int = SZO;
mod0int.Beta(1) = 0;
mod1int = CZO;
mod1int.Beta(1) = 0;
mod2int = SPM;
mod2int.Beta(1) = 0;

%drop-altitude models
mod0alt = SZO;
mod0alt.Beta(2) = 0;
mod1alt = CZO;
mod1alt.Beta(2) = 0;
mod2alt = SPM;
mod2alt.Beta(2) = 0;

%drop-lambda models
mod0lam = SZO;
mod0lam.Lambda = 0;
mod1lam = CZO;
mod1lam.Lambda = 0;
mod2lam = SPM;
mod2lam.Lambda = 0;

%draw samples
tic
parfor i = 1:nsamp
    
    samp0_dropint(:,i) = mod0int.PerfectSample2;
    samp1_dropint(:,i) = mod1int.PerfectSample2;
    samp2_dropint(:,i) = mod2int.PerfectSample2;
    samp0_dropalt(:,i) = mod0alt.PerfectSample2;
    samp1_dropalt(:,i) = mod1alt.PerfectSample2;
    samp2_dropalt(:,i) = mod2alt.PerfectSample2;
    samp0_droplam(:,i) = mod0lam.PerfectSample2;
    samp1_droplam(:,i) = mod1lam.PerfectSample2;
    samp2_droplam(:,i) = mod2lam.PerfectSample2;

end
toc

%Do it for logistic regression
modLR = SZO;
modLR.Lambda=0;
modLR.Beta = fit.Coefficients.Estimate;

modLR_dropint = modLR;
modLR_dropint.Beta(1) = 0;
modLR_dropalt = modLR;
modLR_dropalt.Beta(2) = 0;

sampLR_dropint = zeros(n,nsamp);
sampLR_dropalt = zeros(n,nsamp);
for i = 1:nsamp
    disp(i)
    sampLR_dropint(:,i) = modLR_dropint.PerfectSample2;
    sampLR_dropalt(:,i) = modLR_dropalt.PerfectSample2;
end

% save('hydroSamples2_drop','samp0_dropint','samp1_dropint','samp2_dropint', ...
%      'samp0_dropalt','samp1_dropalt','samp2_dropalt', ...
%      'samp0_droplam','samp1_droplam','samp2_droplam', ...
%      'sampLR_dropint', 'sampLR_dropalt')


%% Get covariate impact measures
load('hydroSamples2_drop')

pred0_dropint = sum(samp0_dropint==SZO.Coding(2), 2)/nsamp;
pred1_dropint = sum(samp1_dropint==CZO.Coding(2), 2)/nsamp;
pred2_dropint = sum(samp2_dropint==SPM.Coding(2), 2)/nsamp;
pred0_dropalt = sum(samp0_dropalt==SZO.Coding(2), 2)/nsamp;
pred1_dropalt = sum(samp1_dropalt==CZO.Coding(2), 2)/nsamp;
pred2_dropalt = sum(samp2_dropalt==SPM.Coding(2), 2)/nsamp;
pred0_droplam = sum(samp0_droplam==SZO.Coding(2), 2)/nsamp;
pred1_droplam = sum(samp1_droplam==CZO.Coding(2), 2)/nsamp;
pred2_droplam = sum(samp2_droplam==SPM.Coding(2), 2)/nsamp;

%is it best to make these absolute values, or would signs be better?
useabs = false;
myabs = @(x) useabs*abs(x) + ~useabs*x;

impact0int = sum(myabs(pred0 - pred0_dropint))/n;
impact1int = sum(myabs(pred1 - pred1_dropint))/n;
impact2int = sum(myabs(pred2 - pred2_dropint))/n;

disp(' ')
disp([impact0int impact1int impact2int])

impact0alt = sum(myabs(pred0 - pred0_dropalt))/n;
impact1alt = sum(myabs(pred1 - pred1_dropalt))/n;
impact2alt = sum(myabs(pred2 - pred2_dropalt))/n;

disp([impact0alt impact1alt impact2alt])

impact0lam = sum(myabs(pred0 - pred0_droplam))/n;
impact1lam = sum(myabs(pred1 - pred1_droplam))/n;
impact2lam = sum(myabs(pred2 - pred2_droplam))/n;

disp([impact0lam impact1lam impact2lam])

%logistic regression
predLR_dropint = sum(sampLR_dropint==modLR.Coding(2), 2)/nsamp;
predLR_dropalt = sum(sampLR_dropalt==modLR.Coding(2), 2)/nsamp;

impactLRint = sum(myabs(pred00 - predLR_dropint))/n;
impactLRalt = sum(myabs(pred00 - predLR_dropalt))/n;

disp(' ')
disp([impactLRint; impactLRalt])


%% ***********************************************************************

%% *** Red deer data ***
clear all 
clc

% Load the data in
load('reddeer')  %-loads table reddeer.
[n, p] = size(reddeer);

% Plot the observations
figure(1)
clf
obshigh = reddeer.obs==1;
obslow=~obshigh;
scatter(reddeer{obslow,1},reddeer{obslow,2})
hold on
scatter(reddeer{obshigh,1},reddeer{obshigh,2},'red')
xlabel('Eastings')
ylabel('Northings')

% Make adjacency matrices.
% A4: 4-nearest neighbor.  A8: 8-nearest neighbor.
% Just do this by dumb looping for simplicity.
A4 = zeros(n);
A8 = zeros(n);
for i = 1:n
    E = reddeer.east(i);
    N = reddeer.north(i);
    leftright = abs(reddeer.east-E)==1 & reddeer.north==N;
    updown = abs(reddeer.north-N)==1 & reddeer.east==E;
    corner = abs(reddeer.east-E)==1 & abs(reddeer.north-N)==1;
    A4(i,:) = leftright | updown;
    A8(i,:) = leftright | updown | corner;
end

% Make graphs
G4 = graph(A4);
G8 = graph(A8);

% Plot graphs
figure(2)
pl4 = plot(G4,'XData',reddeer.east,'YData',reddeer.north);
pl4.NodeCData = 1*obshigh;
colormap([[0 0 0]; [1 0 0]])
figure(3)
pl8 = plot(G8,'XData',reddeer.east,'YData',reddeer.north);
pl8.NodeCData = 1*obshigh;
colormap([[0 0 0]; [1 0 0]])



%% Make ALR models
% NOTE: If we include all predictors it seems the two models behave similarly (need
% detailed checking to be sure, but no clear signs of difference).  
%  (the following observation was due to CZO max PL estimate stopped at a crappy
%  local optimum):
    % BUT, if we drop the east and north predictors, big differences emerge. The centered
    % model predicts presence (high) all over the place (84% of the nodes have marginal
    % prob > 1/2, and the median proportion of high-valued nodes in random draws was
    % 67%). The standard plus/minus model had zero marginal probs over 1/2 (the max was
    % 0.39), and the median percent of high vertices was 13.4%.  (the observed response
    % had 14.9% high)
%
% NOTE: AMB96JAE used only alt^2 and pine as their "reduced" model.  So let's use
% that too!

% Turn the reddeer table into an X matrix. Note that references used altitude^2, not
% the altitude directly, as a predictor.  Also standardize the matrix at the end.
X = reddeer{:,1:5};
X(:,3) = X(:,3).^2;
X = [ones(n,1) standardize(X)];
X = X(:,[1 4 5]);  %-Keep only intercept, alt^2, and pine
p = 3;  %-update number of predictors to match reduced X.


% Make SPM model (model 2). ***Choose 4- or 8-neighbor graph***
SPM = ALRsimple();
SPM.Y = reddeer.obs;
SPM.X = X;
SPM.Graph = G4;
%SPM.Graph = G8;
SPM.Beta = [1 0 0]';
SPM.Coding = [-1 1]/2;  %-plus/minus one half to make magnitudes more comparable.

% Make the corresponding CZO model (model 1).
CZO = SPM;
CZO.Coding = [0 1];
CZO.Centered = 1;

% And the corresponding SZO model (model 0).
SZO = CZO;
SZO.Centered = 0;


%% Try max PSlik for the ALR model.

OF0 = @(theta) PL(theta,SZO);
OF1 = @(theta) PL(theta,CZO);
OF2 = @(theta) PL(theta,SPM);

options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on', ...
          'TolFun',1e-6,'TolX',1e-6,'MaxFunEvals',10^5);
% startval = [1 0 0 0 0 0 0 0]';
startval = 2*rand(p+1,1) - 1;
[out0, fval0, exitflag0, output0] = fminunc(OF0,startval,options);
[out1, fval1, exitflag1, output1] = fminunc(OF1,startval,options);
[out2, fval2, exitflag2, output2] = fminunc(OF2,startval,options);

disp(' ')
disp([out0 out1 out2])
disp(' ')
disp([fval0 fval1 fval2])

SZO.Beta = out0(1:p);
SZO.Lambda = out0(end);
CZO.Beta = out1(1:p);
CZO.Lambda = out1(end);
SPM.Beta = out2(1:p);
SPM.Lambda = out2(end);


%% Draw a bunch of perfect samples from the fitted models for future use

nsamp = 25;
samp0 = zeros(n,nsamp);
samp1 = zeros(n,nsamp);
samp2 = zeros(n,nsamp);
tic
for i = 1:nsamp
    disp(i)
    samp0(:,i) = SZO.PerfectSample2;
    samp1(:,i) = CZO.PerfectSample2;
    samp2(:,i) = SPM.PerfectSample2;
end
toc

%save('RedDeerSamples','reddeer','SZO','CZO','SPM','samp0','samp1','samp2')


%% Compare marginal prediction maps
pred0 = sum(samp0==SZO.Coding(2), 2)/nsamp;
pred1 = sum(samp1==CZO.Coding(2), 2)/nsamp;
pred2 = sum(samp2==SPM.Coding(2), 2)/nsamp;

figure(5)
clf
for k = 1:3
    subplot(1,3,k)
    hold on
    for i = 1:n
        switch k
            case 1
                col = pred0(i)*[1 1 1];
            case 2
                col = pred1(i)*[1 1 1];
            case 3
                col = pred2(i)*[1 1 1];
        end
        centx = reddeer.east(i);
        centy = reddeer.north(i);
        ptch = patch(centx+[-1 1 1 -1]/2, centy+[-1 -1 1 1]/2, col);
%         if obshigh(i)
%             ptch.EdgeColor = [1 1 1];
%             ptch.EdgeAlpha = 0.5;
%         else
            ptch.EdgeColor = 'none';
%         end
    end
    axis square; axis tight; grid on
end


%% Plot a sampling result

pick = 21;

figure(6)
clf
subplot(1,3,1)
hold on
for i = 1:n
    centx = reddeer.east(i);
    centy = reddeer.north(i);
    col = (samp0(i,pick)==SZO.Coding(2))*[1 1 1];
    ptch = patch(centx+[-1 1 1 -1]/2, centy+[-1 -1 1 1]/2, col);
    ptch.EdgeColor = 'none';
end
axis square; axis tight; grid on

subplot(1,3,2)
hold on
for i = 1:n
    centx = reddeer.east(i);
    centy = reddeer.north(i);
    col = (samp1(i,pick)==CZO.Coding(2))*[1 1 1];
    ptch = patch(centx+[-1 1 1 -1]/2, centy+[-1 -1 1 1]/2, col);
    ptch.EdgeColor = 'none';
end
axis square; axis tight; grid on

subplot(1,3,3)
hold on
for i = 1:n
    centx = reddeer.east(i);
    centy = reddeer.north(i);
    col = (samp2(i,pick)==SPM.Coding(2))*[1 1 1];
    ptch = patch(centx+[-1 1 1 -1]/2, centy+[-1 -1 1 1]/2, col);
    ptch.EdgeColor = 'none';
end
axis square; axis tight; grid on


%% Try to measure "covariate impact" as in BGW15ArXiVb

dropcov = 2;  %column of X to drop out.

samp0_drop = zeros(n,nsamp);
samp1_drop = zeros(n,nsamp);
samp2_drop = zeros(n,nsamp);
mod0 = SZO;
mod0.Beta(dropcov) = 0;
mod1 = CZO;
mod1.Beta(dropcov) = 0;
mod2 = SPM;
mod2.Beta(dropcov) = 0;
tic
for i = 1:nsamp
    samp0_drop(:,i) = mod0.PerfectSample2;
    samp1_drop(:,i) = mod1.PerfectSample2;
    samp2_drop(:,i) = mod2.PerfectSample2;
end
toc

pred0_drop = sum(samp0_drop==SZO.Coding(2), 2)/nsamp;
pred1_drop = sum(samp1_drop==CZO.Coding(2), 2)/nsamp;
pred2_drop = sum(samp2_drop==SPM.Coding(2), 2)/nsamp;

%is it best to make these absolute values, or would signs be better?
impact0 = sum(abs(pred0 - pred0_drop))/n;
impact1 = sum(abs(pred1 - pred1_drop))/n;
impact2 = sum(abs(pred2 - pred2_drop))/n;

disp([impact0 impact1 impact2])




%% Use assortativity to check GOF
nsim = nsamp;
ass1 = zeros(nsim,1);
ass2 = zeros(nsim,1);
mod1 = CZO;
mod2 = SPM;
tic
for i = 1:nsim
    mod1.Y = samp1(:,i);
    mod2.Y = samp2(:,i);
    ass1(i) = assortativity(mod1);
    ass2(i) = assortativity(mod2);
end
toc
disp(assortativity(SPM))
disp(median(ass1))
disp(median(ass2))



%% ****************************************************************************

%% *** ipsCC data ***
clear all
clc

% Load the data in
%load('ipsCC') %-with disconnected nodes
load('ipsCCgc')  %-only the "giant connected" part.

n = length(Y);

% Make the ALR objects.  Standardize X in hopes of improving numerics.
SPM = ALRsimple();
SPM.Y = Y;
SPM.X = [ones(n,1) standardize(X)];
SPM.Graph = graph(Adj);
SPM.Beta = [1 0 0 0 0 0 0]';
SPM.Coding = [-1 1]/2;  %-plus/minus one half to make magnitudes more comparable.

CZO = SPM;
CZO.Coding = [0 1];
CZO.Centered = 1;

SZO = CZO;
SZO.Centered = 0;


%% Try max PSlik for the ALR model.

OF0 = @(theta) PL(theta,SZO);
OF1 = @(theta) PL(theta,CZO);
OF2 = @(theta) PL(theta,SPM);
P = SPM.P;

options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on', ...
          'TolFun',1e-6,'TolX',1e-6,'MaxFunEvals',10^5);
% startval = [1 0 0 0 0 0 0 0]';
startval = 2*rand(P+1,1) - 1;
[out0, fval0, exitflag0, output0] = fminunc(OF0,startval,options);
[out1, fval1, exitflag1, output1] = fminunc(OF1,startval,options);
[out2, fval2, exitflag2, output2] = fminunc(OF2,startval,options);

disp(' ')
disp([out0 out1 out2])
disp(' ')
disp([fval0 fval1 fval2])

SZO.Beta = out0(1:P);
SZO.Lambda = out0(end);
CZO.Beta = out1(1:P);
CZO.Lambda = out1(end);
SPM.Beta = out2(1:P);
SPM.Lambda = out2(end);


%% Use assortativity to check GOF
nsim = 500;
ass1 = zeros(nsim,1);
ass2 = zeros(nsim,1);
mod1 = CZO;
mod2 = SPM;
tic
for i = 1:nsim
    mod1.Y = CZO.PerfectSample2;
    mod2.Y = SPM.PerfectSample2;
    ass1(i) = assortativity(mod1);
    ass2(i) = assortativity(mod2);
end
toc

%%
lows = CZO.Y==CZO.Coding(1);
CondProbHigh = CZO.ConditionalProbability;
CondProbHigh(lows) = 1 - CondProbHigh(lows);
pred = CondProbHigh > 0.5;
actual = CZO.Y==CZO.Coding(2);
errrate = 1 - sum(pred==actual)/n;

disp(' ')
disp(errrate)

%%

lows = SPM.Y==SPM.Coding(1);
CondProbHigh = SPM.ConditionalProbability;
CondProbHigh(lows) = 1 - CondProbHigh(lows);
pred = CondProbHigh > 0.5;
actual = SPM.Y==SPM.Coding(2);
errrate = 1 - sum(pred==actual)/n;
disp(errrate)








