%% *** Clearing workspace *** 
clc;
clear all;
clf;

load('MSFT.mat')

%% *** Standard GARCH Model ***

SelectedPeriod     = [1:1:2529];
% SelectedPeriod   = [2530:1:4539];
% SelectedPeriod   = [4540:1:6553];

% Select JP Morgan Data
[Row Col]          = size(SelectedPeriod);
NumberObs          = Col;
MSFT               = MSFT(SelectedPeriod);

% Plot time series data 
% plot(MSFT);

% Compute Square Return
MSFTsqrt            = MSFT.* MSFT;

% Run GARCH with p and q varying from 1 to 4, compute log-likelihood
% LogL contains the value of p in the first column, q in the second,
% and we'll store the log likelihood in the third
LogLVal            = NaN(16, 4);
LogLVal(:, 1)      = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4]';
LogLVal(:, 2)      = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4]';

% Find number degree of freedom
kurcoeff           = kurtosis(MSFT);
DoF                = round(4+(6/(kurcoeff-3)));
tdist = struct('Name','t','DoF',DoF);

% Calculate Loglikelihood and store the data
for i=1:16
    p              = LogLVal(i, 1);
    q              = LogLVal(i, 2);
    GARCHModel     = garch(p, q);
    GARCHModeltdist= garch('Offset',NaN,'GARCHLags',p,'ARCHLags',q,...
                     'Distribution',tdist)
    [~, ~, logL, ~] = estimate(GARCHModel, MSFT);
    [~, ~, logLt, ~]= estimate(GARCHModeltdist, MSFT);
    LogLVal(i,3) = logL;
    LogLVal(i,4) = logLt;
end

% Perform likelihood ratio test for GARCH and GARCH-t
LogL_UNC         = LogLVal(2,3);
LogL_CONST       = LogLVal(1,3);
LogL_UNCt        = LogLVal(2,4);
LogL_CONSTt      = LogLVal(1,4);
LR               = -2*(LogL_CONST - LogL_UNC);
LRt              = -2*(LogL_CONSTt - LogL_UNCt);
PVal             = 1 - chi2cdf(LR, 1)
PValt            = 1 - chi2cdf(LRt, 1)

% Calculate Information Criteria (AIC and BIC)
NumPar           = [2;3;4;5;3;4;5;6;4;5;6;7;5;6;7;8];
[aic,bic]        = aicbic(LogLVal(:,3), NumPar, NumberObs*ones(16,1));
[minAIC,IndexAIC]= min(aic)
[minBIC,IndexBIC]= min(bic)
[aict,bict]        = aicbic(LogLVal(:,4), NumPar, NumberObs*ones(16,1));
[minAICt,IndexAICt]= min(aict)
[minBICt,IndexBICt]= min(bict)

% Display the fitted GARCH model
pAIC             = LogLVal(IndexAIC, 1); 
qAIC             = LogLVal(IndexAIC, 2);
pBIC             = LogLVal(IndexBIC, 1); 
qBIC             = LogLVal(IndexBIC, 2);
model            = 'GARCH';
disp(model)
disp(pAIC)
disp(qAIC)
disp(model)
disp(pBIC)
disp(qBIC)

% Display the fitted GARCH-t model
pAICt            = LogLVal(IndexAICt, 1); 
qAICt            = LogLVal(IndexAICt, 2);
pBICt            = LogLVal(IndexBICt, 1); 
qBICt            = LogLVal(IndexBICt, 2);
model            = 'GARCH-t';
disp(model)
disp(pAICt)
disp(qAICt)
disp(model)
disp(pBICt)
disp(qBICt)

%%%%%......... Backtesting.............%%%%%

% Create a model first
GARCHModel       = garch(1, 1);
GARCHModeltdist  = garch('Offset',NaN,'GARCHLags',1,'ARCHLags',1,...
                     'Distribution',tdist);

% Estimate GARCH for each series, save information for computing estimates
GARCHInfo         = NaN(4,2);

EstimatedModel    = estimate(GARCHModel,MSFT);
GARCHInfo(1, 1)   = EstimatedModel.Constant;
GARCHInfo(2, 1)   = EstimatedModel.ARCH{1};
GARCHInfo(3, 1)   = EstimatedModel.GARCH{1};
GARCHInfo(4, 1)   = EstimatedModel.UnconditionalVariance;

EstimatedModel    = estimate(GARCHModeltdist,MSFT);
GARCHInfo(1, 2)   = EstimatedModel.Constant;
GARCHInfo(2, 2)   = EstimatedModel.ARCH{1};
GARCHInfo(3, 2)   = EstimatedModel.GARCH{1};
GARCHInfo(4, 2)   = EstimatedModel.UnconditionalVariance;


% With the coefficient estimates, we can compute GARCH variances
GARCHVar          = NaN(NumberObs ,2);
GARCHVar(1,1)     = GARCHInfo(4,1);
GARCHVar(1,2)     = GARCHInfo(4,2);

for j=1:2
    for i=2:NumberObs  
       GARCHVar(i, j)= GARCHInfo(1, j) + ...
                         GARCHInfo(2, j) .* MSFTsqrt(i - 1, 1) + ...
                         GARCHInfo(3, j) .* GARCHVar(i - 1, j);
    end
end

% Finally, convert to volatility, eliminate first year of observations
GARCHVol            = sqrt(GARCHVar);

% Calculate Residual for each time period
Resid               = NaN(NumberObs,2);
Resid(:,1)          = MSFT./GARCHVol(:,1);
Resid(:,2)          = MSFT./GARCHVol(:,2);

% Diagnose distribution of Residual by QQ plots
% QQ plots versus Theoretical Normal Standard Distribution 
qqplot(Resid)
% QQ plots versus Theoretical t-student Distribution
nu                  = DoF ;
Tdist               = trnd(nu,NumberObs,2);
qqplot(Resid,Tdist)



% Produce a chart with log returns and +/- 2 Std. Dev. ***
% Chart options
FigNo                 = 1;
ChartData             = [MSFT -2*GARCHVol(:,1) -2*GARCHVol(:,2) 2*GARCHVol(:,1) 2*GARCHVol(:,2)];
StartDate             = 1990;
EndDate               = 2000; 
TimeLabels            = linspace(StartDate, EndDate, NumberObs);

% Charting the series
figure(FigNo)
titleStr = ['Daily Returns for MSFT'];
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate -0.15 0.15]);
grid on;
hold on;
% Note the use of handles (h1, h2, and so on) to limit chart legend to 
% only a few series
h1 = plot(TimeLabels, ChartData(:, 1), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
h2 = plot(TimeLabels, ChartData(:, 2), 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
h3 = plot(TimeLabels, ChartData(:, 4), 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red'); 
h4 = plot(TimeLabels, ChartData(:, 3), 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'green');
h5 = plot(TimeLabels, ChartData(:, 5), 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'green'); 

set(gca, 'YTickMode', 'manual');
set(gca, 'YTickLabel', num2str(100 .* get(gca, 'YTick')', '%1.0f%%'));              
legend1 = legend([h1 h3 h5], {'Log Returns', 'GARCH Vol', 'GARCH-t Vol'}, ...
                     'location','best');            
set(gcf, 'color', 'white');                  
FigNo   = FigNo + 1;

