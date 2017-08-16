%% London School of Economics
%  FM442 - Michaelmas Term 2016
%  Author: Domingos Romualdo
%  Date: 25 October 2016
%  Assignment 3 - Suggested Solution

%% *** Clearing workspace *** 
clc;
clear all;
clf;

%% *** Part 1: Loading daily data ***
% Note: the import file below is the same as that for Homework 2.
% The code is also nearly identical in parts 1-3
Data_D           = importdata('FM442 - HW2 - Daily Data.csv');
Tickers          = ['MSFT'; 'XOM '; 'GE  '; 'JPM '; 'INTC'; 'C   '; ...
                       'SPX '];

%% *** Part 2 - Cleaning up and organizing the data ***
% Retrieving information about data
[NRows, NCols]   = size(Data_D.data);

% Preallocating space for matrices with PermNo and Dates
PermNo           = NaN(NRows, 1);
Dates            = NaN(NRows, 1);

for i = 1:NRows
    % Note: the str2num function converts text to data
    PermNo(i,1)  = str2num(Data_D.textdata{i+1, 1});
    Dates( i,1)  = str2num(Data_D.textdata{i+1, 2});
end

% Collecting data in a single matrix
DataNum          = [Dates, PermNo, Data_D.data];

% Getting unique values of stock identifiers and dates
PermNoUnique     = unique(PermNo);
DatesUnique      = unique(Dates);

% Computing number of stocks and dates
NStocks          = length(PermNoUnique);
NDates           = length(DatesUnique);

%% *** Part 3 - Computing total returns and total returns index ***

% Allocate memory for total returns for each stock
TR_Stocks        = NaN(NDates, NStocks);

for i = 1:NStocks
   % Extract the rows that correspond to that stock
   Indices       = PermNo==PermNoUnique(i);
   StockData     = DataNum(Indices, :);
   
   % Sort the data according to date
   StockData     = sortrows(StockData);

   % Store total returns for this stock in matrix
   TR_Stocks(:, i) = StockData(:, 3);
   
   % Store total returns for S&P 500
   if i == 1
       TR_SPX      = StockData(:, 4);
   end
end

% Compiling single matrix with all total returns data
TR                 = [TR_Stocks, TR_SPX];  

% Computing log returns
LogReturns         = log(1 + TR);

%% *** Part 4: Compute a EWMA volatility for each series ***
Lambda             = 0.97;
Window             = 253;
LogReturnsSq       = LogReturns .* LogReturns;

% Will initialize the computation using two methods, so construct two
% matrices for the estimates: EWMA contains unweighted initial estimate,
% EWMAWgt contains the weighted one
% Note that initially the matrices will contain variances; will convert
% into volatilities at the last step
EWMA               = NaN(NDates - Window, NStocks + 1);
EWMAWgt            = NaN(NDates - Window, NStocks + 1);

% Initial unweighted estimate
EWMA(1, :)         = var(LogReturns(1:Window, :));

% Initial weighted estimate
% A direct computation is possible, but the var function contains an
% option for weighting the observations
Weights            = exp((Window:(-1):1)' * log(Lambda));
Weights            = Weights / sum(Weights);
EWMAWgt(1, :)      = var(LogReturns(1:Window, :), Weights);

% Loop to compute the estimates for the entire series
for i=2:(NDates - Window)
    EWMA(i, :)     = Lambda * EWMA(i-1, :) + (1-Lambda) * ...
                        LogReturnsSq(i + Window - 1, :);
    EWMAWgt(i, :)  = Lambda * EWMAWgt(i-1, :) + (1-Lambda) * ...
                        LogReturnsSq(i + Window - 1, :);
end

% Finally, convert into volatility estimates#
EWMAVol            = sqrt(EWMA);
EWMAVolWgt         = sqrt(EWMAWgt);

%% *** Part 5: Estimate a GARCH(1, 1) model for each series ***

% Create a model first
GARCHModel         = garch(1, 1);

% Estimate GARCH for each series, save information for computing estimates
GARCHInfo          = NaN(4, NStocks + 1);

for i=1:(NStocks + 1)
   EstimatedModel  = estimate(GARCHModel, LogReturns(:, i));
   GARCHInfo(1, i) = EstimatedModel.Constant;
   GARCHInfo(2, i) = EstimatedModel.ARCH{1};
   GARCHInfo(3, i) = EstimatedModel.GARCH{1};
   GARCHInfo(4, i) = EstimatedModel.UnconditionalVariance;
end

% When we inspect the results (or, alternatively, look at the error
% messages), we notice that there are inaccuracies for JPM and C.
% For those stocks, we re-estimate using ugarch
% Note that ugarch inverts the meaning of alpha and beta relative to
% our usage
Indices         = [4 6];
for i = Indices
    [Omega , Alpha , Beta] = ugarch(LogReturns(:, i) , 1 , 1); 
    GARCHInfo(1, i) = Omega;
    GARCHInfo(2, i) = Beta;                 % ARCH term in ugarch function
    GARCHInfo(3, i) = Alpha;                % GARCH term in ugarch function
    GARCHInfo(4, i) = Omega / (1 - Alpha - Beta);
end

% With the coefficient estimates, we can compute GARCH variances
GARCHVar            = NaN(NDates, NStocks + 1);
GARCHVar(1, :)      = GARCHInfo(4, :);

for i=2:NDates
    GARCHVar(i, :)  = GARCHInfo(1, :) + ...
                         GARCHInfo(2, :) .* LogReturnsSq(i - 1, :) + ...
                         GARCHInfo(3, :) .* GARCHVar(i - 1, :);
end

% Finally, convert to volatility, eliminate first year of observations
GARCHVol            = sqrt(GARCHVar((Window+1):end, :));

%% *** Part 6: Produce a chart with log returns and +/- 2 Std. Dev. ***
%  ***         for the various models                               ***
LogReturnsAux       = LogReturns((Window+1):end, :);

% Choose security to chart
i                   = 4;

% Chart options
FigNo               = 1;
ChartData             = [LogReturnsAux(:, i) -2*EWMAVolWgt(:, i) ...
                         2*EWMAVolWgt(:, i) -2*GARCHVol(:, i) ...
                         2*GARCHVol(:, i)];
StartDate             = 1990;
EndDate               = 2015; 
TimeLabels            = linspace(StartDate, EndDate, NDates - Window);

% Charting the series
figure(FigNo)
titleStr = ['Daily Returns for ' Tickers(i, :)];
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
                 'LineWidth', 1, 'Color', 'green');
h3 = plot(TimeLabels, ChartData(:, 3), 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'green'); 
h4 = plot(TimeLabels, ChartData(:, 4), 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
h5 = plot(TimeLabels, ChartData(:, 5), 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');  
set(gca, 'YTickMode', 'manual');
set(gca, 'YTickLabel', num2str(100 .* get(gca, 'YTick')', '%1.0f%%'));              
legend1 = legend([h1 h3 h5], {'Log Returns', 'EWMA Vol', 'GARCH Vol'}, ...
                     'location','best');            
set(gcf, 'color', 'white');                  
FigNo   = FigNo + 1;

%% *** Part 7: Exploring the effect of varying p and q ***
% Will work primarily with the S&P 500
SPXReturns       = LogReturns(:, 7);

% Run GARCH with p and q varying from 1 to 4, compute log-likelihood
% LogL contains the value of p in the first column, q in the second,
% and we'll store the log likelihood in the third
LogLVal          = NaN(16, 3);
LogLVal(:, 1)    = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4]';
LogLVal(:, 2)    = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4]';

for i=1:16
    p            = LogLVal(i, 1);
    q            = LogLVal(i, 2);
    GARCHModel   = garch(p, q);
    [~, ~, logL, ~] = estimate(GARCHModel, SPXReturns);
    LogLVal(i,3) = logL;
end

% Look at some individual models
GARCH_44         = garch(4, 4);
[EstMdl, EstParamCov, logL, info] = estimate(GARCH_44, SPXReturns);

GARCH_12         = garch(1, 2);
[~, ~, logL, ~]  = estimate(GARCH_12, JPM);
LogL_UNC         = logL;

GARCH_11         = garch(1, 1);
[~, ~, logL, ~]  = estimate(GARCH_11, JPM);
LogL_CONST       = logL;

% Perform likelihood ratio test
LR               = -2*(LogL_CONST - LogL_UNC);
PVal             = 1 - chi2cdf(LR, 1);

%% *** Part 8: Estimating GJR-GARCH Model ***
GJRGARCH         = NaN(4, NStocks + 1);

for i=1:(NStocks + 1)
   % Note: parameters are, in order, omega, alpha, gamma and beta
   parameters    = tarch(LogReturns(:, i), 1, 1, 1);
   GJRGARCH(:,i) = parameters';
end

% Extracting parameters
Omega            = GJRGARCH(1, :);
Alpha            = GJRGARCH(2, :);
Gamma            = GJRGARCH(3, :);
Beta             = GJRGARCH(4, :);

% Computing volatility estimates for each series
% First, allocate memory for variance estimates
GJRGARCHVar      = NaN(NDates, NStocks + 1);

% Will use the first year as a burn-in period, so volatility estimates
% are not sensitive to how we estimate the first-period variance

% Can use squared returns or GARCH-like estimates (our choice)
GJRGARCHVar(1, :)= Omega ./ (ones(1, 7) - Alpha - Beta);

% Creating matrix with the sign of returns
ReturnsIsNeg     = ones(size(LogReturns));
ReturnsIsNeg(LogReturns > 0) = 0;

% Loop to compute variance estimates
for i=2:NDates
    GJRGARCHVar(i, :) = Omega + (Alpha + Gamma .* ReturnsIsNeg(i-1, :)) ...
                        .* LogReturnsSq(i-1, :) + Beta .* ...
                        GJRGARCHVar(i-1, :);
end

% Compute volatility, discard first 253 observations
GJRGARCHVol      = sqrt(GJRGARCHVar);
GJRGARCHVol(1:253, :) = [];

% Produce chart for a specific security
i                = 7;
ChartData        = [LogReturnsAux(:, i) -2*GJRGARCHVol(:, i) ...
                         2*GJRGARCHVol(:, i)];
                     
% Charting the series
figure(FigNo)
titleStr = ['Daily Returns for ' Tickers(i, :)];
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
h3 = plot(TimeLabels, ChartData(:, 3), 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');            
legend1 = legend([h1 h2], {'Log Returns', '+/- 2 GJR-GARCH Vol'}, ...
                     'location','best'); 
set(gca, 'YTickMode', 'manual');
set(gca, 'YTickLabel', num2str(100 .* get(gca, 'YTick')', '%1.0f%%'))                               
set(gcf, 'color', 'white');  