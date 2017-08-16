%% London School of Economics
%  FM442 - Michaelmas Term 2016
%  Author: Domingos Romualdo
%  Date: 19 November 2016
%  Assignment 6 - Suggested Solution

%% *** Part 0 - Preliminaries ***
clc;
clear all;
clf;
close all;
FigNo            = 1;

%% *** Part 1 - Loading daily data ***
Data_D = importdata('FM442 - HW2 - Daily Data.csv');
Tickers          = ['MSFT'; 'XOM '; 'GE  '; 'JPM '; 'INTC'; 'C   '; ...
                       'SPX '];
                   
%% *** Part 2 - Cleaning up and organizing the data ***
% Retrieving information about data
[NRows, NCols] = size(Data_D.data);

% Preallocating space for matrices with PermNo and Dates
PermNo = NaN(NRows, 1);
Dates  = NaN(NRows, 1);

for i = 1:NRows
    % Note: the str2num function converts text to data
    PermNo(i,1) = str2num(Data_D.textdata{i+1, 1});
    Dates( i,1) = str2num(Data_D.textdata{i+1, 2});
end

% Collecting data in a single matrix
DataNum       = [Dates, PermNo, Data_D.data];

% Getting unique values of stock identifiers and dates
PermNoUnique  = unique(PermNo);
DatesUnique   = unique(Dates);

% Computing number of stocks and dates
NStocks       = length(PermNoUnique);
NDates        = length(DatesUnique);

%% *** Part 3 - Computing total returns and total returns index ***

% Allocate memory for total returns for each stock
TR_Stocks        = NaN(NDates, NStocks);

for i = 1:NStocks
   % Extract the rows that correspond to that stock
   Indices       = find(PermNo==PermNoUnique(i));
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

% Computing total returns indices
TRIndices          = NaN(NDates + 1, NStocks + 1);
TRIndices(1, :)    = ones(1, NStocks + 1);

for i=1:NDates
    TRIndices(i+1, :) = TRIndices(i, :) .* (1 + TR(i, :));
end

% Computing log returns
LogReturns         = log(1 + TR);

%% *** Part 4 - Computing the dates for which we will estimate model ***
% We will update GARCH(1, 1) coefficients at the end of every month
% We start by computing a vector with the indices of the dates when we
% we will reestimate the model

% Computing day, month and year for each date in sample
Date_Day           = mod(DatesUnique, 100);
AuxVec             = (DatesUnique - Date_Day)/100;
Date_Month         = mod(AuxVec, 100);
Date_Year          = (AuxVec - Date_Month)/100;

% Determining on which dates we will reestimate the model
Reestimate         = ones(NDates, 1);   % Initialize vector

% Set burn-in period length; don't estimate model before it passes
Window             = 253;               % Burn-in period
Reestimate(1:(Window-1), 1) = 0;

% Note: if daily reestimation is desired, skip the rest of this section

% The following line determines which rows do not correspond to 
% end-of-month dates; we won't estimate the model in those dates
NotMonthEndIndices = (Date_Month(1:(end-1),:) == Date_Month(2:end,:));
Reestimate(NotMonthEndIndices, :) = 0;

% Also, there is no need to reestimate on the last day
Reestimate(end, 1) = 0;

% If we want to re-estimate the model only quarterly, run the next line;
% otherwise, reestimation will be monthly
Reestimate(mod(Date_Month, 3) ~= 0, 1) = 0;

% Checking estimation dates
EstimationIndices  = find(Reestimate == 1);
EstimationDates    = DatesUnique(EstimationIndices);

%% *** Part 5 - Estimating GARCH model with increasing samples ***
% Choosing stock, setting parameters
StockIndex         = 5;               % Using INTC
StockData          = LogReturns(:, StockIndex);
StockDataSq        = StockData .* StockData;
ConditionalVar     = NaN(NDates, 1);  % To store conditional variance
                                      %    estimates
GARCHModel         = garch(1, 1);
EstDateCount       = size(EstimationIndices, 1);
GARCHCoeffs        = NaN(EstDateCount, 4);  % To store coefficients for
                                            % later examination

% Loop to reestimate the model and compute conditional variances
for i=1:EstDateCount
    % Extract the data for model estimation
    SampleEnd      = EstimationIndices(i,1);
    GARCHData      = StockData(1:SampleEnd, :);
    
    % Estimate model
    EstimatedModel  = estimate(GARCHModel, GARCHData);
        
    % Extract coefficients
    Omega           = EstimatedModel.Constant;
    Alpha           = EstimatedModel.ARCH{1};
    Beta            = EstimatedModel.GARCH{1};
    UncondVar       = EstimatedModel.UnconditionalVariance;
    
    % Store coefficients
    GARCHCoeffs(i,1)= Omega;
    GARCHCoeffs(i,2)= Alpha;
    GARCHCoeffs(i,3)= Beta;
    GARCHCoeffs(i,4)= UncondVar;
    
    % Determine last date for which we will use forecasts from this
    % estimate
    if (i < size(EstimationIndices, 1))
        LastDate    = EstimationIndices(i+1, 1);
    else
       LastDate     = NDates; 
    end
    
    % Produce conditional variance estimates
    VarianceEst     = NaN(LastDate, 1);
    VarianceEst(1, 1) = UncondVar;
    
    for t=2:LastDate
        VarianceEst(t, 1) = Omega + Alpha * StockDataSq(t - 1, 1) ...
                               + Beta * VarianceEst(t - 1, 1); 
    end
    
    % We will use the variance estimates from this estimation from 
    % SampleEnd + 1 to LastDate
    ConditionalVar((SampleEnd+1):LastDate, :) = ...
        VarianceEst((SampleEnd+1):LastDate, :);
end

% Deleting the observations within the window, computing volatility
ConditionalVol      = sqrt(ConditionalVar((Window+1):end, :));

% Code to save and to load 'ConditionalVol' and 'GARCHCoeffs'
save('ConditionalVol.mat', 'ConditionalVol');
save('GARCHCoeffs.mat', 'GARCHCoeffs');
load('ConditionalVol');
load('GARCHCoeffs');

%% *** Part 6 - Computing VaR estimates for the chosen stock, testing ***
p                   = 0.01;            % VaR level
InvCDF              = norminv(p);
ValueAtRiskEst      = - InvCDF * ConditionalVol;
Breach = (StockData((Window+1):end, :) < - ValueAtRiskEst);
LenBreach           = size(Breach, 1);
pHat                = sum(Breach)/size(Breach, 1);

% Unconditional coverage test
V1                  = sum(Breach);
V0                  = size(Breach, 1) - V1;

LogLConst           = V1 * log(p) + V0 * log(1-p);
LogLUnconst         = V1 * log(pHat) + V0 * log(1-pHat);
TestStatisticUC     = -2 * (LogLConst - LogLUnconst);
SignificanceUC      = 1 - chi2cdf(TestStatisticUC, 1);

% Conditional coverage test
% Estimating the transition probabilities
Aux01               = (ones(LenBreach - 1, 1) - Breach(1:(end -1), :)) ...
                         .* Breach(2:end, :);
p01                 = sum(Aux01) / (LenBreach - 1 ...
                         - sum(Breach(1:(end -1), :)));
                     
Aux11               = Breach(1:(end -1), :) .* Breach(2:end, :);
p11                 = sum(Aux11) / sum(Breach(1:(end -1), :));

% Computing log-likelihood; note that the constrained version here is 
% equal to the unconstrained version above
LogLConst           = LogLUnconst;

% Running loop to compute the unconstrained log-likelihood
LogLUnconst         = Breach(1, 1) * log(pHat) + ...
                         (1 - Breach(1, 1)) * log(1 - pHat);
for i=2:size(Breach, 1)
    if     (Breach(i-1, 1) == 0 && Breach(i, 1) == 0)
        LogLUnconst = LogLUnconst + log(1 - p01);
    elseif (Breach(i-1, 1) == 0 && Breach(i, 1) == 1)
        LogLUnconst = LogLUnconst + log(p01);
    elseif (Breach(i-1, 1) == 1 && Breach(i, 1) == 0)
        LogLUnconst = LogLUnconst + log(1 - p11);
    elseif (Breach(i-1, 1) == 1 && Breach(i, 1) == 1)
        LogLUnconst = LogLUnconst + log(p11);
    end
end

TestStatisticCC     = -2 * (LogLConst - LogLUnconst);
SignificanceCC      = 1 - chi2cdf(TestStatisticCC, 1);

% Extra item: computing expected shortfall (in standard deviation units)
ExpectedES          = normpdf(norminv(p))/p;

% Sample value of expected shortfall
ActualShortfall     = - Breach .* StockData((Window+1):end, :) ...
                         ./ ConditionalVol;
SampleES            = sum(ActualShortfall) / sum(Breach);

%% *** Part 7 - Choosing stocks for DCC, computing univariate GARCH ***
Indices             = [1 5];    % Choose MSFT and INTC
GARCHData           = LogReturns(:, Indices);
GARCHDataSq         = GARCHData .* GARCHData;
GARCHCoeffs         = NaN(4, 2);

for i=1:2
    % Estimate model
    EstimatedModel  = estimate(GARCHModel, GARCHData(:, i));
        
    % Extract coefficients
    GARCHCoeffs(1,i)= EstimatedModel.Constant;
    GARCHCoeffs(2,i)= EstimatedModel.ARCH{1};
    GARCHCoeffs(3,i)= EstimatedModel.GARCH{1};
    GARCHCoeffs(4,i)= EstimatedModel.UnconditionalVariance;
end

% Compute conditional volatility estimates
GARCHVar            = NaN(NDates, 2);
GARCHVar(1, :)      = GARCHCoeffs(4, :);

for t=2:NDates
   GARCHVar(t, :)   = GARCHCoeffs(1, :) ...
                         + GARCHCoeffs(2,i) .* GARCHDataSq(t-1, :) ...
                         + GARCHCoeffs(3,i) .* GARCHVar(t-1, :);
end

GARCHVol            = sqrt(GARCHVar((Window+1):end, :));
MSFTGARCHVol        = GARCHVol(:, 1);
INTCGARCHVol        = GARCHVol(:, 2);

% Computing rolling 20-day correlations
UnivCorr            = NaN(NDates, 1);

for t=(Window+1):NDates
    Data1           = GARCHData((t-20):(t-1), 1);
    Data2           = GARCHData((t-20):(t-1), 2);
    CorrMat         = corrcoef(Data1, Data2);
    UnivCorr(t, 1)  = CorrMat(1, 2);
end

UnivCorr(1:Window,:)= [];

%% *** Part 8 - Running DCC model for two stocks iteratively ***
DCCBacktestVar      = NaN(2, 2, NDates);
DCCData             = LogReturns(:, Indices);
DCCDataSq           = DCCData .* DCCData;

% Loop to reestimate the DCC model and compute conditional variance
% matrices.  The code will combine the loop structure above with the 
% DCC code from Homework 4

% Note: when running the loop below with monthly reestimation, it took
% about one hour and ten minutes on my home PC 
for i=1:size(EstimationIndices, 1)
    % Extract the data for model estimation
    SampleEnd       = EstimationIndices(i,1);
    if (i < size(EstimationIndices, 1))
        FcstEnd     = EstimationIndices(i+1, 1);
    else
        FcstEnd    = NDates; 
    end

    MVGARCHData     = DCCData(1:SampleEnd, :);
    FcstData        = DCCData(1:FcstEnd, :);
    
    % DCC estimation, followed by parameter extraction
    [DCCParms]      = dcc(MVGARCHData, [], 1, 0, 1);

    % Construct DCC model matrices
    DCCOmega        = [DCCParms(1) DCCParms(4)];
    DCCAlpha        = [DCCParms(2) DCCParms(5)];
    DCCBeta         = [DCCParms(3) DCCParms(6)];
    DCCRBar         = [1 DCCParms(7); DCCParms(7) 1];
    DCC_a           = DCCParms(8);
    DCC_b           = DCCParms(9);

    % Computing GARCH volatilities for each asset and standardized residuals
    DCCVar          = NaN(FcstEnd, 2);
    DCCVar(1, :)    = DCCOmega ./ (1 - DCCAlpha - DCCBeta);

    for t=2:FcstEnd
        DCCVar(t, :) = DCCOmega + ...
                             DCCAlpha .* DCCDataSq(t - 1, :) + ...
                             DCCBeta .* DCCVar(t - 1, :);
    end

    DCCVol           = sqrt(DCCVar);
    DCCEps           = FcstData ./ DCCVol;

    % Compute value for DCC Qt and Rt, using recursive relation; 
    % initialize with estimate for RBar
    DCCQt            = NaN(2, 2, FcstEnd);
    DCCRt            = NaN(2, 2, FcstEnd);
    DCCQt(:, :, 1)   = DCCRBar;
    DCCRt(:, :, 1)   = DCCRBar;

    for t=2:FcstEnd
        DCCQt(:, :, t) = (1 - DCC_a - DCC_b) * DCCRBar ...
                           + DCC_a * DCCEps(t-1, :)' * DCCEps(t-1, :) ...
                           + DCC_b * DCCQt(:, :, t-1);
        AuxMat         = [sqrt(DCCQt(1, 1, t))  sqrt(DCCQt(2, 2, t))];
        DCCRt(:, :, t) = DCCQt(:, :, t) ./ (AuxMat' * AuxMat);
        
        % For dates past the end of the sample, store the variance matrix
        if (t > SampleEnd)
            DMat          = [DCCVol(t, 1) 0; 0 DCCVol(t, 2)];
            DCCBacktestVar(:, :, t) = DMat * DCCRt(:, :, t) * DMat;
        end
    end
end

% Code to save and load DCCBacktestVar
save('DCCBacktestVar.mat', 'DCCBacktestVar');
load('DCCBacktestVar');

% Compute historical volatilities and correlations
DCCVolMSFT         = sqrt(DCCBacktestVar(1, 1, (Window+1):end));
DCCVolINTC         = sqrt(DCCBacktestVar(2, 2, (Window+1):end));
DCCCovariance      = DCCBacktestVar(1, 2, (Window+1):end);
DCCCorrel          = DCCCovariance ./ (DCCVolMSFT .* DCCVolINTC);

% Reshape vectors
DCCVolMSFT         = squeeze(DCCVolMSFT);
DCCVolINTC         = squeeze(DCCVolINTC);
DCCCorrel          = squeeze(DCCCorrel);

%% *** Part 9 - Comparing out-of-sample DCC backtest with  ***
%  *** in-sample univariate GARCH and rolling correlations ***                            
% Produce charts comparing output of BEKK with GARCH information
% Creating time labels
StartDate             = 1991;
EndDate               = 2015; 
TimeLabels            = linspace(StartDate, EndDate, NDates - Window);

% Produce figure
figure(FigNo)
subplot(3, 1, 1);
titleStr = 'Volatility Estimates for MSFT';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate 0 1]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * MSFTGARCHVol, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, sqrt(252) * DCCVolMSFT, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('Univariate GARCH', 'Out-of-Sample DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');

subplot(3, 1, 2);
titleStr = 'Volatility Estimates for INTC';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate 0 1.5]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * INTCGARCHVol, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, sqrt(252) * DCCVolINTC, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('Univariate GARCH', 'Out-of-Sample DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');

subplot(3, 1, 3);
titleStr = 'Correlation Estimates';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate -1.2 1.2]);
grid on;
hold on;
plot(TimeLabels, UnivCorr, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, DCCCorrel, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('20-Day', 'Out-of-Sample DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Correlation');
FigNo = FigNo + 1;

%% *** Part 10 - Computing portfolio volatility, charting it ***
Weights             = [0.5; 0.5];
PortRet             = DCCData * Weights;
PortLogRet          = log(1 + PortRet);

% Computing conditional portfolio variance
PortCondVar         = NaN(NDates - Window, 1);
for i=(Window+1):NDates
   CurrentVarMat    = DCCBacktestVar(:, :, i);
   PortCondVar(i,:) = Weights' * CurrentVarMat * Weights;
end

% Computing conditional standard deviation, dropping initial window
ConditionalVol      = sqrt(PortCondVar((Window+1):end, :));

% Producing charts
figure(FigNo)
titleStr = 'Portfolio Returns +/- 2 Std. Dev. (Out-of-Sample DCC)';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate -0.15 0.15]);
grid on;
hold on;
h1 = plot(TimeLabels, PortLogRet((Window+1):end, :), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
h2 = plot(TimeLabels, 2 * ConditionalVol, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
h3 = plot(TimeLabels, (-2) * ConditionalVol, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');  
legend([h1 h3], {'Log Returns', '+/- 2 DCC Vol'}, ...
     'location','best'); 
set(gcf, 'color', 'white'); 

%% *** Part 11 - Performing VaR analysis with portfolio ***
% Parameters for VaR backtest
p                   = 0.01;            % VaR level
InvCDF              = norminv(p);
ValueAtRiskEst      = - InvCDF * ConditionalVol;

% Finding breaches
Breach              = (PortLogRet((Window+1):end, :) < - ValueAtRiskEst);
LenBreach           = size(Breach, 1);
pHat                = sum(Breach)/size(Breach, 1);

% Unconditional coverage test
V1                  = sum(Breach);
V0                  = size(Breach, 1) - V1;

LogLConst           = V1 * log(p) + V0 * log(1-p);
LogLUnconst         = V1 * log(pHat) + V0 * log(1-pHat);
TestStatisticUC     = -2 * (LogLConst - LogLUnconst);
SignificanceUC      = 1 - chi2cdf(TestStatisticUC, 1);

% Conditional coverage test
% Estimating the transition probabilities
Aux01               = (ones(LenBreach - 1, 1) - Breach(1:(end - 1), :)) ...
                         .* Breach(2:end, :);
p01                 = sum(Aux01) / (LenBreach - 1 ...
                         - sum(Breach(1:(end -1), :)));
                     
Aux11               = Breach(1:(end - 1), :) .* Breach(2:end, :);
p11                 = sum(Aux11) / sum(Breach(1:(end -1), :));

% Computing log-likelihood; note that the constrained version here is 
% equal to the unconstrained version above
LogLConst           = LogLUnconst;

% Running loop to compute the unconstrained log-likelihood
LogLUnconst         = Breach(1, 1) * log(pHat) + ...
                         (1 - Breach(1, 1)) * log(1 - pHat);
for i=2:size(Breach, 1)
    if     (Breach(i-1, 1) == 0 && Breach(i, 1) == 0)
        LogLUnconst = LogLUnconst + log(1 - p01);
    elseif (Breach(i-1, 1) == 0 && Breach(i, 1) == 1)
        LogLUnconst = LogLUnconst + log(p01);
    elseif (Breach(i-1, 1) == 1 && Breach(i, 1) == 0)
        LogLUnconst = LogLUnconst + log(1 - p11);
    elseif (Breach(i-1, 1) == 1 && Breach(i, 1) == 1)
        LogLUnconst = LogLUnconst + log(p11);
    end
end

TestStatisticCC     = -2 * (LogLConst - LogLUnconst);
SignificanceCC      = 1 - chi2cdf(TestStatisticCC, 1);

% Extra item: computing expected shortfall (in standard deviation units)
ExpectedES          = normpdf(norminv(p))/p;

% Sample value of expected shortfall
ActualShortfall     = - Breach .* PortLogRet((Window+1):end, :) ...
                         ./ ConditionalVol;
SampleES            = sum(ActualShortfall) / sum(Breach);