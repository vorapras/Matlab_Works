%% London School of Economics
%  FM442 - Michaelmas Term 2016
%  Author: Domingos Romualdo
%  Date: 8 November 2016
%  Assignment 4 - Suggested Solution

% Preliminaries
clc;
clear all;
clf;
close all;

%% *** Part 1: Loading daily data ***
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

%% *** Part 4: Choose stocks, estimate univariate GARCH as benchmark ***
% Choose stocks to estimate the model
% Using MSFT and INTC
Indices          = [1 5];
MVGARCHData      = LogReturns(:, Indices);

% Estimate GARCH for each series, save information for computing estimates
GARCHModel         = garch(1, 1);
GARCHInfo          = NaN(4, 2);

for i=1:2
   EstimatedModel  = estimate(GARCHModel, MVGARCHData(:, i));
   GARCHInfo(1, i) = EstimatedModel.Constant;
   GARCHInfo(2, i) = EstimatedModel.ARCH{1};
   GARCHInfo(3, i) = EstimatedModel.GARCH{1};
   GARCHInfo(4, i) = EstimatedModel.UnconditionalVariance;
end

% With the coefficient estimates, we can compute GARCH variances
GARCHVar            = NaN(NDates, 2);
GARCHVar(1, :)      = GARCHInfo(4, :);
LogReturnsSq        = MVGARCHData .* MVGARCHData;

for i=2:NDates
    GARCHVar(i, :)  = GARCHInfo(1, :) + ...
                         GARCHInfo(2, :) .* LogReturnsSq(i - 1, :) + ...
                         GARCHInfo(3, :) .* GARCHVar(i - 1, :);
end

% Convert to volatility, eliminate first year of observations for later
% comparison and charting
Window              = 253;
GARCHVol            = sqrt(GARCHVar((Window+1):end, :));
GARCHVolMSFT        = GARCHVol(:, 1);
GARCHVolINTC        = GARCHVol(:, 2);

% Compute rolling correlations for comparison with models, dropping
% first year
CorrWindow          = 20;
MWCorr              = NaN(NDates - Window, 1);

for i=1:(NDates-Window)
    Indices         = (Window+i-1-CorrWindow):(Window+i-2);
    Data            = MVGARCHData(Indices,:);
    Temp            = corrcoef(Data);
    MWCorr(i, 1)    = Temp(1, 2);
end

%% *** Part 5: BEKK with stock data ***
% Estimate BEKK(1, 1) - about 10 minutes on my personal machine
[BEKKParms]      = bekk(MVGARCHData, [], 1, 0, 1, 'Full');

% Section to save and retrieve BEKKParms, if needed
save('BEKKParms.mat', 'BEKKParms');
load('BEKKParms');

% Convert parameters to usual format
[BEKK_C, BEKK_A, ~, BEKK_B] = bekk_parameter_transform(BEKKParms, ...
                      1, 0, 1, 2, 3);

% Compute variance matrix estimates at each point in time -
% use first 253 observations to initialize computation.

% The following matrix will hold the estimates; first row is for the
% initial period
BEKKVarMat         = NaN(2, 2, NDates);
BEKKVarMat(:,:, 1) = cov(MVGARCHData(1:Window, :));

for i=2:NDates
    LogRet              = MVGARCHData(i - 1, :);
    BEKKVarMat(:,:,i) = BEKK_C ...
          + BEKK_A' * (LogRet' * LogRet) * BEKK_A ...
          + BEKK_B' * BEKKVarMat(:,:,i-1) * BEKK_B;
end

% Delete first year of estimates
BEKKVarMat         = BEKKVarMat(:, :, (Window+1):end);

% Extract variance estimates for each stock, as well as correlations
BEKKVolMSFT        = sqrt(BEKKVarMat(1, 1, :));
BEKKVolINTC        = sqrt(BEKKVarMat(2, 2, :));
BEKKCorr           = BEKKVarMat(1, 2, :) ./ (BEKKVolMSFT .* BEKKVolINTC);

% Reshape matrices into vectors
BEKKVolMSFT        = reshape(BEKKVolMSFT, NDates - Window, 1);
BEKKVolINTC        = reshape(BEKKVolINTC, NDates - Window, 1);
BEKKCorr           = reshape(BEKKCorr, NDates - Window, 1);

% Produce charts comparing output of BEKK with GARCH information
% Creating time labels
StartDate             = 1991;
EndDate               = 2015; 
TimeLabels            = linspace(StartDate, EndDate, NDates - Window);

% Produce figure
FigNo                 = 1;

figure(FigNo)
subplot(3, 1, 1);
titleStr = 'Volatility Estimates for MSFT';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate 0 1]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * GARCHVolMSFT, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, sqrt(252) * BEKKVolMSFT, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('GARCH', 'BEKK', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');

subplot(3, 1, 2);
titleStr = 'Volatility Estimates for INTC';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate 0 1.5]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * GARCHVolINTC, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, sqrt(252) * BEKKVolINTC, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('GARCH', 'BEKK', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');

subplot(3, 1, 3);
titleStr = 'Correlation Estimates';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate -1.2 1.2]);
grid on;
hold on;
plot(TimeLabels, MWCorr, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, BEKKCorr, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('20-Day', 'BEKK', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');
FigNo = FigNo + 1;

%% *** Part 6: CCC with stock data ***
% Estimate CCC model
[CCCParms,LL,HT,VCV,SCORES]  = ccc_mvgarch(MVGARCHData, [], 1, 0, 1);
CCCOmega         = [CCCParms(1) CCCParms(4)];
CCCAlpha         = [CCCParms(2) CCCParms(5)];
CCCBeta          = [CCCParms(3) CCCParms(6)];
CCCRho           = CCCParms(7);

% Constructing volatility estimates from GARCH Parameters
CCCVar           = NaN(NDates, 2);
CCCVar(1, :)     = CCCOmega ./ (1 - CCCAlpha - CCCBeta);

for i=2:NDates
    CCCVar(i, :) = CCCOmega + ...
                         CCCAlpha .* LogReturnsSq(i - 1, :) + ...
                         CCCBeta .* GARCHVar(i - 1, :);
end

% Convert to volatility, drop first year
CCCVol           = sqrt(CCCVar((Window+1):end,:));
CCCVolMSFT       = CCCVol(:, 1);
CCCVolINTC       = CCCVol(:, 2);

% Chart volatilities
figure(FigNo)
subplot(2, 1, 1);
titleStr = 'Volatility Estimates for MSFT';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate 0 1]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * GARCHVolMSFT, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, sqrt(252) * CCCVolMSFT, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('GARCH', 'CCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');

subplot(2, 1, 2);
titleStr = 'Volatility Estimates for INTC';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate 0 1.5]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * GARCHVolINTC, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, sqrt(252) * CCCVolINTC, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('GARCH', 'CCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');
FigNo            = FigNo + 1;

%% *** Part 7: DCC with stock data ***
% Estimate DCC model
[DCCParms]       = dcc(MVGARCHData, [], 1, 0, 1);

% Construct DCC model matrices
DCCOmega         = [DCCParms(1) DCCParms(4)];
DCCAlpha         = [DCCParms(2) DCCParms(5)];
DCCBeta          = [DCCParms(3) DCCParms(6)];
DCCRBar          = [1 DCCParms(7); DCCParms(7) 1];
DCC_a            = DCCParms(8);
DCC_b            = DCCParms(9);

% Computing GARCH Volatilities for each asset, and standardized residuals
DCCVar           = NaN(NDates, 2);
DCCVar(1, :)     = DCCOmega ./ (1 - DCCAlpha - DCCBeta);

for i=2:NDates
    DCCVar(i, :) = DCCOmega + ...
                         DCCAlpha .* LogReturnsSq(i - 1, :) + ...
                         DCCBeta .* GARCHVar(i - 1, :);
end

DCCVol           = sqrt(DCCVar);
DCCEps           = MVGARCHData ./ DCCVol;

% Compute value for DCC Qt and Rt, using recursive relation; 
% initialize with estimate for RBar
DCCQt            = NaN(2, 2, NDates);
DCCRt            = NaN(2, 2, NDates);
DCCQt(:, :, 1)   = DCCRBar;
DCCRt(:, :, 1)   = DCCRBar;

for i=2:NDates
    DCCQt(:, :, i) = (1 - DCC_a - DCC_b) * DCCRBar ...
                       + DCC_a * DCCEps(i-1, :)' * DCCEps(i-1, :) ...
                       + DCC_b * DCCQt(:, :, i-1);
    AuxMat         = [sqrt(DCCQt(1, 1, i))  sqrt(DCCQt(2, 2, i))];
    DCCRt(:, :, i) = DCCQt(:, :, i) ./ (AuxMat' * AuxMat);
end

% Extracting correlation information, deleting first year
DCCCorr            = DCCRt(1, 2, :);
DCCCorr            = reshape(DCCCorr, NDates, 1);

% Putting together information for chart
DCCVolMSFT         = DCCVol((Window+1):end, 1);
DCCVolINTC         = DCCVol((Window+1):end, 2);
DCCCorr            = DCCCorr((Window+1):end, :);

figure(FigNo)
subplot(3, 1, 1);
titleStr = 'Volatility Estimates for MSFT';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate 0 1]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * GARCHVolMSFT, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, sqrt(252) * DCCVolMSFT, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('GARCH', 'DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');

subplot(3, 1, 2);
titleStr = 'Volatility Estimates for INTC';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate 0 1.5]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * GARCHVolINTC, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, sqrt(252) * DCCVolINTC, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('GARCH', 'DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');

subplot(3, 1, 3);
titleStr = 'Correlation Estimates';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate -1.2 1.2]);
grid on;
hold on;
plot(TimeLabels, MWCorr, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, DCCCorr, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('20-Day', 'DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');
FigNo = FigNo + 1;


%% *** Part 8: Principal components analysis with stock data ***
% Import currency excess returns data
Data_Ccy         = importdata('Currency Excess Returns.csv');
Tickers          = Data_Ccy.textdata(1, :);
Tickers(:, 1)    = [];

% Extracting currency excess returns data from the import
CcyRet           = Data_Ccy.data(:, 2:end);
NumMonths        = size(CcyRet, 1);
CcyLogRet        = log(1 + CcyRet);

% Performing principal components analysis
[Weights, Princomp, Eigenvals] = pca(CcyLogRet);
FactorVariance   = sum(Eigenvals);
TotalVar         = sum(FactorVariance);
FracVarExplained = FactorVariance / TotalVar;

% Displaying principal component weights on each currency
Labels           = num2cell(Tickers);
DisplayMat1      = [Labels' num2cell(Weights')];
DisplayMat2      = [num2str('Factor Variance') num2cell(FactorVariance)];
DisplayMat3      = [num2str('Fraction Total') num2cell(FracVarExplained)];
DisplayMat       = [DisplayMat1; DisplayMat2; DisplayMat3];

% Loading interest rate data
Data_IntRate     = importdata('Interest Rates.csv');
IntRates         = Data_IntRate.data;
IntRates(:, 1)   = [];
IntRateDiff      = IntRates(:,1:(end-1)) - repmat(IntRates(:, end), 1, 9);
AvgIRDiff        = mean(IntRateDiff);

% Display second PC along interest rate differentials
AuxMat           = [Weights(2, :)' AvgIRDiff'];
Correl           = corrcoef(AuxMat);

% Note that the sign of the correlation is immaterial here, only absolute
% value matters