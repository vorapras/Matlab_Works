%% London School of Economics
%  FM442 - Michaelmas Term 2016
%  Author: Domingos Romualdo
%  Date: 13 October 2016
%  Assignment 2 - Suggested Solution

%% *** Clearing workspace *** 
clc;
clear all;
clf;

%% *** Part 1: Loading daily data ***
Data_D = importdata('FM442 - HW2 - Daily Data.csv');

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

%% *** Part 4: Data Work with All Securities ***

% Computing statistics for series
% Note the conversion of mean and standard deviation into annual units
StockStats            = NaN(6, NStocks + 1);
StockStats(1, :)      = 252 * mean(LogReturns);
StockStats(2, :)      = sqrt(252) * std(LogReturns);
StockStats(3, :)      = skewness(LogReturns);
StockStats(4, :)      = kurtosis(LogReturns);

% Loop to compute Jarque-Bera statistics
for i=1:(NStocks+1)
   [H, P, JBSTAT]     = jbtest(LogReturns(:, i));
   StockStats(5, i)   = JBSTAT;
   StockStats(6, i)   = P;
end

% Computing covariance and correlation matrices of returns
CovRet                = cov(LogReturns);
CorrRet               = corrcoef(LogReturns);

% Computing covariance and correlation matrices of squared returns
LogReturnsSq          = LogReturns .* LogReturns;
CovRetSq              = cov(LogReturnsSq);
CorrRetSq             = corrcoef(LogReturnsSq);

%% *** Part 5: Miscellaneous Charts ***

% Choosing stocks to work with (MSFT, JPM, INTC, SPX)
Universe              = [1  4  5  7];
Tickers               = ['MSFT'; 'JPM '; 'INTC'; 'SPX '];
LogReturns2           = LogReturns(:, Universe);

% Creating histogram for each of the series
FigNo                 = 1;
for i=1:4
    figure(FigNo);
    histogram(LogReturns2(:, i));
    TitleStr          = ['Histogram of daily returns for ' Tickers(i, :)];
    title(TitleStr);
    set(gca, 'XTickMode', 'manual');
    set(gca, 'XTickLabel', num2str(100 .* get(gca, 'XTick')', '%g%%'))    
    FigNo             = FigNo + 1;
end

% Alternatively, producing single chart with all histograms
for i=1:4
    subplot(2, 2, i);
    histogram(LogReturns2(:, i));
    TitleStr          = ['Histogram of daily returns for ' Tickers(i, :)];
    title(TitleStr);
    ylabel('Historical Sample');
    set(gca, 'XTickMode', 'manual');
    set(gca, 'XTickLabel', num2str(100 .* get(gca, 'XTick')', '%g%%'))    
end
FigNo                 = FigNo + 1;

% Creating a QQ-plot of stock returns against standard normal distribution
for i=1:4
    subplot(2, 2, i);
    qqplot(LogReturns2(:, i));
    TitleStr          = ['QQ-Plot for ' Tickers(i, :)];
    title(TitleStr);
    ylabel('Historical Sample');
end
FigNo                 = FigNo + 1;

% Creating a QQ-plot of stock returns against t-4 distribution - must first
% create the distribution to plot against returns
tdist                 = makedist('tLocationScale', 'mu', 0, ...
                           'sigma',1, 'nu', 4);
for i=1:4
    subplot(2, 2, i);
    qqplot(LogReturns2(:, i), tdist);
    TitleStr          = ['QQ-Plot for ' Tickers(i, :)];
    title(TitleStr);
    ylabel('Historical Sample');
    xlabel('Quantiles of t(4) distribution');
end
FigNo                 = FigNo + 1;

% Plotting the autocorrelation function for each of the series
for i=1:4
    subplot(2, 2, i);
    autocorr(LogReturns2(:, i));
    TitleStr          = ['Sample autocorrelation for ' Tickers(i, :)];
    title(TitleStr);
    ylabel('Historical Sample');
end
FigNo                 = FigNo + 1;

% Plotting the autocorrelation function for each of the series squared
LogReturns2Sq         = LogReturns2 .* LogReturns2;
for i=1:4
    subplot(2, 2, i);
    autocorr(LogReturns2Sq(:, i));
    TitleStr          = ['Sample autocorrelation for ' Tickers(i, :) ...
                      ' squared'];
    title(TitleStr);
    ylabel('Historical Sample');
end
FigNo                 = FigNo + 1;

%% *** Part 6: Computing moving average variance ***
SPXReturns            = LogReturns(:, 7);
WindowLen             = 50;
MovAvStd              = NaN(NDates - WindowLen, 1);

for i=1:(NDates - WindowLen)
    MovAvStd(i, 1)    = std(SPXReturns(i:(i + WindowLen - 1), :));
end

% Putting together information for chart
ChartData             = [SPXReturns((WindowLen+1):end, :) -2*MovAvStd ...
                         2*MovAvStd];
StartDate             = 1989;
EndDate               = 2015; 
TimeLabels            = linspace(StartDate, EndDate, NDates - WindowLen);

% Charting the series
figure(FigNo)
titleStr = 'S&P 500 Returns and 50-Day Moving Average +/- 2 Std. Dev.';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate -0.15 0.15]);
grid on;
hold on;
plot(TimeLabels, ChartData(:, 1), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, ChartData(:, 2), 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
plot(TimeLabels, ChartData(:, 3), 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');             
set(gcf, 'color', 'white');                  

%% *** Part 7: Portfolio Computations ***
% Note: when computing portfolio returns, must work
% with original returns, not log ones
StockReturns           = TR_Stocks(:, Universe(:, 1:3));
AuxVec                 = ones(3, 1)/3;
PortfolioRet           = StockReturns * AuxVec;
PortfolioLogRet        = log(1 + PortfolioRet);

% Computing portfolio statistics
PortfolioStats         = NaN(6, 1);
PortfolioStats(1, :)   = 252 * mean(PortfolioLogRet);
PortfolioStats(2, :)   = sqrt(252) * std(PortfolioLogRet);
PortfolioStats(3, :)   = skewness(PortfolioLogRet);
PortfolioStats(4, :)   = kurtosis(PortfolioLogRet);
[H, P, JBSTAT]         = jbtest(PortfolioLogRet);
PortfolioStats(5, 1)   = JBSTAT;
PortfolioStats(6, 1)   = P;